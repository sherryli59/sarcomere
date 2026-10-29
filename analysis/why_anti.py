"""Why is am_head anti-aligning?

MECHANISM: a single off-centre head pull TILTS the actin. The head grips at a
material coordinate c along the actin and pulls it toward the myosin; the torque
about the actin centre is r x F with r = c*u, so one grip swings that END toward
the thick filament, i.e. ACROSS the myosin axis -> angle increases. Only two
grips at different axial stations press the actin flat (the span couple).

PREDICTION: tb_am_head's projection onto the aligning direction should be
NEGATIVE for actins gripped by ONE crown and POSITIVE for actins gripped by two
or more. Most pairs are single-crown (86-93%), so the population mean comes out
negative.

Also stratify by |c| (how off-centre the grip is): the tilt torque scales with
the moment arm, so anti-alignment should grow with |c| for single-crown actins.
"""
import os, glob, collections
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
import numpy as np, h5py
BARE, BEAR = 0.10, 0.5*(1.0-0.10)
P = "/scratch/users/xsli/sarcomere/barbed_capsule_normal_full_20261002/data/*kaa3_seed*_5000ms_capsule_normal.h5"
def anchor_axial_table(hph, Lm, spokes):
    spokes = max(1, int(spokes)); cc = max(1, (hph+spokes-1)//spokes)
    out = np.empty(2*hph)
    for hi in (0, 1):
        sgn = -1.0 if hi == 0 else 1.0
        for h in range(hph):
            out[hi*hph+h] = sgn*(0.5*BARE*Lm + ((h//spokes + 0.5)/cc)*BEAR*Lm)
    return out
by_crown = collections.defaultdict(list); by_c = collections.defaultdict(list)
for p in sorted(glob.glob(P)):
    with h5py.File(p, "r", locking=False) as f:
        a = f["meta"].attrs
        hph = int(a["am_heads_per_half"]); spk = int(float(a.get("am_spokes_per_crown", 3)))
        anc = anchor_axial_table(hph, float(a["myosin_length"]), spk)
        cf = int(np.asarray(f["meta/committed_frames"]).reshape(-1)[-1])
        cf = min(cf, f["actin/direction"].shape[0], f["state/current_step"].shape[0])
        st = np.asarray(f["state/current_step"][:cf]).reshape(-1)
        ks = [int(np.searchsorted(st, t*100, "right"))-1 for t in (3000, 4000, 5000)]
        ks = [k for k in ks if float(st[k]) > 0]
        for k in ks:
            u = np.asarray(f["actin/direction"][k], float); u /= np.linalg.norm(u, axis=1, keepdims=True)
            Q = 1.5*np.einsum("ia,ib->ab", u, u)/u.shape[0] - 0.5*np.eye(3)
            d = np.linalg.eigh(Q)[1][:, -1]
            sgn = np.sign(u @ d); sgn[sgn == 0] = 1.0
            dd = d[None, :]*sgn[:, None]
            perp = dd - (np.einsum("ik,ik->i", dd, u))[:, None]*u
            nrm = np.linalg.norm(perp, axis=1, keepdims=True)
            e = np.divide(perp, np.where(nrm == 0, 1, nrm))
            t = np.stack([np.asarray(f["actin/tb_am_head_%s" % c][k], float).ravel() for c in "xyz"], -1)
            proj = np.einsum("ik,ik->i", t, e)
            hs = np.asarray(f["state/am_head_state_current"][k])
            ha = np.asarray(f["state/am_head_bound_actin_current"][k]).astype(int)
            tg = np.asarray(f["state/am_head_target_axial_current"][k], float)
            crowns = collections.defaultdict(set); coords = collections.defaultdict(list)
            mi, hh = np.nonzero((hs >= 1) & (ha >= 0))
            for m_, h_ in zip(mi.tolist(), hh.tolist()):
                ai = int(ha[m_][h_])
                crowns[ai].add(round(float(anc[h_]), 9)); coords[ai].append(abs(float(tg[m_][h_])))
            for ai, cr in crowns.items():
                if proj[ai] == 0: continue
                by_crown[min(len(cr), 3)].append(proj[ai])
                if len(cr) == 1:
                    by_c[min(4, int(np.mean(coords[ai])/0.125))].append(proj[ai])
print("am_head projection onto the aligning direction, by DISTINCT CROWNS gripping the actin")
print("  (positive = aligning, negative = anti-aligning)\n")
print("  crowns      n       mean proj     median")
for kk in sorted(by_crown):
    v = np.array(by_crown[kk])
    print("  %-6d %8d   %+10.4f   %+10.4f" % (kk, v.size, v.mean(), np.median(v)))
print("\nsingle-crown actins only, by mean |grip offset| from the actin centre (um)")
print("  the tilt torque scales with the moment arm, so this should get MORE negative\n")
print("  |c| bin        n       mean proj")
for kk in sorted(by_c):
    v = np.array(by_c[kk])
    print("  %.3f-%.3f %8d   %+10.4f" % (kk*0.125, (kk+1)*0.125, v.size, v.mean()))
