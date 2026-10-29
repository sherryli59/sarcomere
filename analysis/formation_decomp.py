"""Decompose ACTN2 formation rate = (eligible candidate pairs) x (per-pair rate).

Eligibility, from the run's own knobs:
  - barbed tips within aa_cutoff = aa_preferred_centerline + aa_capture_margin
    (0.079 + 0.0396 = 0.1186 um). Barbed end is the MINUS tip.
  - antiparallel: dot(u_i, u_j) < aa_association_max_direction_dot = -0.5
  - both actins below the max_aa_bonds_per_actin = 4 cap

If eligible pairs fall ~37% in arm C, the cause is candidate availability
(geometry / slot jamming). If eligible pairs are flat and the per-pair rate
falls, the cause is a per-attempt gate instead.
"""
import os, glob, collections
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
import numpy as np, h5py
T0, T1 = 1000, 1500
FULL = "/scratch/users/xsli/sarcomere/barbed_capsule_normal_full_20261002/data"
LAD  = "/scratch/users/xsli/sarcomere/hph_ladder_20261005/data"
ARMS = [("A hph9/spokes3 ", f"{FULL}/*kaa3_seed*_5000ms_capsule_normal.h5"),
        ("B hph9/spokes1 ", f"{LAD}/*kaa3_seed*_5000ms_spokes1.h5"),
        ("C hph27/spokes3", f"{LAD}/*kaa3_seed*_5000ms_hph27.h5")]
res = {}
for tag, pat in ARMS:
    rows = []
    for p in sorted(glob.glob(pat)):
        with h5py.File(p, "r", locking=False) as f:
            a = f["meta"].attrs
            cut = float(a["aa_preferred_centerline"]) + float(a["aa_capture_margin"])
            maxdot = float(a["aa_association_max_direction_dot"])
            maxb = int(a["max_aa_bonds_per_actin"]); L = float(a["actin_length"])
            box = np.array([a["box_x"], a["box_y"], a["box_z"]], float)
            dt = float(a["dt"])
            cf = int(np.asarray(f["meta/committed_frames"]).reshape(-1)[-1])
            cf = min(cf, f["actin/center"].shape[0], f["state/current_step"].shape[0])
            st = np.asarray(f["state/current_step"][:cf]).reshape(-1)
            i0 = int(np.searchsorted(st, T0*100, "left")); i1 = min(int(np.searchsorted(st, T1*100, "right")), cf)
            fr = list(range(i0, i1, 25))
            ac = np.asarray(f["actin/center"][fr], float)
            ad = np.asarray(f["actin/direction"][fr], float)
            u = ad/np.linalg.norm(ad, axis=2, keepdims=True)
            barbed = ac - 0.5*L*u                     # MINUS tip
            # bonds per actin from the live link identity
            ident = np.asarray(f["catch_bond/local_link_identity"])
            off = np.asarray(f["catch_bond/local_link_frame_offsets"])
            elig = []; elig_slot = []; sat = []
            for t, fi in enumerate(fr):
                d = barbed[t][:, None, :] - barbed[t][None, :, :]
                d -= box*np.rint(d/box)
                dist = np.linalg.norm(d, axis=2)
                dot = np.einsum("ik,jk->ij", u[t], u[t])
                np.fill_diagonal(dist, 1e9)
                ok = (dist < cut) & (dot < maxdot)
                elig.append(int(ok.sum())//2)
                # bond counts at this frame
                j = int(np.searchsorted(off[:, 0], st[fi], "left"))
                cnt = np.zeros(u.shape[1], int)
                if j < off.shape[0] and int(off[j, 0]) == int(st[fi]):
                    s0 = int(off[j, 1]); n = int(off[j, 2])
                    if n > 0:
                        ii = ident[s0:s0+n, 2].astype(int); jj = ident[s0:s0+n, 3].astype(int)
                        for x in np.concatenate([ii, jj]):
                            if 0 <= x < cnt.size: cnt[x] += 1
                sat.append(float((cnt >= maxb).mean()))
                free = cnt < maxb
                ok2 = ok & free[:, None] & free[None, :]
                elig_slot.append(int(ok2.sum())//2)
            rows.append(dict(elig=float(np.mean(elig)), elig_slot=float(np.mean(elig_slot)),
                             sat=float(np.mean(sat)), cut=cut))
    g = lambda k: float(np.mean([r[k] for r in rows]))
    res[tag[0]] = dict(elig=g("elig"), elig_slot=g("elig_slot"), sat=g("sat"))
    print("\n=== %s  (n=%d seeds, capture cutoff %.4f um, antiparallel dot < %.1f)"
          % (tag, len(rows), rows[0]["cut"], -0.5))
    print("    eligible pairs (geometry only)            %9.1f" % g("elig"))
    print("    eligible pairs (geometry AND free slots)  %9.1f" % g("elig_slot"))
    print("    actins at the 4-bond cap                  %9.1f%%" % (100*g("sat")))
FORM = {"A": 94643.0, "B": 88738.5, "C": 59943.0}   # measured /s, job 46694162
print("\n=== decomposition  (formation rates from job 46694162)")
print("  arm   formation/s   eligible   per-pair rate /s   elig(+slots)   per-pair(+slots) /s")
for k in ("A", "B", "C"):
    r = res[k]
    print("   %s     %9.1f   %8.1f   %14.2f   %12.1f   %18.2f"
          % (k, FORM[k], r["elig"], FORM[k]/r["elig"], r["elig_slot"], FORM[k]/max(1e-9, r["elig_slot"])))
if "A" in res and "C" in res:
    A, C = res["A"], res["C"]
    print("\n  C/A ratios:")
    print("    formation rate                  %.3f   <- the thing to explain" % (FORM["C"]/FORM["A"]))
    print("    eligible pairs (geometry)       %.3f" % (C["elig"]/A["elig"]))
    print("    eligible pairs (+ free slots)   %.3f" % (C["elig_slot"]/A["elig_slot"]))
    print("    per-pair rate (geometry)        %.3f" % ((FORM["C"]/C["elig"])/(FORM["A"]/A["elig"])))
    print("    per-pair rate (+ free slots)    %.3f" % ((FORM["C"]/C["elig_slot"])/(FORM["A"]/A["elig_slot"])))
