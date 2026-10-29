"""Full actin torque budget by MAGNITUDE, so am_head can be compared to the
other sources in like units. The 2-9% figure on record is a PROJECTED restoring
component (tau . e_align); magnitudes and projections differ by ~50x here, so
they must not be mixed.

Also reports the projection and the alignment efficiency
(projected / magnitude) per source -- a source can be large but useless if its
torque points the wrong way.
"""
import os, glob
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
import numpy as np, h5py
P = "/scratch/users/xsli/sarcomere/barbed_capsule_normal_full_20261002/data/*kaa3_seed*_5000ms_capsule_normal.h5"
rows = {}
for p in sorted(glob.glob(P)):
    with h5py.File(p, "r", locking=False) as f:
        cf = int(np.asarray(f["meta/committed_frames"]).reshape(-1)[-1])
        cf = min(cf, f["actin/direction"].shape[0], f["state/current_step"].shape[0])
        st = np.asarray(f["state/current_step"][:cf]).reshape(-1)
        k = int(np.searchsorted(st, 5000*100, "right")) - 1
        if float(st[k]) < 5000*100*0.999: continue
        srcs = sorted({n[3:-2] for n in f["actin"].keys() if n.startswith("tb_") and n.endswith("_x")})
        u = np.asarray(f["actin/direction"][k], float); u /= np.linalg.norm(u, axis=1, keepdims=True)
        # aligning direction: rotate u toward the global director
        Q = 1.5*np.einsum("ia,ib->ab", u, u)/u.shape[0] - 0.5*np.eye(3)
        d = np.linalg.eigh(Q)[1][:, -1]
        dot = np.clip(np.abs(u @ d), 0, 1)
        sgn = np.sign(u @ d); sgn[sgn == 0] = 1.0
        dd = d[None, :]*sgn[:, None]
        perp = dd - (np.einsum("ik,ik->i", dd, u))[:, None]*u
        nrm = np.linalg.norm(perp, axis=1, keepdims=True)
        e_align = np.divide(perp, np.where(nrm == 0, 1, nrm))
        for s in srcs:
            try:
                # tb_* datasets are (frames, nActin, 1); ravel the trailing axis
                t = np.stack([np.asarray(f["actin/tb_%s_%s" % (s, c)][k], float).ravel()
                              for c in "xyz"], -1)
                assert t.shape == (u.shape[0], 3), (s, t.shape)
            except KeyError:
                continue
            mag = np.linalg.norm(t, axis=1)
            proj = np.einsum("ik,ik->i", t, e_align)
            nz = mag > 0
            if nz.sum() == 0: continue
            rows.setdefault(s, []).append((float(mag[nz].mean()), float(proj[nz].mean()),
                                           float(nz.mean())))
print("actin torque budget at 5000 ms, baseline arm A, n=%d seeds\n" % len(next(iter(rows.values()))))
print("  %-18s %12s %12s %10s %9s" % ("source", "|tau| mean", "proj mean", "proj/|tau|", "active%"))
print("  " + "-"*66)
out = []
for s, v in rows.items():
    mag = float(np.mean([x[0] for x in v])); pr = float(np.mean([x[1] for x in v]))
    act = float(np.mean([x[2] for x in v]))
    out.append((mag, s, pr, act))
for mag, s, pr, act in sorted(out, reverse=True):
    print("  %-18s %12.4f %12.4f %9.1f%% %8.1f%%" % (s, mag, pr, 100*pr/mag if mag else 0, 100*act))
tot = {s: m for m, s, p, a in out}
if "am_head" in tot and "actn2_align" in tot:
    print("\n  am_head |tau| is %.1f%% of actn2_align |tau|  (the 2-9%% on record is the PROJECTION)"
          % (100*tot["am_head"]/tot["actn2_align"]))
