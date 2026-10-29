"""S(t) and ACTN2 link count at matched times for every full_20261002 run,
plus the k_soff=27 arm, so the seed-to-seed scatter can be compared against
the between-arm differences. S = largest eigenvalue of the nematic Q."""
import os, glob
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
import numpy as np, h5py
TS = [1000, 2000, 3000, 4200, 5000]

def course(p):
    with h5py.File(p, "r", locking=False) as f:
        cf = int(np.asarray(f["meta/committed_frames"]).reshape(-1)[-1])
        cf = min(cf, f["actin/direction"].shape[0], f["state/current_step"].shape[0])
        st = np.asarray(f["state/current_step"][:cf]).reshape(-1)
        res = {}
        for tms in TS:
            k = int(np.searchsorted(st, tms*100, side="right")) - 1
            if k < 0 or k >= cf:
                res[tms] = (np.nan, np.nan); continue
            ad = np.asarray(f["actin/direction"][k:k+1], float)
            u = ad/np.linalg.norm(ad, axis=2, keepdims=True)
            Q = 1.5*np.einsum("ia,ib->ab", u[0], u[0])/u.shape[1] - 0.5*np.eye(3)
            S = float(np.linalg.eigvalsh(Q).max())
            lk = float(np.asarray(f["state/aa_pairs_current_count"][k]).reshape(-1)[-1])
            res[tms] = (S, lk)
        return res, cf

groups = [("kaa3 capsule(k_soff=160)", "/scratch/users/xsli/sarcomere/barbed_capsule_normal_full_20261002/data/*kaa3_seed*_5000ms_capsule_normal.h5"),
          ("kaa0 capsule(k_soff=160)", "/scratch/users/xsli/sarcomere/barbed_capsule_normal_full_20261002/data/*kaa0_seed*_5000ms_capsule_normal.h5"),
          ("kaa3 k_soff=27",            "/scratch/users/xsli/sarcomere/ksoff27_prod_20261003/data/*kaa3_seed*_5000ms_ksoff27.h5")]
for tag, pat in groups:
    print("\n=== %s" % tag)
    print("  %-10s %6s | %s" % ("seed", "frames", " | ".join("S@%-5d links" % t for t in TS)))
    Sat = {t: [] for t in TS}
    for p in sorted(glob.glob(pat)):
        seed = p.split("seed")[1][:8]
        res, cf = course(p)
        row = " | ".join(("  nan    nan " if np.isnan(res[t][0]) else "%.4f %6.0f" % res[t]) for t in TS)
        print("  %-10s %6d | %s" % (seed, cf, row))
        for t in TS:
            if not np.isnan(res[t][0]): Sat[t].append(res[t][0])
    for t in TS:
        v = Sat[t]
        if len(v) >= 2:
            print("    S@%-5d n=%d  mean=%.4f  min=%.4f  max=%.4f  spread=%.4f  sd=%.4f"
                  % (t, len(v), np.mean(v), min(v), max(v), max(v)-min(v), np.std(v, ddof=1)))
