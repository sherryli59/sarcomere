"""Does spokes=1 actually improve actin-myosin ALIGNMENT? Paired, per seed.

theta_AM = arccos|u_actin . w_myosin| per distinct (actin, myosin) pair carrying
>=1 STRONG head -- the same definition as analysis/am_angles.py. Reported as the
per-seed MEDIAN (the target is a median < 10 deg), paired A vs B, at matched
time. Isotropic null for reference: mean 57.3 deg, median 60 deg.
"""
import os, glob, collections
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
import numpy as np, h5py
D = "/scratch/users/xsli/sarcomere"
A_PAT = f"{D}/barbed_capsule_normal_full_20261002/data/*kaa3_seed*_5000ms_capsule_normal.h5"
B_PAT = f"{D}/hph_ladder_20261005/data/*kaa3_seed*_5000ms_spokes1.h5"
NFR = 5

def med_angle(p, t_ms):
    with h5py.File(p, "r", locking=False) as f:
        cf = int(np.asarray(f["meta/committed_frames"]).reshape(-1)[-1])
        cf = min(cf, f["actin/direction"].shape[0], f["state/current_step"].shape[0],
                 f["state/am_head_state_current"].shape[0])
        st = np.asarray(f["state/current_step"][:cf]).reshape(-1)
        k = int(np.searchsorted(st, t_ms*100, "right")) - 1
        if k < NFR or float(st[k]) < t_ms*100*0.999: return None
        fr = list(range(k-NFR+1, k+1))
        u = np.asarray(f["actin/direction"][fr], float); u /= np.linalg.norm(u, axis=2, keepdims=True)
        w = np.asarray(f["myosin/direction"][fr], float); w /= np.linalg.norm(w, axis=2, keepdims=True)
        hs = np.asarray(f["state/am_head_state_current"][fr])
        ha = np.asarray(f["state/am_head_bound_actin_current"][fr]).astype(int)
        th = []
        for t in range(len(fr)):
            mi, hh = np.nonzero((hs[t] == 2) & (ha[t] >= 0))
            if mi.size == 0: continue
            pairs = np.unique(np.stack([ha[t][mi, hh], mi], 1), axis=0)
            d = np.abs(np.einsum("pk,pk->p", u[t][pairs[:, 0]], w[t][pairs[:, 1]]))
            th.append(np.degrees(np.arccos(np.clip(d, 0, 1))))
        if not th: return None
        v = np.concatenate(th)
        return dict(med=float(np.median(v)), mean=float(v.mean()),
                    p10=float(np.percentile(v, 10)), frac_lt10=float((v < 10).mean()), n=v.size)

for T in (1500, 5000):
    A = {p.split("seed")[1][:8]: med_angle(p, T) for p in sorted(glob.glob(A_PAT))}
    B = {p.split("seed")[1][:8]: med_angle(p, T) for p in sorted(glob.glob(B_PAT))}
    seeds = sorted(s for s in set(A) & set(B) if A[s] and B[s])
    print("\n######## matched t = %d ms   (paired seeds: %s)" % (T, ",".join(seeds) or "none"))
    if not seeds: print("   no seed pair has reached this time"); continue
    print("   seed      A median   B median    delta  |  A frac<10deg  B frac<10deg")
    for s in seeds:
        print("   %-9s %8.2f  %8.2f  %+7.2f  | %11.1f%% %12.1f%%"
              % (s, A[s]["med"], B[s]["med"], B[s]["med"]-A[s]["med"],
                 100*A[s]["frac_lt10"], 100*B[s]["frac_lt10"]))
    d = np.array([B[s]["med"]-A[s]["med"] for s in seeds])
    se = d.std(ddof=1)/np.sqrt(d.size) if d.size > 1 else float("nan")
    print("   MEAN A %.2f   MEAN B %.2f   paired delta %+.2f (SE %.2f, t=%.2f, n=%d)"
          % (np.mean([A[s]["med"] for s in seeds]), np.mean([B[s]["med"] for s in seeds]),
             d.mean(), se, d.mean()/se if se and se > 0 else float("nan"), d.size))
    print("   TARGET median < 10 deg: A %s   B %s"
          % ("MET" if np.mean([A[s]["med"] for s in seeds]) < 10 else "NOT MET",
             "MET" if np.mean([B[s]["med"] for s in seeds]) < 10 else "NOT MET"))
print("\n   isotropic null for scale: mean 57.3 deg, median 60.0 deg")
