"""AM angle distribution for pairs carrying at least one STRONG head.

theta_AM = arccos(|u_actin . w_myosin|), ONE entry per distinct (actin, myosin)
pair that has >=1 strong head.  0 deg = actin parallel to the thick-filament
axis.  Polarity-blind (absolute value), so the range is [0, 90].  This is the
same definition as angles_extract.py part (B), 2026-10-02.

Three arms at matched time.  TARGET=2000 ms is inside the committed window of
every run including the truncated k_soff=27 seeds (shortest is 508 frames
~ 2540 ms), so no arm contributes a repeated last frame.
"""
import os, glob, json, collections
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
import numpy as np, h5py

FULL = "/scratch/users/xsli/sarcomere/barbed_capsule_normal_full_20261002/data"
K27  = "/scratch/users/xsli/sarcomere/ksoff27_prod_20261003/data"
ARMS = [("kaa3_ksoff160", f"{FULL}/*kaa3_seed*_5000ms_capsule_normal.h5"),
        ("kaa0_ksoff160", f"{FULL}/*kaa0_seed*_5000ms_capsule_normal.h5"),
        ("kaa3_ksoff27",  f"{K27}/*kaa3_seed*_5000ms_ksoff27.h5")]
TARGET, NFR, STRONG = 2000, 5, 2
BINS = np.arange(0, 91, 1.0)
out = {"meta": {"target_ms": TARGET, "nframes_per_seed": NFR,
                "bins": BINS.tolist(), "definition":
                "theta_AM = arccos(|u_actin . w_myosin|) per distinct (actin,myosin) pair with >=1 strong head; polarity-blind; deg"},
       "hist": {}, "stats": {}}

for arm, pat in ARMS:
    vals = []; hpp = []; seeds = []; frames_used = 0
    for p in sorted(glob.glob(pat)):
        with h5py.File(p, "r", locking=False) as f:
            cf = int(np.asarray(f["meta/committed_frames"]).reshape(-1)[-1])
            cf = min(cf, f["actin/direction"].shape[0], f["state/current_step"].shape[0])
            st = np.asarray(f["state/current_step"][:cf]).reshape(-1)
            k = int(np.searchsorted(st, TARGET*100, side="right")) - 1
            if k < NFR:
                print("  SKIP %s (only %d frames)" % (os.path.basename(p), cf)); continue
            seeds.append(p.split("seed")[1][:8])
            fr = list(range(k-NFR+1, k+1)); frames_used += len(fr)
            ad = np.asarray(f["actin/direction"][fr], float)
            u = ad/np.linalg.norm(ad, axis=2, keepdims=True)
            md = np.asarray(f["myosin/direction"][fr], float)
            w = md/np.linalg.norm(md, axis=2, keepdims=True)
            hs = np.asarray(f["state/am_head_state_current"][fr])
            ha = np.asarray(f["state/am_head_bound_actin_current"][fr])
            for t in range(len(fr)):
                strong = (hs[t] == STRONG) & (ha[t] >= 0)
                mi, hi = np.nonzero(strong)
                if mi.size == 0: continue
                ai = ha[t][mi, hi]
                stack = np.stack([ai, mi], 1)
                pairs = np.unique(stack, axis=0)
                hpp.extend(collections.Counter(map(tuple, stack.tolist())).values())
                dot = np.abs(np.einsum("pk,pk->p", u[t][pairs[:, 0]], w[t][pairs[:, 1]]))
                vals.append(np.degrees(np.arccos(np.clip(dot, 0, 1))))
    v = np.concatenate(vals) if vals else np.array([])
    h = np.histogram(v, bins=BINS)[0]
    # nothing may be clipped: the histogram must account for every sample
    assert int(h.sum()) == v.size, "clipped %d of %d in %s" % (v.size-int(h.sum()), v.size, arm)
    out["hist"][arm] = h.tolist()
    out["stats"][arm] = dict(
        n_pairs=int(v.size), seeds=seeds, frames_used=frames_used,
        mean=float(v.mean()), median=float(np.median(v)),
        p10=float(np.percentile(v, 10)), p90=float(np.percentile(v, 90)),
        sd=float(v.std(ddof=1)), heads_per_pair=float(np.mean(hpp)))
    s = out["stats"][arm]
    print("%-16s n=%7d seeds=%s  mean=%.2f  median=%.2f  p10=%.2f  p90=%.2f  sd=%.2f  heads/pair=%.3f"
          % (arm, s["n_pairs"], ",".join(s["seeds"]), s["mean"], s["median"],
             s["p10"], s["p90"], s["sd"], s["heads_per_pair"]))
    print("                 histogram sum=%d == n_pairs (nothing clipped)" % int(h.sum()))

dest = "/scratch/users/xsli/sarcomere/ksoff27_prod_20261003/analysis/am_angles.json"
with open(dest, "w") as fh: json.dump(out, fh)
print("wrote", dest)
