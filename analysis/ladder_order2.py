"""M order, Z order and S for the ladder arms, AVERAGED OVER A WINDOW.

The single-frame version (job 46663086) is not trustworthy: band_order is a
global coherence measure whose frame-to-frame scatter is comparable to the
between-arm differences. Average over the 1000-1500 ms window (~100 frames at
save_every=500 -> 5 ms/frame) and report BOTH the across-seed sd and the
within-seed temporal sd, so the two noise sources are separable.
"""
import os, glob
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
import numpy as np, h5py
T0_MS, T1_MS = 1000, 1500

def band_order(x, box_x):
    wrapped = np.asarray(x, dtype=float) % box_x
    harmonics = np.arange(1, 5, dtype=float)
    phase = np.exp(2j*np.pi*wrapped[..., None]*harmonics[None, None, :]/box_x).mean(axis=1)
    return np.max(np.abs(phase), axis=1)

FULL = "/scratch/users/xsli/sarcomere/barbed_capsule_normal_full_20261002/data"
LAD  = "/scratch/users/xsli/sarcomere/hph_ladder_20261005/data"
ARMS = [("A", f"{FULL}/*kaa3_seed*_5000ms_capsule_normal.h5"),
        ("B", f"{LAD}/*kaa3_seed*_5000ms_spokes1.h5"),
        ("C", f"{LAD}/*kaa3_seed*_5000ms_hph27.h5")]
res = {}
for tag, pat in ARMS:
    rows = {}
    for p in sorted(glob.glob(pat)):
        seed = p.split("seed")[1][:8]
        try:
            with h5py.File(p, "r", locking=False) as f:
                a = f["meta"].attrs
                box = float(a["box_x"]); La = float(a["actin_length"])
                cf = int(np.asarray(f["meta/committed_frames"]).reshape(-1)[-1])
                cf = min(cf, f["actin/center"].shape[0], f["state/current_step"].shape[0])
                st = np.asarray(f["state/current_step"][:cf]).reshape(-1)
                i0 = int(np.searchsorted(st, T0_MS*100, side="left"))
                i1 = int(np.searchsorted(st, T1_MS*100, side="right"))
                i1 = min(i1, cf)
                if i1 - i0 < 20: print("  SKIP %s (%d frames in window)" % (seed, i1-i0)); continue
                mc = np.asarray(f["myosin/center"][i0:i1], float)
                ac = np.asarray(f["actin/center"][i0:i1], float)
                ad = np.asarray(f["actin/direction"][i0:i1], float)
                u = ad/np.linalg.norm(ad, axis=2, keepdims=True)
                M = band_order(mc[:, :, 0], box)
                Z = np.maximum(band_order(ac[:, :, 0] - 0.5*La*u[:, :, 0], box),
                               band_order(ac[:, :, 0] + 0.5*La*u[:, :, 0], box))
                Q = 1.5*np.einsum("fia,fib->fab", u, u)/u.shape[1] - 0.5*np.eye(3)[None]
                S = np.linalg.eigvalsh(Q)[:, -1]
                rows[seed] = dict(M=M.mean(), Z=Z.mean(), S=S.mean(), n=i1-i0,
                                  Msd=M.std(ddof=1), Zsd=Z.std(ddof=1), Ssd=S.std(ddof=1))
        except Exception as e:
            print("  ERR %s: %s" % (seed, e))
    res[tag] = rows
    if not rows: continue
    print("\n=== arm %s   n=%d seeds, %d frames each (%d-%d ms)"
          % (tag, len(rows), list(rows.values())[0]["n"], T0_MS, T1_MS))
    for s, r in sorted(rows.items()):
        print("    %-8s M=%.4f(+-%.4f t)  Z=%.4f(+-%.4f t)  S=%.4f(+-%.4f t)"
              % (s, r["M"], r["Msd"], r["Z"], r["Zsd"], r["S"], r["Ssd"]))
    for key in ("M", "Z", "S"):
        v = np.array([r[key] for r in rows.values()])
        tsd = np.mean([r[key+"sd"] for r in rows.values()])
        print("    MEAN %s = %.4f   across-seed sd %.4f   mean within-seed temporal sd %.4f"
              % (key, v.mean(), v.std(ddof=1), tsd))
print("\n=== paired differences over the window")
for t1, t2 in (("A", "B"), ("B", "C")):
    ks = sorted(set(res[t1]) & set(res[t2]))
    for key, nm in (("M", "M order"), ("Z", "Z order"), ("S", "nematic S")):
        d = np.array([res[t2][s][key] - res[t1][s][key] for s in ks])
        se = d.std(ddof=1)/np.sqrt(len(d))
        print("  %s-%s  %-10s %+.4f  (sd %.4f, SE %.4f, t=%6.2f)  improved %d/%d"
              % (t2, t1, nm, d.mean(), d.std(ddof=1), se,
                 d.mean()/se if se else float("inf"), (d > 0).sum(), len(d)))
print("\nfloors: M (N=320) 0.078 mean / 0.136 p99 ; Z (N=1600) 0.035 mean / 0.061 p99")
