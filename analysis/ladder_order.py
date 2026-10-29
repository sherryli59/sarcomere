"""M order and Z order for the three ladder arms at matched time.

band_order is VISUALIZATION_PIPELINE.md verbatim. M order = band_order over
myosin centre x. Z order = max of band_order over the TWO actin-tip
populations (centre -/+ La/2 * u_a). Both are POSITIONAL banding measures, not
the nematic S.

All three arms have identical N (320 myosins, 1600 actins), so the N-dependent
noise floor is the same for all and they ARE comparable. Floors from a 20000-
trial Monte Carlo of uniform random positions: N=320 -> mean 0.078 (p99 0.136);
N=1600 -> mean 0.035 (p99 0.061).
"""
import os, glob
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
import numpy as np, h5py
TARGET_MS = 1500

def band_order(x, box_x):
    wrapped = np.asarray(x, dtype=float) % box_x
    harmonics = np.arange(1, 5, dtype=float)
    phase = np.exp(2j*np.pi*wrapped[..., None]*harmonics[None, None, :]/box_x).mean(axis=1)
    return np.max(np.abs(phase), axis=1)

FULL = "/scratch/users/xsli/sarcomere/barbed_capsule_normal_full_20261002/data"
LAD  = "/scratch/users/xsli/sarcomere/hph_ladder_20261005/data"
ARMS = [("A baseline hph9/spokes3", f"{FULL}/*kaa3_seed*_5000ms_capsule_normal.h5"),
        ("B spokes1  hph9/spokes1", f"{LAD}/*kaa3_seed*_5000ms_spokes1.h5"),
        ("C hph27    hph27/spokes3", f"{LAD}/*kaa3_seed*_5000ms_hph27.h5")]
res = {}
for tag, pat in ARMS:
    rows = []
    for p in sorted(glob.glob(pat)):
        seed = p.split("seed")[1][:8]
        try:
            with h5py.File(p, "r", locking=False) as f:
                a = f["meta"].attrs
                box = float(a["box_x"]); La = float(a["actin_length"])
                cf = int(np.asarray(f["meta/committed_frames"]).reshape(-1)[-1])
                cf = min(cf, f["actin/center"].shape[0], f["state/current_step"].shape[0])
                st = np.asarray(f["state/current_step"][:cf]).reshape(-1)
                k = int(np.searchsorted(st, TARGET_MS*100, side="right"))-1
                if k < 0: print("  SKIP", seed); continue
                mc = np.asarray(f["myosin/center"][k:k+1], float)
                ac = np.asarray(f["actin/center"][k:k+1], float)
                ad = np.asarray(f["actin/direction"][k:k+1], float)
                u = ad/np.linalg.norm(ad, axis=2, keepdims=True)
                m = float(band_order(mc[:, :, 0], box)[0])
                z = float(np.maximum(band_order(ac[:, :, 0] - 0.5*La*u[:, :, 0], box),
                                     band_order(ac[:, :, 0] + 0.5*La*u[:, :, 0], box))[0])
                rows.append((seed, m, z, int(st[k])))
        except Exception as e:
            print("  ERR %s: %s" % (seed, e))
    if not rows: print("\n=== %-26s NO DATA" % tag); continue
    res[tag] = rows
    M = np.array([r[1] for r in rows]); Z = np.array([r[2] for r in rows])
    print("\n=== %-26s n=%d" % (tag, len(rows)))
    for s, m, z, stp in rows:
        print("    seed %-8s  M=%.4f (%.1fx floor)   Z=%.4f (%.1fx floor)   step=%d"
              % (s, m, m/0.078, z, z/0.035, stp))
    print("    MEAN         M=%.4f +- %.4f        Z=%.4f +- %.4f"
          % (M.mean(), M.std(ddof=1), Z.mean(), Z.std(ddof=1)))
print("\n=== paired differences (same seeds)")
def paired(t1, t2, idx, name):
    a = {s: v[idx] for s, *v in [(r[0], r[1], r[2]) for r in res[t1]]}
    b = {s: v[idx] for s, *v in [(r[0], r[1], r[2]) for r in res[t2]]}
    ks = sorted(set(a) & set(b))
    d = np.array([b[s]-a[s] for s in ks])
    se = d.std(ddof=1)/np.sqrt(len(d))
    print("  %-10s %-8s %+.4f  (sd %.4f, SE %.4f, t=%6.2f)  improved %d/%d"
          % (t2.split()[0]+"-"+t1.split()[0], name, d.mean(), d.std(ddof=1), se,
             d.mean()/se if se else float('inf'), (d > 0).sum(), len(d)))
names = list(res)
if len(names) == 3:
    for idx, nm in ((0, "M order"), (1, "Z order")):
        paired(names[0], names[1], idx, nm)
    for idx, nm in ((0, "M order"), (1, "Z order")):
        paired(names[1], names[2], idx, nm)
print("\nnoise floors (uniform random, 20000 trials): M (N=320) 0.078 mean / 0.136 p99;"
      "  Z (N=1600) 0.035 mean / 0.061 p99")
