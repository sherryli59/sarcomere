"""The formation gate is a THRESHOLD, not a linear throttle.

  contour_lower = -L/2
  contour_upper = min(-L/2 + aa_barbed_zone, -L/2 + ratio*L)
  and formation is refused entirely unless contour_upper > contour_lower.

With aa_barbed_zone=0.2 and L=1.0 the ratio term binds only when ratio < 0.2,
and kills formation at ratio == 0. So the arm-mean ratio (0.68 vs 0.42) is
IRRELEVANT -- both are in the inert region. What matters is the mass below 0.2,
and especially the mass at exactly 0.
"""
import os, glob
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
    allr = []
    for p in sorted(glob.glob(pat)):
        with h5py.File(p, "r", locking=False) as f:
            a = f["meta"].attrs
            bz = float(a["aa_barbed_zone"]); L = float(a["actin_length"])
            cf = int(np.asarray(f["meta/committed_frames"]).reshape(-1)[-1])
            cf = min(cf, f["actin/crosslink_ratio"].shape[0], f["state/current_step"].shape[0])
            st = np.asarray(f["state/current_step"][:cf]).reshape(-1)
            i0 = int(np.searchsorted(st, T0*100, "left")); i1 = min(int(np.searchsorted(st, T1*100, "right")), cf)
            allr.append(np.asarray(f["actin/crosslink_ratio"][i0:i1:25], float).ravel())
    r = np.concatenate(allr); thr = bz/L
    # effective window length the gate permits, per actin
    win = np.minimum(bz, np.clip(r, 0, 1)*L)
    res[tag[0]] = dict(r=r, win=win, thr=thr, bz=bz)
    print("\n=== %s   (n=%d actin-frames, gate threshold ratio = %.2f)" % (tag, r.size, thr))
    print("    ratio mean %.4f   median %.4f" % (r.mean(), np.median(r)))
    print("    mass at exactly 0 (NO formation possible): %6.3f%%" % (100*(r <= 1e-12).mean()))
    for t in (0.05, 0.10, thr, 0.30, 0.50):
        print("    fraction below %.2f%s: %6.3f%%" % (t, "  <-- GATE BITES BELOW THIS" if abs(t-thr) < 1e-9 else "      ",
              100*(r < t).mean()))
    print("    MEAN PERMITTED WINDOW: %.4f um of the %.2f um barbed zone  (%.1f%% of full)"
          % (win.mean(), bz, 100*win.mean()/bz))
if "A" in res and "C" in res:
    A, C = res["A"], res["C"]
    print("\n=== C vs A")
    print("  mean permitted window ratio C/A = %.4f" % (C["win"].mean()/A["win"].mean()))
    print("  observed formation rate ratio    = 0.631")
    print("  mass at ratio 0:  A %.3f%%  ->  C %.3f%%" % (100*(A["r"] <= 1e-12).mean(), 100*(C["r"] <= 1e-12).mean()))
    print("  mass below 0.2 :  A %.3f%%  ->  C %.3f%%" % (100*(A["r"] < A["thr"]).mean(), 100*(C["r"] < C["thr"]).mean()))
