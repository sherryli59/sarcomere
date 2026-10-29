"""What is wrong with the ACTN2 formation gate?

The gate (sarcomere_actn2.cpp:295) sizes the attachment window as
  contour_lower + crosslink_ratio * actin.length
and crosslink_ratio is the run of actin from its BARBED TIP to the first point
within the capture radius R of any thick filament, minimised over myosins.

R = actin_radius 0.015 + am_surface_gap 0.015 + am_cutoff 0.050 = 0.080 um, and
am_cutoff is myosin HEAD REACH. So the gate penalises wherever a head COULD
reach, not where one IS, and not penetration.

THE TEST: an actin lying alongside a thick filament -- the well-assembled
geometry we want -- has most of its contour inside that 80 nm envelope. If
crosslink_ratio FALLS as the actin becomes better aligned, the gate actively
penalises the target state. Measure crosslink_ratio vs angle to nearest myosin.
"""
import os, glob
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
import numpy as np, h5py
T0, T1 = 1000, 1500
FULL = "/scratch/users/xsli/sarcomere/barbed_capsule_normal_full_20261002/data"
LAD  = "/scratch/users/xsli/sarcomere/hph_ladder_20261005/data"
ARMS = [("A hph9/spokes3 ", f"{FULL}/*kaa3_seed*_5000ms_capsule_normal.h5"),
        ("C hph27/spokes3", f"{LAD}/*kaa3_seed*_5000ms_hph27.h5")]
BINS = [(0,10),(10,20),(20,30),(30,45),(45,60),(60,75),(75,90)]
for tag, pat in ARMS:
    acc = {b: [] for b in BINS}; allr = []
    for p in sorted(glob.glob(pat)):
        with h5py.File(p, "r", locking=False) as f:
            a = f["meta"].attrs
            cf = int(np.asarray(f["meta/committed_frames"]).reshape(-1)[-1])
            cf = min(cf, f["actin/direction"].shape[0], f["state/current_step"].shape[0])
            st = np.asarray(f["state/current_step"][:cf]).reshape(-1)
            i0 = int(np.searchsorted(st, T0*100, "left")); i1 = min(int(np.searchsorted(st, T1*100, "right")), cf)
            fr = list(range(i0, i1, 25))
            u = np.asarray(f["actin/direction"][fr], float); u /= np.linalg.norm(u, axis=2, keepdims=True)
            w = np.asarray(f["myosin/direction"][fr], float); w /= np.linalg.norm(w, axis=2, keepdims=True)
            ac = np.asarray(f["actin/center"][fr], float); mc = np.asarray(f["myosin/center"][fr], float)
            cr = np.asarray(f["actin/crosslink_ratio"][fr], float)
            box = np.array([a["box_x"], a["box_y"], a["box_z"]], float)
            for t in range(len(fr)):
                d = ac[t][:, None, :] - mc[t][None, :, :]
                d -= box*np.rint(d/box)
                near = np.linalg.norm(d, axis=2).argmin(axis=1)
                th = np.degrees(np.arccos(np.clip(np.abs(np.einsum("ak,ak->a", u[t], w[t][near])), 0, 1)))
                allr.append(cr[t])
                for lo, hi in BINS:
                    m = (th >= lo) & (th < hi)
                    if m.any(): acc[(lo, hi)].append(cr[t][m])
    R = np.concatenate(allr)
    print("\n=== %s   overall crosslink_ratio mean %.4f" % (tag, R.mean()))
    print("    angle to nearest myosin  ->  crosslink_ratio (the formation window multiplier)")
    for b in BINS:
        if not acc[b]: continue
        v = np.concatenate(acc[b])
        bar = "#"*int(round(40*v.mean()))
        print("      theta [%2d,%2d): ratio %.4f  n=%7d  %s" % (b[0], b[1], v.mean(), v.size, bar))
    lo = np.concatenate(acc[(0,10)]) if acc[(0,10)] else None
    hi = np.concatenate(acc[(75,90)]) if acc[(75,90)] else None
    if lo is not None and hi is not None:
        print("    WELL-ALIGNED (0-10 deg) ratio %.4f  vs  PERPENDICULAR (75-90) %.4f  -> penalty factor %.2fx"
              % (lo.mean(), hi.mean(), hi.mean()/lo.mean() if lo.mean() else float("nan")))
