"""Do heads actually invade the ACTN2 barbed zone, and is arm C's ACTN2 loss
from more DELETIONS or fewer FORMATIONS?

The barbed end is the MINUS tip, so with actin_length=1.0 the material
coordinate runs -0.5 (barbed) to +0.5 and aa_barbed_zone=0.2 puts the ACTN2
capture zone at c in [-0.5, -0.3]. am_head_target_axial_current is the head's
attachment coordinate on the actin, clipped to +-0.5.
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
EDGES = np.array([-0.5, -0.4, -0.3, -0.1, 0.1, 0.3, 0.5])

for tag, pat in ARMS:
    coords = []; inzone = []; ratio = []; pen = []
    reasons = collections.Counter(); nev = 0; legend = None
    for p in sorted(glob.glob(pat)):
        with h5py.File(p, "r", locking=False) as f:
            a = f["meta"].attrs
            legend = legend or a.get("aa_local_link_detachment_reason_codes")
            cf = int(np.asarray(f["meta/committed_frames"]).reshape(-1)[-1])
            cf = min(cf, f["state/am_head_state_current"].shape[0])
            st = np.asarray(f["state/current_step"][:cf]).reshape(-1)
            i0 = int(np.searchsorted(st, T0*100, "left")); i1 = min(int(np.searchsorted(st, T1*100, "right")), cf)
            fr = list(range(i0, i1, 20))
            hs = np.asarray(f["state/am_head_state_current"][fr])
            tg = np.asarray(f["state/am_head_target_axial_current"][fr], float)
            m = hs >= 1
            c = tg[m]
            coords.append(c)
            inzone.append(float(((c >= -0.5) & (c <= -0.3)).mean()))
            if "crosslink_ratio" in f["actin"]:
                ratio.append(float(np.asarray(f["actin/crosslink_ratio"][fr], float).mean()))
            if "aa_barbed_tip_penetration" in f["actin"]:
                pen.append(float(np.asarray(f["actin/aa_barbed_tip_penetration"][fr], float).mean()))
            # detachment reasons over the whole run so far
            ev = f["catch_bond/local_link_events_v3"]
            cols = a.get("aa_local_link_events_v3_columns")
            cols = cols.decode().split(",") if cols is not None else None
            if cols and "detachment_reason" in cols:
                j = cols.index("detachment_reason")
                r = np.asarray(ev[:, j], int)
                nev += r.size
                for k, v in zip(*np.unique(r, return_counts=True)):
                    reasons[int(k)] += int(v)
    c = np.concatenate(coords)
    h, _ = np.histogram(c, bins=EDGES)
    print("\n=== %s   (%d bound head-instances)" % (tag, c.size))
    print("    head attachment coordinate along actin (-0.5 = BARBED tip, +0.5 = pointed):")
    for k in range(len(EDGES)-1):
        bar = "#"*int(round(40*h[k]/h.sum()))
        print("      [%+.1f,%+.1f) %5.1f%%  %s%s" % (EDGES[k], EDGES[k+1], 100*h[k]/h.sum(), bar,
              "   <-- ACTN2 zone" if EDGES[k+1] <= -0.3 else ""))
    print("    fraction inside the ACTN2 barbed zone [-0.5,-0.3]: %.1f%% (per-seed mean %.1f%%)"
          % (100*float(((c >= -0.5) & (c <= -0.3)).mean()), 100*float(np.mean(inzone))))
    if ratio: print("    mean actin/crosslink_ratio: %.4f" % float(np.mean(ratio)))
    if pen:   print("    mean aa_barbed_tip_penetration: %.6f um" % float(np.mean(pen)))
    if nev:
        print("    ACTN2 detachment reasons over the run (n=%d events):" % nev)
        for k in sorted(reasons):
            print("      %d: %7d  %5.1f%%" % (k, reasons[k], 100*reasons[k]/nev))
        print("    legend:", legend)
