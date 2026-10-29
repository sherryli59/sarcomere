"""Summarise torque_an npz: mean drive toward the bound partner (kT), by span and angle."""
import numpy as np, sys, glob, os, json
D = "/scratch/users/xsli/sarcomere/gaphunt_20261023/torque_out"
SRC = ["am_head", "RAD", "AX", "LIN", "OTHERPART", "actn2_align", "actn2_lever", "actn2_perp", "actn2_axial",
       "actn2_contact", "aa_lateral", "aa_steric", "am_steric", "TOT", "MTOT"]
def m(x): return (np.mean(x), np.std(x) / np.sqrt(max(1, len(x)))) if len(x) else (np.nan, np.nan)
tags = sys.argv[1:] or [os.path.basename(p)[7:-4] for p in sorted(glob.glob(D + "/torque_*.npz"))]
for tag in tags:
    z = dict(np.load("%s/torque_%s.npz" % (D, tag)))
    print("\n=== %s  (%d pair-frames)" % (tag, len(z["th"])))
    for name, sel0 in (("ALL pairs", np.ones(len(z["th"]), bool)), ("SINGLE-partner actins", z["npart"] == 1)):
        print("-- %s: n %d; median partners/actin %.0f, actins/myosin %.0f" % (name, sel0.sum(), np.median(z["npart"][sel0]), np.median(z["mpart"][sel0])))
        hdr = "%-22s %7s %6s %6s" % ("span | theta", "n", "nHead", "strad") + "".join("%10s" % s for s in SRC) + "%9s%9s%9s%9s" % ("w_a", "Dw_pred", "w_m", "dth/dt")
        print(hdr)
        for lo, hi in ((0, 1), (1, 200), (200, 400), (400, 2000)):
            for tlo, thi in ((5, 15), (15, 30), (30, 60)):
                sel = sel0 & (z["span"] >= lo) & (z["span"] < hi) & (z["th"] >= tlo) & (z["th"] < thi)
                if sel.sum() < 30: continue
                strad = np.mean((z["npos"][sel] > 0) & (z["nneg"][sel] > 0))
                row = "%-22s %7d %6.1f %6.2f" % ("%d-%d nm | %d-%d" % (lo, hi, tlo, thi), sel.sum(), z["n"][sel].mean(), strad)
                row += "".join("%10.1f" % m(z[s][sel])[0] for s in SRC)
                row += "%9.2f%9.2f%9.2f%9.1f" % (z["wa"][sel].mean(), 0.101 * z["TOT"][sel].mean(), z["wm"][sel].mean(), z["dth"][sel].mean())
                print(row)
        # SEMs for the key columns in one representative bin
    sel = (z["span"] >= 200) & (z["th"] >= 15) & (z["th"] < 30)
    print("SEM check (span>=200, 15-30 deg, n %d): " % sel.sum() + ", ".join("%s %.1f+-%.1f" % ((s,) + m(z[s][sel])) for s in ("am_head", "RAD", "TOT", "MTOT", "aa_lateral", "actn2_align")))
