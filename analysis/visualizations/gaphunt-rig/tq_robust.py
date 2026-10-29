"""Robust torque closure per span bin: drop actin-frames hit by the aa_steric
near-coincidence singularity (|aa_steric| > 5000 kT), then report means; check
that the listed sources sum to the recorded total, and that the recorded total
predicts the realised rotation (D*A/kT)."""
import numpy as np, sys
D = "/scratch/users/xsli/sarcomere/gaphunt_20261023/torque_out"
SRC = ["am_head", "actn2_align", "actn2_lever", "actn2_perp", "actn2_axial", "actn2_contact", "aa_lateral", "aa_steric", "am_steric"]
for tag in sys.argv[1:]:
    z = dict(np.load("%s/torque_%s.npz" % (D, tag)))
    ok = np.abs(z["aa_steric"]) < 5000
    ok2 = ok & (np.abs(z["TOT"]) < 5000) & (np.abs(z["MTOT"]) < 5000)
    print("\n=== %s: singular-free %.1f%% of pair-frames; |TOT|,|MTOT|<5000 too: %.1f%%" % (tag, 100 * ok.mean(), 100 * ok2.mean()))
    print("%-12s %7s" % ("span", "n") + "".join("%9s" % s[:9] for s in SRC) + "%9s%9s%9s%9s%9s%9s%9s" % ("SUM", "TOT", "w_a", "D*TOT", "MTOT", "w_m", "Dm*MTOT"))
    for lo, hi in ((0, 1), (1, 200), (200, 400), (400, 2000)):
        for th in ((5, 60),):
            s = ok2 & (z["span"] >= lo) & (z["span"] < hi) & (z["th"] >= th[0]) & (z["th"] < th[1])
            if s.sum() < 30: continue
            v = [z[k][s].mean() for k in SRC]
            print("%-12s %7d" % ("%d-%d" % (lo, hi), s.sum()) + "".join("%9.1f" % x for x in v) +
                  "%9.1f%9.1f%9.2f%9.2f%9.1f%9.2f%9.2f" % (sum(v), z["TOT"][s].mean(), z["wa"][s].mean(), 0.101 * z["TOT"][s].mean(),
                                                     z["MTOT"][s].mean(), z["wm"][s].mean(), 0.059 * z["MTOT"][s].mean()))
    s = ok2
    print("pooled: median TOT %.1f, mean TOT %.1f+-%.1f, w_a %.3f, MTOT %.1f+-%.1f, w_m %.3f" % (
        np.median(z["TOT"][s]), z["TOT"][s].mean(), z["TOT"][s].std() / np.sqrt(s.sum()), z["wa"][s].mean(),
        z["MTOT"][s].mean(), z["MTOT"][s].std() / np.sqrt(s.sum()), z["wm"][s].mean()))
    print("fraction of pair-frames by span: " + ", ".join("%s %.0f%%" % (l, 100 * np.mean((z["span"] >= a) & (z["span"] < b))) for l, a, b in (("0", 0, 1), ("1-200", 1, 200), ("200-400", 200, 400), (">=400", 400, 5000))))
