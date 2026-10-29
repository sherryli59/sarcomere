"""WHY does spokes=1 cost Z-disc banding?

HYPOTHESIS (untested until now): with 3 crowns 225 nm apart the heads act as a
coarse positional GRATING -- an actin gripped by a crown is pinned at one of few
discrete axial offsets, so actin tips cluster at periodic positions. With 9
crowns 75 nm apart the allowed offsets are 3x finer, the grating washes out, and
longitudinal register is lost.

TEST 1: distribution of the actin barbed tip's axial offset from its bound
myosin's centre, projected on the myosin axis. The grating hypothesis predicts
DISCRETE structure at the crown spacing -- peaks ~225 nm apart in spokes=3 and
~75 nm apart in spokes=1 -- and therefore a flatter, less modulated distribution
in spokes=1. Quantify with the harmonic power at the crown spacing.

TEST 2: the barbed-tip axial histogram in the box, which is what Z order sees.

COMPLICATION the simple hypothesis does not explain: M order also fell 17%, so
the MYOSINS are less registered too. Measure myosin centre spacing as well.
"""
import os, glob, collections
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
import numpy as np, h5py
T0, T1 = 3000, 5000
BARE, BEAR = 0.10, 0.5*(1.0-0.10)
D = "/scratch/users/xsli/sarcomere"
ARMS = [("A spokes=3 (225 nm)", f"{D}/barbed_capsule_normal_full_20261002/data/*kaa3_seed*_5000ms_capsule_normal.h5"),
        ("B spokes=1 ( 75 nm)", f"{D}/hph_ladder_20261005/data/*kaa3_seed*_5000ms_spokes1.h5")]
for tag, pat in ARMS:
    off_all = []; tipx = []; myox = []; spacing = None
    for p in sorted(glob.glob(pat)):
        with h5py.File(p, "r", locking=False) as f:
            a = f["meta"].attrs
            hph = int(a["am_heads_per_half"]); spk = int(float(a.get("am_spokes_per_crown", 3)))
            Lm = float(a["myosin_length"]); L = float(a["actin_length"]); box = float(a["box_x"])
            cc = max(1, (hph + spk - 1)//spk); spacing = BEAR*Lm/cc
            cf = int(np.asarray(f["meta/committed_frames"]).reshape(-1)[-1])
            cf = min(cf, f["actin/center"].shape[0], f["state/current_step"].shape[0])
            st = np.asarray(f["state/current_step"][:cf]).reshape(-1)
            i0 = int(np.searchsorted(st, T0*100, "left")); i1 = min(int(np.searchsorted(st, T1*100, "right")), cf)
            if i1 - i0 < 5: continue
            fr = list(range(i0, i1, 40))
            u = np.asarray(f["actin/direction"][fr], float); u /= np.linalg.norm(u, axis=2, keepdims=True)
            w = np.asarray(f["myosin/direction"][fr], float); w /= np.linalg.norm(w, axis=2, keepdims=True)
            ac = np.asarray(f["actin/center"][fr], float); mc = np.asarray(f["myosin/center"][fr], float)
            hs = np.asarray(f["state/am_head_state_current"][fr])
            ha = np.asarray(f["state/am_head_bound_actin_current"][fr]).astype(int)
            barbed = ac - 0.5*L*u
            for t in range(len(fr)):
                tipx.append(barbed[t][:, 0] % box)
                myox.append(mc[t][:, 0] % box)
                mi, hh = np.nonzero((hs[t] >= 1) & (ha[t] >= 0))
                for m_, h_ in zip(mi.tolist(), hh.tolist()):
                    ai = int(ha[t][m_][h_])
                    d = barbed[t][ai] - mc[t][m_]
                    d -= box*np.round(d[0]/box)*np.array([1.0, 0, 0])
                    off_all.append(float(np.dot(d, w[t][m_])))
    o = np.array(off_all); tx = np.concatenate(tipx); mx = np.concatenate(myox)
    print("\n=== %s   crown spacing %.4f um   (n=%d bound pairs, %d-%d ms)" % (tag, spacing, o.size, T0, T1))
    # harmonic power of the offset distribution AT the crown spacing
    for lam, nm in ((spacing, "crown spacing"), (0.225, "225 nm"), (0.075, "75 nm")):
        amp = abs(np.mean(np.exp(2j*np.pi*o/lam)))
        print("    modulation of the AXIAL OFFSET at wavelength %.3f um (%-13s): %.4f" % (lam, nm, amp))
    print("    offset spread: sd %.4f um   p10 %.3f   p90 %.3f" % (o.std(), np.percentile(o,10), np.percentile(o,90)))
    # the box-scale banding that Z order measures, harmonics 1-4 of the 6 um box
    print("    barbed-tip box harmonics |mean exp(2 pi i h x / L)|:")
    for h in (1,2,3,4):
        print("        h=%d (lambda=%.2f um): tips %.4f   myosin centres %.4f"
              % (h, 6.0/h, abs(np.mean(np.exp(2j*np.pi*h*tx/6.0))), abs(np.mean(np.exp(2j*np.pi*h*mx/6.0)))))
