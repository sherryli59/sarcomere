"""Why does hph=27 repair a 60 deg actin in the two-body probe but make
production WORSE? Three hypotheses, measured on the existing data.

H1 COUPLE ABSENT: the head couple simply does not exist at production density.
   Test: am_align_drive_minus/am_align_head_count_minus per myosin, binned by
   that myosin's mean bound-actin angle. The probe gave +21.4 kT at theta
   [50,70). If production reads ~0 there, the couple is absent, not overwhelmed.

H2 FRUSTRATION: each actin is pulled by SEVERAL myosins whose axes disagree, so
   the couples cancel. The probe has exactly one myosin, so this cannot happen
   there. Test: distinct myosin partners per bound actin, and the angular spread
   between those partners' directors.

H3 COMPETITION: ACTN2 is the real aligner (it is absent from the probe, k_aa=0)
   and extra heads pin against it. Test: ACTN2 links per actin alongside the
   above.
"""
import os, glob, collections, itertools
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
import numpy as np, h5py
KT = 4.114e-3
T0, T1 = 1000, 1500
BINS = [(0, 10), (10, 30), (30, 50), (50, 70), (70, 90)]
FULL = "/scratch/users/xsli/sarcomere/barbed_capsule_normal_full_20261002/data"
LAD  = "/scratch/users/xsli/sarcomere/hph_ladder_20261005/data"
ARMS = [("A hph9/spokes3 ", f"{FULL}/*kaa3_seed*_5000ms_capsule_normal.h5"),
        ("B hph9/spokes1 ", f"{LAD}/*kaa3_seed*_5000ms_spokes1.h5"),
        ("C hph27/spokes3", f"{LAD}/*kaa3_seed*_5000ms_hph27.h5")]

for tag, pat in ARMS:
    binvals = {b: [] for b in BINS}; part = []; spread = []; lpa = []
    for p in sorted(glob.glob(pat)):
        with h5py.File(p, "r", locking=False) as f:
            a = f["meta"].attrs
            cf = int(np.asarray(f["meta/committed_frames"]).reshape(-1)[-1])
            cf = min(cf, f["actin/direction"].shape[0], f["state/current_step"].shape[0])
            st = np.asarray(f["state/current_step"][:cf]).reshape(-1)
            i0 = int(np.searchsorted(st, T0*100, "left")); i1 = min(int(np.searchsorted(st, T1*100, "right")), cf)
            fr = list(range(i0, i1, 10))            # every 10th frame in the window
            if len(fr) < 5: print("  SKIP", p.split("seed")[1][:8]); continue
            u = np.asarray(f["actin/direction"][fr], float)
            u /= np.linalg.norm(u, axis=2, keepdims=True)
            w = np.asarray(f["myosin/direction"][fr], float)
            w /= np.linalg.norm(w, axis=2, keepdims=True)
            hs = np.asarray(f["state/am_head_state_current"][fr])
            ha = np.asarray(f["state/am_head_bound_actin_current"][fr]).astype(int)
            d = np.asarray(f["myosin/am_align_drive_minus"][fr], float).reshape(len(fr), -1)
            n = np.asarray(f["myosin/am_align_head_count_minus"][fr], float).reshape(len(fr), -1)
            links = np.asarray(f["state/aa_pairs_current_count"][fr]).reshape(len(fr), -1)[:, -1]
            lpa.append(float(links.mean())/u.shape[1])
            assert hs.shape[1] == u.shape[1] or hs.shape[1] == w.shape[1], hs.shape
            nM = w.shape[1]
            assert n.shape[1] == nM, (n.shape, nM)
            for t in range(len(fr)):
                bound = (hs[t] >= 1) & (ha[t] >= 0)   # (nM, 2*hph), already per-frame
                assert bound.shape == (nM, hs.shape[2]), bound.shape
                # ---- H1: drive/head per myosin, binned by that myosin's mean bound angle
                for m in np.nonzero(n[t] > 0)[0]:
                    ai = np.unique(ha[t][m][bound[m]])
                    if ai.size == 0: continue
                    th = np.degrees(np.arccos(np.clip(np.abs(u[t][ai] @ w[t][m]), 0, 1))).mean()
                    for lo, hi in BINS:
                        if lo <= th < hi: binvals[(lo, hi)].append(d[t][m]/n[t][m]/KT); break
                # ---- H2: partners per actin and the spread of their axes
                pm = collections.defaultdict(set)
                mi, hi_ = np.nonzero(bound)
                for m_, h_ in zip(mi.tolist(), hi_.tolist()):
                    pm[int(ha[t][m_][h_])].add(m_)
                for ai_, ms in pm.items():
                    part.append(len(ms))
                    if len(ms) >= 2:
                        ms = sorted(ms)
                        ang = [np.degrees(np.arccos(np.clip(abs(float(w[t][x] @ w[t][y])), 0, 1)))
                               for x, y in itertools.combinations(ms, 2)]
                        spread.append(float(np.mean(ang)))
    print("\n=== %s" % tag)
    print("  H1 head couple, drive/head (kT), binned by the myosin's mean bound-actin angle:")
    for b in BINS:
        v = np.array(binvals[b])
        print("     theta [%2d,%2d): %s" % (b[0], b[1],
              ("%+7.2f  (n=%6d myosin-frames)" % (v.mean(), v.size)) if v.size else "    --        (none)"))
    pa = np.array(part)
    print("  H2 distinct myosin partners per bound actin: mean %.2f   >=2 partners %.1f%%   >=3 %.1f%%"
          % (pa.mean(), 100*(pa >= 2).mean(), 100*(pa >= 3).mean()))
    sp = np.array(spread)
    if sp.size:
        print("     angular spread between those partners' axes: mean %.1f deg  median %.1f deg  (n=%d actins)"
              % (sp.mean(), np.median(sp), sp.size))
    print("  H3 ACTN2 links per actin: %.3f" % float(np.mean(lpa)))
