"""Is the tug-of-war real, or just a correlation between partner count and bad
alignment? Two predictions that frustration makes and a mere correlation does not.

P1 WITHIN-ARM STRATIFICATION. If multi-master frustration is the harm, then
   inside ONE arm, actins with 1 partner should be well aligned and alignment
   should degrade monotonically with partner count. And arm C's SINGLE-partner
   actins should look like arm A's single-partner actins -- the extra heads
   should only hurt the actins that acquired extra masters.

P2 BISECTOR TEST. A frustrated actin cannot sit on any one myosin's axis; the
   compromise minimising total angular penalty is near the BISECTOR of its
   partners' axes. So for >=2-partner actins, compare the actin's angle to its
   NEAREST partner against its angle to the bisector. If it hugs the bisector
   rather than a single master, it is being pulled by both.
   Null if frustration is absent: the actin should sit on ONE master's axis
   (angle-to-nearest small) and be far from the bisector.

Also locates the single worst concrete example for visualisation.
"""
import os, glob, collections, itertools, json
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
import numpy as np, h5py
T0, T1 = 1000, 1500
FULL = "/scratch/users/xsli/sarcomere/barbed_capsule_normal_full_20261002/data"
LAD  = "/scratch/users/xsli/sarcomere/hph_ladder_20261005/data"
ARMS = [("A hph9/spokes3 ", f"{FULL}/*kaa3_seed20260805*_5000ms_capsule_normal.h5"),
        ("B hph9/spokes1 ", f"{LAD}/*kaa3_seed20260805*_5000ms_spokes1.h5"),
        ("C hph27/spokes3", f"{LAD}/*kaa3_seed20260805*_5000ms_hph27.h5")]

def ang(a, b):
    return float(np.degrees(np.arccos(np.clip(abs(float(np.dot(a, b))), 0, 1))))

best = None
for tag, pat in ARMS:
    byn = collections.defaultdict(list)        # partner count -> angle to nearest
    bis = []                                   # (angle to nearest, angle to bisector)
    for p in sorted(glob.glob(pat)):
        with h5py.File(p, "r", locking=False) as f:
            cf = int(np.asarray(f["meta/committed_frames"]).reshape(-1)[-1])
            cf = min(cf, f["actin/direction"].shape[0], f["state/current_step"].shape[0])
            st = np.asarray(f["state/current_step"][:cf]).reshape(-1)
            i0 = int(np.searchsorted(st, T0*100, "left")); i1 = min(int(np.searchsorted(st, T1*100, "right")), cf)
            fr = list(range(i0, i1, 20))
            u = np.asarray(f["actin/direction"][fr], float); u /= np.linalg.norm(u, axis=2, keepdims=True)
            w = np.asarray(f["myosin/direction"][fr], float); w /= np.linalg.norm(w, axis=2, keepdims=True)
            ac = np.asarray(f["actin/center"][fr], float)
            mc = np.asarray(f["myosin/center"][fr], float)
            hs = np.asarray(f["state/am_head_state_current"][fr])
            ha = np.asarray(f["state/am_head_bound_actin_current"][fr]).astype(int)
            for t in range(len(fr)):
                bound = (hs[t] >= 1) & (ha[t] >= 0)
                pm = collections.defaultdict(set)
                mi, hh = np.nonzero(bound)
                for m_, h_ in zip(mi.tolist(), hh.tolist()):
                    pm[int(ha[t][m_][h_])].add(m_)
                for ai, ms in pm.items():
                    ms = sorted(ms); ua = u[t][ai]
                    angs = [ang(ua, w[t][m]) for m in ms]
                    byn[min(len(ms), 4)].append(min(angs))
                    if len(ms) >= 2:
                        # sign-align the partner axes to a common hemisphere, then bisect
                        ref = w[t][ms[0]]
                        axes = [w[t][m] if np.dot(w[t][m], ref) >= 0 else -w[t][m] for m in ms]
                        b = np.mean(axes, axis=0); nb = np.linalg.norm(b)
                        if nb < 1e-9: continue
                        b /= nb
                        bis.append((min(angs), ang(ua, b)))
                        if len(ms) >= 3:
                            spread = max(ang(w[t][x], w[t][y]) for x, y in itertools.combinations(ms, 2))
                            score = spread * min(angs)
                            if tag.startswith("C") and (best is None or score > best["score"]):
                                best = dict(score=float(score), arm=tag.strip(), frame=int(fr[t]),
                                            step=int(st[fr[t]]), actin=int(ai), myosins=[int(m) for m in ms],
                                            angles_to_myosins=[round(a, 2) for a in angs],
                                            max_pair_spread=round(spread, 2),
                                            angle_to_bisector=round(ang(ua, b), 2),
                                            actin_dir=[round(float(x), 5) for x in ua],
                                            actin_center=[round(float(x), 5) for x in ac[t][ai]],
                                            myosin_dirs=[[round(float(x), 5) for x in w[t][m]] for m in ms],
                                            myosin_centers=[[round(float(x), 5) for x in mc[t][m]] for m in ms])
    print("\n=== %s" % tag)
    print("  P1 angle to NEAREST bound myosin, stratified by partner count:")
    for k in sorted(byn):
        v = np.array(byn[k])
        print("     %d partner%s: n=%7d  mean %6.2f  median %6.2f" % (k, " " if k == 1 else "s", v.size, v.mean(), np.median(v)))
    if bis:
        b = np.array(bis)
        print("  P2 >=2 partners (n=%d): angle to NEAREST master %6.2f   angle to BISECTOR %6.2f"
              % (b.shape[0], b[:, 0].mean(), b[:, 1].mean()))
        print("     closer to the bisector than to any single master: %.1f%% of cases"
              % (100*(b[:, 1] < b[:, 0]).mean()))
if best:
    print("\n=== worst concrete example in arm C (max  pair-spread x angle-to-nearest)")
    for k, v in best.items(): print("   %-20s %s" % (k, v))
    with open("/scratch/users/xsli/sarcomere/hph_ladder_20261005/analysis/frustration_example.json", "w") as fh:
        json.dump(best, fh, indent=1)
    print("   written frustration_example.json")
