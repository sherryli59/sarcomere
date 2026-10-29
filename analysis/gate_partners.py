"""Completes the causal chain. crosslink_ratio is a MIN over engaged myosin
partners, so each extra partner can only lower it. And partner count correlates
with BETTER alignment. So the gate penalises well-aligned actins through the
partner-count channel even though it rewards them through the angle channel.

Measure crosslink_ratio vs number of engaged myosin partners.
"""
import os, glob, collections
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
import numpy as np, h5py
T0, T1 = 1000, 1500
FULL = "/scratch/users/xsli/sarcomere/barbed_capsule_normal_full_20261002/data"
LAD  = "/scratch/users/xsli/sarcomere/hph_ladder_20261005/data"
ARMS = [("A hph9/spokes3 ", f"{FULL}/*kaa3_seed*_5000ms_capsule_normal.h5"),
        ("C hph27/spokes3", f"{LAD}/*kaa3_seed*_5000ms_hph27.h5")]
for tag, pat in ARMS:
    byn = collections.defaultdict(list)
    for p in sorted(glob.glob(pat)):
        with h5py.File(p, "r", locking=False) as f:
            cf = int(np.asarray(f["meta/committed_frames"]).reshape(-1)[-1])
            cf = min(cf, f["actin/direction"].shape[0], f["state/current_step"].shape[0])
            st = np.asarray(f["state/current_step"][:cf]).reshape(-1)
            i0 = int(np.searchsorted(st, T0*100, "left")); i1 = min(int(np.searchsorted(st, T1*100, "right")), cf)
            fr = list(range(i0, i1, 25))
            cr = np.asarray(f["actin/crosslink_ratio"][fr], float)
            hs = np.asarray(f["state/am_head_state_current"][fr])
            ha = np.asarray(f["state/am_head_bound_actin_current"][fr]).astype(int)
            for t in range(len(fr)):
                bound = (hs[t] >= 1) & (ha[t] >= 0)
                pm = collections.defaultdict(set)
                mi, hh = np.nonzero(bound)
                for m_, h_ in zip(mi.tolist(), hh.tolist()):
                    pm[int(ha[t][m_][h_])].add(m_)
                nA = cr.shape[1]
                cnt = np.zeros(nA, int)
                for ai, ms in pm.items(): cnt[ai] = len(ms)
                for k in range(0, 6):
                    m = cnt == k
                    if m.any(): byn[k].append(cr[t][m])
    print("\n=== %s" % tag)
    print("    engaged myosin partners -> crosslink_ratio (formation window multiplier)")
    prev = None
    for k in sorted(byn):
        v = np.concatenate(byn[k])
        arrow = "" if prev is None else ("  (%+.4f)" % (v.mean()-prev))
        print("      %d partner%s: ratio %.4f  n=%8d%s" % (k, " " if k == 1 else "s", v.mean(), v.size, arrow))
        prev = v.mean()
    z = np.concatenate(byn[0]).mean() if 0 in byn else float("nan")
    hi = max(k for k in byn if k >= 1)
    print("    unengaged (0 partners) %.4f  ->  %d partners %.4f   = %.2fx penalty"
          % (z, hi, np.concatenate(byn[hi]).mean(), np.concatenate(byn[hi]).mean()/z if z else float("nan")))
