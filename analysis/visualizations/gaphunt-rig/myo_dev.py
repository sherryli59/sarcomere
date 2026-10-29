"""Per myosin: deviation from the head-weighted average of its own bound actins,
versus candidate causes: bound heads, distinct actins, axial span of its bound
heads, disagreement between its two halves' actin averages, titin links, and
whether its ends are near a Z-disc (barbed-end density).  Spearman correlations
and binned medians."""
import h5py, numpy as np, os, sys
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
CAMP = "/scratch/users/xsli/sarcomere/gaphunt_20261023/data"
def U(d): return d / np.linalg.norm(d, axis=-1, keepdims=True)
def lead(Q): w, v = np.linalg.eigh(Q); return v[..., -1]
ang = lambda a, b: np.degrees(np.arccos(np.clip(np.abs(np.sum(a * b, -1)), 0, 1)))
def rank(x): return np.argsort(np.argsort(x)).astype(float)
def spear(x, y): return float(np.corrcoef(rank(x), rank(y))[0, 1])
for fn in ("hph27m24", "hph81m24"):
    f = h5py.File("%s/%s_seed20260805.h5" % (CAMP, fn), "r", locking=False); a = dict(f["meta"].attrs); hph = int(a["am_heads_per_half"])
    Lm = float(a["myosin_length"]); cc = (hph + 2) // 3
    anc = np.array([(-1 if h < hph else 1) * (0.075 + (((h % hph) // 3 + 0.5) / cc) * 0.675) for h in range(2 * hph)])
    rows = []
    for t in range(360, 401, 10):
        u = U(f["actin/direction"][t]); m = U(f["myosin/direction"][t])
        hs = f["state/am_head_state_current"][t]; hb = f["state/am_head_bound_actin_current"][t].astype(int)
        tl = f["state/titin_bound_actin_current"][t] if "state/titin_bound_actin_current" in f else None
        mi, hi = np.nonzero((hs == 1) | (hs == 2)); ai = hb[mi, hi]; half = (hi >= hph).astype(int)
        Q = np.zeros((len(m), 2, 3, 3)); c = np.zeros((len(m), 2)); np.add.at(Q, (mi, half), u[ai][:, :, None] * u[ai][:, None, :]); np.add.at(c, (mi, half), 1)
        for j in range(len(m)):
            if c[j].sum() < 2: continue
            sel = mi == j; avg = lead(Q[j].sum(0) / c[j].sum())
            hd = ang(lead(Q[j, 0] / c[j, 0]), lead(Q[j, 1] / c[j, 1])) if (c[j] > 0).all() else np.nan
            rows.append((float(ang(m[j], avg)), c[j].sum(), len(set(ai[sel])), 1000 * (anc[hi[sel]].max() - anc[hi[sel]].min()),
                         hd, min(c[j]) / max(c[j]) if max(c[j]) > 0 else 0))
    R = np.array(rows, float); dev = R[:, 0]
    names = ["bound heads", "distinct actins", "head span on myosin (nm)", "half-average disagreement (deg)", "half balance (min/max heads)"]
    print("== %s: %d myosin-frames, deviation from own-actin average: median %.1f deg" % (fn, len(R), np.median(dev)))
    for i, nme in enumerate(names, start=1):
        x = R[:, i]; ok = np.isfinite(x)
        qs = np.nanpercentile(x[ok], [0, 25, 50, 75, 100])
        bins = [np.median(dev[ok][(x[ok] >= lo) & (x[ok] <= hi)]) for lo, hi in zip(qs[:-1], qs[1:])]
        print("   %-32s spearman %+.2f | deviation by quartile: %s  (quartile edges %s)" % (nme, spear(x[ok], dev[ok]), " ".join("%.1f" % b for b in bins), " ".join("%.0f" % q for q in qs)))
