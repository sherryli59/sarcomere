"""Does ACTN2 link survival already depend on the A-A angle or the A-M angle?
Production control r27ctrl, 2.0 -> 2.4 s.  For every live link at frame t (5 ms apart), loss =
link id absent at t+1.  Binned by
  * A-A angle: unsigned angle between the two linked actins (antiparallel = 0)
  * A-M angle: for each linked actin, angle to the head-weighted axis of the myosins holding
    STRONG heads on it; pair value = the larger of the two (NaN if an actin has no strong head)
and by status.  Also: the A-A / A-M angles of links that FORM (id new at t+1)."""
import h5py, numpy as np, os, json, collections
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
f = h5py.File("/scratch/users/xsli/sarcomere/gaphunt_20261023/data/r27ctrl_seed20260805.h5", "r", locking=False)
U = lambda d: d / np.linalg.norm(d, axis=-1, keepdims=True)
ang = lambda a, b: np.degrees(np.arccos(np.clip(np.abs(np.sum(a * b, -1)), 0, 1)))
off = f["catch_bond/local_link_frame_offsets"]; ident = f["catch_bond/local_link_identity"]
def links(t):
    o, c = int(off[t][1]), int(off[t][2]); r = ident[o:o + c]
    return {(int(x[0]), int(x[1])): (int(x[2]), int(x[3]), int(x[4])) for x in r}
def am_axes(t, u):
    m = U(f["myosin/direction"][t]); hs = f["state/am_head_state_current"][t]; hb = f["state/am_head_bound_actin_current"][t].astype(int)
    mi, hi = np.nonzero(hs == 2); ai = hb[mi, hi]
    Q = np.zeros((len(u), 3, 3)); n = np.zeros(len(u))
    np.add.at(Q, ai, m[mi][:, :, None] * m[mi][:, None, :]); np.add.at(n, ai, 1)
    th = np.full(len(u), np.nan); ok = n > 0
    ax = np.linalg.eigh(Q[ok] / n[ok][:, None, None])[1][..., -1]
    th[ok] = ang(u[ok], ax); return th
rows, born = [], []
L0 = links(400)
for t in range(400, 480):
    L1 = links(t + 1)
    u = U(f["actin/direction"][t]); th = am_axes(t, u)
    for k, (i, j, st) in L0.items():
        rows.append((ang(u[i], u[j]), np.nanmax([th[i], th[j]]) if not (np.isnan(th[i]) and np.isnan(th[j])) else np.nan,
                     float(np.isnan(th[i]) or np.isnan(th[j])), st, float(k not in L1)))
    u1 = U(f["actin/direction"][t + 1]); th1 = am_axes(t + 1, u1)
    for k, (i, j, st) in L1.items():
        if k not in L0:
            born.append((ang(u1[i], u1[j]), np.nanmax([th1[i], th1[j]]) if not (np.isnan(th1[i]) and np.isnan(th1[j])) else np.nan))
    L0 = L1
R = np.array(rows); B = np.array(born)
aa, am, nostrong, st, lost = R.T
def h(sel): return (100 * lost[sel].mean(), int(sel.sum())) if sel.sum() else (np.nan, 0)
print("link-frames %d; mean loss per 5 ms %.2f%% (status-1 %.2f%%, status-2 %.2f%%); links formed %d" % (len(R), 100 * lost.mean(), 100 * lost[st == 1].mean(), 100 * lost[st == 2].mean(), len(B)))
print("\nloss per 5 ms by A-A angle (status 1 | status 2):")
for lo, hi in ((0, 5), (5, 10), (10, 20), (20, 40), (40, 91)):
    s = (aa >= lo) & (aa < hi)
    print("   A-A %2d-%2d deg: share %4.1f%% | %5.2f%% (n %6d) | %5.2f%% (n %6d)" % (lo, hi, 100 * s.mean(), *h(s & (st == 1)), *h(s & (st == 2))))
print("\nloss per 5 ms by A-M angle (larger of the two actins) (status 1 | status 2):")
for lo, hi in ((0, 10), (10, 20), (20, 40), (40, 91)):
    s = (am >= lo) & (am < hi)
    print("   A-M %2d-%2d deg: share %4.1f%% | %5.2f%% (n %6d) | %5.2f%% (n %6d)" % (lo, hi, 100 * s.mean(), *h(s & (st == 1)), *h(s & (st == 2))))
s = nostrong > 0
print("   an actin with no strong head: share %4.1f%% | %5.2f%% | %5.2f%%" % (100 * s.mean(), h(s & (st == 1))[0], h(s & (st == 2))[0]))
print("\nangles of links at formation vs the standing population (median, p90):")
print("   A-A: new %.1f / %.1f deg vs existing %.1f / %.1f" % (np.median(B[:, 0]), np.percentile(B[:, 0], 90), np.median(aa), np.percentile(aa, 90)))
print("   A-M: new %.1f / %.1f deg vs existing %.1f / %.1f" % (np.nanmedian(B[:, 1]), np.nanpercentile(B[:, 1], 90), np.nanmedian(am), np.nanpercentile(am, 90)))
json.dump(dict(rows=R.tolist()[:0]), open("/dev/null", "w"))
