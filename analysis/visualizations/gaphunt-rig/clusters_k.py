"""Component sizes when an actin-myosin edge needs >= k strong heads on that pair."""
import h5py, numpy as np, os, collections, sys
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from clusters_an import comps
f = h5py.File("/scratch/users/xsli/sarcomere/gaphunt_20261023/data/r27ctrl_seed20260805.h5", "r", locking=False)
nA = f["actin/center"].shape[1]; nM = f["myosin/center"].shape[1]
for t in (400,):
    hs = f["state/am_head_state_current"][t]; hb = f["state/am_head_bound_actin_current"][t].astype(int)
    mi, hi = np.nonzero(hs == 2); pc = collections.Counter(zip(hb[mi, hi].tolist(), mi.tolist()))
    print("pairs by strong-head count:", sorted(collections.Counter(min(v, 8) for v in pc.values()).items()))
    for k in (1, 2, 3, 4):
        e = [(a, m + nA) for (a, m), v in pc.items() if v >= k]
        r = comps(nA + nM, e); C = collections.Counter(r)
        cl = []
        for root, n in C.items():
            mem = np.flatnonzero(r == root); na = int((mem < nA).sum()); nm = n - na
            if na and nm: cl.append((nm, na))
        cl.sort(reverse=True)
        print("k>=%d: %d clusters; myosins in clusters %d; top (myosins,actins) %s; by myosin count %s"
              % (k, len(cl), sum(c[0] for c in cl), cl[:6], dict(sorted(collections.Counter(min(c[0], 12) for c in cl).items()))))
