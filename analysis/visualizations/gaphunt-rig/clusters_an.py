"""Connected components of the actin-myosin graph, edges = STRONG heads (state 2)
(and, for reference, any bound head).  Nodes: actins 0..nA-1, myosins nA..nA+nM-1.
Prints the component-size distribution at a few frames of the resume arms."""
import h5py, numpy as np, os, sys, collections
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
CAMP = "/scratch/users/xsli/sarcomere/gaphunt_20261023/data"
def comps(n, edges):
    p = np.arange(n)
    def find(x):
        while p[x] != x: p[x] = p[p[x]]; x = p[x]
        return x
    for a, b in edges:
        ra, rb = find(a), find(b)
        if ra != rb: p[rb] = ra
    return np.array([find(i) for i in range(n)])
for arm in sys.argv[1:] or ["r27ctrl", "r27nolat"]:
    f = h5py.File("%s/%s_seed20260805.h5" % (CAMP, arm), "r", locking=False)
    nA = f["actin/center"].shape[1]; nM = f["myosin/center"].shape[1]
    for t in (400, 440, 480):
        hs = f["state/am_head_state_current"][t]; hb = f["state/am_head_bound_actin_current"][t].astype(int)
        for lab, mask in (("strong", hs == 2), ("any bound", (hs == 1) | (hs == 2))):
            mi, hi = np.nonzero(mask); ai = hb[mi, hi]
            e = set(zip(ai.tolist(), (mi + nA).tolist()))
            r = comps(nA + nM, e); C = collections.Counter(r)
            sizes = []
            for root, n in C.items():
                mem = np.flatnonzero(r == root); na = int((mem < nA).sum()); nm = n - na
                if nm >= 1 and na >= 1: sizes.append((n, na, nm))
            sizes.sort(reverse=True)
            free_a = sum(1 for root, n in C.items() if n == 1 and root < nA); free_m = sum(1 for root, n in C.items() if n == 1 and root >= nA)
            hist = collections.Counter(min(s[2], 10) for s in sizes)
            print("%s t=%.2fs %-9s: %d clusters; largest (actins,myosins) %s; unattached actins %d myosins %d; clusters by myosin count %s"
                  % (arm, t * 0.005, lab, len(sizes), [(s[1], s[2]) for s in sizes[:5]], free_a, free_m,
                     dict(sorted(hist.items()))), flush=True)
