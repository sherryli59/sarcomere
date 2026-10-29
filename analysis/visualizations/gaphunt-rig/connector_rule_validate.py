"""Validate connector_rule_impact.py's replication.
(1) Old soft rule: the set of myosin tips (myosin, end) with at least one eligible link must equal
    the tips the model recorded with force > 0 (myosin/aa_zdisc_tip_force_{minus,plus}) in the same frame.
(2) The measured case pinned in status2_connector_gate_test.cpp: r27zt0np frame 475, link 179352,
    myosin 244 must be inside the 130 nm admission cylinder."""
import h5py, numpy as np, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
from connector_rule_impact import slab_hit, U
D = "/scratch/users/xsli/sarcomere/gaphunt_20261023/data/"
def tips_old(path, t):
    f = h5py.File(path, "r", locking=False); a = dict(f["meta"].attrs)
    box = np.array([a["box_x"], a["box_y"], a["box_z"]], float); L = float(a["actin_length"]); Lm = 1.5
    mi = lambda v: v - box * np.round(v / box)
    off = f["catch_bond/local_link_frame_offsets"]; o, c = int(off[t][1]), int(off[t][2])
    idr = f["catch_bond/local_link_identity"][o:o + c]; co = f["catch_bond/local_link_coordinates"][o:o + c]
    i, j = idr[:, 2].astype(int), idr[:, 3].astype(int)
    C = f["actin/center"][t]; Ua = U(f["actin/direction"][t]); Cm = f["myosin/center"][t]; Um = U(f["myosin/direction"][t])
    p = C[i] + (co[:, 0] * L)[:, None] * Ua[i]; q = C[j] + (co[:, 1] * L)[:, None] * Ua[j]
    mid = p + 0.5 * mi(q - p)
    mine = set()
    for m in range(len(Cm)):
        rel = mi(mid - Cm[m]); z = rel @ Um[m]; rho = np.linalg.norm(rel - z[:, None] * Um[m], axis=1)
        for sign, name in ((-1, "minus"), (1, "plus")):
            if np.any((sign * z >= 0) & (sign * z - 0.75 < 0.05) & (rho < 0.13)): mine.add((m, name))
    rec = set((m, s) for s in ("minus", "plus") for m in np.flatnonzero(f["myosin/aa_zdisc_tip_force_" + s][t][:, 0] > 0))
    return mine, rec, f, idr, co, C, Ua, Cm, Um, box, L
for fn, t in (("hph27m24_seed20260805.h5", 400), ("r27zt0_seed20260805.h5", 480)):
    mine, rec = tips_old(D + fn, t)[:2]
    print("%s frame %d: tips with an eligible link (old midpoint rule) %d | recorded tips with force > 0 %d | both %d, mine only %d, recorded only %d"
          % (fn, t, len(mine), len(rec), len(mine & rec), len(mine - rec), len(rec - mine)))
mine, rec, f, idr, co, C, Ua, Cm, Um, box, L = tips_old(D + "r27zt0np_seed20260805.h5", 475)
mi = lambda v: v - box * np.round(v / box)
k = np.flatnonzero(idr[:, 1] == 179352)
print("r27zt0np frame 475 tips: mine %d, recorded %d, both %d" % (len(mine), len(rec), len(mine & rec)))
if len(k):
    k = k[0]; i, j = int(idr[k, 2]), int(idr[k, 3])
    p = C[i] + co[k, 0] * L * Ua[i]; q = C[j] + co[k, 1] * L * Ua[j]; conn = mi(q - p)
    rel_mid = mi(p + 0.5 * conn - Cm[244]); first = rel_mid - 0.5 * conn
    hit, rho = slab_hit(first[None], conn[None], Um[244][None], 0.75, 0.130)
    print("link 179352 (actins %d, %d; test says 1316, 233) vs myosin 244: inside 130 nm cylinder = %s, min radius in slab %.1f nm" % (i, j, bool(hit[0]), 1000 * rho[0]))
else:
    print("link 179352 not present in frame 475")
