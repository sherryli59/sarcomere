"""If the production state is resumed with the 2026-10-28 binary and the hard projection ON (the
production default), which links does the projection act on at once under the NEW full-connector
tip geometry, and by how much?  Mirrors the new zdisc_myosin_tip_gap: per link, myosin and tip,
clip the connector to the tip's own half (0 <= z < L/2 + onset along the outward axis), take the
point of minimum radial distance, eligible if rho < 0.130; the projection fires if the axial gap
there (z - L/2) < 0 and displaces the pair by -gap along the outward axis.  Nearest periodic image
(reach < half of every box side).  Old geometry (midpoint) is held at gap >= 0 by the projection."""
import h5py, numpy as np, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
D = "/scratch/users/xsli/sarcomere/gaphunt_20261023/data/"
U = lambda d: d / np.linalg.norm(d, axis=-1, keepdims=True)
for fn, t in (("hph27m24_seed20260805.h5", 400), ("r27zt0_seed20260805.h5", 480)):
    f = h5py.File(D + fn, "r", locking=False); a = dict(f["meta"].attrs)
    box = np.array([a["box_x"], a["box_y"], a["box_z"]], float); L = float(a["actin_length"]); h = 0.75; onset = 0.05; R = 0.13
    mi = lambda v: v - box * np.round(v / box)
    off = f["catch_bond/local_link_frame_offsets"]; o, c = int(off[t][1]), int(off[t][2])
    idr = f["catch_bond/local_link_identity"][o:o + c]; co = f["catch_bond/local_link_coordinates"][o:o + c]
    st = idr[:, 4].astype(int); i, j = idr[:, 2].astype(int), idr[:, 3].astype(int)
    C = f["actin/center"][t]; Ua = U(f["actin/direction"][t]); Cm = f["myosin/center"][t]; Um = U(f["myosin/direction"][t])
    p = C[i] + (co[:, 0] * L)[:, None] * Ua[i]; q = C[j] + (co[:, 1] * L)[:, None] * Ua[j]
    conn = mi(q - p); mid = p + 0.5 * conn
    deepest = np.zeros(len(p)); nproj = np.zeros(len(p), int)
    for m in range(len(Cm)):
        rel = mi(mid - Cm[m]); first = rel - 0.5 * conn
        for sign in (-1.0, 1.0):
            out = sign * Um[m]
            z0 = first @ out; dz = conn @ out
            with np.errstate(divide="ignore", invalid="ignore"):
                e1 = -z0 / dz; e2 = (h + onset - z0) / dz
            lo = np.where(dz != 0, np.maximum(0, np.minimum(e1, e2)), 0.0); hi = np.where(dz != 0, np.minimum(1, np.maximum(e1, e2)), 1.0)
            ok = np.where(dz != 0, lo <= hi, (z0 >= 0) & (z0 < h + onset))
            rf = first - z0[:, None] * out; rd = conn - dz[:, None] * out; A2 = np.sum(rd * rd, 1)
            with np.errstate(divide="ignore", invalid="ignore"):
                ts = np.where(A2 > 0, -np.sum(rf * rd, 1) / A2, lo)
            tt = np.clip(ts, lo, hi); rho = np.linalg.norm(rf + tt[:, None] * rd, axis=1); gap = z0 + tt * dz - h
            fire = ok & (rho < R) & (gap < -1e-9)
            nproj += fire; deepest = np.where(fire, np.maximum(deepest, -gap), deepest)
    hit = nproj > 0; dpt = 1000 * deepest[hit]
    print("%s frame %d: links the projection would push at once under the new geometry: %d (mature %d, new %d) of %d"
          % (fn, t, hit.sum(), (hit & (st == 2)).sum(), (hit & (st == 1)).sum(), len(p)))
    if hit.any():
        print("   push per link (deepest violation): median %.0f nm, p90 %.0f nm, max %.0f nm; > 50 nm: %d links, > 120 nm: %d links"
              % (np.median(dpt), np.percentile(dpt, 90), dpt.max(), (dpt > 50).sum(), (dpt > 120).sum()))
