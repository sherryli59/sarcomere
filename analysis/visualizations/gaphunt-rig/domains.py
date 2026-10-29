"""Local domains vs global order, production r27ctrl at 2.0 / 2.4 s (and r27nolat_titin0 at 2.4 s).
(a) nematic order of the filaments within radius R of each filament (same species), median over
    filaments, vs the global S; median angle between local director and global director.
(b) sarcomeric order per y-z column (cs x cs um, full x): band order of myosin centres (M) and of
    barbed ends (Z) within the column (columns with >= 12 members), vs the global value, and the
    spread of the band phase (position) across columns."""
import h5py, numpy as np, os, sys
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
D = "/scratch/users/xsli/sarcomere/gaphunt_20261023/data/"
U = lambda d: d / np.linalg.norm(d, axis=-1, keepdims=True)
def S_dir(d):
    Q = 1.5 * (d[:, :, None] * d[:, None, :]).mean(0) - 0.5 * np.eye(3); w, v = np.linalg.eigh(Q); return w[-1], v[:, -1]
def band(x, L, k=None):
    w = np.asarray(x) % L; h = np.arange(1, 5) if k is None else np.array([k])
    z = np.exp(2j * np.pi * w[:, None] * h[None, :] / L).mean(0); i = int(np.argmax(np.abs(z))); return float(np.abs(z[i])), int(h[i]), float(np.angle(z[i]))
for tag, fn, t in (("control 2.0 s", "r27ctrl", 400), ("control 2.4 s", "r27ctrl", 480), ("lateral+titin off 2.4 s", "r27nolat_titin0", 480)):
    f = h5py.File(D + fn + "_seed20260805.h5", "r", locking=False); a = dict(f["meta"].attrs)
    box = np.array([a["box_x"], a["box_y"], a["box_z"]]); La = float(a["actin_length"])
    print("== %s" % tag)
    for sp in ("actin", "myosin"):
        c = f[sp + "/center"][t]; d = U(f[sp + "/direction"][t]); Sg, ng = S_dir(d)
        dc = c[:, None, :] - c[None, :, :]; dc -= box * np.round(dc / box); r = np.linalg.norm(dc, axis=-1)
        row = []
        for R in (0.25, 0.5, 1.0):
            Ss, ang = [], []
            for i in range(len(c)):
                nb = np.flatnonzero(r[i] < R)
                if len(nb) < 5: continue
                s, n = S_dir(d[nb]); Ss.append(s); ang.append(np.degrees(np.arccos(min(1, abs(n @ ng)))))
            row.append("R=%.2f um: local S %.2f, local director off global %.1f deg (n=%d)" % (R, np.median(Ss), np.median(ang), len(Ss)))
        print("   %-6s global S %.3f | %s" % (sp, Sg, " | ".join(row)))
    cm = f["myosin/center"][t]; ca = f["actin/center"][t]; ua = U(f["actin/direction"][t]); bz = ca - 0.5 * La * ua
    for lab, pts in (("M (myosin centres)", cm), ("Z (barbed ends)", bz)):
        g, kg, phg = band(pts[:, 0], box[0])
        for cs in (0.5, 1.0):
            iy = np.floor((pts[:, 1] % box[1]) / cs).astype(int); iz = np.floor((pts[:, 2] % box[2]) / cs).astype(int)
            loc, phs = [], []
            for yy in range(int(round(box[1] / cs))):
                for zz in range(int(round(box[2] / cs))):
                    m = (iy == yy) & (iz == zz)
                    if m.sum() < 12: continue
                    o, k, ph = band(pts[m, 0], box[0], kg); loc.append(o); phs.append(ph)
            phs = np.array(phs); R_ph = float(np.abs(np.exp(1j * phs).mean())) if len(phs) else np.nan
            spread_nm = float(np.degrees(np.angle(np.exp(1j * (phs - np.angle(np.exp(1j * phs).mean())))).std()) / 360 * box[0] / kg * 1000) if len(phs) else np.nan
            print("   %-18s global band order %.3f (harmonic %d, period %.2f um) | %.1f um columns: local order median %.3f (n=%d), phase agreement across columns %.2f, phase spread %.0f nm"
                  % (lab, g, kg, box[0] / kg, cs, np.median(loc), len(loc), R_ph, spread_nm))
