"""Myosins ARE torqued by heads (sarcomere.cpp:4315) and are 2.5x more
rotationally mobile than actins (crowding 64 vs 160). So why does myosin S stay
at its initial value?

Two possibilities:
  FROZEN      -- they barely rotate at all. Then the limit is effective mobility
                 or torque magnitude.
  DIFFUSING   -- they rotate plenty but incoherently, so S does not rise. Then
                 the limit is that the torque has no collective preferred
                 direction (no myosin-myosin coupling).

Measure the angular displacement of each filament from t=0, and the S(t)
course, for both species.
"""
import os, glob
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
import numpy as np, h5py
P = "/scratch/users/xsli/sarcomere/barbed_capsule_normal_full_20261002/data/*kaa3_seed*_5000ms_capsule_normal.h5"
def S(v):
    v = v/np.linalg.norm(v, axis=1, keepdims=True)
    Q = 1.5*np.einsum("ia,ib->ab", v, v)/v.shape[0] - 0.5*np.eye(3)
    return float(np.linalg.eigvalsh(Q).max())
rows = []
for p in sorted(glob.glob(P)):
    with h5py.File(p, "r", locking=False) as f:
        cf = int(np.asarray(f["meta/committed_frames"]).reshape(-1)[-1])
        cf = min(cf, f["actin/direction"].shape[0], f["state/current_step"].shape[0])
        st = np.asarray(f["state/current_step"][:cf]).reshape(-1)
        if float(st[cf-1]) < 5000*100*0.999: continue
        idx = [0] + [int(np.searchsorted(st, t*100, "right"))-1 for t in (100, 500, 1500, 3000, 5000)]
        u = np.asarray(f["actin/direction"][idx], float); u /= np.linalg.norm(u, axis=2, keepdims=True)
        w = np.asarray(f["myosin/direction"][idx], float); w /= np.linalg.norm(w, axis=2, keepdims=True)
        ra = [float(np.median(np.degrees(np.arccos(np.clip(np.abs(np.einsum("ik,ik->i", u[k], u[0])), 0, 1))))) for k in range(len(idx))]
        rm = [float(np.median(np.degrees(np.arccos(np.clip(np.abs(np.einsum("ik,ik->i", w[k], w[0])), 0, 1))))) for k in range(len(idx))]
        Sa = [S(u[k]) for k in range(len(idx))]
        Sm = [S(w[k]) for k in range(len(idx))]
        rows.append((ra, rm, Sa, Sm, [float(st[i])*1e-5*1e3 for i in idx]))
if rows:
    ra = np.mean([r[0] for r in rows], 0); rm = np.mean([r[1] for r in rows], 0)
    Sa = np.mean([r[2] for r in rows], 0); Sm = np.mean([r[3] for r in rows], 0)
    ts = rows[0][4]
    print("baseline arm A, n=%d seeds\n" % len(rows))
    print("   t (ms)   median rotation from t=0      nematic S")
    print("            actin      myosin        actin    myosin")
    for k in range(len(ts)):
        print("   %6.0f   %6.2f deg  %6.2f deg     %.4f   %.4f" % (ts[k], ra[k], rm[k], Sa[k], Sm[k]))
    print("\n   actin  rotated %.1f deg and S went %.3f -> %.3f  (ORDERED)" % (ra[-1], Sa[0], Sa[-1]))
    print("   myosin rotated %.1f deg and S went %.3f -> %.3f  (%s)"
          % (rm[-1], Sm[0], Sm[-1], "DIFFUSING, not ordering" if rm[-1] > 10 else "FROZEN"))
