"""Is the high actin order under am_alignment_k REAL assembly, or just actins
being slaved to a myosin field that the initial condition already supplied?

initial_structure=myosin_x_noise starts myosins near the x axis with
myosin_direction_noise=0.40. If the alignment term simply slaves actins to the
myosin director field, actin S should converge to MYOSIN S -- and myosin S is
an initial-condition property, not an assembly result. Compare actin S and
myosin S at frame 0 and at 80 ms, in all four cells.
"""
import os, glob
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
import numpy as np, h5py
FULL = "/scratch/users/xsli/sarcomere/barbed_capsule_normal_full_20261002/data"
LAD  = "/scratch/users/xsli/sarcomere/hph_ladder_20261005/data"
AL   = "/scratch/users/xsli/sarcomere/amalign_20261005/data"
ARMS = [("A  hph9  k=0 ", f"{FULL}/*kaa3_seed*_5000ms_capsule_normal.h5"),
        ("A9 hph9  k=30", f"{AL}/*kaa3_seed*_5000ms_align9.h5"),
        ("C  hph27 k=0 ", f"{LAD}/*kaa3_seed*_5000ms_hph27.h5"),
        ("C27 hph27 k=30", f"{AL}/*kaa3_seed*_5000ms_align27.h5")]
def S(v):
    v = v/np.linalg.norm(v, axis=1, keepdims=True)
    Q = 1.5*np.einsum("ia,ib->ab", v, v)/v.shape[0] - 0.5*np.eye(3)
    return float(np.linalg.eigvalsh(Q).max())
def director(v):
    v = v/np.linalg.norm(v, axis=1, keepdims=True)
    Q = 1.5*np.einsum("ia,ib->ab", v, v)/v.shape[0] - 0.5*np.eye(3)
    w, V = np.linalg.eigh(Q)
    return V[:, -1]
for tag, pat in ARMS:
    rows = []
    for p in sorted(glob.glob(pat)):
        with h5py.File(p, "r", locking=False) as f:
            cf = int(np.asarray(f["meta/committed_frames"]).reshape(-1)[-1])
            cf = min(cf, f["actin/direction"].shape[0], f["state/current_step"].shape[0])
            st = np.asarray(f["state/current_step"][:cf]).reshape(-1)
            k = int(np.searchsorted(st, 80*100, "right")) - 1
            if k < 1: continue
            a0 = np.asarray(f["actin/direction"][0], float); m0 = np.asarray(f["myosin/direction"][0], float)
            ak = np.asarray(f["actin/direction"][k], float); mk = np.asarray(f["myosin/direction"][k], float)
            da, dm = director(ak), director(mk)
            rows.append((S(a0), S(m0), S(ak), S(mk),
                         float(np.degrees(np.arccos(min(1.0, abs(float(np.dot(da, dm))))))) ))
    if not rows: print("\n%s NO DATA" % tag); continue
    r = np.array(rows)
    print("\n=== %s  n=%d" % (tag, len(rows)))
    print("    frame 0 :  actin S %.4f   myosin S %.4f" % (r[:,0].mean(), r[:,1].mean()))
    print("    80 ms   :  actin S %.4f   myosin S %.4f   actin-vs-myosin director angle %.2f deg"
          % (r[:,2].mean(), r[:,3].mean(), r[:,4].mean()))
    print("    actin S / myosin S at 80 ms = %.3f   (1.0 => actins fully slaved to the myosin field)"
          % (r[:,2].mean()/r[:,3].mean()))
