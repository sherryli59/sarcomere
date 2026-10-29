"""Why do things align in the two-body probe but not production?

HYPOTHESIS: the actins ARE well aligned -- to each other. What is not aligned is
the MYOSIN field. theta_AM then measures the mismatch between a well-ordered
actin director and individually scattered myosin axes, and no amount of head
couple can fix that, because the target itself is scattered.

The probe cannot show this: with ONE myosin, "aligned to my myosin" and
"aligned to the myosin field" are identical.

TEST: at 5000 ms compare
  (a) actin S and myosin S
  (b) the angle from each MYOSIN axis to the global ACTIN director
  (c) the measured theta_AM for strong-head pairs
If (b) ~ (c), myosin scatter sets the floor on theta_AM and the heads are
irrelevant to it.
"""
import os, glob
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
import numpy as np, h5py
D = "/scratch/users/xsli/sarcomere"
ARMS = [("A hph9/spk3  k=0", f"{D}/barbed_capsule_normal_full_20261002/data/*kaa3_seed*_5000ms_capsule_normal.h5"),
        ("B hph9/spk1  k=0", f"{D}/hph_ladder_20261005/data/*kaa3_seed*_5000ms_spokes1.h5"),
        ("C hph27/spk3 k=0", f"{D}/hph_ladder_20261005/data/*kaa3_seed*_5000ms_hph27.h5")]
def Sdir(v):
    v = v/np.linalg.norm(v, axis=1, keepdims=True)
    Q = 1.5*np.einsum("ia,ib->ab", v, v)/v.shape[0] - 0.5*np.eye(3)
    w, V = np.linalg.eigh(Q)
    return float(w[-1]), V[:, -1]
ang = lambda a, b: float(np.degrees(np.arccos(min(1.0, abs(float(np.dot(a, b)))))))
for tag, pat in ARMS:
    rows = []
    for p in sorted(glob.glob(pat)):
        with h5py.File(p, "r", locking=False) as f:
            cf = int(np.asarray(f["meta/committed_frames"]).reshape(-1)[-1])
            cf = min(cf, f["actin/direction"].shape[0], f["state/current_step"].shape[0],
                     f["state/am_head_state_current"].shape[0])
            st = np.asarray(f["state/current_step"][:cf]).reshape(-1)
            k = int(np.searchsorted(st, 5000*100, "right")) - 1
            if float(st[k]) < 5000*100*0.999: continue
            u = np.asarray(f["actin/direction"][k], float); u /= np.linalg.norm(u, axis=1, keepdims=True)
            w = np.asarray(f["myosin/direction"][k], float); w /= np.linalg.norm(w, axis=1, keepdims=True)
            Sa, da = Sdir(u); Sm, dm = Sdir(w)
            # angle from each myosin axis to the GLOBAL ACTIN director
            m_to_adir = np.degrees(np.arccos(np.clip(np.abs(w @ da), 0, 1)))
            # angle from each actin to the global actin director
            a_to_adir = np.degrees(np.arccos(np.clip(np.abs(u @ da), 0, 1)))
            # measured theta_AM over strong-head pairs
            hs = np.asarray(f["state/am_head_state_current"][k])
            ha = np.asarray(f["state/am_head_bound_actin_current"][k]).astype(int)
            mi, hh = np.nonzero((hs == 2) & (ha >= 0))
            pr = np.unique(np.stack([ha[mi, hh], mi], 1), axis=0)
            th = np.degrees(np.arccos(np.clip(np.abs(np.einsum("pk,pk->p", u[pr[:, 0]], w[pr[:, 1]])), 0, 1)))
            rows.append(dict(Sa=Sa, Sm=Sm, dirmis=ang(da, dm),
                             m_med=float(np.median(m_to_adir)), a_med=float(np.median(a_to_adir)),
                             th_med=float(np.median(th))))
    if not rows: print("\n=== %s  no seed at 5000 ms" % tag); continue
    g = lambda k: float(np.mean([r[k] for r in rows]))
    print("\n=== %s   n=%d seeds at 5000 ms" % (tag, len(rows)))
    print("    actin  S = %.4f    median actin-to-actin-director angle  = %5.2f deg" % (g("Sa"), g("a_med")))
    print("    myosin S = %.4f    median MYOSIN-to-actin-director angle = %5.2f deg" % (g("Sm"), g("m_med")))
    print("    actin director vs myosin director: %.2f deg" % g("dirmis"))
    print("    MEASURED theta_AM median (strong pairs)   = %5.2f deg" % g("th_med"))
    print("    -> myosin scatter explains %.0f%% of theta_AM" % (100*g("m_med")/g("th_med")))
