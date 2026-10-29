"""Spread (sd) of each source's instantaneous drive toward the partner, singular
frames excluded, and the rotation-cap event rate.  sd ~ sqrt(kappa)*kT for a
harmonic orientational well, so sd^2 (in kT) estimates the source's stiffness."""
import numpy as np, sys, h5py
D = "/scratch/users/xsli/sarcomere/gaphunt_20261023/torque_out"; C = "/scratch/users/xsli/sarcomere/gaphunt_20261023/data/"
SRC = ["am_head", "RAD", "actn2_align", "actn2_lever", "actn2_axial", "aa_lateral", "am_steric", "TOT"]
FILES = dict(hph27="hph27m24_seed20260805.h5", hph81="hph81m24_seed20260805.h5", tq30_seed20260805="tq30_seed20260805.h5")
for tag in sys.argv[1:]:
    z = dict(np.load("%s/torque_%s.npz" % (D, tag)))
    ok = (np.abs(z["aa_steric"]) < 5000) & (np.abs(z["TOT"]) < 5000)
    print("%-20s" % tag + "  ".join("%s sd %.0f (k~%.0e)" % (s, z[s][ok].std(), z[s][ok].var()) for s in SRC))
    f = h5py.File(C + FILES[tag], "r", locking=False); n = int(f["meta/committed_frames"][0])
    st = f["state/current_step"][:n, 0]; t0 = max(0, n - 101)
    rc = f["actin/cumulative_rotation_cap_events"][n - 1, :, 0] - f["actin/cumulative_rotation_cap_events"][t0, :, 0]
    tc = f["actin/cumulative_torque_cap_events"][n - 1, :, 0] - f["actin/cumulative_torque_cap_events"][t0, :, 0]
    steps = st[n - 1] - st[t0]
    print("   rotation-cap events per actin-step %.2e, torque-cap %.2e  (over %d steps, %d actins)" % (rc.sum() / steps / rc.size, tc.sum() / steps / tc.size, steps, rc.size))
