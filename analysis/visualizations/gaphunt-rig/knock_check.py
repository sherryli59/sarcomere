"""Physics checks that each knockout took effect, 2.0 -> 2.2 s, vs control:
noise0  - frame-to-frame actin and myosin rotation (rad per 5 ms) drops without Brownian kicks
zgate   - Z-disc tip contact events / torque (tb_actn2_contact) on unlinked tips
titin0  - titin force ~ 0
motor0  - strong-head axial strain energy contribution ~ 0 (am_motor_k = 0)"""
import h5py, numpy as np, os
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
D = "/scratch/users/xsli/sarcomere/gaphunt_20261023/data/"
KT = 1.380649e-23 * 310 / 1e-18
for arm in ("r27ctrl", "r27noise0", "r27zgate", "r27titin0", "r27motor0"):
    f = h5py.File(D + arm + "_seed20260805.h5", "r", locking=False)
    ts = range(402, 441, 4)
    rot_a = []; rot_m = []; tc = []; tf = []; se = []
    for t in ts:
        u0 = f["actin/direction"][t - 1]; u1 = f["actin/direction"][t]
        m0 = f["myosin/direction"][t - 1]; m1 = f["myosin/direction"][t]
        n = lambda x: x / np.linalg.norm(x, axis=1, keepdims=True)
        rot_a.append(np.median(np.arccos(np.clip(np.abs(np.sum(n(u0) * n(u1), 1)), 0, 1))))
        rot_m.append(np.median(np.arccos(np.clip(np.abs(np.sum(n(m0) * n(m1), 1)), 0, 1))))
        c = np.linalg.norm(np.stack([f["actin/tb_actn2_contact_%s" % k][t, :, 0] for k in "xyz"], 1), axis=1) / KT
        tc.append(c.mean())
        if "state/titin_force_current" in f: tf.append(np.abs(f["state/titin_force_current"][t]).mean())
        hs = f["state/am_head_state_current"][t]; e = f["state/am_head_strain_energy_current"][t]
        se.append(e[hs == 2].mean() if (hs == 2).any() else 0)
    ev = f["actin/cumulative_aa_zdisc_tip_projection_events"] if "actin/cumulative_aa_zdisc_tip_projection_events" in f else None
    evr = float((ev[440, :, 0] - ev[400, :, 0]).sum()) if ev is not None else float("nan")
    print("%-10s actin rot/frame %.4f rad  myosin rot/frame %.4f rad | Z-tip contact torque %.1f kT/actin, tip projection events 2.0-2.2 s %.0f | titin |F| %.3g pN | strong-head strain energy %.4f pN.um"
          % (arm, np.mean(rot_a), np.mean(rot_m), np.mean(tc), evr, np.mean(tf) if tf else float("nan"), np.mean(se)))
