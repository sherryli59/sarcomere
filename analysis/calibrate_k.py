"""What value of am_alignment_k is actually SMALL?

The term is compute_am_nematic_flat_bottom_tether with
  stiffness = interaction_ramp * am_alignment_k * angular_occupancy
  angular_occupancy = min(1, bound_head_count / am_tether_saturation_heads)
  |drive| = stiffness * (sin(theta) - sin(slack)) * |cos(theta)|
so the per-pair torque is k * occupancy * excess * |dot|, in pN.um.

Calibrate against the torque that ACTUALLY aligns filaments in THESE runs:
tb_actn2_align_{x,y,z}, measured per actin. Then "small" means comparable to it.
"""
import os, glob, collections
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
import numpy as np, h5py
T0, T1 = 1000, 1500
SLACK = np.sin(np.radians(2.0))
FULL = "/scratch/users/xsli/sarcomere/barbed_capsule_normal_full_20261002/data"
LAD  = "/scratch/users/xsli/sarcomere/hph_ladder_20261005/data"
ARMS = [("A hph9 ", f"{FULL}/*kaa3_seed20260805*_5000ms_capsule_normal.h5"),
        ("C hph27", f"{LAD}/*kaa3_seed20260805*_5000ms_hph27.h5")]
for tag, pat in ARMS:
    for p in sorted(glob.glob(pat)):
        with h5py.File(p, "r", locking=False) as f:
            a = f["meta"].attrs
            sath = int(a["am_tether_saturation_heads"])
            cf = int(np.asarray(f["meta/committed_frames"]).reshape(-1)[-1])
            cf = min(cf, f["actin/direction"].shape[0], f["state/current_step"].shape[0])
            st = np.asarray(f["state/current_step"][:cf]).reshape(-1)
            i0 = int(np.searchsorted(st, T0*100, "left")); i1 = min(int(np.searchsorted(st, T1*100, "right")), cf)
            fr = list(range(i0, i1, 25))
            u = np.asarray(f["actin/direction"][fr], float); u /= np.linalg.norm(u, axis=2, keepdims=True)
            w = np.asarray(f["myosin/direction"][fr], float); w /= np.linalg.norm(w, axis=2, keepdims=True)
            hs = np.asarray(f["state/am_head_state_current"][fr])
            ha = np.asarray(f["state/am_head_bound_actin_current"][fr]).astype(int)
            # measured ACTN2 aligning torque magnitude per actin
            tb = np.stack([np.asarray(f["actin/tb_actn2_align_%s" % c][fr], float) for c in "xyz"], -1)
            tbm = np.linalg.norm(tb, axis=-1)
            nz = tbm[tbm > 0]
            occ = []; pred = collections.defaultdict(list)
            for t in range(len(fr)):
                bound = (hs[t] >= 1) & (ha[t] >= 0)
                cnt = collections.Counter()
                mi, hh = np.nonzero(bound)
                for m_, h_ in zip(mi.tolist(), hh.tolist()):
                    cnt[(int(ha[t][m_][h_]), m_)] += 1
                for (ai, mj), nb in cnt.items():
                    o = min(1.0, nb/max(1, sath))
                    dot = abs(float(np.dot(u[t][ai], w[t][mj])))
                    sine = np.sqrt(max(0.0, 1-dot*dot))
                    if sine <= SLACK: continue
                    base = o*(sine-SLACK)*dot          # torque per unit k
                    occ.append(o); pred[1].append(base)
            occ = np.array(occ); base = np.array(pred[1])
            print("\n=== %s  (am_tether_saturation_heads=%d, %d engaged pair-frames)" % (tag, sath, base.size))
            print("    bound heads/pair -> angular_occupancy: mean %.4f  median %.4f" % (occ.mean(), np.median(occ)))
            print("    MEASURED ACTN2 aligning torque |tb_actn2_align|: mean %.4f  median %.4f pN.um (nonzero n=%d)"
                  % (nz.mean(), np.median(nz), nz.size))
            print("    predicted |am_align| per pair at each k (pN.um), and ratio to the ACTN2 torque:")
            for k in (1, 3, 10, 30, 120):
                v = k*base
                print("      k=%-4d  mean %8.4f   = %7.2fx the ACTN2 aligner" % (k, v.mean(), v.mean()/nz.mean()))
