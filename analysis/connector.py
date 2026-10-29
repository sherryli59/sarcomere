"""Are arm C's ACTN2 links being mechanically RIPPED APART?

Reason 1 "geometric_contour" is a misnomer: catch_bond_kinetics.h:102 sets it
from `connector_length > rupture_length` only, i.e. mechanical overextension
past aa_link_rupture_length = 0.12 um. The contour-occlusion branches in
sarcomere_actn2.cpp are all gated on aa_arm_state_model, which production sets
to FALSE, so they never fire. So reason 1 = mechanical rupture.

Prediction if "more heads rip the crosslinks": arm C's live ACTN2 connectors
should be LONGER / closer to 0.12 um than A's.
"""
import os, glob
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
import numpy as np, h5py
T0, T1 = 1000, 1500
FULL = "/scratch/users/xsli/sarcomere/barbed_capsule_normal_full_20261002/data"
LAD  = "/scratch/users/xsli/sarcomere/hph_ladder_20261005/data"
ARMS = [("A hph9/spokes3 ", f"{FULL}/*kaa3_seed*_5000ms_capsule_normal.h5"),
        ("B hph9/spokes1 ", f"{LAD}/*kaa3_seed*_5000ms_spokes1.h5"),
        ("C hph27/spokes3", f"{LAD}/*kaa3_seed*_5000ms_hph27.h5")]
for tag, pat in ARMS:
    allL = []; rupt = None
    for p in sorted(glob.glob(pat)):
        with h5py.File(p, "r", locking=False) as f:
            a = f["meta"].attrs
            rupt = float(a["aa_link_rupture_length"])
            cols = a["aa_local_link_mechanics_v2_columns"].decode().split(",")
            ocols = a["aa_local_link_mechanics_v2_frame_offsets_columns"].decode().split(",")
            jlen = cols.index("length")
            off = np.asarray(f["catch_bond/local_link_mechanics_v2_frame_offsets"])
            oi = {c: i for i, c in enumerate(ocols)}
            cf = int(np.asarray(f["meta/committed_frames"]).reshape(-1)[-1])
            st = np.asarray(f["state/current_step"][:cf]).reshape(-1)
            osteps = off[:, oi["step"]]
            M = f["catch_bond/local_link_mechanics_v2"]
            i0 = int(np.searchsorted(st, T0*100, "left")); i1 = min(int(np.searchsorted(st, T1*100, "right")), cf)
            for fr in range(i0, i1, 20):
                j = int(np.searchsorted(osteps, st[fr], "left"))
                if j >= off.shape[0] or int(off[j, oi["step"]]) != int(st[fr]): continue
                s0 = int(off[j, oi["mechanics_v2_start"]]); n = int(off[j, oi["count"]])
                if n <= 0: continue
                allL.append(np.asarray(M[s0:s0+n, jlen], float))
    L = np.concatenate(allL)
    print("\n=== %s   (%d live link-instances, rupture length %.3f um)" % (tag, L.size, rupt))
    print("    connector length: mean %.4f  median %.4f  p90 %.4f  p99 %.4f  max %.4f um"
          % (L.mean(), np.median(L), np.percentile(L, 90), np.percentile(L, 99), L.max()))
    for thr in (0.08, 0.10, 0.11, 0.115):
        print("    fraction above %.3f um (%.0f%% of rupture): %6.2f%%" % (thr, 100*thr/rupt, 100*(L > thr).mean()))
    print("    mean strain vs rupture: %.1f%%" % (100*L.mean()/rupt))
