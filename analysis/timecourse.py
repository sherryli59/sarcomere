"""Time courses for all five arms, each ending at its own latest committed frame.

A side-by-side comparison at "latest frame" would hide the fact that the arms
are at 340-5000 ms. A time course makes the differing progress explicit and
lets trends be read instead of incommensurate endpoints.

band_order is VISUALIZATION_PIPELINE.md verbatim. Z order = max over the two
actin-tip populations. Both are POSITIONAL; neither is nematic S.
"""
import os, glob, json, collections
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
import numpy as np, h5py
NPTS = 200
def band_order(x, box_x):
    wrapped = np.asarray(x, dtype=float) % box_x
    h = np.arange(1, 5, dtype=float)
    ph = np.exp(2j*np.pi*wrapped[..., None]*h[None, None, :]/box_x).mean(axis=1)
    return np.max(np.abs(ph), axis=1)
S = 20260805
RUNS = [("A  hph9/spokes3 k=0",   f"/scratch/users/xsli/sarcomere/barbed_capsule_normal_full_20261002/data/A1600_M320_Lx6_intermediate_sig0.40_kaa3_seed{S}_5000ms_capsule_normal.h5"),
        ("B  hph9/spokes1 k=0",   f"/scratch/users/xsli/sarcomere/hph_ladder_20261005/data/A1600_M320_Lx6_intermediate_sig0.40_kaa3_seed{S}_5000ms_spokes1.h5"),
        ("C  hph27 k=0",          f"/scratch/users/xsli/sarcomere/hph_ladder_20261005/data/A1600_M320_Lx6_intermediate_sig0.40_kaa3_seed{S}_5000ms_hph27.h5"),
        ("A9 hph9 k=30",          f"/scratch/users/xsli/sarcomere/amalign_20261005/data/A1600_M320_Lx6_intermediate_sig0.40_kaa3_seed{S}_5000ms_align9.h5"),
        ("C27 hph27 k=30",        f"/scratch/users/xsli/sarcomere/amalign_20261005/data/A1600_M320_Lx6_intermediate_sig0.40_kaa3_seed{S}_5000ms_align27.h5")]
out = {"seed": S, "series": {}, "meta": {}}
for lab, p in RUNS:
    with h5py.File(p, "r", locking=False) as f:
        a = f["meta"].attrs
        dt = float(a["dt"]); L = float(a["actin_length"]); box = float(a["box_x"])
        cf = int(np.asarray(f["meta/committed_frames"]).reshape(-1)[-1])
        cf = min(cf, f["actin/direction"].shape[0], f["state/current_step"].shape[0],
                 f["state/am_head_state_current"].shape[0])
        st = np.asarray(f["state/current_step"][:cf]).reshape(-1)
        stride = max(1, cf//NPTS)
        idx = list(range(0, cf, stride))
        u = np.asarray(f["actin/direction"][idx], float); u /= np.linalg.norm(u, axis=2, keepdims=True)
        w = np.asarray(f["myosin/direction"][idx], float); w /= np.linalg.norm(w, axis=2, keepdims=True)
        ac = np.asarray(f["actin/center"][idx], float); mc = np.asarray(f["myosin/center"][idx], float)
        hs = np.asarray(f["state/am_head_state_current"][idx])
        links = np.asarray(f["state/aa_pairs_current_count"][idx]).reshape(len(idx), -1)[:, -1]
        Q = 1.5*np.einsum("fia,fib->fab", u, u)/u.shape[1] - 0.5*np.eye(3)[None]
        Sv = np.linalg.eigvalsh(Q)[:, -1]
        M = band_order(mc[:, :, 0], box)
        Z = np.maximum(band_order(ac[:, :, 0] - 0.5*L*u[:, :, 0], box),
                       band_order(ac[:, :, 0] + 0.5*L*u[:, :, 0], box))
        strong = (hs == 2).sum(axis=(1, 2)).astype(float)
        # status-2 ACTN2 from the live link identity
        ident = np.asarray(f["catch_bond/local_link_identity"])
        off = np.asarray(f["catch_bond/local_link_frame_offsets"])
        s2 = []
        for k in idx:
            j = int(np.searchsorted(off[:, 0], st[k], "left"))
            if j < off.shape[0] and int(off[j, 0]) == int(st[k]):
                s0 = int(off[j, 1]); n = int(off[j, 2])
                s2.append(float((ident[s0:s0+n, 4].astype(int) == 2).sum()) if n > 0 else 0.0)
            else:
                s2.append(float("nan"))
        # AM angle to nearest myosin, mean over actins
        th = []
        for t in range(len(idx)):
            d = ac[t][:, None, :] - mc[t][None, :, :]
            d -= np.array([box, float(a["box_y"]), float(a["box_z"])])*np.rint(
                d/np.array([box, float(a["box_y"]), float(a["box_z"])]))
            near = np.linalg.norm(d, axis=2).argmin(axis=1)
            th.append(float(np.degrees(np.arccos(np.clip(np.abs(
                np.einsum("ak,ak->a", u[t], w[t][near])), 0, 1))).mean()))
        out["series"][lab] = dict(
            t=[round(float(st[k])*dt*1e3, 3) for k in idx],
            S=[round(float(v), 5) for v in Sv],
            links=[round(float(v), 1) for v in links],
            s2=[None if np.isnan(v) else round(float(v), 1) for v in s2],
            M=[round(float(v), 5) for v in M],
            Z=[round(float(v), 5) for v in Z],
            strong=[round(float(v), 1) for v in strong],
            theta=[round(float(v), 3) for v in th])
        out["meta"][lab] = dict(frames=cf, t_max=round(float(st[cf-1])*dt*1e3, 1),
                                hph=int(a["am_heads_per_half"]),
                                spokes=int(float(a.get("am_spokes_per_crown", 3))),
                                align_k=float(a["am_alignment_k"]), n=len(idx))
        print("  %-22s frames=%-5d t_max=%7.1f ms  points=%-4d hph=%-3d spokes=%d k=%.0f"
              % (lab, cf, out["meta"][lab]["t_max"], len(idx), a["am_heads_per_half"],
                 out["meta"][lab]["spokes"], a["am_alignment_k"]))
dest = "/scratch/users/xsli/sarcomere/amalign_20261005/analysis/timecourse.json"
json.dump(out, open(dest, "w"))
print("wrote", dest)
