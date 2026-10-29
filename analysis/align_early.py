"""EARLY read of the am_alignment_k x head-count 2x2 at a matched 80 ms.

80 ms is ~1.6% of the run. Standing caution on this rig: status-2 maturation
does not begin until ~25-35 ms and a 20 ms window once certified two production
runs containing no registered structure. So S and link counts here are
DIRECTIONAL ONLY. What IS meaningful this early is whether the alignment torque
activated at production scale and whether it is moving the eligible-pair count,
which is the mechanism predicted to matter (0.842 of arm C's formation deficit).
"""
import os, glob, collections
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
import numpy as np, h5py
T_MS = 80
FULL = "/scratch/users/xsli/sarcomere/barbed_capsule_normal_full_20261002/data"
LAD  = "/scratch/users/xsli/sarcomere/hph_ladder_20261005/data"
AL   = "/scratch/users/xsli/sarcomere/amalign_20261005/data"
ARMS = [("A  hph9  k=0 ", f"{FULL}/*kaa3_seed*_5000ms_capsule_normal.h5"),
        ("A9 hph9  k=30", f"{AL}/*kaa3_seed*_5000ms_align9.h5"),
        ("C  hph27 k=0 ", f"{LAD}/*kaa3_seed*_5000ms_hph27.h5"),
        ("C27 hph27 k=30", f"{AL}/*kaa3_seed*_5000ms_align27.h5")]
out = {}
for tag, pat in ARMS:
    rows = []
    for p in sorted(glob.glob(pat)):
        try:
            with h5py.File(p, "r", locking=False) as f:
                a = f["meta"].attrs
                cut = float(a["aa_preferred_centerline"]) + float(a["aa_capture_margin"])
                maxdot = float(a["aa_association_max_direction_dot"])
                L = float(a["actin_length"])
                box = np.array([a["box_x"], a["box_y"], a["box_z"]], float)
                cf = int(np.asarray(f["meta/committed_frames"]).reshape(-1)[-1])
                cf = min(cf, f["actin/center"].shape[0], f["state/current_step"].shape[0])
                st = np.asarray(f["state/current_step"][:cf]).reshape(-1)
                k = int(np.searchsorted(st, T_MS*100, "right")) - 1
                if k < 1: continue
                u = np.asarray(f["actin/direction"][k:k+1], float); u /= np.linalg.norm(u, axis=2, keepdims=True)
                w = np.asarray(f["myosin/direction"][k:k+1], float); w /= np.linalg.norm(w, axis=2, keepdims=True)
                ac = np.asarray(f["actin/center"][k:k+1], float)
                hs = np.asarray(f["state/am_head_state_current"][k:k+1])
                ha = np.asarray(f["state/am_head_bound_actin_current"][k:k+1]).astype(int)
                tb = np.stack([np.asarray(f["actin/tb_am_align_" + c][k:k+1], float) for c in "xyz"], -1)
                m = np.linalg.norm(tb, axis=-1); nz = m[m > 0]
                Q = 1.5*np.einsum("ia,ib->ab", u[0], u[0])/u.shape[1] - 0.5*np.eye(3)
                S = float(np.linalg.eigvalsh(Q).max())
                links = float(np.asarray(f["state/aa_pairs_current_count"][k]).reshape(-1)[-1])
                barbed = ac - 0.5*L*u
                d = barbed[0][:, None, :] - barbed[0][None, :, :]
                d -= box*np.rint(d/box)
                dist = np.linalg.norm(d, axis=2); np.fill_diagonal(dist, 1e9)
                dot = np.einsum("ik,jk->ij", u[0], u[0])
                elig = int(((dist < cut) & (dot < maxdot)).sum())//2
                # AM angle for strong-head pairs
                th = []
                mi, hh = np.nonzero((hs[0] == 2) & (ha[0] >= 0))
                for m_, h_ in zip(mi.tolist(), hh.tolist()):
                    ai = int(ha[0][m_][h_])
                    th.append(np.degrees(np.arccos(min(1.0, abs(float(np.dot(u[0][ai], w[0][m_])))))))
                rows.append(dict(S=S, links=links, elig=elig, amal=float(nz.mean()) if nz.size else 0.0,
                                 nz=int(nz.size), th=float(np.median(th)) if th else float("nan"),
                                 step=int(st[k])))
        except Exception as e:
            print("  ERR", os.path.basename(p), e)
    if not rows: print("\n%s  NO DATA" % tag); continue
    g = lambda key: float(np.mean([r[key] for r in rows]))
    out[tag.split()[0]] = {k2: g(k2) for k2 in ("S", "links", "elig", "amal", "th")}
    print("\n=== %s  n=%d seeds, step=%d (~%d ms)" % (tag, len(rows), rows[0]["step"], T_MS))
    print("    |tb_am_align| mean %7.4f pN.um  (nonzero actins %d)" % (g("amal"), int(np.mean([r["nz"] for r in rows]))))
    print("    eligible antiparallel barbed pairs %7.1f" % g("elig"))
    print("    ACTN2 links %7.1f    median AM angle %6.2f deg    S %.4f" % (g("links"), g("th"), g("S")))
print("\n=== the prediction: does k=30 raise eligible pairs, and MORE at hph=27?")
if {"A", "A9", "C", "C27"} <= set(out):
    for lo, hi, lab in (("A", "A9", "hph9 "), ("C", "C27", "hph27")):
        for key, nm in (("elig", "eligible pairs"), ("th", "AM angle"), ("links", "ACTN2 links"), ("S", "S")):
            d = out[hi][key] - out[lo][key]
            print("  %s  %-15s %8.2f -> %8.2f   (%+.2f, %+.1f%%)"
                  % (lab, nm, out[lo][key], out[hi][key], d, 100*d/out[lo][key] if out[lo][key] else float("nan")))
        print()
