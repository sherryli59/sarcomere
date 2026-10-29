"""Resume A/B from the 2.0 s state of the production hph runs (frames 0-400 are
the copied source; 401+ are the arm).  Per frame from 1.9 s on: hph_an's
frame_stats (pair angle, S, heads/pair, ACTN2 links ...), plus the actin angle to
the myosin director, the torque-cap event rate per actin-step since the previous
frame, mean |aa_lateral| torque, and each source's mean drive toward the myosin
director (kT).  Frames 0-400 of every arm are bit-identical to the source by
construction (copied file), which the script asserts on frame 400."""
import h5py, numpy as np, json, os, sys
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hph_an import frame_stats, anchor_table
CAMP = "/scratch/users/xsli/sarcomere/gaphunt_20261023/data"
KT = 1.380649e-23 * 310 / 1e-18
ARMS = [("r27ctrl", "hph27m24", "control"), ("r27nolat", "hph27m24", "lateral repulsion off"),
        ("r27lat4", "hph27m24", "lateral cap 40 to 4 pN"), ("r27soff16", "hph27m24", "strong dwell x10"),
        ("r81ctrl", "hph81m24", "control"), ("r81nolat", "hph81m24", "lateral repulsion off"),
        ("r27noise0", "hph27m24", "no thermal noise"), ("r27zgate", "hph27m24", "Z-stop gated on ACTN2"),
        ("r27titin0", "hph27m24", "titin off"), ("r27motor0", "hph27m24", "no motor force"),
        ("r27nolat_titin0", "hph27m24", "lateral off + titin off"), ("r27kaa30", "hph27m24", "kappa_aa 3 -> 30"),
        ("r27zt0", "hph27m24", "titin off + Z-stop gated"), ("r27zt0np", "hph27m24", "titin off + Z-stop gated + hard projection off"),
        ("r27zt0c2", "hph27m24", "titin off + Z-stop gated + lifetime_coeff 2"),
        ("r27zt0npg130", "hph27m24", "projection off + full-connector status-2 gate 130 nm")]
ONLY = set(sys.argv[2].split(",")) if len(sys.argv) > 2 else None
SRC = ["am_head", "aa_lateral", "actn2_align", "actn2_lever", "actn2_axial", "am_steric"]
def director(m):
    m = m / np.linalg.norm(m, axis=1, keepdims=True); w, v = np.linalg.eigh((m[:, :, None] * m[:, None, :]).mean(0)); return v[:, -1]
out = {}
for arm, src, label in ARMS:
    p = "%s/%s_seed20260805.h5" % (CAMP, arm)
    if not os.path.exists(p) or (ONLY and arm not in ONLY): continue
    f = h5py.File(p, "r", locking=False); a = dict(f["meta"].attrs)
    n = int(min(f["meta/committed_frames"][0], f["actin/center"].shape[0]))
    if n <= 401: print(arm, "has not advanced past the copied frames (n=%d)" % n); continue
    g = h5py.File("%s/%s_seed20260805.h5" % (CAMP, src), "r", locking=False)
    assert np.array_equal(f["actin/center"][400], g["actin/center"][400]), arm + ": frame 400 differs from source"
    anc, crn = anchor_table(int(a["am_heads_per_half"]), float(a["myosin_length"]), int(float(a.get("am_spokes_per_crown", 3))))
    step = f["state/current_step"][:n, 0]; dt = float(a["dt"])
    rows = []
    for t in range(380, n, 2):
        r = frame_stats(f, t, a, anc, crn); r["t_ms"] = float(step[t] * dt * 1e3)
        u = f["actin/direction"][t]; u = u / np.linalg.norm(u, axis=1, keepdims=True); n0 = director(f["myosin/direction"][t])
        nm = np.sign(u @ n0)[:, None] * n0[None, :]; tg = nm - u * np.sum(u * nm, 1)[:, None]
        tg /= np.maximum(np.linalg.norm(tg, axis=1, keepdims=True), 1e-12)
        r["actin_to_director_median"] = float(np.median(np.degrees(np.arccos(np.clip(np.abs(u @ n0), 0, 1)))))
        for s in SRC:
            v = np.stack([f["actin/tb_%s_%s" % (s, c)][t, :, 0] for c in "xyz"], 1)
            pr = np.sum(v * tg, 1) / KT; r["drive_" + s] = float(pr[np.abs(pr) < 5000].mean())
        m_ = f["myosin/direction"][t]; m_ = m_ / np.linalg.norm(m_, axis=1, keepdims=True)
        nc = director(np.vstack([u, m_]))
        r["actin_scatter_median"] = float(np.median(np.degrees(np.arccos(np.clip(np.abs(u @ nc), 0, 1)))))
        r["myosin_scatter_median"] = float(np.median(np.degrees(np.arccos(np.clip(np.abs(m_ @ nc), 0, 1)))))
        na = director(u); r["director_gap"] = float(np.degrees(np.arccos(min(1.0, abs(float(na @ n0))))))
        hs_ = f["state/am_head_state_current"][t]; hb_ = f["state/am_head_bound_actin_current"][t].astype(int)
        mi_, hi_ = np.nonzero((hs_ == 1) | (hs_ == 2)); kk = np.unique(hb_[mi_, hi_].astype(np.int64) * 1000 + mi_)
        pa, pm = kk // 1000, kk % 1000
        r["pair_angle_shuffled_median"] = float(np.median(np.degrees(np.arccos(np.clip(np.abs(np.sum(u[pa] * m_[np.random.default_rng(t).permutation(pm)], 1)), 0, 1)))))
        dsum = sum(f["myosin/am_align_drive_%s" % h][t, :, 0] for h in ("minus", "plus")); csum = sum(f["myosin/am_align_head_count_%s" % h][t, :, 0] for h in ("minus", "plus"))
        r["myosin_drive_per_head"] = float(dsum.sum() / max(csum.sum(), 1) / KT)
        hph_ = int(a["am_heads_per_half"]); half_ = (hi_ >= hph_).astype(int); ai_ = hb_[mi_, hi_]
        Qm = np.zeros((len(m_), 2, 3, 3)); cm_ = np.zeros((len(m_), 2))
        np.add.at(Qm, (mi_, half_), u[ai_][:, :, None] * u[ai_][:, None, :]); np.add.at(cm_, (mi_, half_), 1)
        eng = cm_.sum(1) > 0; both = (cm_ > 0).all(1)
        lead = lambda Q: np.linalg.eigh(Q)[1][..., -1]
        avg = lead(Qm.sum(1)[eng] / cm_.sum(1)[eng][:, None, None])
        r["myosin_to_own_actins_median"] = float(np.median(np.degrees(np.arccos(np.clip(np.abs(np.sum(m_[eng] * avg, 1)), 0, 1)))))
        h0 = lead(Qm[both, 0] / cm_[both, 0][:, None, None]); h1 = lead(Qm[both, 1] / cm_[both, 1][:, None, None])
        r["half_disagreement_median"] = float(np.median(np.degrees(np.arccos(np.clip(np.abs(np.sum(h0 * h1, 1)), 0, 1)))))
        r["own_actins_to_director_median"] = float(np.median(np.degrees(np.arccos(np.clip(np.abs(avg @ nc), 0, 1)))))
        lat = np.linalg.norm(np.stack([f["actin/tb_aa_lateral_%s" % c][t, :, 0] for c in "xyz"], 1), axis=1) / KT
        r["lateral_mag_mean"] = float(lat.mean())
        if t > 0:
            ds = step[t] - step[t - 1]
            r["torque_cap_rate"] = float((f["actin/cumulative_torque_cap_events"][t, :, 0] - f["actin/cumulative_torque_cap_events"][t - 1, :, 0]).sum() / ds / u.shape[0])
            r["rotation_cap_rate"] = float((f["actin/cumulative_rotation_cap_events"][t, :, 0] - f["actin/cumulative_rotation_cap_events"][t - 1, :, 0]).sum() / ds / u.shape[0])
        rows.append(r)
    out[arm] = dict(label=label, source=src, hph=int(a["am_heads_per_half"]), frames=n, rows=rows)
    L = rows[-1]; B = [r for r in rows if r["t_ms"] <= 2000.5][-1]
    print("%-10s %-24s t %.0f ms | pair angle %.1f -> %.1f | actin-director %.1f -> %.1f | S_a %.3f -> %.3f S_m %.3f -> %.3f | links %.2f -> %.2f | |lat| %.0f -> %.0f kT | cap %.3f -> %.3f"
          % (arm, label, L["t_ms"], B["pair_angle_median"], L["pair_angle_median"], B["actin_to_director_median"], L["actin_to_director_median"],
             B["S_actin"], L["S_actin"], B["S_myosin"], L["S_myosin"], B["links_per_actin"], L["links_per_actin"], B["lateral_mag_mean"], L["lateral_mag_mean"],
             B.get("torque_cap_rate", float("nan")), L.get("torque_cap_rate", float("nan"))), flush=True)
json.dump(out, open(sys.argv[1] if len(sys.argv) > 1 else "resab_an.json", "w"))
