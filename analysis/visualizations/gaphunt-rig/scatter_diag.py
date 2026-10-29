"""What sets the residual scatter?  For a frame:
- common director n0 (actins + myosins, Q-tensor leading eigenvector)
- per engaged myosin: head-weighted Q-average of its bound actins' axes, per half and
  both halves; angle myosin <-> that average; angle between the two halves' averages
- per bound actin: angle to the head-weighted average of its myosins, and to its
  ACTN2 partners' average
- orientation correlation <P2(cos)> vs centre distance, myosin-myosin and actin-actin
- S of actins / myosins over the whole trajectory (source run, every 20 frames)"""
import h5py, numpy as np, json, os, sys, collections
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
sys.path.insert(0, "/scratch/users/xsli/sarcomere/actn2_rest74_lateral_20260831/source_snapshot")
CAMP = "/scratch/users/xsli/sarcomere/gaphunt_20261023/data"
def U(d): return d / np.linalg.norm(d, axis=-1, keepdims=True)
def lead(Q): w, v = np.linalg.eigh(Q); return v[..., -1], w[..., -1]
def nem(d): d = U(d); Q = 1.5 * (d[:, :, None] * d[:, None, :]).mean(0) - 0.5 * np.eye(3); return float(np.linalg.eigvalsh(Q)[-1])
ang = lambda a, b: np.degrees(np.arccos(np.clip(np.abs(np.sum(a * b, -1)), 0, 1)))
def g2(c, d, box, rmax=2.0, nb=10, sample=None):
    n = len(c); idx = np.arange(n) if sample is None else sample
    dc = c[idx][:, None, :] - c[None, :, :]; dc -= box * np.round(dc / box); r = np.linalg.norm(dc, axis=-1)
    p2 = 1.5 * np.sum(d[idx][:, None, :] * d[None, :, :], -1) ** 2 - 0.5
    edges = np.linspace(0, rmax, nb + 1); out = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (r > lo) & (r <= hi) & (r > 1e-9); out.append((round(float(0.5 * (lo + hi)), 2), float(p2[m].mean()) if m.any() else None, int(m.sum())))
    return out
def frame(f, t):
    a = dict(f["meta"].attrs); box = np.array([a["box_x"], a["box_y"], a["box_z"]]); hph = int(a["am_heads_per_half"])
    u = U(f["actin/direction"][t]); m = U(f["myosin/direction"][t]); ca = f["actin/center"][t]; cm = f["myosin/center"][t]
    n0, _ = lead((np.vstack([u, m])[:, :, None] * np.vstack([u, m])[:, None, :]).mean(0))
    hs = f["state/am_head_state_current"][t]; hb = f["state/am_head_bound_actin_current"][t].astype(int)
    mi, hi = np.nonzero((hs == 1) | (hs == 2)); ai = hb[mi, hi]; half = (hi >= hph).astype(int)
    Qm = np.zeros((len(m), 2, 3, 3)); cnt = np.zeros((len(m), 2))
    np.add.at(Qm, (mi, half), u[ai][:, :, None] * u[ai][:, None, :]); np.add.at(cnt, (mi, half), 1)
    Qa = np.zeros((len(u), 3, 3)); ca_ = np.zeros(len(u)); np.add.at(Qa, ai, m[mi][:, :, None] * m[mi][:, None, :]); np.add.at(ca_, ai, 1)
    eng = cnt.sum(1) > 0; both = (cnt[:, 0] > 0) & (cnt[:, 1] > 0)
    avg_all, _ = lead(Qm.sum(1)[eng] / cnt.sum(1)[eng][:, None, None])
    avg_h0, _ = lead(Qm[both, 0] / cnt[both, 0][:, None, None]); avg_h1, _ = lead(Qm[both, 1] / cnt[both, 1][:, None, None])
    ba = ca_ > 0; avg_m, _ = lead(Qa[ba] / ca_[ba][:, None, None])
    # ACTN2 partners
    off = f["catch_bond/local_link_frame_offsets"][t]; ident = f["catch_bond/local_link_identity"][int(off[1]):int(off[1]) + int(off[2])]
    Qp = np.zeros((len(u), 3, 3)); cp = np.zeros(len(u))
    if len(ident):
        i_, j_ = ident[:, 2].astype(int), ident[:, 3].astype(int)
        np.add.at(Qp, i_, u[j_][:, :, None] * u[j_][:, None, :]); np.add.at(Qp, j_, u[i_][:, :, None] * u[i_][:, None, :]); np.add.at(cp, i_, 1); np.add.at(cp, j_, 1)
    la = cp > 0; avg_p, _ = lead(Qp[la] / cp[la][:, None, None])
    med = lambda x: float(np.median(x))
    rs = np.random.default_rng(0).choice(len(u), 400, replace=False)
    return dict(
        S_actin=nem(u), S_myosin=nem(m),
        actin_to_n0=med(ang(u, n0[None])), myosin_to_n0=med(ang(m, n0[None])),
        engaged_myosins=int(eng.sum()), actins_per_engaged_myosin=float(np.mean([len(set(ai[mi == j])) for j in np.flatnonzero(eng)])),
        myosin_to_own_actin_avg=med(ang(m[eng], avg_all)), own_actin_avg_to_n0=med(ang(avg_all, n0[None])),
        myosins_with_both_halves=int(both.sum()), half_avg_disagreement=med(ang(avg_h0, avg_h1)),
        myosin_to_half0_avg=med(ang(m[both], avg_h0)), myosin_to_half1_avg=med(ang(m[both], avg_h1)),
        actin_to_own_myosin_avg=med(ang(u[ba], avg_m)), own_myosin_avg_to_n0=med(ang(avg_m, n0[None])),
        actin_to_actn2_partner_avg=med(ang(u[la], avg_p)), linked_actins=int(la.sum()),
        g2_myosin=g2(cm, m, box), g2_actin=g2(ca, u, box, sample=rs))
out = {}
for tag, fn, ts in (("hph27 control 2.0s", "hph27m24", [400]), ("hph27 control 2.4s", "r27ctrl", [480]), ("hph27 lateral off 2.4s", "r27nolat", [480]),
                    ("hph81 control 2.4s", "r81ctrl", [480]), ("hph81 lateral off 2.4s", "r81nolat", [480])):
    f = h5py.File("%s/%s_seed20260805.h5" % (CAMP, fn), "r", locking=False)
    r = frame(f, ts[0]); out[tag] = r
    print("== %s" % tag); [print("   %-28s %s" % (k, ("%.3f" % v) if isinstance(v, float) else v)) for k, v in r.items() if not k.startswith("g2")]
    print("   g2 myosin (r um, <P2>, n):", [(x, round(y, 3) if y is not None else None) for x, y, _ in r["g2_myosin"]])
    print("   g2 actin  (r um, <P2>, n):", [(x, round(y, 3) if y is not None else None) for x, y, _ in r["g2_actin"]], flush=True)
# order over the whole trajectory
for fn in ("hph27m24", "hph81m24"):
    f = h5py.File("%s/%s_seed20260805.h5" % (CAMP, fn), "r", locking=False); n = int(f["meta/committed_frames"][0])
    ts = list(range(0, n, 20)); rows = [(t * 5, nem(f["actin/direction"][t]), nem(f["myosin/direction"][t])) for t in ts]
    out[fn + "_S_time"] = rows
    print(fn, "S over time (ms, actin, myosin):", [(t, round(a_, 3), round(m_, 3)) for t, a_, m_ in rows])
json.dump(out, open("/scratch/users/xsli/sarcomere/gaphunt_20261023/torque_out/scatter_diag.json", "w"))
