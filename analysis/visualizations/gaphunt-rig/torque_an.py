"""Torque balance on bound A-M pairs: why does the head couple turn the actin in
the rigs but not in production?

Per bound (actin a, myosin m) pair and frame:
  t_am = unit tangent at u_a toward +-u_m  (A = drive . t_am > 0 rotates a TOWARD m)
  head drive rebuilt per head exactly as compute_am_two_spring_head
  (interaction.cpp:572, radial line mode, saturating tanh cap) and split into
    RAD  radial channel (the span couple)        AX  axial channel (motor / weak spring)
    LIN  the radial channel with the cap removed (k=600 linear) -- what saturation costs
  every recorded tb_<source> projected on the same t_am.
Reconstruction is validated per actin against tb_am_head (sum over all partners).
Units: drive in pN um, reported in kT (310 K).  dtheta/dt = D_tumble * A / kT."""
import h5py, numpy as np, json, os, sys, collections
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
CAMP = "/scratch/users/xsli/sarcomere/gaphunt_20261023/data"
KT = 1.380649e-23 * 310 / 1e-18          # pN um
SRC = ["am_head", "am_steric", "am_align", "am_legacy", "cb_bias", "actn2_align",
       "actn2_lever", "actn2_perp", "actn2_axial", "actn2_contact", "aa_lateral", "aa_steric"]

def unit(v):
    n = np.linalg.norm(v, axis=-1, keepdims=True); return v / np.where(n > 1e-12, n, 1.0)
def soft(x, s): return np.where(x > s, x - s, np.where(x < -s, x + s, 0.0))
def mimg(d, box): return d - box * np.round(d / box)
def anchors(hph, Lm, spokes=3):
    cc = max(1, (hph + spokes - 1) // spokes); out = np.empty(2 * hph)
    for half in (0, 1):
        sg = -1.0 if half == 0 else 1.0
        for h in range(hph): out[half * hph + h] = sg * (0.05 * Lm + ((h // spokes + 0.5) / cc) * 0.45 * Lm)
    return out

def heads(f, t, A, anc, box, La):
    hs = f["state/am_head_state_current"][t]; ha = f["state/am_head_bound_actin_current"][t].astype(int)
    tg = f["state/am_head_target_axial_current"][t]; sp = f["state/am_head_stroke_progress_current"][t]
    mi, hi = np.nonzero((hs == 1) | (hs == 2))
    if len(mi) == 0: return None
    ai = ha[mi, hi]; strong = hs[mi, hi] == 2
    ac = f["actin/center"][t]; ad = unit(f["actin/direction"][t]); mc = f["myosin/center"][t]; md = unit(f["myosin/direction"][t])
    ua = ad[ai]; um = md[mi]
    m_att = mc[mi] + anc[hi][:, None] * um
    s = np.clip(tg[mi, hi], -0.5 * La, 0.5 * La)
    conn = mimg(mimg(ac[ai] - m_att, box) + s[:, None] * ua, box)
    z = np.einsum("kj,kj->k", conn, ua); tr = conn - z[:, None] * ua; rho = np.linalg.norm(tr, axis=1)
    loaded = soft(rho - A["rest"], A["slack"])
    sc = A["cap"] / A["krad"]
    dU = A["cap"] * np.tanh(loaded / sc); dUlin = A["krad"] * loaded
    # pair head counts (weak+strong) for the motor scale
    key = ai * 100000 + mi; uk, inv, cnt = np.unique(key, return_inverse=True, return_counts=True)
    nb = cnt[inv]; msc = 1.0 / np.maximum(1, np.minimum(nb, A["msat"]))
    kax = A["motor_k"] * np.where(strong, 1.0, A["weak_ratio"]) * msc
    stall = A["stall"] * msc
    strain = sp[mi, hi] - z
    Fs = np.where(strain > 0, np.minimum(kax * strain, stall), 0.0)
    dUz = np.where(strong, -Fs, kax * z)
    def drive(dUr, dUzz):
        gh = (dUr / np.maximum(rho, 1e-12))[:, None] * tr + dUzz[:, None] * ua
        ex = -(dUr * z / np.maximum(rho, 1e-12))[:, None] * conn + dUzz[:, None] * conn
        g = s[:, None] * gh + ex
        tau = -np.cross(ua, g)
        return np.cross(tau, ua)
    return dict(ai=ai, mi=mi, hi=hi, strong=strong, key=key, loaded=loaded, sc=sc, anc=anc[hi], s=s,
                rad=drive(dU, 0 * dU), ax=drive(0 * dU, dUz), lin=drive(dUlin, 0 * dU))

def run(path, frames, A_over=None):
    f = h5py.File(path, "r", locking=False); a = dict(f["meta"].attrs)
    box = np.array([a["box_x"], a["box_y"], a["box_z"]], float); La = float(a["actin_length"])
    A = dict(rest=float(a["am_preferred_centerline"]), slack=float(a["am_tether_slack"]),
             cap=float(a["am_head_radial_force_max"]),
             krad=float(a["am_head_transverse_k"]) + float(a["am_tether_k"]) / float(a["am_tether_saturation_heads"]),
             msat=int(a["am_motor_saturation_heads"]), motor_k=float(a["am_motor_k"]), weak_ratio=0.1,
             stall=float(a["am_stall_force"]))
    anc = anchors(int(a["am_heads_per_half"]), float(a["myosin_length"]))
    Dt = float(a["actin_D_tumble"]); Dm = float(a["myosin_D_tumble"])
    nfr = int(min(f["meta/committed_frames"][0], f["actin/center"].shape[0]))
    step = f["state/current_step"][:nfr, 0].astype(float); dtf = float(a["dt"]) * float(np.median(np.diff(step)))
    rows = []; val = []
    for t in frames:
        if t + 1 >= nfr: continue
        H = heads(f, t, A, anc, box, La)
        if H is None: continue
        ad = unit(f["actin/direction"][t]); md = unit(f["myosin/direction"][t])
        tb = {s: np.stack([f["actin/tb_%s_%s" % (s, c)][t, :, 0] for c in "xyz"], 1) for s in SRC if "actin/tb_%s_x" % s in f}
        tot = f["actin/torque"][t]; mtq = f["myosin/torque"][t]
        ad1 = unit(f["actin/direction"][t + 1]); md1 = unit(f["myosin/direction"][t + 1])
        mpart = collections.Counter(np.unique(H["key"]) % 100000)
        # validation: per actin, sum of rebuilt (rad+ax) vs tb_am_head
        rec = np.zeros_like(tb["am_head"]); np.add.at(rec, H["ai"], H["rad"] + H["ax"])
        has = np.unique(H["ai"])
        val.append((np.linalg.norm(rec[has] - tb["am_head"][has], axis=1), np.linalg.norm(tb["am_head"][has], axis=1)))
        for k in np.unique(H["key"]):
            sel = H["key"] == k; ai = int(H["ai"][sel][0]); mi = int(H["mi"][sel][0])
            u = ad[ai]; m = md[mi]; c = float(np.dot(u, m)); mm = np.sign(c) * m if c != 0 else m
            tg = mm - u * np.dot(u, mm); sn = np.linalg.norm(tg)
            if sn < 1e-6: continue
            tg = tg / sn; th = float(np.degrees(np.arccos(min(1.0, abs(c)))))
            ld = H["loaded"][sel]
            r = dict(t=int(t), a=ai, m=mi, th=th, n=int(sel.sum()), ns=int(H["strong"][sel].sum()),
                     span=1000 * float(H["anc"][sel].max() - H["anc"][sel].min()),
                     npos=int((ld > 0).sum()), nneg=int((ld < 0).sum()),
                     sat=float(np.mean(np.abs(ld) / H["sc"] > 1.5)), meanld=1000 * float(ld.mean()),
                     smean=float(H["s"][sel].mean()),
                     RAD=float(H["rad"][sel].sum(0) @ tg) / KT, AX=float(H["ax"][sel].sum(0) @ tg) / KT,
                     LIN=float(H["lin"][sel].sum(0) @ tg) / KT,
                     npart=int(len(np.unique(H["mi"][H["ai"] == ai]))),
                     TOT=float(tot[ai] @ tg) / KT)
            # myosin side: tangent at u_m toward +-u_a
            uu = np.sign(c) * u if c != 0 else u; tm = uu - m * np.dot(m, uu); tm = tm / max(np.linalg.norm(tm), 1e-12)
            r["MTOT"] = float(mtq[mi] @ tm) / KT
            r["mpart"] = int(mpart[mi])
            # realized rotation over one frame (rad/s), sign: + = toward the partner
            r["wa"] = float((ad1[ai] - u) @ tg) / dtf; r["wm"] = float((md1[mi] - m) @ tm) / dtf
            c1 = float(np.dot(ad1[ai], md1[mi])); r["dth"] = (float(np.degrees(np.arccos(min(1.0, abs(c1))))) - th) / dtf
            for s_, v in tb.items(): r[s_] = float(v[ai] @ tg) / KT
            r["OTHERPART"] = r["am_head"] - r["RAD"] - r["AX"]
            rows.append(r)
    e = np.concatenate([v[0] for v in val]); n = np.concatenate([v[1] for v in val])
    return rows, dict(dtf=dtf, Dm=Dm, rel_med=float(np.median(e / np.maximum(n, 1e-12))), rel_p95=float(np.percentile(e / np.maximum(n, 1e-12), 95)),
                      abs_max=float(e.max()), scale=float(np.median(n)), Dt=Dt, cap=A["cap"], krad=A["krad"])

if __name__ == "__main__":
    OUTD = sys.argv[1] if len(sys.argv) > 1 else "."
    os.makedirs(OUTD, exist_ok=True)
    out = {}
    jobs = []
    for tag, fn in (("hph27", "hph27m24_seed20260805.h5"), ("hph81", "hph81m24_seed20260805.h5"), ("hph54", "hph54m24_seed20260815.h5")):
        p = CAMP + "/" + fn; nfr = int(h5py.File(p, "r", locking=False)["meta/committed_frames"][0])
        jobs.append((tag, p, range(nfr - 100, nfr)))
    for tag in ("tq30_seed20260805", "tq30_seed20260815", "tq30_seed20260825", "tris45_seed20260805", "tris30_seed20260805", "tri30_seed20260805", "tri45_seed20260805"):
        p = CAMP + "/" + tag + ".h5"
        if os.path.exists(p): jobs.append((tag, p, range(0, 2001, 2)))
    for tag, p, fr in jobs:
        rows, v = run(p, fr)
        cols = {k: np.array([r[k] for r in rows]) for k in rows[0]} if rows else {}
        np.savez_compressed("%s/torque_%s.npz" % (OUTD, tag), **cols)
        out[tag] = v; print(tag, "pairs", len(rows), "meta", v, flush=True)
        json.dump(out, open("%s/torque_meta.json" % OUTD, "w"), indent=1)
