"""Case study: Cluster E (2 myosins, 8 actins), production control r27ctrl, 2.0 -> 2.4 s, every 5 ms.

Per member actin:
  * angle to its primary myosin (cluster myosin with most bound heads over the window) and to the
    actin-only local director (mean of the cluster's 8 actins), heads bound per frame
  * recorded per-source torque (tb_*, kT) projected on the tangent that turns the actin toward
    its primary myosin; mean +- sem over frames
  * ACTN2 links: count, status, partner in/out of cluster, partner angle
  * orientational MSD vs lag against free tumbling (<dtheta^2> = 4 D t)
  * geometric lock test at 2.0 / 2.2 / 2.4 s: rotate the actin toward its primary myosin by
    2/5/10/20 deg about its centre and about its barbed end, and count what that creates:
    actin-actin capsule overlaps (< 30 nm), actin-myosin overlaps (< 65 nm), ACTN2 links pulled
    past rupture (0.12 um), barbed-end neighbours inside the 80 nm lateral range (within +-250 nm
    axially).  Also counts at the unrotated pose for reference."""
import os, sys, json, collections, numpy as np, h5py
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import build_cluster_traj as B
KT = 1.380649e-23 * 310 / 1e-18
SRC = ["am_head", "actn2_align", "actn2_lever", "actn2_perp", "actn2_axial", "actn2_contact", "aa_lateral", "aa_steric", "am_steric"]
U = lambda d: d / np.linalg.norm(d, axis=-1, keepdims=True)
def mimg(d, box): return d - box * np.round(d / box)
def ang(a, b): return float(np.degrees(np.arccos(min(1.0, abs(float(np.dot(a, b)))))))

def seg_dist(p0, p1, Q0, Q1):
    """Min distance between segment p0-p1 and many segments Q0-Q1 (already min-imaged near p)."""
    d1 = p1 - p0; d2 = Q1 - Q0; r = p0 - Q0
    a = d1 @ d1; e = np.sum(d2 * d2, 1); f = np.sum(d2 * r, 1); c = r @ d1; b = d2 @ d1
    den = a * e - b * b
    s = np.where(den > 1e-12, np.clip((b * f - c * e) / np.where(den > 1e-12, den, 1), 0, 1), 0.0)
    t = (b * s + f) / np.where(e > 1e-12, e, 1)
    s = np.where(t < 0, np.clip(-c / a, 0, 1), np.where(t > 1, np.clip((b - c) / a, 0, 1), s)); t = np.clip(t, 0, 1)
    cp = p0 + s[:, None] * d1; cq = Q0 + t[:, None] * d2
    return np.linalg.norm(cp - cq, axis=1)

def rotate_toward(u, target, deg):
    """Rotate unit u toward +-target by deg (in the plane of u and target)."""
    tg = np.sign(u @ target) * target; t = tg - u * (u @ tg); n = np.linalg.norm(t)
    if n < 1e-9: return u
    t /= n; a = np.radians(deg); return np.cos(a) * u + np.sin(a) * t

def lock_test(f, t, i, mprim, box, La, Lm, Ra, Rm, rupt, lat_range, lat_win):
    ca = f["actin/center"][t]; ua = U(f["actin/direction"][t]); cm = f["myosin/center"][t]; um = U(f["myosin/direction"][t])
    off = f["catch_bond/local_link_frame_offsets"][t]; o, cnt = int(off[1]), int(off[2])
    ids = f["catch_bond/local_link_identity"][o:o + cnt]; cs = f["catch_bond/local_link_coordinates"][o:o + cnt]
    mine = [(int(r[3]) if int(r[2]) == i else int(r[2]), float(c[1]) if int(r[2]) == i else float(c[0]), float(c[0]) if int(r[2]) == i else float(c[1]))
            for r, c in zip(ids, cs) if i in (int(r[2]), int(r[3]))]          # (partner, partner coord, my coord)
    others = np.array([k for k in range(len(ca)) if k != i])
    def count(ci, ui):
        p0, p1 = ci - 0.5 * La * ui, ci + 0.5 * La * ui
        dA = mimg(ca[others] - ci, box); qa = ci + dA
        nAA = int((seg_dist(p0, p1, qa - 0.5 * La * ua[others], qa + 0.5 * La * ua[others]) < 2 * Ra).sum())
        dM = mimg(cm - ci, box); qm = ci + dM
        nAM = int((seg_dist(p0, p1, qm - 0.5 * Lm * um, qm + 0.5 * Lm * um) < Rm + Ra).sum())
        nL = 0
        for j, sj, si in mine:
            pj = ci + mimg(ca[j] - ci, box) + sj * ua[j]; pi_ = ci + si * ui
            if np.linalg.norm(pj - pi_) > rupt: nL += 1
        bi = ci - 0.5 * La * ui; bj = ci + mimg(ca[others] - ci, box) - 0.5 * La * ua[others]
        sep = bj - bi; ax = U(ui + np.where((ua[others] @ ui)[:, None] >= 0, ua[others], -ua[others]))
        axo = np.sum(sep * ax, 1); tr = np.linalg.norm(sep - axo[:, None] * ax, axis=1)
        nB = int(((np.abs(axo) <= lat_win) & (tr < lat_range)).sum())
        return nAA, nAM, nL, nB
    base = count(ca[i], ua[i]); res = {"base": base}
    for pivot in ("centre", "barbed"):
        for deg in (2, 5, 10, 20):
            un = U(rotate_toward(ua[i], um[mprim], deg))
            ci = ca[i] if pivot == "centre" else (ca[i] - 0.5 * La * ua[i]) + 0.5 * La * un
            res["%s_%d" % (pivot, deg)] = count(ci, un)
    res["links"] = len(mine)
    return res

f = h5py.File("%s/data/r27ctrl_seed20260805.h5" % B.CAMP, "r", locking=False); a = dict(f["meta"].attrs)
box = np.array([a["box_x"], a["box_y"], a["box_z"]]); La = float(a["actin_length"]); Lm = float(a["myosin_length"])
Ra = float(a.get("actin_radius", 0.015)); Rm = float(a.get("myosin_radius", 0.05)); rupt = float(a.get("aa_link_rupture_length", 0.12))
lat_range = float(a["aa_barbed_lateral_range"]); lat_win = float(a["aa_barbed_lateral_axial_window"]); Dt = float(a["actin_D_tumble"])
clusters, adj = B.find_clusters(f, f["actin/center"].shape[1], f["myosin/center"].shape[1])
c = B.pick(clusters)[4]; A = c["actins"]; M = c["myosins"]
assert (len(M), len(A)) == (2, 8)
frames = list(range(400, 481))
hs_all = f["state/am_head_state_current"]; hb_all = f["state/am_head_bound_actin_current"]
# primary myosin per actin = cluster myosin with most bound-head-frames
cnt = collections.Counter()
for t in frames:
    hs = hs_all[t]; hb = hb_all[t]
    for m in M:
        for h in np.flatnonzero((hs[m] == 1) | (hs[m] == 2)):
            if int(hb[m, h]) in A: cnt[(int(hb[m, h]), m)] += 1
prim = {i: max(M, key=lambda m: cnt[(i, m)]) for i in A}
out = dict(actins=A, myosins=M, primary=prim, per_actin={}, frames_ms=[t * 5.0 for t in frames])
series = {i: collections.defaultdict(list) for i in A}
for t in frames:
    ua = U(f["actin/direction"][t]); um = U(f["myosin/direction"][t])
    hs = hs_all[t]; hb = hb_all[t]
    nloc = np.linalg.eigh((ua[A][:, :, None] * ua[A][:, None, :]).mean(0))[1][:, -1]
    tb = {s: np.stack([f["actin/tb_%s_%s" % (s, k)][t, :, 0] for k in "xyz"], 1) for s in SRC}
    tot = f["actin/torque"][t]
    off = f["catch_bond/local_link_frame_offsets"][t]; o, n_ = int(off[1]), int(off[2]); ids = f["catch_bond/local_link_identity"][o:o + n_]
    for i in A:
        m = prim[i]; tgt = np.sign(ua[i] @ um[m]) * um[m]; tg = tgt - ua[i] * (ua[i] @ tgt); tg /= max(np.linalg.norm(tg), 1e-12)
        S = series[i]
        S["theta"].append(ang(ua[i], um[m])); S["theta_other"].append(ang(ua[i], um[[x for x in M if x != m][0]]))
        S["to_local"].append(ang(ua[i], nloc))
        S["heads"].append(int(sum(((hs[mm] == 1) | (hs[mm] == 2))[hb[mm] == i].sum() for mm in M)))
        S["strong"].append(int(sum((hs[mm] == 2)[hb[mm] == i].sum() for mm in M)))
        for s in SRC: S[s].append(float(tb[s][i] @ tg) / KT)
        S["TOT"].append(float(tot[i] @ tg) / KT)
        lk = [(int(r[3]) if int(r[2]) == i else int(r[2]), int(r[4])) for r in ids if i in (int(r[2]), int(r[3]))]
        S["links"].append(len(lk)); S["links2"].append(sum(1 for _, st in lk if st == 2))
        S["links_out"].append(sum(1 for j, _ in lk if j not in A))
        S["partner_angle"].append(float(np.median([ang(ua[i], ua[j]) for j, _ in lk])) if lk else float("nan"))
cap = f["actin/cumulative_torque_cap_events"]; rcap = f["actin/cumulative_rotation_cap_events"]; st = f["state/current_step"]
steps = float(st[480, 0] - st[400, 0])
for i in A:
    S = series[i]; th = np.array(S["theta"])
    ua_t = U(f["actin/direction"][400:481, i])
    msd = {}
    for lag in (1, 2, 4, 8, 16, 32):
        d = np.arccos(np.clip(np.abs(np.sum(ua_t[lag:] * ua_t[:-lag], 1)), 0, 1)); msd[lag] = float(np.mean(d ** 2))
    locks = {int(t): lock_test(f, t, i, prim[i], box, La, Lm, Ra, Rm, rupt, lat_range, lat_win) for t in (400, 440, 480)}
    pa = dict(
        primary_myosin=int(prim[i]), theta_first=float(th[0]), theta_last=float(th[-1]), theta_mean=float(th.mean()),
        theta_other_mean=float(np.mean(S["theta_other"])), to_local_mean=float(np.mean(S["to_local"])),
        heads_mean=float(np.mean(S["heads"])), strong_mean=float(np.mean(S["strong"])), frac_bound=float(np.mean(np.array(S["heads"]) > 0)),
        links_mean=float(np.mean(S["links"])), links2_mean=float(np.mean(S["links2"])), links_out_mean=float(np.mean(S["links_out"])),
        partner_angle_median=float(np.nanmedian(S["partner_angle"])),
        drive={s: (float(np.mean(S[s])), float(np.std(S[s]) / np.sqrt(len(S[s])))) for s in SRC + ["TOT"]},
        drive_sd={s: float(np.std(S[s])) for s in SRC + ["TOT"]},
        msd=msd, msd_free={lag: 4 * Dt * lag * 0.005 for lag in msd},
        torque_cap_rate=float((cap[480, i, 0] - cap[400, i, 0]) / steps), rotation_cap_rate=float((rcap[480, i, 0] - rcap[400, i, 0]) / steps),
        lock=locks, series={k: [round(x, 3) if isinstance(x, float) else x for x in v] for k, v in S.items()})
    out["per_actin"][int(i)] = pa
    L0 = locks[400]
    print("actin %4d  myo %3d  theta %5.1f->%5.1f (mean %5.1f; other myo %5.1f; to local actin dir %5.1f) heads %.1f (bound %.0f%%) | ACTN2 links %.1f (status-2 %.1f, to outside %.1f), partner angle %.1f | cap %.2f rot-cap %.3f"
          % (i, prim[i], th[0], th[-1], th.mean(), pa["theta_other_mean"], pa["to_local_mean"], pa["heads_mean"], 100 * pa["frac_bound"],
             pa["links_mean"], pa["links2_mean"], pa["links_out_mean"], pa["partner_angle_median"], pa["torque_cap_rate"], pa["rotation_cap_rate"]))
    print("     drive toward myosin (kT, mean): " + ", ".join("%s %+.0f" % (s, pa["drive"][s][0]) for s in SRC + ["TOT"]))
    print("     MSD(lag 1,4,16,32 frames) rad^2: %s vs free %s" % ([round(msd[l], 4) for l in (1, 4, 16, 32)], [round(4 * Dt * l * 0.005, 4) for l in (1, 4, 16, 32)]))
    print("     lock @2.0s (AA overlaps, AM overlaps, links past rupture, barbed nbrs<80nm): base %s | centre 5/10/20: %s %s %s | barbed 5/10/20: %s %s %s"
          % (L0["base"], L0["centre_5"], L0["centre_10"], L0["centre_20"], L0["barbed_5"], L0["barbed_10"], L0["barbed_20"]))
# myosins
for m in M:
    um_t = U(f["myosin/direction"][400:481, m]); print("myosin %d: rotation over window %.1f deg" % (m, ang(um_t[0], um_t[-1])))
json.dump(out, open("/scratch/users/xsli/sarcomere/gaphunt_20261023/torque_out/case_cluster_e.json", "w"))
