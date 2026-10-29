"""Offline impact of the 2026-10-28 connector rules on saved production states.
Replicates, per live ACTN2 link and per myosin (nearest periodic image of the connector midpoint;
the bounding reach <= 0.80 um is below half of every box side, so no other image can qualify):
  DELETE (mature links): some connector point within 0.120 um of the myosin axis AND within the
          axial slab |z| <= 0.75 - 0.010 um, at the same point   (actn2_maturation_geometry.h)
  GATE   (maturation):   same test with 0.130 um and |z| <= 0.75 um
  SOFT new: some connector point on the tip's own half (z >= 0 along the outward axis), z < L/2 + onset
          and rho < 0.130 at the same point                       (new zdisc_myosin_tip_gap)
  SOFT old: the connector MIDPOINT on the tip's own half, gap to the tip < onset, rho < 0.130
Site convention validated earlier: site = centre + s*L*u, connector = shortest image site_j - site_i."""
import h5py, numpy as np, os, sys, json
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
D = "/scratch/users/xsli/sarcomere/gaphunt_20261023/data/"
U = lambda d: d / np.linalg.norm(d, axis=-1, keepdims=True)
def slab_hit(first, conn, a, h, R, strict=False, zlo=None, zhi=None):
    """first, conn: (P,3) relative to the myosin centre; a: (P,3) unit axis.  Axial interval
    [-h, h] (or [zlo, zhi]); returns same-point hit and the min radius on the clipped interval."""
    z0 = np.sum(first * a, 1); dz = np.sum(conn * a, 1)
    lo_z = -h if zlo is None else zlo; hi_z = h if zhi is None else zhi
    with np.errstate(divide="ignore", invalid="ignore"):
        e1 = (lo_z - z0) / dz; e2 = (hi_z - z0) / dz
    entry = np.where(dz != 0, np.minimum(e1, e2), -np.inf); exit_ = np.where(dz != 0, np.maximum(e1, e2), np.inf)
    lo = np.maximum(0.0, entry); hi = np.minimum(1.0, exit_)
    ok = np.where(dz != 0, lo <= hi, (z0 >= lo_z) & (z0 <= hi_z))
    lo = np.where(dz != 0, lo, 0.0); hi = np.where(dz != 0, hi, 1.0)
    rf = first - z0[:, None] * a; rd = conn - dz[:, None] * a; aa = np.sum(rd * rd, 1)
    with np.errstate(divide="ignore", invalid="ignore"):
        tstar = np.where(aa > 0, -np.sum(rf * rd, 1) / aa, lo)
    t = np.clip(tstar, lo, hi)
    rho = np.linalg.norm(rf + t[:, None] * rd, axis=1)
    hit = ok & ((rho < R) if strict else (rho <= R))
    return hit, np.where(ok, rho, np.inf)
def frame_impact(path, t):
    f = h5py.File(path, "r", locking=False); a = dict(f["meta"].attrs)
    box = np.array([a["box_x"], a["box_y"], a["box_z"]], float); L = float(a["actin_length"]); Lm = float(a.get("myosin_length", 1.5))
    mi = lambda v: v - box * np.round(v / box)
    off = f["catch_bond/local_link_frame_offsets"]; o, c = int(off[t][1]), int(off[t][2])
    idr = f["catch_bond/local_link_identity"][o:o + c]; co = f["catch_bond/local_link_coordinates"][o:o + c]
    st = idr[:, 4].astype(int); i, j = idr[:, 2].astype(int), idr[:, 3].astype(int)
    C = f["actin/center"][t]; Ua = U(f["actin/direction"][t]); Cm = f["myosin/center"][t]; Um = U(f["myosin/direction"][t])
    p = C[i] + (co[:, 0] * L)[:, None] * Ua[i]; q = C[j] + (co[:, 1] * L)[:, None] * Ua[j]
    conn = mi(q - p); mid = p + 0.5 * conn
    nL, nM = len(p), len(Cm)
    res = {k: np.zeros(nL, bool) for k in ("del", "gate", "soft_new", "soft_old")}
    minrho = np.full(nL, np.inf)
    for m0 in range(0, nM, 40):                      # chunk myosins
        ms = np.arange(m0, min(nM, m0 + 40)); K = len(ms)
        rel_mid = mi(mid[:, None, :] - Cm[None, ms, :]).reshape(-1, 3)
        cn = np.repeat(conn, K, axis=0) if False else np.broadcast_to(conn[:, None, :], (nL, K, 3)).reshape(-1, 3)
        ax = np.broadcast_to(Um[None, ms, :], (nL, K, 3)).reshape(-1, 3)
        first = rel_mid - 0.5 * cn
        h_del, _ = slab_hit(first, cn, ax, 0.75 * Lm / 1.5 - 0.010, 0.120)
        h_gate, r_gate = slab_hit(first, cn, ax, 0.5 * Lm, 0.130)
        sn = np.zeros(len(first), bool); so = np.zeros(len(first), bool)
        zmid = np.sum(rel_mid * ax, 1); rmid = np.linalg.norm(rel_mid - zmid[:, None] * ax, axis=1)
        for sign in (-1.0, 1.0):
            out = sign * ax
            hs, _ = slab_hit(first, cn, out, None, 0.130, strict=True, zlo=0.0, zhi=0.5 * Lm + 0.05 - 1e-12)
            sn |= hs
            zs = sign * zmid
            so |= (zs >= 0) & (zs - 0.5 * Lm < 0.05) & (rmid < 0.130)
        for k, v in (("del", h_del), ("gate", h_gate), ("soft_new", sn), ("soft_old", so)):
            res[k] |= v.reshape(nL, K).any(1)
        minrho = np.minimum(minrho, r_gate.reshape(nL, K).min(1))
    return dict(n=int(nL), n_s2=int((st == 2).sum()), n_s1=int((st == 1).sum()),
                del_s2=int((res["del"] & (st == 2)).sum()), gate_s1=int((res["gate"] & (st == 1)).sum()),
                gate_s2=int((res["gate"] & (st == 2)).sum()),
                soft_new=int(res["soft_new"].sum()), soft_old=int(res["soft_old"].sum()),
                soft_new_only=int((res["soft_new"] & ~res["soft_old"]).sum()),
                minrho_inslab_med_nm=float(1000 * np.median(minrho[np.isfinite(minrho)])) if np.isfinite(minrho).any() else None)
out = {}
for tag, fn, t in (("production 2.0 s (hph27m24 frame 400, projection on)", "hph27m24_seed20260805.h5", 400),
                   ("r27zt0 2.4 s (projection on)", "r27zt0_seed20260805.h5", 480),
                   ("r27zt0np 2.4 s (projection off)", "r27zt0np_seed20260805.h5", 480)):
    r = frame_impact(D + fn, t); out[tag] = r
    print("%-52s links %d (mature %d, new %d)" % (tag, r["n"], r["n_s2"], r["n_s1"]))
    print("   DELETE now (mature, 120 nm, 10 nm inset): %d = %.1f%% of mature links" % (r["del_s2"], 100 * r["del_s2"] / r["n_s2"]))
    print("   GATE (130 nm): new links that could not mature here %d = %.1f%% of new links; mature links inside it %d (%.1f%%)"
          % (r["gate_s1"], 100 * r["gate_s1"] / r["n_s1"], r["gate_s2"], 100 * r["gate_s2"] / r["n_s2"]))
    print("   SOFT tip force eligible: old midpoint rule %d links, new full-connector rule %d links (%d only under the new rule)"
          % (r["soft_old"], r["soft_new"], r["soft_new_only"]))
json.dump(out, open("/scratch/users/xsli/sarcomere/gaphunt_20261023/torque_out/connector_rule_impact.json", "w"))
