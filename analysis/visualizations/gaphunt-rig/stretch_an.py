"""What stretches an ACTN2 link past the 120 nm rupture length?  Production control r27ctrl,
overstretch (geometric_contour) detachments 1.0 -> 2.4 s.
For each one, the two material sites (s_i, s_j recorded on the detachment event) are placed on
the two actins in the saved frame before (t0) and after (t1) the death step:
  site = centre + s * L * u           (validated against the recorded connector below)
and the change of the connector, site_j - site_i, is split exactly into
  slide    : each actin's centre moving ALONG its own axis
  lateral  : each actin's centre moving PERPENDICULAR to its axis
  rotation : each actin turning (s * L * du, lever about the centre)
each projected on the connector direction (contribution to the change in length).
Control: every link alive in both of a pair of consecutive frames (same split).
Also per actin at t0: motor force along the axis (am_active_force), the saved total force,
and t0->t1 increments of the Z-disc tip projection / force / torque / rotation cap counters."""
import h5py, numpy as np, os, json
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
P = "/scratch/users/xsli/sarcomere/gaphunt_20261023/data/r27ctrl_seed20260805.h5"
f = h5py.File(P, "r", locking=False); a = dict(f["meta"].attrs); dt = float(a["dt"]); L = float(a["actin_length"])
box = np.array([a["box_x"], a["box_y"], a["box_z"]], float)
mi = lambda v: v - box * np.round(v / box)
U = lambda d: d / np.linalg.norm(d, axis=-1, keepdims=True)
EC = {c: i for i, c in enumerate(a["aa_local_link_events_v3_columns"].decode().split(","))}
MC = {c: i for i, c in enumerate(a["aa_local_link_mechanics_columns"].decode().split(","))}
nF = int(f["meta/committed_frames"][0]); fsteps = f["state/current_step"][:nF, 0]
S0 = 100000; F0 = int(np.searchsorted(fsteps, S0)) - 1
C = f["actin/center"][F0:nF]; D = U(f["actin/direction"][F0:nF])
fr = lambda t: t - F0
off = f["catch_bond/local_link_frame_offsets"]; ident = f["catch_bond/local_link_identity"]
coord = f["catch_bond/local_link_coordinates"]; mech = f["catch_bond/local_link_mechanics"]
def sites(t, i, j, si, sj):
    return C[fr(t)][i] + (si * L)[:, None] * D[fr(t)][i], C[fr(t)][j] + (sj * L)[:, None] * D[fr(t)][j]
# --- convention check against the recorded connector
t = 400; o, c = int(off[t][1]), int(off[t][2]); idr = ident[o:o + c]; co = coord[o:o + c]; mr = mech[o:o + c]
pi, pj = sites(t, idr[:, 2].astype(int), idr[:, 3].astype(int), co[:, 0], co[:, 1])
rec = mr[:, [MC["connector_x"], MC["connector_y"], MC["connector_z"]]]; mine = mi(pj - pi)
err_pos, err_neg = np.max(np.abs(mine - rec)), np.max(np.abs(-mine - rec))
print("convention check: max |rebuilt - recorded connector| = %.2e um (as j - i), %.2e (as i - j); length err %.2e"
      % (err_pos, err_neg, np.max(np.abs(np.linalg.norm(mine, axis=1) - mr[:, MC["length"]]))))
assert min(err_pos, err_neg) < 1e-6, "site convention does not reproduce the connector"
def split(t0, i, j, si, sj):
    """exact split of the connector change t0 -> t0+1 into slide / lateral / rotation, projected on the mid connector"""
    t1 = t0 + 1
    c0 = mi(np.diff(np.stack(sites(t0, i, j, si, sj)), axis=0)[0]); c1 = mi(np.diff(np.stack(sites(t1, i, j, si, sj)), axis=0)[0])
    n = U(c0 + c1)
    out = {}
    for name, k, s, sg in (("i", i, si, -1.0), ("j", j, sj, 1.0)):
        u0 = D[fr(t0)][k]; dc = mi(C[fr(t1)][k] - C[fr(t0)][k]); du = D[fr(t1)][k] - u0
        ax = np.sum(dc * u0, 1)[:, None] * u0
        out[name] = dict(slide=sg * np.sum(ax * n, 1), lateral=sg * np.sum((dc - ax) * n, 1), rotation=sg * np.sum((s * L)[:, None] * du * n, 1),
                         move=np.linalg.norm(dc + (s * L)[:, None] * du, axis=1))
    tot = {m: out["i"][m] + out["j"][m] for m in ("slide", "lateral", "rotation")}
    return np.linalg.norm(c0, axis=1), np.linalg.norm(c1, axis=1), tot, out
# --- overstretch deaths
ev = f["catch_bond/local_link_events_v3"][:]
d = ev[(ev[:, EC["new_status"]] == 0) & (ev[:, EC["step"]] >= S0) & (ev[:, EC["detachment_reason"]] == 1)]
t0 = np.searchsorted(fsteps, d[:, EC["step"]], side="right") - 1; ok = t0 + 1 < nF; d, t0 = d[ok], t0[ok]
i, j = d[:, EC["actin_i"]].astype(int), d[:, EC["actin_j"]].astype(int); si, sj = d[:, EC["s_i"]], d[:, EC["s_j"]]
st = d[:, EC["old_status"]].astype(int)
R = dict(len0=np.zeros(len(d)), len1=np.zeros(len(d)), slide=np.zeros(len(d)), lateral=np.zeros(len(d)), rotation=np.zeros(len(d)),
         mi_=np.zeros(len(d)), mj_=np.zeros(len(d)))
for t in np.unique(t0):
    q = t0 == t
    l0, l1, tot, per = split(t, i[q], j[q], si[q], sj[q])
    R["len0"][q], R["len1"][q] = l0, l1
    for m in ("slide", "lateral", "rotation"): R[m][q] = tot[m]
    R["mi_"][q], R["mj_"][q] = per["i"]["move"], per["j"]["move"]
# --- control: all links alive in consecutive frame pairs (every 4th frame pair to keep it light)
K = dict(len0=[], len1=[], slide=[], lateral=[], rotation=[])
for t in range(F0 + 1, nF - 1, 4):
    o, c = int(off[t][1]), int(off[t][2]); idr = ident[o:o + c]; co = coord[o:o + c]
    o1, c1 = int(off[t + 1][1]), int(off[t + 1][2]); alive = np.isin(idr[:, 1], ident[o1:o1 + c1][:, 1])
    l0, l1, tot, _ = split(t, idr[alive, 2].astype(int), idr[alive, 3].astype(int), co[alive, 0], co[alive, 1])
    K["len0"].append(l0); K["len1"].append(l1)
    for m in ("slide", "lateral", "rotation"): K[m].append(tot[m])
K = {k: np.concatenate(v) for k, v in K.items()}
def summ(X, sel, lab):
    dl = X["len1"][sel] - X["len0"][sel]
    r = dict(n=int(sel.sum()), len0_med_nm=float(1000 * np.median(X["len0"][sel])), len1_med_nm=float(1000 * np.median(X["len1"][sel])),
             dlen_med_nm=float(1000 * np.median(dl)))
    for m in ("slide", "lateral", "rotation"):
        r[m + "_med_nm"] = float(1000 * np.median(X[m][sel])); r[m + "_mean_nm"] = float(1000 * np.mean(X[m][sel]))
        r[m + "_dominant"] = float(np.mean(np.argmax(np.stack([X[k][sel] for k in ("slide", "lateral", "rotation")]), 0) == ("slide", "lateral", "rotation").index(m)))
    print("%-34s n %6d | length %5.1f -> %5.1f nm (median change %+5.1f) | contribution to the stretch, mean (median): slide %+5.1f (%+5.1f), lateral %+5.1f (%+5.1f), rotation %+5.1f (%+5.1f) nm"
          " | largest term: slide %3.0f%% lateral %3.0f%% rotation %3.0f%%"
          % (lab, r["n"], r["len0_med_nm"], r["len1_med_nm"], r["dlen_med_nm"], r["slide_mean_nm"], r["slide_med_nm"], r["lateral_mean_nm"], r["lateral_med_nm"],
             r["rotation_mean_nm"], r["rotation_med_nm"], 100 * r["slide_dominant"], 100 * r["lateral_dominant"], 100 * r["rotation_dominant"]))
    return r
out = dict(source=P, n_overstretch=int(len(d)))
print()
out["over_s1"] = summ(R, st == 1, "overstretch, status 1")
out["over_s2"] = summ(R, st == 2, "overstretch, status 2")
out["control"] = summ(K, np.ones(len(K["len0"]), bool), "control: surviving links")
dk = K["len1"] - K["len0"]; top = dk >= np.percentile(dk, 99)
out["control_top1"] = summ(K, top, "control: top 1% stretch, survivors")
# who moves: the larger-moving actin's share
big = np.maximum(R["mi_"], R["mj_"]); small = np.minimum(R["mi_"], R["mj_"])
print("\noverstretch: site displacement over the frame, larger-moving actin %.1f nm vs the other %.1f nm (medians)" % (1000 * np.median(big), 1000 * np.median(small)))
out["move_big_small_med_nm"] = [float(1000 * np.median(big)), float(1000 * np.median(small))]
# --- forces on the two actins at t0 vs control actins in the same frames
def actin_feats(t, k):
    u = D[fr(t)][k]
    mf = np.stack([f["actin/am_active_force_" + x][t][k, 0] for x in "xyz"], 1)
    tf = f["actin/force"][t][k]
    out = dict(motor_axial=np.abs(np.sum(mf * u, 1)), total_axial=np.abs(np.sum(tf * u, 1)),
               total_perp=np.linalg.norm(tf - np.sum(tf * u, 1)[:, None] * u, axis=1),
               endstop=f["actin/aa_barbed_tip_endstop_force"][t][k, 0], zdisc=f["actin/aa_zdisc_tip_force"][t][k, 0],
               strong=f["actin/am_strong_count"][t][k, 0])
    for cnt in ("cumulative_aa_zdisc_tip_projection_events", "cumulative_force_cap_events", "cumulative_torque_cap_events", "cumulative_rotation_cap_events"):
        out[cnt.replace("cumulative_", "d_")] = f["actin/" + cnt][t + 1][k, 0] - f["actin/" + cnt][t][k, 0]
    return out
rng = np.random.default_rng(0); feats_o, feats_c = [], []
for t in np.unique(t0)[::3]:
    q = t0 == t; ks = np.unique(np.concatenate([i[q], j[q]]))
    feats_o.append(actin_feats(t, ks)); feats_c.append(actin_feats(t, rng.choice(1600, 200, replace=False)))
cat = lambda L_, k: np.concatenate([x[k] for x in L_])
print("\nactins at t0: overstretch-link actins vs random actins (same frames)")
out["forces"] = {}
for k in feats_o[0]:
    ov, cv = cat(feats_o, k), cat(feats_c, k)
    out["forces"][k] = dict(over_mean=float(np.mean(ov)), ctrl_mean=float(np.mean(cv)), over_frac_pos=float(np.mean(ov > 0)), ctrl_frac_pos=float(np.mean(cv > 0)))
    print("   %-40s overstretch mean %9.3f (nonzero %4.1f%%) | random %9.3f (nonzero %4.1f%%)" % (k, np.mean(ov), 100 * np.mean(ov > 0), np.mean(cv), 100 * np.mean(cv > 0)))
json.dump(out, open("/scratch/users/xsli/sarcomere/gaphunt_20261023/torque_out/stretch_an.json", "w"))
