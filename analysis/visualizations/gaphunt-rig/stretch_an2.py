"""Follow-up to stretch_an.py (same data, same site convention).
1. The actin that moves more over the frame containing the overstretch ("mover"): its mode
   (slide along own axis / lateral / rotation), whether it had Z-disc tip projection events in
   that frame, whether its barbed tip was in Z-disc contact, and its displacement with and without
   projection. Control: the larger-moving actin of the top-1%-stretch surviving links, and random actins.
2. Connector length at formation for links that later overstretch vs links that die otherwise."""
import h5py, numpy as np, os, json
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
P = "/scratch/users/xsli/sarcomere/gaphunt_20261023/data/r27ctrl_seed20260805.h5"
f = h5py.File(P, "r", locking=False); a = dict(f["meta"].attrs); L = float(a["actin_length"])
box = np.array([a["box_x"], a["box_y"], a["box_z"]], float); mi = lambda v: v - box * np.round(v / box)
U = lambda d: d / np.linalg.norm(d, axis=-1, keepdims=True)
EC = {c: i for i, c in enumerate(a["aa_local_link_events_v3_columns"].decode().split(","))}
nF = int(f["meta/committed_frames"][0]); fsteps = f["state/current_step"][:nF, 0]
S0 = 100000; F0 = int(np.searchsorted(fsteps, S0)) - 1
C = f["actin/center"][F0:nF]; D = U(f["actin/direction"][F0:nF])
PJ = f["actin/cumulative_aa_zdisc_tip_projection_events"][F0:nF, :, 0]; PC = f["actin/cumulative_aa_zdisc_tip_projection_correction"][F0:nF, :, 0]
ZC = f["actin/aa_zdisc_tip_contact_count"][F0:nF, :, 0]
fr = lambda t: t - F0
def motion(t, k, s):
    u0 = D[fr(t)][k]; dc = mi(C[fr(t + 1)][k] - C[fr(t)][k]); du = D[fr(t + 1)][k] - u0
    sl = np.abs(np.sum(dc * u0, 1)); lat = np.linalg.norm(dc - np.sum(dc * u0, 1)[:, None] * u0, axis=1); rot = np.linalg.norm((s * L)[:, None] * du, axis=1)
    return np.linalg.norm(dc + (s * L)[:, None] * du, axis=1), sl, lat, rot
ev = f["catch_bond/local_link_events_v3"][:]
d = ev[(ev[:, EC["new_status"]] == 0) & (ev[:, EC["step"]] >= S0)]
rsn = d[:, EC["detachment_reason"]].astype(int)
t0 = np.searchsorted(fsteps, d[:, EC["step"]], side="right") - 1; ok = (t0 + 1 < nF) & (t0 >= F0); d, t0, rsn = d[ok], t0[ok], rsn[ok]
i, j = d[:, EC["actin_i"]].astype(int), d[:, EC["actin_j"]].astype(int); si, sj = d[:, EC["s_i"]], d[:, EC["s_j"]]
mv_i, mv_j = np.zeros(len(d)), np.zeros(len(d)); feats = np.zeros((len(d), 8)); OTH = np.zeros(len(d))
for t in np.unique(t0):
    q = np.flatnonzero(t0 == t)
    Mi = motion(t, i[q], si[q]); Mj = motion(t, j[q], sj[q])
    pick_i = Mi[0] >= Mj[0]; k = np.where(pick_i, i[q], j[q]); other = np.where(pick_i, j[q], i[q])
    M = [np.where(pick_i, x, y) for x, y in zip(Mi, Mj)]; Mo = np.where(pick_i, Mj[0], Mi[0])
    dpj = PJ[fr(t + 1)][k] - PJ[fr(t)][k]; dpc = PC[fr(t + 1)][k] - PC[fr(t)][k]; dpo = PJ[fr(t + 1)][other] - PJ[fr(t)][other]
    feats[q] = np.stack([M[0], M[1], M[2], M[3], Mo, dpj, dpc, ZC[fr(t)][k] > 0], 1); OTH[q] = dpo
names = ["mover_disp", "mover_slide", "mover_lateral", "mover_rotation", "other_disp", "mover_proj_events", "mover_proj_correction", "mover_tip_in_contact"]
# random-actin baseline for the same quantities
rng = np.random.default_rng(1); base = []
for t in range(F0, nF - 1, 3):
    k = rng.choice(1600, 300, replace=False); s = np.full(300, -0.45)
    m = motion(t, k, s); base.append(np.stack([m[0], m[1], m[2], m[3], PJ[fr(t + 1)][k] - PJ[fr(t)][k], ZC[fr(t)][k] > 0], 1))
base = np.concatenate(base)
out = dict(source=P)
def row(sel, lab):
    X = feats[sel]; pj = X[:, 5] > 0
    r = dict(n=int(sel.sum()), mover_disp_med_nm=float(1000 * np.median(X[:, 0])), other_disp_med_nm=float(1000 * np.median(X[:, 4])),
             slide_med_nm=float(1000 * np.median(X[:, 1])), lateral_med_nm=float(1000 * np.median(X[:, 2])), rotation_med_nm=float(1000 * np.median(X[:, 3])),
             mover_projected=float(np.mean(pj)), mover_tip_contact=float(np.mean(X[:, 7] > 0)),
             disp_if_projected_nm=float(1000 * np.median(X[pj, 0])) if pj.any() else None, disp_if_not_nm=float(1000 * np.median(X[~pj, 0])),
             proj_correction_med_nm=float(1000 * np.median(X[pj, 6])) if pj.any() else None)
    po = X[:, 5] > 0; oo = OTH[sel] > 0
    r["other_projected"] = float(np.mean(oo)); r["mover_only_projected"] = float(np.mean(po & ~oo)); r["both_projected"] = float(np.mean(po & oo))
    print("   %-26s other actin projected %3.0f%% | mover only %3.0f%%, both %3.0f%%" % (lab, 100 * np.mean(oo), 100 * np.mean(po & ~oo), 100 * np.mean(po & oo)))
    out[lab] = r
    print("%-26s n %6d | mover %5.1f nm vs other %5.1f nm | mover slide %5.1f, lateral %5.1f, rotation %5.1f nm | mover projected %3.0f%% (tip in contact %3.0f%%) | mover disp if projected %5.1f vs not %5.1f nm | projection correction %s nm"
          % (lab, r["n"], r["mover_disp_med_nm"], r["other_disp_med_nm"], r["slide_med_nm"], r["lateral_med_nm"], r["rotation_med_nm"], 100 * r["mover_projected"],
             100 * r["mover_tip_contact"], r["disp_if_projected_nm"] or np.nan, r["disp_if_not_nm"], "%.1f" % r["proj_correction_med_nm"] if r["proj_correction_med_nm"] is not None else "-"))
st = d[:, EC["old_status"]].astype(int)
print("frame containing the detachment, the actin that moved more ('mover'), medians:")
row((rsn == 1) & (st == 1), "overstretch, status 1"); row((rsn == 1) & (st == 2), "overstretch, status 2")
row((rsn == 2) & (st == 1), "random draw, status 1"); row((rsn == 2) & (st == 2), "random draw, status 2")
pjb = base[:, 4] > 0
print("%-26s n %6d | site disp %5.1f nm | slide %5.1f, lateral %5.1f, rotation %5.1f nm | projected %3.0f%% (tip in contact %3.0f%%) | disp if projected %5.1f vs not %5.1f nm"
      % ("random actins", len(base), 1000 * np.median(base[:, 0]), 1000 * np.median(base[:, 1]), 1000 * np.median(base[:, 2]), 1000 * np.median(base[:, 3]),
         100 * np.mean(pjb), 100 * np.mean(base[:, 5] > 0), 1000 * np.median(base[pjb, 0]), 1000 * np.median(base[~pjb, 0])))
out["random_actins"] = dict(n=int(len(base)), disp_med_nm=float(1000 * np.median(base[:, 0])), projected=float(np.mean(pjb)),
                            disp_if_projected_nm=float(1000 * np.median(base[pjb, 0])), disp_if_not_nm=float(1000 * np.median(base[~pjb, 0])))
# formation length
fm = ev[(ev[:, EC["old_status"]] == 0) & (ev[:, EC["new_status"]] == 1)]
flen = {int(r[EC["link_id"]]): r[EC["connector_length"]] for r in fm}
print("\nconnector length at formation (median, p90):")
for code, lab in ((1, "later overstretch"), (2, "later random draw"), (3, "later burial")):
    v = np.array([flen.get(int(l), np.nan) for l in d[rsn == code][:, EC["link_id"]]]); v = v[~np.isnan(v)]
    print("   %-18s %5.1f nm, p90 %5.1f nm (n %d)" % (lab, 1000 * np.median(v), 1000 * np.percentile(v, 90), len(v)))
    out["formation_len_" + lab.replace(" ", "_")] = [float(1000 * np.median(v)), float(1000 * np.percentile(v, 90))]
json.dump(out, open("/scratch/users/xsli/sarcomere/gaphunt_20261023/torque_out/stretch_an2.json", "w"))
