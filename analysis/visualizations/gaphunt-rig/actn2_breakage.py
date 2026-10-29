"""ACTN2 link breakage reasons, production control r27ctrl (hph27, seed 20260805), links that
detach between 1.0 and 2.4 s (source run + its unchanged resume).
Per detachment event (catch_bond/local_link_events_v3, new_status 0):
  reason (code precedence burial > collision > contour > stochastic), status at death,
  age (formation step from the 0->1 event; links formed before the log are dropped from age stats),
  load and assigned lifetime at that step, connector length, A-A angle, registry error,
  and whether the contour and stochastic tests BOTH fired that step.
Geometry of the stretch: the link's last saved frame (<= 5 ms before death) from
catch_bond/local_link_mechanics: connector split into the component along the two actins'
mean axis (axial, sliding) and perpendicular to it (radial, separation)."""
import h5py, numpy as np, os, json, collections
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
P = "/scratch/users/xsli/sarcomere/gaphunt_20261023/data/r27ctrl_seed20260805.h5"
f = h5py.File(P, "r", locking=False); a = dict(f["meta"].attrs); dt = float(a["dt"])
EC = {c: i for i, c in enumerate(a["aa_local_link_events_v3_columns"].decode().split(","))}
MC = {c: i for i, c in enumerate(a["aa_local_link_mechanics_columns"].decode().split(","))}
RN = dict(x.split("=") for x in a["aa_local_link_detachment_reason_codes"].decode().split(","))
ev = f["catch_bond/local_link_events_v3"][:]
S0 = 100000
born = {}
for r in ev[(ev[:, EC["old_status"]] == 0) & (ev[:, EC["new_status"]] == 1)]:
    born.setdefault(int(r[EC["link_id"]]), r[EC["step"]])
d = ev[(ev[:, EC["new_status"]] == 0) & (ev[:, EC["step"]] >= S0)]
span_s = (ev[:, EC["step"]].max() - S0) * dt
col = lambda k: d[:, EC[k]]
rsn = col("detachment_reason").astype(int); st = col("old_status").astype(int); lid = col("link_id").astype(int)
age = np.array([(s - born[l]) * dt if l in born else np.nan for l, s in zip(lid, col("step"))])
aa = np.degrees(np.arccos(np.clip(np.abs(col("direction_dot")), 0, 1)))
both = (col("geometric_detachment") > 0) & (col("stochastic_detachment") > 0)
# last saved frame of each dying link
want = {int(l): int(s) for l, s in zip(lid, col("step"))}
fsteps = f["state/current_step"][:, 0]; nF = int(f["meta/committed_frames"][0])
off = f["catch_bond/local_link_frame_offsets"]; ident = f["catch_bond/local_link_identity"]; mech = f["catch_bond/local_link_mechanics"]
U = lambda x: x / np.linalg.norm(x, axis=-1, keepdims=True)
last = {}
for t in range(int(np.searchsorted(fsteps[:nF], S0)) - 1, nF):
    o, c = int(off[t][1]), int(off[t][2]); idr = ident[o:o + c]; mr = mech[o:o + c]
    ids = idr[:, 1].astype(int); sel = np.array([want.get(k, -1) >= fsteps[t] for k in ids])
    if not sel.any(): continue
    u = U(f["actin/direction"][t]); i, j = idr[sel, 2].astype(int), idr[sel, 3].astype(int)
    n = U(u[i] - u[j] * np.sign(np.sum(u[i] * u[j], 1))[:, None])      # sign-aligned mean axis
    cv = mr[sel][:, [MC["connector_x"], MC["connector_y"], MC["connector_z"]]]
    ax = np.abs(np.sum(cv * n, 1)); rad = np.linalg.norm(cv - np.sum(cv * n, 1)[:, None] * n, axis=1)
    co = f["catch_bond/local_link_coordinates"][o:o + c][sel]
    for k, A_, R_, L_, (si, sj) in zip(ids[sel], ax, rad, mr[sel][:, MC["length"]], co):
        last[int(k)] = (A_, R_, L_, (want[int(k)] - fsteps[t]) * dt, si, sj)
geo = np.array([last.get(int(l), (np.nan,) * 6) for l in lid])
out = dict(source=P, window_s=[S0 * dt, S0 * dt + span_s], rupture_length=float(a["aa_link_rupture_length"]),
           rest_centerline=float(a["aa_preferred_centerline"]), reason_names=RN, by_status={})
print("detachments %d in %.2f s (%.0f per s); rupture length %.3f um, centre-line rest %.3f um"
      % (len(d), span_s, len(d) / span_s, a["aa_link_rupture_length"], a["aa_preferred_centerline"]))
for sname, ssel in (("all", np.ones(len(d), bool)), ("status 1", st == 1), ("status 2", st == 2)):
    rows = []
    for code in (2, 1, 3, 4):
        s = ssel & (rsn == code)
        if not s.any(): continue
        g = geo[s]
        rows.append(dict(reason=RN[str(code)], n=int(s.sum()), share=float(s.sum() / ssel.sum()), rate_per_s=float(s.sum() / span_s),
                         age_med_ms=float(1000 * np.nanmedian(age[s])), age_mean_ms=float(1000 * np.nanmean(age[s])),
                         load_zero=float(np.mean(col("productive_support")[s] <= 0)), load_mean=float(np.mean(col("productive_support")[s])),
                         length_med_nm=float(1000 * np.median(col("connector_length")[s])), aa_med_deg=float(np.median(aa[s])),
                         registry_med_nm=float(1000 * np.nanmedian(col("registry_error")[s])), both_tests=float(np.mean(both[s])),
                         last_axial_med_nm=float(1000 * np.nanmedian(g[:, 0])), last_radial_med_nm=float(1000 * np.nanmedian(g[:, 1])),
                         last_frame_lag_med_ms=float(1000 * np.nanmedian(g[:, 3]))))
    out["by_status"][sname] = dict(n=int(ssel.sum()), rows=rows)
    print("\n%s: %d detachments" % (sname, ssel.sum()))
    for r in rows:
        print("  %-32s %6d  %4.1f%%  %6.0f/s | age median %5.1f ms (mean %5.1f) | load 0 in %3.0f%%, mean %.2f pN | length %5.1f nm | A-A %4.1f deg | registry %5.1f nm"
              " | last frame: axial %5.1f nm, radial %5.1f nm (%.1f ms before) | contour+stochastic same step %.1f%%"
              % (r["reason"], r["n"], 100 * r["share"], r["rate_per_s"], r["age_med_ms"], r["age_mean_ms"], 100 * r["load_zero"], r["load_mean"],
                 r["length_med_nm"], r["aa_med_deg"], r["registry_med_nm"], r["last_axial_med_nm"], r["last_radial_med_nm"], r["last_frame_lag_med_ms"], 100 * r["both_tests"]))
# how abrupt is the stretch: last saved length vs length at death, and site movement since that frame
ds = np.abs(geo[:, 4] - col("s_i")) + np.abs(geo[:, 5] - col("s_j"))
out["abrupt"] = {}
print("\nlast saved frame (<= 5 ms before death) vs the death step:")
for sname, ssel in (("status 1", st == 1), ("status 2", st == 2)):
    for code in (1, 2):
        q = ssel & (rsn == code); h = q & ~np.isnan(geo[:, 2])
        L = 1000 * geo[h, 2]; pr = np.percentile(L, [10, 50, 90])
        out["abrupt"]["%s %s" % (sname, RN[str(code)])] = dict(n=int(q.sum()), with_frame=float(h.sum() / q.sum()), last_len_p10_p50_p90_nm=pr.tolist(),
            frac_last_below_100nm=float(np.mean(L < 100)), site_shift_med_nm=float(1000 * np.nanmedian(ds[h])))
        print("  %-8s %-18s n %6d, %3.0f%% have a saved frame | last length p10/50/90 %5.1f / %5.1f / %5.1f nm, below 100 nm %3.0f%% | sites moved %4.1f nm (median, sum of both)"
              % (sname, RN[str(code)], q.sum(), 100 * h.sum() / q.sum(), *pr, 100 * np.mean(L < 100), 1000 * np.nanmedian(ds[h])))
# survivors for reference: live-link geometry in the last frame
o, c = int(off[nF - 1][1]), int(off[nF - 1][2]); mr = mech[o:o + c]
out["live_length_med_nm"] = float(1000 * np.median(mr[:, MC["length"]]))
print("\nlive links at the last frame: connector length median %.1f nm" % out["live_length_med_nm"])
json.dump(out, open("/scratch/users/xsli/sarcomere/gaphunt_20261023/torque_out/actn2_breakage.json", "w"))
