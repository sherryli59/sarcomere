"""ACTN2 lifetime against A-M (and A-A) angle for resume arms, 2.0 -> 2.4 s.
Usage: lifetime_am_arm.py ARM[,ARM...]   (files data/<ARM>_seed20260805.h5)
 1. Measured: links formed after the resume (step > 200000); lifetime = detachment - formation
    (events v3), Kaplan-Meier restricted mean to 300 ms (the window is 400 ms), binned by the
    A-M angle at formation (larger of the two actins' angles to the sign-aligned axis sum of the
    myosins holding strong heads on them, nearest saved frame; both actins must hold one) and
    by the A-A angle at formation.
 2. Law-side: mature link-frames in frames 401-480, effective lifetime 1/<1/lifetime> using the
    lifetime the arm's law assigned (catch_bond/local_link_mechanics), and the unloaded share,
    by A-M angle."""
import h5py, numpy as np, os, json, sys
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
D = "/scratch/users/xsli/sarcomere/gaphunt_20261023/"
U = lambda d: d / np.linalg.norm(d, axis=-1, keepdims=True)
def am_angles(f, t, u):
    m = U(f["myosin/direction"][t]); hs = f["state/am_head_state_current"][t]; hb = f["state/am_head_bound_actin_current"][t].astype(int)
    mi, hi = np.nonzero(hs == 2); ai = hb[mi, hi]
    v = np.zeros_like(u); np.add.at(v, ai, m[mi] * np.sign(np.sum(m[mi] * u[ai], 1))[:, None])
    nv = np.linalg.norm(v, axis=1); th = np.full(len(u), np.nan); ok = nv > 1e-12
    th[ok] = np.degrees(np.arccos(np.clip(np.abs(np.sum(u[ok] * v[ok], 1)) / nv[ok], 0, 1))); return th
def km_mean(t, c, tmax=0.3):
    o = np.argsort(t); t, c = t[o], c[o]; n = len(t); S = 1.0; xs = [0.0]; ys = [1.0]
    for k, (ti, ci) in enumerate(zip(t, c)):
        if ci == 0: S *= 1 - 1.0 / (n - k)
        xs.append(ti); ys.append(S)
    xs, ys = np.array(xs), np.array(ys); g = np.linspace(0, tmax, 301)
    Sg = ys[np.searchsorted(xs, g, side="right") - 1]; return float(np.trapezoid(Sg, g))
out = {}
for arm in sys.argv[1].split(","):
    f = h5py.File(D + "data/%s_seed20260805.h5" % arm, "r", locking=False); a = dict(f["meta"].attrs); dt = float(a["dt"])
    EC = {c: i for i, c in enumerate(a["aa_local_link_events_v3_columns"].decode().split(","))}
    MC = {c: i for i, c in enumerate(a["aa_local_link_mechanics_columns"].decode().split(","))}
    nF = int(f["meta/committed_frames"][0]); fs = f["state/current_step"][:nF, 0]; T0 = int(np.searchsorted(fs, 200000))
    ev = f["catch_bond/local_link_events_v3"][:]; ev = ev[ev[:, EC["step"]] > 200000]; s_end = fs[nF - 1]
    AM = {t: am_angles(f, t, U(f["actin/direction"][t])) for t in range(T0, nF)}
    form, det = {}, {}
    for r in ev:
        lid = int(r[EC["link_id"]]); o_, n_ = int(r[EC["old_status"]]), int(r[EC["new_status"]])
        if o_ == 0 and n_ >= 1 and lid not in form:
            form[lid] = (r[EC["step"]], int(r[EC["actin_i"]]), int(r[EC["actin_j"]]), float(np.degrees(np.arccos(min(1.0, abs(r[EC["direction_dot"]]))))))
        elif n_ == 0 and lid in form and lid not in det: det[lid] = (r[EC["step"]], o_)
    rec = []
    for lid, (s0, i, j, aa0) in form.items():
        th = AM[int(np.clip(np.searchsorted(fs, s0), T0, nF - 1))]
        am0 = max(th[i], th[j]) if not (np.isnan(th[i]) or np.isnan(th[j])) else np.nan
        s1, cens, last = (det[lid][0], 0, det[lid][1]) if lid in det else (s_end, 1, -1)
        rec.append(((s1 - s0) * dt, cens, aa0, am0, last))
    R = np.array(rec, float); life, cens, aa0, am0, last = R.T
    res = dict(n_links=int(len(R)), measured_am=[], measured_aa=[], law_am=[])
    for key, x, bins in (("measured_am", am0, [(0, 10), (10, 20), (20, 40), (40, 90.01)]), ("measured_aa", aa0, [(0, 5), (5, 10), (10, 20), (20, 40), (40, 90.01)])):
        for lo, hi in bins:
            q = (x >= lo) & (x < hi)
            if q.sum() < 50: continue
            res[key].append(dict(lo=lo, hi=hi, n=int(q.sum()), mean_ms=1000 * km_mean(life[q], cens[q]),
                                 matured=float(np.mean(last[q & (cens == 0)] == 2)) if (q & (cens == 0)).any() else None))
    # law side
    off = f["catch_bond/local_link_frame_offsets"]; idn = f["catch_bond/local_link_identity"]; mech = f["catch_bond/local_link_mechanics"]
    LT, LD, AMf = [], [], []
    for t in range(T0 + 1, nF):
        o, c = int(off[t][1]), int(off[t][2]); idr = idn[o:o + c]; mr = mech[o:o + c]; s2 = idr[:, 4] == 2
        th = AM[t]; i, j = idr[s2, 2].astype(int), idr[s2, 3].astype(int)
        both = ~np.isnan(th[i]) & ~np.isnan(th[j])
        LT.append(mr[s2, MC["lifetime"]]); LD.append(mr[s2, MC["productive_support"]]); AMf.append(np.where(both, np.fmax(th[i], th[j]), np.nan))
    LT, LD, AMf = map(np.concatenate, (LT, LD, AMf))
    for lo, hi in ((0, 10), (10, 20), (20, 40), (40, 60), (60, 90.01)):
        q = (AMf >= lo) & (AMf < hi)
        res["law_am"].append(dict(lo=lo, hi=hi, n=int(q.sum()), eff_s=float(1 / np.mean(1 / LT[q])), unloaded=float(np.mean(LD[q] <= 0))))
    q = np.isnan(AMf); res["law_no_strong"] = dict(n=int(q.sum()), eff_s=float(1 / np.mean(1 / LT[q])))
    out[arm] = res
    print("\n== %s: %d links formed after 2.0 s" % (arm, len(R)))
    print("  measured lifetime by A-M at formation (KM mean to 300 ms): " + "  ".join("%d-%d: %.0f ms (n %d, matured %.0f%%)" % (b["lo"], min(90, b["hi"]), b["mean_ms"], b["n"], 100 * (b["matured"] or 0)) for b in res["measured_am"]))
    print("  measured lifetime by A-A at formation: " + "  ".join("%d-%d: %.0f ms" % (b["lo"], min(90, b["hi"]), b["mean_ms"]) for b in res["measured_aa"]))
    print("  law: mature effective lifetime by A-M: " + "  ".join("%d-%d: %.2f s (unloaded %.0f%%)" % (b["lo"], min(90, b["hi"]), b["eff_s"], 100 * b["unloaded"]) for b in res["law_am"]) + " | no strong head %.2f s" % res["law_no_strong"]["eff_s"])
json.dump(out, open(D + "torque_out/lifetime_am_arm.json", "w"))
