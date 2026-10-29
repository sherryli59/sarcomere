"""ACTN2 link load vs A-A and A-M angle, production control r27ctrl, 2.0 -> 2.4 s every 25 ms.
Per live link and frame (catch_bond/local_link_mechanics, rows aligned with local_link_identity):
  load      = productive_support (the productive opposed load; the effective load in the
              lifetime law for productive_register mode)
  lifetime  = the lifetime the law assigned that step
  A-A angle = unsigned angle between the two linked actins
  A-M angle = each actin's angle to the sign-aligned sum of the axes of the myosins holding
              strong heads on it; pair value reported as the larger and as the smaller
The law is applied every step as a hazard: P(break) = dt / lifetime with the lifetime of that step
(catch_bond_kinetics.h break_probability), so the lifetime a link actually experiences is the
harmonic mean of its per-step lifetimes, 1 / <1/lifetime>, reported per bin as 'effective lifetime'.
Also: detachment reasons and the load at detachment (events v3, links formed >= 1.0 s), and the
40-step angle-lifetime build check (anglelife_20261027/check): link counts per frame."""
import h5py, numpy as np, os, json
from detachment_reasons import reason_names
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
P = "/scratch/users/xsli/sarcomere/gaphunt_20261023/data/r27ctrl_seed20260805.h5"
f = h5py.File(P, "r", locking=False); a = dict(f["meta"].attrs)
mcols = a["aa_local_link_mechanics_columns"].decode().split(","); MC = {c: i for i, c in enumerate(mcols)}
U = lambda d: d / np.linalg.norm(d, axis=-1, keepdims=True)
off = f["catch_bond/local_link_frame_offsets"]; ident = f["catch_bond/local_link_identity"]; mech = f["catch_bond/local_link_mechanics"]
assert mech.shape[0] == ident.shape[0], "mechanics rows are not aligned with identity rows"
def am_angles(t, u):
    m = U(f["myosin/direction"][t]); hs = f["state/am_head_state_current"][t]; hb = f["state/am_head_bound_actin_current"][t].astype(int)
    mi, hi = np.nonzero(hs == 2); ai = hb[mi, hi]
    v = np.zeros_like(u); np.add.at(v, ai, m[mi] * np.sign(np.sum(m[mi] * u[ai], 1))[:, None])
    nv = np.linalg.norm(v, axis=1); th = np.full(len(u), np.nan); ok = nv > 1e-12
    th[ok] = np.degrees(np.arccos(np.clip(np.abs(np.sum(u[ok] * v[ok], 1)) / nv[ok], 0, 1))); return th
rows = []
for t in range(400, 481, 5):
    o, c = int(off[t][1]), int(off[t][2]); idr = ident[o:o + c]; mr = mech[o:o + c]
    u = U(f["actin/direction"][t]); i, j = idr[:, 2].astype(int), idr[:, 3].astype(int); st = idr[:, 4].astype(int)
    aa = np.degrees(np.arccos(np.clip(np.abs(np.sum(u[i] * u[j], 1)), 0, 1)))
    th = am_angles(t, u); ami, amj = th[i], th[j]
    with np.errstate(invalid="ignore"):
        am_max = np.fmax(ami, amj); am_min = np.fmin(ami, amj)
    both = (~np.isnan(ami)) & (~np.isnan(amj))
    for k in range(c):
        rows.append((aa[k], am_max[k], am_min[k], float(both[k]), st[k], mr[k, MC["productive_support"]], mr[k, MC["lifetime"]], mr[k, MC["bilateral_load"]]))
R = np.array(rows, float); aa, amx, amn, both, st, load, life, bil = R.T
print("link-frames %d; status-2 share %.0f%%; load > 0 in %.0f%% (status-1 %.0f%%, status-2 %.0f%%); median load where > 0: %.2f pN"
      % (len(R), 100 * np.mean(st == 2), 100 * np.mean(load > 0), 100 * np.mean(load[st == 1] > 0), 100 * np.mean(load[st == 2] > 0), np.median(load[load > 0])))
# law check: status-2 lifetime = 0.15 + 0.4 * load
s2 = st == 2; pred = 0.1 + 0.05 + 0.4 * load[s2]
print("lifetime law check (status 2): max |recorded - (0.15 + 0.4 load)| = %.2e s; status-1 lifetimes unique: %s" % (np.max(np.abs(life[s2] - pred)), np.unique(np.round(life[st == 1], 6))[:5]))
out = dict(n=int(len(R)), groups={})
def binned(x, edges, sel_extra, lab):
    res = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        for sname, ssel in (("status 1", st == 1), ("status 2", st == 2)):
            s = (x >= lo) & (x < hi) & ssel & sel_extra
            if s.sum() < 30: continue
            res.append(dict(lo=float(lo), hi=float(hi), status=sname, n=int(s.sum()), mean_load=float(load[s].mean()),
                            median_load=float(np.median(load[s])), frac_loaded=float(np.mean(load[s] > 0)),
                            mean_load_if_loaded=float(load[s][load[s] > 0].mean()) if (load[s] > 0).any() else 0.0,
                            mean_lifetime=float(life[s].mean()), eff_lifetime=float(1.0 / np.mean(1.0 / life[s]))))
    out["groups"][lab] = res
    print("\n%s" % lab)
    for r in res:
        print("   %4.0f-%4.0f deg %-8s n %6d | load mean %.2f pN (median %.2f), loaded %3.0f%%, mean load when loaded %.2f pN | assigned lifetime mean %.3f s, effective %.3f s"
              % (r["lo"], r["hi"], r["status"], r["n"], r["mean_load"], r["median_load"], 100 * r["frac_loaded"], r["mean_load_if_loaded"], r["mean_lifetime"], r["eff_lifetime"]))
edges = np.array([0, 5, 10, 20, 30, 40, 60, 90.01])
binned(aa, edges, np.ones(len(R), bool), "A-A angle")
binned(amx, edges, both > 0, "A-M angle, larger of the two actins (both strong-bound)")
binned(amn, edges, both > 0, "A-M angle, smaller of the two actins (both strong-bound)")
s = both == 0
print("\nat least one actin with no strong head: n %d, load > 0 in %.0f%%, mean load %.2f pN" % (s.sum(), 100 * np.mean(load[s] > 0), load[s].mean()))
out["no_strong"] = dict(n=int(s.sum()), frac_loaded=float(np.mean(load[s] > 0)), mean_load=float(load[s].mean()))
# fine 5-deg curves for the page (status 2, all)
for lab, x, sel in (("aa", aa, np.ones(len(R), bool)), ("am_max", amx, both > 0), ("am_min", amn, both > 0)):
    pts = []
    for lo in range(0, 90, 5):
        for sname, ssel in (("s1", st == 1), ("s2", st == 2)):
            q = (x >= lo) & (x < lo + 5) & sel & ssel
            if q.sum() >= 50: pts.append(dict(mid=lo + 2.5, status=sname, n=int(q.sum()), mean_load=float(load[q].mean()), frac_loaded=float(np.mean(load[q] > 0)),
                                                 mean_lifetime=float(life[q].mean()), eff_lifetime=float(1.0 / np.mean(1.0 / life[q]))))
    out["fine_" + lab] = pts
out["overall"] = dict(frac_s2=float(np.mean(st == 2)), frac_loaded=float(np.mean(load > 0)), frac_loaded_s1=float(np.mean(load[st == 1] > 0)),
                      frac_loaded_s2=float(np.mean(load[st == 2] > 0)), median_nonzero_load=float(np.median(load[load > 0])),
                      law_max_err_s2=float(np.max(np.abs(life[s2] - pred))), eff_lifetime_s2=float(1.0 / np.mean(1.0 / life[s2])),
                      mean_lifetime_s2=float(life[s2].mean()), frames=[400, 480, 5], source=P)
# detachments: reason and the lifetime the law had assigned at that step (0.15 s for status 2 = zero load)
ecols = a["aa_local_link_events_v3_columns"].decode().split(","); EC = {c: i for i, c in enumerate(ecols)}
reasons = reason_names(a)
ev = f["catch_bond/local_link_events_v3"][:]; ev = ev[(ev[:, EC["step"]] >= 100000) & (ev[:, EC["new_status"]] == 0)]
out["detach"] = {}
for sv in (1, 2):
    d = ev[ev[:, EC["old_status"]] == sv]; r = d[:, EC["reason"] if "reason" in EC else EC["detachment_reason"]].astype(int); lt = d[:, EC["lifetime"]]
    sto = r == 2
    out["detach"]["s%d" % sv] = dict(n=int(len(d)), reasons={reasons[str(k)]: float(np.mean(r == k)) for k in (1, 2, 3, 4, 5)},
                                     stochastic_at_zero_load=float(np.mean(lt[sto] <= (0.1 if sv == 1 else 0.15) + 1e-9)) if sto.any() else None,
                                     median_lifetime_at_detach=float(np.median(lt)))
    print("status-%d detachments %d: %s; stochastic ones at zero load %.0f%%" % (sv, len(d), {k: round(v, 3) for k, v in out["detach"]["s%d" % sv]["reasons"].items()},
          100 * (out["detach"]["s%d" % sv]["stochastic_at_zero_load"] or 0)))
# the 40-step angle-lifetime check
C = "/scratch/users/xsli/sarcomere/anglelife_20261027/check/"
print("\n40-step build check (resume at 2.0 s, 10-step frames, 1 thread):")
for tag in ("ref", "new_default", "new_on"):
    g = h5py.File(C + tag + ".h5", "r", locking=False); oo = g["catch_bond/local_link_frame_offsets"]
    cnt = [int(oo[t][2]) for t in range(400, 405)]
    idx = [g["catch_bond/local_link_identity"][int(oo[t][1]):int(oo[t][1]) + int(oo[t][2])][:, 4] for t in (400, 404)]
    print("   %-12s live links at steps 200000..200040: %s | status-2 %d -> %d" % (tag, cnt, int((idx[0] == 2).sum()), int((idx[1] == 2).sum())))
    # a resumed file keeps the source run's meta, so the knob is NOT in meta: read it from the launch script
    import re
    sb = open(C + "../build_check.sbatch").read(); m = re.search(r"run \$\w+ " + tag + r"((?: --\S+)*) &", sb)
    kn = dict(x.lstrip("-").split("=") for x in m.group(1).split()) if m else {}
    out.setdefault("check", {})[tag] = dict(live=cnt, s2=[int((idx[0] == 2).sum()), int((idx[1] == 2).sum())], steps=[200000, 200040],
                                            aa_scale_deg=float(kn.get("aa_lifetime_aa_angle_scale_deg", 0)),
                                            am_scale_deg=float(kn.get("aa_lifetime_am_angle_scale_deg", 0)),
                                            meta_records_knob="aa_lifetime_aa_angle_scale_deg" in g["meta"].attrs)
json.dump(out, open("/scratch/users/xsli/sarcomere/gaphunt_20261023/torque_out/load_angle.json", "w"))
