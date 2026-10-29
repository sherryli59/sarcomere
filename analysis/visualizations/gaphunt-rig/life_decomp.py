"""Take the measured ACTN2 lifetime apart, straight from catch_bond/local_link_events_v3.
For every link formed after S0 (formation = 0->1 event):
  phase 1  = formation -> maturation (1->2) or -> death as status 1 (1->0)
  phase 2  = maturation -> death (2->0), right-censored at the end of the run
A-M angle (larger of the two actins' angles to their strong-head myosins' axis, both actins
must hold one) is taken at formation for phase 1 and AT MATURATION for phase 2, from the
nearest following saved frame.  Means are Kaplan-Meier restricted to TMAX.
Also checks that link ids are never reused (else formation/death pairing is wrong)."""
import h5py, numpy as np, os, sys, json, collections
from detachment_reasons import reason_names
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
D = "/scratch/users/xsli/sarcomere/gaphunt_20261023/"
arm = sys.argv[1]; S0 = int(sys.argv[2]); TMAX = float(sys.argv[3])
f = h5py.File(D + "data/%s_seed20260805.h5" % arm, "r", locking=False); a = dict(f["meta"].attrs); dt = float(a["dt"])
EC = {c: i for i, c in enumerate(a["aa_local_link_events_v3_columns"].decode().split(","))}
RN = reason_names(a)
nF = int(f["meta/committed_frames"][0]); fs = f["state/current_step"][:nF, 0]; T0 = int(np.searchsorted(fs, S0))
U = lambda d: d / np.linalg.norm(d, axis=-1, keepdims=True)
cache = {}
def am(t):
    if t not in cache:
        u = U(f["actin/direction"][t]); m = U(f["myosin/direction"][t]); hs = f["state/am_head_state_current"][t]; hb = f["state/am_head_bound_actin_current"][t].astype(int)
        mi, hi = np.nonzero(hs == 2); ai = hb[mi, hi]
        v = np.zeros_like(u); np.add.at(v, ai, m[mi] * np.sign(np.sum(m[mi] * u[ai], 1))[:, None])
        nv = np.linalg.norm(v, axis=1); th = np.full(len(u), np.nan); ok = nv > 1e-12
        th[ok] = np.degrees(np.arccos(np.clip(np.abs(np.sum(u[ok] * v[ok], 1)) / nv[ok], 0, 1))); cache[t] = th
    return cache[t]
def pair_am(step, i, j):
    th = am(int(np.clip(np.searchsorted(fs, step), T0, nF - 1)))
    return max(th[i], th[j]) if not (np.isnan(th[i]) or np.isnan(th[j])) else np.nan
ev = f["catch_bond/local_link_events_v3"][:]
allform = ev[(ev[:, EC["old_status"]] == 0) & (ev[:, EC["new_status"]] == 1)]
cnt = collections.Counter(allform[:, EC["link_id"]].astype(int))
print("%s: formation events %d, distinct link ids %d, ids formed more than once %d" % (arm, len(allform), len(cnt), sum(1 for v in cnt.values() if v > 1)))
ev = ev[ev[:, EC["step"]] > S0]; s_end = float(fs[nF - 1])
L = {}
for r in ev:
    lid = int(r[EC["link_id"]]); o_, n_, s = int(r[EC["old_status"]]), int(r[EC["new_status"]]), r[EC["step"]]
    if o_ == 0 and n_ == 1: L[lid] = dict(f=s, i=int(r[EC["actin_i"]]), j=int(r[EC["actin_j"]]))
    elif lid in L:
        if o_ == 1 and n_ == 2: L[lid]["m"] = s
        elif n_ == 0: L[lid]["d"] = s; L[lid]["r"] = int(r[EC["detachment_reason"]]); L[lid]["ds"] = o_
def km(t, c):
    if len(t) == 0: return float("nan"), float("nan")
    o = np.argsort(t); t, c = t[o], c[o]; n = len(t); S = 1.0; xs = [0.0]; ys = [1.0]
    for k, (ti, ci) in enumerate(zip(t, c)):
        if ci == 0: S *= 1 - 1.0 / (n - k)
        xs.append(ti); ys.append(S)
    xs, ys = np.array(xs), np.array(ys); g = np.linspace(0, TMAX, 401)
    Sg = ys[np.searchsorted(xs, g, side="right") - 1]
    med = float(xs[np.argmax(ys <= 0.5)]) if (ys <= 0.5).any() else np.nan
    return float(np.trapezoid(Sg, g)), med
p1, p2 = [], []
for lid, x in L.items():
    a_f = pair_am(x["f"], x["i"], x["j"])
    if "m" in x:
        p1.append(((x["m"] - x["f"]) * dt, 0, 1, a_f, 0))
        a_m = pair_am(x["m"], x["i"], x["j"])
        if "d" in x: p2.append(((x["d"] - x["m"]) * dt, 0, a_m, x["r"]))
        else: p2.append(((s_end - x["m"]) * dt, 1, a_m, 0))
    elif "d" in x: p1.append(((x["d"] - x["f"]) * dt, 0, 0, a_f, x["r"]))
    else: p1.append(((s_end - x["f"]) * dt, 1, 0, a_f, 0))
p1, p2 = np.array(p1, float).reshape(-1, 5), np.array(p2, float).reshape(-1, 4)
BINS = [(0, 10), (10, 20), (20, 40), (40, 60), (60, 90.01)]
print("window %.2f-%.2f s; links formed %d; restricted means to %.0f ms" % (S0 * dt, s_end * dt, len(L), 1000 * TMAX))
print("\nPHASE 1 (new link): by A-M at formation")
for lo, hi in BINS + [(None, None)]:
    q = np.isnan(p1[:, 3]) if lo is None else (p1[:, 3] >= lo) & (p1[:, 3] < hi)
    if q.sum() < 30: continue
    mat = p1[q, 2] == 1; died = (p1[q, 2] == 0) & (p1[q, 1] == 0); rr = p1[q][died][:, 4].astype(int)
    print("  %-14s n %6d | matured %4.1f%% (median time to mature %5.1f ms) | died as status 1 %4.1f%% (contour %3.0f%%, random %3.0f%%, burial %3.0f%%) | still status 1 %4.1f%%"
          % ("no strong head" if lo is None else "%d-%d deg" % (lo, min(90, hi)), q.sum(), 100 * mat.mean(), 1000 * np.median(p1[q][mat][:, 0]) if mat.any() else np.nan,
             100 * died.mean(), 100 * np.mean(rr == 1), 100 * np.mean(rr == 2), 100 * np.mean(rr == 3), 100 * np.mean((p1[q, 2] == 0) & (p1[q, 1] == 1))))
print("\nPHASE 2 (mature link): by A-M at maturation")
for lo, hi in BINS + [(None, None)]:
    q = np.isnan(p2[:, 2]) if lo is None else (p2[:, 2] >= lo) & (p2[:, 2] < hi)
    if q.sum() < 30: continue
    mean, med = km(p2[q, 0], p2[q, 1]); d = (p2[q, 1] == 0); rr = p2[q][d][:, 3].astype(int)
    print("  %-14s n %6d | mature life: KM mean %5.1f ms, median %5.1f ms | censored %4.1f%% | deaths: %s"
          % ("no strong head" if lo is None else "%d-%d deg" % (lo, min(90, hi)), q.sum(), 1000 * mean, 1000 * med, 100 * (1 - d.mean()),
             ", ".join("%s %.1f%%" % (RN[str(k)],100*np.mean(rr==k) if len(rr) else 0) for k in (1,2,3,4,5))))
mean, med = km(p2[:, 0], p2[:, 1]); print("  all mature     n %6d | KM mean %5.1f ms, median %5.1f ms" % (len(p2), 1000 * mean, 1000 * med))
