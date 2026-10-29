"""Measured ACTN2 link lifetime vs A-A and A-M angle, production control (hph27, seed 805):
the source run 1.19 -> 2.0 s plus its unchanged resume 2.0 -> 2.4 s (r27ctrl file, which holds both).

Lifetime = detachment step - formation step (catch_bond/local_link_events_v3: formation =
new_status 1 from 0, detachment = new_status 0).  Links formed before the log starts are
excluded; links alive at the last step are right-censored (Kaplan-Meier).
Angles:
  A-A at formation  = acos|direction_dot| recorded on the formation event
  A-A over life     = mean over the saved 5 ms frames while the link is present
  A-M at formation  = larger of the two actins' angles to the sign-aligned sum of the axes of
                      the myosins holding STRONG heads on them, in the nearest saved frame
                      (an actin without a strong head: NaN; the pair value uses the other)
Outputs JSON for the plot page."""
import h5py, numpy as np, os, json, collections
from detachment_reasons import reason_names
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
P = "/scratch/users/xsli/sarcomere/gaphunt_20261023/data/r27ctrl_seed20260805.h5"
f = h5py.File(P, "r", locking=False); a = dict(f["meta"].attrs); dt = float(a["dt"])
cols = a["aa_local_link_events_v3_columns"].decode().split(","); C = {c: i for i, c in enumerate(cols)}
reasons = reason_names(a)
ev = f["catch_bond/local_link_events_v3"][:]
step1 = ev[:, C["step"]].max(); step0 = 100000.0     # links formed after 1.0 s: the settled network, not assembly
ev = ev[ev[:, C["step"]] >= step0]
U = lambda d: d / np.linalg.norm(d, axis=-1, keepdims=True)
nF = int(f["meta/committed_frames"][0]); fsteps = f["state/current_step"][:nF, 0]
F0 = int(np.searchsorted(fsteps, step0))
# per-frame A-M angle of every actin
def am_angles(t):
    u = U(f["actin/direction"][t]); m = U(f["myosin/direction"][t])
    hs = f["state/am_head_state_current"][t]; hb = f["state/am_head_bound_actin_current"][t].astype(int)
    mi, hi = np.nonzero(hs == 2); ai = hb[mi, hi]
    v = np.zeros_like(u); np.add.at(v, ai, m[mi] * np.sign(np.sum(m[mi] * u[ai], 1))[:, None])
    nv = np.linalg.norm(v, axis=1); th = np.full(len(u), np.nan); ok = nv > 1e-12
    th[ok] = np.degrees(np.arccos(np.clip(np.abs(np.sum(u[ok] * v[ok], 1)) / nv[ok], 0, 1))); return th
AM = {t: am_angles(t) for t in range(F0, nF)}
# per-frame A-A angle of present links (by id) for the over-life average
off = f["catch_bond/local_link_frame_offsets"]; ident = f["catch_bond/local_link_identity"]
life_aa = collections.defaultdict(list)
for t in range(F0, nF):
    o, c = int(off[t][1]), int(off[t][2]); r = ident[o:o + c]
    u = U(f["actin/direction"][t]); i, j = r[:, 2].astype(int), r[:, 3].astype(int)
    angs = np.degrees(np.arccos(np.clip(np.abs(np.sum(u[i] * u[j], 1)), 0, 1)))
    for lid, g in zip(r[:, 1].astype(int), angs): life_aa[lid].append(g)
form = {}; det = {}
for row in ev:
    lid = int(row[C["link_id"]]); old, new = int(row[C["old_status"]]), int(row[C["new_status"]])
    if old == 0 and new >= 1 and lid not in form:
        form[lid] = (row[C["step"]], int(row[C["actin_i"]]), int(row[C["actin_j"]]), float(np.degrees(np.arccos(min(1.0, abs(row[C["direction_dot"]]))))))
    elif new == 0 and lid in form and lid not in det:
        det[lid] = (row[C["step"]], int(row[C["detachment_reason"]]), old)
recs = []
for lid, (s0, i, j, aa0) in form.items():
    fr = int(np.clip(np.searchsorted(fsteps, s0), F0, nF - 1))
    th = AM[fr]; am0 = np.nanmax([th[i], th[j]]) if not (np.isnan(th[i]) and np.isnan(th[j])) else np.nan
    if lid in det: s1, rsn, last_status = det[lid]; cens = 0
    else: s1, rsn, last_status = step1, 0, -1; cens = 1
    life = (s1 - s0) * dt
    recs.append((life, cens, aa0, am0, float(np.mean(life_aa[lid])) if life_aa.get(lid) else np.nan, rsn, last_status))
R = np.array(recs, float); life, cens, aa0, am0, aal, rsn, lst = R.T
print("links formed %d (from %.3f to %.3f s), detached %d, censored %d; median lifetime (detached) %.1f ms, mean %.1f ms"
      % (len(R), step0 * dt, step1 * dt, int((cens == 0).sum()), int(cens.sum()), 1000 * np.median(life[cens == 0]), 1000 * life[cens == 0].mean()))
def km(t, c):
    o = np.argsort(t); t, c = t[o], c[o]; n = len(t); S = 1.0; xs = [0.0]; ys = [1.0]
    for k, (ti, ci) in enumerate(zip(t, c)):
        if ci == 0: S *= 1 - 1.0 / (n - k)
        xs.append(ti); ys.append(S)
    xs, ys = np.array(xs), np.array(ys)
    grid = np.linspace(0, 0.6, 121); S_g = np.array([ys[np.searchsorted(xs, g, side="right") - 1] for g in grid])
    med = float(xs[np.argmax(ys <= 0.5)]) if (ys <= 0.5).any() else None
    mean = float(np.trapezoid(S_g, grid))           # restricted mean lifetime up to 0.6 s
    return grid.tolist(), S_g.tolist(), med, mean
BINS = {"aa0": [(0, 5), (5, 10), (10, 20), (20, 40), (40, 90.01)], "am0": [(0, 10), (10, 20), (20, 40), (40, 90.01)], "aal": [(0, 5), (5, 10), (10, 20), (20, 40), (40, 90.01)]}
VAL = {"aa0": aa0, "am0": am0, "aal": aal}
LAB = {"aa0": "A-A angle at formation", "am0": "A-M angle at formation", "aal": "A-A angle averaged over the link's life"}
out = dict(meta=dict(source=P, window_s=[step0 * dt, step1 * dt], n_links=int(len(R)), n_detached=int((cens == 0).sum()),
                     reason_codes=reasons), groups={})
for key, bins in BINS.items():
    g = []
    for lo, hi in bins:
        s = (VAL[key] >= lo) & (VAL[key] < hi)
        if s.sum() < 50: continue
        grid, S, med, mean = km(life[s], cens[s])
        # bootstrap CI on the restricted mean
        rng = np.random.default_rng(0); idx = np.flatnonzero(s); bs = []
        for _ in range(100):
            b = rng.choice(idx, len(idx)); bs.append(km(life[b], cens[b])[3])
        rs = rsn[s & (cens == 0)]
        g.append(dict(lo=lo, hi=hi, n=int(s.sum()), share=float(s.mean()), median_ms=1000 * med if med else None,
                      mean_ms=1000 * mean, mean_ci_ms=[1000 * float(np.percentile(bs, 2.5)), 1000 * float(np.percentile(bs, 97.5))],
                      matured=float(np.mean(lst[s & (cens == 0)] == 2)) if (s & (cens == 0)).any() else None,
                      reasons={reasons[str(k)]: float(np.mean(rs == k)) for k in (1, 2, 3, 4, 5)},
                      km_t=grid, km_S=[round(x, 4) for x in S]))
        print("%-40s %5.1f-%5.1f deg: n %6d (%4.1f%%)  median %5.1f ms  mean %5.1f ms [%5.1f, %5.1f]  matured %.0f%%  stochastic %.0f%% geometric %.0f%% burial %.0f%%"
              % (LAB[key], lo, min(hi, 90), s.sum(), 100 * s.mean(), g[-1]["median_ms"] or np.nan, g[-1]["mean_ms"], *g[-1]["mean_ci_ms"],
                 100 * (g[-1]["matured"] or 0), 100 * g[-1]["reasons"].get("stochastic", 0), 100 * g[-1]["reasons"].get("geometric_contour", 0),
                 100 * g[-1]["reasons"].get("barbed_tip_capsule_penetration", 0)))
    out["groups"][key] = dict(label=LAB[key], bins=g)
# fine-grained scatter-ready curve: restricted mean vs angle in 5-deg bins for both angles
for key in ("aa0", "am0"):
    pts = []
    for lo in range(0, 90, 5):
        s = (VAL[key] >= lo) & (VAL[key] < lo + 5)
        if s.sum() < 100: continue
        pts.append(dict(mid=lo + 2.5, n=int(s.sum()), mean_ms=1000 * km(life[s], cens[s])[3]))
    out["groups"][key]["fine"] = pts
json.dump(out, open("/scratch/users/xsli/sarcomere/gaphunt_20261023/torque_out/lifetime_angle.json", "w"))
print("written")
