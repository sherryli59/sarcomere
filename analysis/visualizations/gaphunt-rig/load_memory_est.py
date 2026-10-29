"""Design estimate for a stronger A-M selection through the ACTN2 catch-bond load.
Production control r27ctrl, every saved frame 2.0 -> 2.4 s (5 ms), mature (status-2) links.
1. Is zero load blinking or persistent?  Run lengths of zero-load stretches per link.
2. Offline estimate of the lifetime a mature link lives with, 1/<1/tau>, by A-M angle, if
   tau = 0.15 s + c * L_hold, with L_hold a peak-hold of the recorded load decaying with time
   constant tau_mem (tau_mem = 0 is today's instantaneous load):
       L_hold(t) = max(L(t), L_hold(t - 5 ms) * exp(-5 ms / tau_mem))
   Frame resolution (5 ms) cannot see sub-frame blinking; the tau_mem = 0, c = 0.4 row must
   reproduce the measured 0.64 / 0.22 s and is the check on the binning."""
import h5py, numpy as np, os, json, collections
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
P = "/scratch/users/xsli/sarcomere/gaphunt_20261023/data/r27ctrl_seed20260805.h5"
f = h5py.File(P, "r", locking=False); a = dict(f["meta"].attrs)
MC = {c: i for i, c in enumerate(a["aa_local_link_mechanics_columns"].decode().split(","))}
off = f["catch_bond/local_link_frame_offsets"]; ident = f["catch_bond/local_link_identity"]; mech = f["catch_bond/local_link_mechanics"]
U = lambda d: d / np.linalg.norm(d, axis=-1, keepdims=True)
def am_angles(t, u):
    m = U(f["myosin/direction"][t]); hs = f["state/am_head_state_current"][t]; hb = f["state/am_head_bound_actin_current"][t].astype(int)
    mi, hi = np.nonzero(hs == 2); ai = hb[mi, hi]
    v = np.zeros_like(u); np.add.at(v, ai, m[mi] * np.sign(np.sum(m[mi] * u[ai], 1))[:, None])
    nv = np.linalg.norm(v, axis=1); th = np.full(len(u), np.nan); ok = nv > 1e-12
    th[ok] = np.degrees(np.arccos(np.clip(np.abs(np.sum(u[ok] * v[ok], 1)) / nv[ok], 0, 1))); return th
series = collections.defaultdict(list)          # link id -> [(frame, status, load, am_max)]
for t in range(400, 481):
    o, c = int(off[t][1]), int(off[t][2]); idr = ident[o:o + c]; mr = mech[o:o + c]
    u = U(f["actin/direction"][t]); th = am_angles(t, u); i, j = idr[:, 2].astype(int), idr[:, 3].astype(int)
    both = ~np.isnan(th[i]) & ~np.isnan(th[j]); amx = np.where(both, np.fmax(th[i], th[j]), np.nan)
    for k, s_, L_, A_ in zip(idr[:, 1].astype(int), idr[:, 4].astype(int), mr[:, MC["productive_support"]], amx):
        series[k].append((t, s_, L_, A_))
# 1. zero-load run lengths for status-2 links (consecutive frames only)
runs = []; censored = 0
for k, s in series.items():
    s = [x for x in s if x[1] == 2]
    if len(s) < 2: continue
    run = 0; prev = None
    for (t, st, L, A) in s:
        if prev is not None and t != prev + 1:
            run = 0
        if L <= 0: run += 1
        else:
            if run: runs.append(run)
            run = 0
        prev = t
    if run: censored += 1
runs = np.array(runs)
print("mature links, zero-load stretches that ended within the window: n %d; length in 5 ms frames: 1 frame %.0f%%, 2 %.0f%%, 3-4 %.0f%%, 5+ %.0f%% (median %d); still unloaded at the end or at a gap: %d"
      % (len(runs), 100 * np.mean(runs == 1), 100 * np.mean(runs == 2), 100 * np.mean((runs >= 3) & (runs <= 4)), 100 * np.mean(runs >= 5), np.median(runs), censored))
# 2. offline estimate
EDG = [(0, 10), (10, 20), (20, 40), (40, 60), (60, 90.01)]
out = dict(source=P, runs=dict(n=int(len(runs)), frac1=float(np.mean(runs == 1)), frac2=float(np.mean(runs == 2)), frac5p=float(np.mean(runs >= 5))), grid=[])
print("\neffective lifetime of mature links (s) by A-M angle (more misaligned actin, both strong-bound) | links with an actin lacking strong heads | ratio 0-10 / 60-90")
for c in (0.4, 1.0, 2.0):
    for tm in (0.0, 0.010, 0.025, 0.050):
        tau_all, am_all = [], []
        for k, s in series.items():
            hold = 0.0; prev = None
            for (t, st, L, A) in s:
                if st != 2: prev = None; continue
                if prev is None or t != prev + 1 or tm == 0: hold = L
                else: hold = max(L, hold * np.exp(-0.005 / tm))
                prev = t; tau_all.append(0.15 + c * hold); am_all.append(A)
        tau_all, am_all = np.array(tau_all), np.array(am_all)
        row = dict(c=c, tau_mem=tm, bins=[], no_strong=float(1 / np.mean(1 / tau_all[np.isnan(am_all)])))
        for lo, hi in EDG:
            q = (am_all >= lo) & (am_all < hi); row["bins"].append(dict(lo=lo, hi=hi, n=int(q.sum()), eff=float(1 / np.mean(1 / tau_all[q]))))
        out["grid"].append(row)
        print("  c %.1f s/pN, memory %2.0f ms: %s | no strong head %.2f | ratio %.1f"
              % (c, 1000 * tm, "  ".join("%d-%d: %.2f" % (b["lo"], min(90, b["hi"]), b["eff"]) for b in row["bins"]), row["no_strong"], row["bins"][0]["eff"] / row["bins"][-1]["eff"]))
json.dump(out, open("/scratch/users/xsli/sarcomere/gaphunt_20261023/torque_out/load_memory_est.json", "w"))
