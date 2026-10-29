"""How much A-M contrast can any lifetime(load) expression give, with today's instantaneous load?
Each step a link breaks with probability dt / tau(L(t)), so the lifetime it lives with is
1 / <1 / tau(L)> over its frames (snapshots sample time without bias, so this is exact in
expectation for an instantaneous law).  Mature (status-2) link-frames, r27ctrl 2.0 -> 2.4 s,
binned by A-M angle of the more misaligned actin (both strong-bound).
Forms (tau(0) = 0.15 s in all):  linear tau0 + c L (today c = 0.4 s/pN);  exponential catch
tau0 exp(L / Fc);  and the ceiling tau0 / (unloaded fraction) that any form approaches as the
catch gets infinitely strong."""
import h5py, numpy as np, os, json
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
Ls, As = [], []
for t in range(400, 481):
    o, c = int(off[t][1]), int(off[t][2]); idr = ident[o:o + c]; mr = mech[o:o + c]; s2 = idr[:, 4] == 2
    u = U(f["actin/direction"][t]); th = am_angles(t, u); i, j = idr[s2, 2].astype(int), idr[s2, 3].astype(int)
    both = ~np.isnan(th[i]) & ~np.isnan(th[j]); Ls.append(mr[s2, MC["productive_support"]]); As.append(np.where(both, np.fmax(th[i], th[j]), np.nan))
L = np.concatenate(Ls); A = np.concatenate(As)
EDG = [(0, 10), (10, 20), (20, 40), (40, 60), (60, 90.01)]
sel = [(A >= lo) & (A < hi) for lo, hi in EDG] + [np.isnan(A)]
lab = ["0-10", "10-20", "20-40", "40-60", "60-90", "no strong head"]
forms = [("linear c 0.4 (today)", lambda x: 0.15 + 0.4 * x), ("linear c 2", lambda x: 0.15 + 2 * x), ("linear c 10", lambda x: 0.15 + 10 * x),
         ("exp Fc 2 pN", lambda x: 0.15 * np.exp(np.minimum(x / 2, 50))), ("exp Fc 1 pN", lambda x: 0.15 * np.exp(np.minimum(x / 1, 50))),
         ("exp Fc 0.5 pN", lambda x: 0.15 * np.exp(np.minimum(x / 0.5, 50))),
         ("base 10 ms, c 0.4", lambda x: 0.010 + 0.4 * x), ("base 1 ms, c 0.4", lambda x: 0.001 + 0.4 * x)]
out = dict(source=P, n=int(len(L)), bins=lab, unloaded=[float(np.mean(L[q] <= 0)) for q in sel], forms={})
print("mature link-frames %d; unloaded fraction by bin: %s" % (len(L), "  ".join("%s %.2f" % (l, u_) for l, u_ in zip(lab, out["unloaded"]))))
print("\neffective lifetime (s) by A-M bin | 0-10/20-40 | 0-10/60-90")
for name, fn in forms:
    tau = fn(L); eff = [float(1 / np.mean(1 / tau[q])) for q in sel]; out["forms"][name] = eff
    print("  %-22s %s | %.1fx | %.1fx" % (name, "  ".join("%s %.4f" % (l, e) for l, e in zip(lab, eff)), eff[0] / eff[2], eff[0] / eff[4]))
cap = [0.15 / max(u_, 1e-9) for u_ in out["unloaded"]]; out["forms"]["ceiling (infinitely strong catch)"] = cap
print("  %-22s %s | %.1fx | %.1fx" % ("ceiling", "  ".join("%s %.2f" % (l, e) for l, e in zip(lab, cap)), cap[0] / cap[2], cap[0] / cap[4]))
json.dump(out, open("/scratch/users/xsli/sarcomere/gaphunt_20261023/torque_out/load_form_est.json", "w"))
