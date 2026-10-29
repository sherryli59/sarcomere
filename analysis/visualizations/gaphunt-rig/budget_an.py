"""Data for the torque-budget page; every number the page shows comes from here.
1) per span bin: mean drive toward the bound partner (kT) by source, production
   hph27/54/81 (last 0.5 s, seed 805/815/805) vs the sticky rig tq30 x3 seeds and
   the production-kinetics rig k160_30 x3 (aa_steric-singular frames dropped);
2) per actin: drive toward the MYOSIN DIRECTOR by source (last 50 frames, every 5th);
3) torque-cap and rotation-cap event rates per actin-step over the last 0.5 s;
4) angle vs time for the head-kinetics rig ladder (every 4th frame = 2 ms)."""
import h5py, numpy as np, json, os, sys
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import torque_an as T
D = "/scratch/users/xsli/sarcomere/gaphunt_20261023/torque_out"; C = T.CAMP + "/"
KT = T.KT
BINS = [(0, 1, "one crown (0 nm)"), (1, 200, "1-200 nm"), (200, 400, "200-400 nm"), (400, 5000, "400 nm or more")]
ONLY_KIN = "--kin-only" in sys.argv
out = json.load(open("%s/budget.json" % D)) if ONLY_KIN else dict(span=[], director=[], caps=[], kin=[])
if ONLY_KIN: out["kin"] = []
# rigs at production kinetics need the npz too
for arm in (() if ONLY_KIN else ("k160_30",)):
    for s in ("20260805", "20260815", "20260825"):
        tag = "%s_seed%s" % (arm, s)
        if not os.path.exists("%s/torque_%s.npz" % (D, tag)):
            rows, v = T.run(C + tag + ".h5", range(0, 2001, 2))
            np.savez_compressed("%s/torque_%s.npz" % (D, tag), **{k: np.array([r[k] for r in rows]) for k in rows[0]})
GROUPS = [("hph27", ["hph27"], "production, 27 heads/half"), ("hph54", ["hph54"], "production, 54 heads/half"),
          ("hph81", ["hph81"], "production, 81 heads/half"),
          ("rig_sticky", ["tq30_seed20260805", "tq30_seed20260815", "tq30_seed20260825"], "1 actin + 1 myosin, sticky heads"),
          ("rig_k160", ["k160_30_seed20260805", "k160_30_seed20260815", "k160_30_seed20260825"], "1 actin + 1 myosin, production head kinetics")]
for key, tags, label in ([] if ONLY_KIN else GROUPS):
    z = {}
    for tg in tags:
        d = dict(np.load("%s/torque_%s.npz" % (D, tg)))
        for k, v in d.items(): z.setdefault(k, []).append(v)
    z = {k: np.concatenate(v) for k, v in z.items()}
    ok = (np.abs(z["aa_steric"]) < 5000) & (np.abs(z["TOT"]) < 5000) & (np.abs(z["MTOT"]) < 5000)
    g = dict(key=key, label=label, n=int(ok.sum()), bins=[])
    for lo, hi, bl in BINS:
        s = ok & (z["span"] >= lo) & (z["span"] < hi) & (z["th"] >= 5) & (z["th"] < 60)
        if s.sum() < 20: g["bins"].append(dict(bin=bl, n=int(s.sum()))); continue
        e = lambda x: (float(x.mean()), float(x.std() / np.sqrt(len(x))))
        actn2 = z["actn2_align"] + z["actn2_lever"] + z["actn2_perp"] + z["actn2_axial"] + z["actn2_contact"]
        g["bins"].append(dict(bin=bl, n=int(s.sum()), share=float(s.sum() / (ok & (z["th"] >= 5) & (z["th"] < 60)).sum()),
            pair_radial=e(z["RAD"][s]), heads_all=e(z["am_head"][s]), lateral=e(z["aa_lateral"][s]), actn2=e(actn2[s]),
            am_steric=e(z["am_steric"][s]), total=e(z["TOT"][s]), myosin_total=e(z["MTOT"][s]),
            w_actin=e(z["wa"][s]), pred_actin=float(float(np.mean(z["TOT"][s])) * 0.10096018660404828),
            heads=float(z["n"][s].mean())))
    out["span"].append(g)
SRC = [("am_head", "all myosin heads"), ("aa_lateral", "barbed-end lateral repulsion"), ("actn2_align", "ACTN2 angular (kappa_aa)"),
       ("actn2_lever", "ACTN2 lever"), ("actn2_axial", "ACTN2 axial register"), ("am_steric", "actin-myosin sterics")]
for tag, fn in (() if ONLY_KIN else (("hph27", "hph27m24_seed20260805.h5"), ("hph54", "hph54m24_seed20260815.h5"), ("hph81", "hph81m24_seed20260805.h5"))):
    f = h5py.File(C + fn, "r", locking=False); n = int(f["meta/committed_frames"][0])
    acc = {k: [] for k, _ in SRC}; th = []; lat = []
    for t in range(n - 50, n, 5):
        u = f["actin/direction"][t]; u = u / np.linalg.norm(u, axis=1, keepdims=True)
        m = f["myosin/direction"][t]; m = m / np.linalg.norm(m, axis=1, keepdims=True)
        w, v = np.linalg.eigh((m[:, :, None] * m[:, None, :]).mean(0)); n0 = v[:, -1]
        nm = np.sign(u @ n0)[:, None] * n0[None, :]; tg = nm - u * np.sum(u * nm, 1)[:, None]; tg /= np.maximum(np.linalg.norm(tg, axis=1, keepdims=True), 1e-12)
        th.append(np.degrees(np.arccos(np.clip(np.abs(u @ n0), 0, 1))))
        for k, _ in SRC:
            vv = np.stack([f["actin/tb_%s_%s" % (k, c)][t, :, 0] for c in "xyz"], 1); p = np.sum(vv * tg, 1) / KT; acc[k].append(p[np.abs(p) < 5000])
        lat.append(np.linalg.norm(np.stack([f["actin/tb_aa_lateral_%s" % c][t, :, 0] for c in "xyz"], 1), axis=1) / KT)
    lat = np.concatenate(lat)
    out["director"].append(dict(key=tag, angle_median=float(np.median(np.concatenate(th))),
        sources=[dict(src=k, label=l, mean=float(np.concatenate(acc[k]).mean()), sem=float(np.concatenate(acc[k]).std() / np.sqrt(len(np.concatenate(acc[k]))))) for k, l in SRC],
        lateral_any=float(np.mean(lat > 1e-9)), lateral_median=float(np.median(lat[lat > 1e-9])), lateral_ge4000=float(np.mean(lat >= 4000))))
    st = f["state/current_step"][:n, 0]; t0 = n - 101; steps = float(st[n - 1] - st[t0])
    rc = f["actin/cumulative_rotation_cap_events"][n - 1, :, 0] - f["actin/cumulative_rotation_cap_events"][t0, :, 0]
    tc = f["actin/cumulative_torque_cap_events"][n - 1, :, 0] - f["actin/cumulative_torque_cap_events"][t0, :, 0]
    out["caps"].append(dict(key=tag, torque_cap=float(tc.sum() / steps / tc.size), rotation_cap=float(rc.sum() / steps / rc.size)))
for tag in (() if ONLY_KIN else ("tq30_seed20260805", "k160_30_seed20260805")):
    f = h5py.File(C + tag + ".h5", "r", locking=False); n = int(f["meta/committed_frames"][0]); st = f["state/current_step"][:n, 0]
    rc = f["actin/cumulative_rotation_cap_events"][n - 1, 0, 0] - f["actin/cumulative_rotation_cap_events"][0, 0, 0]
    tc = f["actin/cumulative_torque_cap_events"][n - 1, 0, 0] - f["actin/cumulative_torque_cap_events"][0, 0, 0]
    out["caps"].append(dict(key=tag, torque_cap=float(tc / (st[n - 1] - st[0])), rotation_cap=float(rc / (st[n - 1] - st[0]))))
for arm, lab in (("tq30", "sticky heads"), ("k160_30", "production kinetics"), ("k160trap_30", "production kinetics, actin trapped"),
                 ("k160h81trap_30", "production kinetics, 81 heads/half, trapped")):
    for s in ("20260805", "20260815", "20260825"):
        f = h5py.File(C + "%s_seed%s.h5" % (arm, s), "r", locking=False); n = int(f["meta/committed_frames"][0])
        ua = T.unit(f["actin/direction"][:n, 0]); um = T.unit(f["myosin/direction"][:n, 0])
        th = np.degrees(np.arccos(np.clip(np.abs(np.sum(ua * um, 1)), 0, 1)))
        hs = f["state/am_head_state_current"][:n, 0]; nb = ((hs == 1) | (hs == 2)).sum(1)
        dtf = float(f["meta"].attrs["dt"]) * float(np.median(np.diff(f["state/current_step"][:n, 0])))
        idx = np.concatenate([np.arange(0, min(n, int(0.3 / dtf) + 1)), np.arange(int(0.3 / dtf) + 4, n, 4)])  # full resolution to 300 ms
        dsum = sum(f["myosin/am_align_drive_%s" % h][:n, :, 0] for h in ("minus", "plus")); csum = sum(f["myosin/am_align_head_count_%s" % h][:n, :, 0] for h in ("minus", "plus"))
        out["kin"].append(dict(arm=arm, label=lab, seed=s[-3:], myosin_drive_per_head=float(dsum.sum() / max(csum.sum(), 1) / KT),
                               t_ms=[round(float(i * dtf * 1e3), 2) for i in idx],
                               theta=[round(float(th[i]), 2) for i in idx], heads=[int(nb[i]) for i in idx]))
json.dump(out, open("%s/budget.json" % D, "w"))
for g in out["span"]:
    print(g["label"], [(b["bin"], b["n"], round(b.get("pair_radial", (0,))[0], 1), round(b.get("lateral", (0,))[0], 1)) for b in g["bins"]])
for d in out["director"]: print(d["key"], d["angle_median"], [(s["src"], round(s["mean"])) for s in d["sources"]])
print(out["caps"])
