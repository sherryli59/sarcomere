"""Hard-projection test (resab5, job 49878921): r27zt0 (titin off, Z-stop gated on ACTN2,
projection on) vs r27zt0np (same, projection off), both resumed from the hph27 production state
at 2.0 s (step 200000) to 2.4 s.  Readouts per arm, frames 401-480 (5 ms):
 1. projection events per step (actin counter / 2: each projection increments two actins)
 2. ACTN2-patch penetration of myosin tips: the model's recorded per-tip axial gap of the
    selected (largest-force) patch, myosin/aa_zdisc_tip_gap_{minus,plus}, where the tip force > 0.
    Negative = the selected contact point is past the myosin tip. Its geometry
    depends on the executable version (midpoint or full connector), so use
    connector_metrics.py for consistent cross-version geometric overshoot.
 3. barbed-tip penetration into myosin capsules of actins carrying ACTN2 links
 4. ACTN2 detachment reasons by status, 2.0 -> 2.4 s
 5. live / mature link counts
 6. wall clock from the job logs"""
import h5py, numpy as np, os, json, re, glob, sys
from detachment_reasons import reason_names
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
D = "/scratch/users/xsli/sarcomere/gaphunt_20261023/"
ARMS = sys.argv[1].split(",") if len(sys.argv) > 1 else ["r27zt0", "r27zt0np"]
OUTPUT = sys.argv[2] if len(sys.argv) > 2 else D + "torque_out/proj_test.json"
out = {}
for arm in ARMS:
    P = D + "data/%s_seed20260805.h5" % arm
    f = h5py.File(P, "r", locking=False); a = dict(f["meta"].attrs); dt = float(a["dt"])
    nF = int(f["meta/committed_frames"][0]); fs = f["state/current_step"][:nF, 0]
    T0 = int(np.searchsorted(fs, 200000)); T1 = nF - 1
    R = dict(note=a.get("resume_ab_note", b"").decode() if isinstance(a.get("resume_ab_note", b""), bytes) else str(a.get("resume_ab_note")),
             frames=[T0, T1], t_end_s=float(fs[T1] * dt))
    pe = f["actin/cumulative_aa_zdisc_tip_projection_events"]
    R["projections_per_step"] = float((pe[T1][:, 0].sum() - pe[T0][:, 0].sum()) / 2 / max(1, fs[T1] - fs[T0]))
    # 2. patch penetration
    G, T = [], []
    per_frame = []
    for t in range(T0 + 1, T1 + 1):
        gs = []
        for s in ("minus", "plus"):
            g = f["myosin/aa_zdisc_tip_gap_" + s][t][:, 0]; fo = f["myosin/aa_zdisc_tip_force_" + s][t][:, 0]
            gs.append(g[fo > 0])
        gs = np.concatenate(gs); G.append(gs)
        per_frame.append((float(fs[t] * dt), int(len(gs)), int((gs < -1e-3).sum()), float(-1000 * gs.min()) if len(gs) else 0.0,
                          float(-1000 * np.median(gs[gs < -1e-3])) if (gs < -1e-3).any() else 0.0))
    G = np.concatenate(G)
    pen = -1000 * G[G < -1e-3]
    R["tips_in_contact_per_frame"] = float(len(G) / (T1 - T0))
    R["penetrating_per_frame"] = float(len(pen) / (T1 - T0))
    R["penetration_depth_nm_p50_p90_max"] = [float(np.median(pen)), float(np.percentile(pen, 90)), float(pen.max())] if len(pen) else [0, 0, 0]
    R["per_frame"] = per_frame
    # 3. barbed tip penetration of linked actins
    bt = []; kp = []
    for t in range(T0 + 1, T1 + 1, 4):
        lc = f["actin/aa_local_link_count"][t][:, 0] > 0
        bt.append(f["actin/aa_barbed_tip_penetration"][t][lc, 0]); kp.append(f["actin/aa_barbed_tip_kill_penetration"][t][lc, 0])
    bt = np.concatenate(bt); kp = np.concatenate(kp)
    R["linked_tip_penetration"] = dict(frac_pos=float(np.mean(bt > 0)), p90_nm=float(1000 * np.percentile(bt, 90)), max_nm=float(1000 * bt.max()),
                                       kill_frac_pos=float(np.mean(kp > 0)), kill_max_nm=float(1000 * kp.max()))
    # 4. detachment reasons
    EC = {c: i for i, c in enumerate(a["aa_local_link_events_v3_columns"].decode().split(","))}
    RN = reason_names(a)
    ev = f["catch_bond/local_link_events_v3"][:]
    span = (fs[T1] - 200000) * dt
    d = ev[(ev[:, EC["new_status"]] == 0) & (ev[:, EC["step"]] > 200000)]
    formed = int(((ev[:, EC["old_status"]] == 0) & (ev[:, EC["new_status"]] == 1) & (ev[:, EC["step"]] > 200000)).sum())
    R["formed_per_s"] = formed / span; R["detached_per_s"] = len(d) / span
    R["reasons"] = {}
    for sv in (1, 2):
        q = d[d[:, EC["old_status"]] == sv]; r = q[:, EC["detachment_reason"]].astype(int)
        R["reasons"]["status %d" % sv] = dict(n=int(len(q)), per_s=float(len(q) / span), share={RN[str(k)]: float(np.mean(r == k)) if len(r) else 0.0 for k in (1, 2, 3, 4, 5)},
                                              per_s_by={RN[str(k)]: float(np.sum(r == k) / span) for k in (1, 2, 3, 4, 5)})
    R["selected_tip_diagnostic_note"] = "Recorded gap/rho use the simulated winner and its runtime geometry; use geometry_overshoot for cross-version comparisons."
    # 5. links
    off = f["catch_bond/local_link_frame_offsets"]; idn = f["catch_bond/local_link_identity"]
    def counts(t):
        o, c = int(off[t][1]), int(off[t][2]); s = idn[o:o + c][:, 4]; return int(c), int((s == 2).sum())
    R["links_start"], R["links_end"] = counts(T0), counts(T1)
    # 6. wall clock
    log_pattern = {
        "r27zt0": "resab5-*_0.out",
        "r27zt0np": "resab5-*_1.out",
        "r27zt0c2": "resab6-*.out",
    }.get(arm, arm + "_seed20260805.simulation.out")
    logs = sorted(glob.glob(D + "logs/" + log_pattern))
    if logs:
        txt = open(logs[-1]).read(); walls = re.findall(r"^exit=.*?\bwall_s=(\d+(?:\.\d+)?)", txt, re.M); h = re.search(r"ARM \S+ from \S+ on (\S+) (.*)", txt)
        R["wall_s"] = float(walls[-1]) if walls else None; R["host"] = h.group(1) if h else None; R["cpu"] = h.group(2).strip() if h else None
        R["profile_line"] = next((l for l in txt.splitlines() if l.startswith("actn2_tip_projection profile")), None)
    out[arm] = R
    print("\n== %s  (to %.3f s)  %s" % (arm, R["t_end_s"], R["note"][:160]))
    print("projections per step %.2f | tips in patch contact per frame %.1f, penetrating %.1f; depth p50/p90/max %.1f / %.1f / %.1f nm"
          % (R["projections_per_step"], R["tips_in_contact_per_frame"], R["penetrating_per_frame"], *R["penetration_depth_nm_p50_p90_max"]))
    print("linked actins: barbed tip inside a myosin %.1f%% (p90 %.1f nm, max %.1f nm); kill-penetration > 0 %.2f%% (max %.1f nm)"
          % (100 * R["linked_tip_penetration"]["frac_pos"], R["linked_tip_penetration"]["p90_nm"], R["linked_tip_penetration"]["max_nm"],
             100 * R["linked_tip_penetration"]["kill_frac_pos"], R["linked_tip_penetration"]["kill_max_nm"]))
    print("links formed %.0f/s, detached %.0f/s | live links %d -> %d, mature %d -> %d" % (R["formed_per_s"], R["detached_per_s"], R["links_start"][0], R["links_end"][0], R["links_start"][1], R["links_end"][1]))
    for k, v in R["reasons"].items():
        print("   %s detachments %.0f/s: %s" % (k, v["per_s"], ", ".join("%s %.1f%% (%.0f/s)" % (n, 100 * v["share"][n], v["per_s_by"][n]) for n in v["share"])))
    print("wall %s s on %s (%s)" % (R.get("wall_s"), R.get("host"), R.get("cpu")))
    if len(per_frame):
        print("penetration over time (t, tips in contact, penetrating, deepest nm, median depth nm):", [tuple(round(x, 1) if isinstance(x, float) else x for x in p) for p in per_frame[::10]])
json.dump(out, open(OUTPUT, "w"))
