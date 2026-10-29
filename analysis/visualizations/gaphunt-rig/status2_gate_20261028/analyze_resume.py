"""Matched readout using the established projection, Z-order, alignment and life analyses."""
import contextlib
import hashlib
import json
import os
from pathlib import Path
import re
import runpy
import sys

import h5py
import numpy as np
from connector_metrics import analyze as geometry_overshoot

HERE = Path(__file__).resolve().parent
RIG = HERE.parent
CAMP = Path(os.environ["SCRATCH"]) / "sarcomere/gaphunt_20261023"
OUT = CAMP / "torque_out/status2_gate_20261028"
ARMS = ["r27zt0", "r27zt0np", "r27zt0npg130"]
NEW = ARMS[-1]


def clean(value):
    if isinstance(value, dict): return {str(k): clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)): return [clean(v) for v in value]
    if isinstance(value, np.ndarray): return clean(value.tolist())
    if isinstance(value, np.integer): return int(value)
    if isinstance(value, (float, np.floating)):
        return float(value) if np.isfinite(value) else None
    return value


def invoke(name, args, log_name):
    previous = sys.argv
    sys.argv = [str(RIG / name)] + list(map(str, args))
    print("ANALYSIS", name, *args, flush=True)
    try:
        with open(OUT / log_name, "w") as log, contextlib.redirect_stdout(log):
            return runpy.run_path(str(RIG / name), run_name="__main__")
    finally:
        sys.argv = previous


def life_group(ns, rows):
    if len(rows) == 0:
        return dict(n=0)
    mean, median = ns["km"](rows[:, 0], rows[:, 1])
    death = rows[:, 1] == 0
    return dict(n=len(rows), km_restricted_mean_ms=1000 * mean,
                km_median_ms=1000 * median, censored_fraction=float(np.mean(~death)),
                deaths_by_reason={str(k): int(np.sum(death & (rows[:, 3] == k)))
                                  for k in (1, 2, 3, 4, 5)})


def event_readout(arm):
    path = CAMP / f"data/{arm}_seed20260805.h5"
    with h5py.File(path, "r", locking=False) as f:
        a = f["meta"].attrs
        cols = a["aa_local_link_events_v3_columns"].decode().split(",")
        ec = {s: i for i, s in enumerate(cols)}
        events = f["catch_bond/local_link_events_v3"][:]
        post = events[events[:, ec["step"]] > 200000]
        form = (post[:, ec["old_status"]] == 0) & (post[:, ec["new_status"]] == 1)
        mat = (post[:, ec["old_status"]] == 1) & (post[:, ec["new_status"]] == 2)
        born = set(post[form, ec["link_id"]].astype(int))
        matured = set(post[mat, ec["link_id"]].astype(int))
        off = f["catch_bond/local_link_frame_offsets"]
        identities = f["catch_bond/local_link_identity"]
        inherited_rows = identities[int(off[400, 1]):int(off[400, 1] + off[400, 2])]
        inherited = set(inherited_rows[inherited_rows[:, 4] == 2, 1].astype(int))
        # Cohort membership is determined by link history, not status at death:
        # geometric deletion can coincide with a 2->1 transition.
        death2 = post[(post[:, ec["new_status"]] == 0)]
        deaths_inherited = np.isin(death2[:, ec["link_id"]].astype(int), list(inherited))
        dead_born = np.isin(death2[:, ec["link_id"]].astype(int), list(born & matured))
        end = int(f["meta/committed_frames"][0]) - 1
        ids_end = identities[int(off[end, 1]):int(off[end, 1] + off[end, 2])]
        return dict(
            source=str(path), post_resume_formed=int(form.sum()),
            post_resume_maturation_events=int(mat.sum()),
            maturation_events_per_s=float(mat.sum() / .4),
            new_link_fraction_matured=len(born & matured) / max(1, len(born)),
            inherited_mature_at_start=len(inherited),
            inherited_mature_at_end=int(np.sum(
                np.isin(ids_end[:, 1], list(inherited)) & (ids_end[:, 4] == 2))),
            inherited_mature_deaths_by_reason={
                str(k): int(np.sum(deaths_inherited & (death2[:, ec["detachment_reason"]] == k)))
                for k in (1, 2, 3, 4, 5)},
            post_resume_born_mature_deaths_by_reason={
                str(k): int(np.sum(dead_born & (death2[:, ec["detachment_reason"]] == k)))
                for k in (1, 2, 3, 4, 5)},
        )


def force_budget(arm):
    path = CAMP / f"data/{arm}_seed20260805.h5"
    result = {}
    with h5py.File(path, "r", locking=False) as f:
        for species in ("actin", "myosin"):
            cap = f[species + "/cumulative_force_cap_events"]
            n = f[species + "/center"].shape[1]
            values = cap[400:481, :, 0]
            changes = np.diff(values, axis=0)
            decreases = changes < 0
            reset_frames = np.flatnonzero(np.any(decreases, axis=1))
            # Legacy resumes reset myosin feature counters globally. Some
            # arms were also preempted and resumed again inside this window.
            # A simple endpoint subtraction therefore gives negative rates.
            if species == "myosin":
                changes[reset_frames] = values[1:][reset_frames]
            else:
                changes = np.where(decreases, values[1:], changes)
            fraction = float(np.sum(changes) / (40000 * n))
            assert 0 <= fraction <= 1, (arm, species, fraction)
            result[species + "_force_cap_fraction"] = fraction
            result[species + "_counter_reset_frames"] = (reset_frames + 401).tolist()
    return result


def load_angle_readout(ns):
    """Same 25 ms live-link sampling as load_angle.py; report loads, not
    a frame-weighted reciprocal hazard mislabeled as measured lifetime."""
    f = ns["f"]
    raw_cols = f["meta"].attrs["aa_local_link_mechanics_columns"]
    cols = (raw_cols.decode() if isinstance(raw_cols, bytes) else raw_cols).split(",")
    load_col = cols.index("productive_support")
    records = []
    for frame in range(400, 481, 5):
        offset = f["catch_bond/local_link_frame_offsets"][frame]
        start, count = int(offset[1]), int(offset[2])
        ident = f["catch_bond/local_link_identity"][start:start + count]
        mech = f["catch_bond/local_link_mechanics"][start:start + count]
        mask = ident[:, 4] == 2
        i, j = ident[mask, 2], ident[mask, 3]
        u = f["actin/direction"][frame]
        u /= np.linalg.norm(u, axis=1, keepdims=True)
        aa = np.degrees(np.arccos(np.clip(np.abs(np.sum(u[i] * u[j], axis=1)), 0, 1)))
        th = ns["am"](frame)
        am = np.maximum(th[i], th[j])
        records.append(np.column_stack((aa, am, mech[mask, load_col])))
    rows = np.concatenate(records)
    def summarize(q):
        v = rows[q, 2]
        return dict(n_link_frames=len(v), mean_load_pN=float(np.mean(v)) if len(v) else None,
                    zero_load_fraction=float(np.mean(v <= 0)) if len(v) else None)
    result = dict(all=summarize(np.ones(len(rows), dtype=bool)),
                  am_angle_definition="larger A-M angle; both linked actins require strong heads",
                  frame_stride_ms=25, by_am_angle={}, by_aa_angle={})
    for lo, hi in ns["BINS"]:
        key = f"{lo}-{min(90, hi):g}"
        result["by_am_angle"][key] = summarize((rows[:, 1] >= lo) & (rows[:, 1] < hi))
        result["by_aa_angle"][key] = summarize((rows[:, 0] >= lo) & (rows[:, 0] < hi))
    result["by_am_angle"]["no_strong_head"] = summarize(np.isnan(rows[:, 1]))
    return result


def main():
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Run analysis in Slurm")
    OUT.mkdir(parents=True, exist_ok=True)
    sys.path.insert(0, str(RIG))
    for arm in ARMS:
        with h5py.File(CAMP / f"data/{arm}_seed20260805.h5", "r", locking=False) as f:
            n = int(f["meta/committed_frames"][0])
            assert n == 481 and int(f["state/current_step"][480, 0]) == 240000
    run_manifest = json.loads((OUT / "run_manifest.json").read_text())
    assert run_manifest.get("complete")
    invoke("proj_test_an.py", [",".join(ARMS), OUT / "projection.json"], "projection.txt")
    invoke("zcheck.py", [",".join(ARMS), OUT / "zcheck.json"], "zcheck.txt")
    invoke("resab_an.py", [OUT / "alignment.json", ",".join(ARMS)], "alignment.txt")
    projection = json.loads((OUT / "projection.json").read_text())
    zdisc = json.loads((OUT / "zcheck.json").read_text())
    alignment = json.loads((OUT / "alignment.json").read_text())
    summary = dict(
        window_s=[2.0, 2.4], reference_arm="r27zt0np", new_arm=NEW,
        lifetime_restriction_s=.3, run_manifest=run_manifest, arms={},
        caveats=[
            "One seed and a 0.4 s continuation; inherited status-2 links are not retroactively gated.",
            "A-M angles at maturation use the next saved frame (up to 5 ms later), as in life_decomp.py.",
            "Different CPU models and interrupted control execution prevent a total-wall-time speedup claim.",
            "Force-cap counts are post-summation cap frequency, not per-contact delivered force; legacy myosin counter resets are corrected at saved-frame boundaries.",
        ])
    for arm in ARMS:
        ns = invoke("life_decomp.py", [arm, "200000", "0.3"], arm + ".life.txt")
        groups = {"all": life_group(ns, ns["p2"])}
        for lo, hi in ns["BINS"]:
            p2 = ns["p2"]
            q = (p2[:, 2] >= lo) & (p2[:, 2] < hi)
            groups[f"{lo}-{min(90, hi):g}"] = life_group(ns, p2[q])
        groups["no_strong_head"] = life_group(ns, ns["p2"][np.isnan(ns["p2"][:, 2])])
        loads = load_angle_readout(ns)
        ns["f"].close()
        final = alignment[arm]["rows"][-1]
        zfinal = zdisc[arm][-1]
        summary["arms"][arm] = dict(
            projection=projection[arm], mature_lifetime=groups,
            geometry_overshoot=geometry_overshoot(CAMP / f"data/{arm}_seed20260805.h5",200000,1),
            mature_load_vs_angle=loads,
            event_readout=event_readout(arm), force_budget=force_budget(arm),
            alignment_final=final, zdisc_final=zfinal)
    # Preserve timing inputs and final cumulative phase profiles.
    profiles = {}
    for arm, pattern in (
        ("r27zt0", "resab5-*_0.out"), ("r27zt0np", "resab5-*_1.out"),
        (NEW, NEW + "_seed20260805.simulation.out"),
    ):
        paths = sorted((CAMP / "logs").glob(pattern))
        logs = []
        for path in paths:
            lines = path.read_text().splitlines()
            logs.append(dict(path=str(path), final_runtime_profile=next(
                (line for line in reversed(lines) if line.startswith("runtime_profile cumulative:")), None)))
        profiles[arm] = logs
    summary["timing_profiles"] = profiles
    summary["analysis_sources_sha256"] = {
        str(path): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in [Path(__file__), RIG / "proj_test_an.py", RIG / "zcheck.py",
                     RIG / "resab_an.py", RIG / "hph_an.py", RIG / "life_decomp.py"]}
    (OUT / "summary.json").write_text(json.dumps(clean(summary), indent=2, allow_nan=False) + "\n")
    report = ["# Status-2 connector maturation gate: matched production resume",
              "", "All arms: seed 20260805, 27 heads/half, 2.0–2.4 s. "
              "New arm uses projection off and 130 nm full-connector maturation exclusion.",
              "", "| Quantity | Projection on | Projection off, old gate | Projection off, new gate |",
              "|---|---:|---:|---:|"]
    specs = [
        ("Mature links at 2.4 s", lambda r: r["projection"]["links_end"][1], ".0f"),
        ("New links matured (%)", lambda r: 100*r["event_readout"]["new_link_fraction_matured"], ".1f"),
        ("Mature overstretch deaths/s", lambda r: r["projection"]["reasons"]["status 2"]["per_s_by"]["geometric_contour"], ".0f"),
        ("New-link overstretch deaths/s", lambda r: r["projection"]["reasons"]["status 1"]["per_s_by"]["geometric_contour"], ".0f"),
        ("Mature life, KM restricted mean (ms)", lambda r: r["mature_lifetime"]["all"]["km_restricted_mean_ms"], ".1f"),
        ("Mature-link mean load (pN)", lambda r: r["mature_load_vs_angle"]["all"]["mean_load_pN"], ".2f"),
        ("Mature link-frames at zero load (%)", lambda r: 100*r["mature_load_vs_angle"]["all"]["zero_load_fraction"], ".1f"),
        ("Selected-tip overshoots/frame", lambda r: r["projection"]["penetrating_per_frame"], ".1f"),
        ("Full-connector overlapping link-tip pairs/frame", lambda r: r["geometry_overshoot"]["summary"]["all"]["mean_penetrating_link_tip_pairs_per_frame"], ".1f"),
        ("Full-connector overshoot median (nm)", lambda r: (r["geometry_overshoot"]["summary"]["all"]["depth_nm_p50_p90_max"] or [0,0,0])[0], ".1f"),
        ("Full-connector overshoot p90 (nm)", lambda r: (r["geometry_overshoot"]["summary"]["all"]["depth_nm_p50_p90_max"] or [0,0,0])[1], ".1f"),
        ("Full-connector overshoot maximum (nm)", lambda r: (r["geometry_overshoot"]["summary"]["all"]["depth_nm_p50_p90_max"] or [0,0,0])[2], ".1f"),
        ("Z capture", lambda r: r["zdisc_final"]["zcap"], ".3f"),
        ("Z-band order", lambda r: r["zdisc_final"]["z_order"], ".3f"),
        ("Barbed NN spacing (nm)", lambda r: r["zdisc_final"]["nn_median_nm"], ".1f"),
        ("Pair angle median (deg)", lambda r: r["alignment_final"]["pair_angle_median"], ".1f"),
        ("Actin order S", lambda r: r["alignment_final"]["S_actin"], ".3f"),
        ("Myosin order S", lambda r: r["alignment_final"]["S_myosin"], ".3f"),
    ]
    for label, getter, fmt in specs:
        values = [format(getter(summary["arms"][arm]), fmt) for arm in ARMS]
        report.append("| " + label + " | " + " | ".join(values) + " |")
    report += ["", "Mature lifetimes above use links born after the resume, "
               "with Kaplan–Meier means restricted to 300 ms.",
               "", "See summary.json for angle bins, death reasons, inherited-link cohorts, "
               "force-cap frequencies, timing profiles and input provenance."]
    (OUT / "summary.md").write_text("\n".join(report) + "\n")
    print("\n".join(report), flush=True)
    # Small derived results have a durable copy alongside the campaign scripts.
    durable = HERE / "results"
    durable.mkdir(exist_ok=True)
    for name in ("summary.json", "summary.md", "run_manifest.json"):
        (durable / name).write_text((OUT / name).read_text())
    print("ANALYSIS_COMPLETE", OUT, flush=True)


if __name__ == "__main__":
    main()
