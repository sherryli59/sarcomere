#!/usr/bin/env python3
"""Embed matched pre-equilibrated trajectories for an inline comparison."""

import csv
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent
RUNS = {
    "h35": (
        ROOT / "diagnostic_runs/force_attribution/"
        "corrected_register_titin_h35_pre50_release300_full.csv",
        ROOT / "diagnostic_runs/force_attribution/"
        "corrected_register_titin_h35_pre50_release300_full_heads.csv",
    ),
    "h141": (
        ROOT / "diagnostic_runs/force_attribution/"
        "corrected_register_titin_h141_pre50_release300_full.csv",
        ROOT / "diagnostic_runs/force_attribution/"
        "corrected_register_titin_h141_pre50_release300_full_heads.csv",
    ),
}
TEMPLATE = ROOT / "preequil-head-density.template.html"
DESTINATION = Path(
    "/home/users/xsli/.codex/visualizations/2026/07/29/"
    "019faf8b-ec1d-74d1-8391-d59895a6c70a/"
    "corrected-balanced-heads.html"
)

META = {
    "keys": ["h35", "h141"],
    "labels": ["35 heads / myosin half", "141 heads / myosin half"],
    "trace_labels": ["35 heads/half", "141"],
    "aria": (
        "Synchronized animation comparing a force-balanced sarcomere lattice "
        "with 35 and 141 myosin heads per half after fixed-lattice chemical "
        "pre-equilibration."
    ),
}
AUDIT_MODE = False

if len(sys.argv) == 2 and sys.argv[1] == "--aa-fix":
    AUDIT_MODE = True
    RUNS = {
        "buggy": (
            ROOT / "diagnostic_runs/aa_bug_audit/baseline_v2.csv",
            ROOT / "diagnostic_runs/aa_bug_audit/baseline_v2_heads.csv",
        ),
        "fixed": (
            ROOT / "diagnostic_runs/aa_bug_audit/consistent_v2.csv",
            ROOT / "diagnostic_runs/aa_bug_audit/consistent_v2_heads.csv",
        ),
    }
    META = {
        "keys": ["buggy", "fixed"],
        "labels": [
            "Buggy: endpoint force only",
            "Fixed: force + endpoint torque",
        ],
        "trace_labels": ["buggy", "fixed"],
        "aria": (
            "Synchronized deterministic animation comparing the incomplete "
            "A-A endpoint register with the conservative force-and-torque "
            "correction under the same seed."
        ),
    }
    DESTINATION = Path(
        "/home/users/xsli/.codex/visualizations/2026/07/29/"
        "019faf8b-ec1d-74d1-8391-d59895a6c70a/"
        "conservative-aa-register-fix.html"
    )
elif len(sys.argv) != 1:
    raise SystemExit("usage: make_preequil_compare_visual.py [--aa-fix]")


def number(value):
    return round(float(value), 4)


def load_run(geometry_path, heads_path):
    with geometry_path.open(newline="") as handle:
        geometry_rows = list(csv.DictReader(handle))

    selected_steps = {
        int(row["step"])
        for index, row in enumerate(geometry_rows)
        if index % 4 == 0 or index == len(geometry_rows) - 1
    }
    heads_by_step = {step: [] for step in selected_steps}
    pairs_by_step = {
        step: [[0, 0], [0, 0], [0, 0], [0, 0]]
        for step in selected_steps
    }
    pair_slot = {(0, 0): 0, (0, 2): 1, (1, 1): 2, (1, 3): 3}
    with heads_path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            step = int(row["step"])
            if step not in selected_steps:
                continue
            state = int(row["state"])
            anchor_r = number(row["anchor_y"]) + 0.35 * number(row["anchor_z"])
            head = [
                state,
                number(row["anchor_x"]),
                round(anchor_r, 4),
            ]
            if int(row["bound_actin"]) >= 0 and state in (1, 2):
                pair = (int(row["myosin"]), int(row["bound_actin"]))
                if pair in pair_slot:
                    pairs_by_step[step][pair_slot[pair]][state - 1] += 1
                binding_r = (
                    number(row["binding_y"]) +
                    0.35 * number(row["binding_z"])
                )
                head.extend([
                    number(row["binding_x"]),
                    round(binding_r, 4),
                    number(row["extension"]),
                ])
            heads_by_step[step].append(head)

    frames = []
    for row in geometry_rows:
        step = int(row["step"])
        if step not in selected_steps:
            continue
        actins = []
        for filament in range(4):
            actins.append([
                number(row[f"a{filament}_x"]),
                round(
                    number(row[f"a{filament}_y"]) +
                    0.35 * number(row[f"a{filament}_z"]),
                    4,
                ),
                number(row[f"a{filament}_ux"]),
                round(
                    number(row[f"a{filament}_uy"]) +
                    0.35 * number(row[f"a{filament}_uz"]),
                    4,
                ),
            ])
        myosins = []
        for filament in range(2):
            myosins.append([
                number(row[f"m{filament}_x"]),
                round(
                    number(row[f"m{filament}_y"]) +
                    0.35 * number(row[f"m{filament}_z"]),
                    4,
                ),
                number(row[f"m{filament}_ux"]),
                round(
                    number(row[f"m{filament}_uy"]) +
                    0.35 * number(row[f"m{filament}_uz"]),
                    4,
                ),
            ])
        frames.append({
            "t": number(row["time"]),
            "L": number(row["Lx"]),
            "aa": number(row["aa_overlap"]),
            "am": number(row["am_overlap"]),
            "q": number(row["thin_q"]),
            "w": int(float(row["weak_heads"])),
            "s": int(float(row["strong_heads"])),
            "a": actins,
            "m": myosins,
            "h": heads_by_step[step],
            "b": pairs_by_step[step],
        })
    return frames


payload = {
    name: load_run(*paths)
    for name, paths in RUNS.items()
}
payload["_meta"] = META
fragment = TEMPLATE.read_text().replace(
    "__RUN_DATA__",
    json.dumps(payload, separators=(",", ":"), allow_nan=False),
)
if AUDIT_MODE:
    fragment = fragment.replace(
        "preequil-head-density", "aa-register-fix"
    ).replace("phd-", "aar-")
DESTINATION.parent.mkdir(parents=True, exist_ok=True)
DESTINATION.write_text(fragment)
print(DESTINATION)
print(
    f"{DESTINATION.stat().st_size} bytes, "
    f"{len(payload[META['keys'][0]])} synchronized frames"
)
