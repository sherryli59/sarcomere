#!/usr/bin/env python3
"""Embed matched ACTN2 failure/correction trajectories."""

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
RUNS = {
    "failed": ROOT / "diagnostic_runs/aa_bug_audit/legacy_k1x_300_full",
    "cohesive": ROOT / "diagnostic_runs/aa_bug_audit/productive_k100x_300_full",
}
TEMPLATE = ROOT / "actn2-cohesion.template.html"
DESTINATION = Path(
    "/home/users/xsli/.codex/visualizations/2026/07/29/"
    "019faf8b-ec1d-74d1-8391-d59895a6c70a/"
    "actual-actn2-cohesion.html"
)


def number(value):
    return round(float(value), 4)


def load_run(prefix):
    with (Path(str(prefix) + ".csv")).open(newline="") as handle:
        geometry_rows = list(csv.DictReader(handle))
    selected_steps = {
        int(row["step"])
        for index, row in enumerate(geometry_rows)
        if index % 4 == 0 or index == len(geometry_rows) - 1
    }

    heads = {step: [] for step in selected_steps}
    with (Path(str(prefix) + "_heads.csv")).open(newline="") as handle:
        for row in csv.DictReader(handle):
            step = int(row["step"])
            if step not in selected_steps:
                continue
            state = int(row["state"])
            record = [
                state,
                number(row["anchor_x"]),
                number(row["anchor_y"]),
                number(row["anchor_z"]),
            ]
            if int(row["bound_actin"]) >= 0 and state in (1, 2):
                record.extend([
                    number(row["binding_x"]),
                    number(row["binding_y"]),
                    number(row["binding_z"]),
                    number(row["extension"]),
                    int(row["bound_actin"]),
                ])
            heads[step].append(record)

    links = {step: [] for step in selected_steps}
    with (Path(str(prefix) + "_links.csv")).open(newline="") as handle:
        for row in csv.DictReader(handle):
            step = int(row["step"])
            if step not in selected_steps:
                continue
            links[step].append([
                int(row["status"]),
                int(row["actin_i"]),
                int(row["actin_j"]),
                number(row["first_x"]),
                number(row["first_y"]),
                number(row["first_z"]),
                number(row["second_x"]),
                number(row["second_y"]),
                number(row["second_z"]),
                number(row["length"]),
            ])

    frames = []
    for row in geometry_rows:
        step = int(row["step"])
        if step not in selected_steps:
            continue
        actins = [[
            number(row[f"a{i}_x"]),
            number(row[f"a{i}_y"]),
            number(row[f"a{i}_z"]),
            number(row[f"a{i}_ux"]),
            number(row[f"a{i}_uy"]),
            number(row[f"a{i}_uz"]),
        ] for i in range(4)]
        myosins = [[
            number(row[f"m{i}_x"]),
            number(row[f"m{i}_y"]),
            number(row[f"m{i}_z"]),
            number(row[f"m{i}_ux"]),
            number(row[f"m{i}_uy"]),
            number(row[f"m{i}_uz"]),
        ] for i in range(2)]
        frames.append({
            "t": number(row["time"]),
            "L": number(row["Lx"]),
            "aa01": [
                number(row["aa01_axial"]),
                number(row["aa01_transverse"]),
                number(row["aa01_distance"]),
                int(row["aa01_status1"]),
                int(row["aa01_status2"]),
            ],
            "aa23": [
                number(row["aa23_axial"]),
                number(row["aa23_transverse"]),
                number(row["aa23_distance"]),
                int(row["aa23_status1"]),
                int(row["aa23_status2"]),
            ],
            "a": actins,
            "m": myosins,
            "h": heads[step],
            "links": links[step],
        })
    return frames


payload = {name: load_run(prefix) for name, prefix in RUNS.items()}
fragment = TEMPLATE.read_text().replace(
    "__RUN_DATA__",
    json.dumps(payload, separators=(",", ":"), allow_nan=False),
)
DESTINATION.parent.mkdir(parents=True, exist_ok=True)
DESTINATION.write_text(fragment)
print(DESTINATION)
print(
    f"{DESTINATION.stat().st_size} bytes, "
    f"{len(payload['failed'])} synchronized frames"
)
