#!/usr/bin/env python3
"""Embed one representative full-simulator trajectory in the inline visual."""

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
GEOMETRY = ROOT / "results/representative/register_full_capped.csv"
HEADS = ROOT / "results/representative/register_full_capped_heads.csv"
TEMPLATE = ROOT / "balanced-overlap-diagnostic.template.html"
DESTINATION = Path(
    "/home/users/xsli/.codex/visualizations/2026/07/29/"
    "019faf8b-ec1d-74d1-8391-d59895a6c70a/"
    "balanced-overlap-diagnostic.html"
)


def number(value):
    return round(float(value), 5)


with GEOMETRY.open(newline="") as handle:
    geometry_rows = list(csv.DictReader(handle))

selected_steps = {
    int(row["step"])
    for index, row in enumerate(geometry_rows)
    if index % 2 == 0 or index == len(geometry_rows) - 1
}

heads_by_step = {step: [] for step in selected_steps}
with HEADS.open(newline="") as handle:
    for row in csv.DictReader(handle):
        step = int(row["step"])
        if step not in selected_steps:
            continue
        state = int(row["state"])
        bound = int(row["bound_actin"])
        anchor_r = number(row["anchor_y"]) + 0.35 * number(row["anchor_z"])
        if bound >= 0 and state in (1, 2):
            binding_x = number(row["binding_x"])
            binding_r = (
                number(row["binding_y"]) +
                0.35 * number(row["binding_z"])
            )
        else:
            binding_x = None
            binding_r = None
        heads_by_step[step].append([
            state, bound, number(row["anchor_x"]), round(anchor_r, 5),
            binding_x, None if binding_r is None else round(binding_r, 5),
            number(row["extension"]), number(row["active_force"]),
        ])

frames = []
for index, row in enumerate(geometry_rows):
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
                5,
            ),
            number(row[f"a{filament}_ux"]),
            round(
                number(row[f"a{filament}_uy"]) +
                0.35 * number(row[f"a{filament}_uz"]),
                5,
            ),
        ])
    myosins = []
    for filament in range(2):
        myosins.append([
            number(row[f"m{filament}_x"]),
            round(
                number(row[f"m{filament}_y"]) +
                0.35 * number(row[f"m{filament}_z"]),
                5,
            ),
            number(row[f"m{filament}_ux"]),
            round(
                number(row[f"m{filament}_uy"]) +
                0.35 * number(row[f"m{filament}_uz"]),
                5,
            ),
        ])
    frames.append({
        "t": number(row["time"]),
        "L": number(row["Lx"]),
        "aa": number(row["aa_overlap"]),
        "am": number(row["am_overlap"]),
        "weak": number(row["weak_heads"]),
        "strong": number(row["strong_heads"]),
        "wext": number(row["weak_extension"]),
        "sext": number(row["strong_extension"]),
        "s1": number(row["aa_status1"]),
        "s2": number(row["aa_status2"]),
        "a": actins,
        "m": myosins,
        "h": heads_by_step[step],
    })

payload = json.dumps(frames, separators=(",", ":"), allow_nan=False)
fragment = TEMPLATE.read_text().replace("__FRAME_DATA__", payload)
DESTINATION.parent.mkdir(parents=True, exist_ok=True)
DESTINATION.write_text(fragment)
print(DESTINATION)
print(f"{DESTINATION.stat().st_size} bytes, {len(frames)} frames")
