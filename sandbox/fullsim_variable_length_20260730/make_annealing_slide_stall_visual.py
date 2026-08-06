#!/usr/bin/env python3
"""Embed the ACTN2-only and centered 1 s annealing challenges."""

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
RUNS = {
    "unretained": ROOT / (
        "diagnostic_runs/anneal_stall_fixed/"
        "aa_register_all_l3_AA016_1s_full"
    ),
    "centered": ROOT / (
        "diagnostic_runs/anneal_stall_fixed/"
        "titin_register_all_l3_AA016_1s_full_v2"
    ),
}
TEMPLATE = ROOT / "actn2-cohesion.template.html"
DESTINATION = Path(
    "/home/users/xsli/.codex/visualizations/2026/07/29/"
    "019faf8b-ec1d-74d1-8391-d59895a6c70a/"
    "annealing-slide-stall.html"
)


def number(value):
    return round(float(value), 4)


def wrap(value, length):
    return value - length * round(value / length)


def dot(first, second):
    return sum(a * b for a, b in zip(first, second))


def am_overlap(actin, myosin, length):
    displacement = [
        wrap(actin[0] - myosin[0], length),
        actin[1] - myosin[1],
        actin[2] - myosin[2],
    ]
    axis = myosin[3:6]
    center = dot(displacement, axis)
    half = 0.5 * dot(actin[3:6], axis)
    lower, upper = sorted((center - half, center + half))
    return max(0.0, min(upper, 0.75) - max(lower, -0.75))


def load_run(prefix):
    with Path(str(prefix) + ".csv").open(newline="") as handle:
        geometry_rows = list(csv.DictReader(handle))
    selected_steps = {
        int(row["step"])
        for index, row in enumerate(geometry_rows)
        if index % 10 == 0 or index == len(geometry_rows) - 1
    }

    heads = {step: [] for step in selected_steps}
    with Path(str(prefix) + "_heads.csv").open(newline="") as handle:
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
    with Path(str(prefix) + "_links.csv").open(newline="") as handle:
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
        length = number(row["Lx"])
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
        owning_myosin = (0, 1, 0, 1)
        overlaps = [
            am_overlap(actins[i], myosins[owning_myosin[i]], length)
            for i in range(4)
        ]
        frames.append({
            "t": number(row["time"]),
            "L": length,
            "am": round(sum(overlaps) / 4.0, 4),
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
fragment = TEMPLATE.read_text()
replacements = {
    "actual-actn2-cohesion": "annealing-slide-stall",
    "aac-": "aas-",
    'const keys = ["failed", "cohesive"];':
        'const keys = ["unretained", "centered"];',
    '    "Original: legacy state + 0.01 pN/nm total",\n'
    '    "Corrected: productive state + 1 pN/nm total",':
        '    "ACTN2 only: anneals, then loses bilateral register",\n'
        '    "+ reciprocal A–M centering: slide, then stall",',
    "  const frames = runs.failed;": "  const frames = runs.unretained;",
    "mean 3D A–A endpoint separation": "mean A–M axial overlap",
    "Math.max(0.65, ...keys.flatMap((key) =>\n"
    "      runs[key].flatMap((frame) => [frame.aa01[2], frame.aa23[2]])\n"
    "    ))":
        "Math.max(0.75, ...keys.flatMap((key) =>\n"
        "      runs[key].map((frame) => frame.am)\n"
        "    ))",
    "const value = 0.5 * (frame.aa01[2] + frame.aa23[2]);":
        "const value = frame.am;",
    'context.fillText("300 ms", right, bottom + 18);':
        'context.fillText("1,000 ms", right, bottom + 18);',
    "      `At ${Math.round(time * 1000)} milliseconds, the original A0-A1 ` +\n"
    "      `pair has ${runs.failed[index].aa01[3] + runs.failed[index].aa01[4]} ` +\n"
    "      `live ACTN2 links and ${runs.failed[index].aa01[2].toFixed(3)} ` +\n"
    "      `micrometers full endpoint separation; the corrected pair has ` +\n"
    "      `${runs.cohesive[index].aa01[3] + runs.cohesive[index].aa01[4]} ` +\n"
    "      `links and ${runs.cohesive[index].aa01[2].toFixed(3)} micrometers.`;":
        "      `At ${Math.round(time * 1000)} milliseconds, mean A-M overlap is ` +\n"
        "      `${runs.unretained[index].am.toFixed(3)} micrometers without ` +\n"
        "      `centering and ${runs.centered[index].am.toFixed(3)} micrometers ` +\n"
        "      `with reciprocal centering.`;",
}
for old, new in replacements.items():
    if old not in fragment:
        raise RuntimeError(f"template text not found: {old[:80]}")
    fragment = fragment.replace(old, new)
fragment = fragment.replace(
    "__RUN_DATA__",
    json.dumps(payload, separators=(",", ":"), allow_nan=False),
)
DESTINATION.parent.mkdir(parents=True, exist_ok=True)
DESTINATION.write_text(fragment)
print(DESTINATION)
print(
    f"{DESTINATION.stat().st_size} bytes, "
    f"{len(payload['unretained'])} synchronized frames"
)
