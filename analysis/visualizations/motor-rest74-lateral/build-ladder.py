#!/usr/bin/env python3
from __future__ import annotations

import base64
import importlib.util
import json
import os
from pathlib import Path
import subprocess

os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

import h5py
import numpy as np


ROOT = Path("/scratch/users/xsli/sarcomere")
VIS_ROOT = Path(
    "/home/groups/rotskoff/sherry/sarcomere/sarcomere/analysis/"
    "visualizations/motor-rest74-lateral"
)
# build-simulations.py supplies committed_frames and clip_periodic_segments and
# is copied alongside this builder so the viewer does not depend on .codex.
TEMPLATE = VIS_ROOT / "motor-rest74-lateral-template.html"
OUTPUT = VIS_ROOT / "ladder-latest-snapshots.html"
PIPELINE = VIS_ROOT / "build-simulations.py"
TARGET_TIME_MS = 5000.0

RUNS = [
    {
        "key": "k4",
        "label": "k=4 \u00b7 seed 20260825 \u00b7 COMPLETE \u00b7 S=0.25 \u00b7 NEVER ORDERS",
        "config": "A1600 M320 Lx=6 \u00b7 azimuth-OFF radial \u00b7 kappa_aa 3 \u00b7 rest 79 nm \u00b7 titin ON (k_off_unsup 10) \u00b7 weak-alignment ladder, LIVE",
        "jobId": "live",
        "path": ROOT / "weakalign_5s_20260922/data/" "A1600_M320_Lx6_intermediate_alignk4_seed20260825_5000ms.h5",
    },
    {
        "key": "k12-ordered",
        "label": "k=12 seed 20260805 \u00b7 933/1001 \u00b7 S=0.933 \u00b7 BIMODAL: ordered branch",
        "config": "A1600 M320 Lx=6 \u00b7 azimuth-OFF radial \u00b7 kappa_aa 3 \u00b7 rest 79 nm \u00b7 titin ON (k_off_unsup 10) \u00b7 weak-alignment ladder, LIVE",
        "jobId": "live",
        "path": ROOT / "weakalign_5s_20260922/data/" "A1600_M320_Lx6_intermediate_alignk12_seed20260805_5000ms.h5",
    },
    {
        "key": "k12-disordered",
        "label": "k=12 seed 20260835 \u00b7 289/1001 \u00b7 S=0.352 \u00b7 BIMODAL: disordered branch",
        "config": "A1600 M320 Lx=6 \u00b7 azimuth-OFF radial \u00b7 kappa_aa 3 \u00b7 rest 79 nm \u00b7 titin ON (k_off_unsup 10) \u00b7 weak-alignment ladder, LIVE",
        "jobId": "live",
        "path": ROOT / "weakalign_5s_20260922/data/" "A1600_M320_Lx6_intermediate_alignk12_seed20260835_5000ms.h5",
    },
    {
        "key": "k36",
        "label": "k=36 \u00b7 seed 20260845 \u00b7 694/1001 \u00b7 S=0.980 \u00b7 ORDERS IN 6/6 SEEDS",
        "config": "A1600 M320 Lx=6 \u00b7 azimuth-OFF radial \u00b7 kappa_aa 3 \u00b7 rest 79 nm \u00b7 titin ON (k_off_unsup 10) \u00b7 weak-alignment ladder, LIVE",
        "jobId": "live",
        "path": ROOT / "weakalign_5s_20260922/data/" "A1600_M320_Lx6_intermediate_alignk36_seed20260845_5000ms.h5",
    },
    {
        "key": "prealigned",
        "label": "PRE-ALIGNED k=0 \u00b7 sigma 0.125 \u00b7 seed 20260815 \u00b7 corrected, links now building",
        "config": "A1600 M320 Lx=6 \u00b7 am_alignment_k = 0 \u00b7 PRE-ALIGNED, polar_fraction 0.5 (CORRECTED: mixed polarity, antiparallel pairs restored 0.0000 -> 0.5007) \u00b7 LIVE",
        "jobId": "live",
        "path": ROOT / "prealigned_k0_5s_20260922/data/" "A1600_M320_Lx6_intermediate_sig0.125_seed20260815_5000ms.h5",
    },
    {
        "key": "ref-k120",
        "label": "REFERENCE k=120 \u00b7 seed 20260805 \u00b7 COMPLETE \u00b7 S=0.978",
        "config": "A1600 M320 Lx=6 \u00b7 am_alignment_k = 120 \u00b7 COMPLETE reference \u00b7 S=0.978 but network retention only 0.640",
        "jobId": "live",
        "path": ROOT / "production_titin_5s_20260918/data/" "A1600_M320_Lx6_intermediate_titinon_seed20260805_5000ms.h5",
    },
]






def scheduler_states() -> dict[str, dict[str, str]]:
    ids = ",".join(run["jobId"] for run in RUNS)
    result = subprocess.run(
        ["squeue", "-h", "-j", ids, "-o", "%i|%T|%M|%R"],
        capture_output=True,
        text=True,
        check=False,
    )
    states: dict[str, dict[str, str]] = {}
    for line in result.stdout.splitlines():
        job_id, state, runtime, reason = line.strip().split("|", 3)
        states[job_id] = {
            "stateLabel": state.title(),
            "runTime": runtime,
            "reason": reason.strip("()"),
        }
    return states


def encode(values: np.ndarray) -> str:
    packed = np.ascontiguousarray(values, dtype=np.float32)
    return base64.b64encode(packed.tobytes()).decode("ascii")


def order(values: np.ndarray, box_x: float) -> np.ndarray:
    wrapped = np.asarray(values, dtype=float) % box_x
    harmonics = np.arange(1, 5, dtype=float)
    phase = np.exp(
        2j * np.pi * wrapped[..., None] * harmonics[None, None, :] / box_x
    ).mean(axis=1)
    return np.max(np.abs(phase), axis=1)


spec = importlib.util.spec_from_file_location("snapshot_pipeline", PIPELINE)
pipeline = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(pipeline)


def snapshot(path: Path) -> dict | None:
    if not path.exists():
        return None
    with h5py.File(path, "r", locking=False) as handle:
        usable = pipeline.committed_frames(handle)
        if usable < 1:
            return None
        frame = usable - 1
        first_recent = max(0, frame - 20)
        frames = slice(first_recent, frame + 1)
        meta = handle["meta"].attrs
        box = np.asarray([meta["box_x"], meta["box_y"], meta["box_z"]], float)
        dt = float(meta["dt"])
        actin_length = float(meta["actin_length"])
        myosin_length = float(meta["myosin_length"])
        steps = np.asarray(handle["state/current_step"][frames, 0], dtype=int)
        centers_a_recent = np.asarray(handle["actin/center"][frames], dtype=float)
        directions_a_recent = np.asarray(
            handle["actin/direction"][frames], dtype=float
        )
        centers_m_recent = np.asarray(handle["myosin/center"][frames], dtype=float)
        m_values = order(centers_m_recent[:, :, 0], box[0])
        plus = centers_a_recent[:, :, 0] + 0.5 * actin_length * directions_a_recent[:, :, 0]
        minus = centers_a_recent[:, :, 0] - 0.5 * actin_length * directions_a_recent[:, :, 0]
        z_values = np.maximum(order(plus, box[0]), order(minus, box[0]))

        centers_a = centers_a_recent[-1]
        directions_a = directions_a_recent[-1]
        centers_m = centers_m_recent[-1]
        directions_m = np.asarray(handle["myosin/direction"][frame], dtype=float)
        statuses_a = np.asarray(handle["actin/cb_status"][frame, :, 0], dtype=int)
        a_half = 0.5 * actin_length * directions_a
        m_half = 0.5 * myosin_length * directions_m
        groups: dict[str, np.ndarray] = {}
        for status, key in ((0, "actin0"), (1, "actin1"), (2, "actin2")):
            keep = statuses_a == status
            groups[key] = np.asarray(
                pipeline.clip_periodic_segments(
                    (centers_a - a_half)[keep], (centers_a + a_half)[keep], box
                ),
                dtype=np.float32,
            )
        groups["myosin"] = np.asarray(
            pipeline.clip_periodic_segments(
                centers_m - m_half, centers_m + m_half, box
            ),
            dtype=np.float32,
        )
        barbed = centers_a - a_half
        groups["barbed"] = ((barbed + 0.5 * box) % box - 0.5 * box).astype(
            np.float32
        )

        link_counts = {"1": 0, "2": 0}
        groups["aa1"] = np.empty((0, 6), dtype=np.float32)
        groups["aa2"] = np.empty((0, 6), dtype=np.float32)
        offset_path = "catch_bond/local_link_frame_offsets"
        identity_path = "catch_bond/local_link_identity"
        coordinate_path = "catch_bond/local_link_coordinates"
        if all(name in handle for name in (offset_path, identity_path, coordinate_path)):
            _, offset, count = np.asarray(handle[offset_path][frame], dtype=int)
            identity = np.asarray(
                handle[identity_path][offset : offset + count], dtype=int
            )
            coordinates = np.asarray(
                handle[coordinate_path][offset : offset + count], dtype=float
            )
            for status, key in ((1, "aa1"), (2, "aa2")):
                keep = identity[:, 4] == status
                selected_identity = identity[keep]
                selected_coordinates = coordinates[keep]
                link_counts[str(status)] = int(np.count_nonzero(keep))
                if len(selected_identity):
                    first_ids = selected_identity[:, 2]
                    second_ids = selected_identity[:, 3]
                    starts = centers_a[first_ids] + selected_coordinates[:, :1] * directions_a[first_ids]
                    raw_ends = centers_a[second_ids] + selected_coordinates[:, 1:] * directions_a[second_ids]
                    ends = starts + (
                        raw_ends
                        - starts
                        - box[None, :] * np.rint((raw_ends - starts) / box[None, :])
                    )
                    groups[key] = np.asarray(
                        pipeline.clip_periodic_segments(starts, ends, box),
                        dtype=np.float32,
                    )

        time_ms = float(steps[-1] * dt * 1e3)
        return {
            "frame": int(frame),
            "step": int(steps[-1]),
            "timeMs": time_ms,
            "progressPct": 100.0 * time_ms / TARGET_TIME_MS,
            "box": np.round(box, 4).tolist(),
            "actinStatusCounts": {
                str(status): int(np.count_nonzero(statuses_a == status))
                for status in (0, 1, 2)
            },
            "linkStatusCounts": link_counts,
            "myosinCount": int(len(centers_m)),
            "mOrder": float(m_values[-1]),
            "zOrder": float(z_values[-1]),
            "mDelta100": float(m_values[-1] - m_values[0]),
            "zDelta100": float(z_values[-1] - z_values[0]),
            "layers": {key: encode(values) for key, values in groups.items()},
            "layerCounts": {
                key: int(len(values)) for key, values in groups.items()
            },
        }


states = scheduler_states()
payload = []
for definition in RUNS:
    state = states.get(
        definition["jobId"],
        {"stateLabel": "Not queued", "runTime": "0:00", "reason": ""},
    )
    payload.append(
        {
            "key": definition["key"],
            "label": definition["label"],
            "config": definition["config"],
            "jobId": definition["jobId"],
            "stateLabel": state["stateLabel"],
            "runTime": state["runTime"],
            "reason": state["reason"],
            "targetTimeMs": TARGET_TIME_MS,
            "snapshot": snapshot(definition["path"]),
        }
    )

template = TEMPLATE.read_text()
marker = "/*__RUNS__*/"
if template.count(marker) != 1:
    raise RuntimeError("template must contain exactly one data marker")
OUTPUT.write_text(
    template.replace(marker, json.dumps(payload, separators=(",", ":")))
)
print(
    json.dumps(
        {
            "output": str(OUTPUT),
            "bytes": OUTPUT.stat().st_size,
            "runs": [
                {
                    "key": run["key"],
                    "state": run["stateLabel"],
                    "runTime": run["runTime"],
                    "reason": run["reason"],
                    "snapshot": None
                    if run["snapshot"] is None
                    else {
                        key: value
                        for key, value in run["snapshot"].items()
                        if key not in ("layers",)
                    },
                }
                for run in payload
            ],
        },
        indent=2,
    )
)
