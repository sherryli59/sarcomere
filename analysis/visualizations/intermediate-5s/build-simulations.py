#!/usr/bin/env python3
"""Extract committed sarcomere observables and inject them into the viewer."""

from __future__ import annotations

import argparse
from datetime import datetime
import json
import os
from itertools import product
from pathlib import Path

os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

import h5py
import numpy as np


CORE = (
    "state/current_step",
    "actin/center",
    "actin/direction",
    "myosin/center",
    "myosin/direction",
)


def band_order(x: np.ndarray, box_x: float) -> np.ndarray:
    wrapped = np.asarray(x, dtype=float) % box_x
    harmonics = np.arange(1, 5, dtype=float)
    phase = np.exp(
        2j * np.pi * wrapped[..., None] * harmonics[None, None, :] / box_x
    ).mean(axis=1)
    return np.max(np.abs(phase), axis=1)


def clip_periodic_segments(
    starts: np.ndarray, ends: np.ndarray, box: np.ndarray
) -> list[list[float]]:
    lower, upper = -0.5 * box, 0.5 * box
    output: list[list[float]] = []
    for start, end in zip(np.asarray(starts), np.asarray(ends)):
        choices: list[list[float]] = []
        for axis in range(3):
            offsets = [0.0]
            if min(start[axis], end[axis]) < lower[axis]:
                offsets.append(float(box[axis]))
            if max(start[axis], end[axis]) > upper[axis]:
                offsets.append(float(-box[axis]))
            choices.append(offsets)
        for offset in product(*choices):
            first = start + np.asarray(offset)
            delta = end - start
            enter, leave = 0.0, 1.0
            valid = True
            for axis in range(3):
                if abs(delta[axis]) < 1e-14:
                    if first[axis] < lower[axis] or first[axis] > upper[axis]:
                        valid = False
                    continue
                t0 = (lower[axis] - first[axis]) / delta[axis]
                t1 = (upper[axis] - first[axis]) / delta[axis]
                enter = max(enter, min(t0, t1))
                leave = min(leave, max(t0, t1))
            if valid and leave > enter + 1e-12:
                clipped_start = first + enter * delta
                clipped_end = first + leave * delta
                row = np.concatenate((clipped_start, clipped_end))
                output.append(np.round(row, 4).tolist())
    return output


def committed_frames(handle: h5py.File) -> int:
    common = min(handle[name].shape[0] for name in CORE)
    if "meta/committed_frames" in handle:
        common = min(common, int(handle["meta/committed_frames"][0]))
    return common


def read_run(path: Path, key: str, label: str, state: str) -> tuple[dict, dict | None]:
    with h5py.File(path, "r", locking=False) as handle:
        frame_count = committed_frames(handle)
        if frame_count < 1:
            raise RuntimeError(f"no committed frames in {path}")
        meta = handle["meta"].attrs
        box = np.array(
            [meta["box_x"], meta["box_y"], meta["box_z"]], dtype=float
        )
        dt = float(meta["dt"])
        actin_length = float(meta["actin_length"])
        myosin_length = float(meta["myosin_length"])
        steps = np.asarray(handle["state/current_step"][:frame_count, 0], dtype=int)
        actin_center = np.asarray(handle["actin/center"][:frame_count], dtype=float)
        actin_direction = np.asarray(
            handle["actin/direction"][:frame_count], dtype=float
        )
        myosin_center = np.asarray(
            handle["myosin/center"][:frame_count], dtype=float
        )

        m_order = band_order(myosin_center[:, :, 0], box[0])
        actin_norm = np.linalg.norm(actin_direction, axis=2)
        unit_actin = actin_direction / np.maximum(
            actin_norm[:, :, None], np.finfo(float).tiny
        )
        plus = actin_center[:, :, 0] + 0.5 * actin_length * unit_actin[:, :, 0]
        minus = actin_center[:, :, 0] - 0.5 * actin_length * unit_actin[:, :, 0]
        z_order = np.maximum(band_order(plus, box[0]), band_order(minus, box[0]))
        values = [
            {
                "x": round(float(step * dt * 1e3), 3),
                "m": round(float(m_value), 6),
                "z": round(float(z_value), 6),
            }
            for step, m_value, z_value in zip(steps, m_order, z_order)
        ]
        run = {
            "key": key,
            "label": label,
            "state": state,
            "values": values,
            "latestStep": int(steps[-1]),
            "latestTimeMs": round(float(steps[-1] * dt * 1e3), 3),
            "latestM": round(float(m_order[-1]), 6),
            "latestZ": round(float(z_order[-1]), 6),
            "frames": int(frame_count),
        }

        if key != "m320":
            return run, None

        frame = frame_count - 1
        centers_a = actin_center[frame]
        directions_a = actin_direction[frame]
        centers_m = myosin_center[frame]
        directions_m = np.asarray(handle["myosin/direction"][frame], dtype=float)
        statuses_a = np.asarray(handle["actin/cb_status"][frame, :, 0], dtype=int)

        a_half = 0.5 * actin_length * directions_a
        m_half = 0.5 * myosin_length * directions_m
        segment_groups: dict[str, list[list[float]]] = {}
        for status, group_name in ((0, "actin0"), (1, "actin1"), (2, "actin2")):
            mask = statuses_a == status
            segment_groups[group_name] = clip_periodic_segments(
                (centers_a - a_half)[mask], (centers_a + a_half)[mask], box
            )
        segment_groups["myosin"] = clip_periodic_segments(
            centers_m - m_half, centers_m + m_half, box
        )

        link_counts = {"1": 0, "2": 0}
        for status, group_name in ((1, "aa1"), (2, "aa2")):
            segment_groups[group_name] = []
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
            for status, group_name in ((1, "aa1"), (2, "aa2")):
                keep = identity[:, 4] == status
                selected_identity = identity[keep]
                selected_coordinates = coordinates[keep]
                link_counts[str(status)] = int(np.count_nonzero(keep))
                if len(selected_identity):
                    first_ids = selected_identity[:, 2]
                    second_ids = selected_identity[:, 3]
                    starts = (
                        centers_a[first_ids]
                        + selected_coordinates[:, :1] * directions_a[first_ids]
                    )
                    raw_ends = (
                        centers_a[second_ids]
                        + selected_coordinates[:, 1:] * directions_a[second_ids]
                    )
                    ends = starts + (
                        raw_ends - starts
                        - box[None, :] * np.rint((raw_ends - starts) / box[None, :])
                    )
                    segment_groups[group_name] = clip_periodic_segments(
                        starts, ends, box
                    )

        status_counts = {
            str(status): int(np.count_nonzero(statuses_a == status))
            for status in (0, 1, 2)
        }
        snapshot = {
            "frame": int(frame),
            "step": int(steps[-1]),
            "timeMs": round(float(steps[-1] * dt * 1e3), 3),
            "box": np.round(box, 4).tolist(),
            "actinStatusCounts": status_counts,
            "linkStatusCounts": link_counts,
            "myosinCount": int(len(centers_m)),
            "segments": segment_groups,
        }
        return run, snapshot


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--m320", type=Path, required=True)
    parser.add_argument("--m640", type=Path, required=True)
    parser.add_argument("--template", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    run_320, snapshot = read_run(
        args.m320, "m320", "A1600 / M320", "running"
    )
    run_640, _ = read_run(
        args.m640, "m640", "A1600 / M640", "cancelled"
    )
    payload = {
        "checkedLocal": datetime.now().astimezone().strftime("%Y-%m-%d %H:%M %Z"),
        "runs": [run_320, run_640],
        "snapshot": snapshot,
    }
    template = args.template.read_text()
    marker = "/*__SIM_DATA__*/"
    if template.count(marker) != 1:
        raise RuntimeError("template must contain exactly one data marker")
    rendered = template.replace(
        marker,
        "const SIM_DATA = " + json.dumps(payload, separators=(",", ":")) + ";",
    )
    args.output.write_text(rendered)
    print(json.dumps({
        "output": str(args.output),
        "bytes": args.output.stat().st_size,
        "m320_latest": run_320,
        "m640_latest": run_640,
        "snapshot_segment_counts": {
            key: len(value) for key, value in snapshot["segments"].items()
        },
    }, indent=2))


if __name__ == "__main__":
    main()
