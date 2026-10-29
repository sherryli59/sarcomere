#!/usr/bin/env python3
import collections
import glob
import json
import math
import os
import re

import h5py
import numpy as np

ROOT = os.environ.get(
    "BARBED_DATA",
    "/scratch/users/xsli/sarcomere/barbed_steric080_20260930/data")


def pair_angle_deg(u, v):
    dot = np.sum(u * v, axis=-1)
    dot = np.clip(dot, -1.0, 1.0)
    return 180.0 - np.degrees(np.arccos(dot))


def summarize(path):
    match = re.search(r"kappa(\d+)_seed(\d+)_5s", os.path.basename(path))
    kappa, seed = int(match.group(1)), int(match.group(2))
    with h5py.File(path, "r") as h5:
        offsets = np.asarray(h5["catch_bond/local_link_frame_offsets"])
        identities = h5["catch_bond/local_link_identity"]
        counts = offsets[:, 2].astype(int)
        nframes = len(counts)
        pair01 = np.zeros(nframes, dtype=int)
        pair23 = np.zeros(nframes, dtype=int)
        status2 = np.zeros(nframes, dtype=int)
        for frame, (_, start, count) in enumerate(offsets.astype(int)):
            rows = np.asarray(identities[start:start + count])
            if not len(rows):
                continue
            pairs = np.sort(rows[:, 2:4], axis=1)
            pair01[frame] = np.count_nonzero(np.all(pairs == [0, 1], axis=1))
            pair23[frame] = np.count_nonzero(np.all(pairs == [2, 3], axis=1))
            status2[frame] = np.count_nonzero(rows[:, 4] == 2)

        directions = np.asarray(h5["actin/direction"])
        angle01 = pair_angle_deg(directions[:, 0], directions[:, 1])
        angle23 = pair_angle_deg(directions[:, 2], directions[:, 3])
        penetration = np.asarray(
            h5["actin/aa_barbed_tip_penetration"][:, :, 0])
        clearance = np.asarray(
            h5["actin/aa_barbed_tip_signed_clearance"][:, :, 0])
        force = np.asarray(
            h5["actin/aa_barbed_tip_endstop_force"][:, :, 0])
        events = np.asarray(h5["catch_bond/local_link_events_v3"])

    tail = slice(max(1, nframes - 500), nframes)
    reasons = collections.Counter(events[:, 30].astype(int))
    finite_clearance = clearance[1:][np.isfinite(clearance[1:])]
    result = {
        "kappa": kappa,
        "seed": seed,
        "frames": nframes,
        "final_links": int(counts[-1]),
        "tail_mean_links": float(np.mean(counts[tail])),
        "tail_any_link_fraction": float(np.mean(counts[tail] > 0)),
        "tail_pair01_mean_links": float(np.mean(pair01[tail])),
        "tail_pair23_mean_links": float(np.mean(pair23[tail])),
        "tail_pair01_alive_fraction": float(np.mean(pair01[tail] > 0)),
        "tail_pair23_alive_fraction": float(np.mean(pair23[tail] > 0)),
        "tail_both_pairs_alive_fraction": float(np.mean(
            (pair01[tail] > 0) & (pair23[tail] > 0))),
        "tail_status2_mean": float(np.mean(status2[tail])),
        "final_angle01_deg": float(angle01[-1]),
        "final_angle23_deg": float(angle23[-1]),
        "tail_angle01_deg": float(np.mean(angle01[tail])),
        "tail_angle23_deg": float(np.mean(angle23[tail])),
        "min_clearance_nm": float(np.min(finite_clearance) * 1000.0),
        "max_penetration_nm": float(np.max(penetration[1:]) * 1000.0),
        "frames_any_penetration": int(np.count_nonzero(
            np.max(penetration[1:], axis=1) > 0.0)),
        "frames_rupture_penetration": int(np.count_nonzero(
            np.max(penetration[1:], axis=1) >= 0.010)),
        "max_endstop_force_pn": float(np.max(force[1:])),
        "frames_endstop_active": int(np.count_nonzero(
            np.max(force[1:], axis=1) > 0.0)),
        "reason1": int(reasons[1]),
        "reason2": int(reasons[2]),
        "reason3": int(reasons[3]),
        "reason4": int(reasons[4]),
    }
    return result


rows = [summarize(path) for path in sorted(glob.glob(os.path.join(ROOT, "*.h5")))]
print("PER_RUN")
for row in rows:
    print(json.dumps(row, sort_keys=True))

print("AGGREGATE")
metrics = [
    "final_links", "tail_mean_links", "tail_any_link_fraction",
    "tail_pair01_mean_links", "tail_pair23_mean_links", "tail_status2_mean",
    "tail_pair01_alive_fraction", "tail_pair23_alive_fraction",
    "tail_both_pairs_alive_fraction",
    "final_angle01_deg", "final_angle23_deg", "tail_angle01_deg",
    "tail_angle23_deg", "min_clearance_nm", "max_penetration_nm",
    "frames_any_penetration", "frames_rupture_penetration",
    "max_endstop_force_pn", "frames_endstop_active",
    "reason1", "reason2", "reason3", "reason4",
]
for kappa in (3, 0):
    subset = [row for row in rows if row["kappa"] == kappa]
    out = {"kappa": kappa, "n": len(subset)}
    for metric in metrics:
        values = np.asarray([row[metric] for row in subset], dtype=float)
        out[metric] = float(np.mean(values))
        out[metric + "_sd"] = float(np.std(values, ddof=1)) if len(values) > 1 else 0.0
    print(json.dumps(out, sort_keys=True))
