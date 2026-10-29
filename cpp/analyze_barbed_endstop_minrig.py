#!/usr/bin/env python3
import collections
import h5py
import numpy as np
import sys

path = sys.argv[1]
with h5py.File(path, "r") as h5:
    events = np.asarray(h5["catch_bond/local_link_events_v3"])
    clearance = np.asarray(
        h5["actin/aa_barbed_tip_signed_clearance"][:, :, 0])
    penetration = np.asarray(
        h5["actin/aa_barbed_tip_penetration"][:, :, 0])
    force = np.asarray(h5["actin/aa_barbed_tip_endstop_force"][:, :, 0])
    chosen = np.asarray(h5["actin/aa_barbed_tip_myosin_index"][:, :, 0])
    crosslink = np.asarray(h5["actin/crosslink_ratio"][:, :, 0])

reasons = events[:, 30].astype(int)
event_codes = events[:, 4].astype(int)
print("event_count", len(events))
print("event_codes", dict(sorted(collections.Counter(event_codes).items())))
print("detachment_reasons", dict(sorted(collections.Counter(reasons).items())))
for reason in (1, 2, 3, 4):
    rows = events[reasons == reason]
    print("reason", reason, "count", len(rows))
    for row in rows[:10]:
        print(" event", {
            "link": int(row[0]), "actins": (int(row[1]), int(row[2])),
            "step": int(row[3]), "event": int(row[4]),
            "old_status": int(row[5]), "new_status": int(row[6]),
            "connector_length": float(row[21]),
            "probability": float(row[31]), "draw": float(row[32]),
        })

finite = np.isfinite(clearance)
min_flat = np.argmin(np.where(finite, clearance, np.inf))
frame, actin = np.unravel_index(min_flat, clearance.shape)
print("minimum_clearance", {
    "frame": int(frame), "step": int(frame * 10), "actin": int(actin),
    "clearance_um": float(clearance[frame, actin]),
    "penetration_um": float(penetration[frame, actin]),
    "force_pN": float(force[frame, actin]),
    "myosin": int(chosen[frame, actin]),
})
print("max_penetration_um", float(np.max(penetration)))
print("frames_any_penetration", int(np.count_nonzero(np.max(penetration, axis=1) > 0)))
print("frames_rupture_depth", int(np.count_nonzero(np.max(penetration, axis=1) >= 0.010)))
print("max_force_pN", float(np.max(force)))
print("frames_zero_crosslink_ratio", int(np.count_nonzero(np.min(crosslink, axis=1) <= 1e-12)))

for step in (125290, 125300):
    frame = step // 10
    print("snapshot", step, [
        {
            "actin": i,
            "clearance_um": float(clearance[frame, i]),
            "penetration_um": float(penetration[frame, i]),
            "force_pN": float(force[frame, i]),
            "myosin": int(chosen[frame, i]),
            "crosslink_ratio": float(crosslink[frame, i]),
        }
        for i in range(clearance.shape[1])
    ])

window = events[(events[:, 3] >= 125000) & (events[:, 3] <= 125600)]
print("events_125000_125600", [
    {"link": int(r[0]), "actins": (int(r[1]), int(r[2])),
     "step": int(r[3]), "event": int(r[4]), "reason": int(r[30])}
    for r in window
])
