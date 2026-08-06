#!/usr/bin/env python3
"""Audit the shared 4-thick/8-thin force-balanced lattice diagnostic."""
import argparse
import os
from pathlib import Path

os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/sarcomere_matplotlib")

import h5py
import matplotlib.pyplot as plt
import numpy as np


def pbc_delta(delta, box, periodic):
    result = np.array(delta, dtype=float, copy=True)
    for axis in range(3):
        if periodic[axis]:
            result[..., axis] -= box[axis] * np.round(
                result[..., axis] / box[axis])
    return result


def nematic_order(directions):
    q = 1.5 * np.einsum("...i,...j->...ij", directions, directions)
    q -= 0.5 * np.eye(3)
    return np.linalg.eigvalsh(q.mean(axis=-3))[..., -1]


def segment_distance(c1, u1, c2, u2, length, box, periodic):
    """Finite-segment distance after minimum-imaging the segment centers."""
    d = pbc_delta(c1 - c2, box, periodic)
    w = d
    a = float(u1 @ u1)
    b = float(u1 @ u2)
    c = float(u2 @ u2)
    d1 = float(u1 @ w)
    e = float(u2 @ w)
    denominator = a * c - b * b
    if denominator > 1e-14:
        s = np.clip((b * e - c * d1) / denominator,
                    -0.5 * length, 0.5 * length)
    else:
        s = np.clip(-d1 / max(a, 1e-14),
                    -0.5 * length, 0.5 * length)
    t = np.clip((b * s + e) / max(c, 1e-14),
                -0.5 * length, 0.5 * length)
    s = np.clip((b * t - d1) / max(a, 1e-14),
                -0.5 * length, 0.5 * length)
    return np.linalg.norm(w + s * u1 - t * u2)


def initial_graph(actin_center, myosin_center, target, box, periodic):
    am = []
    for half in range(2):
        for a in range(8 * half, 8 * (half + 1)):
            for m in range(4):
                d = pbc_delta(actin_center[a] - myosin_center[m],
                              box, periodic)
                transverse = np.linalg.norm(d[1:])
                if abs(transverse - target) < 1e-5:
                    am.append((a, m, half))
    aa = []
    for half in range(2):
        offset = 8 * half
        for i in range(offset, offset + 8):
            for j in range(i + 1, offset + 8):
                d = pbc_delta(actin_center[i] - actin_center[j],
                              box, periodic)
                if abs(np.linalg.norm(d[1:]) - target) < 1e-5:
                    aa.append((i, j))
    return am, aa


def feature(handle, species, name, frames, count):
    path = f"{species}/{name}"
    if path not in handle:
        return np.zeros((frames, count))
    return np.asarray(handle[path], dtype=float).reshape(frames, count)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("trajectory", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    with h5py.File(args.trajectory, "r", locking=False) as h:
        ac = np.asarray(h["actin/center"])
        au = np.asarray(h["actin/direction"])
        mc = np.asarray(h["myosin/center"])
        mu = np.asarray(h["myosin/direction"])
        steps = np.asarray(h["state/current_step"])[:, 0]
        meta = h["meta"].attrs
        dt = float(meta["dt"])
        box = np.array([meta["box_x"], meta["box_y"], meta["box_z"]],
                       dtype=float)
        periodic = np.array(
            [meta["periodic_x"], meta["periodic_y"], meta["periodic_z"]],
            dtype=bool)
        target_am = float(meta["am_preferred_centerline"])
        target_aa = float(meta["aa_preferred_centerline"])
        actin_radius = float(meta["actin_radius"])
        actin_length = float(meta["actin_length"])
        states = np.asarray(h["state/am_head_state_current"])
        bound = np.asarray(h["state/am_head_bound_actin_current"])
        pair_count = np.asarray(h["state/aa_pairs_current_count"])[:, 0]
        pairs = np.asarray(h["state/aa_pairs_current"])
        pair_status = np.asarray(h["state/aa_status_current"])[:, :, 0]
        frames = ac.shape[0]
        actin_force_caps = feature(
            h, "actin", "cumulative_force_cap_events", frames, 16)
        actin_torque_caps = feature(
            h, "actin", "cumulative_torque_cap_events", frames, 16)
        actin_rotation_caps = feature(
            h, "actin", "cumulative_rotation_cap_events", frames, 16)
        myosin_force_caps = feature(
            h, "myosin", "cumulative_force_cap_events", frames, 4)
        myosin_torque_caps = feature(
            h, "myosin", "cumulative_torque_cap_events", frames, 4)
        myosin_rotation_caps = feature(
            h, "myosin", "cumulative_rotation_cap_events", frames, 4)
        aa_s1_molecular = 0.5 * feature(
            h, "actin", "aa_local_status1_count", frames, 16).sum(axis=1)
        aa_s2_molecular = 0.5 * feature(
            h, "actin", "aa_local_status2_count", frames, 16).sum(axis=1)

    time = steps * dt
    am_edges, aa_edges = initial_graph(
        ac[0], mc[0], target_am, box, periodic)
    expected_aa = {tuple(sorted(edge)) for edge in aa_edges}

    am_distance = np.empty((frames, len(am_edges)))
    psi6 = np.empty((frames, 8))
    for f in range(frames):
        for edge, (a, m, _) in enumerate(am_edges):
            d = pbc_delta(ac[f, a] - mc[f, m], box, periodic)
            axis = mu[f, m] / np.linalg.norm(mu[f, m])
            am_distance[f, edge] = np.linalg.norm(d - axis * (d @ axis))
        for half in range(2):
            for m in range(4):
                actins = [a for a, owner, h in am_edges
                          if owner == m and h == half]
                axis = mu[f, m] / np.linalg.norm(mu[f, m])
                reference = np.array([0.0, 0.0, 1.0])
                if abs(reference @ axis) > 0.9:
                    reference = np.array([0.0, 1.0, 0.0])
                e1 = np.cross(axis, reference)
                e1 /= np.linalg.norm(e1)
                e2 = np.cross(axis, e1)
                angles = []
                for a in actins:
                    d = pbc_delta(ac[f, a] - mc[f, m], box, periodic)
                    d -= axis * (d @ axis)
                    angles.append(np.arctan2(d @ e2, d @ e1))
                psi6[f, 4 * half + m] = abs(
                    np.mean(np.exp(6j * np.asarray(angles))))

    occupied_expected = np.zeros(frames)
    unexpected_edges = np.zeros(frames)
    aa_status1_pairs = np.zeros(frames)
    aa_status2_pairs = np.zeros(frames)
    for f in range(frames):
        present = set()
        for slot in range(int(pair_count[f])):
            i, j = map(int, pairs[f, slot])
            if i < 0 or j < 0:
                continue
            edge = tuple(sorted((i, j)))
            present.add(edge)
            if pair_status[f, slot] == 1:
                aa_status1_pairs[f] += 1
            elif pair_status[f, slot] == 2:
                aa_status2_pairs[f] += 1
        occupied_expected[f] = len(present & expected_aa) / len(expected_aa)
        unexpected_edges[f] = len(present - expected_aa)

    minimum_aa = np.full(frames, np.inf)
    overlap_pairs = np.zeros(frames)
    severe_overlap_pairs = np.zeros(frames)
    for f in range(frames):
        for i in range(16):
            for j in range(i + 1, 16):
                distance = segment_distance(
                    ac[f, i], au[f, i], ac[f, j], au[f, j],
                    actin_length, box, periodic)
                minimum_aa[f] = min(minimum_aa[f], distance)
                overlap_pairs[f] += distance < 2.0 * actin_radius - 1e-8
                severe_overlap_pairs[f] += distance < 1.8 * actin_radius

    weak = np.count_nonzero((states == 1) & (bound >= 0), axis=(1, 2))
    strong = np.count_nonzero((states == 2) & (bound >= 0), axis=(1, 2))
    cap_total = (
        actin_force_caps.sum(1) + actin_torque_caps.sum(1) +
        actin_rotation_caps.sum(1) + myosin_force_caps.sum(1) +
        myosin_torque_caps.sum(1) + myosin_rotation_caps.sum(1))
    cap_increment = np.diff(cap_total, prepend=cap_total[0])

    columns = np.c_[
        time, am_distance.mean(1), am_distance.std(1), psi6.mean(1),
        occupied_expected, unexpected_edges, weak, strong,
        aa_status1_pairs, aa_status2_pairs,
        aa_s1_molecular, aa_s2_molecular,
        nematic_order(au), nematic_order(mu),
        1e3 * minimum_aa, overlap_pairs, severe_overlap_pairs,
        cap_total, cap_increment]
    names = [
        "time_s", "am_mean_um", "am_std_um", "mean_psi6",
        "expected_aa_edge_occupancy", "unexpected_aa_edges",
        "weak_am_heads", "strong_am_heads",
        "aa_status1_pairs", "aa_status2_pairs",
        "aa_status1_molecules", "aa_status2_molecules",
        "actin_nematic", "myosin_nematic",
        "minimum_aa_segment_distance_nm", "soft_core_contact_count",
        "severe_overlap_pair_count",
        "cumulative_cap_events", "new_cap_events"]
    np.savetxt(args.output / "timeseries.csv", columns, delimiter=",",
               header=",".join(names), comments="")

    fig, axes = plt.subplots(6, 1, figsize=(10, 15), sharex=True)
    axes[0].plot(time, 1e3 * am_distance.mean(1), color="#0072B2")
    axes[0].fill_between(
        time, 1e3 * (am_distance.mean(1) - am_distance.std(1)),
        1e3 * (am_distance.mean(1) + am_distance.std(1)),
        color="#0072B2", alpha=0.18)
    axes[0].axhline(1e3 * target_am, color="black", ls="--")
    axes[0].set_ylabel("AM spacing (nm)\nmean ± SD")
    axes[1].plot(time, psi6.mean(1), color="#CC79A7", label=r"$|\psi_6|$")
    axes[1].plot(time, nematic_order(au), color="#009E73",
                 label="actin nematic")
    axes[1].plot(time, nematic_order(mu), color="#E69F00",
                 label="myosin nematic")
    axes[1].set_ylim(0, 1.03)
    axes[1].set_ylabel("order")
    axes[1].legend(ncol=3)
    axes[2].plot(time, weak, color="#56B4E9", label="weak AM")
    axes[2].plot(time, strong, color="#D55E00", label="strong AM")
    axes[2].set_ylabel("AM heads")
    axes[2].legend()
    axes[3].plot(time, aa_s1_molecular, color="#E69F00",
                 label="status-1 molecules")
    axes[3].plot(time, aa_s2_molecular, color="#009E73",
                 label="status-2 molecules")
    axes[3].set_ylabel("ACTN2")
    axes[3].legend()
    axes[4].plot(time, occupied_expected, color="#5E3C99",
                 label="intended edge occupancy")
    axes[4].plot(time, unexpected_edges, color="#B2182B",
                 label="unexpected edges")
    axes[4].set_ylabel("AA topology")
    axes[4].legend()
    axes[5].plot(time, 1e3 * minimum_aa, color="#1B9E77",
                 label="minimum AA distance")
    axes[5].axhline(2e3 * actin_radius, color="black", ls="--",
                    label="hard-core diameter")
    cap_axis = axes[5].twinx()
    cap_axis.plot(time, cap_total, color="#D95F02", alpha=0.7,
                  label="cumulative caps")
    axes[5].set_ylabel("distance (nm)")
    cap_axis.set_ylabel("cap events")
    axes[5].set_xlabel("time (s)")
    axes[5].legend(loc="upper left")
    cap_axis.legend(loc="upper right")
    for axis in axes:
        axis.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(args.output / "shared_lattice_audit.png", dpi=200)

    tail = slice(max(0, int(0.75 * frames)), None)
    print(f"frames={frames}, duration={time[-1]:.6g} s")
    print(f"AM spacing={1e3*np.mean(am_distance.mean(1)[tail]):.2f} ± "
          f"{1e3*np.mean(am_distance.std(1)[tail]):.2f} nm")
    print(f"psi6={np.mean(psi6.mean(1)[tail]):.4f}; "
          f"S_actin={np.mean(nematic_order(au)[tail]):.4f}; "
          f"S_myosin={np.mean(nematic_order(mu)[tail]):.4f}")
    print(f"expected AA occupancy={np.mean(occupied_expected[tail]):.4f}; "
          f"unexpected edges={np.mean(unexpected_edges[tail]):.3f}")
    print(f"weak={np.mean(weak[tail]):.2f}; strong={np.mean(strong[tail]):.2f}; "
          f"AA1={np.mean(aa_s1_molecular[tail]):.2f}; "
          f"AA2={np.mean(aa_s2_molecular[tail]):.2f}")
    print(f"minimum AA distance={minimum_aa.min()*1e3:.3f} nm; "
          f"max soft-core contacts={int(overlap_pairs.max())}; "
          f"max severe overlaps={int(severe_overlap_pairs.max())}; "
          f"cap events={int(cap_total[-1] - cap_total[0])}")


if __name__ == "__main__":
    main()
