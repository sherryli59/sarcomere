#!/usr/bin/env python3
"""Summarize connectivity, periodic winding, M/Z order, and motor state."""

from __future__ import annotations

import argparse
from collections import Counter, deque
import csv
from pathlib import Path
import re

import h5py
import numpy as np


def components(n: int, pairs: set[tuple[int, int]]):
    parent = np.arange(n)
    size = np.ones(n, dtype=int)
    degree = np.zeros(n, dtype=int)

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = int(parent[x])
        return x

    for a, b in pairs:
        degree[a] += 1; degree[b] += 1
        ra, rb = find(a), find(b)
        if ra == rb:
            continue
        if size[ra] < size[rb]:
            ra, rb = rb, ra
        parent[rb] = ra; size[ra] += size[rb]
    roots = np.array([find(i) for i in range(n)], dtype=int)
    counts = np.bincount(roots, minlength=n)
    return roots, counts, degree


def winding(centers: np.ndarray, box: np.ndarray,
            pairs: set[tuple[int, int]], roots: np.ndarray) -> np.ndarray:
    adjacency = [[] for _ in range(len(centers))]
    for a, b in pairs:
        delta = centers[b] - centers[a]
        delta -= box * np.rint(delta / box)
        adjacency[a].append((b, delta)); adjacency[b].append((a, -delta))
    assigned = np.full((len(centers), 3), np.nan)
    axes = np.zeros(3, dtype=bool)
    for seed in range(len(centers)):
        if np.isfinite(assigned[seed, 0]) or not adjacency[seed]:
            continue
        assigned[seed] = centers[seed]
        queue = deque([seed])
        while queue:
            a = queue.popleft()
            for b, delta in adjacency[a]:
                proposal = assigned[a] + delta
                if not np.isfinite(assigned[b, 0]):
                    assigned[b] = proposal; queue.append(b)
                else:
                    residual = proposal - assigned[b]
                    integer = np.rint(residual / box).astype(int)
                    if np.max(np.abs(residual - integer * box)) < 1e-6:
                        axes |= integer != 0
    return axes


def band_order(x: np.ndarray, box_x: float) -> float:
    """Match the established four-run M/Z analysis (best harmonic 1--4)."""
    if len(x) == 0:
        return 0.0
    wrapped = np.asarray(x, dtype=float) % box_x
    harmonics = np.arange(1, 5, dtype=float)
    phase = np.exp(2j * np.pi * wrapped[:, None] * harmonics[None, :] / box_x)
    return float(np.max(np.abs(phase.mean(axis=0))))


def local_links(handle: h5py.File, frame: int):
    offset = np.asarray(handle["catch_bond/local_link_frame_offsets"][frame], int)
    identity = np.asarray(handle["catch_bond/local_link_identity"][
        offset[1]: offset[1] + offset[2]], int)
    raw = [tuple(sorted((int(r[2]), int(r[3])))) for r in identity]
    all_pairs = set(raw)
    mature = {tuple(sorted((int(r[2]), int(r[3])))) for r in identity if r[4] == 2}
    status1 = int(np.count_nonzero(identity[:, 4] == 1)) if len(identity) else 0
    status2 = int(np.count_nonzero(identity[:, 4] == 2)) if len(identity) else 0
    return raw, all_pairs, mature, status1, status2


def zdisc_capture(centers, directions, box, actin_length, pairs, roots, counts):
    """Largest component within each Z-disc, as a fraction of that Z-disc.

    Z-discs are located from the dominant harmonic of the barbed-end axial
    distribution (the barbed end is the MINUS tip: centre - L/2 * direction).
    Returns (n_zdiscs, mean_capture, max_capture, span_y, span_z) where the
    spans are the transverse extent of the largest component as a fraction of
    the box.
    """
    n = len(centers)
    barbed = (centers[:, 0] - 0.5 * actin_length * directions[:, 0]) % box[0]
    harmonics = np.arange(1, 9, dtype=float)
    amps = np.abs(np.exp(2j * np.pi * barbed[:, None] * harmonics[None, :]
                         / box[0]).mean(axis=0))
    k = int(harmonics[int(np.argmax(amps))])
    phase = np.angle(np.exp(2j * np.pi * k * barbed / box[0]).mean())
    zid = np.floor(((barbed / box[0] * k) - phase / (2 * np.pi) + 0.5) % k).astype(int)

    captures = []
    for z in range(k):
        member = zid == z
        idx = np.flatnonzero(member)
        if len(idx) < 10:
            continue
        remap = {int(j): q for q, j in enumerate(idx)}
        sub = {(remap[a], remap[b]) for a, b in pairs if member[a] and member[b]}
        _, sub_counts, _ = components(len(idx), sub)
        captures.append(float(sub_counts.max()) / len(idx))
    if not captures:
        return k, 0.0, 0.0, 0.0, 0.0

    big_root = max(set(roots.tolist()), key=lambda r: counts[r])
    big = np.array([roots[i] == big_root for i in range(n)])
    spans = []
    for axis in (1, 2):
        if not big.any():
            spans.append(0.0); continue
        wrapped = centers[big, axis] % box[axis]
        hist, _ = np.histogram(wrapped, bins=20, range=(0.0, box[axis]))
        spans.append(float((hist > 0).sum()) / 20.0)
    return (k, float(np.mean(captures)), float(np.max(captures)),
            spans[0], spans[1])


def analyze(path: Path, tail_fraction: float) -> dict[str, object]:
    with h5py.File(path, "r", locking=False) as h:
        meta = h["meta"].attrs
        committed = int(h["meta/committed_frames"][0])
        keys = ["state/current_step", "actin/center", "actin/direction", "myosin/center",
                "actin/cb_status", "catch_bond/local_link_frame_offsets"]
        frames = min(committed, *(h[k].shape[0] for k in keys))
        start = max(0, int((1.0 - tail_fraction) * frames))
        box = np.array([meta[k] for k in ("box_x", "box_y", "box_z")], float)
        actin_length = float(meta["actin_length"])
        records = []
        for frame in range(frames):
            centers = np.asarray(h["actin/center"][frame], float)
            directions = np.asarray(h["actin/direction"][frame], float)
            myosin = np.asarray(h["myosin/center"][frame], float)
            raw, pairs, mature, status1, status2 = local_links(h, frame)
            roots, counts, degree = components(len(centers), pairs)
            roots2, counts2, _ = components(len(centers), mature)
            axes = winding(centers, box, pairs, roots)
            axes2 = winding(centers, box, mature, roots2)
            plus_tips = centers[:, 0] + 0.5 * actin_length * directions[:, 0]
            minus_tips = centers[:, 0] - 0.5 * actin_length * directions[:, 0]
            records.append({
                "m_order": band_order(myosin[:, 0], box[0]),
                "z_order": max(band_order(plus_tips, box[0]),
                               band_order(minus_tips, box[0])),
                "lcc": float(np.max(counts) / len(centers)),
                "status2_lcc": float(np.max(counts2) / len(centers)),
                "perp_any": bool(axes[1] or axes[2]),
                "perp_both": bool(axes[1] and axes[2]),
                "status2_perp_any": bool(axes2[1] or axes2[2]),
                "status2_perp_both": bool(axes2[1] and axes2[2]),
                "unique_pairs": len(pairs), "molecular_links": len(raw),
                "molecules_per_unique_edge": len(raw) / max(1, len(pairs)),
                "duplicate_fraction": 1.0 - len(pairs) / len(raw) if raw else 0.0,
                "status1_links": status1,
                "status2_links": status2,
                "mean_degree": float(np.mean(degree)),
                **dict(zip(("zdisc_count", "zdisc_capture_mean",
                            "zdisc_capture_max", "lcc_span_y", "lcc_span_z"),
                           zdisc_capture(centers, directions, box,
                                         actin_length, pairs, roots, counts))),
            })
        tail = records[start:]
        mean = lambda key: float(np.mean([r[key] for r in tail]))
        steps = np.asarray(h["state/current_step"][:frames, 0], float)
        time = steps * float(meta["dt"])
        median = lambda key: float(np.median([r[key] for r in tail]))

        def slope(key: str) -> float:
            values = np.asarray([r[key] for r in tail], float)
            if len(values) < 2 or time[-1] == time[start]:
                return float("nan")
            return float(np.polyfit(time[start:], values, 1)[0])

        def first_time(key: str) -> float:
            for index, record in enumerate(records):
                if record[key]: return float(time[index])
            return float("nan")
        tail_duration = max(float(time[-1] - time[start]), float(meta["dt"]))
        strong = np.asarray(h["actin/am_strong_count"][:frames], float).reshape(frames, -1).sum(1)
        on = (np.asarray(h["myosin/am_on_count_minus"][:frames], float).reshape(frames, -1).sum(1) +
              np.asarray(h["myosin/am_on_count_plus"][:frames], float).reshape(frames, -1).sum(1))
        off = (np.asarray(h["myosin/am_off_count_minus"][:frames], float).reshape(frames, -1).sum(1) +
               np.asarray(h["myosin/am_off_count_plus"][:frames], float).reshape(frames, -1).sum(1))
        q = np.asarray(h["actin/am_q"][:frames], float).reshape(frames, -1).mean(1)
        weak = np.asarray(h["actin/am_weak_count"][:frames], float).reshape(frames, -1).sum(1)
        tension = np.concatenate([
            np.asarray(h["myosin/am_tension_minus"][:frames], float).reshape(frames, -1),
            np.asarray(h["myosin/am_tension_plus"][:frames], float).reshape(frames, -1),
        ], axis=1)
        tail_tension = tension[start:].reshape(-1)

        def cumulative_rate(species: str, feature: str) -> float:
            key = f"{species}/{feature}"
            if key not in h: return float("nan")
            values = np.asarray(h[key][:frames], float)
            return float((values[-1].sum() - values[start].sum()) / tail_duration)
        turnover = np.asarray(h["actin/turnover_events"], float) if "actin/turnover_events" in h else np.empty((0,20))
        label = path.stem
        # Naming-agnostic: the arm is the token immediately before _seed.
        # Works for A800_M160_<arm>_seed... and A1600_M320_Lx6_<arm>_seed...
        match = re.search(r"_([^_]+)_seed\d+_", label)
        arm = match.group(1) if match else label
        return {
            "label": label, "arm": arm, "path": str(path), "frames": frames,
            "final_time_s": float(time[-1]),
            "turnover_rate": float(meta.get("actin_turnover_rate", 0.0)),
            "thick_feedback": int(meta.get("am_thick_feedback_enabled", 0)),
            "thick_recruit_rate": float(meta.get("am_thick_recruit_rate", 0.0)),
            "thick_fold_rate": float(meta.get("am_thick_fold_rate", 0.0)),
            "thick_tension_scale": float(meta.get("am_thick_tension_scale", 0.0)),
            "tail_m_order": mean("m_order"), "tail_z_order": mean("z_order"),
            "tail_m_order_median": median("m_order"), "tail_z_order_median": median("z_order"),
            "tail_m_order_slope_per_s": slope("m_order"),
            "tail_z_order_slope_per_s": slope("z_order"),
            "tail_window_duration_s": tail_duration,
            "tail_lcc": mean("lcc"), "tail_status2_lcc": mean("status2_lcc"),
            "tail_perp_any_persistence": mean("perp_any"),
            "tail_perp_both_persistence": mean("perp_both"),
            "tail_status2_perp_any_persistence": mean("status2_perp_any"),
            "tail_status2_perp_both_persistence": mean("status2_perp_both"),
            "first_perp_any_s": first_time("perp_any"),
            "first_perp_both_s": first_time("perp_both"),
            "first_status2_perp_any_s": first_time("status2_perp_any"),
            "first_status2_perp_both_s": first_time("status2_perp_both"),
            "tail_unique_pairs": mean("unique_pairs"),
            "tail_molecules_per_unique_edge": mean("molecules_per_unique_edge"),
            "tail_status1_links": mean("status1_links"),
            "tail_status2_links": mean("status2_links"),
            "tail_duplicate_fraction": mean("duplicate_fraction"),
            "tail_mean_unique_degree": mean("mean_degree"),
            "tail_zdisc_count": median("zdisc_count"),
            "tail_zdisc_capture_mean": mean("zdisc_capture_mean"),
            "tail_zdisc_capture_max": mean("zdisc_capture_max"),
            "tail_lcc_span_y": mean("lcc_span_y"),
            "tail_lcc_span_z": mean("lcc_span_z"),
            "tail_mean_q": float(np.mean(q[start:])),
            "tail_strong_heads": float(np.mean(strong[start:])),
            "tail_on_fraction": float(np.mean(on[start:] / np.maximum(1.0, on[start:] + off[start:]))),
            "tail_weak_heads": float(np.mean(weak[start:])),
            "tail_tension_p50": float(np.quantile(tail_tension, 0.50)),
            "tail_tension_p90": float(np.quantile(tail_tension, 0.90)),
            "tail_actin_force_cap_events_per_s": cumulative_rate("actin", "cumulative_force_cap_events"),
            "tail_myosin_force_cap_events_per_s": cumulative_rate("myosin", "cumulative_force_cap_events"),
            "tail_hard_contact_events_per_s": cumulative_rate("myosin", "cumulative_hard_contact_events"),
            "tail_tip_projection_events_per_s": cumulative_rate("myosin", "cumulative_aa_zdisc_tip_projection_events"),
            "successful_turnovers": int(np.count_nonzero(turnover[:, 3] > 0.5)),
            "failed_turnovers": int(np.count_nonzero(turnover[:, 3] <= 0.5)),
        }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("inputs", nargs="+", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--tail-fraction", type=float, default=0.2)
    args = parser.parse_args()
    rows = [analyze(path, args.tail_fraction) for path in args.inputs]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    print(args.output)


if __name__ == "__main__":
    main()
