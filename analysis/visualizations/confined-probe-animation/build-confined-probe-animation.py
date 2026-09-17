#!/usr/bin/env python3
"""Build the interactive alignment-probe trajectory animation.

Follows VISUALIZATION_PIPELINE.md: inline HTML, base64 float32 payload, manual
decoder in the viewer, species-hue layer contract with artifact-owned colour
fallbacks.  Extends the canonical snapshot viewer with a full-trajectory
payload plus playback controls, because a single committed frame cannot show a
binding episode.

Motor-arm geometry is reconstructed exactly as compute_am_two_spring_head()
(interaction.cpp:555) does and validated against
state/am_head_strain_energy_current before anything is written.
"""
from __future__ import annotations
import argparse, base64, json, os, re, sys
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
import numpy as np, h5py

DATA = "/scratch/users/xsli/sarcomere/confined_probe_20260910/data"
JOB_ID = "42753933"
BARE_FRAC = 0.10
BEAR_FRAC = 0.5 * (1.0 - BARE_FRAC)
DETACHED, WEAK, STRONG, OFF = 0, 1, 2, 3
NSCAL = 8                      # must match NSCAL in the template


def encode(values: np.ndarray, dtype=np.float32) -> str:
    packed = np.ascontiguousarray(values, dtype=dtype)
    return base64.b64encode(packed.tobytes()).decode("ascii")


def unit(v):
    n = np.linalg.norm(v, axis=-1, keepdims=True)
    return v / np.where(n > 1e-12, n, 1.0)


def soft_threshold(x, slack):
    return np.where(x > slack, x - slack, np.where(x < -slack, x + slack, 0.0))


def anchor_axial_table(heads_per_half: int, myosin_length: float) -> np.ndarray:
    """head.anchor_axial exactly as _initialize_am_head_pools, sarcomere.cpp:2935."""
    crown_count = max(1, (heads_per_half + 2) // 3)
    out = np.empty(2 * heads_per_half)
    for half_index in (0, 1):
        half_sign = -1.0 if half_index == 0 else 1.0
        for h in range(heads_per_half):
            fraction = (h // 3 + 0.5) / crown_count
            out[half_index * heads_per_half + h] = half_sign * (
                0.5 * BARE_FRAC * myosin_length
                + fraction * BEAR_FRAC * myosin_length)
    return out


def min_image(d, box, periodic):
    return d - np.where(periodic, box * np.round(d / box), 0.0)


def build_run(path: str, stride: int = 1) -> dict:
    with h5py.File(path, "r", locking=False) as handle:
        attrs = dict(handle["meta"].attrs)
        # Mandatory rule: never read past the last committed frame.
        frames = int(min(handle["meta/committed_frames"][0],
                         handle["actin/center"].shape[0],
                         handle["state/am_head_state_current"].shape[0]))
        box = np.array([attrs["box_x"], attrs["box_y"], attrs["box_z"]], float)
        periodic = np.array([attrs["periodic_x"], attrs["periodic_y"],
                             attrs["periodic_z"]], bool)
        dt = float(attrs["dt"])
        actin_length = float(attrs["actin_length"])
        myosin_length = float(attrs["myosin_length"])
        myosin_radius = float(attrs["myosin_radius"])
        heads_per_half = int(attrs["am_heads_per_half"])
        rho_rest = float(attrs["actin_radius"] + attrs["am_surface_gap"])
        slack = float(attrs["am_tether_slack"])
        k_rho = float(attrs["am_head_transverse_k"])
        k_axial_weak = float(attrs["am_motor_k"]) * 0.1
        reach_rho = rho_rest + float(attrs["am_rupture_distance"])
        cap_half = int(attrs["am_max_heads_per_actin_half"])
        alignment_k = float(attrs["am_alignment_k"])

        # Uniform decimation, applied to every per-frame read below so the
        # frame axis stays consistent.  frameMs is recomputed from the
        # decimated step column, so the viewer clock stays correct.
        sl = slice(0, frames, max(1, int(stride)))
        step = handle["state/current_step"][sl, 0].astype(float)
        frames = len(step)
        actin_centre = handle["actin/center"][sl, 0, :]
        actin_dir = unit(handle["actin/direction"][sl, 0, :])
        myosin_centre = handle["myosin/center"][sl, 0, :]
        myosin_dir = unit(handle["myosin/direction"][sl, 0, :])
        e1_raw = np.concatenate(
            [handle[f"myosin/material_e1_{c}"][sl, 0, :] for c in "xyz"], -1)
        e2_raw = np.concatenate(
            [handle[f"myosin/material_e2_{c}"][sl, 0, :] for c in "xyz"], -1)
        head_state = handle["state/am_head_state_current"][sl, 0, :]
        head_actin = handle["state/am_head_bound_actin_current"][sl, 0, :]
        head_azimuth = handle["state/am_head_binding_azimuth_current"][sl, 0, :]
        head_target = handle["state/am_head_target_axial_current"][sl, 0, :]
        head_stroke = handle["state/am_head_stroke_progress_current"][sl, 0, :]
        head_energy = handle["state/am_head_strain_energy_current"][sl, 0, :]
        weak_agg = handle["actin/am_weak_count"][sl, 0, 0]
        strong_agg = handle["actin/am_strong_count"][sl, 0, 0]

    # Head state codes are calibrated against the aggregate counts, not assumed.
    assert np.array_equal((head_state == WEAK).sum(axis=1), weak_agg.astype(int)), \
        f"{path}: weak-count mismatch"
    assert np.array_equal((head_state == STRONG).sum(axis=1), strong_agg.astype(int)), \
        f"{path}: strong-count mismatch"

    anchor_axial = anchor_axial_table(heads_per_half, myosin_length)
    # Gram-Schmidt material frame, interaction.cpp:571-580.
    e1 = unit(e1_raw - myosin_dir * np.einsum("fk,fk->f", e1_raw, myosin_dir)[:, None])
    e2 = unit(np.cross(myosin_dir, e1))
    e2 = np.where(np.einsum("fk,fk->f", e2, e2_raw)[:, None] < 0.0, -e2, e2)

    # Everything is drawn relative to the myosin centre: the myosin is held at
    # the origin with its lab-frame orientation, so the actin is what moves.
    d_centre = min_image(actin_centre - myosin_centre, box, periodic)
    axial = np.einsum("fk,fk->f", d_centre, myosin_dir)
    perp = np.linalg.norm(d_centre - axial[:, None] * myosin_dir, axis=1)
    theta = np.degrees(np.arccos(np.clip(
        np.abs(np.einsum("fk,fk->f", actin_dir, myosin_dir)), 0.0, 1.0)))

    half_a = 0.5 * actin_length
    half_m = 0.5 * myosin_length
    actin_seg = np.concatenate(
        [d_centre - half_a * actin_dir, d_centre + half_a * actin_dir], axis=1)
    myosin_seg = np.concatenate(
        [-half_m * myosin_dir, half_m * myosin_dir], axis=1)
    barbed = d_centre - half_a * actin_dir      # barbed end is the MINUS tip

    bound = (head_state == WEAK) | (head_state == STRONG)
    fi, hi = np.nonzero(bound)
    radial_anchor = myosin_radius * (
        np.cos(head_azimuth[fi, hi])[:, None] * e1[fi]
        + np.sin(head_azimuth[fi, hi])[:, None] * e2[fi])
    myosin_attachment = anchor_axial[hi][:, None] * myosin_dir[fi] + radial_anchor
    # The C++ min-images the actin centre against the attachment, not the centre.
    d_att = min_image(actin_centre[fi] - (myosin_centre[fi] + myosin_attachment),
                      box, periodic)
    coordinate = np.clip(head_target[fi, hi], -half_a, half_a)
    connector = min_image(d_att + coordinate[:, None] * actin_dir[fi], box, periodic)
    actin_attachment = myosin_attachment + connector

    z = np.einsum("kj,kj->k", connector, actin_dir[fi])
    rho = np.linalg.norm(connector - z[:, None] * actin_dir[fi], axis=1)
    arm_length = np.linalg.norm(connector, axis=1)

    # Validate the whole reconstruction against the C++ strain energy.
    loaded = soft_threshold(rho - rho_rest, slack)
    radial_energy = 0.5 * k_rho * loaded ** 2
    is_weak = head_state[fi, hi] == WEAK
    # The strong-head penalty uses am_motor_k * am_weak_stiffness_ratio, NOT
    # am_head_transverse_k (sarcomere.cpp:3473).  The two are both 100 pN/um in
    # most configs, which hid this until the notransverse arm set the transverse
    # stiffness to zero and the validation below failed.
    mine = np.where(is_weak,
                    radial_energy + 0.5 * k_axial_weak * z ** 2,
                    radial_energy + 0.5 * k_axial_weak
                    * np.maximum(0.0, z - head_stroke[fi, hi]) ** 2)
    reference = head_energy[fi, hi]
    if len(mine):
        error = np.abs(mine - reference) / np.maximum(1e-9, np.abs(reference))
        assert error.max() < 1e-6, (
            f"{path}: arm reconstruction disagrees with "
            f"state/am_head_strain_energy_current, max rel error {error.max():.3e}")

    order = np.argsort(fi, kind="stable")
    fi, hi = fi[order], hi[order]
    arms = np.concatenate([myosin_attachment[order], actin_attachment[order]], axis=1)
    arm_state = head_state[fi, hi].astype(np.int32)
    arm_offset = np.zeros(frames + 1, dtype=np.int32)
    arm_offset[1:] = np.cumsum(np.bincount(fi, minlength=frames))

    n_weak = (head_state == WEAK).sum(axis=1).astype(float)
    n_strong = (head_state == STRONG).sum(axis=1).astype(float)
    mean_len = np.zeros(frames)
    mean_rho = np.zeros(frames)
    counts = np.bincount(fi, minlength=frames)
    np.add.at(mean_len, fi, arm_length[order] * 1e3)
    np.add.at(mean_rho, fi, rho[order] * 1e3)
    mean_len /= np.maximum(counts, 1)
    mean_rho /= np.maximum(counts, 1)

    scalars = np.stack([theta, perp, axial, n_weak, n_strong,
                        mean_len, mean_rho, step], axis=1)
    assert scalars.shape[1] == NSCAL

    # True per-axis bounding box over every drawn point, so the viewer fits the
    # data rather than a cube.  `engaged` is the same box restricted to frames
    # with at least one bound head; it is the default camera because the escaped
    # excursions are 5-6x larger and would shrink the pair to nothing.
    def bounding_box(frame_mask):
        points = [actin_seg[frame_mask, :3], actin_seg[frame_mask, 3:],
                  myosin_seg[frame_mask, :3], myosin_seg[frame_mask, 3:]]
        stacked = np.concatenate([p for p in points if len(p)], axis=0)
        lo = stacked.min(axis=0)
        hi = stacked.max(axis=0)
        pad = 0.06 * np.maximum(hi - lo, 0.1)
        return [(lo - pad).tolist(), (hi + pad).tolist()]

    engaged = (n_weak + n_strong) > 0
    bounds_all = bounding_box(np.ones(frames, bool))
    bounds_engaged = bounding_box(engaged) if engaged.any() else bounds_all

    name = os.path.basename(path)
    match = re.match(r"(.+?)_seed(\d+)\.h5", name)
    assert match, f"unparsable run name: {name}"
    config, seed = match.group(1), match.group(2)
    return dict(
        key=f"{config}|{seed}", config=config, seed=seed, jobId=JOB_ID,
        label=f"{config.replace('free_', '')} · s{seed}",
        frames=frames, frameMs=float(step[1] - step[0]) * dt * 1e3 if frames > 1 else 0.0,
        capHalf=cap_half, alignmentK=alignment_k,
        myosinLength=myosin_length, myosinRadius=myosin_radius, reachRho=reach_rho,
        bounds=bounds_all, boundsEngaged=bounds_engaged,
        actin=encode(actin_seg), myosin=encode(myosin_seg), barbed=encode(barbed),
        scalars=encode(scalars), arms=encode(arms),
        armState=encode(arm_state, np.int32), armOffset=encode(arm_offset, np.int32),
        _armCount=int(len(arm_state)),
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--runs", nargs="+", required=True,
                        help="basenames under the data directory")
    parser.add_argument("--template", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--title", required=True)
    parser.add_argument("--data-dir", default=DATA)
    parser.add_argument("--stride", type=int, default=1,
                        help="keep every Nth committed frame")
    args = parser.parse_args()

    payload = []
    for name in args.runs:
        run = build_run(os.path.join(args.data_dir, name), args.stride)
        print(f"  {name:38s} frames={run['frames']:5d} "
              f"arm-records={run['_armCount']:6d} "
              f"frame={run['frameMs']:.2f} ms  cap_half={run['capHalf']} "
              f"align_k={run['alignmentK']}")
        del run["_armCount"]
        payload.append(run)
    assert payload, "no runs built"

    # The template carries UTF-8 glyphs in its labels; only the *payload* has
    # to be ASCII-safe, which base64 guarantees and this asserts.
    for run in payload:
        for field in ("actin", "myosin", "barbed", "scalars", "arms",
                      "armState", "armOffset"):
            assert re.fullmatch(r"[A-Za-z0-9+/]*={0,2}", run[field]), (
                f"{run['key']}: layer {field} is not base64-clean")

    with open(args.template, "r", encoding="utf-8") as handle:
        template = handle.read()
    assert template.count("/*__RUNS__*/") == 1, "expected exactly one RUNS marker"
    assert template.count("__TITLE__") == 1, "expected exactly one TITLE marker"
    html = (template.replace("/*__RUNS__*/", json.dumps(payload))
                    .replace("__TITLE__", args.title))
    with open(args.output, "w", encoding="utf-8") as handle:
        handle.write(html)
    print(f"wrote {args.output}  ({os.path.getsize(args.output)/1e6:.2f} MB, "
          f"{len(payload)} runs)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
