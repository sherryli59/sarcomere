#!/usr/bin/env python3
"""Standalone 3-D actin--myosin alignment demo with explicit stochastic heads.

The bookkeeping follows Lamson et al. (2021): every motor is an independent
crosslink with one material coordinate on each filament; filament force and
torque are sums over the explicitly bound motors.  The motor mechanics are
adapted to muscle: a weak/strong head has an extensional spring and a
state-dependent angular spring.  A separate flat-bottom proximity tether acts
through filament centers and therefore supplies translation but no torque.

Units are micrometres, seconds, picoNewtons, and pN um.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/sarcomere_matplotlib")

import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter
import numpy as np


DETACHED, WEAK, STRONG = 0, 1, 2
STATE_COLORS = {WEAK: "#8e6bbf", STRONG: "#d7308f"}


@dataclass
class Parameters:
    duration: float = 2.0
    dt: float = 2.0e-4
    save_every: int = 5
    seed: int = 2407
    n_heads: int = 48
    actin_length: float = 1.0
    myosin_length: float = 1.5
    initial_angle_deg: float = 20.0
    initial_out_of_plane_deg: float = 5.0

    # Production rod mobilities used by the C++ model.
    actin_D_parallel: float = 0.3002634333442037
    actin_D_perpendicular: float = 0.1854881703179433
    actin_D_tumble: float = 0.10096018660404828
    myosin_D_parallel: float = 0.08114427500387118
    myosin_D_perpendicular: float = 0.05223912612667491
    myosin_D_tumble: float = 0.058983847004612185
    temperature: float = 310.0
    thermal_noise_scale: float = 1.0
    actin_fixed: bool = False

    # Two-spring muscle head.  The angular rest state changes on the power
    # stroke; this replaces a processive walk along both filaments.
    head_rest_length: float = 0.055
    head_length_k_weak: float = 15.0   # pN / um
    head_length_k_strong: float = 180.0
    head_angle_k_weak: float = 0.001   # pN um / rad^2
    head_angle_k_strong: float = 0.030
    weak_angle_deg: float = 90.0
    strong_angle_deg: float = 68.0

    # Explicit-head KMC.
    capture_distance: float = 0.115
    rupture_distance: float = 0.145
    k_on: float = 80.0
    k_weak_to_strong: float = 45.0
    k_weak_off: float = 10.0
    k_strong_off: float = 3.0

    # Translational-only pair containment.  This never contributes torque.
    tether_rest_distance: float = 0.055
    tether_slack: float = 0.010
    tether_k: float = 35.0

    finite_difference_eps: float = 2.0e-6

    @property
    def beta(self) -> float:
        return 1.0 / (1.380649e-5 * self.temperature)


def normalize(v: np.ndarray) -> np.ndarray:
    n = np.linalg.norm(v)
    if n < 1e-14:
        raise ValueError("cannot normalize a zero vector")
    return v / n


def rotate_vector(v: np.ndarray, rotation_vector: np.ndarray) -> np.ndarray:
    """Rodrigues rotation using an axis-angle rotation vector."""
    angle = np.linalg.norm(rotation_vector)
    if angle < 1e-14:
        return normalize(v + np.cross(rotation_vector, v))
    axis = rotation_vector / angle
    return normalize(
        v * np.cos(angle)
        + np.cross(axis, v) * np.sin(angle)
        + axis * np.dot(axis, v) * (1.0 - np.cos(angle))
    )


def nematic_angle_deg(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.degrees(np.arccos(np.clip(abs(np.dot(a, b)), 0.0, 1.0))))


def head_geometry(
    rm: np.ndarray,
    um: np.ndarray,
    ra: np.ndarray,
    ua: np.ndarray,
    sm: np.ndarray,
    sa: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    bases = rm + sm[:, None] * um
    sites = ra + sa[:, None] * ua
    connector = sites - bases
    length = np.linalg.norm(connector, axis=1)
    unit = connector / np.maximum(length[:, None], 1e-14)
    theta = np.arccos(np.clip(unit @ um, -1.0, 1.0))
    return bases, sites, length, theta


def individual_head_energy(
    p: Parameters,
    rm: np.ndarray,
    um: np.ndarray,
    ra: np.ndarray,
    ua: np.ndarray,
    sm: np.ndarray,
    sa: np.ndarray,
    state: np.ndarray,
) -> np.ndarray:
    _, _, length, theta = head_geometry(rm, um, ra, ua, sm, sa)
    target = np.where(
        state == STRONG,
        np.radians(p.strong_angle_deg),
        np.radians(p.weak_angle_deg),
    )
    k_angle = np.where(
        state == STRONG, p.head_angle_k_strong, p.head_angle_k_weak
    )
    k_length = np.where(
        state == STRONG, p.head_length_k_strong, p.head_length_k_weak
    )
    energy = (
        0.5 * k_length * (length - p.head_rest_length) ** 2
        + 0.5 * k_angle * (theta - target) ** 2
    )
    return np.where(state != DETACHED, energy, 0.0)


def total_head_energy(
    p: Parameters,
    rm: np.ndarray,
    um: np.ndarray,
    ra: np.ndarray,
    ua: np.ndarray,
    sm: np.ndarray,
    sa: np.ndarray,
    state: np.ndarray,
) -> float:
    return float(np.sum(individual_head_energy(p, rm, um, ra, ua, sm, sa, state)))


def head_forces_and_torques(
    p: Parameters,
    rm: np.ndarray,
    um: np.ndarray,
    ra: np.ndarray,
    ua: np.ndarray,
    sm: np.ndarray,
    sa: np.ndarray,
    state: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Numerical energy gradient for exact reciprocal forces and torques."""
    eps = p.finite_difference_eps
    fm = np.zeros(3)
    tm = np.zeros(3)
    ta = np.zeros(3)
    axes = np.eye(3)

    for k in range(3):
        shift = eps * axes[k]
        ep = total_head_energy(p, rm + shift, um, ra, ua, sm, sa, state)
        em = total_head_energy(p, rm - shift, um, ra, ua, sm, sa, state)
        fm[k] = -(ep - em) / (2.0 * eps)

        up = rotate_vector(um, eps * axes[k])
        un = rotate_vector(um, -eps * axes[k])
        ep = total_head_energy(p, rm, up, ra, ua, sm, sa, state)
        em = total_head_energy(p, rm, un, ra, ua, sm, sa, state)
        tm[k] = -(ep - em) / (2.0 * eps)

        up = rotate_vector(ua, eps * axes[k])
        un = rotate_vector(ua, -eps * axes[k])
        ep = total_head_energy(p, rm, um, ra, up, sm, sa, state)
        em = total_head_energy(p, rm, um, ra, un, sm, sa, state)
        ta[k] = -(ep - em) / (2.0 * eps)

    # Translation invariance gives the equal-and-opposite force exactly.
    return fm, -fm, tm, ta


def translational_tether_force(
    p: Parameters, rm: np.ndarray, um: np.ndarray, ra: np.ndarray, ua: np.ndarray
) -> np.ndarray:
    """Force on myosin from a center-applied flat-bottom proximity tether."""
    target = ua if np.dot(ua, um) >= 0.0 else -ua
    common = normalize(um + target)
    separation = ra - rm
    radial = separation - np.dot(separation, common) * common
    distance = np.linalg.norm(radial)
    if distance < 1e-14:
        return np.zeros(3)
    extension = distance - p.tether_rest_distance
    loaded = 0.0
    if extension > p.tether_slack:
        loaded = extension - p.tether_slack
    elif extension < -p.tether_slack:
        loaded = extension + p.tether_slack
    return p.tether_k * loaded * radial / distance


def brownian_translate(
    rng: np.random.Generator,
    p: Parameters,
    center: np.ndarray,
    axis: np.ndarray,
    force: np.ndarray,
    d_parallel: float,
    d_perpendicular: float,
) -> np.ndarray:
    f_parallel = np.dot(force, axis) * axis
    f_perpendicular = force - f_parallel
    drift = p.beta * p.dt * (
        d_parallel * f_parallel + d_perpendicular * f_perpendicular
    )
    gaussian = rng.normal(size=3)
    g_parallel = np.dot(gaussian, axis) * axis
    g_perpendicular = gaussian - g_parallel
    noise = p.thermal_noise_scale * np.sqrt(2.0 * p.dt) * (
        np.sqrt(d_parallel) * g_parallel
        + np.sqrt(d_perpendicular) * g_perpendicular
    )
    return center + drift + noise


def brownian_rotate(
    rng: np.random.Generator,
    p: Parameters,
    axis: np.ndarray,
    torque: np.ndarray,
    d_tumble: float,
) -> np.ndarray:
    torque_perp = torque - np.dot(torque, axis) * axis
    gaussian = rng.normal(size=3)
    gaussian -= np.dot(gaussian, axis) * axis
    rotation = (
        p.beta * d_tumble * p.dt * torque_perp
        + p.thermal_noise_scale * np.sqrt(2.0 * d_tumble * p.dt) * gaussian
    )
    return rotate_vector(axis, rotation)


def attachment_candidates(
    p: Parameters,
    rm: np.ndarray,
    um: np.ndarray,
    ra: np.ndarray,
    ua: np.ndarray,
    sm: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    bases = rm + sm[:, None] * um
    sa = np.clip((bases - ra) @ ua, -0.5 * p.actin_length, 0.5 * p.actin_length)
    sites = ra + sa[:, None] * ua
    distance = np.linalg.norm(sites - bases, axis=1)
    return sa, distance


def kmc_step(
    rng: np.random.Generator,
    p: Parameters,
    rm: np.ndarray,
    um: np.ndarray,
    ra: np.ndarray,
    ua: np.ndarray,
    sm: np.ndarray,
    sa: np.ndarray,
    state: np.ndarray,
) -> None:
    candidate_sa, candidate_distance = attachment_candidates(p, rm, um, ra, ua, sm)
    detached = state == DETACHED
    trial_state = np.full_like(state, WEAK)
    trial_energy = individual_head_energy(
        p, rm, um, ra, ua, sm, candidate_sa, trial_state
    )
    # Lamson-style Boltzmann gating for formation of the second attachment.
    on_rate = p.k_on * np.exp(-np.minimum(p.beta * trial_energy, 50.0))
    attach = (
        detached
        & (candidate_distance <= p.capture_distance)
        & (rng.random(len(state)) < -np.expm1(-on_rate * p.dt))
    )
    state[attach] = WEAK
    sa[attach] = candidate_sa[attach]

    weak_indices = np.flatnonzero(state == WEAK)
    if len(weak_indices):
        total = p.k_weak_to_strong + p.k_weak_off
        event = rng.random(len(weak_indices)) < -np.expm1(-total * p.dt)
        choice = rng.random(len(weak_indices))
        selected = weak_indices[event]
        selected_choice = choice[event]
        to_strong = selected[selected_choice < p.k_weak_to_strong / total]
        to_off = selected[selected_choice >= p.k_weak_to_strong / total]
        state[to_strong] = STRONG
        state[to_off] = DETACHED

    strong = state == STRONG
    detach_strong = strong & (
        rng.random(len(state)) < -np.expm1(-p.k_strong_off * p.dt)
    )
    state[detach_strong] = DETACHED

    # A bound head cannot stretch without limit.
    _, _, length, _ = head_geometry(rm, um, ra, ua, sm, sa)
    state[(state != DETACHED) & (length > p.rupture_distance)] = DETACHED


def simulate(p: Parameters) -> dict[str, np.ndarray]:
    rng = np.random.default_rng(p.seed)
    ua = np.array([1.0, 0.0, 0.0])
    angle = np.radians(p.initial_angle_deg)
    out = np.radians(p.initial_out_of_plane_deg)
    um = normalize(np.array([np.cos(angle) * np.cos(out), np.sin(angle) * np.cos(out), np.sin(out)]))
    rm = np.array([0.0, 0.0, 0.0])
    ra = np.array([-0.28, p.tether_rest_distance, 0.0])

    bare_half = 0.05 * p.myosin_length
    sm = -(bare_half + (np.arange(p.n_heads) + 0.5) / p.n_heads *
           (0.5 * p.myosin_length - bare_half))
    sa = np.zeros(p.n_heads)
    state = np.zeros(p.n_heads, dtype=np.int8)

    n_steps = int(round(p.duration / p.dt))
    n_frames = n_steps // p.save_every + 1
    result = {
        "time": np.zeros(n_frames),
        "myosin_center": np.zeros((n_frames, 3)),
        "myosin_axis": np.zeros((n_frames, 3)),
        "actin_center": np.zeros((n_frames, 3)),
        "actin_axis": np.zeros((n_frames, 3)),
        "state": np.zeros((n_frames, p.n_heads), dtype=np.int8),
        "sm": np.repeat(sm[None, :], n_frames, axis=0),
        "sa": np.zeros((n_frames, p.n_heads)),
        "force_myosin": np.zeros((n_frames, 3)),
        "torque_myosin": np.zeros((n_frames, 3)),
        "tether_force": np.zeros((n_frames, 3)),
        "angle_deg": np.zeros(n_frames),
    }

    def save(frame: int, time: float, fm: np.ndarray, tm: np.ndarray, ft: np.ndarray) -> None:
        result["time"][frame] = time
        result["myosin_center"][frame] = rm
        result["myosin_axis"][frame] = um
        result["actin_center"][frame] = ra
        result["actin_axis"][frame] = ua
        result["state"][frame] = state
        result["sa"][frame] = sa
        result["force_myosin"][frame] = fm
        result["torque_myosin"][frame] = tm
        result["tether_force"][frame] = ft
        result["angle_deg"][frame] = nematic_angle_deg(um, ua)

    zeros = np.zeros(3)
    save(0, 0.0, zeros, zeros, zeros)
    frame = 1
    for step in range(1, n_steps + 1):
        kmc_step(rng, p, rm, um, ra, ua, sm, sa, state)
        fm, fa, tm, ta = head_forces_and_torques(p, rm, um, ra, ua, sm, sa, state)
        ft = translational_tether_force(p, rm, um, ra, ua)
        fm += ft
        fa -= ft

        rm = brownian_translate(
            rng, p, rm, um, fm, p.myosin_D_parallel, p.myosin_D_perpendicular
        )
        um = brownian_rotate(rng, p, um, tm, p.myosin_D_tumble)
        if not p.actin_fixed:
            ra = brownian_translate(
                rng, p, ra, ua, fa, p.actin_D_parallel, p.actin_D_perpendicular
            )
            ua = brownian_rotate(rng, p, ua, ta, p.actin_D_tumble)

        if step % p.save_every == 0:
            save(frame, step * p.dt, fm, tm, ft)
            frame += 1
    return result


def rod_endpoints(center: np.ndarray, axis: np.ndarray, length: float) -> np.ndarray:
    return np.vstack((center - 0.5 * length * axis, center + 0.5 * length * axis))


def render_animation(p: Parameters, data: dict[str, np.ndarray], output: Path, fps: int) -> None:
    fig = plt.figure(figsize=(13.5, 6.2), facecolor="white")
    writer = FFMpegWriter(
        fps=fps, codec="libx264",
        extra_args=("-pix_fmt", "yuv420p", "-crf", "18", "-threads", "1"),
    )
    frames = len(data["time"])
    stride = max(1, int(np.ceil(frames / 600)))
    render_frames = list(range(0, frames, stride))
    if render_frames[-1] != frames - 1:
        render_frames.append(frames - 1)

    output.parent.mkdir(parents=True, exist_ok=True)
    with writer.saving(fig, str(output), dpi=120):
        for number, frame in enumerate(render_frames):
            fig.clear()
            ax = fig.add_subplot(1, 2, 1, projection="3d")
            metrics = fig.add_subplot(1, 2, 2)
            rm, um = data["myosin_center"][frame], data["myosin_axis"][frame]
            ra, ua = data["actin_center"][frame], data["actin_axis"][frame]
            mp = rod_endpoints(rm, um, p.myosin_length)
            ap = rod_endpoints(ra, ua, p.actin_length)
            ax.plot(*mp.T, lw=12, color="#2171b5", solid_capstyle="round")
            ax.plot(*ap.T, lw=9, color="#ef5a29", solid_capstyle="round")

            bound = np.flatnonzero(data["state"][frame] != DETACHED)
            for h in bound:
                base = rm + data["sm"][frame, h] * um
                site = ra + data["sa"][frame, h] * ua
                state = int(data["state"][frame, h])
                ax.plot(*np.vstack((base, site)).T, color=STATE_COLORS[state],
                        lw=2.2 if state == STRONG else 1.4, alpha=0.9)
                ax.scatter(*base, color=STATE_COLORS[state], s=10, depthshade=False)

            all_centers = np.vstack((data["myosin_center"], data["actin_center"]))
            midpoint = np.mean(all_centers, axis=0)
            ax.set_xlim(midpoint[0] - 0.9, midpoint[0] + 0.9)
            ax.set_ylim(midpoint[1] - 0.55, midpoint[1] + 0.55)
            ax.set_zlim(midpoint[2] - 0.45, midpoint[2] + 0.45)
            ax.set_box_aspect((1.8, 1.1, 0.9))
            ax.set_proj_type("ortho")
            ax.view_init(elev=24, azim=-67)
            ax.set_xlabel("x (µm)")
            ax.set_ylabel("y (µm)")
            ax.set_zlabel("z (µm)")
            weak = np.count_nonzero(data["state"][frame] == WEAK)
            strong = np.count_nonzero(data["state"][frame] == STRONG)
            ax.set_title(
                f"Explicit two-spring myosin heads\nweak={weak}, strong={strong}"
            )

            t = data["time"][: frame + 1]
            metrics.plot(t, data["angle_deg"][: frame + 1], color="#1f2937",
                         lw=2.4, label="nematic misalignment")
            counts_weak = np.sum(data["state"][: frame + 1] == WEAK, axis=1)
            counts_strong = np.sum(data["state"][: frame + 1] == STRONG, axis=1)
            count_axis = metrics.twinx()
            count_axis.plot(t, counts_weak, color=STATE_COLORS[WEAK], lw=1.4,
                            label="weak heads")
            count_axis.plot(t, counts_strong, color=STATE_COLORS[STRONG], lw=1.8,
                            label="strong heads")
            metrics.set_xlim(0.0, p.duration)
            metrics.set_ylim(0.0, max(p.initial_angle_deg * 1.35,
                                      np.max(data["angle_deg"]) * 1.08))
            count_axis.set_ylim(0, p.n_heads)
            metrics.set_xlabel("time (s)")
            metrics.set_ylabel("misalignment (degrees)")
            count_axis.set_ylabel("bound heads")
            metrics.grid(alpha=0.22)
            lines = metrics.lines + count_axis.lines
            metrics.legend(lines, [line.get_label() for line in lines], loc="upper right")
            fig.suptitle(f"t = {data['time'][frame]:.3f} s", fontsize=15, weight="bold")
            fig.subplots_adjust(left=0.03, right=0.93, bottom=0.10, top=0.88, wspace=0.18)
            writer.grab_frame(facecolor="white")
            if number % 50 == 0:
                print(f"rendered {number + 1}/{len(render_frames)}", flush=True)
    plt.close(fig)


def plot_summary(p: Parameters, data: dict[str, np.ndarray], output: Path) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(8, 6.5), sharex=True)
    axes[0].plot(data["time"], data["angle_deg"], color="#1f2937", lw=1.8)
    axes[0].set_ylabel("misalignment (degrees)")
    axes[0].grid(alpha=0.25)
    axes[1].plot(data["time"], np.sum(data["state"] == WEAK, axis=1),
                 color=STATE_COLORS[WEAK], label="weak")
    axes[1].plot(data["time"], np.sum(data["state"] == STRONG, axis=1),
                 color=STATE_COLORS[STRONG], label="strong")
    axes[1].set_xlabel("time (s)")
    axes[1].set_ylabel("bound heads")
    axes[1].legend()
    axes[1].grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path,
                        default=Path("analysis/standalone_am_alignment_demo"))
    parser.add_argument("--seed", type=int, default=2407)
    parser.add_argument("--duration", type=float, default=2.0)
    parser.add_argument("--n-heads", type=int, default=48)
    parser.add_argument(
        "--fix-actin",
        action="store_true",
        help="Immobilize actin as a diagnostic control; both rods move by default.",
    )
    parser.add_argument("--noise-scale", type=float, default=1.0)
    parser.add_argument("--fps", type=int, default=20)
    parser.add_argument("--no-animation", action="store_true")
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()

    p = Parameters(
        duration=0.02 if args.smoke else args.duration,
        seed=args.seed,
        n_heads=args.n_heads,
        actin_fixed=args.fix_actin,
        thermal_noise_scale=args.noise_scale,
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    data = simulate(p)
    np.savez_compressed(args.output_dir / "trajectory.npz", **data)
    with (args.output_dir / "parameters.json").open("w") as handle:
        json.dump(asdict(p), handle, indent=2)
    plot_summary(p, data, args.output_dir / "summary.png")
    if not args.no_animation:
        render_animation(p, data, args.output_dir / "alignment_3d.mp4", args.fps)

    state = data["state"]
    print(f"initial angle: {data['angle_deg'][0]:.3f} deg")
    print(f"final angle:   {data['angle_deg'][-1]:.3f} deg")
    print(f"minimum angle: {np.min(data['angle_deg']):.3f} deg")
    print(f"mean weak heads:   {np.mean(np.sum(state == WEAK, axis=1)):.2f}")
    print(f"mean strong heads: {np.mean(np.sum(state == STRONG, axis=1)):.2f}")


if __name__ == "__main__":
    main()
