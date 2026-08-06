#!/usr/bin/env python3
"""Build collaborator-facing intermediate-result figures from frozen runs."""

from __future__ import annotations

import csv
from pathlib import Path

import h5py
import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt


OUT = Path(__file__).resolve().parent
DATA = Path("/scratch/users/xsli/sarcomere")

PRODUCTION = {
    "seed 1": DATA / "prod_titin_serial_20260802/data/PRODALIGN120_A1280_M640_h35_100ms_seed2080201.h5",
    "seed 2": DATA / "prod_titin_serial_20260802/data/PRODALIGN120_A1280_M640_h35_100ms_seed2080202.h5",
}

HEAD_AB = {
    r"coupled transverse escape ($\theta_\perp=0.1$)":
        DATA / "prod_head_tether_ab_20260804/data/A1280_M640_qeq_coupled_seed20260804_head_tether_k100_kalign120.h5",
    r"uncoupled transverse escape ($\theta_\perp=0$)":
        DATA / "prod_split_detachment_k50_20260805/data/A1280_M640_qeq_coupled_seed20260804_split_k100_thetaperp0_kalign120_100ms.h5",
}

CURRENT_FAILURE = (
    DATA / "prod_productive_register_qeq_ramp_20260805/data/"
    "A1280_M640_productive_register_qeq_ramp_seed20260805_100ms.h5"
)

MYOMESIN = {
    0: DATA / "prod_titin_serial_20260802/data/ALIGN120_A640_M320_h35_100ms_seed2080301.h5",
    5: DATA / "prod_titin_serial_20260802/data/MMB5_A640_M320_h35_a120_100ms_seed2080301.h5",
    15: DATA / "prod_titin_serial_20260802/data/MMB15_A640_M320_h35_a120_100ms_seed2080301.h5",
    30: DATA / "prod_titin_serial_20260802/data/MMB30_A640_M320_h35_a120_100ms_seed2080301.h5",
    50: DATA / "prod_titin_serial_20260802/data/MMB50_A640_M320_h35_a120_100ms_seed2080301.h5",
    200: DATA / "prod_titin_serial_20260802/data/MMB200_A640_M320_h35_a120_100ms_seed2080301.h5",
}

BLUE = "#3366CC"
ORANGE = "#E67E22"
GREEN = "#2CA25F"
PURPLE = "#756BB1"
RED = "#C0392B"
INK = "#202124"
MUTED = "#6B7280"


def setup_style() -> None:
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "axes.edgecolor": "#9CA3AF",
        "axes.labelcolor": INK,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "text.color": INK,
        "legend.frameon": False,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "figure.dpi": 150,
    })


def save(fig: plt.Figure, stem: str) -> None:
    fig.savefig(OUT / f"{stem}.png", dpi=220, bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def nematic_series(direction: np.ndarray) -> np.ndarray:
    result = []
    for frame in direction:
        u = frame / np.linalg.norm(frame, axis=1, keepdims=True)
        q = (3.0 * u.T @ u / len(u) - np.eye(3)) / 2.0
        result.append(np.linalg.eigvalsh(q)[-1])
    return np.asarray(result)


def active_link_series(events: np.ndarray, steps: np.ndarray) -> np.ndarray:
    if len(events) == 0:
        return np.zeros_like(steps)
    order = np.argsort(events[:, 3], kind="stable")
    e = events[order]
    delta = np.where(e[:, 4] == 1, 1, np.where(e[:, 4] == 2, -1, 0))
    cumulative = np.cumsum(delta)
    index = np.searchsorted(e[:, 3], steps, side="right") - 1
    return np.where(index >= 0, cumulative[np.maximum(index, 0)], 0)


def production_series(path: Path) -> dict[str, np.ndarray]:
    with h5py.File(path, "r", locking=False) as h:
        dt = float(h["meta"].attrs["dt"])
        steps = np.asarray(h["state/current_step"][:, 0], dtype=np.int64)
        cb = np.asarray(h["actin/cb_status"][:, :, 0], dtype=int)
        strong = np.asarray(h["actin/am_strong_count"][:, :, 0], dtype=float)
        direction = np.asarray(h["actin/direction"], dtype=float)
        events = np.asarray(h["catch_bond/local_link_events"], dtype=float)
    return {
        "ms": steps * dt * 1e3,
        "active": active_link_series(events, steps),
        "status2": np.sum(cb == 2, axis=1),
        "nematic": nematic_series(direction),
        "strong_fraction": np.mean(strong > 0, axis=1),
    }


def make_scaleup_figure() -> list[dict[str, float]]:
    series = {label: production_series(path) for label, path in PRODUCTION.items()}
    fig, axes = plt.subplots(2, 2, figsize=(10.8, 7.2), sharex=True)
    metrics = [
        ("active", "Active ACTN2 links", "links"),
        ("status2", "Mature ACTN2-linked actins", "status-2 actins"),
        ("nematic", "Actin nematic order", r"$S$"),
        ("strong_fraction", "Actins carrying a strong head", "fraction"),
    ]
    colors = [BLUE, ORANGE]
    summaries: list[dict[str, float]] = []
    for ax, (key, title, ylabel) in zip(axes.flat, metrics):
        for color, (label, data) in zip(colors, series.items()):
            ax.plot(data["ms"], data[key], color=color, lw=2, label=label)
            ax.scatter(data["ms"][-1], data[key][-1], s=20, color=color, zorder=3)
        ax.set_title(title, loc="left")
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.18, lw=0.7)
    axes[1, 0].set_xlabel("time (ms)")
    axes[1, 1].set_xlabel("time (ms)")
    axes[0, 0].legend(loc="upper left")
    axes[1, 0].set_ylim(0, 1)
    axes[1, 1].set_ylim(0, 1)
    fig.suptitle("Production-scale assembly reproducibility (1,280 actins, model 19)", x=0.08, ha="left", fontsize=13)
    fig.text(0.08, 0.015, "Two independent seeds; older weak-titin configuration. This is an intermediate result, not the final calibrated model.", color=MUTED)
    fig.tight_layout(rect=(0, 0.035, 1, 0.95))
    save(fig, "01-production-scale-assembly")
    for label, data in series.items():
        summaries.append({
            "run": label,
            "active_actn2_links": float(data["active"][-1]),
            "status2_actins": float(data["status2"][-1]),
            "actin_nematic": float(data["nematic"][-1]),
            "strong_actin_fraction": float(data["strong_fraction"][-1]),
        })
    return summaries


def wrap(displacement: np.ndarray, box: np.ndarray, periodic: np.ndarray) -> np.ndarray:
    result = displacement.copy()
    for axis in range(3):
        if periodic[axis] and box[axis] > 0:
            result[axis] -= box[axis] * np.round(result[axis] / box[axis])
    return result


def final_overlap_samples(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    with h5py.File(path, "r", locking=False) as h:
        meta = h["meta"].attrs
        box = np.asarray([meta.get("box_x", 0), meta.get("box_y", 0), meta.get("box_z", 0)], dtype=float)
        periodic = np.asarray([bool(meta.get("periodic_x", 1)), bool(meta.get("periodic_y", 1)), bool(meta.get("periodic_z", 1))])
        actin_length = float(meta.get("actin_length", 1.0))
        myosin_length = float(meta.get("myosin_length", 1.5))
        ac = np.asarray(h["actin/center"][-1], dtype=float)
        au = np.asarray(h["actin/direction"][-1], dtype=float)
        au /= np.linalg.norm(au, axis=1, keepdims=True)
        mc = np.asarray(h["myosin/center"][-1], dtype=float)
        mu = np.asarray(h["myosin/direction"][-1], dtype=float)
        mu /= np.linalg.norm(mu, axis=1, keepdims=True)
        cb = np.asarray(h["actin/cb_status"][-1, :, 0], dtype=int)
        bound = np.asarray(h["state/am_head_bound_actin_current"][-1], dtype=int)
        offset = np.asarray(h["catch_bond/local_link_frame_offsets"][-1], dtype=int)
        start, count = int(offset[1]), int(offset[2])
        identity = np.asarray(h["catch_bond/local_link_identity"][start:start + count], dtype=int)

    aa = []
    for row in identity:
        i, j, status = int(row[2]), int(row[3]), int(row[4])
        if status != 2:
            continue
        ui, uj = au[i], au[j].copy()
        if ui.dot(uj) < 0:
            uj = -uj
        axis = ui + uj
        if np.linalg.norm(axis) < 1e-10:
            axis = ui.copy()
        axis /= np.linalg.norm(axis)
        d = wrap(ac[j] - ac[i], box, periodic)
        pj = d.dot(axis)
        hi = 0.5 * actin_length * abs(au[i].dot(axis))
        hj = 0.5 * actin_length * abs(au[j].dot(axis))
        overlap = min(hi, pj + hj) - max(-hi, pj - hj)
        aa.append(max(0.0, overlap))

    # One A-M sample per mature actin. Using every myosin touched by an actin
    # overcounts crowded secondary neighbours; select the myosin carrying the
    # largest number of material heads on that actin.
    partners: dict[int, list[tuple[int, int]]] = {}
    for j in range(len(mc)):
        actins, counts = np.unique(bound[j][bound[j] >= 0], return_counts=True)
        for i, count in zip(actins.astype(int), counts.astype(int)):
            partners.setdefault(i, []).append((count, j))

    productive_am = []
    barbed_distance = []
    half_m = 0.5 * myosin_length
    bare_half = 0.05 * myosin_length
    for i in np.where(cb == 2)[0]:
        if i not in partners:
            continue
        _, j = max(partners[i])
        displacement = wrap(ac[i] - mc[j], box, periodic)
        center_projection = displacement.dot(mu[j])
        projected_half_length = 0.5 * actin_length * abs(au[i].dot(mu[j]))
        amin = center_projection - projected_half_length
        amax = center_projection + projected_half_length
        if au[i].dot(mu[j]) > 0:
            head_min, head_max = -half_m, -bare_half
        else:
            head_min, head_max = bare_half, half_m
        overlap = min(amax, head_max) - max(amin, head_min)
        productive_am.append(max(0.0, overlap))
        barbed = displacement - 0.5 * actin_length * au[i]
        barbed_distance.append(abs(barbed.dot(mu[j])))
    return np.asarray(aa), np.asarray(productive_am), np.asarray(barbed_distance)


def make_overlap_figure() -> None:
    aa_all, am_all, barbed_all = [], [], []
    for path in PRODUCTION.values():
        aa, am, barbed = final_overlap_samples(path)
        aa_all.append(aa)
        am_all.append(am)
        barbed_all.append(barbed)
    aa = np.concatenate(aa_all)
    am = np.concatenate(am_all)
    barbed = np.concatenate(barbed_all)
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.1))
    axes[0].hist(aa, bins=np.linspace(0, 0.25, 26), color=GREEN, alpha=0.85, edgecolor="white")
    axes[0].axvline(0.10, color=RED, lw=2, ls="--", label=r"preferred A--A overlap = 0.10 $\mu$m")
    axes[0].axvline(np.median(aa), color=INK, lw=1.5, label=rf"median = {np.median(aa):.3f} $\mu$m")
    axes[0].set_title("Mature A–A overlap", loc="left")
    axes[0].set_xlabel(r"geometric overlap ($\mu$m)")
    axes[0].set_ylabel("status-2 links")
    axes[0].legend()
    axes[1].hist(am, bins=np.linspace(0, 0.70, 29), color=PURPLE, alpha=0.85, edgecolor="white")
    axes[1].axvline(0.675, color=RED, lw=2, ls="--", label=r"full head-bearing half = 0.675 $\mu$m")
    axes[1].axvline(np.median(am), color=INK, lw=1.5, label=rf"median = {np.median(am):.3f} $\mu$m")
    axes[1].set_title("Productive A–M overlap", loc="left")
    axes[1].set_xlabel(r"primary head-bearing half ($\mu$m)")
    axes[1].set_ylabel("status-2 actins")
    axes[1].legend()
    axes[2].hist(barbed, bins=np.linspace(0, 1.8, 31), color=ORANGE, alpha=0.85, edgecolor="white")
    axes[2].axvline(0.75, color=RED, lw=2, ls="--", label=r"thick-filament tip = 0.75 $\mu$m")
    inside = np.mean(barbed < 0.75)
    axes[2].set_title("A–M placement remains imperfect", loc="left")
    axes[2].set_xlabel(r"barbed-end distance from M-line ($\mu$m)")
    axes[2].set_ylabel("status-2 actins")
    axes[2].legend(title=f"{inside:.0%} barbed ends inside thick filament")
    for ax in axes:
        ax.grid(alpha=0.18, lw=0.7)
    fig.suptitle("Final overlap diagnostics at 100 ms (two production seeds pooled)", x=0.06, ha="left", fontsize=13)
    fig.text(0.06, 0.015, "A–M uses each mature actin's primary material-head partner. The old run reached the A–A target but did not fully solve sarcomeric A–M placement.", color=MUTED)
    fig.tight_layout(rect=(0, 0.05, 1, 0.92))
    save(fig, "02-overlap-distributions")


def head_ab_series(path: Path) -> dict[str, np.ndarray]:
    with h5py.File(path, "r", locking=False) as h:
        dt = float(h["meta"].attrs["dt"])
        steps = np.asarray(h["state/current_step"][:, 0], dtype=np.int64)
        events = np.asarray(h["catch_bond/local_link_events"], dtype=float)
        strong = np.asarray(h["actin/am_strong_count"][:, :, 0], dtype=float)
    return {
        "ms": steps * dt * 1e3,
        "active": active_link_series(events, steps),
        "strong_fraction": np.mean(strong > 0, axis=1),
    }


def make_detachment_ab_figure() -> None:
    series = {label: head_ab_series(path) for label, path in HEAD_AB.items()}
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.9))
    colors = [BLUE, RED]
    for color, (label, data) in zip(colors, series.items()):
        mask = data["ms"] <= 20.01
        axes[0].plot(data["ms"][mask], data["active"][mask], color=color, lw=2, label=label)
        axes[1].plot(data["ms"][mask], data["strong_fraction"][mask], color=color, lw=2)
    axes[0].set_title("ACTN2 retention", loc="left")
    axes[0].set_ylabel("active links")
    axes[0].legend(fontsize=8)
    axes[1].set_title("Motor-head accumulation", loc="left")
    axes[1].set_ylabel("actins with ≥1 strong head")
    axes[1].set_ylim(0, 1)
    with h5py.File(CURRENT_FAILURE, "r", locking=False) as h:
        dt = float(h["meta"].attrs["dt"])
        events = np.asarray(h["catch_bond/local_link_events_v3"], dtype=float)
    det = events[(events[:, 4] == 2) & (events[:, 3] * dt * 1e3 <= 36.0)]
    geometric = int(np.sum(det[:, 33] > 0.5))
    stochastic = int(np.sum(det[:, 34] > 0.5))
    axes[2].bar(["geometric\ncontour", "stochastic"], [geometric, stochastic], color=[RED, BLUE], width=0.65)
    axes[2].set_title("Why current links disappeared", loc="left")
    axes[2].set_ylabel("detachments by 36 ms")
    for index, value in enumerate([geometric, stochastic]):
        axes[2].text(index, value + max(geometric, stochastic) * 0.025, str(value), ha="center")
    for ax in axes:
        ax.set_xlabel("time (ms)" if ax is not axes[2] else "event class")
        ax.grid(axis="y", alpha=0.18, lw=0.7)
    fig.suptitle("Resolved regression: transverse force without transverse strain-sensitive escape", x=0.07, ha="left", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.91))
    save(fig, "03-head-detachment-regression")


def banding_noise_multiple(x: np.ndarray, length: float, maximum_mode: int = 30) -> float:
    modes = np.arange(1, maximum_mode + 1)
    amplitudes = np.abs(np.exp(2j * np.pi * np.outer(x / length, modes)).mean(axis=0))
    return float(amplitudes.max() * np.sqrt(len(x)))


def myomesin_metrics(path: Path) -> dict[str, float]:
    with h5py.File(path, "r", locking=False) as h:
        meta = h["meta"].attrs
        length = float(meta["box_x"])
        cb = np.asarray(h["actin/cb_status"][-1, :, 0], dtype=int)
        actin_direction = np.asarray(h["actin/direction"][-1], dtype=float)
        myosin_center = np.asarray(h["myosin/center"][-1], dtype=float)
        count = int(h["catch_bond/local_link_frame_offsets"][-1, 2])
    return {
        "status2": float(np.sum(cb == 2)),
        "links": float(count),
        "banding": banding_noise_multiple(myosin_center[:, 0], length),
        "nematic": float(nematic_series(actin_direction[None, ...])[0]),
    }


def make_myomesin_figure() -> list[dict[str, float]]:
    values = [(k, myomesin_metrics(path)) for k, path in MYOMESIN.items()]
    labels = [str(k) for k, _ in values]
    x = np.arange(len(values))
    fig, axes = plt.subplots(1, 3, figsize=(12.2, 4.0))
    status2 = np.asarray([v["status2"] for _, v in values])
    links = np.asarray([v["links"] for _, v in values])
    banding = np.asarray([v["banding"] for _, v in values])
    nematic = np.asarray([v["nematic"] for _, v in values])
    axes[0].bar(x - 0.18, links, width=0.36, color=BLUE, label="all ACTN2 links")
    axes[0].bar(x + 0.18, status2, width=0.36, color=GREEN, label="status-2 actins")
    axes[0].set_title("More encounters, less maturation", loc="left")
    axes[0].set_ylabel("count at 100 ms")
    axes[0].legend(fontsize=8)
    axes[1].plot(x, banding, marker="o", color=PURPLE, lw=2)
    axes[1].axhline(1, color=MUTED, lw=1, ls="--")
    axes[1].set_title("Myosin centres form M-bands", loc="left")
    axes[1].set_ylabel("peak Fourier amplitude / noise")
    axes[2].plot(x, nematic, marker="o", color=ORANGE, lw=2)
    axes[2].set_ylim(0, 1)
    axes[2].set_title("Strong myomesin disrupts actin order", loc="left")
    axes[2].set_ylabel(r"actin nematic $S$")
    for ax in axes:
        ax.set_xticks(x, labels)
        ax.set_xlabel(r"$k_{\rm mm}$")
        ax.grid(axis="y", alpha=0.18, lw=0.7)
    fig.suptitle("Myomesin scan: M-band formation trades against mature ACTN2", x=0.07, ha="left", fontsize=13)
    fig.text(0.07, 0.01, "640 actins, one seed. kmm = 5 preserves maturation but produces only modest M-band order; ≥15 strongly bands myosin while suppressing status 2.", color=MUTED)
    fig.tight_layout(rect=(0, 0.04, 1, 0.91))
    save(fig, "04-myomesin-tradeoff")
    return [{"k_mm": float(k), **v} for k, v in values]


def write_metrics(production: list[dict[str, float]], myomesin: list[dict[str, float]]) -> None:
    with (OUT / "production_summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(production[0]))
        writer.writeheader()
        writer.writerows(production)
    with (OUT / "myomesin_summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(myomesin[0]))
        writer.writeheader()
        writer.writerows(myomesin)


def main() -> None:
    setup_style()
    production = make_scaleup_figure()
    make_overlap_figure()
    make_detachment_ab_figure()
    myomesin = make_myomesin_figure()
    write_metrics(production, myomesin)


if __name__ == "__main__":
    main()
