#!/usr/bin/env python3
"""Time-resolved lattice diagnostics for two six-actin myosin half-rings."""
import argparse
from pathlib import Path
import h5py
import matplotlib.pyplot as plt
import numpy as np


def main():
    p = argparse.ArgumentParser()
    p.add_argument("trajectory", type=Path)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    with h5py.File(args.trajectory, "r") as h:
        ac = np.asarray(h["actin/center"], float)
        mc = np.asarray(h["myosin/center"], float)[:, 0]
        mu = np.asarray(h["myosin/direction"], float)[:, 0]
        step = np.asarray(h["state/current_step"], int)[:, 0]
        dt = float(h["meta"].attrs["dt"])
        target = float(h["meta"].attrs["am_preferred_centerline"])
        state = np.asarray(h["state/am_head_state_current"], int)[:, 0]
        owner = np.asarray(
            h["state/am_head_bound_actin_current"], int)[:, 0]
        aa_n = np.asarray(h["state/aa_pairs_current_count"], int)[:, 0]
        aa_s = np.asarray(h["state/aa_status_current"], int)[:, :, 0]
    time = step * dt
    d = ac - mc[:, None]
    axial = np.sum(d * mu[:, None], axis=2)
    transverse = d - axial[:, :, None] * mu[:, None]
    radius = np.linalg.norm(transverse, axis=2)
    # Any orthonormal transverse coordinates give invariant |psi6| and gap CV.
    ref = np.tile(np.array([0., 0., 1.]), (len(mu), 1))
    bad = np.abs(np.sum(ref * mu, axis=1)) > .9
    ref[bad] = np.array([0., 1., 0.])
    e1 = np.cross(mu, ref); e1 /= np.linalg.norm(e1, axis=1)[:, None]
    e2 = np.cross(mu, e1)
    theta = np.arctan2(np.sum(transverse * e2[:, None], axis=2),
                       np.sum(transverse * e1[:, None], axis=2))
    values = []
    for lo in (0, 6):
        r = radius[:, lo:lo + 6]
        th = theta[:, lo:lo + 6]
        order = np.abs(np.mean(np.exp(6j * th), axis=1))
        sort = np.sort(np.mod(th, 2*np.pi), axis=1)
        gaps = np.diff(np.c_[sort, sort[:, 0] + 2*np.pi], axis=1)
        values.append((r.mean(1), r.std(1), order,
                       gaps.std(1) / gaps.mean(1)))
    weak = np.sum(state == 1, axis=tuple(range(1, state.ndim)))
    strong = np.sum(state == 2, axis=tuple(range(1, state.ndim)))
    weak_half = np.stack([
        np.sum((state == 1) & (owner >= lo) & (owner < lo + 6), axis=1)
        for lo in (0, 6)], axis=1)
    strong_half = np.stack([
        np.sum((state == 2) & (owner >= lo) & (owner < lo + 6), axis=1)
        for lo in (0, 6)], axis=1)
    axial_mean = np.stack([
        np.mean(axial[:, lo:lo + 6], axis=1) for lo in (0, 6)], axis=1)
    s1 = np.array([np.count_nonzero(aa_s[i, :n] == 1)
                   for i, n in enumerate(aa_n)])
    s2 = np.array([np.count_nonzero(aa_s[i, :n] == 2)
                   for i, n in enumerate(aa_n)])
    args.output.mkdir(parents=True, exist_ok=True)
    columns = [time]
    names = ["time_s"]
    for half, (mean, std, psi, gap) in enumerate(values):
        columns += [mean, std, psi, gap]
        names += [f"h{half}_am_mean_um", f"h{half}_am_std_um",
                  f"h{half}_psi6", f"h{half}_gap_cv"]
    columns += [weak, strong, s1, s2]
    names += ["weak_heads", "strong_heads", "aa_status1", "aa_status2"]
    for half in range(2):
        columns += [weak_half[:, half], strong_half[:, half],
                    axial_mean[:, half]]
        names += [f"h{half}_weak_heads", f"h{half}_strong_heads",
                  f"h{half}_axial_mean_um"]
    np.savetxt(args.output/"timeseries.csv", np.c_[*columns], delimiter=",",
               header=",".join(names), comments="")
    fig, ax = plt.subplots(5, 1, figsize=(9, 12), sharex=True)
    colors = ("#0072B2", "#D55E00")
    for half, (mean, std, psi, gap) in enumerate(values):
        ax[0].plot(time, 1000*mean, color=colors[half],
                   label=f"half {half}")
        ax[0].fill_between(time, 1000*(mean-std), 1000*(mean+std),
                           color=colors[half], alpha=.2)
        ax[1].plot(time, psi, color=colors[half], label=f"half {half}")
        ax[2].plot(time, gap, color=colors[half], label=f"half {half}")
    ax[0].axhline(1000*target, color="k", ls="--", label=r"$d^*_{AM}$")
    ax[0].set_ylabel("AM spacing (nm)\nmean ± SD"); ax[0].legend(ncol=3)
    ax[1].set_ylabel(r"$|\psi_6|$"); ax[1].set_ylim(0, 1.03)
    ax[2].set_ylabel("angular-gap CV"); ax[2].set_ylim(bottom=0)
    for half in range(2):
        ax[3].plot(time, 1000*axial_mean[:, half], color=colors[half],
                   label=f"half {half}")
    ax[3].axhline(0, color="k", ls="--", lw=1)
    ax[3].set_ylabel("axial position\nrelative to M-line (nm)")
    ax[3].legend(ncol=2)
    ax[4].plot(time, strong_half[:, 0], color=colors[0],
               label="strong AM, half 0")
    ax[4].plot(time, strong_half[:, 1], color=colors[1],
               label="strong AM, half 1")
    ax[4].plot(time, s1, color="#009E73", label="AA status 1")
    ax[4].plot(time, s2, color="#CC79A7", label="AA status 2")
    ax[4].set_ylabel("count"); ax[4].set_xlabel("time (s)")
    ax[4].legend(ncol=2)
    fig.tight_layout(); fig.savefig(args.output/"bilateral_lattice.png", dpi=180)
    tail = slice(int(.75*len(time)), None)
    for half, (mean, std, psi, gap) in enumerate(values):
        print(f"half {half}: AM={1000*np.mean(mean[tail]):.2f} ± "
              f"{1000*np.mean(std[tail]):.2f} nm, "
              f"psi6={np.mean(psi[tail]):.3f}, gapCV={np.mean(gap[tail]):.3f}")
    print(f"tail strong={np.mean(strong[tail]):.2f}, "
          f"AA1={np.mean(s1[tail]):.2f}, AA2={np.mean(s2[tail]):.2f}")

if __name__ == "__main__":
    main()
