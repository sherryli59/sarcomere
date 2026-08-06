#!/usr/bin/env python3
"""Diagnostics for the 2-myosin, 4-half, 24-actin periodic lattice cell."""
import argparse
from pathlib import Path
import h5py
import matplotlib.pyplot as plt
import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("trajectory", type=Path)
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    with h5py.File(a.trajectory, "r") as h:
        ac = h["actin/center"][:]
        mc = h["myosin/center"][:]
        mu = h["myosin/direction"][:]
        step = h["state/current_step"][:, 0]
        dt = float(h["meta"].attrs["dt"])
        target = float(h["meta"].attrs["am_preferred_centerline"])
        hs = h["state/am_head_state_current"][:]
        ho = h["state/am_head_bound_actin_current"][:]
        aa_n = h["state/aa_pairs_current_count"][:, 0]
        aa_s = h["state/aa_status_current"][:, :, 0]
    time = step * dt
    owners = (0, 1, 0, 1)
    labels = ("M0 −half", "M1 +half", "M0 +half", "M1 −half")
    colors = ("#0072B2", "#D55E00", "#56B4E9", "#E69F00")
    metrics = []
    for ring, owner in enumerate(owners):
        sl = slice(6*ring, 6*ring+6)
        d = ac[:, sl] - mc[:, owner, None]
        axis = mu[:, owner, None]
        axial = np.sum(d * axis, axis=2)
        trans = d - axial[:, :, None] * axis
        radius = np.linalg.norm(trans, axis=2)
        ref = np.tile([0., 0., 1.], (len(time), 1))
        bad = np.abs(np.sum(ref * mu[:, owner], axis=1)) > .9
        ref[bad] = [0., 1., 0.]
        e1 = np.cross(mu[:, owner], ref)
        e1 /= np.linalg.norm(e1, axis=1)[:, None]
        e2 = np.cross(mu[:, owner], e1)
        theta = np.arctan2(np.sum(trans*e2[:, None], axis=2),
                           np.sum(trans*e1[:, None], axis=2))
        order = np.abs(np.mean(np.exp(6j*theta), axis=1))
        sort = np.sort(np.mod(theta, 2*np.pi), axis=1)
        gaps = np.diff(np.c_[sort, sort[:, 0]+2*np.pi], axis=1)
        strong = np.sum((hs == 2) & (ho >= 6*ring) &
                        (ho < 6*ring+6), axis=(1, 2))
        metrics.append((radius.mean(1), radius.std(1), order,
                        gaps.std(1)/gaps.mean(1), axial.mean(1), strong))
    s1 = np.array([np.count_nonzero(aa_s[i, :n] == 1)
                   for i, n in enumerate(aa_n)])
    s2 = np.array([np.count_nonzero(aa_s[i, :n] == 2)
                   for i, n in enumerate(aa_n)])
    a.output.mkdir(parents=True, exist_ok=True)
    cols, names = [time], ["time_s"]
    for ring, m in enumerate(metrics):
        for name, value in zip(("am_mean_um", "am_std_um", "psi6",
                                "gap_cv", "axial_um", "strong"), m):
            cols.append(value); names.append(f"ring{ring}_{name}")
    cols += [s1, s2]; names += ["aa_status1", "aa_status2"]
    np.savetxt(a.output/"timeseries.csv", np.c_[*cols], delimiter=",",
               header=",".join(names), comments="")
    fig, ax = plt.subplots(5, 1, figsize=(10, 13), sharex=True)
    for ring, (mean, std, psi, gap, axial, strong) in enumerate(metrics):
        ax[0].plot(time, 1e3*mean, color=colors[ring], label=labels[ring])
        ax[0].fill_between(time, 1e3*(mean-std), 1e3*(mean+std),
                           color=colors[ring], alpha=.12)
        ax[1].plot(time, psi, color=colors[ring])
        ax[2].plot(time, gap, color=colors[ring])
        ax[3].plot(time, 1e3*axial, color=colors[ring])
        ax[4].plot(time, strong, color=colors[ring], label=labels[ring])
    ax[0].axhline(1e3*target, color="k", ls="--", label=r"$d^*_{AM}$")
    ax[0].set_ylabel("AM spacing (nm)\nmean ± SD"); ax[0].legend(ncol=3)
    ax[1].set_ylabel(r"$|\psi_6|$"); ax[1].set_ylim(0, 1.03)
    ax[2].set_ylabel("angular-gap CV"); ax[2].set_ylim(bottom=0)
    ax[3].set_ylabel("ring axial position\nrelative to M-line (nm)")
    ax[3].axhline(0, color="k", ls="--", lw=1)
    aa = ax[4].twinx()
    aa.plot(time, s1, color="#009E73", ls="--", label="AA status 1")
    aa.plot(time, s2, color="#CC79A7", ls="--", label="AA status 2")
    ax[4].set_ylabel("strong heads"); aa.set_ylabel("ACTN2 links")
    ax[4].set_xlabel("time (s)")
    lines = ax[4].lines + aa.lines
    ax[4].legend(lines, [x.get_label() for x in lines], ncol=3)
    fig.tight_layout(); fig.savefig(a.output/"force_balanced_lattice.png", dpi=180)
    tail = slice(int(.75*len(time)), None)
    for ring, (mean, std, psi, gap, axial, strong) in enumerate(metrics):
        print(f"{labels[ring]}: AM={1e3*np.mean(mean[tail]):.2f}±"
              f"{1e3*np.mean(std[tail]):.2f} nm psi6={np.mean(psi[tail]):.3f} "
              f"gapCV={np.mean(gap[tail]):.3f} "
              f"axial={1e3*np.mean(axial[tail]):.1f} nm "
              f"strong={np.mean(strong[tail]):.1f}")
    print(f"tail AA1={np.mean(s1[tail]):.2f} AA2={np.mean(s2[tail]):.2f}")

if __name__ == "__main__":
    main()
