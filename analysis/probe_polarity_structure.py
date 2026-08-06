#!/usr/bin/env python3
"""Probe 2: polarity structure factor S_p(q) along the sarcomere (x) axis.

Probe 1 showed nondir load-bearing bonds are majority PARALLEL. If that produces
parallel BUNDLES rather than an antiparallel sarcomeric REGISTER, it should show
up in how filament polarity is arranged in space.

Let p_k = sign(dir_x_k) (filament polarity) and x_k its center. Define

    S_p(q) = | sum_k p_k exp(i q x_k) |^2 / N ,   q = 2*pi*n / Lx

Three regimes are distinguishable by where S_p(q) peaks:
  * peak at FINITE q  -> polarity alternates with a period  = sarcomeric register
                         (expect directional; q should match the myosin phi_3d mode)
  * peak at q -> 0    -> filaments phase-separate into large same-polarity domains
                         = parallel bundling (the nondir prediction)
  * flat / no peak    -> polarity randomly mixed

Reported per (campaign, lambda): net polarization S_p(0)/N, the dominant nonzero
mode n* and its S_p, averaged over late frames of completed runs.
"""
from __future__ import annotations

import glob
import os
import re
from collections import defaultdict

os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

import h5py
import numpy as np

SCRATCH = "/scratch/users/xsli/sarcomere"
CAMPAIGNS = [
    ("directional", "data_a2000_status2_lambda"),
    ("nondir", "data_a2000_status2_nondir_lambda"),
]
LX = 12.0
NMAX = 30


def stepof(h):
    for k in ("/resume/current_step", "/resume_prev/current_step"):
        if k in h:
            return int(np.asarray(h[k]).reshape(-1)[-1])
    return int(np.asarray(h["/state/current_step"])[-1, 0])


def lt_of(name):
    m = re.search(r"_lt([0-9e.\-]+?)_(?:nodir_)?r\d", name)
    return m.group(1) if m else "?"


def ltf(s):
    if s in ("0", "lt0"):
        return 0.0
    try:
        mant, exp = s.split("e")
        return float(mant) * 10 ** float(exp)
    except ValueError:
        return float(s)


def sp_q(path, nframes=20, minstep=290000):
    qs = 2 * np.pi * np.arange(0, NMAX + 1) / LX
    with h5py.File(path, "r", locking=False) as h:
        if stepof(h) < minstep:
            return None
        C, D = h["/actin/center"], h["/actin/direction"]
        n = min(C.shape[0], D.shape[0])
        acc = np.zeros(NMAX + 1)
        cnt = 0
        for fr in range(max(0, n - nframes), n):
            x = np.asarray(C[fr])[:, 0]
            p = np.sign(np.asarray(D[fr])[:, 0])
            N = len(p)
            phase = np.exp(1j * np.outer(qs, x))       # (NMAX+1, N)
            s = np.abs(phase @ p) ** 2 / N
            acc += s
            cnt += 1
    return acc / cnt if cnt else None


def main():
    for tag, dirn in CAMPAIGNS:
        by = defaultdict(list)
        for p in sorted(glob.glob(f"{SCRATCH}/{dirn}/*.h5")):
            s = sp_q(p)
            if s is not None:
                by[lt_of(os.path.basename(p))].append(s)
        print(f"\n=== {tag}: polarity structure factor S_p(q) along x (Lx={LX}) ===")
        print(f"  {'lambda':>8} {'Sp(0)=netpol':>12} {'peak_n*':>8} "
              f"{'q*':>7} {'Sp(n*)':>9} {'period=Lx/n*':>12}")
        for lt in sorted(by, key=ltf):
            s = np.mean(by[lt], axis=0)
            nstar = int(1 + np.argmax(s[1:]))          # dominant nonzero mode
            qstar = 2 * np.pi * nstar / LX
            print(f"  {ltf(lt):>8.3g} {s[0]:>12.1f} {nstar:>8d} "
                  f"{qstar:>7.2f} {s[nstar]:>9.1f} {LX / nstar:>12.3f}")


if __name__ == "__main__":
    main()
