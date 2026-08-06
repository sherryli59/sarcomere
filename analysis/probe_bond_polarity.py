#!/usr/bin/env python3
"""Probe 1: relative polarity of LOAD-BEARING (status-2) actin-actin bonds.

The directional kinetics gate a status-2 catch bond on cos_angle < 0 (the two
actin filaments must be ANTIPARALLEL); the non-directional kinetics drop that
gate. If the "indiscriminate glue" explanation for why nondir fails to assemble
is right, nondir should spend a large fraction of its long-lived, load-bearing
bonds on PARALLEL (cos > 0) pairs, whereas directional is all cos < 0 by
construction (which also serves as a positive control on the field reads).

cos_angle = dir_i . dir_j from /actin/direction; bonds from /actin/bonds (which
is status-2-filtered at write time, so every listed pair is load-bearing).
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


def bond_cos(path, nframes=20, minstep=290000):
    with h5py.File(path, "r", locking=False) as h:
        if stepof(h) < minstep:
            return None
        B, D = h["/actin/bonds"], h["/actin/direction"]
        na = D.shape[1]
        n = min(B.shape[0], D.shape[0])
        coss = []
        for fr in range(max(0, n - nframes), n):
            b = np.asarray(B[fr]).reshape(-1, 2)
            dirs = np.asarray(D[fr])
            i, j = b[:, 0].astype(int), b[:, 1].astype(int)
            m = (i >= 0) & (j >= 0) & (i < na) & (j < na) & (i != j)
            i, j = i[m], j[m]
            if len(i) == 0:
                continue
            di, dj = dirs[i], dirs[j]
            c = np.sum(di * dj, axis=1) / (
                np.linalg.norm(di, axis=1) * np.linalg.norm(dj, axis=1) + 1e-12)
            coss.append(c)
        return np.concatenate(coss) if coss else None


def main():
    for tag, dirn in CAMPAIGNS:
        by = defaultdict(list)
        for p in sorted(glob.glob(f"{SCRATCH}/{dirn}/*.h5")):
            c = bond_cos(p)
            if c is not None and len(c):
                by[lt_of(os.path.basename(p))].append(c)
        print(f"\n=== {tag}: relative polarity of status-2 (load-bearing) bonds ===")
        print(f"  {'lambda':>8} {'n_bonds':>10} {'mean_cos':>9} "
              f"{'frac_parallel':>14} {'frac_anti':>10}")
        for lt in sorted(by, key=ltf):
            c = np.concatenate(by[lt])
            print(f"  {ltf(lt):>8.3g} {len(c):>10d} {c.mean():>9.3f} "
                  f"{np.mean(c > 0):>14.3f} {np.mean(c < 0):>10.3f}")


if __name__ == "__main__":
    main()
