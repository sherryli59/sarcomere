"""
Transverse coherence: smectic (dephased columns) vs crystal (registered Z-sheets),
measured WITHOUT projecting out y,z. See reports/transverse_coherence_spec.md.

Headline number (replaces the amplitude-blind smectic_fraction):

  R = A_coh / A_col  in [0,1]
      A_coh = |<exp(-i q* x)>_all|              (transverse-coherent amplitude)
      A_col = sum_c n_c |psi_c| / N             (mean per-column band amplitude)
      psi_c = <exp(-i q* x)>_{filaments in transverse cell c}

  R -> 1 : all columns share phase -> registered crystal (Z-sheets)
  R -> 0 : columns banded but phases random -> dephased smectic (bands cancel)

Also: per-column phasors psi_c (visualize alignment), R_sigma on the crosslink
tension field, and R(t). q* taken from the cb2 longitudinal winning mode.

NOTE: data0430 box is 12x2x2 -- thin transversely. R is a global ratio and is
reliable here; the transverse correlation length needs a wider box (e.g. 12x8x8).
"""
import os
import sys
import argparse
import numpy as np
import h5py

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from coarse_grain_fields import (  # noqa: E402
    parse_box, parse_lt, parse_actin_count, n_frames,
    qdom_timeseries, periodic_midpoint,
)


def steady_qstar(ts, Lx):
    """Steady longitudinal wavevector from the second-half median winning nx*."""
    nx = ts["nx_star"]
    tail = nx[len(nx) // 2:]
    tail = tail[tail > 0]
    if tail.size == 0:
        return 0.0, 0
    nx_star = int(np.median(tail))
    return 2.0 * np.pi * nx_star / Lx, nx_star


def transverse_cell_index(y, z, Ly, Lz, nyb, nzb):
    iy = np.clip((np.mod(y, Ly) / Ly * nyb).astype(int), 0, nyb - 1)
    iz = np.clip((np.mod(z, Lz) / Lz * nzb).astype(int), 0, nzb - 1)
    return iy * nzb + iz


def registration_from_points(x, y, z, qstar, Ly, Lz, nyb, nzb):
    """R, A_coh, A_col, and the per-cell phasor map psi_c for one frame/population."""
    N = x.shape[0]
    if N < 2 or qstar <= 0:
        return np.nan, np.nan, np.nan, np.full(nyb * nzb, np.nan, dtype=complex), \
               np.zeros(nyb * nzb, dtype=int)
    ph = np.exp(-1j * qstar * x)
    A_coh = np.abs(ph.mean())

    cell = transverse_cell_index(y, z, Ly, Lz, nyb, nzb)
    ncell = nyb * nzb
    psi_sum = np.zeros(ncell, dtype=complex)
    counts = np.zeros(ncell, dtype=int)
    np.add.at(psi_sum, cell, ph)
    np.add.at(counts, cell, 1)
    occupied = counts > 0
    psi_c = np.full(ncell, np.nan, dtype=complex)
    psi_c[occupied] = psi_sum[occupied] / counts[occupied]
    # per-column band amplitude, weighted by occupancy
    A_col = np.sum(counts[occupied] * np.abs(psi_c[occupied])) / N
    R = A_coh / A_col if A_col > 0 else np.nan
    return R, A_coh, A_col, psi_c, counts


def process_file(h5path, out_npz, stride=5, nyb=4, nzb=4, population="cb2"):
    name = os.path.basename(h5path)
    Lx, Ly, Lz = parse_box(name)

    with h5py.File(h5path, "r") as fh:
        nf = n_frames(fh)
        frame_idx = list(range(0, nf, stride))
        ts = qdom_timeseries(fh, Lx, Ly, Lz, frame_idx, population=population)
        qstar, nx_star = steady_qstar(ts, Lx)

        # bulk read once
        centers = np.asarray(fh["/actin/center"][:], dtype=float)[np.asarray(frame_idx)]
        cb = np.asarray(fh["/actin/cb_status"][:])[np.asarray(frame_idx)].reshape(len(frame_idx), -1)

        R_t = np.full(len(frame_idx), np.nan)
        Acoh_t = np.full(len(frame_idx), np.nan)
        Acol_t = np.full(len(frame_idx), np.nan)
        psi_last = None
        counts_last = None
        for j in range(len(frame_idx)):
            mask = cb[j] == 2 if population == "cb2" else np.ones(centers.shape[1], bool)
            pts = centers[j][mask]
            R, ac, al, psi_c, counts = registration_from_points(
                pts[:, 0], pts[:, 1], pts[:, 2], qstar, Ly, Lz, nyb, nzb)
            R_t[j], Acoh_t[j], Acol_t[j] = R, ac, al
            psi_last, counts_last = psi_c, counts

        # sigma (crosslink tension) registration: aa bonds, phase at bond x-midpoint,
        # binned transversely at the bond's (y,z) midpoint
        R_sigma_t = np.full(len(frame_idx), np.nan)
        bonds_ds = fh["/actin/bonds"]
        ax_ds = fh["/actin/center"]
        for jj, idx in enumerate(frame_idx):
            bonds = np.asarray(bonds_ds[idx]).reshape(-1, 2)
            valid = (bonds[:, 0] >= 0) & (bonds[:, 1] >= 0)
            if valid.sum() < 2 or qstar <= 0:
                continue
            ac_centers = np.asarray(ax_ds[idx], dtype=float)
            bi = bonds[valid, 0].astype(int)
            bj = bonds[valid, 1].astype(int)
            xm = periodic_midpoint(ac_centers[bi, 0], ac_centers[bj, 0], Lx)
            ym = periodic_midpoint(ac_centers[bi, 1], ac_centers[bj, 1], Ly)
            zm = periodic_midpoint(ac_centers[bi, 2], ac_centers[bj, 2], Lz)
            Rs, *_ = registration_from_points(xm, ym, zm, qstar, Ly, Lz, nyb, nzb)
            R_sigma_t[jj] = Rs

    def tail_mean(a):
        a = a[len(a) // 2:]
        a = a[np.isfinite(a)]
        return float(np.mean(a)) if a.size else np.nan

    payload = {
        "fname": name, "Lx": Lx, "Ly": Ly, "Lz": Lz,
        "lt": parse_lt(name), "n_actins": parse_actin_count(name),
        "nyb": nyb, "nzb": nzb, "qstar": qstar, "nx_star": nx_star,
        "frame": np.asarray(frame_idx),
        "R_t": R_t, "A_coh_t": Acoh_t, "A_col_t": Acol_t, "R_sigma_t": R_sigma_t,
        "R_steady": tail_mean(R_t), "R_sigma_steady": tail_mean(R_sigma_t),
        "psi_last": psi_last if psi_last is not None else np.array([]),
        "counts_last": counts_last if counts_last is not None else np.array([]),
    }
    np.savez_compressed(out_npz, **payload)

    print(f"{name}")
    print(f"  nx*={nx_star}  q*={qstar:.3f}  ({nyb}x{nzb} transverse grid)")
    print(f"  R_steady       = {payload['R_steady']:.3f}"
          f"   (->1 registered crystal / ->0 dephased smectic)")
    print(f"  R_sigma_steady = {payload['R_sigma_steady']:.3f}   (tension-field registration)")
    print(f"  saved {out_npz}")
    return payload


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("h5", nargs="+")
    ap.add_argument("--out_dir", default=".")
    ap.add_argument("--stride", type=int, default=5)
    ap.add_argument("--nyb", type=int, default=4, help="transverse y bins")
    ap.add_argument("--nzb", type=int, default=4, help="transverse z bins")
    ap.add_argument("--population", default="cb2", choices=("cb2", "all"))
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    for path in args.h5:
        stem = os.path.splitext(os.path.basename(path))[0]
        out = os.path.join(args.out_dir, f"tc_{stem}.npz")
        try:
            process_file(path, out, stride=args.stride, nyb=args.nyb,
                         nzb=args.nzb, population=args.population)
        except Exception as e:
            print(f"{os.path.basename(path)}  ERROR  {e}")
