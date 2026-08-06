import argparse
import csv
import os
import re
from collections import defaultdict
from pathlib import Path

os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

import h5py
import numpy as np


DEFAULT_BOX = (12.0, 2.0, 2.0)
DEFAULT_GRID = (120, 20, 20)


def minimal_image(dx, box):
    return dx - box * np.round(dx / box)


def parse_lifetime_coeff(lt_str):
    if lt_str == "0":
        return 0.0
    try:
        return float(lt_str)
    except ValueError:
        base, exp = lt_str.split("e")
        return float(base) * (10 ** float(exp))


def parse_params_from_name(path):
    match = re.search(r"A(\d+)_M\d+_lt([^_]+)", Path(path).name)
    if not match:
        return None, None
    return int(match.group(1)), match.group(2)


def myosin_center_dataset(fh):
    if "myosin/center" in fh:
        return fh["myosin/center"]
    if "myosin/centers" in fh:
        return fh["myosin/centers"]
    raise KeyError("Expected myosin/center or myosin/centers")


def analyze_tilt_invariant_order(x, y, z, box_sizes, grid_sizes, sigma_xy, global_k_star=None):
    Lx, Ly, Lz = box_sizes
    Nx, Ny, Nz = grid_sizes

    X = np.linspace(-Lx / 2, Lx / 2, Nx, endpoint=False)
    Y = np.linspace(-Ly / 2, Ly / 2, Ny, endpoint=False)
    Z = np.linspace(-Lz / 2, Lz / 2, Nz, endpoint=False)

    rho = np.zeros((Nx, Ny, Nz), dtype=float)
    box = np.asarray(box_sizes, dtype=float)
    for j in range(len(x)):
        dx = minimal_image(X - x[j], box[0])
        dy = minimal_image(Y - y[j], box[1])
        dz = minimal_image(Z - z[j], box[2])

        gx = np.exp(-(dx**2) / (2 * sigma_xy**2))
        gy = np.exp(-(dy**2) / (2 * sigma_xy**2))
        gz = np.exp(-(dz**2) / (2 * sigma_xy**2))
        rho += gx[:, None, None] * gy[None, :, None] * gz[None, None, :]

    rho_sum = rho.sum()
    if rho_sum > 0:
        rho /= rho_sum

    rho_centered = rho - rho.mean(axis=0, keepdims=True)
    F_1d = np.fft.rfft(rho_centered, axis=0)

    if global_k_star is None:
        rho_1d = rho.sum(axis=(1, 2))
        rho_1d -= rho_1d.mean()
        F_global = np.fft.rfft(rho_1d)
        global_k_star = int(np.argmax(np.abs(F_global[1:])) + 1) if len(F_global) > 1 else 0

    F_yz = F_1d[global_k_star, :, :]
    F_transverse = np.fft.fft2(F_yz)
    amplitude = np.abs(F_transverse)
    max_order = float(np.max(amplitude))

    ky_idx, kz_idx = np.unravel_index(int(np.argmax(amplitude)), amplitude.shape)
    ky_freq = float(np.fft.fftfreq(Ny)[ky_idx])
    kz_freq = float(np.fft.fftfreq(Nz)[kz_idx])
    tilt_mag = float(np.sqrt(ky_freq**2 + kz_freq**2))

    return max_order, ky_freq, kz_freq, tilt_mag, int(global_k_star)


def convergence_summary(values, tail_count=None, rel_tol=0.05, abs_tol=0.01, min_frames=4):
    arr = np.asarray(values, dtype=float)
    finite = arr[np.isfinite(arr)]
    if finite.size < min_frames:
        return False, np.nan, np.nan, np.nan, np.nan, finite.size

    if tail_count is None:
        tail_count = max(finite.size // 2, 1)

    tail = finite[-min(tail_count, finite.size):]
    tail_mean = float(np.mean(tail))
    tail_std = float(np.std(tail, ddof=1)) if tail.size > 1 else 0.0

    previous = finite[-2 * tail.size : -tail.size]
    if previous.size == 0:
        previous = finite[:tail.size]
    previous_mean = float(np.mean(previous))
    tail_increase = tail_mean - previous_mean

    if tail.size > 1:
        tail_slope = float(np.polyfit(np.arange(tail.size, dtype=float), tail, 1)[0])
    else:
        tail_slope = 0.0

    allowed = max(abs_tol, rel_tol * max(abs(tail_mean), 1e-12))
    per_sample_allowed = allowed / max(tail.size - 1, 1)

    converged = bool(
        tail_increase <= allowed
        and tail_slope <= per_sample_allowed
        and tail_std <= 2.0 * allowed
    )
    return converged, tail_mean, tail_std, tail_increase, tail_slope, finite.size


def analyze_file(
    h5_path,
    last_frames=500,
    stride=50,
    box_sizes=DEFAULT_BOX,
    grid_sizes=DEFAULT_GRID,
    sigma_xy=0.1,
    convergence_tail_count=None,
    convergence_rel_tol=0.05,
    convergence_abs_tol=0.01,
):
    h5_path = Path(h5_path)
    A_val, lt_str = parse_params_from_name(h5_path)
    records = []

    with h5py.File(h5_path, "r") as fh:
        centers_ds = myosin_center_dataset(fh)
        n_frames = int(centers_ds.shape[0])
        start = max(0, n_frames - last_frames)
        frame_indices = list(range(start, n_frames, stride))

        for frame in frame_indices:
            centers = np.asarray(centers_ds[frame], dtype=float)
            order, ky, kz, tilt_mag, k_star = analyze_tilt_invariant_order(
                centers[:, 0],
                centers[:, 1],
                centers[:, 2],
                box_sizes=box_sizes,
                grid_sizes=grid_sizes,
                sigma_xy=sigma_xy,
                global_k_star=None,
            )
            records.append(
                {
                    "file": h5_path.name,
                    "path": str(h5_path),
                    "A_val": A_val,
                    "lt_str": lt_str,
                    "frame": frame,
                    "n_frames": n_frames,
                    "k_star": k_star,
                    "tilt_order": order,
                    "ky_freq": ky,
                    "kz_freq": kz,
                    "tilt_magnitude": tilt_mag,
                }
            )

    orders = [row["tilt_order"] for row in records]
    tilts = [row["tilt_magnitude"] for row in records]
    converged, mean_order, std_order, order_tail_increase, order_tail_slope, n_eval = convergence_summary(
        orders,
        tail_count=convergence_tail_count,
        rel_tol=convergence_rel_tol,
        abs_tol=convergence_abs_tol,
    )
    _, mean_tilt, std_tilt, tilt_tail_increase, tilt_tail_slope, _ = convergence_summary(
        tilts,
        tail_count=convergence_tail_count,
        rel_tol=convergence_rel_tol,
        abs_tol=convergence_abs_tol,
    )

    summary = {
        "file": h5_path.name,
        "path": str(h5_path),
        "A_val": A_val,
        "lt_str": lt_str,
        "n_frames": records[0]["n_frames"] if records else 0,
        "n_eval_frames": n_eval,
        "first_eval_frame": records[0]["frame"] if records else "",
        "last_eval_frame": records[-1]["frame"] if records else "",
        "order_mean_tail": mean_order,
        "order_std_tail": std_order,
        "order_tail_increase": order_tail_increase,
        "order_tail_slope": order_tail_slope,
        "order_converged": converged,
        "tilt_mean_tail": mean_tilt,
        "tilt_std_tail": std_tilt,
        "tilt_tail_increase": tilt_tail_increase,
        "tilt_tail_slope": tilt_tail_slope,
        "converged": bool(converged),
    }
    return records, summary


def discover_h5_files(data_dir):
    return sorted(Path(data_dir).glob("*.h5"))


def write_csv(path, rows, fieldnames):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def build_phase_grid(summaries):
    grouped = defaultdict(list)
    for row in summaries:
        if row["A_val"] is None or row["lt_str"] is None:
            continue
        grouped[(row["A_val"], row["lt_str"])].append(row)

    A_vals = sorted({key[0] for key in grouped})
    lt_vals = sorted({key[1] for key in grouped}, key=parse_lifetime_coeff)
    order_mean = np.full((len(A_vals), len(lt_vals)), np.nan)
    order_std = np.full_like(order_mean, np.nan)
    tilt_mean = np.full_like(order_mean, np.nan)
    tilt_std = np.full_like(order_mean, np.nan)
    convergence_fraction = np.full_like(order_mean, np.nan)
    n = np.zeros_like(order_mean, dtype=int)

    for i, A_val in enumerate(A_vals):
        for j, lt_str in enumerate(lt_vals):
            rows = grouped.get((A_val, lt_str), [])
            if not rows:
                continue
            orders = np.asarray([r["order_mean_tail"] for r in rows], dtype=float)
            tilts = np.asarray([r["tilt_mean_tail"] for r in rows], dtype=float)
            conv = np.asarray([r["converged"] for r in rows], dtype=float)
            ddof = 1 if len(rows) > 1 else 0
            order_mean[i, j] = np.nanmean(orders)
            order_std[i, j] = np.nanstd(orders, ddof=ddof)
            tilt_mean[i, j] = np.nanmean(tilts)
            tilt_std[i, j] = np.nanstd(tilts, ddof=ddof)
            convergence_fraction[i, j] = np.nanmean(conv)
            n[i, j] = len(rows)

    return A_vals, lt_vals, order_mean, order_std, tilt_mean, tilt_std, convergence_fraction, n


def main():
    parser = argparse.ArgumentParser(
        description="Compute tilt-invariant structural order and convergence for sarcomere trajectories."
    )
    parser.add_argument("--data-dir", default="/scratch/users/xsli/sarcomere/data0416")
    parser.add_argument("--output-dir", default="tilt_invariant_order_results")
    parser.add_argument("--last-frames", type=int, default=500)
    parser.add_argument("--stride", type=int, default=50)
    parser.add_argument("--sigma-xy", type=float, default=0.1)
    parser.add_argument("--grid", type=int, nargs=3, default=DEFAULT_GRID, metavar=("NX", "NY", "NZ"))
    parser.add_argument("--box", type=float, nargs=3, default=DEFAULT_BOX, metavar=("LX", "LY", "LZ"))
    parser.add_argument(
        "--tail-count",
        type=int,
        default=None,
        help="Number of sampled points in the final convergence window. Default: half of sampled points.",
    )
    parser.add_argument("--rel-tol", type=float, default=0.05)
    parser.add_argument("--abs-tol", type=float, default=0.01)
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    files = discover_h5_files(args.data_dir)
    if not files:
        raise FileNotFoundError(f"No .h5 files found in {args.data_dir}")

    all_records = []
    summaries = []
    for idx, h5_path in enumerate(files, start=1):
        print(f"[{idx}/{len(files)}] {h5_path.name}")
        try:
            records, summary = analyze_file(
                h5_path,
                last_frames=args.last_frames,
                stride=args.stride,
                box_sizes=tuple(args.box),
                grid_sizes=tuple(args.grid),
                sigma_xy=args.sigma_xy,
                convergence_tail_count=args.tail_count,
                convergence_rel_tol=args.rel_tol,
                convergence_abs_tol=args.abs_tol,
            )
            all_records.extend(records)
            summaries.append(summary)
            status = "CONVERGED" if summary["converged"] else "NOT CONVERGED"
            print(
                f"  {status}: order={summary['order_mean_tail']:.6g} "
                f"+/- {summary['order_std_tail']:.3g}, "
                f"tail_increase={summary['order_tail_increase']:.3g}, "
                f"tail_slope={summary['order_tail_slope']:.3g}, "
                f"tilt={summary['tilt_mean_tail']:.6g} +/- {summary['tilt_std_tail']:.3g}"
            )
        except Exception as exc:
            print(f"  ERROR: {exc}")
            summaries.append(
                {
                    "file": h5_path.name,
                    "path": str(h5_path),
                    "A_val": None,
                    "lt_str": None,
                    "n_frames": 0,
                    "n_eval_frames": 0,
                    "first_eval_frame": "",
                    "last_eval_frame": "",
                    "order_mean_tail": np.nan,
                    "order_std_tail": np.nan,
                    "order_tail_increase": np.nan,
                    "order_tail_slope": np.nan,
                    "order_converged": False,
                    "tilt_mean_tail": np.nan,
                    "tilt_std_tail": np.nan,
                    "tilt_tail_increase": np.nan,
                    "tilt_tail_slope": np.nan,
                    "converged": False,
                }
            )

    record_fields = [
        "file",
        "path",
        "A_val",
        "lt_str",
        "frame",
        "n_frames",
        "k_star",
        "tilt_order",
        "ky_freq",
        "kz_freq",
        "tilt_magnitude",
    ]
    summary_fields = [
        "file",
        "path",
        "A_val",
        "lt_str",
        "n_frames",
        "n_eval_frames",
        "first_eval_frame",
        "last_eval_frame",
        "order_mean_tail",
        "order_std_tail",
        "order_tail_increase",
        "order_tail_slope",
        "order_converged",
        "tilt_mean_tail",
        "tilt_std_tail",
        "tilt_tail_increase",
        "tilt_tail_slope",
        "converged",
    ]
    write_csv(output_dir / "tilt_invariant_order_timeseries.csv", all_records, record_fields)
    write_csv(output_dir / "tilt_invariant_order_summary.csv", summaries, summary_fields)

    A_vals, lt_vals, order_mean, order_std, tilt_mean, tilt_std, conv_frac, n = build_phase_grid(summaries)
    np.savez(
        output_dir / "tilt_invariant_order_phase_grid.npz",
        A_vals=np.asarray(A_vals),
        lt_vals=np.asarray(lt_vals, dtype=str),
        order_mean=order_mean,
        order_std=order_std,
        tilt_mean=tilt_mean,
        tilt_std=tilt_std,
        convergence_fraction=conv_frac,
        n=n,
    )
    print(f"Wrote results to {output_dir}")


if __name__ == "__main__":
    main()
