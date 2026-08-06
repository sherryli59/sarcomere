"""Aggregate per-frame sarcomeric metrics by replica and plot A x lambda maps."""

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.patheffects as path_effects
import matplotlib.pyplot as plt
import numpy as np


DERIVED_METRICS = {
    "shared_coherence": (
        "shared_coherence_yz_myosin",
        "shared_coherence_yz_actin_cb2",
    ),
}

PLOT_METRICS = {
    "q_sarcomere": ("Tilt-invariant sarcomeric score", "magma", 0.0, 1.0),
    "q_sarcomere_axial": ("Axial-only sarcomeric score", "magma", 0.0, 1.0),
    "q_gain": ("Tilt-invariant minus axial score", "coolwarm", None, None),
    "shared_joint_band_amplitude": (
        "Joint corrected band amplitude", "viridis", 0.0, 1.0
    ),
    "shared_coherence": ("Shared-mode transverse coherence", "viridis", 0.0, 1.0),
    "shared_registration_factor": (
        "M/Z antiphase registration factor", "viridis", 0.0, 1.0
    ),
    "phi_3d_myosin": ("Myosin $\\Phi_{3D}$", "viridis", 0.0, 1.0),
    "phi_3d_actin_cb2": ("CB2 actin $\\Phi_{3D}$", "viridis", 0.0, 1.0),
    "shared_theta_3d_deg": ("Shared-mode tilt (degrees)", "plasma", 0.0, None),
}


def lambda_to_float(label):
    label = str(label).strip()
    if "e" in label.lower():
        mantissa, exponent = label.lower().split("e", 1)
        return float(mantissa) * 10.0 ** float(exponent)
    return float(label)


def parse_filename(filename, default_mmb):
    patterns = {
        "A": r"(?:^|_)A(\d+)",
        "M": r"(?:^|_)M(\d+)",
        "lambda_label": r"_lt([^_]+)_",
        "replicate": r"_(r\d+)_",
        "mmb": r"_mmb(\d+)_",
        "box": r"_L([^_]+)\.h5$",
    }
    matches = {key: re.search(pattern, filename) for key, pattern in patterns.items()}
    required = ("A", "M", "lambda_label", "replicate", "box")
    missing = [key for key in required if matches[key] is None]
    if missing:
        raise ValueError(f"Cannot parse {missing} from {filename}")
    box_match = re.fullmatch(
        r"([0-9.]+)x([0-9.]+)x([0-9.]+)", matches["box"].group(1)
    )
    if box_match is None:
        raise ValueError(f"Cannot parse box dimensions from {filename}")
    lambda_label = matches["lambda_label"].group(1)
    return {
        "A": int(matches["A"].group(1)),
        "M": int(matches["M"].group(1)),
        "lambda": lambda_to_float(lambda_label),
        "lambda_label": lambda_label,
        "replicate": matches["replicate"].group(1),
        "mmb": (
            int(matches["mmb"].group(1))
            if matches["mmb"] is not None
            else default_mmb
        ),
        "Lx": float(box_match.group(1)),
        "Ly": float(box_match.group(2)),
        "Lz": float(box_match.group(3)),
    }


def finite_float(row, key):
    value = float(row[key])
    return value if np.isfinite(value) else np.nan


def load_by_replicate(path, default_mmb):
    runs = defaultdict(list)
    with Path(path).open(newline="") as stream:
        for row in csv.DictReader(stream):
            record = dict(row)
            record["frame"] = int(row["frame"])
            for metric in PLOT_METRICS:
                if metric in row:
                    record[metric] = finite_float(row, metric)
            record["shared_coherence"] = float(
                np.sqrt(
                    finite_float(row, DERIVED_METRICS["shared_coherence"][0])
                    * finite_float(row, DERIVED_METRICS["shared_coherence"][1])
                )
            )
            registration_pi = finite_float(row, "shared_registration_pi")
            record["shared_registration_factor"] = 0.5 * (1.0 + registration_pi)
            record["q_gain"] = (
                finite_float(row, "q_sarcomere")
                - finite_float(row, "q_sarcomere_axial")
            )
            runs[row["filename"]].append(record)

    parsed = {}
    for filename, rows in runs.items():
        rows.sort(key=lambda row: row["frame"])
        parsed[filename] = (parse_filename(filename, default_mmb), rows)
    return parsed


def mean_finite(values):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    return float(np.mean(values)) if values.size else np.nan


def summarize_replicates(runs, tail_samples):
    summaries = []
    for filename, (params, rows) in sorted(runs.items()):
        tail = rows[-min(tail_samples, len(rows)):]
        summary = {
            "filename": filename,
            **params,
            "n_tail_samples": len(tail),
            "first_tail_frame": tail[0]["frame"],
            "last_frame": tail[-1]["frame"],
        }
        for metric in PLOT_METRICS:
            summary[metric] = mean_finite([row[metric] for row in tail])
        summaries.append(summary)
    return summaries


def aggregate_parameters(replicates):
    cells = defaultdict(list)
    for record in replicates:
        key = (
            record["A"], record["M"], record["lambda"], record["lambda_label"],
            record["mmb"], record["Lx"], record["Ly"], record["Lz"],
        )
        cells[key].append(record)

    aggregates = []
    for key, records in sorted(cells.items()):
        A, M, lambda_value, lambda_label, mmb, Lx, Ly, Lz = key
        result = {
            "A": A, "M": M, "lambda": lambda_value,
            "lambda_label": lambda_label, "mmb": mmb,
            "Lx": Lx, "Ly": Ly, "Lz": Lz,
            "n_replicates": len(records),
        }
        for metric in PLOT_METRICS:
            values = np.asarray([record[metric] for record in records], dtype=float)
            values = values[np.isfinite(values)]
            result[f"{metric}_mean"] = float(np.mean(values)) if values.size else np.nan
            result[f"{metric}_std"] = (
                float(np.std(values, ddof=1)) if values.size > 1 else 0.0
            )
        aggregates.append(result)
    return aggregates


def write_csv(path, records):
    if not records:
        raise ValueError(f"No records available for {path}")
    with Path(path).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    print(f"Wrote {len(records)} rows to {path}")


def validate_fixed_parameters(records, varying_key):
    fixed_keys = ["A", "M", "mmb", "Lx", "Ly", "Lz"]
    fixed_keys.remove(varying_key)
    for key in fixed_keys:
        values = sorted({record[key] for record in records})
        if len(values) != 1:
            raise ValueError(
                f"Cannot make one {varying_key} x lambda diagram: "
                f"{key} has values {values}"
            )


def phase_grid(records, metric, varying_key):
    x_values = sorted({record[varying_key] for record in records})
    lambda_values = sorted({record["lambda"] for record in records})
    labels = {}
    grid = np.full((len(lambda_values), len(x_values)), np.nan)
    counts = np.zeros_like(grid, dtype=int)
    for record in records:
        i = lambda_values.index(record["lambda"])
        j = x_values.index(record[varying_key])
        grid[i, j] = record[f"{metric}_mean"]
        counts[i, j] = record["n_replicates"]
        labels[record["lambda"]] = record["lambda_label"]
    return x_values, lambda_values, [labels[value] for value in lambda_values], grid, counts


def color_limits(metric, grid, configured_min, configured_max):
    finite = grid[np.isfinite(grid)]
    if metric == "q_gain":
        limit = max(float(np.max(np.abs(finite))) if finite.size else 0.0, 1e-12)
        return -limit, limit
    vmin = configured_min
    vmax = configured_max
    if vmax is None:
        vmax = float(np.max(finite)) if finite.size else 1.0
    return vmin, vmax


def draw_phase_panel(ax, records, metric, varying_key, include_title=True):
    title, cmap_name, configured_min, configured_max = PLOT_METRICS[metric]
    x_values, _, lambda_labels, grid, counts = phase_grid(
        records, metric, varying_key
    )
    masked = np.ma.masked_invalid(grid)
    cmap = plt.get_cmap(cmap_name).copy()
    cmap.set_bad("#eeeeee")
    vmin, vmax = color_limits(metric, grid, configured_min, configured_max)
    image = ax.imshow(
        masked, origin="lower", aspect="auto", interpolation="none",
        cmap=cmap, vmin=vmin, vmax=vmax,
    )
    ax.set_xticks(range(len(x_values)), labels=x_values)
    ax.set_yticks(range(len(lambda_labels)), labels=lambda_labels)
    ax.set_xlabel("Actin count A" if varying_key == "A" else "Myosin count M")
    ax.set_ylabel("Catch-bond parameter $\\lambda$")
    if include_title:
        ax.set_title(title)
    for i in range(grid.shape[0]):
        for j in range(grid.shape[1]):
            if not np.isfinite(grid[i, j]):
                continue
            annotation = ax.text(
                j, i, f"{grid[i, j]:.3g}\nn={counts[i, j]}",
                ha="center", va="center", fontsize=8, color="black",
            )
            annotation.set_path_effects(
                [path_effects.withStroke(linewidth=2, foreground="white")]
            )
    return image, title


def plot_phase_diagrams(records, output_dir, varying_key, context):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    validate_fixed_parameters(records, varying_key)
    for metric in PLOT_METRICS:
        fig, ax = plt.subplots(figsize=(6.5, 5.0))
        image, title = draw_phase_panel(ax, records, metric, varying_key)
        ax.set_title(f"{title}\n{context}")
        fig.colorbar(image, ax=ax, label=title)
        fig.tight_layout()
        for suffix in ("png", "svg"):
            path = output_dir / f"phase_{metric}.{suffix}"
            fig.savefig(path, dpi=200)
            print(f"Wrote {path}")
        plt.close(fig)

    fig, axes = plt.subplots(3, 3, figsize=(16, 13), constrained_layout=True)
    for ax, metric in zip(axes.flat, PLOT_METRICS):
        image, title = draw_phase_panel(ax, records, metric, varying_key)
        fig.colorbar(image, ax=ax, shrink=0.85, label=title)
    fig.suptitle(context, fontsize=18)
    combined = output_dir / "phase_sarcomeric_components.png"
    fig.savefig(combined, dpi=200)
    fig.savefig(combined.with_suffix(".svg"))
    plt.close(fig)
    print(f"Wrote {combined} and {combined.with_suffix('.svg')}")

    metric_rows = list(PLOT_METRICS)
    for row_index in range(3):
        row_metrics = metric_rows[3 * row_index:3 * (row_index + 1)]
        fig, axes = plt.subplots(
            1, 3, figsize=(16, 5.2), constrained_layout=True
        )
        for ax, metric in zip(axes, row_metrics):
            image, title = draw_phase_panel(ax, records, metric, varying_key)
            fig.colorbar(image, ax=ax, shrink=0.9, label=title)
        fig.suptitle(context, fontsize=17)
        row_path = output_dir / f"phase_sarcomeric_components_row{row_index + 1}.png"
        fig.savefig(row_path, dpi=200)
        fig.savefig(row_path.with_suffix(".svg"))
        plt.close(fig)
        print(f"Wrote {row_path} and {row_path.with_suffix('.svg')}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_csv", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--tail-samples", type=int, default=20,
        help="Number of final sampled frames averaged within each replica.",
    )
    parser.add_argument(
        "--default-mmb", type=int, default=6,
        help="max_myosin_bonds used when the filename omits _mmbN.",
    )
    args = parser.parse_args()
    if args.tail_samples < 1:
        parser.error("--tail-samples must be at least 1")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    runs = load_by_replicate(args.input_csv, args.default_mmb)
    replicates = summarize_replicates(runs, args.tail_samples)
    aggregates = aggregate_parameters(replicates)
    write_csv(args.output_dir / "replicate_steady_state.csv", replicates)
    write_csv(args.output_dir / "parameter_averages.csv", aggregates)

    M_values = sorted({record["M"] for record in aggregates})
    A_values = sorted({record["A"] for record in aggregates})
    plotted = []
    for M in M_values:
        records = [record for record in aggregates if record["M"] == M]
        if len({record["A"] for record in records}) < 2:
            continue
        output = (
            args.output_dir if len(M_values) == 1
            else args.output_dir / f"A_lambda_M{M}"
        )
        plot_phase_diagrams(records, output, "A", f"M={M}")
        plotted.append(str(output))
    for A in A_values:
        records = [record for record in aggregates if record["A"] == A]
        if len({record["M"] for record in records}) < 2:
            continue
        output = args.output_dir / f"M_lambda_A{A}"
        plot_phase_diagrams(records, output, "M", f"A={A}")
        plotted.append(str(output))
    if not plotted:
        raise ValueError("No two-dimensional A x lambda or M x lambda slice found")
    print(
        f"Aggregated {len(replicates)} replicas into {len(aggregates)} "
        f"parameter cells; plotted {', '.join(plotted)}."
    )


if __name__ == "__main__":
    main()
