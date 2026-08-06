"""Phase diagrams of entropy-production rates from entropy_transduction_results CSV."""

from __future__ import annotations

import csv
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np

CSV_PATH = Path(__file__).parent / "entropy_transduction_results_data0528.csv"
OUT_DIR = Path(__file__).parent

# ── helpers ──────────────────────────────────────────────────────────────────

def _parse_lt(lt_str: str) -> float:
    """Convert lt_str like '1e-3', '1e-1.5', or '0' to float."""
    if lt_str == "0":
        return 0.0
    m = re.fullmatch(r"1e([+-]?\d+(?:\.\d+)?)", lt_str)
    if m:
        return 10.0 ** float(m.group(1))
    return float(lt_str)


def _lt_label(lt_str: str) -> str:
    if lt_str == "0":
        return "0"
    m = re.fullmatch(r"1e([+-]?\d+(?:\.\d+)?)", lt_str)
    if m:
        exp = m.group(1)
        return f"$10^{{{exp}}}$"
    return lt_str


def load_data(csv_path: Path) -> dict:
    with open(csv_path, newline="") as f:
        rows = list(csv.DictReader(f))

    a_vals = sorted({int(r["A_val"]) for r in rows})

    # Order lt values: '0' first, then by actual float value ascending
    lt_strs_raw = sorted({r["lt_str"] for r in rows},
                         key=lambda s: (_parse_lt(s), s))
    # Build lookup (a_val, lt_str) → row
    lookup: dict[tuple[int, str], dict] = {}
    for r in rows:
        lookup[(int(r["A_val"]), r["lt_str"])] = r

    return {"a_vals": a_vals, "lt_strs": lt_strs_raw, "lookup": lookup}


def build_grid(data: dict, field: str) -> np.ndarray:
    a_vals = data["a_vals"]
    lt_strs = data["lt_strs"]
    lookup = data["lookup"]
    grid = np.full((len(a_vals), len(lt_strs)), np.nan)
    for i, a in enumerate(a_vals):
        for j, lt in enumerate(lt_strs):
            row = lookup.get((a, lt))
            if row is not None:
                try:
                    v = float(row[field])
                    grid[i, j] = v if np.isfinite(v) else np.nan
                except (ValueError, KeyError):
                    pass
    return grid


# ── plotting ─────────────────────────────────────────────────────────────────

PANELS = [
    ("Mech_Rate",       r"Mech. entropy rate $\dot{S}_\mathrm{mech}$"),
    ("Chem_Rate",       r"Chem. entropy rate $\dot{S}_\mathrm{chem}$"),
    ("Total_Rate",      r"Total entropy rate $\dot{S}_\mathrm{tot}$"),
    ("Active_Work_Rate",r"Active work rate $\dot{W}_\mathrm{active}$"),
]

PANELS_PER_ACTIN = [
    ("Mech_Rate_per_actin",        r"$\dot{S}_\mathrm{mech}$ / actin"),
    ("Chem_Rate_per_actin",        r"$\dot{S}_\mathrm{chem}$ / actin"),
    ("Total_Rate_per_actin",       r"$\dot{S}_\mathrm{tot}$ / actin"),
    ("Active_Work_Rate_per_actin", r"$\dot{W}_\mathrm{active}$ / actin"),
]


def plot_entropy_phase_diagrams(data: dict, out_path: Path) -> None:
    a_vals  = data["a_vals"]
    lt_strs = data["lt_strs"]
    x_idx   = np.arange(len(lt_strs))
    y_idx   = np.arange(len(a_vals))
    x_labels = [_lt_label(s) for s in lt_strs]
    y_labels = [str(a) for a in a_vals]

    fig, axes = plt.subplots(2, 2, figsize=(13, 9), constrained_layout=True)
    fig.suptitle("Entropy production rate phase diagrams", fontsize=14)

    for ax, (field, title) in zip(axes.flat, PANELS):
        grid = build_grid(data, field)

        # symmetric log normalisation so zero is visible
        vmax = np.nanmax(np.abs(grid)) if not np.all(np.isnan(grid)) else 1.0
        vmin = np.nanmin(grid[grid > 0]) if np.any(grid > 0) else 1e-6

        norm = mcolors.SymLogNorm(linthresh=vmin, vmin=0, vmax=vmax, base=10)
        im = ax.pcolormesh(
            np.arange(len(lt_strs) + 1) - 0.5,
            np.arange(len(a_vals) + 1) - 0.5,
            grid,
            norm=norm,
            cmap="plasma",
            shading="flat",
        )
        cb = fig.colorbar(im, ax=ax, pad=0.02, fraction=0.046)
        cb.set_label("rate (kBT / time)", fontsize=8)

        # annotate cells with values
        for i in range(len(a_vals)):
            for j in range(len(lt_strs)):
                v = grid[i, j]
                if not np.isnan(v):
                    txt = f"{v:.0f}" if abs(v) >= 1 else f"{v:.2g}"
                    brightness = im.norm(v) if not np.isnan(im.norm(v)) else 0
                    color = "white" if brightness < 0.55 else "black"
                    ax.text(j, i, txt, ha="center", va="center",
                            fontsize=5.5, color=color)

        ax.set_xticks(x_idx)
        ax.set_xticklabels(x_labels, rotation=45, ha="right", fontsize=7)
        ax.set_yticks(y_idx)
        ax.set_yticklabels(y_labels, fontsize=8)
        ax.set_xlabel(r"Lifetime coeff $\lambda_t$", fontsize=9)
        ax.set_ylabel("Actin amplitude A", fontsize=9)
        ax.set_title(title, fontsize=10)

    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Saved: {out_path}")


def plot_log_ratio_diagram(data: dict, out_path: Path) -> None:
    """Log-ratio Chem/Mech rate — shows which channel dominates."""
    a_vals  = data["a_vals"]
    lt_strs = data["lt_strs"]
    x_idx   = np.arange(len(lt_strs))
    y_idx   = np.arange(len(a_vals))

    chem  = build_grid(data, "Chem_Rate")
    mech  = build_grid(data, "Mech_Rate")
    with np.errstate(divide="ignore", invalid="ignore"):
        # sign(chem) * log10(|chem|/mech): negative chem → negative ratio (blue),
        # so anti-entropic chemical production stays on the mech-dominated side.
        ratio = np.where(
            mech > 0,
            np.sign(chem) * np.log10(np.abs(chem) / mech),
            np.nan,
        )

    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
    vabs = np.nanmax(np.abs(ratio[np.isfinite(ratio)])) if np.any(np.isfinite(ratio)) else 1.0
    norm = mcolors.TwoSlopeNorm(vmin=-vabs, vcenter=0, vmax=vabs)
    im = ax.pcolormesh(
        np.arange(len(lt_strs) + 1) - 0.5,
        np.arange(len(a_vals) + 1) - 0.5,
        ratio,
        norm=norm,
        cmap="RdBu_r",
        shading="flat",
    )
    cb = fig.colorbar(im, ax=ax, pad=0.02)
    cb.set_label(r"$\log_{10}(\dot{S}_\mathrm{chem}/\dot{S}_\mathrm{mech})$", fontsize=9)

    for i in range(len(a_vals)):
        for j in range(len(lt_strs)):
            v = ratio[i, j]
            if np.isfinite(v):
                brightness = norm(v)
                color = "white" if 0.3 < brightness < 0.7 else "black"
                ax.text(j, i, f"{v:.1f}", ha="center", va="center",
                        fontsize=6, color=color)

    ax.set_xticks(x_idx)
    ax.set_xticklabels([_lt_label(s) for s in lt_strs], rotation=45, ha="right", fontsize=7)
    ax.set_yticks(y_idx)
    ax.set_yticklabels([str(a) for a in a_vals], fontsize=8)
    ax.set_xlabel(r"Lifetime coeff $\lambda_t$", fontsize=9)
    ax.set_ylabel("Actin amplitude A", fontsize=9)
    ax.set_title(r"$\mathrm{sgn}(\dot{S}_\mathrm{chem})\,\log_{10}(|\dot{S}_\mathrm{chem}|/\dot{S}_\mathrm{mech})$"
                 "\nblue: mech dominates or $\\dot{S}_\\mathrm{chem}<0$;  red: chem dominates", fontsize=9)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Saved: {out_path}")


def plot_individual_panels(data: dict, out_dir: Path) -> None:
    """Save each entropy quantity as its own SVG phase diagram."""
    a_vals  = data["a_vals"]
    lt_strs = data["lt_strs"]
    x_idx   = np.arange(len(lt_strs))
    y_idx   = np.arange(len(a_vals))
    x_labels = [_lt_label(s) for s in lt_strs]

    fields = [
        ("Mech_Rate",        r"$\dot{S}_\mathrm{mech}$",  "plasma"),
        ("Chem_Rate",        r"$\dot{S}_\mathrm{chem}$",  "viridis"),
        ("Total_Rate",       r"$\dot{S}_\mathrm{tot}$",   "inferno"),
        ("Active_Work_Rate", r"$\dot{W}_\mathrm{active}$","magma"),
    ]
    for field, label, cmap in fields:
        grid = build_grid(data, field)
        fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
        vmax = np.nanmax(np.abs(grid)) if not np.all(np.isnan(grid)) else 1.0
        vmin_pos = np.nanmin(grid[grid > 0]) if np.any(grid > 0) else 1e-6
        norm = mcolors.SymLogNorm(linthresh=vmin_pos, vmin=0, vmax=vmax, base=10)
        im = ax.pcolormesh(
            np.arange(len(lt_strs) + 1) - 0.5,
            np.arange(len(a_vals) + 1) - 0.5,
            grid, norm=norm, cmap=cmap, shading="flat",
        )
        cb = fig.colorbar(im, ax=ax, pad=0.02)
        cb.set_label(f"{label}  (kBT / time)", fontsize=9)
        for i in range(len(a_vals)):
            for j in range(len(lt_strs)):
                v = grid[i, j]
                if not np.isnan(v):
                    txt = f"{v:.0f}" if abs(v) >= 1 else f"{v:.2g}"
                    brightness = float(im.norm(v))
                    color = "white" if brightness < 0.5 else "black"
                    ax.text(j, i, txt, ha="center", va="center",
                            fontsize=6, color=color)
        ax.set_xticks(x_idx)
        ax.set_xticklabels(x_labels, rotation=45, ha="right", fontsize=7)
        ax.set_yticks(y_idx)
        ax.set_yticklabels([str(a) for a in a_vals], fontsize=8)
        ax.set_xlabel(r"Lifetime coeff $\lambda_t$", fontsize=9)
        ax.set_ylabel("Actin amplitude A", fontsize=9)
        ax.set_title(f"Phase diagram: {label}", fontsize=11)
        svg = out_dir / f"PhaseDiagram_Entropy_{field}.svg"
        fig.savefig(svg, bbox_inches="tight")
        png = out_dir / f"PhaseDiagram_Entropy_{field}.png"
        fig.savefig(png, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved: {svg}")


def _plot_panel_grid(
    data: dict,
    panels: list[tuple[str, str]],
    out_dir: Path,
    *,
    title_prefix: str = "Phase diagram",
    cmaps: list[str] | None = None,
    combined_out: Path | None = None,
    combined_title: str = "Entropy production rate phase diagrams (per actin)",
) -> None:
    """Generic helper: one SVG per panel + optional combined 2×2 figure."""
    default_cmaps = ["plasma", "viridis", "inferno", "magma"]
    if cmaps is None:
        cmaps = default_cmaps

    a_vals  = data["a_vals"]
    lt_strs = data["lt_strs"]
    x_idx   = np.arange(len(lt_strs))
    y_idx   = np.arange(len(a_vals))
    x_labels = [_lt_label(s) for s in lt_strs]

    # combined figure
    if combined_out is not None:
        fig_c, axes_c = plt.subplots(2, 2, figsize=(13, 9), constrained_layout=True)
        fig_c.suptitle(combined_title, fontsize=14)
        combined_axes = list(axes_c.flat)
    else:
        combined_axes = [None] * len(panels)

    for (field, label), cmap, cax in zip(panels, cmaps, combined_axes):
        grid = build_grid(data, field)

        vmax = np.nanmax(np.abs(grid)) if not np.all(np.isnan(grid)) else 1.0
        vmin_pos = np.nanmin(grid[grid > 0]) if np.any(grid > 0) else 1e-6
        norm = mcolors.SymLogNorm(linthresh=vmin_pos, vmin=0, vmax=vmax, base=10)

        for ax in ([cax] if cax is not None else []) + [None]:
            if ax is None:
                # individual figure
                fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
                standalone = True
            else:
                standalone = False

            im = ax.pcolormesh(
                np.arange(len(lt_strs) + 1) - 0.5,
                np.arange(len(a_vals) + 1) - 0.5,
                grid, norm=norm, cmap=cmap, shading="flat",
            )
            cb = (fig_c if not standalone else fig).colorbar(im, ax=ax, pad=0.02)
            cb.set_label(f"{label}  (kBT / time)", fontsize=8 if not standalone else 9)
            for i in range(len(a_vals)):
                for j in range(len(lt_strs)):
                    v = grid[i, j]
                    if not np.isnan(v):
                        txt = f"{v:.1f}" if abs(v) >= 1 else f"{v:.2g}"
                        brightness = float(im.norm(v))
                        color = "white" if brightness < 0.5 else "black"
                        ax.text(j, i, txt, ha="center", va="center",
                                fontsize=5.5 if not standalone else 6, color=color)
            ax.set_xticks(x_idx)
            ax.set_xticklabels(x_labels, rotation=45, ha="right",
                               fontsize=7)
            ax.set_yticks(y_idx)
            ax.set_yticklabels([str(a) for a in a_vals], fontsize=8)
            ax.set_xlabel(r"Lifetime coeff $\lambda_t$", fontsize=9)
            ax.set_ylabel("Actin amplitude A", fontsize=9)
            ax.set_title(f"{title_prefix}: {label}", fontsize=10 if not standalone else 11)

            if standalone:
                clean = field.replace("/", "_")
                svg = out_dir / f"PhaseDiagram_Entropy_{clean}.svg"
                png = out_dir / f"PhaseDiagram_Entropy_{clean}.png"
                fig.savefig(svg, bbox_inches="tight")
                fig.savefig(png, dpi=150, bbox_inches="tight")
                plt.close(fig)
                print(f"Saved: {svg}")

    if combined_out is not None:
        fig_c.savefig(combined_out, dpi=150, bbox_inches="tight")
        print(f"Saved: {combined_out}")
        plt.close(fig_c)


def _signed_log_ratio_grid(data: dict, num_field: str, den_field: str) -> np.ndarray:
    """sign(den) * log10(|num/den|), NaN where den==0 or either is NaN."""
    num = build_grid(data, num_field)
    den = build_grid(data, den_field)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(
            (den != 0) & np.isfinite(num) & np.isfinite(den),
            np.sign(den) * np.log10(np.abs(num) / np.where(den != 0, np.abs(den), np.nan)),
            np.nan,
        )
    return ratio


def _annotate_grid(ax, grid, im, fontsize=6):
    for i in range(grid.shape[0]):
        for j in range(grid.shape[1]):
            v = grid[i, j]
            if np.isfinite(v):
                txt = f"{v:.1f}" if abs(v) >= 1 else f"{v:.2g}"
                brightness = float(im.norm(v))
                color = "white" if brightness < 0.5 else "black"
                ax.text(j, i, txt, ha="center", va="center", fontsize=fontsize, color=color)


def _phase_ax(ax, fig, data, field, title, cmap, *, colorbar_label="kBT / time",
              ratio_mode=False, fontsize_annot=6):
    """Draw one phase diagram panel onto ax."""
    a_vals  = data["a_vals"]
    lt_strs = data["lt_strs"]
    grid = build_grid(data, field) if not ratio_mode else field  # pass grid directly in ratio mode

    if ratio_mode:
        vabs = np.nanmax(np.abs(grid[np.isfinite(grid)])) if np.any(np.isfinite(grid)) else 1.0
        norm = mcolors.TwoSlopeNorm(vmin=-vabs, vcenter=0, vmax=vabs)
        cmap_use = "RdBu_r"
    else:
        vmax = np.nanmax(np.abs(grid)) if not np.all(np.isnan(grid)) else 1.0
        vmin_pos = np.nanmin(grid[grid > 0]) if np.any(grid > 0) else 1e-6
        norm = mcolors.SymLogNorm(linthresh=vmin_pos, vmin=0, vmax=vmax, base=10)
        cmap_use = cmap

    im = ax.pcolormesh(
        np.arange(len(lt_strs) + 1) - 0.5,
        np.arange(len(a_vals) + 1) - 0.5,
        grid, norm=norm, cmap=cmap_use, shading="flat",
    )
    cb = fig.colorbar(im, ax=ax, pad=0.02)
    cb.set_label(colorbar_label, fontsize=8)
    _annotate_grid(ax, grid, im, fontsize=fontsize_annot)
    ax.set_xticks(np.arange(len(lt_strs)))
    ax.set_xticklabels([_lt_label(s) for s in lt_strs], rotation=45, ha="right", fontsize=7)
    ax.set_yticks(np.arange(len(a_vals)))
    ax.set_yticklabels([str(a) for a in a_vals], fontsize=8)
    ax.set_xlabel(r"Lifetime coeff $\lambda_t$", fontsize=9)
    ax.set_ylabel("Actin amplitude A", fontsize=9)
    ax.set_title(title, fontsize=9)


def plot_steadystate_figure(data: dict, out_path: Path) -> None:
    """Figure 1: steady-state rates (tail-300-frame mean) per filament.

    Panels: per-actin mech, per-actin chem, sign*log(mech/chem) per actin, per-myosin mech.
    Current tail mean uses every frame in the last 300; regenerate CSV with --stride 10
    to subsample every 10th frame (30 quasi-independent samples).
    """
    fig, axes = plt.subplots(2, 2, figsize=(14, 9), constrained_layout=True)
    fig.suptitle(
        r"Steady-state entropy production rates (last-300-frame mean, per filament)",
        fontsize=12,
    )

    mech_grid = build_grid(data, "SS_Actin_Mech_Rate_per_actin")
    chem_grid = build_grid(data, "SS_Chem_Rate_per_actin")
    myosin_grid = build_grid(data, "SS_Myosin_Mech_Rate_per_myosin")
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio_grid = np.where(
            (chem_grid != 0) & np.isfinite(chem_grid) & np.isfinite(mech_grid),
            np.sign(chem_grid) * np.log10(np.abs(mech_grid) / np.where(chem_grid != 0, np.abs(chem_grid), np.nan)),
            np.nan,
        )

    _phase_ax(axes[0, 0], fig, data, "SS_Actin_Mech_Rate_per_actin",
              r"$\dot{S}_\mathrm{mech}^\mathrm{actin}$ per actin  (kBT/time)", "plasma")
    _phase_ax(axes[0, 1], fig, data, "SS_Chem_Rate_per_actin",
              r"$\dot{S}_\mathrm{chem}$ per actin  (kBT/time)", "viridis")
    _phase_ax(axes[1, 0], fig, data, ratio_grid,
              r"$\mathrm{sgn}(\dot{S}_\mathrm{chem})\,\log_{10}(|\dot{S}_\mathrm{mech}^\mathrm{actin}/\dot{S}_\mathrm{chem}|)$",
              None, colorbar_label="log ratio", ratio_mode=True)
    _phase_ax(axes[1, 1], fig, data, "SS_Myosin_Mech_Rate_per_myosin",
              r"$\dot{S}_\mathrm{mech}^\mathrm{myosin}$ per myosin  (kBT/time)", "magma")

    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    svg_path = out_path.with_suffix(".svg")
    fig.savefig(svg_path, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_path}")
    print(f"Saved: {svg_path}")


def plot_totaltime_figure(data: dict, out_path: Path) -> None:
    """Figure 2: total entropy production / total time, per filament.

    Panels: per-actin mech, per-actin chem, sign*log(mech/chem), per-myosin mech.
    """
    fig, axes = plt.subplots(2, 2, figsize=(14, 9), constrained_layout=True)
    fig.suptitle(
        r"Total entropy production / total time (per filament)",
        fontsize=12,
    )

    mech_grid = build_grid(data, "TT_Actin_Mech_Rate_per_actin")
    chem_grid = build_grid(data, "TT_Chem_Rate_per_actin")
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio_grid = np.where(
            (chem_grid != 0) & np.isfinite(chem_grid) & np.isfinite(mech_grid),
            np.sign(chem_grid) * np.log10(np.abs(mech_grid) / np.where(chem_grid != 0, np.abs(chem_grid), np.nan)),
            np.nan,
        )

    _phase_ax(axes[0, 0], fig, data, "TT_Actin_Mech_Rate_per_actin",
              r"$\dot{S}_\mathrm{mech}^\mathrm{actin}$ per actin  (kBT/time)", "plasma")
    _phase_ax(axes[0, 1], fig, data, "TT_Chem_Rate_per_actin",
              r"$\dot{S}_\mathrm{chem}$ per actin  (kBT/time)", "viridis")
    _phase_ax(axes[1, 0], fig, data, ratio_grid,
              r"$\mathrm{sgn}(\dot{S}_\mathrm{chem})\,\log_{10}(|\dot{S}_\mathrm{mech}^\mathrm{actin}/\dot{S}_\mathrm{chem}|)$",
              None, colorbar_label="log ratio", ratio_mode=True)
    _phase_ax(axes[1, 1], fig, data, "TT_Myosin_Mech_Rate_per_myosin",
              r"$\dot{S}_\mathrm{mech}^\mathrm{myosin}$ per myosin  (kBT/time)", "magma")

    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    svg_path = out_path.with_suffix(".svg")
    fig.savefig(svg_path, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_path}")
    print(f"Saved: {svg_path}")


def _ratio_grid(data, mech_field, chem_field):
    mech = build_grid(data, mech_field)
    chem = build_grid(data, chem_field)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(
            (chem != 0) & np.isfinite(chem) & np.isfinite(mech),
            np.sign(chem) * np.log10(np.abs(mech) / np.where(chem != 0, np.abs(chem), np.nan)),
            np.nan,
        )


def _make_4panel(title, panels_cmaps, ratio_grid):
    """Build a 2×2 figure given (field, label, cmap) for panels [0,1,3] and a pre-built ratio grid."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 9), constrained_layout=True)
    fig.suptitle(title, fontsize=12)
    (f0, l0, c0), (f1, l1, c1), (f3, l3, c3) = panels_cmaps
    return fig, axes, [(f0, l0, c0), (f1, l1, c1), (f3, l3, c3)]


def plot_engaged_steadystate_figure(data: dict, out_path: Path) -> None:
    """Steady-state rates normalised by total n_actins (A_val)."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 9), constrained_layout=True)
    fig.suptitle(
        r"Steady-state entropy production rates per actin / per myosin",
        fontsize=12,
    )
    rg = _ratio_grid(data, "SS_Actin_Mech_Rate_per_actin", "SS_Chem_Rate_per_actin")
    _phase_ax(axes[0, 0], fig, data, "SS_Actin_Mech_Rate_per_actin",
              r"$\dot{S}_\mathrm{mech}^\mathrm{actin}$ / actin  (kBT/time)", "plasma")
    _phase_ax(axes[0, 1], fig, data, "SS_Chem_Rate_per_actin",
              r"$\dot{S}_\mathrm{chem}$ / actin  (kBT/time)", "viridis")
    _phase_ax(axes[1, 0], fig, data, rg,
              r"$\mathrm{sgn}(\dot{S}_\mathrm{chem})\,\log_{10}(|\dot{S}_\mathrm{mech}/\dot{S}_\mathrm{chem}|)$  per actin",
              None, colorbar_label="log ratio", ratio_mode=True)
    _phase_ax(axes[1, 1], fig, data, "SS_Myosin_Mech_Rate_per_myosin",
              r"$\dot{S}_\mathrm{mech}^\mathrm{myosin}$ / myosin  (kBT/time)", "magma")
    for ext in (".png", ".svg"):
        p = out_path.with_suffix(ext)
        fig.savefig(p, dpi=150 if ext == ".png" else None, bbox_inches="tight")
        print(f"Saved: {p}")
    plt.close(fig)


def plot_engaged_totaltime_figure(data: dict, out_path: Path) -> None:
    """Total entropy / time normalised by total n_actins (A_val)."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 9), constrained_layout=True)
    fig.suptitle(
        r"Total entropy production / time per actin / per myosin",
        fontsize=12,
    )
    rg = _ratio_grid(data, "TT_Actin_Mech_Rate_per_actin", "TT_Chem_Rate_per_actin")
    _phase_ax(axes[0, 0], fig, data, "TT_Actin_Mech_Rate_per_actin",
              r"$\dot{S}_\mathrm{mech}^\mathrm{actin}$ / actin  (kBT/time)", "plasma")
    _phase_ax(axes[0, 1], fig, data, "TT_Chem_Rate_per_actin",
              r"$\dot{S}_\mathrm{chem}$ / actin  (kBT/time)", "viridis")
    _phase_ax(axes[1, 0], fig, data, rg,
              r"$\mathrm{sgn}(\dot{S}_\mathrm{chem})\,\log_{10}(|\dot{S}_\mathrm{mech}/\dot{S}_\mathrm{chem}|)$  per actin",
              None, colorbar_label="log ratio", ratio_mode=True)
    _phase_ax(axes[1, 1], fig, data, "TT_Myosin_Mech_Rate_per_myosin",
              r"$\dot{S}_\mathrm{mech}^\mathrm{myosin}$ / myosin  (kBT/time)", "magma")
    for ext in (".png", ".svg"):
        p = out_path.with_suffix(ext)
        fig.savefig(p, dpi=150 if ext == ".png" else None, bbox_inches="tight")
        print(f"Saved: {p}")
    plt.close(fig)


if __name__ == "__main__":
    data = load_data(CSV_PATH)
    plot_entropy_phase_diagrams(data, OUT_DIR / "Fig_EntropyRate_Component_PhaseDiagrams.svg")
    plot_log_ratio_diagram(data, OUT_DIR / "PhaseDiagram_Entropy_LogRatio.svg")
    plot_individual_panels(data, OUT_DIR)
    # per-actin normalized versions
    _plot_panel_grid(
        data, PANELS_PER_ACTIN, OUT_DIR,
        title_prefix="Per-actin phase diagram",
        combined_out=OUT_DIR / "Fig_EntropyRate_PerActin_PhaseDiagrams.svg",
        combined_title="Entropy production rate phase diagrams (per actin filament)",
    )
    # focused per-filament breakdown figures
    plot_steadystate_figure(data, OUT_DIR / "Fig_SteadyState_PerFilament.png")
    plot_totaltime_figure(data, OUT_DIR / "Fig_TotalTime_PerFilament.png")
    plot_engaged_steadystate_figure(data, OUT_DIR / "Fig_SteadyState_PerEngaged.png")
    plot_engaged_totaltime_figure(data, OUT_DIR / "Fig_TotalTime_PerEngaged.png")
    print("Done.")
