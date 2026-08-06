"""
Aggregate cg_*.npz (from coarse_grain_fields.py) over the full A x lt x replica
grid into phase diagrams:

  table_phase.csv          one row per run: A, lt, rep, steady_nx, period,
                           smectic_fraction, phi_final
  fig_period_vs_density    period vs A, one line per lt (replicas avg +/- sd)
                           -> the emergent-wavelength vs packing-law test
  fig_smectic_heatmap      smectic fraction over (lt rows, A cols)
  fig_nx_heatmap           steady nx* over (lt rows, A cols)
  fig_period_vs_lt         period vs lt, one line per A

Robust lt parsing: labels include non-integer exponents (1e-0.5, 1e-2.5, ...),
which float() cannot parse -> handled as mantissa * 10**float(exp).
"""
import os
import re
import glob
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from collections import defaultdict


def lt_to_float(lt):
    lt = str(lt).strip()
    if lt in ("0", "lt0", "None", ""):
        return 0.0
    if "e" in lt:
        mant, exp = lt.split("e")
        return float(mant) * 10.0 ** float(exp)   # exp may be non-integer
    return float(lt)


def rep_from_name(name):
    m = re.search(r"_(r\d+)_", name)
    return m.group(1) if m else "r??"


def steady_nx(ts_nx):
    tail = ts_nx[len(ts_nx) // 2:]
    tail = tail[tail > 0]
    return float(np.median(tail)) if tail.size else 0.0


def load_records(npz_dir):
    recs = []
    for path in sorted(glob.glob(os.path.join(npz_dir, "cg_*.npz"))):
        d = np.load(path, allow_pickle=True)
        name = str(d["fname"])
        A = int(d["n_actins"])
        lt = str(d["lt"])
        Lx = float(d["Lx"])
        nx = steady_nx(d["ts_nx_star"])
        phi = d["ts_phi"]
        recs.append({
            "name": name, "A": A, "lt": lt, "lt_val": lt_to_float(lt),
            "rep": rep_from_name(name),
            "steady_nx": nx,
            "period": (Lx / nx) if nx > 0 else np.nan,
            "smectic_fraction": float(d["smectic_fraction"]),
            "phi_final": float(np.mean(phi[len(phi) // 2:])),
        })
    return recs


def write_csv(recs, out):
    import csv
    fields = ["name", "A", "lt", "lt_val", "rep", "steady_nx", "period",
              "smectic_fraction", "phi_final"]
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(recs)
    print(f"wrote {out}  ({len(recs)} runs)")


def _grid(recs, value_key):
    """Return (A_sorted, lt_sorted, lt_labels, mean_grid) averaging over reps."""
    As = sorted({r["A"] for r in recs})
    lts = sorted({r["lt_val"] for r in recs})
    lt_label = {}
    for r in recs:
        lt_label[r["lt_val"]] = r["lt"]
    cells = defaultdict(list)
    for r in recs:
        v = r[value_key]
        if v is not None and np.isfinite(v):
            cells[(r["lt_val"], r["A"])].append(v)
    G = np.full((len(lts), len(As)), np.nan)
    for i, lv in enumerate(lts):
        for j, A in enumerate(As):
            vals = cells.get((lv, A))
            if vals:
                G[i, j] = np.mean(vals)
    return As, lts, [lt_label[lv] for lv in lts], G


def fig_period_vs_density(recs, out):
    As, lts, labels, _ = _grid(recs, "period")
    fig, ax = plt.subplots(figsize=(7, 5))
    cmap = plt.cm.viridis(np.linspace(0, 1, len(lts)))
    for i, (lv, lab) in enumerate(zip(lts, labels)):
        xs, ys, es = [], [], []
        for A in As:
            vals = [r["period"] for r in recs
                    if r["A"] == A and r["lt_val"] == lv
                    and np.isfinite(r["period"])]
            if vals:
                xs.append(A); ys.append(np.mean(vals)); es.append(np.std(vals))
        if xs:
            ax.errorbar(xs, ys, yerr=es, marker="o", color=cmap[i],
                        label=f"lt={lab}", capsize=2, lw=1)
    ax.set_xlabel("actin count A"); ax.set_ylabel("steady period $L_x/n_x$")
    ax.set_title("period vs density (per lt)\nflat => steric length; "
                 "$\\propto 1/A$ => fixed cluster; steeper => emergent")
    ax.legend(fontsize=7, ncol=2)
    fig.tight_layout(); fig.savefig(out, dpi=150); plt.close(fig)
    print(f"wrote {out}")


def fig_period_vs_lt(recs, out):
    As = sorted({r["A"] for r in recs})
    fig, ax = plt.subplots(figsize=(7, 5))
    cmap = plt.cm.plasma(np.linspace(0, 1, len(As)))
    for k, A in enumerate(As):
        pts = sorted([(r["lt_val"], r["period"]) for r in recs
                      if r["A"] == A and np.isfinite(r["period"])])
        if not pts:
            continue
        agg = defaultdict(list)
        for lv, p in pts:
            agg[lv].append(p)
        xs = sorted(agg)
        ys = [np.mean(agg[x]) for x in xs]
        ax.plot([x + 1e-7 for x in xs], ys, marker="o", color=cmap[k],
                label=f"A{A}", lw=1)
    ax.set_xscale("symlog", linthresh=1e-5)
    ax.set_xlabel("catch-bond strength lt"); ax.set_ylabel("steady period $L_x/n_x$")
    ax.set_title("period vs catch-bond strength (per density)")
    ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(out, dpi=150); plt.close(fig)
    print(f"wrote {out}")


def fig_heatmap(recs, value_key, title, out):
    As, lts, labels, G = _grid(recs, value_key)
    fig, ax = plt.subplots(figsize=(7, 6))
    im = ax.imshow(G, aspect="auto", origin="lower", cmap="magma")
    ax.set_xticks(range(len(As))); ax.set_xticklabels(As)
    ax.set_yticks(range(len(labels))); ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel("actin count A"); ax.set_ylabel("catch-bond strength lt")
    ax.set_title(title)
    fig.colorbar(im, ax=ax)
    fig.tight_layout(); fig.savefig(out, dpi=150); plt.close(fig)
    print(f"wrote {out}")


def load_tc_records(npz_dir):
    """Records from transverse_coherence.py outputs (tc_*.npz)."""
    recs = []
    for path in sorted(glob.glob(os.path.join(npz_dir, "tc_*.npz"))):
        d = np.load(path, allow_pickle=True)
        name = str(d["fname"])
        recs.append({
            "name": name, "A": int(d["n_actins"]), "lt": str(d["lt"]),
            "lt_val": lt_to_float(str(d["lt"])), "rep": rep_from_name(name),
            "R_steady": float(d["R_steady"]),
            "R_sigma_steady": float(d["R_sigma_steady"]),
            "nx_star": int(d["nx_star"]),
        })
    return recs


def fig_R_vs_lt(recs, out):
    As = sorted({r["A"] for r in recs})
    fig, ax = plt.subplots(figsize=(7, 5))
    cmap = plt.cm.plasma(np.linspace(0, 1, len(As)))
    for k, A in enumerate(As):
        agg = defaultdict(list)
        for r in recs:
            if r["A"] == A and np.isfinite(r["R_steady"]):
                agg[r["lt_val"]].append(r["R_steady"])
        if not agg:
            continue
        xs = sorted(agg)
        ys = [np.mean(agg[x]) for x in xs]
        es = [np.std(agg[x]) for x in xs]
        ax.errorbar([x + 1e-7 for x in xs], ys, yerr=es, marker="o",
                    color=cmap[k], label=f"A{A}", capsize=2, lw=1)
    ax.set_xscale("symlog", linthresh=1e-5)
    ax.set_ylim(0, 1)
    ax.set_xlabel("catch-bond strength lt")
    ax.set_ylabel("registration $R$")
    ax.set_title("registration vs catch-bond strength (replicas avg)\n"
                 "->1 registered crystal / ->0 dephased smectic")
    ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(out, dpi=150); plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("npz_dir")
    ap.add_argument("--out_dir", default=None)
    args = ap.parse_args()
    out_dir = args.out_dir or args.npz_dir
    os.makedirs(out_dir, exist_ok=True)

    recs = load_records(args.npz_dir)
    if not recs:
        raise SystemExit(f"no cg_*.npz in {args.npz_dir}")
    print(f"loaded {len(recs)} runs across "
          f"{len({r['A'] for r in recs})} densities x "
          f"{len({r['lt'] for r in recs})} lt x "
          f"{len({r['rep'] for r in recs})} reps")

    write_csv(recs, os.path.join(out_dir, "table_phase.csv"))
    fig_period_vs_density(recs, os.path.join(out_dir, "fig_period_vs_density.png"))
    fig_period_vs_lt(recs, os.path.join(out_dir, "fig_period_vs_lt.png"))
    fig_heatmap(recs, "smectic_fraction",
                "smectic fraction (->1 layered / quasi-1D)",
                os.path.join(out_dir, "fig_smectic_heatmap.png"))
    fig_heatmap(recs, "steady_nx", "steady winning $n_x$ (period $=L_x/n_x$)",
                os.path.join(out_dir, "fig_nx_heatmap.png"))

    # transverse registration R (if transverse_coherence.py was run)
    tc = load_tc_records(args.npz_dir)
    if tc:
        print(f"loaded {len(tc)} transverse-coherence runs")
        import csv as _csv
        with open(os.path.join(out_dir, "table_transverse.csv"), "w", newline="") as f:
            w = _csv.DictWriter(f, fieldnames=["name", "A", "lt", "lt_val", "rep",
                                               "nx_star", "R_steady", "R_sigma_steady"])
            w.writeheader(); w.writerows(tc)
        fig_heatmap(tc, "R_steady",
                    "registration $R$ (->1 crystal / ->0 smectic, reps avg)",
                    os.path.join(out_dir, "fig_R_heatmap.png"))
        fig_heatmap(tc, "R_sigma_steady",
                    "$\\sigma$-field registration $R_\\sigma$ (reps avg)",
                    os.path.join(out_dir, "fig_Rsigma_heatmap.png"))
        fig_R_vs_lt(tc, os.path.join(out_dir, "fig_R_vs_lt.png"))
