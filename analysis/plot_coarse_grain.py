"""
Turn cg_*.npz (from coarse_grain_fields.py) into the figures that decide what
the field theory must reproduce:

  fig 1  q_dom(t) and nx*(t) per run          -> selection vs coarsening
  fig 2  steady period vs lt, and A1200 vs A1500 -> steric length vs emergent wavelength
  fig 3  coarse-grained fields rho,p,m,sigma at early/mid/late snapshots
  (smectic fraction is printed per run         -> smectic vs crystal)

Pure matplotlib + numpy; reads the .npz files written by the job.
"""
import os
import glob
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def load_all(npz_dir):
    runs = []
    for path in sorted(glob.glob(os.path.join(npz_dir, "cg_*.npz"))):
        d = dict(np.load(path, allow_pickle=True))
        d["_path"] = path
        runs.append(d)
    return runs


def lt_to_float(lt):
    lt = str(lt)
    if lt in ("0", "lt0", "None"):
        return 0.0
    try:
        return float(lt.replace("e", "e"))
    except ValueError:
        return np.nan


def fig_qdom(runs, out):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))
    for d in runs:
        f = d["ts_frame"]
        nx = d["ts_nx_star"]
        phi = d["ts_phi"]
        lbl = f"A{int(d['n_actins'])} lt{d['lt']}"
        good = nx > 0
        ax1.plot(f[good], nx[good], lw=1, label=lbl)
        ax2.plot(f, phi, lw=1, label=lbl)
    ax1.set_xlabel("frame"); ax1.set_ylabel("winning $n_x$ (period $=L_x/n_x$)")
    ax1.set_title("nx*(t): decreasing => coarsening, flat => selection")
    ax2.set_xlabel("frame"); ax2.set_ylabel(r"$\Phi$ (order parameter)")
    ax2.set_title(r"$\Phi(t)$: pattern amplitude")
    ax2.legend(fontsize=7, ncol=2)
    fig.tight_layout(); fig.savefig(out, dpi=150); plt.close(fig)
    print(f"wrote {out}")


def fig_period(runs, out):
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    for d in runs:
        nx = d["ts_nx_star"]
        nx_tail = nx[len(nx) // 2:]
        nx_tail = nx_tail[nx_tail > 0]
        if nx_tail.size == 0:
            continue
        period = d["Lx"] / np.median(nx_tail)
        ax.scatter(lt_to_float(d["lt"]) + 1e-6, period,
                   s=60, label=f"A{int(d['n_actins'])}")
        ax.annotate(f"A{int(d['n_actins'])}", (lt_to_float(d['lt']) + 1e-6, period),
                    fontsize=7, xytext=(3, 3), textcoords="offset points")
    ax.set_xscale("symlog", linthresh=1e-3)
    ax.set_xlabel("catch-bond strength lt (= lifetime_coeff)")
    ax.set_ylabel("steady band period $L_x/n_x$")
    ax.set_title("period vs catch-bond strength & density\n"
                 "(flat in A => steric/overlap length; moves => emergent)")
    fig.tight_layout(); fig.savefig(out, dpi=150); plt.close(fig)
    print(f"wrote {out}")


def fig_fields(d, out):
    x = d["x"]
    snap_frames = d["snap_frames"]
    pick = [0, len(snap_frames) // 2, len(snap_frames) - 1]
    fig, axes = plt.subplots(4, 1, figsize=(7, 9), sharex=True)
    panels = [("snap_rho", r"$\rho$ (total)"),
              ("snap_p", r"$p$ (polarity $\rho_+-\rho_-$)"),
              ("snap_m_myo", r"$m$ (bound myosin, deduped)"),
              ("snap_sigma", r"$\sigma$ (crosslink tension)")]
    for ax, (key, ylabel) in zip(axes, panels):
        arr = d[key]
        for pi in pick:
            ax.plot(x, arr[pi], lw=1, label=f"frame {int(snap_frames[pi])}")
        ax.set_ylabel(ylabel)
    axes[0].legend(fontsize=7, ncol=3)
    axes[-1].set_xlabel("x (sarcomere axis)")
    axes[0].set_title(f"{str(d['fname'])}  (smectic frac={float(d['smectic_fraction']):.2f})")
    fig.tight_layout(); fig.savefig(out, dpi=150); plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("npz_dir")
    ap.add_argument("--out_dir", default=None)
    args = ap.parse_args()
    out_dir = args.out_dir or args.npz_dir
    os.makedirs(out_dir, exist_ok=True)

    runs = load_all(args.npz_dir)
    if not runs:
        raise SystemExit(f"no cg_*.npz in {args.npz_dir}")
    print(f"loaded {len(runs)} runs")
    for d in runs:
        print(f"  {str(d['fname']):50s} smectic_frac={float(d['smectic_fraction']):.2f}")

    fig_qdom(runs, os.path.join(out_dir, "fig1_qdom_timeseries.png"))
    fig_period(runs, os.path.join(out_dir, "fig2_period_vs_lt_density.png"))
    for d in runs:
        stem = os.path.splitext(os.path.basename(str(d["_path"])))[0]
        fig_fields(d, os.path.join(out_dir, f"fig3_fields_{stem}.png"))
