"""
Decouple PATTERN EXISTENCE from PATTERN SPACING, then ask whether the band
period is steric (flat) or interaction-set (moves) -- separately vs density and
vs catch-bond strength lt, RESTRICTED to where a pattern actually exists.

Motivation: period = Lx/nx is meaningless where no pattern forms. phi_final is a
BAD existence gate -- at low lt there are few cb2 actins and the Fourier amplitude
of a handful of points is spuriously high (small-N), giving a fake "strong pattern"
at a large spurious period. So we gate on a robust LONGITUDINAL existence metric:

  nx_stability = fraction of steady-window frames whose winning nx is within +/-1
                 of the modal nx  (a real pattern sits at a stable band number;
                 noise jumps around)
  n_cb2        = median number of doubly-engaged actins in the window
                 (need enough particles for the mode to be meaningful)

Existence gate = (nx_stability >= STAB_MIN) and (n_cb2 >= NCB2_MIN).
Spacing        = Lx / modal_nx, computed only for gated runs.

Outputs (per gate): period vs density (per lt) and vs lt (per A), replicas
averaged, GATED; an existence map over (A,lt); and log-log slope fits
  d log(period)/d log(A)   and   d log(period)/d log(lt)
in the gated regime. slope ~ 0 -> steric; slope != 0 -> interaction-set.
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
    if lt in ("0", "None", ""):
        return 0.0
    if "e" in lt:
        m, e = lt.split("e")
        return float(m) * 10.0 ** float(e)
    return float(lt)


def run_metrics(path):
    d = np.load(path, allow_pickle=True)
    nx = np.asarray(d["ts_nx_star"])
    npart = np.asarray(d["ts_n_particles"]) if "ts_n_particles" in d.files else np.full_like(nx, 999)
    Lx = float(d["Lx"])
    tail = slice(len(nx) // 2, None)          # steady window
    nx_t = nx[tail]
    good = nx_t[nx_t > 0]
    if good.size == 0:
        modal, stab, period = 0, 0.0, np.nan
    else:
        modal = int(np.bincount(good).argmax())
        stab = float(np.mean(np.abs(good - modal) <= 1))
        period = Lx / modal if modal > 0 else np.nan
    return {
        "A": int(d["n_actins"]), "lt": str(d["lt"]), "lt_val": lt_to_float(str(d["lt"])),
        "rep": (re.search(r"_(r\d+)_", str(d["fname"])) or [None, "r??"])[1]
               if re.search(r"_(r\d+)_", str(d["fname"])) else "r??",
        "modal_nx": modal, "nx_stability": stab, "n_cb2": float(np.median(npart[tail])),
        "period": period, "phi_final": float(np.mean(np.asarray(d["ts_phi"])[tail])),
        "smectic": float(d["smectic_fraction"]),
    }


def loglog_slope(xs, ys):
    xs, ys = np.asarray(xs, float), np.asarray(ys, float)
    ok = (xs > 0) & (ys > 0) & np.isfinite(xs) & np.isfinite(ys)
    if ok.sum() < 3:
        return np.nan, np.nan
    lx, ly = np.log(xs[ok]), np.log(ys[ok])
    A = np.vstack([lx, np.ones_like(lx)]).T
    slope, intc = np.linalg.lstsq(A, ly, rcond=None)[0]
    resid = ly - (slope * lx + intc)
    r2 = 1 - np.sum(resid**2) / (np.sum((ly - ly.mean())**2) + 1e-30)
    return float(slope), float(r2)


def gated_mean(recs, key_x, x_vals, fixed_key, fixed_val):
    """mean period +/- std over gated reps, per x_val, at fixed other control."""
    xs, ys, es, ns = [], [], [], []
    for xv in x_vals:
        vals = [r["period"] for r in recs
                if r[key_x] == xv and r[fixed_key] == fixed_val
                and r["_gated"] and np.isfinite(r["period"])]
        if vals:
            xs.append(xv); ys.append(np.mean(vals)); es.append(np.std(vals)); ns.append(len(vals))
    return xs, ys, es, ns


def main(npz_dir, out_dir, stab_min, ncb2_min):
    recs = [run_metrics(p) for p in sorted(glob.glob(os.path.join(npz_dir, "cg_*.npz")))]
    for r in recs:
        r["_gated"] = (r["nx_stability"] >= stab_min) and (r["n_cb2"] >= ncb2_min)
    As = sorted({r["A"] for r in recs})
    lts = sorted({r["lt_val"] for r in recs})
    lt_label = {r["lt_val"]: r["lt"] for r in recs}
    n_gated = sum(r["_gated"] for r in recs)
    print(f"{len(recs)} runs, {n_gated} pass existence gate "
          f"(nx_stability>={stab_min}, n_cb2>={ncb2_min})")

    # ---- existence map (fraction of reps gated) over (A, lt) ----
    fig, ax = plt.subplots(figsize=(7, 6))
    G = np.full((len(lts), len(As)), np.nan)
    for i, lv in enumerate(lts):
        for j, A in enumerate(As):
            grp = [r for r in recs if r["A"] == A and r["lt_val"] == lv]
            if grp:
                G[i, j] = np.mean([r["_gated"] for r in grp])
    im = ax.imshow(G, aspect="auto", origin="lower", cmap="viridis", vmin=0, vmax=1)
    ax.set_xticks(range(len(As))); ax.set_xticklabels(As)
    ax.set_yticks(range(len(lts))); ax.set_yticklabels([lt_label[v] for v in lts], fontsize=7)
    ax.set_xlabel("actin count A"); ax.set_ylabel("catch-bond strength lt")
    ax.set_title("PATTERN EXISTENCE: fraction of reps passing gate")
    fig.colorbar(im, ax=ax)
    fig.tight_layout(); fig.savefig(os.path.join(out_dir, "fig_existence_map.png"), dpi=150)
    plt.close(fig)

    # ---- period vs density (per lt), GATED ----
    fig, ax = plt.subplots(figsize=(7.5, 5))
    cmap = plt.cm.viridis(np.linspace(0, 1, len(lts)))
    slopes_A = {}
    for i, lv in enumerate(lts):
        xs, ys, es, ns = gated_mean(recs, "A", As, "lt_val", lv)
        if len(xs) >= 3:
            s, r2 = loglog_slope(xs, ys); slopes_A[lt_label[lv]] = (s, r2)
            ax.errorbar(xs, ys, yerr=es, marker="o", color=cmap[i],
                        label=f"lt={lt_label[lv]} (slope={s:+.2f})", capsize=2, lw=1)
    ax.set_xlabel("actin count A"); ax.set_ylabel("band period Lx/nx (GATED)")
    ax.set_title("period vs density, where patterns exist\nslope~0 steric | slope<0 interaction-set")
    ax.legend(fontsize=6, ncol=2)
    fig.tight_layout(); fig.savefig(os.path.join(out_dir, "fig_period_vs_density_gated.png"), dpi=150)
    plt.close(fig)

    # ---- period vs lt (per A), GATED ----
    fig, ax = plt.subplots(figsize=(7.5, 5))
    cmap = plt.cm.plasma(np.linspace(0, 1, len(As)))
    slopes_lt = {}
    for k, A in enumerate(As):
        xs, ys, es, ns = gated_mean(recs, "lt_val", lts, "A", A)
        xs2 = [x for x in xs if x > 0]; ys2 = [y for x, y in zip(xs, ys) if x > 0]
        es2 = [e for x, e in zip(xs, es) if x > 0]
        if len(xs2) >= 3:
            s, r2 = loglog_slope(xs2, ys2); slopes_lt[A] = (s, r2)
            ax.errorbar([x + 1e-7 for x in xs2], ys2, yerr=es2, marker="o", color=cmap[k],
                        label=f"A{A} (slope={s:+.2f})", capsize=2, lw=1)
    ax.set_xscale("log")
    ax.set_xlabel("catch-bond strength lt"); ax.set_ylabel("band period Lx/nx (GATED)")
    ax.set_title("period vs catch-bond strength, where patterns exist\n"
                 "slope~0 lt only ENABLES | slope!=0 lt TUNES the scale")
    ax.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(os.path.join(out_dir, "fig_period_vs_lt_gated.png"), dpi=150)
    plt.close(fig)

    # ---- verdict ----
    print("\n-- log-log slopes d log(period)/d log(A) per lt (gated) --")
    for lt, (s, r2) in sorted(slopes_A.items(), key=lambda kv: lt_to_float(kv[0])):
        print(f"  lt={lt:8s} slope={s:+.2f}  (fit R2={r2:.2f})")
    print("-- log-log slopes d log(period)/d log(lt) per A (gated) --")
    for A, (s, r2) in sorted(slopes_lt.items()):
        print(f"  A{A:<5d} slope={s:+.2f}  (fit R2={r2:.2f})")
    med_sA = np.nanmedian([s for s, _ in slopes_A.values()]) if slopes_A else np.nan
    med_slt = np.nanmedian([s for s, _ in slopes_lt.values()]) if slopes_lt else np.nan
    print(f"\nVERDICT (gated regime): median slope vs density = {med_sA:+.2f}, "
          f"vs lt = {med_slt:+.2f}")
    print("  |slope|<~0.1 => STERIC (flat); |slope|>~0.2 => interaction-set (emergent).")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("npz_dir")
    ap.add_argument("--out_dir", default=None)
    ap.add_argument("--stab_min", type=float, default=0.5, help="min nx-stability to count as a pattern")
    ap.add_argument("--ncb2_min", type=float, default=100.0, help="min median cb2 count to trust the mode")
    args = ap.parse_args()
    main(args.npz_dir, args.out_dir or args.npz_dir, args.stab_min, args.ncb2_min)
