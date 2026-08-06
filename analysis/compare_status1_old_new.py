"""Compare status-1 engagement between the OLD sweep and the NEW A2000 sweeps.

Why this comparison exists
--------------------------
The old sweep (data0430, cpp @ 30829b0, built -O0/Debug) and the new A2000
sweeps (current working tree, built -O2/RelWithDebInfo) differ by a large mixed
changeset: the status-1/status-2 kinetics gate AND a pile of performance work
(std::map -> struct, std::vector -> std::array, printf removal, -O0 -> -O2).
That is a confound: any old/new difference could be kinetics OR speedups.

The `s1s2` sweep breaks the confound. Its parameters are:
    base_lifetime = 0, directional_base_lifetime = 0, stabilize_status1 = true
which reproduces the OLD kinetics law exactly

    lifetime(status 1) = lifetime(status 2) = lambda * f_load

but running the NEW optimized code.  So:

    s1s2  ~=  old   ->  speedups are physics-neutral; differences are kinetics.
    s1s2  !=  old   ->  the refactor changed the physics.

The `s2only` sweep (stabilize_status1 = false) has lifetime(status 1) = 0
exactly -- status-1 contacts cannot survive a single step.

Matched-window comparison
-------------------------
The runs have very different lengths (old ~4007 frames; s1s2 ~3-825).  A
"last 25%" statistic therefore samples different physical times.  This script
compares a COMMON frame window [--win_lo, --win_hi] across all sweeps and drops
replicas too short to cover it, so old/new are compared like for like.
"""

import argparse
import os
import re
from collections import defaultdict

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

_LT_RE = re.compile(r"_lt([^_]+?)_r\d+_")
_A_RE = re.compile(r"A(\d+)_")
_REP_RE = re.compile(r"_(r\d+)_")


def lt_to_float(lt):
    lt = str(lt).strip()
    if lt in ("0", "lt0", "None", ""):
        return 0.0
    if "e" in lt:
        mant, exp = lt.split("e")
        return float(mant) * 10.0 ** float(exp)
    return float(lt)


def load(npz_path, want_A=None):
    """-> {(A, lt): [n1_series, ...]} keyed by config, one entry per replica."""
    z = np.load(npz_path)
    out = defaultdict(list)
    for key in z.files:
        if not key.endswith("__n1"):
            continue
        stem = key[: -len("__n1")]
        a = _A_RE.search(stem)
        lt = _LT_RE.search(stem + "_")
        if not (a and lt):
            continue
        A = int(a.group(1))
        if want_A is not None and A not in want_A:
            continue
        n1 = np.asarray(z[key], dtype=float)
        n2 = np.asarray(z[stem + "__n2"], dtype=float)
        out[(A, lt.group(1))].append((stem, n1, n2))
    return out


def window_stat(series, lo, hi):
    """Mean over frames [lo,hi) for replicas long enough; None if none qualify."""
    vals1, vals2, used = [], [], 0
    for _stem, n1, n2 in series:
        if len(n1) < hi:
            continue
        vals1.append(n1[lo:hi].mean())
        vals2.append(n2[lo:hi].mean())
        used += 1
    if not used:
        return None
    return float(np.mean(vals1)), float(np.mean(vals2)), used


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--old", default="status1_vs_lambda_old0430_timeseries.npz")
    p.add_argument("--s1s2", default="status1_new_data_a2000_s1s2_lambda_timeseries.npz")
    p.add_argument("--s2only", default="status1_new_data_a2000_status2_lambda_timeseries.npz")
    p.add_argument("--win_lo", type=int, default=200)
    p.add_argument("--win_hi", type=int, default=450)
    p.add_argument("--out", default="compare_status1_old_new.png")
    args = p.parse_args()

    old = load(args.old, want_A={1800, 2100})
    s1s2 = load(args.s1s2)
    s2only = load(args.s2only)

    lts = sorted({lt for (_A, lt) in list(s1s2) + list(s2only)},
                 key=lt_to_float)

    lo, hi = args.win_lo, args.win_hi
    print(f"matched frame window [{lo}, {hi})  "
          f"(replicas shorter than {hi} frames are dropped)\n")
    header = (f"{'lambda':>10} | {'OLD A1800':>18} {'OLD A2100':>18} | "
              f"{'s1s2 A2000':>18} | {'s2only A2000':>18}")
    print(header)
    print("-" * len(header))

    rows = []
    for lt in lts:
        lam = lt_to_float(lt)
        cells, rec = [], {"lt": lt, "lam": lam}
        for label, src, A in (("old1800", old, 1800), ("old2100", old, 2100),
                              ("s1s2", s1s2, 2000), ("s2only", s2only, 2000)):
            st = window_stat(src.get((A, lt), []), lo, hi)
            if st is None:
                cells.append(f"{'--':>18}")
                rec[label] = None
            else:
                n1, n2, k = st
                cells.append(f"{n1:>7.1f}/{n2:<7.1f}({k})")
                rec[label] = (n1, n2, k)
        print(f"{lam:>10.3g} | {cells[0]} {cells[1]} | {cells[2]} | {cells[3]}")
        rows.append(rec)

    print("\ncells are  n1/n2 (n_replicas);  n1 = actins with cb_status==1, "
          "n2 = cb_status==2;  '--' = no replica long enough")

    # ---------------- figure ----------------
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8))
    series = [
        ("OLD A1800 (old code, -O0)", "old1800", "tab:blue", "o", "-"),
        ("OLD A2100 (old code, -O0)", "old2100", "tab:cyan", "^", "-"),
        ("NEW s1s2 A2000 (old kinetics, new code)", "s1s2", "tab:red", "s", "--"),
        ("NEW s2only A2000 (tau1 = 0)", "s2only", "tab:green", "D", ":"),
    ]
    for ax, idx, name in ((axes[0], 0, "n1  (status-1-only actins)"),
                          (axes[1], 1, "n2  (status-2 actins)")):
        for label, key, color, mk, ls in series:
            xs = [r["lam"] + 1e-7 for r in rows if r[key] is not None]
            ys = [r[key][idx] for r in rows if r[key] is not None]
            if xs:
                ax.plot(xs, ys, marker=mk, ls=ls, color=color, ms=5, label=label)
        ax.set_xscale("symlog", linthresh=1e-5)
        ax.set_xlabel(r"$\lambda$ (lifetime_coeff)")
        ax.set_ylabel(name)
        ax.grid(alpha=0.3)
    axes[0].set_title(f"Status-1 population, matched frames [{lo},{hi})")
    axes[1].set_title("Status-2 population (same window)")
    axes[0].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(args.out, dpi=150)
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
