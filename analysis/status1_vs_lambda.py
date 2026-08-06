"""Count status-1 catch-bond engagement vs the catch-bond strength lambda (= lt).

Motivation
----------
In the ORIGINAL (pre-split) kinetics every bond with status >= 1 received the
load-dependent lifetime bonus

    lifetime = base_lifetime + [directional_base_lifetime if status==2]
                             + lambda * f_load          <-- status 1 AND 2

so a merely "loaded" status-1 contact could be stabilized indefinitely.  The
new kinetics (cpp/include/catch_bond_kinetics.h) gates that last term behind
`stabilize_status1`, defaulting to status-2-only.  This script quantifies how
much status-1 population the OLD runs actually carried as a function of lambda.

What is measurable (IMPORTANT)
------------------------------
`/actin/bonds` is NOT a full bond list.  `Sarcomere::_extract_bonded_pairs`
(cpp/src/sarcomere.cpp) filters to status-2 pairs at write time:

    if (actin_actin_bonds[i][j] == 1 &&
        actin_actin_status[i][j] == 2 &&
        actin.cb_status[i] == 2 && actin.cb_status[j] == 2)

Status-1 pairs are therefore never serialized, and a per-BOND status-1 count
cannot be recovered from the trajectories.  The only surviving trace of status
1 is the per-ACTIN field `/actin/cb_status`, which `_cb_decide` assigns as the
MAX over that actin's accepted bonds:

    cb_status[i] == 0  ->  actin i has no accepted contact
    cb_status[i] == 1  ->  actin i has >=1 accepted contact and NONE is status 2
    cb_status[i] == 2  ->  actin i participates in >=1 status-2 bond

So `n1 = count(cb_status == 1)` counts ACTINS ENGAGED ONLY IN STATUS-1
CONTACTS.  It is a lower bound on status-1 participation: an actin holding both
a status-2 and a status-1 bond registers as 2 and is invisible here.  Read the
reported n1 as "actins whose only engagement is a status-1 contact", not as a
bond count.

Outputs
-------
  <prefix>_per_traj.csv    one row per trajectory (steady-state + peak stats)
  <prefix>_by_lambda.csv   aggregated over replicas, grouped by (A, lt)
  <prefix>_timeseries.npz  per-trajectory n0/n1/n2 time series
  <prefix>_n1_vs_lambda.png / _fractions_vs_lambda.png
"""

import os

os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

import argparse
import csv
import glob
import re
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed

import h5py
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

_A_RE = re.compile(r"A(\d+)_")
_LT_RE = re.compile(r"_lt([^_]+?)_r\d+_")
_REP_RE = re.compile(r"_(r\d+)_")


# --------------------------------------------------------------------------- #
# filename / metadata helpers (conventions shared with aggregate_phase.py)
# --------------------------------------------------------------------------- #
def lt_to_float(lt):
    """lt label -> lambda. Exponents may be non-integer ('1e-2.5')."""
    lt = str(lt).strip()
    if lt in ("0", "lt0", "None", ""):
        return 0.0
    if "e" in lt:
        mant, exp = lt.split("e")
        return float(mant) * 10.0 ** float(exp)
    return float(lt)


def parse_name(path):
    name = os.path.basename(path)
    a = _A_RE.search(name)
    lt = _LT_RE.search(name)
    rep = _REP_RE.search(name)
    return {
        "file": name,
        "A": int(a.group(1)) if a else -1,
        "lt": lt.group(1) if lt else "?",
        "lt_val": lt_to_float(lt.group(1)) if lt else float("nan"),
        "rep": rep.group(1) if rep else "r??",
    }


def n_frames(fh):
    """Number of committed frames (meta/committed_frames if present)."""
    if "/meta/committed_frames" in fh:
        c = int(np.asarray(fh["/meta/committed_frames"])[0])
        if c > 0:
            return c
    return fh["/actin/cb_status"].shape[0]


# --------------------------------------------------------------------------- #
# per-trajectory measurement
# --------------------------------------------------------------------------- #
def process(path, stride, tail_frac):
    meta = parse_name(path)
    try:
        with h5py.File(path, "r") as fh:
            if "/actin/cb_status" not in fh:
                return {**meta, "error": "no /actin/cb_status"}
            ds = fh["/actin/cb_status"]
            nf = min(n_frames(fh), ds.shape[0])
            if nf <= 0:
                return {**meta, "error": "no committed frames"}
            # cb_status is {T, Na, 1} int -> a few tens of MB, safe to slurp.
            st = np.asarray(ds[0:nf:stride]).reshape(-1, ds.shape[1])
    except Exception as exc:  # noqa: BLE001 - report, do not abort the sweep
        return {**meta, "error": f"{type(exc).__name__}: {exc}"}

    n_actins = st.shape[1]
    n0 = (st == 0).sum(axis=1)
    n1 = (st == 1).sum(axis=1)
    n2 = (st == 2).sum(axis=1)

    tail = max(1, int(round(tail_frac * len(n1))))
    s1, s2, s0 = n1[-tail:], n2[-tail:], n0[-tail:]

    return {
        **meta,
        "error": "",
        "n_actins": n_actins,
        "n_frames": nf,
        "n_sampled": len(n1),
        # steady state = last tail_frac of the trajectory
        "n1_mean": float(s1.mean()),
        "n1_std": float(s1.std()),
        "n2_mean": float(s2.mean()),
        "n2_std": float(s2.std()),
        "n0_mean": float(s0.mean()),
        "f1_mean": float(s1.mean() / n_actins),
        "f2_mean": float(s2.mean() / n_actins),
        # transient behaviour: status 1 is the putative nucleation precursor
        "n1_peak": float(n1.max()),
        "n1_peak_frame": int(np.argmax(n1) * stride),
        "n1_traj_mean": float(n1.mean()),
        "_ts": (n0.astype(np.int32), n1.astype(np.int32), n2.astype(np.int32)),
    }


# --------------------------------------------------------------------------- #
# aggregation over replicas
# --------------------------------------------------------------------------- #
def aggregate(rows):
    groups = defaultdict(list)
    for r in rows:
        groups[(r["A"], r["lt"])].append(r)

    out = []
    for (A, lt), rs in groups.items():
        get = lambda k: np.array([r[k] for r in rs], dtype=float)  # noqa: E731
        out.append(
            {
                "A": A,
                "lt": lt,
                "lt_val": lt_to_float(lt),
                "n_reps": len(rs),
                "n_actins": rs[0]["n_actins"],
                "n1_mean": float(get("n1_mean").mean()),
                "n1_sem": float(get("n1_mean").std() / max(1, np.sqrt(len(rs)))),
                "n2_mean": float(get("n2_mean").mean()),
                "n2_sem": float(get("n2_mean").std() / max(1, np.sqrt(len(rs)))),
                "n0_mean": float(get("n0_mean").mean()),
                "f1_mean": float(get("f1_mean").mean()),
                "f2_mean": float(get("f2_mean").mean()),
                "n1_peak_mean": float(get("n1_peak").mean()),
                "n1_traj_mean": float(get("n1_traj_mean").mean()),
            }
        )
    out.sort(key=lambda d: (d["A"], d["lt_val"]))
    return out


def write_csv(rows, path, fields):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"wrote {path}  ({len(rows)} rows)")


# --------------------------------------------------------------------------- #
# plots
# --------------------------------------------------------------------------- #
def _xs(agg):
    # symlog with a small offset so lt=0 is plottable alongside the decades
    return [d["lt_val"] + 1e-7 for d in agg]


def plot_counts(agg, out):
    by_A = defaultdict(list)
    for d in agg:
        by_A[d["A"]].append(d)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), sharex=True)
    cmap = plt.get_cmap("viridis")
    keys = sorted(by_A)
    for k, A in enumerate(keys):
        ds = sorted(by_A[A], key=lambda d: d["lt_val"])
        c = cmap(k / max(1, len(keys) - 1))
        axes[0].errorbar(
            _xs(ds), [d["n1_mean"] for d in ds], yerr=[d["n1_sem"] for d in ds],
            marker="o", ms=4, color=c, label=f"A={A}", capsize=2,
        )
        axes[1].errorbar(
            _xs(ds), [d["n2_mean"] for d in ds], yerr=[d["n2_sem"] for d in ds],
            marker="s", ms=4, color=c, label=f"A={A}", capsize=2,
        )

    axes[0].set_ylabel("# actins with cb_status == 1")
    axes[0].set_title("Status-1-only actins (steady state)")
    axes[1].set_ylabel("# actins with cb_status == 2")
    axes[1].set_title("Status-2 actins (steady state)")
    for ax in axes:
        ax.set_xscale("symlog", linthresh=1e-5)
        ax.set_xlabel(r"$\lambda$  (lifetime_coeff, lt)")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    print(f"wrote {out}")


def plot_fractions(agg, out):
    by_A = defaultdict(list)
    for d in agg:
        by_A[d["A"]].append(d)

    fig, ax = plt.subplots(figsize=(7, 4.5))
    cmap = plt.get_cmap("viridis")
    keys = sorted(by_A)
    for k, A in enumerate(keys):
        ds = sorted(by_A[A], key=lambda d: d["lt_val"])
        c = cmap(k / max(1, len(keys) - 1))
        ax.plot(_xs(ds), [d["f1_mean"] for d in ds], marker="o", ms=4,
                color=c, ls="-", label=f"A={A} status 1")
        ax.plot(_xs(ds), [d["f2_mean"] for d in ds], marker="s", ms=4,
                color=c, ls="--", alpha=0.6, label=f"A={A} status 2")
    ax.set_xscale("symlog", linthresh=1e-5)
    ax.set_xlabel(r"$\lambda$  (lifetime_coeff, lt)")
    ax.set_ylabel("fraction of actins")
    ax.set_title("Status-1-only (solid) vs status-2 (dashed) fraction")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=6, ncol=2)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    print(f"wrote {out}")


# --------------------------------------------------------------------------- #
def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--data_dir", default="/scratch/users/xsli/sarcomere/data0430")
    p.add_argument("--glob", default="*.h5",
                   help="filename glob within data_dir")
    p.add_argument("--stride", type=int, default=1,
                   help="frame stride (cb_status is small; 1 is affordable)")
    p.add_argument("--tail_frac", type=float, default=0.25,
                   help="trailing fraction of frames treated as steady state")
    p.add_argument("--workers", type=int,
                   default=int(os.environ.get("SLURM_CPUS_PER_TASK", 4)))
    p.add_argument("--prefix", default="status1_vs_lambda")
    args = p.parse_args()

    paths = sorted(glob.glob(os.path.join(args.data_dir, args.glob)))
    paths = [p_ for p_ in paths if p_.endswith(".h5")]
    if not paths:
        raise SystemExit(f"no .h5 matched {args.data_dir}/{args.glob}")
    print(f"{len(paths)} trajectories from {args.data_dir}", flush=True)

    rows, ts_store, failures = [], {}, []
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(process, p_, args.stride, args.tail_frac): p_
                for p_ in paths}
        for i, fut in enumerate(as_completed(futs), 1):
            r = fut.result()
            if r.get("error"):
                failures.append((r["file"], r["error"]))
                print(f"[{i}/{len(paths)}] SKIP {r['file']}: {r['error']}",
                      flush=True)
                continue
            ts = r.pop("_ts")
            # npz keys become filenames inside the zip: keep them separator-free
            stem = r["file"][:-3] if r["file"].endswith(".h5") else r["file"]
            ts_store[stem + "__n0"] = ts[0]
            ts_store[stem + "__n1"] = ts[1]
            ts_store[stem + "__n2"] = ts[2]
            rows.append(r)
            print(f"[{i}/{len(paths)}] {r['file']}  lt={r['lt']} "
                  f"n1={r['n1_mean']:.1f} n2={r['n2_mean']:.1f} "
                  f"(peak n1={r['n1_peak']:.0f} @ frame {r['n1_peak_frame']})",
                  flush=True)

    if not rows:
        raise SystemExit("no trajectory produced data")

    rows.sort(key=lambda r: (r["A"], r["lt_val"], r["rep"]))
    write_csv(rows, f"{args.prefix}_per_traj.csv",
              ["file", "A", "lt", "lt_val", "rep", "n_actins", "n_frames",
               "n_sampled", "n0_mean", "n1_mean", "n1_std", "n2_mean",
               "n2_std", "f1_mean", "f2_mean", "n1_peak", "n1_peak_frame",
               "n1_traj_mean"])

    agg = aggregate(rows)
    write_csv(agg, f"{args.prefix}_by_lambda.csv",
              ["A", "lt", "lt_val", "n_reps", "n_actins", "n0_mean", "n1_mean",
               "n1_sem", "n2_mean", "n2_sem", "f1_mean", "f2_mean",
               "n1_peak_mean", "n1_traj_mean"])

    np.savez_compressed(f"{args.prefix}_timeseries.npz", stride=args.stride,
                        **ts_store)
    print(f"wrote {args.prefix}_timeseries.npz")

    plot_counts(agg, f"{args.prefix}_n1_vs_lambda.png")
    plot_fractions(agg, f"{args.prefix}_fractions_vs_lambda.png")

    print("\n=== status-1-only actins vs lambda (steady state) ===")
    print(f"{'A':>6} {'lt':>8} {'lambda':>10} {'n1':>9} {'n2':>9} "
          f"{'f1':>7} {'f2':>7}")
    for d in agg:
        print(f"{d['A']:>6} {d['lt']:>8} {d['lt_val']:>10.3g} "
              f"{d['n1_mean']:>9.1f} {d['n2_mean']:>9.1f} "
              f"{d['f1_mean']:>7.3f} {d['f2_mean']:>7.3f}")

    if failures:
        print(f"\n{len(failures)} file(s) skipped:")
        for f, e in failures:
            print(f"  {f}: {e}")


if __name__ == "__main__":
    main()
