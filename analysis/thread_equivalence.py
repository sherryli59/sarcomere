#!/usr/bin/env python3
"""Decide whether a code change is behaviourally equivalent under multithreading.

Bitwise comparison only works at one thread.  Above that, force reductions run
in team order and two identical runs decorrelate, so equivalence has to be
argued statistically.  The construction is a paired null:

  * thermal noise is drawn serially from one GSL stream before the parallel
    region, so with a fixed seed the ONLY source of run-to-run difference is
    reduction order;
  * therefore, if `treat` is equivalent, d(ctrl_a, treat) and d(ctrl_a, ctrl_b)
    are draws from the same distribution, seed by seed;
  * the statistic is the paired log-ratio  r_s = log10 d(ctrl_a, treat)_s
                                                - log10 d(ctrl_a, ctrl_b)_s,
    which is centred on zero under equivalence.  Pairing removes the seed-to-
    seed spread, which is far larger than the effect being looked for.

The pert_* arms are a calibration ladder, not controls: they scale
actin_steric_k by 1+eps.  The smallest eps that separates from the null is the
rig's detection threshold, and it is reported alongside the verdict so that a
"no difference" answer carries a known sensitivity.

Usage:
    thread_equivalence.py [--root DIR] [--out-prefix PATH] [--alpha 0.05]
"""
from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

import h5py
import numpy as np

CONTROL_A = "ctrl_a"
CONTROL_B = "ctrl_b"
FLOOR = 1e-18       # below this a divergence is treated as exactly zero


# --------------------------------------------------------------------------
# loading
# --------------------------------------------------------------------------
def committed(handle: h5py.File, *keys: str) -> int:
    n = int(handle["meta/committed_frames"][0])
    return min([n] + [handle[k].shape[0] for k in keys if k in handle])


def load(path: Path) -> dict:
    with h5py.File(path, "r", locking=False) as h:
        meta = h["meta"].attrs
        needed = ["actin/center", "actin/direction", "actin/cb_status",
                  "myosin/center"]
        missing = [k for k in needed if k not in h]
        if missing:
            raise KeyError(f"missing {missing}")
        n = committed(h, *needed, "catch_bond/local_link_frame_offsets")
        if n < 2:
            raise KeyError(f"only {n} committed frames")
        out = {
            "box": np.array([meta["box_x"], meta["box_y"], meta["box_z"]], float),
            "frames": n,
            "center": np.asarray(h["actin/center"][:n], float),
            "direction": np.asarray(h["actin/direction"][:n], float),
            "cb_status": np.asarray(h["actin/cb_status"][:n]),
            "myosin": np.asarray(h["myosin/center"][:n], float),
        }
        key = "catch_bond/local_link_frame_offsets"
        if key in h and h[key].shape[0] >= n:
            offsets = h[key]
            out["n_links"] = np.array([int(offsets[f][2]) for f in range(n)], int)
        else:
            out["n_links"] = np.zeros(n, int)
    return out


# --------------------------------------------------------------------------
# divergence between two runs that share a seed
# --------------------------------------------------------------------------
def divergence(a: dict, b: dict) -> dict:
    """Per-frame RMS separation of matched filaments, plus discrete disagreement.

    Actin turnover is off in the rig, so index i is the same filament in both
    runs and the comparison needs no matching step.
    """
    n = min(a["frames"], b["frames"])
    box = a["box"]
    d = a["center"][:n] - b["center"][:n]
    d -= box * np.round(d / box)
    position = np.sqrt((d * d).sum(-1).mean(-1))            # um, per frame

    dot = (a["direction"][:n] * b["direction"][:n]).sum(-1)
    orientation = np.sqrt((np.arccos(np.clip(dot, -1.0, 1.0)) ** 2).mean(-1))

    status = (a["cb_status"][:n] != b["cb_status"][:n]).mean(-1)
    links = np.abs(a["n_links"][:n] - b["n_links"][:n]).astype(float)
    return {"frames": n, "position": position, "orientation": orientation,
            "status": status, "links": links}


def terminal(series: np.ndarray) -> float:
    """Mean over the last quarter of the trajectory, to damp frame-to-frame noise."""
    tail = series[max(1, int(0.75 * len(series))):]
    return float(np.mean(tail)) if len(tail) else float("nan")


# A terminal statistic alone would be blind here.  The model decorrelates in
# ~20 ms, so by the end of a 2000-step run both the null and a genuinely
# different arm have saturated at the same chaotic plateau and their ratio
# returns to 1.  All the discriminating information is in how fast divergence
# grows BEFORE saturation, so the comparison is run frame by frame and the
# verdict is taken from the frames where the null is still small.
def frame_resolved(null_curves: dict, test_curves: dict, metric: str,
                   alpha: float) -> dict:
    seeds = sorted(set(null_curves) & set(test_curves))
    if not seeds:
        return {"frames": 0}
    n = min(min(null_curves[s]["frames"] for s in seeds),
            min(test_curves[s]["frames"] for s in seeds))
    per_frame = []
    for f in range(n):
        null_f = {s: float(null_curves[s][metric][f]) for s in seeds}
        test_f = {s: float(test_curves[s][metric][f]) for s in seeds}
        stats = compare(null_f, test_f)
        stats["frame"] = f
        stats["null_median"] = float(np.median(list(null_f.values())))
        stats["test_median"] = float(np.median(list(test_f.values())))
        per_frame.append(stats)
    flagged = [p["frame"] for p in per_frame
               if p["n_usable"] >= 6 and not math.isnan(p["sign_p"])
               and p["sign_p"] <= alpha and p["median_log10_ratio"] > 0]
    return {"frames": n, "per_frame": per_frame,
            "first_flagged_frame": flagged[0] if flagged else None,
            "n_flagged_frames": len(flagged)}


def windowed(null_curves: dict, test_curves: dict, metric: str,
             last_frame: int) -> dict:
    """One statistic per seed, so the verdict costs exactly one test.

    Testing every frame separately invites multiplicity: at 81 frames and
    alpha=0.05 about four frames flag by chance, so "one flagged frame" is
    noise, not detection.  Instead each seed contributes the mean log10 ratio
    over the pre-saturation window, and a single sign test decides.
    """
    seeds = sorted(set(null_curves) & set(test_curves))
    per_seed = {}
    for s in seeds:
        n = min(null_curves[s]["frames"], test_curves[s]["frames"], last_frame + 1)
        # Frame 0 is included: it follows the pre-phase, so it already carries
        # 500 steps of amplification and is the single most informative frame.
        a = np.asarray(null_curves[s][metric][0:n], float)
        b = np.asarray(test_curves[s][metric][0:n], float)
        keep = (a > FLOOR) & (b > FLOOR)
        if not keep.any():
            continue
        per_seed[s] = float(np.mean(np.log10(b[keep]) - np.log10(a[keep])))
    values = list(per_seed.values())
    positive = sum(1 for v in values if v > 0)
    return {"n_seeds": len(values),
            "median_log10_ratio": float(np.median(values)) if values else float("nan"),
            "iqr_log10_ratio": (float(np.percentile(values, 75) -
                                      np.percentile(values, 25))
                                if len(values) > 3 else float("nan")),
            "n_positive": positive,
            "sign_p": sign_test_p(positive, len(values)),
            "per_seed": per_seed}


def informative_window(null_curves: dict, metric: str,
                       fraction: float = 0.10) -> int | None:
    """Last frame at which the null is still well below its decorrelated level.

    Measured behaviour: the null reaches 90% of its final level only near the
    end of the run, so a "not yet saturated" cutoff admits almost every frame
    and buries the signal.  The separation actually lives in the steep early
    growth -- at 2000 steps the perturbation ladder is ordered monotonically at
    the first saved frame and has merged into the null by the fourth.  So the
    window ends once the null has climbed to `fraction` of its final value,
    beyond which both runs are simply following the same chaotic envelope.
    """
    seeds = sorted(null_curves)
    if not seeds:
        return None
    n = min(null_curves[s]["frames"] for s in seeds)
    med = np.median(np.vstack([null_curves[s][metric][:n] for s in seeds]), axis=0)
    if not np.isfinite(med).any() or med[-1] <= FLOOR:
        return None
    beyond = np.flatnonzero(med >= fraction * med[-1])
    last = int(beyond[0]) if beyond.size else n - 1
    return max(0, last)


# --------------------------------------------------------------------------
# the paired test
# --------------------------------------------------------------------------
def sign_test_p(positive: int, total: int) -> float:
    """Two-sided exact binomial p under H0: median log-ratio is zero."""
    if total == 0:
        return float("nan")
    k = min(positive, total - positive)
    tail = sum(math.comb(total, i) for i in range(0, k + 1)) / (2.0 ** total)
    return min(1.0, 2.0 * tail)


def compare(null: dict[str, float], test: dict[str, float]) -> dict:
    """Paired log10 ratio of test divergence to null divergence, seed by seed."""
    seeds = sorted(set(null) & set(test))
    ratios, degenerate = [], 0
    for s in seeds:
        a, b = null[s], test[s]
        if a <= FLOOR and b <= FLOOR:
            degenerate += 1
            continue
        if a <= FLOOR:                      # null is zero but test is not
            ratios.append(float("inf"))
            continue
        ratios.append(math.log10(max(b, FLOOR) / a))
    finite = [r for r in ratios if math.isfinite(r)]
    positive = sum(1 for r in ratios if r > 0)
    return {
        "n_seeds": len(seeds),
        "n_usable": len(ratios),
        "n_degenerate": degenerate,
        "median_log10_ratio": float(np.median(finite)) if finite else float("nan"),
        "iqr_log10_ratio": (float(np.percentile(finite, 75) -
                                  np.percentile(finite, 25))
                            if len(finite) > 3 else float("nan")),
        "n_positive": positive,
        "sign_p": sign_test_p(positive, len(ratios)),
        "per_seed_log10_ratio": {s: (ratios[i] if i < len(ratios) else None)
                                 for i, s in enumerate(seeds)},
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path,
                    default=Path("/scratch/users/xsli/sarcomere/thread_equiv"))
    ap.add_argument("--out-prefix", type=Path, default=None)
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--window-fraction", type=float, default=0.10,
                    help="window ends once the null reaches this fraction of "
                         "its final divergence")
    ap.add_argument("--metric", default="position",
                    choices=["position", "orientation", "status", "links"])
    args = ap.parse_args()

    runs: dict[tuple[str, str], Path] = {}
    pattern = re.compile(r"^(?P<arm>.+)_seed(?P<seed>\d+)\.h5$")
    for path in sorted((args.root / "data").glob("*.h5")):
        m = pattern.match(path.name)
        if m:
            runs[(m["arm"], m["seed"])] = path

    arms = sorted({a for a, _ in runs})
    seeds = sorted({s for _, s in runs})
    if CONTROL_A not in arms or CONTROL_B not in arms:
        raise SystemExit(f"need {CONTROL_A} and {CONTROL_B}; found {arms}")
    print(f"root {args.root}")
    print(f"{len(runs)} runs, {len(arms)} arms, {len(seeds)} seeds\n")

    cache: dict[tuple[str, str], dict] = {}

    def get(arm: str, seed: str):
        key = (arm, seed)
        if key not in runs:
            return None
        if key not in cache:
            try:
                cache[key] = load(runs[key])
            except (OSError, KeyError) as exc:
                # A run still being written has partial groups; skipping keeps a
                # partial sweep readable instead of aborting on the first one.
                print(f"  skipping {runs[key].name}: {exc}")
                cache[key] = None
        return cache[key]

    # divergence of every arm from ctrl_a, seed by seed
    curves: dict[str, dict[str, dict]] = {}
    for arm in arms:
        if arm == CONTROL_A:
            continue
        per_seed = {}
        for seed in seeds:
            a, b = get(CONTROL_A, seed), get(arm, seed)
            if a is None or b is None:
                continue
            per_seed[seed] = divergence(a, b)
        if per_seed:
            curves[arm] = per_seed

    scalar = {arm: {s: terminal(d[args.metric]) for s, d in per.items()}
              for arm, per in curves.items()}

    null = scalar[CONTROL_B]
    null_values = np.array([v for v in null.values() if v > FLOOR])
    print(f"NULL  (ctrl_a vs ctrl_b), metric = {args.metric}")
    if null_values.size == 0:
        print("  every control pair is identical: the model is deterministic at "
              "this thread count, so use the bitwise test instead.\n")
    else:
        print(f"  median {np.median(null_values):.3e}   "
              f"range [{null_values.min():.3e}, {null_values.max():.3e}]   "
              f"n={null_values.size}\n")

    sat = informative_window(curves[CONTROL_B], args.metric, args.window_fraction)
    n_all = min(c["frames"] for c in curves[CONTROL_B].values())
    print(f"  informative window ends at frame {sat} of {n_all} "
          f"(null reaches {100*args.window_fraction:.0f}% of its decorrelated "
          f"level there)\n" if sat is not None else "  no informative window\n")

    window_end = sat if sat is not None else n_all - 1
    print(f"  verdict window = frames 0..{window_end} "
          f"(one statistic per seed, one test per arm)\n")

    rows, results = [], {}
    for arm in sorted(curves):
        if arm == CONTROL_B:
            continue
        win = windowed(curves[CONTROL_B], curves[arm], args.metric, window_end)
        fr = frame_resolved(curves[CONTROL_B], curves[arm], args.metric, args.alpha)
        different = ((not math.isnan(win["sign_p"])) and
                     win["sign_p"] <= args.alpha and
                     win["median_log10_ratio"] > 0)
        stats = {**win,
                 "terminal_log10_ratio": compare(null, scalar[arm])["median_log10_ratio"],
                 "n_flagged_frames": fr.get("n_flagged_frames", 0),
                 "verdict": "DIFFERENT" if different else "EQUIVALENT"}
        results[arm] = stats
        rows.append((arm, stats))

    print(f"{'arm':12s}{'window log10 ratio':>20s}{'IQR':>9s}{'n+/n':>9s}"
          f"{'sign p':>10s}{'terminal':>10s}   verdict")
    for arm, st in rows:
        print(f"{arm:12s}{st['median_log10_ratio']:>20.4f}"
              f"{st['iqr_log10_ratio']:>9.4f}"
              f"{st['n_positive']:>5d}/{st['n_seeds']:<3d}"
              f"{st['sign_p']:>10.5f}{st['terminal_log10_ratio']:>10.3f}"
              f"   {st['verdict']}")

    # calibration: the smallest deliberate perturbation the rig separates
    ladder = [(arm, st) for arm, st in rows if arm.startswith("pert")]
    detected = [arm for arm, st in ladder if st["verdict"] == "DIFFERENT"]
    print()
    if ladder:
        print("CALIBRATION (deliberate actin_steric_k perturbations)")
        for arm, st in ladder:
            print(f"  {arm:10s} {st['verdict']:10s} "
                  f"window log10 ratio {st['median_log10_ratio']:+.4f}  "
                  f"sign p {st['sign_p']:.5f}")
        print(f"  -> rig detects: {', '.join(detected) if detected else 'NONE'}")
        if not detected:
            print("  -> WARNING: the rig failed to detect any deliberate change, "
                  "so a null result from it means nothing.  Raise the seed count "
                  "or the step count before trusting any verdict.")
    if "treat" in results:
        print(f"\nTREAT VERDICT: {results['treat']['verdict']}")
        print("  Interpret only against the calibration line above: an "
              "EQUIVALENT verdict bounds the change below the smallest "
              "perturbation the rig can see, not below zero.")

    if args.out_prefix:
        payload = {"metric": args.metric, "alpha": args.alpha,
                   "saturation_frame": sat,
                   "arms": arms, "seeds": seeds,
                   "null_median": (float(np.median(null_values))
                                   if null_values.size else None),
                   "results": results,
                   "detected_perturbations": detected}
        Path(f"{args.out_prefix}.json").write_text(json.dumps(payload, indent=1))
        _plot(curves, args.metric, f"{args.out_prefix}.png")
        print(f"\nwrote {args.out_prefix}.json and {args.out_prefix}.png")


def _plot(curves, metric, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    order = [CONTROL_B, "treat", "pert1em8", "pert1em6", "pert1em4"]
    colours = {CONTROL_B: "#2a78d6", "treat": "#eb6834", "pert1em8": "#8c8c8c",
               "pert1em6": "#5c5c5c", "pert1em4": "#111111"}
    fig, ax = plt.subplots(figsize=(7.6, 4.8))
    for arm in order:
        if arm not in curves:
            continue
        stack = [c[metric] for c in curves[arm].values()]
        n = min(len(s) for s in stack)
        band = np.vstack([s[:n] for s in stack])
        t = np.arange(n)
        med = np.median(band, axis=0)
        ax.fill_between(t, np.percentile(band, 25, axis=0),
                        np.percentile(band, 75, axis=0),
                        color=colours.get(arm, "0.5"), alpha=0.18, lw=0)
        ax.plot(t, med, color=colours.get(arm, "0.5"), lw=2,
                label=f"{arm} (n={len(stack)})")
    ax.set_yscale("log")
    ax.set_xlabel("saved frame")
    ax.set_ylabel(f"RMS {metric} difference from ctrl_a")
    ax.set_title("Thread-equivalence rig: divergence against the paired null\n"
                 "band = seed IQR; grey = deliberate perturbation ladder",
                 fontsize=10.5)
    ax.grid(alpha=0.25, lw=0.6)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, fontsize=8.5)
    fig.tight_layout()
    fig.savefig(path, dpi=160)


if __name__ == "__main__":
    main()
