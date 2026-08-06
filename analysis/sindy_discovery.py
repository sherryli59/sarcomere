"""
Data-driven PDE discovery (strong-form STLS / SINDy) for the coarse-grained
sarcomere fields {rho, p, sigma, m_myo}. The field-theory structure is the
OUTPUT: a rich unbiased library spans plausible terms; sparse regression selects.

Runs as a script or opens as a notebook (cells delimited by `# %%`).

Reads cg_*.npz from coarse_grain_fields.py. For quality temporal derivatives you
want ~20-30 snapshots concentrated in a window with real dynamics (assembly
transient), i.e. regenerate with:
    coarse_grain_fields.py <files> --n_snapshots 30 --start_frame 300 --smooth 0.05
With the default 6 snapshots this still runs, but is LOW-TEMPORAL-RESOLUTION
(only a few central-difference points) -- a proof of concept.

CAVEATS baked into the reporting:
 - Strong form: spatial derivs are spectral; the 4th derivative (Swift-Hohenberg
   term) amplifies high-k binning noise ~k^4. Smoothing (--smooth) is essential;
   a weak-form inner product is the robust upgrade (build_library is written so
   a test-function projection can replace pointwise evaluation).
 - Pure steady state has d_t u ~ 0 -> weak identifiability. Prefer a window that
   includes the transient.
"""
# %% imports & config -------------------------------------------------------
import os
import re
import glob
import json
import argparse
import numpy as np

DT_PER_FRAME = 1e-3            # data0430: save_every=100 * dt=1e-5 (confirm per dataset)
FIELD_KEYS = {"rho": "snap_rho", "p": "snap_p",
              "sigma": "snap_sigma", "m": "snap_m_myo",
              "rho_cb2": "snap_rho_cb2"}   # loaded always; used as a library field on demand
TARGETS = ["rho", "p", "sigma", "m"]
THRESHOLD_LEVELS = [1e-3, 1e-2, 1e-1]   # x median|c| ; report stability across these


# %% spectral derivative & lt parsing ---------------------------------------
def spec_deriv(u, Lx, order):
    """d^order u / dx^order for a periodic 1D field via FFT."""
    n = u.shape[-1]
    k = 2.0 * np.pi * np.fft.rfftfreq(n, d=Lx / n)
    return np.fft.irfft((1j * k) ** order * np.fft.rfft(u), n=n)


def lt_to_float(lt):
    lt = str(lt).strip()
    if lt in ("0", "None", ""):
        return 0.0
    if "e" in lt:
        m, e = lt.split("e")
        return float(m) * 10.0 ** float(e)     # exponent may be non-integer
    return float(lt)


def rep_from_name(name):
    m = re.search(r"_(r\d+)_", str(name))
    return m.group(1) if m else "r??"


# %% load one run -----------------------------------------------------------
def load_run(path):
    d = np.load(path, allow_pickle=True)
    fields = {k: np.asarray(d[key], dtype=float) for k, key in FIELD_KEYS.items()}
    snap_frames = np.asarray(d["snap_frames"], dtype=float)
    return {
        "fields": fields,                       # each (n_snap, nbins)
        "t": snap_frames * DT_PER_FRAME,        # snapshot times
        "Lx": float(d["Lx"]),
        "A": int(d["n_actins"]),
        "lt": str(d["lt"]),
        "lt_val": lt_to_float(str(d["lt"])),
        "rep": rep_from_name(str(d["fname"])),
        "n_snap": snap_frames.shape[0],
    }


# %% library builder --------------------------------------------------------
def build_terms(fld):
    """Given a dict of the 4 fields at ONE snapshot (each shape (nbins,)) plus a
    box length, return an ordered dict {label: column(nbins)}. Written so a weak
    form can later replace pointwise columns with test-function inner products."""
    r, p, s, m, Lx = fld["rho"], fld["p"], fld["sigma"], fld["m"], fld["Lx"]
    dx = lambda u, o: spec_deriv(u, Lx, o)
    r_x, r_xx, r_xxxx = dx(r, 1), dx(r, 2), dx(r, 4)
    p_x, p_xx = dx(p, 1), dx(p, 2)
    s_x, s_xx = dx(s, 1), dx(s, 2)
    m_x, m_xx = dx(m, 1), dx(m, 2)
    one = np.ones_like(r)
    return {
        # zeroth-order linear
        "rho": r, "p": p, "sigma": s, "m": m,
        # first-order linear
        "rho_x": r_x, "p_x": p_x, "sigma_x": s_x, "m_x": m_x,
        # second-order linear
        "rho_xx": r_xx, "p_xx": p_xx, "sigma_xx": s_xx, "m_xx": m_xx,
        # fourth-order (Swift-Hohenberg / PFC)
        "rho_xxxx": r_xxxx,
        # quadratic densities
        "rho^2": r * r, "rho*p": r * p, "rho*sigma": r * s,
        "p^2": p * p, "p*sigma": p * s, "sigma^2": s * s,
        # advective / flux
        "rho*rho_x": r * r_x, "p*p_x": p * p_x, "rho*p_x": r * p_x, "p*rho_x": p * r_x,
        # mixed gradient
        "rho*sigma_xx": r * s_xx, "sigma*rho_xx": s * r_xx,
        "rho*p_xx": r * p_xx, "p*rho_xx": p * r_xx,
        # motor-density
        "rho*m_xx": r * m_xx, "m*rho_xx": m * r_xx, "rho*m": r * m,
        # bias
        "1": one,
    }


LIBRARY_LABELS = list(build_terms(
    {k: np.zeros(8) for k in ("rho", "p", "sigma", "m")} | {"Lx": 1.0}).keys())


def assemble(runs, interior=None):
    """Stack library rows and per-target d_t over all interior snapshots and all
    runs in the group. Returns Theta (Npts, Nterms), targets {field: b}, labels."""
    rows, dudt = [], {t: [] for t in TARGETS}
    for run in runs:
        F, tt, Lx, ns = run["fields"], run["t"], run["Lx"], run["n_snap"]
        idxs = interior if interior is not None else range(1, ns - 1)
        for s in idxs:
            if s <= 0 or s >= ns - 1:
                continue                        # need both-sided neighbour
            fld = {k: F[k][s] for k in ("rho", "p", "sigma", "m")}
            fld["Lx"] = Lx
            terms = build_terms(fld)
            rows.append(np.stack([terms[l] for l in LIBRARY_LABELS], axis=1))  # (nbins,Nterms)
            dt = tt[s + 1] - tt[s - 1]
            for t in TARGETS:
                dudt[t].append((F[t][s + 1] - F[t][s - 1]) / dt)
    Theta = np.concatenate(rows, axis=0)
    b = {t: np.concatenate(dudt[t]) for t in TARGETS}
    return Theta, b, LIBRARY_LABELS


# %% STLS solver ------------------------------------------------------------
def stls(Theta, b, rel_threshold, max_iter=20):
    """Sequentially thresholded least squares with column normalization.
    rel_threshold is a fraction of the median |normalized coefficient|."""
    scales = np.linalg.norm(Theta, axis=0)
    scales[scales == 0] = 1.0
    Tn = Theta / scales
    c = np.linalg.lstsq(Tn, b, rcond=None)[0]
    # threshold relative to the LARGEST normalized coefficient (max-relative is
    # the standard sparsifying choice; median-relative barely prunes here).
    thr = rel_threshold * np.max(np.abs(c)) if np.any(c != 0) else 0.0
    keep = np.abs(c) >= thr
    for _ in range(max_iter):
        c = np.zeros(Theta.shape[1])
        if keep.sum() == 0:
            break
        c[keep] = np.linalg.lstsq(Tn[:, keep], b, rcond=None)[0]
        new_keep = np.abs(c) >= thr
        if np.array_equal(new_keep, keep):
            break
        keep = new_keep
    pred = Tn @ c
    ss_res = np.sum((b - pred) ** 2)
    ss_tot = np.sum((b - b.mean()) ** 2) + 1e-30
    r2 = 1.0 - ss_res / ss_tot
    return c / scales, keep, float(r2)          # coef in ORIGINAL units


# %% discovery over one (A, lt) group ---------------------------------------
def discover_group(runs):
    Theta, b, labels = assemble(runs)
    out = {"A": runs[0]["A"], "lt": runs[0]["lt"], "lt_val": runs[0]["lt_val"],
           "n_reps": len(runs), "n_points": Theta.shape[0], "equations": {}}
    for t in TARGETS:
        per_thr = {}
        for lvl in THRESHOLD_LEVELS:
            coef, keep, r2 = stls(Theta, b[t], lvl)
            per_thr[lvl] = {"terms": {labels[i]: float(coef[i])
                                      for i in np.flatnonzero(keep)}, "r2": r2}
        # structural terms = survive ALL thresholds
        survived = set.intersection(*[set(per_thr[l]["terms"]) for l in THRESHOLD_LEVELS]) \
            if per_thr else set()
        out["equations"][t] = {"per_threshold": per_thr, "structural_terms": sorted(survived)}
    return out


def crossval_r2(runs, target, rel_threshold=1e-2):
    """Fit on all but one rep, score R^2 on the held-out rep."""
    if len(runs) < 2:
        return None
    Th_tr, b_tr, labels = assemble(runs[:-1])
    coef, keep, _ = stls(Th_tr, b_tr[target], rel_threshold)
    Th_te, b_te, _ = assemble(runs[-1:])
    pred = Th_te @ coef
    ss_res = np.sum((b_te[target] - pred) ** 2)
    ss_tot = np.sum((b_te[target] - b_te[target].mean()) ** 2) + 1e-30
    return float(1.0 - ss_res / ss_tot)


# %% human-readable equation ------------------------------------------------
def format_equation(target, term_dict, r2):
    if not term_dict:
        return f"  d_t {target} = 0   (no terms; R2={r2:.2f})"
    parts = [f"{c:+.3g}*{lbl}" for lbl, c in sorted(term_dict.items(),
                                                    key=lambda kv: -abs(kv[1]))]
    return f"  d_t {target} = " + " ".join(parts) + f"   (R2={r2:.2f})"


# %% main -------------------------------------------------------------------
def group_runs(npz_dir, A_filter=None, lt_filter=None):
    runs = [load_run(p) for p in sorted(glob.glob(os.path.join(npz_dir, "cg_*.npz")))]
    groups = {}
    for r in runs:
        if A_filter and r["A"] not in A_filter:
            continue
        if lt_filter and r["lt"] not in lt_filter:
            continue
        groups.setdefault((r["A"], r["lt"]), []).append(r)
    return groups


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("npz_dir")
    ap.add_argument("--out_dir", default=None)
    ap.add_argument("--A", type=int, nargs="+", default=None, help="restrict to these actin counts")
    ap.add_argument("--lt", nargs="+", default=None, help="restrict to these lt labels")
    ap.add_argument("--report_threshold", type=float, default=1e-2,
                    help="which threshold level to print equations for")
    args = ap.parse_args()
    out_dir = args.out_dir or args.npz_dir
    os.makedirs(out_dir, exist_ok=True)

    groups = group_runs(args.npz_dir, args.A, args.lt)
    if not groups:
        raise SystemExit("no matching cg_*.npz")

    all_results = []
    for (A, lt), runs in sorted(groups.items(), key=lambda kv: (kv[0][0], lt_to_float(kv[0][1]))):
        ns = runs[0]["n_snap"]
        res = discover_group(runs)
        all_results.append(res)
        print(f"\n=== A{A} lt{lt}  ({len(runs)} reps, {ns} snaps, "
              f"{res['n_points']} points){'  [LOW-RES: <8 snaps]' if ns < 8 else ''} ===")
        for t in TARGETS:
            eq = res["equations"][t]
            thr_res = eq["per_threshold"][args.report_threshold]
            print(format_equation(t, thr_res["terms"], thr_res["r2"]))
            xv = crossval_r2(runs, t)
            print(f"      structural (all thresholds): {eq['structural_terms']}"
                  + (f"   | held-out R2={xv:.2f}" if xv is not None else ""))

    with open(os.path.join(out_dir, "sindy_results.json"), "w") as f:
        json.dump(all_results, f, indent=2)
    print(f"\nwrote {os.path.join(out_dir, 'sindy_results.json')}")
