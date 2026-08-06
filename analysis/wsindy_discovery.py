"""
Weak-form PDE discovery (WSINDy) for the CG sarcomere fields {rho, p, sigma, m}.

Fixes the strong-form failure (near-zero, non-generalizing fits): here we NEVER
differentiate the noisy data. Library terms are d_x^alpha[P(u)] (spatial deriv of
a monomial); the weak inner product against a smooth compact test function phi
moves all derivatives onto phi (analytic):

    <phi, d_t u>          = -<d_t phi, u>
    <phi, d_x^alpha P(u)> = (-1)^alpha <d_x^alpha phi, P(u)>

So every integral is (data monomial) x (analytic test-function derivative) -- no
numerical derivative of data, in space or time.

Three fixes fused:
  1. weak form                      -> the inner products above
  2. tight dynamic window           -> regenerate CG with --start_frame/--end_frame
  3. conservation for rho, p        -> restrict their library to alpha >= 1
                                       (RHS must be d_x(flux); no zeroth-order term)

Test functions: separable phi(x,t) = G((x-xq)/Hx) H((t-tq)/Ht), G=H=(1-s^2)^beta
compactly supported on |s|<1; periodic wrap in x. Derivatives of (1-s^2)^beta are
exact polynomials (computed via numpy.polynomial), so d_x^alpha phi is exact.

Reuses load_run / stls / format_equation / crossval from sindy_discovery.
"""
# %% imports ----------------------------------------------------------------
import os
import glob
import json
import argparse
import itertools
import numpy as np
from numpy.polynomial import Polynomial

import sindy_discovery as sd     # load_run, stls, format_equation, TARGETS, lt_to_float

BUMP_BETA = 6                    # (1-s^2)^6: >=4 continuous derivatives vanishing at |s|=1
FIELD_ORDER = ["rho", "p", "sigma", "m"]
CONSERVED = {"rho", "p"}
ODD_FIELDS = {"p"}               # parity under P:(x->-x, +<->-) -- only polarity is odd
                                 # (rho, sigma, m, rho_cb2 are all even)

_trapz = np.trapezoid if hasattr(np, "trapezoid") else np.trapz   # numpy 2.x rename


# %% test-function derivative polynomials -----------------------------------
def bump_polys(beta=BUMP_BETA, max_order=4):
    """Exact polynomial derivatives of G(s)=(1-s^2)^beta, orders 0..max_order."""
    base = Polynomial([1.0, 0.0, -1.0]) ** beta        # (1 - s^2)^beta
    return [base.deriv(k) if k > 0 else base for k in range(max_order + 1)]


def eval_bump(poly, s, mask):
    out = poly(s)
    out[~mask] = 0.0
    return out


# %% monomial library (degree <= 2 over the 4 fields) -----------------------
def monomial_exponents(max_degree=2):
    exps = []
    for combo in itertools.product(range(max_degree + 1), repeat=len(FIELD_ORDER)):
        if sum(combo) <= max_degree:
            exps.append(combo)
    return exps                                        # deg<=2 monomials over FIELD_ORDER


def monomial_label(exps):
    if sum(exps) == 0:
        return "1"
    parts = []
    for e, f in zip(exps, FIELD_ORDER):
        if e == 1:
            parts.append(f)
        elif e > 1:
            parts.append(f"{f}^{e}")
    return "*".join(parts)


def eval_monomial(stacks, exps):
    out = np.ones_like(stacks["rho"])
    for e, f in zip(exps, FIELD_ORDER):
        if e:
            out = out * stacks[f] ** e
    return out


def term_parity(alpha, exps):
    """Sign of d_x^alpha[monomial] under P:(x->-x, +<->-). Each d_x contributes -1;
    each odd-field factor contributes -1. Only fields in ODD_FIELDS are odd, so the
    sign depends on alpha + (total power of the odd fields)."""
    n_odd = sum(e for e, f in zip(exps, FIELD_ORDER) if f in ODD_FIELDS)
    return (-1) ** ((alpha + n_odd) % 2)


def parity_ok(target, alpha, exps):
    """A term may appear in d_t[target] only if it transforms with the SAME parity as
    the target field (d_t target has the target's parity). Half of the raw monomial x
    derivative library is symmetry-forbidden; dropping it here stops STLS from selecting
    parity-violating columns to fit fluctuations. E.g. for the even conserved rho, the
    tension current enters as d2[rho*sigma] (flux ~ rho d_x sigma, odd) -- kept -- while
    d1[rho*sigma] (flux ~ rho*sigma, even) is forbidden -- dropped."""
    target_parity = -1 if target in ODD_FIELDS else 1
    return term_parity(alpha, exps) == target_parity


def term_list(target):
    """(alpha, exps, label) for the target. Conserved -> alpha>=1 (flux form);
    non-conserved -> alpha in {0,1,2}. alpha=4 (Swift-Hohenberg) only on the pure
    field monomials rho, p. Parity-forbidden terms are excluded (see parity_ok)."""
    mons = monomial_exponents(2)
    alphas = [1, 2] if target in CONSERVED else [0, 1, 2]
    terms = []
    for a in alphas:
        for exps in mons:
            if sum(exps) == 0 and a >= 1:
                continue        # d^a[const]=0 -> a zero column (spurious huge coef)
            if not parity_ok(target, a, exps):
                continue        # symmetry-forbidden -> cannot be in a covariant law
            terms.append((a, exps, f"d{a}[{monomial_label(exps)}]"))
    # 4th-order (curvature/SH) on the pure fields only
    for f_i, f in enumerate(FIELD_ORDER):
        if f in ("rho", "p"):
            exps = tuple(1 if j == f_i else 0 for j in range(len(FIELD_ORDER)))
            if parity_ok(target, 4, exps):
                terms.append((4, exps, f"d4[{f}]"))
    return terms


# %% weak system for one (A,lt) group ---------------------------------------
def build_weak_system(runs, target, Hx, Ht, nqx=24, nqt=12, max_alpha=4):
    polys = bump_polys(BUMP_BETA, max_alpha)
    terms = term_list(target)
    alphas_needed = sorted({a for a, _, _ in terms} | {0})
    rowsT, rowsb = [], []

    for run in runs:
        F, t, Lx = run["fields"], run["t"], run["Lx"]
        nbins = F["rho"].shape[1]
        dx = Lx / nbins
        x = np.arange(nbins) * dx
        stacks = {f: F[f] for f in FIELD_ORDER}
        mon_cache = {exps: eval_monomial(stacks, exps) for _, exps, _ in terms}
        u = F[target]

        xqs = np.linspace(0.0, Lx, nqx, endpoint=False)
        tqs = np.linspace(t[0] + Ht, t[-1] - Ht, nqt)
        for tq in tqs:
            st = (t - tq) / Ht
            tmask = np.abs(st) < 1.0
            if tmask.sum() < 4:
                continue
            h0 = eval_bump(polys[0], st, tmask)
            h1 = eval_bump(polys[1], st, tmask) / Ht
            for xq in xqs:
                sx = ((x - xq + Lx / 2.0) % Lx) - Lx / 2.0
                ux = sx / Hx
                xmask = np.abs(ux) < 1.0
                if xmask.sum() < 6:
                    continue
                g = {a: eval_bump(polys[a], ux, xmask) / Hx ** a for a in alphas_needed}
                # b = -<d_t phi, u> = -int_t h1 * (int_x g0 * u dx) dt
                Iu0 = (g[0][None, :] * u).sum(axis=1) * dx
                b_val = -_trapz(h1 * Iu0, t)
                # Theta[.,term] = (-1)^alpha int_t h0 * (int_x g_alpha * P dx) dt
                row = []
                for a, exps, _ in terms:
                    Ix = (g[a][None, :] * mon_cache[exps]).sum(axis=1) * dx
                    row.append(((-1.0) ** a) * _trapz(h0 * Ix, t))
                rowsT.append(row)
                rowsb.append(b_val)
    return np.asarray(rowsT), np.asarray(rowsb), [l for _, _, l in terms]


# %% discovery over one group -----------------------------------------------
def discover_group(runs, Hx, Ht):
    out = {"A": runs[0]["A"], "lt": runs[0]["lt"], "lt_val": runs[0]["lt_val"],
           "n_reps": len(runs), "equations": {}}
    for tgt in sd.TARGETS:
        Theta, b, labels = build_weak_system(runs, tgt, Hx, Ht)
        out["n_query"] = Theta.shape[0]
        per_thr = {}
        for lvl in sd.THRESHOLD_LEVELS:
            coef, keep, r2 = sd.stls(Theta, b, lvl)
            per_thr[lvl] = {"terms": {labels[i]: float(coef[i]) for i in np.flatnonzero(keep)},
                            "r2": r2}
        # held-out replicate R^2 (fit on all but last rep, score on last)
        xv = None
        if len(runs) >= 2:
            Th_tr, b_tr, _ = build_weak_system(runs[:-1], tgt, Hx, Ht)
            Th_te, b_te, _ = build_weak_system(runs[-1:], tgt, Hx, Ht)
            coef, _, _ = sd.stls(Th_tr, b_tr, 1e-1)
            pred = Th_te @ coef
            ss = np.sum((b_te - pred) ** 2)
            tot = np.sum((b_te - b_te.mean()) ** 2) + 1e-30
            xv = float(1.0 - ss / tot)
        survived = set.intersection(*[set(per_thr[l]["terms"]) for l in sd.THRESHOLD_LEVELS])
        out["equations"][tgt] = {"per_threshold": per_thr,
                                 "structural_terms": sorted(survived),
                                 "heldout_r2": xv}
    return out


# %% main -------------------------------------------------------------------
if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("npz_dir")
    ap.add_argument("--out_dir", default=None)
    ap.add_argument("--Hx", type=float, default=1.2, help="x test-function half-width (box units)")
    ap.add_argument("--Ht", type=float, default=0.4, help="t test-function half-width (sim time)")
    ap.add_argument("--report_threshold", type=float, default=1e-1)
    ap.add_argument("--include_cb2", action="store_true",
                    help="add rho_cb2 as a 5th library field (closure-enrichment A/B test)")
    args = ap.parse_args()
    out_dir = args.out_dir or args.npz_dir
    os.makedirs(out_dir, exist_ok=True)

    if args.include_cb2 and "rho_cb2" not in FIELD_ORDER:
        FIELD_ORDER.append("rho_cb2")     # monomials/derivatives pick this up via len(FIELD_ORDER)
    print(f"library fields: {FIELD_ORDER}")

    runs = [sd.load_run(p) for p in sorted(glob.glob(os.path.join(args.npz_dir, "cg_*.npz")))]
    groups = {}
    for r in runs:
        groups.setdefault((r["A"], r["lt"]), []).append(r)

    results = []
    for (A, lt), grp in sorted(groups.items(), key=lambda kv: (kv[0][0], sd.lt_to_float(kv[0][1]))):
        res = discover_group(grp, args.Hx, args.Ht)
        results.append(res)
        print(f"\n=== A{A} lt{lt}  ({len(grp)} reps, {grp[0]['n_snap']} snaps, "
              f"{res['n_query']} weak queries) ===")
        for tgt in sd.TARGETS:
            eq = res["equations"][tgt]
            thr = eq["per_threshold"][args.report_threshold]
            print(sd.format_equation(tgt, thr["terms"], thr["r2"]))
            print(f"      structural: {eq['structural_terms']}"
                  + (f"   | held-out R2={eq['heldout_r2']:.2f}" if eq["heldout_r2"] is not None else ""))

    with open(os.path.join(out_dir, "wsindy_results.json"), "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nwrote {os.path.join(out_dir, 'wsindy_results.json')}")
