# Weak-form SINDy scaffold — spec

**Goal.** From the coarse-grained trajectory fields, *discover* the 1D PDEs that
govern $\partial_t\rho,\ \partial_t p,\ \partial_t m$ (and test whether $m,\sigma$
are slaved or dynamical), and **cross-check every discovered term against the
analytic coarse-graining (Route 1)**. The deliverable is not "a PDE" — it is a
ranked term list per field with coefficients, plus the falsification tests that
say which terms are real (microscopic origin + systematic dependence on `lt`)
vs fit noise.

This is Route 2 of the plan; see `reports/` field-theory notes and the project
memory. Data semantics: axis = x, polarity $p=\rho_+-\rho_-$, `lt`=catch-bond
strength (`lifetime_coeff`), $\ell_0$=`am_optimal`=0.06 (fixed), catch-bond law
$k_\mathrm{off}\propto|\cos\theta|/(b+\texttt{lt}\,f_\mathrm{load})$.

---

## 0. Prerequisite — data contract (BLOCKS the scaffold)

The current `coarse_grain_fields.py` saves only **6 field snapshots** per run
(`snap_*`) — fine for eyeballing, **insufficient for SINDy**, which needs the
field *densely in time* to form $\partial_t U$.

**Required change:** add `--dump_fields` to `coarse_grain_fields.py` that writes,
per run, the full strided field stack on a fixed grid:

```
U[field] : shape (n_t, nbins)   for field in {rho, p, m, sigma, n_bond, rho_cb2}
x        : (nbins,)             x grid, dx = Lx/nbins
t        : (n_t,)               real time = frame_idx * save_every * dt
```

- `dx` ≈ `smooth` ≈ a fraction of $\ell_0=0.06$ → `nbins=240` (dx=0.05) is fine
  but **note dx > $\ell_0$**; if the $\ell_0$ structure matters, go `nbins≥480`.
- `t`: needs `save_every` and `dt` (defaults 200, 1e-5 → frame spacing 2e-3;
  **confirm against the sweep sbatch** `SAVE_EVERY`/`nsteps`, pass via CLI).
- **Include the assembly transient.** At steady state $\partial_t U\approx0$ and
  nothing is identifiable; the dynamic information is in the ramp-up. Sample
  from frame 0, not just the tail.

---

## 1. Weak (integral) formulation — why, and the mechanics

Differentiating particle-binned fields amplifies shot noise; **weak-form SINDy**
avoids it. For each field $U$ and each compactly-supported test function
$\phi_q(x,t)$, multiply the candidate PDE $\partial_t U=\sum_j\xi_j\,\Theta_j[U]$
by $\phi_q$ and integrate by parts over its support, moving derivatives onto
$\phi$ (known analytically):

$$-\!\int\!\!\int (\partial_t\phi_q)\,U \;=\; \sum_j \xi_j \!\int\!\!\int \big(\mathcal{D}_j^{\!*}\phi_q\big)\,\Pi_j[U]$$

where $\mathcal{D}_j^*$ is the adjoint (sign-flipped) spatial derivative operator
and $\Pi_j$ the *non-derivative* nonlinearity (e.g. $\rho^2$, $m\,p$). No numerical
derivative of data is ever taken.

- **Test functions:** separable $\phi_q(x,t)=g\!\big(\tfrac{x-x_q}{H_x}\big)\,
  h\!\big(\tfrac{t-t_q}{H_t}\big)$, with $g,h=(1-s^2)^\beta$ bumps ($\beta\!\ge\!2$),
  $g$ **periodic** in $x$ (box is PBC). Query grid $(x_q,t_q)$; widths
  $H_x\sim$ a few $\ell_0$, $H_t\sim$ a few frames. Each query point = one row of
  the linear system; each candidate term = one column.

---

## 2. Candidate library — constrained by the microscopics (the crux)

Symmetry under $x\to-x$: $\rho,m,\sigma$ even, $p$ odd, currents odd. Build each
field's library to respect its conservation law.

**$\partial_t\rho$ — conserved (Model B): every term is $\partial_x J_\rho$.**
$$J_\rho \in \{\,\partial_x\rho,\ \rho\,\partial_x\rho,\ \partial_x^3\rho,\
\rho\,p,\ p\,\partial_x p,\ m\,p,\ \sigma\,\partial_x p,\ \rho\,\partial_x m\,\}$$
- $\rho p$ / $mp$ = polarity-advective current $J\sim\rho v$, $v\propto p$ (motors).
- $\partial_x^3\rho$ = surface tension $\kappa$ (sign decides cohesion).

**$\partial_t p$ — polarity (not conserved; sorting creates it).** Model-A +
current. Novelty lives here.
$$\Theta\in\{\,\partial_x^2 p,\ p,\ p^3,\ \rho p,\ \partial_x\rho,\
\underbrace{m\,\partial_x\rho}_{\text{alignment}},\ \underbrace{\sigma\,\partial_x p,\ m\,p}_{\text{catch-bond coupling}},\ \partial_x\sigma\,\}$$
- The **directional catch-bond** enters as the $\sigma$–$p$ / $m$–$p$ couplings.
  Their coefficients **must scale with `lt`** and **vanish as `lt`→0** — built-in
  falsification using the existing sweep.

**$\partial_t m$ — reaction-dominated (binding/unbinding).**
$$\Theta\in\{\,\rho \,(\text{source }s_0),\ -m\,(k_\text{off}),\ \sigma m,\
\rho-m,\ \partial_x^2 m\,\}$$
- Catch-bond off-rate $k_\mathrm{off}\propto 1/(b+\texttt{lt}\,\sigma)$ is
  **rational**; include the linearization $-k_0 m + k_1\,\sigma m$ in the
  polynomial library, and offer a constrained rational fit
  $-m/(b+\texttt{lt}\,\sigma)$ as a second pass.
- **Adiabatic-elimination test:** if the $\partial_t m$ fit returns
  $m\approx s_0\rho/k_\mathrm{off}$ with tiny $\partial_t m$ residual, $m$ is
  slaved → elimination justified (what the field theory has been assuming).
  If not, $m$ stays dynamical. Same diagnostic for $\sigma$.

**Nonlocality / $\ell_0$.** A polynomial-gradient library captures a steric
period only by piling up high-order $\partial_x^n$ terms. **If the regression
demands $\partial_x^4,\partial_x^6,\dots$ to fit $\rho$, that is the signal we are
in the DFT/packing regime** (the document's caveat) — then augment with a
nonlocal convolution column $K_{\ell_0}\!*\rho$ ($K$ a kernel of range
$\ell_0$=0.06) rather than forcing a gradient expansion.

---

## 3. Regression

- **STLSQ** (sequential thresholded least squares) with ridge; sweep sparsity
  threshold $\lambda$ → **Pareto front** (#terms vs residual). Pick the knee.
- **Column normalization** before thresholding (terms have wildly different
  scales); report **dimensional** coefficients after un-normalizing.
- One independent regression per field ($\rho,p,m$).
- Self-contained NumPy implementation (weak-form library + STLSQ is ~few hundred
  lines) to avoid env churn; `pysindy.WeakPDELibrary` is an acceptable drop-in
  **only if already in the `sarcomere` env** (do not add a dependency blindly).

---

## 4. Acceptance / falsification (tie to Route 1)

A term is trusted only if it passes **all**:
1. **Hold-out predictive residual** — fit on a subset of query points, test on held-out $(x_q,t_q)$.
2. **Replica stability** — coefficient consistent across `r00–r02`.
3. **Systematic `lt` dependence** — the catch-bond coupling coefficients scale with `lt` and → 0 as `lt`→0; a noise term shows no such trend. *(This is the headline physics plot.)*
4. **Microscopic origin gate** — the term also falls out of the Dean–Kawasaki coarse-graining (Route 1), and free coefficients match known inputs: diffusion ≈ `actin_diff_coeff_trans=0.2`, source ≈ `k_on`, etc.

A term present in the regression but with **no microscopic origin and no `lt`
trend is fitting noise**, not a discovery.

---

## 5. Outputs

- `sindy_<run>.json` — per field: ranked terms, dimensional coefficients, Pareto, residuals.
- `fig_coeff_vs_lt.png` — discovered coupling coefficients vs `lt` (the central result).
- `table_discovered_vs_micro.csv` — discovered coefficient vs Route-1 / known input.
- `fig_pareto_<field>.png` — sparsity/accuracy trade-off per field.

---

## 6. Risks / scope

- **1D projection** integrates out $y,z$; valid only where dynamics are quasi-1D.
  **Gate on `smectic_fraction` high** (from the grounding pass) before trusting it.
- **Identifiability needs the transient** — steady-state-only data is degenerate.
- **Low per-bin counts** make $p$ noisy; tie smoothing width to $\ell_0$ and
  prefer the weak-form integral scale $H_x$ over pointwise values.
- **Discovery ≠ truth.** SINDy always returns *some* sparse PDE; §4 gates 1–4 are
  the discipline. Cross-validation against Route 1 is non-negotiable.

---

## 7. Build order

1. `--dump_fields` in `coarse_grain_fields.py` (the §0 data contract). *Blocks everything.*
2. `analysis/sindy_fields.py` — weak-form library builder + STLSQ + CV (NumPy).
3. Run across the `lt` sweep × replicas; produce §5 outputs.
4. Route-1 comparison table (needs the Dean–Kawasaki notes).
