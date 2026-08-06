# Coarse-grained field-theory investigation of the sarcomere KMC model

**Data:** `/scratch/users/xsli/sarcomere/data0430` — full grid A ∈ {600,900,1200,
1500,1800,2100} × lt ∈ {0, 1e-5.5 … 1e0} × 3 replicas = 234 runs, box 12×2×2,
M1000, ~4007 frames. (`lt` = `--lifetime_coeff` = directional catch-bond strength;
overlap length `am_optimal`=`aa_optimal`=0.06 fixed.)
**Outputs:** `/scratch/users/xsli/sarcomere/cg_fields0430/` (phase diagrams,
tc/cg npz), `cg_wsindy0430/` (discovery). **Tooling:** `analysis/`.

Goal: go particle-sim → continuum field theory, *grounded in trajectory data*
(not postulate-and-fit). This report consolidates what the data does and does
not support.

---

## 1. Longitudinal structure — period, threshold, coarsening

Measured by the 3D Fourier order parameter (winning nx → q_dom, period = Lx/nx)
run over the FULL trajectory (`coarse_grain_fields.py`), aggregated across the
grid (`aggregate_phase.py`, replicas averaged).

- **Band period is a density-controlled CROSSOVER, NOT steric everywhere.**
  [CORRECTED — an earlier read called it "steric, inherited not derived"; that
  holds ONLY in the dense/strong-bond corner.] Measured properly by DECOUPLING
  pattern existence from spacing (`period_vs_control.py`): gate on longitudinal
  existence (nx-stability + cb2 count), NOT phi_final (which small-N inflates at
  low lt), then fit log-log slopes.
    - HIGH density (A≥1500) OR HIGH lt (≥1e-1.5): period flat ≈1.0 (≈actin length),
      slope ~0 → STERIC packing (`fig_period_vs_lt_gated.png`).
    - LOW density + intermediate lt: period LARGE (6–12), → 1.0 as A↑, slope
      d log(λ)/d log(A) ≈ −1.8 to −2.7 (fit R²=0.66–0.80) → DENSITY-SET, not
      steric. Physically: nx→1 = a single condensate filling the box, splitting
      into more bands as density rises until it hits the steric limit.
      (`fig_period_vs_density_gated.png`, `fig_existence_map.png` confirms these
      are REAL patterns: existence=1, ≥100 cb2, stable nx — not small-N noise.)
    - lt ≲ 1e-4: no pattern (existence gate → 0).
  Catch bond is primarily an EXISTENCE GATE (sharp threshold ~1e-4); it tunes the
  scale only weakly and only at low density (d log(λ)/d log(lt) ≈ −0.4 at A600,
  ≈0 at A≥1500). CAVEATS: near the crossover it's bistable (reps split nx≈1 vs
  nx≈10 → huge error bars = real coexistence); the low-density arm is PINNED at
  box size (nx→1, λ→Lx=12) so cannot distinguish a genuine interaction-set
  wavelength from finite-size single-condensate — the wider box would settle it.
  Implication: the model DOES set its own length over much of the phase diagram
  (partial wavelength selection), steric only asymptotically.
- **Sharp catch-bond assembly threshold** (`fig_smectic_heatmap.png`): lt ≲ 1e-4
  (and lt=0) → zero order; above ~1e-3.5 layering switches on, strongest at high
  density × intermediate lt (1e-2.5–1e-1, smectic fraction 0.8–0.99); degrades at
  lt1e0 (kinetic trapping). The directional catch bond is NECESSARY, with a
  quantitative onset and an optimal intermediate strength.
- **Assembly = coarsening → arrest:** Φ(t) grows monotonically; nx*(t) starts
  high (~20–28 bands) and settles to nx≈10–11 (A1200)/15 (A1500). Fields are
  physically correct: sharp ρ bands, co-located bound-myosin peaks, oscillating
  polarity p, near-binary crosslink tension σ.

## 2. Transverse registration — smectic vs crystal

The heart of the crosslinker hypothesis (registered Z-sheets vs dephased columns)
lives in the y,z directions, which a 1D projection discards. The old
`smectic_fraction` is amplitude-blind (registered crystal and dephased smectic
both win at (ny,nz)=(0,0)). Replaced by the amplitude ratio (`transverse_coherence.py`):

    R = A_coh / A_col ∈ [0,1],   A_coh=|<e^{-i q* x}>_all|,  A_col=Σ_c n_c|ψ_c|/N
    R→1 registered crystal (Z-sheets) ; R→0 dephased smectic

Full-grid R sweep (replicas averaged, `fig_R_heatmap.png`, `fig_R_vs_lt.png`):
- **R_σ ≈ R everywhere** — tension and density register as one structure.
- **Strong-bond + low-density dephasing** (lt1e0: R→0.3–0.55 at A600/900/1200
  while A1800/2100 hold ~0.7) — kinetic trapping, confirmed across replicas.
- **KEY LIMIT: 12×2×2 cannot resolve absolute registration.** R floors ~0.3
  (never →0 even in disorder) — finite-size-inflated; "high" R ≠ true crystal.
  Only relative trends are trustworthy. A single-replica A1200 dip did NOT survive
  averaging (retracted).
- **Mechanism note:** `k_bundle` (=0.1 here) is a force toward the global transverse
  COM — mere confinement, NOT band-phase registration. The registration mechanisms
  in the model are myosin multi-binding (`max_myosin_bonds`) and the aa crosslink
  network. So the mechanism test = sweep `max_myosin_bonds`, NOT k_bundle.

## 3. Field-theory discovery (WSINDy) — structure yes, coefficients no

Strong-form STLS (`sindy_discovery.py`) FAILED: spectral-derivative noise, R²≈0,
negative held-out R², unphysical non-divergence terms in the conserved d_t ρ.

Weak-form WSINDy (`wsindy_discovery.py`) — library terms d_x^α[monomial], weak
inner product moves all derivatives onto a compact test function (never
differentiates noisy data); conserved ρ,p restricted to α≥1 (flux form):
- **Discovers consistent, physical STRUCTURE across all (A,lt):**
  - `d_t σ ≈ −k σ(1−σ)` — logistic self-saturation (σ near-binary, slaved;
    its ~static evolution is correctly discovered as ≈0 dynamics).
  - `d_t ρ, d_t p` dominated by `d1[σ], d1[σ²]` — density/polarity transported by
    crosslink-**tension gradients**, ∂_x[f(σ)]. Conservation respected.
- **But coefficients are NOT reliably identified.** In-sample R² 0.05–0.23;
  held-out-replicate R² mostly ≤0. Every lever was tested and did NOT rescue it:
  - window: fields genuinely evolve in [300,1800] (not dead time);
  - test-function width Hx∈{1.2,2.5,4}×Ht∈{.15,.3}: mean held-out stays negative;
    only a minority of slow-mode equations discoverable (max held-out 0.66);
  - closure enrichment with rho_cb2: does NOT help (mean up but max collapses
    0.66→0.16 — collinear with ρ,σ; "missing local field" refuted);
  - nonlocal ℓ0 kernel: sub-grid (dx≈ℓ0), not separable.
- **Conclusion: the bottleneck is STATISTICAL.** Single-realization ∂_t in a
  12×2×2 box with 3 replicates is fluctuation-dominated; the deterministic reduced
  PDE is identifiable in FORM but not in COEFFICIENT from this data.

## 4. Discovered reduced-model structure (the takeaway)

Even though coefficients are not pinned, the *form* that survives is a coherent
reduced field theory:

    ∂_t σ = −k σ(1−σ)                         (fast, slaved: σ saturates → binary)
    ∂_t ρ = ∂_x[ f(σ) + ... ]                 (density flux driven by tension gradients)
    ∂_t p = ∂_x[ g(σ,ρ) + ... ]               (polarity flux, tension-gradient coupled)
    ∂_t m : couples to σ and m                (motor engagement tracks tension)

with the longitudinal length fixed by steric packing (§1) and registration set by
myosin multi-binding / crosslink network (§2), not by k_bundle.

## 5. What needs new data — the wider-box simulation

Two coefficient-level questions both require escaping the thin 12×2×2 box:
- **Registration:** absolute smectic-vs-crystal needs a wider cross-section to
  resolve the transverse coherence length ξ_perp (R floors in 2×2).
- **PDE coefficients:** more spatial statistics + smaller relative fluctuations
  are what the WSINDy closure needs to move from form to reliable coefficients.

Drafted (NOT submitted): `cpp/wider_box_registration.sbatch` — 12×6×6 (9× area,
density-matched A13500/M9000), lt1e-2, `max_myosin_bonds` ∈ {6,1} (the mechanism
contrast) × replicas. Heavy: ~9× particles, days of runtime, ~60GB/file. Recommend
a 1-task pilot first. One sim advances BOTH open threads.

---

## Tooling index (`analysis/`)
- `coarse_grain_fields.py` — CG fields + q_dom(t)/smectic (has --start_frame/--end_frame)
- `transverse_coherence.py` — registration R, R_σ, ψ(y,z)
- `aggregate_phase.py` — A×lt phase diagrams (period, smectic, R), replicas averaged
- `plot_coarse_grain.py` — per-run field plots
- `sindy_discovery.py` — strong-form STLS (kept for comparison; fails)
- `wsindy_discovery.py` — weak-form WSINDy (the working method)
- Specs: `reports/sindy_scaffold_spec.md`, `reports/transverse_coherence_spec.md`
