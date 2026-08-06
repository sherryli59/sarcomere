# Sarcomere Self-Assembly Model — Code Review

**Date:** 2026-06-16
**Branch:** `am_kmc`
**Scope:** actin/myosin/crosslink mechanics, Langevin dynamics, catch-bond KMC,
entropy/phase-diagram formulation, physical validity, and a path toward a minimal
load-dependent self-assembly model.

**Files reviewed (primary):**
- `cpp/src/langevin.cpp`, `cpp/include/langevin.h`
- `cpp/src/sarcomere.cpp`, `cpp/include/sarcomere.h`
- `cpp/src/interaction.cpp`, `cpp/include/interaction.h`
- `cpp/include/components.h`, `cpp/include/geometry.h`
- `cpp/src/main.cpp`
- `analysis/entropy_transduction.py`

---

## 1. Langevin Dynamics — `langevin.cpp:sample_step`

**Formulation (valid).** The translational update is overdamped Euler–Maruyama
(Einstein–Smoluchowski):

```
Δx = β·D·F·dt  +  v·dt  +  √(2·D·dt)·η
```

Noise variance `2·D·dt = 2·μ·kBT·dt` satisfies the fluctuation–dissipation
theorem (`μ = βD`). The active velocity `v` enters additively as a prescribed
drift, which is correct for a motor-driven overdamped system.

**Rotational update is an approximate S² step.** `langevin.cpp:332`:

```cpp
vec new_dir = static_cast<vec>(model.myosin.direction[i]) + delta_u;
new_dir.normalize();
```

Adding the stochastic increment directly to the direction vector and then
normalizing is an Euler step on the sphere. It does not project the increment
onto the tangent plane before stepping, so the effective rotational diffusion
picks up O(dt) corrections proportional to `D_rot`, and detailed balance on
orientation is broken at finite dt. Fine for small `dt`/`D_rot`, but it will give
slightly wrong orientation marginals if you ever test orientation against a
Boltzmann weight. Standard fix: geodesic (Lie-group) update,
`u → normalize(rotate(u, dt·τ×u + √(2 D_rot dt)·ξ_⊥))`.

**Inconsistent capping of translation vs. rotation.** Translations are capped
with a 30% slack (`displacement_slack = 1.3`, `langevin.cpp:278`), but the
myosin rotation cap block is stubbed out (`langevin.cpp:317–321`) and the actin
one computes magnitudes without acting on them (`langevin.cpp:392–399`). A large
torque can therefore still flip a filament in a single step. Recommend either
capping rotational increments symmetrically or logging when
`delta_u.norm() >> D_rot·dt`.

---

## 2. Motor Velocity Model — `sarcomere.cpp:_calc_am_force_velocity`

**Valid structure.** Linear force–velocity `V = V₀·(1 − f_load)` is the simplest
valid motor model. Real myosin-II is convex (Hill), but linear is standard in
coarse-grained filament work.

**Two distinct load variables.** `f_load` (linear, drives velocity) and
`f_load_cb` (exponential saturation, drives catch-bond kinetics) are computed
separately (`sarcomere.cpp:729`, `:744`). The split is physically motivated
(binding geometry affects sliding and bond stability differently) but couples
two non-trivial nonlinearities. For a minimal model both should collapse to a
single order parameter.

**Discontinuous actin velocity at bond formation** (`sarcomere.cpp:958–1004`):

- `cb_status == 2`: `v = (1 − f_load)·v_am·d̂`
- bound, `cb_status != 2`: `v = v_am·d̂` (full speed, load-independent)
- unbound: `v = 0`

The middle branch slides at full speed irrespective of force until a catch bond
establishes, producing a discrete jump in velocity when `cb_status` reaches 2.
Expect numerical transients around bond formation; worth smoothing in any
reduced model.

---

## 3. Catch-Bond KMC — `sarcomere.cpp:_cb_decide`, `determine_cb_status`

**Lifetime model (valid, leading order).**

```cpp
lifetime  = base_lifetime + directional_base_lifetime + lifetime_coeff·f_load;
k_off_adj = dt / lifetime;
```

This is a *linear* catch bond, τ(f) = τ₀ + λ·f — the leading Taylor term of any
catch bond. It omits the catch→slip turnover seen in real bonds (selectin,
α-actinin). Adequate for sweeps; replace with a two-pathway Bell model for
quantitative protein-specific comparison.

**Status semantics are physically correct.** Status 2 requires both actins to be
within `aa_cutoff`, under tension, antiparallel (directional case), and engaged
to *different* myosin filaments. That is precisely the sarcomeric topology (two
oppositely-oriented actins bridged by a common thick filament region).

**Load combination rule** (`sarcomere.cpp:1342`):

```cpp
double f_load = abs_cos_angle * std::min(f_load_i, f_load_j);
```

Taking the minimum makes the bond as strong as the weaker partner — conservative
and defensible (bond fails at the weaker tether), but a lightly loaded actin can
suppress catch bonding with a heavily loaded one. A weighted average may be more
physical.

---

## 4. BUG — Myomesin spring forces computed but never applied

`sarcomere.cpp:1058–1096`, `_apply_myomesin_spring`:

```cpp
vec force_on_i = force_scalar * unit_axial;
vec force_on_j = -force_on_i;
// forces are NEVER written to local_myosin_forces[i] / [j]
```

Verified: within the function body the only reference to `local_myosin_forces`
is the parameter declaration. The M-line (myomesin) coupling between thick
filaments is inert — any run with `k_mm > 0` (`--k_mm` / `--kmm`) experienced
zero myomesin force. **Any sweep using `k_mm > 0` should be rerun after fixing.**

Fix:

```cpp
local_myosin_forces[i] += force_on_i;
local_myosin_forces[j] += force_on_j;
```

**Severity: High** (silently invalidates a parameter axis).

---

## 5. BUG — Data race on diagnostic bond counters

`sarcomere.cpp:1427–1431`, inside `_set_cb(i, j, status)`:

```cpp
actin_n_bonds[i] += 1;
actin_n_bonds[j] += 1;       // j touched by multiple threads
if (status == 2){
    actin_strong_cb_count[i] += 1;
    actin_strong_cb_count[j] += 1;
}
```

`_process_catch_bonds` runs the outer loop over `i` in parallel; different
threads handling different `i` can write the same `actin_n_bonds[j]` /
`actin_strong_cb_count[j]`. These counters are not on the force critical path,
so trajectories are likely unaffected, but the recorded counts can be corrupted.
Fix with `std::atomic` or thread-local accumulation + reduction.

**Severity: Low** (diagnostics only — confirm they are not used downstream).

---

## 6. Physical concern — always-on axial-gap energy

`interaction.cpp:120–123`, `aa_energy`:

```cpp
real ds_pos = (ds - 0.1*length1)/(0.1*length1);
U_ax += 0.5 * k_aa * (ds_pos) * (ds_pos);
real energy = 0.5 * (k_aa * dist * dist) + angle_energy + U_ax;
```

The quadratic term `0.5·k_aa·ds_pos²` is added unconditionally (the early
`raw_dist > cutoff` return only skips the radial term, not this). It grows with
any axial offset, so it acts as an effectively infinite-range positional
potential pulling any neighbor-listed actin pair toward ~10% barbed-end stagger
(the "optimal overlap" from commit `6f69152`). This is likely stronger and
longer-ranged than intended. Recommend gating `U_ax` inside the
`raw_dist < cutoff` branch.

**Severity: Medium** (affects spatial registration / spacing results).

---

## 7. Entropy Production — `analysis/entropy_transduction.py`

**Mechanical (Sekimoto motor-friction) estimator (valid).**

```
ΔS/kB = β·γ|V|²·Δt = |V|²/D·Δt,   V = (1 − f_load)·v_am
```

Correct for a prescribed-velocity overdamped system. The documented exclusion of
the conservative spring term is sound: ⟨F·dx⟩_spring → 0 at steady state, and the
F·V cross-term is load-resisting (anti-correlated) and dominates before
converging, producing spurious negatives.

**Chemical (transition-affinity) entropy (valid).**

```
affinity = log(k_on / k_off) = log(k_on · τ(f_load))
```

Formed bond → `+affinity`, broken bond → `−affinity`. Textbook stochastic
thermodynamics for a two-state system with load-broken detailed balance.
Engine-side cumulative accumulation at event time is preferred over the
saved-frame fallback (which misses sub-frame bind/unbind cycles).

**Attribution caveat.** `_record_chem_transition` (`sarcomere.cpp:1311–1323`)
splits each pair event 50/50 across both actins; the Python sum over actins then
recovers the full affinity once — correct for *system totals*. For *per-filament*
entropy maps (e.g., distinguishing Z-disc actins), the symmetric split is not
physically grounded (ownership depends on polarity) and can mislead.

**IFT consistency check.** `⟨exp(−ΔS_segment)⟩` (`:654`) is the integral
fluctuation theorem proxy; ≈1 near equilibrium (small λ), dropping with
irreversibility. Good built-in validity check to plot vs. `lifetime_coeff`.

---

## 8. Phase-Diagram Formulation

Sweep axes are (# actins A, `lifetime_coeff` λ); order parameters include
coherence `phi_coh`, density-wave amplitude `phi_amp`, and entropy-rate
components. The transition tracks closure of the catch-bond feedback loop:

```
sliding → f_load↑ → τ(f)↑ → P(bound)↑ → antiparallel registration → sliding (amplified)
```

Entropy rate peaking near (not at) the boundary is consistent with maximum
dissipation at the transition.

**Recommendation.** The natural dimensionless control parameter is
`λ·f_char / τ₀`, not λ alone. Because density A shifts the characteristic
f_load, the boundary will be non-vertical in raw (A, λ) space — fine to study,
but worth stating, and collapsing onto the dimensionless group would make the
boundary sharper and more transferable.

---

## 9. Toward a Minimal Load-Dependent Self-Assembly Model

The essential physics is the catch-bond positive feedback loop (Sec. 8). Three
reduction candidates:

**A. 1D driven lattice with internal state.** Each actin → position + polarity
bit (±1); myosin → spring linking two actins; bond state binary with
`k_on` (proximity-gated) and `k_off(f) = 1/(τ₀ + λf)`, `f` from relative
velocity. Yields a driven two-species lattice; NESS analyzable by master
equation; transition = condensation of antiparallel pairs.

**B. Mean-field self-consistency.** Let `p` = pair catch-bond probability,
`f = f(⟨p⟩)`:

```
dp/dt = k_on·(1−p) − k_off(f(p))·p = 0
```

The 1→3 fixed-point bifurcation in (A, λ) analytically predicts the phase
boundary — a cheap first map to compare against the full simulation.

**C. Agent-based lattice reduction.** Keep KMC kinetics, discretize positions to
binding sites; polarity constraint → Ising-like coupling; motor drift breaks
detailed balance. Tractable for thousands of filaments.

**Observables already available for validating any reduction:**
- completed-lifetime distribution (`aa_completed_lifetimes`)
- f_load at formation/breaking (`cb_breakage_events`)
- IFT proxy `⟨exp(−ΔS)⟩` (`entropy_transduction.py`)

---

## Summary Table

| # | Issue | Severity | Location |
|---|---|---|---|
| 4 | Myomesin spring forces computed but never applied | **High** | `sarcomere.cpp:1058–1096` |
| 6 | `U_ax` always-on (effectively infinite-range) | Medium | `interaction.cpp:120–123` |
| 1 | Rotation update not geodesic on S² | Approximation | `langevin.cpp:332–339` |
| 1 | Rotation cap stubbed out / inert | Numerical safety | `langevin.cpp:317–321, 392–399` |
| 5 | Race on `actin_n_bonds` / `strong_cb_count` | Low (diag) | `sarcomere.cpp:1427–1431` |
| 2 | Actin velocity discontinuity at cb_status==2 | Model choice | `sarcomere.cpp:958–1004` |
| 3 | `min(f_load_i, f_load_j)` bond load rule | Model choice | `sarcomere.cpp:1342` |
| 7 | Chemical entropy 50/50 per-actin attribution | Analysis caveat | `sarcomere.cpp:1311–1323` |

**Top priority:** fix the myomesin spring (#4) and rerun any `k_mm > 0` sweeps.
All other items are known approximations or low-impact diagnostics.
