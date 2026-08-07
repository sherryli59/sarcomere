# Accelerated-dynamics parameter fitting

Written 2026-08-06. Purpose: choose a parameter set that reaches sarcomere
assembly inside ~500k steps, without pretending the numbers are physiological,
but keeping the *proportions* that decide the outcome.

All numbers below are measured from this session's runs unless marked
"nominal" (i.e. read off the CLI).

---

## 1. The obstruction: step count is a dynamic range

For explicit overdamped integration,

```
dt      <=  0.2 * tau_fast          tau_fast = 1 / (k_max * mu)
N_steps  =  T_window / dt           T_window >= tau_slow
=>  N_steps >= 5 * tau_slow / tau_fast
```

**Multiplying every rate by a constant buys nothing** — `tau_slow` and
`tau_fast` shrink together and `N_steps` is unchanged. Acceleration requires
*compressing the ratio*, which necessarily breaks some dimensionless group.
The design question is which group is safe to break.

Measured today:

| quantity | value | source |
|---|---|---|
| `mu_actin` (parallel) | 70.2 um/(pN*s) | D*beta, D=0.3003, beta=233.644 |
| `tau_fast` | **4.75 us** | 1/(am_tether_k=3000 * mu) |
| `tau_slow` (rotational diffusion) | **4.95 s** | 1/(2*D_tumble), D_tumble=0.10096 |
| dynamic range | **~1e6** | => ~5e6 steps needed |
| budget | 5e5 steps | **10x short** |

So the model as parameterised is an order of magnitude outside the budget, and
that is *before* asking whether assembly needs more than one rotational time.

---

## 2. What is free and what is not

This is the single most useful distinction found today.

**Free to accelerate (no dt cost):** anything whose rate does not enter
`tau_fast`. Chemical rates, power-stroke size, activation rates. Driven
sliding speed is `stroke x cycle_rate`, and neither term touches dt.

**Not free:** diffusion. With `dt_max ~ 1/(k*mu)`,

```
<x^2> per step = 2*D*dt = 2*(mu*kT)*dt  ~  mu / (k*mu) = 1/k
```

The mobility **cancels**. Raising `mu` to "speed up diffusion" is exactly
self-defeating: dt shrinks in proportion. **Diffusive exploration per step is
set by stiffness alone.** The only lever is softer springs.

Consequence for design: *prefer mechanisms that are driven, not diffusive*.
A model that waits on diffusive search cannot be accelerated cheaply; a model
where motors and alignment torques do the searching can.

---

## 3. Fix these in physical units — they define "assembled"

Lengths are the one thing that must not be rescaled, because the target
structure is defined by them.

| parameter | value (um) | why fixed |
|---|---|---|
| `actin_length` | 1.0 | sets sarcomere geometry |
| `myosin_length` | 1.5 | " |
| `actin_radius` / `myosin_radius` | 0.015 / 0.025 | lattice spacing |
| `am_optimal`, `am_surface_gap` | 0.055, 0.015 | AM centreline |
| `aa_preferred_centerline` | 0.055 | Z-disc spacing |
| `aa_status2_registry_on/off` | 0.040 | **the register window** |
| `am_sites_per_actin` | 182 (5.49 nm) | actin helical repeat |
| box / density | 41.45 actins/um^3 | crowding |

Everything else is fair game.

---

## 4. Dimensionless groups to preserve

These are the proportions that decide whether assembly happens. Keep them
physical even when the absolute units are not.

| # | group | meaning | current value | target |
|---|---|---|---|---|
| G1 | `dt / tau_fast` | integration accuracy | 0.11 | **<= 0.2** |
| G2 | `duty = k_ps/(k_ps+k_weak_off)` etc. | fraction of heads bound | ~0.39 | keep |
| G3 | `F_motor / F_rupture_AA` | can motors tear crosslinks | ~13.5 / 21.6 = 0.6 | **< 1** |
| G4 | `v * tau_AA / L_register` | sliding per crosslink lifetime | 0.52*0.0087/0.04 = **0.11** | ~1 |
| G5 | `(L_rupture - L_pref) / L_register` | rupture headroom vs register | (0.12-0.055)/0.04 = **1.6** | **>= 10** |
| G6 | `sqrt(kT/k) / L_register` | thermal jitter vs register | 0.0012/0.04 = 0.03 | << 1 |
| G7 | `Pe = v * L_actin / D` | drive vs diffusion | 0.52*1/0.3 = 1.7 | >= 1 |
| G8 | `T_window / tau_AA` | crosslink turnovers per run | **57** | **3 - 10** |

**G5 and G8 are the broken ones.** G5 = 1.6 means a crosslink snaps after
1.6 register windows of relative motion, while registration requires the
filaments to search over many windows. G8 = 57 is the same defect seen from
the other end.

Measured cause (v3 event stream, 790 deaths):

```
geometric rupture (length > aa_link_rupture_length):  717   90.8%
stochastic (kinetic draw):                             73    9.2%
implied kinetic-only lifetime:                      100.2 ms  (= base_lifetime, correct)
actual mean lifetime:                                 8.7 ms
```

The kinetics are fine. The geometry is not.

---

## 5. Full parameter inventory

### A. Integration
| parameter | current | class | note |
|---|---|---|---|
| `dt` | 5e-7 | **fit** | set from G1 once `k_max` is chosen |

Measured dt ladder (20 ms physical, A160/M80): torque-cap events are
**exactly 0** at dt <= 1e-6 and switch on at 2e-6; force caps jump 750x
between 2e-6 and 4e-6. **Nothing ever NaN'd** — the model clamps silently
rather than blowing up, so cap counters are the only honest stability signal.

### B. Mobility (do not use as an acceleration knob — see section 2)
`physical_mobility`, `temperature_kelvin`, `solvent_viscosity`,
`actin_hydrodynamic_diameter`, `myosin_hydrodynamic_diameter`.

### C. Stiffness — sets dt, so this is where acceleration is bought
| parameter | current | class | note |
|---|---|---|---|
| `am_tether_k` | **3000** | **fit** | *currently sets dt.* Effective value is `3000 * min(1, heads/6)` ~ 1000 at 2 bound heads. Measured mean force 26 pN vs motor 0.34 pN |
| `actin_steric_k` | 2000 | fit | second stiffest |
| `am_motor_k` | 1000 | constrain | with `am_power_stroke` sets force per stroke |
| `myosin_steric_k`, `am_steric_k` | 1000 | fit | |
| `k_aa` | 1000 | constrain | /`aa_stiffness_reference_links`=3 -> 333 per link; sets G3 |
| `kappa_aa`, `aa_perpendicular_k`, `aa_status1_radial_k` | 100, 50, 75 | fit | |
| `am_alignment_k` | 120 | **fit** | the only thing that aligns filaments; scaled by `min(1, heads/am_tether_saturation_heads)` |
| `am_head_transverse_k` | 100 | fit | |
| `titin_k` | 100 | fit | |
| `k_mm` | **0** | **fit** | M-band spring; implemented and momentum-conserving but never switched on |

### D. Motor kinetics — free to accelerate
| parameter | current | class |
|---|---|---|
| `am_k_on`, `am_k_weak_off`, `am_k_powerstroke`, `am_k_strong_off` | 100/200/200/160 | **fit (scale together to hold G2)** |
| `am_power_stroke` | 0.005 | **fit** — best speed knob, v ∝ stroke at constant duty |
| `am_unloaded_velocity` | 5 | derived — must exceed `stroke / t_attached` or the stroke cannot finish |
| `am_stall_force` | 4.5 | constrain (G3) |
| `am_strain_detachment_theta` | 0.1 | **fit** — sets the force-velocity ceiling |
| `am_heads_per_half` | 9 | constrain — realized heads/actin is emergent from stoichiometry |
| `am_max_heads_per_actin`, `am_tether_saturation_heads` | 24, 6 | fit |

Measured speed vs rate scaling factor f (all four rates + `unloaded_velocity`):

| f | 1 | 2 | 4 | 8 | 16 |
|---|---|---|---|---|---|
| v (um/s) | 0.52 | 1.27 | 2.09 | 2.33 | 1.79 |
| strong heads | 1.2 | 0.4 | 0.2 | 0.0 | 0.0 |

**Saturates near f=4 and reverses by f=16.** Cause: `am_strain_dependent_detachment`
couples detachment to mechanical strain, which does not rescale with the
chemical rates, so at high f heads are thrown off before the filament responds.
Speed past f~4 is bought by spending force. Use `am_power_stroke` instead.

### E. Crosslink (ACTN2) — where the defect is
| parameter | current | class | note |
|---|---|---|---|
| `aa_link_rupture_length` | **0.12** | **fit — primary** | G5 = 1.6, needs >= 10 -> **0.4-0.5** |
| `base_lifetime` | 0.1 | constrain | already correct (100.2 ms verified) |
| `directional_base_lifetime`, `lifetime_coeff` | 0.05, 0.4 | constrain | status-2 only |
| `k_on` / `k_off` | 2500 / 1 | fit | attach is fast, thermal detach is irrelevant |
| `aa_status1_site_diffusion` | 0.010 | **fit** | only 13 nm of site migration per lifetime — 5x too slow to relieve 65 nm of strain |
| `aa_capture_margin`, `aa_overlap_capture`, `aa_barbed_zone` | 0.045, 0.10, 0.20 | geometry | keep |
| `max_aa_bonds_per_actin` | 2 | constrain | |

### F. Activation / titin
`am_q_calcium_open_rate`, `q_close_rate` (20, 20 -> tau 50 ms): scale with the
motor rates. `titin_k_off` = 0.01 (tau = 100 s) is effectively irreversible at
any feasible window — titin never releases once captured. Decide whether that
is intended.

---

## 6. Proposed fitting procedure

Ordered, each step measurable.

**Step 1 — fix geometry** (section 3). Never rescale.

**Step 2 — repair the broken groups before accelerating anything.**
Set `aa_link_rupture_length` from G5:
```
aa_link_rupture_length = aa_preferred_centerline + 10 * L_register
                       = 0.055 + 10*0.040 = 0.455 um
```
Then check G3 stays < 1: rupture force = `(k_aa/3) * (L_rup - L_pref)`
= 333 * 0.40 = 133 pN, vs peak motor 13.5 pN -> G3 = 0.10. Safe.
Raise `aa_status1_site_diffusion` until site migration per lifetime exceeds the
strain it must relieve.

**Step 3 — choose the acceleration factor from the budget.**
```
Lambda = tau_slow_target / (N_budget * dt)
```
with `tau_slow` = whatever assembly actually waits on. Today that is *not*
diffusion (G7 = 1.7, drive already beats diffusion) but rotational alignment,
which is motor-driven, so it is in the "free" class.

**Step 4 — accelerate the free things first.** Scale the four AM rates and
`am_unloaded_velocity` by f, and raise `am_power_stroke`. Hold G2 by scaling
all four together. Stop at f~4 (measured saturation); take further speed from
stroke size.

**Step 5 — buy dt with stiffness, cheapest spring first.** Sort springs by
`k`, and reduce only those not pinned by a dimensionless group. `am_tether_k`
is the current binding constraint at 3000; its *effective* value is already
~1000. Reducing to ~300 gives `tau_fast` = 47 us and **dt = 1e-5**, and the
geometry still holds: peak motor force 13.5 pN stretches a 300 pN/um tether by
45 nm, inside `am_cutoff` = 100 nm. Below ~300 heads stretch out of binding range.

**Step 6 — set dt = 0.2 * tau_fast, then verify with cap counters, not eyeballs.**
Force / torque / rotation cap events must stay at the dt=5e-7 baseline. A run
that clamps is not integrating.

**Step 7 — re-derive the derived quantities.** `am_unloaded_velocity` >=
`stroke / t_attached`; `save_every` for the frame budget; `interaction_ramp_steps`
should be a startup transient (~10% of the run), NOT the whole run as in
job 37757726 where the ramp never reached 1.0 until the final step.

---

## 7. Proposed starting point

| parameter | from | to | rationale |
|---|---|---|---|
| `aa_link_rupture_length` | 0.12 | **0.455** | G5: 1.6 -> 10 |
| `aa_status1_site_diffusion` | 0.010 | **0.10** | site can track 40 nm/lifetime |
| `am_tether_k` | 3000 | **300** | tau_fast 4.75 -> 47.5 us |
| `actin_steric_k` | 2000 | 300 | keep below the tether |
| `am_steric_k`, `myosin_steric_k` | 1000 | 300 | " |
| `k_aa` | 1000 | 1000 | pinned by G3 |
| `dt` | 5e-7 | **1e-5** | 0.2 * tau_fast |
| AM rates (4x) | 100/200/200/160 | 400/800/800/640 | f=4, G2 held |
| `am_unloaded_velocity` | 5 | 20 | stroke must complete |
| `am_power_stroke` | 0.005 | 0.020 | 4x sliding at constant duty |
| `k_mm` | 0 | **scan** | never been switched on |

Resulting window: 500k * 1e-5 = **5 s**, versus 50 ms today — a 100x increase
in physical time at the same cost.

Predicted at 5 s: axial diffusion RMS 1.73 um (1.7 actin lengths), driven
sliding at f=4 with 4x stroke ~ 40 um (unbounded — likely needs a throttle),
crosslink turnovers ~ 5-10 (G8 target), rotational diffusion still only 1 rad
so alignment stays motor-driven.

---

## 8. Validation gates

Do not accept a parameter set that fails any of these.

1. **Cap counters at baseline.** force/torque/rotation cap events per
   particle-step no worse than the dt=5e-7 reference. Silent clamping is the
   failure mode, not NaN.
2. **ACTN2 death channel.** geometric rupture < 50% of deaths (currently 91%).
3. **Duty ratio preserved.** strong-head fraction within 2x of the f=1 value —
   catches the "speed bought with force" trap.
4. **G3 < 1.** motor force must not exceed crosslink rupture force.
5. **Seeds.** >= 12 seeds, report P(assembled), never a 3-seed mean — this
   system is bistable per seed.
6. **Same-node A/B.** Sherlock nodes vary ~25%; never compare across jobs.

---

## 9. Open questions this does not settle

- Is `am_tether_k` behaving as intended? It is engaged 99.9% of the time,
  stretched ~20x past its 2 nm slack band, carrying a mean force 75x the
  motor's. A "flat-bottom containment" spring should be idle most of the time.
- Should `titin_k_off` = 0.01 (irreversible at any window) be reversible?
- `k_mm` has never run in any production job. The M-band is the second
  registration plane and is currently absent.
- Does assembly need long-range diffusion at all? If not, section 2 says the
  model is much cheaper to accelerate than assumed.
