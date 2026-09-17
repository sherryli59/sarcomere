# Dynamic actin-myosin model

The physiological motor model is the command-line default (`--dynamic_am=true`).
Setting `--dynamic_am=false` selects the original instantaneous geometric AM
bonds and prescribed sliding velocity.

Each myosin half owns `am_heads_per_half` spatially resolved head slots. A slot
can bind at most one actin site in the physical overlap of that half and follows

```text
detached --k_on a(t) g--> weak --k_powerstroke--> strong --k_strong_off(F)--> detached
                         \--k_weak_off-------------------------------> detached
```

Here `g` is the eligible overlap fraction and `a(t)` is either tonic
(`period <= 0`) or a normalized difference-of-exponentials heartbeat
transient. The transient is evaluated only during
`am_activation_active_fraction * am_cycle_period`; the remainder of the cycle
is exactly `am_activation_floor`, providing a true diastolic interval.
Transition probabilities use `1-exp(-k dt)`.

Only the polarity-compatible half is productive. A head on the opposite half
may form a weak stereochemically nonproductive encounter when a thin filament
truly crosses the M-line, but its weak-to-strong rate is exactly zero. It still
occupies a real head and actin site and can therefore compete with productive
heads. The central bare zone contains no head anchors.

The conserved pool is not duplicated for every neighboring actin. In addition,
one thin filament cannot monopolize all 141 heads on a half-filament. Its
simultaneous weak-plus-strong capacity is

```text
N_capacity(i,m,half) = ceil[
    min(am_heads_per_half, am_max_heads_per_actin_half)
    * productive_overlap / head_bearing_half_length]
```

clamped between zero and the configured full-overlap capacity. The default
`am_max_heads_per_actin_half=24` corresponds to sharing 141 heads among the
approximately six thin-filament neighbors of a cardiac thick filament. Strong
heads reserve capacity before weak heads. If sliding lowers capacity, the most
strained excess weak heads leave through the capacity-limit pathway; strong
heads are never forcibly removed and retain their load-dependent KMC turnover.
Setting the cap to zero disables neighbor sharing for a control.

Actin capacity is also global, not recreated for every neighboring myosin.
Each 1-um thin filament has `am_sites_per_actin=182` material sites (about
5.5-nm axial spacing). One site has exactly one weak or strong owner across all
myosins, with strong chemistry and then lower connector strain determining the
winner if a resumed/legacy state contains a collision. Successful parallel KMC
associations atomically reserve a free site. The separate default
`am_max_heads_per_actin=24` limits total simultaneous occupancy across all
sites and all myosins touching one actin.

An attached slot stores two material quantities that remain fixed until
detachment: its crown anchor `(a_m, phi_m)` in the transported myosin material
frame and its bound actin coordinate `s_a`. The weak state has a soft axial
spring and the same flat-bottom radial head spring; it has no active stroke.

On the weak-to-strong transition, a separate internal stroke coordinate `q`
starts at zero and advances no farther than `am_power_stroke`, with the linear
force-velocity law

```text
v_stroke(F) = am_unloaded_velocity * max(0, 1 - F/am_stall_force)
```

Let `p_a=r_a+s_a u_a`, let `p_m` be the transported myosin crown anchor, and
define `h=p_a-p_m`, `z=h.u_a`, and `rho=|h-z u_a|`. One objective two-spring
energy generates all head forces and torques:

```text
U_head = 1/2 k_radial flat_bottom(rho-rho0)^2 + U_axial
U_axial(weak)   = 1/2 k_weak z^2
dU_axial(strong)/d(q-z) = min(k_motor max(q-z,0), F_stall)
```

Thus a strong motor exerts the tension-only, stall-capped axial force

```text
F_axial = min(am_motor_k * max(q - z, 0),
              am_stall_force)
```

This is a single finite stroke per attachment, not an indefinitely processive
walk. A freely moving actin can still unload the head; an isometric or
network-balanced configuration develops force and slows the remaining stroke.

along the actin polarity. Because `s_a` no longer follows the closest point,
relative actin-myosin sliding stretches or unloads the attached head rather
than silently moving its binding site.

Resisting force slows strong-state detachment:

```text
k_det(F) = am_k_strong_off *
           [f_floor + (1-f_floor)*exp(-max(F,0)/am_force_scale)]
```

The default `f_floor=1` represents a load-insensitive ATP-driven turnover
pathway: productive load cannot freeze a head across diastole. With the
default zero-load rate of 160/s, the 50-ms diagnostic heartbeat relaxes to
only a few residual strong heads during its 25-ms diastole.

Nonproductive elastic distortion provides a separate, continuous exit path:

```text
k_det,total = k_det,chemical *
              exp(min(theta * beta * E_nonproductive, exponent_cap))
```

For a strong head, productive tensile strain remains catch-stabilizing and is
not included in `E_nonproductive`; transverse distortion and negative axial
strain are. Weak-head energy includes both passive spring components. A hard
rupture is retained only as an emergency contour-extension guard measured from
the relaxed connector. Bound heads are evaluated directly from their stored
material coordinates, independent of association candidate or neighbor lists.

Head anchors are distributed through the productive region outside the central
10% bare zone. Each axial crown contains three sites separated by 120 degrees,
and successive crowns rotate by 40 degrees. The 141-site default gives 47
crowns per half, approximately 14.3 nm apart for a 1.5-um thick filament.
Forces and the physical torque
on each rod are negative derivatives of `U_head`; they are not assembled from
separate endpoint forces. Consequently every isolated head has zero net force
and zero net torque about an arbitrary origin. Summed head moments rotate the
thick filament toward the local actin bundle. The explicit path does not use
the legacy direct `kappa_am` angular potential.

There is no status-2-only bare-zone wall or `k_short` force. Excess penetration
is handled only by finite head-bearing overlap, shared-site competition,
explicit connector mechanics, and network forces (ACTN2/titin where present).
Consequently an artificially pre-compressed state is not guaranteed to reverse:
tension-only titin is slack under compression and cannot supply an outward
force. This is an explicit physical limitation rather than a hidden penalty.

Status-2 actin crosslinks are deliberately absent from the AM transition rules.
They change AM kinetics mechanically: a crosslinked actin relaxes its axial motor
strain less readily, so the force entering the rates persists. Only strong AM
states count when classifying an actin-actin bond as status 2.

The command-line default is `--aa_status_mode=arm_load_productive`. It treats
status 2 as a mechanically loaded state of an already-existing AA crosslink.
For each actin, the productive strong-head force is allocated across all ACTN2
arms in proportion to their favorable raw transmitted load.  For link `l` on
actin `i`,

```text
allocated_motor_i,l = productive_motor_i * raw_arm_i,l / sum_l(raw_arm_i,l)
effective_arm_i,l = min(raw_arm_i,l, allocated_motor_i,l)
catch_load_l = min(effective_arm_i,l, effective_arm_j,l)
```

Thus the spring force supplies a transmission capacity and load-sharing weight,
not an independent source of catch stabilization. The sum of allocated loads
cannot exceed the productive myosin force on an actin, and a slack link cannot
become status 2 merely because motors are active nearby. Both arms must be
supported and the local link must carry bilateral load. In schematic form,

```text
status2 = bilateral_ACTN2_load && productive_opposed_motor_support
```

Per-frame local-link diagnostics store the effective KMC arm loads first and
append the raw spring-arm forces and allocated motor budgets separately.

The legacy `k_aa` and `kappa_aa` values describe the former one-edge-per-pair
interaction. Explicit ACTN2 molecules therefore use fixed per-link values

```text
k_single = k_aa / aa_stiffness_reference_links
kappa_single = kappa_aa / aa_stiffness_reference_links
```

with `aa_stiffness_reference_links=3` by default. The divisor never uses the
instantaneous occupancy: a surviving molecule does not stiffen when a neighbor
detaches. Consequently pair stiffness rises and falls naturally with the
number and geometry of occupied explicit links.

`opposed_and_transmitted`, `aa_force_projection`, `motor_only`,
`transmitted_only`, and the old `legacy_motor_identity` rule remain available
as diagnostic controls. All load-based modes use hysteretic on/off thresholds,
and all newly formed AA bonds begin in status 1.

Because individual strong attachments turn over faster than an AA crosslink,
`am_status2_memory` can retain status 2 for an already-established AA bond after
its mechanical support disappears. Mechanical load is still required to
establish or re-establish status 2. The memory cannot create a new status-2
bond, and it clears on AA rupture, geometric loss, or expiry. It applies to both
the legacy and load-based classifiers. The default is zero, which gives
instantaneous classification. This is state memory only: it does not invent
residual force, so the catch-lifetime term still uses the instantaneous load.

Parallel head slots already provide an emergent memory: status-2 support is
lost only when the relevant actins simultaneously lose sufficient strong-head
force. In the fixed two-actin diagnostic, status-2 occupancy increased from
about 0.62 with one slot per AM pair to 0.992 with eight and 0.998 with 35,
without explicit memory.

## Default dynamic parameters

| option | default | meaning |
|---|---:|---|
| `am_k_on` | 100 | weak attachment rate |
| `am_heads_per_half` | 141 | 47 three-site crowns owned by each myosin half |
| `am_sites_per_actin` | 182 | global material sites (~5.5-nm spacing for a 1-um actin) |
| `am_max_heads_per_actin` | 24 | total weak+strong occupancy across all neighboring myosins |
| `am_weak_binding_floor` | 0.1 | basal weak-only attachment gate at zero activation |
| `am_k_weak_off` | 200 | weak detachment rate |
| `am_k_powerstroke` | 200 | weak-to-strong rate |
| `am_k_strong_off` | 160 | zero-load strong detachment rate |
| `am_force_scale` | 1 | coarse-grained force scale for load sensitivity |
| `am_motor_k` | 1000 | axial spring constant (pN/um; 1 pN/nm) |
| `am_power_stroke` | 0.005 | stroke distance (5 nm if one length unit is 1 micrometre) |
| `am_stall_force` | 4.5 | maximum force per strong head (pN) |
| `am_unloaded_velocity` | 5 | zero-load working-stroke speed (um/s) |
| `am_weak_stiffness_ratio` | 0.1 | weak/strong passive stiffness ratio |
| `am_tether_saturation_heads` | 6 | occupied heads needed for mature pair-level transverse stiffness |
| `am_strain_dependent_detachment` | true | enable continuous elastic-strain detachment |
| `am_strain_detachment_theta` | 0.1 | barrier fraction coupled to nonproductive energy |
| `am_strong_off_floor_fraction` | 1 | residual ATP-turnover fraction at high load |
| `am_rupture_distance` | 0.1 | emergency excess connector extension |
| `am_cycle_period` | 0 | tonic activation; set positive for heartbeat pulses |
| `am_activation_rise` | 0.02 | pulse rise time |
| `am_activation_decay` | 0.2 | pulse decay time |
| `am_activation_floor` | 0 | minimum activation |
| `am_activation_active_fraction` | 0.5 | systolic calcium-active fraction of each cycle |
| `am_status2_memory` | 0 | status-2 retention time after strong AM support is lost |

Example:

```bash
sarcomere --dynamic_am --am_cycle_period 1 \
  --am_activation_rise 0.02 --am_activation_decay 0.2 \
  --am_activation_active_fraction 0.5 \
  --am_k_strong_off 160 --am_power_stroke 0.005
```

Per-frame diagnostics are stored under `/actin/am_weak_count`,
`/actin/am_strong_count`, `/actin/am_motor_force`, and `/actin/am_activation`.
The weak and strong counts are physical head counts, not AM-pair counts.
Compact aggregate pair state plus every head state, fixed actin material
coordinate, binding azimuth, stroke progress, nonproductive strain energy and
instantaneous detachment hazard are stored under `/state`; resumable mechanical
state remains under `/resume`. Dynamic AM model version 22 rejects older
checkpoints because the crown layout, retention law, ATP-limited strong-head
off-rate, global site ownership, opposite-half weak competition, per-link
motor-load allocation, and finite-systole heartbeat gating changed. An occupied
AM pair also has a weak, occupancy-scaled nematic
compliance outside a force-free angular band (defaults:
`am_alignment_k=30 pN um/rad^2`, `am_alignment_slack_deg=2`). This represents
passive resistance of the unresolved head ensemble; explicit strong-head
forces remain the source of active alignment torque. Every dynamic motor
transition is stored under
`/actin_myo/head_transitions` with
actin/myosin IDs, step, old/new state, transition reason, activation, strain,
trajectory frame are preserved in `/resume/am_transition_events_pending`.
Schema version 3 appends `head_id`; the aggregate pair state remains for legacy
analysis compatibility.

Summarize completed weak/strong dwells and cause-specific exits with:

```bash
python analysis/summarize_am_transitions.py trajectory.h5 --dt 1e-5
```

The current implementation is a minimal coarse-grained crossbridge cycle, not a
full ATP hydrolysis scheme. HCM perturbations can independently target motor
stiffness/stroke, strong-state detachment, power-stroke rate, attachment rate, or
heartbeat activation without changing the actin-actin catch-bond law.

The staged validation workflow is in `DYNAMIC_AM_DIAGNOSTIC_PLAN.md`.
