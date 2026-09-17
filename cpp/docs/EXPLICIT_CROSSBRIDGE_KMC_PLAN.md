# Explicit Actin-Myosin Cross-Bridge KMC Plan

## Objective

Replace the current conflation of actin-myosin proximity, active motor force,
and ACTN2 load classification with three distinct pieces:

1. A passive AM tether keeps an attached motor close to its actin binding
   region.
2. An explicit detached/weak/strong cross-bridge cycle generates active force
   with a polarity-fixed power stroke.
3. ACTN2 status 2 reports load carried by an existing ACTN2 crosslink. It is
   not created directly by motor occupancy and does not itself pin an actin.

The central mechanical requirement is that stationarity emerge from force
balance. A free actin must move and unload an attached motor. A constrained or
oppositely loaded actin can remain nearly stationary while the motor carries
near-isometric force.

## Modeling boundaries

This is an explicit, nonprocessive cross-bridge model. Each geometrically
eligible AM filament pair carries a sparse ensemble of independent head slots.
Each attached slot represents one coarse-grained head, performs one power
stroke per weak-to-strong transition, and then detaches. Do not advance the
same strong head continuously along actin as if it were a processive motor.
`am_heads_per_pair` is the potential capacity at full productive overlap; it is
not the instantaneous number of attached heads.

The muscle force-velocity relation should emerge from repeated cross-bridge
cycles, elastic strain, and load-dependent transition rates. A prescribed
`v_am` remains available only as a legacy control.

## 1. Separate AM mechanics

### 1.1 Passive proximity tether

For every head attachment, store:

- actin index and myosin index;
- the actin material coordinate or binding-region coordinate;
- the myosin half and head-anchor coordinate;
- the productive polarity sign;
- the attachment step and current KMC state.

The passive tether should maintain transverse proximity without functioning as
the motor or the ACTN2 load sensor. A minimal energy is

```text
U_tether = 0.5 * k_perp * (d_perp - d0)^2
```

where `d_perp` is the shortest transverse distance between the stored AM
contact regions. The tether should have no intentional axial sliding
stiffness. The legacy direct `kappa_am` angular potential remains available
through the legacy AM path for regression comparisons, but is not combined
with the geometric torque in the explicit-KMC path.

Required behavior:

- attached AM pairs remain within the proximity cutoff;
- longitudinal actin-myosin sliding is not resisted by the proximity tether;
- loss of overlap or excessive transverse separation causes geometric
  detachment;
- tether forces obey action-reaction;
- tether force magnitude and direction are never used as motor load or ACTN2
  catch-bond load.

Both passive and active cross-bridge forces act at the head attachment point,
not at the thick-filament center. For each attached head,

```text
tau_myosin = (r_head - r_myosin_center) x F_on_myosin
```

and analogously on actin. The current unit-vector orientation integrator uses
the tangent drive `tau x filament_axis`. Directional attachments restrict
`r_head` to the productive outer third of the appropriate half of the bipolar
thick filament. Summed moments from multiple attached actins therefore rotate
myosin toward the local bundle without imposing an always-on angular bias.

The existing `compute_am_force_and_energy` routine should be split rather than
silently reinterpreted. The legacy full segment spring must remain selectable
for regression comparisons.

### 1.2 Active cross-bridge elastic element

Maintain a distinct active axial coordinate for each attached head. Let
`x` be the current axial coordinate and `x0` the active rest coordinate.

On detached-to-weak attachment:

```text
x0 = x
active_force = 0
```

On weak-to-strong power stroke:

```text
x0 <- x0 + polarity_sign * d_ps
```

The productive active strain is defined with one tested sign convention:

```text
epsilon = polarity_sign * (x0 - x)
f_active = k_motor * max(epsilon, 0)
```

The active force is

```text
F_active_on_actin = polarity_sign * f_active * actin_axis
F_active_on_myosin = -F_active_on_actin
```

This is a tension-only active element. Overshoot or assisting strain sets the
active contribution to zero rather than reversing the motor. A separate
assisting-load detachment rule can be added later if data require it.

The passive AM tether remains active in weak and strong states. Only the
strong state generates the full active force. A weak state may optionally use
a small active stiffness, but the first implementation should use zero so the
roles remain identifiable.

### 1.3 Force cap and force-velocity behavior

Add a motor stall-force parameter `am_stall_force`. Initially cap the active
force:

```text
f_active = min(k_motor * max(epsilon, 0), am_stall_force)
```

For explicit nonprocessive heads, do not impose a force-velocity curve on each
strong state as continuous walking. Instead, obtain the ensemble
force-velocity relation from:

- the power-stroke displacement;
- active spring stiffness and stall cap;
- load-dependent weak-to-strong transition;
- load-dependent strong-state detachment;
- repeated attachment cycles.

If this ensemble does not recover a reasonable force-velocity curve, add a
coarse-grained ensemble mode later. Do not mix that correction into the first
explicit-head implementation.

## 2. KMC states and rates

Retain the current minimal state graph:

```text
detached --k_on(a, geometry)--> weak
weak --k_weak_off--> detached
weak --k_ps(f_h)--> strong
strong --k_strong_off(f_h)--> detached
```

Use the nonnegative hindering active load

```text
f_h = f_active
```

for load-dependent rates. Do not use the passive tether force.

The existing exponential form may remain for the first implementation:

```text
k_ps(f) = k_ps0 * exp(-f / F_ps)
k_strong_off(f) = k_off0 * exp(-f / F_off)
```

Use exact discrete-time transition probabilities:

```text
p = 1 - exp(-k * dt)
```

and competing hazards for weak detachment versus power stroke. All rate scales
must be explicit parameters rather than sharing an undocumented force scale.
Head transitions carry an explicit `head_id`, and weak/strong diagnostics count
heads rather than filament pairs. Checkpoint and resume must preserve every
head's state and axial target.

## 3. Productive AM orientation

Determine the motor direction at attachment and store it. It must not change
because a passive spring changes direction.

The productive-orientation gate should include:

- actin polarity;
- myosin half identity;
- actin-myosin overlap;
- an optional minimum axis alignment.

Reuse the existing myosin-half and productive-polarity helpers, but add unit
tests for every combination of myosin half and actin polarity. The test should
state explicitly whether `actin.direction` is barbed-to-pointed or
pointed-to-barbed. No force implementation should proceed until that semantic
is fixed in a comment and test.

## 4. ACTN2 mechanics and status

### 4.1 Existing crosslink first

ACTN2 binding remains independent of strong AM occupancy. A newly accepted AA
crosslink begins in status 1. Status 2 is a mechanical state of that existing
crosslink.

### 4.2 Separate motor drive from transmitted load

Motor occupancy alone is insufficient: two motors can pull two crosslinked
actins in the same direction and translate the pair without substantially
loading the crosslink. Passive AA prestress alone is also insufficient.

For each existing AA bond, compute two nonnegative scalars.

First, the opposed active drive across the linked actins:

```text
L_drive = 0.5 * norm(F_AM_active_i - F_AM_active_j)
```

where each force is the sum over strong, productively oriented AM bonds on the
actin. A refined version can project the force difference onto the mechanically
relevant AA shear or separation mode, but the initial version should avoid a
signed projection along either actin axis.

Second, the magnitude of load actually transmitted by the AA connector:

```text
L_transmitted = norm(F_AA_on_i)
```

This is magnitude-only. The old requirement that both signed AA arm
projections be positive is removed.

The candidate catch-bond load is

```text
L_catch = min(L_drive, L_transmitted)
```

This construction requires both active motor drive and mechanical transmission
through the crosslink. It rejects passive AA prestress, same-direction
translation, and strong-motor occupancy with no force balance.

Before adopting it as the production rule, compare it against two alternatives
in the minimal tests:

1. `motor_only = min(|F_AM_i|, |F_AM_j|)`;
2. `transmitted_only = |F_AA|`;
3. `opposed_and_transmitted = min(L_drive, L_transmitted)`.

The third is the default hypothesis, not a foregone conclusion.

### 4.3 Status hysteresis and lifetime

For an existing bond:

```text
status 1 -> status 2 when L_catch > F_on
status 2 -> status 1 when L_catch < F_off
```

with `F_off < F_on`. A short confirmation interval may be added only after the
raw force signal is characterized. Do not use temporal filtering to hide an
incorrect load definition.

Only status 2 receives the catch-bond lifetime increment:

```text
tau_1 = base_lifetime
tau_2 = base_lifetime + directional_bonus + lambda * L_catch
```

Status changes do not create or remove the AA bond. Bond loss occurs through
the AA KMC off event, geometric rupture, or an explicitly recorded capacity
event.

## 5. Integration and force balance

The overdamped actin update must use the total physical force:

```text
gamma_a * v_a = F_AM_active + F_AM_tether + F_AA + F_steric + F_other + noise
```

There is no status-dependent velocity assignment. In particular, status 2
must not set `v_a = 0` or multiply a prescribed velocity by a catch-bond load
factor.

Near-stationarity is diagnosed using relative motion, not absolute coordinates:

- actin velocity relative to its attached myosin;
- velocity difference between the two AA-linked actins;
- AA-link extension rate;
- center-of-mass velocity of the minimal system.

A force-balanced pair may be internally stationary while its center of mass
drifts. If absolute stationarity is required, represent an external anchor or
remove global center-of-mass drift; do not pin every status-2 actin.

## 6. Code changes

### 6.1 New runtime modes

Add explicit compatibility switches:

```text
--am_mechanics=legacy_velocity|explicit_crossbridge
--am_tether=legacy_full_spring|transverse_proximity
--aa_load_mode=legacy_identity|signed_aa_projection|motor_only|
               transmitted_only|opposed_and_transmitted
--am_unidirectional_active_force=true
```

Existing trajectories must remain reproducible when legacy modes are selected.

### 6.2 Data structures

Per active AM pair, store:

- KMC state;
- attachment actin coordinate;
- myosin half/anchor coordinate;
- productive polarity sign;
- active rest coordinate `x0`;
- active strain;
- active force magnitude;
- passive tether force magnitude;
- last transition step.

Avoid adding another dense `n_actin x n_myosin` double matrix if the active-pair
representation can hold the new continuous state. If dense storage is retained
for checkpoint compatibility, document its memory cost.

### 6.3 Force routines

Split the current AM force path into:

```text
compute_am_tether_force(...)
compute_am_active_crossbridge_force(...)
```

Add helpers for:

```text
productive_motor_sign(...)
am_hindering_load(...)
aa_opposed_motor_drive(...)
aa_transmitted_load_magnitude(...)
aa_catch_load(...)
```

Keep the helpers unit-testable and free of global state where possible.

### 6.4 HDF5 diagnostics

At saved frames, write per-actin:

- summed productive AM force vector;
- summed active force magnitude;
- passive AM tether force magnitude;
- weak and strong AM counts;
- actin velocity and net-force components.

For AM transitions, record:

- old/new state and reason;
- active strain and force before the event;
- passive tether force;
- productive sign and myosin half;
- local relative velocity.

For AA transitions, record:

- `L_drive`;
- `L_transmitted`;
- selected `L_catch`;
- the three load-mode candidates;
- active force vectors on both actins;
- explicit exit reason: KMC, geometry, capacity, or status reclassification.

This exit-reason field is required to validate AA KMC hazards independently of
geometric rupture.

### 6.5 Checkpointing

Checkpoint every new per-AM-pair state needed for exact continuation, including
active rest coordinate and productive sign. Extend the checkpoint round-trip
test so resumed and uninterrupted deterministic trajectories match event by
event.

## 7. Unit tests

### 7.1 Geometry and tether

1. Zero tether force at the preferred transverse distance.
2. Equal and opposite tether forces.
3. Axial translation at fixed transverse distance generates negligible tether
   force.
4. Transverse displacement generates the expected restoring force.
5. Myosin-half identity remains stable during permitted axial sliding.

### 7.2 Active crossbridge

1. Weak attachment starts at zero active strain.
2. Power stroke creates positive force in the productive direction for every
   polarity/half combination.
3. Active force never becomes negative.
4. Overshoot produces zero active force, not reversed force.
5. Active action-reaction residual is below numerical tolerance.
6. Stall-force cap is exact.

### 7.3 KMC

1. Exact `1-exp(-k*dt)` probabilities.
2. Weak competing-hazard event rate and branch fraction.
3. Strong dwell-time hazard in several fixed-force bins.
4. Load-slowed power-stroke hazard.
5. Deterministic reproducibility across thread counts where supported.

### 7.4 ACTN2 load

1. Motors off plus passive AA prestress gives `L_catch = 0` in the combined
   rule.
2. Equal same-direction motor forces give near-zero opposed drive.
3. Equal opposite motor forces give the expected opposed drive.
4. Opposed motors with no AA bond cannot establish status 2.
5. Existing AA bond with opposed drive but zero transmitted force remains
   status 1.
6. Existing AA bond with opposed and transmitted load becomes status 2.
7. Reversing the passive AA spring vector without changing its magnitude does
   not change the magnitude-based load.

## 8. Minimal test simulations

All minimal tests use deterministic initial geometry, one thread, fixed seeds,
high-frequency output, and no bundling, steric contacts, or unrelated springs
unless explicitly listed.

### Test A: passive tether only

System:

- one actin and one fixed myosin;
- AM pair permanently weak;
- no active force;
- no AA crosslinks.

Experiments:

1. Initialize at the preferred transverse separation.
2. Initialize with transverse displacement.
3. Translate actin axially across the myosin while retaining overlap.

Expected:

- transverse separation relaxes to `d0`;
- there is no systematic axial drift from the tether;
- axial tether force is negligible compared with transverse force;
- the bond ruptures only on overlap loss or excessive separation.

### Test B: one free actin, no force balance

System:

- one free actin;
- one fixed myosin;
- explicit AM KMC or a forced weak-to-strong event;
- no AA crosslink and no external load.

Expected sequence:

1. Power stroke creates positive active strain and force.
2. The actin moves in the productive direction.
3. Active strain and force relax toward zero.
4. The actin is not stationary while unloaded.

Primary metrics:

- displacement after each power stroke;
- peak and relaxation time of active force;
- mean shortening velocity;
- fraction of time with force above the ACTN2 threshold.

This is the critical negative control. A strong AM attachment alone must not be
interpreted as sustained load.

### Test C: one isometrically constrained actin

System is identical to Test B, but fix the actin position or connect it to an
explicit very stiff external anchor.

Expected:

- actin velocity is zero or near zero;
- power-stroke strain persists;
- active force approaches the spring/stall limit;
- strong-state lifetime follows the load-dependent rate;
- the result is independent of the passive tether stiffness over its stable
  range.

This is the positive isometric control.

### Test D: external force-clamp force-velocity curve

System:

- one actin and one fixed myosin ensemble, or enough statistically independent
  AM heads to produce steady motion;
- opposing external force fractions `0, 0.25, 0.5, 0.75, 1.0, 1.25` of the
  measured stall force;
- at least 20 replicates per condition.

Measure steady actin velocity, active force, attachment duty ratio, and state
lifetimes. Expected:

- shortening velocity decreases monotonically with hindering load;
- velocity is near zero at the measured stall force;
- super-stall load produces no productive shortening;
- KMC hazards remain consistent with instantaneous force-bin exposures.

Do not prescribe the target force-velocity curve in this explicit-head test.
Use it to determine whether the selected cross-bridge rates reproduce one.

### Test E: minimal two-actin ACTN2 system

System:

- two actins with controlled polarity;
- one existing ACTN2 bond;
- one independently controlled myosin/head ensemble on each actin;
- fixed myosin backbones so motor reactions have a defined mechanical sink.

Run the following conditions:

| Condition | Motor configuration | AA crosslink | Expected behavior |
|---|---|---:|---|
| E0 | motors off | on | no active load; status 1 |
| E1 | one motor active | on | actin moves; transient AA load; little sustained status 2 |
| E2 | two motors, same-direction forces | on | pair translates; low opposed drive; status 1 |
| E3 | equal opposing motors | off | actins move independently; status 2 impossible |
| E4 | equal opposing motors | on | high transmitted load; low relative velocity; status 2 |
| E5 | equal opposing motors, both actins pinned | on | isometric positive control; status 2 |
| E6 | unequal opposing motors | on | drift plus partial load; intermediate status-2 duty |

For E1-E6, sweep:

- passive AM tether stiffness;
- ACTN2 stiffness;
- motor stall force or power-stroke size;
- ACTN2 `F_on/F_off` thresholds.

The purpose is to distinguish stationarity caused by real force balance from
stationarity imposed by a code path.

### Test F: orientation controls

Repeat E2-E4 for:

- antiparallel actins with correct myosin halves;
- parallel actins;
- one reversed actin polarity;
- correct axes but incorrect myosin-half assignment;
- misaligned actins over a controlled angular sweep.

Record the active force direction, opposed drive, transmitted load, and status.
These tests determine which orientation gate is needed independently of the
passive spring direction.

## 9. Stationarity diagnostics

For every minimal run, report:

```text
v_actin_parallel
v_actin_relative_to_myosin
v_i - v_j for the AA-linked pair
AA extension rate
system center-of-mass velocity
active AM force
passive AM tether force
AA transmitted force
net force residual
status-1 and status-2 occupancy and dwell times
```

Define a dimensionless stationarity score relative to the unloaded control:

```text
S_stationary = 1 - RMS(v_relative) / RMS(v_free)
```

and report it together with force, rather than classifying a trajectory as
stationary from displacement alone.

For the two-actin system, separate:

- internal stationarity: small relative velocity and AA extension rate;
- absolute stationarity: small center-of-mass velocity.

E2 should expose the distinction: a pair may preserve its internal geometry
while translating as a unit.

## 10. Statistical validation

For AM KMC, compute exposure-integrated expected events:

```text
N_expected = sum_over_steps_and_bonds p_event(f, dt)
```

Compare with observed counts by state, load bin, and transition type. Also test
survival curves against the cumulative hazard rather than comparing only mean
dwell times.

For AA KMC, use explicit exit reasons to separate:

- stochastic KMC unbinding;
- geometric rupture;
- capacity removal;
- status 2-to-1 reclassification.

Required outputs:

- observed/expected hazard ratios with confidence intervals;
- weak-state branch fractions;
- force-binned strong-state survival;
- force-binned AA status-2 survival;
- censoring counts and reasons.

## 11. Numerical validation

Repeat representative minimal conditions at:

```text
dt = 0.5e-5, 1e-5, 2e-5
```

Check convergence of velocity, force, occupancy, dwell statistics, and
stationarity. Sweep tether and active stiffness to identify stable values. The
tether must not require a stiffness that makes the Langevin integrator unstable.

Run deterministic resume tests across a checkpoint inside a strong AM dwell and
inside a status-2 AA dwell.

## 12. Acceptance criteria

Before a full A=2000, M=1000 run:

1. All force action-reaction tests pass to numerical tolerance.
2. The active AM force never points in the nonproductive direction.
3. The proximity tether produces no meaningful axial propulsion or braking.
4. Test B moves and unloads; it does not remain falsely stationary.
5. Test C is stationary and force-bearing.
6. The force-clamp velocity decreases monotonically and reaches a reproducible
   stall regime.
7. E0-E3 do not create sustained false-positive status 2.
8. E4-E5 create sustained force-dependent status 2.
9. Same-direction pair translation is distinguished from opposed internal
   loading.
10. AM and AA KMC observed/expected event ratios are statistically consistent
    with one.
11. Results are stable under a factor-of-two time-step change.
12. Resume and uninterrupted deterministic runs agree.

## 13. Implementation sequence

### Phase 1: mechanics separation

- document polarity conventions;
- split passive tether and active force routines;
- add unidirectional active strain and stall cap;
- add unit tests for forces and orientation;
- implement Tests A-C.

### Phase 2: KMC and force-velocity validation

- route KMC load dependence through active force only;
- add exposure-based diagnostics;
- implement Test D;
- tune only motor parameters needed to recover sensible force-velocity behavior.

### Phase 3: ACTN2 load definition

- implement all three candidate load modes;
- add explicit AA exit reasons;
- implement Tests E-F;
- choose the load definition from false-positive and force-balance behavior.

### Phase 4: integration

- complete checkpoint/HDF5 support;
- run small stochastic many-particle tests;
- run A=2000, M=1000 only after the minimal acceptance criteria pass;
- compare against legacy `v_am`, legacy motor-identity status, and signed-AA
  projection controls.

## 14. Full-system first comparison

Use identical initial conditions and seeds for:

1. legacy velocity plus legacy status;
2. current explicit AM plus signed AA projection;
3. separated tether/crossbridge plus motor-only ACTN2 load;
4. separated tether/crossbridge plus transmitted-only load;
5. separated tether/crossbridge plus opposed-and-transmitted load.

Track:

- myosin `Phi_3D`;
- actin and myosin nematic order;
- status-1/status-2 bond counts and dwell times;
- actin velocity conditioned on AA status;
- AM active force and AA transmitted load;
- status-2 false positives in same-direction or unbalanced configurations;
- assembly, turnover, and force-velocity statistics.

This comparison should determine whether force-stabilized assembly arises from
true load-bearing constraints rather than motor occupancy, passive spring sign,
or an imposed velocity slowdown.

## 15. Multi-head validation completed 2026-07-20

- The fixed one-pair KMC test used 35 independent slots for three runs of two
  time units. Measured/expected ratios were 1.0016 for attachment, 1.0024 for
  total weak exit, 0.9909 for the power-stroke branch, and 1.0076 for strong
  detachment.
- In the tilted mechanics control, eight slots reduced the 20-degree initial
  AM angle to 5.55 +/- 0.72 degrees at time 0.2; the 35-slot torque-off control
  remained at 19.71 +/- 2.64 degrees. The unconstrained noisy sweep is retained
  alongside the controlled probe.
- With an intentionally persistent AA scaffold and no explicit memory, the
  status-2 fraction was 0.622, 0.992, and 0.998 for 1, 8, and 35 slots,
  respectively. Parallel heads therefore provide substantial emergent
  persistence.
- A persistence-only `am_status2_memory` bug in the load-based classifier was
  found and fixed. Five milliseconds of state memory raised the one-slot
  status-2 fraction to 0.749 but changed little for 8 or 35 slots.
- With the real catch-lifetime law restored, exposure/break lifetime estimates
  were 0.00279, 0.0535, and 0.149 for 1, 8, and 35 slots. State memory alone did
  not alter those lifetimes because it correctly does not create residual load.
