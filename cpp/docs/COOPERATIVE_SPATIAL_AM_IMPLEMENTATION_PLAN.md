# Cooperative, Spatially Resolved Actin-Myosin Implementation Plan

## Objective

Build the next actin-myosin model in a sequence that keeps the mechanisms
identifiable:

1. Refactor head ownership from actin-myosin-pair-owned to myosin-half-owned.
2. Give every head a spatial anchor and a single actin target, and validate
   steric organization without an imposed lattice potential.
3. Add a per-thin-filament activation variable `q_i` and feedback from strongly
   bound heads.
4. Test calcium pulses with thin-filament cooperativity alone.
5. Add tension-dependent thick-filament `OFF/SRX <-> ON` recruitment.
6. Combine thin- and thick-filament feedback and test their synergy.
7. Only after those layers pass their isolated tests, study assembly and ACTN2
   catch-bond persistence.

This version deliberately excludes M-band, Z-disc, titin, and MyBP-C scaffolds.
It also excludes a direct sixfold-ordering potential. Any transverse order must
come from spatially resolved AM attachments, conserved heads, filament sterics,
and motor mechanics. This makes the test stringent: motors may promote local
hexagonal coordination, but the model is not guaranteed to reproduce the
long-range mature lattice without the omitted structural components.

## Core modeling rules

- A thick filament may bind several actin filaments simultaneously.
- An individual myosin head may bind at most one actin at a time.
- The number of available heads belongs to the thick filament and is
  independent of how many actins happen to be nearby.
- Calcium/thin activation controls access to actin binding sites; it does not
  directly apply force.
- Strongly bound heads can further activate their thin filament.
- Thick-filament tension recruits heads from `OFF/SRX` to `ON`; passive AM
  tether force and geometric overlap alone do not count as thick tension.
- Strong heads generate active force at spatial head anchors. Passive AM
  tethers maintain proximity but do not act as motors.
- ACTN2 binding remains independent of AM occupancy. Status 2 remains a loaded
  state of an already-existing ACTN2 crosslink.

## Scope and explicit non-goals

For all tests in this plan, set `k_mm=0`, `titin_k=0`, and
`k_bundle_max=0`. Do not add MyBP-C-like tethers, M/Z anchors, lattice sites,
or a term proportional to `-|psi_6|^2`.

The pair-owned dynamic-head representation will be replaced rather than kept
as a second production implementation. Pre-refactor dynamic-AM checkpoints do
not need to resume into the new head layout. The loader should fail clearly on
the old schema instead of silently inventing or duplicating heads. Legacy files
remain readable by analysis tools where practical.

## Stage 1: myosin-half-owned, spatially resolved heads

### 1.1 State ownership

Replace

```cpp
am_head_ensembles[actin_index][myosin_index][pair_head_index]
```

with one fixed-capacity pool per myosin half:

```cpp
myosin_heads[myosin_index][half][head_index]
```

Suggested head state:

```cpp
struct AMHeadState {
    AMChemicalState state;       // OFF, detached_ON, weak, strong
    int8_t half;                 // -1 or +1
    int bound_actin;             // -1 when unbound
    double anchor_axial;         // material coordinate on thick filament
    double anchor_azimuth;       // material angle around thick filament
    double actin_material;       // attachment coordinate on bound actin
    double target_axial;         // active spring target after power stroke
    uint64_t attach_step;
};
```

Use `am_heads_per_half`, not `am_heads_per_pair`. Head conservation is an exact
invariant:

```text
N_OFF + N_detached_ON + N_weak + N_strong = 2*N_myosin*N_heads_per_half.
```

Changing the number of neighboring actins must not change this total.

### 1.2 A material frame for each thick filament

The current axial unit vector cannot represent rotation about the filament's
own axis. Store an orthonormal material frame
`(e1_m, e2_m, u_m)` or equivalently `u_m` plus a transported azimuthal phase.
Head-root positions are

```text
r_h = r_m + s_h*u_m
            + R_m*(cos(phi_h)*e1_m + sin(phi_h)*e2_m).
```

Update the frame under the same finite rotation used for the thick filament.
Include the torque component about `u_m`; otherwise azimuthally resolved heads
cannot rotate around the thick-filament backbone. Re-orthonormalize after each
update and test that repeated rotations do not drift the frame.

The initial crown layout should be parameterized rather than embedded in the
force routine. A suitable first coarse graining is three azimuthally spaced
head pairs per crown, with successive crowns staggered in azimuth and placed
along the productive region of each half. Exact crown spacing and rotation can
be calibrated later without changing ownership or KMC logic.

### 1.3 Candidate selection and exclusive binding

For each detached, ON head, build candidates from actins satisfying:

- the correct myosin half and productive actin polarity;
- sufficient axial overlap;
- radial distance below the capture cutoff;
- angular accessibility from the head's azimuthal anchor;
- `q_i > 0`.

Use a rate weight

```text
k_attach(h,i) = k_on * q_i * p_ON(m,half)
                * g_overlap * g_polarity * g_distance * g_azimuth.
```

Because ON/OFF is represented explicitly, `p_ON` is either implicit in the
number of `detached_ON` heads or used only in a mean-field diagnostic, not
multiplied a second time. Choose among candidates with a competing-hazard draw.
One head can accept one actin. Several different heads can select the same
actin, and one thick filament can distribute its heads among multiple actins.

Iterate KMC by myosin owner, not by actin-myosin pair. This makes exclusivity
natural and avoids OpenMP races. Accumulate actin forces and torques in
thread-local arrays as the current mechanics already does.

### 1.4 Head mechanics

Store the actin material coordinate at attachment. The weak and strong passive
tethers connect the spatial head root to that material point and provide only
the intended proximity elasticity. The strong state additionally supplies the
polarity-fixed active axial spring and stall cap already defined in the
explicit-crossbridge model.

Apply equal and opposite forces at both attachment points. For each head,

```text
tau_m = (r_h - r_m) x F_on_myosin
tau_a = (r_attach - r_a) x F_on_actin.
```

The full torque, including the component around the thick-filament axis, must
be retained. Geometric rupture detaches only that head, not every head between
the same two filaments.

### 1.5 Sterics

Retain finite-radius rod-rod repulsion for actin-actin, actin-myosin, and
myosin-myosin pairs. Before using it for ordering, verify:

- force is equal and opposite;
- overlap energy/force is continuous near the cutoff;
- parallel rods do not pass through each other at the production timestep;
- repulsion is not disabled merely because a filament has an ACTN2 bond;
- AM attraction cannot overwhelm the force cap and create numerical overlap.

No angular sixfold bias is added. A noisy ring of six actins around one thick
filament should be an emergent test of radial AM attraction plus actin sterics,
not an imposed target.

### Stage-1 tests and gate

Add unit or minimal-system tests for:

1. Exact head conservation with 0, 1, 2, and 6 candidate actins.
2. Exclusive binding: no head has two actin targets.
3. Correct half/polarity eligibility for every sign combination.
4. Correct head-root position after translation, tilt, and axial rotation.
5. Action-reaction and off-center torque for a single spatial head.
6. Analytical D/W/S dwell and branching statistics at fixed `q_i=1` and all
   heads forced ON; retain the existing 5% statistical tolerance.
7. Checkpoint/resume equality with bound heads distributed across actins.
8. One myosin plus six initially mis-spaced parallel actins, comparing sterics
   on/off and spatial anchors versus centerline anchors.

Do not start thin-feedback implementation until conservation, single binding,
and KMC statistics pass.

## Stage 2: per-thin-filament activation and strong-head feedback

### 2.1 Thin activation variable

Add one bounded activation variable `q_i in [0,1]` per actin filament. The
global calcium waveform becomes an input `C(t)`, not the attachment gate by
itself. A minimal deterministic coarse-grained equation is

```text
dq_i/dt = k_C * H_C(C(t)) * (1-q_i)
          + k_SB * H_SB(s_i) * (1-q_i)
          - k_close * q_i,
```

where

```text
H_C(C)  = C^n_C/(K_C^n_C + C^n_C)
H_SB(s) = s^n_SB/(K_SB^n_SB + s^n_SB)
s_i     = N_strong_on_actin_i / N_accessible_heads_i.
```

Use a stable bounded update and clamp only as a final numerical safeguard.
`N_accessible_heads_i` is a geometric normalization, so feedback is not larger
just because global myosin density is larger. Strong occupancy, rather than
passive tether force, drives thin feedback. Initially keep feedback local to
the same thin filament; do not propagate `q` through ACTN2 bonds.

`q_i` gates detached-ON to weak attachment. A falling `q_i` should stop most
new attachment but should not forcibly erase already-weak or strong heads;
those decay through their physical detachment rates.

### 2.2 Calcium pulse protocol

Replace the current single global `am_activation` diagnostic with a recorded
calcium drive `C(t)` plus per-actin `q_i(t)`. Support a prescribed waveform and
a tabulated waveform; retain the current difference-of-exponentials pulse as
the first protocol.

For the thin-only phase, force every thick-filament head to remain ON and turn
off tension-dependent recruitment. Run identical calcium pulses with:

- strong-head feedback on;
- strong-head feedback off (`k_SB=0`);
- no strong heads (`k_powerstroke=0`), keeping the calcium response intact;
- no calcium, to test spontaneous activation;
- repeated pulses followed by a long diastolic interval.

Required behavior:

- `q_i` remains in `[0,1]`;
- calcium alone produces a finite rise and relaxation;
- strong heads amplify and/or prolong activation when feedback is enabled;
- removing strong heads removes most motor-dependent cooperativity;
- at zero calcium, activation and strong occupancy return reproducibly toward
  baseline and do not ratchet upward over successive beats.

The thin-only gate passes only when pulse statistics and diastolic decay are
stable under `dt`, `dt/2`, and `dt/4`.

## Stage 3: tension-dependent thick-filament recruitment

### 3.1 OFF/SRX and ON states

Extend the head state graph to

```text
OFF/SRX <-> detached_ON -> weak -> strong -> detached_ON.
```

Only detached heads may return to OFF in the first implementation. This avoids
destroying a force-bearing attachment merely because the backbone recruitment
variable changed.

Use exact KMC probabilities. A minimal tension-sensitive rate pair is

```text
k_OFF_to_ON(T_h) = k_recruit_0 * exp(+T_h/F_recruit)
k_ON_to_OFF(T_h) = k_fold_0    * exp(-T_h/F_fold),
```

with exponent caps or a bounded alternative to prevent numerical overflow.
The zero-tension rates set the resting SRX/ON balance.

### 3.2 Thick-half tension

Compute tension separately for the two halves. It must be the summed
nonnegative force transmitted toward the bare zone by productive strong heads,
not:

- the maximum load of one AM pair;
- net force on the whole freely translating thick filament;
- passive transverse tether force;
- overlap or head occupancy alone.

A first definition is

```text
T_mh = sum_{strong heads on half h}
       max(0, hindering axial force transmitted toward the bare zone).
```

Record both this scalar and the signed force sum used to construct it. Verify
with hand-built geometries that assisting or unloaded heads do not recruit the
thick filament.

### 3.3 Force-clamp and load tests

Hold calcium and `q_i` fixed and disable thin feedback. Test:

- prescribed half-filament tension levels;
- an isometric actin that retains motor strain;
- a freely sliding/low-load actin that unloads the heads;
- zero active force with passive tether force present;
- release from high load to zero load.

The measured OFF-to-ON and ON-to-OFF hazards must match the input rates within
5% in prescribed-tension tests. ON fraction and OFF-to-ON recruitment must rise
monotonically with thick-filament tension and relax after load release.

## Stage 4: combine the two feedback loops

The combined loop is

```text
Ca -> q_i -> AM attachment -> strong heads -> q_i
                                  |
                                  v
                         thick-half tension
                                  |
                                  v
                          OFF/SRX -> ON heads
                                  |
                                  +----> more AM attachment.
```

Run the factorial comparison under the same calcium waveform and matched
zero-load baseline availability:

| Thin strong-head feedback | Thick tension feedback | Label |
|---|---|---|
| off | off | neither |
| on  | off | thin only |
| off | on  | thick only |
| on  | on  | combined |

Repeat under isometric and low-load conditions. Report peak and integrated
`q`, ON fraction, weak/strong occupancy, active force, and relaxation time.
Quantify interaction rather than calling any increase synergy by inspection:

```text
Delta_synergy = R_both - R_thin_only - R_thick_only + R_neither,
```

for each response `R`. Positive feedback may sharpen or prolong activation,
but the combined system must still return to baseline during diastole.

### Decisive controls

The required controls are:

1. Same calcium pulse, synergy on versus off.
2. Same calcium pulse under isometric versus low-load conditions.
3. Thin-feedback only versus thick-feedback only.
4. No strong heads: thin activation should lose most motor-dependent
   cooperativity.
5. Force clamp: `OFF -> ON` recruitment should rise with thick-filament
   tension.
6. Diastole: activation must reliably decay rather than latch.

No assembly conclusion is allowed before all six controls have been run.

## Stage 5: scaffold-free transverse organization

Use the combined model with ACTN2 disabled first.

### 5.1 Minimal systems

1. One thick filament and six parallel thin filaments initialized as a noisy
   ring, with small axial and angular perturbations.
2. Three thick filaments and one shared thin filament to test the reciprocal
   coordination geometry.
3. A small periodic transverse patch at approximately two thin filaments per
   thick filament.

Compare spatial heads plus sterics against:

- spatial heads with sterics disabled;
- centerline head anchors;
- the same total head pool artificially restricted to one actin;
- motors present but active force disabled;
- thin feedback only, thick feedback only, both, and neither.

Measure local transverse vectors after projecting perpendicular to the local
myosin axis. For thick filament `m`,

```text
psi6_m = (1/N_m) * sum_a exp(6*i*theta_ma).
```

Also report:

- actin coordination around each thick filament;
- thick coordination around each actin;
- angular-gap and radial-distance distributions;
- species-resolved transverse structure factors;
- head occupancy by actin neighbor and by myosin half;
- axial and transverse force balance;
- order persistence after motor activation is switched off.

`|psi6|` alone is insufficient: six actins on a ring around one myosin is not a
long-range double-hexagonal lattice. Coordination and structure-factor results
must accompany it.

The expected result without structural scaffolds may be transient local order
rather than a persistent global lattice. Treat that as a mechanistic result,
not an implementation failure, provided the isolated mechanics and KMC tests
pass.

## Stage 6: assembly and ACTN2 catch-bond persistence

Only after Stages 1-5 pass, enable AA binding and the status-2 load classifier.
Do not allow `q_i`, head availability, or AM attachment to read ACTN2 status
directly. ACTN2 affects motors only through mechanics: constrained actins retain
strain and therefore alter strong-head force, lifetime, and thick tension.

Use independent, geometrically compatible AM and AA centerline targets:

```text
d_AM* = r_M + r_A + am_surface_gap
d_AA* = rho_AA_AM * d_AM*
```

The physiological default is `rho_AA_AM=0.92`, with focused controls at
`0.8`, `0.9`, and `1.0`. A direct `aa_preferred_centerline` override is allowed
for calibration, but its effective ratio must be recorded. The AA binding
cutoff is `d_AA* + aa_capture_margin`; it must not be smaller than the spring
minimum. This spacing is an assembly compatibility hypothesis, not a claim
that ACTN2 permanently crosslinks the mature A-band hexagon.

Run in this order:

1. Combined AM model, ACTN2 disabled.
2. ACTN2 with constant lifetime.
3. Status-2 catch law with `lifetime_coeff=0`.
4. Status-2 catch law across the focused lifetime-coefficient range.
5. Matched occupancy/residence/turnover constant-lifetime controls.
6. Pulsatile versus matched-mean tonic activation.

For each condition measure assembly order, time to assembly, status-1/status-2
counts, AA residence and turnover, load distributions, status-2 interruptions,
AM strong-head residence, `q_i`, thick ON fraction, and transverse coordination.
During diastole, ask separately whether:

- thin and thick activation switch off;
- passive/catch-stabilized structure persists;
- status 2 decays when mechanical load disappears;
- the network progressively latches across beats.

Catch-bond persistence is supported only if matched constant-lifetime controls
cannot reproduce the same history dependence and if activation itself still
decays normally.

## Diagnostics and HDF5 schema

Increment the state schema and record enough information to reconstruct every
feedback loop.

Per head/event:

```text
myosin_id, half, head_id, bound_actin, old_state, new_state, reason,
anchor_axial, anchor_azimuth, actin_material, strain, active_force,
q_actin, half_tension, calcium, step.
```

Per actin per frame:

```text
q_i, calcium, accessible_head_count, weak_count, strong_count,
active_AM_force, passive_AM_force, ACTN2 status counts.
```

Per myosin half per frame:

```text
N_OFF, N_detached_ON, N_weak, N_strong, half_tension,
signed_axial_force, number_of_distinct_bound_actins.
```

Per myosin per frame, also store the material frame or axial phase needed to
reconstruct head roots. Checkpoint/resume must preserve head targets, bindings,
chemical states, frames, `q_i`, and pending transition events exactly.

Add automated audits for:

- head conservation and exclusive targets;
- measured versus dictated KMC hazards;
- thin-activation rise/decay rates;
- OFF/ON recruitment versus force clamp;
- force/torque action-reaction residuals;
- pulse-to-pulse baseline drift;
- coordination, `psi6`, and structure factor.

## Expected code changes

- `include/sarcomere.h`: replace pair-owned containers; add head ownership,
  `q_i`, material frames, half tension, and OFF/ON parameters.
- `include/am_motor_kinetics.h`: extend states/rates and add bounded thin- and
  thick-regulation helpers suitable for unit testing.
- `src/sarcomere.cpp`: change AM updates to myosin-owned KMC, candidate draws,
  spatial mechanics, feedback updates, tension accumulation, save/resume.
- `src/interaction.cpp` and `include/interaction.h`: spatial head-root tether,
  material-coordinate attachment, torque, and steric validation helpers.
- `src/components.cpp`: initialize and update thick-filament material frames.
- `src/main.cpp`: replace `am_heads_per_pair` with `am_heads_per_half`; add
  calcium, thin-feedback, and thick-feedback parameters and explicit toggles.
- `src/h5_utils.cpp` and `include/h5_utils.h`: new checkpoint and diagnostic
  schema with a hard version check.
- `tests/`: ownership, material-frame, spatial-force, thin-feedback,
  force-clamp recruitment, and checkpoint tests.
- `analysis/`: pulse audit, head-conservation audit, load-recruitment curves,
  synergy factorial analysis, and transverse-order analysis.

## Implementation order and commit gates

1. **Ownership/schema:** conserved myosin-half pools and checkpoint format.
2. **Spatial mechanics:** material frames, exclusive candidates, spatial
   forces/torques, and steric tests.
3. **Thin loop:** `C(t) -> q_i`, strong-head feedback, calcium-pulse controls.
4. **Thick loop:** half tension and `OFF/SRX <-> ON`, force-clamp controls.
5. **Combined loops:** full factorial, load dependence, and diastolic decay.
6. **Local organization:** scaffold-free minimal and periodic tests.
7. **Assembly/ACTN2:** focused production runs and matched-lifetime controls.

Each gate requires passing unit tests, numerical convergence, dictated-rate
audits, and a checkpoint/resume test. If a gate fails, fix that isolated layer
before adding the next feedback loop.
