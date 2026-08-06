# AM–ACTN2 mechanics and heartbeat correction plan

Date: 2026-07-22

## Objective

Make the assembly model distinguish four coupled but different quantities:

1. calcium/thin-filament permission for myosin attachment;
2. explicit weak and strong actomyosin (AM) head states;
3. force actually transmitted through each ACTN2–actin binding arm;
4. status-1 and load-stabilized status-2 ACTN2 links.

The corrected model must allow actins and myosins to move. A loaded ACTN2 link
may translate with a bundle or creep under sub-stall motor force; neither zero
velocity nor motor occupancy alone defines status 2.

## Implementation status (2026-07-22)

Implemented:

- `arm_load_productive` is now the production-default ACTN2 classifier.
- Every local ACTN2 link is classified from its two signed mechanical arm
  loads plus productive, polarity-correct strong-head support.
- Status-2 memory can bridge a short support interruption only while the link
  continues to carry bilateral arm load; it cannot create status 2.
- Material-head association and retention now use exact 3-D distance to the
  finite actin segment. Pair-level overlap is no longer multiplied into the
  rate of an already spatially eligible head.
- Forced AM rupture events have a backward-compatible v2 stream containing the
  exact head connector/surface distance and association/retention thresholds.
- Sparse ACTN2 trajectories have aligned per-link connector, arm-load,
  productive-support, gate, lifetime, and legacy-tension records, plus an
  extended transition-event stream.
- Minimal and large-system renderers now include synchronized activation,
  weak/strong AM, status-1/status-2 ACTN2, and order/mechanics dashboards. The
  myosin-head renderer uses the trajectory bare-zone fraction rather than the
  former one-third hard-coded geometry.
- Production-like minimal variants cover perfect, annealing, one/no motor,
  parallel, and wrong-polarity configurations. The ideal configuration seeds
  six geometrically reachable weak heads on each arm and uses two spatial
  ACTN2 links; controls rupture invalid heads through ordinary geometry logic.

Validation completed:

- all 17 CTest targets pass;
- the exposure/event estimate of the weak chemical timescale is 2.46 ms versus
  2.5 ms expected from `k_weak_off+k_powerstroke=400/s`; geometric ruptures are
  treated as censored rather than misclassified as chemical exits;
- the aggregate weak-to-strong outcome fraction is 0.495 versus 0.5 expected;
- 83 ACTN2 detachments were observed versus 79.7 expected by integrating the
  saved instantaneous off-hazard (ratio 1.04);
- no NaNs or force/torque caps occurred in the 18-run minimal factorial;
- one/no-motor, parallel, and wrong-polarity controls never entered status 2;
- a 5 ms support-memory control increased perfect-config status-2 occupancy
  from about 5% to 14%, but did not robustly maintain or anneal the freely
  diffusing two-actin geometry.

The last result is a model finding, not a numerical failure. In the current
scaffold-free minimal system, freely mobile thick filaments can follow their
actins and unload, and transient loss of all ACTN2 links permits rotational
diffusion. Increasing the ACTN2 perpendicular stiffness from 50 to 300 did not
rescue the behavior. Do not hide this by immobilizing myosin. The next physical
decision is whether the minimal load-bearing unit should include thick-filament
network connectivity, a mature multi-step ACTN2 state, or calibrated longer
baseline ACTN2 residence; these controls should precede a new large pilot.

Throughout the code, plots, and documentation, use these unambiguous names:

- `AM weak head` and `AM strong head` for cross-bridge-cycle states;
- `AA status 1` and `AA status 2` for ACTN2-link states.

Do not call both of them merely “weak/strong catch bonds.”

## Required physical definitions

### Spatially accessible AM heads

Each conserved head has a material anchor on one myosin half. For head `h` and
actin `i`, compute the exact distance from the head's surface anchor to the
finite actin segment. A new attachment is eligible only when:

- the head lies on the polarity-correct myosin half;
- its anchor is outside the bare zone;
- a finite actin point is within the 3-D capture reach;
- the binding azimuth is reachable;
- the actin is thin-filament activated.

The available head population must emerge from the explicit anchors. In the
parallel continuum limit it should satisfy

```text
E[N_accessible] ~= N_half * L_effective_overlap / L_head_bearing_half.
```

Remove the pair-level overlap factor from each already spatially eligible
head's rate, or demonstrate analytically that it is not being counted twice.
A head chooses among all reachable actins by competing rates and can own only
one attachment.

Retention uses the same finite-segment geometry with a larger rupture reach
than association reach. A head must detach after losing its own productive
axial overlap; overlap elsewhere on the same half is insufficient.

### AM force–velocity mechanics

Retain a finite internal working-stroke coordinate for each strong head:

```text
epsilon_h = target_h - current_AM_coordinate
F_h       = min(k_motor * max(epsilon_h, 0), F_stall)
d target_h/dt = v0 * max(0, 1 - F_h/F_stall)
```

This is the internal stroke velocity, not prescribed filament velocity.
Actual actin and myosin velocities must continue to arise from reciprocal
forces, anisotropic mobilities, all other mechanical constraints, and thermal
noise.

Calibrate the emergent force–velocity curve with force-clamp simulations. Do
not assume that the internal linear law guarantees the desired ensemble curve.
Check for double counting between load-slowed internal stroke kinetics and
load-slowed detachment.

### Per-link ACTN2 mechanics

An ACTN2 link stores fixed material coordinates `(s_i, s_j)` while bound:

```text
p_i = r_i + s_i u_i
p_j = r_j + s_j u_j
e   = minimum_image(p_j - p_i)
```

The existing extensional and perpendicularity energies continue to generate
equal-and-opposite forces and torques. They resist relative axial sliding but
must not prevent common bundle translation. Detachment and rebinding may pick
new material sites and thereby permit annealing.

### Arm-resolved catch load

For each physical ACTN2 link, expose the total interface force on its first
actin, `F_AA_i`; the second receives `-F_AA_i`. Project the interface reactions
onto the experimentally favorable directions along the two actins, using one
explicitly documented force-on-actin versus force-on-ACTN2 convention:

```text
f_arm_i = max(0, favorable_projection_i(F_AA_i))
f_arm_j = max(0, favorable_projection_j(-F_AA_i))
```

Unit tests, rather than comments alone, must establish the signs for both
actin polarities.

Keep the scalar extensional tension

```text
T_extension = max(0, k_AA * (link_length - preferred_length))
```

as a structural diagnostic. It is not by itself the directional catch load.

### Status-2 gate and lifetime

Define productive opposed motor support using strong heads on the correct
myosin half for both actins. The motor condition identifies the intended
sarcomeric source of load; it does not substitute its magnitude for the force
borne by an individual ACTN2 arm.

For the first implementation, use

```text
f_effective = min(f_arm_i, f_arm_j)

status 1 -> status 2 iff
    f_arm_i > f_arm_on
    and f_arm_j > f_arm_on
    and productive_opposed_motor_support

status 2 -> status 1 iff
    f_arm_i < f_arm_off
    or f_arm_j < f_arm_off
    or productive support is absent beyond the configured memory interval.
```

Evaluate the existing catch-lifetime law with `f_effective`, not total myosin
force and not connector extension alone. Preserve on/off hysteresis and allow
status-2 memory only to bridge a short interruption; memory cannot establish
status 2.

The later, more physical option is an explicit two-arm link:

```text
k_link_off = k_arm_off(f_arm_i) + k_arm_off(f_arm_j),
```

because loss of either ACTN2–actin interface destroys the crosslink. Bank this
as a model refinement after the coarse bilateral-load version is validated.

For multiple ACTN2 links between the same actins, compute arm forces separately
from each link's mechanics. Do not give every link the entire pair-level motor
force. Mechanical force balance should produce the load partition naturally.

## Heartbeat model and required outputs

### Activation sequence

For `am_cycle_period > 0`, retain the global calcium-like transient and local
thin-filament state `q_i`. The expected causal sequence is:

```text
calcium pulse
  -> q_i rises
  -> AM weak occupancy rises
  -> AM strong occupancy rises after weak-to-strong conversion
  -> motor and ACTN2 arm loads rise
  -> AA status-2 occupancy rises
  -> calcium and q_i decay
  -> new AM attachment is suppressed
  -> strong AM heads detach
  -> arm load and AA status-2 occupancy decay.
```

Heads remain independent KMC objects. “Collective” means correlated ensemble
occupancy driven by shared calcium, thin-filament feedback, load-dependent duty
ratio, and optional thick-filament recruitment; it does not mean synchronous
state transitions.

Existing production scripts currently use tonic activation
`am_cycle_period=0`. Heartbeat diagnostics and later heartbeat pilots must set
an explicit nonzero period.

### Diastolic correction audit

Quantify whether load-slowed strong-head detachment retains heads into the next
cycle. If diastolic occupancy does not converge close to baseline, compare:

1. the present load-slowed detachment law;
2. a calibrated minimum strong-head off-rate;
3. an explicit nucleotide/post-stroke exit state;
4. a calcium-permissive attachment gate plus a separate biochemical
   relaxation pathway for already strong heads.

Do not choose among these by visual appearance. Require multi-cycle occupancy,
dwell-time, and force–velocity agreement.

### HDF5 diagnostics

Retain the existing exact AM head-state arrays and add backward-compatible
per-link mechanical fields. Each saved local ACTN2 link should provide:

- link identity, actin indices, status, and material coordinates;
- connector vector and length;
- extensional tension;
- favorable arm loads `f_arm_i` and `f_arm_j`;
- `f_effective` used by the lifetime law;
- productive opposed motor-support value and Boolean gate;
- instantaneous lifetime/off-rate and status-transition reason.

Global time series derived from state, not independently maintained counters:

- calcium activation `C(t)`;
- mean and distribution of `q_i(t)`;
- `N_AM_weak(t)` and `N_AM_strong(t)`;
- detached-ON and OFF/SRX head counts when thick feedback is enabled;
- `N_AA_status1(t)` and `N_AA_status2(t)` counted once per physical link;
- summed and quantile ACTN2 arm loads;
- mean productive AM force and mean relative AM sliding velocity;
- actin and myosin nematic order for larger systems.

Use transition-event logs for KMC rate and dwell-time tests; saved-frame counts
alone cannot resolve millisecond weak-head lifetimes unless output cadence is
correspondingly fine.

## Common animation dashboard

Refactor `animate_minimal_mechanics.py`, `animate_sarcomere_3d.py`, and the
standalone AM demo to use a shared time-series/dashboard helper. Every new
animation must include a current-time cursor and at least:

1. 3-D geometry:
   - actin polarity;
   - myosin head-bearing halves and bare zone;
   - weak and strong AM heads with different colors;
   - AA status-1 and status-2 links with different colors.
2. Activation panel:
   - calcium `C(t)`;
   - mean `q_i(t)`.
3. AM occupancy panel:
   - `N_AM_weak(t)`;
   - `N_AM_strong(t)`.
4. ACTN2 occupancy panel:
   - `N_AA_status1(t)`;
   - `N_AA_status2(t)`.

Minimal mechanics movies should add relative sliding, bundle center-of-mass
motion, and selected-link arm loads. Large-system movies should add actin and
myosin nematic order when space permits.

Compute time series once before rendering. Do not reread or rescan the entire
HDF5 file for every frame. Open HDF5 with file locking disabled and `locking=False`.
The renderer must read the simulator's bare-zone fraction and filament lengths
from metadata.

## Diagnostic test systems

### 1. ACTN2 force/sign unit tests

Two parallel and two antiparallel actins, no motors and no diffusion.

- Rest-length perpendicular link: zero arm load and zero restoring force.
- Transverse extension: positive extensional tension and reciprocal forces.
- Positive and negative axial slip: force reverses and always opposes slip.
- Swap link ownership: energy invariant, forces reverse, arm-load labels swap.
- Both actin polarities: favorable arm projections match the documented truth
  table.
- Compression: no catch-stabilizing load.

Extend `local_aa_link_test.cpp` to assert axial restoring-force signs, not just
nonzero energy and torque.

### 2. Status-2 truth table

One pre-existing ACTN2 link with mechanics and motor support independently
clamped.

| Arm load | Productive opposed support | Expected |
|---|---|---|
| low | absent | status 1 |
| high | absent | status 1 |
| low | present | status 1 |
| high on both arms | present | status 2 |
| high on one arm only | present | status 1 |
| high | same-direction or wrong-half motors | status 1 |

Sweep around `f_arm_on` and `f_arm_off` to verify hysteresis, timestep
independence, and memory that preserves but never creates status 2.

### 3. Relative sliding versus common translation

Two mobile actins with zero, one, and two ACTN2 links.

- Equal common forces: bundle translates, relative displacement stays small,
  and a rest-length link remains status 1.
- Equal opposed forces: relative sliding initially rises, ACTN2 restoring force
  develops, and status 2 appears only after bilateral arm load crosses threshold.
- Unequal opposed forces: simultaneous bundle translation and status-2 load.
- One-sided motor force: ACTN2 resists and drags the second actin, but the
  productive-opposed gate keeps the assembly-specific status at 1.
- Remove one motor after status 2: arm load relaxes and status returns to 1
  after the configured memory interval without deleting the ACTN2 link.

Record `v_1`, `v_2`, bundle velocity, relative velocity, link geometry, both arm
loads, motor forces, and status.

### 4. Per-head AM-overlap sweep

Translate both actin polarities through both myosin halves, the central bare
zone, and beyond the filament.

- Count exact spatially accessible heads versus effective overlap.
- Compare with `N_half * L_effective/L_head_bearing_half`.
- Verify the correct staircase for discrete anchors.
- Require zero accessible heads in the bare zone and on the wrong half.
- Require new attachment to vanish at capture loss and existing attachment to
  survive only within the larger retention reach.
- Repeat across timestep and neighbor-cell boundaries.

With several actins around one myosin, verify conservation: one head can bind
only one actin, and candidates compete without duplicating the head pool.

### 5. Force–velocity clamps

One AM pair with both particles mobile, then with one fixed only as an explicit
control. Use one head and multiple heads.

- Apply a series of opposing force clamps from zero to above stall.
- Measure actual relative filament sliding velocity, internal stroke velocity,
  force, strain, dwell time, and completed stroke fraction.
- Require maximum sliding at zero load, decreasing velocity with load, and no
  productive backward force after overshoot.
- Compare no ACTN2, one ACTN2 link, and a two-link bundle.
- Test whether the current stroke law plus force-dependent detachment
  double-counts load sensitivity.

### 6. Chemical heartbeat without mechanics

Fixed productive geometry only for this rate-validation test. Run at least ten
cycles with fine event logging.

- thin feedback off/on;
- thick feedback off/on;
- synergy off/on;
- no strong heads;
- tonic activation control;
- several `dt` values.

Compare measured attachment, weak-to-strong, weak-off, strong-off, OFF/SRX
recruitment, and folding hazards with implemented rates. Measure phase lag and
cross-correlation among calcium, `q`, weak occupancy, and strong occupancy.

### 7. Mobile minimal heartbeat catch system

Two myosins, two oppositely polarized mobile actins, and one pre-existing ACTN2
link, all with physical diffusion.

Render several heartbeat cycles and plot:

- calcium and `q`;
- AM weak and strong counts;
- AA status 1 and status 2 counts;
- each ACTN2 arm load;
- relative actin sliding and bundle translation.

Controls:

- no motors;
- motors but no ACTN2 catch stabilization;
- one productive motor side only;
- tonic calcium;
- catch law with productive gate removed;
- isometric versus freely mobile mechanics.

Expected result: periodic weak occupancy followed by strong occupancy, followed
by delayed arm loading and status 2. Diastole must reliably return strong AM and
AA status-2 occupancy toward baseline rather than latch.

### 7A. Production-faithful minimal catch matrix: persistence and annealing

Repeat the original minimal catch configuration after all corrections above are
enabled. This is not a reduced chemical toy: use the same feature switches and
parameter path as the next production pilot, with only particle number and box
size reduced.

Production-faithful means:

- two mobile actins and mobile myosins with physical anisotropic mobility and
  full thermal noise;
- dynamic detached/weak/strong AM heads;
- explicit spatial head anchors, transported myosin material frames, binding-
  time azimuth selection, off-center active forces, passive flat-bottom tether,
  finite force–velocity stroke, and load-dependent detachment;
- polarity-correct, per-head finite-overlap eligibility, including the shared
  bare-zone constant;
- thin-filament activation and strong-head feedback;
- tension-dependent thick-filament OFF/SRX recruitment in the combined-feedback
  condition, with thin-only and thick-only controls retained;
- local multi-ACTN2 links with fixed material sites, perpendicularity mechanics,
  site exclusion, production association law, bilateral arm-load catch law,
  productive-opposed support gate, hysteresis, and configured status-2 memory;
- filament sterics, production interaction ramp, production timestep, force and
  torque caps, output cadence, and HDF5 event logging;
- a nonzero heartbeat period for the heartbeat series, plus a tonic-activation
  comparison using the otherwise identical setup.

Do not immobilize either actin. A fixed-particle run is permitted only as an
explicit mechanical sign control and is not evidence of persistence.

#### Geometry and motor-presence factorial

Use the following core matrix. “Correct” means the actin overlaps the
head-bearing half that recognizes its polarity; bare-zone and wrong-half
overlap do not count.

| Case | Actin relation | Myosin support | Expected status-2 behavior |
|---|---|---|---|
| A | antiparallel, correct | two correct myosins | recurrent/persistent status 2 under bilateral transmitted load |
| B1 | antiparallel, correct | remove myosin 1 | status 1 only after memory; one-sided drag may remain |
| B2 | antiparallel, correct | remove myosin 2 | status 1 only after memory; mirror of B1 |
| C | antiparallel, correct | remove both myosins | status 1 prebundling/turnover only |
| D+ | parallel, both `+` polarity | two geometrically nearby myosins | never status 2 |
| D- | parallel, both `-` polarity | two geometrically nearby myosins | never status 2 |
| E1 | antiparallel | one actin assigned to the wrong myosin half | never bilateral status 2 |
| E2 | antiparallel | both actins on wrong halves | no productive strong heads and no status 2 |
| F | antiparallel | one or both actins wholly in the bare zone | inaccessible heads contribute zero; no status 2 |
| G | antiparallel | correct half but zero effective axial overlap | no new AM attachment and no status 2 |

Cases B1 and B2 must be statistically symmetric. Cases D–G are truth-table
failures if even one new status-2 episode is established; short memory may only
finish an episode established before the perturbation.

Also test one shared bipolar myosin positioned between the two actins as a
separate geometry. Its two halves may support both actins only if each explicit
head independently satisfies the correct-half, correct-polarity, and finite-
overlap test. This must not be treated as equivalent to simply multiplying the
force from one AM pair.

#### Persistence protocol

Initialize the ideal antiparallel geometry at preferred AM and AA spacing with
pre-existing ACTN2 link(s). Begin every link in status 1; status 2 must be
established dynamically by bilateral favorable arm load.

Run at least ten heartbeat cycles after relaxation. Report both molecular-link
persistence and pair-level structural persistence:

- survival time of every ACTN2 link identity;
- fraction of time each link spends in status 1 and status 2;
- number and duration of temporary complete-unbinding gaps;
- fraction of time the actin pair has at least one ACTN2 link;
- time to irreversible loss of pair association, using a predeclared gap and
  distance criterion;
- AM weak/strong occupancy and AA status-1/status-2 occupancy by cycle;
- both ACTN2 arm loads, motor-support gate, extension, obliquity, and lifetime;
- relative sliding, common translation, AM overlap, spacing error, and nematic
  misalignment;
- cycle-to-cycle recovery and diastolic convergence.

Temporary replacement of one ACTN2 molecule by another does not count as loss
of the structural configuration. Conversely, a continuously bound but unloaded
status-1 link does not count as evidence of catch stabilization.

Required controls use the same Brownian-noise streams where possible:

- catch lifetime enabled versus force-independent matched-lifetime control;
- productive-support gate enabled versus recorded-but-not-applied gate;
- heartbeat versus tonic activation;
- thin-only, thick-only, combined feedback, and all synergy off;
- one, two, and zero allowed ACTN2 links per pair where applicable.

#### Annealing protocol

Start without a mature status-2 link, or with status-1 prebundle links only.
Perturb one coordinate at a time and then selected combinations:

- AM/AA radial spacing at `0.8`, `0.9`, `1.0`, `1.1`, and `1.2` times the
  preferred value;
- actin nematic misalignment of `5`, `10`, `20`, and `30` degrees, including an
  out-of-plane component;
- axial material-site mismatch and relative axial offset;
- unequal AM overlap on the two actins;
- small displacement across a periodic boundary;
- a link initially oblique but within capture reach.

Define annealing before inspecting results. Success requires recovery into the
same tolerance basin as the perfect configuration for spacing, nematic angle,
and productive AM overlap, followed by dynamically established status 2.
Measure:

- first-passage time into the correct basin;
- probability of annealing within a fixed number of heartbeats;
- number of ACTN2 detach/rebind events required;
- whether the pair becomes kinetically trapped in a wrong material-site
  registration;
- changes in relative sliding and arm load before and after annealing;
- persistence after the first successful annealing event.

The parallel and wrong-polarity systems may bundle transiently through status-1
ACTN2 links, but must not be reported as successfully annealed unless they
actually reach the target antiparallel productive configuration. The current
rod model cannot reverse an actin's intrinsic polarity; a parallel pair should
therefore remain a negative control rather than being expected to flip polarity.

#### Replication and convergence

For each deterministic sign/geometry test, use zero noise first. For stochastic
persistence and annealing, use at least 20 seeds for the compact minimal system,
then repeat the key A, B1/B2, D, and E cases at three timesteps and at least two
neighbor-cell arrangements. Report survival and annealing distributions, not
only representative trajectories.

#### Required movies

Render representative full-duration movies for:

- perfect two-myosin persistence;
- successful annealing from angular plus axial perturbation;
- one-myosin removal after an established status-2 episode;
- zero-myosin turnover;
- parallel-actin negative control;
- wrong-polarity/wrong-half negative control.

Every movie uses the common dashboard: calcium and `q`, AM weak/strong counts,
AA status-1/status-2 counts, arm loads, relative sliding, and a current-time
cursor. Use identical colors, camera, axis limits, and time scaling across the
factorial so visual comparisons are meaningful.

### 8. Scaffold-free six-actin lattice heartbeat

One mobile myosin with six surrounding actins at preferred AM/AA spacing. Start
from both perfect and perturbed geometries. Run tonic and pulsed activation.

Measure AM accessibility by crown/half, radial and azimuthal occupancy,
hexagonal order, actin/myosin nematic order, relative bundle slip, ACTN2 arm
loads, and cycle-resolved status populations. Confirm that heartbeat turnover
does not destroy a perfect lattice and that imperfect configurations can anneal.

### 9. Full-scale pilot

Only after tests 1–8 pass, run a short 2000-actin/1000-myosin pilot with a
multi-cycle heartbeat and three replicates. Acceptance checks:

- no NaNs or exploding forces;
- force/torque cap fraction below the established tolerance;
- exact head conservation and overlap-scaled accessibility;
- KMC hazards and dwell distributions match rates;
- both AM weak/strong and AA status-1/status-2 populations are present;
- reproducible cycle phase relationships;
- diastolic convergence without progressive latching;
- plausible memory, output size, and runtime;
- no loss of nematic or lattice order solely from visualization/output cadence.

## Implementation order

1. Add arm-load helpers, sign tests, and per-link HDF5 diagnostics without
   changing status behavior.
2. Implement the bilateral arm-load plus productive-support status classifier;
   preserve the old modes as explicit compatibility controls.
3. Route the catch lifetime through `f_effective` and validate per-link KMC.
4. Implement per-head finite-overlap eligibility and correct the thin-feedback
   accessible-head denominator.
5. Add force-clamp and sliding diagnostics; calibrate before changing kinetic
   parameters.
6. Add heartbeat aggregation and diastolic convergence tests.
7. Refactor both animation pipelines around the common dashboard and regenerate
   all minimal movies with AM and AA population plots.
8. Run chemical, mobile minimal, and six-actin heartbeat systems.
9. Launch the full-scale pilot only after the preceding acceptance tests pass.

Each step should land with unit tests and a compatibility flag until its paired
diagnostic demonstrates the expected behavior. Do not tune thresholds or rates
using full-scale assembly order before the minimal force and KMC tests pass.

## Rejected perfect-configuration stabilization attempt (2026-07-22)

The first production-physics minimal test exposed two implementation defects
and one missing mechanical degree of freedom:

1. The explicit-head AM path passed `am_cutoff` to a centerline-distance
   geometry routine. This omitted the actin radius, myosin radius, and preferred
   surface gap, sometimes set the crosslinkable segment to zero despite visible
   AM occupancy, and prevented replacement ACTN2 links from associating.
2. Configurations with more than two local ACTN2 molecules reused the same two
   axial sites. Local sites are now distinct and evenly distributed over the
   instantaneous actin overlap.
3. Connector perpendicularity alone leaves an in-plane rotational hinge in
   3D. Mature local links now inherit the legacy nematic AA energy
   `0.5*kappa_aa*(1-(u_i.u_j)^2)`. The pair-level `kappa_aa` is divided over
   simultaneous molecular links so changing occupancy does not multiply the
   total intended alignment stiffness. The encounter-state association energy
   still receives zero nematic penalty.

The minimal diagnostic now reconstructs weak/strong AM populations from the
conserved per-head state and ownership arrays, including frame zero. The old
derived feature counters are zero before the first force pass and had made the
initial movie frame incorrectly appear unoccupied.

The attempted perfect-system diagnostic parameters were four distinct effective
ACTN2 sites, `base_lifetime=0.2 s`, `am_status2_memory=5 ms`, and the existing
`kappa_aa=50`. These were useful for isolating AA rotational mechanics but are
not accepted as stable defaults. `run_minimal_catch_factorial.sh` has been
returned to its prior two-link, 50 ms status-1, and 0.5 ms-memory defaults.

Across three 200 ms (four accelerated heartbeat) replicas:

- the minimum simultaneous ACTN2 population was two in every replica;
- zero-link frame fraction was exactly zero;
- individual links still detached 2--3 times and replacement links formed;
- status 2 appeared in 37.4--52.2% of saved frames;
- mean strong-head population was 12.4--20.6;
- final actin--actin misalignment was 0.304--0.624 degrees and maximum
  misalignment was 2.13--2.90 degrees;
- force- and torque-cap counts were zero and every saved value was finite.

Results and the dashboard movie are in
`/scratch/users/xsli/sarcomere/minimal_catch_diagnostics_20260722/perfect_stable_life200ms`.
These AA-only metrics did not establish persistence of the full configuration.
Inspection of the movie and full coordinates showed that the myosin separation
collapsed from 2.66 to 0.44 micrometres: each free myosin was loaded on only one
half and therefore translated toward the central crosslinked actins. One
intended AM pair ultimately lost productive overlap. The actins stayed linked
and nearly parallel, and status 2 was transiently nonzero, but the four-body
AM--AA configuration fell apart. This result is rejected.

The corrected acceptance test must include myosin--actin axial registration,
productive overlap on each intended half, myosin separation, pair identity,
and absence of half-crossing in addition to AA occupancy and angle. With all
particles mobile, the full-physics minimal cell must also contain the opposing
load that a real bipolar thick filament receives from its other half. Use a
periodic two-Z-disc unit (two central actins of interest plus counter-actins on
the unused halves) or an explicitly labeled isometric clamp control; do not
call the unbalanced two-actin/two-myosin cell stationary.

## Whole-cell follow-up: transient weak-head friction

`minimal_catch_balanced` now supplies the two periodic counter-actins needed
to occupy both halves of each freely mobile thick filament. No rod is fixed;
production anisotropic mobility and thermal noise remain enabled. The explicit
head test also required `dt=0.5 us`: the previous `1 us` step triggered caps
during peak recruitment, whereas the smaller step produced zero force and
torque caps without changing physical rates.

The remaining failure was diastolic decoupling. The flat-bottom AM tether is
radial and intentionally permits axial sliding. Each weak head now also stores
its attachment-time axial coordinate and applies the reversible passive force

`F_weak = k_weak * (x_bind - x_current) * u_actin`,

where `k_weak = am_weak_stiffness_ratio * am_motor_k`. KMC turnover releases
and renews the coordinate, so this is transient sliding friction rather than a
permanent preferred overlap. A configurable basal weak-association gate allows
weak encounters at `q=0`; weak-to-strong conversion is multiplied by `q`, so
basal weak binding alone cannot produce motor force or status 2.

Three one-heartbeat replicas with basal weak access 0.1 passed the expanded
checks: minimum intended overlap 0.195--0.296 um, maximum myosin-separation
change 0.031--0.079 um, final actin misalignment 0.49--0.83 degrees, strong
heads cleared in diastole, and zero force/torque caps. A three-heartbeat run at
0.1 nevertheless lost two overlaps in the third cycle. Basal access 0.2 kept
all overlaps for three cycles (minimum 0.091 um, zero caps), but deep-diastolic
weak occupancy was about 13% and overlap drift continued monotonically.

Therefore the present scaffold-free cell is **metastable, not indefinitely
stable**. Do not choose basal weak access solely to preserve the lattice.
Long-time heartbeat stability requires a separately justified counterconstraint
(for example M-line/titin mechanics) or experimentally calibrated slower
crowded mobility. That missing physics must not be hidden by fixed myosins or
strong-head latching.
