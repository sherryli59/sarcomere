# ACTN2 directional loading and filament-alignment design

## Purpose

This document specifies how the particle model should represent the coupling
among actin polarity, bipolar myosin geometry, actin--myosin (AM) overlap,
motor-generated force, actin--actin (AA) crosslinks, and directional ACTN2
catch-bond stabilization.

The central modeling decision is to keep **crosslink existence** separate from
**crosslink loading state**:

- An ACTN2 crosslink can bind actin without myosin being strongly attached.
- Myosin activity can load an existing crosslink, but does not chemically
  create the crosslink.
- Directional catch stabilization occurs only when the force transmitted
  through both ACTN2 arms has the experimentally favored sign.

This document began as a design and implementation plan. The first isolated
implementation slice described below is now available behind an opt-in runtime
mode; the legacy classifier remains the default.

## Implementation status

Implemented:

- `actin.direction` was verified to point from the modeled barbed end
  (`left_end`) to the pointed end (`right_end`).
- Unit-tested helpers define signed ACTN2 arm load, bilateral load, hysteretic
  status classification, nematic AM alignment, myosin half, and productive
  polarity score.
- `--aa_status_mode=aa_force_projection` classifies only an existing AA bond
  from signed force on both arms. New AA bonds always start as status 1.
- `--aa_status2_force_on` and `--aa_status2_force_off` set the hysteretic
  transition thresholds.
- Because the current AA potential contains passive geometric prestress, the
  first implementation limits each favorable raw AA arm reaction by the motor
  load available on that actin. Dynamic AM uses summed strong-head spring
  strain; legacy AM uses its existing overlap/load proxy. This prevents passive
  AA relaxation from creating status 2 in the absence of motor engagement.
- Force-projection mode uses favorable bilateral AA load in the status-2
  lifetime law and does not use the fixed AM status-memory window.
- The existing `legacy_motor_identity` mode is unchanged and remains default.
- Per-actin overlap/orientation and AA-load summaries are written as trajectory
  features. AA bond/status changes receive a separate
  `/catch_bond/load_transitions` force audit.
- New trajectories store classifier mode and thresholds under `/meta`; resume
  refuses mismatched settings.

Not yet implemented:

- overlap-dependent rates beyond the existing geometric binding-ratio gate;
- per-AM-pair overlap datasets;
- relaxation of the current directional/antiparallel formation geometry;
- overlap-scaled angular alignment;
- pN calibration or an experimentally fitted ACTN2 kinetic law;
- continuum closure.

## Scope and non-goals

The immediate target is a minimal mechanistic model of thin-filament sorting,
Z-body/Z-disc consolidation, and force-dependent maturation. It is not a full
model of de novo sarcomerogenesis. In particular, the first implementation will
not explicitly model actin polymerization, FHOD3, nonmuscle myosin II, titin,
troponin/tropomyosin biochemistry, myosin-head super-relaxed states, or
individual ACTN2 molecules within a multi-protein Z-disc.

The design must support two later reductions:

1. A reduced particle model in which motor load is estimated from overlap,
   duty fraction, and relative sliding.
2. A continuum model in which local actin polarity and myosin density generate
   active stress, and ACTN2 stabilization depends on projected tensile stress.

## Biological interpretation of the discrete states

AA bond existence and AA loading status must remain distinct variables.

| Bond flag | Status | Interpretation |
|---|---:|---|
| 0 | 0 | No ACTN2-mediated AA crosslink exists. |
| 1 | 1 | An ACTN2 crosslink exists but is slack, compressed, weakly loaded, loaded in the unfavorable direction, or favorably loaded on only one arm. |
| 1 | 2 | The same crosslink exists and both ACTN2 arms carry sufficient load in the favorable direction. |

A `1 -> 2` or `2 -> 1` transition is a change in mechanical/kinetic state, not
binding or unbinding. Only AA association and rupture change the bond flag.

Strong AM states may generate the forces responsible for a status transition,
but neither strong motor occupancy nor the identity of the attached myosins is
part of the molecular definition of status 2.

## Coordinate and polarity conventions

### Actin

Define one documented actin polarity vector everywhere:

\[
\mathbf p_i = \text{unit vector from the barbed }(+)\text{ end toward the
pointed }(-)\text{ end of actin }i.
\]

With this convention:

- the expected force of myosin on actin is approximately in the
  `+p_i` direction;
- favorable ACTN2 loading is toward `+p_i`;
- two actins correctly arranged across a Z-disc satisfy
  `p_i dot p_j approximately -1`;
- their barbed-end regions are close to the crosslink, while their pointed ends
  extend away in opposite directions.

Before behavioral changes are made, the existing meaning of
`actin.direction`, `left_end`, and `right_end` must be verified. If the stored
direction has the opposite convention, either change it globally with a
careful migration or introduce explicit helpers such as
`actin_pointed_direction(i)` and `actin_barbed_end(i)`. Do not scatter sign
corrections through interaction code.

### Bipolar myosin

Treat a myosin filament axis `m_m` as an unoriented/nematic axis with two
head-bearing halves. For a contact point `r_c` and myosin center `r_m`, define

\[
u_{im}=(\mathbf r_c-\mathbf r_m)\cdot\mathbf m_m,
\qquad h_{im}=\operatorname{sign}(u_{im}).
\]

For an ideally assembled sarcomere,

\[
\mathbf p_i \simeq -h_{im}\mathbf m_m.
\]

The corresponding polarity diagnostic is

\[
q_{im}^{\rm productive}=-h_{im}\,\mathbf p_i\cdot\mathbf m_m.
\]

`q approximately 1` indicates the final productive polarity for that half,
`q approximately -1` reversed polarity, and `q approximately 0` a transverse
contact. This score should initially be diagnostic rather than a hard binding
gate: incorrectly oriented actins must be able to interact with motors if
motor activity is to sort them dynamically.

## Actin--myosin overlap

### What overlap represents

AM overlap measures the number of motor heads that could engage an actin. It is
a proxy for force-generating **capacity**, not for actual load. Large overlap
can carry little load during free sliding, while a small stalled overlap can
carry substantial load.

### Geometric overlap

For each geometrically neighboring actin--myosin pair, compute a head-accessible
overlap length `L_geom`. The computation should:

1. Respect periodic boundaries.
2. Restrict the myosin interval to the appropriate head-bearing regions.
3. Project the actin segment onto the local myosin axis.
4. Intersect the projected actin interval with the myosin head-bearing
   interval or intervals.
5. Return zero when lateral separation is beyond the interaction range.

The first version may use a single axial projection if the filaments are close
to parallel. A later exact version can integrate the length over which the
distance between the two finite segments is within head reach.

Define an effective accessible overlap

\[
L_{im}^{\rm eff}=L_{im}^{\rm geom}
g_{\rm radial}(d_\perp)
g_{\rm axis}(|\mathbf p_i\cdot\mathbf m_m|).
\]

Use a smooth nematic alignment factor, for example

\[
g_{\rm axis}=\exp\left[-\frac{(1-|\mathbf p_i\cdot\mathbf m_m|)^2}
{2\sigma_{AM}^2}\right].
\]

The absolute dot product is deliberate: attachment compatibility depends on
axial alignment, while polarity selection should emerge from bipolar motor
geometry and directed stepping.

### Overlap-dependent AM attachment

The weak attachment rate should eventually scale with accessible overlap:

\[
k_{\rm on,im}=k_{\rm on}^0 a(t)
\frac{L_{im}^{\rm eff}}{L_{\rm ref}},
\]

with an explicit upper bound if necessary for timestep stability. Equivalently,

\[
N_{im}^{\rm available}=\rho_{\rm head}L_{im}^{\rm eff}
\]

can set the number of available coarse-grained heads or the maximum AM bond
occupancy.

The implementation should first record overlap without changing rates. The
overlap-dependent on-rate is a separate behavioral change with its own toggle
and matched-occupancy controls.

## From motor engagement to active force

### Full particle model

In the mechanistic particle model, actual motor force should come from the sum
of strains in strong AM attachments:

\[
\mathbf F_i^{AM}=\sum_{m\in\mathrm{strong}(i)}
k_{\rm head}\Delta x_{im}\,\mathbf p_i
\]

plus any explicitly modeled transverse spring contributions. Overlap controls
how many heads can bind; attachment state and strain determine actual force.

### Reduced particle proxy

If explicit head strain is removed, estimate active load using overlap, strong
duty fraction, and a force--velocity factor:

\[
F_{im}^{\rm proxy}=
\rho_{\rm head}L_{im}^{\rm eff}
d_{im}^{\rm strong}f_{\rm head}
\left(1-\frac{v_{im}}{v_0}\right)_+.
\]

This distinguishes freely sliding and mechanically stalled overlap. A bare
`force proportional to overlap` rule is permitted only as an explicitly named
minimal control.

### Motor alignment and torque

A single motor attachment should behave primarily as a point spring. Multiple
attachments distributed over a longer overlap naturally create a torque. If an
effective angular AM interaction is retained, its strength should scale with
engaged overlap and occupancy, for example

\[
\kappa_{im}^{AM}=\kappa_0^{AM}
\frac{L_{im}^{\rm eff}}{L_{\rm ref}}N_{im}^{\rm strong}.
\]

It should favor nematic axial alignment
`abs(p_i dot m_m) -> 1`, not directly impose the final polarity sign. A direct,
strong polarity-aligning potential would encode the desired sarcomere rather
than allow it to self-assemble.

## Actin--actin crosslink geometry and alignment

### Crosslink formation

An AA crosslink should form according to ACTN2-compatible geometry:

- proximity of the two crosslinkable/barbed-end regions;
- finite ACTN2 span;
- steric accessibility;
- optionally, a broad and independently controlled angular compatibility.

Initial formation should not require strong myosin attachment or strict
antiparallel orientation. This allows parallel and misoriented short-lived
contacts to act as substrates for force-dependent selection.

### AA forces and torques

The AA interaction routine should expose the connector force at each ACTN2
attachment point before distributing it into center-of-mass force and torque.
Let

\[
\mathbf F_{i\leftarrow AA}^{(ij)},\qquad
\mathbf F_{j\leftarrow AA}^{(ij)}
\]

be the forces applied by crosslink `ij` to the two actins. Their translational
parts should obey Newton's third law, while torques follow from the attachment
lever arms and any explicitly modeled ACTN2 orientational elasticity.

A very strong direct antiparallel AA angular potential should not be the
baseline. Antiparallel order is mechanistically interesting only if it persists
when direct angular preference is absent or modest. Required controls are
`kappa_aa = 0`, weak structural preference, and the current strong preference.

## Force-based status-2 classification

### Signed arm loads

For every existing AA bond, calculate the force transmitted through each ACTN2
arm and project it onto the corresponding pointed-end direction. The exact sign
depends on whether the routine exposes force on actin or force on ACTN2, so the
final helper must be validated against an explicit two-filament pulling
geometry.

After that convention is fixed, define favorable signed loads

\[
f_i^- = \mathbf F_{i}^{\rm arm}\cdot\mathbf p_i,
\qquad
f_j^- = \mathbf F_{j}^{\rm arm}\cdot\mathbf p_j.
\]

Use only the force attributable to this AA connector. Do not use total actin
force, scalar whole-filament tension, motor count, or myosin identity.

For the current fixed-rest-geometry AA potential, keep both a raw connector
reaction and a motor-limited effective load:

\[
f_{i,\mathrm{eff}}^-=\min\left[\max(0,f_{i,\mathrm{raw}}^-),
F_i^{AM}\right].
\]

The same definition applies to arm `j`. This temporary cap is necessary because
the existing axial AA energy can produce passive prestress even with motor
attachment disabled. Dynamic AM supplies `F_i^{AM}` from summed favorable
strong-head spring strain; legacy AM supplies its documented overlap proxy. A
future attachment-specific AA rest configuration should replace this cap if
external, non-myosin loading needs to stabilize ACTN2.

The bilateral load is

\[
f_{ij}^{\rm bilateral}=\min(f_i^-,f_j^-).
\]

The minimum is appropriate because both ACTN2 arms must be favorably loaded.

### Minimal classifier with hysteresis

For an existing bond:

\[
1\rightarrow2 \quad\text{if}\quad
f_i^->f_{\rm on}\;\text{and}\;f_j^->f_{\rm on},
\]

\[
2\rightarrow1 \quad\text{if}\quad
f_i^-<f_{\rm off}\;\text{or}\;f_j^-<f_{\rm off},
\]

with `f_off < f_on`. Hysteresis suppresses timestep-scale force noise without
inventing a fixed temporal memory. Thresholds must eventually be calibrated in
pN.

Antiparallel orientation may be retained as a diagnostic and, initially, as a
safety gate. The target implementation should make it redundant: a genuinely
bilateral favorable load through a barbed-end crosslink should itself select the
correct geometry.

### Later kinetic classifier

Replace the hard status threshold with a force-dependent internal state only
after the minimal classifier is validated:

```text
bonded-unloaded <-> bonded-load-stabilized -> unbound
```

with rates

\[
k_{1\to2}(f_i^-,f_j^-),\quad
k_{2\to1}(f_i^-,f_j^-),\quad
k_{\rm off,1},\quad
k_{\rm off,2}(f_i^-,f_j^-).
\]

This is the appropriate place to incorporate a fitted two-state or
multi-pathway ACTN2 catch-bond model.

## Catch-bond lifetime law

The first force-based implementation may retain the current phenomenological
law for controlled comparison:

\[
\tau_1=\tau_{\rm base},
\qquad
\tau_2=\tau_{\rm base}+\lambda
\max(0,f_{ij}^{\rm bilateral}).
\]

Status 1 remains short-lived and receives no load-dependent term. The old
`stabilize_status1` behavior remains available only as a legacy control.

The linear law must be labeled phenomenological and bounded to its calibrated
force range. It should later be replaced by direction-specific experimental
off-rates. A constant `directional_base_lifetime` should not be assumed unless
zero-load ACTN2 data justify it.

## Update order and circular-dependency resolution

Force affects bond lifetime, but a bond must exist to transmit force. Use this
update order:

1. Restore the existing AA bond set and previous internal states.
2. Compute AM states, motor strains, and motor forces.
3. For every existing AA bond, compute its trial connector force from the
   current geometry.
4. Project the two arm forces and update status 1/2 using hysteresis or kinetic
   transitions.
5. Evaluate the status-dependent AA off-rate and remove ruptured bonds.
6. Apply forces from surviving AA bonds.
7. Evaluate geometrically eligible unbound AA pairs and form new bonds.
8. Initialize every newly formed bond as status 1.
9. Allow a new bond to become status 2 only on a later step after it has carried
   measurable load.

If force evaluation and integration ordering require it, use the previous
step's connector force for status classification and document the one-step lag.
Do not use strong myosin occupancy to avoid the lag.

## Role of status-2 memory

The biological baseline should use zero fixed status memory. The existing
`am_status2_memory` option may remain as a numerical sensitivity control while
the force-based classifier is developed.

If hysteresis is insufficient, prefer an explicit ACTN2 internal-state rate or
a documented low-pass force measurement over an arbitrary grace time. Any
filter time must be reported and swept.

## Continuum reduction

The particle-level measurements should determine the closure used in the
continuum model. Candidate fields are:

- actin densities `rho_plus(x)` and `rho_minus(x)` or total density `rho_a` and
  polarity `p(x)`;
- myosin density `rho_m(x)` and nematic orientation `m(x)`;
- an AM overlap/engagement field `o(x)`;
- a crosslink density `c(x)`;
- a stabilized-crosslink fraction `s(x)`;
- velocity and stress fields.

A minimal active stress is

\[
\boldsymbol\sigma^{\rm active}
=\zeta\,o(\mathbf x)\,d_{\rm strong}(\mathbf x)
\left(\mathbf p\mathbf p-\frac{1}{d}\mathbf I\right),
\]

or a two-polarity generalization. Mechanical force balance converts overlap and
activity into load. Catch stabilization then depends on the tensile stress
projected along each polarity, rather than directly on overlap:

\[
\partial_t s
=k_{1\to2}(\sigma_i^-,\sigma_j^-)(c-s)
-k_{2\to1}(\sigma_i^-,\sigma_j^-)s.
\]

The continuum closure should be fitted or checked against particle-level
conditional measurements such as
`E[AA load | AM overlap, duty fraction, relative velocity, orientation]`.

## Implementation phases

### Phase 0: convention audit and frozen legacy reference

1. Document the exact existing meanings of actin endpoints and direction.
2. Construct fixed two-actin and actin--myosin geometries that establish all
   force signs.
3. Freeze paired-seed legacy trajectories and summary statistics.
4. Add no behavioral changes in this phase.

Exit criterion: polarity, motor stepping, AA connector force, and optical-trap
loading conventions are mutually consistent and covered by tests.

### Phase 1: overlap and orientation diagnostics only

1. Implement a pure AM overlap geometry helper.
2. Compute `L_geom`, `L_eff`, myosin half, axial alignment, productive polarity
   score, strong occupancy, relative sliding velocity, and strain-based force.
3. Save these quantities without changing any transition or force rule.
4. Compare overlap with actual motor force to determine when overlap is an
   adequate proxy.

Exit criterion: overlap is continuous under small displacements, symmetric
under exchange/reflection where appropriate, periodic-boundary safe, and
predictive of engagement but not incorrectly labeled as load.

### Phase 2: overlap-dependent AM attachment

1. Add an opt-in attachment mode in which `k_on` scales with `L_eff`.
2. Preserve the current constant-on-rate mode.
3. Match mean AM occupancy between the two modes before comparing assembly.
4. Test convergence with timestep and overlap smoothing width.

Exit criterion: expected attachment counts and dwell statistics agree with
analytic fixed-geometry predictions.

### Phase 3: expose AA arm force and add force-based status mode

1. Refactor the AA interaction so connector force is available independently
   of force accumulation.
2. Add signed force projections for both arms.
3. Add a runtime status-classification mode, initially alongside the legacy
   strong-distinct-myosin classifier.
4. Initialize new bonds as status 1.
5. Add configurable `f_on` and `f_off`; default the new experimental mode to
   zero fixed memory.
6. Retain the current lifetime formula initially so classification is the only
   changed mechanism.

Suggested explicit modes:

```text
legacy_motor_identity
aa_force_projection
aa_force_projection_with_hysteresis
```

Exit criterion: the controlled status-classification tests below all pass, and
strong motor occupancy has no direct effect when AA connector load is held
fixed.

### Phase 4: alignment-mechanism audit

1. Separate AM nematic alignment from polarity sorting.
2. Make any effective AM angular strength scale with engaged overlap/occupancy.
3. Test AA angular interaction at zero, weak, and current strength.
4. Relax strict antiparallel formation gating so incorrectly oriented status-1
   crosslinks can form and turn over.
5. Quantify whether directional loading, rather than a prescribed angular
   potential, selects antiparallel organization.

Exit criterion: the intended ordering is not produced trivially by a dominant
alignment potential and is lost or altered under directional-catch controls.

### Phase 5: calibration and kinetic catch law

1. Map simulation length, time, stiffness, and force units to micrometers,
   seconds, pN/um, and pN.
2. Fit or tabulate direction-specific ACTN2 off-rates over the measured force
   range.
3. Replace the unbounded linear lifetime with the fitted kinetic law.
4. Fit AM zero-load rates and load sensitivities independently for power-stroke
   and detachment transitions.
5. Propagate uncertainty through ensemble simulations.

Exit criterion: single-bond force-clamp simulations reproduce target lifetime
distributions, not merely their means.

### Phase 6: reduced and continuum closures

1. Regress actual particle-level AM force against overlap, strong duty,
   relative velocity, and orientation.
2. Select the smallest closure that predicts conditional mean and variance.
3. Derive active stress and crosslink-state equations.
4. Compare particle and continuum phase diagrams using matched dimensionless
   groups and history-switch protocols.

## Required diagnostics and stored data

At minimum, record per active AM pair or suitable sparse summaries of:

- geometric and effective overlap;
- myosin half and contact coordinate;
- `abs(p dot m)` and productive polarity score;
- AM state, axial strain, force, and relative sliding velocity.

For every active AA bond, record:

- attachment positions and barbed-end distances;
- `p_i dot p_j`;
- connector force vector;
- the two signed arm loads;
- bilateral load;
- status before and after classification;
- status-transition reason;
- rupture probability and rupture reason;
- nearby AM overlap, strong occupancy, and motor force for each actin.

Per-frame summaries should include status-1/status-2 counts and fractions,
parallel/antiparallel AA counts, bilateral-load distributions, AM overlap and
force distributions, and correlations with `phi_3d`.

Checkpoint/resume data must preserve any hysteresis state, filtered force, or
kinetic ACTN2 state exactly.

## Controlled geometry tests

### AM overlap and polarity

1. Zero axial overlap gives zero accessible overlap.
2. Full parallel overlap gives the expected segment length.
3. Reversing actin polarity preserves nematic overlap but reverses productive
   polarity score.
4. Moving a contact across the myosin center changes the half sign correctly.
5. Lateral separation smoothly suppresses effective overlap.
6. Periodic images give the same result as equivalent unwrapped geometry.

### AA status classifier

| Configuration | Expected status |
|---|---:|
| Existing antiparallel bond, both arms favorably loaded | 2 |
| Existing antiparallel bond, slack | 1 |
| Existing antiparallel bond, compressed | 1 |
| Only one arm favorably loaded | 1 |
| Parallel pair under a common-mode force | 1 |
| Two strong distinct myosins but no transmitted AA tension | 1 |
| Favorable externally generated AA load without myosin | 2 |
| Newly formed AA bond | 1 |
| Existing status-2 bond unloaded below `f_off` | 1 without forced unbinding |
| Geometric AA bond loss | 0 and memory/internal state cleared |

### Stochastic kinetics

1. Fixed-force AA dwell times follow the configured exponential or fitted
   multi-state distribution.
2. Status-1 lifetime is independent of load in status-2-only mode.
3. Status-2 lifetime responds to the bilateral load, not the larger single-arm
   load.
4. Results converge with timestep under `1-exp(-k dt)` transitions.
5. Resume and uninterrupted trajectories agree statistically and preserve
   deterministic behavior where requested.

## Validation experiments

Use paired seeds and matched occupancy where possible.

1. **Classifier comparison:** legacy motor-identity status versus AA-force
   status, with identical motor and lifetime parameters.
2. **Overlap ablation:** constant AM on-rate versus overlap-dependent on-rate,
   matched for total AM occupancy.
3. **Load proxy comparison:** actual head strain versus overlap-only versus
   overlap-times-stall proxy.
4. **Alignment audit:** `kappa_aa = 0`, weak, and current; constant versus
   overlap-scaled AM alignment.
5. **Direction controls:** signed bilateral load, unsigned load, reversed
   polarity sign, and nondirectional stabilization.
6. **Memory controls:** zero, hysteresis only, and prior fixed-memory values.
7. **History switches:** assemble at intermediate catch strength and switch to
   high strength; compare with disordered high-strength starts.
8. **Constant-lifetime controls:** match occupancy, mean residence time, and
   turnover rate, not only nominal lifetime.
9. **Force calibration:** single-AA force clamps across the experimental force
   range before full assembly sweeps.

Primary outcome measures are `phi_3d`, polarity registration, Z-disc spatial
coherence, AA status composition, bond turnover, AM occupancy and force, and
the fraction of status 2 events that truly carry bilateral favorable load.

## Acceptance criteria for the mechanistic claim

The model supports a claim of directional, force-stabilized self-assembly only
if all of the following hold:

1. Status 2 is predicted from local AA force, not directly from motor identity.
2. Favorable load on both arms is necessary for stabilization.
3. Overlap controls available motor engagement but does not create stabilization
   in a freely sliding, unloaded configuration.
4. Antiparallel order persists with weak or absent imposed AA angular alignment.
5. Reversing or removing load directionality reduces or qualitatively changes
   assembly.
6. Matched constant-lifetime controls do not reproduce all structure and
   turnover results solely through residence time.
7. Results are robust to timestep, force threshold/hysteresis choices within a
   calibrated range, and paired initial conditions.

## Backward compatibility and migration

- Keep legacy classification and constant AM attachment available through
  explicit runtime modes during validation.
- Write all selected modes and new parameters into trajectory metadata.
- Never infer a new mode while resuming an old file; resume must recover the
  mode stored in the checkpoint or fail clearly.
- Use new output directories/tags for force-based classification.
- Do not overwrite or reinterpret existing status datasets without versioned
  metadata.

## Main risks and mitigations

### Force-sign errors

Mitigation: fixed geometries with analytically known pulling directions and
explicit diagrams/tests before production runs.

### Circular force/status coupling

Mitigation: classify only existing bonds from trial or previous-step connector
force; initialize new bonds as status 1.

### Noisy threshold switching

Mitigation: calibrated hysteresis first, explicit state kinetics later; retain
zero-memory baseline.

### Assembly imposed by alignment potentials

Mitigation: zero/weak angular controls and torque decomposition by interaction
type.

### Overlap mistaken for load

Mitigation: save overlap, velocity, strain, and force separately; compare their
conditional relationships.

### Too many new mechanisms changed simultaneously

Mitigation: phase the work so overlap diagnostics, overlap-dependent binding,
force-based status, alignment changes, and fitted kinetics are independent
toggles and comparisons.

## Recommended first implementation slice

The smallest scientifically interpretable change is:

1. Audit and expose polarity helpers.
2. Expose per-bond AA connector forces.
3. Save signed arm-load diagnostics without changing behavior.
4. Add an opt-in force-projection classifier with hysteresis.
5. Keep the existing lifetime law and all AM rules fixed.
6. Compare it against the legacy strong-distinct-myosin classifier in controlled
   geometries and paired assembly runs.

Only after this slice is validated should AM overlap alter attachment rates or
alignment strengths. This ordering isolates the central conceptual correction:
ACTN2 status 2 is a mechanically loaded state of an existing crosslink, while
AM overlap controls the capacity to generate that load.
