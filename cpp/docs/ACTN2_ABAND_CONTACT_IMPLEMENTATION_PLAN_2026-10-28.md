# ACTN2 / A-band contact: consolidated implementation plan

Updated: 2026-10-28 (America/Los_Angeles).

**Status: proposed simulator changes, not implemented production physics.**
The geometry gates, buffered instantaneous deletion, and old full-connector
force detection are implemented. The conservative energy, revised weak-site
boundaries, energetic association, and contact-driven failure below are not.
The isotropic (`alpha=1`) predecessor has been evaluated in standalone
numerical visualizations. The weighted candidate below has not yet been
implemented or plotted.

This document supersedes the earlier energy proposal with a 100–130 nm
radial taper, and the proposal to reflect weak sites at both ends of the
200 nm barbed zone. It also supersedes the unweighted shortest-exit energy
as the axial-stop candidate: its geometry favors radial escape at the
approximately 80 nm actin shell as soon as a tip overshoots. The current
candidate weights radial escape by `alpha`, initially 3. Implementation must
distinguish accepted structure from parameters awaiting calibration.

## 1. Objective and model decisions

Prevent mature Z-disc connectors from being placed inside an A-band, allow
existing overlaps to resolve through mechanics, and distinguish molecular
failure from geometric contact. Determine whether the resulting network can
resist motor drive without excessive lateral rearrangement or link loss.

| Mechanism | Target behavior |
|---|---|
| Maturation into status 2 | Exact full-connector gate: 130 nm radius, physical myosin end planes |
| Buffered 120 nm / 10 nm region | Diagnostic occupancy/entry region; no instantaneous deletion |
| ACTN2 barbed-tip burial deletion | Retire as a separate A-band geometric death criterion |
| Soft connector–myosin contact | Gradient of one bounded anisotropic conservative energy, summed across distinct actin pairs; radial escape weighted to favor axial resistance near tips |
| Redundant connectors on one actin pair | One energy budget per pair and myosin; select maximum energy, not maximum force |
| Weak-site sliding | Energy-gradient drift plus thermal noise; absorbing physical endpoints and soft inner barbed-zone edge |
| Mature attachment sites | Fixed material coordinates until detachment, preserving current binding physics |
| Failure | Existing baseline turnover/contour rupture first; contact-specific failure added and calibrated separately |
| Hard projection | Explicitly off for every new-model comparison |
| Independent actin–myosin mechanics | Preserve and document existing settings; do not remove as a side effect of ACTN2 policy changes |

Do not add mechanism on/off switches for admission, connector geometry,
buffered-region handling, or the replacement contact law. Use pinned,
versioned implementations for staged experiments. Physical parameters and
numerical tolerances are recorded parameters, not a new switch family.

The 130 nm radius is a **structural exclusion scale**, not the myosin's
50 nm physical radius. It may cover much of a dense A-band, but ideal lattice
coverage is not proof of coverage in disordered geometry or near staggered
filament ends. Preserve separate physical-core sterics.

## 2. Current evidence and source provenance

- Active source: `$SCRATCH/sarcomere/zdisclink_20261020/cpp`.
- Durable repository: `$GROUP_HOME/sherry/sarcomere/sarcomere`.
- The current shell's older `sarcomere_am_kmc` checkout is not the target.
- Both active and durable source trees have been compiled and tested:
  44/44 active tests and 40/40 durable-repository tests passed.
- Current code already has the exact 130 nm maturation gate, immediate
  buffered mature deletion, full-connector detection with the old capped
  axial force, and parallelized deletion geometry.
- Current soft force still chooses one strongest connector per tip, is
  axial only, saturates at 20 pN, and splits its reaction equally. Changing
  detection geometry did not replace those mechanics.
- Production uses `aa_arm_state_model=false`; do not silently enable the
  arm-state model to obtain load-dependent failure.
- Existing cap levels depend on body/binding state. A myosin's roughly
  26.4 pN cap is not a bound on its raw motor drive.

Validated current executable:

```text
$SCRATCH/sarcomere/connector_parallel_20261028/sarcomere
SHA256 8167e9a5e933912125ca1cd5f684e45a74b7fcb20527356b773821f9dbd3c52d
```

Preserve that binary and the historical production binaries. Do not
overwrite trajectories, and do not copy one source tree wholesale onto the
other: their feature sets differ. Existing unrelated dirty changes must be
preserved.

Performance and comparable-geometry evidence is recorded in:

```text
$GROUP_HOME/sherry/sarcomere/sarcomere/analysis/visualizations/gaphunt-rig/
  status2_gate_20261028/results/performance_diagnostics.md
  status2_gate_20261028/results/parallel_step_benchmark.json
  status2_gate_20261028/results/comparable_connector_metrics.json
```

The paired benchmark matches death counts after reason-code remapping, but
poses and event records are not bitwise identical. Same-binary repeats also
diverged. Exact replay remains unresolved; kernel equivalence and matching
counts do not establish trajectory equivalence. Investigate this before
attributing small model differences to the new physics.

## 3. Shared full-connector geometry

Use actual material attachment sites. If `v_i` is the simulator's stored
actin direction and `s_i` its attachment coordinate in micrometers:

```math
p_i = c_i + s_i v_i,\qquad p_j = c_j + s_j v_j.
```

Do not multiply `s` by filament length. Do not silently renormalize stored
actin directions when reconstructing existing binding sites. Use the
existing coordinate convention and corresponding derivatives consistently.

Construct the short periodic connector:

```math
d = \operatorname{minimumImage}(p_j-p_i),\qquad
x(t)=p_i+t d,\quad 0\le t\le1.
```

For each relevant image of myosin center `c_m`, normalized axis `u_m`, and
physical half-length `h=L_m/2`:

```math
z(t)=(x(t)-c_m)\cdot u_m,\qquad
r(t)=x(t)-c_m-z(t)u_m,\qquad \rho(t)=|r(t)|.
```

Formation and buffered diagnostics must test radial and axial conditions
at the **same connector point**:

```math
\text{maturation excluded}\iff
\exists t:\rho(t)\le0.130\ \text{and}\ |z(t)|\le h,
```

```math
\text{buffered occupancy}\iff
\exists t:\rho(t)\le0.120\ \text{and}\ |z(t)|\le h-0.010.
```

Distances are in micrometers. Clip the segment analytically to the axial
slab, then minimize its radial quadratic on the surviving interval. This
gate needs no numerical integration. Keep the existing periodic edge-case
tests and enumerate all relevant images, not only the closest center image.

Use a conservative segment/cylinder bounding-sphere rejection first, with
the connector half-length included. A support bound for mechanical contact
uses the extended half-length `H=h+a`. Visit all potentially contacting
myosins; endpoint neighbor lists alone can miss a connector interior.

## 4. Revised conservative energy

### 4.1 Weighted penetration into a finite cylinder

Use mechanical radius `R=0.130 um` and half-length `H=h+a`. The initial
visualization uses an approach distance `a=0.050 um`; this is uncalibrated.
The exact maturation gate continues to use `h`, not `H`.

Use an anisotropic penetration coordinate, with `alpha>=1`:

```math
p_\alpha(x)=\max[0,\min(\alpha(R-\rho),H-|z|)].
```

For `alpha>1` this is a weighted escape coordinate, not the Euclidean signed
distance. The support remains within the same finite cylinder, but radial
escape costs more energy and the tip face can provide axial resistance
across a larger part of the actin shell.

The reason for this revision follows directly from geometry. At a point
buried `d>=0` inside the physical tip, the axial distance to the mechanical
end plane is `a+d`. The unsmoothed axial branch is selected when

```math
\alpha(R-\rho)>a+d,\qquad
\rho<R-\frac{a+d}{\alpha}.
```

With `R=130 nm`, `a=50 nm`, and `alpha=1`, that becomes `rho<80 nm-d`.
Thus points at approximately 80 nm or farther outward already favor radial
escape after overshoot. This failure of the unweighted candidate near the
actin shell is geometrically established, not an experimental uncertainty.

Use `alpha=3` as the first weighted candidate. Its approximate axial-branch
coverage (before smoothing) is:

| Burial d | Axial branch for rho below |
|---|---:|
| 0 nm | 113.3 nm |
| 50 nm | 96.7 nm |
| 100 nm | 80.0 nm |
| 150 nm | 63.3 nm |

It therefore extends axial capture at the 80 nm shell to about 100 nm of
burial. At rho=100 nm it extends only to 40 nm of burial, and at rho=110 nm
only to 10 nm. It does not make the whole 130 nm cylinder an axial wall,
nor guarantee reversal of the deepest inherited overlaps. Finite connectors
integrate contributions from all their points; the branch condition is
local, not a midpoint rule.

Specify a target radius and target burial before selecting alpha:

```math
\alpha(R-\rho_{target})>a+d_{target}+\text{smoothing margin}.
```

Do not simply increase alpha until it covers every inherited overshoot:
radial forces and local stiffness increase with it. The network's ability
to stop a driven myosin within the chosen axial capture region still needs
controlled loading tests.

Use the following smoothed version as the first candidate:

```math
\rho_\epsilon=\sqrt{|r|^2+\epsilon^2},\quad
z_\epsilon=\sqrt{z^2+\epsilon^2},\quad
A=\alpha(R-\rho_\epsilon),\quad B=H-z_\epsilon,
```

```math
q=\frac{A+B-\sqrt{(A-B)^2+\epsilon^2}}{2},\qquad
p_\epsilon=\max(0,q).
```

Start numerical evaluation with `epsilon=0.005 um`, `alpha=3`, and vary both
in sensitivity tests. Smoothing erodes the mechanical region slightly,
particularly near the side/end transition. It blends both force directions;
the unsmoothed branch boundary is not an exact force-direction boundary.
In particular the radial contribution is multiplied by alpha. Verify the
actual axial/radial force ratio across the target shell. Do not apply this
weighting or smoothing erosion to the exact maturation gate.

### 4.2 Bound and report the anisotropic point force

```math
p_* = F_*/k,\qquad
\Phi(p)=
\begin{cases}
\tfrac12 k p^2,&p\le p_*\\
F_*(p-\tfrac12p_*),&p>p_*.
\end{cases}
```

Point energy is `u(x)=Phi(p_epsilon(x))`. Let

```math
D=\sqrt{(A-B)^2+\epsilon^2},\quad
w_A=\tfrac12[1-(A-B)/D],\quad w_B=1-w_A.
```

For `q>0`, the force on a connector point is

```math
f(x)=\min(kp_\epsilon,F_*)
\left(w_B\frac{z}{z_\epsilon}u_m+
\alpha w_A\frac{r}{\rho_\epsilon}\right).
```

Set force to zero for `q<=0`. The weighted coordinate has
`|grad p_epsilon|<=alpha`, so the correct bounds are

```math
|f_{axial}|\le F_*,\qquad |f_{radial}|\le\alpha F_*,\qquad
|f|\le\alpha F_*.
```

Here `F_*` is the saturation slope of Phi and an axial force bound. It is
not the whole-force bound when alpha>1. For example, `alpha=3`, `F_*=30 pN`
allows up to 90 pN radially. In the unsmoothed radial quadratic branch the
effective radial stiffness is `alpha^2*k`, or 9000 pN/um at k=1000 and alpha=3;
additional curvature arises in the smooth transition. Re-evaluate timestep
and inherited-overlap response accordingly.

Do not normalize or clip the final vector f to preserve the earlier 30 pN
claim: that would generally stop it being the gradient of this energy.
If a whole-force budget F_budget is required, use `F_*=F_budget/alpha`
in Phi and explicitly accept the corresponding smaller axial force bound.
Calibrate the axial and radial budgets together.

Do not substitute `W(rho)*Phi(axial_depth)` with only the axial force capped:
that still leaves radial force proportional to an energy that grows with
axial depth. The older 30 nm taper and its large radial curvature are
superseded by this candidate.

The existing interactive example uses `alpha=1`, `k=1000 pN/um`,
`F_*=30 pN`, and `epsilon=5 nm`. It does not visualize the alpha=3 revision.
These are illustrative, **not selected production values**.
Smooth symmetric geometry can have zero gradient at symmetric interior
points; this is not evidence of a strict constraint or an impenetrable wall.

### 4.3 Full-connector and pair energy

For each live two-arm connector and myosin:

```math
U_{\ell m}=\int_0^1 \Phi[p_\epsilon(x(t))]dt.
```

This is a normalized integral over the material connector, not an integral
multiplied by its length. Its point limit is `Phi`; a single connector's net
translation force is bounded by `alpha*F_*` for one connector–myosin
interaction; its axial component is bounded by `F_*`. These are not bounds
on the sum over distinct actin pairs or different myosins.

For unordered actin pair `P`, deduplicate redundant links by maximum energy:

```math
U_{Pm}=\max_{\ell\in P}U_{\ell m},\qquad
U_c=\sum_{P,m}U_{Pm}.
```

Every distinct pair contributes. Do not select one pair per tip. Do not
sum the redundant links or choose the winner by force magnitude. Handle
periodic images without duplicate counting.

For the selected connector, differentiate its energy:

```math
F_i=\int_0^1(1-t)f(t)dt,\qquad
F_j=\int_0^1t f(t)dt,\qquad F_m=-(F_i+F_j),
```

```math
\tau_i=(p_i-c_i)\times F_i,\qquad
\tau_j=(p_j-c_j)\times F_j,\qquad
\tau_m=-\int_0^1[x(t)-c_m]\times f(t)dt.
```

Use consistent unwrapped lever arms. Convert physical torques through the
existing orientation representation. Finite-difference orientation tests
must verify this conversion.

Maximum-energy selection is conservative almost everywhere but nonsmooth
at winner changes. Use deterministic ties, test switching explicitly, and
check work error under timestep refinement. Do not average arbitrary
winner forces: they must remain derivatives of the chosen pair energy.

### 4.4 Integration and implementation shape

```text
for each eligible two-arm connector:
    construct actual periodic material segment
    reject distant myosins conservatively
    for each candidate myosin and relevant image:
        clip segment to a conservative possible-contact interval
        integrate energy, endpoint forces, and myosin torque
        save candidate result under (unordered actin pair, myosin)

for each (pair, myosin), in a deterministic order:
    select greatest candidate energy with deterministic tie handling
    scatter selected derivatives to the two actins and myosin
```

Start with 8-point versus 16-point Gauss integration and refine intervals
where needed. Split or refine around force-onset and quadratic/linear
transitions. Initial candidate tolerances: absolute energy `1e-6 pN um`,
force `1e-3 pN`, plus relative `1e-5`; test tighter tolerances. Never silently
accept a refinement-limit failure. The earlier claim of exact quadratic
integration in a flat radial region does not directly apply to this new
smoothed distance energy.

Keep geometry and energy pure/testable. Cache only with explicit invalidation
on changes to poses, sites, topology, periodic images, ramp strength, alpha,
or other parameters. Recompute after kinetic commits before applying body
forces.

## 5. Mechanics, timestep, calibration, and inherited overlaps

First measure raw motor drive and the forces actually used by the integrator
on overshooting/contacting myosins. Record axial distributions and extremes,
noncontact force before/after caps, contact-pair counts, torques, and lateral
motion. First confirm analytic axial capture coverage for the target shell
and burial, including finite-connector force ratios. Select `alpha`, `k`,
`F_*`, approach distance, and smoothing from those coverage requirements and
controlled loading/timestep tests. Count projected pair contributions rather
than assuming every pair supplies `F_*` axially. Include alpha^2 radial
stiffness scaling and the larger radial force budget in stability estimates.

Preserve the intended noncontact cap policy initially, then add contact
gradients without subsequent independent bodywise clipping:

```text
body_force = cap_existing_noncontact_drive(noncontact_force)
             + conservative_contact_force
```

Audit torque/rotation caps as well; clipping contact-induced rotation can
also change the intended dynamics. Use stable integration rather than
clipping the new contact response. This does not make the rest of the
active, capped simulator conservative; the energy-gradient claim applies
to this contact term at fixed topology and fixed ramp strength.

For translations, rotations, and sliding coordinates together, inspect the
mobility-weighted Hessian

```math
\mathcal A=M^{1/2}(\nabla^2U)M^{1/2}.
```

In locally stable configurations, use `dt*lambda_max(A) <= 0.1` as an initial
accuracy target, not a universal stability certificate. Include existing
mechanical stiffness and changes of redundant-link winner. Verify relevant
observables at `dt`, `dt/2`, and, where needed, `dt/4`. Keep noise consistent
with mobility and timestep. Initially reduce the complete integration step;
do not assume force-only substeps with frozen chemistry are equivalent.

For resumes:

```text
load the checkpoint without geometric removal of inherited overlaps
ramp the contact energy from zero to full strength during startup
multiply energy, forces, and torques by the same ramp factor
use a stable timestep throughout, including at full strength
record external ramp work and changes in topology/overlap
separate the ramp transient from subsequent comparison statistics
```

Test a 5 ms ramp against a longer ramp. A ramp does not fix an unstable
final timestep. Do not claim it preserves the inherited state distribution.
Do not use projection to prepare or repair the new-model initial state.

## 6. Weak attachment motion and the corrected boundary rules

Use `b` as distance inward from the barbed end, with `b=0` at that physical
endpoint, `b=L` at the pointed endpoint, and `w=0.200 um` at the structural
inner barbed-zone edge. Implement conversion using the existing polarity
and material-coordinate helpers and test both actin orientations.

### 6.1 Sliding and physical end loss

Only weak sites slide in this change. Mature sites remain fixed, but their
actins translate and rotate under the same contact reactions.

```math
db=-\mu_b\frac{\partial U_{total}}{\partial b}dt
   +\sqrt{2D_bdt}\,\xi,\qquad D_b=\mu_b k_BT.
```

Use the total applicable attachment energy, including the contact term.
The chain rule must match the stored coordinate/direction convention.

- Remove reflecting walls at both ends of the **barbed zone**.
- Absorb a weak attachment when its path crosses an actual filament
  endpoint. In the current two-arm weak-link model, remove the weak link
  once when either attachment is lost and record endpoint loss.
- If the arm-state model is separately used, map the event to the relevant
  arm opening; do not silently change the production state model.
- Absorption is a coarse-grained end-loss assumption, not a unique molecular
  boundary law. Finite-rate end loss is a different future hypothesis.

Crossing detection must include excursions between step endpoints. For a
locally constant-drift one-dimensional Brownian substep with positive start
and end distances `b0,b1` from an absorbing endpoint, the bridge crossing
probability is `exp(-b0*b1/(D_b*dt))`. Treat this as a substep approximation
when forces vary; verify first-passage statistics and timestep convergence.
Handle either filament endpoint, correlated event ordering, and links that
already lie at an absorbing boundary explicitly. Do not reflect, clamp,
or reject a crossing as a replacement for end loss.

### 6.2 Soft structural edge at 200 nm

Weak sites may slide to `b>w`; no mechanical mirror exists there. Add a
continuous weak-link detachment hazard:

```math
\lambda_{edge}=k_{edge}\sum_{a=i,j}
\left(\frac{[b_a-w]_+}{\ell_{edge}}\right)^2,
\qquad [x]_+=\max(x,0).
```

Calibrate `k_edge` and `ell_edge` from survival versus distance and weak-site
distributions. These are not yet assigned production values. A larger
off-rate reduces binding occupancy outside the zone; do not additionally
introduce an uncalibrated mechanical soft-edge potential in the same change.

Birth still samples the allowed barbed zone, `0<=b<=w`. Maturation requires
both attachment sites inside that zone plus existing registry/support and
the exact 130 nm connector gate. Prevent a weak link that wandered outside
from maturing through a bypass path. Soft-edge loss is separate from
contact-driven molecular failure.

```text
for each weak link:
    propose site motion from total energy gradients and thermal noise
    if either path crosses a physical endpoint:
        remove link once; record weak endpoint loss
        continue
    accept sites, including positions beyond the 200 nm edge
    evaluate baseline + soft-edge + applicable contact hazards
    sample competing detachment channels
    otherwise permit maturation only if all placement conditions hold
```

## 7. Energetic association and kinetic-route audit

For a proposed weak link, compute the change of the **pair-budget contact
energy**, including all affected myosins:

```math
\Delta U_c=U_c^{after}-U_c^{before},\qquad
k_{on,new}=k_{on,existing}\exp[-\beta\max(\Delta U_c,0)].
```

This avoids charging identical redundant links repeatedly. Retain a
separately justified physical-core veto. Do not retain the old implicit
`axial_tip_contact => reject weak birth` as an additional soft-contact veto.
Likewise replace one-way weak diffusion penetration vetoes with the energy
gradient and the corrected endpoint/soft-edge rules.

At the illustrative energies this Boltzmann factor can be nearly zero:
report actual penalties, birth attempts, accepted births, distinct link
identities, and sensitivity to approach distance. This is a proposed
modification of association kinetics, not a claim that the whole existing
active kinetic model satisfies detailed balance.

Audit every route, including association, 1-to-2 maturation, memory/rebinding,
arm rebinding, weak diffusion, prebinding, status transitions, and resume.
The full-tip geometry helper is currently shared with association/diffusion
and projection-related queries; replacing mechanics must not accidentally
leave those implicit old gates active.

Inventory coarse actin/barbed-tip A-band proxies, including local-site and
whole-contour occlusion checks. Remove those that implement the retired
ACTN2 geometric policy, while preserving physical core sterics, actual
binding-coordinate constraints, and independent actin–myosin interactions.
Record the exact bare-tip endstop setting; the matched projection-off
configuration currently has the slab endstop disabled. Do not turn it on
implicitly as a replacement for ACTN2 burial deletion.

## 8. Contact-driven failure: separate calibrated stage

Stage A/B preserve the current baseline lifetime law and contour rupture.
Do not introduce a slip term based on the link's entire existing traction:
that would change lifetime behavior throughout the network.

For stage C, evaluate this candidate added channel:

```math
\lambda_{contact,\ell}=k_c
\left[\exp(T_{c,\ell}/F_b)-1\right],\qquad
\lambda_{total}=\lambda_{existing}+\lambda_{edge}+\lambda_{contact}.
```

`T_c` is contact-induced attachment traction from the same endpoint energy
derivatives. An initial explicit scalar choice is the larger of the two
endpoint contact-force magnitudes, after vector accumulation over relevant
myosins. Test directional alternatives if molecular evidence warrants them.
Never copy a whole-actin force to every attached link. Under the maximum
energy budget, the selected connector receives that contact derivative;
unselected redundant links retain their existing baseline mechanics.

The scalar failure law and `k_c,F_b` require controlled loading calibration;
they are not derived from the conservative potential. Compute loads before
the kinetic decision. Use stable exponential/hazard arithmetic and refine
steps when rates become too fast; do not hide this with an arbitrary force
cap. Do not add a parallel depth-triggered death hazard in this experiment.

For hazards evaluated over a substep:

```math
P_{detach}=1-\exp(-\lambda_{total}dt).
```

Attribute a stochastic death to its competing channel consistently with
the sampled hazards. Keep deterministic contour rupture distinct. Endpoint
first-passage versus chemical-event ordering must converge with timestep.

## 9. Step ordering, code targets, and event provenance

```text
refresh current geometry / contact energy at the current ramp strength
advance weak sites when that stage includes revised site motion
commit endpoint losses; recompute geometry/loads after accepted site changes
evaluate baseline, soft-edge and stage-C contact kinetics
check maturation placement; commit topology in canonical order
invalidate affected pair winners and cached contact derivatives
recompute contact forces for surviving/new connectors before body motion
combine noncontact drive with contact gradients using audited cap ordering
advance bodies with validated mobility, noise and timestep
record post-motion geometry and exact-step entry/exit histories
```

Preserve RNG/mutation order where practical, but new kinetics and path
crossing draws intentionally change trajectories. Statistical comparisons
cannot assume identical stochastic histories.

Proposed implementation locations (names may be adapted to existing layout):

- `cpp/include/actn2_contact_energy.h`: pure point/segment energy, derivatives,
  quadrature, conservative bounds, and result types.
- Existing `actn2_maturation_geometry.h`: retain the exact gate; reuse its
  buffered query for diagnostics, not a deletion decision.
- `src/sarcomere.cpp` / `src/interaction.cpp`: all-pair contact accumulation,
  endpoint forces and torques, body integration/cap ordering, ramp handling.
- `src/sarcomere_actn2.cpp`: weak diffusion/boundaries, kinetic-route audit,
  energetic association, stage-C failure and canonical commits.
- `src/sarcomere_io.cpp`, headers/main: versioned metadata, parameter
  persistence, resume validation and explicitly projection-off run manifests.
- `cpp/tests`: gradient/geometry, periodic, pair-budget, boundary first-passage,
  controlled mechanical loading and kinetic-channel regressions.
- Existing analysis/readers: new provenance, comparable geometry, cohorts,
  exact-step excursions and force/work decomposition.

Keep reason codes 0–5 historically interpretable. In particular, code 5
remains `mature_connector_exclusion` for old runs even after that mechanism
is retired. Proposed additions, after checking all reader/storage schemas:
6 `weak_endpoint_loss`, 7 `weak_outside_zone`, 8 `contact_load_failure`.
Do not reuse code 3 (legacy burial) for these events. Update array bounds,
metadata dictionaries, cohort analysis, and old-file fallbacks together.

Record status/cohort at the event, link identity, attachment coordinates,
force/load channel, and model version. Keep existing event rows readable;
use a versioned supplemental table for diagnostics that do not fit. Batch
and compress records rather than writing one file per event.

## 10. Validation and staged experiments

### Kernel and controlled-system acceptance

1. Finite-difference point, endpoint, translation, and orientation gradients.
2. Zero net contact force and torque for the complete interacting system.
3. Correct anisotropic bounds per connector–myosin interaction: axial
   `F_*`, radial/net `alpha*F_*`. Verify axial/radial force ratios across
   the intended capture region, including finite outward-running connectors.
4. Closed-loop work convergence at fixed topology/ramp strength, including
   redundant-link winner changes. Bookkeep energy changes at molecular
   events and during ramping separately.
5. Periodic invariance, relevant non-nearest images, endpoint/interior
   contacts, degenerate short connectors, cylinder corners, axis/midplane,
   force onset and quadratic/linear transition.
6. Identical redundant-link invariance and distinct-pair additivity.
7. Quadrature refinement and timestep convergence including rotations.
8. Weak-site free diffusion and energy drift; physical endpoint first-passage;
   no artificial reflection at 200 nm; survival outside the zone; inability
   to mature outside allowed placement; both filament polarities.
9. Controlled myosin-against-connector and myosin-against-network loading:
   force/displacement, lateral escape, rupture, and work balance.
10. Zero-contact kinetic behavior preserved; competing failure reasons
    attributed correctly; no stale forces after link removal.

### Stage definitions

| Version | Change relative to the previous stage |
|---|---|
| Reference | Pinned `r27zt0np` model and trajectory; old projection-off baseline |
| A | New conservative all-pair contact; remove instantaneous connector/burial deletion; preserve old weak birth/sliding rules for attribution |
| B | A plus energetic weak association, gradient-consistent weak sliding, absorbing physical endpoints, soft 200 nm edge, and complete maturation checks |
| C | B plus calibrated contact-induced failure |

Stage A deliberately retains the old weak-site restrictions as an interim
comparison; it is not the target final model. State that limitation in its
metadata. Versions are pinned source/binary milestones, not runtime switches.

Use the same 2.0 s checkpoint and the established 1600-actin/320-myosin
configuration. First run representative short continuations and startup
tests; then a matched 0.4 s continuation per validated stage. Compare common
physical times, not step counts when timesteps differ. If ramps differ,
compare equivalent post-ramp windows and show the omitted transient.
Use additional seeds for claims about lifetime, network integrity or
overshoot, rather than treating a single paired stochastic path as proof.

Record for every stage:

- Model version, git/source diff, binary hash, arguments, Slurm resources,
  timestep, seed, checkpoint, ramp history, alpha, target capture region,
  force budgets, and all other new parameters.
- Full-connector geometric overlap and buffered occupancy, independent of
  force winners; unique links/tips as well as link-tip pairs.
- Current mature, inherited mature, newly mature, and all-link cohorts.
- Exact-step entries/exits, depth, survival, recovery, demotion, and actual
  cause of death, including competing events.
- Raw and used motor/noncontact forces, contact endpoint forces/torques,
  pair counts, winner changes, axial/lateral displacement and work.
- Link lifetimes, connectivity, contour rupture, baseline death, endpoint
  loss, soft-edge loss, and contact failure; weak-site distributions.
- Association attempts versus distinct vetoed links, energy penalties,
  site-diffusion excursions, numerical refinement counts and wall time.

Saved 5 ms frames cannot establish sub-5 ms excursion lifetimes. First-entry
counts in old trajectories are counterfactual geometric exposures, not
additive predictions or lower bounds on *additional* deaths under a changed
network. Use exact-step histories for those claims.

## 11. Computational cost and profiling plan

Typical current configuration: about 6200 live connectors and 320 myosins,
or about two million cheap broad-phase checks per step. With `C` surviving
connector/image contacts and `Q` quadrature evaluations, aim for

```text
O(L*M + Q*C + C*log(C))
```

using a deterministic sort/group implementation of pair budgets. A hash
index may reduce grouping cost, but deterministic iteration must be retained.
Do not pay quadrature cost on all `L*M` pairs. Profile association and site
energy queries separately; they add candidate evaluations beyond mechanics.

Existing four-CPU, 200-step benchmark (job 50013712):

| Quantity | Measured |
|---|---:|
| Optimized wall time including startup/checkpoint I/O | 41.890 s |
| Effective wall time per step | 209.45 ms |
| Current tip-contact phase | about 27.3 ms/step |
| Remaining serial classification phase | about 3.2 ms/step |

This baseline **already includes** parallelized deletion geometry. Do not
credit the old 100-to-3 ms serial-phase saving again. Some geometry work
moved into the parallel mechanics phase; retiring deletion may save some
of that work, but the amount must be measured.

A provisional mechanics-stage scenario is

```math
t_{step}\simeq182+27f\ \mathrm{ms},
```

where `f` is replacement contact-phase cost relative to the current phase.
For 0.4 s simulated time and original `dt=10 us`:

| Contact-phase multiplier | Original dt | dt/2 | dt/4 |
|---|---:|---:|---:|
| 2x | 2.6 h | 5.2 h | 10.5 h |
| 4x | 3.2 h | 6.4 h | 12.9 h |
| 8x | 4.4 h | 8.8 h | 17.7 h |

These are scenarios on four CPUs, not measured new-kernel throughput or a
prediction for a 16-CPU allocation. Chemistry/association changes in B/C,
adaptive refinement and network evolution may change them. Multiplying a
radial distance by alpha adds little arithmetic, but the alpha^2 radial
stiffness and sharper force-direction transitions may require a smaller
timestep or more refinement. That may dominate integration cost.

Benchmark the pure kernel first on representative saved geometries, then
the full short continuation at the validated timestep and allocation.
Report candidate counts, quadrature work, pair grouping, mechanics, site
updates, kinetics and I/O. Preserve immutable job scripts and refuse to
overwrite existing benchmark result files. Read current Sherlock guidance
before submitting jobs; compile, numerics and plotting belong in Slurm.

## 12. Existing visualization and next concrete milestone

The existing standalone visualization of the isotropic predecessor is:

```text
$GROUP_HOME/sherry/sarcomere/sarcomere/analysis/visualizations/gaphunt-rig/
  status2_gate_20261028/results/interactive_energy/connector_energy_3d.html
```

It uses alpha=1 and **does not show the weighted alpha=3 candidate**. Preserve
it as the isotropic reference. Before candidate validation, regenerate the
visualization with alpha=1 versus alpha=3 comparisons and correctly labeled
axial/radial/net force bounds. It includes point and 40 nm finite-connector modes, equal-energy surfaces,
eight independent example connectors, custom positions, proportional force
arrows, endpoint contributions and a myosin reaction toggle. Its numerical
gradient checks and independent Python/JavaScript readout comparisons passed.
Actual headless Chromium rendering was not verified on Sherlock because
the browser requires a newer glibc. Download and open the self-contained HTML
in a modern local browser. These visualization checks are not a simulator
integration or production validation.

**Next implementation milestone:** update the candidate contours and force
ratios for alpha=3, then build the pure weighted-energy kernel, its
gradient/torque and periodic tests, and its representative cost benchmark.
Measure motor drive and coupled stiffness before selecting production
parameters or launching the staged 0.4 s comparisons.

Open calibration choices are `alpha`, `k`, `F_*`, approach distance, smoothing,
stable timestep, ramp duration, soft-edge length/rate, and the contact-failure
law/parameters. The alpha=1 geometric deficiency at the approximately 80 nm
actin shell is established. The weighted candidate's axial coverage must
first satisfy an explicit geometric target; sustained stopping, network
integrity, and failure under motor loading then require simulation tests.
