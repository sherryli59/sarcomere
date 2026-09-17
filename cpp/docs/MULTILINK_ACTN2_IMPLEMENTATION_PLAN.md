# Sparse local ACTN2 crosslinks: implementation plan

## Objective

Replace the current single binary AA edge with a small number of sparse,
local, KMC ACTN2 links. The model must support transient anywhere
prebundling, load-dependent catch stabilization, alignment from multiple
approximately perpendicular links, annealing by turnover, and collective
translation of an actin bundle relative to myosin.

This is an implicit-molecule model: an active ACTN2 link has explicit
attachment coordinates and a KMC state, but ACTN2 does not become a freely
diffusing rigid particle.

## Non-negotiable state semantics

Three concepts must remain separate.

1. **Occupancy:** whether an ACTN2 link exists and where it is attached.
2. **Catch state:** status 1 is an ordinary short-lived link; status 2 is the
   same existing link under favorable bilateral load.
3. **Cluster maturity:** whether several links have organized into a
   Z-disc-like cluster. This is not status 2 and will not be included in the
   first implementation phase.

A new link always enters status 1. Strong myosin occupancy never creates an AA
link and never directly assigns status 2.

## Phase-1 scope

Implement only local prebundling and catch stabilization:

- ACTN2 can associate anywhere within the axial overlap of two nearby actins.
- At most `aa_prebundle_links_per_pair = 2` links may exist between one pair.
- The number of links is also limited by a finite global ACTN2 pool when that
  option is enabled.
- Each link stores its own material attachment coordinates.
- Each link prefers an effective centerline separation and an orientation
  approximately perpendicular to both actin axes.
- Status-1 turnover permits annealing.
- Favorable transmitted load can promote an existing link to status 2 and
  extend its lifetime.
- No individual local link applies the mature 100--130-nm barbed-end registry
  potential.
- The current pair-level overlap potential remains available only as a
  temporary comparison control, then is removed from the default path.

## Local link geometry

For link `a` between actins `i` and `j`, store material contour coordinates
`s_i^a` and `s_j^a`, measured from each actin center:

```text
x_i^a = r_i + s_i^a u_i
x_j^a = r_j + s_j^a u_j
l_a   = minimum_image(x_j^a - x_i^a)
```

The coordinates remain fixed to their respective actin while the link is
bound. The first implementation uses an effective centerline connector with
rest length `d_AA*`; explicit actin-surface azimuths are deferred.

Use

```text
U_a = 1/2 k_l (|l_a| - d_AA*)^2
    + 1/2 k_perp [ (lhat_a.u_i)^2 + (lhat_a.u_j)^2 ].
```

The first term enforces preferred spacing. The second favors a connector
perpendicular to both filaments. Forces and orientation drives are obtained
from the full gradient with respect to both centers and both unit vectors.

Do not add the old global nematic AA term initially. Two separated local links
must be allowed to generate alignment mechanically. A weak direct angular term
may be restored only if a controlled two-link test shows that molecular
perpendicularity alone is insufficient.

The interaction must satisfy

```text
F_i^a + F_j^a = 0
```

to machine precision. Consequently AA links constrain internal modes but do
not pin the bundle center of mass.

## Selecting attachment coordinates

For each neighbor pair:

1. Construct the common axial overlap interval using the existing segment
   geometry and polarity-aware coordinate mapping.
2. Require the closest achievable connector length to lie within the soft
   capture reach.
3. Draw a candidate position within the common overlap rather than always
   selecting the segment midpoint or closest-point degeneracy.
4. Map that position to `s_i` and `s_j` so the initial connector is as nearly
   perpendicular as the geometry permits.
5. Reject the candidate if either attachment lies within
   `aa_site_exclusion` of an occupied attachment coordinate on that actin.

Initial parameter:

```text
aa_site_exclusion = 0.0192 um
```

This is a steric/site-exclusion scale, not permission to occupy every 19.2-nm
site. The explicit pair cap remains two during prebundling.

## Association capacity and rates

Do not make the maximum occupancy proportional to a one-micron overlap. Use

```text
N_cap,ij = min(N_prebundle,
               floor(max(0, ell_overlap - 2 delta_end)/aa_site_exclusion) + 1)
N_prebundle = 2.
```

Overlap should increase encounter probability only until it saturates:

```text
g_overlap(ell) = 1 - exp(-ell/ell_capture)
R_on,ij = (N_cap,ij - n_ij) k_on g_overlap(ell) c_free.
```

Initial value:

```text
ell_capture = 0.1 um.
```

`c_free=1` gives the current reservoir limit. The finite-pool option uses

```text
c_free = N_free / N_ACTN2,total.
```

The event probability is

```text
P_on = 1 - exp(-R_on,ij dt).
```

Retain the soft encounter barrier with `theta=0.1`, using only capture-scale
radial and perpendicular costs. The mature bound-link energy must not prevent
an initially imperfect encounter.

## Catch-state mechanics and load sharing

For one local link, let `F_i^a` be the force the link applies to actin `i`.
With actin axes pointing barbed to pointed, define raw favorable arm loads

```text
f_i,raw^a = max(0, -F_i^a.u_i)
f_j,raw^a = max(0,  F_i^a.u_j)
f_bilat,raw^a = min(f_i,raw^a, f_j,raw^a).
```

Passive prestress must not establish a catch state. Let

```text
f_opp = 1/2 |F_motor,i - F_motor,j|.
```

Distribute the available active drive among the links in proportion to their
transmitted force:

```text
w_a = |F_a| / sum_b |F_b|
f_active,budget^a = w_a f_opp
f_eff^a = min(f_bilat,raw^a, f_active,budget^a).
```

If the denominator is zero, all weights and effective loads are zero.

Use the existing hysteresis first:

```text
status 1 -> 2 when f_eff > f_on
status 2 -> 1 when f_eff < f_off
f_off < f_on.
```

The optional status-2 memory may preserve an already established catch state
briefly, but cannot establish one. A later control can replace deterministic
hysteresis with explicit force-dependent 1-to-2 and 2-to-1 KMC rates.

Lifetimes remain

```text
tau_1 = tau_0
tau_2 = tau_0 + tau_directional + lambda f_eff
P_off = min(1, dt/tau_status).
```

For `n_ij` independent links, evaluate one detachment draw per link. Do not
give every link the full pair force.

## Why bundles can move relative to myosin

AA forces are internal. For one connected bundle,

```text
sum_i F_AA,i = 0.
```

The common translation is therefore governed by net AM force and total drag:

```text
zeta_bundle V_bundle = sum_i F_AM,i + thermal noise.
```

Local AA links suppress actin-actin relative motion, while transient AM heads
allow the entire bundle to slide relative to the thick filament. No
isometric constraint or prescribed AM velocity is used in production.

Relative actin-actin sliding requires ACTN2 compliance and turnover. A bound
local link keeps its stored material coordinates; after it detaches, a new
link binds at the current geometry. Repeated turnover therefore relaxes the
bundle registry without requiring the whole pair to become completely
unlinked at once.

## Sparse state representation

Do not allocate link slots for all `N_actin^2` pairs.

```cpp
struct AACrosslink {
    std::uint64_t id;
    int actin_i;
    int actin_j;
    double s_i;
    double s_j;
    std::uint8_t status;
    std::uint32_t attach_step;
    std::uint32_t status2_memory_until;
};
```

Maintain:

- `aa_links`: a slot pool of active records;
- `aa_free_slots`: recycled slots;
- `aa_active_ids`: compact active-link iteration order;
- `aa_links_by_actin[i]`: incident link IDs;
- `aa_links_by_pair[pair_key(i,j)]`: zero, one, or two active IDs;
- `aa_degree[i]`: number of distinct AA neighbor filaments;
- `aa_link_count[i]`: total incident ACTN2 links.

Use a canonical 64-bit pair key from `min(i,j), max(i,j)`. With two links per
pair and approximately two AA neighbors per actin, a 2000-actin simulation has
only a few thousand records. This is smaller than the existing dense AA state
matrices.

## Deterministic parallel update

Split each step into explicit stages.

1. Freeze the previous active-link state.
2. Evaluate local-link forces in parallel into thread-local actin force and
   torque arrays.
3. Evaluate status and detachment proposals independently for active links.
4. Generate association proposals from the actin neighbor list in parallel.
5. Sort proposals by canonical pair key, attachment coordinate and deterministic
   random key.
6. Commit proposals in one ordered pass while reserving pair capacity, actin
   sites and the global ACTN2 pool.
7. Rebuild or incrementally update sparse adjacency indices.

No force loop may mutate KMC occupancy. No OpenMP worker may reserve an ACTN2
site directly.

## HDF5 and checkpointing

Bump the state schema version. Avoid a padded
`frames x maximum_links x fields` trajectory.

Use ragged CSR-style snapshots:

```text
/aa_links/records
    [id, i, j, s_i, s_j, status, attach_step, memory_until]
/aa_links/frame_offsets
    cumulative record offset for each committed frame
/aa_links/events
    [step, id, i, j, old_state, new_state, reason, s_i, s_j, f_eff]
```

The resume checkpoint stores only the current active pool, free-list state,
next stable ID and RNG state. Preserve the two-generation checkpoint
transaction already used by the AM implementation.

Continue writing per-actin aggregate diagnostics:

- total AA link count;
- status-1 count;
- status-2 count;
- distinct AA neighbor count;
- maximum and mean effective link load;
- AA force and torque;
- association and detachment counters.

During the temporary validation period, derive the old pair-level bond and
status summaries from the sparse records at save time. Do not keep dense
matrices as the source of truth.

## Implementation sequence

### Step 1: Pure geometry and mechanics

- Add `AACrosslink` and local-link energy/gradient functions.
- Add overlap-coordinate mapping and site-exclusion helpers.
- Implement forces for manually constructed links with KMC disabled.
- Verify translational invariance, pair-order symmetry and PBC behavior.

### Step 2: Sparse state container

- Add stable IDs, slot recycling, per-actin adjacency and per-pair lookup.
- Add deterministic insertion/removal and consistency checks.
- Do not connect this container to the production update yet.

### Step 3: KMC occupancy

- Add capped association and independent detachment.
- Add finite-pool accounting behind a runtime option.
- Implement link events and restart serialization.
- Compare occupancy statistics with analytic birth-death predictions.

### Step 4: Per-link catch classification

- Compute per-link raw arm loads.
- Distribute the opposed active-force budget across simultaneously occupied
  links.
- Apply status hysteresis and status-2-only lifetime stabilization.
- Validate rates under force clamps.

### Step 5: Replace the production AA path

- Drive mechanics and KMC from sparse links.
- Derive compatibility diagnostics from sparse state.
- Remove the pair-level barbed-overlap force from ordinary links.
- Keep a temporary runtime comparison flag only until validation completes.

### Step 6: Full diagnostics

- Run local mechanics, KMC, annealing and bundle-sliding tests below.
- Run a reduced assembly pilot before 2000-actin production.
- Delete the temporary legacy comparison path once results are archived.

### Step 7: Mature Z-cluster model, only after Phase 1 passes

Add a separate pair-level maturity variable `m_ij`. Candidate maturation
criteria to compare are:

- at least two sustained local links;
- antiparallel polarity;
- sustained favorable load;
- barbed-end enrichment;
- a minimum maturation dwell time.

Only a mature cluster may increase its capacity toward four to six layers or
acquire a collective 100--130-nm registry potential. Run a no-maturation
control and do not assume which maturation criterion is necessary before the
prebundling model is tested.

## Unit tests

1. A perpendicular link at its rest length has zero force and torque.
2. Axial or angular perturbations generate restoring forces of the correct
   sign.
3. Swapping actin indices exchanges forces and torques without changing
   energy.
4. Translating both linked actins together leaves energy unchanged.
5. Internal AA forces sum to zero at every tested geometry.
6. Two separated links align initially tilted parallel and antiparallel rods.
7. One link permits the expected residual rotational degree of freedom.
8. Site exclusion rejects attachments that are too close on either actin.
9. Pair occupancy never exceeds two in prebundle mode.
10. Global occupancy never exceeds a finite ACTN2 pool.
11. Sparse insertion/removal preserves all adjacency invariants.
12. Checkpoint/resume reproduces an uninterrupted deterministic trajectory.

## Minimal stochastic tests

### Birth-death statistics

For fixed geometry and load-free identical sites, verify the stationary
occupancy distribution against the finite-site birth-death model. Check
association counts, detachment counts and dwell-time survival curves directly,
not only mean occupancy.

### Catch activation

Run no-load, favorable-load and unfavorable-load clamps. Require:

- all new links begin in status 1;
- passive AA prestress alone cannot establish status 2;
- favorable bilateral motor load produces status 2;
- removing load returns links to status 1 after hysteresis/memory;
- measured status-specific off rates match the implemented lifetimes.

### Annealing

Start two actins with imperfect axial and angular registry. Compare:

- permanent links;
- status-1 turnover;
- status-2-only stabilization;
- status-1 stabilization as a trapping control.

Require turnover to improve alignment/registry without requiring simultaneous
loss of all pair connectivity.

### Bundle translation relative to myosin

Crosslink two or more actins, attach transient AM heads and apply:

- balanced opposing motor forces;
- common-direction motor force;
- an external force clamp;
- no AA control.

Measure bundle center-of-mass displacement, internal relative displacement,
AM head cycling and AA strain. Under common-direction drive, AA must reduce
internal displacement while allowing bundle translation relative to myosin.
Balanced opposing drive should instead generate tension with little center-of-
mass motion.

### Alignment from link geometry

Use one versus two local links at controlled axial separation. Measure angular
relaxation and verify that alignment strength grows with link number and lever
arm separation. Compare against a direct nematic-potential control.

## Reduced assembly tests

Before full scale, run matched replicates with:

1. AA disabled;
2. one local link per pair;
3. two local links per pair;
4. two links with catch stabilization disabled;
5. two links with status-2-only stabilization;
6. two links with status-1 stabilization as the kinetic-trapping control.

Track actin and myosin nematic order, legacy `phi_3d`, bundle-size
distribution, status populations, link occupancy, turnover, relative sliding,
force/torque caps and AM cycling.

## Acceptance gates for full-scale production

- All geometry, sparse-state and restart tests pass.
- No NaNs or invalid link IDs.
- Force and torque cap fractions below `1e-4` outside documented fixed-particle
  diagnostics.
- Measured KMC statistics agree with rates within sampling uncertainty.
- A perfect local bundle remains stable.
- An imperfect bundle anneals with status-1 turnover.
- A linked actin bundle translates under common AM drive instead of becoming
  artificially isometric.
- Parallel status-1 prebundles remain transient.
- Favorably loaded status-2 links persist longer than status-1 links.
- Sparse AA memory and output size scale with active links, not `N_actin^2`.
- Reduced assembly shows improvement over AA-disabled controls without
  pathological arrest.

Only after these gates pass should the 2000-actin/1000-myosin production array
be submitted with its babysitter and completion analysis dependency.

## Implementation status (2026-07-21)

Implemented as the production default through `--aa_local_multilink=true`.
The legacy one-edge-per-pair model remains available explicitly with
`--aa_local_multilink=false`:

- sparse persistent link instances with stable IDs, pair/actin indexes,
  free-slot reuse, local contour coordinates, status, attachment step and
  status-2 memory;
- up to two transient links per actin pair, distinct-neighbor valence, an
  optional global ACTN2 pool, overlap-saturating association, pair-local site
  exclusion and optional transition-state energy gating;
- local connector mechanics with a preferred centerline distance and a
  dimensional penalty on connector displacement along both actin axes;
- per-link tensile status-2 classification, hysteresis/memory, shared opposed
  motor load and status-dependent KMC detachment;
- sparse ragged HDF5 trajectory records, exact checkpoint/resume state, link
  event logs, metadata validation and per-actin status/link/tension fields;
- local-link prebinding for the six-actin ring and derived dense pair summaries
  for existing analysis.

Validation completed:

- all 13 compiled tests pass, including local rest geometry, restoring force,
  pair-order symmetry, oblique-link torque and periodic-boundary invariance;
- a fixed-geometry KMC control produced 1,984 completed lifetimes with mean
  1.995 ms versus the prescribed 2.000 ms (SEM 0.043 ms), while pair and pool
  caps held;
- a perfect six-actin lattice remained at the 0.0506 um AA spacing with zero
  force/torque caps while the second link per edge populated;
- a perturbed two-link lattice reduced nearest-neighbor spacing RMSE from
  0.0145 to 0.00679 um and actin-axis tilt RMS from 10.09 to 6.12 degrees
  without force/torque caps.

The single-link control explains why two separated links are important: one
local linker tethers a pair but only weakly determines their relative axes.
Link turnover plus two separated attachment sites supplies the intended
alignment and annealing mechanism.

## Coupled motor-load controls (2026-07-21)

Matched five-replicate minimal-mechanics controls (ten replicates for the
high-catch condition) passed the two outstanding coupling checks:

- Under common-direction motor drive, a two-link actin bundle translated by
  0.1661 +/- 0.0040 um while accumulating only 0.00368 +/- 0.00024 um of
  relative path.  Without AA links, relative path was 0.1153 +/- 0.0039 um.
  Thus local links reduced internal slip about 31-fold without pinning bundle
  translation.
- Under opposing fixed-myosin drive, bundle center-of-mass displacement was
  statistically near zero while relative motion and AA tension developed.
  Status 2 occupied 6.2% of link time at catch coefficient 0.1 and 9.3% at
  coefficient 1.0.  Free-myosin and no-powerstroke controls produced no
  status-2 transitions, demonstrating that proximity alone does not activate
  the catch state.
- At catch coefficient 1.0, the pooled status-1 detachment rate was 104.8/s
  versus the prescribed 100/s.  The status-2 detachment rate was 34.48/s,
  consistent with 35.27/s inferred independently from the event-recorded
  force-dependent lifetimes.
- All 35 trajectories were finite and both force- and torque-cap fractions
  were zero.

The reduced-system motor-load and bundle-translation acceptance gates are
therefore satisfied.  Full-scale assembly still requires its own cap,
occupancy, turnover, runtime and ordering audit.
