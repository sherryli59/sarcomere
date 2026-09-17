# Dynamic titin-like capture for large sarcomere assembly

Date: 2026-07-31

## Objective

Extend the validated six-filament centering mechanism to the established
large assembly (`1280` actins, `640` myosins, `141` heads per myosin half)
without preassigning actin-to-myosin ownership.

The model must answer whether local A-M engagement, turnover-capable ACTN2,
and a slowly remodeling titin/Z-disc linkage can organize finite A-M overlap
and approximately `0.10 um` A-A register in a fully Brownian active assembly.

The existing `balanced_owning` mode remains a small-system positive control.
It is not a large-assembly model and must not be labeled production.

The current implementations of both `legacy_crosslinked` and
`balanced_owning` are bidirectional and can generate force under compression.
They must therefore be removed from the physical titin namespace. Preserve
their enum values only for checkpoint compatibility, report them as legacy
effective-centering modes, and add a separate tension-only `balanced_titin`
regression mode. Previously generated `balanced_owning` trajectories are
effective-centering controls, not physical titin results.

## Physical interpretation

Myosin heads dynamically recognize a geometrically compatible actin. Sustained
head occupancy can nucleate a titin-like structural connection, but the
connection is not another myosin head:

- Head binding/unbinding occurs on the fast motor-cycle timescale.
- Titin-like capture requires persistent A-M engagement and local Z-disc
  evidence.
- Once captured, the link survives ordinary head turnover.
- Titin-like release and rebinding occur on a slower, force-dependent
  timescale.
- A released link may recapture only a nearby compatible barbed-end site; it
  cannot teleport to a distant actin.

Every mode named titin will be tension-only. The compression-bearing Hookean
spring that produced the six-filament positive control will remain available
only under an explicitly named `effective_centering` diagnostic mode. If
tension-only bilateral links do not generate a finite-overlap stall, the
result will be reported as a missing physical interaction rather than quietly
treating compressed titin as literal titin.

## Coarse-grained topology

Each myosin has two independently capturable halves. Each half can own one
coarse-grained titin-like link in the first implementation.

A link connects:

- the myosin center/M-line proxy; and
- the barbed endpoint of one actin that belongs to an ACTN2-supported local
  Z-disc environment.

The two half-links must lie on opposite axial sides of the myosin. Ownership
therefore emerges as assembly proceeds rather than from particle indices.

Initial constraints:

- At most one active titin-like link per myosin half.
- At most one titin-like owner per actin.
- The candidate actin must match the myosin half used by its engaged heads.
- The barbed end must lie on the corresponding side of the myosin center.
- The actin must have a live barbed-zone ACTN2 link or a retained ACTN2
  neighbor slot.

This produces at most `1280` links in the `1280/640` production system and
keeps the mechanics linear in system size.

## State and data structures

Add a sparse state layer to `Sarcomere`:

```cpp
enum class TitinCaptureState : uint8_t {
    unbound = 0,
    candidate = 1,
    bound = 2,
    refractory = 3,
};

struct TitinLink {
    uint64_t id;
    int myosin;
    int half;
    int actin;
    TitinCaptureState state;
    size_t candidate_since;
    size_t attach_step;
    size_t refractory_until;
    double support_ema;
    double rest_length;
    bool active;
};
```

Maintain:

- two slots per myosin (`2 * n_myosins`);
- an `actin -> titin slot` ownership vector;
- thread-local candidate accumulators;
- event buffers and per-frame mechanics buffers;
- checkpointed candidate dwell and refractory state.

Do not use a dense `n_actin * n_myosin` titin matrix.

## Dynamic recognition algorithm

Evaluate capture at a configurable interval, initially every `100` steps
(`0.1 ms` at `dt=1 us`), rather than adding a second all-pairs search every
integration step.

For every myosin half:

1. Traverse its existing explicit head array and aggregate weak and strong
   occupied heads by bound actin.
2. Reject an actin unless all geometry gates pass:
   - productive A-M overlap exceeds `titin_capture_overlap_min`;
   - the actin is assigned to this myosin half by existing head geometry;
   - its barbed end is on the correct side of the myosin center;
   - its barbed end lies inside a finite capture shell around the target
     center-to-Z distance;
   - it has live or retained barbed-zone ACTN2 support;
   - neither the half nor actin is already owned by another active link.
3. Rank valid candidates deterministically by sustained support EMA,
   strong-head count, total occupied-head count, productive overlap, smaller
   elastic capture cost, and actin index as the final tie-breaker.
4. Continue candidate dwell only while the same actin remains the winner.
   Switching candidates resets dwell rather than transferring accumulated
   recognition time.
5. Form a link after the dwell threshold and a counter-based stochastic
   association test. Include an elastic-energy gate so capture cannot inject
   a large discontinuous spring energy.

Initial mesoscopic test values, not final calibration:

- capture interval: `0.1 ms`;
- minimum occupied heads: `3`;
- minimum strong heads: `1`;
- support-EMA timescale: `2 ms`;
- continuous dwell: `5 ms`;
- minimum productive overlap: `0.10 um`;
- capture shell: `+/-0.15 um` around the nominal rest length;
- one link per myosin half.

## Link mechanics

### Physical mode: `dynamic_capture_tension_only`

Use a flat-bottom or WLC-like tensile law:

```text
extension = distance - rest_length
force = 0 exactly                              if extension <= 0
force = min(k_titin * extension, force_cap)    if extension > 0
```

Apply equal and opposite forces to conserve linear momentum. Because the
Z-side force acts at an actin barbed endpoint, apply the corresponding
`r x F` torque to the rigid actin. The myosin anchor is initially its center,
so it receives no off-center titin torque.

Two links on opposite myosin halves provide centering through tension
imbalance: displacement toward one Z environment slackens that side and
stretches the other. A one-sided link is permitted during assembly but must
be reported separately because it is not a centered sarcomere.

### Diagnostic mode: `dynamic_capture_effective_centering`

Use the validated bidirectional Hookean law only as an effective centering
network. Store this mode distinctly in metadata and figures. Never label its
compression-bearing force as literal titin.

### Release and local recapture

Use a slow slip-bond release law initially:

```text
k_off(F) = k_off_0 * exp(max(F, 0) / F_detach)
```

After release:

- clear ownership immediately;
- preserve the former local partner as a preferred candidate during a short
  refractory/local-memory interval;
- require normal geometry and A-M recognition gates before recapture;
- allow a different actin only after local memory expires.

Loss of all ACTN2 support starts a grace timer; it does not instantly delete
the titin link during molecular ACTN2 turnover.

## CLI and metadata

Add explicit options while preserving current defaults:

```text
--titin_mode=off|balanced_tension_only|dynamic_capture_tension_only
--effective_centering_mode=off|legacy_crosslinked|balanced_owning|
                            dynamic_capture
--titin_k
--titin_rest_length
--titin_force_cap
--titin_capture_interval
--titin_capture_heads
--titin_capture_strong_heads
--titin_capture_dwell
--titin_capture_support_tau
--titin_capture_overlap_min
--titin_capture_shell
--titin_require_aa_support
--titin_links_per_half
--titin_k_on
--titin_k_off
--titin_detach_force
--titin_refractory_time
--titin_aa_loss_grace
```

Write all values to `/meta`. Increment a titin schema version and reject
incompatible resumes.

Add HDF5 datasets:

- `/titin/link_identity`;
- `/titin/link_mechanics`;
- `/titin/frame_offsets`;
- `/titin/events`;
- checkpoint copies of active links, candidate dwell, ownership, refractory
  state, next link ID, and support EMA.

Event reasons must distinguish capture, stochastic release, force-assisted
release, lost-ACTN2 release, geometry failure, and local recapture.

## Code changes

1. Add pure capture/rate helpers in `cpp/include/titin_capture.h`.
2. Add link state, parameters, sparse indices, and event buffers in
   `cpp/include/sarcomere.h`.
3. Add candidate aggregation, deterministic commit, force/torque application,
   release kinetics, and checkpoint handling in `cpp/src/sarcomere.cpp`.
4. Add CLI parsing, validation, and mode names in `cpp/src/main.cpp`.
5. Keep `minimal_catch_balanced` as the exact positive-control topology in
   `cpp/src/initialization.cpp`.
6. Add analysis for titin occupancy, bilateral fraction, lifetimes, force
   distribution, partner switching, and capture distances.
7. Create a new large-production SLURM manifest derived from
   `cpp/production_productive_register_20260726.sbatch`; do not alter or
   resume the prior control trajectory.

## Verification tests

### Unit tests

- Deterministic candidate ranking and tie-breaking.
- Correct myosin-half and barbed-end polarity gates.
- Candidate dwell resets on partner switching.
- Capacity: one owner per actin and configured links per half.
- Tension-only spring is exactly slack under compression.
- No mode whose metadata name contains `titin` produces compressive force;
  the legacy bidirectional mode is recorded as `effective_centering`.
- Equal-and-opposite force and correct endpoint torque.
- Force-dependent release statistics.
- Counter-RNG invariance across OpenMP thread counts.
- Serialization/resume equality for candidate and bound states.

### Small-system regressions

- `4 actin / 2 myosin`: dynamic capture must recover the intended two-sided
  topology without index ownership.
- Head-turnover ablation: captured titin survives temporary zero AM occupancy.
- ACTN2-turnover challenge: brief molecular loss does not erase the Z
  association; prolonged cluster loss can release it.
- No-ACTN2 control: A-M occupancy alone cannot permanently create a Z tether.
- No-motor control: titin capture cannot occur without recognition.
- Brownian and deterministic variants, with rotational diffusion and endpoint
  torque enabled.

### Mesoscopic gates

Run `40 actin / 20 myosin`, then `128 actin / 64 myosin`, for at least four
seeds and compare:

1. no titin;
2. fixed-owner positive control where applicable;
3. dynamic tension-only capture;
4. dynamic effective-centering capture.

Required diagnostics:

- fraction of myosin halves with zero, one, or bilateral titin capture;
- A-A and A-M overlap distributions, not only their means;
- capture/release/recapture lifetimes;
- rate of partner switching and maximum switch distance;
- weak/strong head occupancy;
- ACTN2 live-link and neighbor-slot occupancy;
- net titin force and torque conservation;
- force/torque/rotation cap counts;
- cluster size and sarcomeric-order metrics.

## Production rollout

The production configuration will retain the previous assembly scale and
baseline:

```text
n_actins=1280
n_myosins=640
am_heads_per_half=141
box=7.72 x 2 x 2 um
pbc_mask=111
initial_structure=myosin_x_noise
dt=1 us
physical translational and rotational rod mobility at 310 K
```

Overlay the promoted ACTN2 settings:

```text
aa_status_mode=productive_register
aa_register_all_statuses=true
aa_antiparallel_overlap=0.10 um
aa_barbed_zone=0.10 um
aa_prebundle_links_per_pair=3
aa_stiffness_reference_links=3
```

Do not copy the six-filament `k_aa=1000 pN/um` blindly into the crowded
production assembly. Compare it with the prior production value
`k_aa=3000 pN/um` at the mesoscopic gate and select by per-link force and
turnover, not aggregate stiffness alone.

Rollout sequence:

1. `20 ms` large-system initialization/capture smoke test.
2. `100 ms` large-system pilot with frequent checkpoints.
3. Inspect capture topology, numerical caps, and overlap distributions.
4. Launch the `1 s`, single-seed `1280/640` production run only after those
   gates pass.
5. Add three more seeds only if the single-seed trajectory shows bilateral
   capture, ACTN2 turnover, and a finite-overlap plateau.

## Acceptance criteria

Implementation is complete only when:

- no index-based ownership is used in dynamic mode;
- capture depends on sustained, geometrically valid A-M engagement;
- links persist independently of fast head turnover;
- ACTN2 support is required but molecular ACTN2 turnover remains possible;
- no actin or myosin is multiply owned beyond configured capacity;
- forces are reciprocal and endpoint torques are accounted for;
- no long-distance partner teleportation occurs;
- checkpoint/resume reproduces the uninterrupted trajectory;
- the large-system pilot remains finite and does not rely on numerical caps;
- overlap stall is measured from a late-time velocity confidence interval;
- tension-only and compression-bearing results are labeled separately.

## Decision rule

If the tension-only mode produces bilateral capture but does not stall deeper
A-M sliding, do not increase titin stiffness until it behaves like a
compression spring. That outcome means literal titin is insufficient for the
requested finite-overlap selector. The next physical candidates would be an
explicit Z-disc/M-line network, lattice confinement, or another passive
overlap-dependent interaction.
