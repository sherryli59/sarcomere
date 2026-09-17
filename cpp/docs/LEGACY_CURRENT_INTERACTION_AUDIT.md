# Legacy versus current interaction audit

Date: 2026-07-21

## Scope and provenance

This audit compares the code that produced the high-order `data0430` sweep with
the current explicit-crossbridge development tree. It compares executed force,
velocity, geometry, and KMC paths rather than option names alone.

- Legacy source: git commit `30829b0` (the revision identified in
  `analysis/compare_status1_old_new.py` for `data0430`).
- Legacy clean reference checkout: `/tmp/sarcomere_legacy_30829b0.uAuzGj`.
- Legacy reference executable SHA-256:
  `36e214918345e695583df5a155349e2f86c995ea47f1c931cd943f6779474632`.
- Current tracked base: git commit `ef6e4d0c63b6bf2e6f419813435f622ccb5c7f47`,
  plus the uncommitted explicit-AM, ACTN2, mobility, sterics, and diagnostics
  changes in the working tree.
- Current `cpp/build/sarcomere` SHA-256 at audit time:
  `2717514e6dd55721261ef0d88a59eb95a111448e845efea7fc4526ae078a9b80`.

The old executable was rebuilt successfully from the exact legacy source in a
separate clone. The temporary checkout is evidence, not a permanent archive;
the commit hash is the durable provenance.

## Executive conclusion

The old high-`phi_3d` behavior cannot be attributed to the catch-lifetime law
alone. It emerged from a composite coarse-grained assembly mechanism:

1. every geometrically eligible AM pair was continuously attached;
2. every such AM pair received a strong, global nematic alignment energy;
3. actins were given a prescribed polarity-dependent sliding velocity;
4. AA links imposed a very stiff and asymmetric axial registry in addition to
   radial and angular alignment;
5. a weak system-wide transverse bundling field was active; and
6. the pre-phase removed myosin-myosin overlaps but left dense AA and AM
   contacts largely intact.

The current model simultaneously replaces (1)-(3), changes (4), usually
disables (5), and reverses (6). Those changes are physically motivated, but
they make the new full-scale runs a different model rather than a clean
mechanistic refinement of the old one.

The highest-priority missing organizer is the legacy AM nematic torque. Under
`dynamic_am=true`, `kappa_am` does not enter the AM force path. The explicit
model must generate alignment only through spatial head forces and their
off-center lever arms. Therefore setting `--kappa_am=50` in an explicit-AM
submission does not restore the old aligning interaction. Similarly,
`--k_am=3000` does not control explicit-AM proximity; `am_tether_k` does.

## Reconstructed legacy high-order operating point

The following is the relevant `data0430` regime, reconstructed from the exact
legacy defaults and the matched `cpp/catch_knockout.sbatch` manifest.

| Quantity | Legacy high-order value |
|---|---:|
| Box and PBC | `12 x 2 x 2`, `111` |
| Myosins | 1000 |
| Actins | swept; best documented group includes 2100 |
| Time step | `1e-5` s |
| Save cadence | 100 steps |
| Initial condition | isotropic actins, noisy x-directed myosins |
| Sterics-only pre-phase | 500 steps |
| Actin scalar diffusion | translational 0.2, rotational 0.1 |
| Myosin scalar diffusion | translational 0.05, rotational 0.05 |
| `beta` | 241 |
| Actin length | 1 um |
| Myosin length | 1.5 um |
| Myosin radius | 0.058 um |
| AM/AA capture cutoff | 0.08 um |
| AM/AA preferred distance | 0.06 um |
| `k_am`, `k_aa` | 10, 10 in normalized-distance potentials |
| `kappa_am`, `kappa_aa` | 100, 100 |
| Prescribed motor speed `v_am` | 10 um/s |
| Maximum actins per myosin | 6, split by polarity group |
| AA `k_on` | 1000 1/s |
| Base AA lifetime | 0 in the pure-catch sweep |
| Catch coefficient | lambda sweep |
| Transverse bundling | `k_bundle_max=0.1`, ramp 5000 steps |
| Myomesin and titin | off |

Representative group means already extracted from the old trajectories are:

| Actins | lambda | mean `phi_3d` | replicate SD | n |
|---:|---:|---:|---:|---:|
| 2100 | `1e-3` | 0.922 | 0.004 | 2 |
| 2100 | `1e-2` | 0.851 | 0.066 | 3 |
| 1800 | `10^-2.5` | 0.835 | 0.034 | 3 |
| 1800 | `1e-3` | 0.811 | 0.003 | 2 |
| 2100 | `10^-2.5` | 0.797 | 0.176 | 3 |

These are finite-time outcomes of the entire legacy interaction set; they are
not measurements of the catch law in isolation.

## Interaction-by-interaction comparison

| Mechanism | Exact legacy behavior | Current behavior | Assessment |
|---|---|---|---|
| AM occupancy | A geometrically eligible actin-myosin pair is effectively continuously bound. There is no detached/weak/strong crossbridge cycle. | A conserved set of heads on each myosin half undergoes detached/OFF, weak, and strong KMC transitions. Pair occupancy is intermittent. | Intentional physical refinement; large loss of duty ratio and organizing time. |
| AM radial proximity | For each eligible AM pair, an overlap-scaled spring acts on segment distance. | With dynamic AM, a flat-bottom pair tether is distributed over its occupied heads. Default `am_tether_k=1000`, slack 0.005 um. | Replacement exists. `k_am` is not the dynamic-AM stiffness. |
| AM angular alignment | Every geometric AM pair has `0.5*kappa_am*(1-(a.m)^2)`, even outside the radial spring cutoff. | Dynamic AM does not call the legacy AM potential. Alignment comes from off-center tether plus active head forces. | **Major lost organizer.** `kappa_am` is inert for dynamic AM. |
| AM active motion | Bound actins get prescribed velocity `v_am*a`; status-2 load reduces it by `(1-f_load)`. Status-2 reaction drift is assigned to myosin. | Strong heads use an axial motor spring, power stroke, stall force, and force-velocity evolution. Particle motion follows force and mobility. | Intentional replacement; old high order may partly be velocity-driven sorting. |
| Myosin polarity | Myosin is nematic, while actin polarity selects motion and per-half occupancy groups. | Myosin remains nematic. Explicit halves and a motorless central bare zone are represented; head ownership uses half and productive polarity. | Retained and made more explicit. |
| M-line/end-zone handling | A nonphysical axial endstop acts on the actin right tip and can push it outward after deep penetration. | The endstop was removed. Productive half-zone overlap, bare-zone geometry, head retention, and dissociation are intended to handle exit. | Do not restore the exact force. Preserve its sorting function with correct half-zone rules. |
| AA representation | At most one implicit pair bond; one matrix element stores its status. | Pair AA remains available. Sparse local multi-ACTN2 with multiple spatial links is implemented behind `aa_local_multilink`. | New model is richer, but local multilink is off by default and absent from most full-scale jobs. |
| AA capture | Distance between crosslinkable actin subsegments must be below the cutoff. | Capture may use pair geometry or local overlap sites, with encounter energy gating and site exclusion. | Refined; more restrictive in current physical pilots. |
| AA radial spring | `0.5*k_aa*((d-d*)/d*)^2`; with `k_aa=10`, `d*=0.06`, the small-extension stiffness is about 2778 in dimensional units. It ceases outside the cutoff. | Pair AA uses `0.5*k_aa*(d-d*)^2` and remains mechanically active until KMC detachment. Physical runs use about 3000. Local AA has its own dimensional connector spring. | Current `k_aa=3000` roughly restores radial curvature and fixes silent mechanical rupture. |
| AA angular alignment | Nematic `0.5*kappa_aa*(1-(a_i.a_j)^2)`. | Pair AA retains this. Local AA instead penalizes connector projections along each rod, driving an approximately perpendicular ACTN2 geometry. | Pair term retained; local term is more physical. |
| AA axial registry | Strong asymmetric barbed-end energy: a tenfold penalty for one sign of axial offset plus a target near 0.1 actin length. | Pair AA uses a symmetric, polarity-weighted barbed-end registry and configurable antiparallel overlap. Local AA sites can bind anywhere in overlap and slide/rebind. | Major semantic change. Legacy registry was a strong assembly bias and pair-order-sensitive. |
| AA bond persistence | Off-rate lifetime is `base + directional_base(status2) + lambda*f_load`. Both status 1 and status 2 receive `lambda*f_load`. | Default has short-lived status 1; catch stabilization is status-2-only unless `stabilize_status1=true`. | Intentional mechanistic intervention, not a legacy match. |
| Meaning of AA load | `f_load` is a bounded AM overlap/alignment proxy multiplied by AA alignment; it is not mechanical force. | Default `productive_opposed` uses strong heads with productive half/polarity and opposed active motor force. Other classifiers exist. | Intentional correction; status-2 occupancy will be much lower. |
| Meaning of status 2 | Antiparallel AA geometry, nonzero AM proxy on both actins, and attachment to different myosin identities. | Existing AA link plus productive strong-head forces on both actins with opposed pulling, hysteretic on/off thresholds, and optional memory. | More physical but far stricter than legacy. |
| Refractory time | `tau_rec=0.05`, but the legacy Boolean combines prior bonding and recovery before testing recovery, making the rejection condition unreachable. | Reassociation is actually blocked until recovery expires. | **Legacy bug fixed; major kinetic difference.** A 50-ms cooldown can suppress dense rebinding. |
| AA valency | Pair matrix plus geometric neighbors; the effective limitation is indirect. | `max_aa_bonds_per_actin` constrains pair bonds. Local links also have per-pair link count, site exclusion, and global pool. | Current pilots can be substantially more valency-limited. |
| `max_strong_actin_bonds` | Old 2-D paths implement an enforcement routine; the exact 3-D high-order path does not rely on it. | The field and declaration remain, but no active 3-D call/definition was found. | Configuration trap; do not treat it as an operative 3-D constraint. |
| Steric pre-phase | Calls only myosin-myosin exclusion. AA and AM contacts in the dense random initial state are retained. | Calls full volume exclusion, separating AA, AM, and MM overlaps before chemistry. | **Major initial-condition change.** Longer current relaxation erases the dense contact network that seeded the old run. |
| Production sterics | MM repulsion; AA repulsion is conditional on crosslink status; no general AM hard core. | AA, AM, and MM rod sterics are applied independently of bonding. | More physical; can strongly slow capture in a dense start. |
| Weak transverse bundling | `k_bundle_max=0.1`, ramped over 5000 steps in matched high-order runs. | Capability remains, but newer physical/scaffold-free pilots set it to zero. | A real legacy organizer and an uncontrolled difference. |
| Mobility | One scalar translational and one scalar rotational diffusion coefficient per species. | Legacy scalar mobility remains available. Physical pilots use rod `D_parallel`, `D_perp`, and `D_rot` from hydrodynamic drag with effective viscosity/crowding. | Intentional refinement; time scales cannot be compared without matching. |
| Thin/thick cooperativity | Absent. | Thin-filament activation, strong-head feedback, thick OFF/SRX recruitment, calcium protocols, and material-frame head azimuths are implemented/configurable. | New biology; should be disabled in a legacy bridge test. |
| Myomesin/titin | Code paths exist but were zero in `data0430`. | Code paths remain; current scaffold-free pilots also set them to zero. | Not an explanation for old high order. |
| Force/torque caps | Force caps follow maximum deterministic displacement. Angular input is capped in update paths, with limited diagnostics. | Caps are instrumented and should be monitored; dense mature-spring activation previously caused pervasive capping. | Retained with improved observability. |

## Parameter names that are misleading in explicit-AM runs

When `dynamic_am=true`:

- `k_am` does not set the head tether stiffness; use `am_tether_k`.
- `kappa_am` does not supply the legacy nematic AM alignment potential.
- `v_am` does not prescribe actin sliding; use the explicit motor parameters
  (`am_unloaded_velocity`, `am_stall_force`, `am_motor_k`, power stroke and KMC
  rates).
- `k_off` is parsed, but AA detachment is controlled by the implemented
  lifetime law; it should not be used as an assumed independent AA off-rate.
- `max_strong_actin_bonds` should not be assumed to limit current 3-D AA
  valency. Use the explicitly exercised pair/local AA controls.

Submission manifests should stop carrying inactive legacy parameters without a
comment, because their presence currently suggests a match that does not exist.

## Features implemented but not yet present in the main full-scale setup

Sparse spatial multi-ACTN2 is implemented and is the appropriate representation
for prebundling plus load-matured links, but `aa_local_multilink` defaults to
false. The principal full-scale explicit-AM manifests do not all enable it. The
dedicated `local_aa_assembly_pilot.sbatch` does enable it, together with:

- two candidate prebundle links per pair;
- site exclusion and a finite global ACTN2 pool;
- energy-gated association with `theta=0.1`;
- physical rod mobility;
- a two-AA-bond-per-actin limit;
- productive-opposed status 2; and
- thick-filament feedback.

That combination changes too many factors at once to identify why assembly is
lost. The local-link feature itself should be tested in a staged bridge rather
than only inside the full biological stack.

## Bugs and unresolved risks found during the audit

1. **Legacy cooldown was ineffective.** This is a genuine old-code bug, but
   fixing it changes kinetics enough that `tau_rec=0.05` must not be silently
   inherited in a legacy-match run.
2. **Explicit AM still has a candidate-lifecycle failure.** In the one-actin,
   one-myosin axial-exit diagnostic, all attachments disappear together near
   205 ms even though radial distance remains inside the tightened tether band
   and productive axial overlap is still nonzero. This indicates that an
   ownership/candidate retention gate remains inconsistent with the intended
   physical detachment rule. It must be fixed before interpreting full-scale AM
   occupancy or catch activation.
3. **Current documentation is stale in places.** For example, the ACTN2 design
   document describes a legacy default classifier, whereas the current CLI
   defaults to `productive_opposed`.
4. **Pre-phase time bookkeeping merits a test.** The model step counter is reset
   after pre-relaxation, while integrator start-step state may retain the end of
   the pre-phase in some paths. This could affect ramps or checkpoint metadata.
5. **Pair and local AA are different models.** Results must record
   `aa_local_multilink`; otherwise a trajectory cannot be interpreted from
   `k_aa`, `kappa_aa`, and lambda alone.

## What should be restored versus what should not

### Restore or emulate in controlled tests

1. A configurable **weak passive AM nematic organizer**, separate from active
   crossbridge torque, so its contribution can be swept down to zero. This is
   the cleanest bridge to the old model and directly tests the myosin-nematic
   deficit. It should act only within a biologically defensible local AM
   encounter/overlap zone, not globally outside capture as the old term did.
2. An explicit choice of **initial-contact protocol**: legacy MM-only prephase,
   full-sterics prephase, and gradual release/capture. This determines whether
   the old result is contact-seeded kinetic assembly.
3. The weak `k_bundle=0.1` field as an explicit positive-control arm, not as a
   hidden default.
4. A switchable **legacy AA registry control** beside the symmetric/current
   registry and local anywhere-binding model. This measures how much old order
   was imposed by the axial potential.
5. A **high-duty-ratio explicit-head bridge** that approaches continuous AM
   occupancy while retaining explicit head ownership and force bookkeeping.
6. A no-cooldown arm (`tau_rec=0`) beside the physical cooldown to quantify the
   legacy bug's contribution.

### Do not directly restore

1. The legacy AM endstop formula. It mixes absolute projection, a hard-coded
   `0.4*Lm` scale, and a propulsion-like axial force. Replace its function with
   correct head-bearing half zones, a motorless bare zone, polarity-aware head
   ownership, and loss of attachment after productive overlap ends.
2. Load-dependent stabilization of all status-1 links as the sole default. Keep
   it as a legacy/diagnostic toggle because it tests the proposed parallel-link
   trapping mechanism.
3. The asymmetric, pair-order-sensitive AA axial wall as the physical final
   model. Preserve it only as a legacy control.
4. MM-only relaxation as the sole production initial condition. It is useful as
   a mechanistic control, not a physical conclusion.

## Decisive staged comparison

Use identical seeds, particle counts, box, output cadence, and initial generated
coordinates. Measure actin and myosin nematic order, `phi_3d`, AM occupancy and
turnover, AA status populations, local-link count, cap fraction, and contact
network evolution.

### Stage A: executable reproduction

- Run the clean `30829b0` executable at A=2100, M=1000, lambda=`1e-3` for at
  least the documented seeds.
- Confirm the historical order trajectory, not only its endpoint.
- Archive executable hash, command line, seed, and initial frame.

### Stage B: exact coarse-model bridge in the current executable

- `dynamic_am=false`;
- legacy scalar mobilities;
- legacy AM/AA spacing and normalized-equivalent stiffness;
- `stabilize_status1=true`;
- `tau_rec=0`;
- `k_bundle=0.1`;
- thin/thick feedback off;
- pair AA, not local multilink;
- legacy MM-only versus current full-sterics prephase as the only factorial.

This establishes whether code evolution outside the intended physics changed
the legacy result.

### Stage C: one-factor migration ladder

Starting from the reproduced bridge, change exactly one row at a time:

1. full AA/AM/MM steric prephase;
2. physical production sterics;
3. symmetric current pair-AA registry;
4. real refractory time;
5. status-2-only catch stabilization;
6. explicit AM KMC at matched high duty ratio;
7. remove passive legacy AM nematic torque;
8. remove transverse bundling;
9. enable physical mobility;
10. replace pair AA with sparse local AA;
11. add thin feedback;
12. add thick feedback.

The order of steps 6-8 should also be tested as a small factorial because AM
intermittency, passive torque, and bundling can compensate for one another.

### Stage D: physical explicit model

Only after the migration ladder identifies the order-generating terms should
the full spatial-head, local-ACTN2, physical-mobility model be interpreted. Its
minimum acceptance criteria are:

- no false AM geometric detachment while productive overlap remains;
- radial tether retention without pathological force capping;
- KMC transition statistics matching specified rates;
- both AA status populations with status-2 episode durations long enough for
  catch stabilization to act;
- myosin alignment drive exceeding rotational diffusion on the assembly time
  scale; and
- lattice maintenance plus annealing from perturbed configurations.

## Bottom line

The legacy result is valuable because it demonstrates that this geometry can
assemble, but its mechanism is a coupled coarse-grained sorting machine. The
current model should build from it through a controlled migration ladder, not
attempt to recover `phi_3d > 0.8` by tuning lambda inside a simultaneously
changed AM, AA, steric, mobility, and initialization model.
