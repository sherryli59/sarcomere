# Actin-turnover and myosin-activation assembly plan

Date: 2026-08-28

## Objective

Create two independent, interpretable controls over sarcomere assembly:

1. actin turnover as a topology-remodeling/annealing control;
2. myosin activation as a motor-availability control spanning healthy-like to
   HCM-like behavior.

The desired output is not just one endpoint score.  Each condition should map
to an assembly phenotype containing:

- transverse ACTN2 percolation;
- M-order;
- Z-order;
- network size, persistence, and link redundancy;
- motor and ACTN2 state variables that explain the phenotype.

No production-scale factorial should be launched until the two mechanisms pass
their separate pilot gates.

## What the current control means

The A1600/M320 control uses tonic calcium (`am_cycle_period=0`) and thin
activation, but `am_thick_feedback=false`.  In the implementation this converts
all OFF/SRX heads to detached-ON, so all 9 heads per half are available.  It is
therefore a useful full-availability upper bound, not a healthy-like motor
state.

At the latest inspected control frame, mean `q_i` is about 0.65 and the mean
strong-head count is about 0.46 per actin.  Half-filament tension is very broad:
the median is zero, the upper quartile is about 2.9 pN, and the 90th percentile
is about 4.5 pN.  The default thick-filament tension scale of 0.1 pN would
therefore make recruitment nearly switch-like and saturated under an ordinary
single-head load.  It must be calibrated before a healthy/HCM comparison.

## Assembly observables

Measure the following at every saved frame and summarize both onset and the
last 20% of each trajectory.

### Primary outputs

- `P_perp_any`: an ACTN2 component has a periodic winding in y or z;
- `P_perp_both`: the same component winds in both y and z (the strong
  transverse-percolation criterion);
- winding persistence: fraction of tail frames satisfying each criterion;
- largest-component fraction and transverse extent;
- M-order and Z-order, including tail median and tail slope.

Compute graph quantities twice: once using every live ACTN2 link and once using
only mature/status-2 links.  Collapse multilinks to unique actin pairs for
topology, while also reporting molecules per unique edge and the fraction of
ACTN2 molecules placed on duplicate edges.

### Mechanistic outputs

- mean and distribution of thin activation `q_i`;
- OFF/SRX, detached-ON, weak, and strong myosin-head fractions;
- strong heads per actin and half-filament tension distribution;
- ACTN2 status-1/status-2 counts, formation/promotion/loss rates, and lifetimes;
- actin renewal-event count, successful-placement count, and links/heads lost
  per renewal event;
- force-cap rates, hard-contact corrections, minimum separations, and all
  finite-state checks.

Do not reduce these outputs to a single arbitrary "assembly score" during the
pilot.  First build a response surface.  A scalar tuning coordinate can be fit
later only if connectivity, M-order, and Z-order move together reproducibly;
otherwise report a Pareto map or target boxes.

## Axis 1: minimal actin turnover

There is currently no actin-filament turnover mechanism.  The first
implementation should be fixed-population, whole-filament renewal rather than
polymerization or length dynamics.

For each actin, use an independent Poisson hazard

```text
p_renew(dt) = 1 - exp(-k_turn dt).
```

On an event:

1. detach every AM head, ACTN2 molecule, and other explicit link involving that
   actin and release all associated site ownership;
2. clear pair status, recovery/memory, reserved-neighbor state, cached
   constraints, and local force/history fields involving that actin;
3. keep the actin count and length fixed, but redraw its center and orientation
   from the original initialization distribution, preserving the population's
   polarity balance;
4. reject insertions with hard steric overlap and record failed attempts;
5. reset `q_i` to its calcium-equilibrium value and reset the filament age;
6. rebuild affected broad-phase/neighborhood data before forces are evaluated.

This is an accelerated remodeling surrogate, not a claim that a full cardiac
thin filament physically disappears in a second.  Report it using both
`k_turn` and the realized renewed fraction
`1-exp(-k_turn*T_active)` so model-time comparisons remain explicit.

The hazard must initially be uniform and independent of degree, load, status,
or order.  Preferentially replacing underconnected filaments would bake the
desired percolation result into the rule.  Stress- or age-dependent turnover is
a later biological refinement.

Suggested controls:

```text
--actin_turnover=false
--actin_turnover_rate=<s^-1>
--actin_turnover_start_time=<s>
--actin_turnover_max_placement_attempts=<integer>
```

Checkpoint/resume must preserve actin ages and the renewal RNG stream exactly.
HDF5 should record renewal events with actin ID, step, old/new center and
direction, number of AM/ACTN2 links released, and placement attempts.

### Turnover verification gate

Before assembly runs, test a small periodic system and require:

- constant actin count, length, and polarity counts;
- observed event counts consistent with the dictated Poisson hazard;
- no stale ACTN2 pair state, AM site ownership, or hard-contact constraints;
- exact checkpoint/resume reproduction;
- no overlap spike after accepted insertion;
- `k_turn=0` is trajectory-identical to the current code.

### Turnover assembly pilot

Use the density-preserving A800/M160, 3 x 2 x 2 box as a screening system; its
transverse dimensions match the A1600 box.  First verify one no-turnover arm
against the early A1600/M320 behavior so that the shorter axial box is not
silently changing the conclusion.

After a 0.2-s common burn-in, test `k_turn = 0, 0.1, 0.35, 1.0 s^-1` for a
further 0.8-1.0 s, using three paired seeds.  These levels renew approximately
0%, 10%, 30%, and 63% of filaments per second.  Use the current full-availability
motor setting in this pilot so turnover is the only changed mechanism.

Expected qualitative behavior is potentially nonmonotonic: no turnover may
trap redundant local bundles; intermediate turnover may allow partner exchange
and transverse annealing; excessive turnover should erase mature links and
reduce M/Z order.  A nonmonotonic result is informative rather than a failed
tuning knob.

## Axis 2: healthy-like versus HCM-like myosin activation

Do not represent HCM by increasing myosin count or stall force.  Those changes
confound motor availability with density/mechanical strength, and the existing
M640 run already shows that more myosin can mainly increase force clipping.

The primary HCM axis should be reduced OFF/SRX sequestration: change the basal
OFF-to-ON and ON-to-OFF rates at fixed head count, stall force, attachment
chemistry, calcium drive, and thin-filament parameters.

Parameterize the zero-tension ON fraction as

```text
p_ON,0 = k_recruit,0 / (k_recruit,0 + k_fold,0).
```

Hold the total zero-tension switching rate near 20 s^-1 while changing the
ratio, so availability changes without also changing the relaxation timescale:

| label | p_ON,0 | k_recruit,0 | k_fold,0 |
|---|---:|---:|---:|
| healthy-like | 0.10 | 2 s^-1 | 18 s^-1 |
| intermediate | 0.30 | 6 s^-1 | 14 s^-1 |
| HCM-like | 0.50 | 10 s^-1 | 10 s^-1 |
| full-ON bound | 1.00 | thick feedback disabled | all heads available |

Initialize each arm at its zero-load equilibrium fraction.  These labels are
model phenotypes, not mutation-specific claims.  A later mutation-specific
study can separately change strong-state detachment/duty ratio, but it should
not be mixed into this first activation axis.

Keep tonic calcium for the first assembly pilot.  Pulsatile calcium is a later
validation because heartbeat duty cycle and diastolic relaxation add a second
activation variable.

### Motor pilot 1: cheap kinetic calibration

Use the existing force-clamp recruitment rig and
`analysis/check_thick_recruitment.py` before any network run.

1. Measure the dictated and observed ON fraction and transition hazards at
   clamp tensions 0, 1, 3, and 6 pN.
2. Compare tension scales 1, 3, and 6 pN.  Do not use the current 0.1-pN default
   for production assembly unless calibration shows it does not saturate.
3. Require measured hazards/ON fractions to agree with theory within the
   existing statistical tolerance and require ON fraction to relax after load
   release.
4. Select a scale that spans recruitment across the observed 0-5-pN range
   rather than becoming all-ON after one loaded head.

### Motor pilot 2: density-preserving assembly screen

On A800/M160, compare healthy-like, intermediate, HCM-like, and the current
full-ON upper bound, with turnover off.  Use the calibrated tension scale,
three paired seeds, a 0.2-s burn-in, and 0.8-1.0 s of measured assembly.

The pilot passes if the conditions produce separated ON/strong-head occupancy
distributions without numerical failure or pathological force clipping.  Only
then ask whether activation changes percolation or M/Z order.  If the ON states
do not separate, fix recruitment kinetics before interpreting structural null
results.

If availability changes but structure does not, run one secondary duty-ratio
pilot at fixed healthy-like availability by lowering `am_k_strong_off` from
160 to 80 s^-1.  Treat that as a distinct mechanism, not part of the initial
healthy/HCM label.

## Interaction pilot

After selecting one useful intermediate turnover rate and a calibrated motor
scale, run a small factorial on A800/M160:

| turnover | motor availability |
|---|---|
| off | healthy-like, intermediate, HCM-like |
| selected intermediate | healthy-like, intermediate, HCM-like |
| high/destructive control | healthy-like, HCM-like |

Use three paired seeds.  This tests the central interaction hypothesis:
activation creates and loads candidate structure, while turnover determines
whether it remains kinetically trapped or can anneal into a transversely
percolated network.

Analyze time-to-first winding, winding persistence, M/Z tail values, unique
edge creation, duplicate-edge fraction, and microscopic motor state.  A useful
control regime should change assembly without simply driving force caps or
continually erasing status-2 links.

## Production validation

Promote only three phenotypes from the pilots: low, intermediate, and high
assembly.  Validate them in the full A1600/M320, 6 x 2 x 2 system for 5 s with
at least 6 independent seeds (preferably 8-12 for percolation probability).

Use the spatial broad-phase executable, identical ACTN2 parameters
(`links/pair=3`, stiffness reference 3, encounter rate 10), identical output
cadence, and paired seed sets.  The canceled divergent A2400 original-one-link
trajectory is not a control for this study.

Predeclare success as a reproducible response surface that spans:

- low to high probability of persistent transverse winding;
- a meaningful range of M-order and Z-order;
- stable particle geometry and acceptable force/contact diagnostics;
- mechanism-consistent changes in renewal, head availability, strong
  occupancy, and unique ACTN2 edges.

## Recommended execution order

1. Add and validate the actin-renewal mechanism with default off.
2. Run the cheap thick-filament force-clamp calibration and choose the tension
   scale.
3. Run the separate A800 turnover and activation screens.
4. Inspect the response curves before choosing the interaction matrix.
5. Run the small interaction pilot.
6. Promote only selected low/mid/high phenotypes to replicated A1600 5-s runs.


## Addendum 2026-08-28: renewal clearance had to be calibrated

The first end-to-end smoke run of the `turn1` arm fired one renewal event and
that event **failed placement after exhausting all 1000 attempts**.  Left alone
this would have made axis 1 a silent null: `k_turn` would have been swept over
four levels with essentially no filament ever being renewed.

The cause is not crowding, it is an inconsistent acceptance rule.
`_actin_renewal_position_is_clear` demanded the full nominal hard core
(`2*r_actin` between actins, `r_actin+r_myosin` to myosin), but the steric
potential is soft and force-capped, so the equilibrium population itself sits
inside that distance.  Measured on the actual configurations:

| configuration | incumbent actins passing the renewal test | median min a-a separation / nominal | proposal acceptance |
|---|---:|---:|---:|
| A800/M160 3x2x2 | 8.7% | 0.885 | 0/1500 |
| A1600/M320 6x2x2, 200 ms | 23.7% | 0.943 | 0.60% |

A candidate was therefore being held to a standard that ~90% of the existing
filaments already violate.

Fix: `--actin_turnover_clearance_ratio` (default 0.80) scales both clearances.
The default admits candidates that are no more overlapped than the incumbent
10th-25th percentile, and `--actin_turnover_max_placement_attempts` was raised
1000 -> 5000 because the mean accepted proposal now takes ~410 attempts and the
observed maximum was 2346.

Verification after the change, at an inflated `k_turn=200 s^-1` so the statistics
are meaningful in a short run:

- 322 successful renewals against 320 expected, Poisson z = 0.11;
- 0 failed placements;
- 97 AM heads, 65 ACTN2 molecules and 44 titins released and rebuilt;
- packing *improved* rather than spiking (median min a-a separation 0.0268 ->
  0.0297), so the "no overlap spike" gate passes;
- +13% wall clock at 200 renewals/s, i.e. negligible at the campaign's
  `k_turn <= 1`;
- `k_turn=0` remains **bitwise identical** to the pre-change binary across
  actin/myosin centres, directions and crossbridge status.

The ratio is a calibration of the acceptance rule to the model's own steady
state, not a physical parameter.  It is written to the HDF5 attributes and
required on resume (legacy trajectories promote to 1.0).  If the production
A1600 box turns out to be looser, the ratio can be raised; the campaign figure's
renewal panel annotates placement failures whenever they occur, so a regression
to the silent-null mode is visible rather than hidden.

## Addendum 2026-08-28 (2): Lx=6 is now mandatory

Standing instruction: **every run uses Lx=6**.  At the campaign density
(66.7 actin/um^3, 13.3 myosin/um^3) that means **A1600/M320 in 6x2x2**.  There is
no cheaper screening box, and density matching does not rescue a shorter one.

Reason: the model selects a 2.00 um sarcomere period.  In Lx=6 the myosin band
spectrum is sharply peaked at harmonic k=3 (0.478 at 1 s, 0.639 at 3.2 s).  A
2.00 um period in Lx=3 requires k=1.5, which PBC cannot host, so the A800/M160
3x2x2 box shows a flat spectrum (max 0.166 at k=1) and 3-5x lower axial order
than A1600 at matched time and density:

| | lcc | M order | Z order | links/actin |
|---|---:|---:|---:|---:|
| A800 3x2x2 @ 1 s | 0.149 | 0.167 | 0.084 | 2.53 |
| A1600 6x2x2 @ 1 s | 0.099 | 0.479 | 0.427 | 3.39 |
| A1600 6x2x2 @ 3.2 s | 0.141 | 0.631 | 0.537 | 4.00 |

Consequence for the 2026-08-28 campaign (array 41061060): the M/Z-order and
axial-percolation readouts of all 21 arms are void.  Connectivity readouts are
**not** void — largest-component fraction is far less sensitive to axial
commensurability, and the within-box motor comparison is valid.

Use `submit_assembly_control_pilots_Lx6_20260828.sbatch`.

## Addendum 2026-08-28 (3): the success criterion has to change

`P_perp_both` is **zero in every run measured** — three box sizes, two link
stoichiometries, out to 5 s.  A predeclared binary criterion no arm can reach
makes every arm read as a null regardless of mechanism.  Two candidate causes
are already excluded:

- molecular redundancy: `aa_prebundle_links_per_pair=1` takes the duplicate
  fraction 0.63 -> 0.000 and `lcc` only 0.122 -> 0.133;
- link range: measured capture geometry is 59 nm transverse / 104 nm axial,
  100% antiparallel, i.e. correct alpha-actinin physics that must not change.

The blocker is barbed-end packing uniformity.  At A1600 3.2 s each actin has
4.30 antiparallel partners available inside the real capture envelope but only
2.90 linked, with **31% having no partner in reach** and **64% saturated at
`max_aa_bonds_per_actin=4`**.  Voids next to jammed pockets, where a real Z-disc
percolates because its lattice is uniform.  Raising density makes it worse
(A2400: 1.89 available, 46% with none), so the lever is registration.

Replace the primary outcome with graded measures that actually vary:
largest-component fraction, its transverse span, and the gap to the
random-placement expectation.  Untested levers, cheapest first:
`max_aa_bonds_per_actin` 4 -> 6, then axial registration of barbed ends.

Note the caveat on the ~4.5 mean-degree percolation threshold: that is the 2D
random-geometric result and the Z-disc here is a finite-thickness slab, so treat
it as "near the transition", not a precise number.

## Addendum 2026-08-31: barbed-end lateral repulsion

### Why this term exists

ACTN2 is a purely attractive crosslinker with no excluded volume of its own.
Its rest length (`aa_preferred_centerline`, 55 nm) is a *bond* property, so it
constrains only pairs that are already linked; per barbed end there are on
average 1.34 linked partners inside that distance and **2.72 non-partners**,
which nothing pushes away.  The only all-pairs term is the 30 nm actin steric,
which is soft and force-capped.  The result is progressive lateral collapse,
measured within a single run as links accumulate:

| ACTN2 links | barbed-end NN spacing | condensation vs mean |
|---:|---:|---:|
| 0 | 41.2 nm | 4.3x |
| 417 | 38.3 nm | 5.3x |
| 1647 | 30.5 nm | 7.2x |
| 2495 | 26.1 nm | 9.8x |

The collapse is self-defeating: each link shortens the reach of the next, so a
crosslinked cluster contracts instead of spreading.  Tiling one Z-disc face at
this filament count (603 barbed ends in 2x2 um) needs ~81 nm spacing, and
ACTN2's transverse reach is 93 nm, so even spacing *would* connect.

Offline verification on the real final frame (relaxing only transverse
barbed-end coordinates, per Z-disc, then rebuilding the link graph with the
measured envelope):

| condition | <k> | lcc | fraction of the 1/n_Zdiscs ceiling |
|---|---:|---:|---:|
| as-is (26 nm) | 5.79 | 0.194 | 58% |
| 70 nm spacing | 2.22 | 0.296 | 89% |
| 80 nm spacing | 2.20 | 0.306 | 92% |
| 80 nm + polarity ordering | 2.60 | 0.328 | 98% |

### What was ruled out first

- **Bond cap.** Removing `max_aa_bonds_per_actin` entirely and doubling the
  capture radius (<k> 3.1 -> 13.6) left `lcc` at 0.194.  Link budget is not the
  constraint.
- **`aa_zdisc_tip_*`.** That is the myosin-tip/Z-patch axial end-stop, not a
  filament-spacing term; median contact count 0 and max force 7.2 pN against a
  20 pN cap, i.e. working and not saturated.
- **Myosin crowding.** Actin *centres* sit 67 nm from the nearest myosin axis,
  but the *barbed ends* sit 337 nm away, out in the open Z-disc gap.  Myosin
  does not confine them.

### The term

`compute_barbed_end_lateral_repulsion` in `interaction.cpp`.  Linear soft
repulsion on the transverse separation of two barbed ends
(`centre - L/2 * direction`; `direction` points barbed -> pointed):

- applies to **every** barbed-end pair, both polarities, because the collapse is
  driven by unlinked neighbours;
- force acts **only** in the plane transverse to the pair's common axis, so
  axial registration is left to the register terms;
- gated on axial proximity (`aa_barbed_lateral_axial_window`, default 250 nm) so
  ends in different Z-discs do not push on one another;
- reciprocal, with torques applied at the material barbed-end caps via the
  existing `steric_contact_orientation_drive`;
- additive to the capsule sterics rather than replacing them.

Controls, all defaulting to inert:

```text
--aa_barbed_lateral_range=<um>          # 0 disables; try 0.070-0.080
--aa_barbed_lateral_k=<pN/um>           # 0 disables
--aa_barbed_lateral_force_max=<pN>      # <=0 uncapped
--aa_barbed_lateral_axial_window=<um>   # default 0.250
```

Written to the HDF5 attributes and required on resume; legacy trajectories
promote to the inert defaults.  Unit test `barbed_lateral_test` covers
inertness when unconfigured, transverse-only direction, reciprocity, the range
and axial gates, force-cap saturation, coincident-end degeneracy, and inclusion
of parallel pairs.

### Expected outcome and the known risk

Predicted from the offline check: per-Z-disc capture should rise from ~0.48
toward ~0.9, and global `lcc` from 0.19 toward the 0.33 ceiling.  It is **not**
expected to deliver `P_perp_both`: at 80 nm spacing the graph is marginal
(<k> ~ 2.6) and too tree-like to wrap two axes reliably.

Known risk, untestable offline: a repulsion strong enough to hold 80 nm at the
barbed ends fights the AM binding that pins actin *centres* 80 nm from myosin
axes.  That tension may surface as force-cap clipping rather than clean
spacing.  Watch `cumulative_force_cap_events` and the realised barbed-end
spacing together - if spacing does not open while clipping rises, the two terms
are fighting and the stiffness must come down.
