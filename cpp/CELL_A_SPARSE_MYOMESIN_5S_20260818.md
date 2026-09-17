# Cell A plus sparse myomesin comparison (2026-08-18)

## Question

Does the sparse, local, finite-valency myomesin module improve M-band
registration when added to the established Cell A assembly model, while
retaining Cell A's spontaneous Z-band assembly?

This comparison deliberately does not use the altered ACTN2 and steric
settings of the 750-ms six-run pilot. It reproduces Cell A at its original
`10 us` integration step and changes only the sparse myomesin interaction in
the treatment arm.

## Matched arms

Both arms use 800 actins, 320 myosins, the same `4.825 x 2 x 2 um` periodic
box, seed `20260805`, a 5-ms sterics-only prephase, a 20-ms interaction ramp,
and 5 seconds of simulated time.

| Arm | `k_mm` | Sparse myomesin model |
| --- | ---: | --- |
| Cell A control | 0 pN/um | disabled |
| Cell A + myomesin | 240 pN/um | local dynamic modules, 40-pN force limit |

The treatment module has transverse stiffness `300 pN/um`, radial rest
spacing `110 nm`, axial capture mismatch at most `250 nm`, radial capture
window `110 +/- 40 nm`, and two sites per myosin filament. One simulated edge
is a coarse-grained myomesin module, not one molecule.

## Cell A settings restored

Relative to the currently running 750-ms production pilot, this comparison
restores the original Cell A values:

- `k_aa`: 3000 to 1000 pN/um;
- myosin steric stiffness: 300 to 1000 pN/um;
- effective M-M steric contact: explicit 95.3 nm to the natural hard-core
  distance (`myosin_exclusion_distance=0`);
- ACTN2 maturation: `motor_transmitted_register` to Cell A's
  `productive_register`;
- ACTN2 force hysteresis: 1.0/0.5 to 0.25/0.10 pN;
- ACTN2 association energy gate: enabled to disabled.

The numerical schedule now matches Cell A exactly: `dt=10 us`, 500,000 total
steps, a 500-step (5-ms) sterics-only prephase, a 2,000-step (20-ms)
interaction ramp, and frames every 500 steps (5 ms). Each trajectory therefore
contains 1,001 frames from physical step 0 through step 500,000.

An earlier `dt=2.5 us` smoke attempt revealed that this frozen executable
retains the hidden prephase offset in the integrator's internal save counter
even after resetting the physical step to zero. A save stride that does not
divide the prephase can therefore omit the intended final physical state. At
`dt=10 us`, the 500-step save stride exactly divides the 500-step prephase. The
launch script rejects misaligned alternatives.

Each long arm requests 32 CPUs, 64 GB, and a 48-hour limit. The original Cell A
job took about 5.6 hours for the same 500,000 steps on 32 CPUs. The treatment
also pays for sparse-edge search and mechanics, so the limit retains generous
headroom.

## Validation and decision readout

The two-task smoke array uses the full particle counts and both mechanical
arms for 25 ms before the long jobs are submitted. Every final trajectory is
checked for complete monotonic frames, exact model metadata, finite values,
and numerical-force-cap pathologies.

The live checker was also corrected for Cell A's classifier: the
`local_link_motor_support` table is required only for status mode 8
(`motor_transmitted_register`). Cell A uses status mode 7
(`productive_register`), which legitimately has link-mechanics offsets but no
source-resolved motor-support table. Tests cover both the allowed mode-7
absence and the required mode-8 presence.

At `dt=10 us`, original Cell A is intentionally force-cap dominated: its
historical maximum per-interval cap fractions are 0.6505 for actin and 0.4288
for myosin. The first 10-us smoke reproduced that profile closely (control:
0.6689 and 0.2598 at 25 ms) but failed the 0.05 thresholds designed for the
2.5-us convergence runs. The Cell A comparison therefore uses explicit
reference-envelope health limits of 0.75 for actin and 0.50 for myosin. These
limits do not claim timestep convergence; they detect excursions beyond the
known Cell A numerical regime while honoring the requested original timestep.

`analysis/audit_cell_a_sparse_myomesin.py` keeps numerical integrity separate
from the scientific outcome. It reports:

1. M-, Z-, and A-band order in each arm and their paired differences;
2. module mismatch reduction relative to each bridge's formation mismatch;
3. finite valency, duplicate-edge, and connected-component safeguards;
4. retention of mature ACTN2 links and preservation of Cell A's Z-band.

The sparse module is supported only if its bridges anneal mismatch, final
M-band order exceeds the matched control, Z-band order is not materially
damaged, and the topology remains sparse and local.

### Historical Cell A replication gate

Cell A replication is now an explicit prerequisite rather than an inference
from matching command-line parameters. Historical job `38011216` ended with
Z-band order `0.708`; the median over its final 20% was `0.699`. The current
auditor therefore requires the no-myomesin control to have both final Z-band
order and final-20% median Z-band order at least `0.60`. This requires a
strong, persistent periodic pattern while allowing a different stochastic
multithreaded trajectory.

`analysis/check_production_myomesin_pilot.py` reports the endpoint, tail
median, and tail minimum. `analysis/audit_cell_a_sparse_myomesin.py` exposes
the pass/fail result as `control_reproduces_cell_a_z_band`.

The standalone final verifier `analysis/validate_cell_a_replication.py`
compares historical Cell A directly with the corrected control, checks the
complete 1,001-frame/500,000-step trajectory and Cell A metadata, and writes
both a JSON decision record and aligned Z-order time series. Slurm job
`39688239` is dependent on the corrected production/recovery chain and will
run this gate automatically when that chain exits.

Paired treatment audit job `39687591` waits for both tasks of array
`39662780`, recomputes the control and myomesin summaries with the updated
tail metrics, and writes `cellA_sparse_myomesin_audit_v2.json`. This keeps the
control-replication decision separate from the treatment readout: the report
records global M-order, local mismatch annealing, finite valency, duplicate
edges, component size, mature ACTN2 retention, and Z-band preservation even if
the myomesin arm does not improve global M-order.

### Backward-compatibility check

`cpp/cell_a_backward_compat_20260818.sbatch` recovers the authoritative Cell A
argument vector directly from the historical job log, substitutes only the
output path and shortened duration/save cadence, and runs the historical and
current binaries sequentially with one OpenMP thread. The comparator
`analysis/compare_cell_a_backward_compat.py` requires exact equality for 24
coordinate, force, chemical-state, motor-head, ACTN2, and titin arrays. When
the treatment is active it additionally requires exact equality of myomesin
bridge rows, association/release events, frame offsets, and per-frame bridge
IDs (28 datasets total).

The frozen treatment-capable binary passed all 24 arrays exactly through the
complete 500-step prephase and all five saved states from physical steps 0 to
2,000 (20 ms). This includes stochastic binding and force evolution, and rules
out a change in initialization, active Cell A force laws, or random-number
state. The rebuilt no-barrier binary also passed all 24 arrays exactly at all
five saved states through physical step 2,000. Thus both the feature addition
and the scheduling-compatibility patch have end-to-end state-equivalence
evidence, not only unit-level zero-force evidence.

The audit did identify one execution-level compatibility leak: the sparse
myomesin implementation entered an unconditional OpenMP `single` region even
when `k_mm=0`. Its callee returned without changing model state, but the
region's implicit team barrier changed multithread scheduling. The call is now
guarded outside the OpenMP region by `k_mm > 0`, restoring Cell A's historical
zero-myomesin execution graph while leaving the treatment path unchanged. The
rebuilt executable is:

`/scratch/users/xsli/sarcomere/cellA_backward_compat_20260818/bin/sarcomere_kmm0_nobarrier_20260818`

SHA-256:

`5ce40655c5349dc5de079713ea29b32984ef5b4a803264d9e77201a63427f51e`

As a scheduling diagnostic (not an exact multithread reproducibility claim),
the maximum actin-coordinate difference from historical Cell A immediately
after pre-equilibration fell from `2.59e-5 um` with the unconditional barrier
to `7.95e-6 um` without it; the corresponding myosin difference fell from
`3.15e-6 um` to `9.70e-7 um`. Chemical states and bound-head identities were
identical in all three files at that state. This is the direction expected
when removing the execution-only perturbation.

The active-treatment preservation probe used the complete 176-argument
Cell-A-plus-module command with `k_mm=240`, one OpenMP thread, and the same
2,000-step physical interval. The pre-patch treatment binary and rebuilt
no-barrier binary matched exactly in all 28 datasets, including every sparse
myomesin bridge and event (`myomesin_retained_nobarrier/comparison.txt`;
interactive Slurm allocation `39693614`). Thus the compatibility guard changes
only the `k_mm=0` scheduling path; it does not merely leave treatment code
present, but preserves its realized stochastic trajectory exactly.

### Why local M-M mismatch did not immediately become global M-order

At 2.005 s the treatment had 312 bridges covering 319/320 myosins, so the
failure to propagate was not low bridge occupancy. Valency two was nearly
saturated (305 myosins had degree two), but it partitioned the M-M graph into
54 disconnected chains/cycles: median component size 5 and largest component
21. Each component was internally almost perfectly phase-registered
(`0.997` size-weighted order), while the component phases had order `0.154`,
matching the global M-order of `0.150`. The local 7.3-nm bridge mismatch was
therefore real but confined within many independently phased islands.

Myosins sharing an actin were coupled more strongly than random pairs but were
not forced onto one M-plane. At the same frame, 2,205 myosin pairs shared at
least one attached actin and had median M-phase mismatch 220 nm; the 1,047
pairs sharing a strongly bound actin had median mismatch 182 nm, versus 300 nm
over all pairs. Only 24% (any attached) and 29% (strong) were within 25 nm.
This is physical: two motors can bind different axial sites on a 1-um actin,
so shared actin transmits load without imposing equal thick-filament centers.

By 2.675 s the bridge graph had only coarsened from 54 to 49 components (the
largest remained 21), but component-phase order rose to `0.302` and global
M-order to `0.295`; within-component order remained `0.991`. Propagation is
therefore occurring slowly through the broader actomyosin/ACTN2 network, while
the finite-valency M-band graph remains the principal connectivity bottleneck.
Increasing `k_mm` or adding parallel springs to the same pair would stiffen an
already registered island without connecting islands. A physical connectivity
test should instead give the coarse-grained bare-zone module three or four
angularly distinct, local neighbor sites (with one edge per sector), while
keeping finite valency, local radial/axial capture, and maturation locking.

All 28 C++ tests and 28 targeted Python tests pass after this change. A new
one-seed, 5-second control using the rebuilt binary is job `39684207`; its
output is isolated under `full5s_nobarrier`. The original paired array remains
available as a before-fix scheduling comparison, and its myomesin arm remains
valid because the active-module execution path is unchanged.

The corrected run's observed throughput remained inside its 12-hour limit but
slowed as the interaction graph grew. Sherlock does not permit this user to
extend a running allocation. Recovery job `39688235_0` therefore waits for
`39684207_0`, opens the same frozen executable and configuration with
`--resume=true`, and targets the same absolute 500,000-step endpoint. If the
first job completes, the successor is a zero-step validation/postprocessing
pass; if it times out or is preempted, the successor continues from the latest
committed resume snapshot. Final validator `39688239` waits for the successor.

Final completion-audit job `39697351` waits for both the corrected-control
validator (`39688239`) and paired treatment audit (`39687591`). It writes
`cell_a_goal_completion.json` and exits successfully only if all of the
following are simultaneously proven: complete/sustained Cell A Z-order,
healthy corrected HDF5, frozen-binary path and SHA-256, paired-run integrity,
the exact 5-s treatment configuration, healthy finite-valency/unique/local
annealing topology, and the 28/28 exact active-treatment comparison across
the scheduling patch. Global treatment M-order is recorded as a scientific
outcome but is not substituted for the explicit treatment-retention gates.

## Reproducibility

Launch script: `cpp/cell_a_sparse_myomesin_5s_20260818.sbatch`

Frozen executable:

`/scratch/users/xsli/sarcomere/prod_sparse_myomesin_module_A800_M320_750ms_dt2p5us_300k_20260818/bin/sarcomere_sparse_myomesin_module_20260818`

SHA-256:

`2a133f81b13fcbc89f869eb1067f13984ab8be6aa03544244b1d776dbe85e9b6`

Output root:

`/scratch/users/xsli/sarcomere/cellA_sparse_myomesin_A800_M320_5s_dt10us_20260818`

The launch configuration, documentation, final auditor, HDF5 monitor, and
their tests are copied into `source_snapshot/` under this root before the long
array is submitted. Long-running jobs call the frozen snapshot rather than
mutable working-tree analysis files.

Initial owners-queue smoke array `39645960` was cancelled while still pending;
it never launched and created no output. Smoke array `39646773` used the
group-owned `rotskoff` partition and exposed the save-alignment issue above:
both simulator processes exited normally after healthy intermediate states,
then the post-run final-step/frame-count gate correctly failed both tasks.
Aligned replacement array `39648650` was cancelled while pending because no
group node had 32 free CPUs. The identical smoke gate was resubmitted with 8
CPUs per task as Slurm job `39648775`; the long jobs retain 32 CPUs each.

Aligned smoke `39648775` passed both arms with exit code 0, exact final step
10,000, six monotonic frames, finite geometry, and low force-cap fractions.
At 25 ms the module arm had 22 bridges, maximum degree 2, no duplicate pair,
and largest connected component 3/320. For bridges formed with at least 20 nm
initial mismatch, median reduction from formation was 85.4%; all 22 evaluated
bridges improved and 95.5% improved by at least 10%.

The `dt=2.5 us` full array `39650838` and its audit `39650841` were cancelled
after 42 minutes when the numerical schedule was changed to match original
Cell A. Their partial files remain preserved under the old `dt2p5us` root.

The first `dt=10 us` smoke array `39660922` reproduced the original Cell A
force-cap profile but failed the stricter 2.5-us monitoring threshold after
both simulator processes completed normally. Reference-envelope smoke array
`39662478` tests the documented Cell A-specific limits end to end.

The full 5-second array `39662780` is staged with
`afterok:39662478_*`, so it cannot start unless both smoke arms pass. Final
audit `39662789` is staged with `afterany:39662780`.
