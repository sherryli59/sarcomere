# Sparse myomesin production run: 300,000 steps

## Final run design

The production comparison retains three matched seeds at full A800/M320 scale:

| Tasks | Arm | Seeds | `k_mm` | Module force limit |
| --- | --- | --- | ---: | ---: |
| 0--2 | no-myomesin control | 831704--831706 | 0 pN/um | inactive |
| 3--5 | sparse myomesin module | 831704--831706 | 240 pN/um | 40 pN |

Each trajectory uses:

- 300,000 integration steps at the converged `dt = 2.5e-6 s`;
- 750 ms of simulated time;
- output every 2,000 steps (5 ms), giving 151 frames including step zero and
  the final frame;
- a 2,000-step (5 ms) pre-phase and 8,000-step (20 ms) interaction ramp;
- 8 CPUs and 16 GB per task, a 12-hour wall-time limit, and Slurm
  requeue/resume support;
- at most two concurrent tasks (16 CPUs total).

The only matched-arm difference is the active sparse myomesin module. All other
physical and kinetic parameters are identical to the calibrated production
pilot documented in `SPARSE_MYOMESIN_PRODUCTION_PILOT_20260818.md`.

Active Slurm array: `39580183`, tasks `0-5%2` on `owners`.

Independent final audit: `39580184`, submitted with an `afterany`
dependency on the entire array.

## Why the timestep changed

The earlier `dt = 1e-5 s` setup was not numerically trustworthy. Its 5 nm
Brownian trust radius imposed an effective actin-force clipping threshold of
about 7.1 pN, and at 25 ms it clipped 510.5 of 800 actins per integration step
(63.8%) and 24.7% of myosins per step. Those counts represented actual capped
dynamics, not corrupt HDF5 output.

Matched 25 ms convergence tests gave:

| `dt` | Actin cap fraction/step at 25 ms | Myosin cap fraction/step | Result |
| ---: | ---: | ---: | --- |
| 10 us | 63.82% | 24.74% | rejected |
| 5 us | 21.47% | not selected | rejected |
| 2.5 us | 0.370% | 0.075% | selected |

The 2.5 us test remained finite, retained unit direction vectors to machine
precision, and passed the structural HDF5 checker. Parameter ablations showed
that the clipping came from the combined passive actin, ACTN2/titin, and A--M
force budget rather than one defective interaction. Raising only the
displacement cap reduced clipping but permitted larger excursions, so it was
not used as a substitute for timestep convergence.

Keeping the user-requested 300,000 steps at the converged timestep yields 750
ms of physical time. This is longer than the original 500 ms pilot, although
shorter than the rejected 3-second numerical setup.

## Output and reproducibility

Run root:

`/scratch/users/xsli/sarcomere/prod_sparse_myomesin_module_A800_M320_750ms_dt2p5us_300k_20260818`

The run root contains the frozen executable, its SHA-256 check, source/config
snapshot, parameter manifests, trajectories, analysis records, and Slurm logs.
The batch task refuses to overwrite an existing trajectory unless Slurm is
restarting the task or resume is explicitly authorized.

## Live integrity monitoring

Each task launches `analysis/monitor_live_production_hdf5.py`, which first
checks after one minute and then every five minutes. Compact records are
appended to `analysis/*.live.jsonl`; the final strict result is written to
`analysis/*.live_final.json`. An unhealthy final state fails the task.

The checker verifies finite coordinates, unit filament directions, monotonic
committed steps, consistent core dataset lengths, ACTN2 status-2 retention,
M/A order, and sparse-module endpoint validity, valency, pair uniqueness,
mismatch, radial spacing, configured force capacity, and frame alignment. An
active module arm must contain the myomesin offset/state stream. The checker
cross-checks the committed-frame marker
and ACTN2 identity/mechanics/motor-support offset tables against the saved
steps and their target row counts. It also computes the latest per-step force
cap fractions and scans every saved post-ramp interval. Any post-ramp interval
above 5% for either actin or myosin is a hard integrity failure, even if a
later interval recovers.

Concrete stop/fix/resubmit conditions are a simulator or metadata-validation
failure, NaN/Inf state, non-monotonic or inconsistent committed frames, invalid
module endpoints, duplicate M--M pair edges, valency above two, or the force-cap
gate. Scientific outcomes such as weak assembly or poor mismatch reduction are
reported without silently retuning the matched comparison.

The independent audit requires all six trajectories and final health records,
requires the last committed step to equal 300,000, checks strictly increasing
saved steps, verifies the 800/320 counts and arm-specific myomesin/ACTN2
parameters, and writes `analysis/production_300k_audit.json`. Because its
dependency is `afterany`, it will still run and report missing or failed outputs
if any array task does not succeed.

`analysis/test_production_hdf5_monitoring.py` supplies nine synthetic-HDF5 unit
tests for the force-cap pass/fail boundary, the six-run final audit's exact
151-frame requirement, committed-frame consistency, and rejection of a
mislabeled arm parameter. It also verifies that ACTN2 detachments without a
held motor identity are not counted as identity releases and that a historical
post-ramp cap spike fails even after the latest interval recovers. It also
verifies that an active myomesin configuration cannot pass without its
frame-offset stream. All nine pass.
The rebuilt `actn2_load_kinetics` and `myomesin_mband` CTest targets also pass;
they cover weak-state partner exchange/commitment/identity locking and local
capture, reciprocal force, independent axial/radial mechanics, force capacity,
turnover, and finite valency.
The complete 28-test CTest suite also exits successfully, covering the broader
initialization, A--M, ACTN2, titin, steric, neighbor-list, and mobility paths.

The final per-run analyzer additionally reports ACTN2 motor-pair selection,
switch, dropout, maturation-lock, and release event counts; final selected and
locked identity counts; active challenger counts, residence, and support; and
historical saved-row challenger occupancy, residence, and support; and the
first/max status-2 population. It decomposes the maturation gate into rows that
meet the configured residence time, filtered-support threshold, and both. These
are scientific observables rather than integrity gates, so a negative annealing
outcome remains visible instead of causing silent retuning.

The frozen production executable emits event-code 5 notifications for every
ACTN2 detachment, including rows whose old motor identity is invalid. This is a
logging-only issue: dynamics and checkpoints are unaffected. The frozen
analyzer counts code 5 as an identity release only when both old myosin indices
are valid and reports empty notifications separately. The working C++ source
has also been patched to suppress empty code-5 rows in future builds; the
production executable itself was deliberately not changed mid-comparison.

## Superseded submissions

- `39564537`: 50,000-step array, canceled pending when the requested length
  changed to 300,000 steps.
- `39572149`: 300,000-step array, canceled pending to embed the persistent HDF5
  sidecar.
- `39572811`: 10 us timestep array; tasks 0 and 1 were stopped at committed
  checkpoints after the force-cap pathology was identified.
- `39573007`: audit paired with `39572811`, canceled with that array.
- `39573805` and `39573813`: held replacement and audit, canceled without
  running after timestep convergence showed the dynamics were over-clipped.
