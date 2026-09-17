# Sparse myomesin module production pilot (2026-08-18)

## Purpose

This is a production-scale pilot, not a parameter sweep. It asks whether the
calibrated sparse myomesin module improves M-band registration in the full
A800/M320 assembly model without creating a system-spanning myosin network or
damaging the actin/ACTN2 structure.

The tested myomesin interaction is one sparse, local, dynamic, finite-valency
M--M edge. Each edge represents a coarse-grained myomesin module rather than a
single protein. The selected calibration is:

- axial stiffness `k_mm = 240 pN/um`;
- transverse stiffness `myomesin_radial_k = 300 pN/um`;
- module force limit `myomesin_force_max = 40 pN`;
- radial rest spacing `110 nm`;
- axial capture mismatch at most `250 nm`;
- radial capture window `110 +/- 40 nm`;
- at most two module sites per myosin filament;
- local association rate `5000 /s`, zero-load release rate `1 /s`, and Bell
  distance `0.5 nm`.

The no-myomesin control uses the same geometry, sterics, ACTN2 chemistry,
motors, titin, seeds, and inactive myomesin parameters; its only intended
mechanical difference is `k_mm = 0`.

## Matched design

All trajectories use 800 actins, 320 myosins, a periodic
`4.825 x 2 x 2 um` box, `dt = 1e-5 s`, 50,000 steps (500 ms), and frames every
500 steps (5 ms). The array mapping is:

| Task | Arm | Seed | Axial module stiffness | Force limit |
| ---: | --- | ---: | ---: | ---: |
| 0 | control | 831704 | 0 pN/um | inactive |
| 1 | control | 831705 | 0 pN/um | inactive |
| 2 | control | 831706 | 0 pN/um | inactive |
| 3 | module40 | 831704 | 240 pN/um | 40 pN |
| 4 | module40 | 831705 | 240 pN/um | 40 pN |
| 5 | module40 | 831706 | 240 pN/um | 40 pN |

The full array is throttled to two simultaneous tasks: at most 64 requested
CPUs and 128 GB across this pilot. Each task has an eight-hour wall-time limit.

## Staged submission

The two-arm smoke array runs the full A800/M320 initialization and interaction
pathway for ten steps with one control and one module task. The full array can
be released only if both smoke tasks exit successfully.

- Initial smoke array: Slurm job `39555444`, tasks `0,3%1`, 20-minute limit;
  failed before simulator launch as described below.
- Initial full pilot: Slurm job `39555445`, tasks `0-5%2`; cancelled
  automatically when the smoke dependency failed.

The first submission stopped safely at the smoke gate: tasks `39555444_0` and
`39555444_3` exited before simulator launch because the Sherlock Bash version
treated expansion of an empty resume-option array as an unbound variable under
`set -u`. Consequently, dependent array `39555445` was cancelled without
running. The script now represents the resume option with a scalar Boolean and
appends `--resume=true` only inside an explicit conditional. No simulation
parameter or executable changed. The replacement job identifiers are recorded
below:

- Replacement smoke array: Slurm job `39564535`, tasks `0,3%1`.
- Replacement full pilot: Slurm job `39564537`, tasks `0-5%2`, with
  `afterok:39564535_*`.

This dependency checks executable/configuration integrity. It is not evidence
that the myomesin hypothesis passes; that decision is made from the completed
500-ms trajectories.

## Reproducibility and outputs

The frozen executable is:

`/scratch/users/xsli/sarcomere/prod_sparse_myomesin_module_A800_M320_500ms_20260818/bin/sarcomere_sparse_myomesin_module_20260818`

SHA-256:

`2a133f81b13fcbc89f869eb1067f13984ab8be6aa03544244b1d776dbe85e9b6`

The run root is:

`/scratch/users/xsli/sarcomere/prod_sparse_myomesin_module_A800_M320_500ms_20260818`

It contains `data/` trajectories, `logs/`, per-run `params/`, per-run JSON
summaries in `analysis/`, and the source/configuration snapshot in
`source_snapshot/`. The smoke trajectories use the separate root:

`/scratch/users/xsli/sarcomere/prod_sparse_myomesin_module_A800_M320_smoke_20260818`

The submitted script refuses to overwrite an existing trajectory. Resume is
allowed only when Slurm marks the task as a requeue of a preempted run. It also
checks the frozen executable hash, A800/M320 size and density, and output HDF5
metadata.

## Analysis and decision rule

Each completed trajectory is summarized by
`analysis/check_production_myomesin_pilot.py`. For the module arm it records
association and release events, exact module mismatch relative to bridge
formation, radial spacing, module force, valency, duplicate pair edges, and
connected-component size. For both arms it reports final M-, Z-, and A-band
order, retained mature ACTN2 links, and global numerical force-cap events.

The module is safe enough to consider for a longer production campaign only if:

1. every module trajectory respects valency two and has no duplicate pair
   edges;
2. no module network becomes system spanning (diagnostic guard: largest
   component no more than 45% of myosin filaments);
3. modules formed with at least 20 nm initial axial mismatch show a positive
   median reduction, and at least half improve rather than worsen;
4. module-arm M-band order improves reproducibly over its seed-matched control,
   rather than being driven by one seed;
5. final Z-/A-band order and mature ACTN2 retention are not materially worse
   than the matched controls;
6. no unexpected global numerical force caps occur, and module forces remain
   consistent with the configured 40 pN limit.

The analyzer also reports the fraction of modules reaching at least 10%
mismatch reduction. That stricter statistic is retained for interpretation,
but the production decision should use the full matched-seed picture because
new modules can form near zero mismatch and have little room for a 10% change.

Useful status command:

```bash
squeue -j 39555444,39555445 \
  -o '%.18i %.9P %.20j %.2t %.10M %.10l %.4C %R'
```

Submitted configuration:

`cpp/production_sparse_myomesin_module_pilot_20260818.sbatch`

## Single-seed scope change (2026-08-18)

The production decision was changed to one matched seed per experiment. For
array `39580183`, seed `831704` is the retained pair: control task 0 completed,
and module task 3 remains active. Pending control task 2, active module task 4,
and pending module task 5 were cancelled. Control task 1 had already completed
before this change; its file is retained as historical output but is not part
of the requested one-seed decision. No trajectory was deleted.

The obsolete six-run aggregate audit `39580184` was cancelled. The newer Cell
A comparison already consists of one seed (`20260805`) with matched control
and myomesin arms, so no change to that array was required.
