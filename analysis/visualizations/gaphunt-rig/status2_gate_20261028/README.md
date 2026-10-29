# Matched production resume: status-2 connector gate

Submitted 2026-10-28 for the requested production comparison.

- Production job: **49948644**, `production.sbatch`.
- Readout job: **49948755**, `analysis.sbatch`, dependency
  `afterok:49948644`.
- New arm: `r27zt0npg130`.
- Primary comparison: `r27zt0np` (projection off, previous maturation rule).
- Additional reference: `r27zt0` (projection on, previous maturation rule).

The production job was pending at the last inspection. These are queued jobs,
not completed results. The readout executes automatically after successful
production completion; it does not schedule another agent session.

## Protocol

Use the exact source `hph27m24_seed20260805.h5` and original saved argv,
resuming step 200000 (2.0 s) to 240000 (2.4 s). Keep seed 20260805,
27 heads/half, Lx=6, lifetime_coeff=0.4, k_mm=0, titin_k=1e-6,
requires_link=true and hard projection=false. The only intended model change
relative to `r27zt0np` is the new status 1 -> 2 connector exclusion.
The initial 3881 mature links retain their existing retention semantics.

First stage an independent source copy for a 500-step smoke test. Verify
resume advancement, the unchanged initial saved pose and checkpoint gate
metadata. Only after that succeeds, stage another original-source copy for
the full 40000-step continuation. Smoke output is retained separately.

Pinned binary: `$SCRATCH/sarcomere/gaphunt_20261023/bin/sarcomere_status2_gate_20261028`,
copied inside the production job from the validated 44/44-test build.
Required SHA256:
`dfa63b085df368246e93f8c27f121f4884b6f1e08fa0fc3f77379b65bb84c47b`.

The preparation script refuses duplicate output paths and records the source
HDF5 hash, binary hash, original and actual argv, host/CPU, runtime, model
settings and resume history in `run_manifest.json`.

## Outputs

Campaign root: `$SCRATCH/sarcomere/gaphunt_20261023`.

- Trajectory: `data/r27zt0npg130_seed20260805.h5`.
- Batch log: `logs/status2g130-49948644.out`.
- Simulation log: `logs/r27zt0npg130_seed20260805.simulation.out`.
- Analysis log: `logs/status2readout-49948755.out`.
- Measurements: `torque_out/status2_gate_20261028/`.
- Small completed summary/provenance files are also copied to `results/`
  beside this README after analysis finishes.

Readout covers projection counts, selected-patch axial overshoots, barbed-tip
penetration, detachment reasons and rates by status, live/mature populations,
maturation fraction, mature lifetimes by A-M angle (KM means restricted to
300 ms), mature load versus A-A/A-M angle, zero-load fraction, Z capture,
barbed spacing, Z-band order, filament/pair alignment and nematic order,
force-cap frequency and runtime phase profiles. Inherited mature-link
survival and death reasons are reported separately.

The analysis pipeline was exercised on both completed reference arms in a
compute allocation before the new trajectory became available. It exposed
legacy myosin counter resets at resume/preemption; force-cap frequencies
are reconstructed across those resets rather than subtracting endpoint
counters (which would produce negative frequencies).

The established projection analysis now accepts an explicit output path
and parses only the final `exit=... wall_s=...` footer for wall time.
The alignment analysis recognizes the new arm. Lifetime analysis handles
an empty mature cohort without a shape error.

## Interpretation and continuation

This is one seed and 0.4 s of resumed dynamics. Lifetime results concern links
born after the resume, and maturation angles use the next 5-ms saved pose,
following the established analysis. The gate does not retroactively reject
inherited mature links or prevent later movement of a mature connector into
the excluded region. Total runtimes on different CPU models are not a
controlled performance comparison.

Before continuing on Sherlock, read `/etc/agents/AGENTS.md` and the relevant
topic files. Use jobs for substantial computation and I/O. Scheduler status
checks must be at least 60 seconds apart. Cancel only exact job IDs.
