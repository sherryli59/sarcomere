# Connector deletion optimization and comparable diagnostics

Admission, buffered deletion and full-connector repulsion remain always on.
No new runtime switches were introduced.

## Performance and correctness

The deletion query now rejects distant myosin centers before exact periodic
image enumeration, computing the conservative reach once per link. Mature
connector checks run in the existing parallel mechanics pass. Classification,
deletion, slot recycling and event insertion retain canonical order.

On a saved production frame, the reference and optimized kernels both delete
32 of 4,062 mature links against 320 myosins. All 100,000 additional randomized
periodic/nonperiodic geometry cases agree exactly. The serial pass takes
86.12 ms before and 33.20 ms after the distance rejection (2.59x faster).

Final whole-simulation timing: Slurm job **50013712**, host
`sh03-17n14.int`, four CPUs, seed 20260805, explicit projection off, 200 steps
from the same 2.0 s checkpoint. Both executables already include admission,
buffered deletion and full-connector repulsion.

| Measurement | Before optimization | After optimization |
|---|---:|---:|
| Wall time including startup and checkpoint I/O | 59.735 s | 41.890 s |
| Serial link-classification phase, cumulative | 19.990 s | 0.645 s |
| Parallel mechanics phase, cumulative | 4.936 s | 6.793 s |

The complete short run is **1.426x faster**, or **29.9% less wall time**.
The mechanics phase now includes the parallel geometry work. This is one
paired short benchmark on four CPUs, not a production-throughput estimate
for a different allocation.

All detachment counts match after remapping old mature-deletion code 3 to
new code 5: 34 contour, 92 stochastic and 116 mature-connector deletions.
The final poses and event values are not bitwise identical. Separate repeats
of an unchanged binary also diverged during the earlier timing runs; this
work does not certify exact trajectory replay or resolve that existing issue.

Validation: **44/44 active-tree tests**, **40/40 shared-repository tests**,
and the connector regression run directly with four threads all pass, with
assertions enabled. The repository build also exposed and fixed an earlier
redundant dwell-counter reset that referenced an active-tree-only field.

Validated executable:
`$SCRATCH/sarcomere/connector_parallel_20261028/sarcomere`

SHA256: `8167e9a5e933912125ca1cd5f684e45a74b7fcb20527356b773821f9dbd3c52d`.

## Detachment reasons

Code **5** now identifies mature full-connector exclusion, independent of
status at death. Code **3** retains legacy barbed-tip burial. New-file and
resume metadata identify the current rule. Analysis handles code 5 even
when a resumed file retains an older reason dictionary in `/meta`; mature
cohorts follow link identity and maturation history.

Historical ambiguous code-3 events cannot be classified retrospectively.

## Existing trajectories remeasured with fixed geometry

Every live two-arm connector/myosin/tip pair is considered, independently of
the selected force contact. The radius is fixed at 130 nm, with physical
myosin end planes and a 1 nm reporting threshold. Depth is the greatest
axial burial of any eligible connector point on that tip's half. Counts below
are mature link-tip pairs per saved frame; depths are pooled pair-frame
percentiles. All arms use the same 80 frames at 5 ms spacing, 2.0–2.4 s.

| Existing arm | Overlapping mature pairs/frame | Median depth | p90 depth |
|---|---:|---:|---:|
| Projection on, old gate (`r27zt0`) | 101.5 | 26.5 nm | 187.8 nm |
| Projection off, old gate (`r27zt0np`) | 271.2 | 82.4 nm | 462.3 nm |
| Projection off, admission gate (`r27zt0npg130`) | 234.1 | 63.3 nm | 302.0 nm |

Maximum depth is 750 nm for all three under this full-connector definition.
These are remeasurements of existing trajectories, not a new 0.4 s run with
deletion and repulsion changes. Recorded selected-tip diagnostics are retained
separately and are not substituted for the common geometry metric.

The final buffered-region mature-link occupancies are 32, 223 and 151,
respectively. The JSON also includes unique-link counts, inherited and
post-resume mature cohorts, midpoint-only references and per-frame readouts.
The batched implementation agrees with the scalar reference on randomized
periodic cases and is independently checked against dense point sampling.

## Saved evidence

- `parallel_step_benchmark.json`: final timing, binary hashes, phase profiles,
  reason counts and short-run geometric readouts.
- `deletion_kernel_benchmark.json`: reference-kernel equivalence and timing.
- `comparable_connector_metrics.json`: identical geometric analysis for all arms.
- `performance_active_ctest.log`, `performance_repo_ctest.log`, and
  `performance_parallel_regression.log`: validation logs.

An earlier timing job script was edited while running, causing the driver to
repeat and overwrite its first JSON summary. That job was stopped; the final
results above come from a fresh paired run with overwrite protection.
Production trajectories and pinned production binaries were not overwritten.
