# Audit of newly matured ACTN2 links — 2026-10-28

The cohort is every link with a logged status 1 -> 2 transition AFTER step
200000, excluding all links already mature at the resume. It includes both
previously weak links and links whose initial 0 -> 1 birth occurred after
the resume.

## Full production trajectory, independent saved-pose geometry

| Quantity | Old gate, projection off | New gate, projection off |
|---|---:|---:|
| Post-resume maturation events | 9497 | 9200 |
| Newly mature links observed in a saved frame | 8983 | 8740 |
| Died before any saved mature pose | 514 | 460 |
| Excluded at first saved mature pose | 644 (7.17%) | 420 (4.81%) |
| First-pose depth >50 nm AND radius <120 nm | 276 | 50 |
| First saved mature pose clear, later pose excluded | 913 | 932 |
| Newly mature overstretch deaths | 1328 | 1307 |
| Newly mature burial deaths | 185 | 84 |

The independent Python audit clips the full connector to each finite myosin
axial slab, then minimizes radial distance on that clipped connector. It
checks every myosin and distinguishes weak, inherited mature and newly
mature cohorts. The observed production connector lengths and box dimensions
allow a single nearest midpoint image; the code asserts that conservative
bounding-sphere condition for every connector rather than assuming it.

The production data PROVE later entry for the 932 links first observed clear
and subsequently observed inside the cylinder. The 420 first-pose overlaps
do NOT prove admission-time violations: the saved pose follows integration,
and is 0–4.99 ms after the actual event. The median first-sample delay is
2.45 ms.

## Exact-step diagnostic continuation

Slurm job 49979927 built an isolated diagnostic translation unit with a logger
immediately after the final maturation gate. It recorded the connector
endpoints and ALL myosin positions/directions at each admitted transition,
before later filament motion. The same model settings and source checkpoint
were used for 2000 steps (200000–202000), with 4 threads permitted by the dev
partition's per-user CPU limit.

- 554 admitted promotions were checked independently in Python.
- **Zero** admitted geometries violated the 130 nm + finite axial overlap rule.
- 5046 prospective promotions were vetoed.

Example from this diagnostic continuation: link 164164, actins 1168/687,
relative to myosin 260:

| Instant | Minimum connector radius within axial overlap | State |
|---|---:|---|
| Maturation, step 200227 | 136.733 nm | Admission allowed |
| Next saved pose, step 200500 (2.73 ms later) | 114.021 nm | Still status 2 |

At that later pose the closest in-slab connector point is 92.566 nm behind
an end plane. This is direct evidence that a newly mature link can pass
the gate and then move into the excluded region while retaining status 2.

### Repeatability limitation

The diagnostic is NOT an exact reproduction of the original production
trajectory: its post-resume event set and poses diverged, and the script's
repeatability assertion correctly failed. The gate geometry audit itself
passed all 554 admissions. The original production smoke test and full run,
both using the original binary and 16 threads, also have different poses by
step 200500; therefore changing thread count is not the only repeatability
limitation observed. The cause of that existing repeatability issue has not
been established.

Do not use this diagnostic's matching numerical link IDs as proof of exact
geometry for those IDs in the original production run. The example table
above uses the diagnostic's own event-time trace AND its own saved pose.
The full production cohort statistics above use only the production files.
Exact admission geometry was not saved for all 9200 production maturations,
so those events cannot all be retrospectively certified from 5-ms snapshots.

## Conclusion

No bypass of the implemented maturation gate was found in the source paths
or the independent exact-step diagnostic. The admission condition is being
enforced in the checked continuation. The source explicitly applies it only
when old status != 2 and proposed status == 2. After admission, the same
130 nm cylinder rule is not enforced as a retention/deletion condition.

Consequently the gate improves the early geometry and burial-death counts,
but later entry persists among newly mature links and their overstretch
deaths barely change. This cannot be explained solely by inherited mature
links. A further change would need to address continued exclusion/retention
and the mechanical causes of overstretch; this audit makes no such model
change.

## Inputs and artifacts

Production files, seed 20260805:

- `$SCRATCH/sarcomere/gaphunt_20261023/data/r27zt0np_seed20260805.h5`
- `$SCRATCH/sarcomere/gaphunt_20261023/data/r27zt0npg130_seed20260805.h5`

Full saved-pose audit JSON:
`$SCRATCH/sarcomere/gaphunt_20261023/torque_out/status2_gate_20261028/maturation_audit/`.
Audit script: `../audit_new_mature.py`.

Exact-step diagnostic:
`$SCRATCH/sarcomere/status2_exact_audit_20261028.4gQ5CB/` contains
`admission.jsonl`, `exact_summary.json`, `example_exact_geometry.json`,
`argv.json`, `replay.h5`, the logger-only translation unit and executable,
and the simulation log. Build/run driver: `../replay_exact_maturation.py`;
logger: `../exact_maturation_trace.h`.

On Sherlock continuation, read `/etc/agents/AGENTS.md` and relevant topic
files anew after compaction. Significant computation and I/O belong in Slurm.
