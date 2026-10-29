# Connector deletion cost and diagnostics

The 130 nm admission gate, buffered 120 nm deletion cylinder (10 nm inset
end planes), and full-connector tip contact remain always enabled. No new
model switches were added.

## Deletion cost

`_aa_status2_connector_excluded` now computes its segment bounding sphere
once per link and rejects distant myosin centers by minimum-image distance.
The exact periodic cylinder kernel also rejects by this sphere before image
enumeration and axis normalization. This is conservative: the closest image
of the cylinder center being outside the sphere rules out every image.
The exact same-point radial/axial test and non-nearest-image support remain.

Established mature-link geometry is evaluated in the existing parallel
mechanics pass, after site diffusion. A byte per slot records the result;
the existing barrier makes it available to status classification. Link
destruction, neighbor-slot release, RNG decisions and event insertion remain
in canonical slot order. A link that loses status 2 in the classification
step still uses its stored geometric result. A prospective new mature link
retains the direct check after admission.

## Detachment provenance

Reason codes 0–4 retain their original meanings. New code **5** means
`mature_connector_exclusion`. It has its own decision flag, participates in
the existing geometric-rupture cleanup, and records probability 1 without
drawing RNG. Code 3 is reserved for the legacy barbed-tip capsule rule.
The event row layout remains unchanged. New-file reason metadata and resume
checkpoints identify code 5 and the version-2 deletion rule.

Analysis readers recognize code 5 even when a resumed file's original
metadata lists only 0–4. Mature cohorts use maturation/link identity history,
rather than relying on status at death. Historical ambiguous reason-3 rows
cannot be separated retrospectively.

## Comparable overshoot

`analysis/visualizations/gaphunt-rig/status2_gate_20261028/connector_metrics.py`
recomputes the same fixed geometric metric on every trajectory:

- Enumerate all live connector / myosin / tip pairs and relevant periodic
  images, independently of the simulated force winner and force cutoff.
- Within a fixed 130 nm radial cylinder on a tip's physical half, report the
  greatest axial burial of any connector point. Radial and axial conditions
  must hold at the same point. Report overlaps deeper than 1 nm.
- Separate link-tip pair counts, unique links and unique myosin tips. Report
  all, current mature, inherited mature and post-resume mature cohorts.
- Retain a midpoint-only geometric reference using the same all-pairs
  population; never interpret recorded winner gap/rho as a comparable metric
  when the simulation's contact geometry changes.

The 120 nm buffered-region occupancy is also reported. Saved geometry is
after integration and can overlap briefly before the next deletion check.
The matched analysis now includes these readouts, while retaining recorded
selected-tip diagnostics with an explicit cross-version caveat.

Validation scripts, baseline reference kernel, and timing/analysis entry
points are in `analysis/visualizations/gaphunt-rig/status2_gate_20261028`.
Short timing continuations use disposable copies and explicit projection-off
arguments. They do not overwrite production trajectories or pinned binaries.

Final validation and measured results:
`analysis/visualizations/gaphunt-rig/status2_gate_20261028/results/performance_diagnostics.md`.
The active tree passed 44/44 tests and the shared repository passed 40/40.
The serial geometry pass is 2.59x faster; a paired four-CPU, 200-step resume
took 59.735 s before and 41.890 s after (29.9% less wall time).
