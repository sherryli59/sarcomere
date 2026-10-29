# Status-2 connector exclusion gate — 2026-10-28

## Agreed rule

Block an ACTN2 status 1 -> 2 transition if any point of its complete material
connector lies within 0.130 um of any myosin axis AND between that myosin's
two physical end planes. Both conditions must hold at the same connector point.
The radial boundary and end planes count as excluded. There is no additional
50 nm approach margin and no spherocylinder end cap.

The radius is fixed independently of the soft contact cutoff, force amplitude,
head occupancy and hard projection setting. Existing load/registry conditions
are still required for maturation.

## Implementation

- `include/actn2_maturation_geometry.h`: clip the connector to the axial slab,
  then minimize its radial distance on the clipped interval. Use ACTN2's shortest
  periodic connector and enumerate relevant periodic myosin images, including
  non-nearest images for long or tilted cylinders.
- `src/sarcomere_actn2.cpp`: apply the veto after status classification and
  status memory, before committing mature motor identity and arm state.
  A veto leaves the link at status 1 and clears a proposed status-memory lock.
  Every myosin is considered, with a conservative sphere prefilter, so force
  neighbor lists cannot silently omit a blocker.
- `src/sarcomere_io.cpp`: write the rule and radius in new-file metadata and
  every saved resume snapshot. Checkpoint attributes also record the active
  rule when resuming a trajectory whose original metadata predates this gate.

The original change affected admission only. The subsequent buffered deletion
rule below also enforces retention of mature links.

The changes were applied to both the durable repository at
`/home/groups/rotskoff/sherry/sarcomere/sarcomere/cpp` and the active source at
`$SCRATCH/sarcomere/zdisclink_20261020/cpp`, preserving their pre-existing
differences. Validation used the latter source tree.

## Validation

Added `tests/status2_connector_gate_test.cpp` covering:

- midpoint outside 130 nm while a connector end is inside;
- radial and axial boundaries, axial clipping, zero-length and crossing
  connectors, and simultaneous radial/axial eligibility at one point;
- rotated axes, periodic connectors and relevant non-nearest myosin images;
- rejection of the measured geometry from
  `gaphunt_20261023/data/r27zt0np_seed20260805.h5`, frame 475, step 237500,
  link 179352, actins 1316 and 233, against myosin 244;
- blocked maturation despite sufficient filtered productive support, successful
  maturation just outside the actual end plane, established mature retention,
  status-memory veto, and metadata/checkpoint provenance.

Final Release build and CTest run in Sherlock Slurm job 49947901:
**44/44 tests passed**, with assertions enabled for the regression executable.
Initial test-fixture issues (a nonperiodic connector that actually crossed the
cylinder and a registry reference inconsistent with the fixture's full actin
overlap) were corrected before this final run.

Build directory and logs:
`$SCRATCH/sarcomere/status2_gate_20261028.983jQo`
(`build.log`, `ctest.log`, `Testing/Temporary/LastTest.log`).

Validated executable:
`$SCRATCH/sarcomere/status2_gate_20261028.983jQo/sarcomere`

SHA256:
`dfa63b085df368246e93f8c27f121f4884b6f1e08fa0fc3f77379b65bb84c47b`

No simulation campaign was launched, and pinned campaign binaries were not
replaced. Production throughput of the new prospective-transition scan has
not yet been measured.

## Sherlock continuation reminder

Read `/etc/agents/AGENTS.md` at the start of the session and after context
compaction, then the topic files relevant to the next task. Compile and run
substantial analysis in a Slurm allocation; preserve inputs and model version
when preparing any follow-on campaign.

## Buffered mature-link deletion (subsequent change, 2026-10-28)

Mature ACTN2 links now detach immediately (reason 5) if any point on the full
connector is both within 120 nm of any myosin axis and between end planes
inset by 10 nm from the physical tips. This allows 10 nm penetration relative
to the formation cylinder (130 nm, actual end planes). Both conditions must
hold at the same connector point, including boundary equality. This is a
geometric tolerance, not an expansion of the deletion region.

The shared segment/cylinder kernel handles admission and deletion, including
periodic images and every myosin independently of force neighbor lists. The
deletion check uses old or current status 2, so a simultaneous loss of support
does not let an existing mature link evade deletion. It applies to inherited
mature links on their first update and to both arm-state models. It removes
the complete link through the existing geometric-detachment bookkeeping.

For mature links this replaces the barbed-tip capsule reason-3 hazard/dwell;
it does not depend on those legacy settings or the repulsion switch. Weak
links retain their previous burial rule. Contour rupture, collision rupture,
and stochastic detachment remain available. The force law is unchanged.

New-file metadata and every resume snapshot record the deletion radius,
axial inset, and immediate reason-5 rule. Original metadata in resumed files
remains historical; inspect the current resume group's rule for provenance.
Fixed constants are in `include/actn2_maturation_geometry.h`.

Validation: Release build and **44/44 CTest tests passed** in Slurm allocation
50007466, including full-connector midpoint misses, radial and axial buffer
boundaries, simultaneous point eligibility, periodic images, mature/inherited
link deletion, support loss, arm-state integration, reason-3 provenance, and
metadata/checkpoint fields. The weak-link physical-collision fixture explicitly
uses status 1 to keep testing reason 4 independently of the earlier mature-link
exclusion. Assertions are enabled in the regression executables.

Validated build: `$SCRATCH/sarcomere/connector_deletion_20261028`.
The test log is also saved under
`analysis/visualizations/gaphunt-rig/status2_gate_20261028/results/deletion_ctest.log`.
The pinned production executable and existing trajectories have not been
replaced. The original 0.4 s production comparison has not been repeated with
this deletion rule. The initial deletion implementation used reason 3; it is
now superseded by dedicated reason 5 and a faster geometry pass. See
`CONNECTOR_PERFORMANCE_DIAGNOSTICS_2026-10-28.md` for the current implementation
and short-continuation validation.
