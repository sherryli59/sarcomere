# Full-connector myosin-tip repulsion detection

The tip contact previously tested only the ACTN2 connector midpoint. A bridge
could therefore enter the 130 nm radial support with an endpoint or interior
section while requesting zero repulsion. The contact kernel now tests the
entire material connector, using its shortest periodic image.

For each tip, clip the connector to that tip's own myosin half and the axial
onset plane (the existing 50 nm approach margin). Minimize radial distance on
that interval, then evaluate the existing axial force law at the same point.
If the minimum is on the zero-force onset plane while other points qualify,
use an interior point of the simultaneous radial/axial intersection. Equal
radial distances prefer greater penetration. Relevant periodic myosin images
are enumerated, including non-nearest images for tilted or long filaments.

The preliminary sphere test now includes half the connector length. The
force scan checks every myosin with that conservative test, rather than
depending on the union of endpoint actin neighbor lists. Existing canonical
link order and strongest-connector-per-tip selection are preserved.

The 100–130 nm smoothstep, 20 pN maximum, 50 nm axial onset, axial-only force
direction, equal force split between the actins, opposite myosin reaction,
orientation drives and total-force caps retain their existing meanings. The
smoothstep uses an algebraically equivalent expression that stays positive
just inside its outer cutoff. Formation and buffered mature-link deletion
retain their separate, already implemented rules.

The lightweight gap query shares this geometry with the force kernel. Thus
hard projection, when enabled, also queries full-connector contacts; it is not
enabled by this change. Association/diffusion accessibility queries calling
the same contact kernel also use its new geometry. Recorded tip gap and rho
now refer to the selected connector contact point rather than its midpoint.
New-file and resume metadata record
`aa_zdisc_tip_contact_geometry=full_connector_min_radius_same_point_axial_v1`.

Changes are applied to both durable `cpp` and the active source at
`$SCRATCH/sarcomere/zdisclink_20261020/cpp`, preserving prior differences.
The existing production executable and trajectories are retained.

Regression coverage includes midpoint misses, endpoint and interior contact,
same-point eligibility, opposite halves, onset/cutoff boundaries, degenerate
connectors, rotated axes, periodic translation and non-nearest images,
conservative sphere coverage, the measured link 179352 / myosin 244 pose,
3,000 connectors checked against independent dense point sampling, and force
application with both endpoint actins' myosin neighbor lists empty.

Validation: **44/44 CTest tests passed**, assertions enabled, in Slurm
allocation 50009605. The gap/force parity test checks 300,000 randomized
cases; the full-connector dense-sampling test checks another 3,000 segments.
The first neighbor-coverage fixture used endpoints still inside the broad
neighbor sphere; increasing its bridge length isolated the intended omission
without changing production geometry or disabling production checks.

Validated Release build: `$SCRATCH/sarcomere/connector_force_20261028`.
Executable SHA256:
`c23ae96b5e6024983ee333b6613dced1f34143a669ec51a1e3f28cf39596685d`.
The durable test log is
`analysis/visualizations/gaphunt-rig/status2_gate_20261028/results/full_connector_force_ctest.log`.
No production continuation has been run with this contact geometry yet.
