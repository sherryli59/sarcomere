# Shared thick/thin lattice and ACTN2 neighbor identity

## Geometry represented now

`--initial_structure=shared_thick_thin_lattice` is a periodic,
force-balanced A-band diagnostic. It contains four thick filaments and eight
shared thin-filament positions on each axial half (16 actins total). The
orthorhombic transverse unit cell is

\[
L_y=2\sqrt{3}\,d_{\rm AM},\qquad L_z=3d_{\rm AM}.
\]

Under transverse periodic boundaries:

- every thick filament has six thin-filament neighbors;
- every thin filament is shared by three thick filaments;
- the two axial halves have opposite actin polarity;
- the two myosin halves have equal transverse coordination, so no particle is
  fixed to obtain axial force balance;
- opposing pointed ends have a one-actin-diameter M-line gap and are not
  incorrectly treated as a Z-disc overlap.

The initializer can prebind the 3-regular same-polarity AA graph on each half
and seed six weak AM heads on each productive thick/thin pair.

## AA topology versus molecular occupancy

`aa_neighbor_identity::Registry` separates a retained filament-neighbor slot
from the zero, one, or two explicit ACTN2 molecules currently occupying it.

- A status-1 slot survives a complete molecular vacancy for
  `aa_neighbor_slot_memory` (default 20 ms), allowing ordinary turnover.
- A slot observed in status 2 survives for `aa_status2_slot_memory` (default
  200 ms), so a catch-supported lattice edge has priority without deleting a
  competing physical bond.
- Pinned diagnostic lattice edges do not anneal away.
- `max_aa_bonds_per_actin` now caps retained distinct neighbors (default four),
  not the number of explicit ACTN2 molecules.
- `aa_prebundle_links_per_pair` separately caps molecular redundancy per
  retained edge (default two).
- New neighbors wait for an available slot; the model never forcibly ruptures
  an occupied ACTN2 link to admit another candidate.

This gives status 2 retention priority while permitting status-1 topology to
anneal after its short vacancy memory.

## Steric correction

Actin hard-core forces are now always equal and opposite, independent of
catch-bond status. Exactly coincident parallel/antiparallel rods use a
deterministic pair-symmetric transverse regularization rather than silently
receiving zero repulsion. The numerical hard-core penalty is
`10000 pN/um`; it is an excluded-volume penalty, not an ACTN2 stiffness.

## Scope boundary: Z-disc

The current diagnostic tests the shared A-band hexagonal lattice. It does not
pretend that the same transverse graph is the cardiac Z-disc. A subsequent
Z-disc diagnostic should introduce the experimentally observed staggered
small-square/basketweave opposite-polarity lattice at the axial boundary, with
its own neighbor sectors and spacing calibration. Reusing the old six-private-
actin ring for this purpose creates coaxial overlaps and conflates two
different lattices.

## Acceptance diagnostics

`analysis/analyze_shared_lattice.py` reports:

- AM centerline spacing and dispersion;
- local sixfold order, actin nematic order, and myosin nematic order;
- weak/strong AM head populations;
- molecular status-1/status-2 ACTN2 populations;
- occupancy of the intended AA graph and unexpected edges;
- minimum finite-segment AA distance, soft-core contacts, and severe overlaps;
- force, torque, and rotation cap events.

`analysis/animate_shared_lattice.py` renders the 3D cell, both transverse
halves, live AA edges, and synchronized AM/ACTN2 KMC counts.
