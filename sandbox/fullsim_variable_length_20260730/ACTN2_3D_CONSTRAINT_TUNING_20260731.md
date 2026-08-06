# ACTN2 3D constraint tuning (2026-07-31)

Scope: the 2-myosin/12-actin `explicit_titin_hex_fixed` sandbox diagnostic only. No production run was launched.

## Failure found

The old baseline used a per-link distance stiffness of only `10 / 3 = 3.33 pN/um`. Opposing-triplet ACTN2 links reached 463--583 nm despite a 55 nm preferred length and persisted because their lifetime was at least 1 s. The fixed-index A--A diagnostic also became invalid after dynamic neighbor exchange.

## Tuned model

- `actn2_tuned` is now the default force variant for `explicit_titin_hex_fixed`.
- Per-link distance stiffness: `1000 / 3 = 333.3 pN/um`.
- Axial register ratio: 3; saturating force ceiling: 36 pN.
- New-association polarity gate: director dot product must be below -0.5.
- Finite ACTN2 contour rupture: 0.12 um. A ruptured link releases its remembered neighbor slot immediately so the network can reanneal.
- Torque limiting now acts on the director-tangent component that actually rotates actin, rather than on the irrelevant director-parallel gradient.
- CSV output now includes `aa_linked_overlap`, the overlap of the dynamically connected antiparallel neighbor pairs. The legacy `aa_overlap` field is retained for fixed-index comparisons.

These are local molecular rules. They do not introduce a Z-node, centroid force, prescribed ring, or pinned neighbor list.

## Four-seed, 200 ms, full-diffusion validation

| Seed | Linked A--A overlap | A--M overlap | Mean/max actin tilt | Maximum saved ACTN2 length | Cap frames |
|---:|---:|---:|---:|---:|---:|
| 1001 | 90.7 nm | 566.6 nm | 4.4 / 5.6 deg | 94.3 nm | 0 |
| 1002 | 87.0 nm | 600.9 nm | 4.2 / 6.5 deg | 95.0 nm | 0 |
| 1003 | 85.9 nm | 567.5 nm | 6.0 / 8.2 deg | 95.2 nm | 0 |
| 1004 | 84.3 nm | 574.2 nm | 5.2 / 6.4 deg | 97.1 nm | 0 |

All links remained below the 120 nm contour limit; the 463--583 nm crosslinks and large fan-out did not recur.

## Remaining limitation

The catastrophic 3D crosslink-distance failure is fixed. Complete sixfold transverse order is not yet robust (`z_radius_cv` remains seed-dependent), so these runs validate a stable local ACTN2 constraint and finite A--M overlap, not a fully annealed hexagonal lattice.
