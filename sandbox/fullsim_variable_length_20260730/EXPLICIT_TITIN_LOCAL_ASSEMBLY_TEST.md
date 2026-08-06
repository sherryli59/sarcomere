# Explicit local titin: minimal assembly test

Date: 2026-07-31

Scope: sandbox only. No production run was launched.

## Mechanism

- Removed occupancy-triggered dynamic titin.
- Removed ACTN2-component and Z-node discovery from titin mechanics.
- Removed Z-centroid centering, Z-ring radial shape, and two-myosin backbone
  constraints from the axial-cell sandbox.
- Added independent coarse-grained titin molecules. Each has one material
  anchor on a myosin half and can bind at most one individual actin barbed end.
- Capture uses only local half/polarity geometry, a finite capture distance,
  and a stochastic on-rate. It does not query myosin-head occupancy or ACTN2.
- A bound titin applies a reciprocal 3D tension-only endpoint force and the
  corresponding actin and myosin orientation drives.
- Titin turnover uses its own optional force-dependent off-rate.

The 12-actin diagnostic uses three titins per myosin half. It starts from a
near-contact, jittered arrangement with no ACTN2 bonds and no seeded myosin
heads. The population and filament polarities are supplied, but no Z-disc,
neighbor graph, ring radius, or ring azimuth is supplied.

## Verification

All 24 sandbox tests pass. Focused assertions show that:

- titin captures individual actins when every ACTN2 link is inactive and every
  myosin head is detached;
- one titin pulls only its own actin rather than sharing force across a group;
- compressed titin generates exactly zero force;
- stretched titin gives equal-and-opposite actin/myosin forces;
- off-center endpoint forces generate orientation drive;
- force-dependent or force-independent turnover remains local to each titin.

## Full-diffusion smoke run

Configuration:

- mechanism: `explicit_titin_hex_fixed`
- 2 myosins, 12 actins, 3 titins per myosin half
- full translational and rotational Brownian noise
- seed 1001
- 50 ms
- no preequilibration, no prebound ACTN2, no seeded heads
- baseline local ACTN2 mechanics

Output prefix:

`diagnostics/explicit_titin_local_20260731/seed1001_50ms`

Observed:

- all 12 titins captured individual actins within the first millisecond and
  remained bound through 50 ms;
- zero force-cap, torque-cap, or rotation-cap events;
- 16 ACTN2 links across 10 distinct actin pairs were active at 50 ms;
- myosin heads attached and cycled dynamically (274 attachments, 100 power
  strokes, 2 geometric ruptures);
- A-A overlap changed from 164 nm to 173 nm after a fast transient;
- A-M overlap changed from 492 nm to 574 nm;
- transverse hexagonal order did not emerge: final A-M radius CV was 67%,
  Z-radius CV 53%, Z-edge CV 117%, and maximum angular-gap error 115 degrees.

## Interpretation

This is a successful mechanics/causality test and a negative assembly result.
Local titin now works before any Z-disc-like component exists, but explicit
titin plus the current ACTN2 and motor models do not self-organize the
transverse hexagonal lattice on this timescale. The prior apparent stability
therefore depended on the removed global scaffold terms.

The initializer is still a near-contact minimal assay, not a de novo assembly
from components randomized throughout the full box. That broader test should
wait until local ACTN2 bending/binding-domain mechanics and local M-line
crosslinking are specified.
