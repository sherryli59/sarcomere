# Three-actin-per-half hexagonal lattice scaffolded control (superseded)

> This configuration used collective Z-centroid, ring-shape, and thick-lattice
> constraints. It is retained only as a positive-control record. The emergent
> local-titin implementation is documented in
> `EXPLICIT_TITIN_LOCAL_ASSEMBLY_TEST.md`.

Date: 2026-07-31

Scope: minimal sandbox only. No production run was launched.

## Topology

- 2 bipolar myosin filaments and 12 actin filaments.
- Each of the four myosin halves is surrounded by one three-actin population
  at 120-degree spacing.
- At each Z-disc, the two opposite-polarity triplets are offset by 60 degrees,
  making one six-actin ring.
- Each actin has two ACTN2 neighbor identities. Each ring has six unique edges
  and two turnover-capable links per edge (24 active links total).
- Four dynamic Z-ring--myosin-half links are recognized from aggregate A--M
  occupancy. Each stores the full six-actin connected component.

## Stabilization terms tested

The perimeter-only ACTN2 ring was insufficient: an equilateral hexagon has
internal shear modes. The stable diagnostic therefore uses four distinct
mechanical roles:

1. ACTN2 edge links set the edge length and axial A--A register.
2. A weak Z-ring centroid--thick-centerline term removes rigid transverse
   drift without fixing ring azimuth.
3. An internal thick-lattice backbone opposes relative transverse displacement
   and relative tilt of neighboring myosins without pinning the common frame.
4. A rotationally invariant Z-disc shape term keeps the six barbed ends at a
   common radius about their instantaneous centroid. It supplies the shear
   rigidity missing from perimeter springs but does not choose angular phase.

The axial titin term remains separate, tension-only, and distributed over the
whole Z-ring. Compressed titin generates no force.

For a regular hexagon, edge length equals circumradius. The clean diagnostic
therefore uses commensurate 55-nm AA and AM preferred spacings rather than the
frustrated 50.6/55-nm pair.

## Brownian validation

Command configuration:

- mechanism: `dynamic_titin_hex_fixed`
- initial cell length: 4.7 um
- initial A--A overlap: 160 nm
- seed: 1001
- 5-ms fixed-geometry pre-equilibration
- 50-ms three-dimensional Brownian trajectory with rotational diffusion
- force variant: `aa_k100x_register_all_l2`

Result file prefix:

`diagnostics/hex_three_membered_20260731/zshape_k500_full50ms`

Results:

- all 24 ACTN2 links remained active;
- all four dynamic Z--M links captured six-actin rings (three at 1 ms, the
  initially unoccupied fourth half at 4 ms);
- zero geometric ruptures, neighbor losses, or force/torque/rotation caps;
- final mean AM radius: 54.2 nm;
- final AM-radius CV: 5.0%;
- final triplet-chord CV: 8.7%;
- final maximum Z-ring angular-gap error: 12.2 degrees;
- over 10--50 ms: mean AM-radius CV 8.2%, mean triplet-chord CV 8.6%, and mean
  maximum angular-gap error 9.5 degrees;
- A--A overlap moved from 160 nm to 115.9 nm;
- mean A--M overlap moved from 495 nm to 516.8 nm.

## Still unresolved

- The transverse hexagonal structure is stable for this 50-ms Brownian test,
  but the A--A overlap has not yet reached the 100-nm target. A longer sandbox
  run is needed to distinguish slow continued annealing from a true stall near
  116 nm.
- The stiffnesses for Z-disc shear, thick-lattice registration, and transverse
  centering are diagnostic values, not experimentally calibrated parameters.
- Turnover is enabled but deliberately slow (0.5--1 s link lifetimes); a longer
  run is required to demonstrate repeated edge loss and recovery while the
  hexagon remains stable.
- The 12-actin topology and stabilizers have not been promoted to the
  hundreds-of-filaments production assembly.
