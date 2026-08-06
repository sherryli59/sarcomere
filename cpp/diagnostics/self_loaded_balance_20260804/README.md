# Endogenous ACTN2 load validation (2026-08-04)

This diagnostic removes the external force probe and disables titin.  It uses
the same motor and ACTN2 parameter set in three minimal systems.

## Coordinates and forces

For the balanced cell, the force-bearing sliding coordinate is

`q_ext(t) = overlap_AA(0) - overlap_AA(t)`.

It increases while opposed motors pull the two antiparallel actins out of their
initial 160 nm interdigitation.  The ACTN2 load is computed from the actual 3D
connector forces, projected onto each actin axis:

`L_AA = max(0, -sum_i F_AA,i dot u_i / N_actin)`.

No applied force, overlap proxy, titin force, centroid force, or prescribed
velocity enters this value.  Two-millisecond windows provide paired
`(L_AA, dq_ext/dt)` force-velocity samples.

Status-1 ACTN2 binding sites diffuse along the 200 nm barbed-end association
zone.  A status-2 molecule locks its two material sites until stochastic
dissociation or contour rupture.  Its lifetime remains load dependent.

## Three-scenario result, seed 2080206

The 100 ms/300 ms/100 ms suite in `suite/` uses 35 heads per myosin half and
full translational and rotational Brownian noise.

- Free actin: 0.322 um directed travel in 100 ms; terminal velocity
  3.30 um/s.
- One-sided crosslinked pair: the motor-bound actin traveled 0.541 um and the
  unbound actin followed by 0.530 um.  Nine unsupported ACTN2 molecules turned
  over, with 43 ms observed mean dwell; none reached status 2.
- Balanced cell: A-A overlap annealed from 0.160 to 0.099 um.  The unlocked
  extension velocity was 6.03 um/s and the locked velocity was 0.018 um/s.
  The mean late ACTN2 load was 1.09 pN/actin and the deterministic force
  residual was 0.14%.

Thus the motor, unsupported-turnover, endogenous-load, and load-bearing stall
mechanisms all work.  The old external-load diagnostic is no longer needed for
this check.

## Remaining topology limitation

The same run also falsifies the stronger claim that ACTN2 alone can maintain
every individual A-M overlap indefinitely.  A status-2 ACTN2 link has energy
`U_AA(r_i - r_j, u_i, u_j)`.  Translating both linked actins by the same vector
leaves `U_AA` unchanged, so ACTN2 supplies exactly zero restoring force for the
pair centroid.  In the 100 ms 35-head run, one owner-specific A-M overlap falls
to zero even though A-A register remains stalled.  At 70 heads/half all four
overlaps survive 100 ms, but one falls to 30 nm by 300 ms.  At 140 heads/half a
different actin reaches zero by 300 ms.  Increasing molecular number therefore
reduces short-time noise but does not remove the zero mode.

A stable finite A-M overlap for every actin needs one additional physical
ingredient that couples a Z-complex centroid to the thick-filament lattice
(for example explicit tension-only titin), or an explicit unit-cell symmetry
constraint.  The latter is a boundary condition, not emergent assembly.
Neither was silently added here.

## Reproduction

Run:

`cpp/run_motor_validation_suite.sh`

The script returns failure if any scenario-specific gate fails.  Its optional
sixth argument sets the balanced duration; using 300 ms activates the stricter
long-time exposure of the centroid zero mode.
