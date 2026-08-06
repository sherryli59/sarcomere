# Bilateral explicit-titin perturbation gate (2026-08-05)

## Question

Can local, tension-only explicit titin reverse a deliberately buried thick
filament without introducing a compression force or a pre-existing centroid/Z
node?

## Correct initializer

1. Assemble the centered 2-myosin/4-actin cell.
2. During fixed-geometry pre-equilibration, capture one explicit titin molecule
   on each actin barbed end (one partner on each thick-filament half).
3. Preserve those four material bonds.
4. Displace myosin 0 by -0.55 um and myosin 1 by +0.55 um.
5. Release mechanical motion and chemistry together.

This ordering is essential. If the thick filaments are displaced before
capture, the far-side molecule is outside the local capture range and cannot
provide the restoring tension.

## Local observable

For each myosin, project the two sets of titin-bound actin barbed ends onto the
instantaneous myosin axis. The translation-invariant centering error is

    e_m = abs(0.5 * (mean(s_minus) + mean(s_plus))).

The reported error is `max_m e_m`. It is unchanged by rigid translation of the
whole local assembly; the previous lab-frame myosin-center metric was not.

## Results

- Immediately after the perturbation: 2 titins stretched, 2 titins slack,
  4/4 actins titin-bound, and maximum local midpoint error = 0.550 um.
- Fresh 50 ms thermal run: error 0.550 -> 0.218 um (39.6% remaining),
  4/4 bonds retained, AA overlap 0.102 um, mean AM overlap 0.371 um: PASS.
- Previously saved 200 ms deterministic trajectory, recomputed with the local
  metric: 0.550 -> 0.087 um.
- Previously saved 200 ms thermal trajectory, recomputed with the local metric:
  0.550 -> 0.180 um.
- Free-actin 100 ms control still passes directed translocation.
- Full CTest suite: 28/28 passed.

The equilibrium is a motor-titin force balance, not necessarily zero geometric
error. The gate therefore requires more than 55% reversal and a final local
error below 0.25 um, while retaining bilateral ownership, AA register, and a
finite AM overlap.

## Mechanical regression

The explicit-titin restoring test now checks the buried configuration directly:
after a 0.55 um shift, the near-side titin has exactly zero force, the opposite
titin bears tension, the net molecular force is reciprocal, and all original
owners remain bound. No short-length penalty or compression-bearing titin force
is used.

## Files

- Fresh log: `local_midpoint_gate/seed20260825_50ms_gate.log`
- Fresh CSV: `local_midpoint_gate/seed20260825_50ms_gate.csv`
- Free control: `controls/free_100ms.log`
- Saved trajectories used for the 200 ms recomputation:
  `k100_t0_200ms.h5`, `k100_t1_200ms.h5`

