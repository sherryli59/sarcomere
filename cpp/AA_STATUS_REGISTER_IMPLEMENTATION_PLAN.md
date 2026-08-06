# Minimal AA status/register implementation plan

## Model rules

1. Status 1 is a local ACTN2 connector at fixed material coordinates. It has
   radial, perpendicular, and nematic connector mechanics, but no barbed-end
   axial-register energy or force.
2. Status-1 association is extensive in available axial overlap:
   `k_pair = (free pair slots) * k_on * overlap / overlap_reference`.
   The reference length only fixes the units/calibration of `k_on`; the factor
   is not capped or exponentially saturated.
3. A 1 -> 2 transition requires opposed productive motor force above
   `aa_status2_force_on`, antiparallel actin polarity, and barbed-end registry
   error below `aa_status2_registry_on`. Per-link spring-load allocation does
   not enter this state decision. An established status-2 link persists until
   stochastic detachment or loss of `aa_status2_registry_off`; transient motor
   force loss does not instantaneously downgrade it.
4. Only status 2 receives an axial-register potential centered at
   `aa_antiparallel_overlap`. Each scalar register-error contribution uses

       U(e) = 0.5 k sigma^2 [1 - exp(-e^2 / (2 sigma^2))].

   It is quadratic near the target, has a finite maximum energy and force, and
   contributes both translational forces and orientational torques from the
   gradient of the same energy.

## Implementation sequence

1. Add unit-testable helpers for overlap-extensive association and registry
   hysteresis.
2. Extend the local-link mechanical result with registry diagnostics and
   replace the unbounded, force-only axial term with the bounded conservative
   potential.
3. Evaluate status classification from status-1 connector mechanics, apply the
   registry gate, then recompute final mechanics with the register enabled only
   for links classified as status 2.
4. Persist and resume-check the new registry parameters.
5. Add focused unit tests and the existing minimal force-balanced fixture.
6. Build, run tests and a short smoke trajectory, and fix inconsistent results.
