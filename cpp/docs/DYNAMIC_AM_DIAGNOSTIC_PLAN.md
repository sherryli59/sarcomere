# Dynamic AM diagnostic and validation plan

## Goal

Determine whether actin crosslinks stabilize dynamic actin-myosin attachments
through retained mechanical strain, whether instantaneous status-2 classification
flickers too rapidly, and which particle-level variables are needed in a minimal
continuum model.

## Phase 1: event-level instrumentation

Record every AM transition with

```text
actin_id, myosin_id, step, old_state, new_state, transition_reason,
activation, axial_strain, motor_force, preceding_actin_cb_status
```

Track AA bond presence separately from its status, including continuous bond
episodes, status-1/status-2 episodes, motor force on both actins, and strong-motor
counts. Produce AM dwell distributions, status-2 interruption distributions,
AA survival curves, and event-triggered trajectories.

## Phase 2: isolated motor verification

Use one fixed AM pair under tonic activation. Verify

```text
mean weak dwell = 1 / (k_weak_off + k_powerstroke)
P(weak -> strong) = k_powerstroke / (k_weak_off + k_powerstroke)
mean zero-load strong dwell = 1 / k_strong_off
```

Require measured values to agree with the analytical three-state model within
5%. Repeat under prescribed loads `F/F* = 0, 0.5, 1, 2` and verify the intended
linear dependence of log transition rate on force.

## Phase 3: minimal crosslink assay

Use two antiparallel actins and fixed flanking myosins. Run at least 20 short
replicates for:

1. Free actins without an AA crosslink.
2. Actins joined by a persistent AA crosslink.
3. Fixed actins without an AA crosslink.
4. Fixed actins joined by a persistent AA crosslink.

Compare strong-AM residence, weak-to-strong probability, strain relaxation,
motor force, AA survival, and status-2 duty fraction. The difference between
conditions 1 and 2 contains the crosslink's mechanical effect. Conditions 3 and
4 must have indistinguishable AM kinetics because fixing the actins removes the
crosslink's mechanical effect; any residual difference would reveal a hidden
direct status-to-AM coupling.

## Phase 4: status-2 flicker audit

For each continuous AA bond, calculate status-2 episode duration, embedded
status-1 gap duration, the number of `2 -> 1 -> 2` interruptions, survival through
each gap, and delay from strong-AM loss to AA rupture. Compare

```text
P(status-1 gap > tau1)
```

across parameters. Keep the instantaneous rule if gaps are much shorter than
`tau1`. Otherwise replace it with actual-force classification plus hysteresis or
a fitted load-relaxation timescale. AA bond presence must remain distinct from
its force-bearing status.

## Phase 5: reduced-network causal controls

Before full-size sweeps, compare three replicates of:

| AM model | activation | AA lifetime law |
|---|---|---|
| legacy geometric | tonic | status-2 catch |
| dynamic, force-dependent | tonic | status-2 catch |
| dynamic, force-independent | tonic | status-2 catch |
| dynamic, matched constant dwell | tonic | status-2 catch |
| dynamic, force-dependent | pulsatile | status-2 catch |
| dynamic, force-dependent | matched-mean tonic | status-2 catch |
| dynamic, force-dependent | pulsatile | matched constant AA lifetime |

Match mean AM occupancy, strong duty ratio, AM residence, and AA occupancy before
comparing `phi_3d`, assembly time, force distributions, turnover, and conditional
residence times.

## Phase 6: causal timing

Compute event-triggered averages around weak-to-strong transitions, status
`1 -> 2` and `2 -> 1` transitions, and AA attachment/rupture. Test the sequence

```text
power stroke -> AM force -> loaded AA bond -> slower strain relaxation
             -> longer strong-AM residence
```

Measure lagged correlations among motor force, status-2 fraction, and `phi_3d`.

## Phase 7: numerical robustness

Repeat representative tests at `dt`, `dt/2`, and `dt/4`. Require convergence of
transition rates, dwell distributions, force distributions, status flicker, AA
survival, and `phi_3d`. Confirm exact uninterrupted versus resumed trajectories.

## Phase 8: focused full-size runs

At `N_actin=2100`, `N_myosin=1000`, run three replicates initially for legacy,
dynamic load-dependent tonic, dynamic load-independent tonic, matched-dwell,
dynamic pulsatile, and dynamic pulsatile with matched constant AA lifetime.
Increase to five replicates near qualitative transitions.

## Phase 9: continuum reduction

Coarse-grain detached/weak/strong motor fractions, crosslink density,
force-bearing crosslink fraction, active stress, orientational order, and
`phi_3d`. Fit the smallest closed model capable of reproducing measured response
and memory times. Add a separate load-bearing crosslink variable only if the
flicker and event-triggered tests show that instantaneous stress is insufficient.

## Immediate implementation order

1. AM transition event table and checkpoint-safe buffering.
2. AM dwell/hazard summarizer.
3. Continuous AA bond/status episode table.
4. Single-pair and force-clamp drivers.
5. Minimal two-actin crosslink driver.
6. Reduced-network manifest and aggregation scripts.

## Implementation status

- Complete: checkpoint-safe AM transition table and dwell/cause analysis.
- Complete: checkpoint-safe AA bond/status transition table and flicker analysis.
- Complete: fixed single-motor force-clamp initializer, array, aggregation, and watchdog.
- Complete: four-arm minimal-crosslink array, aggregation, and watchdog.
- Complete: 20-replicate force clamp. Measured/expected strong off-rate ratios
  are 0.990, 1.057, 0.933, and 1.007 at `F/F*=0,0.5,1,2`.
- Complete: 20-replicate crosslink assay. Fixed no-AA and persistent-AA motor
  off-rates are statistically indistinguishable (difference 1.40/s, approximate
  95% CI -1.11 to 3.92), ruling out a direct status-to-AM rate multiplier.
- Complete: flicker audit. For fixed persistent AA pairs, the median internal
  status-1 gap is 3.63 ms and 84.8% exceed `tau1=1 ms`; instantaneous status is
  therefore too brittle for short-lived status-1 bonds.
- Complete: opt-in, checkpoint-safe `am_status2_memory`; default zero preserves
  legacy behavior and strong AM support remains necessary to establish status 2.
- Complete: 20-replicate status-memory sweep at 0, 1, 5, and 20 ms with
  `tau1=1 ms`. In fixed geometry, AA occupancy rises from 0.818 to 0.852,
  0.930, and 0.979; status-2 duty rises from 0.791 to 0.838, 0.933, and 0.988,
  while the AM off-rate remains 76--77/s. This confirms that memory acts only
  through AA classification and lifetime, not a hidden AM-rate multiplier.
- Decision: use 5 ms as the provisional reduced-network memory. It is close to
  the measured 3.63 ms median interruption and recovers most status-2 occupancy
  without the near-locking behavior of 20 ms. Keep 0, 1, and 20 ms as controls.
- Next: construct the reduced-network force-dependence, activation-history, and
  matched-dwell controls before launching full-size runs.
