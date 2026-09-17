# Four-myosin ACTN2 partner-annealing test

Date: 2026-08-17

## Question

Can a force-balanced periodic cell with two myosins bound to every actin use
weak ACTN2 partner exchange to select one coherent transverse thick-filament
lane, then improve M-band registration without `k_mm`?

## Test geometry

`minimal_catch_balanced_double_myosin_anneal` contains four actins and four
myosins. Every actin initially has heads from two myosins. The two transverse
myosin neighbors on each axial side start 150 nm out of register. The cell is
periodic in x and inversion symmetric, so the mismatch is not a center-of-mass
kick. Each antiparallel actin pair has one explicit ACTN2 molecule.

Key controls:

```text
initial_structure=minimal_catch_balanced_double_myosin_anneal
n_actins=4, n_myosins=4, Lx=5.32 um, pbc_mask=100
dt=1e-6 s, nsteps=50000, save_every=250
k_aa=3000 pN/um, aa_stiffness_reference_links=1
k_mm=0
aa_status_mode=motor_transmitted_register
aa_motor_pair_support_window=0.005 s
aa_motor_pair_commit_time=0.010 s
aa_motor_pair_exchange_time=0.001 s
aa_motor_pair_dropout_time=0.005 s
aa_motor_pair_switch_ratio=1.5
aa_status2_force_on=1 pN
physical_mobility=true, temperature_kelvin=310, thermal_noise_scale=1
am_activation_floor=1, am_thin_activation=false, am_thick_feedback=false
seed=831704, deterministic, OMP_NUM_THREADS=1
```

The sustained activation is deliberate: this is a partner-competition test,
not a calcium-pulse test.

## Diagnostic iteration

The first event-instrumented run produced 1,155 switches in 50 ms. An
incumbent with zero instantaneous support made the ratio threshold zero, so any
transient challenger won and the commitment clock never advanced. Adding
dropout grace alone still produced 995 switches because positive but rapidly
fluctuating source forces alternated the strongest pair.

The final implementation requires the same qualifying challenger to remain
best continuously for 1 ms. It produced exactly four identity events:

```text
step 204:   link 2 selects myosin pair (2,3)
step 425:   link 3 selects myosin pair (2,3)
step 10204: link 2 locks pair (2,3)
step 10425: link 3 locks pair (2,3)
switches: 0
dropouts: 0
```

Both links reached status 2. Candidate support EMA at locking was 6.74 and
7.79 pN. The identity state machine is therefore stable and the maturation
gate is active.

## Physical result: failed annealing

The transverse M-band registration errors did not improve:

```text
pair       initial    final     median over final 10 ms
M0--M2     150.0 nm   198.1 nm  195.6 nm
M1--M3     150.0 nm   191.3 nm  190.5 nm
```

The selected `(2,3)` lane locks, but myosins 0 and 1 remain AM-bound to the same
actins. Candidate identity currently affects ACTN2 state/provenance only; it
does not penalize competing AM partners or create a mechanical coupling between
the selected thick filaments. Consequently identity locking alone cannot
produce M-band annealing.

This test also does not exercise source-budget oversubscription because there
is only one ACTN2 arm per actin. Budget conservation is covered by the
deterministic unit test; a separate multi-link test is required for its physical
validation.

## Remaining implementation gaps

1. **No selected-partner feedback to AM kinetics.** Heads from unselected
   myosins attach, power stroke, and detach exactly as before. A physical next
   step needs a weak-state, local rule that makes mechanically incompatible
   competitors easier to release; hard deletion of all nonselected heads would
   be unphysical and order dependent.
2. **No transverse lane compatibility gate.** Distinct myosin IDs are required,
   but the selector has no explicit material-azimuth/lane test. The present run
   chose a coherent lane, but that is not guaranteed generally.
3. **No global mature-partner capacity.** Separate ACTN2 molecules can lock
   mutually incompatible uses of the same myosin because locks are local.
4. **No M-band mechanical interaction.** `k_mm=0` is intentional here. If
   selected identities are supposed to register thick filaments mechanically,
   that coupling has not yet been implemented.

## Runtime and outputs

The final 50-ms run took 38.71 s wall time and 14.3 MB peak resident memory on
one CPU thread while saving 201 frames. Inside 50,000 local-link updates,
source construction used 1.058 s of worker time and candidate scans used
0.0167 s. HDF5/frame output dominates this small test; the identity scan is
negligible here.

Trajectory and rendered replay:

```text
/home/users/xsli/.codex/visualizations/2026/08/13/019ffcf1-524c-7e43-bf23-de3548be8297/motor-transmitted-register/multi-myosin-annealing/k3000_50ms_challenger_dwell.h5
/home/users/xsli/.codex/visualizations/2026/08/13/019ffcf1-524c-7e43-bf23-de3548be8297/motor-transmitted-register/multi-myosin-annealing/k3000_50ms_challenger_dwell.mp4
```
