# Identity-resolved motor-transmitted ACTN2 maturation

## Purpose

`--aa_status_mode=motor_transmitted_register` is an opt-in alternative to
`productive_register`. It prevents status-2 ACTN2 maturation from being
triggered by an aggregate motor vector or by passive ACTN2 prestress alone.
The default remains `productive_register`, so existing runs do not change.

This mode is available only with `dynamic_am=true` and
`aa_local_multilink=true`.

## Implemented rule

For every actin, productive strong-head forces are accumulated separately by
`(myosin_id, half)`. For a local ACTN2 link between actins `i` and `j`, let
`C_i` be the force the link applies to actin `i`; by construction
`C_j = -C_i`. For a motor contribution `M` on the same actin, the supported
part of that connector force is

```text
S(M, C) = max(0, -M dot C / |C|).
```

For a candidate source pair `(m_i, m_j)`, the evidence is

```text
L_pair = min(favorable ACTN2 arm i,
             favorable ACTN2 arm j,
             S(M_i, C_i),
             S(M_j, C_j)).
```

The two source myosins must be distinct. Thus a strained AA spring, an AA-AA
network force, or a steric collision cannot create maturation evidence unless
two explicit productive motor sources also balance this particular connector
on its two actins. Steric and other non-ACTN2 forces are never added to
`C_i`.

Each source vector is also a finite budget. If the same `(myosin_id, half)`
opposes several ACTN2 arms, its raw projected demands are summed first and the
source magnitude is distributed among those arms in proportion to demand.
The full source force therefore cannot be credited independently to every
incident link.

### Weak-state partner exchange

While the ACTN2 molecule is status 1:

1. The strongest source pair is selected with deterministic identity
   tie-breaking.
2. An alternative qualifies only when its instantaneous support exceeds
   `aa_motor_pair_switch_ratio` times the current pair's support.
3. The same qualifying alternative must remain the strongest challenger
   continuously for `aa_motor_pair_exchange_time`; a momentary force spike is
   not an exchange.
4. A zero-support gap cannot be exploited by a challenger. The incumbent is
   retained through `aa_motor_pair_dropout_time`, then dropped and required to
   earn a fresh commitment interval if it returns.
5. The selected pair's support is exponentially averaged over
   `aa_motor_pair_support_window`.
6. Every identity switch or dropout resets the commitment clock and support
   average.
7. Maturation requires the same selected identities to reside for
   `aa_motor_pair_commit_time`, the averaged support to exceed
   `aa_status2_force_on`, and the narrow registry gate to pass.

This is a coarse-grained annealing rule, not a claim that ACTN2 directly reads
a myosin filament ID. The support and registry tests are mechanical; the
identity residence clock represents unresolved repeated head turnover from a
persistent thick-filament partner. Exchange is hysteretic/deterministic rather
than a separately calibrated chemical KMC step.

### Mature-state identity lock

On the status-1 to status-2 transition, the selected myosin IDs and halves are
locked while the mature conformation persists. Head turnover cannot relabel
the mature link, but its support is still filtered so sustained unloading can
be detected. With a positive `aa_status2_unlock_rate`, a mature link can relax
to status 1 without ACTN2 detachment; that releases the motor identities for
weak-state partner exchange. A zero unlock rate preserves the historical lock.
Its catch-bond lifetime uses the current favorable bilateral force carried by
the ACTN2 connector, not the historical maturation evidence.

## Difference from `productive_register`

| Behavior | `productive_register` | `motor_transmitted_register` |
|---|---|---|
| Motor accounting | One aggregate productive vector per actin | Separate vector per myosin ID and half |
| Source identity | Not retained | One distinct myosin pair selected per ACTN2 molecule |
| Transmission test | Opposed aggregate vectors | Each motor source must oppose this link's connector reaction |
| Multiple myosins on one cross-linked actin | Contributions can cancel or combine before classification | Candidate combinations are evaluated explicitly |
| Passive AA prestress | Can coincide with aggregate motor evidence without source attribution | Cannot mature the link without bilateral motor-to-connector support |
| Weak-state annealing | Aggregate support EMA only | Hysteretic partner exchange plus pair-specific EMA and dwell time |
| Mature identity | None | IDs and halves frozen until mature-state relaxation or molecular detachment |

## Controls

- `--aa_motor_pair_support_window=0.010`: seconds used by the selected-pair
  support EMA. Zero uses instantaneous support.
- `--aa_motor_pair_commit_time=0.010`: minimum continuous candidate residence
  before maturation. Zero disables the dwell-time requirement.
- `--aa_motor_pair_switch_ratio=1.5`: alternative/current support ratio needed
  to switch candidates; it must be at least one.
- `--aa_motor_pair_exchange_time=0.001`: the same stronger challenger must
  remain best continuously for this many seconds before exchange.
- `--aa_motor_pair_dropout_time=0.010`: zero-support grace interval for the
  incumbent identity. At timeout the identity and residence clock are cleared.
- Existing `aa_status2_force_on`, `aa_status2_registry_on`, and
  `aa_status2_registry_off` retain their usual units and meanings.
- The optional `aa_status2_unlock_rate`,
  `aa_status2_unlock_force_scale`, and
  `aa_status2_unlock_registry_scale` controls are documented in
  `REVERSIBLE_ACTN2_MATURATION.md`.

The 10 ms defaults are initial coarse-grained choices, not fitted constants.
They should be swept against registration quality, component-size growth, and
head residence statistics.

## Diagnostics and restart state

Trajectory metadata records all three controls and the force-source semantics.
`/catch_bond/local_link_motor_support` is row-aligned with
`/catch_bond/local_link_identity` and contains:

```text
candidate_myosin_i, candidate_myosin_j,
candidate_half_i, candidate_half_j,
support_ema, support_instantaneous, candidate_age,
identity_locked, candidate_since_step, candidate_last_seen_step,
raw_source_support_i, raw_source_support_j,
allocated_source_support_i, allocated_source_support_j,
challenger_myosin_i, challenger_myosin_j,
challenger_half_i, challenger_half_j,
challenger_support_ema, challenger_age
```

`/catch_bond/motor_pair_events` records exact selections, completed exchanges,
dropouts, locks, and ACTN2 releases. Candidate/challenger state and pending
events are preserved in `/resume`, so residence and exchange clocks survive a
restart exactly.

The same identity, timing, lock, and evidence fields are written into the
atomic `/resume` snapshot. A run in this mode refuses a resume snapshot that
lacks them.

At normal process exit the executable prints cumulative profiling fields:

- local ACTN2 update wall time;
- source-building worker-thread time;
- candidate-scan worker-thread time;
- number of distinct candidate combinations evaluated.

Worker-thread times can exceed wall time in multithreaded runs and should not
be added to wall time.

## Tests

`actn2_load_kinetics_test` covers mode parsing, connector opposition signs,
bilateral limiting, finite source-budget conservation, rejection of
same-myosin and missing sources, EMA behavior, switch hysteresis, continuous
challenger dwell, dropout grace, commitment gating, identity lock, registry
formation/retention hysteresis, and mature-state persistence after transient
motor loss.

## Current mechanistic boundary

Identity selection is presently an ACTN2 maturation state and provenance
record. It does **not** alter AM attachment/detachment rates, reserve a myosin
for one lane globally, or apply a force directly between selected thick
filaments. The four-myosin test in `MULTI_MYOSIN_ANNEALING_TEST.md` demonstrates
the consequence: a coherent pair can lock while competing myosins remain bound
to the same actins and drift farther out of M-band register. This is a known
implementation gap, not a successful annealing result.

Validation performed on 2026-08-17:

```text
cmake --build cpp/build -j4
ctest --test-dir cpp/build --output-on-failure -j4
28/28 tests passed
```

A 5-repeat, single-thread timing comparison used the two-sarcomere unit with
4 actins, 2 myosins, 1,000 fixed-geometry chemistry steps, and 10,000 dynamic
steps per repeat. The only changed option was the status mode:

```text
productive_register:          mean wall 2.490 s (range 2.41-2.54 s)
motor_transmitted_register:   mean wall 2.400 s (range 2.33-2.44 s)
observed wall difference:     -3.6% (no resolved slowdown)
```

In the last new-mode repeat, 11,811 candidate combinations were evaluated.
The internal counters reported 0.0663 s of source-building thread time and
0.00121 s of candidate-scan thread time across 11,000 profiled updates. This
small benchmark is a regression check, not a large-system scaling result. The
candidate scan scales approximately with the product of source counts on the
two actins incident on each active ACTN2 link.
