# Collapsing ACTN2 to one link per actin pair — implementation plan

Written 2026-08-06. Supersedes the multi-link arrangement described in
`MULTILINK_ACTN2_IMPLEMENTATION_PLAN.md` for the production path.

## Why this is safe

Traced 2026-08-06 against `aa_status_mode=productive_register` (what production
runs). Nothing that gates status or lifetime depends on link count:

| quantity | count-dependent? | at N=1 |
|---|---|---|
| `effective_load` (drives lifetime + status gates) | **no** — equals `productive_support` in this mode | unchanged |
| `productive_support` | **no** — function of the two filament directions and cached motor force | unchanged |
| `weight = link_tension / pair_tension_sum` | share | self-normalises to 1.0 |
| `allocated_motor_*` | share | self-normalises; branch only active in `arm_load_productive` mode |
| `pair_tension_sum` | sum over pair | equals the single link's tension |
| `occupied` -> association rate | **yes** | see k_on decision below |
| `aa_local_link_count`, `aa_local_status1_count` | **yes** | reporting only, drops ~3x |

So the status-2 thresholds (`aa_status2_force_on=0.25`) compare against a
count-free quantity and carry over unchanged.

A single link is also not mechanically impoverished. `compute_local_aa_link`
already carries three independent constraints:

```
energy = 0.5*radial_k*extension^2                       // separation
       + 0.5*perpendicular_k*(proj_i^2 + proj_j^2)      // resists axial sliding
       + 0.5*alignment_k*(1 - direction_dot^2)          // resists relative tilt
```

The usual reason for multiple links -- a lever arm against sliding and tilt --
is already covered by the perpendicular and alignment terms.

## Phase 1 — parameter-only collapse (NO code change)

The machinery is already parameterised by N.

```
--aa_prebundle_links_per_pair=1
--aa_stiffness_reference_links=1
```

- `_evaluate_aa_association`: `if (occupied >= 1) return` -> one link per pair
- `_local_aa_attachment_site`: `site_count=1`, `site_index=0`,
  `fraction = 1/(1+1) = 1/2` -> midpoint of the accessible overlap
- `aa_site_exclusion` becomes unreachable (only one link per pair)

**`aa_stiffness_reference_links=1` is mandatory, not cosmetic.** It is the
divisor in

```cpp
inline double fixed_single_link_stiffness(double legacy_pair_stiffness, int reference_links) {
    return std::max(0.0, legacy_pair_stiffness) / std::max(1, reference_links);
}
```

Leave it at 3 and the single link inherits `k_aa/3`, i.e. a pair bond 3x too
soft. That would be a silent physics change masquerading as a coarse-graining.

### Decision: does k_on need compensating?

Association rate is `(N - occupied) * k_on * overlap_extensive_factor(...)`.
For an unbonded pair that is `3*k_on` at N=3 and `1*k_on` at N=1, so the rate
at which a *pair* first becomes crosslinked drops 3x.

Two defensible choices:

- **A (matched pair-formation rate):** `k_on = 2500 -> 7500`. Compares the
  coarse-graining at equal pair kinetics. Recommended for the first A/B, since
  it isolates "one bond vs three" from "fewer formation attempts".
- **B (unchanged k_on):** accepts a 3x slower pair formation as part of the
  coarse-graining.

Run A first. If A matches the control, B tells you whether formation rate
mattered at all.

## Phase 2 — validation against control 37876203

Control: A800/M320, 9 heads/half, 50 ms, pull-through on, N=3.
Arm: identical + Phase 1 parameters.

Gates:

1. **Pair-edge stiffness preserved.** `pair_tension_sum` distribution should
   match the control within noise. If it is ~3x lower, `aa_stiffness_reference_links`
   was not set.
2. **Distinct crosslinked pairs**, not link count. Compare
   `_local_aa_distinct_neighbors` / `aa_neighbor_slots.degree()`, since raw
   `aa_local_status1_count` legitimately drops ~3x and will otherwise read as a
   collapse.
3. **status-2 count and 1->2 rate** at matched model time.
4. **Rupture fraction** (`length > aa_link_rupture_length` share of deaths);
   currently 90.8%.
5. **Nematic S** unchanged within seed scatter.
6. **Wall time** — expect the `_update_local_aa_links` ordered loop to shrink
   ~3x; this is the top barrier site in the profile.
7. **>= 12 seeds, report P(assembled).** This system is bistable per seed;
   3-seed means are meaningless.

## Phase 3 — register-minimising site selection (code change)

Only worth doing once N=1, because with N>1 the sites must be distinct and
evenly spaced.

Today the site is the midpoint of the accessible interval and **registration
plays no part in choosing it**:

```cpp
const double fraction = (site_index + 1) / (site_count + 1);
const double target   = geometry.lower + fraction * overlap;
```

while promotion later tests

```cpp
registry_error = sqrt(0.5 * (axial_error_first^2 + axial_error_second^2))
axial_error_first  = ds_first  - target_overlap
axial_error_second = ds_second + target_overlap
target_overlap     = aa_antiparallel_overlap * 0.5*(1 - direction_dot)
```

So registration is a rejection test on a configuration the model never tries to
achieve. Replace `target` with the value minimising `registry_error`, clamped
to `[lower, upper]` (the accessible interval, already barbed-zone limited).
Guard behind a new flag, default off, so it is A/B-able:

```
--aa_site_minimises_register=true
```

Fall back to the midpoint when the minimiser falls outside `[lower, upper]`.

## Phase 4 — optional: axial register at status 1

`aa_register_all_statuses` exists and is `false`. With it off,
`axial_stiffness_ratio = 0` for status 1, so a nascent link exerts **no axial
force toward register** -- it can only wait for sliding to deliver the pair.
Enabling a weak status-1 register term makes links actively anneal toward
register. Test after Phase 3; they are alternative routes to the same end and
should be compared, not stacked blindly.

## Phase 5 — code simplification (only after physics is validated)

Not required for correctness; this is the performance and clarity payoff.

- delete the ordinal probe loop in `_evaluate_aa_association`
- delete `_local_aa_site_available` and `aa_site_exclusion`
- `aa_links_by_pair` vector -> single optional slot
- keep `aa_prebundle_links_per_pair` as a compatibility no-op or remove from CLI

Expect ~3x fewer iterations in the `#pragma omp for ordered schedule(static,1)`
loop, which is 96/384 stack samples and the single largest barrier site.

## Risks

1. **No redundancy.** A pair now dies on its first rupture. Mostly illusory:
   `aa_neighbor_slots.erase(i, j)` is already *pair-level*, so one geometric
   rupture destroyed pair identity even with three links surviving. But the
   mechanical bond now vanishes at the same moment.
   Mitigation: run this campaign together with `aa_link_rupture_length=0`
   (the model **default** -- production explicitly enables it at 0.12), which
   removes 91% of deaths and is the single highest-value one-parameter test
   outstanding.
2. **Reporting discontinuity.** Every status-1 count quoted in today's analysis
   is in units of links, not pairs. Recompute in pairs before comparing.
3. **Site exclusion silently gone.** Intended, but note it in the params file so
   the run is self-describing.

## Suggested first campaign (4 arms x >= 12 seeds)

| arm | N | reference_links | k_on | rupture_length |
|---|---|---|---|---|
| control | 3 | 3 | 2500 | 0.12 |
| single | 1 | 1 | 7500 | 0.12 |
| no-rupture | 3 | 3 | 2500 | **0** |
| single + no-rupture | 1 | 1 | 7500 | **0** |

This separates the coarse-graining from the rupture channel, which are
independent and both suspected.
