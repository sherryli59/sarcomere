# Porting the legacy (am_kmc) AA bond into the current model

Written 2026-08-06. Source: `sarcomere_am_kmc` @ ef6e4d0 (2026-05-28),
`cpp/src/interaction.cpp::aa_energy`.

## 1. What the legacy bond actually is

```cpp
raw_dist = segment_segment_distance(a, b, c, d, box, periodic);   // SHORTEST filament-filament
angle_energy = 0.5 * kappa_aa * (1 - dot^2);
if (cutoff > 0 && raw_dist > cutoff) return angle_energy;         // outside cutoff: alignment ONLY

dist = (raw_dist - optimal) / optimal;                            // NORMALISED radial extension

ds  = (c - a) . dir1;                                             // barbed-end axial offset
U_ax = 0;
if (ds < 0) U_ax  = 0.5 * (10*k_aa) * (ds/optimal)^2;             // one-sided anti-crossing wall
ds_pos = (ds - 0.1*length1) / (0.1*length1);
U_ax += 0.5 * k_aa * ds_pos^2;                                    // harmonic toward ds = 0.1*L

energy = 0.5*k_aa*dist^2 + angle_energy + U_ax;
```

Three terms, not two. Legacy values: `k_aa=10`, `kappa_aa=100`,
`aa_optimal = am_optimal = 0.03`, `aa_cutoff = am_cutoff = 0.05`.

**No attachment sites anywhere.** No `s_i`/`s_j`, no site exclusion, no
ordinals, no barbed-zone interval, no per-pair link multiplicity. The bond is a
property of the filament pair.

## 2. Two corrections to earlier analysis in this session

**(a) Legacy springs are NOT soft.** The radial extension is normalised by
`optimal`, so `k_aa` is dimensionless and

```
k_eff = k_aa / optimal^2 = 10 / 0.03^2 = 11,111 pN/um
```

which is 11x the current `k_aa=1000` and 3.7x `am_tether_k=3000`. The `ds<0`
wall is 10x stiffer again, ~111,000 pN/um. Earlier in this session I claimed
legacy was ~300x softer and that this was why dt=1e-5 worked. That was wrong.

**(b) Legacy dt=1e-5 was cap-dominated, not stable.** With mu_actin = 70.2
um/(pN*s), tau = 1/(k_eff*mu) = 1.28 us, so dt=1e-5 is **7.8 tau** -- four times
past the 2*tau hard limit. It did not blow up because
`max_actin_displacement=0.005` clamps every step. So legacy's dt is not
evidence that dt=1e-5 is safe here; it is evidence that legacy was integrating
through a limiter. Cf. the dt ladder in this session, where nothing NaN'd even
at 1.7 tau because the same caps absorbed it.

Implication: do NOT port legacy's dt along with the bond. Set dt from the
cap-counter ladder as usual.

## 3. Mapping legacy terms onto the current model

| legacy term | current counterpart | difference |
|---|---|---|
| `segment_segment_distance` | `connector` between chosen sites `s_i`,`s_j` | **current is site-based; legacy is site-free** |
| `0.5*k_aa*((d-opt)/opt)^2` | `0.5*radial_k*(len - preferred)^2` | normalised vs absolute; different units for k |
| `0.5*kappa_aa*(1-dot^2)` | `0.5*alignment_k*(1-dot^2)` | identical |
| `U_ax` harmonic to `ds=0.1*L` | status-2 saturating log-cosh register to `aa_antiparallel_overlap=0.10` | **legacy always on; current status-2 only.** Same target value |
| `U_ax` one-sided `ds<0` wall | *(none)* | current has no anti-crossing term |
| cutoff = switch-off | `aa_cutoff` capture + separate `aa_link_rupture_length` rupture | legacy: one reversible threshold. Current: two, and rupture erases pair identity |
| *(none)* | `perpendicular_k * (proj_i^2 + proj_j^2)` | current-only; the term to drop at status 1 |

## 4. Design

Add `--aa_bond_model={site,legacy}`, default `site`, so both run in one gate.

`legacy` path in `compute_local_aa_link`:

1. `raw_dist = segment_segment_distance(...)` replaces the site connector.
   `s_i`/`s_j` become unused (keep them in the struct for the site path and for
   diagnostics; do not delete).
2. Radial: keep the **current** absolute form `0.5*k*(d - preferred)^2` rather
   than legacy's normalised one, so `k_aa` stays in pN/um and the value is
   comparable across the two paths. Document that legacy `k_aa=10` maps to
   `k_eff = 11,111`, and do NOT copy 10 across.
3. Alignment: unchanged.
4. Axial: reuse the existing saturating register term (bounded is strictly
   better than legacy's unbounded harmonic; the current code comment says the
   force-only unbounded version was already removed once for causing runaway
   prestress). Gate it by status as today, with `aa_register_all_statuses` to
   test legacy's always-on behaviour.
5. Anti-crossing wall: port as an **optional** one-sided term
   `aa_barbed_crossing_k` (default 0 = off). Legacy used `10*k_aa`; that should
   be re-derived, not copied.
6. Force application: legacy differentiates the energy w.r.t. centre and
   direction, so torque falls out of the closest-approach geometry rather than
   a site lever arm `-s_first * connector_gradient`. The legacy path must do the
   same or momentum/torque balance will not hold. `actn2_motor_removal_test`
   already asserts `total_force.norm() < 1e-8`; add a torque-balance assertion.

## 5. Open questions to settle before coding

- **Status gating.** Legacy has no status-dependent mechanics at all -- one bond
  law, with status only affecting lifetime. Do we keep the status-1/2 split
  (radial+alignment, then +axial) or go fully legacy (all terms always)?
  Recommend keeping the split; it is the only mechanism that distinguishes a
  nascent encounter from a mature Z-disc link.
- **Cutoff semantics.** Legacy's cutoff is a reversible switch-off: beyond it the
  radial term vanishes but the pair survives and re-forms on return. Current
  treats exceeding `aa_link_rupture_length` as an identity-destroying event
  (`aa_neighbor_slots.erase`). Recommend porting the reversible semantics --
  this is the third option not tested in the 96-run gate:
  1. rupture=0.12 (current), 2. rupture=0 (no geometric channel),
  3. **reversible switch-off at the capture cutoff, identity retained.**
- **Multiplicity.** With a site-free bond, `aa_prebundle_links_per_pair` has no
  meaning -- there is one bond per pair by construction. The gate showed
  links=1 loses to links=3 *only when rupture is enabled* (6.3 vs 11.3); with
  rupture off the gap nearly closes (10.7 vs 11.7). Since the legacy path
  implies multiplicity 1, it should be tested together with reversible cutoff,
  not against rupture=0.12.

## 6. Validation gates

1. `compute_local_aa_link` site path bit-identical when `aa_bond_model=site`.
2. Force AND torque balance on the legacy path (extend
   `actn2_motor_removal_test`).
3. Cap counters at the dt=5e-7 baseline -- especially given section 2(b), do not
   let the legacy path be quietly stabilised by the displacement clamp.
4. Connector-length distribution: with a reversible cutoff, lengths should pile
   up below `aa_cutoff` rather than at a rupture threshold.
5. Gate sweep {site, legacy} x {rupture 0.12, 0, reversible} x >= 6 seeds, on a
   rig large enough that `max_aa_bonds_per_actin` actually binds -- the 4-actin
   balanced rig never reaches degree 2, so that factor was untestable in the
   96-run sweep.
6. Animation for every run (standing instruction from 2026-08-06).
