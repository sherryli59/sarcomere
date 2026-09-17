# ACTN2 two-arm annealing and A-M boundary contact

Date: 2026-08-19

## Outcome

The simulator now has an opt-in, molecularly interpretable ACTN2 arm-state
model. It replaces continuous weak-link coordinate diffusion and an infinitely
locked mature link with a finite sequence:

1. weak two-arm bridge;
2. mature two-arm bridge;
3. one arm opens while the other remains at its material actin site;
4. the free arm either closes onto a locally accessible actin site or the
   remaining arm dissociates.

The feature is disabled by default (`--aa_arm_state_model=false`), so existing
inputs retain schema-4 behavior and can still resume. Enabling it uses local
ACTN2 schema 5 and requires schema-5 checkpoint state on resume.

## What changed relative to the previous implementation

| Process | Previous model | Two-arm model |
|---|---|---|
| Weak mobility | Both attachment coordinates diffuse continuously along a prescribed barbed-end zone. | Binding coordinates are material sites. Mobility occurs through a discrete opening and local reclosure event. |
| Maturation | Status 2 stops coordinate diffusion and, at zero unlock rate, remains locked until whole-link detachment. | A positive mature-to-weak rate permits conformational relaxation without detachment; the arm-state model additionally gives each mature arm a finite force-dependent opening rate. |
| Partner exchange | A link either keeps both actin identities or disappears into the free pool. | One identity remains physically tethered while the free arm integrates over the continuously reachable, exposed contour of local candidate partners. |
| A-M boundary | Whole-link rupture occurred only when a crosslinkable fraction vanished; a stored attachment could otherwise remain behind the moving boundary. | New binding is excluded from A-M overlap. Existing arms respond to explicit local ACTN2-myosin contact; deterministic opening is reserved for complete loss of crosslinkable contour. |
| Thick-filament crossing | No ACTN2-connector/thick-filament contact. | The finite ACTN2 connector has WCA contact with thick-filament capsules, with equal-and-opposite forces and torques. |
| Lifetime load | Status/catch-bond paths could use aggregate productive motor support. | Arm opening uses only the link's own tensile connector load plus direct ACTN2-myosin contact on that arm. Actin-actin springs, generic volume exclusion, and aggregate motor force are not reinterpreted as arm tension. |

`productive_register`, `arm_load_productive`, and
`motor_transmitted_register` maturation semantics are otherwise unchanged. In
particular, the source-resolved transmitted-register option remains the way to
require a consistent bilateral myosin source. The new arm kinetics controls how
an already formed bridge opens and anneals; it does not invent a new maturation
signal.

## Physical state and kinetics

States are stored in `actn2_arm_kinetics::State`:

- `1`: weak two-arm bridge
- `2`: mature two-arm bridge
- `3`: tethered to stored actin `i`
- `4`: tethered to stored actin `j`

Each bound arm opens through a Bell slip pathway

`k_open(F) = k_open(0) exp(beta F delta)`

using nonnegative direct arm load. The two arms are competing hazards, so at
most one opens per timestep. Weak and mature bridges have separate zero-force
rates. "Mature/locked" now means that the two stereospecific material sites do
not diffuse; it does not mean that their identities can never change.

A tether's free arm enumerates nearby antiparallel actins. For each candidate,
the implementation intersects its exposed barbed-end contour with the sphere
reachable from the surviving material attachment, then subtracts contour
already occupied by another ACTN2 molecule on the same pair-facing surface.
This is a continuous local capture interval, not the canonical pair sites used
for de-novo association. The interval is integrated by midpoint bins no wider
than 10 nm (or half the site-exclusion length when smaller).

The total closing hazard for a candidate contour is

`k_rebind = (k_attempt/L_ref) integral_A exp(-beta theta U_local(s)) ds`.

`L_ref` is `aa_overlap_capture`. `U_local` contains the deformation needed at
the molecular capture transition: connector extension and local
perpendicularity. It deliberately excludes whole-filament alignment and axial
registry energies. Those mechanical terms act after the second arm binds; one
free molecular domain is not required to rotate two entire filaments before it
can close.

Candidate hazards are summed and a candidate is selected in proportion to its
integrated rate. A material coordinate is then sampled inside the selected
contour bin. The surviving bound coordinate is never moved. A separate
first-order rate removes the remaining arm into the free molecular pool.

Finite neighbor valency is checked before a candidate contributes to the KMC
hazard, and a successful closure is inserted into the live pair/site index
immediately. Thus an impossible valency event is not silently consumed and two
tethers cannot select the same pair-facing contour in one timestep.

## A-M boundary semantics

The scalar A-M overlap fraction is an accessibility field, not a resolved
molecular collision. It limits de-novo association and the free arm's capture
contour. Once an ACTN2 arm is bound, partial boundary passage alone no longer
assigns an infinite opening rate. The explicit ACTN2-myosin contact supplies
the local force used by the Bell arm-opening pathway, and deep contact retains
its separate geometric-rupture threshold. If an actin's crosslinkable contour
falls to zero, its arm still opens deterministically. The same rule removes a
one-arm tether whose last bound actin has no exposed contour.

## Molecular conservation

`aa_total_link_pool > 0` represents a finite number of coarse-grained ACTN2
modules in the periodic cell. A bridge and a one-arm tether each consume one
module. Only complete detachment returns it to the free pool; intramolecular
reclosure does not require a free module. The finite pool now also caps
deterministic prebinding. `aa_total_link_pool=0` remains the backward-compatible
reservoir behavior.

All arm events use counter-based random streams keyed by timestep and physical
link ID, preserving OpenMP scheduling invariance.

## Finite ACTN2-thick-filament contact

The segment joining the two material ACTN2 sites is treated as a coarse-grained
connector of radius `aa_actn2_contact_radius`. Its centerline exclusion distance
from a thick filament is

`myosin.radius + aa_actn2_contact_radius`.

The WCA cutoff is exactly this distance. The connector reaction is distributed
between the two actins by the closest-point position along the connector, and
the thick filament receives the equal-and-opposite reaction. Off-center
contacts also apply orientation drives to both actins and the myosin. A force
cap is an explicit coarse-grained load-sharing/numerical regularizer; deep
penetration has a separate geometric opening/rupture threshold.

An optional soft corona can begin outside that unchanged hard core. It
coarse-grains steric/crowding from the myosin-head layer as a harmonic shoulder:

`U_shell(d) = 0.5 k_shell (d_core + width - d)^2`,

with `k_shell = 2 E_shell / width^2`, so `E_shell` is the shoulder energy at
the hard-core surface. The shell continues smoothly inside the core and the
original WCA term is added there. Crucially, geometric rupture still uses only
`max(0, d_core - d)`: increasing the shell width does not move the rupture
surface outward.

When the corona is enabled, a proposed new bridge or tether reclosure that
crosses the hard core is rejected. A proposal lying only within the corona is
allowed with the conservative factor `exp(-beta U_shell)`. This prevents a new
ACTN2 connector from appearing through a thick filament while retaining
thermally accessible configurations in its soft outer layer. Width zero
disables both the shoulder and proposal gating and reproduces the original
contact path.

Contact geometry is evaluated once per pre-existing bridge per step and cached
for force application. Newly associated bridges are evaluated once during the
final force pass. This avoids the original implementation draft's duplicate
segment-distance calculation.

## Runtime parameters and current defaults

| CLI option | Default | Meaning |
|---|---:|---|
| `aa_arm_state_model` | `false` | Enable the complete model. |
| `aa_arm_open_rate_weak` | 20 s^-1 | Zero-force opening rate of each weak arm. |
| `aa_arm_open_rate_mature` | 1 s^-1 | Zero-force opening rate of each mature arm. |
| `aa_arm_bell_distance` | 0.001 um | Bell distance (1 nm). |
| `aa_tether_off_rate` | 5 s^-1 | Loss rate of the remaining arm. |
| `aa_tether_rebind_rate` | 100 s^-1 | Free-arm closing attempt rate per `aa_overlap_capture` of reachable contour. |
| `aa_tether_capture_distance` | 0.10 um | Maximum candidate connector length. |
| `aa_arm_rebind_theta` | 0.5 | Transition-state fraction for elastic gating. |
| `aa_actn2_contact_radius` | 0.005 um | Effective connector radius. |
| `aa_myosin_contact_epsilon` | 0.00428 pN um | WCA energy scale (approximately kBT at 310 K). |
| `aa_myosin_contact_force_cap` | 40 pN | Contact cap; nonpositive means uncapped. |
| `aa_myosin_contact_shell_width` | 0 um | Soft crowded-head-corona width; zero disables it. Recommended first test: 0.015 um. |
| `aa_myosin_contact_shell_energy` | 0 pN um | Corona energy at the hard-core surface. Recommended first test: 0.0214 pN um (5 kBT at 310 K). |
| `aa_occlusion_rupture_penetration` | 0.010 um | Deep-contact geometric threshold; zero disables it. |

These are mechanically safe starting values, not fitted ACTN2 kinetic
constants. The weak/mature opening rates, Bell distance, and tether closing
rate should be swept before interpreting M-order quantitatively.

## HDF5 and restart

Metadata stores every parameter above and the state/event schemas.

- `/catch_bond/local_link_arm_state`: one row per sparse link identity with
  state, state age, last partner, bound-arm flags, arm loads, contact force, and
  penetration.
- `/catch_bond/arm_events`: opening, reclosure, and final tether-detachment
  events.
- `/resume/aa_local_link_arm_state` and
  `/resume/aa_local_link_arm_loads`: exact arm identity and diagnostic state.
- `/resume/aa_arm_events_pending`: unflushed event continuity.

The legacy `local_link_identity` dataset remains present. A tether retains the
old/free partner in that row for provenance; `local_link_arm_state` is the
authoritative statement of which arm is physically bound.

For tether-detachment rows in `/catch_bond/arm_events`, `event_rate=inf`
identifies complete loss of crosslinkable contour; a finite value is the
stochastic remaining-arm off-rate. This makes geometric and molecular loss
separable without changing the event table width.

## Verification

New tests:

- `actn2_arm_kinetics_test`: state helpers, exact event probability, Bell-load
  monotonicity, competing arm hazards, elastic closing suppression, and the
  distinction between partial and complete A-M occlusion.
- `actn2_myosin_contact_test`: legacy hard-core contact, force cap, force
  conservation, endpoint load sharing, the 15-nm/5-kBT corona, unchanged
  hard-core penetration under a wider shell, and the separated case.
- `actn2_arm_state_integration_test`: mature opening to exactly one tether,
  absence of a two-actin spring while tethered, exact checkpoint/resume of the
  bound arm, continuous local reclosure when every canonical site is out of
  reach, exclusion of whole-filament alignment from the capture barrier, and
  enforcement of the finite pool during deterministic prebinding.

The full Release CTest suite passed: 32/32 tests, including the legacy
A-M-occlusion regression. The manual timing harness is
`actn2_arm_benchmark`. Five paired 50,000-step, one-thread runs with one active
bridge and a real capped ACTN2-myosin contact gave:

- legacy mean: 39.92 us/step;
- arm/contact mean: 41.27 us/step;
- mean overhead: 3.4%;
- median overhead: 2.1%.

This is a single-link microbenchmark, not a production scaling claim. Cost
scales with active ACTN2 bridges times nearby thick-filament contact candidates;
the motor-head and global neighbor-list work is unchanged. A production-sized
paired timing should still be done before a long array submission.

The corrected benchmark also has a `tether` mode that performs the continuous
reachable-contour scan every step. Five interleaved 10,000-step, one-thread
runs on the shared login node gave median times of 300.4 us/step (legacy),
309.5 us/step (two-arm bridge), and 291.8 us/step (scanning tether). The host
noise is larger than the differences, so the microbenchmark detects no new
local-capture overhead. The A800/M320 pilot wall time remains the relevant
production-scale measurement.

## Corrected A800/M320 pilot

The matched 100-ms, 10-us-step, seed-20260805 myomesin pilot completed as
Sherlock job `39820550_1` with 101 committed frames, finite coordinates, and a
healthy live-check result. It used 800 actins, 320 thick filaments,
`myosin_radius=0.05 um`, three myomesin sites per thick filament, and a finite
pool of 360 ACTN2 modules. Relative to the earlier arm-state pilot:

| Metric at 100 ms | Previous | Corrected |
|---|---:|---:|
| free-arm reclosures | 3 | 32 |
| reclosures onto a different partner | 1 | 11 |
| total arm openings | 2090 | 1397 |
| final one-arm tethers | 249 | 185 |
| final two-arm bridges | 87 | 82 |
| final active ACTN2 modules | 336 | 267 |
| maximum active modules | 359 | 297 (pool limit 360) |
| M-band order | 0.0762 | 0.0678 |
| Z-band order | 0.0675 | 0.0732 |

Thus the implementation fix produces real partner exchange, reduces spurious
opening by one third, and lowers the tether population by 26%. Because the
finite pool reduces the total active population by 21%, the final absolute
bridge count is nearly preserved and the bridge fraction rises from 25.9% to
30.7%. The 100-ms structural orders do not show improved global M-band order;
that is a longer-time assembly question, not a failure of the now-observed
local annealing pathway.

The corrected job used 16 CPUs and completed its simulation in 7:37, versus
10:06 on 32 CPUs for the prior pilot. This does not indicate a slowdown, though
the finite pool and differing active-link population prevent treating it as a
strict kernel benchmark. Force clipping remained inherited from the matched
mechanical setup: the post-ramp maxima were 70.2% of actin steps and 39.6% of
myosin steps, within the pilot's configured health thresholds.

## Known coarse-graining limits

- ACTN2 has no explicit azimuth or internal hinge; axial site exclusion is a
  pair-facing approximation.
- The contact connector is a single straight segment rather than resolved
  ACTN2 domains.
- Force-dependent opening is a Bell slip pathway only; no fitted multi-pathway
  catch/slip law is assumed.
- The A-M accessibility boundary remains geometrical and does not resolve
  individual actin monomers. Partial-overlap release is therefore delegated to
  the resolved local connector/thick-filament contact rather than inferred
  directly from the scalar boundary.
- Parameter sweeps must establish whether local partner exchange improves
  global Z/M order; the implementation supplies the pathway but does not
  guarantee network-scale propagation.
