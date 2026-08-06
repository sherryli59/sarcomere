# Full-simulator variable-length overlap diagnostics

Date: 2026-07-30  
Scope: isolated six-filament sandbox only; **no production run was launched**

## Critical diagnostic erratum

The original sandbox driver called `set_am_radial_geometry()` but omitted the
production initialization call to `set_aa_spacing_geometry()`. Consequently,
`am_preferred_centerline` remained zero. `compute_am_head_tether()` treats a
non-positive target as disabled, so the nominal 3000 pN/µm occupied-pair A–M
radial tether was silently absent from the campaign and from the original
animations.

This invalidates the campaign's conclusions about A–M occupancy loss and
geometric rupture. It does not constitute a production-code defect: the
production driver sets the preferred centerline correctly.

The sandbox now follows the production path and sets
`am_preferred_centerline = 0.055 µm`. In matched 141-head/half controls:

| condition | geometric ruptures, first 30 ms | mean occupied heads, first 30 ms |
|---|---:|---:|
| old deterministic diagnostic | 76 | 8.36 |
| corrected deterministic diagnostic | 0 | 84.5 |
| corrected full Brownian diagnostic | 4 | 84.8 |

Over 300 ms, the corrected deterministic and full-Brownian trajectories retain
mean terminal-window occupancies of approximately 82 and 87 heads,
respectively. They still develop large axial/orientational distortions, so
stable A–A and A–M overlaps remain unresolved; however, wholesale motor
detachment is no longer the limiting failure.

All original campaign tables below are retained as an audit trail but must not
be used to infer motor-retention behavior.

## A–A rotation bug found after the corrected-retention animation

The native turnover-capable ACTN2 local-link mechanics are **not**
sign-reversed. Their translational and orientational forces pass
finite-difference checks against the implemented link energy, and the
first/second filament forces and torques are scattered reciprocally.

The bug is in the separate diagnostic A–A register positive control. Its
coordinate is the signed separation of the two actin **left endpoints**,

`delta = (center_i - L u_i/2)_x - (center_j - L u_j/2)_x`,

but the old implementation applied only the corresponding center forces. It
omitted the two orientational gradients associated with `u_i` and `u_j`.
Consequently, off-center A–M head torque could change the register coordinate
by rotating the rods without receiving any restoring reaction from the
register spring. The spring was nonconservative for rotating rods.

An opt-in conservative correction now applies both the center force and the
matching endpoint-orientation drive. Matched deterministic seed-7 diagnostics
(141 heads/half, 50 ms fixed-lattice pre-equilibration, 150 ms release) give:

| diagnostic | final A–A overlap (µm) | actin tilt (rad) |
|---|---:|---:|
| old endpoint force only | 0.581 | 0.894 |
| force plus endpoint torque | 0.050 | 0.250 |
| corrected register, A–M head torque disabled | 0.073 | 0.013 |
| corrected register, native ACTN2 mechanics disabled | 0.053 | 0.235 |

The ACTN2 knockout changes little in this short positive-control diagnostic,
which directly rules out a backwards native A–A torque as the cause. The
off-center A–M head moments are the dominant rotational drive; their
`r × F` sign and conversion to director drive are mechanically correct, but
the incomplete register previously provided no conjugate rotational
resistance. Increasing the register cap from 12 to 100 pN has no effect here,
because the corrected trajectory never reaches the original cap.

The corrected A–A overlap settles below the 0.10 µm target rather than exactly
at it because active motor forces continue to oppose the finite 180 pN/µm
register spring. This is a force-balanced displacement, not a sign error.

## Deterministic force/torque attribution

At the pre-equilibrated symmetric release geometry, the native interaction
calculation has no transverse force on any filament. Its actin forces are
axial (33.7, -58.3, -46.7, and 52.6 pN for A0--A3). The director-relevant
perpendicular actin torques are 1.93--2.58 pN µm; these generate later tilt but
not the initial transverse separation. The much larger ±66.7 pN µm axial
components are projected out of actin director tumbling.

The force-source audit is:

| source | release contribution | role in original rupture |
|---|---|---|
| explicit weak/strong heads | axial force plus small off-center torque | secondary rotation/axial motion |
| occupied-pair radial tether | **zero because target centerline was unset** | missing restoring force; primary setup defect |
| occupied-pair angular tether | zero inside the initial 2° slack band | acts only after tilt |
| ACTN2 links | internal A–A mechanics | not the initial A–M radial drive |
| A–M/A–A sterics | zero at the initialized non-overlapping distances | not causal |
| diagnostic titin | full 3-D center force, no endpoint torque | primary outward radial drive on M0 side |
| diagnostic A–A register | zero at the 0.10 µm target and axial-only | not causal |
| axial-cell update | affine x-only rescaling | not causal for radial escape |
| Brownian translation/rotation | zero in the deterministic control | not required |

The diagnostic titin itself is asymmetric at this topology. For A0/A2--M0,
the endpoint distance is 1.0764 µm, shorter than the 1.125 µm rest length.
The compression-bearing Hookean spring therefore pushes each actin away from
M0 with a -0.546 pN transverse component while pushing M0 oppositely with the
sum from both springs. The initial relative radial speed predicted from the
configured rod mobilities is approximately 37 µm/s, enough to traverse the
45 nm from the initialized 55 nm centerline spacing to the 100 nm rupture
threshold in about 1.2 ms. The observed first deterministic rupture is at
0.840 ms.

For A1/A3--M1, the endpoint distance is 1.1763 µm, so the same springs are in
tension and initially pull that side together. Thus the old deterministic
failure begins on M0, not symmetrically.

With the radial target corrected, the 3000 pN/µm occupied-pair tether activates
as soon as the pair leaves its 0.055 ± 0.002 µm flat-bottom band. It eliminates
all deterministic geometric rupture over the first 30 ms despite the
compression-bearing titin. This causal ablation identifies the combination of
an unset radial target and the diagnostic titin's compressive 3-D force as the
source of the apparent spontaneous detachment.

## What was implemented

All changes are under this sandbox; the original source tree is untouched.

- An opt-in axial cell with fixed, zero-stress, or series-spring boundary
  modes. Fixed mode is the default and a strict no-op.
- Pairwise minimum-image A–M and A–A generalized-force bookkeeping.
- Affine axial cell rescaling and neighbor-list rebuilding.
- An always-connected diagnostic Z–M/titin spring.
- A persistent A–A register positive control.
- A one-sided opposing-side/cross-midpoint force surrogate, with restoring
  and assisting signs.
- A native bidirectional-crossbridge comparison.
- Independent control of rotational Brownian noise for the diagnostic
  ablation.
- Full output of all 140 myosin head anchors, weak/strong binding state,
  binding site, signed connector extension, stroke progress, and active force.

The full simulator retains explicit weak and strong heads, turnover-capable
ACTN2, 3-D anisotropic rod diffusion, off-center head torque, sterics, and
periodic force-balanced 4-actin/2-myosin topology.

## Verification

All **23/23** C++ tests pass. The new tests cover:

- fixed-mode neutrality;
- minimum-image virial invariance;
- zero net added particle force;
- cell integration and clamping signs;
- affine fractional-coordinate and orientation preservation;
- the exact initial A–A and A–M overlaps;
- compressed-titin force sign;
- A–A register restoring force;
- exact conservative endpoint-orientation gradient for the corrected register.

The capped and uncapped versions of the representative full-diffusion
trajectory are identical. Its recorded force, torque, and rotation-cap counts
are all zero, so its failure is not a clamping artifact.

## Main campaign

The screen comprised:

- 6 mechanisms;
- initial A–A overlap = 0.04, 0.10, or 0.16 µm at `Lx = 4.3 µm`;
- initial `Lx = 3.9, 4.3, or 4.7 µm` at A–A = 0.10 µm;
- 4 seeds;
- 0.3 s per trajectory at `dt = 1 µs`;
- 120 trajectories and 36 million integrated steps.

Every arm used the same zero-stress cell boundary. The diagnostic titin rest
length was 1.125 µm, corresponding to A–M = 0.625 µm in the symmetric
A–A = 0.10 µm geometry.

Tail statistics below pool the five initial conditions and four seeds.

| Mechanism | mean A–A (µm) | A–A spread | mean A–M (µm) | A–M spread | `Lx` spread (µm) | mean strong heads | result |
|---|---:|---:|---:|---:|---:|---:|---|
| turnover ACTN2 | 0.144 | 0.425 | 0.614 | 0.219 | 0.800 | 0.08 | fail |
| titin | 0.146 | 0.323 | 0.688 | 0.130 | 0.729 | 0.23 | fail |
| register + titin | **0.0966** | **0.0169** | 0.636 | 0.387 | 0.454 | 0.18 | fail |
| opposing-side restoring surrogate + titin | 0.217 | 0.970 | 0.680 | 0.461 | 0.746 | 0.14 | fail/runaway |
| opposing-side assisting sign control + titin | 0.538 | 2.267 | 0.500 | 0.691 | 0.685 | 0.11 | fail/runaway |
| bidirectional heads + titin | 0.144 | 0.292 | 0.686 | 0.129 | 0.727 | 0.27 | fail |

All mechanisms also fail the terminal-slope criteria. No arm maintained an
average of one strong head across the campaign. Native status-2 ACTN2 never
formed; the retained links were status 1.

## Longer positive control

Register+titin was extended to 1 s from all five initial geometries for one
seed.

- A–A tail means converged to 0.096–0.100 µm.
- `Lx` tail means remained separated from 3.905 to 4.250 µm.
- A–M tail means remained separated from 0.119 to 0.677 µm.
- Four of five trajectories had essentially no bound strong heads in the
  terminal window.
- The long-repeat trajectory reached a mean actin tilt of 1.41 radians and
  lost A–M overlap despite the axial targets.

This rules out slow axial-cell relaxation as the sole explanation.

## Rotational-noise ablation

At the reference initial geometry (`Lx = 4.3 µm`, A–A = 0.10 µm),
register+titin was repeated for four seeds:

| Noise mode | actin tilt (rad) | weak heads | strong heads | A–M (µm) |
|---|---:|---:|---:|---:|
| full translation + rotation | 0.194 | 0.45 | 0.37 | 0.644 |
| translation, no rotational Brownian noise | 0.063 | 2.55 | 1.09 | 0.635 |
| deterministic mechanics | 0.044 | 0.20 | 0.12 | 0.628 |

Rotational Brownian noise materially accelerates disengagement, but deleting
it is not a physical solution. Deterministic mechanics also loses the motor
population because Brownian encounters are part of rebinding. The needed
change is an interaction that retains the 3-D lattice while leaving axial
coordinates free.

## Fixed-lattice chemistry pre-equilibration

The reference `register+titin` geometry was held completely fixed for 50 ms
while the existing calcium/thin-filament and weak/strong-head kinetics were
advanced. The filaments were then released for 300 ms with full translational
and rotational Brownian diffusion. One matched seed compared 35 versus 141
heads per myosin half. Pre-equilibration events and numerical cap counters
were cleared at release.

| heads per half | release thin `q` | release weak + strong | first sampled zero occupancy | mean occupancy, 0–50 ms | mean occupancy, 200–300 ms | final occupancy |
|---:|---:|---:|---:|---:|---:|---:|
| 35 | 0.832 | 19 + 20 | 24.5 ms | 9.65 | 0.239 | 0 |
| 141 | 0.838 | 36 + 48 | 56.0 ms | 29.76 | 0.353 | 0 |

The additional heads delay the first complete disengagement and raise early
occupancy, but neither condition retains an engaged A–M lattice. In the
terminal 100 ms, mean strong-head occupancy is 0.000 for 35 heads per half and
0.050 for 141. Final A–M overlaps are 0.633 and 0.629 µm, respectively, but
this apparent axial similarity is supplied by the diagnostic target forces,
not sustained motor engagement.

Post-release transition accounting closes exactly:

- 35 heads/half: 39 release-bound + 179 attachments =
  93 weak detachments + 61 strong detachments + 64 geometric ruptures.
- 141 heads/half: 84 release-bound + 588 attachments =
  298 weak detachments + 150 strong detachments + 216 geometric ruptures +
  8 capacity-limit removals.

Thus the calcium-initiation transient is not the root cause. Removing it
provides a larger initially strong population, and the physiological head
density extends engagement by tens of milliseconds, but rotation/transverse
drift still closes the association gate and drives occupancy to zero.

## ACTN2 cohesion audit and sandbox correction

The apparent A0–A1 gap had two independent causes, plus one visualization
error:

1. The direct diagnostic constructed `Sarcomere` without selecting the
   production driver's `productive_register` ACTN2 status classifier. It
   silently retained `legacy_motor_identity`. At 23.883 ms, the first A0–A1
   link therefore detached as a short-lived status-1 link (0.1 s assigned
   lifetime), even though 8.118 pN of opposed productive motor support was
   available. The pair was already outside the 0.0706 µm capture cutoff, so
   the lost link could not immediately rebind.
2. The original aggregate `k_aa = 10 pN/µm` gives only 3.33 pN/µm =
   0.00333 pN/nm per link with three reference links. With the correct
   classifier, both links matured and survived, but stretched to about
   225 nm and supplied only about 0.7 pN per link against motor forces of
   tens of pN. `aa_perpendicular_k` controls connector/rod geometry; it is
   not a radial cohesion stiffness.
3. The old movie drew the intended pair as a dashed ACTN2 connector even
   after no link existed. Its A–A trace also used only axial endpoint
   separation, so a pair could look axially registered while separating by
   hundreds of nanometres in y/z.

The diagnostic now records full pair-resolved 3-D geometry, all live ACTN2
attachments and their material attachment points, connector force/torque,
status, and every bind/mature/detach event. It also explicitly selects
`productive_register`; `legacy_aa_status` remains available as a control.

A sandbox-only candidate correction combines the production classifier with
100× ACTN2 radial stiffness (`0.333 pN/nm` per link, `1 pN/nm` aggregate).
This does not alter a production default.

| matched 300 ms full-Brownian run | A0–A1 full 3-D separation (µm) | A2–A3 full 3-D separation (µm) | final live links |
|---|---:|---:|---:|
| legacy classifier, original stiffness | 0.213 | 0.596 | 0 + 0 |
| productive classifier, 100× stiffness | 0.109 | 0.108 | 2 + 2 |

The corrected run retains cohesion through 300 ms with both translational and
rotational diffusion. It still permits turnover: one mature A0–A1 link
detached at 276.929 ms, rebound 1.012 ms later, and rematured 1 µs after
rebinding. A matched deterministic run ends at 0.107 and 0.109 µm, showing
that the correction is not merely Brownian trapping. Head torque is the
dominant deterministic driver of the original transverse opening; disabling
it with original stiffness leaves A0–A1 at 0.108 µm after 150 ms.

This is a successful cohesion diagnostic, not yet a calibrated ACTN2 model.
The 100× value must be checked against rupture-force/load-dependent turnover,
multi-seed behavior, and timestep convergence before any production use.

## Annealing and slide–stall challenge

The corrected cohesion run was **not** an annealing test. Its diagnostic
initialization explicitly left `aa_axial_register_ratio = 0`, so local ACTN2
could maintain radial proximity and resist instantaneous material-site shear,
but it had no force that selected the 0.10 µm A–A register after turnover.
Relaxed rebinding at the current coordinates erases rather than restores an
axial target.

An off-register fixed-cell challenge started at A–A = 0.16 µm and mean A–M =
0.495 µm. Enabling the local barbed-end register only for status-2 links
initially annealed A–A to about 0.10 µm. It was not durable: brief loss of
productive motor support dematured links to status 1 and simultaneously
removed their axial mechanics. This creates a destabilizing feedback precisely
near stall, where motor support fluctuates.

The sandbox candidate therefore separates mechanics from lifetime state:
every attached ACTN2 link supplies the barbed-end register force, while status
1/2 still controls its turnover lifetime. Three sites per pair provide
redundancy without changing the aggregate stiffness. This is default-off and
does not alter production behavior.

That correction exposes a more fundamental result. ACTN2 alone can anneal the
internal A–A coordinate, but it cannot center the common A–A pair relative to
the two thick filaments. Its reciprocal internal forces sum to zero on the pair
center of mass. In the 1 s full-Brownian ACTN2-only run, the system showed a
transient A–M plateau, then a head-occupancy asymmetry translated one Z pair
toward one myosin and away from the other; replacement ACTN2 eventually became
geometrically inaccessible and A2–A3 lost cohesion.

A reciprocal titin-like A–M centering spring was then used as a positive
control while keeping the periodic cell length fixed:

| full-Brownian positive control | release | 100 ms | 1 s |
|---|---:|---:|---:|
| mean A–A overlap (µm) | 0.160 | 0.091 | 0.089 |
| mean A–M overlap (µm) | 0.495 | 0.532 | 0.532 |
| A0–A1 / A2–A3 full 3-D separation (µm) | 0.168 / 0.168 | 0.099 / 0.101 | 0.097 / 0.104 |
| live ACTN2 links | 5 | 6 | 6 |

Thus this run does display the requested sequence: rapid A–M sliding over the
first approximately 100 ms, followed by a 900 ms plateau while A–A anneals
toward 0.10 µm. It also retains real turnover. An A2–A3 link detached at
26.682 ms, rebound 11 µs later, and rematured 1 µs later; an A0–A1 link lost at
365.534 ms rebound at new material coordinates at 888.778 ms and rematured the
next microsecond.

This is not evidence that motors plus ACTN2 alone select finite A–M overlap.
It is a positive control demonstrating that the requested behavior appears
once the otherwise neutral A–M centering coordinate is supplied. The physical
identity and calibration of that centering interaction remain open.

The challenge also uncovered and corrected an A–M diagnostic error: periodic
wrapping had been applied independently to both actin endpoints, which could
make a seam-crossing 1 µm actin appear to overlap 1.5 µm of myosin. The metric
now wraps the actin center displacement once and reconstructs both endpoints
from the same periodic image.

## Interpretation of the candidate mechanisms

1. **Turnover ACTN2 alone:** transient elastic links do not select either
   overlap.
2. **Titin alone:** improves local A–M spacing but leaves A–A free and does not
   preserve motor engagement.
3. **Register+titin:** proves the A–A force sign and target are correct, but
   remains a prescribed two-target control, not an emergent mechanism.
4. **One-sided opposing-side surrogate:** fails the full topology and has a
   runaway branch. It should not be promoted to production.
5. **Bidirectional existing heads:** does not reproduce actual other-half
   capture and does not improve convergence.

The cross-midpoint biological idea therefore remains untested at the explicit
head-ownership level. A valid implementation must allow a head from the
opposite thick-filament half to bind only after a geometrically reachable
midpoint crossing, while preserving reciprocal force/torque and continuous
kinetics.

## Still unresolved

1. What physical interaction retains transverse lattice spacing and
   orientation in this minimal force-balanced unit: neighboring thick
   filaments/myomesin, Z-disc lattice coupling, titin geometry, or explicit
   radial confinement?
2. Can that interaction preserve at least one active strong head and
   turnover-capable ACTN2 for at least 1 s under full rotational diffusion?
3. Can explicit opposite-half head capture then provide the missing
   overlap-dependent active feedback without a discontinuity or runaway?
4. Does one parameter set still support both free-actin sliding/unbinding and
   balanced partial-overlap arrest?
5. What physical external compliance should set the global cell mobility and
   zero-stress/series boundary?
6. The diagnostic cell virial still needs native steric and any future
   M-line/Z-disc interaction components before it can become a production
   integrator.

## Recommended next tests

1. Add a default-off transverse/orientational lattice-retention interaction
   with no axial component or preferred axial register.
2. Require survival and force/torque reciprocity before repeating any overlap
   sweep.
3. Implement actual opposite-half head eligibility and retention as a separate
   switch; do not reuse the failed phenomenological surrogate.
4. Repeat the five-initial-condition/four-seed test under full diffusion.
5. Only after a mechanism passes, test cell mobility and time-step convergence,
   then run the F–V sweep. Do not launch production before those gates pass.
