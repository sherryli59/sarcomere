# Axial myosin-tip / ACTN2-patch contact

## Purpose

The existing ACTN2--myosin corona is a radial capsule contact. At ordinary
actin--myosin spacing it supplies little axial resistance when actin sliding
brings a myosin tip into the ACTN2-rich Z-disc patch. This optional interaction
adds that missing axial steric resistance without making the entire myosin
filament artificially thicker.

## Geometry and force

For each bilateral ACTN2 bridge, the midpoint of its two material attachment
sites is the local Z-patch position. For a candidate myosin end, let `g` be the
signed gap from the tip to that patch along the outward myosin axis, and `rho`
the transverse separation. A patch is considered only on the same half of the
myosin. The transverse weight is one for `rho <= R_flat`, a smooth cubic taper
between `R_flat` and `R_cut`, and zero beyond `R_cut`.

For `g < g0`, the axial compression and force are

```
delta = g0 - g
F = S(rho) Fmax tanh(k delta / Fmax).
```

Thus `k` controls the small-compression stiffness and `Fmax` prevents an
unbounded force during a rare deep overlap. Equal half-forces act on the two
ACTN2-bound actins and the exactly opposite force acts at the myosin tip.

The radial gate selects the spatial extent of the coarse-grained patch. The
implemented force itself is purely axial; it intentionally does not add a
second transverse attraction or repulsion in the taper region.

## Finite valency

Several ACTN2 molecules may describe the same local Z-disc patch. The simulator
therefore selects only the strongest candidate for each myosin tip on each
step. A left and right tip can each carry at most one contact, and each contact
is bounded by `Fmax`. This avoids multiplying the end-stop by the number of
nearby ACTN2 links.

## Physical accessibility and nonpenetration

Steric accessibility belongs to the thick-filament body, so it is independent
of whether a motor head happens to be bound. The same all-tip geometry is used
for mechanics and kinetics. A new weak bilateral bridge, including tether-arm
closure, is rejected if any nearby myosin tip occupies its axial contact zone.

A status-1 diffusive step cannot enter a new tip zone. If a site is already in
a tip corona, it may only reduce the penetration of that same physical
`(myosin, tip)` contact. This per-tip no-flux rule prevents a move from escaping
one tip by silently entering another. It preserves contour annealing without
allowing diffusion through a thick-filament end.

After every overdamped position update, a mobility-weighted projected
Gauss--Seidel solve enforces `g >= 0` for every active bilateral bridge and
nearby tip. Both linked actins and a mobile myosin share the displacement in
proportion to their translational mobility; a fixed myosin contributes no
motion. A sub-picometer positive numerical buffer prevents roundoff in one
shared-filament constraint from reopening another. This is the hard-contact
limit of steric repulsion, not a link-deletion rule. Projection retains its
strict `0.001 pm` solver tolerance. The final audit permits a `1 pm` residual;
larger residuals are written to the scheduler error log without terminating
the trajectory, and the following step attempts projection again.

Existing status-1 and status-2 bridges are therefore not deleted merely because
a tip approaches. They remain subject to the existing physical detachment,
overstretch, and deep connector-clash pathways. Complete loss of crosslinkable
actin contour also retains the existing detachment rule.

## Options and defaults

The feature is opt-in for backward compatibility:

```
--aa_zdisc_tip_contact=true
--aa_zdisc_tip_onset=0.050
--aa_zdisc_tip_stiffness=1000
--aa_zdisc_tip_force_max=20
--aa_zdisc_tip_transverse_flat=0.100
--aa_zdisc_tip_transverse_cutoff=0.130
```

Lengths are in micrometers, stiffness in pN/um, and force in pN. All settings
are stored in HDF5 metadata and checked on resume. A checkpoint created with
the earlier motor-occupancy-scoped accessibility rule is rejected when this
feature is enabled, because it may already violate the all-tip invariant.

## Diagnostics

Actin frames store `aa_zdisc_tip_force` and
`aa_zdisc_tip_contact_count`. Myosin frames store force, gap, transverse
distance, and gate for each `minus` and `plus` tip. These fields distinguish
successful pre-occlusion contact from subsequent emergency rupture.

`/catch_bond/tip_accessibility_events` audits every accepted association or
tether closure. Its 13 columns are link id, step, event code, both actin ids,
both contour coordinates, nearest myosin id, tip sign, nearest axial gap,
transverse distance, tip-zone occupancy, and motor engagement. Event code 1 is
association and 2 is tether closure. Accepted events must have tip-zone
occupancy zero. Cumulative projection-event and projection-correction features
on actin and myosin report how often the hard constraint was active and the
total corrected penetration.

## Scope

This term is a steric end-stop. It does not contribute to ACTN2 maturation,
motor-support classification, arm opening load, or the existing deep-clash
rupture criterion. In particular, observing `g < 0` does not itself detach a
link; the position constraint prevents that state from being committed.
