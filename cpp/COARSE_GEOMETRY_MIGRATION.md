# Coarse-geometry baseline for the explicit-head model

## Reference implementation

- Working coarse-grained baseline: commit `ef6e4d0` (`am_kmc`).
- Preferred-overlap/regular-spacing change: commit `6f69152`.
- Representative assembled sweeps used `myosin_radius=0.058`,
  `am_cutoff=0.08`, `am_optimal=0.06`, `max_myosin_bonds=6`, and a weak
  transverse confinement `k_bundle_max=0.1` ramped over 5000 steps.

## Geometry retained

1. `actin.direction` points from the barbed (`left_end`) toward the pointed
   (`right_end`) end.
2. AM accessibility and polarity continue to come from the coarse model's
   `geometry::analyze_am` overlap and partial-overlap construction.
3. AM and AA radial interactions retain preferred centerline separations.
4. Angular alignment is nematic (`1-(u_i dot u_j)^2`), so myosin has no
   artificial plus/minus polarity.
5. AA axial registry is defined by the two barbed-end regions. Parallel pairs
   register their barbed ends; antiparallel pairs target a finite 0.1-um
   interdigitation, recovering the old 0.1-actin-length overlap in a
   pair-symmetric form.

### Meaning of the axial-overlap target

The pair-level AA edge coarse-grains an ACTN2-rich Z-disc overlap zone; it is
not one literal ACTN2 dimer. Cardiac Z-bands contain roughly four to six
longitudinal ACTN2 layers at about 19.2-nm spacing, and six-layer cardiac
Z-bands are about 130 nm wide. Thus the historical 0.1-actin-length target is a
reasonable first overlap scale for a 1-um actin. The present binary AA edge
still does not represent the layer occupancy explicitly. A later refinement
should give each eligible filament pair a finite site capacity or occupancy
`n_ij` while conserving the ACTN2 pool; the current overlap target should be
retained as the mature cluster geometry.

For barbed ends `b_i`, filament axes `u_i`, minimum-image displacement
`delta=b_j-b_i`, and `q=u_i.u_j`, define

```
W_ij = W_AA (1-q)/2
s_i = delta.u_i
s_j = delta.u_j
U_ax = k_aa/4 [(s_i-W_ij)^2 + (s_j+W_ij)^2].
```

This gives zero axial offset to parallel neighbors and `W_AA` overlap to
antiparallel neighbors, without dependence on pair ordering.

## Mechanics deliberately replaced

1. The prescribed coarse AM velocity and overlap-dependent velocity shutdown
   are replaced by explicit detached/weak/strong heads and a finite
   force-velocity working stroke.
2. The passive AM axial end-stop is not restored; axial force belongs to the
   motor heads.
3. Status 2 remains a load state of an existing AA bond, not a switch that
   creates the crosslink.
4. Global x-alignment bias is not active. Weak transverse confinement remains
   a diagnostic control rather than an assumed scaffold-free mechanism.

## Regression found and corrected

The first explicit-head lattice version registered AA-linked filament centers.
That makes an ideal antiparallel Z-disc pair appear approximately one actin
length out of register and generated force caps under isometric conditions.
The AA energy now registers the `left_end` barbed-end overlap zone, following
the working coarse model, while preserving pair symmetry and physical radial
spring units.

## Validation sequence

1. Unit tests: parallel and antiparallel barbed-end rest states, pair-order
   symmetry, radial lattice stability, finite working-stroke bounds.
2. Minimal two-actin test: free/isometric crossed with no-AA/persistent-AA;
   require no pervasive force caps and reduced relative axial motion.
3. Six-actin lattice: verify preferred AM/AA spacing and annealing.
4. Matched assembly pilots: coarse model, explicit heads with migrated
   geometry, and explicit heads plus the historical weak confinement.
5. Only after these pass, launch the 2000-actin/1000-myosin production run.
