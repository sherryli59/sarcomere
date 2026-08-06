# Variable-length validation plan for the full simulator

## Why this is required

`two_sarcomere_unit` is force-balanced internally, but its axial periodic box
length is fixed.  In the symmetric geometry,

```text
Z--M distance = Lx / 4 + A--A overlap / 2
```

so changing `Lx` directly changes the mean A--M overlap.  Particles can
redistribute inside the box, but no internal interaction can make simulations
started at `Lx = 3.9, 4.3, 4.7 um` converge to one common `Lx` while those
three values remain immutable.  A fixed-box comparison therefore tests local
force balance, not selection of sarcomere length.

## Safest next implementation

Add one opt-in axial cell degree of freedom, disabled by default:

```text
axial_cell_mode = fixed | stress_clamp | series_spring
cell_mobility
external_axial_force
series_k
series_rest_length
cell_thermal_scale
```

At each integration step:

1. Integrate all filament translations and rotations normally.
2. Accumulate the pairwise axial virial using minimum-image connectors.
3. Form the generalized cell force

   ```text
   Qx = virial_xx / Lx
        + external_axial_force
        - series_k * (Lx - series_rest_length)
   ```

4. Integrate `dLx = cell_mobility * Qx * dt + cell noise`.
5. Affinely rescale all axial coordinates by `Lx_new / Lx_old`, wrap them, and
   rebuild the neighbor list.

The cell noise should be off for the first mechanical diagnostic.  If enabled
later, its amplitude must satisfy the fluctuation-dissipation relation for the
chosen cell mobility.

## Force bookkeeping required

Record these separately in HDF5:

- instantaneous `Lx` and `dLx/dt`;
- active A--M virial;
- ACTN2 elastic virial;
- any ACTN2 overlap-free-energy/entropic virial;
- titin virial;
- myomesin/M-line virial;
- external or series-spring force;
- mean and per-pair A--A overlap;
- mean and per-half-sarcomere A--M overlap;
- Z-disc--M-line distances.

The sum of the recorded components must reproduce the integrated cell force.

## Interaction variants to port, one at a time

1. **Control:** existing turnover-capable elastic ACTN2.
2. **Johann adhesion:** `U = -gamma * y` within the admissible A--A overlap
   zone.  This is a constant overlap-increasing force and is not by itself a
   preferred-length interaction.
3. **Fixed-pool entropic pressure:** `F = N kBT / y`, with a documented finite
   force cap near zero overlap.  This assumes crosslinkers are longitudinally
   mobile and end-confined; the current discrete ACTN2 springs do not have
   that property.
4. **Persistent barbed-end register:** a symmetric potential centered at the
   desired A--A overlap.  This is a positive control because it explicitly
   supplies the missing A--A restoring curvature.
5. **Titin-like Z--M elasticity:** expose `titin_rest_length` as an input rather
   than computing it silently.  Test both always-connected and the current
   bond-state-dependent implementation.
6. **Series elasticity:** a shared external/structural spring acting on the
   repeat length, not separately on each A--M pair.  This preserves local
   motor--ACTN2 competition and is mechanically different from titin directly
   canceling each motor pair.

Do not combine mechanisms until each single-mechanism force curve and sign have
been verified.

## Unit tests before trajectories

- each reciprocal pair has zero net particle force;
- pair virial is invariant to wrapping either member across the boundary;
- active A--M virial is contractile in the reference geometry;
- compressed titin is expansile and stretched titin is contractile;
- a zero-force cell does not drift;
- affine rescaling preserves fractional coordinates and rod orientations;
- the cell integrator is deterministic for a fixed seed;
- passive spring-only equilibrium is independent of time step;
- disabling the new mode produces bitwise-identical fixed-box behavior.

## Decisive diagnostic campaign

Use the same external boundary condition for every initial state:

```text
A--A initial overlap: 0.04, 0.10, 0.16 um
Lx initial length:    3.9, 4.3, 4.7 um
seeds:                at least 4
```

Do not call a mechanism successful merely because its velocity is small.
Require:

- all initial lengths converge to the same final `Lx`;
- all A--A initial conditions converge to the same overlap near `0.1 um`;
- the common A--M overlap is below the maximum/bare-zone limit;
- terminal slopes of A--A, A--M, and `Lx` are statistically indistinguishable
  from zero;
- small perturbations on either side generate restoring, not merely opposing,
  forces;
- the result persists when the cell mobility and integration time step are
  changed.

Only after this diagnostic passes should the mechanism be ported into a larger
3-D or production sweep.

