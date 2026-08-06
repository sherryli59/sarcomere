# Physical Rod Mobility and Lattice-Stability Plan

## Objective

Replace the phenomenological isotropic diffusion constants with a thermally
consistent rigid-rod mobility model, repair AA reassociation semantics, and
require the interaction model to preserve and recover a mechanically compatible
six-actin lattice at physiological temperature. Diffusion may set the rate of
relaxation; it must not be used to hide an unstable force law.

The initial scope remains scaffold-free: no M-line, Z-disc plate, titin, or
MyBP-C lattice template is added. Failure of the scaffold-free system after the
passive mechanics are made internally consistent is a biological result that
identifies missing organization, rather than a reason to lower temperature.

## Evidence and identifiability

### What is measured

- Dilute, in-vitro F-actin translational diffusion has been measured by dynamic
  light scattering. For polydisperse samples with number-average lengths
  approximately 0.47--0.69 um and weight-average lengths 0.81--1.25 um, reported
  orientation-averaged values are approximately 1.56--2.67 um^2/s at 25 C.
- Single fluorescent F-actin studies have measured length-dependent
  translational anisotropy and rotational diffusion in confined aqueous layers.
- Cardiac live imaging and FRAP establish that actin, ACTN2, and myosin
  components exchange and assemble dynamically, but these measurements usually
  report molecular recovery or structural incorporation, not the Brownian
  motion of an intact 1--1.6 um filament.
- Intact cardiac thick-filament molecular exchange occurs on minute scales and
  is spatially restricted within a sarcomere. This is not a whole-filament
  diffusion coefficient.

### What is not directly measured

There is no sufficiently direct experimental value for all of

```text
D_parallel, D_perpendicular, D_tumble, D_roll
```

for intact cardiac thin and thick filaments during early myofibrillogenesis.
The model should therefore distinguish:

1. a hydrodynamic free-solution prediction;
2. an intracellular effective-drag correction;
3. binding and scaffold constraints, represented mechanically rather than by
   silently changing the intrinsic diffusion coefficient.

## Reference hydrodynamic model

For a cylinder of length `L`, hydrodynamic diameter `d_h`, aspect ratio
`p=L/d_h`, viscosity `eta`, and temperature `T`, use the finite-cylinder
corrections of Tirado, Martinez, and Garcia de la Torre:

```text
D_parallel = kBT [ln(p) + nu_parallel] / (2 pi eta L)
D_perp     = kBT [ln(p) + nu_perp]     / (4 pi eta L)
D_tumble   = 3 kBT [ln(p) + nu_r]      / (pi eta L^3)

nu_parallel = -0.207 + 0.980/p - 0.133/p^2
nu_perp     =  0.839 + 0.185/p + 0.233/p^2
nu_r        = -0.662 + 0.917/p - 0.050/p^2
```

These formulas apply to a dilute rigid cylinder. They are an upper-bound
reference for an assembling cardiomyocyte, not the final intracellular values.
F-actin has a persistence length much larger than the 1 um simulated thin
filament, so the rigid-rod approximation is suitable at this coarse-graining.

The hydrodynamic diameter must be separate from the steric/capture envelope.
Candidate physical inputs are:

| species | L (um) | hydrodynamic diameter (um) | comment |
|---|---:|---:|---|
| thin filament | 1.0 | 0.008--0.012 | F-actin plus regulatory envelope |
| thick filament | 1.5--1.6 | 0.030--0.040 | cardiac thick-filament backbone |

At 310 K and `eta=1 mPa s`, representative midpoint estimates are:

| species | D_parallel (um^2/s) | D_perp (um^2/s) | D_tumble (1/s) |
|---|---:|---:|---:|
| thin, L=1.0 um, d=0.010 um | 3.00 | 1.85 | 16.2 |
| thick, L=1.6 um, d=0.035 um | 1.55 | 0.99 | 3.17 |

The thin-filament orientation average `(D_parallel+2 D_perp)/3 = 2.24
um^2/s` is consistent with the approximate 1.56--2.67 um^2/s dilute F-actin
measurements. This is a useful validation of the starting calculation.

Because spatial myosin heads use a transported material frame, axial roll is a
real degree of freedom. Add `D_roll` separately from tumbling. Initially derive
it from cylinder axial rotational drag and expose it as a recorded parameter.
Do not reuse `D_tumble` without a sensitivity test; rotation about a cylinder's
long axis has different drag from end-over-end rotation.

## Stage 1: establish physical units and temperature

1. Declare the simulation units in one shared header and in every HDF5 file:

   ```text
   length = um
   time = s
   force = pN
   energy = pN um
   stiffness = pN/um
   torque = pN um
   ```

2. Set the physiological default temperature to 310 K and calculate rather
   than independently enter beta:

   ```text
   kBT = 0.00428 pN um
   beta = 1/kBT = approximately 234 (pN um)^-1
   ```

3. Retain `--beta` only as an explicit diagnostic override. Emit a warning and
   metadata field when beta does not match the requested temperature.
4. Remove `beta=1e6` from any test described as physiological. Keep it only in
   a named `athermal_mechanics` test family.
5. Unit-test the conversion at 298 K and 310 K.

Acceptance: the default beta is generated from temperature, and no production
submission can silently inherit the athermal lattice-test beta.

## Stage 2: introduce rod mobility tensors

For rod axis `u`, construct the translational diffusion and mobility tensors:

```text
D_t = D_perp I + (D_parallel - D_perp) u u^T
M_t = beta D_t
```

Replace the scalar translational update by

```text
Delta r = M_t F dt + sqrt(2 dt) B xi,
B B^T = D_t.
```

An efficient square root is

```text
B = sqrt(D_perp) (I-u u^T) + sqrt(D_parallel) u u^T.
```

For orientation, project torque and noise onto the plane perpendicular to `u`,
use `D_tumble`, and update the direction with an exponential-map or quaternion
rotation rather than adding a vector and renormalizing. Apply the axial torque
component to the transported material frame using `D_roll`.

New per-species parameters:

```text
hydrodynamic_length
hydrodynamic_diameter
solvent_viscosity
crowding_parallel
crowding_perp
crowding_tumble
crowding_roll
```

Derived `D` and mobility values must be written to HDF5. Direct `D` overrides
may remain for diagnostics but must be marked as overrides.

Unit tests:

1. With no force, measured parallel/perpendicular MSD slopes match `2D`.
2. Orientation correlation decays at the theoretical rigid-rod rate.
3. With a constant parallel or perpendicular force, drift/noise satisfy the
   Einstein relation independently in each channel.
4. Rotating the complete initial system leaves all measured statistics
   invariant.
5. Axial material-frame roll does not spuriously change the rod axis.

## Stage 3: calibrate free-solution and intracellular mobility

Implement a small reference calculator that takes `L`, `d_h`, `T`, and `eta`
and reports all four diffusion channels. Validate it against published dilute
F-actin diffusion data before using it for cardiac filaments.

Represent intracellular effects using dimensionless drag multipliers:

```text
D_parallel,eff = D_parallel,water / chi_parallel
D_perp,eff     = D_perp,water     / chi_perp
D_tumble,eff   = D_tumble,water   / chi_tumble
D_roll,eff     = D_roll,water     / chi_roll
```

Do not assume the same `chi` for translation and rotation. A long filament can
be rotationally caged while retaining appreciable axial sliding.

Initial uncertainty ranges, used as priors rather than conclusions:

```text
chi_parallel = 3--30
chi_perp     = 5--100
chi_tumble   = 10--300
chi_roll     = 3--100
```

Run sensitivity sweeps in log space. Report results versus the actual derived
relaxation times, not only versus `chi`.

Experimental route for narrowing these ranges:

1. Reconstituted thin and thick filaments in buffer: sparse fluorescence or
   SHG tracking to measure center and orientation MSDs.
2. Repeat in cytoplasmic extract or a controlled crowding series.
3. In assembling iPSC cardiomyocytes, sparsely label nascent filaments and fit
   short-lag parallel/perpendicular displacement and orientation correlations
   before stable incorporation.
4. Treat FRAP of ACTN2 or myosin molecules as turnover data, not filament
   diffusion data.
5. Fit drag factors hierarchically, retaining cell-to-cell and assembly-stage
   variation.

Acceptance: the dilute calculation reproduces in-vitro F-actin measurements,
and every intracellular value has a documented drag factor and uncertainty.

## Stage 4: enforce fluctuation-dissipation and numerical convergence

The same mobility tensor must multiply deterministic force and set Brownian
noise covariance. Never tune force mobility and noise independently in a run
described as thermal.

1. Audit every translational and rotational update for `M=beta D`.
2. Apply displacement caps only as emergency numerical guards. Record every
   cap event; accepted runs require a negligible cap fraction.
3. Perform timestep convergence at `dt`, `dt/2`, and `dt/4` for:
   - free-rod MSD and orientation correlation;
   - harmonic AM radial confinement;
   - one AA-linked pair;
   - the perfect six-actin lattice.
4. Verify the Boltzmann distribution in passive harmonic tests, including
   equipartition and covariance predicted from the Hessian.
5. Confirm that changing mobility while holding the energy model fixed changes
   relaxation time but not the passive equilibrium distribution.

Acceptance: passive equilibrium statistics agree across timesteps and mobility
rescaling. A lattice cannot be declared stable merely because `D` is tiny.

## Stage 5: mechanically couple crosslinked filaments without changing D

Keep each rod's intrinsic mobility tensor unchanged by AA or AM status. A
crosslink produces equal-and-opposite forces and torques; correlated movement
then emerges from the mechanical connection.

For compliant ACTN2 and AM attachments, spring coupling is the primary model:

```text
v_i = M_i F_i
v_j = M_j F_j
F_i = -F_j = F_link
```

If a link is intentionally treated as a rigid constraint, project deterministic
and stochastic increments through a constraint solve using the joint mobility
matrix. Do not set one particle's `D` equal to the other's and do not set either
to zero.

Hydrodynamic interactions between separate rods are a later extension. They
are not required for the first physically consistent release, but the API
should allow off-diagonal mobility blocks in the future.

Repair AA reassociation semantics before interpreting persistence:

1. Fix the currently unreachable recovery/cooldown guard.
2. Separate `aa_refractory_until` from `aa_association_memory_until`.
3. Define whether the coarse AA edge represents one ACTN2 dimer or an ensemble.
4. If it represents an ensemble, use an occupancy or birth-death variable rather
   than relying on instantaneous binary rebinding to create macro-persistence.
5. Record molecular dwell, off-gap, edge occupancy, and macro-association
   lifetime separately.
6. Require newly bound edges to begin in status 1; only transmitted favorable
   load may establish status 2.

Acceptance: bond persistence remains finite after the cooldown fix, and
crosslinking changes correlated motion through force transmission rather than
through a hidden diffusion switch.

## Stage 6: make the perfect lattice a mechanics acceptance test

Define two related reference geometries.

1. **Unstrained geometric hexagon:** six parallel actins equally spaced around
   one myosin. Adjacent AA distance equals AM radius, so exact simultaneous AM
   and AA rest-length compatibility requires `rho=1`.
2. **Mechanically equilibrated hexagon:** for `rho != 1`, minimize the symmetric
   AM+AA energy first. With six equivalent AM and six adjacent AA springs, a
   simple radial estimate is

   ```text
   R_eq = (k_AM d_AM + k_AA d_AA) / (k_AM + k_AA),
   ```

   before angular and finite-overlap corrections. The initializer should use
   the actual numerical minimizer, not assume `R=d_AM`.

Test layers:

### 6A. Force and Hessian audit

- Fix the intended six-edge AA topology.
- Disable KMC, motors, and noise.
- Verify net force and torque are zero at the minimized geometry.
- Numerically differentiate the energy/force and compute the transverse
  Hessian.
- Apart from global translation and rotation modes, all eigenvalues must be
  non-negative and the lattice modes intended to be stable must be positive.
- Check every AM and AA force with finite differences of the stored energy.

### 6B. Passive thermal stability

- Use `T=310 K` and calibrated mobility tensors.
- Keep intended AA edges persistently occupied.
- Require no secular radius drift, no topology change, and stationary
  distributions consistent with the Hessian.
- Initial target: tail `|psi6| >= 0.9`, angular-gap CV `<=0.15`, and bounded
  radial variance for at least 100 slowest-mode relaxation times.

### 6C. KMC turnover stability

- Enable AA binding/unbinding with the corrected recovery logic.
- Temporary molecular exchange may occur, but the six-neighbor lattice should
  remain structurally stable.
- Limit AA valency or use geometry/energy-dependent association so an arbitrary
  diagonal chord cannot replace a correct nearest-neighbor edge merely because
  it entered the capture radius.
- Check detailed balance in the passive constant-lifetime control before adding
  catch stabilization.

### 6D. Active-cycle stability

- Enable explicit heads, physiological beta, and a calcium cycle.
- Permit reversible periodic strain, but require recovery in diastole and no
  cycle-to-cycle loss of `psi6`.
- Audit axial and transverse AM force, myosin torque, material-frame roll, and
  AA strain to identify the term that destabilizes the ring.
- A perfect lattice must survive with synergy off and on. If synergy destroys
  it, fix motor force geometry rather than reducing diffusion.

Acceptance: the model maintains a mechanically equilibrated perfect lattice at
physiological temperature before it is used for assembly claims.

## Stage 7: demonstrate annealing and then return to production assembly

Stability does not prove annealing. Starting from matched perturbations, test:

```text
angular noise: 5, 10, 20, 30 degrees
radial noise: 0.02, 0.05, 0.10 times R_eq
missing AA edges: 0, 1, 2
incorrect diagonal edge: absent/present
rho: 0.8, 0.9, 0.92, 1.0
```

Use at least 10 matched replicates for the final comparison. Measure:

- first-passage time to `|psi6| >= 0.9`;
- probability of reaching and remaining in the target basin;
- angular-gap CV and radial CV;
- correct-neighbor precision and recall;
- AA molecular and macro-association persistence;
- AM occupancy, force balance, and material-frame torque;
- recovery after one deliberately removed or incorrect AA edge;
- diastolic convergence and absence of latching.

Required controls:

1. AA off;
2. AM active force off but passive tether retained;
3. constant AA lifetime matched for occupancy and turnover;
4. catch stabilization off;
5. thin feedback only, thick feedback only, combined, and neither;
6. physiological mobility versus uniformly rescaled mobility;
7. athermal mechanics, clearly labeled diagnostic only.

Acceptance criteria:

- perfect-lattice stability passes first;
- perturbed states recover more often and faster than AA-off controls;
- final passive order is insensitive to uniform mobility rescaling after enough
  relaxation time;
- the preferred `rho` is supported by stability and annealing, not merely by
  increased bond count;
- full 2000-actin/1000-myosin simulations reproduce the local coordination
  result without new force caps or KMC inconsistencies.

Only after these criteria pass should the calibrated mobility model replace the
current production defaults.

## Implementation order and deliverables

1. Unit/temperature header and HDF5 metadata.
2. Standalone finite-cylinder mobility calculator plus reference tests.
3. Tensor translational integrator.
4. Tumbling and material-frame roll integrator.
5. Fluctuation-dissipation and timestep tests.
6. AA cooldown fix and explicit molecular-versus-edge persistence diagnostics.
7. Force-gradient/Hessian lattice audit.
8. Passive perfect-lattice stability.
9. KMC-turnover and active-cycle stability.
10. Perturbation/annealing factorial.
11. Calibrated full assembly runs.

## Implementation status (2026-07-20)

Implemented in the active C++ build:

- exact `kBT`/`beta(T)` conversion, declared units, HDF5 thermal/mobility
  metadata, and an explicit-beta diagnostic override warning;
- a standalone finite-cylinder calculator and per-species
  `D_parallel`, `D_perp`, `D_tumble`, and `D_roll` with independent drag
  multipliers;
- tensor-consistent deterministic drift and Brownian noise, tangent-plane
  tumbling, sphere exponential-map updates, and transported-frame axial roll;
- cumulative force/torque cap counters in every trajectory;
- corrected AA refractory semantics, separation of capture from occupied-bond
  persistence, pair-order-symmetric axial registry, and optional diagnostic AA
  valency;
- configurable transition-state splitting of AA elastic energy,
  `k_on = k_on0 exp(-beta theta DeltaE_link)`, plus a passive-control mode with
  `k_off = k_off0 exp(+beta (1-theta) DeltaE_link)` and exact finite-step event
  probabilities. The passive mode preserves the local rate ratio while the
  active catch-bond kinetics remain a separate model;
- physical radial spring units (`pN/um`), removal of the passive AM axial
  endstop, an equilibrated six-ring initializer, optional intended-edge
  prebinding, and a true isometric axial constraint for controls;
- force-symmetry, rest-force, transverse Hessian, mobility covariance,
  Einstein-drift, refractory, and existing motor/KMC regression tests;
- lattice analysis for topology precision/recall, molecular versus
  gap-bridged macro occupancy, first passage, basin occupancy, and cap rates.

Current validation conclusion: the mechanical Hessian is stable and a perfect
prebound ring can be maintained at 310 K without numerical caps. With
`k_AM=k_AA=3000 pN/um`, `kappa_AA=50 pN um`, and `dt=1e-6 s`, three matched
runs gave tail `|psi6|=0.9562 +/- 0.0049`, angular-gap CV `0.0563`, complete
intended topology, and zero force/torque-cap fraction. The old `k=10` value is
far too soft when interpreted as `pN/um`; its thermal bond fluctuation is about
20.7 nm at 310 K.

The first perturbed-ring controls also exposed a real kinetic defect rather
than a mobility defect: uniform association lets strained diagonal AA chords
bind and long-lived status-1 links then trap them. Fast uniform turnover
improves motion but also erodes a perfect ring. The energy-gated association
law above is implemented and unit-tested as the detailed-balance-compatible
control. A matched three-replicate passive sweep over
`theta={0,0.1,0.25,0.5,1}` found a kinetic tradeoff: `theta=0` annealed most
strongly (`|psi6|=0.924+/-0.045`) but retained a perfect ring at only
`0.870+/-0.076`, whereas `theta=1` retained the perfect ring
(`0.952+/-0.012`) but annealed poorly (`0.790+/-0.136`). No tested theta passed
both criteria over 0.2 s. This indicates that barrier splitting alone does not
solve structural persistence and annealing simultaneously; binding free energy,
capture-state mechanics, or a multistep ACTN2 encounter complex must be tested.
The full perturbation factorial in Stage 7 therefore remains unaccepted.
`--aa_energy_gated_association` therefore remains an
explicit opt-in until that validation is complete. Low-load loss of AM overlap
is retained as a control and must not be hidden by changing an intrinsic
diffusion coefficient.

### Provisional AA kinetic choice and deferred two-arm model

For current opt-in simulations, use `aa_transition_state_theta=0.1` as a
provisional finite-time compromise. In the three-replicate, 0.2 s passive
six-ring test it gave:

```text
perturbed/annealing tail |psi6| = 0.873 +/- 0.105
initially perfect tail |psi6|   = 0.896 +/- 0.073
```

This is not an equilibrium calibration or evidence that `theta=0.1` is a
molecular constant. It is a working kinetic choice and energy-gated AA remains
off by default.

TODO after the current assembly study: replace the binary AA chemical state by
an explicit two-arm ACTN2 scheme,

```text
unbound <-> singly bound to one actin <-> doubly bound AA bridge.
```

Only the doubly bound state should exert interfilament AA mechanics. The free
arm of a singly bound ACTN2 should search locally and bind its second actin with
a rate determined by connector reach, axial overlap, and binding-site
orientation. A doubly bound bridge may be mechanically status 1 or status 2;
chemical occupancy and mechanical load status must remain separate. Loss of
one arm should return the molecule to the singly bound state, enabling rapid
local rebinding and macro-crosslink persistence without inventing a weak
interfilament spring. Parameterize first-arm kinetics from isolated ABD data,
second-arm capture from effective local concentration/reach, and bridge
detachment from passive or catch-bond measurements. Before production use,
test passive thermodynamic consistency, perfect-lattice persistence, perturbed
annealing, and active motor-driven assembly against the binary-theta control.

All fourteen current unit/regression tests pass, including empirical tensor-noise
covariance, Einstein drift, AA refractory timing, association energetics,
pair-order symmetry, and the six-ring force/Hessian audit.

Each stage must produce a machine-readable CSV, a compact plot, and a pass/fail
summary. Failed stages block later biological interpretation but do not block
diagnostic work on the cause.

## TODO: calibrate effective filament mobility

Status: unresolved; the finite-cylinder hydrodynamic formulas are principled,
but the environmental drag multipliers used by the production defaults are
still phenomenological:

```text
actin translation:   10x drag
actin rotation:     160x drag
myosin translation:  20x drag
myosin rotation:      64x drag
```

Do not select four independent multipliers by maximizing assembly order. Keep
the parallel/perpendicular mobility ratio fixed by finite-cylinder
hydrodynamics. Calibrate at most one translational slowdown per filament
species from trajectories of unbound filaments using

\[
\langle \Delta r_\parallel^2(t)\rangle=2D_\parallel t,
\qquad
\langle \Delta r_\perp^2(t)\rangle=4D_\perp t.
\]

Calibrate one rotational slowdown per species from the second-Legendre
orientation autocorrelation,

\[
\left\langle
P_2\!\left[\hat{\mathbf u}(t)\cdot\hat{\mathbf u}(0)\right]
\right\rangle
\sim e^{-6D_{\mathrm{tumble}}t}.
\]

Use experimentally measured filament trajectories when available. In their
absence, report an explicit sensitivity sweep over a shared effective-viscosity
or slowdown scale rather than presenting one fitted multiplier set as a
physical calibration. Separate localization error and short-time tracking
noise from the diffusive fitting window.

The calibrated values are intrinsic single-filament mobilities. Do not change
them when AM or AA bonds form. Crosslinks, motors, sterics, and other mechanical
connections must reduce motion by coupling forces and producing collective
drag, not by switching a filament's diffusion coefficient.

Acceptance deliverables:

1. Unbound-actin and unbound-myosin trajectory ensembles with no interactions.
2. Body-frame parallel/perpendicular MSD curves and identified linear fitting
   windows.
3. Orientation autocorrelation curves and fitted tumbling diffusion constants.
4. Recovery of the imposed coefficients in synthetic-data and simulator-only
   controls.
5. A parameter file containing the calibrated slowdowns, confidence intervals,
   provenance, temperature, viscosity convention, filament dimensions, and
   fitting window.
6. A production sensitivity check showing which assembly conclusions survive
   the mobility uncertainty.

## Primary references

- Tirado MM, Garcia de la Torre J. Translational friction coefficients of rigid,
  symmetric top macromolecules. J Chem Phys 71, 2581 (1979).
  https://doi.org/10.1063/1.438613
- Tirado MM, Garcia de la Torre J. Rotational dynamics of rigid, symmetric top
  macromolecules. J Chem Phys 73, 1986 (1980).
  https://doi.org/10.1063/1.440288
- Tirado MM, Martinez CL, Garcia de la Torre J. Comparison of theories for the
  translational and rotational diffusion coefficients of rod-like
  macromolecules. J Chem Phys 81, 2047 (1984).
  https://doi.org/10.1063/1.447827
- Hydrodynamic and Polyelectrolyte Properties of Actin Filaments: Theory and
  Experiments. Polymers 14, 2438 (2022).
  https://pmc.ncbi.nlm.nih.gov/articles/PMC9230757/
- Li G, Tang JX. Diffusion of actin filaments within a thin layer between two
  walls. Phys Rev E 69, 061921 (2004).
  https://doi.org/10.1103/PhysRevE.69.061921
- Myofibril Assembly Visualized by Imaging N-RAP, Alpha-Actinin, and Actin in
  Living Cardiomyocytes. Exp Cell Res 315, 2126 (2009).
  https://pmc.ncbi.nlm.nih.gov/articles/PMC2742992/
- Visualization of cardiac thick filament dynamics in ex vivo heart
  preparations. J Mol Cell Cardiol 185, 88 (2023).
  https://doi.org/10.1016/j.yjmcc.2023.10.013
