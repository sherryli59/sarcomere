# A–A / A–M overlap mechanism diagnostic

Date: 2026-07-30  
Scope: sandbox and short diagnostics only; **no production run was launched**

## Bottom line

The present fixed-periodic setup cannot answer whether the system selects a
finite mean A–M overlap.  In the exact `two_sarcomere_unit` geometry, mean
Z-disc--M-line spacing is constrained by the prescribed axial box length.
Changing the starting "A–M overlap" by changing `Lx` creates three different
immutable system lengths; no internal particle interaction can make those
lengths converge.

The tests do answer several narrower questions:

1. Turnover-capable elastic ACTN2 links anneal their rest lengths and creep.
   They can slow sliding, but they do not generate a persistent preferred
   A–A overlap by themselves.
2. A Johann-like overlap adhesion or fixed-pool entropic pressure can oppose
   motor-driven loss of A–A overlap, but neither supplies the missing global
   A–M restoring direction in a fixed cell.
3. A titin-like Z--M spring controls local A–M spacing, but it does not select
   A–A overlap.
4. A persistent A–A register plus titin stabilizes both coordinates in the
   reduced free-length test.  This is a positive control, not an emergent
   explanation: both target lengths are encoded.
5. A fixed-pool ACTN2 pressure plus local titin is **not automatically**
   equivalent to the Johann mechanism.  In the force-balanced network, titin
   can cancel the motor on each A–M pair, leaving the one-sided crosslinker
   pressure unopposed in the A–A coordinate.
6. A shared series elastic load is mechanically different: it stabilizes the
   global length while preserving local motor--crosslinker competition.
   However, the tested entropic-plus-series system still has an unstable
   A–A/A–M mode and fails the full stochastic convergence test.

The minimal defensible conclusion is therefore: **the current ingredients do
not yet demonstrate emergent selection of both overlaps.**  The missing test is
a variable-length cell or open-ended chain run under one common external load.

## What was implemented

Everything is isolated under this directory.  Existing C++ source files were
not changed.

### 1. Local four-body stochastic model

`overlap_model.py` contains two antiparallel actins meeting at one Z-disc and
two neighboring myosin centers.  It includes reciprocal motor, ACTN2, titin,
register, and optional M-line forces plus translational Brownian diffusion.

This model is useful for force signs and creep, but its local topology omits the
second Z-disc needed to close the periodic force-balanced network.

### 2. Exact fixed-periodic topology

`periodic_model.py` reproduces the axial connectivity of C++
`two_sarcomere_unit`:

- four actins;
- two bipolar myosins;
- two A–A pairs;
- every myosin engaged from both sides;
- minimum-image periodic geometry;
- reciprocal force accounting for every interaction.

### 3. Dynamic-repeat topology

`variable_cell_model.py` adds one overdamped axial repeat-length coordinate.
The cell responds to the pairwise minimum-image axial virial, an optional
applied tension, and an optional series spring.  Axial coordinates deform
affinely when the cell length changes.

This is a diagnostic barostat-like model, not a production implementation.  It
exists to determine whether a candidate interaction supplies restoring force
in the global-length direction that a fixed cell removes.

## Candidate interactions

| Label | Physical interpretation | What it can constrain |
|---|---|---|
| `baseline_turnover` | Elastic ACTN2 unbinds and rebinds stress-free | transient A–A resistance |
| `permanent_actn2` | Non-turning-over elastic links | A–A memory, not annealing |
| `johann_adhesion` | `U(y) = -gamma y` | constant force favoring larger A–A overlap |
| `johann_entropic` | `F(y) = N kBT / y` | overlap-dependent A–A pressure |
| `titin` | Z-disc--M-line spring | local A–M/Z--M spacing |
| `register_plus_titin` | harmonic A–A target plus titin | prescribed A–A and A–M targets |
| `entropic_plus_titin` | fixed-pool pressure plus local titin | two passive terms, but not necessarily two independent restoring modes |
| `entropic_plus_series` | fixed-pool pressure plus global series elasticity | motor–crosslinker competition plus global length control |

The Johann fixed-pool force assumes longitudinally mobile, end-confined
crosslinkers.  The production ACTN2 implementation consists of discrete
turnover-capable elastic links and therefore does not currently implement that
free energy.

## Verification

### C++ regression baseline

A clean out-of-tree build was made under `/tmp`, and all **22/22 C++ tests
passed**.  A short 20,000-step baseline smoke trajectory was also completed.
The longer C++ HDF5 matrix was stopped because output overhead dominated and it
could not resolve global-length selection in a fixed box.

### Sandbox unit tests

All **15/15 Python tests passed**, including:

- exact initial A–A symmetry and equal A–M overlaps;
- zero net force for every reciprocal internal interaction;
- correct motor and titin cell-force signs;
- minimum-image invariance;
- affine cell-coordinate preservation;
- deterministic stochastic trajectories for a fixed seed;
- turnover erasing stored ACTN2 strain;
- correct motor, adhesion, and titin force directions.

### Campaign design

The fixed-periodic campaign used:

- seven mechanisms;
- five initial conditions;
- six seeds;
- 4 s per trajectory;
- 210 trajectories total.

The dynamic-repeat campaign used six candidate/control mechanisms, the same
five initial conditions, and three seeds (90 trajectories).  Initial
`Lx = 3.9, 4.3, 4.7 um` states were evolved under the **same** external
boundary condition.  Trajectories lasted 12 s, or 15 s for the slower series
spring.  Terminal means, initial-condition spreads, and terminal slopes were
recorded rather than judging stability from a final frame.

## Fixed-periodic results

| Mechanism | mean A–A (um) | A–A IC spread | mean A–M (um) | A–M IC spread | max terminal speed (um/s) |
|---|---:|---:|---:|---:|---:|
| baseline turnover | -0.239 | 0.804 | 0.602 | 0.142 | 0.309 |
| permanent ACTN2 | 0.091 | 0.122 | 0.599 | 0.150 | 0.022 |
| Johann adhesion | 0.043 | 0.251 | 0.586 | 0.160 | 0.058 |
| Johann entropic | 0.093 | **0.020** | 0.593 | 0.162 | 0.030 |
| titin | 0.066 | 0.239 | 0.639 | 0.065 | 0.041 |
| entropic + titin | 0.087 | 0.233 | 0.628 | 0.069 | 0.035 |
| register + titin | 0.088 | 0.125 | 0.624 | 0.104 | **0.009** |

No mechanism passes the nominal two-coordinate convergence rule in a fixed
cell.  The important exception is interpretive: the A–M initial-condition
spread cannot go to zero because the "deep" and "shallow" groups have different
fixed box lengths.  Those failures must not be used as parameter-rejection
evidence.

The fixed-pool entropic interaction is the only non-prescriptive candidate
that narrows the A–A initial-condition spread below `0.025 um`.  It still
drifts and does not control global A–M.

## Dynamic-repeat results

The complete table is generated in `results/variable_cell_summary.csv`.

| Mechanism | mean A–A (um) | A–A IC spread | mean A–M (um) | A–M IC spread | `Lx` IC spread (um) | max terminal speed (um/s) | pass |
|---|---:|---:|---:|---:|---:|---:|---|
| baseline turnover | -0.901 | 0.906 | 0.169 | 0.321 | 0.825 | 0.335 | no |
| Johann entropic | 0.097 | 0.089 | 0.572 | 0.461 | 0.885 | 0.074 | no |
| titin | 0.087 | 0.103 | 0.649 | **0.0003** | 0.206 | 0.0009 | no |
| entropic + titin | 0.223 | 0.035 | 0.643 | **0.0012** | 0.075 | 0.0246 | no |
| register + titin | **0.101** | **0.0006** | **0.649** | **0.0003** | **0.0021** | **0.0007** | **yes** |
| entropic + series | 0.160 | 0.030 | 0.264 | 0.145 | **0.0090** | 0.0320 | no |

Titin makes the A–M overlap nearly independent of initial conditions, but the
A–A coordinate remains free.  In the entropic-plus-titin run, A–A overlap
keeps increasing while the repeat contracts; these two motions compensate so
that the titin-controlled A–M overlap looks stationary.  A stationary A–M trace
therefore does not imply a stationary sarcomere.

The strong series spring does the complementary job: repeat length converges,
but A–A continues to increase and A–M continues to decrease.  It has a
candidate fixed point in the exactly symmetric instantaneous phase portrait,
but the full stochastic six-body trajectory leaves that manifold.  This is why
testing the complete force-balanced topology is essential.

The dynamic-repeat model is a topology/force-rank diagnostic.  Passing it is
necessary evidence that a mechanism can restore both coordinates; it is not
sufficient evidence that the parameter values or coarse-grained cell mobility
are biological.

## Visual diagnostics

- `results/periodic_timeseries.png`: fixed-periodic mechanism comparison.
- `results/variable_cell_timeseries.png`: A–A, A–M, and repeat-length traces.
- `results/variable_cell_length_convergence.png`: long/short repeat comparison
  for entropic+titin versus register+titin.
- `results/variable_cell_phase_portraits.png`: instantaneous two-coordinate
  vector fields and nullclines.
- `results/periodic_*.gif`: exact fixed-periodic topology animations.
- `results/variable_register_plus_titin_from_4.7.gif`: successful positive
  control contracting from a long repeat.
- `results/variable_entropic_plus_series_from_4.7.gif`: series-loaded failure
  mode.

These are rod/topology animations.  They do **not** show explicit weak/strong
myosin heads, head strain, stereospecific binding sites, or rotational
diffusion; those variables do not exist in the reduced sandbox.  A clean
head-state animation remains a separate unresolved deliverable for the full
simulator.

## Why the Johann papers overlap with this model—but do not answer it directly

Johann et al. demonstrate stable overlaps when active motors and passive
crosslinker forces act on the same overlap coordinate.  That supports the
general idea that turnover-capable passive crosslinkers can participate in an
active steady state without becoming permanent tethers.

The sarcomere question has at least two independent axial coordinates:

```text
A--A registry y
A--M / Z--M spacing x
```

An interaction that produces curvature or a restoring force only in `y`
leaves `x` free.  Conversely, titin can control `x` while leaving `y` free.
Furthermore, attaching titin and motors to the same actin--myosin pair changes
the force paths: once those forces balance locally, the motor may no longer be
available to oppose an ACTN2 pressure at the Z-disc.  The microtubule overlap
mechanism therefore cannot be transferred by force-sign analogy alone; the
network connectivity and rank of the restoring-force matrix matter.

There is also a crucial sign difference.  In the symmetric periodic geometry,

```text
Z--M distance = Lx / 4 + y / 2 .
```

At fixed `Lx`, increasing A–A overlap `y` increases Z--M distance and generally
*decreases* A–M overlap and the number of engaged myosin heads.  Thus the motor
force opposing ACTN2 pressure becomes weaker as `y` grows.  Johann's stable
microtubule fixed point instead relies on motor force growing with the very
same overlap coordinate that the passive crosslinkers favor.  The naive
mapping therefore has the wrong feedback sign in this sarcomere geometry.

Chew et al. reinforce the broader point that bipolar organization depends on
opposing active processes and boundary/filament architecture.  That is useful
conceptually, but their opposing-motor system does not establish that
myosin + ordinary turnover elastic ACTN2 alone selects both sarcomeric
coordinates.

## Production-code implications

The existing C++ implementation has two relevant limitations for this
question:

1. `two_sarcomere_unit` fixes `box[0]`, so mean sarcomere length is imposed.
2. `titin_rest_length` is computed in `main.cpp` as
   `0.5 * myosin_length + (2/3) * actin_length` rather than exposed as a
   command-line parameter.  Titin application is also conditional on current
   A–A and A–M bond state, which is not identical to the always-connected
   reduced spring used as a clean control here.

The recommended opt-in production design and validation matrix are in
`VARIABLE_LENGTH_TEST_PLAN.md`.  The fixed-cell behavior must remain the
default and should be bitwise unchanged when the new mode is disabled.

## Remaining unresolved questions

1. Is ACTN2 longitudinally mobile/end-confined enough for a Johann-like
   entropic pressure to be biologically meaningful, or is a barbed-end
   register/binding-zone potential more appropriate?
2. What physical structure supplies global series elasticity in this
   preparation: titin, myomesin/M-line compliance, Z-disc compliance,
   neighboring sarcomeres, or the external boundary?
3. Does the same mechanism remain stable with explicit weak/strong myosin
   heads, load-dependent kinetics, rotational diffusion, and torque?
4. Is the observed fixed point stable to small perturbations on both sides,
   or merely slow because bonds turn over?
5. Does it persist after time-step, cell-mobility, crosslinker-pool, and
   system-size convergence tests?
6. Can one parameter set reproduce both free-actin sliding/unbinding and
   force-balanced partial-overlap stalling?
7. Can cross-midpoint/other-side myosin capture supply the missing
   overlap-increasing motor feedback without producing an unphysical
   discontinuity?  This mechanism remains deferred in the production model.
8. How should weak/strong heads, head strain, binding sites, filament rotation,
   torque, and rotational diffusion be rendered without hiding the axial
   overlap dynamics?

## Decision

Do not launch a production sweep yet.  The simple Johann mapping is rejected
for the current network, not merely under-tuned.  The next full-simulator
diagnostic should first add a sandboxed/opt-in variable-length or open-boundary
degree of freedom and expose titin rest length.  The two most defensible next
mechanisms are:

1. a physical barbed-end/Z-disc register interaction plus Z--M elasticity
   (the positive control that passed here); and
2. a separately justified opposing-side/cross-midpoint motor feedback that
   makes active force grow in the A–A restoring direction, tested under a
   shared compliant boundary.

Then repeat the common-load convergence test from both shorter and longer
initial lengths before any production-scale sweep.
