# From the legacy `am_kmc` sweep model to the current sarcomere model

Date: 2026-07-27

## Scope and provenance

This document summarizes changes between:

- the legacy three-dimensional `am_kmc` model used for the large `data0430`
  and related actin-number/catch-lifetime parameter sweeps; and
- the current development model, whose defaults use explicit spatial
  myosin heads, sparse local ACTN2 links, productive-register status-2
  kinetics, physical rod mobility, and full filament sterics.

The legacy reference is git commit `30829b0`, previously reconstructed and
rebuilt for the sweep audit. The present reference is the working tree as of
2026-07-27. Because the current tree contains substantial uncommitted
development, current results should record the executable SHA-256 in addition
to the git commit.

The central conclusion is that these are not two parameterizations of the
same model. The legacy model was a coarse-grained assembly machine with
continuous geometric AM occupancy, an always-available nematic AM torque, and
prescribed actin sliding. The current model attempts to obtain those behaviors
from individual cross-bridge kinetics, elastic forces, force balance, and
thermal rod motion.

## Executive comparison

| Mechanism | Legacy sweep model | Current default |
|---|---|---|
| AM state | Every eligible AM pair is immediately and continuously bonded | Conserved heads cycle through OFF/SRX, detached, weak, and strong states |
| Active motion | Prescribed actin velocity `v_am * actin.direction` | Motion follows forces through overdamped rod mobility |
| AM force | Pair-level radial spring plus direct nematic potential | Per-head radial tether plus weak pair compliance and a finite active axial stroke |
| Force-velocity law | Overlap proxy throttles the prescribed velocity | Strong-head stroke speed decreases with its tensile force |
| AM torque | Direct pair energy proportional to `kappa_am` | Off-center head forces plus a weaker occupied-pair flat-bottom compliance |
| Thick-filament geometry | Coarse pair overlap | Two head-bearing halves, explicit central bare zone, material crowns and head azimuths |
| AM occupancy capacity | Up to a configured number of actins per myosin | One conserved pool per half; per-neighbor capacity scales with physical half-overlap; material actin sites and a global actin cap are shared across all myosins |
| Heartbeat | Absent; tonic motor activity | Optional finite systole/diastole, thin activation and strong-head feedback |
| AA representation | At most one implicit bond per actin pair | Several sparse implicit ACTN2 molecules with fixed material attachment coordinates |
| Status 1 | Bound pair; load can extend its lifetime in the legacy law | Ordinary ACTN2 link; no load-dependent stabilization by default |
| Status 2 | Geometry plus AM-overlap proxies and distinct myosin identities | Existing, registered antiparallel link established by opposed productive strong-head force |
| Status-2 retention | Reclassified from the AM proxy every step | Current `productive_register` retains status 2 through transient motor loss while registry persists |
| AA axial organization | Strong asymmetric pair-level barbed-end potential acts on any bond | Local status-1 links bind anywhere; bounded barbed-end registry acts only in status 2 |
| Mobility | One translational and one rotational scalar per species | Parallel, perpendicular, tumbling and rolling mobilities from finite-cylinder hydrodynamics |
| Sterics | Limited and bond-dependent; pre-phase was principally MM relaxation | Independent AA, AM and MM rod exclusion; full steric pre-relaxation |
| Randomness | Thread-owned mutable GSL streams | Object/timestep-keyed KMC streams; scheduling does not change an object's random draw |

## 1. Legacy actin-myosin mechanics

### 1.1 Geometric occupancy

For every neighboring actin \(i\) and myosin \(m\), the legacy code computed
coarse overlap fractions. If the pair had nonzero binding and productive
partial overlap, it immediately set

\[
B_{im}=1.
\]

There was no detached/weak/strong cross-bridge cycle. A pair remained occupied
as long as the geometric test succeeded. Consequently the effective AM duty
ratio was close to one for every eligible pair.

### 1.2 Pair-level AM energy

Let \(d_{im}\) be the shortest distance between the finite rods, \(d_{\rm
AM}^*\) the preferred distance, and \(\mathbf u_i,\mathbf u_m\) their axes.
Inside the cutoff, the legacy energy was

\[
U_{\rm AM}^{\rm legacy}
=
\frac{k_{\rm AM}}{2}
\left|\mathbf u_i\!\cdot\!\mathbf u_m\right|
\left(\frac{d_{im}-d_{\rm AM}^*}{d_{\rm AM}^*}\right)^2
+
\frac{\kappa_{\rm AM}}{2}
\left[1-(\mathbf u_i\!\cdot\!\mathbf u_m)^2\right].
\]

Outside the radial cutoff, the angular term remained active for a pair that
reached this force path. Thus every eligible pair supplied a strong direct
nematic organizer. With the historical \(\kappa_{\rm AM}=100\), this was a
major source of actin-myosin alignment.

The legacy code also included an axial end-stop based on an actin tip and a
hard-coded \(0.4L_m\) scale. It could generate propulsion-like axial force
independently of the prescribed motor velocity. This has been removed rather
than interpreted as a physical cross-bridge force.

### 1.3 Prescribed motor sliding

The legacy motor did not compute velocity from an elastic head force. It first
assigned

\[
\mathbf v_i^{(0)}=v_{\rm AM}\mathbf u_i.
\]

For a status-2 actin, a bounded overlap proxy \(f_i\in[0,1]\) reduced this to

\[
\mathbf v_i=(1-f_i)\,v_{\rm AM}\mathbf u_i.
\]

Reaction drift was assigned to associated myosin filaments. This imposed the
direction and scale of relative sliding even when the mechanical forces were
unbalanced. It therefore acted as an active sorting rule, not an overdamped
cross-bridge model.

### 1.4 Legacy AM “load”

For each AM neighbor, define a normalized productive-overlap contribution

\[
c_{im}=3\min(g_{im},1/3),
\qquad
a_{im}=|\mathbf u_i\!\cdot\!\mathbf u_m|,
\]

where \(g_{im}\) is the partial binding ratio. The actin-level load proxy was
approximately

\[
f_i^{\rm proxy}
=
\frac{1-\exp[-2\max_m(c_{im}a_{im})]}{1-\exp(-2)}.
\]

This quantity is dimensionless. It measures favorable overlap and alignment;
it is not tensile force exerted by a myosin head.

## 2. Legacy ACTN2 model

### 2.1 One edge per actin pair

The state was a dense pair matrix. Each actin pair could have at most one
implicit AA bond and one status value.

The pair energy combined radial spacing, nematic alignment, and a strong
barbed-end registry:

\[
U_{\rm AA}^{\rm legacy}
=
\frac{k_{\rm AA}}{2}
\left(\frac{d_{ij}-d_{\rm AA}^*}{d_{\rm AA}^*}\right)^2
+
\frac{\kappa_{\rm AA}}{2}
\left[1-(\mathbf u_i\!\cdot\!\mathbf u_j)^2\right]
+
U_{\rm axial}^{\rm legacy}.
\]

For a signed barbed-end offset \(s\), the axial term contained both a target
near \(0.1L_a\) and a tenfold one-sided wall:

\[
U_{\rm axial}^{\rm legacy}
=
\frac{k_{\rm AA}}{2}
\left(\frac{s-0.1L_a}{0.1L_a}\right)^2
+
\begin{cases}
5k_{\rm AA}\left(s/d_{\rm AA}^*\right)^2,&s<0,\\
0,&s\ge 0.
\end{cases}
\]

This potential strongly selected a particular antiparallel axial
configuration. It was also pair-order-sensitive.

### 2.2 Legacy status classifier

An eligible pair was status 1 if it was within AA capture distance. It became
status 2 only when:

1. both actins had nonzero AM alignment/overlap proxies;
2. the actins were antiparallel in directional mode;
3. both had geometrically assigned myosin bonds; and
4. at least one myosin identity on actin \(i\) differed from one on actin
   \(j\).

This did not require strong heads, actual opposed forces, or load transmitted
through ACTN2. A pair of aligned AM overlaps was enough.

### 2.3 Legacy lifetime law

For an existing bond, the code formed

\[
f_{ij}^{\rm legacy}
=
|\mathbf u_i\!\cdot\!\mathbf u_j|
\min(f_i^{\rm proxy},f_j^{\rm proxy}).
\]

Its lifetime was

\[
\tau_{ij}^{\rm legacy}
=
\tau_0
+
\mathbb 1_{s=2}\tau_{\rm directional}
+
\lambda f_{ij}^{\rm legacy},
\]

and the per-step break probability was approximately

\[
P_{\rm off}^{\rm legacy}
=
\min\!\left(1,\frac{\Delta t}{\tau_{ij}^{\rm legacy}}\right).
\]

The important point is that \(\lambda f\) stabilized both status 1 and status
2. Status 1 was therefore not necessarily short-lived in the historical
sweeps.

New association used \(P_{\rm on}\approx k_{\rm on}\Delta t\). The intended
refractory period contained a Boolean error that made its rejection condition
unreachable, so the historical model effectively allowed immediate
reassociation after KMC rupture.

## 3. Current explicit myosin-head model

### 3.1 Conserved, half-filament-owned heads

Each myosin half owns a fixed pool of spatial head slots. A head follows

\[
\mathrm{OFF/SRX}
\rightleftarrows
\mathrm{detached\ ON}
\rightarrow
\mathrm{weak}
\rightarrow
\mathrm{strong}
\rightarrow
\mathrm{detached\ ON},
\]

with a competing weak-to-detached transition. Thick-filament OFF/SRX
recruitment is optional and currently disabled by default.

For a transition with rate \(k\), the current implementation uses

\[
P(\text{event in }\Delta t)=1-e^{-k\Delta t},
\]

and uses competing hazards when multiple exits are possible.

The weak attachment rate to actin \(i\) is schematically

\[
k_{\rm attach}^{h\to i}
=
k_{\rm on}^{\rm AM}\,
g_{\rm overlap}^{h i}\,
\left[b+(1-b)q_i\right],
\]

where \(b\) is the weak-binding floor and \(q_i\) is thin-filament activation.
The polarity-compatible myosin half is the only productive overlap. The central
bare zone has no head anchors. A true opposite-half overlap may bind weakly and
consume a head and a globally shared actin site, but it has zero power-stroke
hazard and cannot become a strong motor.

### 3.2 Overlap-scaled capacity

The 141 default slots on a half-filament are shared among neighboring thin
filaments. For one actin,

\[
N_{\rm cap}
=
\left\lceil
\min(N_{\rm half},N_{\rm actin,max})
\frac{\ell_{\rm productive}}{\ell_{\rm head-bearing}}
\right\rceil,
\]

clamped to the configured full-overlap capacity. The current default is
\(N_{\rm actin,max}=24\). Strong heads reserve capacity before weak heads.

This pair/half capacity is not a new site pool. Each thin filament carries 182
material axial sites, shared by every neighboring myosin, and at most 24 total
weak-plus-strong heads may occupy the actin at once. Parallel association
events atomically claim a free site; detachment releases it immediately.

This replaces the legacy “up to six actin partners per myosin” rule with an
explicit conserved molecular pool and overlap-dependent sharing.

### 3.3 Spatial material heads

A bound head stores:

- its myosin half and crown;
- a transported axial and azimuthal material anchor;
- a fixed actin material coordinate \(s_a\);
- its weak/strong state;
- stroke progress \(q_h\); and
- its binding-time swivel angle.

For actin attachment point \(\mathbf p_a\), myosin anchor \(\mathbf p_m\), and
actin axis \(\mathbf u_a\), define

\[
\mathbf h=\mathbf p_a-\mathbf p_m,\qquad
z=\mathbf h\!\cdot\!\mathbf u_a,\qquad
\boldsymbol\rho=\mathbf h-z\mathbf u_a.
\]

The head uses a two-spring energy

\[
U_h
=
\frac{k_\perp}{2}
\operatorname{flatbottom}(|\boldsymbol\rho|-\rho_0;\delta_\rho)^2
+
U_z.
\]

For a weak head,

\[
U_z^{\rm weak}=\frac{k_{\rm weak}}2z^2,
\qquad
k_{\rm weak}=r_{\rm weak}k_{\rm motor}.
\]

For a strong tension-only head, the productive strain and force are

\[
\epsilon_h=q_h-z,
\qquad
F_h=\min\!\left(k_{\rm motor}\max(\epsilon_h,0),F_{\rm stall}\right).
\]

The default remains tension-only. A bidirectional elastic strong-head control
exists but is not the default.

### 3.4 Finite force-velocity stroke

A strong head performs one finite working stroke rather than walking
processively:

\[
\dot q_h
=
v_0\max\!\left(0,1-\frac{\max(F_h,0)}{F_{\rm stall}}\right),
\qquad
0\le q_h\le d_{\rm ps}.
\]

Thus unloaded stroke motion is fastest, and an isometric load stalls the
remaining internal stroke. Filament velocity itself is not assigned from this
formula; it follows the summed physical forces and rod mobility.

### 3.5 AM detachment

The chemical strong-head off-rate is

\[
k_{\rm strong,off}(F)
=
k_{\rm strong,off}^0
\left[
f_{\rm floor}
+
(1-f_{\rm floor})e^{-\max(F,0)/F_*}
\right].
\]

With the current default \(f_{\rm floor}=1\), this reduces to the
load-insensitive ATP-turnover rate \(k_{\rm strong,off}^0\). Productive load
therefore cannot freeze a head indefinitely.

Nonproductive stored energy accelerates detachment:

\[
k_{\rm off,total}
=
k_{\rm chemical}
\exp\!\left[
\min(\theta\beta E_{\rm nonproductive},x_{\max})
\right].
\]

The implementation also records capacity exits, geometric exits, and an
emergency contour-extension rupture separately.

### 3.6 Force and torque balance

Head forces are action-reaction pairs:

\[
\mathbf F_a^{h}=-\mathbf F_m^{h}.
\]

They act at off-center attachment points, so a thick filament receives

\[
\boldsymbol\tau_m
=
\sum_h(\mathbf p_m^h-\mathbf r_m)\times\mathbf F_m^h.
\]

This is now the active source of myosin alignment. The current model also has
a weaker occupied-pair flat-bottom nematic compliance,

\[
U_{\rm align}
=
\frac{k_{\rm align}}2
\left[
\max\!\left(0,\sin\theta-\sin\theta_{\rm slack}\right)
\right]^2,
\]

but the legacy global \(\kappa_{\rm AM}\) potential is not used when
`dynamic_am=true`.

## 4. Current heartbeat and cooperative activation

For heartbeat period \(T\), the calcium-like drive during systole is a
normalized difference of exponentials,

\[
C(t)\propto
e^{-\phi/\tau_{\rm decay}}
-
e^{-\phi/\tau_{\rm rise}},
\qquad
\phi=t\bmod T.
\]

It is set exactly to the activation floor during the configured diastolic
fraction.

Each actin carries \(q_i\in[0,1]\):

\[
\frac{dq_i}{dt}
=
k_{{\rm open},i}(1-q_i)-k_{\rm close}q_i,
\]

where

\[
k_{{\rm open},i}
=
k_C H(C)
+
k_{\rm strong}H(C)H(n_{{\rm strong},i}).
\]

The calcium factor on both terms prevents strong-head feedback from latching
the thin filament in diastole. This entire activation layer was absent from
the legacy sweeps. A nonpositive heartbeat period gives tonic activation.

## 5. Current spatial ACTN2 model

### 5.1 Multiple implicit molecules

`aa_local_multilink=true` is now the default. Each ACTN2 link stores

\[
(i,j,s_i,s_j,\mathrm{status},t_{\rm attach}),
\]

where \(s_i,s_j\) are fixed material coordinates while the molecule is bound.
ACTN2 is explicit as a sparse connector state, but not as a freely diffusing
particle.

The default transient capacity is two links per actin pair, subject to site
exclusion, actin neighbor valency, optional finite total pool, and retained
lattice-neighbor identity.

### 5.2 Status-1 connector mechanics

For

\[
\mathbf x_i=\mathbf r_i+s_i\mathbf u_i,\qquad
\mathbf x_j=\mathbf r_j+s_j\mathbf u_j,\qquad
\boldsymbol\ell=\operatorname{minimumimage}(\mathbf x_j-\mathbf x_i),
\]

the ordinary local-link energy is

\[
U_1
=
\frac{k_\ell}{2}(|\boldsymbol\ell|-d_{\rm AA}^*)^2
+
\frac{k_\perp}{2}
\left[
(\boldsymbol\ell\!\cdot\!\mathbf u_i)^2
+
(\boldsymbol\ell\!\cdot\!\mathbf u_j)^2
\right]
+
\frac{k_{\rm align}}2
\left[1-(\mathbf u_i\!\cdot\!\mathbf u_j)^2\right].
\]

The first term sets AA spacing; the second favors an approximately
perpendicular ACTN2 connector; separated links plus the weak nematic term
align the rods. All forces and torques are derivatives of the same energy, and
internal AA force sums to zero.

The preferred distances are now separated by molecular geometry:

\[
d_{\rm AM}^*=r_M+r_A+\ell_{\rm head},
\qquad
d_{\rm AA}^*=2r_A+\ell_{\rm ACTN2,projected},
\]

with the current default ratio \(d_{\rm AA}^*/d_{\rm AM}^*=0.92\).

### 5.3 Association

Local association is extensive in accessible axial overlap:

\[
R_{{\rm on},ij}
=
(N_{{\rm cap},ij}-n_{ij})k_{\rm on}^{\rm AA}
\frac{\ell_{{\rm overlap},ij}}{\ell_{\rm ref}}c_{\rm free}.
\]

When transition-state gating is enabled,

\[
k_{\rm on,eff}
=
k_{\rm on}
\exp(-\theta\beta\Delta E_{\rm encounter}).
\]

The encounter energy uses a soft radial stiffness and little or no angular
penalty; the mature link energy acts after association. The default CLI has
energy gating disabled, although production manifests may enable
\(\theta=0.1\).

### 5.4 Productive-register status 2

The current default classifier is `productive_register`. An existing status-1
link becomes status 2 only if:

1. the actins are antiparallel;
2. their barbed-end registry error is below
   \(\delta_{\rm on}\);
3. both actins receive productive strong-head force from the correct myosin
   halves; and
4. the opposed productive support exceeds \(F_{\rm on}\).

The opposed support is

\[
F_{\rm opposed}
=
\chi_{\rm force}
\min\!\left(
|\mathbf F^{\rm productive}_i|,
|\mathbf F^{\rm productive}_j|
\right),
\]

where

\[
\chi_{\rm force}
=
\max\!\left(
0,
-
\frac{\mathbf F_i^{\rm productive}\!\cdot\!
      \mathbf F_j^{\rm productive}}
     {|\mathbf F_i^{\rm productive}|
      |\mathbf F_j^{\rm productive}|}
\right).
\]

In code this is the force-opposition factor times the smaller productive-force
magnitude.

Once established, status 2 persists while the wider registry gate
\(\delta_{\rm off}\) remains satisfied. Transient loss of strong AM support
does not demote it:

\[
s=2
\;\Longrightarrow\;
s'=2
\quad\text{if}\quad
\delta_{\rm registry}\le\delta_{\rm off}.
\]

Therefore `aa_status2_force_off` does not control demotion in the current
`productive_register` mode. This is a major change from both the intermediate
arm-load model and the desired heartbeat-synchronous status-2 diagnostic.

### 5.5 Status-2-only registry mechanics

Status 1 has no collective barbed-end registry force. Status 2 adds a bounded
registry centered at the antiparallel overlap target. For each scalar registry
error \(e\), the default Gaussian-attenuated branch is

\[
U_{\rm reg}(e)
=
\frac12 k_{\rm reg}\sigma^2
\left[
1-e^{-e^2/(2\sigma^2)}
\right].
\]

Near the target this is quadratic; far away its energy and restoring force are
bounded. A separate optional `tanh` force-saturating branch is available.

### 5.6 Current lifetime law

The helper now implements

\[
\tau_1=\tau_0,
\]

\[
\tau_2
=
\tau_0+\tau_{\rm directional}
+
\lambda F_{\rm opposed}
\qquad\text{in `productive_register` mode}.
\]

The break probability remains

\[
P_{\rm off}
=
\min\!\left(1,\frac{\Delta t}{\tau_s}\right).
\]

`stabilize_status1=true` restores the legacy rule
\(\tau_1=\tau_0+\lambda F\) as a diagnostic toggle, but the “load” is still the
current physical classifier's load, not the old dimensionless overlap proxy.

## 6. Mobility and integration

The legacy scalar overdamped update was schematically

\[
\Delta\mathbf r
=
\beta D\mathbf F\Delta t
+
\mathbf v_{\rm prescribed}\Delta t
+
\sqrt{2D\Delta t}\,\boldsymbol\xi.
\]

The current model removes the prescribed velocity and uses the anisotropic rod
mobility tensor

\[
\mathbf D
=
D_\parallel\mathbf u\mathbf u^\mathsf T
+
D_\perp(\mathbf I-\mathbf u\mathbf u^\mathsf T).
\]

Thus

\[
\Delta\mathbf r
=
\beta\mathbf D\mathbf F\Delta t
+
\sqrt{2\Delta t}\,\mathbf D^{1/2}\boldsymbol\xi.
\]

For aspect ratio \(p=L/d\), the finite-cylinder estimates are

\[
D_\parallel
=
\frac{k_BT[\ln p+\nu_\parallel(p)]}
     {2\pi\eta L\,c_\parallel},
\qquad
D_\perp
=
\frac{k_BT[\ln p+\nu_\perp(p)]}
     {4\pi\eta L\,c_\perp},
\]

\[
D_{\rm tumble}
=
\frac{3k_BT[\ln p+\nu_r(p)]}
     {\pi\eta L^3c_r}.
\]

The \(c\)'s are phenomenological crowding multipliers. Crosslinking does not
arbitrarily replace actin diffusion by myosin diffusion; bonds slow relative
and collective motion through their mechanical forces.

## 7. Sterics and initialization

The historical sterics-only pre-phase principally removed myosin-myosin
overlaps. Dense AA and AM contacts survived into the active phase. In current
runs:

- AA, AM, and MM steric repulsions are independent of bond state;
- full steric relaxation is normally applied;
- mature AM/AA stiffness can be ramped gradually after relaxation; and
- force, torque, displacement, and rotation caps are instrumented.

This changes the initial contact network substantially. The old high-order
state may therefore partly reflect contact-seeded kinetic assembly.

The historical weak transverse organizer
`k_bundle_max=0.1` was enabled in important high-order sweeps. The current
scaffold-free default is zero. Titin, myomesin, MyBP-C, explicit M-lines, and
explicit Z-discs remain absent unless separately enabled.

## 8. Randomness, parallelism, and checkpointing

The legacy KMC used mutable RNG engines owned by OpenMP threads. Changing
thread count or scheduling changed which random numbers an object received.

Current AM-head, AA-association, and AA-detachment draws are keyed by

\[
(\mathrm{seed},\mathrm{step},\mathrm{process},\mathrm{object\ IDs},
\mathrm{local\ draw}),
\]

so the same object receives the same stochastic stream independently of
OpenMP scheduling. Current cross-thread regression tests give bitwise-identical
HDF5 trajectories at 1, 2, 4, and 8 threads for the tested cases.

The current checkpoint includes individual material head states, stroke
progress, transported material frame, sparse ACTN2 links, recovery state,
pending events, and schema checks. It is not compatible with old checkpoints.

## 9. Parameter names whose meanings changed

| Parameter | Legacy meaning | Meaning with current defaults |
|---|---|---|
| `v_am` | Directly sets prescribed actin speed | Inactive when `dynamic_am=true`; use `am_unloaded_velocity` |
| `k_am` | Pair radial AM stiffness | Legacy AM path only; explicit proximity uses `am_tether_k` and per-head radial stiffness |
| `kappa_am` | Strong always-available pair nematic torque | Legacy AM path only; current uses head torque and `am_alignment_k` |
| `am_optimal` | Pair centerline distance | Legacy path; current preferred AM centerline follows radii plus surface gap |
| `lifetime_coeff` | Seconds per unit dimensionless overlap proxy | Lifetime increment per current physical load variable |
| `aa_status2_force_off` | Not present as a physical threshold | Used in load-hysteresis modes, but not for current `productive_register` demotion |
| `k_aa`, `kappa_aa` | Whole pair-edge stiffnesses | Divided by `aa_stiffness_reference_links` to obtain fixed per-link values |
| `max_myosin_bonds` | Maximum geometric actin partners | Compatibility/pair bookkeeping; explicit head capacity is controlled by half pools and overlap |
| scalar diffusion options | Production mobility | Used only when `physical_mobility=false` |

Merely copying the old numerical values into a current command does not create
a matched model.

## 10. Current defaults most relevant to comparison

| Quantity | Current default |
|---|---:|
| `dynamic_am` | true |
| `am_heads_per_half` | 141 |
| `am_max_heads_per_actin_half` | 24 |
| AM weak attachment / weak off / power stroke | 100 / 200 / 200 s\(^{-1}\) |
| AM strong off | 160 s\(^{-1}\) |
| Motor stiffness / stroke / stall force | 1000 pN/µm / 0.005 µm / 4.5 pN |
| Unloaded internal stroke velocity | 5 µm/s |
| Aggregate AM tether stiffness / slack | 3000 pN/µm / 0.002 µm |
| Occupied-pair alignment stiffness / slack | 30 pN µm rad\(^{-2}\) / 2° |
| `physical_mobility` | true |
| `aa_local_multilink` | true |
| local links per pair | 2 |
| AA site exclusion | 0.0192 µm |
| `aa_status_mode` | `productive_register` |
| status-2 formation threshold | 1 pN |
| registry on / off tolerances | 0.10 / 0.15 µm |
| status-1 / directional lifetime | 0.10 / 0.05 s |
| lifetime coefficient | 0.4 |
| `stabilize_status1` | false |
| AA neighbor-slot memory, status 1 / status 2 | 0.02 / 0.20 s |
| `sterics` | true |
| `k_bundle_max`, titin, myomesin | 0 |

Production manifests often override these defaults, so trajectory metadata and
the submitted command remain authoritative.

## 11. What the old parameter sweeps demonstrated

The legacy sweeps demonstrated that the filament geometry can assemble under a
coupled set of coarse organizing rules. They did not isolate a molecular
ACTN2 catch-bond law. In particular, the nonmonotonic optimum at intermediate
\(\lambda\) reflects the combination of:

- nearly continuous AM occupancy;
- strong direct AM nematic alignment;
- imposed polarity-dependent sliding;
- a strong AA axial registry;
- load stabilization of status-1 as well as status-2 bonds;
- ineffective post-rupture cooldown;
- weak transverse bundling in important runs; and
- a dense, contact-preserving initial relaxation.

At small \(\lambda\), bonds turn over and can anneal but do not retain
structure. At intermediate \(\lambda\), the imposed sorting and alignment can
be captured long enough to accumulate order. At large \(\lambda\), numerous
geometrically loaded links—including status-1 links—can trap incompatible
local arrangements. That mechanism remains a strong interpretation of the old
nonmonotonic curve, but it is a property of the entire legacy feedback loop.

## 12. What can and cannot be compared directly

Reasonably comparable observables:

- actin and myosin nematic order;
- `phi_3d`, provided its exact historical definition is used;
- AM-eligible overlap;
- AA occupancy and turnover;
- final structural snapshots.

Not directly comparable without remapping:

- the numerical value of “AM load”;
- the numerical value of \(\lambda\);
- AM bond count versus explicit head count;
- status-2 fraction;
- time measured in steps when mobility and \(\Delta t\) differ;
- force or torque magnitudes;
- axial registry energy;
- residence times when refractory kinetics differ.

The closest legacy bridge in the current executable is:

```text
dynamic_am=false
physical_mobility=false
aa_local_multilink=false
stabilize_status1=true
tau_rec=0
thin/thick feedback disabled
legacy scalar mobilities
legacy bundling and pre-phase protocol restored explicitly
```

Even this is a mechanistic bridge rather than guaranteed bitwise reproduction:
sterics, pair symmetry, bug fixes, RNG ownership, and checkpoint logic have
changed. Exact reproduction should use the archived `30829b0` executable.

## 13. Current validation status and open implications

The current explicit-head implementation passes focused geometry, ownership,
mobility, cross-bridge, ACTN2, lattice, checkpoint, and OpenMP exactness tests.
In the force-balanced four-cycle diagnostic:

- strong heads still pulse with the heartbeat;
- AM detachment statistics match integrated hazards;
- AM misalignment remains only a few degrees;
- no force, torque, or rotation caps occur; and
- no forced AM ruptures occur.

However, the current `productive_register` rule makes status 2 kinetically
persistent once formed. In the latest 0.2-s test, strong-head heartbeat
correlation was about 0.756, whereas status-2 correlation fell to about 0.166
and seven to eight links remained status 2 through most diastoles. This is
consistent with the implemented classifier, not a stochastic anomaly.

Therefore two distinct biological questions should not be conflated:

1. Does ACTN2 molecular occupancy persist while its mechanical catch state
   follows each heartbeat?
2. Does a matured, registered Z-disc link retain a long-lived state through
   transient unloading?

The legacy model did not distinguish these questions. The current model now
does so structurally, but its single status-2 variable presently represents
the second interpretation. If heartbeat-periodic catch state is required, a
separate instantaneous loaded flag should be recorded rather than forcing
molecular occupancy or registry maturity to turn off every diastole.

## Bottom line

The legacy `am_kmc` model achieved assembly through continuous geometric
attachment, imposed sliding, strong direct nematic alignment, and catch-like
retention of a coarse AA edge. The current model replaces those rules with
explicit intermittent heads, finite elastic strokes, action-reaction forces,
off-center torques, spatial ACTN2 molecules, anisotropic mobility, and
productive-register maturation.

The current model is more mechanistically interpretable, but it should not be
expected to reproduce the old phase diagram by sweeping the same numerical
\(\lambda\). The old sweep is best used as a target for a one-factor migration
ladder: first reproduce the coarse model, then replace occupancy, alignment,
sliding, AA registry, status-1 stabilization, mobility, sterics, and
initialization one at a time.
