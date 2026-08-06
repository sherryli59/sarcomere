# Transverse coherence — spec (smectic vs crystal, done right)

**Problem with what we have.** The `fig3` field plots are 1D projections (bin
along x, integrate y,z) and the `smectic_fraction` metric is **amplitude-blind**:
it asks only *which* transverse mode wins $(n_y,n_z)$, not its amplitude. A
registered 3D crystal and a dephased smectic both win at $(n_y,n_z)=(0,0)$ — the
difference is whether the coherent amplitude $|\langle e^{i(q^*,0,0)\cdot r}\rangle|$
is **large** (phases locked across the cross-section) or **collapses** (columns
dephase). That distinction — registered Z-sheets vs dephased columns — is the
heart of the crosslinker hypothesis and is invisible in the projection.

**Goal.** Turn "smectic frac = 0.91" into a sharp, amplitude-resolved statement
about *how registered* the lattice is, using the y,z information in the 3D
trajectory we already have.

---

## 1. Per-column phase order parameter

Use the longitudinal winning wavevector $q^*=2\pi n_x^*/L_x$ (already measured;
take it from the cb2 population). For each filament $i$ at $(x_i,y_i,z_i)$ it
carries a longitudinal phase $e^{-iq^* x_i}$. Bin filaments into a transverse
grid of cells $c$ over $[0,L_y)\times[0,L_z)$ (PBC). Per cell:

$$\psi_c \;=\; \frac{1}{n_c}\sum_{i\in c} e^{-i q^* x_i}$$

- $|\psi_c|$ = how strongly *that column* is banded along x (longitudinal order).
- $\arg\psi_c$ = the **phase/offset** of that column's bands. Registration is
  whether $\arg\psi_c$ is **locked across cells**.

## 2. Registration order parameter (the amplitude ratio)

$$A_\mathrm{coh}=\Big|\tfrac{1}{N}\sum_i e^{-iq^* x_i}\Big|
\quad\text{(transverse-coherent amplitude} = S(q^*,0,0)),$$
$$A_\mathrm{col}=\frac{1}{N}\sum_c n_c\,|\psi_c|\quad\text{(mean per-column band amplitude)},$$
$$\boxed{\,R = A_\mathrm{coh}/A_\mathrm{col}\in[0,1]\,}$$

- $R\to1$: every column shares the same phase → **registered crystal** (Z-sheets).
- $R\to0$: columns strongly banded individually but phases random → **dephased
  smectic**; the bands cancel transversely.

$R$ is a *global* ratio, so it is computable even in the thin 2×2 box (it only
asks whether the few columns are in phase, not a decay length). This is the
headline number that replaces `smectic_fraction`.

## 3. Transverse correlation function (needs width — caveat)

$$C_\perp(\Delta) = \frac{\langle \psi_c\,\psi_{c'}^* \rangle_{|r_c-r_{c'}|=\Delta}}
{\langle |\psi_c|^2\rangle}$$

Decay length $\xi_\perp$ = transverse coherence length. **In 2×2 this is barely
resolvable** (≈ a few lateral spacings across the box); report it with the
caveat and treat $R$ as primary. $\xi_\perp \gtrsim L_\perp$ ⇒ registered;
$\xi_\perp \ll L_\perp$ ⇒ dephased.

## 4. Transverse coherence of $\sigma$ (the cleanest field)

$\sigma$ (crosslink tension) is near-binary, so it's the best field for this.
Bin aa-bond loads at their PBC midpoints into 3D cells; at each x-band, ask
whether a tension plateau **spans the cross-section** (registered crosslinking)
or is **patchy in y,z**. Compute the same $R_\sigma$ and $C_\perp^\sigma$ on the
$\sigma$ field. If $\sigma$ registers but $\rho$ doesn't (or vice versa), that
localizes where registration is set.

---

## 5. What's measurable now vs needs new sims

- **Now (data0430, 12×2×2, k_bundle=0.1):** $R$, $R_\sigma$ across the full
  A×lt grid — *how registered* the lattice is at fixed coupling, and how that
  depends on density and catch-bond strength.
- **`k_bundle` is the WRONG knob for the mechanism test.** It pulls all particles
  to the transverse COM (confinement); it does not couple band phase across y,z.
  The lateral-registration mechanisms in the model are **myosin multi-binding**
  (`max_myosin_bonds`) and the **aa crosslink network** bridging transversely-
  offset filaments. The registration-mechanism test = sweep `max_myosin_bonds`
  (and/or crosslink connectivity), *not* `k_bundle`.
- **Wider cross-section sim (recommended one new run):** same axial length, much
  wider transverse box, e.g. **L 12×8×8**, to actually resolve $\xi_\perp$ and
  $C_\perp$. 2×2 can measure $R$ but is too thin to map the dephasing length.

---

## 6. Outputs / build

- `analysis/transverse_coherence.py` — per run: $q^*$, $R$, $R_\sigma$, $|\psi_c|$
  map, $C_\perp(\Delta)$, time series of $R(t)$. Writes `tc_<run>.npz`.
- Aggregate $R$ over the A×lt grid (replicas averaged) → heatmap, alongside the
  existing period/smectic phase diagrams.
- Plot: $\psi_c$ in the complex plane per frame (are the column phasors aligned?)
  — the most direct visual of registration.

Transverse grid: with ~1200 actins in 12×2×2, a 4×4 transverse grid gives ~75
filaments/cell — enough for $\psi_c$; finer grids run out of statistics in 2×2.
