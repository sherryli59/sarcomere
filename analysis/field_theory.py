"""
Preliminary continuum field theory for sarcomere self-assembly (1D, along x).

This is a MINIMAL, RUNNABLE model to test ideas -- NOT the final grounded theory.
Every term is labeled by provenance:

  [D] DERIVED from symmetry (holds without any fitting):
      - flux-form conservation for rho, p:  d_t(.) = -d_x J(.)
      - parity-forced period-doubling: p^2 is the lowest even-in-p coupling to rho,
        so density must carry the 2q harmonic of the polarity modulation.

  [W] DISCOVERED by WSINDy -- but the v2 PARITY filter splits the sigma coupling
      by target, and the split matters (see the PARITY WARNING on the engine below):
      - POLARITY d_t p (odd, conserved): carries d1[sigma], d1[sigma^2] (alpha=1 ->
        odd -> LEGAL for the odd target p). Survived v2. BUT with sigma=f(p^2) (even)
        these live at 2q, so they cannot sustain the fundamental-q polarity band.
      - DENSITY d_t rho (even, conserved): carries d2[sigma], d2[rho*sigma] (alpha=2
        -> even -> LEGAL for rho). This is a SAME-SIGN chemotactic drift (both
        polarities climb the same sigma-gradient), NOT an opposite-sign sorter.
      - d1[sigma] in d_t rho is PARITY-FORBIDDEN (odd term, even target) and was
        DELETED by the v2 filter. Do NOT attribute density transport to it.
      - sigma fast, near-binary; NOT slaved to rho alone (corr(sigma,rho)~0).

  [P] PLACEHOLDER (added for well-posedness / because the real thing is sub-grid):
      - steric entropic log-barrier  f ~ rho+ ln rho+ + rho- ln rho-  (bounds |p|<=rho).
      - Swift-Hohenberg kappa*(d_x^2 + q0^2)^2 for wavelength selection. The review
        showed the length is STERIC (nonlocal l0 kernel); SH is a surrogate stand-in
        (l0 ~ dx is sub-grid, so the true kernel can't yet be resolved). FLAGGED.

State: rho_+(x), rho_-(x) (the two polarity populations, both conserved) and
sigma(x) (tension, fast reaction). Diagnostics use rho = rho_+ + rho_-, p = rho_+ - rho_-.

The active engine is an OPPOSITE-SIGN sorting current:
    J_+ ~ +v * rho_+ * d_x sigma ,   J_- ~ -v * rho_- * d_x sigma
opposite polarities drift oppositely in the tension gradient. The coefficient v is
set by catch-bond strength, and v=0 IS the catch-bond knockout (mirrors
cpp/catch_knockout.sbatch) -- so the model is testable against the intervention.

    *** PARITY WARNING (do not treat this engine as WSINDy-derived) ***
    Under P: x->-x, +<->- (with rho even, p ODD, sigma even), this opposite-sign
    current makes d_t rho = -v d_x(p d_x sigma) map to MINUS itself while rho is
    even -> PARITY-VIOLATING. It is the continuum face of d1[sigma] in the even rho
    equation -- exactly the odd-parity column the v2 filter DELETES. The RELATIVE
    MINUS that makes it a sorter is the same thing that breaks parity. The
    parity-LEGAL sigma coupling (v2) is the SAME-SIGN chemotactic drift
    d2[sigma]/d2[rho*sigma], which is not a sorter. A parity-clean sharpening engine
    must be ODD-IN-P, e.g. J_+ ~ +w rho_+ d_x p (gives d_t p = -w d_x(rho d_x p),
    which IS P-covariant and acts at the fundamental q). TODO: replace this engine.

Numerics: pseudo-spectral in x (periodic, FFT derivatives); ETD1 exponential-Euler
for the stiff linear part (diffusion + SH) treated exactly in Fourier, everything
else explicit. sigma relaxes with fast linear kinetics to a bistable setpoint.
"""
import numpy as np
from dataclasses import dataclass, field


# %% parameters --------------------------------------------------------------
@dataclass
class Params:
    # --- geometry / discretization ---
    Lx: float = 12.0            # match the sims (box length)
    N: int = 256               # grid points
    dt: float = 2e-4           # time step
    nsteps: int = 200000       # integration steps
    save_every: int = 1000

    # --- mean density (maps to A = n_actins) ---
    rho_bar: float = 1.0        # mean total density; sets packing fraction

    # --- [W] active non-reciprocal sorting current (maps to catch-bond strength lt) ---
    v: float = 1.0              # sorting velocity;  v = 0  -> CATCH-BOND KNOCKOUT
    sigma_grad_power: int = 1   # 1 -> d_x sigma (leading WSINDy term); 2 -> d_x(sigma^2)
    sorting_sign: float = +1.0  # sign of the sorting current J_+ = sign*v*rho_+*d_x sigma.
                                # +1 (with mode=core) DISSOLVES; -1 AMPLIFIES (filaments
                                # pulled INTO high-tension cores -> sharpening/assembly).

    # --- [P] transport regularization ---
    D: float = 0.02            # bare diffusion of each population
    Gamma: float = 0.02        # mobility for the free-energy (steric) current
    eps_floor: float = 1e-4    # floor inside ln() for the steric log-barrier

    # --- [P] Swift-Hohenberg wavelength selection (SURROGATE for steric l0 kernel) ---
    use_SH: bool = True
    kappa: float = 0.02        # SH stiffness
    n_bands: float = 6.0       # target band count -> q0 = 2*pi*n_bands/Lx

    # --- [W] sigma reaction (fast bistable auxiliary) ---
    tau_sigma: float = 0.02    # relaxation time (fast -> sigma slaves quickly)
    sigma_max: float = 1.0
    sigma_thresh: float = 0.15 # bistable switch threshold on the (p/rho)^2 signal
    sigma_sharp: float = 6.0   # switch sharpness
    sigma_mode: str = "core"   # MODELING CHOICE (flag): where sigma concentrates.
                               #   "core": sigma high where sorting order (p/rho)^2 is
                               #           HIGH -> peaks in polar cores -> d_x sigma points
                               #           inward -> sorting current AMPLIFIES bands
                               #           (assembling regime; positive feedback).
                               #   "seam": sigma high where (p/rho)^2 is LOW -> peaks in
                               #           antiparallel seams (crosslink-tension picture)
                               #           -> tends to DISSOLVE. Flip to test the biology.

    seed: int = 0

    # --- initial condition ---
    seed_amp: float = 0.0      # 0 -> near-uniform+noise (tests spontaneous growth);
                               # >0 -> finite-amplitude banded polarity seed at n_bands
                               # (tests NUCLEATION: does a seeded sarcomere pattern
                               # survive/sharpen? homogeneous state is linearly stable,
                               # so assembly here is nucleation-type, as in the sim).

    def q0(self):
        return 2.0 * np.pi * self.n_bands / self.Lx


# %% spectral helpers --------------------------------------------------------
def _k(P):
    return 2.0 * np.pi * np.fft.fftfreq(P.N, d=P.Lx / P.N)   # wavenumbers

def ddx(fhat_or_f, k, spectral=False):
    """d/dx via FFT. Pass real-space f (spectral=False) or its FFT (spectral=True)."""
    fhat = fhat_or_f if spectral else np.fft.fft(fhat_or_f)
    return np.real(np.fft.ifft(1j * k * fhat))


# %% sigma setpoint (bistable in antiparallel overlap) -----------------------
def sigma_setpoint(rho_p, rho_m, P):
    """[W] sigma is a fast bistable auxiliary that is BAND-SELECTIVE on the local
    sorting order  m2 = (p/rho)^2  in [0,1]  (0 = antiparallel seam, 1 = polar core).
    m2 is parity-EVEN (as sigma must be) and needs both populations (-> not slaved to
    rho alone, matching corr(sigma,rho)~0). A sharp switch makes sigma near-binary.

    sigma_mode selects WHERE sigma concentrates (a modeling choice, see Params):
      "core": signal = m2      -> sigma peaks in cores  -> d_x sigma inward -> AMPLIFY
      "seam": signal = 1 - m2  -> sigma peaks in seams  -> tends to DISSOLVE
    The spatial structure is the point: uniform sigma has no gradient, so the sorting
    current (v * rho_pm * d_x sigma) is idle regardless of v -- that was the null
    catch-bond test in the first band-seeded run."""
    rho = rho_p + rho_m
    p = rho_p - rho_m
    m2 = (p / (rho + P.eps_floor)) ** 2                  # sorting order in [0,1], even in p
    signal = m2 if P.sigma_mode == "core" else (1.0 - m2)
    switch = 0.5 * (1.0 + np.tanh(P.sigma_sharp * (signal - P.sigma_thresh)))
    return P.sigma_max * switch


# %% right-hand side (explicit, nonlinear part) ------------------------------
def nonlinear_rhs(rho_p, rho_m, sigma, P, k):
    """Return the EXPLICIT part of d_t rho_+/-, i.e. -d_x J_explicit, plus the sigma RHS.
    The stiff linear part (diffusion + SH) is handled by the ETD integrator, so it
    is EXCLUDED here."""
    dsig = ddx(sigma, k)
    if P.sigma_grad_power == 2:
        drive = ddx(sigma ** 2, k)                      # d_x(sigma^2)
    else:
        drive = dsig                                    # d_x(sigma)

    # OPPOSITE-SIGN sorting current: opposite polarities drift oppositely.
    # PARITY WARNING (see module docstring): with sigma even, the relative minus
    # makes d_t rho parity-VIOLATING (= the forbidden d1[sigma]-in-rho column). NOT
    # WSINDy-derived. The parity-clean sharpening engine is odd-in-p (rho_+ * d_x p).
    Jsort_p = +P.sorting_sign * P.v * rho_p * drive
    Jsort_m = -P.sorting_sign * P.v * rho_m * drive

    # [P] steric log-barrier current: -Gamma * rho * d_x( mu_steric ), mu = ln(rho+eps)
    mu_p = np.log(np.maximum(rho_p, 0.0) + P.eps_floor)
    mu_m = np.log(np.maximum(rho_m, 0.0) + P.eps_floor)
    Jster_p = P.Gamma * rho_p * ddx(mu_p, k)
    Jster_m = P.Gamma * rho_m * ddx(mu_m, k)

    drho_p = -ddx(Jsort_p + Jster_p, k)
    drho_m = -ddx(Jsort_m + Jster_m, k)

    # [W] sigma fast relaxation to bistable setpoint
    S = sigma_setpoint(rho_p, rho_m, P)
    dsigma = (S - sigma) / P.tau_sigma
    return drho_p, drho_m, dsigma


# %% linear operator (diffusion + SH) in Fourier -----------------------------
def linear_op(P, k):
    """L(k) acting on each population's FFT (applied to rho_+ and rho_- identically).
    CONSERVED form: d_t rho = d_xx mu, so every rate carries an overall -k^2 and
    L(k=0)=0 (mass is preserved). Diffusion: -D k^2. [P] SH free energy
    mu_SH = kappa(d_x^2+q0^2)^2 rho -> conserved rate -kappa k^2 (q0^2-k^2)^2, which
    (unlike a bare SH rate) is ZERO at k=0. Selects the q0 band; stabilizing -> ETD stable.
    NOTE run#1 bug: the SH term lacked the k^2 factor -> drained the mean density."""
    L = -P.D * k ** 2
    if P.use_SH:
        L = L - P.kappa * k ** 2 * (P.q0() ** 2 - k ** 2) ** 2
    return L


# %% ETD1 exponential-Euler step ---------------------------------------------
def etd1_coeffs(L, dt):
    """e^{L dt} and (e^{L dt} - 1)/L, with the L->0 limit handled."""
    E = np.exp(L * dt)
    with np.errstate(divide="ignore", invalid="ignore"):
        phi = np.where(np.abs(L) > 1e-12, (E - 1.0) / L, dt)
    return E, phi


# %% integrator --------------------------------------------------------------
def run(P: Params, verbose=True):
    rng = np.random.default_rng(P.seed)
    x = np.arange(P.N) * (P.Lx / P.N)
    k = _k(P)
    L = linear_op(P, k)
    E, phi = etd1_coeffs(L, P.dt)

    # initial condition
    base = P.rho_bar * (1.0 + 0.01 * rng.standard_normal(P.N))
    if P.seed_amp > 0.0:
        # finite-amplitude banded polarity seed at the target wavenumber (nucleation test)
        pol0 = P.seed_amp * P.rho_bar * np.sin(P.q0() * x)
    else:
        pol0 = 0.01 * rng.standard_normal(P.N)          # spontaneous-growth test
    pol0 = np.clip(pol0, -base + P.eps_floor, base - P.eps_floor)   # keep |p|<=rho
    rho_p = 0.5 * (base + pol0)
    rho_m = 0.5 * (base - pol0)
    sigma = sigma_setpoint(rho_p, rho_m, P)

    snaps = {"t": [], "rho": [], "p": [], "sigma": [], "rho_p": [], "rho_m": []}
    for step in range(P.nsteps + 1):
        if step % P.save_every == 0:
            rho = rho_p + rho_m
            snaps["t"].append(step * P.dt)
            snaps["rho"].append(rho.copy())
            snaps["p"].append((rho_p - rho_m).copy())
            snaps["sigma"].append(sigma.copy())
            snaps["rho_p"].append(rho_p.copy()); snaps["rho_m"].append(rho_m.copy())
            if verbose and step % (P.save_every * 20) == 0:
                nx, r2q = diagnostics_frame(rho_p - rho_m, rho_p + rho_m, P)
                print(f"  step {step:7d}  t={step*P.dt:7.2f}  "
                      f"|p|max={np.abs(rho_p-rho_m).max():.3f}  nx={nx}  "
                      f"2q/q(rho)={r2q:.2f}", flush=True)

        Np, Nm, Ds = nonlinear_rhs(rho_p, rho_m, sigma, P, k)
        # ETD1 for the two conserved fields (linear part exact in Fourier)
        rho_p = np.real(np.fft.ifft(E * np.fft.fft(rho_p) + phi * np.fft.fft(Np)))
        rho_m = np.real(np.fft.ifft(E * np.fft.fft(rho_m) + phi * np.fft.fft(Nm)))
        # positivity safety (log-barrier should prevent this; clip as backstop)
        rho_p = np.maximum(rho_p, 0.0)
        rho_m = np.maximum(rho_m, 0.0)
        sigma = sigma + P.dt * Ds                       # explicit fast relaxation
        sigma = np.clip(sigma, 0.0, P.sigma_max)

        if not np.isfinite(rho_p).all():
            print(f"  BLEW UP at step {step}", flush=True); break

    for kk in snaps:
        snaps[kk] = np.asarray(snaps[kk])
    snaps["x"] = x
    return snaps


# %% diagnostics -------------------------------------------------------------
def diagnostics_frame(p, rho, P):
    """Return (band count nx from |p| spectrum, power ratio 2q/q in rho).
    Period-doubling [D] predicts rho carries the 2q harmonic of p -> ratio ~O(1)+."""
    pk = np.abs(np.fft.rfft(p - p.mean()))
    if pk[1:].max() <= 0:
        return 0, 0.0
    q_idx = 1 + int(np.argmax(pk[1:]))                  # dominant polarity mode
    rk = np.abs(np.fft.rfft(rho - rho.mean()))
    p_at_q = rk[q_idx] if q_idx < len(rk) else 0.0
    p_at_2q = rk[2 * q_idx] if 2 * q_idx < len(rk) else 0.0
    ratio = p_at_2q / (p_at_q + 1e-12)
    return q_idx, ratio


def summarize(snaps, P):
    """Final-state summary: band count, period-doubling ratio, well-posedness check."""
    p, rho, sig = snaps["p"][-1], snaps["rho"][-1], snaps["sigma"][-1]
    nx, r2q = diagnostics_frame(p, rho, P)
    rho_pf, rho_mf = snaps["rho_p"][-1], snaps["rho_m"][-1]
    viol = np.max(np.abs(p) - rho)                      # |p|<=rho must hold ([P] barrier)
    return {
        "nx_bands": int(nx),
        "period_doubling_ratio": float(r2q),
        "p_amplitude": float(np.abs(p).max()),
        "sigma_mean": float(sig.mean()),
        "sigma_binary_frac": float(np.mean((sig < 0.1) | (sig > 0.9))),
        "wellposed_|p|<=rho": bool(viol <= 1e-6),
        "max_violation": float(viol),
    }


# %% main: baseline vs catch-bond knockout -----------------------------------
if __name__ == "__main__":
    import argparse, json, os
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out_dir", default="/scratch/users/xsli/sarcomere/field_theory")
    ap.add_argument("--nsteps", type=int, default=200000)
    ap.add_argument("--N", type=int, default=256)
    ap.add_argument("--n_bands", type=float, default=6.0)
    ap.add_argument("--v", type=float, default=1.0, help="sorting velocity (catch-bond strength)")
    ap.add_argument("--seed_amp", type=float, default=0.3,
                    help="banded polarity seed amplitude (nucleation test); 0 = noise IC")
    ap.add_argument("--no_SH", action="store_true", help="disable Swift-Hohenberg length selection")
    ap.add_argument("--sigma_mode", default="core", choices=("core","seam"),
                    help="where sigma concentrates: core (amplify?) vs seam (dissolve?) -- test both")
    ap.add_argument("--sigma_thresh", type=float, default=0.15,
                    help="switch threshold on (p/rho)^2; MUST sit inside the seed m2 range to band sigma")
    ap.add_argument("--sorting_sign", type=float, default=1.0,
                    help="+1 dissolves (core), -1 amplifies/assembles")
    ap.add_argument("--knockout_only", action="store_true")
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    common = dict(nsteps=args.nsteps, N=args.N, n_bands=args.n_bands,
                  use_SH=not args.no_SH, seed_amp=args.seed_amp, sigma_mode=args.sigma_mode,
                  sigma_thresh=args.sigma_thresh, sorting_sign=args.sorting_sign)

    results = {}
    runs = [] if args.knockout_only else [("baseline", args.v)]
    runs.append(("knockout", 0.0))
    for tag, v in runs:
        print(f"\n=== {tag}  (v={v}) ===", flush=True)
        P = Params(v=v, **common)
        snaps = run(P)
        summ = summarize(snaps, P)
        results[tag] = summ
        print(f"  -> {json.dumps(summ)}", flush=True)
        np.savez_compressed(os.path.join(args.out_dir, f"ft_{tag}.npz"),
                            x=snaps["x"], t=snaps["t"], rho=snaps["rho"],
                            p=snaps["p"], sigma=snaps["sigma"])

    with open(os.path.join(args.out_dir, "ft_summary.json"), "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nwrote {args.out_dir}/ft_summary.json")
    if "baseline" in results and "knockout" in results:
        b, ko = results["baseline"], results["knockout"]
        print(f"\n=== CATCH-BOND TEST (field theory) ===")
        print(f"  baseline: {b['nx_bands']} bands, p_amp={b['p_amplitude']:.3f}")
        print(f"  knockout: {ko['nx_bands']} bands, p_amp={ko['p_amplitude']:.3f}")
        print(f"  -> sorting current {'NECESSARY (bands vanish at v=0)' if ko['p_amplitude'] < 0.3*b['p_amplitude'] else 'NOT decisive'}")
