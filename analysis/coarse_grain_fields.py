"""
Coarse-grain KMC sarcomere trajectories into the continuum fields a field
theory would evolve, and measure -- over the *whole* trajectory, not just the
steady-state tail -- the three things that decide what the field theory has to
reproduce:

  1. q_dom(t)        : selection vs coarsening.  Does the dominant longitudinal
                       wavevector lock in early (pattern *selection*) or drift
                       toward smaller q as t grows (Ostwald *coarsening*)?
  2. period vs density / lt : is the band spacing set by the microscopic overlap
                       length (am_optimal/aa_optimal, fixed at 0.06 here) -- a
                       steric/DFT length, density-independent -- or is it an
                       emergent wavelength that moves with actin number / catch
                       bond strength?
  3. transverse coherence : is the Z-band a smectic (modulated along x only;
                       winning transverse mode (ny,nz)=(0,0), liquid-like in
                       y,z) or a crystal (Bragg peaks at finite k_perp)?

It also bins the continuum fields themselves, which are the inputs a later
weak-form sparse-regression (SINDy) pass would regress against:

  rho_pm(x,t)  +/- actin density, split by instantaneous polarity sign(dir_x)
  rho(x,t)     = rho_+ + rho_-                      (even / total density)
  p(x,t)       = rho_+ - rho_-                      (odd / polarity field)
  m(x,t)       bound actin-myosin density           (the motor field)
  sigma(x,t)   crosslink tension = mean bond_pair_load per bin  (measured, not
               assumed -- this is the "gold" field)
  n_bond(x,t)  crosslink (alpha-actinin) bond density
  rho_cb2(x,t) density of doubly-engaged (cb_status==2) actins -- the Z-disc
               -building load-bearing crosslink

Schema facts this relies on (verified against the HDF5 and the C++):
  axis = x ; box parsed from filename token Lxx<Ly>x<Lz> (e.g. L12x2x2 -> 12,2,2)
  /actin/center      {T,Na,3}    /actin/direction {T,Na,3}  (polarity = sign dir[:,0])
  /actin/cb_status   {T,Na,1}    0 unbound, 1 single, 2 doubly engaged
  /actin/bonds       {T,Nb,2}    actin-actin crosslink pairs; empty slot = (-1,-1)
  /actin/bond_pair_load {T,Nb,1} per-bond load ~ catch-bond f_load
  /actin_myo/bonds   {T,Nb,2}    (actin_idx, myosin_idx); empty slot = (-1,-1)

Catch-bond law (src/sarcomere_2d.cpp:997):
  k_off ~ dt * |cos angle| / (base_lifetime + lt * f_load)
  -> lt (filename token 'lt' = --lifetime_coeff) is the catch-bond strength.

Reuses the single-mode Fourier order parameter from extract_structure_3d.py so
q_dom is measured exactly the way the steady-state phase diagram already is.
"""

import os
import re
import sys
import argparse
import numpy as np
import h5py

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from extract_structure_3d import (  # noqa: E402
    build_k_matrix,
    build_mode_matrix,
    fourier_order_3d_frames,
    resolve_compute_device,
)

# --------------------------------------------------------------------------- #
# box / file helpers
# --------------------------------------------------------------------------- #
_BOX_RE = re.compile(r"_L(\d+(?:\.\d+)?)x(\d+(?:\.\d+)?)x(\d+(?:\.\d+)?)")
_LT_RE = re.compile(r"_lt([^_]+)")
_A_RE = re.compile(r"A(\d+)_M")


def parse_box(name):
    m = _BOX_RE.search(name)
    if not m:
        raise ValueError(f"cannot parse box LxXxYxZ from {name!r}")
    return tuple(float(v) for v in m.groups())


def parse_lt(name):
    m = _LT_RE.search(name)
    return m.group(1) if m else None


def parse_actin_count(name):
    m = _A_RE.search(name)
    return int(m.group(1)) if m else None


def n_frames(fh):
    """Number of committed frames (meta/committed_frames if present, else shape)."""
    if "/meta/committed_frames" in fh:
        c = int(np.asarray(fh["/meta/committed_frames"])[0])
        if c > 0:
            return c
    return fh["/actin/center"].shape[0]


# --------------------------------------------------------------------------- #
# periodic binning utilities
# --------------------------------------------------------------------------- #
def periodic_histogram(x, L, nbins, weights=None):
    """Histogram of x in [0,L) with PBC wrap; returns density per unit length."""
    xw = np.mod(x, L)
    counts, edges = np.histogram(xw, bins=nbins, range=(0.0, L), weights=weights)
    dx = L / nbins
    return counts / dx, 0.5 * (edges[:-1] + edges[1:])


def periodic_gaussian_smooth(field, L, sigma):
    """Smooth a periodic 1D field with a Gaussian of width sigma (same units as L)."""
    if sigma <= 0:
        return field
    n = field.shape[-1]
    dx = L / n
    k = 2.0 * np.pi * np.fft.rfftfreq(n, d=dx)
    ker = np.exp(-0.5 * (k * sigma) ** 2)
    return np.fft.irfft(np.fft.rfft(field) * ker, n=n)


def periodic_midpoint(xa, xb, L):
    """PBC-aware midpoint of two x-coordinates via circular mean."""
    ta = 2 * np.pi * np.mod(xa, L) / L
    tb = 2 * np.pi * np.mod(xb, L) / L
    s = np.sin(ta) + np.sin(tb)
    c = np.cos(ta) + np.cos(tb)
    tm = np.mod(np.arctan2(s, c), 2 * np.pi)
    return tm * L / (2 * np.pi)


# --------------------------------------------------------------------------- #
# field coarse-graining for a single frame
# --------------------------------------------------------------------------- #
def coarse_grain_frame(fh, idx, Lx, nbins=240, smooth=0.0):
    """Return dict of 1D fields along x for one frame."""
    centers = np.asarray(fh["/actin/center"][idx], dtype=float)        # (Na,3)
    dirs = np.asarray(fh["/actin/direction"][idx], dtype=float)        # (Na,3)
    cb_status = np.asarray(fh["/actin/cb_status"][idx]).reshape(-1)    # (Na,)

    xa = centers[:, 0]
    plus = dirs[:, 0] > 0
    minus = dirs[:, 0] < 0

    rho_p, xc = periodic_histogram(xa[plus], Lx, nbins)
    rho_m, _ = periodic_histogram(xa[minus], Lx, nbins)
    rho_cb2, _ = periodic_histogram(xa[cb_status == 2], Lx, nbins)

    # motor fields from actin-myosin bonds (NOT the same as the aa catch bond):
    #   m_myo : distinct bound myosins, binned at the MYOSIN center -> bound-myosin
    #           density (deduped; a myosin bound to k actins counts once).
    #   m_cb2 : am attachments whose actin has cb_status==2, binned at the actin
    #           center -> motor engagement on the Z-disc-building filaments
    #           (population-consistent with the q_dom / R structure analysis).
    am = np.asarray(fh["/actin_myo/bonds"][idx]).reshape(-1, 2)
    am_valid = (am[:, 0] >= 0) & (am[:, 1] >= 0)
    a_bound = am[am_valid, 0].astype(int)
    m_bound = am[am_valid, 1].astype(int)
    myo_centers = np.asarray(fh["/myosin/center"][idx], dtype=float)
    m_myo, _ = periodic_histogram(myo_centers[np.unique(m_bound), 0], Lx, nbins)
    cb2_attach = a_bound[cb_status[a_bound] == 2]
    m_cb2, _ = periodic_histogram(xa[cb2_attach], Lx, nbins)

    # crosslink tension: aa bonds with their load, binned at the PBC bond midpoint
    bonds = np.asarray(fh["/actin/bonds"][idx]).reshape(-1, 2)
    load = np.asarray(fh["/actin/bond_pair_load"][idx]).reshape(-1)
    b_valid = (bonds[:, 0] >= 0) & (bonds[:, 1] >= 0)
    bi = bonds[b_valid, 0].astype(int)
    bj = bonds[b_valid, 1].astype(int)
    bload = load[b_valid]
    xmid = periodic_midpoint(xa[bi], xa[bj], Lx)
    n_bond, _ = periodic_histogram(xmid, Lx, nbins)
    load_sum, _ = periodic_histogram(xmid, Lx, nbins, weights=bload)
    # mean load per unit length / bond density -> mean tension where bonds exist
    with np.errstate(invalid="ignore", divide="ignore"):
        sigma = np.where(n_bond > 0, load_sum / n_bond, 0.0)

    out = {
        "x": xc,
        "rho_plus": rho_p,
        "rho_minus": rho_m,
        "rho": rho_p + rho_m,
        "p": rho_p - rho_m,
        "m_myo": m_myo,
        "m_cb2": m_cb2,
        "sigma": sigma,
        "n_bond": n_bond,
        "rho_cb2": rho_cb2,
        "n_aa_bonds": int(b_valid.sum()),
        "n_am_bonds": int(am_valid.sum()),
        "n_cb2": int((cb_status == 2).sum()),
    }
    if smooth > 0:
        for key in ("rho_plus", "rho_minus", "rho", "p", "m_myo", "m_cb2", "sigma", "n_bond", "rho_cb2"):
            out[key] = periodic_gaussian_smooth(out[key], Lx, smooth)
    return out


# --------------------------------------------------------------------------- #
# q_dom(t): selection vs coarsening
# --------------------------------------------------------------------------- #
def qdom_timeseries(fh, Lx, Ly, Lz, frame_idx, population="cb2",
                    n_max_x=30, n_max_perp=5):
    """
    Winning longitudinal mode nx*, amplitude Phi, transverse modes ny*,nz* and
    tilt theta, per sampled frame. Single pass; the mode grid is built once and
    the needed arrays are bulk-read once (datasets are uncompressed, chunked
    {10,Na,3}, so one sequential read beats thousands of strided reads), then we
    score each frame against the prebuilt k-matrix. population:
      'cb2'  -> doubly-engaged actins (Z-disc builders; matches phase diagram)
      'all'  -> all actins
      'plus'/'minus' -> one polarity population
    """
    fr = np.asarray(frame_idx, dtype=int)
    nx_vals = np.arange(1, n_max_x + 1)
    ny_vals = np.arange(-n_max_perp, n_max_perp + 1)
    nz_vals = np.arange(-n_max_perp, n_max_perp + 1)
    k_matrix = build_k_matrix(Lx, Ly, Lz, nx_values=nx_vals,
                              ny_values=ny_vals, nz_values=nz_vals)        # (K,3)
    mode_matrix = build_mode_matrix(nx_values=nx_vals, ny_values=ny_vals,
                                    nz_values=nz_vals)                     # (K,3)
    kT = np.ascontiguousarray(k_matrix.T)                                 # (3,K)

    # bulk read: full dataset once, then index the strided frames in memory
    centers = np.asarray(fh["/actin/center"][:], dtype=float)[fr]         # (F,Na,3)
    if population in ("plus", "minus"):
        dx = np.asarray(fh["/actin/direction"][:], dtype=float)[fr, :, 0]
        masks = dx > 0 if population == "plus" else dx < 0
    elif population == "cb2":
        masks = np.asarray(fh["/actin/cb_status"][:])[fr].reshape(len(fr), -1) == 2
    else:
        masks = np.ones(centers.shape[:2], dtype=bool)

    F = len(fr)
    nx_star = np.zeros(F, dtype=int)
    ny_star = np.zeros(F, dtype=int)
    nz_star = np.zeros(F, dtype=int)
    phi = np.zeros(F)
    theta = np.zeros(F)
    npart = np.zeros(F, dtype=int)

    for j in range(F):
        coords = centers[j][masks[j]]
        npart[j] = coords.shape[0]
        if coords.shape[0] < 2:
            continue
        phases = coords @ kT                       # (n,K)
        real = np.cos(phases).mean(axis=0)         # (K,)
        imag = np.sin(phases).mean(axis=0)
        amp = np.sqrt(real * real + imag * imag)
        idx = int(np.argmax(amp))
        nx_star[j], ny_star[j], nz_star[j] = (int(mode_matrix[idx, 0]),
                                              int(mode_matrix[idx, 1]),
                                              int(mode_matrix[idx, 2]))
        phi[j] = float(amp[idx])
        kperp = np.hypot(k_matrix[idx, 1], k_matrix[idx, 2])
        theta[j] = float(np.degrees(np.arctan2(kperp, abs(k_matrix[idx, 0]))))

    with np.errstate(divide="ignore", invalid="ignore"):
        period = np.where(nx_star > 0, Lx / nx_star, np.nan)
    qdom = 2.0 * np.pi * nx_star / Lx
    return {
        "frame": fr,
        "nx_star": nx_star,
        "ny_star": ny_star,
        "nz_star": nz_star,
        "q_dom": qdom,
        "period": period,
        "phi": phi,
        "theta_deg": theta,
        "n_particles": npart,
    }


def smectic_fraction_from_ts(ts):
    """Fraction of ordered frames whose winning transverse mode is (0,0):
    ->1 means pure longitudinal layering (smectic); <1 means transverse order."""
    valid = ts["nx_star"] > 0
    is_smectic = valid & (ts["ny_star"] == 0) & (ts["nz_star"] == 0)
    return float(is_smectic.sum() / max(1, int(valid.sum()))), int(valid.sum())


# --------------------------------------------------------------------------- #
# driver for one file
# --------------------------------------------------------------------------- #
def process_file(h5path, out_npz, stride=5, nbins=240, smooth=0.0,
                 n_snapshots=6, population="cb2", device="numpy",
                 start_frame=0, end_frame=-1):
    name = os.path.basename(h5path)
    Lx, Ly, Lz = parse_box(name)

    with h5py.File(h5path, "r") as fh:
        nf = n_frames(fh)
        frame_idx = list(range(0, nf, stride))

        ts = qdom_timeseries(fh, Lx, Ly, Lz, frame_idx, population=population)
        smectic_frac, _ = smectic_fraction_from_ts(ts)

        # binned fields at n_snapshots evenly spaced snapshots in [start,end].
        # For PDE discovery, concentrate these in a window with real dynamics
        # (assembly transient), not pure steady state where d_t u ~ 0.
        lo = max(0, start_frame)
        hi = (nf - 1) if end_frame < 0 else min(end_frame, nf - 1)
        snap_frames = np.unique(np.linspace(lo, hi, n_snapshots).astype(int))
        snaps = {}
        for sf in snap_frames:
            cg = coarse_grain_frame(fh, int(sf), Lx, nbins=nbins, smooth=smooth)
            for k, v in cg.items():
                snaps.setdefault(k, []).append(v if np.ndim(v) else v)
        x_axis = snaps.pop("x")[0]

    payload = {
        "fname": name,  # not "file": collides with np.savez's positional arg
        "Lx": Lx, "Ly": Ly, "Lz": Lz,
        "lt": parse_lt(name),
        "n_actins": parse_actin_count(name),
        "n_frames": nf,
        "stride": stride,
        "x": x_axis,
        "snap_frames": snap_frames,
        "smectic_fraction": smectic_frac,
    }
    for k, v in ts.items():
        payload[f"ts_{k}"] = v
    for k, v in snaps.items():
        # scalar per-frame quantities stay 1D; field arrays become (n_snap, nbins)
        payload[f"snap_{k}"] = np.asarray(v)

    np.savez_compressed(out_npz, **payload)

    # console summary -- the headline numbers
    nx_tail = ts["nx_star"][len(ts["nx_star"]) // 2:]
    nx_tail = nx_tail[nx_tail > 0]
    print(f"{name}")
    print(f"  frames={nf}  lt={payload['lt']}  A={payload['n_actins']}  box={Lx}x{Ly}x{Lz}")
    if nx_tail.size:
        print(f"  steady nx* (median, second half) = {int(np.median(nx_tail))}"
              f"  -> period = {Lx/np.median(nx_tail):.3f}")
    print(f"  nx* early->late: "
          f"{_first_nonzero(ts['nx_star'])} -> {_last_nonzero(ts['nx_star'])}"
          f"   (decreasing => coarsening; flat => selection)")
    print(f"  smectic fraction = {smectic_frac:.2f}"
          f"   (->1 smectic / layered; <1 transverse order)")
    print(f"  saved {out_npz}")
    return payload


def _first_nonzero(a):
    nz = np.flatnonzero(a > 0)
    return int(a[nz[0]]) if nz.size else 0


def _last_nonzero(a):
    nz = np.flatnonzero(a > 0)
    return int(a[nz[-1]]) if nz.size else 0


# --------------------------------------------------------------------------- #
if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("h5", nargs="+", help="trajectory HDF5 file(s)")
    ap.add_argument("--out_dir", default=".", help="directory for .npz outputs")
    ap.add_argument("--stride", type=int, default=5,
                    help="frame stride for the q_dom(t) time series")
    ap.add_argument("--nbins", type=int, default=240,
                    help="x bins for the coarse-grained fields (dx=Lx/nbins)")
    ap.add_argument("--smooth", type=float, default=0.0,
                    help="Gaussian smoothing width for fields (in box units, 0=off)")
    ap.add_argument("--n_snapshots", type=int, default=6,
                    help="number of evenly spaced field snapshots in [start,end]")
    ap.add_argument("--start_frame", type=int, default=0,
                    help="first frame for the snapshot window (for PDE discovery, "
                         "start in the assembly transient, e.g. 300)")
    ap.add_argument("--end_frame", type=int, default=-1,
                    help="last frame for the snapshot window (-1 = last frame)")
    ap.add_argument("--population", default="cb2",
                    choices=("cb2", "all", "plus", "minus"),
                    help="which actin population defines q_dom (cb2 matches phase diagram)")
    ap.add_argument("--device", default="numpy", choices=("numpy", "cpu", "cuda", "auto"))
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    for path in args.h5:
        stem = os.path.splitext(os.path.basename(path))[0]
        out = os.path.join(args.out_dir, f"cg_{stem}.npz")
        try:
            process_file(path, out, stride=args.stride, nbins=args.nbins,
                         smooth=args.smooth, n_snapshots=args.n_snapshots,
                         population=args.population, device=args.device,
                         start_frame=args.start_frame, end_frame=args.end_frame)
        except Exception as e:  # keep going across a batch of files
            print(f"{os.path.basename(path)}  ERROR  {e}")
