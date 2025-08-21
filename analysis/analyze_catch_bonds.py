import os
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

import argparse
import h5py
import numpy as np
import matplotlib.pyplot as plt


def compute_pair_fload(
    dir_i: np.ndarray,
    dir_j: np.ndarray,
    f_i: float,
    f_j: float,
) -> tuple[float, float]:
    """Return (f_load, angle_deg) for a bonded actin pair.

    Parameters
    ----------
    dir_i, dir_j : ndarray
        Unit direction vectors of actins i and j.
    f_i, f_j : float
        Scalar load values for actins i and j.
    """
    # Ensure unit length (robust to small numeric drift)
    di = dir_i / max(np.linalg.norm(dir_i), 1e-12)
    dj = dir_j / max(np.linalg.norm(dir_j), 1e-12)

    cos_val = float(np.clip(np.dot(di, dj), -1.0, 1.0))
    # Only anti-parallel contributes (cos<0); use the smaller f as the pair load
    f_load = abs(cos_val) * min(f_i, f_j) * (cos_val < 0)
    angle = np.degrees(np.arccos(cos_val))
    return f_load, angle


def analyze_catch_bonds(h5file: str, dt: float = 1.0, prefix: str = "analysis") -> None:
    """Analyze actin catch bonds and actin--myosin connectivity in a trajectory file.

    Parameters
    ----------
    h5file : str
        Path to the HDF5 trajectory produced by the simulation.
    dt : float, optional
        Time between stored frames (converts lifetimes from frames to time).
    prefix : str, optional
        Prefix for output files.
    """
    with h5py.File(h5file, "r") as fh:
        bonds_ds = fh["/actin/bonds"]
        dirs_ds = fh["/actin/direction"]
        fload_ds = fh["/actin/f_load"]
        cb_status_actin = fh["/actin/cb_status"]

        actin_vel_ds = fh["/actin/velocity"]
        myosin_vel_ds = fh["/myosin/velocity"] if "/myosin/velocity" in fh else None
        myosin_fload_ds = fh["/myosin/f_load"] if "/myosin/f_load" in fh else None

        # Optional myosin datasets
        myosin_bonds_ds = fh["/myosin/bonds"] if "/myosin/bonds" in fh else None
        am_bonds_ds = fh["/actin_myo/bonds"] if "/actin_myo/bonds" in fh else None
        n_frames = bonds_ds.shape[0]

        active: dict[tuple[int, int], dict] = {}
        lifetimes: list[float] = []
        mean_floads: list[float] = []
        angles_pairwise: list[float] = []

        # Collect all actin directions across frames for a global distribution
        all_dirs = []

        # Ratios of filaments engaged in catch bonds per frame
        ratio_actin_cb: list[float] = []
        ratio_myosin_cb: list[float] = []

        # Distributions of actin--myosin connectivity
        myosins_per_actin: list[int] = []
        actins_per_myosin: list[int] = []

        actin_speeds: list[float] = []
        actin_load_vals: list[float] = []
        myosin_speeds: list[float] = []
        myosin_load_vals: list[float] = []

        for frame in range(n_frames):
            bonds = bonds_ds[frame]
            dirs = np.asarray(dirs_ds[frame])  # (N,3)
            f_load = fload_ds[frame, :, 0]     # (N,)
            cb_strength_frame = cb_status_actin[frame, :, 0]

            actin_speed = np.linalg.norm(actin_vel_ds[frame], axis=1)
            actin_speeds.extend(actin_speed)
            actin_load_vals.extend(f_load)
            if myosin_vel_ds is not None and myosin_fload_ds is not None:
                myosin_speed = np.linalg.norm(myosin_vel_ds[frame], axis=1)
                myosin_speeds.extend(myosin_speed)
                myosin_load_vals.extend(myosin_fload_ds[frame, :, 0])

            # Accumulate raw directions for global distribution
            all_dirs.append(dirs)

            bonded_actins: set[int] = set()
            current_pairs: set[tuple[int, int]] = set()
            for pair in bonds:
                a, b = int(pair[0]), int(pair[1])
                if a < 0 or b < 0:
                    continue
                if a > b:
                    a, b = b, a
                current_pairs.add((a, b))
                bonded_actins.update([a, b])

                pair_fload, angle = compute_pair_fload(
                    dirs[a], dirs[b], float(f_load[a]), float(f_load[b])
                )
                angles_pairwise.append(angle)

                if (a, b) in active:
                    entry = active[(a, b)]
                    entry["last"] = frame
                    entry["sum_fload"] += pair_fload
                    entry["count"] += 1
                else:
                    active[(a, b)] = {
                        "start": frame,
                        "last": frame,
                        "sum_fload": pair_fload,
                        "count": 1,
                    }

            # Ratio of actins engaged in catch bonds (cb_status > 0)
            n_actins = cb_strength_frame.shape[0]
            catch_actins = [i for i in bonded_actins if cb_strength_frame[i] > 0]
            ratio_actin_cb.append(len(catch_actins) / max(n_actins, 1))

            # Actin--myosin connectivity for this frame
            if am_bonds_ds is not None:
                am_pairs = am_bonds_ds[frame]
                a2m: dict[int, set[int]] = {}
                m2a: dict[int, set[int]] = {}
                for pair in am_pairs:
                    a, m = int(pair[0]), int(pair[1])
                    if a < 0 or m < 0:
                        continue
                    a2m.setdefault(a, set()).add(m)
                    m2a.setdefault(m, set()).add(a)

                myosins_per_actin.extend(len(v) for v in a2m.values())
                actins_per_myosin.extend(len(v) for v in m2a.values())

            # Ratio of myosins engaged in catch bond
            if myosin_bonds_ds is not None:
                myo_bonds = myosin_bonds_ds[frame]
                n_myosins = 100
                bonded_myosins: set[int] = set()
                for pair in myo_bonds:
                    m, n = int(pair[0]), int(pair[1])
                    if m < 0 or n < 0:
                        continue
                    bonded_myosins.update([m, n])
                catch_myosins = [i for i in bonded_myosins]
                ratio_myosin_cb.append(len(catch_myosins) / max(n_myosins, 1))

            # Close out bonds that ended this frame
            ended = [p for p in active if p not in current_pairs]
            for p in ended:
                entry = active.pop(p)
                lifetime = (entry["last"] - entry["start"]) * dt
                mean_fload = entry["sum_fload"] / max(entry["count"], 1)
                lifetimes.append(lifetime)
                mean_floads.append(mean_fload)

        # Close out bonds that persist to the final frame
        for entry in active.values():
            lifetime = (entry["last"] - entry["start"] + 1) * dt
            mean_fload = entry["sum_fload"] / max(entry["count"], 1)
            lifetimes.append(lifetime)
            mean_floads.append(mean_fload)

    if actin_speeds:
        plt.figure()
        plt.scatter(actin_load_vals, actin_speeds, s=10, alpha=0.7)
        plt.xlabel("f_load")
        plt.ylabel("Actin speed")
        plt.tight_layout()
        plt.savefig(f"{prefix}_actin_velocity_vs_load.png", dpi=300)
        plt.close()

    if myosin_speeds:
        plt.figure()
        plt.scatter(myosin_load_vals, myosin_speeds, s=10, alpha=0.7)
        plt.xlabel("f_load")
        plt.ylabel("Myosin speed")
        plt.tight_layout()
        plt.savefig(f"{prefix}_myosin_velocity_vs_load.png", dpi=300)
        plt.close()

    # Plot lifetime vs load for bonded pairs
    if lifetimes:
        plt.figure()
        plt.scatter(mean_floads, lifetimes, s=10, alpha=0.7)
        plt.xlabel("f_load")
        plt.ylabel("Lifetime")
        plt.tight_layout()
        plt.savefig(f"{prefix}_lifetime_vs_fload.png", dpi=300)
        plt.close()

    # Pairwise angle distribution (between bonded actin directions)
    if angles_pairwise:
        plt.figure()
        plt.hist(angles_pairwise, bins=50, density=True)
        plt.xlabel("Angle between bonded actins (degrees)")
        plt.ylabel("Probability density")
        plt.tight_layout()
        plt.savefig(f"{prefix}_bonded_pair_angle_distribution.png", dpi=300)
        plt.close()

    # Global actin direction distribution (angle to +x axis across all frames)
    if all_dirs:
        all_dirs = np.vstack(all_dirs)  # (total_actins, 3)
        # Normalize to be safe
        norms = np.linalg.norm(all_dirs, axis=1, keepdims=True)
        norms[norms == 0] = 1.0
        U = all_dirs / norms
        # Angle to +x axis
        cos_x = np.clip(U[:, 0], -1.0, 1.0)
        angles_to_x = np.degrees(np.arccos(cos_x))  # 0°=+x, 180°=-x

        plt.figure()
        plt.hist(angles_to_x, bins=72, density=True)  # 2.5° bins over [0,180]
        plt.xlabel("Actin direction: angle to +x (degrees)")
        plt.ylabel("Probability density")
        plt.tight_layout()
        plt.savefig(f"{prefix}_actin_direction_angle_to_x.png", dpi=300)
        plt.close()

    # Time series of catch bond engagement ratios
    if ratio_actin_cb:
        times = np.arange(len(ratio_actin_cb)) * dt
        plt.figure()
        plt.plot(times, ratio_actin_cb)
        plt.xlabel("Time")
        plt.ylabel("Actins in catch bond (ratio)")
        plt.tight_layout()
        plt.savefig(f"{prefix}_actin_catch_bond_ratio.png", dpi=300)
        plt.close()

    if ratio_myosin_cb:
        times = np.arange(len(ratio_myosin_cb)) * dt
        plt.figure()
        plt.plot(times, ratio_myosin_cb)
        plt.xlabel("Time")
        plt.ylabel("Myosins in catch bond (ratio)")
        plt.tight_layout()
        plt.savefig(f"{prefix}_myosin_catch_bond_ratio.png", dpi=300)
        plt.close()

    # Distributions of actin--myosin connections
    if myosins_per_actin:
        plt.figure()
        bins = np.arange(1, max(myosins_per_actin) + 2) - 0.5
        plt.hist(myosins_per_actin, bins=bins, density=True)
        plt.xlabel("Myosins bound per actin")
        plt.ylabel("Probability density")
        plt.tight_layout()
        plt.savefig(f"{prefix}_myosins_per_actin_distribution.png", dpi=300)
        plt.close()

    if actins_per_myosin:
        plt.figure()
        bins = np.arange(1, max(actins_per_myosin) + 2) - 0.5
        plt.hist(actins_per_myosin, bins=bins, density=True)
        plt.xlabel("Actins bound per myosin")
        plt.ylabel("Probability density")
        plt.tight_layout()
        plt.savefig(f"{prefix}_actins_per_myosin_distribution.png", dpi=300)
        plt.close()


def main() -> None:
    parser = argparse.ArgumentParser(description="Analyze actin catch bonds from trajectory file")
    parser.add_argument("h5file", help="Path to HDF5 trajectory")
    parser.add_argument("--dt", type=float, default=0.02, help="Time between frames")
    parser.add_argument("--prefix", default="analysis", help="Prefix for output files")
    args = parser.parse_args()
    analyze_catch_bonds(args.h5file, dt=args.dt, prefix=args.prefix)


if __name__ == "__main__":
    main()
