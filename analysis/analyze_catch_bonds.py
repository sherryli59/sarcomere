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
    """Analyze actin catch bonds in a trajectory file.

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
        n_frames = bonds_ds.shape[0]

        active: dict[tuple[int, int], dict] = {}
        lifetimes: list[float] = []
        mean_floads: list[float] = []
        angles_pairwise: list[float] = []

        # Collect all actin directions across frames for a global distribution
        all_dirs = []

        for frame in range(n_frames):
            bonds = bonds_ds[frame]
            dirs = np.asarray(dirs_ds[frame])  # (N,3)
            f_load = fload_ds[frame, :, 0]     # (N,)

            # Accumulate raw directions for global distribution
            all_dirs.append(dirs)

            current_pairs: set[tuple[int, int]] = set()
            for pair in bonds:
                a, b = int(pair[0]), int(pair[1])
                if a < 0 or b < 0:
                    continue
                if a > b:
                    a, b = b, a
                current_pairs.add((a, b))

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


def main() -> None:
    parser = argparse.ArgumentParser(description="Analyze actin catch bonds from trajectory file")
    parser.add_argument("h5file", help="Path to HDF5 trajectory")
    parser.add_argument("--dt", type=float, default=0.02, help="Time between frames")
    parser.add_argument("--prefix", default="analysis", help="Prefix for output files")
    args = parser.parse_args()
    analyze_catch_bonds(args.h5file, dt=args.dt, prefix=args.prefix)


if __name__ == "__main__":
    main()
