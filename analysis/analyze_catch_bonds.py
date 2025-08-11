import os
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

import argparse
import h5py
import numpy as np
import matplotlib.pyplot as plt


def compute_pair_fload(dir_i: np.ndarray, dir_j: np.ndarray, f_i: float, f_j: float) -> tuple[float, float]:
    """Return (f_load, cos_angle) for a bonded actin pair.

    Parameters
    ----------
    dir_i, dir_j : ndarray
        Unit direction vectors of actins *i* and *j*.
    f_i, f_j : float
        Scalar load values for actins *i* and *j*.
    """
    cos_angle = np.abs(np.dot(dir_i, dir_j))
    f_load = cos_angle * min(f_i, f_j)
    return f_load, cos_angle


def analyze_catch_bonds(h5file: str, dt: float = 1.0, prefix: str = "analysis") -> None:
    """Analyze actin catch bonds in a trajectory file.

    Parameters
    ----------
    h5file : str
        Path to the HDF5 trajectory produced by the simulation.
    dt : float, optional
        Time between stored frames.  Used to convert bond lifetimes from
        frames to physical time units.  Defaults to 1.0.
    prefix : str, optional
        Prefix for output files.  Default is ``"analysis"``.
    """
    with h5py.File(h5file, "r") as fh:
        bonds_ds = fh["/actin/bonds"]
        dirs_ds = fh["/actin/direction"]
        fload_ds = fh["/actin/f_load"]
        n_frames = bonds_ds.shape[0]

        active: dict[tuple[int, int], dict] = {}
        lifetimes: list[float] = []
        mean_floads: list[float] = []
        angles: list[float] = []

        for frame in range(n_frames):
            bonds = bonds_ds[frame]
            dirs = dirs_ds[frame]
            f_load = fload_ds[frame, :, 0]

            current_pairs: set[tuple[int, int]] = set()
            for pair in bonds:
                a, b = int(pair[0]), int(pair[1])
                if a < 0 or b < 0:
                    continue
                if a > b:
                    a, b = b, a
                current_pairs.add((a, b))

                pair_fload, cos_angle = compute_pair_fload(dirs[a], dirs[b], f_load[a], f_load[b])
                angle = np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0)))
                angles.append(angle)

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

            ended = [p for p in active if p not in current_pairs]
            for p in ended:
                entry = active.pop(p)
                lifetime = (entry["last"] - entry["start"] + 1) * dt
                mean_fload = entry["sum_fload"] / entry["count"]
                lifetimes.append(lifetime)
                mean_floads.append(mean_fload)

        for entry in active.values():
            lifetime = (entry["last"] - entry["start"] + 1) * dt
            mean_fload = entry["sum_fload"] / entry["count"]
            lifetimes.append(lifetime)
            mean_floads.append(mean_fload)

    if lifetimes:
        plt.figure()
        plt.scatter(mean_floads, lifetimes, s=10, alpha=0.7)
        plt.xlabel("f_load")
        plt.ylabel("Lifetime")
        plt.tight_layout()
        plt.savefig(f"{prefix}_lifetime_vs_fload.png", dpi=300)
        plt.close()

    if angles:
        plt.figure()
        plt.hist(angles, bins=50, density=True)
        plt.xlabel("Angle between actins (degrees)")
        plt.ylabel("Probability density")
        plt.tight_layout()
        plt.savefig(f"{prefix}_angle_distribution.png", dpi=300)
        plt.close()

    if lifetimes:
        arr = np.column_stack([mean_floads, lifetimes])
        np.savetxt(f"{prefix}_lifetime_vs_fload.dat", arr, header="f_load lifetime")

    if angles:
        np.savetxt(f"{prefix}_angles.dat", angles, header="angle_degrees")


def main() -> None:
    parser = argparse.ArgumentParser(description="Analyze actin catch bonds from trajectory file")
    parser.add_argument("h5file", help="Path to HDF5 trajectory")
    parser.add_argument("--dt", type=float, default=1.0, help="Time between frames")
    parser.add_argument("--prefix", default="analysis", help="Prefix for output files")
    args = parser.parse_args()
    analyze_catch_bonds(args.h5file, dt=args.dt, prefix=args.prefix)


if __name__ == "__main__":
    main()
