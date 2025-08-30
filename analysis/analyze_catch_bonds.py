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


def plot_breakage_events(h5file: str, dt: float = 1.0, prefix: str = "analysis") -> None:
    """Read recorded catch-bond breakage events and plot distance, angle,
    tension and myosin-attachment metrics vs time."""
    with h5py.File(h5file, "r") as fh:
        if "/catch_bond/breakage" not in fh:
            print("No catch-bond breakage data found in file")
            return
        data = np.asarray(fh["/catch_bond/breakage"])
        act_vel_ds = fh["/actin/velocity"] if "/actin/velocity" in fh else None
        myo_vel_ds = fh["/myosin/velocity"] if "/myosin/velocity" in fh else None

        act_speeds: list[float] = []
        myo_speeds: list[float] = []
        if data.size:
            if act_vel_ds is not None:
                for row in data:
                    step_idx = int(row[2])
                    i = int(row[0])
                    j = int(row[1])
                    frame_vel = act_vel_ds[min(step_idx, act_vel_ds.shape[0] - 1)]
                    act_speeds.extend(
                        [
                            float(np.linalg.norm(frame_vel[i])),
                            float(np.linalg.norm(frame_vel[j])),
                        ]
                    )
            if myo_vel_ds is not None:
                max_myo = (data.shape[1] - 9) // 2
                for row in data:
                    step_idx = int(row[2])
                    frame_myo_vel = myo_vel_ds[min(step_idx, myo_vel_ds.shape[0] - 1)]
                    for idx in row[9:9 + max_myo]:
                        mi = int(idx)
                        if mi >= 0:
                            myo_speeds.append(float(np.linalg.norm(frame_myo_vel[mi])))
                    for idx in row[9 + max_myo:9 + 2 * max_myo]:
                        mi = int(idx)
                        if mi >= 0:
                            myo_speeds.append(float(np.linalg.norm(frame_myo_vel[mi])))

    if data.size == 0:
        print("No catch-bond breakage events recorded")
        return

    steps = data[:, 2]
    times = steps * dt
    distances = data[:, 3]
    angles = np.degrees(np.arccos(np.clip(data[:, 4], -1.0, 1.0)))
    tension_i = data[:, 5]
    tension_j = data[:, 6]
    min_tension = np.minimum(tension_i, tension_j)
    count_i = data[:, 7]
    count_j = data[:, 8]
    total_myo = count_i + count_j

    plt.figure()
    plt.scatter(times, distances, s=10, alpha=0.7)
    plt.xlabel("Time")
    plt.ylabel("Segment distance at break")
    plt.tight_layout()
    plt.savefig(f"{prefix}_cb_break_distance_vs_time.png", dpi=300)
    plt.close()

    plt.figure()
    plt.scatter(times, angles, s=10, alpha=0.7)
    plt.xlabel("Time")
    plt.ylabel("Angle at break (deg)")
    plt.tight_layout()
    plt.savefig(f"{prefix}_cb_break_angle_vs_time.png", dpi=300)
    plt.close()

    plt.figure()
    plt.scatter(times, min_tension, s=10, alpha=0.7)
    plt.xlabel("Time")
    plt.ylabel("Minimum tension at break")
    plt.tight_layout()
    plt.savefig(f"{prefix}_cb_break_min_tension_vs_time.png", dpi=300)
    plt.close()

    plt.figure()
    plt.scatter(times, total_myo, s=10, alpha=0.7)
    plt.xlabel("Time")
    plt.ylabel("Total myosins attached at break")
    plt.tight_layout()
    plt.savefig(f"{prefix}_cb_break_total_myosins_vs_time.png", dpi=300)
    plt.close()

    if act_speeds:
        plt.figure()
        plt.hist(act_speeds, bins=50, density=True)
        plt.xlabel("Actin speed before break")
        plt.ylabel("Probability density")
        plt.tight_layout()
        plt.savefig(f"{prefix}_cb_break_actin_speed_distribution.png", dpi=300)
        plt.close()

    if myo_speeds:
        plt.figure()
        plt.hist(myo_speeds, bins=50, density=True)
        plt.xlabel("Myosin speed before break")
        plt.ylabel("Probability density")
        plt.tight_layout()
        plt.savefig(f"{prefix}_cb_break_myosin_speed_distribution.png", dpi=300)
        plt.close()

    tensionless = np.sum(min_tension < 1e-6)
    detached = np.sum((count_i == 0) | (count_j == 0))
    print(f"{tensionless} of {len(times)} break events occurred with near-zero tension")
    print(f"{detached} events involved an actin with no myosin attachments")


def summarize_limit_removals(h5file: str) -> None:
    """Print a table summarizing limit-enforced catch-bond removals."""
    with h5py.File(h5file, "r") as fh:
        if "/catch_bond/limit_removal" not in fh:
            print("No limit-removal data found in file")
            return
        data = np.asarray(fh["/catch_bond/limit_removal"])

    if data.size == 0:
        print("No limit-removal events recorded")
        return

    header = f"{'step':>10} {'i':>5} {'j':>5} {'bonds_i':>8} {'bonds_j':>8}"
    print("Limit-enforced removals")
    print(header)
    for row in data:
        print(f"{int(row[2]):>10} {int(row[0]):>5} {int(row[1]):>5} {int(row[3]):>8} {int(row[4]):>8}")


def analyze_catch_bonds(h5file: str, dt: float = 1.0,
                        prefix: str = "analysis",
                        start_frame: int = 0) -> None:
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
        partial_binding_ds = fh["/actin/partial_binding_ratio"]

        actin_vel_ds = fh["/actin/velocity"]
        myosin_vel_ds = fh["/myosin/velocity"] if "/myosin/velocity" in fh else None
        myosin_fload_ds = fh["/myosin/f_load"] if "/myosin/f_load" in fh else None

        # Optional datasets
        am_bonds_ds = fh["/actin_myo/bonds"] if "/actin_myo/bonds" in fh else None
        n_myosins_total = (
            myosin_vel_ds.shape[1]
            if myosin_vel_ds is not None
            else (myosin_fload_ds.shape[1] if myosin_fload_ds is not None else 0)
        )
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

        # Nematic order per frame
        nematic_order: list[float] = []

        # Distributions of actin--myosin connectivity
        myosins_per_actin: list[int] = []
        actins_per_myosin: list[int] = []

        actin_speeds: list[float] = []
        actin_load_vals: list[float] = []
        myosin_speeds: list[float] = []
        myosin_load_vals: list[float] = []

        # Track actins that report cb status 2 but zero partial binding ratio
        cb2_zero_partial: list[tuple[int, np.ndarray]] = []

        for frame in range(n_frames):
            if frame < start_frame:
                continue
            bonds = bonds_ds[frame]
            dirs = np.asarray(dirs_ds[frame])  # (N,3)
            f_load = fload_ds[frame, :, 0]     # (N,)
            cb_strength_frame = cb_status_actin[frame, :, 0]
            partial_binding = partial_binding_ds[frame, :, 0]

            cb2_indices = np.where(cb_strength_frame == 2)[0]
            zero_partial = cb2_indices[partial_binding[cb2_indices] <= 0]
            if zero_partial.size:
                cb2_zero_partial.append((frame, zero_partial))

            actin_speed = np.linalg.norm(actin_vel_ds[frame], axis=1)
            actin_speeds.extend(actin_speed)
            actin_load_vals.extend(f_load)
            if myosin_vel_ds is not None and myosin_fload_ds is not None:
                myosin_speed = np.linalg.norm(myosin_vel_ds[frame], axis=1)
                myosin_speeds.extend(myosin_speed)
                myosin_load_vals.extend(myosin_fload_ds[frame, :, 0])

            # Accumulate raw directions for global distribution
            all_dirs.append(dirs)

            # Nematic order parameter for this frame
            if dirs.size:
                norms = np.linalg.norm(dirs, axis=1, keepdims=True)
                norms[norms == 0] = 1.0
                U = dirs / norms
                q = 1.5 * (U.T @ U) / U.shape[0] - 0.5 * np.eye(3)
                eigvals = np.linalg.eigvalsh(q)
                nematic_order.append(float(eigvals[-1]))

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

                # Myosins attached to actins with cb_status == 2 are catch bonded
                cb_actins = np.where(cb_strength_frame == 2)[0]
                catch_myosins: set[int] = set()
                for a in cb_actins:
                    catch_myosins.update(a2m.get(int(a), set()))
                ratio_myosin_cb.append(
                    len(catch_myosins) / max(n_myosins_total, 1)
                )

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

        if cb2_zero_partial:
            for frame_idx, indices in cb2_zero_partial:
                print(
                    f"Frame {frame_idx}: actins {indices.tolist()} have cb_status 2"
                    " but zero partial binding ratio"
                )

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

    # Time series of nematic order parameter
    if nematic_order:
        times = np.arange(len(nematic_order)) * dt
        plt.figure()
        plt.plot(times, nematic_order)
        plt.xlabel("Time")
        plt.ylabel("Nematic order")
        plt.tight_layout()
        plt.savefig(f"{prefix}_nematic_order_vs_time.png", dpi=300)
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
        plt.hist(myosins_per_actin, bins=bins)
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
    parser.add_argument("--start_frame", type=int, default=0,help="First frame to include in analysis (skip earlier frames)")

    args = parser.parse_args()
    analyze_catch_bonds(args.h5file,
                        dt=args.dt,
                        prefix=args.prefix,
                        start_frame=args.start_frame)
    plot_breakage_events(args.h5file, dt=args.dt, prefix=args.prefix)
    summarize_limit_removals(args.h5file)

if __name__ == "__main__":
    main()
