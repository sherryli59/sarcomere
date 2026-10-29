"""Render one measured Z-disc patch / myosin-tip overshoot from a trajectory.

The HDF5 file records the winning (largest-force) ACTN2 patch contact for each
myosin tip, but not the winning link ID.  This script selects a severe recorded
contact, reconstructs every live ACTN2 patch in that frame, and requires the
best matching patch to reproduce the saved values within the sub-step motion
between force evaluation and the stored pose.
"""

import argparse
import json
import os

os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
os.environ.setdefault("MPLCONFIGDIR", f"/tmp/codex-matplotlib-{os.getuid()}")

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def unit(v):
    n = np.linalg.norm(v)
    if n <= 1e-15:
        raise ValueError("cannot normalize a zero vector")
    return v / n


def min_image(d, box, periodic):
    out = np.asarray(d, dtype=float).copy()
    for k in range(3):
        if periodic[k]:
            out[k] -= box[k] * np.round(out[k] / box[k])
    return out


def radial_gate(rho, flat, cutoff):
    if rho <= flat:
        return 1.0
    q = (rho - flat) / (cutoff - flat)
    return 1.0 - 3.0 * q * q + 2.0 * q * q * q


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("input", help="projection-off trajectory HDF5")
    p.add_argument("output_png")
    p.add_argument("--output-json", default=None)
    p.add_argument(
        "--min-force",
        type=float,
        default=1.0,
        help="minimum recorded contact force in pN when choosing the deepest example",
    )
    return p.parse_args()


def main():
    args = parse_args()
    out_json = args.output_json or os.path.splitext(args.output_png)[0] + ".json"

    with h5py.File(args.input, "r", locking=False) as f:
        a = dict(f["meta"].attrs)
        n_frames = int(f["meta/committed_frames"][0])
        dt = float(a["dt"])
        box = np.array([a["box_x"], a["box_y"], a["box_z"]], dtype=float)
        periodic = np.array(
            [a.get("periodic_x", 1), a.get("periodic_y", 1), a.get("periodic_z", 1)],
            dtype=bool,
        )
        la = float(a["actin_length"])
        lm = float(a["myosin_length"])
        ra = float(a.get("actin_radius", 0.015))
        rm = float(a.get("myosin_radius", 0.05))
        onset = float(a["aa_zdisc_tip_onset"])
        stiffness = float(a["aa_zdisc_tip_stiffness"])
        force_max = float(a["aa_zdisc_tip_force_max"])
        transverse_flat = float(a["aa_zdisc_tip_transverse_flat"])
        transverse_cutoff = float(a["aa_zdisc_tip_transverse_cutoff"])

        # Select the deepest saved contact that still carries a meaningful force.
        candidates = []
        all_depths = []
        for end_name, end_sign in (("minus", -1), ("plus", 1)):
            gap = f[f"myosin/aa_zdisc_tip_gap_{end_name}"][:n_frames, :, 0]
            force = f[f"myosin/aa_zdisc_tip_force_{end_name}"][:n_frames, :, 0]
            rho = f[f"myosin/aa_zdisc_tip_rho_{end_name}"][:n_frames, :, 0]
            penetrating = (force > 0.0) & (gap < -1e-3)
            all_depths.append(-1000.0 * gap[penetrating])
            eligible = penetrating & (force >= args.min_force)
            for frame, myosin in np.argwhere(eligible):
                candidates.append(
                    (
                        float(gap[frame, myosin]),
                        end_name,
                        end_sign,
                        int(frame),
                        int(myosin),
                        float(rho[frame, myosin]),
                        float(force[frame, myosin]),
                    )
                )
        if not candidates:
            raise RuntimeError(f"no penetrating contact has force >= {args.min_force:g} pN")
        saved_gap, end_name, end_sign, frame, myosin, saved_rho, saved_force = min(
            candidates, key=lambda row: row[0]
        )
        all_depths = np.concatenate(all_depths)

        step = int(f["state/current_step"][frame, 0])
        ca = f["actin/center"][frame]
        ua = f["actin/direction"][frame]
        ua = ua / np.linalg.norm(ua, axis=1, keepdims=True)
        cm = f["myosin/center"][frame, myosin].astype(float)
        um = unit(f["myosin/direction"][frame, myosin].astype(float))
        outward = end_sign * um
        tip = cm + end_sign * 0.5 * lm * um

        offset = f["catch_bond/local_link_frame_offsets"][frame]
        start, count = int(offset[1]), int(offset[2])
        identities = f["catch_bond/local_link_identity"][start : start + count]
        coordinates = f["catch_bond/local_link_coordinates"][start : start + count]

        reconstructed = []
        for local_index, (identity, coord) in enumerate(zip(identities, coordinates)):
            i, j = int(identity[2]), int(identity[3])
            si, sj = float(coord[0]), float(coord[1])
            first = ca[i] + si * ua[i]
            second_raw = ca[j] + sj * ua[j]
            connector = min_image(second_raw - first, box, periodic)
            second = first + connector
            patch = first + 0.5 * connector
            center_to_patch = min_image(patch - cm, box, periodic)
            if end_sign * float(center_to_patch @ um) < 0.0:
                continue
            tip_to_patch = min_image(patch - tip, box, periodic)
            gap = float(tip_to_patch @ outward)
            transverse = tip_to_patch - gap * outward
            rho = float(np.linalg.norm(transverse))
            if rho >= transverse_cutoff or gap >= onset:
                continue
            gate = radial_gate(rho, transverse_flat, transverse_cutoff)
            force = gate * force_max * np.tanh(stiffness * (onset - gap) / force_max)
            reconstructed.append(
                dict(
                    force=float(force),
                    gap=gap,
                    rho=rho,
                    local_index=local_index,
                    identity=identity.astype(int),
                    i=i,
                    j=j,
                    si=si,
                    sj=sj,
                    first=first,
                    second=second,
                    patch=patch,
                    connector=connector,
                    transverse=transverse,
                )
            )
        if not reconstructed:
            raise RuntimeError("no live link reconstructed the saved contact")
        # The production selection scans only each linked actin's spatial
        # myosin-neighbor list.  That neighbor list is not written to HDF5, so
        # identify the saved winner by its recorded (gap, rho, force) triple.
        # The match is unique at floating-point precision in this frame.
        winner = min(
            reconstructed,
            key=lambda row: max(
                abs(row["gap"] - saved_gap) / 1e-3,
                abs(row["rho"] - saved_rho) / 1e-3,
                abs(row["force"] - saved_force) / 0.2,
            ),
        )
        errors = {
            "gap_um": abs(winner["gap"] - saved_gap),
            "rho_um": abs(winner["rho"] - saved_rho),
            "force_pN": abs(winner["force"] - saved_force),
        }
        if (
            errors["gap_um"] > 1e-3
            or errors["rho_um"] > 1e-3
            or errors["force_pN"] > 0.2
        ):
            raise RuntimeError(
                "reconstructed winner does not match saved contact: " + repr(errors)
            )

        # Unwrap both actins consistently around their connector, then express
        # the measured geometry in a local basis at the selected myosin tip.
        i, j = winner["i"], winner["j"]
        first, second, patch = winner["first"], winner["second"], winner["patch"]
        ci = first - winner["si"] * ua[i]
        cj = second - winner["sj"] * ua[j]
        er = unit(winner["transverse"])
        et = unit(np.cross(outward, er))

        def local(p):
            d = min_image(np.asarray(p) - tip, box, periodic)
            return 1000.0 * np.array([d @ outward, d @ er, d @ et])

        actin_i = np.vstack([ci - 0.5 * la * ua[i], ci + 0.5 * la * ua[i]])
        actin_j = np.vstack([cj - 0.5 * la * ua[j], cj + 0.5 * la * ua[j]])
        ai_local = np.vstack([local(p) for p in actin_i])
        aj_local = np.vstack([local(p) for p in actin_j])
        first_local, second_local, patch_local = local(first), local(second), local(patch)
        myosin_center_local = np.array([-500.0 * lm, 0.0, 0.0])
        myosin_far_local = np.array([-1000.0 * lm, 0.0, 0.0])

        fig = plt.figure(figsize=(17.0, 5.6), constrained_layout=True)
        gs = fig.add_gridspec(1, 3, width_ratios=(1.13, 1.28, 0.92))

        # A: actual local 3-D pose.
        ax3 = fig.add_subplot(gs[0, 0], projection="3d")
        ax3.plot(
            [myosin_far_local[0], 0], [0, 0], [0, 0],
            color="#40464f", lw=9, solid_capstyle="round", label="myosin",
        )
        ax3.plot(*ai_local.T, color="#1f9e89", lw=5, label=f"actin {i}")
        ax3.plot(*aj_local.T, color="#e76f51", lw=5, label=f"actin {j}")
        ax3.plot(
            [first_local[0], second_local[0]],
            [first_local[1], second_local[1]],
            [first_local[2], second_local[2]],
            color="#7b2cbf", lw=4, label="ACTN2 link",
        )
        ax3.scatter(*patch_local, s=95, color="#ffd166", edgecolor="black", zorder=10)
        ax3.scatter(0, 0, 0, s=85, marker="D", color="#d62828", label="myosin tip")
        ax3.plot([saved_gap * 1000, 0], [saved_rho * 1000] * 2, [0, 0], "--", color="#d62828", lw=1.5)
        ax3.set_xlim(-900, 170)
        ax3.set_ylim(-260, 300)
        ax3.set_zlim(-260, 260)
        ax3.set_xlabel("outward axis from tip (nm)")
        ax3.set_ylabel("radial coordinate (nm)")
        ax3.set_zlabel("orthogonal coordinate (nm)")
        ax3.view_init(elev=23, azim=-58)
        ax3.set_title("A  Actual frame geometry", loc="left", fontweight="bold")
        ax3.legend(loc="upper left", fontsize=8, frameon=False)

        # B: longitudinal projection makes the two measured distances explicit.
        ax = fig.add_subplot(gs[0, 1])
        ax.fill_between(
            [-900, 0], -1000 * rm, 1000 * rm,
            color="#59636f", alpha=0.55, label=f"myosin body (radius {1000*rm:.0f} nm)",
        )
        ax.axhspan(-1000 * transverse_cutoff, 1000 * transverse_cutoff, color="#4cc9f0", alpha=0.07)
        ax.axhline(1000 * transverse_flat, color="#4cc9f0", ls=":", lw=1)
        ax.axhline(1000 * transverse_cutoff, color="#4cc9f0", ls="--", lw=1)
        ax.axvline(0, color="#d62828", ls="--", lw=1.5)
        ax.plot(ai_local[:, 0], ai_local[:, 1], color="#1f9e89", lw=5)
        ax.plot(aj_local[:, 0], aj_local[:, 1], color="#e76f51", lw=5)
        ax.plot(
            [first_local[0], second_local[0]], [first_local[1], second_local[1]],
            color="#7b2cbf", lw=4,
        )
        ax.scatter(patch_local[0], patch_local[1], s=105, color="#ffd166", edgecolor="black", zorder=6)
        ax.scatter(0, 0, s=90, marker="D", color="#d62828", zorder=6)
        ax.annotate(
            "", xy=(0, patch_local[1] + 27), xytext=(patch_local[0], patch_local[1] + 27),
            arrowprops=dict(arrowstyle="<->", color="#d62828", lw=2),
        )
        ax.text(
            0.5 * patch_local[0], patch_local[1] + 42,
            f"overshoot = {-1000*saved_gap:.1f} nm", color="#a4161a", ha="center", fontsize=10,
        )
        ax.annotate(
            "", xy=(patch_local[0] - 25, 0), xytext=(patch_local[0] - 25, patch_local[1]),
            arrowprops=dict(arrowstyle="<->", color="#277da1", lw=2),
        )
        ax.text(
            patch_local[0] - 38, 0.5 * patch_local[1], f"$\\rho$ = {1000*saved_rho:.1f} nm",
            color="#1d6a8a", rotation=90, va="center", ha="right", fontsize=10,
        )
        ax.text(8, -37, "tip plane", color="#a4161a", fontsize=9)
        ax.text(-880, 117, "100 nm flat zone", color="#277da1", fontsize=8)
        ax.text(-880, 137, "130 nm cutoff", color="#277da1", fontsize=8)
        ax.set_xlim(-900, 170)
        ax.set_ylim(-230, 265)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("outward axis from selected myosin tip (nm)")
        ax.set_ylabel("radial coordinate (nm)")
        ax.set_title("B  Longitudinal projection", loc="left", fontweight="bold")
        ax.text(
            0.02, 0.02,
            f"Recorded soft force: {saved_force:.1f} pN\n"
            f"Patch is outside the {1000*rm:.0f} nm myosin radius;\n"
            "the overshoot is axial, not radial burial.",
            transform=ax.transAxes, fontsize=9, va="bottom",
            bbox=dict(boxstyle="round,pad=0.35", facecolor="white", alpha=0.9, edgecolor="#bbbbbb"),
        )

        # C: place the example in the measured trajectory-wide distribution.
        ah = fig.add_subplot(gs[0, 2])
        bins = np.linspace(0, max(800.0, np.ceil(all_depths.max() / 50.0) * 50.0), 44)
        ah.hist(all_depths, bins=bins, color="#8ecae6", edgecolor="white", lw=0.5)
        p50, p90, pmax = np.percentile(all_depths, 50), np.percentile(all_depths, 90), all_depths.max()
        for x, color, label in (
            (p50, "#2a9d8f", f"median {p50:.0f} nm"),
            (p90, "#f4a261", f"p90 {p90:.0f} nm"),
            (-1000 * saved_gap, "#d62828", f"shown {(-1000*saved_gap):.0f} nm"),
        ):
            ah.axvline(x, color=color, lw=2, label=label)
        ah.set_yscale("log")
        ah.set_xlabel("patch overshoot past myosin tip (nm)")
        ah.set_ylabel("selected tip-contacts")
        ah.set_title("C  Full 2.0–2.4 s distribution", loc="left", fontweight="bold")
        ah.legend(frameon=False, fontsize=9)
        ah.text(
            0.98, 0.97, f"maximum: {pmax:.1f} nm\n$n$ = {len(all_depths):,}",
            transform=ah.transAxes, ha="right", va="top", fontsize=9,
        )

        fig.suptitle(
            f"Projection-off trajectory: a force-bearing myosin-tip overshoot\n"
            f"r27zt0np · frame {frame} · step {step} · {step*dt:.3f} s · myosin {myosin} {end_name} tip",
            fontsize=15, fontweight="bold",
        )
        os.makedirs(os.path.dirname(os.path.abspath(args.output_png)), exist_ok=True)
        fig.savefig(args.output_png, dpi=220, facecolor="white")
        plt.close(fig)

        manifest = {
            "input": os.path.abspath(args.input),
            "selection": {
                "rule": "deepest penetrating saved tip contact with recorded force >= min_force_pN",
                "min_force_pN": args.min_force,
            },
            "frame": frame,
            "step": step,
            "time_s": step * dt,
            "myosin": myosin,
            "end": end_name,
            "link_identity_row": winner["identity"].tolist(),
            "actins": [i, j],
            "link_status": int(winner["identity"][4]),
            "link_coordinates_um": [winner["si"], winner["sj"]],
            "saved_contact": {
                "axial_gap_um": saved_gap,
                "overshoot_nm": -1000.0 * saved_gap,
                "rho_um": saved_rho,
                "rho_nm": 1000.0 * saved_rho,
                "force_pN": saved_force,
            },
            "reconstruction_errors": errors,
            "reconstruction_tolerance": {
                "gap_um": 1e-3,
                "rho_um": 1e-3,
                "force_pN": 0.2,
                "note": "allows the sub-step motion between contact-force evaluation and the stored pose",
            },
            "distribution_nm": {
                "n": int(len(all_depths)),
                "median": float(p50),
                "p90": float(p90),
                "max": float(pmax),
            },
            "model_geometry_um": {
                "myosin_radius": rm,
                "actin_radius": ra,
                "soft_contact_transverse_flat": transverse_flat,
                "soft_contact_transverse_cutoff": transverse_cutoff,
                "soft_contact_axial_onset": onset,
            },
            "output_png": os.path.abspath(args.output_png),
        }
        with open(out_json, "w", encoding="utf-8") as handle:
            json.dump(manifest, handle, indent=2)
            handle.write("\n")
        print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
