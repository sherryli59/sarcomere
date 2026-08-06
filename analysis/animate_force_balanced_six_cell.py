#!/usr/bin/env python3
"""Clear 3D + four-half transverse animation of the 24-actin periodic cell."""
import argparse
import os
from pathlib import Path
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/sarcomere_matplotlib")
import h5py
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter
from matplotlib.patches import Circle
import numpy as np

COLORS = ("#0072B2", "#D55E00", "#56B4E9", "#E69F00")
LABELS = ("M0 −half", "M1 +half", "M0 +half", "M1 −half")
OWNERS = (0, 1, 0, 1)


def rod(ax, center, direction, length, color, width):
    ends = np.vstack((center-.5*length*direction,
                      center+.5*length*direction))
    ax.plot(*ends.T, color=color, lw=width, solid_capstyle="round")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("trajectory", type=Path)
    p.add_argument("output", type=Path)
    p.add_argument("--stride", type=int, default=8)
    p.add_argument("--fps", type=int, default=15)
    a = p.parse_args()
    with h5py.File(a.trajectory, "r", locking=False) as h:
        ac, ad = h["actin/center"][:], h["actin/direction"][:]
        mc, md = h["myosin/center"][:], h["myosin/direction"][:]
        step = h["state/current_step"][:, 0]
        dt = float(h["meta"].attrs["dt"])
        al = float(h["meta"].attrs["actin_length"])
        ml = float(h["meta"].attrs["myosin_length"])
        mr = float(h["meta"].attrs.get("myosin_radius", .025))
        ar = float(h["meta"].attrs["actin_radius"])
        target = float(h["meta"].attrs["am_preferred_centerline"])
        hs = h["state/am_head_state_current"][:]
        aa_n = h["state/aa_pairs_current_count"][:, 0]
        aa_s = h["state/aa_status_current"][:, :, 0]
    time = step*dt
    strong = np.sum(hs == 2, axis=(1, 2))
    s1 = np.array([np.count_nonzero(aa_s[i, :n] == 1)
                   for i, n in enumerate(aa_n)])
    s2 = np.array([np.count_nonzero(aa_s[i, :n] == 2)
                   for i, n in enumerate(aa_n)])
    frames = range(0, len(time), max(1, a.stride))
    a.output.parent.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(14, 9), facecolor="white")
    writer = FFMpegWriter(fps=a.fps, codec="libx264",
                          extra_args=("-pix_fmt", "yuv420p", "-crf", "18",
                                      "-threads", "1"))
    with writer.saving(fig, str(a.output), dpi=120):
        for f in frames:
            fig.clear()
            grid = fig.add_gridspec(3, 3, height_ratios=(1, 1, .55),
                                    width_ratios=(1.45, 1, 1),
                                    hspace=.36, wspace=.28)
            ax3 = fig.add_subplot(grid[:2, 0], projection="3d")
            ax3.set_proj_type("ortho")
            for m in range(2):
                rod(ax3, mc[f, m], md[f, m], ml, "#263238", 8)
            for ring in range(4):
                for i in range(6*ring, 6*ring+6):
                    rod(ax3, ac[f, i], ad[f, i], al, COLORS[ring], 3)
            ax3.set(xlim=(-2.2, 2.2), ylim=(-.13, .13), zlim=(-.13, .13),
                    xlabel="x (µm)", ylabel="y (µm)", zlabel="z (µm)")
            ax3.set_box_aspect((3.2, 1, 1))
            ax3.view_init(elev=22, azim=-70)
            ax3.set_title("2 myosins × 2 halves × 6 actins = 24 actins")
            for ring, owner in enumerate(OWNERS):
                ax = fig.add_subplot(grid[ring//2, 1+ring%2])
                axis = md[f, owner]
                ref = np.array([0., 0., 1.])
                if abs(ref@axis) > .9: ref = np.array([0., 1., 0.])
                e1 = np.cross(axis, ref); e1 /= np.linalg.norm(e1)
                e2 = np.cross(axis, e1)
                d = ac[f, 6*ring:6*ring+6]-mc[f, owner]
                d -= (d@axis)[:, None]*axis
                yz = 1e3*np.c_[d@e1, d@e2]
                ax.add_patch(Circle((0, 0), 1e3*mr, fc="#455A64",
                                    ec="none", alpha=.75))
                ax.add_patch(Circle((0, 0), 1e3*target, fill=False,
                                    ec="black", ls="--", lw=1.2))
                for point in yz:
                    ax.add_patch(Circle(point, 1e3*ar, fc=COLORS[ring],
                                        ec="white", lw=.7, alpha=.9))
                closed = np.vstack((yz, yz[0]))
                ax.plot(closed[:, 0], closed[:, 1], color=COLORS[ring],
                        lw=1.2, alpha=.65)
                ax.set(xlim=(-100, 100), ylim=(-100, 100),
                       aspect="equal", title=LABELS[ring],
                       xlabel="transverse 1 (nm)",
                       ylabel="transverse 2 (nm)")
                ax.grid(alpha=.2)
            counts = fig.add_subplot(grid[2, :])
            counts.plot(1e3*time, strong, color="#5E3C99",
                        label="strong AM heads")
            counts.plot(1e3*time, s1, color="#E69F00",
                        label="status-1 ACTN2")
            counts.plot(1e3*time, s2, color="#009E73",
                        label="status-2 ACTN2")
            counts.axvline(1e3*time[f], color="black", lw=1)
            counts.set(xlim=(0, 1e3*time[-1]), xlabel="time (ms)",
                       ylabel="count")
            counts.legend(ncol=3, loc="upper right")
            counts.grid(alpha=.2)
            fig.suptitle(f"Force-balanced 24-actin cell — t={1e3*time[f]:.1f} ms",
                         fontsize=15, weight="bold")
            writer.grab_frame(facecolor="white")
    plt.close(fig)
    print(a.output)

if __name__ == "__main__":
    main()
