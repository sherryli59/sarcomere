"""Clean MP4 export of a trajectory, drawn exactly as the interactive viewer draws
it: same orthographic camera (rotate() with azimuth -0.46, elevation 0.30), same
painter's depth order, same species hues and per-layer width/alpha.  Only the
legend is drawn -- no title, time stamp, scale bar or other text.

Input: a list (one per frame) of dicts mapping layer key -> float array of
segments (n, 6) or points (n, 3) in um."""
import os, subprocess, tempfile, shutil, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from matplotlib.colors import to_rgba

# viewer layer style contract (VISUALIZATION_PIPELINE.md / cluster-traj-template.html)
STYLE = {
    "ctx":       dict(color="#1baf7a", width=0.9, alpha=0.22, kind="segments", label=None),
    "actin":     dict(color="#1baf7a", width=1.9, alpha=0.92, kind="segments", label="actin"),
    "myosin":    dict(color="#eda100", width=4.8, alpha=0.86, kind="segments", label="myosin"),
    "barbed":    dict(color="#008300", width=3.2, alpha=0.90, kind="points",   label="barbed end"),
    "armWeak":   dict(color="#7c5cd6", width=1.4, alpha=0.72, kind="segments", label="weak head"),     # viewer 1.2/0.60, bolder for video
    "armStrong": dict(color="#7c5cd6", width=2.6, alpha=1.00, kind="segments", label="strong head"),
    "aa1":       dict(color="#e87ba4", width=1.2, alpha=0.50, kind="segments", label="ACTN2 status-1"),  # viewer 0.9/0.30, bolder for video
    "aa2":       dict(color="#e87ba4", width=1.5, alpha=0.85, kind="segments", label="ACTN2 status-2"),
}
LEGEND = ["actin", "myosin", "barbed", "armWeak", "armStrong", "aa1", "aa2"]

def rotate(p, az=-0.46, el=0.30):
    ca, sa, ce, se = np.cos(az), np.sin(az), np.cos(el), np.sin(el)
    x1 = ca * p[..., 0] - sa * p[..., 1]; y1 = sa * p[..., 0] + ca * p[..., 1]
    return np.stack([x1, ce * y1 - se * p[..., 2], se * y1 + ce * p[..., 2]], -1)

def render(frames, out, fps=12, size=(1920, 1080), az=-0.46, el=0.30, px_scale=1.8, dpi=100, margin=0.05,
           fit_keys=("actin", "myosin", "barbed", "armWeak", "armStrong", "aa1", "aa2"), legend_frac=0.17):
    W, H = size
    # fixed camera over the whole trajectory (as the viewer), fit to the drawn geometry
    pts = []
    for fr in frames:
        for k, a in fr.items():
            if len(a) == 0 or k not in fit_keys: continue
            pts.append(a.reshape(-1, 3))
    P = rotate(np.concatenate(pts), az, el)
    xmin, xmax = P[:, 0].min(), P[:, 0].max(); zmin, zmax = P[:, 2].min(), P[:, 2].max()
    Wd = W * (1 - legend_frac)                       # drawing area; the legend owns the right strip
    scale = min(Wd * (1 - 2 * margin) / max(xmax - xmin, 0.1), H * (1 - 2 * margin) / max(zmax - zmin, 0.1))
    cx = Wd / 2 - scale * (xmin + xmax) / 2; cy = H / 2 + scale * (zmin + zmax) / 2
    lw = lambda px: px * px_scale * 72.0 / dpi
    tmp = tempfile.mkdtemp(prefix="mp4_")
    try:
        for i, fr in enumerate(frames):
            fig = plt.figure(figsize=(W / dpi, H / dpi), dpi=dpi); ax = fig.add_axes([0, 0, 1, 1])
            ax.set_xlim(0, W); ax.set_ylim(H, 0); ax.axis("off"); fig.patch.set_facecolor("white")
            segs, cols, widths, depth = [], [], [], []
            for k, st in STYLE.items():
                a = fr.get(k)
                if a is None or len(a) == 0: continue
                rgba = to_rgba(st["color"], st["alpha"])
                if st["kind"] == "segments":
                    s = a.reshape(-1, 2, 3); r = rotate(s, az, el)
                    xy = np.stack([cx + scale * r[..., 0], cy - scale * r[..., 2]], -1)
                    d = 0.5 * (r[:, 0, 1] + r[:, 1, 1])
                    segs += list(xy); depth += list(d); cols += [rgba] * len(xy); widths += [lw(st["width"])] * len(xy)
                else:   # points: zero-length round-capped segment, diameter = 2 * viewer radius
                    r = rotate(a.reshape(-1, 3), az, el)
                    xy = np.stack([cx + scale * r[:, 0], cy - scale * r[:, 2]], -1)
                    segs += [np.stack([p, p + 1e-3]) for p in xy]; depth += list(r[:, 1])
                    cols += [rgba] * len(xy); widths += [lw(2 * st["width"])] * len(xy)
            order = np.argsort(depth, kind="stable")
            lc = LineCollection([segs[j] for j in order], colors=[cols[j] for j in order],
                                linewidths=[widths[j] for j in order], capstyle="round")
            ax.add_collection(lc)
            from matplotlib.patches import Rectangle
            clip = Rectangle((0, 0), Wd, H, transform=ax.transData); lc.set_clip_path(clip)
            handles = []
            for k in LEGEND:
                st = STYLE[k]
                if st["kind"] == "points":
                    handles.append(Line2D([], [], ls="", marker="o", ms=lw(2 * st["width"]), color=to_rgba(st["color"], st["alpha"]), label=st["label"]))
                else:
                    handles.append(Line2D([], [], color=to_rgba(st["color"], st["alpha"]), lw=lw(st["width"]), solid_capstyle="round", label=st["label"]))
            ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1 - legend_frac + 0.005, 0.5), frameon=False,
                      fontsize=15, handlelength=2.4, borderpad=0.8, labelspacing=0.6)
            fig.savefig(os.path.join(tmp, "f%05d.png" % i), dpi=dpi, facecolor="white"); plt.close(fig)
        subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", str(fps), "-i", os.path.join(tmp, "f%05d.png"),
                        "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "18", "-movflags", "+faststart", out], check=True)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    return out
