#!/usr/bin/env python3
"""Render the interactive viewer's own payload to an MP4.

The frames come from the HTML file itself -- its embedded int16 payload and its
layer table are parsed and replayed -- so the video cannot drift from the
viewer. Camera, colours, widths, alphas, depth ordering and default layer
visibility are all taken from the file rather than restated here.
"""
from __future__ import annotations
import argparse, base64, json, os, re, subprocess
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mplcfg")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

AZ, EL, ZOOM = -0.46, 0.38, 0.93          # the template's defaults
EDGES = [(0,1),(0,2),(0,4),(1,3),(1,5),(2,3),(2,6),(3,7),(4,5),(4,6),(5,7),(6,7)]


def rotate(p, az=AZ, el=EL):
    ca, sa, ce, se = np.cos(az), np.sin(az), np.cos(el), np.sin(el)
    x, y, z = p[..., 0], p[..., 1], p[..., 2]
    x1 = ca*x - sa*y
    y1 = sa*x + ca*y
    return np.stack([x1, ce*y1 - se*z, se*y1 + ce*z], -1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--html", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--workdir", default="/scratch/users/xsli/sarcomere/html_mp4")
    ap.add_argument("--fps", type=int, default=12)
    ap.add_argument("--width", type=int, default=2560)
    ap.add_argument("--height", type=int, default=1440)
    ap.add_argument("--dpi", type=int, default=160,
                    help="figure stays 16x9 in; dpi sets the pixels")
    ap.add_argument("--lw", type=float, default=0.62, help="viewer px -> pt scale")
    ap.add_argument("--azimuth", type=float, default=AZ)
    ap.add_argument("--elevation", type=float, default=EL)
    ap.add_argument("--auto-aspect", action="store_true",
                    help="derive height from the projected box extent")
    ap.add_argument("--only-frame", type=int, default=-1)
    ap.add_argument("--merge-species", action="store_true",
                    help="one style and one legend entry per species")
    ap.add_argument("--all-layers", action="store_true",
                    help="draw hidden layers too; default follows the viewer")
    a = ap.parse_args()

    html = Path(a.html).read_text()
    run = json.loads(re.search(r"const runs = (.+?);\n", html, re.S).group(1))[0]
    layers = [{"key": k, "label": lab, "color": c, "width": float(w),
               "alpha": float(al), "on": on == "true"}
              for k, lab, c, w, al, on in re.findall(
                  r"\{key:'([a-z0-9]+)',label:'([^']*)',color:'(#[0-9a-fA-F]{6})',"
                  r"width:([0-9.]+),alpha:([0-9.]+),kind:'segments',on:(true|false)\}",
                  html)]
    if not layers:
        raise SystemExit("no layer table found in the HTML")
    draw = [l for l in layers if a.all_layers or l["on"]]
    if a.merge_species:
        SPECIES = [("actin",  ["actin0", "actin1", "actin2"], "#08306b", 1.7, 0.80),
                   ("myosin", ["myosin"],                     "#f5a45b", 6.6, 0.92),
                   ("ACTN2",  ["aa1", "aa2"],                 "#e87ba4", 1.35, 0.85)]
        present = {l["key"] for l in draw}
        draw = [{"key": keys, "label": label, "color": c, "width": w, "alpha": al,
                 "on": True}
                for label, keys, c, w, al in SPECIES
                if present.intersection(keys)]
        for d in draw:
            d["keys"] = d["key"] if isinstance(d["key"], list) else [d["key"]]
    else:
        for d in draw:
            d["keys"] = [d["key"]]
    box = np.asarray(run["box"], float)
    half = box / 2.0
    frames = run["frames"]
    print(f"{len(frames)} frames; drawing {[l['key'] for l in draw]}")

    corners = np.array([[sx*half[0], sy*half[1], sz*half[2]]
                        for sx in (-1, 1) for sy in (-1, 1) for sz in (-1, 1)])
    rc = rotate(corners, a.azimuth, a.elevation)
    pad = 0.03*(rc[:, 0].max() - rc[:, 0].min())
    xlim = (rc[:, 0].min()-pad, rc[:, 0].max()+pad)
    ylim = (rc[:, 2].min()-pad, rc[:, 2].max()+pad)
    if a.auto_aspect:
        ratio = (ylim[1]-ylim[0]) / (xlim[1]-xlim[0])
        a.height = int(round(a.width * ratio / 0.90 / 2) * 2)
        print(f"auto aspect: projected box {xlim[1]-xlim[0]:.2f} x "
              f"{ylim[1]-ylim[0]:.2f} um -> {a.width}x{a.height}")
    span = max(xlim[1]-xlim[0], (ylim[1]-ylim[0])*a.width/a.height)
    xlim = (0.5*(xlim[0]+xlim[1]) - span/2/ZOOM, 0.5*(xlim[0]+xlim[1]) + span/2/ZOOM)
    yspan = span*a.height/a.width
    ylim = (0.5*(ylim[0]+ylim[1]) - yspan/2/ZOOM, 0.5*(ylim[0]+ylim[1]) + yspan/2/ZOOM)

    work = Path(a.workdir); work.mkdir(parents=True, exist_ok=True)
    todo = [a.only_frame] if a.only_frame >= 0 else range(len(frames))
    for n, i in enumerate(todo):
        fr = frames[i]
        fig = plt.figure(figsize=(a.width/a.dpi, a.height/a.dpi),
                         dpi=a.dpi)
        # legend strip along the bottom; the axes no longer fill the frame
        ax = fig.add_axes([0.005, 0.085, 0.99, 0.90]); ax.set_axis_off()
        ax.set_xlim(*xlim); ax.set_ylim(*ylim); ax.set_aspect("equal")
        fig.patch.set_facecolor("white")
        for u, v in EDGES:
            p = rc[[u, v]]
            ax.plot(p[:, 0], p[:, 2], color="#c9ccc9", lw=0.6, zorder=0)
        for z, layer in enumerate(draw):
            chunks = [np.frombuffer(base64.b64decode(fr["layers"][k]), np.int16)
                      for k in layer["keys"]
                      if fr["layers"].get(k)]
            if not chunks:
                continue
            q = np.concatenate(chunks).reshape(-1, 2, 3)
            seg = q * half[None, None, :] / 32767.0
            p0 = rotate(seg[:, 0], a.azimuth, a.elevation)
            p1 = rotate(seg[:, 1], a.azimuth, a.elevation)
            order = np.argsort(0.5*(p0[:, 1] + p1[:, 1]))
            ax.add_collection(LineCollection(
                np.stack([np.stack([p0[order, 0], p0[order, 2]], -1),
                          np.stack([p1[order, 0], p1[order, 2]], -1)], 1),
                colors=layer["color"], linewidths=layer["width"]*a.lw,
                alpha=layer["alpha"], capstyle="round", zorder=1+z))
        handles = [plt.Line2D([], [], color=l["color"],
                              lw=max(1.6, l["width"]*a.lw), alpha=max(0.6, l["alpha"]),
                              label=l["label"]) for l in draw]
        ax.legend(handles=handles, loc="upper center",
                  bbox_to_anchor=(0.5, -0.012), ncol=len(handles), frameon=False,
                  fontsize=12, labelcolor="#5a5852", handlelength=2.2,
                  columnspacing=2.4, borderaxespad=0.0)
        fig.savefig(work / f"f{n:04d}.png", facecolor="white")
        plt.close(fig)
        if n % 20 == 0 or n == len(list(todo)) - 1:
            print(f"  {n+1}  t={fr['timeMs']:.0f} ms", flush=True)

    if a.only_frame >= 0:
        print(f"single frame -> {work}/f0000.png"); return
    cmd = ["ffmpeg", "-y", "-framerate", str(a.fps), "-i", str(work/"f%04d.png"),
           "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "16",
           "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2", a.out]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode:
        print(r.stderr[-1200:]); raise SystemExit(1)
    print(f"wrote {a.out}  {Path(a.out).stat().st_size/1048576:.1f} MB  "
          f"{len(frames)} frames @ {a.fps} fps")


if __name__ == "__main__":
    main()
