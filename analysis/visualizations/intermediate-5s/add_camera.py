#!/usr/bin/env python3
"""Expose the camera angle, and frame the aspect to the projected box.

rotate() maps world -> (screen_x, depth, screen_y). With azimuth 0 and
elevation pi/2 that reduces to screen = (x, y) and depth = -z, i.e. the xy
plane seen head-on with z as the view axis.

The projected box is then 6 x 2 um, a 3:1 rectangle, so a 16:9 frame would sit
it in a band with ~40% of the height empty. --auto-aspect sets the output
height from the projected extent instead, which spends the pixels on the
sample rather than on margins.
"""
import pathlib, sys
p = pathlib.Path("/home/groups/rotskoff/sherry/sarcomere/sarcomere/analysis/"
                 "visualizations/intermediate-5s/html_to_mp4.py")
s = p.read_text()
if "--azimuth" in s:
    sys.exit("already patched")

s = s.replace('    ap.add_argument("--only-frame", type=int, default=-1)',
              '    ap.add_argument("--azimuth", type=float, default=AZ)\n'
              '    ap.add_argument("--elevation", type=float, default=EL)\n'
              '    ap.add_argument("--auto-aspect", action="store_true",\n'
              '                    help="derive height from the projected box extent")\n'
              '    ap.add_argument("--only-frame", type=int, default=-1)')

# every rotate() call must use the requested angles
s = s.replace("    rc = rotate(corners)", "    rc = rotate(corners, a.azimuth, a.elevation)")
s = s.replace("        p0, p1 = rotate(seg[:, 0]), rotate(seg[:, 1])",
              "        p0 = rotate(seg[:, 0], a.azimuth, a.elevation)\n"
              "            p1 = rotate(seg[:, 1], a.azimuth, a.elevation)")

old_frame = """    pad = 0.03*(rc[:, 0].max() - rc[:, 0].min())
    xlim = (rc[:, 0].min()-pad, rc[:, 0].max()+pad)
    ylim = (rc[:, 2].min()-pad, rc[:, 2].max()+pad)"""
new_frame = """    pad = 0.03*(rc[:, 0].max() - rc[:, 0].min())
    xlim = (rc[:, 0].min()-pad, rc[:, 0].max()+pad)
    ylim = (rc[:, 2].min()-pad, rc[:, 2].max()+pad)
    if a.auto_aspect:
        ratio = (ylim[1]-ylim[0]) / (xlim[1]-xlim[0])
        a.height = int(round(a.width * ratio / 2) * 2)
        print(f"auto aspect: projected box {xlim[1]-xlim[0]:.2f} x "
              f"{ylim[1]-ylim[0]:.2f} um -> {a.width}x{a.height}")"""
if old_frame not in s:
    sys.exit("framing anchor not found")
s = s.replace(old_frame, new_frame, 1)
p.write_text(s)
print("camera + auto-aspect options added")
