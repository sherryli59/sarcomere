#!/usr/bin/env python3
"""Scale resolution with dpi, not figure size.

figsize = width/100 at dpi=100 grows the canvas in INCHES, so line widths --
which matplotlib measures in points -- cover proportionally fewer pixels and
the image looks soft however many pixels it has. Holding the figure at a fixed
physical size and raising dpi keeps every stroke the same relative weight while
adding real resolution.
"""
import pathlib, sys
p = pathlib.Path("/home/groups/rotskoff/sherry/sarcomere/sarcomere/analysis/"
                 "visualizations/intermediate-5s/html_to_mp4.py")
s = p.read_text()
if "--dpi" in s:
    sys.exit("already patched")
s = s.replace('    ap.add_argument("--width", type=int, default=1600)\n'
              '    ap.add_argument("--height", type=int, default=900)',
              '    ap.add_argument("--width", type=int, default=2560)\n'
              '    ap.add_argument("--height", type=int, default=1440)\n'
              '    ap.add_argument("--dpi", type=int, default=160,\n'
              '                    help="figure stays 16x9 in; dpi sets the pixels")')
s = s.replace("        fig = plt.figure(figsize=(a.width/100, a.height/100), dpi=100)",
              "        fig = plt.figure(figsize=(a.width/a.dpi, a.height/a.dpi),\n"
              "                         dpi=a.dpi)")
s = s.replace('"-crf", "18",', '"-crf", "16",')
p.write_text(s)
print("resolution now driven by dpi; default 2560x1440, crf 16")
