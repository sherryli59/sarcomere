#!/usr/bin/env python3
"""Strip the text, keep the legend; soften actin, thicken myosin.

The text nodes are hidden rather than deleted: updateText() writes into them
every frame, so removing them from the DOM would throw on a null reference and
blank the viewer. Hiding costs nothing and keeps the code path intact.
"""
import pathlib, sys

D = pathlib.Path("/home/groups/rotskoff/sherry/sarcomere/sarcomere/analysis/"
                 "visualizations/intermediate-5s")

t = D / "intermediate-5s-animation-template.html"
s = t.read_text()
if "hide every text block" not in s:
    anchor = "    #a2400-comparison-snapshots .canvas-wrap { position: relative;"
    css = """    /* Figure mode: hide every text block and leave only the layer legend
       and the transport controls. The nodes stay in the DOM because
       updateText() writes to them on every frame. */
    #a2400-comparison-snapshots h3,
    #a2400-comparison-snapshots .snapshot-meta,
    #a2400-comparison-snapshots .snapshot-detail,
    #a2400-comparison-snapshots #a2400-time,
    #a2400-comparison-snapshots #a2400-run-controls,
    #a2400-comparison-snapshots .comparison-header { display: none; }
"""
    if anchor not in s:
        sys.exit("css anchor missing")
    s = s.replace(anchor, css + anchor, 1)
    t.write_text(s)
    print("template: text hidden, legend kept")
else:
    print("template: already in figure mode")

b = D / "build_full.py"
s = b.read_text()
old = """LAYERS = [("actin0", "actin (unlinked)", "#b8bfbc", 0.9, 0.35, "false"),
          ("actin1", "actin (status-1)", "#6baed6", 1.2, 0.55, "false"),
          ("actin2", "actin (status-2)", "#08306b", 1.7, 0.95, "true"),
          ("myosin", "myosin", "#f5a45b", 4.8, 0.90, "true"),
          ("aa1", "ACTN2 (status-1)", "#f2b8cd", 0.8, 0.40, "false"),
          ("aa2", "ACTN2 (mature)", "#e87ba4", 1.35, 0.85, "true")]"""
new = """LAYERS = [("actin0", "actin (unlinked)", "#b8bfbc", 0.9, 0.30, "false"),
          ("actin1", "actin (status-1)", "#6baed6", 1.2, 0.50, "false"),
          ("actin2", "actin (status-2)", "#08306b", 1.7, 0.80, "true"),
          ("myosin", "myosin", "#f5a45b", 6.6, 0.92, "true"),
          ("aa1", "ACTN2 (status-1)", "#f2b8cd", 0.8, 0.40, "false"),
          ("aa2", "ACTN2 (mature)", "#e87ba4", 1.35, 0.85, "true")]"""
if old not in s:
    sys.exit("LAYERS table not found")
b.write_text(s.replace(old, new, 1))
print("builder: actin alpha -> 0.80, myosin width 4.8 -> 6.6")
