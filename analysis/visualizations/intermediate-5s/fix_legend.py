#!/usr/bin/env python3
"""Put the legend in its own strip instead of over the sample.

The axes filled the whole figure, so an in-axes legend necessarily covers
geometry -- and in the head-on view it landed on the filaments. Reserve a band
at the bottom and lay the entries out horizontally there.
"""
import pathlib, sys
p = pathlib.Path("/home/groups/rotskoff/sherry/sarcomere/sarcomere/analysis/"
                 "visualizations/intermediate-5s/html_to_mp4.py")
s = p.read_text()
if "legend strip" in s:
    sys.exit("already patched")

s = s.replace("        ax = fig.add_axes([0, 0, 1, 1]); ax.set_axis_off()",
              "        # legend strip along the bottom; the axes no longer fill the frame\n"
              "        ax = fig.add_axes([0.005, 0.085, 0.99, 0.90]); ax.set_axis_off()")
s = s.replace("""        ax.legend(handles=handles, loc="upper right", frameon=False, fontsize=11,
                  labelcolor="#5a5852", borderaxespad=1.2)""",
              """        ax.legend(handles=handles, loc="upper center",
                  bbox_to_anchor=(0.5, -0.012), ncol=len(handles), frameon=False,
                  fontsize=12, labelcolor="#5a5852", handlelength=2.2,
                  columnspacing=2.4, borderaxespad=0.0)""")
# the reserved strip changes the usable aspect
s = s.replace("        a.height = int(round(a.width * ratio / 2) * 2)",
              "        a.height = int(round(a.width * ratio / 0.90 / 2) * 2)")
p.write_text(s)
print("legend moved to a bottom strip")
