#!/usr/bin/env python3
"""Add an opt-in uniform actin colour.

--actin_display cb already restricts drawing to crosslinked (cb_status == 2)
filaments. What hides most of them is the colour: actin is shaded through
plt.cm.Blues by load, so a crosslinked filament carrying little load renders
near-white and disappears against the background.

This adds --actin_color {load,uniform}; `load` is the existing behaviour and
stays the default, `uniform` puts every drawn filament at the top of the Blues
ramp. Implemented by flattening the metric rather than bypassing the colormap,
so the two paths share one code path.
"""
import pathlib, re, sys
p = pathlib.Path("/home/groups/rotskoff/sherry/sarcomere/sarcomere/analysis/"
                 "visualize_traj.py")
s = p.read_text()
if "--actin_color" in s:
    sys.exit("already added")

# 1. CLI flag
s = s.replace('''    parser.add_argument(
        "--actin_display",''',
'''    parser.add_argument(
        "--actin_color",
        choices=["load", "uniform"],
        default="load",
        help="Shade crosslinked actin by load (default) or draw it all one dark blue.",
    )
    parser.add_argument(
        "--actin_display",''')

# 2. plot_system signature + the metric
s = s.replace('                myosin_radius, myosin_display="all", actin_display="cb",\n',
              '                myosin_radius, myosin_display="all", actin_display="cb",\n'
              '                actin_color="load",\n', 1)
old_metric = '''    load_metric = np.maximum(f_load, 0) * (cb_status == 2)'''
new_metric = '''    if actin_color == "uniform":
        # Flat metric -> every drawn filament lands at the top of the Blues
        # ramp, so crosslinked-but-unloaded actin is visible instead of white.
        load_metric = (cb_status == 2).astype(float)
    else:
        load_metric = np.maximum(f_load, 0) * (cb_status == 2)'''
if old_metric not in s:
    sys.exit("load_metric not found")
s = s.replace(old_metric, new_metric, 1)

# 3. matplotlib fallback signature, so the kwarg does not blow it up
s = s.replace('                           myosin_radius, myosin_display="all", actin_display="cb",\n',
              '                           myosin_radius, myosin_display="all", actin_display="cb",\n'
              '                           actin_color="load",\n', 1)

# 4. thread it from argv into the worker call
s = s.replace("    myosin_display = args.myosin_display",
              "    myosin_display = args.myosin_display\n    actin_color = args.actin_color", 1)
s = s.replace("                      myosin_display=myosin_display,",
              "                      myosin_display=myosin_display,\n"
              "                      actin_color=actin_color,", 1)
p.write_text(s)
print("added --actin_color {load,uniform}; default unchanged")
