#!/usr/bin/env python3
"""One actor per species instead of one per rod.

plot_filaments_3d called plotter.add_mesh() inside the loop, so a frame with
ACTN2 enabled created ~8100 VTK actors (1120 actin + 320 myosin + 6651
bridges, before periodic duplicates). VTK does not degrade linearly in actor
count, and that turned a ~1 minute render into one that had not produced a
frame after 79 minutes.

Merging the cylinders into a single PolyData per call collapses this to one
actor. Per-filament colour survives as a point-data RGB array, so the Blues
load ramp still works; cylinder resolution is now a parameter so thin bridges
need not pay for 24 facets.
"""
import pathlib, sys
p = pathlib.Path("/home/groups/rotskoff/sherry/sarcomere/sarcomere/analysis/"
                 "visualize_traj.py")
s = p.read_text()
if "merged geometry" in s:
    sys.exit("already batched")

old_sig = ("def plot_filaments_3d(plotter, center, direction, radius, l, Lx, Ly, Lz, "
           "color='k', color_spectrum=None, opacity=1.0):")
new_sig = ("def plot_filaments_3d(plotter, center, direction, radius, l, Lx, Ly, Lz, "
           "color='k', color_spectrum=None, opacity=1.0, resolution=24):")
if old_sig not in s:
    sys.exit("signature not found")
s = s.replace(old_sig, new_sig, 1)

start = s.index("    for i in range(center.shape[0]):\n        if l[i] < 0.01:")
end = s.index("def add_vector_arrows(")
body = '''    # Accumulate merged geometry and add it as ONE actor. Adding a mesh per
    # rod is what made ACTN2-enabled frames unrenderable.
    pieces = []

    def _emit(cyl_center, this_color):
        cyl = pv.Cylinder(center=cyl_center, direction=direction[i],
                          height=l[i], radius=radius[i], resolution=resolution)
        rgb = np.array(mcolors.to_rgb(this_color)) * 255.0
        cyl.point_data["RGB"] = np.tile(rgb, (cyl.n_points, 1)).astype(np.uint8)
        pieces.append(cyl)

    for i in range(center.shape[0]):
        if l[i] < 0.01:
            continue
        this_color = color
        if color_values is not None:
            this_color = plt.cm.Blues(norm(color_values[i]))
        _emit(center[i], this_color)

        dx = 0.5 * l[i] * np.array(direction[i])
        ox_set = get_offsets_for_filament(center[i][0], abs(dx[0]), Lx)
        oy_set = get_offsets_for_filament(center[i][1], abs(dx[1]), Ly)
        oz_set = get_offsets_for_filament(center[i][2], abs(dx[2]), Lz)
        for ox, oy, oz in product(ox_set, oy_set, oz_set):
            if ox == oy == oz == 0:
                continue
            _emit(center[i] + np.array([ox, oy, oz]), this_color)

    if not pieces:
        return
    merged = pieces[0] if len(pieces) == 1 else pv.merge(pieces)
    plotter.add_mesh(merged, scalars="RGB", rgb=True, opacity=opacity)


'''
s = s[:start] + body + s[end:]

# thin bridges do not need 24 facets
s = s.replace('''                        plotter=plotter,
                        color="#e87ba4",
                    )''',
              '''                        plotter=plotter,
                        color="#e87ba4",
                        resolution=8,
                    )''')
p.write_text(s)
print("plot_filaments_3d now emits one merged actor per call")
