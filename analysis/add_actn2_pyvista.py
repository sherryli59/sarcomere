#!/usr/bin/env python3
"""Draw ACTN2 bridges in visualize_traj.py.

The script renders actin and myosin but has no ACTN2 layer at all, so the
PyVista path cannot show the crosslink network. Each live bridge is a short rod
between two attachment points on two actin filaments; expressing it as
centre/direction/length lets it reuse plot_filaments_3d, so it inherits the same
shaded cylinders and periodic-image handling as the filaments.
"""
import pathlib, sys
p = pathlib.Path("/home/groups/rotskoff/sherry/sarcomere/sarcomere/analysis/"
                 "visualize_traj.py")
s = p.read_text()
if "show_actn2" in s:
    sys.exit("already added")

# CLI
s = s.replace('''    parser.add_argument(
        "--actin_color",''',
'''    parser.add_argument("--show_actn2", action="store_true",
                        help="Draw live ACTN2 crosslinks as rods between attachment points.")
    parser.add_argument("--actn2_radius", type=float, default=0.006,
                        help="Cylinder radius for ACTN2 bridges.")
    parser.add_argument("--actn2_status", choices=["2", "all"], default="2",
                        help="Draw only mature (status-2) bridges, or every live one.")
    parser.add_argument(
        "--actin_color",''')

# plot_system signature
s = s.replace('                actin_color="load",\n',
              '                actin_color="load", show_actn2=False,\n'
              '                actn2_radius=0.006, actn2_status="2",\n', 1)
s = s.replace('                           actin_color="load",\n',
              '                           actin_color="load", show_actn2=False,\n'
              '                           actn2_radius=0.006, actn2_status="2",\n', 1)

# draw the bridges just before the myosin block
anchor = """    # ------------------------------------------------------------------
    # Myosin filaments
    # ------------------------------------------------------------------"""
draw = '''    # ------------------------------------------------------------------
    # ACTN2 crosslinks
    # ------------------------------------------------------------------
    if show_actn2:
        off_key = "/catch_bond/local_link_frame_offsets"
        id_key = "/catch_bond/local_link_identity"
        co_key = "/catch_bond/local_link_coordinates"
        if all(k in data for k in (off_key, id_key, co_key)):
            _, offset, count = np.asarray(data[off_key][frame], dtype=int)
            identity = np.asarray(data[id_key][offset:offset + count], dtype=int)
            coords = np.asarray(data[co_key][offset:offset + count], dtype=float)
            if actn2_status == "2":
                sel = identity[:, 4] == 2
                identity, coords = identity[sel], coords[sel]
            if len(identity):
                all_center = np.asarray(data["/actin/center"][frame], dtype=float)
                all_dir = np.asarray(data["/actin/direction"][frame], dtype=float)
                first, second = identity[:, 2], identity[:, 3]
                starts = all_center[first] + coords[:, :1] * all_dir[first]
                raw_ends = all_center[second] + coords[:, 1:] * all_dir[second]
                box_vec = np.array([Lx, Ly, Lz], dtype=float)
                # minimum image before measuring the rod, or a bridge across a
                # periodic face becomes a span of the whole box
                ends = starts + (raw_ends - starts
                                 - box_vec * np.rint((raw_ends - starts) / box_vec))
                spans = ends - starts
                lengths = np.linalg.norm(spans, axis=1)
                keep = lengths > 1e-6
                if np.any(keep):
                    print(f"Frame {frame}: {int(np.sum(keep))} ACTN2 bridges drawn")
                    plot_filaments_3d(
                        center=0.5 * (starts + ends)[keep],
                        direction=spans[keep] / lengths[keep, None],
                        radius=actn2_radius,
                        l=lengths[keep],
                        Lx=Lx, Ly=Ly, Lz=Lz,
                        plotter=plotter,
                        color="#e87ba4",
                    )

''' + anchor
if anchor not in s:
    sys.exit("myosin anchor not found")
s = s.replace(anchor, draw, 1)

# thread from argv
s = s.replace("    actin_color = args.actin_color",
              "    actin_color = args.actin_color\n"
              "    show_actn2 = args.show_actn2\n"
              "    actn2_radius = args.actn2_radius\n"
              "    actn2_status = args.actn2_status", 1)
s = s.replace("                      actin_color=actin_color,",
              "                      actin_color=actin_color,\n"
              "                      show_actn2=show_actn2,\n"
              "                      actn2_radius=actn2_radius,\n"
              "                      actn2_status=actn2_status,", 1)
p.write_text(s)
print("ACTN2 layer added to visualize_traj.py")
