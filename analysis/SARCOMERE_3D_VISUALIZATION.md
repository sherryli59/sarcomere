# Sarcomere 3-D visualization workflow

`animate_visualize_traj3d.py` is the primary clean structural renderer. It is
faithful to the supplied `visualize_traj3d` implementation: PyVista cylinders,
orange status-2 actins, light-blue myosins, dark-blue myosins attached to
status-2 actins, a white background, periodic copies and the fixed camera. The
adaptation only streams current HDF5 frames and assembles PNG/MP4 output.

```bash
python analysis/animate_visualize_traj3d.py run.h5 clean.mp4 \
  --box 12 2 2 --stride 2 --fps 15
```

`animate_sarcomere_3d.py` is the primary mechanical/KMC debugging renderer.
It streams one saved frame at a time and batches each object class into a small
number of VTK meshes, so it is suitable for 2000-actin/1000-myosin runs.

The default `--style clean` follows the visual hierarchy of
`visualize_traj3d.py`: thin filaments, a light background and only the loaded
network. Individual ACTN2 and AM-head connectors are hidden unless
`--show-links` or `--show-heads` is requested. `--style debug` restores the
dense all-filament view. `--selection legacy` reproduces the older
`cb_status*f_load >= 0.01` actin filter and bonded-myosin filter.

Use `--selection engaged` for the larger prebundled network (status 1 or any
bound head), and `--selection loaded` for the mechanically selected network
(status 2 or a strong head).

The default status palette is:

- blue: unlinked actin;
- cyan: status-1 actin and local ACTN2 links;
- magenta: status-2 actin and loaded ACTN2 links;
- gold: nonpolar thick filaments;
- yellow: weak AM heads;
- orange-red: strong AM heads.

Each frame reports actin and myosin nematic order and the rendered AA/AM state
counts.  `--show-polarity` marks each actin barbed end.  `--actin-color
polarity` instead colors actins by the sign of their global x polarity.

Examples:

```bash
export HDF5_USE_FILE_LOCKING=FALSE

python analysis/animate_sarcomere_3d.py run.h5 last_frame.png \
  --box 12 2 2

python analysis/animate_sarcomere_3d.py run.h5 debug.mp4 \
  --box 12 2 2 --stride 2 --fps 15 --style debug \
  --show-links --show-heads --show-forces both

python analysis/animate_sarcomere_3d.py run.h5 transverse.mp4 \
  --box 12 2 2 --camera transverse --parallel-projection --stride 5
```

For routine debugging, make three synchronized movies: isometric/status,
transverse/status, and isometric/polarity.  Force arrows should be enabled for
short diagnostic windows rather than entire production runs.  OVITO remains a
reasonable option for interactive exploration, but PyVista directly supports
the ragged local-link records and explicit head anchors without converting the
trajectory to another format.

For dense full-scale runs, use `--crop XMIN XMAX YMIN YMAX ZMIN ZMAX` to make a
local mechanics movie in addition to the whole-box structural view. Rendering
requires a GPU/EGL node, an X display, or OSMesa. The current login node has
none of these, while the existing `gpu_testjob` allocation renders correctly.

`render_sarcomere_debug_views.sbatch` generates the three recommended views on
a GPU node:

```bash
sbatch analysis/render_sarcomere_debug_views.sbatch run.h5 movies 12 2 2 2 15
```
