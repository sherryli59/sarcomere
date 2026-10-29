#!/usr/bin/env python3
"""Run analysis/animate_sarcomere_3d.py unchanged, zooming the camera for the small rig.

The renderer frames the camera on the whole display box, which leaves the
0.5 um rig cross-section tiny in the axial view. This wrapper only zooms the
camera before each screenshot; all geometry, colors and styles are the
renderer's own. Arguments are passed through to animate_sarcomere_3d.py.
"""
import runpy
import sys

import pyvista as pv

RENDERER = '/home/groups/rotskoff/sherry/sarcomere/sarcomere/analysis/animate_sarcomere_3d.py'
ZOOM = {'axial': 3.0, 'transverse': 1.0, 'isometric': 1.0, 'oblique': 1.15}
# 'oblique' views from the -y side and above, so every filament shows its full
# length and +x reads left to right. The renderer only knows its own cameras,
# so it is handed 'isometric' and the position is replaced before each shot.
OBLIQUE = [(0.45, -3.6, 1.5), (0., 0., 0.), (0., 0., 1.)]
camera = sys.argv[sys.argv.index('--camera')+1] if '--camera' in sys.argv else 'isometric'
if camera == 'oblique':
    sys.argv[sys.argv.index('--camera')+1] = 'isometric'
screenshot = pv.Plotter.screenshot


def zoomed(self, *args, **kwargs):
    if camera == 'oblique':
        self.camera_position = OBLIQUE
    self.camera.zoom(ZOOM[camera])
    return screenshot(self, *args, **kwargs)


pv.Plotter.screenshot = zoomed
sys.argv = [RENDERER]+sys.argv[1:]
runpy.run_path(RENDERER, run_name='__main__')
