#!/usr/bin/env python3
"""Write standalone ACTN2 rig trajectories in the animate_sarcomere_3d.py schema.

The rig simulates only connector endpoints and one myosin. Each endpoint is
drawn on an actin filament following the standalone weak-site geometry: the
80 nm site sits on a +x actin and the 120 nm site on its -x antiparallel
partner, each 100 nm from the barbed end. Those actins translate rigidly with
their endpoints and are not simulated. Coordinates are shifted so the scene is
centred in a nonperiodic display box.
"""
import argparse
from pathlib import Path

import h5py
import numpy as np

ACTIN_LENGTH, MYOSIN_LENGTH, BARBED_OFFSET = 1.0, 1.5, .100
SHIFT = np.array([-.775, 0., 0.])  # puts x in about [-1.6, 1.6]
DT, STRIDE = 1.e-7, 100  # rig integration step and save stride
HOME0 = .080*np.stack([np.zeros(4), np.cos(np.arange(4)*np.pi/2), np.sin(np.arange(4)*np.pi/2)], -1)+[.750, 0, 0]
HOME1 = HOME0+.040*(HOME0-[.750, 0, 0])/.080


def export(z, key, realization, removed, path, law):
    c = z[f'{key}__center'][:, realization]
    p0, p1 = z[f'{key}__p0'][:, realization].copy(), z[f'{key}__p1'][:, realization].copy()
    fi, fj = z[f'{key}__fi'][:, realization], z[f'{key}__fj'][:, realization]
    alive = z[f'{key}__alive'][:, realization].copy()
    for j in removed:  # parked far away during the run; show the bare actins at home
        alive[:, j] = False
        p0[:, j], p1[:, j] = HOME0[j], HOME1[j]
    T = len(c)
    # Actin 2j carries p0 of connector j (direction +x); actin 2j+1 carries p1 (direction -x).
    direction = np.tile(np.array([[1., 0, 0], [-1., 0, 0]]), (4, 1))
    sites = np.stack([p0, p1], 2).reshape(T, 8, 3)
    coordinate = .5*ACTIN_LENGTH-BARBED_OFFSET  # site position from actin centre along direction
    actin_center = sites-coordinate*direction+SHIFT
    actin_force = np.stack([fi, fj], 2).reshape(T, 8, 3)*alive.repeat(2, 1)[..., None]
    status = alive.repeat(2, 1).astype(np.int32)*2
    offsets, identity, coords = [], [], []
    for t in range(T):
        live = np.flatnonzero(alive[t])
        offsets.append([t, len(identity), len(live)])
        for j in live:
            identity.append([j, 0, 2*j, 2*j+1, 2])
            coords.append([coordinate, coordinate])
    with h5py.File(path, 'w') as h:
        h['/actin/center'] = actin_center
        h['/actin/direction'] = np.broadcast_to(direction, (T, 8, 3)).copy()
        h['/actin/cb_status'] = status
        h['/actin/force'] = actin_force
        h['/myosin/center'] = (c+SHIFT)[:, None, :]
        h['/myosin/direction'] = np.tile([[[1., 0, 0]]], (T, 1, 1))
        h['/myosin/force'] = -(fi+fj).sum(1)[:, None, :]
        h['/state/current_step'] = (np.arange(T)*STRIDE)[:, None]
        h['/catch_bond/local_link_frame_offsets'] = np.asarray(offsets, dtype=np.int64)
        h['/catch_bond/local_link_identity'] = np.asarray(identity, dtype=np.int64).reshape(-1, 5)
        h['/catch_bond/local_link_coordinates'] = np.asarray(coords, dtype=float).reshape(-1, 2)
        m = h.create_group('meta').attrs
        m['dt'] = DT
        m['actin_hydrodynamic_length'] = ACTIN_LENGTH
        m['myosin_hydrodynamic_length'] = MYOSIN_LENGTH
        m['contact_law'] = law
        m['source'] = f'standalone_actn2_axial_contact trajectories.npz key {key} realization {realization}'
        m['actins_are_display_only'] = True
    print(path, T, 'frames; live connectors at end:', int(alive[-1].sum()))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('trajectories', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    z = np.load(args.trajectories)
    for law in ('weighted', 'axial'):
        export(z, f'{law}_failure', 1, (), args.output/f'escape_run1_{law}.h5', law)
        export(z, f'{law}_remove1', 0, (0,), args.output/f'c1_removed_{law}.h5', law)


if __name__ == '__main__':
    main()
