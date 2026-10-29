"""Pack axial-vs-weighted rig trajectories into compact int16 base64 for the viewer."""
import base64, json, sys
import numpy as np

z = np.load(sys.argv[1]); out = {}
def b64(a, dtype): return base64.b64encode(np.ascontiguousarray(a, dtype='<' + dtype).tobytes()).decode()
theta = np.arange(4) * np.pi / 2
radial = np.stack([np.zeros(4), np.cos(theta), np.sin(theta)], -1)
home0, home1 = .080 * radial + [.750, 0, 0], .120 * radial + [.750, 0, 0]
for case, removed in (('failure', ()), ('remove0', ()), ('remove1', (0,)), ('remove2', (0, 1)), ('remove3', (0, 1, 2))):
    for law in ('weighted', 'axial'):
        k = f'{law}_{case}__'
        c, p0, p1 = z[k + 'center'], z[k + 'p0'].copy(), z[k + 'p1'].copy()
        alive = z[k + 'alive'].copy(); force = z[k + 'fi'] + z[k + 'fj']
        for j in removed:  # parked far away in the run; draw as ghosts at home
            alive[:, :, j] = False; p0[:, :, j] = home0[j]; p1[:, :, j] = home1[j]
        T, R = c.shape[:2]
        pos = np.concatenate([c[:, :, None, :], p0, p1], 2)  # (T,R,9,3) um
        assert np.abs(pos).max() * 1e4 < 32000 and np.abs(force).max() * 10 < 32000
        out[f'{law}_{case}'] = dict(T=T, R=R,
            pos=b64(np.round(pos * 1e4), 'i2'),       # 0.1 nm units
            force=b64(np.round(force * 10), 'i2'),    # 0.1 pN units
            alive=b64(alive, 'u1'))
out['time_ms'] = (z['axial_failure__time'] * 1000).round(4).tolist()
json.dump(out, open(sys.argv[2], 'w'), separators=(',', ':'))
print({k: (v['T'], v['R']) for k, v in out.items() if isinstance(v, dict)})
