import h5py, numpy as np
for t in ['all_soft','bond_hard']:
    h=h5py.File(t+'.h5','r')
    def last(sp,name,n):
        a=np.array(h[sp][name]) if name in h[sp] else None
        return a
    out={}
    for sp,n in [('actin',16),('myosin',4)]:
        g=h[sp]
        for nm in ['cumulative_force_cap_events','cumulative_torque_cap_events','cumulative_rotation_cap_events']:
            if nm in g:
                arr=np.array(g[nm])
                # arr shape (frames, n) or (n,frames); take final frame sum
                tot=arr.reshape(arr.shape[0],-1)[-1].sum() if arr.ndim>1 else arr[-1]
                out[sp+'.'+nm.replace('cumulative_','').replace('_cap_events','')]=int(tot)
    print(t, out)
    h.close()
