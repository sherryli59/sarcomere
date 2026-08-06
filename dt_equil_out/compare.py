import numpy as np, subprocess, os
PY='/home/groups/rotskoff/sherry/mamba_root/envs/torch_env/bin/python'
base='/home/groups/rotskoff/sherry/sarcomere/sarcomere'
def analyze(tag):
    d=base+'/dt_equil_out/'+tag+'.csv'
    subprocess.run([PY,base+'/analysis/analyze_shared_lattice.py',base+'/dt_equil_out/'+tag+'.h5','--output',d],
                   stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    a=np.genfromtxt(d+'/timeseries.csv',delimiter=',',names=True)
    return a
def stat(a,col):
    v=a[col]; h=v[len(v)//2:]  # second half = equilibrated
    return h.mean(), h.std()
cols=['mean_psi6','expected_aa_edge_occupancy','actin_nematic','myosin_nematic','am_std_um','minimum_aa_segment_distance_nm']
ref=analyze('ref_1e6'); test=analyze('test_1e5')
print('%-34s %-18s %-18s'%('observable (2nd-half mean +/- std)','dt=1e-6 ref','dt=1e-5 test'))
for c in cols:
    rm,rs=stat(ref,c); tm,ts=stat(test,c)
    print('%-34s %8.4f+/-%-7.4f %8.4f+/-%-7.4f'%(c,rm,rs,tm,ts))
# cap totals
print('caps: ref=%d  test=%d'%(int(ref['cumulative_cap_events'][-1]),int(test['cumulative_cap_events'][-1])))
