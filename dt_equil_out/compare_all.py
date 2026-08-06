import numpy as np, subprocess, os
PY='/home/groups/rotskoff/sherry/mamba_root/envs/torch_env/bin/python'
base='/home/groups/rotskoff/sherry/sarcomere/sarcomere'
def analyze(tag):
    d=base+'/dt_equil_out/'+tag+'.csv'
    if not os.path.exists(d+'/timeseries.csv'):
        subprocess.run([PY,base+'/analysis/analyze_shared_lattice.py',base+'/dt_equil_out/'+tag+'.h5','--output',d],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    return np.genfromtxt(d+'/timeseries.csv',delimiter=',',names=True)
def stat(a,col):
    v=a[col]; h=v[len(v)//2:]; return h.mean(), h.std()
runs=[('1e-6','ref_1e6'),('2e-6','test_2e6'),('3e-6','test_3e6'),('5e-6','test_5e6'),('1e-5','test_1e5')]
print('%-6s %-15s %-15s %-15s %-11s %s'%('dt','occupancy','psi6','tether_std','minAA_nm','caps'))
for lbl,tag in runs:
    try: a=analyze(tag)
    except Exception: print('%-6s MISSING'%lbl); continue
    om,osd=stat(a,'expected_aa_edge_occupancy'); pm,ps=stat(a,'mean_psi6')
    tm,ts=stat(a,'am_std_um'); mm,ms=stat(a,'minimum_aa_segment_distance_nm')
    caps=int(a['cumulative_cap_events'][-1])
    print('%-6s %.3f+/-%.3f   %.3f+/-%.3f   %.4f+/-%.4f  %.1f+/-%.1f %d'%(lbl,om,osd,pm,ps,tm,ts,mm,ms,caps))
