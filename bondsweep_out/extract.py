import numpy as np, glob, os
print('%-14s %-8s %-6s %-6s %s'%('tag','caps','psi6','occ','minAA_nm'))
for t in ['bond_hard','bond_mid','bond_soft','bond_vsoft','all_soft','steric_only']:
    p=t+'.csv/timeseries.csv'
    if not os.path.exists(p): print('%-14s MISSING'%t); continue
    d=np.genfromtxt(p,delimiter=',',names=True)
    L=d[-1]
    print('%-14s %-8d %-6.3f %-6.2f %.1f'%(t,int(L['cumulative_cap_events']),L['mean_psi6'],L['expected_aa_edge_occupancy'],L['minimum_aa_segment_distance_nm']))
