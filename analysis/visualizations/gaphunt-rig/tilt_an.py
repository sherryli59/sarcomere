"""Tilt scan: do the actin and the free myosin align with each other?"""
import h5py, numpy as np, math, json, os, statistics, sys
PFX=sys.argv[1] if len(sys.argv)>1 else "ts"
OUTJ=sys.argv[2] if len(sys.argv)>2 else "tilt_summary.json"
os.environ.setdefault("HDF5_USE_FILE_LOCKING","FALSE")
CAMP="/scratch/users/xsli/sarcomere/gaphunt_20261023"
out=[]
for T in (0,15,30,45,60,75,85):
    for S in (20260805,20260815,20260825):
        p="%s/data/%s%d_seed%d.h5"%(CAMP,PFX,T,S)
        if not os.path.exists(p): continue
        f=h5py.File(p,"r",locking=False)
        n=int(min(f["meta/committed_frames"][0],f["actin/direction"].shape[0]))
        d=f["actin/direction"][:n,0]; md=f["myosin/direction"][:n,0]
        c=f["actin/center"][:n,0]; mc=f["myosin/center"][:n,0]
        u=d/np.linalg.norm(d,axis=1,keepdims=True); m=md/np.linalg.norm(md,axis=1,keepdims=True)
        am=np.degrees(np.arccos(np.clip(np.abs(np.sum(u*m,axis=1)),0,1)))
        a0=np.degrees(np.arccos(np.clip(np.abs(u@u[0]),0,1)))      # actin's own rotation
        m0=np.degrees(np.arccos(np.clip(np.abs(m@m[0]),0,1)))      # myosin's own rotation
        hs=f["state/am_head_state_current"][:n]
        strong=(hs==2).sum(axis=(1,2)); weak=(hs==1).sum(axis=(1,2))
        box=np.array([f["meta"].attrs["box_x"],f["meta"].attrs["box_y"],f["meta"].attrs["box_z"]],float)
        r=c-mc; r=r-box*np.round(r/box)          # minimum image: the box is periodic
        perp=np.linalg.norm(r-np.sum(r*m,axis=1)[:,None]*m,axis=1)*1000
        ht=f["actin_myo/head_transitions"][:]
        rupt=int(((ht[:,3]==2)&(ht[:,5]==4)).sum()) if len(ht) else 0
        t0=int(n*0.6); dt=5000.0/1e6*50  # ms per frame: save_every 50 x dt 1e-5
        frame_ms=float(f["state/current_step"][1,0]-f["state/current_step"][0,0])*1e-5*1e3 if n>1 else 0
        half=next((k*frame_ms for k in range(n) if am[k]<=am[0]/2),None) if am[0]>1 else None
        lost=next((k for k in range(n) if strong[k]==0 and weak[k]==0),None)
        sc=np.sum(u*m,axis=1); axc=np.sum(r*m,axis=1); pointed=axc+0.5*sc*1.0
        bare=0.5*0.10*1.5
        midzone=next((k*frame_ms for k in range(n) if abs(pointed[k])<bare),None)
        held=bool(((strong+weak)[t0:]>0).mean()>0.9)
        out.append(dict(midzone_ms=midzone,held_to_end=held,heads_end=float((strong+weak)[t0:].mean()),pointed0=float(pointed[0]),t_lost_ms=(lost*frame_ms if lost is not None else None),
             am_at_lost=(float(am[lost]) if lost is not None else None),am_min=float(am.min()),tilt=T,seed=S,frames=n,start=float(am[0]),end=float(am[t0:].mean()),
             half_ms=half,actin_turn=float(a0[-1]),myosin_turn=float(m0[-1]),
             strong=float(strong[t0:].mean()),strong0=int(strong[0]),weak=float(weak[t0:].mean()),
             perp_end=float(perp[-1]),ruptures=rupt))
json.dump(out,open(OUTJ,"w"),indent=1)
hdr=("tilt","seed","start","t_half","min A-M","heads lost at","A-M then","strong@0","ruptures","perp end")
print("%5s %5s | %6s %7s %8s | %13s %9s | %8s %8s %9s"%hdr)
for r in out:
    th="%.0f ms"%r["half_ms"] if r["half_ms"] is not None else "-"
    tl="%.0f ms"%r["t_lost_ms"] if r["t_lost_ms"] is not None else "never"
    at="%.1f"%r["am_at_lost"] if r["am_at_lost"] is not None else "-"
    print("%5d %5s | %6.1f %7s %8.1f | %13s %9s | %8d %8d %7.0fnm"
          %(r["tilt"],str(r["seed"])[-3:],r["start"],th,r["am_min"],tl,at,
            r["strong0"],r["ruptures"],r["perp_end"]))
