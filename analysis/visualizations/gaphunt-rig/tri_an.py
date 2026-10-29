"""One actin shared by three parallel myosins vs one myosin, same actin pose."""
import h5py, numpy as np, math, json, os
os.environ.setdefault("HDF5_USE_FILE_LOCKING","FALSE")
CAMP="/scratch/users/xsli/sarcomere/gaphunt_20261023"; SEED=20260805
def ang(a,b): return np.degrees(np.arccos(np.clip(np.abs(np.sum(a*b,axis=-1)),0,1)))
out=[]
for name in ("tri","mono"):
    for T in (0,15,30,45):
        f=h5py.File("%s/data/%s%d_seed%d.h5"%(CAMP,name,T,SEED),"r",locking=False)
        n=int(f["meta/committed_frames"][0]); box=np.array([f["meta"].attrs[k] for k in ("box_x","box_y","box_z")],float)
        fm=float(f["state/current_step"][1,0]-f["state/current_step"][0,0])*1e-5*1e3
        u=f["actin/direction"][:n,0]; u=u/np.linalg.norm(u,axis=1,keepdims=True)
        md=f["myosin/direction"][:n]; md=md/np.linalg.norm(md,axis=2,keepdims=True); nM=md.shape[1]
        am=ang(u[:,None,:],md)                       # (frames, nM) actin vs each myosin
        amm=am.mean(axis=1)
        hs=f["state/am_head_state_current"][:n]; per=((hs==1)|(hs==2)).sum(axis=2)   # (frames, nM)
        tot=per.sum(axis=1)
        spread=np.zeros(n)
        if nM>1:
            spread=np.median(np.stack([ang(md[:,i],md[:,j]) for i in range(nM) for j in range(i+1,nM)],1),axis=1)
        mbar=md.mean(axis=1); mbar=mbar/np.linalg.norm(mbar,axis=1,keepdims=True)
        c=f["actin/center"][:n,0]; mc=f["myosin/center"][:n].mean(axis=1)
        r=c-mc; r=r-box*np.round(r/box); ax=np.sum(r*mbar,axis=1); pointed=ax+0.5*np.sum(u*mbar,axis=1)
        last=int(np.nonzero(tot>0)[0].max()) if (tot>0).any() else 0
        half=next((k*fm for k in range(n) if amm[k]<=amm[0]/2),None) if amm[0]>1 else None
        t0=int(n*0.6)
        out.append(dict(rig=name,tilt=T,nM=nM,start=float(amm[0]),min_am=float(amm.min()),half_ms=half,
            am_end=float(amm[t0:].mean()),spread_max=float(spread.max()),spread_end=float(spread[-1]),
            heads_release=int(tot[0]),heads_mid=float(tot[min(n-1,int(100/fm))]),
            per_myosin_mid=[int(x) for x in per[min(n-1,int(100/fm))]],
            attached_until=(None if tot[-1]>0 else last*fm), run_ms=(n-1)*fm,
            pointed0=float(pointed[0]),pointed_100=float(pointed[min(n-1,int(100/fm))])))
json.dump(out,open("tri_summary.json","w"),indent=1)
print("%-5s %4s | %6s %8s %7s | %13s | %8s %19s | %10s %16s"
      %("rig","tilt","start","t_half","min A-M","myosin spread","heads@0","heads@100ms (per M)","attached","pointed x 0->100ms"))
for r in out:
    print("%-5s %4d | %6.1f %8s %7.1f | %5.1f max/%4.1f | %8d %8.0f %-10s | %10s %8.2f -> %5.2f"
          %(r["rig"],r["tilt"],r["start"],("%.0f ms"%r["half_ms"]) if r["half_ms"] is not None else "-",r["min_am"],
            r["spread_max"],r["spread_end"],r["heads_release"],r["heads_mid"],str(r["per_myosin_mid"]),
            "whole run" if r["attached_until"] is None else "%.0f ms"%r["attached_until"],r["pointed0"],r["pointed_100"]))
