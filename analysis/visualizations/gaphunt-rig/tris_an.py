"""Shared actin + three sigma-0.4 myosins: geometry at release, then time course."""
import h5py, numpy as np, math, json, os
os.environ.setdefault("HDF5_USE_FILE_LOCKING","FALSE")
CAMP="/scratch/users/xsli/sarcomere/gaphunt_20261023"
def ang(a,b): return np.degrees(np.arccos(np.clip(np.abs(np.sum(a*b,axis=-1)),0,1)))
def segdist(p0,p1,q0,q1):
    d1=p1-p0; d2=q1-q0; r=p0-q0; a=d1@d1; e=d2@d2; f=d2@r; c=d1@r; b=d1@d2; den=a*e-b*b
    s=np.clip((b*f-c*e)/den,0,1) if den>1e-14 else 0.0
    t=np.clip((b*s+f)/e,0,1); s=np.clip((b*t-c)/a,0,1)
    return float(np.linalg.norm((p0+s*d1)-(q0+t*d2)))
drawn=json.load(open(CAMP+"/tris_draw.json"))
TT=(0,25,50,100,200,300,400,500,700,1000); out={"drawn":drawn,"t":TT,"runs":{}}
for T in (0,15,30,45):
    f=h5py.File("%s/data/tris%d_seed20260805.h5"%(CAMP,T),"r",locking=False)
    n=int(f["meta/committed_frames"][0]); box=np.array([f["meta"].attrs[k] for k in ("box_x","box_y","box_z")],float)
    fm=float(f["state/current_step"][1,0]-f["state/current_step"][0,0])*1e-5*1e3
    u=f["actin/direction"][:n,0]; u=u/np.linalg.norm(u,axis=1,keepdims=True)
    md=f["myosin/direction"][:n]; md=md/np.linalg.norm(md,axis=2,keepdims=True)
    c=f["actin/center"][:n,0]; mc=f["myosin/center"][:n]
    am=ang(u[:,None,:],md); sp=np.median(np.stack([ang(md[:,i],md[:,j]) for i in range(3) for j in range(i+1,3)],1),axis=1)
    hs=f["state/am_head_state_current"][:n]; per=((hs==1)|(hs==2)).sum(axis=2); tot=per.sum(1)
    # geometry at release (frame 0 = after the sterics pre-phase and the 100 ms bind)
    perp=[]; clr=[]
    for k in range(3):
        r=c[0]-mc[0,k]; r=r-box*np.round(r/box); perp.append(float(np.linalg.norm(r-(r@md[0,k])*md[0,k])))
        q0=mc[0,k]-0.75*md[0,k]; q1=mc[0,k]+0.75*md[0,k]; mcen=mc[0,k]
        clr.append(segdist(c[0]-0.5*u[0],c[0]+0.5*u[0],q0,q1))
    mm=min(segdist(mc[0,i]-0.75*md[0,i],mc[0,i]+0.75*md[0,i],mc[0,j]-0.75*md[0,j],mc[0,j]+0.75*md[0,j]) for i in range(3) for j in range(i+1,3))
    last=int(np.nonzero(tot>0)[0].max()) if (tot>0).any() else 0
    idx=[min(n-1,int(round(t/fm))) for t in TT]
    out["runs"]["tris%d"%T]=dict(
        release=dict(myosin_to_x=[round(float(math.degrees(math.acos(min(1,abs(v[0]))))),1) for v in md[0]],
                     spread=round(float(sp[0]),1),actin_tilt_x=round(float(math.degrees(math.acos(min(1,abs(u[0,0]))))),1),
                     actin_to_each=[round(float(a),1) for a in am[0]],perp_nm=[round(1000*p) for p in perp],
                     actin_clear_nm=[round(1000*x) for x in clr],myo_pair_clear_nm=round(1000*mm),heads=[int(x) for x in per[0]]),
        am=[round(float(am[i].mean()),1) for i in idx],worst=[round(float(am[i].max()),1) for i in idx],
        spread=[round(float(sp[i]),1) for i in idx],heads=[int(tot[i]) for i in idx],
        myos_engaged=[int((per[i]>0).sum()) for i in idx],
        attached_until=(None if tot[-1]>0 else round(last*fm)))
json.dump(out,open("tris_summary.json","w"),indent=1)
print("drawn scaffold: myosins %s deg to x, pairwise %s"%(drawn["to_axis"],[round(x,1) for x in drawn["pair"]]))
for T in (0,15,30,45):
    r=out["runs"]["tris%d"%T]; g=r["release"]
    print("\nactin tilt %d deg -- AT RELEASE: myosins %s deg to x, spread %.1f, actin %.1f to x, actin-to-each %s,"
          " actin-myosin axis distance %s nm, clearance %s nm, myo-myo clearance %d nm, heads per myosin %s"
          %(T,g["myosin_to_x"],g["spread"],g["actin_tilt_x"],g["actin_to_each"],g["perp_nm"],g["actin_clear_nm"],g["myo_pair_clear_nm"],g["heads"]))
    print("   %-26s"%"t (ms)"+"".join("%7d"%t for t in TT))
    for lab,key in (("actin-myosin, mean of 3","am"),("  worst of 3","worst"),("myosin spread","spread"),("heads bound","heads"),("myosins engaged","myos_engaged")):
        print("   %-26s"%lab+"".join(("%7.1f" if isinstance(v,float) else "%7d")%v for v in r[key]))
    print("   attached until: %s"%("whole run" if r["attached_until"] is None else "%d ms"%r["attached_until"]))
