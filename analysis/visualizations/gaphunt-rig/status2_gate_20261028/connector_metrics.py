"""Comparable saved-pose overshoot, independent of the simulated force winner.

Depth is the greatest axial burial of ANY connector point within a fixed
130 nm radius, on the selected tip's own half of the finite myosin. Every
link/myosin/tip pair is considered. Units internally are um; output is nm.
"""
import argparse
import itertools
import json
import os
from pathlib import Path
import h5py
import numpy as np

REASONS = {1:"geometric_contour",2:"stochastic",3:"barbed_tip_or_historical_occlusion",
           4:"actn2_myosin_collision",5:"mature_connector_exclusion"}


def wrap(x, box, periodic):
    # C++ std::round convention, including half-box ties.
    scaled = x/box
    return x-box*np.copysign(np.floor(np.abs(scaled)+.5),scaled)*periodic


def depths(first, delta, axis, half, radius):
    """Clip to radial cylinder AND axial half-slab at the same segment point."""
    z = np.sum(first*axis,axis=1); dz = np.sum(delta*axis,axis=1)
    p = first-z[:,None]*axis; d = delta-dz[:,None]*axis
    a = np.sum(d*d,axis=1); b = np.sum(p*d,axis=1)
    c = np.sum(p*p,axis=1)-radius*radius
    discriminant = b*b-a*c
    divisor = np.where(a>0,a,1)
    root = np.sqrt(np.maximum(0,discriminant))
    lo = np.where(a>0,(-b-root)/divisor,0)
    hi = np.where(a>0,(-b+root)/divisor,1)
    ok = np.where(a>0,discriminant>=0,c<=0)
    denom = np.where(dz!=0,dz,1)
    x = -z/denom; y = (half-z)/denom
    lo = np.maximum(np.maximum(lo,0),np.where(dz!=0,np.minimum(x,y),0))
    hi = np.minimum(np.minimum(hi,1),np.where(dz!=0,np.maximum(x,y),1))
    ok &= (lo<=hi)&((dz!=0)|((z>=0)&(z<=half)))
    value = half-np.minimum(z+lo*dz,z+hi*dz)
    return np.where(ok,value,-np.inf)


def self_test():
    u=np.array([1.,0,0])
    assert depths(np.array([[.9,.05,0]]),np.array([[-.2,.25,0]]),u,.75,.13)[0] == -np.inf
    assert np.isclose(depths(np.array([[.7,.1,0]]),np.array([[.1,0,0]]),u,.75,.13)[0],.05)
    assert np.isfinite(depths(np.array([[.5,-.2,0]]),np.array([[0,.4,0]]),u,.75,.13)[0])
    random=np.random.default_rng(20261028)
    for _ in range(200):
        p=random.uniform(-1,1,3); d=random.uniform(-1,1,3)
        point=p+np.linspace(0,1,2001)[:,None]*d
        q=(point[:,0]>=0)&(point[:,0]<=.75)&(np.linalg.norm(point[:,1:],axis=1)<=.13)
        actual=depths(p[None,:],d[None,:],u,.75,.13)[0]
        if q.any(): assert actual >= np.max(.75-point[q,0])-1e-10


def text(value):
    return value.decode() if isinstance(value,bytes) else str(value)


def frame_geometry_reference(first, second, centers, axes, box, periodic, length):
    d=wrap(second-first,box,periodic); midpoint=first+.5*d
    reach=np.hypot(.5*length,.13)+.5*np.linalg.norm(d,axis=1)
    max_reach=float(np.max(reach,initial=0))
    ranges=[range(-int(np.floor(max_reach/size+.5)),int(np.floor(max_reach/size+.5))+1)
            if periodic[k] else (0,) for k,size in enumerate(box)]
    pairs=[]; mid_pairs=[]; buffered=np.zeros(len(first),bool)
    for center,axis in zip(centers,axes):
        relative=wrap(midpoint-center,box,periodic)
        best=[np.full(len(first),-np.inf),np.full(len(first),-np.inf)]
        midbest=[np.full(len(first),-np.inf),np.full(len(first),-np.inf)]
        for image in itertools.product(*ranges):
            rel=relative+np.asarray(image)*box
            near=np.flatnonzero(np.sum(rel*rel,axis=1)<=reach*reach)
            if not len(near): continue
            p=rel[near]-.5*d[near]
            for side,sign in enumerate((-1,1)):
                best[side][near]=np.maximum(best[side][near],depths(p,d[near],sign*axis,.5*length,.13))
                midbest[side][near]=np.maximum(midbest[side][near],depths(rel[near],np.zeros_like(p),sign*axis,.5*length,.13))
                buffered[near] |= np.isfinite(depths(p,d[near],sign*axis,.5*length-.01,.12))
        pairs.extend(best); mid_pairs.extend(midbest)
    return np.asarray(pairs),np.asarray(mid_pairs),buffered


def frame_geometry(first, second, centers, axes, box, periodic, length):
    """Batch myosins, retaining the scalar implementation as an oracle."""
    d=wrap(second-first,box,periodic); midpoint=first+.5*d
    reach=np.hypot(.5*length,.13)+.5*np.linalg.norm(d,axis=1)
    max_reach=float(np.max(reach,initial=0))
    ranges=[range(-int(np.floor(max_reach/size+.5)),int(np.floor(max_reach/size+.5))+1)
            if periodic[k] else (0,) for k,size in enumerate(box)]
    pairs=np.full((2*len(centers),len(first)),-np.inf)
    mid_pairs=np.full_like(pairs,-np.inf);buffered=np.zeros(len(first),bool)
    images=list(itertools.product(*ranges))
    for begin in range(0,len(centers),32):
        end=min(begin+32,len(centers))
        relative=wrap(midpoint[:,None,:]-centers[None,begin:end,:],box,periodic)
        for image in images:
            rel=relative+np.asarray(image)*box
            ai,mi=np.nonzero(np.sum(rel*rel,axis=2)<=reach[:,None]**2)
            if not len(ai): continue
            middle=rel[ai,mi];delta=d[ai];p=middle-.5*delta
            for side,sign in enumerate((-1,1)):
                axis=sign*axes[begin+mi];rows=2*(begin+mi)+side
                pairs[rows,ai]=np.maximum(pairs[rows,ai],depths(p,delta,axis,.5*length,.13))
                mid_pairs[rows,ai]=np.maximum(mid_pairs[rows,ai],depths(middle,np.zeros_like(p),axis,.5*length,.13))
                hit=np.isfinite(depths(p,delta,axis,.5*length-.01,.12))
                buffered |= np.bincount(ai,weights=hit,minlength=len(first))>0
    return pairs,mid_pairs,buffered


def batch_self_test():
    random=np.random.default_rng(20261028)
    for _ in range(12):
        box=np.array([2.,1.,1.]);periodic=np.array([True,True,True])
        p=random.uniform(-1,1,(40,3));q=random.uniform(-1,1,(40,3))
        c=random.uniform(-1,1,(4,3));u=random.normal(size=(4,3))
        u/=np.linalg.norm(u,axis=1,keepdims=True)
        x=frame_geometry_reference(p,q,c,u,box,periodic,1.5)
        y=frame_geometry(p,q,c,u,box,periodic,1.5)
        for a,b in zip(x,y): assert np.allclose(a,b,atol=1e-12,rtol=0)


def analyze(path, start_step, stride):
    with h5py.File(path,"r",locking=False) as f:
        a=f["meta"].attrs; dt=float(a["dt"])
        box=np.array([a["box_"+k] for k in "xyz"])
        periodic=np.array([bool(a.get("periodic_"+k,True)) for k in "xyz"])
        n=int(f["meta/committed_frames"][0]); steps=f["state/current_step"][:n,0]
        ec={key:i for i,key in enumerate(text(a["aa_local_link_events_v3_columns"]).split(","))}
        ev=f["catch_bond/local_link_events_v3"][:]
        off=f["catch_bond/local_link_frame_offsets"]
        identities=f["catch_bond/local_link_identity"]
        base=max(0,int(np.searchsorted(steps,start_step,side="right"))-1)
        o,c=map(int,off[base,1:3]); initial=identities[o:o+c]
        inherited=set(initial[initial[:,4]==2,1].astype(int))
        maturation=ev[(ev[:,ec["old_status"]]==1)&(ev[:,ec["new_status"]]==2)&(ev[:,ec["step"]]>start_step)]
        mat_step={}
        for row in maturation:
            lid=int(row[ec["link_id"]]);mat_step[lid]=min(mat_step.get(lid,np.inf),int(row[ec["step"]]))
        frames=[]; samples={name:[] for name in ("all","mature","inherited_mature","post_resume_mature","midpoint_all")}
        for frame in range(base+1,n,stride):
            o,c=map(int,off[frame,1:3]); ids=identities[o:o+c]
            if not c: continue
            sites=f["catch_bond/local_link_coordinates"][o:o+c]
            arm_states=f.get("catch_bond/local_link_arm_state")
            if arm_states is not None and arm_states.shape[0]==identities.shape[0]:
                states=arm_states[o:o+c,0]
                bridge=np.isin(states,[1,2])
                if np.all(np.isin(states,[1,2,3,4])):
                    ids=ids[bridge]; sites=sites[bridge]; c=len(ids)
                    if not c: continue
            elif bool(a.get("aa_arm_state_model",False)):
                raise ValueError("Arm-state trajectory lacks aligned bridge/tether diagnostics")
            ca=f["actin/center"][frame]; ua=f["actin/direction"][frame]
            assert np.all((ids[:,2:4]>=0)&(ids[:,2:4]<len(ca))),"bridge has invalid actin index"
            # Match material-site construction in compute_actn2_patch_geometry:
            # s is in um and multiplies the stored actin direction directly.
            first=ca[ids[:,2]]+sites[:,0,None]*ua[ids[:,2]]
            second=ca[ids[:,3]]+sites[:,1,None]*ua[ids[:,3]]
            mu=f["myosin/direction"][frame];mu/=np.linalg.norm(mu,axis=1,keepdims=True)
            pair,mid,buffered=frame_geometry(first,second,f["myosin/center"][frame],mu,box,periodic,float(a["myosin_length"]))
            mature=ids[:,4]==2
            old=mature&np.isin(ids[:,1],list(inherited))
            new=mature&np.array([mat_step.get(int(lid),np.inf)<=steps[frame] for lid in ids[:,1]])
            record=dict(frame=frame,step=int(steps[frame]),time_s=float(steps[frame]*dt))
            for name,values,mask in (("all",pair,np.ones(c,bool)),("mature",pair,mature),("inherited_mature",pair,old),
                                     ("post_resume_mature",pair,new),("midpoint_all",mid,np.ones(c,bool))):
                q=values[:,mask]; penetration=q>.001 # fixed 1 nm reporting threshold
                positive=q[penetration]*1000
                samples[name].extend(positive.tolist())
                record[name]=dict(live_links=int(mask.sum()),penetrating_link_tip_pairs=int(penetration.sum()),
                    penetrating_links=int(np.any(penetration,axis=0).sum()),penetrating_tips=int(np.any(penetration,axis=1).sum()),
                    buffered_exclusion_links=int(np.sum(buffered&mask)))
            frames.append(record)
        death=ev[(ev[:,ec["new_status"]]==0)&(ev[:,ec["step"]]>start_step)]
        cohort={"all":death,"inherited_mature":death[np.isin(death[:,ec["link_id"]],list(inherited))],
                "post_resume_mature":death[np.array([mat_step.get(int(r[ec["link_id"]]),np.inf)<=r[ec["step"]] for r in death],bool)]}
        reasons={name:{str(k):int(np.sum(rows[:,ec["detachment_reason"]]==k)) for k in REASONS}
                 for name,rows in cohort.items()}
        summary={}
        for name,values in samples.items():
            summary[name]=dict(pair_samples=len(values),
                depth_nm_p50_p90_max=np.percentile(values,[50,90,100]).tolist() if values else None,
                mean_penetrating_links_per_frame=float(np.mean([r[name]["penetrating_links"] for r in frames])),
                mean_penetrating_link_tip_pairs_per_frame=float(np.mean([r[name]["penetrating_link_tip_pairs"] for r in frames])))
        return dict(source=str(path),window_steps=[start_step,int(steps[-1])],frame_stride=stride,
                    summary=summary,death_counts=reasons,frames=frames)


def main():
    parser=argparse.ArgumentParser(); parser.add_argument("paths",nargs="*")
    parser.add_argument("--output",required=True);parser.add_argument("--start-step",type=int,default=200000)
    parser.add_argument("--stride",type=int,default=1);args=parser.parse_args()
    assert os.environ.get("SLURM_JOB_ID"),"Run trajectory analysis in Slurm"
    self_test()
    batch_self_test()
    root=Path(os.environ["SCRATCH"])/"sarcomere/gaphunt_20261023/data"
    paths=args.paths or [root/(arm+"_seed20260805.h5") for arm in ("r27zt0","r27zt0np","r27zt0npg130")]
    out=dict(metric_version="full_connector_all_pairs_fixed_130nm_v1",radius_nm=130,reporting_depth_threshold_nm=1,
        definition="Maximum axial burial over connector points inside a fixed radial cylinder on each tip's own physical half; all link/myosin/tip pairs, independent of force eligibility and winner.",
        caveats=["Saved poses follow integration; instantaneous overlaps may appear before the next deletion check.",
                 "Historical reason-3 events retain their original meaning; old ambiguous mature deletions cannot be retroactively separated.",
                 "The same link may contribute to several myosin-tip pairs; pair counts and unique-link counts are separate."],
        reason_codes=REASONS,arms={})
    for path in paths:
        out["arms"][Path(path).stem]=analyze(path,args.start_step,args.stride)
        print("CONNECTOR_METRICS",Path(path).stem,out["arms"][Path(path).stem]["summary"],flush=True)
    Path(args.output).write_text(json.dumps(out,indent=2,allow_nan=False)+"\n")


if __name__=="__main__":main()
