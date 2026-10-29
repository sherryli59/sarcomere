#!/usr/bin/env python3
"""Standalone ACTN2 proposal: geometry, conservative contact, and failure.

Independent of simulator files and trajectory data. Units: um, pN, seconds.
Requires numpy, scipy, matplotlib. Run numerical work in a Slurm allocation.
All constants are illustrative, not a production calibration. See report.md.
"""
from __future__ import annotations
import argparse
import base64
from dataclasses import asdict, dataclass, replace
from functools import lru_cache
import json
import math
from pathlib import Path
import platform
import time
import os
import hashlib

import numpy as np
import scipy
from numpy.polynomial.legendre import leggauss
from scipy.special import erf
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize


@dataclass(frozen=True)
class Model:
    alpha: float = 3.
    k: float = 1000.
    slope: float = 30.
    epsilon: float = .005
    radius: float = .130
    half: float = .750
    approach: float = .050
    core_radius: float = .055  # 50 nm myosin + 5 nm connector, capsule reference
    barbed_zone: float = .200
    filament_length: float = 1.
    temperature: float = 310.

    @property
    def kbt(self):
        return 1.380649e-5*self.temperature  # pN um


@lru_cache(None)
def gauss(n):
    t,w=leggauss(n)
    return (t+1)/2,w/2


def point_energy(x, model=Model(), center=(0.,0.,0.), axis=(1.,0.,0.)):
    x=np.asarray(x,dtype=float)-np.asarray(center)
    u=np.asarray(axis,dtype=float);u=u/np.linalg.norm(u)
    z=np.sum(x*u,axis=-1)
    r=x-z[...,None]*u
    rho=np.sqrt(np.sum(r*r,axis=-1)+model.epsilon**2)
    za=np.sqrt(z*z+model.epsilon**2)
    aa=model.alpha*(model.radius-rho)
    bb=model.half+model.approach-za
    dd=np.sqrt((aa-bb)**2+model.epsilon**2)
    depth=np.maximum((aa+bb-dd)/2,0.)
    ps=model.slope/model.k
    energy=np.where(depth<=ps,.5*model.k*depth**2,
                    model.slope*(depth-.5*ps))
    wa=.5*(1-(aa-bb)/dd);wb=1-wa
    force=np.minimum(model.k*depth,model.slope)[...,None]*(
        wb[...,None]*(z/za)[...,None]*u+
        model.alpha*wa[...,None]*r/rho[...,None])
    return energy,force


def contact(p0,p1,model=Model(),center=(0.,0.,0.),axis=(1.,0.,0.),n=32):
    """Normalized segment integral and gradients; arbitrary nonperiodic pose."""
    p0,p1=np.broadcast_arrays(np.asarray(p0,dtype=float),np.asarray(p1,dtype=float))
    c=np.broadcast_to(np.asarray(center),p0.shape)
    ts,ws=gauss(n)
    x=p0[...,None,:]+ts[:,None]*(p1-p0)[...,None,:]
    e,f=point_energy(x,model,c[...,None,:],axis)
    fi=np.sum(f*(ws*(1-ts))[:,None],axis=-2)
    fj=np.sum(f*(ws*ts)[:,None],axis=-2)
    torque=-np.sum(np.cross(x-c[...,None,:],f)*ws[:,None],axis=-2)
    return np.sum(e*ws,axis=-1),fi,fj,torque


def cylinder_intersects(p0,p1,radius=.130,half=.750):
    """Exact same-point segment/cylinder intersection in the myosin frame."""
    p0,p1=np.broadcast_arrays(np.asarray(p0,dtype=float),np.asarray(p1,dtype=float))
    d=p1-p0;z=p0[...,0];dz=d[...,0]
    moving=np.abs(dz)>1.e-15
    safe=np.where(moving,dz,1.)
    ta=(-half-z)/safe;tb=(half-z)/safe
    lo=np.where(moving,np.maximum(0,np.minimum(ta,tb)),0.)
    hi=np.where(moving,np.minimum(1,np.maximum(ta,tb)),1.)
    valid=(lo<=hi)&(moving|(np.abs(z)<=half))
    den=np.sum(d[...,1:]**2,axis=-1)
    t=np.clip(-np.sum(p0[...,1:]*d[...,1:],axis=-1)/np.where(den>0,den,1),lo,hi)
    r=p0[...,1:]+t[...,None]*d[...,1:]
    return valid&(np.sum(r*r,axis=-1)<=radius**2)


def capsule_intersects(p0,p1,model=Model()):
    """Physical-core veto: cylinder plus spherical endcaps, connector inflated."""
    p0,p1=np.asarray(p0),np.asarray(p1);d=p1-p0
    hit=cylinder_intersects(p0,p1,model.core_radius,model.half)
    den=np.sum(d*d,axis=-1)
    for sign in (-1,1):
        v=p0-np.array([sign*model.half,0,0])
        t=np.clip(-np.sum(v*d,axis=-1)/np.where(den>0,den,1),0,1)
        hit=hit|(np.sum((v+t[...,None]*d)**2,axis=-1)<=model.core_radius**2)
    return hit


def maturation_allowed(p0,p1,b=(.1,.1),registry=True,model=Model()):
    b=np.asarray(b)
    return (registry & np.all((b>=0)&(b<=model.barbed_zone),axis=-1)
            & ~cylinder_intersects(p0,p1,model.radius,model.half))


def association_factor(p0,p1,b=(.1,.1),old_pair_energy=0.,model=Model()):
    # Maximum pair budget: adding a redundant link changes max(E_old,E_new).
    energy=contact(p0,p1,model,n=64)[0]
    delta=np.maximum(energy-old_pair_energy,0.)
    valid=np.all((np.asarray(b)>=0)&(np.asarray(b)<=model.barbed_zone),axis=-1)
    valid=valid&~capsule_intersects(p0,p1,model)
    return np.where(valid,np.exp(-delta/model.kbt),0.)


def pair_budget(p0,p1,pair_ids,model=Model(),n=64):
    e,fi,fj,tm=contact(p0,p1,model,n=n)
    selected=[]
    for key in sorted(set(pair_ids)):
        ids=np.flatnonzero(np.asarray(pair_ids)==key)
        selected.append(ids[np.argmax(e[ids])])  # stable first-index tie
    selected=np.asarray(selected)
    return dict(energy=float(e[selected].sum()),
        net_force=np.sum(fi[selected]+fj[selected],axis=0),selected=selected)


def edge_hazard(b,model=Model(),rate=100.,length=.020):
    return rate*np.sum((np.maximum(np.asarray(b)-model.barbed_zone,0)/length)**2,axis=-1)


def contact_hazard(fi,fj,rate=50.,scale=5.):
    traction=np.maximum(np.linalg.norm(fi,axis=-1),np.linalg.norm(fj,axis=-1))
    # No arbitrary exponent or force cap; this demo's forces keep exp finite.
    return rate*np.expm1(traction/scale)


def sample_channels(rates,dt,rng):
    """Exact frozen-rate competing exponential clocks within a substep."""
    rates=np.asarray(rates)
    thresholds=rng.exponential(size=rates.shape)
    waits=np.divide(thresholds,rates,out=np.full_like(rates,np.inf),where=rates>0)
    channel=np.argmin(waits,axis=-1)
    fired=np.min(waits,axis=-1)<dt
    return np.where(fired,channel,-1)


def endpoint_crossing(old,new,diffusion,dt,rng,length=1.):
    hit=(old<=0)|(old>=length)|(new<=0)|(new>=length)
    if diffusion>0:
        # Brownian bridge correction for paths whose sampled endpoints survive.
        good=(old>0)&(new>0)&(old<length)&(new<length)
        p0=np.exp(np.minimum(-old*new/(diffusion*dt),0.))
        p1=np.exp(np.minimum(-(length-old)*(length-new)/(diffusion*dt),0.))
        hit|=good&((rng.random(old.shape)<p0)|(rng.random(old.shape)<p1))
    return hit


def geometry_fixtures(model):
    fixtures=[
        ('outside end',[.840,.080,0],[.840,.120,0],True,True),
        ('midpoint misses connector',[.700,.115,0],[.700,.175,0],False,True),
        ('radial and axial at different points',[.700,.200,0],[.800,.100,0],True,True),
        ('outside radial boundary',[.700,.140,0],[.700,.180,0],True,True),
        ('within approach, outside physical end',[.780,.080,0],[.780,.100,0],True,True),
        ('crosses physical core',[0,-.080,0],[0,.080,0],False,False),
        ('exact physical end plane',[.750,.080,0],[.750,.120,0],False,True),
        ('radial tangency point',[.700,.130,0],[.700,.130,0],False,True),
        ('cylinder corner point',[.750,.130,0],[.750,.130,0],False,True),
    ]
    rows=[]
    for label,p0,p1,expected_mature,expected_core in fixtures:
        allowed=bool(maturation_allowed(p0,p1,model=model))
        core_free=not bool(capsule_intersects(p0,p1,model))
        assert allowed==expected_mature,(label,allowed)
        assert core_free==expected_core,(label,core_free)
        e,fi,fj,_=contact(p0,p1,model,n=128)
        rows.append(dict(label=label,p0_nm=(np.array(p0)*1000).tolist(),
            p1_nm=(np.array(p1)*1000).tolist(),maturation_allowed=allowed,
            physical_core_free=core_free,weak_association_factor=float(association_factor(p0,p1,model=model)),
            energy_pNnm=float(1000*e),net_force_pN=(fi+fj).tolist()))
    assert not bool(maturation_allowed(fixtures[0][1],fixtures[0][2],b=(.201,.1),model=model))
    assert not bool(maturation_allowed(fixtures[0][1],fixtures[0][2],registry=False,model=model))
    return rows


def check_model(model,rng):
    checks={};step=1.e-7
    # Endpoint and myosin translation/rotation derivatives, arbitrary pose.
    p0=rng.uniform([.65,.065,-.02],[.83,.145,.02],(150,3))
    p1=p0+rng.uniform([-.015,.015,-.01],[.015,.065,.01],(150,3))
    center=np.array([.008,-.003,.004]);axis=np.array([1,.03,-.02]);axis/=np.linalg.norm(axis)
    e,fi,fj,tm=contact(p0,p1,model,center,axis,n=64)
    endpoint_error=0.;myosin_error=0.;torque_error=0.
    for j in range(3):
        dx=np.eye(3)[j]*step
        fa=-(contact(p0+dx,p1,model,center,axis,64)[0]-contact(p0-dx,p1,model,center,axis,64)[0])/(2*step)
        fb=-(contact(p0,p1+dx,model,center,axis,64)[0]-contact(p0,p1-dx,model,center,axis,64)[0])/(2*step)
        fm=-(contact(p0,p1,model,center+dx,axis,64)[0]-contact(p0,p1,model,center-dx,axis,64)[0])/(2*step)
        endpoint_error=max(endpoint_error,float(np.max(np.abs(fa-fi[:,j]))),float(np.max(np.abs(fb-fj[:,j]))))
        myosin_error=max(myosin_error,float(np.max(np.abs(fm+fi[:,j]+fj[:,j]))))
        da=np.cross(np.eye(3)[j],axis)*1.e-6
        torque_fd=-(contact(p0,p1,model,center,axis+da,64)[0]-contact(p0,p1,model,center,axis-da,64)[0])/(2.e-6)
        torque_error=max(torque_error,float(np.max(np.abs(torque_fd-tm[:,j]))))
    balance=np.cross(p0-center,fi)+np.cross(p1-center,fj)+tm
    assert endpoint_error<1.e-3 and myosin_error<1.e-3 and torque_error<1.e-4
    assert np.max(np.abs(balance))<1.e-10
    checks['gradients']=dict(endpoint_max_error_pN=endpoint_error,myosin_max_error_pN=myosin_error,
        orientation_max_error_pNum=torque_error,torque_balance_max_pNum=float(np.max(np.abs(balance))))
    point=np.array([.730,.080,0.]);pe,pf=point_energy(point,model)
    ce,cf,cg,_=contact(point,point,model,n=64)
    assert np.isclose(pe,ce) and np.allclose(cf,pf/2) and np.allclose(cg,pf/2)
    checks['zero_length_connector']='point energy and equal endpoint force shares passed'
    xx=rng.uniform([-.85,-.18,-.18],[.85,.18,.18],(10000,3))
    for alpha in (1.,3.):
        mm=replace(model,alpha=alpha);_,f=point_energy(xx,mm)
        assert np.max(np.abs(f[:,0]))<=mm.slope+1e-10
        assert np.max(np.linalg.norm(f,axis=1))<=alpha*mm.slope+1e-10
    checks['force_bounds']='passed for 10,000 points at alpha=1 and alpha=3'
    # Dense independent same-point reference; avoid measure-zero tangencies.
    aa=rng.uniform([.60,-.22,-.1],[.90,.22,.1],(1200,3))
    bb=aa+rng.uniform(-.15,.15,(1200,3));ts=np.linspace(0,1,3001)
    dense=aa[:,None,:]+ts[None,:,None]*(bb-aa)[:,None,:]
    sampled=np.any((np.abs(dense[...,0])<=model.half)&
        (np.sum(dense[...,1:]**2,axis=-1)<=model.radius**2),axis=1)
    exact=cylinder_intersects(aa,bb,model.radius,model.half)
    assert np.array_equal(exact,sampled)
    checks['formation_geometry']='1,200 random segments agree with 3,001-point independent reference'
    pp=np.array([[.730,.080,0],[.730,.080,0],[.730,.080,0]])
    qq=pp+[0,.04,0]
    single=pair_budget(pp[:1],qq[:1],[1],model)
    duplicate=pair_budget(pp,qq,[1,1,1],model)
    distinct=pair_budget(pp,qq,[1,2,3],model)
    assert np.isclose(single['energy'],duplicate['energy'])
    assert np.allclose(single['net_force'],duplicate['net_force'])
    assert np.isclose(distinct['energy'],3*single['energy'])
    factor=association_factor(pp[0],qq[0],old_pair_energy=single['energy'],model=model)
    assert np.isclose(factor,1.)
    checks['pair_budget']='redundant-link invariance; distinct-pair additivity; incremental birth penalty passed'
    checks['closed_loop_work_pNnm']={}
    # Translate two redundant connectors around a loop, selecting max energy.
    for count in (128,512,2048):
        corners=np.array([[.69,.075,0],[.84,.075,0],[.84,.15,0],[.69,.15,0],[.69,.075,0]])
        path=np.concatenate([a+(b-a)*np.arange(count)[:,None]/count for a,b in zip(corners[:-1],corners[1:])])
        dp=np.roll(path,-1,axis=0)-path;mid=path+dp/2
        e0,f0,g0,_=contact(mid+[-.008,-.004,0],mid+[-.008,.036,0],model,n=64)
        e1,f1,g1,_=contact(mid+[.008,.004,0],mid+[.008,.044,0],model,n=64)
        force=np.where((e0>=e1)[:,None],f0+g0,f1+g1)
        work=float(np.sum(force*dp)*1000)
        checks['closed_loop_work_pNnm'][str(count)]=work
    assert abs(work)<.02,work
    # Explicitly show that finite outward connectors can still push laterally.
    probes=[]
    for alpha in (1.,3.):
        mm=replace(model,alpha=alpha)
        for ell in (0.,.040,.080):
            e,fa,fb,_=contact([.730,.080,0],[.730,.080+ell,0],mm,n=256)
            probes.append(dict(alpha=alpha,length_nm=ell*1000,burial_nm=20.,
                energy_pNnm=float(1000*e),force_pN=(fa+fb).tolist()))
    checks['force_direction_probes']=probes
    return checks


def rig_initial(replicates):
    theta=np.arange(4)*np.pi/2
    radial=np.stack([np.zeros(4),np.cos(theta),np.sin(theta)],axis=-1)
    p0=.080*radial+[.750,0,0];p1=.120*radial+[.750,0,0]
    return (np.tile(p0,(replicates,1,1)),np.tile(p1,(replicates,1,1)),
            np.tile([.025,0,0],(replicates,1)))


def simulate_rig(model,drive,duration=.0006,dt=1.e-7,replicates=1,
                 failure=False,anchor=1000.,seed=20261028,nquad=16):
    """Four distinct radial pairs; translated endpoints and myosin, no rotations.

    Harmonic anchors represent an explicit external elastic scaffold, not
    the simulated ACTN2 network. Mechanics is deterministic; chemical death
    uses independent exponential integrated-hazard thresholds. No births.
    """
    rng=np.random.default_rng(seed)
    p0,p1,c=rig_initial(replicates);rest0=p0.copy();rest1=p1.copy()
    alive=np.ones((replicates,4),bool);death=np.zeros((replicates,4),int)
    mu_end,mu_m=70.,19.;ks=1000.;rest_length=.040
    hazard=np.zeros((replicates,4,2));threshold=rng.exponential(size=hazard.shape)
    count=round(duration/dt);stride=max(1,count//200)
    traces={k:[] for k in ('time','center','p0','p1','fi','fj','alive','energy','penetration','radial_shift','death',
        'connector_axial_overlap')}
    dis=extwork=0.;balance=[]

    def forces():
        e,fi,fj,tm=contact(p0,p1,model,c[:,None,:],n=nquad)
        fi=fi*alive[...,None];fj=fj*alive[...,None]
        delta=p1-p0;ell=np.linalg.norm(delta,axis=-1)
        tension=ks*(ell-rest_length)*alive
        spring=tension[...,None]*delta/ell[...,None]
        f0=fi+spring-anchor*(p0-rest0)
        f1=fj-spring-anchor*(p1-rest1)
        fm=-np.sum(fi+fj,axis=1)+np.array([drive,0,0])
        energy=np.sum(e*alive+.5*ks*(ell-rest_length)**2*alive+
            .5*anchor*np.sum((p0-rest0)**2+(p1-rest1)**2,axis=-1),axis=1)
        return fi,fj,f0,f1,fm,energy,ell

    fi,fj,f0,f1,fm,energy,ell=forces();energy_initial=float(energy[0])
    for step in range(count+1):
        fi,fj,f0,f1,fm,energy,ell=forces()
        if step%stride==0 or step==count:
            values=dict(time=step*dt,center=c.copy(),p0=p0.copy(),p1=p1.copy(),fi=fi.copy(),fj=fj.copy(),
                alive=alive.copy(),energy=energy.copy(),penetration=c[:,0]+model.half-.750,
                radial_shift=np.max(np.linalg.norm(p1[...,1:],axis=-1)-.120,axis=1),death=death.copy(),
                connector_axial_overlap=c[:,None,0]+model.half-np.minimum(p0[...,0],p1[...,0]))
            for key in traces:traces[key].append(values[key])
        if step==count:break
        if failure:
            rate=np.stack([np.full_like(ell,2.),contact_hazard(fi,fj)],axis=-1)
            hazard+=rate*dt
            crossed=hazard>=threshold
            fired=np.any(crossed,axis=-1)&alive
            # Pick earliest crossing fraction when both clocks cross this step.
            fraction=np.divide(threshold-(hazard-rate*dt),rate*dt,
                out=np.full_like(hazard,np.inf),where=rate>0)
            channel=np.argmin(fraction,axis=-1)+1
            contour=(ell>.120)&alive
            death=np.where(fired,channel,death);death=np.where(contour,3,death)
            alive&=~(fired|contour)
            if np.any(fired|contour):
                fi,fj,f0,f1,fm,energy,ell=forces()  # no stale forces after death
        if not failure and replicates==1:
            extwork+=drive*mu_m*float(fm[0,0])*dt
            dis+=float(mu_end*np.sum(f0*f0+f1*f1)+mu_m*np.sum(fm*fm))*dt
        p0+=mu_end*f0*dt;p1+=mu_end*f1*dt;c+=mu_m*fm*dt
    out={key:np.asarray(value) for key,value in traces.items()}
    summary=dict(alpha=model.alpha,drive_pN=drive,dt_s=dt,duration_s=duration,
        replicates=replicates,chemical_failure=failure,anchor_pN_per_um=anchor,
        endpoint_mobility_um_per_pNs=mu_end,myosin_mobility_um_per_pNs=mu_m,
        final_penetration_nm=(out['penetration'][-1]*1000).tolist(),
        penetration_reference='myosin +tip relative to the original scaffold plane X=750 nm',
        final_connector_axial_overlap_nm=(out['connector_axial_overlap'][-1]*1000).tolist(),
        final_alive_links=out['alive'][-1].sum(axis=1).tolist(),
        final_myosin_lateral_shift_nm=(np.linalg.norm(c[:,1:],axis=1)*1000).tolist(),
        passed_with_live_connectors=int(np.sum((out['penetration'][-1]>0)&alive.any(axis=1))),
        final_max_radial_shift_nm=(out['radial_shift'][-1]*1000).tolist(),
        any_death_fraction=float(np.mean(np.any(~alive,axis=1))),
        all_dead_fraction=float(np.mean(~np.any(alive,axis=1))),
        final_overshoot_fraction=float(np.mean(out['penetration'][-1]>0)),
        death_counts={str(j):int(np.sum(death==j)) for j in (1,2,3)},
        death_channels={'1':'baseline','2':'contact','3':'contour'})
    if not failure and replicates==1:
        summary['work_balance_residual_pNnm']=1000*(float(out['energy'][-1,0])-energy_initial-extwork+dis)
    return out,summary


def weak_sites(model,kind,reflecting=False,replicates=1024,dt=1.e-5,duration=.010,seed=1):
    """Fixed antiparallel actins; true contact-gradient motion along their axes.

    Old reflecting boundaries are a counterfactual reference. The target
    stage-B model has baseline + soft-edge loss, no contact slip channel yet.
    """
    rng=np.random.default_rng(seed);ds=.010;mob=ds/model.kbt
    initial=np.array([.008,.100] if kind=='barbed' else [.100,.190])
    b=np.tile(initial,(replicates,1));orientation=np.array([[1.,0,0],[-1.,0,0]])
    sites=np.array([[.770,.080,0],[.770,.120,0]])
    barbed=sites+initial[:,None]*orientation
    alive=np.ones(replicates,bool);reason=np.zeros(replicates,int)
    nsteps=round(duration/dt);stride=max(1,nsteps//200)
    times=[];survival=[];histories=[];reason_history=[]
    for step in range(nsteps+1):
        if step%stride==0 or step==nsteps:
            times.append(step*dt);survival.append(np.mean(alive));histories.append(b[:12].copy())
            reason_history.append([int(np.sum(reason==j)) for j in (1,2,3)])
        if step==nsteps:break
        ids=np.flatnonzero(alive)
        if len(ids)==0:continue
        positions=barbed-b[ids,:,None]*orientation
        _,fi,fj,_=contact(positions[:,0],positions[:,1],model,n=16)
        ff=np.stack([fi,fj],axis=1)
        # p=barbed-b*u -> dU/db=F_endpoint dot u.
        drift=-mob*np.sum(ff*orientation,axis=-1)
        proposal=b[ids]+drift*dt+np.sqrt(2*ds*dt)*rng.normal(size=(len(ids),2))
        if reflecting:
            folded=np.mod(proposal,2*model.barbed_zone)
            b[ids]=np.where(folded<=model.barbed_zone,folded,2*model.barbed_zone-folded)
            rates=np.stack([np.full(len(ids),10.),np.zeros(len(ids))],axis=-1)
            crossed=np.zeros(len(ids),bool)
        else:
            crossed=np.any(endpoint_crossing(b[ids],proposal,ds,dt,rng,model.filament_length),axis=1)
            reason[ids[crossed]]=1;alive[ids[crossed]]=False
            b[ids[~crossed]]=proposal[~crossed]
            rates=np.stack([np.full(len(ids),10.),edge_hazard(b[ids],model)],axis=-1)
        event=sample_channels(rates,dt,rng)
        fired=(event>=0)&~crossed
        reason[ids[fired]]=np.where(event[fired]==0,2,3);alive[ids[fired]]=False
    summary=dict(kind=kind,reflecting_reference=reflecting,replicates=replicates,
        dt_s=dt,alive_fraction=float(np.mean(alive)),endpoint_loss=int(np.sum(reason==1)),
        baseline_loss=int(np.sum(reason==2)),soft_edge_loss=int(np.sum(reason==3)),
        surviving_outside_zone=int(np.sum(alive&np.any(b>model.barbed_zone,axis=1))),
        final_mean_alive_b_nm=(1000*np.mean(b[alive],axis=0)).tolist() if np.any(alive) else None)
    return dict(time=np.array(times),survival=np.array(survival),history=np.array(histories),
        reason_history=np.array(reason_history),b=b,alive=alive,summary=summary)


def statistical_checks(model,rng):
    data={};checks={}
    # Pure endpoint first-passage test; no chemical death in this analytic fixture.
    n=20000;dt=1.e-5;duration=.005;diff=.010
    b=np.full(n,.010);live=np.ones(n,bool)
    for _ in range(round(duration/dt)):
        ids=np.flatnonzero(live)
        proposal=b[ids]+np.sqrt(2*diff*dt)*rng.normal(size=len(ids))
        dead=endpoint_crossing(b[ids],proposal,diff,dt,rng)
        live[ids[dead]]=False;b[ids[~dead]]=proposal[~dead]
    expected=float(erf(.010/np.sqrt(4*diff*duration)))
    observed=float(np.mean(live));sigma=np.sqrt(expected*(1-expected)/n)
    assert abs(observed-expected)<6*sigma
    checks['absorbing_first_passage']=dict(observed_survival=observed,analytic_survival=expected,
        binomial_standard_error=sigma,replicates=n)
    # The far physical endpoint has the same first-passage law by reflection.
    b=np.full(n,.990);live=np.ones(n,bool)
    for _ in range(round(duration/dt)):
        ids=np.flatnonzero(live)
        proposal=b[ids]+np.sqrt(2*diff*dt)*rng.normal(size=len(ids))
        dead=endpoint_crossing(b[ids],proposal,diff,dt,rng)
        live[ids[dead]]=False;b[ids[~dead]]=proposal[~dead]
    observed=float(np.mean(live))
    assert abs(observed-expected)<6*sigma
    checks['far_endpoint_first_passage']=dict(observed_survival=observed,
        analytic_survival=expected,binomial_standard_error=sigma,replicates=n)
    # Constant endpoint traction clamps: Monte Carlo vs analytic survival.
    times=np.linspace(0,.020,101);tractions=np.array([0.,5.,15.,30.]);rates=2+50*np.expm1(tractions/5.)
    survival=[];theory=[];causes=[]
    for load,lam in zip(tractions,rates):
        waits_base=rng.exponential(1/2,size=n)
        contact_rate=50*np.expm1(load/5)
        waits_contact=(rng.exponential(1/contact_rate,size=n) if contact_rate>0 else np.full(n,np.inf))
        deaths=np.minimum(waits_base,waits_contact)
        measured=np.mean(deaths[:,None]>times[None,:],axis=0);analytic=np.exp(-lam*times)
        se=np.sqrt(analytic*(1-analytic)/n)
        assert np.all(np.abs(measured-analytic)<6*se+2/n)
        survival.append(measured);theory.append(analytic)
        causes.append(float(np.mean(waits_contact<waits_base)))
    data['failure']=dict(time=times,traction=tractions,rates=rates,
        survival=np.array(survival),theory=np.array(theory))
    checks['failure_clamps']=dict(replicates_per_load=n,traction_pN=tractions.tolist(),
        rates_per_s=rates.tolist(),contact_first_fraction=causes,
        max_survival_error=float(np.max(np.abs(np.array(survival)-theory))))
    # Exercise the actual competing-channel sampler used by weak-site kinetics.
    channel_rates=np.array([10.,30.,60.]);channel_dt=.01
    events=sample_channels(np.broadcast_to(channel_rates,(100000,3)),channel_dt,rng)
    probabilities=np.r_[np.exp(-channel_rates.sum()*channel_dt),
        -np.expm1(-channel_rates.sum()*channel_dt)*channel_rates/channel_rates.sum()]
    measured=np.array([np.mean(events==j) for j in (-1,0,1,2)])
    stderr=np.sqrt(probabilities*(1-probabilities)/len(events))
    assert np.all(np.abs(measured-probabilities)<6*stderr)
    assert np.all(sample_channels(np.zeros((100,3)),channel_dt,rng)==-1)
    assert contact_hazard(np.zeros(3),np.zeros(3))==0
    checks['competing_channels']=dict(expected=probabilities.tolist(),observed=measured.tolist(),
        realizations=len(events),zero_rate_and_zero_contact_checks=True)
    assert edge_hazard([.2,.1],model)==0
    assert np.isclose(edge_hazard([.22,.1],model),100.)
    assert edge_hazard([.24,.1],model)>edge_hazard([.22,.1],model)
    # Real weak-birth geometry at the physical tip, moving toward onset.
    centers=np.array([[x,.080,0] for x in (.780,.785,.790,.795,.800,.805,.810,.820)])
    factors=association_factor(centers,centers+[0,.040,0],model=model)
    accepted=np.mean(rng.random((n,len(factors)))<factors,axis=0)
    se=np.sqrt(factors*(1-factors)/n)
    assert np.all(np.abs(accepted-factors)<6*se+2/n)
    checks['weak_birth_sampling']=dict(x_nm=(centers[:,0]*1000).tolist(),
        predicted_factor=factors.tolist(),observed_fraction=accepted.tolist(),replicates_per_point=n)
    return data,checks


def save_figures(output,model,rigs,weak,stats):
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    x=np.linspace(.650,.835,371);y=np.linspace(.050,.165,231)
    xx,yy=np.meshgrid(x,y);p0=np.stack([xx,yy,np.zeros_like(xx)],axis=-1)
    p1=p0+[0,.040,0]
    gate=maturation_allowed(p0,p1,model=model)
    birth=association_factor(p0,p1,model=model)
    fig,axs=plt.subplots(2,2,figsize=(13,9.5),layout='constrained')
    extent=(1000*x[0],1000*x[-1],1000*y[0],1000*y[-1])
    ax=axs[0,0];ax.imshow(gate,origin='lower',extent=extent,aspect='auto',cmap='RdYlGn',vmin=0,vmax=1)
    ax.set_title('A · Maturation: exact full-connector gate\ngreen allowed / red rejected; 40 nm along +Y')
    # Full connector whose midpoint lies outside but near endpoint lies inside.
    ax.plot([700,700],[115,155],color='white',lw=4);ax.scatter([700],[135],c='black',s=25)
    ax.annotate('midpoint outside;\nconnector rejected',(700,135),xytext=(660,148),
        arrowprops={'arrowstyle':'->'},fontsize=9)
    ax=axs[0,1]
    im=ax.imshow(np.log10(np.maximum(birth,1e-12)),origin='lower',extent=extent,aspect='auto',cmap='viridis',vmin=-12,vmax=0)
    ax.set_title('B · Weak birth: incremental energy penalty\nlog10(rate / existing rate), alpha = 3')
    fig.colorbar(im,ax=ax,label='log10 association factor',shrink=.8)
    for alpha,ax in zip((1.,3.),axs[1]):
        mm=replace(model,alpha=alpha);_,fi,fj,_=contact(p0,p1,mm,n=32);ff=fi+fj
        axial=np.abs(ff[...,0]);radial=np.linalg.norm(ff[...,1:],axis=-1)
        fraction=np.divide(axial,axial+radial,out=np.zeros_like(axial),where=(axial+radial)>1e-8)
        im=ax.imshow(np.ma.masked_where(axial+radial<.01,fraction),origin='lower',extent=extent,
            aspect='auto',cmap='coolwarm_r',vmin=0,vmax=1)
        ax.set_title(f'{"C" if alpha==1 else "D"} · Finite-connector force, alpha = {alpha:g}\nblue axial / red radial; arrows proportional to force')
        qx,qy=np.meshgrid(np.linspace(.680,.800,7),np.linspace(.070,.135,5))
        qa=np.stack([qx,qy,np.zeros_like(qx)],axis=-1)
        _,a,b,_=contact(qa,qa+[0,.04,0],mm,n=64);f=a+b
        ax.quiver(1000*qx,1000*qy,.3*f[...,0],.3*f[...,1],angles='xy',scale_units='xy',scale=1,
            color='#172b3e',width=.004)
        ax.text(.02,.025,'Arrow scale: 0.3 nm/pN',transform=ax.transAxes,fontsize=8)
        fig.colorbar(im,ax=ax,label='|Fx| / (|Fx| + |Fradial|)',shrink=.8)
    for ax in axs.flat:
        ax.axvline(750,color='#333',ls=':',lw=1);ax.axvline(800,color='#333',ls='--',lw=1)
        ax.set(xlabel='Connector start X (nm)',ylabel='Connector start Y (nm)')
    fig.suptitle('Standalone proposal: formation and repulsion\n'
        'k = 1000 pN/µm · potential slope = 30 pN · smoothing = 5 nm · values are uncalibrated',fontsize=13)
    fig.savefig(output/'formation_repulsion.png',dpi=170);plt.close(fig)

    fig,axs=plt.subplots(2,3,figsize=(15,8.5),layout='constrained')
    colors={'alpha1':'#d77531','alpha3':'#2588ac','overdrive':'#a477b9','free':'#58a88c'}
    for key in ('alpha1','alpha3','overdrive'):
        v=rigs[key]['data'];label=rigs[key]['label']
        axs[0,0].plot(v['time']*1000,v['penetration'][:,0]*1000,label=label,color=colors[key])
    axs[0,0].axhline(0,color='#777',lw=1,ls=':')
    axs[0,0].set(title='A · Driven myosin, elastic scaffold\npositive = tip past the original scaffold plane',xlabel='Time (ms)',ylabel='Tip penetration (nm)')
    axs[0,0].legend(fontsize=8)
    for key in ('alpha1','alpha3'):
        v=rigs[key]['data'];axs[0,1].plot(v['time']*1000,v['radial_shift'][:,0]*1000,label=rigs[key]['label'],color=colors[key])
    axs[0,1].set(title='B · Sideways scaffold deformation\n4 distinct pairs, spring anchors explicitly present',xlabel='Time (ms)',ylabel='Largest outer-site radial shift (nm)')
    axs[0,1].legend(fontsize=8)
    v=rigs['free']['data'];axs[0,2].plot(v['time']*1000,v['energy'][:,0]*1000,color=colors['free'])
    axs[0,2].set(title='C · Free relaxation, no drive or anchors\ncontact + connector-spring energy',xlabel='Time (ms)',ylabel='Energy (pN nm)')
    failure=stats['failure']
    for j,load in enumerate(failure['traction']):
        col=plt.cm.viridis(j/3)
        axs[1,0].plot(failure['time']*1000,failure['theory'][j],color=col,label=f'{load:g} pN')
        axs[1,0].plot(failure['time'][::10]*1000,failure['survival'][j,::10],'.',color=col)
    axs[1,0].set(title='D · Failure under clamped contact traction\nlines: analytic; dots: 20,000 random lifetimes/load',xlabel='Time (ms)',ylabel='Survival probability',ylim=(-.02,1.04))
    axs[1,0].legend(fontsize=8)
    fail=rigs['failure']['data'];control=rigs['loaded_control']['data'];t=fail['time']*1000
    lo,med,hi=np.percentile(fail['penetration']*1000,[10,50,90],axis=1)
    axs[1,1].fill_between(t,lo,hi,color='#d45b6c',alpha=.2,label='failure: 10–90%')
    axs[1,1].plot(t,med,color='#c74359',label='failure: median')
    axs[1,1].plot(control['time']*1000,control['penetration'][:,0]*1000,color='#2588ac',label='fixed topology')
    axs[1,1].axhline(0,color='#777',lw=1,ls=':')
    axs[1,1].set(title='E · Contact failure coupled to mechanics\nalpha = 3, 40 pN drive, all realizations included',xlabel='Time (ms)',ylabel='Tip penetration (nm)')
    axs[1,1].legend(fontsize=8)
    counts=fail['alive'].sum(axis=-1)
    axs[1,2].plot(t,np.mean(counts,axis=1),color='#c74359',label='mean surviving links')
    axs[1,2].fill_between(t,*np.percentile(counts,[10,90],axis=1),color='#d45b6c',alpha=.2)
    axs[1,2].set(title='F · Molecular loss; no geometric kill\n4 initial links, no rebirth in this loading experiment',xlabel='Time (ms)',ylabel='Surviving links',ylim=(-.1,4.2))
    fig.suptitle('Standalone controlled mechanics and failure · not a production network\n'
        'Failure parameters: baseline 2/s, contact 50/s × [exp(T / 5 pN) − 1]',fontsize=13)
    fig.savefig(output/'mechanics_failure.png',dpi=170);plt.close(fig)

    fig,axs=plt.subplots(1,2,figsize=(12,4.8),layout='constrained')
    lateral=np.linalg.norm(fail['center'][...,1:],axis=-1)*1000
    lo,med,hi=np.percentile(lateral,[10,50,90],axis=1)
    axs[0].fill_between(t,lo,hi,color='#d45b6c',alpha=.2,label='10–90%')
    axs[0].plot(t,med,color='#c74359',label='median')
    axs[0].set(xlabel='Time (ms)',ylabel='Myosin lateral displacement (nm)',
        title='A · Asymmetric loss allows sideways escape')
    axs[0].legend()
    im=axs[1].scatter(fail['penetration'][-1]*1000,lateral[-1],
        c=counts[-1],cmap=plt.get_cmap('viridis',5),vmin=-.5,vmax=4.5,alpha=.75)
    axs[1].axvline(0,color='#777',ls=':')
    axs[1].set(xlabel='Final tip position past original scaffold plane (nm)',
        ylabel='Final myosin lateral displacement (nm)',
        title='B · Passing does not require losing every connector')
    fig.colorbar(im,ax=axs[1],ticks=range(5),label='Surviving connectors')
    fig.suptitle('40 pN drive · alpha = 3 · illustrative failure law · all realizations',fontsize=12)
    fig.savefig(output/'lateral_escape.png',dpi=170);plt.close(fig)

    fig,axs=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    for column,kind in enumerate(('barbed','inner')):
        for reflecting,color,label in ((True,'#d77531','reflecting reference'),(False,'#2588ac','absorbing + soft edge')):
            v=weak[(kind,reflecting)]
            axs[0,column].plot(v['time']*1000,v['survival'],color=color,label=label)
            which=0 if kind=='barbed' else 1
            survivors=v['b'][v['alive'],which]*1000
            if len(survivors):axs[1,column].hist(survivors,bins=35,histtype='step',density=True,color=color,label=label,lw=1.5)
        axs[0,column].set(title=f'{"Barbed endpoint" if kind=="barbed" else "Inner structural edge"}: survival',xlabel='Time (ms)',ylabel='Survival probability',ylim=(-.02,1.03))
        axs[0,column].legend(fontsize=9)
        axs[1,column].axvline(0 if kind=='barbed' else 200,color='#777',ls='--')
        axs[1,column].set(title='Final positions of surviving weak sites',xlabel='Distance inward from barbed tip (nm)',ylabel='Conditional density (1/nm)')
        axs[1,column].legend(fontsize=9)
    fig.suptitle('Weak-site boundary experiment: actual contact-gradient drift + diffusion\n'
        'No reflection in target model · baseline 10/s · soft edge 100/s × ([b − 200 nm]+ / 20 nm)²',fontsize=12)
    fig.savefig(output/'weak_boundaries.png',dpi=170);plt.close(fig)


def json_value(value):
    if isinstance(value,np.ndarray):return value.tolist()
    if isinstance(value,np.generic):return value.item()
    raise TypeError(type(value).__name__)


def interactive_html(output,rigs,summary):
    cases=[]
    for key in ('alpha1','alpha3','overdrive','loaded_control','failure','free'):
        v=rigs[key]['data']
        # All failure realizations retained; display one, never cherry-pick it.
        cases.append(dict(key=key,label=rigs[key]['label'],alpha=rigs[key]['summary']['alpha'],
            drive=rigs[key]['summary']['drive_pN'],time=np.round(v['time']*1000,6).tolist(),
            center=np.round(v['center']*1000,3).tolist(),p0=np.round(v['p0']*1000,3).tolist(),
            p1=np.round(v['p1']*1000,3).tolist(),fi=np.round(v['fi'],3).tolist(),fj=np.round(v['fj'],3).tolist(),
            alive=v['alive'].astype(int).tolist(),energy=np.round(v['energy']*1000,3).tolist()))
    encoded=json.dumps(dict(cases=cases,summary=summary),default=json_value,separators=(',',':'))
    images=''.join(f'<figure><img src="data:image/png;base64,{base64.b64encode((output/name).read_bytes()).decode()}" alt="{name}"></figure>'
        for name in ('formation_repulsion.png','mechanics_failure.png','lateral_escape.png','weak_boundaries.png'))
    template=HTML_TEMPLATE.replace('__DATA__',encoded).replace('__FIGURES__',images)
    (output/'standalone_contact.html').write_text(template)


HTML_TEMPLATE=r'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Standalone ACTN2 contact experiment</title>
<style>body{margin:0;background:#f5f7fb;color:#203047;font:14px system-ui,sans-serif}header,main{max-width:1250px;margin:auto;padding:20px}h1{font-size:25px;margin:0 0 8px}.muted{color:#61718a;line-height:1.6}.card{background:white;border:1px solid #dce4ef;border-radius:10px;padding:16px;margin:14px 0}.controls{display:flex;gap:12px;align-items:center;flex-wrap:wrap}button,select,input{font:inherit}button,select{border:1px solid #bccbdd;border-radius:5px;background:white;padding:7px;cursor:pointer}input[type=range]{accent-color:#277d9b}canvas{width:100%;height:460px;display:block;background:#f8fbff;touch-action:none;border-radius:8px}#time{flex:1;min-width:200px}#readout{font-variant-numeric:tabular-nums;line-height:1.7;padding:10px 0}img{width:100%;height:auto}figure{margin:0}details{line-height:1.7}label{white-space:nowrap}.warning{background:#fff4d8;padding:10px;border-radius:6px}a{color:#246f9d}</style></head>
<body><header><h1>Standalone ACTN2: formation, repulsion, and failure</h1><div class="muted">Controlled experiments for the proposed model. Not a production trajectory or calibrated sarcomere. All failure realizations are available below.</div></header><main>
<div class="card"><div class="controls"><select id="case" aria-label="Experiment"></select><label>Realization <input id="rep" type="number" min="0" max="0" value="0" style="width:65px"></label><button id="play">Play</button><button id="reset">Reset</button><button id="view">Whole / tip</button></div>
<div class="controls" style="margin:12px 0"><input id="time" type="range" min="0" max="200" value="0"><label><input type="checkbox" id="arrows" checked>Contact-force arrows</label><label><input type="checkbox" id="anchors" checked>Elastic anchors</label><label><input type="checkbox" id="dead" checked>Lost connectors</label></div>
<canvas id="scene" aria-label="Rotatable 3D controlled connector loading experiment"></canvas><div id="readout"></div><div class="muted">Drag to rotate, scroll to zoom. Arrows: 2 nm/pN, without normalization. Blue arrows: contact forces on attachment endpoints. Orange arrow: total contact reaction on myosin, drawn at its center. Gray segments: lost links. The four pairs have independent harmonic anchors; body rotations, thermal body motion, births, and network remodeling are omitted in this mechanical rig.</div></div>
<div class="card"><b>What is actually demonstrated?</b><div id="findings" class="muted"></div></div>
__FIGURES__
<div class="card"><details><summary>Model, assumptions, and scope</summary><p>Myosin length 1500 nm; physical radius 50 nm; structural cylinder radius 130 nm; mechanical approach 50 nm; smoothing 5 nm; k=1000 pN/µm; potential slope 30 pN. At alpha=3 the bounds are 30 pN axial and 90 pN radial/net per connector–myosin interaction.</p><p>Energy is the normalized integral along each material connector. Reactions are derivatives of that energy. Each of four distinct pairs contributes; redundant-link budget invariance is tested separately. The exact maturation gate uses the physical end planes, not the extended mechanical approach. Weak birth has a physical-core veto and exp(−positive incremental pair energy / kBT) factor.</p><p>Mechanical loading starts from 25 nm inherited tip penetration. Elastic anchors represent an external scaffold, not a simulated Z-disc network. The fixed-topology cases intentionally isolate mechanics. The failure ensemble includes baseline and contact-specific molecular loss plus contour rupture; no geometric deletion, rebirth, or body thermal noise. The free-relaxation case has no anchors or external drive.</p><p>Weak-site experiments use actual energy-gradient drift along fixed antiparallel actins, diffusion D=0.01 µm²/s, absorbing physical endpoints with Brownian-bridge correction, and an outside-zone hazard. Weak cohorts are held weak to isolate boundaries. The old reflecting case is a counterfactual reference. These experiments test the chosen implementation and consequences; they do not establish molecular calibration, production stability, or a realistic failure threshold.</p><p>Four static figures above contain all realizations/statistical outcomes; this animation displays a user-selected realization starting at index zero. No external libraries or network requests are needed to open this file.</p></details></div>
<script id="data" type="application/json">__DATA__</script><script>
'use strict';const D=JSON.parse(document.getElementById('data').textContent),$=s=>document.getElementById(s);let ci=0,rep=0,frame=0,playing=false,whole=false,yaw=-.28,pitch=.30,zoom=1,last=0;
const canvas=$('scene'),ctx=canvas.getContext('2d');let W=1000,H=460;
const plus=(a,b)=>a.map((v,j)=>v+b[j]),scale=(a,s)=>a.map(v=>v*s),norm=a=>Math.hypot(...a);
function project(p){const shift=whole?[150,0,0]:[750,0,0],q=p.map((v,j)=>v-shift[j]);let a=q[0]*Math.cos(yaw)-q[1]*Math.sin(yaw),b=q[0]*Math.sin(yaw)+q[1]*Math.cos(yaw),c=q[2];const k=(whole?Math.min(W/2100,H/650):Math.min(W/700,H/430))*zoom;return [W/2+k*a,H/2+k*(b*Math.sin(pitch)-c*Math.cos(pitch))];}
function line(points,color,width=2,dashed=false){ctx.beginPath();points.forEach((p,i)=>{const s=project(p);i?ctx.lineTo(...s):ctx.moveTo(...s)});ctx.strokeStyle=color;ctx.lineWidth=width;ctx.setLineDash(dashed?[4,4]:[]);ctx.stroke();ctx.setLineDash([]);}
function dot(p,color,r=4){const [x,y]=project(p);ctx.beginPath();ctx.arc(x,y,r,0,2*Math.PI);ctx.fillStyle=color;ctx.fill();}
function arrow(p,f,color){if(norm(f)<.01)return;const tip=plus(p,scale(f,2)),a=project(p),b=project(tip),dx=b[0]-a[0],dy=b[1]-a[1],l=Math.hypot(dx,dy);if(l<1)return;line([p,tip],color,2);const u=[dx/l,dy/l],head=Math.min(9,.3*l);ctx.beginPath();ctx.moveTo(...b);ctx.lineTo(b[0]-head*u[0]-.4*head*u[1],b[1]-head*u[1]+.4*head*u[0]);ctx.lineTo(b[0]-head*u[0]+.4*head*u[1],b[1]-head*u[1]-.4*head*u[0]);ctx.closePath();ctx.fillStyle=color;ctx.fill();}
function draw(){const v=D.cases[ci];rep=Math.max(0,Math.min(v.center[0].length-1,Number($('rep').value)||0));frame=Math.max(0,Math.min(v.time.length-1,Number($('time').value)||0));ctx.clearRect(0,0,W,H);const c=v.center[frame][rep],p0=v.p0[frame][rep],p1=v.p1[frame][rep],live=v.alive[frame][rep];
 for(let j=0;j<12;j++){const a=2*Math.PI*j/12;line([plus(c,[-750,50*Math.cos(a),50*Math.sin(a)]),plus(c,[750,50*Math.cos(a),50*Math.sin(a)])],'#8795a7',2);}
 for(const x of [-750,750]){let pts=[];for(let j=0;j<=64;j++){const a=2*Math.PI*j/64;pts.push(plus(c,[x,50*Math.cos(a),50*Math.sin(a)]));}line(pts,'#65788e',2);}
 let boundary=[];for(let j=0;j<=64;j++){const a=2*Math.PI*j/64;boundary.push([750,130*Math.cos(a),130*Math.sin(a)]);}line(boundary,'#bccada',1,true);
 let reaction=[0,0,0],count=0;for(let j=0;j<4;j++){if(!live[j]&&!$('dead').checked)continue;const col=live[j]?'#218e93':'#b5bcc8';if(live[j])count++;line([p0[j],p1[j]],col,live[j]?5:2,!live[j]);dot(p0[j],col);dot(p1[j],col);
 if($('anchors').checked&&v.key!=='free'){line([v.p0[0][rep][j],p0[j]],'#acb8ca',1,true);line([v.p1[0][rep][j],p1[j]],'#acb8ca',1,true);dot(v.p0[0][rep][j],'#cbd5e3',2);dot(v.p1[0][rep][j],'#cbd5e3',2);}
 if(live[j]){const a=v.fi[frame][rep][j],b=v.fj[frame][rep][j];reaction=plus(reaction,scale(plus(a,b),-1));if($('arrows').checked){arrow(p0[j],a,'#277cc1');arrow(p1[j],b,'#277cc1');}}}
 if($('arrows').checked)arrow(c,reaction,'#e58a3d');
 const tip=plus(c,[750,0,0]);dot(tip,'#d26058',4);const screen=project(tip);ctx.fillStyle='#b64a43';ctx.font='12px system-ui';ctx.fillText('myosin + tip',screen[0]+7,screen[1]-7);
 ctx.fillStyle='#52677f';ctx.font='13px system-ui';ctx.fillText(`${v.label} · alpha ${v.alpha} · drive ${v.drive} pN`,15,24);ctx.fillText('Contact arrows: 2 nm/pN',15,45);
 $('readout').innerHTML=`<b>t = ${v.time[frame].toFixed(4)} ms</b> · realization ${rep} / ${v.center[0].length-1} · surviving links ${live.reduce((a,b)=>a+b,0)} / 4 · tip past original scaffold plane ${(tip[0]-750).toFixed(2)} nm<br>Contact + spring + anchor energy: ${v.energy[frame][rep].toFixed(2)} pN·nm · myosin contact reaction (${reaction.map(x=>x.toFixed(2)).join(', ')}) pN`;
}
function resize(){const rect=canvas.getBoundingClientRect();W=rect.width;H=rect.height;const dpr=window.devicePixelRatio||1;canvas.width=W*dpr;canvas.height=H*dpr;ctx.setTransform(dpr,0,0,dpr,0,0);draw();}
function experiment(){ci=Number($('case').value);frame=0;rep=0;$('time').value=0;$('rep').value=0;$('rep').max=D.cases[ci].center[0].length-1;$('time').max=D.cases[ci].time.length-1;draw();}
$('case').innerHTML=D.cases.map((v,i)=>`<option value="${i}">${v.label}</option>`).join('');$('case').onchange=experiment;$('rep').onchange=draw;$('time').oninput=draw;['arrows','anchors','dead'].forEach(id=>$(id).onchange=draw);
$('play').onclick=()=>{playing=!playing;$('play').textContent=playing?'Pause':'Play';};$('reset').onclick=()=>{$('time').value=0;draw();};$('view').onclick=()=>{whole=!whole;zoom=1;draw();};
let drag=null;canvas.onpointerdown=e=>{drag=[e.clientX,e.clientY];canvas.setPointerCapture(e.pointerId);};canvas.onpointermove=e=>{if(!drag)return;yaw+=(e.clientX-drag[0])*.006;pitch=Math.max(-1.5,Math.min(1.5,pitch+(e.clientY-drag[1])*.006));drag=[e.clientX,e.clientY];draw();};canvas.onpointerup=()=>drag=null;canvas.onwheel=e=>{e.preventDefault();zoom=Math.max(.25,Math.min(5,zoom*Math.exp(-e.deltaY*.001)));draw();};
function tick(now){if(playing&&now-last>35){last=now;$('time').value=(frame+1)%D.cases[ci].time.length;draw();}requestAnimationFrame(tick);}window.addEventListener('resize',resize);
$('findings').innerHTML=D.summary.findings.map(x=>`<p>${x}</p>`).join('');experiment();resize();requestAnimationFrame(tick);window.standaloneContact={D,project,draw};
</script></main></body></html>'''


def write_report(output,summary):
    rigs=summary['rigs'];checks=summary['checks'];weak=summary['weak_sites']
    grad=checks['gradients'];conv=checks['numerical_convergence']
    lines=[
        '# Standalone ACTN2 formation, repulsion, and failure',
        '',
        'This is an independently runnable implementation of the current proposal, not a production simulator change or a calibrated molecular model.',
        '',
        '## Findings',
        '',
        *['- '+s for s in summary['findings']],
        '',
        '## Implemented rules',
        '',
        '- Exact same-point full-connector maturation gate: radius 130 nm and physical half-length 750 nm; both weak sites must remain in the 200 nm barbed zone and the supplied registry condition must hold.',
        '- Weak-birth physical-core veto using a 55 nm inflated myosin capsule; incremental maximum-pair contact-energy penalty exp(-positive delta U/kBT). This physical-core reference is an explicit standalone approximation, not a port of all production sterics.',
        '- Smoothed weighted coordinate p=max(0,smin(alpha*(R-rho_epsilon),H-z_epsilon)); H=800 nm, epsilon=5 nm. Quadratic-then-linear Phi with k=1000 pN/um and saturation slope 30 pN. At alpha=3, radial/net force is bounded by 90 pN per connector–myosin interaction, axial force by 30 pN.',
        '- Normalized integral over the entire connector, with endpoint derivatives weighted by (1-t) and t. Equal-and-opposite myosin reaction and exact physical torque from the same energy. Maximum-energy budget over redundant links is tested; distinct pairs add.',
        '- Weak diffusion along fixed antiparallel actins, D=0.01 um^2/s and mobility D/kBT; actual contact-gradient drift, absorbing physical endpoints with local Brownian-bridge correction; soft inner-edge hazard 100/s * sum(([b-200 nm]+/20 nm)^2). Baseline weak off-rate is 10/s. Weak cohorts are held weak in these boundary experiments.',
        '- Failure loading experiment: baseline 2/s plus contact hazard 50/s * [exp(T/5 pN)-1], where T is the larger endpoint contact-force magnitude. Contour rupture at 120 nm. Competing integrated-hazard clocks; no deterministic A-band deletion. These failure parameters are deliberately illustrative, not fitted.',
        '',
        '## Numerical and stochastic checks',
        '',
        f'- Maximum endpoint finite-difference force error: {grad["endpoint_max_error_pN"]:.3g} pN; myosin translation error: {grad["myosin_max_error_pN"]:.3g} pN.',
        f'- Maximum orientation-gradient error: {grad["orientation_max_error_pNum"]:.3g} pN um; total contact torque residual: {grad["torque_balance_max_pNum"]:.3g} pN um.',
        '- Formation agrees with an independent 3,001-point segment reference on 1,200 random segments, plus explicit adversarial fixtures.',
        '- Point-force bounds checked on 10,000 random points at both alpha values. Duplicate-link invariance, distinct-pair additivity and incremental redundant birth energy pass.',
        f'- Closed-loop work, including maximum-energy winner changes, at 128 / 512 / 2048 samples per side: {list(checks["closed_loop_work_pNnm"].values())} pN nm.',
        f'- Timestep/quadrature checks: {json.dumps(conv)}.',
        '- Absorbing endpoint survival checked against the analytic first-passage law with 20,000 walkers. Weak-birth Bernoulli acceptance and fixed-contact-load lifetimes checked against their analytic probabilities (20,000 realizations per condition).',
        '',
        '## Controlled mechanical rigs',
        '',
        'Four distinct radial connectors have endpoints initially at 80 and 120 nm radius, on the same X plane, in four azimuthal directions. Myosin initially overshoots by 25 nm. Each endpoint translates with mobility 70 um/(pN s), myosin with 19 um/(pN s); connector extension stiffness 1000 pN/um and rest length 40 nm.',
        '',
        'The supported cases attach each endpoint to an explicit harmonic external scaffold with stiffness 1000 pN/um. Those springs provide the missing surrounding-network support: a stopped myosin in this rig is not evidence that an unconstrained production Z-disc would hold. Body rotations and body Brownian noise are omitted. The free-relaxation case has no anchors or drive.',
        '',
        'Fixed-topology cases isolate contact mechanics. The stochastic ensemble includes molecular death but no new births, remodeling or replacement. All realizations are included in statistics; the interactive view starts with realization zero and permits selecting every realization.',
        '',
        'Penetration is measured relative to the original scaffold plane X=750 nm, not a moving connector. `connector_axial_overlap` separately records the current tip position minus each connector\'s nearest axial endpoint; radial overlap and the live mask must also be considered.',
        '',
        '| Case | Drive (pN) | Alpha | Final scaffold penetration, median (nm) | Any link death | All links dead |',
        '|---|---:|---:|---:|---:|---:|',
    ]
    for key in ('alpha1','alpha3','overdrive','loaded_control','failure','free'):
        s=rigs[key]
        lines.append(f'| {key} | {s["drive_pN"]:g} | {s["alpha"]:g} | {np.median(s["final_penetration_nm"]):.2f} | {s["any_death_fraction"]:.3f} | {s["all_dead_fraction"]:.3f} |')
    lines += ['', '## Weak-site boundary outcomes', '',
        '| Initial condition | Rule | Alive at 10 ms | Endpoint loss | Soft-edge loss | Baseline loss |',
        '|---|---|---:|---:|---:|---:|']
    for s in weak:
        lines.append(f'| {s["kind"]} | {"reflecting reference" if s["reflecting_reference"] else "absorbing + soft edge"} | {s["alive_fraction"]:.3f} | {s["endpoint_loss"]} | {s["soft_edge_loss"]} | {s["baseline_loss"]} |')
    lines += ['', '## Interpretation and limits', '',
        'The unweighted geometry fails the intended near-tip direction test. Alpha=3 repairs point-level axial capture over a finite range, but an outward finite connector can still carry a larger radial than axial net force. Collective cancellation of lateral reactions on a myosin does not eliminate the sideways loading on each actin.',
        '',
        'Energy-weighted weak association is effectively a hard veto through much of the 50 nm approach zone at these illustrative energies. The separate maturation gate can allow a location where weak-birth acceptance is extremely small. This is a substantive consequence to calibrate, not an implementation error to hide.',
        '',
        'Failure tests verify the specified stochastic law and its interaction with this supported toy rig. They do not show that the chosen off-rate matches ACTN2, or that the production model will stop myosins, preserve connectivity, or share loads correctly.',
        '',
        'This script is nonperiodic and uses one myosin. It does not reproduce production motor chemistry, catch-bond lifetimes, rigid-rod rotations, the full force-cap policy, multiple myosin images, site-registry kinetics, arm-state turnover, or network reformation. Those still require C++ kernel/integrator tests and the staged production comparisons in the implementation plan.',
        '',
        '## Files and reproduction', '',
        '- `standalone_actn2_contact.py`: copied standalone source, with no simulator imports.',
        '- `summary.json`: all numerical checks, parameter values, seeds, timings, fixture outputs and experiment summaries.',
        '- `trajectories.npz`: mechanical positions, contact endpoint forces, live/dead masks, local death codes and energies; weak-site histories/final cohorts; statistical clamp curves.',
        '- `formation_repulsion.png`, `mechanics_failure.png`, `lateral_escape.png`, `weak_boundaries.png`: scientific figures.',
        '- `standalone_contact.html`: self-contained rotatable 3D trajectory viewer and all figures; no WebGL, external dependency or internet connection needed.',
        '',
        'From a compute allocation with numpy, scipy and matplotlib:', '',
        '```bash',
        'python standalone_actn2_contact.py --output /path/to/new_empty_results --replicates 64 --seed 20261028',
        '```',
        '',
        'The output directory must be empty; the script refuses to overwrite existing results. Read current Sherlock agent instructions before scheduling work. A one-CPU run aggregates the deterministic, statistical, stochastic-ensemble and rendering work.',
        '',
        f'Wall time: {summary["wall_seconds"]:.2f} s. Python {summary["python"]}, NumPy {summary["numpy"]}. Seed {summary["seed"]}. No production source or trajectories were modified.',
    ]
    (output/'report.md').write_text('\n'.join(lines)+'\n')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--seed',type=int,default=20261028)
    parser.add_argument('--replicates',type=int,default=64)
    args=parser.parse_args()
    if args.replicates<8:parser.error('Use at least 8 mechanical failure realizations.')
    if args.output.exists() and any(args.output.iterdir()):
        raise FileExistsError(f'Refusing to overwrite nonempty output: {args.output}')
    args.output.mkdir(parents=True,exist_ok=True)
    started=time.perf_counter();model=Model();rng=np.random.default_rng(args.seed)
    print('Checking geometry, gradients, force bounds, and pair budgets',flush=True)
    fixtures=geometry_fixtures(model);checks=check_model(model,rng)
    print('Checking endpoint first-passage, association sampling, and failure statistics',flush=True)
    stats,stochastic=statistical_checks(model,rng);checks.update(stochastic)
    rigs={}
    specifications=[('alpha1','alpha 1 · 20 pN drive',replace(model,alpha=1.),20.,.0006,False,1000.,1),
        ('alpha3','alpha 3 · 20 pN drive',model,20.,.0006,False,1000.,1),
        ('overdrive','alpha 3 · 150 pN overdrive',model,150.,.0006,False,1000.,1),
        ('loaded_control','alpha 3 · 40 pN · fixed topology',model,40.,.002,False,1000.,1),
        ('failure','alpha 3 · 40 pN · stochastic failure',model,40.,.002,True,1000.,args.replicates),
        ('free','alpha 3 · free relaxation',model,0.,.0002,False,0.,1)]
    for key,label,mm,drive,duration,failure,anchor,nrep in specifications:
        begin=time.perf_counter();print(f'Mechanical rig: {key} ({nrep} realization(s))',flush=True)
        data,s=simulate_rig(mm,drive,duration=duration,replicates=nrep,failure=failure,anchor=anchor,seed=args.seed)
        s['wall_seconds']=time.perf_counter()-begin
        rigs[key]=dict(label=label,data=data,summary=s)
        print(json.dumps({key:s}),flush=True)
    print('Checking timestep and quadrature refinement',flush=True)
    coarse=rigs['alpha3']['data']
    finer,_=simulate_rig(model,20.,dt=5.e-8,nquad=32)
    quad,_=simulate_rig(model,20.,dt=1.e-7,nquad=32)
    dt_error=float(np.max(np.abs(coarse['center']-finer['center']))*1000)
    quad_error=float(np.max(np.abs(coarse['center']-quad['center']))*1000)
    free=rigs['free']['data']
    assert np.max(np.diff(free['energy'][:,0]))<1e-9
    assert dt_error<.5 and quad_error<.5,(dt_error,quad_error)
    checks['numerical_convergence']=dict(myosin_center_dt_plus_quadrature_change_nm=dt_error,
        myosin_center_quadrature_only_change_nm=quad_error,
        free_relaxation_energy_monotone=True,dt_s=1e-7,refined_dt_s=5e-8,
        fixed_Gauss_nodes=16,refined_nodes=32,
        trajectory_scope='maximum center difference over all 201 common saved times',
        scope='deterministic supported 20 pN rig; does not certify stochastic production integration')
    print('Weak-site boundary experiments',flush=True)
    weak={}
    for kind in ('barbed','inner'):
        for reflecting in (True,False):
            key=(kind,reflecting);weak[key]=weak_sites(model,kind,reflecting,seed=args.seed)
            print(json.dumps(weak[key]['summary']),flush=True)
    small={k:v['summary'] for k,v in rigs.items()}
    a1=small['alpha1']['final_penetration_nm'][0];a3=small['alpha3']['final_penetration_nm'][0]
    high=small['overdrive']['final_penetration_nm'][0]
    behavior=dict(weighted_supported_rig_recovers=a3<0,unweighted_rig_does_not_recover=a1>25,
        above_capacity_drive_passes=high>25,
        barbed_absorption_releases=weak[('barbed',False)]['summary']['endpoint_loss']>0,
        soft_edge_does_not_reflect=weak[('inner',False)]['summary']['surviving_outside_zone']>0,
        molecular_contact_deaths=small['failure']['death_counts']['2']>0)
    assert all(behavior.values()),behavior
    checks['controlled_behavior']=behavior
    findings=[
        f'Formation: full-connector maturation rejects a connector whose midpoint lies outside 130 nm; radial and axial conditions at different points do not cause a false rejection. Weak birth and maturation use distinct rules.',
        f'With an explicit elastic scaffold and 20 pN drive, alpha=1 finishes at {a1:.1f} nm penetration; alpha=3 finishes at {a3:.1f} nm. Positive means the physical tip is past the original scaffold plane. This result is conditional on the scaffold.',
        f'Alpha=3 does not eliminate sideways loading: at 20 nm burial the 40 nm outward connector carries approximately 19.9 pN axial and 30.3 pN radial. A point at its inner endpoint is almost purely axial.',
        f'A 150 pN drive exceeds the four-pair axial capacity (at most 120 pN) and passes, ending at {high:.1f} nm penetration. This is a bounded soft resistance, not an impenetrable wall.',
        f'In {args.replicates} coupled 40 pN failure realizations, {100*small["failure"]["any_death_fraction"]:.1f}% lose at least one link and {100*small["failure"]["all_dead_fraction"]:.1f}% lose all four by 2 ms. No connector is killed merely for entering the buffered region.',
        f'{100*small["failure"]["final_overshoot_fraction"]:.1f}% of the failure ensemble pass the original scaffold plane; {small["failure"]["passed_with_live_connectors"]} realizations do so with live connectors. Asymmetric loss permits lateral myosin escape, shown in lateral_escape.png.',
        f'Weak sites detach at the barbed endpoint or cross the 200 nm structural edge and incur a higher hazard; the reflecting reference instead confines surviving sites. Failure and edge parameters are illustrative.'
    ]
    summary=dict(model=asdict(model),seed=args.seed,python=platform.python_version(),numpy=np.__version__,
        scipy=scipy.__version__,matplotlib=matplotlib.__version__,slurm_job_id=os.environ.get('SLURM_JOB_ID'),
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        fixtures=fixtures,checks=checks,rigs=small,weak_sites=[v['summary'] for v in weak.values()],
        findings=findings,production_model_changed=False,calibrated=False)
    print('Rendering figures and self-contained interactive viewer',flush=True)
    save_figures(args.output,model,rigs,weak,stats)
    summary['wall_seconds']=time.perf_counter()-started
    interactive_html(args.output,rigs,summary)
    arrays={f'{key}_{name}':value for key,rig in rigs.items() for name,value in rig['data'].items()}
    arrays.update({f'weak_{kind}_{"reflecting" if reflecting else "absorbing"}_{name}':value
        for (kind,reflecting),run in weak.items() for name,value in run.items() if name!='summary'})
    arrays.update({f'stat_{key}_{name}':value for key,run in stats.items() for name,value in run.items()})
    np.savez_compressed(args.output/'trajectories.npz',**arrays)
    (args.output/'summary.json').write_text(json.dumps(summary,default=json_value,indent=2))
    write_report(args.output,summary)
    import shutil
    shutil.copy2(__file__,args.output/'standalone_actn2_contact.py')
    print(json.dumps(dict(output=str(args.output),wall_seconds=summary['wall_seconds'],findings=findings)),flush=True)


if __name__=='__main__':main()
