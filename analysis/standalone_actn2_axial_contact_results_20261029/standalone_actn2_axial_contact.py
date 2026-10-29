#!/usr/bin/env python3
"""Standalone ACTN2 axial-stop contact: myosin pushed back only along its axis.

Replaces the weighted-energy contact of standalone_actn2_contact.py with a law
whose force on every connector point is parallel to the myosin axis. Radial
position only gates whether contact applies; it never produces a force.

    U(x)   = g(rho) Phi(depth(z)),   depth = max(0, H - sqrt(z^2 + eps^2))
    f(x)   = -(dU/dz) u              (connector point; myosin gets -f)
    dropped: -Phi g'(rho) rho_hat    (the radial gradient of U)

f is exactly the axial part of -grad U, so at fixed radial positions it is
conservative and U is the work to push the connector out along the axis.
It is not the full gradient of any energy: closed loops that cross the gate
boundary do net work Phi*(g_in - g_out). Checks below quantify that.

Reuses the audited rig integrator and hazard law from standalone_actn2_contact
unchanged by swapping its contact kernel. Units: um, pN, seconds. Illustrative
constants only. Run in a Slurm allocation.
"""
from __future__ import annotations
import argparse
from dataclasses import asdict, replace
import importlib.util
import sys
import json
from pathlib import Path
import platform
import time

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

here=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('weighted',here/'standalone_actn2_contact.py')
weighted=importlib.util.module_from_spec(spec);sys.modules['weighted']=weighted;spec.loader.exec_module(weighted)
Model,gauss=weighted.Model,weighted.gauss
weighted_contact=weighted.contact
GATE_WIDTH=.010  # gate falls 1 -> 0 between R and R + 10 nm


def gate(rho,model,width=GATE_WIDTH):
    s=np.clip((rho-model.radius)/width,0.,1.)
    return 1-s*s*(3-2*s),-6*s*(1-s)/width


def point_axial(x,model=Model(),center=(0.,0.,0.),axis=(1.,0.,0.)):
    """Pseudo-energy U, axial force, and the dropped radial gradient term."""
    x=np.asarray(x,dtype=float)-np.asarray(center)
    u=np.asarray(axis,dtype=float);u=u/np.linalg.norm(u)
    z=np.sum(x*u,axis=-1)
    r=x-z[...,None]*u
    rho=np.sqrt(np.sum(r*r,axis=-1))
    za=np.sqrt(z*z+model.epsilon**2)
    depth=np.maximum(model.half+model.approach-za,0.)
    ps=model.slope/model.k
    phi=np.where(depth<=ps,.5*model.k*depth**2,model.slope*(depth-.5*ps))
    g,dg=gate(rho,model)
    force=(g*np.minimum(model.k*depth,model.slope)*z/za)[...,None]*u
    rhat=np.divide(r,rho[...,None],out=np.zeros_like(r),where=rho[...,None]>0)
    dropped=-(phi*dg)[...,None]*rhat
    return g*phi,force,dropped


def axial_contact(p0,p1,model=Model(),center=(0.,0.,0.),axis=(1.,0.,0.),n=32):
    """Same signature and quadrature as the weighted kernel."""
    p0,p1=np.broadcast_arrays(np.asarray(p0,dtype=float),np.asarray(p1,dtype=float))
    c=np.broadcast_to(np.asarray(center),p0.shape)
    ts,ws=gauss(n)
    x=p0[...,None,:]+ts[:,None]*(p1-p0)[...,None,:]
    e,f,_=point_axial(x,model,c[...,None,:],axis)
    fi=np.sum(f*(ws*(1-ts))[:,None],axis=-2)
    fj=np.sum(f*(ws*ts)[:,None],axis=-2)
    torque=-np.sum(np.cross(x-c[...,None,:],f)*ws[:,None],axis=-2)
    return np.sum(e*ws,axis=-1),fi,fj,torque


def check_axial(model,rng):
    checks={};step=1.e-7
    p0=rng.uniform([.65,.065,-.02],[.83,.145,.02],(150,3))
    p1=p0+rng.uniform([-.015,.015,-.01],[.015,.065,.01],(150,3))
    center=np.array([.008,-.003,.004]);axis=np.array([1,.03,-.02]);axis/=np.linalg.norm(axis)
    e,fi,fj,tm=axial_contact(p0,p1,model,center,axis,n=64)
    # 1. No sideways component on either endpoint or on the myosin.
    side=lambda f:np.linalg.norm(f-np.sum(f*axis,axis=-1)[...,None]*axis,axis=-1)
    lateral=float(max(side(fi).max(),side(fj).max(),side(fi+fj).max()))
    assert lateral<1e-12
    # 2. Axial component equals -dU/dz for a rigid axial shift of the connector.
    du=axis*step
    fd=-(axial_contact(p0+du,p1+du,model,center,axis,64)[0]-axial_contact(p0-du,p1-du,model,center,axis,64)[0])/(2*step)
    axial_error=float(np.max(np.abs(fd-np.sum((fi+fj)*axis,axis=-1))))
    assert axial_error<1e-3
    # 3. Angular momentum: endpoint torques plus myosin reaction torque cancel.
    balance=np.cross(p0-center,fi)+np.cross(p1-center,fj)+tm
    assert np.max(np.abs(balance))<1e-10
    xx=rng.uniform([-.85,-.18,-.18],[.85,.18,.18],(10000,3))
    _,f,_=point_axial(xx,model)
    assert np.max(np.abs(f[:,1:]))==0 and np.max(np.abs(f[:,0]))<=model.slope+1e-10
    checks['axial_law']=dict(max_lateral_force_pN=lateral,axial_fd_error_pN=axial_error,
        torque_balance_max_pNum=float(np.max(np.abs(balance))),
        point_bound='|f| <= 30 pN, f_y = f_z = 0 on 10,000 random points')
    # 4. Non-conservation: loop work of one point around a (z, rho) rectangle.
    loops={}
    for label,(za,zb,ra,rb) in dict(inside_gate=(.770,.790,.080,.120),
            crossing_gate=(.770,.790,.100,.150)).items():
        corners=np.array([[za,ra,0],[zb,ra,0],[zb,rb,0],[za,rb,0],[za,ra,0]])
        out={}
        for count in (512,4096):
            path=np.concatenate([a+(b-a)*np.arange(count)[:,None]/count for a,b in zip(corners[:-1],corners[1:])])
            dp=np.roll(path,-1,axis=0)-path
            out[str(count)]=float(np.sum(point_axial(path+dp/2,model)[1]*dp)*1000)
        phi=lambda z:point_axial(np.array([z,0.,0.]),model)[0]
        g_in,g_out=gate(np.array([ra,rb]),model)[0]
        out['analytic']=float(1000*(phi(za)-phi(zb))*(g_in-g_out))
        loops[label]=out
    assert abs(loops['inside_gate']['4096'])<1e-6
    checks['loop_work_pNnm']=loops
    # 5. The probe that exposed the weighted law: 40 nm outward connector, 20 nm burial.
    probes=[]
    for name,fn,mm in (('weighted alpha=3',weighted_contact,model),('axial stop',axial_contact,model)):
        for ell in (0.,.040,.080):
            _,fa,fb,_=fn([.730,.080,0],[.730,.080+ell,0],mm,n=256);ff=fa+fb
            probes.append(dict(law=name,length_nm=ell*1000,axial_pN=float(ff[0]),
                radial_pN=float(np.linalg.norm(ff[1:]))))
    checks['force_direction_probes']=probes
    return checks


def run_rig(law,drive,removed=(),**kw):
    """Audited integrator; removed connectors are parked at 400-440 nm radius,
    outside both laws' range, so they exert no contact force."""
    original_init,original_contact=weighted.rig_initial,weighted.contact
    def init(replicates):
        p0,p1,c=original_init(replicates)
        for j in removed:
            p0[:,j,1:]*=5.;p1[:,j,1:]*=440/120
        return p0,p1,c
    weighted.rig_initial=init
    weighted.contact=axial_contact if law=='axial' else original_contact
    try:return weighted.simulate_rig(Model(),drive,**kw)
    finally:weighted.rig_initial,weighted.contact=original_init,original_contact


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--seed',type=int,default=20261028)
    parser.add_argument('--replicates',type=int,default=64)
    args=parser.parse_args()
    if args.output.exists() and any(args.output.iterdir()):
        raise FileExistsError(f'Refusing to overwrite nonempty output: {args.output}')
    args.output.mkdir(parents=True,exist_ok=True)
    started=time.perf_counter();model=Model();rng=np.random.default_rng(args.seed)
    checks=check_axial(model,rng);print(json.dumps(checks),flush=True)
    cases={};data={}
    for law in ('weighted','axial'):
        for removed in ((),(0,),(0,1),(0,1,2)):
            key=f'{law}_remove{len(removed)}'
            d,s=run_rig(law,40.,removed,duration=.002,seed=args.seed)
            data[key]=d;cases[key]=s
        d,s=run_rig(law,40.,duration=.002,failure=True,replicates=args.replicates,seed=args.seed)
        data[f'{law}_failure']=d;cases[f'{law}_failure']=s
        d,s=run_rig(law,20.,duration=.0006,seed=args.seed);cases[f'{law}_20pN']=s
        if law=='axial':
            fine,_=run_rig(law,20.,duration=.0006,dt=5.e-8,nquad=32,seed=args.seed)
            checks['axial_dt_quadrature_center_change_nm']=float(np.max(np.abs(d['center']-fine['center']))*1000)
            checks['axial_20pN_work_balance_residual_pNnm']=s['work_balance_residual_pNnm']
    rows=[]
    for key,s in cases.items():
        lat=np.asarray(s['final_myosin_lateral_shift_nm']);pen=np.asarray(s['final_penetration_nm'])
        rows.append(dict(case=key,drive_pN=s['drive_pN'],realizations=s['replicates'],
            median_penetration_nm=float(np.median(pen)),median_lateral_nm=float(np.median(lat)),
            max_lateral_nm=float(lat.max()),passed_fraction=s['final_overshoot_fraction'],
            any_death=s['any_death_fraction'],all_dead=s['all_dead_fraction'],
            passed_with_live=s['passed_with_live_connectors'],
            max_outer_radial_shift_nm=float(np.max(s['final_max_radial_shift_nm']))))
        print(json.dumps(rows[-1]),flush=True)
    # Realization 1 replays the saved sideways-escape instance with identical clocks.
    escape={law:dict(final_lateral_nm=cases[f'{law}_failure']['final_myosin_lateral_shift_nm'][1],
        final_penetration_nm=cases[f'{law}_failure']['final_penetration_nm'][1],
        final_alive=cases[f'{law}_failure']['final_alive_links'][1]) for law in ('weighted','axial')}
    assert all(r['max_lateral_nm']<1e-9 for r in rows if r['case'].startswith('axial'))
    figures(args.output,data,checks)
    np.savez_compressed(args.output/'trajectories.npz',**{f'{k}__{f}':v for k,d in data.items() for f,v in d.items()})
    summary=dict(model=asdict(model),gate_width_um=GATE_WIDTH,seed=args.seed,checks=checks,
        cases=rows,escape_instance_1=escape,python=platform.python_version(),numpy=np.__version__,
        wall_seconds=time.perf_counter()-started)
    (args.output/'summary.json').write_text(json.dumps(summary,indent=1))
    print(json.dumps(dict(escape_instance_1=escape,wall_seconds=summary['wall_seconds'])),flush=True)


def figures(output,data,checks):
    colors={'weighted':'#d77531','axial':'#2588ac'}
    fig,axs=plt.subplots(2,3,figsize=(15,8.5),constrained_layout=True)
    probes=[p for p in checks['force_direction_probes'] if p['length_nm']==40]
    ax=axs[0,0];xs=np.arange(2)
    ax.bar(xs-.18,[p['axial_pN'] for p in probes],.36,label='axial',color='#2588ac')
    ax.bar(xs+.18,[p['radial_pN'] for p in probes],.36,label='sideways (radial)',color='#d45b6c')
    ax.set_xticks(xs,[p['law'] for p in probes]);ax.legend()
    ax.set(title='A · 40 nm outward connector, 20 nm burial\nnet force on connector',ylabel='Force (pN)')
    for col,metric,ylabel in ((1,'pen','Tip past original plane (nm)'),(2,'lat','Myosin sideways displacement (nm)')):
        ax=axs[0,col]
        for law in ('weighted','axial'):
            for k,ls in zip(range(4),('-','--','-.',':')):
                d=data[f'{law}_remove{k}'];t=d['time']*1000
                y=d['penetration'][:,0]*1000 if metric=='pen' else np.linalg.norm(d['center'][:,0,1:],axis=-1)*1000
                ax.plot(t,y,ls,color=colors[law],label=f'{law}, {4-k} connectors')
        ax.set(xlabel='Time (ms)',ylabel=ylabel,
            title=('B · Axial stop vs connectors present\n40 pN, fixed topology' if metric=='pen' else
                   'C · Sideways motion from contact\nasymmetric connector removal'))
        if metric=='pen':ax.set_ylim(-100,600);ax.legend(fontsize=7,ncol=2)
    for col,law in ((0,'weighted'),(1,'axial')):
        d=data[f'{law}_failure'];t=d['time']*1000;lat=np.linalg.norm(d['center'][...,1:],axis=-1)*1000
        ax=axs[1,col]
        for i in range(lat.shape[1]):ax.plot(d['penetration'][:,i]*1000,lat[:,i],color=colors[law],alpha=.25,lw=.8)
        ax.set(xlabel='Tip past original plane (nm)',ylabel='Myosin sideways displacement (nm)',
            title=f'{"D" if col==0 else "E"} · {law} law, 40 pN stochastic failure\nall {lat.shape[1]} realizations',
            xlim=(-100,2000),ylim=(-5,max(150,lat.max()*1.05)))
    ax=axs[1,2]
    for law in ('weighted','axial'):
        d=data[f'{law}_failure'];t=d['time']*1000
        ax.plot(t,np.mean(d['penetration']>0,axis=1),color=colors[law],label=f'{law}: passed plane')
        ax.plot(t,np.mean(d['alive'].sum(-1)<4,axis=1),'--',color=colors[law],label=f'{law}: lost >= 1')
    ax.set(xlabel='Time (ms)',ylabel='Fraction of realizations',title='F · Failure ensemble outcomes');ax.legend(fontsize=8)
    fig.suptitle('Axial-stop contact: radial position gates contact, force acts only along the myosin axis',fontsize=12)
    fig.savefig(output/'axial_vs_weighted.png',dpi=160);plt.close(fig)


if __name__=='__main__':
    main()
