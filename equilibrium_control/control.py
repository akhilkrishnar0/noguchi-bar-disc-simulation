"""Reproducible numerical control: tracers in a smooth fixed potential.

This is NOT a self-gravitating N-body disc or a reconstruction of Noguchi (2001).
It deliberately isolates orbit integration and external bar forcing.
"""
import argparse,json
from pathlib import Path
import numpy as np
MD,MH,MB=.6,.2,.2
A_DISK,B_DISK,B_HALO,B_BULGE=.25,.08,.30,.12
BAR_STRENGTH,BAR_SCALE,OMEGA,BAR_RAMP=.008,.30,1.5,4.
TIME_MYR=75.

def initial(n=3000,seed=2001):
    rng=np.random.default_rng(seed);r=[]
    while len(r)<n:
        v=rng.gamma(2,.25,size=n)
        r.extend(v[v<1.].tolist())
    R=np.array(r[:n]);phi=rng.uniform(0,2*np.pi,n)
    amp=rng.rayleigh(.020,n);phase=rng.uniform(0,2*np.pi,n)
    z=amp*np.sin(phase)
    vc2=(MD*R**2/(R**2+(A_DISK+B_DISK)**2)**1.5+
         MH*R**2/(R**2+B_HALO**2)**1.5+
         MB*R**2/(R**2+B_BULGE**2)**1.5)
    nu2=(MD*(A_DISK+B_DISK)/(B_DISK*(R**2+(A_DISK+B_DISK)**2)**1.5)+
         MH/(R**2+B_HALO**2)**1.5+MB/(R**2+B_BULGE**2)**1.5)
    vp=np.sqrt(vc2)
    x=np.column_stack((R*np.cos(phi),R*np.sin(phi),z))
    v=np.column_stack((-vp*np.sin(phi),vp*np.cos(phi),np.sqrt(nu2)*amp*np.cos(phase)))
    return x,v

def gravity(x,t,forced=False,with_pot=False):
    xx,yy,z=x.T;R2=xx*xx+yy*yy
    B=np.sqrt(z*z+B_DISK**2);D=np.sqrt(R2+(A_DISK+B)**2)
    h=(R2+z*z+B_HALO**2)**1.5;b=(R2+z*z+B_BULGE**2)**1.5
    f=MD/D**3+MH/h+MB/b
    a=-f[:,None]*x
    a[:,2]+=-MD*z*((A_DISK+B)/B-1)/D**3
    phi=-MD/D-MH/np.sqrt(R2+z*z+B_HALO**2)-MB/np.sqrt(R2+z*z+B_BULGE**2)
    if forced:
        angle=OMEGA*t;c=np.cos(angle);s=np.sin(angle)
        X=c*xx+s*yy;Y=-s*xx+c*yy;Q=X*X-Y*Y;den=R2+BAR_SCALE**2
        ramp=0 if t<=0 else (1. if t>=BAR_RAMP else .5*(1-np.cos(np.pi*t/BAR_RAMP)))
        gx=BAR_STRENGTH*ramp*(2*X*den-3*X*Q)/den**2.5
        gy=BAR_STRENGTH*ramp*(-2*Y*den-3*Y*Q)/den**2.5
        a[:,0]+=c*gx-s*gy;a[:,1]+=s*gx+c*gy
        phi+=-BAR_STRENGTH*ramp*Q/den**1.5
    return (a,phi) if with_pot else a

def diagnostics(x,v,t,forced):
    _,potential=gravity(x,t,forced,with_pot=True)
    E=np.mean(.5*np.sum(v*v,axis=1)+potential)
    Lz=np.mean(x[:,0]*v[:,1]-x[:,1]*v[:,0])
    R=np.hypot(x[:,0],x[:,1]);az=np.arctan2(x[:,1],x[:,0]);sel=(R>.05)&(R<.6)
    m2=np.mean(np.exp(2j*az[sel])) if sel.sum()>30 else np.nan+0j
    return [t*TIME_MYR,E,Lz,E-OMEGA*Lz,float(np.sqrt(np.mean(x[:,2]**2))),float(np.median(R)),float(abs(m2)),int(sel.sum())]

def run(n=3000,seed=2001,step_myr=1.,duration_myr=4500.,output_myr=5.,forced=False,out=None):
    steps=round(duration_myr/step_myr);cadence=round(output_myr/step_myr)
    assert steps>0 and cadence>0 and abs(steps*step_myr-duration_myr)<1e-8
    assert abs(cadence*step_myr-output_myr)<1e-8
    x,v=initial(n,seed);dt=step_myr/TIME_MYR
    snapshots=[x.astype(np.float32).copy()];history=[diagnostics(x,v,0,forced)]
    acc=gravity(x,0,forced)
    for step in range(1,steps+1):
        v+=.5*dt*acc;x+=dt*v
        t=step*dt;acc=gravity(x,t,forced);v+=.5*dt*acc
        if step%cadence==0:
            snapshots.append(x.astype(np.float32).copy())
            history.append(diagnostics(x,v,t,forced))
    data=np.array(history)
    result={'mode':'forced analytic rotating m=2' if forced else 'axisymmetric analytic control',
      'self_gravity':False,'n':n,'seed':seed,'dt_myr':step_myr,'output_myr':output_myr,
      't_end_myr':duration_myr,'rms_height_change_pct':float(100*(data[-1,4]/data[0,4]-1)),
      'R50_change_pct':float(100*(data[-1,5]/data[0,5]-1)),
      'energy_change_pct':float(100*(data[-1,1]/data[0,1]-1)),
      'Lz_change_pct':float(100*(data[-1,2]/data[0,2]-1)),
      'A2_initial_final':[float(data[0,6]),float(data[-1,6])],
      'min_A2_count':int(data[:,7].min())}
    if forced:
        post=data[data[:,0]>=BAR_RAMP*TIME_MYR]
        result['jacobi_post_ramp_max_fraction']=float(np.max(np.abs((post[:,3]-post[0,3])/post[0,3])))
    if out:
        out=Path(out);out.mkdir(parents=True,exist_ok=True)
        np.savez_compressed(out/'evolution.npz',x=np.array(snapshots),history=data)
        np.savetxt(out/'diagnostics.csv',data,delimiter=',',header='time_myr,Emean,Lzmean,Jacobi,zrms,R50,A2,N_A2',comments='')
        (out/'summary.json').write_text(json.dumps(result,indent=2))
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--n',type=int,default=3000)
    p.add_argument('--seed',type=int,default=2001)
    p.add_argument('--dt-myr',type=float,default=1.)
    p.add_argument('--duration-myr',type=float,default=4500.)
    p.add_argument('--output-myr',type=float,default=5.)
    p.add_argument('--out',default='equilibrium_control/results/unforced')
    p.add_argument('--forced',action='store_true')
    args=p.parse_args()
    print(json.dumps(run(args.n,args.seed,args.dt_myr,args.duration_myr,args.output_myr,args.forced,args.out),indent=2))
