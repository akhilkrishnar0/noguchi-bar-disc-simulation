"""Low-N experimental galaxy seed, not a reproduction of Noguchi Model A."""
import argparse, json
from pathlib import Path
import numpy as np
p=argparse.ArgumentParser();p.add_argument('--n-disc',type=int,default=750);p.add_argument('--n-halo',type=int,default=300);p.add_argument('--n-bulge',type=int,default=450);p.add_argument('--output',default='lowram/seed.bin');args=p.parse_args()
rng=np.random.default_rng(2001)
N=[args.n_disc,args.n_halo,args.n_bulge]; assert min(N)>0
M=[.6,.2,.2]; soft=[.05,.08,.04]
# Disc exponential with smooth outer taper, approximate Gaussian thickness.
def disc_r(n):
    out=[]
    while len(out)<n:
        r=rng.gamma(2,.25,size=n)
        w=np.where(r<.85,1.,np.where(r<1.15,.5*(1+np.cos(np.pi*np.minimum(r-.85,.30)/.30)),0.))
        out.extend(r[rng.random(n)<w].tolist())
    return np.array(out[:n])
rd=disc_r(N[0]);a=rng.uniform(0,2*np.pi,N[0]);z=rng.normal(0,.06,N[0]);pos=[np.column_stack([rd*np.cos(a),rd*np.sin(a),z])]
# Initial circular speed in rough spherical monopole, corrected during optional validation.
# This is NOT an axisymmetric Jeans distribution function.
md=.6*(1-(1+rd/.25)*np.exp(-rd/.25));mh=.2*rd**2/(rd+.2)**2;mb=.2*rd**2/(rd+.05/1.8153)**2
vphi=np.sqrt((md+mh+mb)*rd**2/(rd**2+.04**2)**1.5)
vR=rng.normal(0,.085,N[0]);vp=vphi+rng.normal(0,.06,N[0]);vz=rng.normal(0,.06,N[0]);vel=[np.column_stack([vR*np.cos(a)-vp*np.sin(a),vR*np.sin(a)+vp*np.cos(a),vz])]
# Hernquist inverse CDF, finite radius cut: assumed spheroid profiles.
for n,scale in [(N[1],.2),(N[2],.05/1.8153)]:
    u=rng.uniform(0,.92,n);r=scale*np.sqrt(u)/(1-np.sqrt(u));dirs=rng.normal(size=(n,3));dirs/=np.linalg.norm(dirs,axis=1)[:,None]
    pos.append(dirs*r[:,None]);vdisp=np.sqrt(.45*(.6*(1-(1+r/.25)*np.exp(-r/.25))+.2*r*r/(r+.2)**2+.2*r*r/(r+.05/1.8153)**2)/np.sqrt(r*r+.04**2))
    vel.append(rng.normal(size=(n,3))*vdisp[:,None])
xyz=np.concatenate(pos);v=np.concatenate(vel);mass=np.concatenate([np.full(n,m/n) for n,m in zip(N,M)]);eps=np.concatenate([np.full(n,s) for n,s in zip(N,soft)]);comp=np.repeat(np.arange(3),N)
xyz-=np.average(xyz,axis=0,weights=mass);v-=np.average(v,axis=0,weights=mass)
out=np.column_stack([xyz,v,mass,eps,comp]).astype('<f8');Path(args.output).parent.mkdir(parents=True,exist_ok=True);out.tofile(args.output)
Path(args.output).with_suffix('.json').write_text(json.dumps({'N':N,'mass_fractions':M,'softening':soft,'seed':2001,'warning':'approximate Jeans/circular velocities; relaxed spheroids and vertical calibration necessary but not sufficient for stationarity'},indent=2))
print('Created',args.output,'particles',sum(N),'bytes',out.nbytes)
