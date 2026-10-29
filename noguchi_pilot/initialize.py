"""Documented approximate Model-A realization; NOT original Noguchi initial conditions."""
import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.special import roots_legendre
import json
from pathlib import Path
out=Path(__file__).parent
rng=np.random.default_rng(2001)
# These reconstruction choices are explicit; published settings are kept separately.
N=[30000,10000,19985]; M=[.6,.2,.2]; eps=[.02,.04,.004]
h=.25; z0=.02; Re=.05; ah=.2
r=np.geomspace(1e-5,4,900)
# Exact Abel deprojection of a Sersic n=4 surface profile, normalization below.
x,w=roots_legendre(160); u=(x+1)*6; wu=w*6
R=r[:,None]*np.cosh(u)
b=7.6692494425
sp=-b/(4*Re)*(R/Re)**(-.75)*np.exp(-b*((R/Re)**.25-1))
rhob=-np.sum(sp*wu,axis=1)/np.pi
# Hernquist halo, smoothly tapered between 0.8 and 1; assumed, not published.
taper=np.where(r<.8,1,np.where(r<1,.5*(1+np.cos(np.pi*(r-.8)/.2)),0))
rhoh=taper/(r*(r+ah)**3)
rhos=[]; cumulative=[]; dms=[]
for rho,m in [(rhoh,.2),(rhob,.2)]:
 c=cumulative_trapezoid(4*np.pi*r*r*rho,r,initial=0);rho*=m/c[-1]
 c=cumulative_trapezoid(4*np.pi*r*r*rho,r,initial=0)
 rhos.append(rho);cumulative.append(c)
 dms.append(4*np.pi*r*r*rho*np.gradient(r))
# Spherical average of finite exponential disc, used only for spheroid Jeans moments.
md=.6*(1-(1+r/h)*np.exp(-r/h))/(1-(1+1/h)*np.exp(-1/h));md[r>=1]=.6
# Spheroid Jeans moments use softened spherical angular integrals for each component,
# with the disc approximated by its spherical monopole. No fake pre-relaxation.
dmd=np.gradient(md,r)*np.gradient(r); dmd[r>=1]=0
def spherical_phi(target_eps):
 p=np.zeros_like(r)
 for dm,ep in zip([dmd,*dms],[.02,.04,.004]):
  e2=(target_eps**2+ep**2)/2
  a=r[:,None];s=r[None,:]
  kernel=2/(np.sqrt((a+s)**2+e2)+np.sqrt((a-s)**2+e2))
  p-=kernel@dm
 return p
# Axisymmetric finite disc force in midplane, softened and thickness-approximated.
rd=np.geomspace(.0002,1.15,450)
x,w=roots_legendre(220); sr=(x+1)/2; sw=w/2
sig_norm=.6/(2*np.pi*h*h*(1-(1+1/h)*np.exp(-1/h)))
Sigma=lambda R:sig_norm*np.exp(-R/h)
ph=(np.arange(256)+.5)*2*np.pi/256
pdisc=np.zeros_like(rd)
for k,Rt in enumerate(rd):
 dist=Rt*Rt+sr[:,None]**2-2*Rt*sr[:,None]*np.cos(ph)+eps[0]**2+z0*z0
 pdisc[k]=-np.sum(2*np.pi*sr*Sigma(sr)*sw*np.mean(1/np.sqrt(dist),axis=1))
psph=np.zeros_like(rd)
for dm,ep in zip(dms,[.04,.004]):
 a=rd[:,None];s=r[None,:];e2=(.02**2+ep**2)/2
 psph-= (2/(np.sqrt((a+s)**2+e2)+np.sqrt((a-s)**2+e2)))@dm
vc2=rd*np.gradient(pdisc+psph,rd)
omega2=vc2/rd**2
kappa2=rd*np.gradient(omega2,rd)+4*omega2
if np.any(kappa2<=0):raise ValueError('unstable circular orbit in initialization grid')
sr2=(3.36*Sigma(rd)*1.5/np.sqrt(kappa2))**2
sp2=sr2*kappa2/(4*omega2)
dln=np.gradient(np.log(Sigma(rd)*sr2),np.log(rd))
vmean2=vc2+sr2*(1+dln)-sp2
# Finite-thickness Jeans approximation; external vertical harmonic restoring force.
nu2=np.maximum(np.gradient(psph,rd)/rd,0)
sz2=np.pi*Sigma(rd)*z0+nu2*z0*z0
cdf=(1-(1+rd/h)*np.exp(-rd/h)); cdf/=1-(1+1/h)*np.exp(-1/h)
rr=np.interp(rng.random(N[0]),cdf,rd); az=rng.uniform(0,2*np.pi,N[0]);zz=z0*np.arctanh(rng.uniform(-1,1,N[0]))
vR=rng.normal(size=N[0])*np.sqrt(np.interp(rr,rd,sr2))
vphi=np.sqrt(np.maximum(np.interp(rr,rd,vmean2),0))+rng.normal(size=N[0])*np.sqrt(np.interp(rr,rd,sp2))
vz=rng.normal(size=N[0])*np.sqrt(np.interp(rr,rd,sz2))
positions=[np.column_stack([rr*np.cos(az),rr*np.sin(az),zz])]
velocities=[np.column_stack([vR*np.cos(az)-vphi*np.sin(az),vR*np.sin(az)+vphi*np.cos(az),vz])]
for n,rho,c,ep in zip(N[1:],rhos,cumulative,eps[1:]):
 rr=np.interp(rng.random(n)*c[-1],c,r)
 dirs=rng.normal(size=(n,3));dirs/=np.linalg.norm(dirs,axis=1)[:,None]
 potential=spherical_phi(ep);g=np.gradient(potential,r)
 integ=cumulative_trapezoid((rho*g)[::-1],r[::-1],initial=0)[::-1]*-1
 s2=np.divide(integ,rho,out=np.zeros_like(r),where=rho>0)
 vel=rng.normal(size=(n,3))*np.sqrt(np.maximum(np.interp(rr,r,s2),0))[:,None]
 positions.append(rr[:,None]*dirs);velocities.append(vel)
p=np.vstack(positions);v=np.vstack(velocities)
m=np.repeat(np.array(M)/N,N);e=np.repeat(eps,N);comp=np.repeat([0,1,2],N)
p-=np.sum(m[:,None]*p,axis=0);v-=np.sum(m[:,None]*v,axis=0)
data=np.column_stack([p,v,m,e,comp]).astype('float64')
data.tofile(out/'initial.bin')
np.savez(out/'initial.npz',data=data)
np.savez(out/'target_profiles.npz',r=rd,vc=np.sqrt(vc2),sigmaR=np.sqrt(sr2),vmean=np.sqrt(np.maximum(vmean2,0)),Sigma=Sigma(rd))
meta={'label':'approximate equilibrium pilot, not exact reproduction','N':N,'masses':M,'softening':eps,'h':h,'Re':Re,'assumed_z0':z0,'assumed_halo_a':ah,'seed':2001,'negative_streaming_grid_points':int(np.sum(vmean2[rd<1]<0)),'initialization':'Jeans moments; no staged relaxation; spherical monopole approximation for spheroids; local vertical approximation for disc; Gaussian velocities; no escape rejection'}
(out/'initialization.json').write_text(json.dumps(meta,indent=2));print(meta)
