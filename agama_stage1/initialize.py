#!/usr/bin/env python3
"""AGAMA action-DF self-consistent 3-component initial model (G=M=L=1).

A new dynamically motivated reconstruction, NOT Noguchi's exact Model A.
No fabricated equilibrium: the live 300 Myr evolution must pass the tests.
Based on Agama's official example_self_consistent_model_bdh.py.
"""
import json, argparse, hashlib, sys
from pathlib import Path
import numpy as np
import agama

def build(out:Path, nd:int=6000, nh:int=2000, nb:int=2000, iterations:int=4):
    out.mkdir(parents=True,exist_ok=True)
    rd=.25;hd=.025; re=.05; rb=re/1.8153
    # Hernquist-like bulge/halo are *assumptions*, not Noguchi's halo law.
    diskpar=dict(type="Disk",surfaceDensity=.6/(2*np.pi*rd**2),scaleRadius=rd,scaleHeight=-hd)
    bulgepar=dict(type="Spheroid",mass=.2,scaleRadius=rb,gamma=1,beta=4,outerCutoffRadius=3)
    halopar=dict(type="Spheroid",mass=.2,scaleRadius=.2,gamma=1,beta=4,outerCutoffRadius=4)
    dd,db,dh=[agama.Density(**p) for p in (diskpar,bulgepar,halopar)]
    sc=agama.SelfConsistentModel(rminSph=.001,rmaxSph=30,sizeRadialSph=40,
        lmaxAngularSph=4,RminCyl=.003,RmaxCyl=5,sizeRadialCyl=28,
        zminCyl=.003,zmaxCyl=3,sizeVerticalCyl=24)
    sc.components=[agama.Component(density=dd,disklike=True),
                   agama.Component(density=db,disklike=False),
                   agama.Component(density=dh,disklike=False)]
    print("AGAMA: initial static density guess",flush=True)
    sc.iterate()
    init_potential=sc.potential
    # Set Q minimum to 1.5 using the *unsoftened* potential and approximate
    # exponential radial dispersion. Recalculate the realized Q later.
    R=np.geomspace(.065,.70,180);xyz=np.column_stack((R,0*R,0*R))
    force,der=init_potential.eval(xyz,acc=True,der=True)
    kappa=np.sqrt(np.maximum(-der[:,0]-3*force[:,0]/R,1.e-12))
    rsigma=2*rd
    sigma0_target=diskpar["surfaceDensity"]
    baseq=np.exp(-R/rsigma)*kappa/(3.36*sigma0_target*np.exp(-R/rd)*agama.G)
    sig0=1.5/np.min(baseq)
    diskDFpar=dict(type="QuasiIsothermal",potential=init_potential,
                   Sigma0=sigma0_target,Rdisk=rd*.96,Hdisk=hd,
                   sigmar0=sig0,Rsigmar=rsigma,sigmamin=.01*sig0)
    dfDisk0=agama.DistributionFunction(**diskDFpar)
    diskDFpar["Sigma0"]*=.6/dfDisk0.totalMass()
    dfDisk=agama.DistributionFunction(**diskDFpar)
    dfBulge=agama.DistributionFunction(type="QuasiSpherical",potential=init_potential,density=db)
    dfHalo=agama.DistributionFunction(type="QuasiSpherical",potential=init_potential,density=dh)
    sc.components=[
        agama.Component(df=dfDisk,disklike=True,RminCyl=.006,RmaxCyl=2.5,sizeRadialCyl=20,
                        zminCyl=.005,zmaxCyl=.5,sizeVerticalCyl=16),
        agama.Component(df=dfBulge,disklike=False,rminSph=.001,rmaxSph=2,
                        sizeRadialSph=26,lmaxAngularSph=4),
        agama.Component(df=dfHalo,disklike=False,rminSph=.004,rmaxSph=15,
                        sizeRadialSph=30,lmaxAngularSph=4)]
    masses=[]
    for k in range(iterations):
        print(f"AGAMA: self-consistency iteration {k+1}/{iterations}",flush=True)
        sc.iterate()
        masses.append([float(c.density.totalMass()) for c in sc.components])
        print("  component masses:",masses[-1],flush=True)
    # Density representation/softening must be checked: AGAMA's equilibrium
    # potential is smooth, whereas the live C++ integrator uses softened forces.
    particles=[]
    for comp,(df,n,eps) in enumerate([(dfDisk,nd,.02),(dfHalo,nh,.04),(dfBulge,nb,.004)]):
        print(f"AGAMA: DF-sampling component {comp}: N={n}",flush=True)
        xv,m=agama.GalaxyModel(potential=sc.potential,df=df,af=sc.af).sample(n)
        xv=np.asarray(xv,dtype=np.float64)
        m=np.asarray(m,dtype=np.float64).reshape(-1)
        if xv.shape!=(n,6) or m.shape!=(n,) or not np.all(np.isfinite(xv)):
            raise RuntimeError(f"Invalid AGAMA sample comp={comp}: {xv.shape} {m.shape}")
        particles.append(np.column_stack((xv,m,np.full(n,eps),np.full(n,comp))))
    data=np.vstack(particles)
    # Preserve phase-space DF while centring global COM; no arbitrary velocity renormalization.
    weight=data[:,6];data[:,:3]-=np.average(data[:,:3],axis=0,weights=weight)
    data[:,3:6]-=np.average(data[:,3:6],axis=0,weights=weight)
    if not np.all(np.isfinite(data)) or np.min(data[:,6])<=0:raise RuntimeError("Bad particle masses/positions")
    data.astype("<f8").tofile(out/"initial.bin")
    sha=hashlib.sha256((out/"initial.bin").read_bytes()).hexdigest()
    # Smooth-potential circular and epicyclic-frequency references for Q diagnostics.
    R=np.geomspace(.04,1.5,300);pts=np.column_stack((R,R*0,R*0))
    force,der=sc.potential.eval(pts,acc=True,der=True)
    vc=np.sqrt(np.maximum(-R*force[:,0],0))
    kappa=np.sqrt(np.maximum(-der[:,0]-3*force[:,0]/R,0))
    np.savez(out/"target.npz",R=R,vc=vc,kappa=kappa)
    meta=dict(status="candidate_IC_not_yet_live_validated",method="AGAMA action-DF SCM",
              N=[nd,nh,nb],softenings=[.02,.04,.004],
              realized_component_mass=[float(data[data[:,8]==i,6].sum()) for i in (0,1,2)],
              df_mass_iteration_history=masses,disk_scale_radius=rd,disk_scale_height=hd,
              requested_Qmin=1.5,approx_sigmar0=float(sig0),
              binary_sha256=sha,agama_version=str(getattr(agama,"__version__","unknown")),
              caveats=["DF uses smooth potential; N-body uses Plummer softening",
                       "small N and finite live halo/bulge may evolve",
                       "new assumed halo profile; not exact Noguchi (2001) Model A"])
    (out/"initial_conditions.json").write_text(json.dumps(meta,indent=2))
    print(json.dumps(meta,indent=2),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out",default="agama_stage1/results")
    p.add_argument("--nd",type=int,default=6000)
    p.add_argument("--nh",type=int,default=2000)
    p.add_argument("--nb",type=int,default=2000)
    p.add_argument("--iterations",type=int,default=4)
    a=p.parse_args();build(Path(a.out),a.nd,a.nh,a.nb,a.iterations)
