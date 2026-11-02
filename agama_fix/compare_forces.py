#!/usr/bin/env python3
"""Direct-pair force vs AGAMA reference for identical 10k particle ICs.

Finite-N force fluctuations are distinct from systematic softening bias;
report component-specific and radial/vertical diagnostics. No calibration of
particle velocities or artificial relaxation is performed here.
"""
import json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parent
OUT=ROOT/"results"
p=np.fromfile(OUT/"initial.bin",dtype="<f8").reshape(-1,9)
ref=np.load(OUT/"agama_force_reference.npz")
ids=ref["indices"];orig=ref["forces"];xyz=p[:,:3];m=p[:,6]
z=p[ids,2];vz=p[ids,5];sample=p[ids]
assert len(ids)==len(orig) and np.allclose(sample[:,:3],ref["positions"])
assert np.all(np.isfinite(orig))
disk=sample[:,8]==0
def pair_force(scale):
    out=np.zeros((len(ids),3),dtype=np.float64)
    eps=p[:,7]*scale
    for st in range(0,len(ids),16):
        target=ids[st:st+16]
        delta=xyz[None,:,:]-xyz[target,None,:]
        dist2=np.einsum('ijk,ijk->ij',delta,delta)
        dist2+=.5*(eps[target,None]**2+eps[None,:]**2)
        dist2[np.arange(len(target)),target]=np.inf
        weights=m[None,:]*dist2**(-1.5)
        out[st:st+len(target)]=np.einsum("ij,ijk->ik",weights,delta,optimize=True)
    return out
def summary(fc, label):
    norm=np.linalg.norm(orig,axis=1)
    rel=np.linalg.norm(fc-orig,axis=1)/np.maximum(norm,1e-10)
    zd=z[disk];vm=vz[disk]
    Kz=float(np.sum(sample[disk,6]*(vm-np.average(vm,weights=sample[disk,6]))**2))
    def virial(force):
        Wz=float(-np.sum(sample[disk,6]*zd*force[disk,2]))
        return Kz/Wz if Wz>0 else None
    # Midplane force disagreement dominated by finite-N shot noise;
    # use median over disc separately from spheroids.
    dx=sample[disk,:2];rr=np.maximum(np.linalg.norm(dx,axis=1),1e-8)
    ar_model=np.sum(orig[disk,:2]*dx,axis=1)/rr
    ar_disc=np.sum(fc[disk,:2]*dx,axis=1)/rr
    return dict(tag=label,median_force_error_pct=float(np.median(rel)*100),
      disc_force_median_error_pct=float(np.median(rel[disk])*100),
      disc_radial_force_median_fractional_bias=float(np.median((ar_disc-ar_model)/np.maximum(np.abs(ar_model),1e-10))),
      force_error_95pct=float(np.percentile(rel,95)*100),
      sampled_disc_vertical_virial_ratio=virial(fc),
      sampled_AGAMA_smooth_vertical_virial_ratio=virial(orig))
results=[]
for fac in (1.0,.5,.25):
    f=pair_force(fac)
    tag=f"eps_{str(fac).replace('.','p')}"
    results.append(summary(f,tag))
    new=p.copy();new[:,7]*=fac
    new.astype("<f8").tofile(OUT/f"initial_{tag}.bin")
    print(results[-1],flush=True)
(OUT/"force_comparison.json").write_text(json.dumps(dict(
     setup="same AGAMA particle realization; exact finite-N pairwise Plummer force vs smooth unsoftened AGAMA potential",
     N=int(len(p)),sample_count=int(len(ids)),cases=results,
     interpretation="No single force-error statistic establishes an equilibrated live galaxy. Shot noise, softening and dynamics require separate tests."),indent=2))
