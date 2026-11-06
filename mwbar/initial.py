"""Approximate stationary epicycle-ensemble disc; not a full action DF."""
import numpy as np
from .model import Galaxy,KMS_TO_KPCMYR

def sample_tracers(model:Galaxy,n:int,seed:int=2026):
    rng=np.random.default_rng(seed)
    # Guiding-centre distribution with exponential surface density.
    rg=rng.gamma(shape=2.,scale=model.tracer_rd,size=n)
    invalid=(rg<.65)|(rg>17.)
    while invalid.any():
        rg[invalid]=rng.gamma(shape=2.,scale=model.tracer_rd,size=invalid.sum())
        invalid=(rg<.65)|(rg>17.)
    omega,kappa,nu=model.frequencies(rg)
    theta_r=rng.uniform(0,2*np.pi,size=n)
    sigmaR=np.clip(model.sigma_R_sun_kms*np.exp((model.solar_radius-rg)/12.),18,85)*KMS_TO_KPCMYR
    # Epicycle amplitudes (Rayleigh), truncated to avoid radius inversion.
    ar=np.minimum(rng.rayleigh(sigmaR/np.maximum(kappa,1e-9)),.45*rg)
    R=rg+ar*np.cos(theta_r)
    vr=-kappa*ar*np.sin(theta_r)
    Lz=rg*model.circular_speed(rg)
    vphi=Lz/R
    angle=rng.uniform(0,2*np.pi,size=n)
    c,s=np.cos(angle),np.sin(angle)
    x,y=R*c,R*s
    # Harmonic vertical action approximation at each guiding centre.
    theta_z=rng.uniform(0,2*np.pi,size=n)
    az=rng.rayleigh(model.tracer_hz,size=n)
    z=az*np.cos(theta_z)
    vz=-nu*az*np.sin(theta_z)
    pos=np.column_stack((x,y,z))
    vel=np.column_stack((vr*c-vphi*s,vr*s+vphi*c,vz))
    if not (np.all(np.isfinite(pos)) and np.all(np.isfinite(vel))):
        raise RuntimeError('Invalid sampled initial positions or velocities')
    return pos,vel
