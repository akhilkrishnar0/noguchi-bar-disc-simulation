"""Inertial Cartesian kick-drift-kick integration of collisionless tracers."""
from dataclasses import dataclass
import numpy as np
from .model import Galaxy
from .initial import sample_tracers

@dataclass(frozen=True)
class Settings:
    n:int=10000
    seed:int=2026
    timestep_myr:float=.25
    warmup_myr:float=500.
    growth_myr:float=300.
    science_myr:float=300.
    snapshot_myr:float=5.

    def validate(self):
        if self.n<1 or self.timestep_myr<=0 or self.science_myr<=0:return False
        if self.warmup_myr<self.growth_myr or self.growth_myr<=0:return False
        for period in (self.warmup_myr,self.growth_myr,self.science_myr,self.snapshot_myr):
            if abs(round(period/self.timestep_myr)*self.timestep_myr-period)>1e-8:return False
        return True

def grown_bar_amplitude(t,s:Settings):
    if t<=-s.warmup_myr:return 0.
    if t>=-s.warmup_myr+s.growth_myr:return 1.
    u=(t+s.warmup_myr)/s.growth_myr
    return 6*u**5-15*u**4+10*u**3

def leapfrog(pos,vel,t,dt,model:Galaxy,barred:bool,amplitude1=1.,amplitude2=1.):
    a,_=model.force(pos,t,barred,amplitude=amplitude1)
    half=vel+.5*dt*a
    pos_new=pos+dt*half
    a1,_=model.force(pos_new,t+dt,barred,amplitude=amplitude2)
    vel_new=half+.5*dt*a1
    return pos_new,vel_new

def integrate_pair(model:Galaxy=Galaxy(),settings:Settings=Settings(),progress=None):
    if not settings.validate():raise ValueError('Times must be multiples of the integration step and growth<=warmup')
    start_pos,start_vel=sample_tracers(model,settings.n,settings.seed)
    outputs={}
    for barred in (False,True):
        pos=start_pos.copy();vel=start_vel.copy()
        dt=settings.timestep_myr
        n_pre=round(settings.warmup_myr/dt)
        # Prior to science t=0, the barred case ramps smoothly, then remains
        # fully established for (warmup-growth) Myr. Control stays axisymmetric.
        for j in range(n_pre):
            t=-settings.warmup_myr+j*dt
            a0=grown_bar_amplitude(t,settings) if barred else 0.
            a1=grown_bar_amplitude(t+dt,settings) if barred else 0.
            pos,vel=leapfrog(pos,vel,t,dt,model,barred,a0,a1)
        prep_pos,prep_vel=pos.copy(),vel.copy()
        cadence=round(settings.snapshot_myr/dt)
        steps=round(settings.science_myr/dt)
        samples_pos=[pos.astype('float32').copy()]
        samples_vel=[vel.astype('float32').copy()]
        times=[0.]
        for j in range(1,steps+1):
            t=(j-1)*dt
            pos,vel=leapfrog(pos,vel,t,dt,model,barred)
            if j%cadence==0:
                times.append(j*dt)
                samples_pos.append(pos.astype('float32').copy())
                samples_vel.append(vel.astype('float32').copy())
        outputs['barred' if barred else 'unbarred']={
            't':np.asarray(times),
            'x':np.asarray(samples_pos),
            'v':np.asarray(samples_vel),
            'phase0_x':prep_pos,'phase0_v':prep_vel}
        if progress:progress('barred' if barred else 'unbarred',len(times))
    return outputs
