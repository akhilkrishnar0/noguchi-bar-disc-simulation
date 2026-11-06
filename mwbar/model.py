"""An analytic Milky Way-like potential; positions in kpc, time in Myr.

All tracer particles have zero self-gravity.  These potentials are an
observationally motivated approximation, NOT a fully fitted MW mass model.
The rotating m=2 bar has constant pattern speed in the science interval.
"""
from dataclasses import dataclass
import numpy as np

KMS_TO_KPCMYR=1.0227121650537077e-3
G=4.30091e-6*KMS_TO_KPCMYR**2  # kpc^3 / (Msun Myr^2)

@dataclass(frozen=True)
class Galaxy:
    solar_radius:float=8.2
    solar_vcirc:float=232.8
    disc_mass:float=5.3e10
    disc_a:float=2.75
    disc_b:float=.30
    bulge_mass:float=.8e10
    bulge_a:float=.55
    halo_rs:float=18.0
    halo_c:float=12.0
    bar_length:float=4.8
    bar_core:float=3.5
    bar_axis_ratio_z:float=.60
    bar_pattern_kms_kpc:float=39.0
    bar_force_fraction:float=.12
    bar_angle_deg:float=0.0
    tracer_rd:float=2.6
    tracer_hz:float=.26
    sigma_R_sun_kms:float=35.0

    @property
    def omega_bar(self):
        return self.bar_pattern_kms_kpc*KMS_TO_KPCMYR
    @property
    def halo_mass(self):
        """Set the NFW normalization to match the fiducial solar rotation speed."""
        r=self.solar_radius
        vd2=G*self.disc_mass*r*r/(r*r+(self.disc_a+self.disc_b)**2)**1.5
        vb2=G*self.bulge_mass*r/(r+self.bulge_a)**2
        h=self.halo_rs
        fx=np.log1p(r/h)-(r/h)/(1+r/h)
        fc=np.log1p(self.halo_c)-self.halo_c/(1+self.halo_c)
        target=(self.solar_vcirc*KMS_TO_KPCMYR)**2
        remaining=target-vd2-vb2
        if remaining<=0:raise ValueError('Disc+bulge exceed requested solar circular speed')
        return remaining*r*fc/(G*fx)

    def axisymmetric(self,xyz):
        """Returns acceleration[kpc/Myr²], potential[kpc²/Myr²]."""
        p=np.asarray(xyz,dtype=float)
        if p.ndim!=2 or p.shape[1]!=3:raise ValueError('xyz must be (N,3)')
        x,y,z=p.T
        R2=x*x+y*y
        rr=np.maximum(np.linalg.norm(p,axis=1),1e-9)
        # Miyamoto-Nagai stellar disc
        b=np.sqrt(z*z+self.disc_b**2)
        q=self.disc_a+b
        D2=R2+q*q
        D=np.sqrt(D2)
        ff=-G*self.disc_mass/(D2*D)
        acc=ff[:,None]*p
        acc[:,2]*=q/b
        potential=-G*self.disc_mass/D
        # Hernquist-like round bulge
        fac=-G*self.bulge_mass/(rr*(rr+self.bulge_a)**2)
        acc+=fac[:,None]*p
        potential+=-G*self.bulge_mass/(rr+self.bulge_a)
        # NFW halo, with mass inside rvir=c*rs normalized to halo_mass
        rs=self.halo_rs
        u=rr/rs
        mass_shape=np.log1p(u)-u/(1+u)
        fc=np.log1p(self.halo_c)-self.halo_c/(1+self.halo_c)
        coeff=G*self.halo_mass/fc
        acc+=-(coeff*mass_shape/rr**3)[:,None]*p
        potential+=-coeff*np.log1p(u)/rr
        return acc,potential

    def circular_frequency(self,radius):
        r=np.asarray(radius,dtype=float)
        v=self.circular_speed(r)
        return v/np.maximum(r,1e-9)

    def circular_speed(self,radius):
        r=np.atleast_1d(radius).astype(float)
        pts=np.column_stack((r,np.zeros_like(r),np.zeros_like(r)))
        a,_=self.axisymmetric(pts)
        return np.sqrt(np.maximum(-r*a[:,0],0))

    def frequencies(self,radius):
        r=np.atleast_1d(radius).astype(float)
        d=np.maximum(.002,.004*r)
        om=self.circular_frequency(r)
        # kappa² = r d(Omega²)/dr + 4 Omega²
        oplus=self.circular_frequency(r+d)**2
        ominus=self.circular_frequency(np.maximum(r-d,1e-4))**2
        k2=r*(oplus-ominus)/(2*d)+4*om**2
        zp=.003
        pts=np.column_stack((r,np.zeros_like(r),np.full_like(r,zp)))
        a,_=self.axisymmetric(pts)
        n2=-a[:,2]/zp
        return om,np.sqrt(np.maximum(k2,0)),np.sqrt(np.maximum(n2,0))

    @property
    def bar_amplitude(self):
        # Calibrate max F_tangential/F_radial over 1..6 kpc at z=0.
        radii=np.linspace(1.,6.,501)
        a,_=self.axisymmetric(np.column_stack((radii,0*radii,0*radii)))
        ft_unit=2*radii/(radii**2+self.bar_core**2)**2.5
        strength=np.max(ft_unit/np.maximum(np.abs(a[:,0]),1e-20))
        return self.bar_force_fraction/strength

    def bar(self,xyz,t,amplitude=1.0):
        """Rotating quadrupole, analytic Cartesian 3D forces.

        Phi=-A*(X²-Y²) / (X²+Y²+(z/q)²+b²)^(5/2)
        with X,Y rotated with respect to the fixed inertial frame.
        The bar component is externally prescribed, not a live Ferrers bar.
        """
        xyz=np.asarray(xyz,dtype=float)
        if amplitude==0:
            return np.zeros_like(xyz),np.zeros(len(xyz))
        a=self.bar_amplitude*amplitude
        phi=self.omega_bar*t+np.deg2rad(self.bar_angle_deg)
        c,s=np.cos(phi),np.sin(phi)
        x,y,z=xyz.T
        X=c*x+s*y
        Y=-s*x+c*y
        Q=X*X-Y*Y
        S=X*X+Y*Y+(z/self.bar_axis_ratio_z)**2+self.bar_core**2
        common=a*S**(-3.5)
        fx=common*X*(2*S-5*Q)
        fy=common*Y*(-2*S-5*Q)
        fz=-common*5*Q*z/self.bar_axis_ratio_z**2
        f=np.column_stack((c*fx-s*fy,s*fx+c*fy,fz))
        pot=-a*Q*S**(-2.5)
        return f,pot

    def force(self,xyz,t,barred=False,amplitude=1.0):
        a,phi=self.axisymmetric(xyz)
        if barred:
            b,p=self.bar(xyz,t,amplitude=amplitude)
            a+=b;phi+=p
        return a,phi

    def resonance_radii(self):
        r=np.linspace(.35,24.,2800)
        om,kappa,_=self.frequencies(r)
        out={}
        for name,f in [('ILR',om-kappa/2),('CR',om),('OLR',om+kappa/2)]:
            v=f-self.omega_bar
            roots=np.flatnonzero(v[:-1]*v[1:]<=0)
            out[name]=[float(np.interp(0,[v[i],v[i+1]][::(-1 if v[i]>v[i+1] else 1)], [r[i],r[i+1]][::(-1 if v[i]>v[i+1] else 1)])) for i in roots]
        return out
