"""Gigayear-scale *tracer-response* observables in a prescribed barred potential.

Equal-mass tracers are assigned a fixed total disc mass and constant r-band M/L.
Intrinsic bar gravity is fixed in the model; derived A2/response extent are not
self-gravitating bar properties or changes in the bar's imposed length.
"""
import numpy as np
from mwbar.measure import profile,scalelength

M_R_SUN_AB=4.65
M_TO_L_R=1.5  # assumed constant Msun/Lsun, r-band
CONVERSION=21.572
FOURIER_EDGES=np.arange(0.5,8.01,0.5)
PHASE_TOL_DEG=20.0


def mu_r_from_sigma(sigma,ml_r=M_TO_L_R):
    """Face-on mu_r proxy [AB mag/arcsec^2] from surface mass density [Msun/pc^2]."""
    sig=np.asarray(sigma,dtype=float)
    if ml_r<=0:raise ValueError('M/L must be positive')
    with np.errstate(divide='ignore',invalid='ignore'):
        return np.where(sig>0,M_R_SUN_AB+CONVERSION-2.5*np.log10(sig/ml_r),np.nan)


def fourier_radial(x,model,t,edges=FOURIER_EDGES):
    """Binned m=2 amplitude/phase relative to the imposed rotating potential."""
    R=np.hypot(x[:,0],x[:,1]);az=np.arctan2(x[:,1],x[:,0]);bangle=model.omega_bar*t+np.deg2rad(model.bar_angle_deg)
    count,_=np.histogram(R,bins=edges)
    real,_=np.histogram(R,bins=edges,weights=np.cos(2*az))
    imag,_=np.histogram(R,bins=edges,weights=np.sin(2*az))
    Z=(real+1j*imag)/np.maximum(count,1)
    amp=np.abs(Z)
    phase_offset_deg=np.rad2deg(.5*np.angle(Z*np.exp(-2j*bangle)))
    centers=.5*(edges[:-1]+edges[1:])
    return centers,amp,phase_offset_deg,count


def response_extent(x,model,t,min_count=70,amp_floor=0.08,phase_tolerance_deg=20.):
    """Fourier m=2 coherent-response extent, not an intrinsic potential bar length.

    Find the inner peak of |A2(R)| in 0.75<R<3.75 kpc, then move outward
    until A2 drops below max(0.08,0.5 * local peak) or the m=2 phase departs
    by >20deg from the imposed bar. No detected coherent region -> NaN.
    Used only as a *noisy diagnostic*, with threshold and fit uncertainty.
    """
    R,A2,delta,n=fourier_radial(x,model,t)
    zone=(R>=.75)&(R<=3.75)&(n>=min_count)
    if not np.any(zone):return np.nan,np.nan,np.nan
    candidates=np.flatnonzero(zone)
    peak=int(candidates[np.argmax(A2[candidates])]);peakamp=float(A2[peak])
    valid=(n>=min_count)&(A2>=max(amp_floor,.5*peakamp))&(np.abs(delta)<=phase_tolerance_deg)
    if peakamp<.10 or not valid[peak]:return np.nan,peakamp,float(R[peak])
    last=peak
    while last+1<len(R) and valid[last+1]:last+=1
    # Demand at least 2 accepted radial bins; otherwise a chance peak is not a coherent bar.
    first=peak
    while first-1>=0 and valid[first-1]:first-=1
    if last-first+1<2:return np.nan,peakamp,float(R[peak])
    return float(FOURIER_EDGES[last+1]),peakamp,float(R[peak])


def outer_fit(model,x,fit_min=4.,fit_max=11.):
    """Weighted log Sigma fit with unconstrained zero point; returns Rd, mu0 proxy."""
    r,sigma,n=profile(model,x)
    good=(r>=fit_min)&(r<=fit_max)&(n>=8)&(sigma>0)
    if np.count_nonzero(good)<4:return np.nan,np.nan
    slope,zero=np.polyfit(r[good],np.log(sigma[good]),1,w=np.sqrt(n[good]))
    if slope>=0:return np.nan,np.nan
    return float(-1./slope),float(mu_r_from_sigma(np.exp(zero)))


def diagnostics_v3(model,runs,ml_r=M_TO_L_R):
    """Compute per-frame bar response, brightness and scalelength with explicit columns."""
    if ml_r != M_TO_L_R:raise ValueError('Version 3 fixed M/L in comparison; explicitly reconfigure to vary it')
    cols=('time_Myr','A2_inner','A2_peak_inner','response_extent_kpc','response_A2_peak_R_kpc',
          'bar_axis_phase_offset_deg_at_peak','Rd_4_11_kpc','Rd_6_12_kpc',
          'mu0_r_extrapolated_mag_arcsec2','R50_kpc','disc_zrms_kpc')
    out={}; profiles={}
    for case,run in runs.items():
        lines=[];prof=[]
        for t,x in zip(run['t'],run['x']):
            r,amp,phase,n=fourier_radial(x,model,t)
            extent,amax,rmax=response_extent(x,model,t)
            good=np.flatnonzero(np.isclose(r,rmax)) if np.isfinite(rmax) else []
            delta=float(phase[good[0]]) if len(good) else np.nan
            radius=np.hypot(x[:,0],x[:,1]);phi=np.arctan2(x[:,1],x[:,0]);inner=(radius>.5)&(radius<6)
            a2_inner=float(np.abs(np.mean(np.exp(2j*phi[inner]))))
            rd,mu0=outer_fit(model,x);rd2=outer_fit(model,x,6,12)[0]
            lines.append([t,a2_inner,amax,extent,rmax,delta,rd,rd2,mu0,
                          float(np.median(radius)),float(np.sqrt(np.mean(x[:,2]**2)))])
            pR,sigma,count=profile(model,x)
            for R,si,ci,mu in zip(pR,sigma,count,mu_r_from_sigma(sigma)):
                prof.append([t,R,si,mu,ci])
        out[case]=(cols,np.asarray(lines,float))
        profiles[case]=(('time_Myr','radius_kpc','tracer_Sigma_Msun_pc2',
                         'mu_r_ABmag_arcsec2','tracer_count'),np.asarray(prof,float))
    return out,profiles
