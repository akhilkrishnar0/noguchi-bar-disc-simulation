"""Guiding-centre, radial profile, and numerical validation measurements."""
import numpy as np
from .model import Galaxy

def guiding_radius(model:Galaxy,x,v):
    Lz=x[:,0]*v[:,1]-x[:,1]*v[:,0]
    grid=np.linspace(.03,60.,12000)
    Lcirc=grid*model.circular_speed(grid)
    rg=np.interp(np.maximum(Lz,0),Lcirc,grid,left=grid[0],right=grid[-1])
    return rg

def profile(model:Galaxy,x,edge=None):
    if edge is None:edge=np.arange(.5,16.01,.5)
    R=np.hypot(x[:,0],x[:,1])
    n,_=np.histogram(R,bins=edge)
    area=np.pi*(edge[1:]**2-edge[:-1]**2)
    sigma=(n*model.disc_mass/len(x))/(area*1.e6) # Msun/pc² weighted tracers
    center=(edge[1:]+edge[:-1])*.5
    return center,sigma,n

def scalelength(model,x,fit_min=4.,fit_max=11.):
    r,sigma,n=profile(model,x)
    fit=(r>=fit_min)&(r<=fit_max)&(n>=8)
    if fit.sum()<4:return np.nan
    slope,_=np.polyfit(r[fit],np.log(sigma[fit]),1,w=np.sqrt(n[fit]))
    return -1/slope if slope<0 else np.nan

def diagnostics(model:Galaxy,run,barred):
    n=len(run['t']); arr=[]
    x0,v0=run['x'][0].astype(float),run['v'][0].astype(float)
    rg0=guiding_radius(model,x0,v0)
    for i in range(n):
        t=run['t'][i];x=run['x'][i].astype(float);v=run['v'][i].astype(float)
        _,phi=model.force(x,t,barred)
        E=.5*np.sum(v*v,axis=1)+phi
        Lz=x[:,0]*v[:,1]-x[:,1]*v[:,0]
        J=E-model.omega_bar*Lz
        R=np.hypot(x[:,0],x[:,1]);az=np.arctan2(x[:,1],x[:,0])
        region=(R>.5)&(R<6.)
        a2=np.abs(np.mean(np.exp(2j*az[region])))
        rg=guiding_radius(model,x,v)
        dRg=rg-rg0
        arr.append((t,np.median(R),np.sqrt(np.mean(x[:,2]**2)),scalelength(model,x),
                    a2,np.mean(E),np.mean(Lz),np.mean(J),
                    np.median(np.abs(dRg)),np.mean(np.abs(dRg)>.5),
                    np.mean(np.abs(dRg)>1.),np.mean(R),np.mean(rg)))
    header='time_myr,R50_kpc,zrms_kpc,Rd_kpc,A2_inner,Emean_kpc2Myr2,Lzmean_kpc2Myr,Jacobi_mean,median_abs_dRg_kpc,frac_abs_dRg_gt_0p5,frac_abs_dRg_gt_1,Rmean_kpc,Rguiding_mean_kpc'
    data=np.asarray(arr)
    a=data[:,5] if not barred else data[:,7]
    change=np.abs((a-a[0])/np.maximum(np.abs(a[0]),1e-12))
    report={
       'mode':'existing rotating external bar' if barred else 'axisymmetric control',
       'R50_change_pct':float(100*(data[-1,1]/data[0,1]-1)),
       'zrms_change_pct':float(100*(data[-1,2]/data[0,2]-1)),
       'Rd_initial_kpc':float(data[0,3]),
       'Rd_final_kpc':float(data[-1,3]),
       'A2_initial':float(data[0,4]),
       'A2_final':float(data[-1,4]),
       'mean_Lz_change_pct':float(100*(data[-1,6]/data[0,6]-1)),
       'integral_name':'Jacobi E-Omega_b Lz' if barred else 'inertial E',
       'max_relative_mean_integral_deviation':float(change.max()),
       'median_abs_guiding_shift_kpc':float(data[-1,8]),
       'frac_abs_guiding_shift_gt_0p5_kpc':float(data[-1,9]),
       'frac_abs_guiding_shift_gt_1_kpc':float(data[-1,10])
    }
    return header,data,report
