#!/usr/bin/env python3
"""Analysis and reject-on-failure gate for the 300 Myr live equilibrium test."""
import json, numpy as np
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

P=Path(__file__).resolve().parent
D=P/"results";S=D/"live"
def load(step):
    v=np.fromfile(S/f"snapshot_{step}.bin",dtype="<f8")
    if v.size%9:raise RuntimeError("Corrupt particle snapshot")
    return v.reshape(-1,9)
def structural(a):
    m=a[:,8]==0;d=a[m];pos=d[:,:3]
    xyz_com=np.average(a[:,:3],weights=a[:,6],axis=0)
    x=pos[:,0]-xyz_com[0];y=pos[:,1]-xyz_com[1];z=pos[:,2]-xyz_com[2]
    r=np.hypot(x,y)
    ring=(r>.05)&(r<.6)
    amp=abs(np.mean(np.exp(2j*np.arctan2(y[ring],x[ring]))))
    return dict(R50=float(np.median(r)),zrms=float(np.sqrt(np.mean(z*z))),A2=float(amp))
def profiles(a,ref):
    d=a[a[:,8]==0];R=np.hypot(d[:,0],d[:,1]);vR=(d[:,0]*d[:,3]+d[:,1]*d[:,4])/np.maximum(R,1e-8)
    edges=np.arange(.075,.826,.05);center=.5*(edges[:-1]+edges[1:])
    values=[]
    for k,r in enumerate(center):
        mask=(R>=edges[k])&(R<edges[k+1]);n=int(mask.sum())
        surf=float(d[mask,6].sum()/(np.pi*(edges[k+1]**2-edges[k]**2)))
        sig=float(np.std(vR[mask],ddof=1)) if n>10 else np.nan
        kappa=float(np.interp(r,ref["R"],ref["kappa"]))
        Q=sig*kappa/(3.36*surf) if surf>0 else np.nan
        values.append([r,n,surf,sig,Q])
    return np.asarray(values)
if __name__=="__main__":
    initial,mid,final=[load(s) for s in (0,600,1200)]
    m=[structural(a) for a in (initial,mid,final)]
    f=np.genfromtxt(S/"diagnostics.csv",delimiter=",",names=True)
    ref=np.load(D/"target.npz")
    ini=profiles(initial,ref);end=profiles(final,ref)
    np.savetxt(D/"profile_initial.csv",ini,delimiter=",",header="R,N_disk,Sigma,sigma_R,Q",comments="")
    np.savetxt(D/"profile_final.csv",end,delimiter=",",header="R,N_disk,Sigma,sigma_R,Q",comments="")
    def exact(path):return float(path.read_text().strip().split(",")[-1])
    Ei=exact(D/"exact_initial.csv");Ef=exact(D/"exact_final.csv")
    de=100*(Ef-Ei)/abs(Ei)
    lz=100*(f["Lz"][-1]/f["Lz"][0]-1)
    dz=100*(m[-1]["zrms"]/m[0]["zrms"]-1)
    dr=100*(m[-1]["R50"]/m[0]["R50"]-1)
    good=(ini[:,1]>=60)&np.isfinite(ini[:,4])&(ini[:,0]>.12)&(ini[:,0]<.65)
    q_init=float(np.median(ini[good,4])) if np.any(good) else float("nan")
    checks=dict(finite=np.isfinite([de,lz,dz,dr,q_init]).all(),
                exact_energy_abs_pct_below_0p1=abs(de)<.1,
                total_Lz_abs_pct_below_0p05=abs(lz)<.05,
                vertical_rms_abs_pct_below_5=abs(dz)<5,
                R50_abs_pct_below_2=abs(dr)<2,
                measured_initial_Q_median_between_1_and_2p5=1<q_init<2.5)
    status="PASS_SHORT_EQUILIBRIUM_SCREEN_NOT_FULL_REPRODUCTION" if all(checks.values()) else "FAILED_EQUILIBRIUM_SCREEN"
    result=dict(status=status,passed=bool(all(checks.values())),N=int(len(initial)),
      duration_Myr=300,integrator_timestep_Myr=.25,snapshot_Myr=5,
      exact_energy_change_pct=float(de),Lz_change_pct=float(lz),
      disc_z_rms_change_pct=float(dz),disc_R50_change_pct=float(dr),
      measured_Q_initial_median=q_init,
      epochs_Myr=[0,150,300],structures=m,checks={k:bool(v) for k,v in checks.items()},
      caveats=["AGAMA initial DF is in smooth potential; simulation applies Plummer softening",
      "Model halo and truncation are assumptions, not original paper's Fig.1 halo",
      "300 Myr stationarity cannot establish t=60 stability",
      "N=10000 is not converged; force, softening, timestep and N study still required",
      "No bar intentionally imposed or established"])
    (D/"validation.json").write_text(json.dumps(result,indent=2,allow_nan=False))
    print(json.dumps(result,indent=2))
    plt.rcParams.update({"font.family":"serif","font.size":9,"xtick.direction":"in","ytick.direction":"in"})
    fig,ax=plt.subplots(2,3,figsize=(12,7),layout="constrained")
    for j,(a,t) in enumerate(zip((initial,mid,final),(0,150,300))):
        d=a[a[:,8]==0];x,y,z=d[:,0]*15,d[:,1]*15,d[:,2]*15
        ax[0,j].hist2d(x,y,bins=105,range=[[-18,18],[-18,18]],cmap="magma",vmin=0,vmax=15)
        ax[1,j].hist2d(x,z,bins=(105,80),range=[[-18,18],[-5,5]],cmap="magma",vmin=0,vmax=20)
        ax[0,j].set(title=f"{t} Myr",xlabel="x [kpc]",ylabel="y [kpc]",aspect="equal")
        ax[1,j].set(xlabel="x [kpc]",ylabel="z [kpc]")
    fig.suptitle(f"Live 10,000-particle AGAMA DF candidate — {status}\nNo bar imposed; not an exact Noguchi (2001) reproduction")
    fig.savefig(D/"three_epochs.png",dpi=160);plt.close(fig)
    fig,ax=plt.subplots(2,2,figsize=(10,7),layout="constrained")
    ax[0,0].plot(f["t"]*75,100*(f["E"]-f["E"][0])/abs(f["E"][0]),label="tree-estimated")
    ax[0,0].scatter([300],[de],label="exact endpoint",color="black",zorder=5)
    ax[0,1].plot(f["t"]*75,100*(f["Lz"]/f["Lz"][0]-1))
    for key,col,ylabel in (("zrms",0,"RMS z [kpc]"),("R50",1,"R50 [kpc]")):
        ax[1,col].plot([0,150,300],[q[key]*15 for q in m],marker="o")
        ax[1,col].set_ylabel(ylabel)
    ax[0,0].set_ylabel("energy drift [%]");ax[0,1].set_ylabel("Lz drift [%]")
    for a in ax.flat:a.set_xlabel("Time [Myr]");a.grid(alpha=.2)
    ax[0,0].legend(fontsize=7)
    fig.suptitle(f"300 Myr live N-body stationarity diagnostics — {status}")
    fig.savefig(D/"stationarity.png",dpi=165);plt.close(fig)
    fig,ax=plt.subplots(1,2,figsize=(10,4),layout="constrained")
    for p,l in ((ini,"initial"),(end,"final")):
        ax[0].semilogy(p[:,0]*15,p[:,2],marker="o",label=l)
        ax[1].plot(p[:,0]*15,p[:,4],marker="o",label=l)
    ax[1].axhline(1.5,color="black",ls="--",label="nominal target Q=1.5")
    for a in ax:a.set_xlabel("R [kpc]");a.legend(fontsize=8);a.grid(alpha=.2)
    ax[0].set_ylabel("Disc surface density [code mass/area]");ax[1].set_ylabel("Measured thin-disc Q proxy")
    fig.suptitle("Sampled density and Q — finite N; no claim of exact Q=1.5")
    fig.savefig(D/"profiles.png",dpi=160);plt.close(fig)
    # Do not make the workflow fail before it commits a *failed* numerical report.
