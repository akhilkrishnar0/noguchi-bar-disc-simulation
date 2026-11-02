#!/usr/bin/env python3
"""Science gate: identical AGAMA DF realization under three numerical configurations."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
P=Path(__file__).resolve().parent/"results"
CASES={"baseline":(1200,600),"refined":(2400,1200),"smaller_eps":(2400,1200)}
def get(case,step):
    a=np.fromfile(P/case/f"snapshot_{step}.bin",dtype="<f8")
    if a.size%9!=0:raise RuntimeError(f"corrupt snapshot {case}/{step}")
    return a.reshape(-1,9)
def metrics(a):
    centre=np.average(a[:,:3],axis=0,weights=a[:,6])
    d=a[a[:,8]==0];xyz=d[:,:3]-centre
    r=np.hypot(xyz[:,0],xyz[:,1])
    selection=(r>.05)&(r<.6)
    a2=np.abs(np.mean(np.exp(2j*np.arctan2(xyz[selection,1],xyz[selection,0]))))
    return dict(R50_code=float(np.median(r)),
                zrms_code=float(np.sqrt(np.mean(xyz[:,2]**2))),A2=float(a2))
def exact(f):
    return float((P/f).read_text().strip().split(",")[-1])
def run(case,steps,half):
    m=[metrics(get(case,j)) for j in (0,half,steps)]
    data=np.genfromtxt(P/case/"diagnostics.csv",delimiter=",",names=True)
    Ei=exact(f"{case}/energy_initial.csv"); Ef=exact(f"{case}/energy_final.csv")
    values=dict(
      epsilon_scale=1.0 if case!="smaller_eps" else .5,
      theta=.5 if case=="baseline" else .35,
      step_myr=.25 if case=="baseline" else .125,
      energy_drift_exact_percent=100*(Ef-Ei)/abs(Ei),
      total_Lz_drift_percent=100*(data["Lz"][-1]/data["Lz"][0]-1),
      disc_zrms_change_percent=100*(m[-1]["zrms_code"]/m[0]["zrms_code"]-1),
      disc_R50_change_percent=100*(m[-1]["R50_code"]/m[0]["R50_code"]-1),
      metrics_0_150_300_Myr=m)
    flags={"energy_abs_lt_0p1pct":abs(values["energy_drift_exact_percent"])<.1,
           "Lz_abs_lt_0p05pct":abs(values["total_Lz_drift_percent"])<.05,
           "vertical_abs_lt_5pct":abs(values["disc_zrms_change_percent"])<5,
           "R50_abs_lt_2pct":abs(values["disc_R50_change_percent"])<2}
    values["checks"]={k:bool(x) for k,x in flags.items()}
    values["status"]="PASS_SHORT_SCREEN_ONLY" if all(flags.values()) else "FAILED_SHORT_SCREEN"
    return values,data
res={}; histories={}
for case,(last,half) in CASES.items():
    res[case],histories[case]=run(case,last,half)
    print(case,json.dumps(res[case],indent=2),flush=True)
force=json.loads((P/"force_comparison.json").read_text())
report=dict(description="300 Myr live-disc/halo/bulge force-softening diagnostic, same AGAMA particle phase-space realization",
   baseline_short_screen_passed=res["baseline"]["status"]=="PASS_SHORT_SCREEN_ONLY",
   refined_short_screen_passed=res["refined"]["status"]=="PASS_SHORT_SCREEN_ONLY",
   smaller_softening_short_screen_passed=res["smaller_eps"]["status"]=="PASS_SHORT_SCREEN_ONLY",
   sample_particle_count=10000, force_comparison=force,
   runs=res,
   interpretive_limits=[
      "Refined case changes both timestep and opening angle, so causes of changes are not separately identified",
      "Smaller-epsilon test changes softened physics, not only numerical precision",
      "Finite particle discreteness and gravitational softening prevent the DF from being exact",
      "A 300 Myr screen does not establish long-duration equilibrium or bar formation",
      "No physically converged N, softening, timestep or random-seed series yet performed"])
(P/"fix_validation.json").write_text(json.dumps(report,indent=2,allow_nan=False))
plt.rcParams.update({"font.size":10,"font.family":"serif","xtick.direction":"in","ytick.direction":"in"})
fig,ax=plt.subplots(2,3,figsize=(12.5,7.5),layout="constrained",sharex="col")
for row in (0,1):
    target="baseline" if row==0 else "smaller_eps"
    last,half=CASES[target]
    for col,step in enumerate((0,half,last)):
        d=get(target,step);d=d[d[:,8]==0]
        x,y,z=d[:,0]*15,d[:,1]*15,d[:,2]*15
        ax[row,col].hist2d(x,y if row==0 else z,bins=(120,90),range=[[-18,18],[-18,18] if row==0 else [-5,5]],cmap="magma",vmin=0,vmax=20)
        ax[row,col].set(title=f"{target}: {col*150} Myr",xlabel="x [kpc]",ylabel="y [kpc]" if row==0 else "z [kpc]")
        if row==0:ax[row,col].set_aspect("equal")
fig.suptitle("Live galaxy 300 Myr softening test — no bar imposed\nSame initial particle phase-space realization, softened force differs")
fig.savefig(P/"three_epochs_comparison.png",dpi=155);plt.close(fig)
fig,axes=plt.subplots(2,2,figsize=(10,7),layout="constrained")
for key,dd in histories.items():
    t=dd["t"]*75
    axes[0,0].plot(t,100*(dd["E"]-dd["E"][0])/abs(dd["E"][0]),label=key)
    axes[0,1].plot(t,100*(dd["Lz"]/dd["Lz"][0]-1),label=key)
    mn=[metrics(get(key,st)) for st in (0,CASES[key][1],CASES[key][0])]
    axes[1,0].plot((0,150,300),[v["zrms_code"]*15 for v in mn],marker="o",label=key)
    axes[1,1].plot((0,150,300),[v["R50_code"]*15 for v in mn],marker="o",label=key)
for a,y in zip(axes.flat,["Tree energy drift [%]","Total Lz drift [%]","Disc RMS height [kpc]","Disc R50 [kpc]"]):
    a.set(xlabel="Time [Myr]",ylabel=y);a.grid(alpha=.25);a.legend(fontsize=8)
fig.suptitle("Numerical/softening sensitivity: exact endpoint energy values are in JSON")
fig.savefig(P/"numerical_diagnostics.png",dpi=155);plt.close(fig)
print("Validation file:",P/"fix_validation.json",flush=True)
