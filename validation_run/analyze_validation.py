import numpy as np,json
from pathlib import Path
from scipy.interpolate import CubicSpline
from scipy.ndimage import gaussian_filter1d
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).parent;old=P.parent/'noguchi_pilot'
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
def load(p):return np.fromfile(p,dtype='<f8').reshape(-1,9)
def measure(p,acc=None):
 com=np.sum(p[:,:3]*p[:,6,None],axis=0);bulk=np.sum(p[:,3:6]*p[:,6,None],axis=0)
 d=p[p[:,8]==0];x=d[:,:3]-com;v=d[:,3:6]-bulk;r=np.hypot(x[:,0],x[:,1]);az=np.arctan2(x[:,1],x[:,0]);s=(r>.05)&(r<.6)
 out={'R50_kpc':float(np.median(r)*15),'z_rms_kpc':float(np.sqrt(np.mean(x[:,2]**2))*15),'A2':float(abs(np.mean(np.exp(2j*az[s])))),'outside_R1_fraction':float(np.mean(r>1)),'Lz_components':[],'component_structure':{}}
 for c in range(3):
  d=p[p[:,8]==c];out['Lz_components'].append(float(np.sum(d[:,6]*(d[:,0]*d[:,4]-d[:,1]*d[:,3]))))
 for c in range(3):
  d=p[p[:,8]==c];xx=d[:,:3]-com;rr=np.linalg.norm(xx,axis=1);core=rr<np.quantile(rr,.9)
  out['component_structure'][str(c)]={'spherical_R50_kpc':float(np.median(rr)*15),'q_inertia_inner90':float(np.sqrt(2*np.mean(xx[core,2]**2)/np.mean(xx[core,0]**2+xx[core,1]**2)))}
 if acc is not None:
  out['virial_tensor_support_ratios']={}
  for c in range(3):
   s=p[:,8]==c;m=p[s,6,None];x=p[s,:3]-com;v=p[s,3:6]-bulk;aa=acc[s]
   tt=np.sum(m*v*v,axis=0);ww=-np.sum(m*x*aa,axis=0)
   out['virial_tensor_support_ratios'][str(c)]=(tt/ww).tolist()
 return out
runs={'prior theta 0.8':(old/'half_step',.0025),'live theta 0.5':(P/'live_tight',.0025),'frozen field':(P/'frozen_field',.0025),'dt 0.00125':(P/'small_step',.00125)}
colors=['#999999','#2476a7','#cb6c32','#4c9852'];summary={};series={}
for label,(folder,dt) in runs.items():
 paths=sorted(folder.glob('snapshot_*.bin'),key=lambda x:int(x.stem.split('_')[-1]));arr=[]
 for path in paths:
  ap=folder/path.name.replace('snapshot','acceleration');a=np.fromfile(ap).reshape(-1,3) if ap.exists() else None
  m=measure(load(path),a);m['t']=int(path.stem.split('_')[-1])*dt;arr.append(m)
 series[label]=arr;summary[label]={'initial':arr[0],'final':arr[-1]}
 for name in ['force_validation.csv','force_final.csv']:
  f=folder/name
  if f.exists():
   v=np.genfromtxt(f,delimiter=',',names=True)['relative_error'];summary[label][name]=np.percentile(v,[50,95,100]).tolist()
 f=np.genfromtxt(folder/'diagnostics.csv',delimiter=',',names=True)
 summary[label]['Lz_change']=float(f['Lz'][-1]-f['Lz'][0]);summary[label]['Lz_drift_percent']=float((f['Lz'][-1]/f['Lz'][0]-1)*100)
 # A frozen nonaxisymmetric external field can exert a real torque; Lz is not conserved there.
 if label!='prior theta 0.8':
  e=np.loadtxt(folder/'exact_final.csv',delimiter=',',usecols=(1,2,3)) if label!='frozen field' else np.loadtxt(folder/'exact_final.csv',delimiter=',')
  e0=np.loadtxt(old/'exact_initial.csv',delimiter=',',usecols=(1,2,3));ref=e0[0]+2*e0[1] if label=='frozen field' else e0[2]
  summary[label]['exact_energy_drift_percent']=float(100*(e[-1]-ref)/abs(ref))
# Equal-time comparisons isolate time step and tree angle.
def at(label,t):return min(series[label],key=lambda a:abs(a['t']-t))
summary['comparisons']={}
for label,t in [('prior theta 0.8',1),('dt 0.00125',.25),('frozen field',1)]:
 a=at('live theta 0.5',t);b=at(label,t)
 summary['comparisons'][label]={'time':t,'relative_z_difference_percent':100*(b['z_rms_kpc']/a['z_rms_kpc']-1),'relative_R50_difference_percent':100*(b['R50_kpc']/a['R50_kpc']-1),'absolute_A2_difference':b['A2']-a['A2']}
# Diagnose classical Q from the directly measured softened midplane force.
profiles={}
for key,binpath in [('initial',old/'initial.bin'),('final',P/'live_tight/snapshot_400.bin')]:
 ring=np.genfromtxt(P/f'probe_{key}.csv',delimiter=',',names=True);rr=ring['R'];vc2=-rr*ring['mean_ar'];spl=CubicSpline(rr,gaussian_filter1d(vc2,1))
 p=load(binpath);com=np.sum(p[:,:3]*p[:,6,None],axis=0);d=p[p[:,8]==0];x=d[:,:3]-com;r=np.hypot(x[:,0],x[:,1]);vR=np.sum(x[:,:2]*d[:,3:5],axis=1)/r;vphi=(x[:,0]*d[:,4]-x[:,1]*d[:,3])/r
 edges=np.linspace(.05,.85,17);mid=(edges[1:]+edges[:-1])/2;Q=[];disp=[];rotation=[];counts=[]
 for lo,hi,R in zip(edges[:-1],edges[1:],mid):
  s=(r>=lo)&(r<hi);sigma=np.std(vR[s],ddof=1);Sigma=d[s,6].sum()/(np.pi*(hi*hi-lo*lo));kap2=spl(R,1)/R+2*spl(R)/R**2
  Q.append(sigma*np.sqrt(max(0,kap2))/(3.36*Sigma));disp.append(sigma);rotation.append(np.mean(vphi[s]));counts.append(int(s.sum()))
 profiles[key]={'r':mid.tolist(),'Q':Q,'sigmaR':disp,'mean_vphi':rotation,'counts':counts}
summary['Q_profiles']=profiles
fig,ax=plt.subplots(2,2,figsize=(11,8),layout='constrained')
for (label,arr),color in zip(series.items(),colors):
 t=np.array([a['t'] for a in arr])*75
 ax[0,0].plot(t,[a['z_rms_kpc'] for a in arr],label=label,color=color)
 ax[0,1].plot(t,[100*(a['R50_kpc']/arr[0]['R50_kpc']-1) for a in arr],color=color)
 ax[1,0].plot(t,[a['A2'] for a in arr],color=color)
ax[0,0].set(xlabel='Time [Myr]',ylabel='Disc RMS height [kpc]');ax[0,0].legend(fontsize=8)
ax[0,1].set(xlabel='Time [Myr]',ylabel='Half-mass radius change [%]')
ax[1,0].set(xlabel='Time [Myr]',ylabel='Global A₂, 0.75 < R < 9 kpc')
for key,style in [('initial','-'),('final','--')]:
 q=profiles[key];ax[1,1].plot(np.array(q['r'])*15,q['Q'],style,label=key)
ax[1,1].axhline(1.5,color='gray',ls=':',label='initial target 1.5');ax[1,1].legend(fontsize=8);ax[1,1].set(xlabel='Radius [kpc]',ylabel='Classical Q diagnostic')
for a in ax.flat:a.grid(alpha=.2)
fig.suptitle('Equilibrium validation: same particles, controlled numerical changes\nFrozen field tests initial orbital balance; it cannot develop a collective bar',fontsize=13)
fig.savefig(P/'validation_results.png',dpi=170);plt.close(fig)
# Vertical phase-space plot illustrates disequilibrium without suggesting a bar.
fig,ax=plt.subplots(1,3,figsize=(12,4),layout='constrained')
for a,(title,path) in zip(ax,[('Initial',old/'initial.bin'),('Live, tighter forces',P/'live_tight/snapshot_400.bin'),('Frozen initial field',P/'frozen_field/snapshot_400.bin')]):
 d=load(path)[:30000];r=np.hypot(d[:,0],d[:,1]);s=(r>.1)&(r<.8)
 a.hist2d(d[s,2]*15,d[s,5]*200,bins=70,range=[[-1.2,1.2],[-100,100]],cmap='magma',cmin=1)
 a.set(title=title,xlabel='z [kpc]',ylabel='v_z [km/s]')
fig.suptitle('Vertical phase space · disc annulus 1.5–12 kpc · final panels at 75 Myr')
fig.savefig(P/'vertical_phase_space.png',dpi=170);plt.close(fig)
(P/'validation_results.json').write_text(json.dumps(summary,indent=2))
print(json.dumps({k:{a:b for a,b in v.items() if a in ['final','exact_energy_drift_percent','Lz_drift_percent']} for k,v in summary.items() if k in runs},indent=2))
