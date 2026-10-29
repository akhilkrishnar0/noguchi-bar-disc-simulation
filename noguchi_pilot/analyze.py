import numpy as np,json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.ndimage import gaussian_filter
from PIL import Image
from pathlib import Path
P=Path(__file__).parent
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.facecolor':'white'})
def load(p):return np.fromfile(p,dtype=np.float64).reshape(-1,9)
def disc(a):return a[a[:,8]==0]
def center(a):
 # Use total center of mass for consistent conserved-coordinate diagnostics.
 return np.sum(a[:,:3]*a[:,6,None],axis=0)/np.sum(a[:,6])
def metrics(a):
 d=disc(a);xyz=d[:,:3]-center(a);r=np.linalg.norm(xyz[:,:2],axis=1);mask=(r>.05)&(r<.6);az=np.arctan2(xyz[:,1],xyz[:,0])
 return {'disc_R50_kpc':float(np.median(r)*15),'disc_z_rms_kpc':float(np.sqrt(np.mean(xyz[:,2]**2))*15),'A2_0p05_to_0p6':float(abs(np.mean(np.exp(2j*az[mask])))),'N_A2':int(mask.sum())}
def density(a,side=False):
 d=disc(a);xyz=(d[:,:3]-center(a))*15
 extent=[-15,15,-3,3] if side else [-15,15,-15,15]
 H,*_=np.histogram2d(xyz[:,0],xyz[:,2 if side else 1],bins=(220,100 if side else 220),range=[extent[:2],extent[2:]],weights=d[:,6]*1.5e11)
 area=30/220*(6/100 if side else 30/220)*1e6
 return gaussian_filter(H.T/area,1),extent
runs={key:np.genfromtxt(P/key/'diagnostics.csv',delimiter=',',names=True) for key in ['baseline','half_step']}
paths=sorted((P/'half_step').glob('snapshot_*.bin'),key=lambda p:int(p.stem.split('_')[-1]))
start=load(paths[0]);end=load(paths[-1]);tfinal=int(paths[-1].stem.split('_')[-1])*.0025
fig,axes=plt.subplots(2,2,figsize=(10,7),layout='constrained')
for row,side in enumerate([False,True]):
 for col,(a,t) in enumerate([(start,0),(end,tfinal)]):
  H,ex=density(a,side);im=axes[row,col].imshow(H,origin='lower',extent=ex,cmap='magma',norm=LogNorm(vmin=.1,vmax=1000),aspect='auto' if side else 'equal');axes[row,col].set(xlabel='x [kpc]',ylabel='z [kpc]' if side else 'y [kpc]',title=f'{"Edge-on" if side else "Face-on"} · t = {t*75:.0f} Myr')
fig.colorbar(im,ax=axes,label='Stellar disc surface density [solar masses / pc²]',shrink=.8)
fig.suptitle('First live N-body pilot · 59,985 particles\nApproximate Model A initialization; half timestep (0.0025)',fontsize=14)
fig.savefig(P/'first_run_snapshots.png',dpi=180);plt.close(fig)
# Physical diagnostics, including numerical convergence rather than claiming a bar.
fig,ax=plt.subplots(2,2,figsize=(10,7),layout='constrained')
colors={'baseline':'#c26330','half_step':'#266b9c'}
summary={}
for key,f in runs.items():
 ax[0,0].plot(f['t']*75,100*(f['E']-f['E'][0])/abs(f['E'][0]),label=key.replace('_',' '),color=colors[key]);ax[0,1].plot(f['t']*75,100*(f['Lz']/f['Lz'][0]-1),color=colors[key],label=key.replace('_',' '))
 snaps=sorted((P/key).glob('snapshot_*.bin'),key=lambda p:int(p.stem.split('_')[-1]));dt=.005 if key=='baseline' else .0025
 mm=[metrics(load(q)) for q in snaps];tt=[int(q.stem.split('_')[-1])*dt*75 for q in snaps]
 ax[1,0].plot(tt,[m['A2_0p05_to_0p6'] for m in mm],color=colors[key],label=key.replace('_',' '))
 summary[key]={'initial':mm[0],'final':mm[-1],'tree_energy_drift_fraction':float((f['E'][-1]-f['E'][0])/abs(f['E'][0])),'Lz_drift_fraction':float(f['Lz'][-1]/f['Lz'][0]-1)}
for a,label in [(start,'initial'),(end,'half-step final')]:
 d=disc(a);r=np.linalg.norm((d[:,:3]-center(a))[:,:2],axis=1);edges=np.linspace(.02,1.2,45);mass,_=np.histogram(r,bins=edges,weights=d[:,6]);sigma=mass/(np.pi*np.diff(edges**2))*1.5e11/(15000**2)
 ax[1,1].semilogy((edges[1:]+edges[:-1])/2*15,sigma,label=label)
ax[0,0].set(xlabel='Time [Myr]',ylabel='Tree-estimated energy change [%]');ax[0,1].set(xlabel='Time [Myr]',ylabel='Angular momentum Lz change [%]');ax[1,0].set(xlabel='Time [Myr]',ylabel='Global m=2 amplitude, 0.75 < R < 9 kpc',ylim=(0,None));ax[1,1].set(xlabel='Radius [kpc]',ylabel='Disc surface density [solar masses / pc²]')
for a in ax.flat:a.legend(fontsize=8);a.grid(alpha=.2)
fig.suptitle('Numerical and structural checks · short pilot, not a bar-growth result',fontsize=13)
fig.savefig(P/'first_run_diagnostics.png',dpi=180);plt.close(fig)
frames=[]
for q in paths:
 a=load(q);t=int(q.stem.split('_')[-1])*.0025
 fig,axes=plt.subplots(1,2,figsize=(9,4.2),layout='constrained')
 for side,ax in zip([False,True],axes):
  H,ex=density(a,side);im=ax.imshow(H,origin='lower',extent=ex,cmap='magma',norm=LogNorm(vmin=.1,vmax=1000),aspect='auto' if side else 'equal');ax.set(xlabel='x [kpc]',ylabel='z [kpc]' if side else 'y [kpc]',title='Edge-on' if side else 'Face-on')
 fig.suptitle(f'Live N-body pilot · t = {t*75:5.1f} Myr\nApproximate initialization · no imposed bar · fixed density colour scale',fontsize=12)
 fig.canvas.draw();frames.append(Image.fromarray(np.asarray(fig.canvas.buffer_rgba()).copy()).convert('RGB'));plt.close(fig)
frames[0].save(P/'first_run.gif',save_all=True,append_images=frames[1:],duration=180,loop=0)
for key in ['baseline','tight_force']:
 f=np.genfromtxt(P/key/'force_validation.csv',delimiter=',',names=True)
 summary.setdefault(key,{})['force_error_fraction_percentiles_50_95_max']=np.percentile(f['relative_error'],[50,95,100]).tolist()

for key in ['initial','baseline','half_step']:
 text=(P/f'exact_{key}.csv').read_text().strip().split(',')
 summary.setdefault(key,{})['exact_energy']=float(text[-1])
 summary[key]['exact_K']=float(text[-3]);summary[key]['exact_U']=float(text[-2])
for key in ['baseline','half_step']:
 summary[key]['exact_energy_drift_fraction']=(summary[key]['exact_energy']-summary['initial']['exact_energy'])/abs(summary['initial']['exact_energy'])
summary['duration_code_units']=tfinal
summary['duration_Myr']=75*tfinal
(P/'results.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
