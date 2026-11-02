"""Create three-epoch plots and 5 Myr-cadence diagnostic GIFs.

All panels use fixed physical axes; external bar overlay denotes imposed
potential orientation, not a measured stellar bar.
"""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from PIL import Image, ImageDraw, ImageFont
P=Path(__file__).resolve().parent
runs=['live_60','forced_60']
def load(name,i):
    a=np.fromfile(P/name/f"snap_{i if i else '0000'}.bin",dtype='<f8').reshape(-1,9)
    return a[a[:,8]==0,:3]*15
fig,ax=plt.subplots(2,3,figsize=(12,7.3),sharex=True,sharey=True,constrained_layout=True)
for row,name in enumerate(runs):
    for col,i in enumerate((0,450,900)):
        xyz=load(name,i);x,y=xyz[:,0],xyz[:,1]
        h,_,_=np.histogram2d(x,y,bins=120,range=[[-65,65],[-65,65]])
        ax[row,col].imshow(h.T+0.1,origin='lower',extent=[-65,65,-65,65],cmap='magma',norm=LogNorm(vmin=.1,vmax=30),interpolation='nearest')
        ax[row,col].set_title(f'{"Unforced" if row==0 else "Externally forced"} · {i*5:,} Myr')
        ax[row,col].set_xlabel('x [kpc]')
        if col==0:ax[row,col].set_ylabel('y [kpc]')
        ax[row,col].set_aspect('equal')
fig.suptitle('Three observed simulation epochs · LOW N · EQUILIBRIUM FAILED\nFixed scale; point-noise and expansion are numerical/initialization-contaminated',fontsize=13)
fig.savefig(P/'three_snapshots.png',dpi=150);plt.close(fig)
fig,axes=plt.subplots(2,2,figsize=(11,7),constrained_layout=True)
summary={}
for name in runs:
    f=np.genfromtxt(P/name/'diagnostics.csv',delimiter=',',names=True)
    label='No imposed bar' if name=='live_60' else 'Forced rotating quadrupole'
    kw={'label':label}
    axes[0,0].plot(f['t']*75,f['zrms']*15,**kw)
    axes[0,1].plot(f['t']*75,f['A2'],**kw)
    axes[1,0].plot(f['t']*75,100*(f['E_self']-f['E_self'][0])/abs(f['E_self'][0]),**kw)
    axes[1,1].plot(f['t']*75,100*(f['Lz']/f['Lz'][0]-1),**kw)
    summary[name]={'N':1500,'duration_code_units':60,'duration_Myr':4500,'output_interval_Myr':5,'E_self_change_percent':float(100*(f['E_self'][-1]-f['E_self'][0])/abs(f['E_self'][0])),'Lz_change_percent':float(100*(f['Lz'][-1]/f['Lz'][0]-1)),'z_rms_kpc_start_end':[float(f['zrms'][0]*15),float(f['zrms'][-1]*15)],'A2_start_end':[float(f['A2'][0]),float(f['A2'][-1])],'stationary_equilibrium':False}
for a,title,y in zip(axes.flat,['Disc RMS height','Global m=2 amplitude','Self-gravity-only energy change','Angular momentum change'],['kpc','A2','percent','percent']):
    a.set(xlabel='Time [Myr]',ylabel=y,title=title);a.grid(alpha=.3);a.legend(fontsize=8)
fig.suptitle('LOW-RESOLUTION DIAGNOSTIC — not a validated bar or GLSB result\nExternal rotating potential supplies energy and angular momentum in forced case')
fig.savefig(P/'long_run_diagnostics.png',dpi=145);plt.close(fig)
font=ImageFont.load_default()
for name in runs:
    frames=[]
    for i in range(901):
        xyz=load(name,i); im=Image.new('RGB',(640,330),(9,12,24));dr=ImageDraw.Draw(im)
        # Fixed physical field: 130 x 130 kpc face-on and edge-on.
        for X,Y,origin in [(xyz[:,0],xyz[:,1],(165,168)),(xyz[:,0],xyz[:,2],(475,168))]:
            px=origin[0]+np.rint(X*2.2).astype(int);py=origin[1]-np.rint(Y*2.2).astype(int)
            for xx,yy in zip(px,py):
                if origin[0]-146<xx<origin[0]+146 and 21<yy<316: dr.point((int(xx),int(yy)),fill=(250,173,93))
        for center in (165,475):dr.rectangle((center-145,22,center+145,313),outline=(70,77,98),width=1)
        dr.text((80,6),'FACE-ON: 130 kpc',font=font,fill='white');dr.text((390,6),'EDGE-ON: 130 kpc',font=font,fill='white')
        if name=='forced_60':
            ang=1.8*(i*5/75);xx=165+int(28*np.cos(ang));yy=168-int(28*np.sin(ang));x2=165-int(28*np.cos(ang));y2=168+int(28*np.sin(ang))
            dr.line((x2,y2,xx,yy),fill=(85,215,255),width=2)
        dr.text((15,317),f'{"LIVE / NO BAR IMPOSED" if name=="live_60" else "IMPOSED ROTATING QUADRUPOLE"}  |  t={i*5:4d} Myr  |  NOT EQUILIBRATED',fill=(230,230,230),font=font)
        frames.append(im.quantize(colors=32,method=Image.Quantize.FASTOCTREE))
    outfile=P/('live_60.gif' if name=='live_60' else 'forced_bar_60.gif')
    frames[0].save(outfile,save_all=True,append_images=frames[1:],optimize=False,duration=45,loop=0,disposal=2)
    print(outfile,round(outfile.stat().st_size/1e6,2),'MB')
(P/'long_run_summary.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
