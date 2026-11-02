"""Render three epochs, diagnostics and pair of low-memory GIFs."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image,ImageDraw
P=Path(__file__).resolve().parent
runs={name:np.load(P/'results'/name/'evolution.npz') for name in ('unforced','forced')}
plt.rcParams.update({'font.size':10,'axes.linewidth':1,'xtick.direction':'in','ytick.direction':'in'})
fig,ax=plt.subplots(2,3,figsize=(12.5,8),sharex=True,sharey=True,layout='constrained')
for row,name in enumerate(runs):
 for col,i in enumerate((0,450,900)):
  xy=runs[name]['x'][i]*15
  ax[row,col].hist2d(xy[:,0],xy[:,1],bins=110,range=[[-18,18],[-18,18]],cmap='inferno',vmin=0,vmax=9)
  ax[row,col].set_aspect('equal')
  ax[row,col].set(title=f'{name.capitalize()} • {i*5:,} Myr',xlabel='x [kpc]',ylabel='y [kpc]')
fig.suptitle('3,000 NON-SELF-GRAVITATING disc tracers in an analytic galactic potential\nExternal m=2 field in lower row — not a Noguchi (2001) reproduction',fontsize=12)
fig.savefig(P/'three_snapshots_control.png',dpi=145);plt.close(fig)
fig,ax=plt.subplots(2,2,figsize=(11,7.3),layout='constrained')
for name,data in runs.items():
 h=data['history'];tt=h[:,0];label='Axisymmetric control' if name=='unforced' else 'External rotating bar'
 for a,k in ((ax[0,0],4),(ax[0,1],5),(ax[1,0],6)):a.plot(tt,h[:,k],label=label)
 j=np.searchsorted(tt,300)
 y=100*(h[:,1]/h[0,1]-1) if name=='unforced' else 100*(h[:,3]-h[j,3])/abs(h[j,3])
 ax[1,1].plot(tt,y,label='Energy drift' if name=='unforced' else 'Jacobi drift (post-ramp)')
for a,l in ((ax[0,0],'RMS vertical height'),(ax[0,1],'Half-mass radius'),(ax[1,0],'m=2 amplitude'),(ax[1,1],'Conservation [%]')):a.set_ylabel(l)
for a in ax.flat:a.set_xlabel('Time [Myr]');a.legend(fontsize=8);a.grid(alpha=.2)
fig.suptitle('Orbital integration control — NOT SELF-GRAVITATING',fontsize=12)
fig.savefig(P/'control_diagnostics.png',dpi=155);plt.close(fig)
for name,data in runs.items():
 frames=[]
 for i in range(0,901,10):  # GIF 50 Myr; underlying snapshots still 5 Myr
  xyz=data['x'][i]*15;im=Image.new('RGB',(760,382),(10,12,21));draw=ImageDraw.Draw(im)
  for panel,(cx,cy) in enumerate(((192,196),(570,196))):
   draw.rectangle((cx-171,26,cx+171,363),outline=(74,83,104),width=1)
   xx=xyz[:,0];yy=xyz[:,1] if panel==0 else xyz[:,2]
   pix_x=np.rint(cx+xx*9).astype(int)
   pix_y=np.rint(cy-yy*(9 if panel==0 else 30)).astype(int)
   good=(pix_x>=cx-169)&(pix_x<=cx+169)&(pix_y>28)&(pix_y<360)
   for x,y in zip(pix_x[good],pix_y[good]):draw.point((int(x),int(y)),fill=(250,171,90))
  if name=='forced':
   angle=1.5*(i*5/75.)
   draw.line((192-32*np.cos(angle),196+32*np.sin(angle),192+32*np.cos(angle),196-32*np.sin(angle)),fill=(70,213,239),width=3)
  draw.text((80,8),'FACE ON',fill='white');draw.text((460,8),'EDGE ON',fill='white')
  draw.text((18,368),f'{name.upper()} | t={i*5:4d} Myr | EXTERNAL FIELD IF FORCED | NO SELF-GRAVITY',fill=(200,203,210))
  frames.append(im.quantize(colors=32,method=Image.Quantize.FASTOCTREE))
 frames[0].save(P/(name+'_control.gif'),save_all=True,append_images=frames[1:],duration=100,loop=0,disposal=2,optimize=False)
