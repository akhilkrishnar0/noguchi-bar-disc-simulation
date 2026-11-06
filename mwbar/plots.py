"""Manuscript-style diagnostic figures and GIF generated from simulated orbits."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image,ImageDraw
from .model import Galaxy
from .measure import profile

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,
                     'axes.linewidth':.8,'xtick.direction':'in',
                     'ytick.direction':'in','savefig.dpi':160})

def figures(model:Galaxy,runs:dict,measured:dict,out:Path):
    out.mkdir(parents=True,exist_ok=True)
    # 2 cases × 4 equally spaced science-time epochs; strictly fixed axes.
    frames=[0,20,40,60]
    fig,ax=plt.subplots(2,4,figsize=(13.3,6.7),sharex=True,sharey=True,layout='constrained')
    for row,(case,run) in enumerate(runs.items()):
        for col,i in enumerate(frames):
            x=run['x'][i];tt=run['t'][i]
            ax[row,col].hist2d(x[:,0],x[:,1],bins=115,
                               range=[[-15,15],[-15,15]],cmap='inferno',vmin=0,vmax=13)
            ax[row,col].set_aspect('equal')
            ax[row,col].set(title=f'{tt:g} Myr',xlabel='x [kpc]',ylabel='y [kpc]')
            if col==0:
                ax[row,col].text(.03,.96,'BARRED' if case=='barred' else 'UNBARRED',
                    transform=ax[row,col].transAxes,ha='left',va='top',color='white',fontsize=11)
    fig.suptitle('Present-day Milky Way analogue — 10,000 collisionless disc tracers\n'
                 'Bar is an imposed rotating potential, not a self-gravitating stellar structure')
    fig.savefig(out/'face_on_epochs.png');plt.close(fig);print('FIG face-on',flush=True)

    fig,axes=plt.subplots(1,2,figsize=(11.5,4.3),sharex=True,sharey=True,layout='constrained')
    colors=plt.rcParams['axes.prop_cycle'].by_key()['color']
    for ax,(case,run) in zip(axes,runs.items()):
        for ci,i in enumerate(frames):
            rr,sigma,n=profile(model,run['x'][i])
            ax.plot(rr,sigma,'o-',lw=1,ms=2.5,label=f'{run["t"][i]:g} Myr',color=colors[ci])
        ax.set(title=case.title(),xlabel='Galactocentric radius [kpc]',
               ylabel=r'Tracer-weighted $\Sigma_\star$ [$M_\odot\,\mathrm{pc}^{-2}$]',yscale='log',
               xlim=(.5,15.5))
        ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.suptitle('Radial stellar-tracer surface-density proxies (not live disc mass)')
    fig.savefig(out/'surface_density.png');plt.close(fig);print('FIG surface density',flush=True)

    fig,ax=plt.subplots(2,2,figsize=(11,7.4),layout='constrained',sharex=True)
    for case,(_,data,_) in measured.items():
        t=data[:,0];label=case.title()
        ax[0,0].plot(t,data[:,3],label=label,lw=1.5)
        ax[0,1].plot(t,data[:,4],label=label,lw=1.5)
        ax[1,0].plot(t,data[:,8],label=label,lw=1.5)
        ax[1,1].plot(t,100*(data[:,6]/data[0,6]-1),label=label,lw=1.5)
    labels=['Fitted outer $R_d$ [kpc]',r'Inner $|A_2|$ ($0.5<R<6$ kpc)',
            r'Median $|\Delta R_g|$ since $t=0$ [kpc]',
            r'Change in mean tracer $L_z$ [%]']
    for a,name in zip(ax.flat,labels):
        a.set(xlabel='Science time [Myr]',ylabel=name);a.grid(alpha=.2);a.legend(fontsize=8)
    fig.suptitle('Response to an EXISTING bar: density profile, $m=2$, guiding-centre shifts and torque')
    fig.savefig(out/'structural_evolution.png');plt.close(fig);print('FIG structural',flush=True)

    fig,ax=plt.subplots(figsize=(8.2,4.7),layout='constrained')
    R=np.linspace(.4,16.,800)
    om,kap,_=model.frequencies(R)
    ax.plot(R,om/model.omega_bar*model.bar_pattern_kms_kpc,label=r'$\Omega$',lw=1.7)
    ax.plot(R,(om-kap/2)/model.omega_bar*model.bar_pattern_kms_kpc,label=r'$\Omega-\kappa/2$',lw=1.7)
    ax.plot(R,(om+kap/2)/model.omega_bar*model.bar_pattern_kms_kpc,label=r'$\Omega+\kappa/2$',lw=1.7)
    ax.axhline(model.bar_pattern_kms_kpc,ls='--',c='black',label=r'$\Omega_\mathrm{bar}$')
    for kind,v in model.resonance_radii().items():
        for rv in v:
            if rv<16:
                ax.axvline(rv,ls=':',alpha=.45)
                ax.text(rv+0.05,95,kind,rotation=90,fontsize=8,va='top')
    ax.set(xlabel='Radius [kpc]',ylabel=r'Pattern / orbital frequency [km s$^{-1}$ kpc$^{-1}$]',
           ylim=(-25,150),xlim=(.4,16.))
    ax.legend(loc='upper right',fontsize=8);ax.grid(alpha=.2)
    fig.suptitle('Resonances in the axisymmetric background (approximate epicycle analysis)')
    fig.savefig(out/'resonances.png');plt.close(fig);print('FIG resonances',flush=True)

    # Small reproducible GIF: a pair of true orbital snapshots every 10 Myr.
    # The two panels share the *same* sky-plane range and equal pixel scale.
    panels=[]
    for i in range(0,len(runs['barred']['t']),2):
        frame=Image.new('RGB',(780,425),(11,12,21));draw=ImageDraw.Draw(frame)
        for case,cx in [('unbarred',195),('barred',585)]:
            x=runs[case]['x'][i]
            ypix=np.rint(222-x[:,1]*11.5).astype(int)
            xpix=np.rint(cx+x[:,0]*11.5).astype(int)
            ok=(xpix>=cx-183)&(xpix<=cx+183)&(ypix>=37)&(ypix<=401)
            draw.rectangle((cx-185,36,cx+185,403),outline=(58,65,88),width=1)
            for px,py in zip(xpix[ok],ypix[ok]):
                draw.point((int(px),int(py)),fill=(255,177,103))
        bar_ang=model.omega_bar*runs['barred']['t'][i]+np.deg2rad(model.bar_angle_deg)
        dx=26*np.cos(bar_ang);dy=26*np.sin(bar_ang)
        draw.line((585-dx,222+dy,585+dx,222-dy),fill=(91,226,228),width=3)
        draw.text((105,16),'UNBARRED CONTROL',fill=(235,237,239))
        draw.text((484,16),'EXISTING BAR',fill=(235,237,239))
        draw.text((12,409),f'Science epoch {runs["barred"]["t"][i]:g} Myr  | 10,000 stellar tracers | fixed galaxy potential',fill=(201,203,212))
        panels.append(frame.quantize(colors=64,method=Image.Quantize.FASTOCTREE))
    print('GIF generated frames',len(panels),flush=True)
    panels[0].save(out/'barred_vs_unbarred.gif',save_all=True,
                   append_images=panels[1:],duration=130,loop=0,disposal=2,optimize=False)

