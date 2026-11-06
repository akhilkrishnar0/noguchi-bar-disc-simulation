"""Measured 2 Gyr tracer diagnostic figures and animations (no generated imagery)."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image,ImageDraw,ImageFont
from mwbar.measure import profile
from .measure import mu_r_from_sigma,M_TO_L_R,M_R_SUN_AB,fourier_radial

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.linewidth':0.8,
                     'xtick.direction':'in','ytick.direction':'in','savefig.dpi':165})


def _smoothed_hist(x,disc_mass,bins=84,extent=14.,smooth_pix=1.4):
    edges=np.linspace(-extent,extent,bins+1)
    counts,_,_=np.histogram2d(x[:,1],x[:,0],bins=(edges,edges))
    area=(edges[1]-edges[0])**2*1e6
    sigma=counts*(disc_mass/len(x))/area
    # 2D circular Gaussian smoothing; mass conserving, constant photometric assumption.
    half=int(np.ceil(4*smooth_pix));axis=np.arange(-half,half+1)
    kernel=np.exp(-.5*(axis/smooth_pix)**2);kernel/=kernel.sum()
    sigma=np.apply_along_axis(lambda a:np.convolve(np.pad(a,(half,half),mode='constant'),kernel,mode='valid'),axis=0,arr=sigma)
    sigma=np.apply_along_axis(lambda a:np.convolve(np.pad(a,(half,half),mode='constant'),kernel,mode='valid'),axis=1,arr=sigma)
    return mu_r_from_sigma(sigma)


def _epoch_indices(times,targets=(0,500,1000,2000)):
    return [int(np.argmin(np.abs(times-t))) for t in targets]


def plot_all(model,runs,diag,profiles,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    times=runs['barred']['t'];indices=_epoch_indices(times)
    colors=plt.rcParams['axes.prop_cycle'].by_key()['color']

    # Figure 1: 2Gyr tracer face-on evolution; scientifically common color scaling.
    fig,ax=plt.subplots(2,4,figsize=(14,7.0),sharex=True,sharey=True,layout='constrained')
    for row,case in enumerate(('unbarred','barred')):
        for col,idx in enumerate(indices):
            x=runs[case]['x'][idx]
            ax[row,col].hist2d(x[:,0],x[:,1],bins=115,range=[[-15,15],[-15,15]],
                                cmap='inferno',vmin=0,vmax=14)
            ax[row,col].set_aspect('equal')
            ax[row,col].set(title=f't={times[idx]/1000:g} Gyr',xlabel='x [kpc]',ylabel='y [kpc]')
            if col==0:ax[row,col].text(.02,.98,case.upper(),va='top',ha='left',transform=ax[row,col].transAxes,color='white')
    fig.suptitle('Milky Way z~0: 10,000 stellar tracers, pre-existing externally prescribed bar\n'
                 'No live galaxy self-gravity, star formation, or bar formation')
    fig.savefig(out/'face_on_2gyr.png');plt.close(fig)

    # Figure 2: force-defined bar Qb versus tracer response A2 and Fourier extent.
    fig,ax=plt.subplots(2,2,figsize=(11,7.2),layout='constrained',sharex=True)
    for case,style in (('unbarred','-'),('barred','--')):
        h,y=diag[case];lab=case.title();time=y[:,0]
        ax[0,0].plot(time,y[:,1],linestyle=style,label=lab,linewidth=1.5)
        ax[0,1].plot(time,y[:,2],linestyle=style,label=lab,linewidth=1.4)
    b=diag['barred'][1]
    ext=b[:,3]
    ax[1,0].plot(b[:,0],ext,label='Coherent m=2 tracer-response extent',lw=1.25)
    ax[1,0].axhline(model.bar_length,color='black',ls=':',label='4.8 kpc reference potential scale')
    ax[1,1].axhline(model.bar_force_fraction,color='black',linewidth=2,label='Imposed peak force fraction Qb (constant)')
    ax[1,1].set(ylim=(0,.15))
    lab=['Inner tracer $|A_2|$ (0.5–6 kpc)', 'Inner peak tracer $|A_2(R)|$',
         'm=2 tracer-response extent [kpc]', 'Prescribed bar force strength $Q_b$']
    for p,l in zip(ax.flat,lab):p.set(xlabel='Science time [Myr]',ylabel=l);p.grid(alpha=.2);p.legend(fontsize=8)
    fig.suptitle('Stellar response vs *fixed* imposed bar strength/scale\n'
                 'The measured extent is a Fourier threshold proxy, not the gravitational bar length')
    fig.savefig(out/'bar_strength_length_2gyr.png');plt.close(fig)

    # Figure 3: brightness radial profiles from constant assumed M/L.
    fig,axes=plt.subplots(1,2,figsize=(12,5.1),layout='constrained',sharex=True,sharey=True)
    for a,case in zip(axes,('unbarred','barred')):
        for j,i in enumerate(indices):
            r,s,n=profile(model,runs[case]['x'][i])
            mu=mu_r_from_sigma(s)
            a.plot(r,mu,'o-',ms=3,lw=1.3,color=colors[j],label=f'{times[i]/1000:g} Gyr')
        a.set(title=case.title(),xlabel='Radius [kpc]',ylabel=r'Illustrative $\mu_r$ [AB mag arcsec$^{-2}$]',
              xlim=(.5,15.5),ylim=(31,19))
        a.grid(alpha=.2);a.legend(fontsize=8)
    fig.suptitle(r'Face-on r-band surface-brightness proxy: fixed $M_\star/L_r=1.5$, $M_{r,\odot}=4.65$'
                 '\nEqual-weight massless tracers; no stellar-population evolution, dust or star formation')
    fig.savefig(out/'brightness_profiles_2gyr.png');plt.close(fig)

    # Figure 4: 2D map of the same proxy and same brightness range, no image synthesis.
    fig,ax=plt.subplots(2,4,figsize=(14,7.1),layout='constrained',sharex=True,sharey=True)
    ims=[]
    for row,case in enumerate(('unbarred','barred')):
        for col,idx in enumerate(indices):
            mu=_smoothed_hist(runs[case]['x'][idx],model.disc_mass)
            # Empty pixels mean no tracer luminosity; render at faint-limit
            # instead of matplotlib's white missing-data colour.
            mu=np.nan_to_num(mu,nan=28.5,posinf=28.5)
            im=ax[row,col].imshow(mu,origin='lower',extent=[-14,14,-14,14],
                 cmap='magma_r',vmin=19.5,vmax=28.5,interpolation='nearest')
            ims.append(im)
            ax[row,col].set_aspect('equal')
            ax[row,col].set(title=f'{times[idx]/1000:g} Gyr',xlabel='x [kpc]',ylabel='y [kpc]')
            if col==0:ax[row,col].text(.03,.96,case.upper(),transform=ax[row,col].transAxes,
                                     va='top',ha='left',color='white',fontsize=9)
    fig.colorbar(ims[-1],ax=ax.ravel().tolist(),shrink=.77,location='right',label='Illustrative r-band [mag arcsec$^{-2}$]')
    fig.suptitle('Smoothed face-on brightness proxy: uniform mass-to-light ratio\n'
                 'Tracer surface density × fixed mass normalization; not observational photometry')
    fig.savefig(out/'brightness_maps_2gyr.png');plt.close(fig)

    # Figure 5: Rd and extrapolated mu0; fit-range sensitivity is shown, no smooth cosmetization.
    fig,axes=plt.subplots(2,2,figsize=(11.3,7.3),layout='constrained',sharex=True)
    for case,ls in (('unbarred','-'),('barred','--')):
        h,y=diag[case];t=y[:,0];lbl=case.title()
        axes[0,0].plot(t,y[:,6],ls,label=lbl)
        axes[0,1].plot(t,y[:,7],ls,label=lbl)
        axes[1,0].plot(t,y[:,8],ls,label=lbl)
        axes[1,1].plot(t,y[:,9],ls,label=lbl)
    for a,l in zip(axes.flat, ['Outer exponential $R_d$ (4–11 kpc) [kpc]',
              'Fit-range check: $R_d$ (6–12 kpc) [kpc]',
              r'Extrapolated outer-disc $\mu_{0,r}$ [AB mag arcsec$^{-2}$]',
              'Median radius $R_{50}$ [kpc]']):
        a.set(xlabel='Science time [Myr]',ylabel=l);a.grid(alpha=.2);a.legend(fontsize=8)
    fig.suptitle('Disc profile evolution in equal-weight tracer population (not self-gravitating)')
    fig.savefig(out/'disc_length_brightness_2gyr.png');plt.close(fig)

    # Figure 6: Guiding shift vs time, and mean Lz (one panel).
    from mwbar.measure import diagnostics as olddiag
    fig,ax=plt.subplots(1,2,figsize=(11,4.5),layout='constrained')
    for case in ('unbarred','barred'):
        _,old,_=olddiag(model,runs[case],case=='barred')
        ax[0].plot(old[:,0],old[:,8],label=case.title())
        ax[1].plot(old[:,0],100*(old[:,6]/old[0,6]-1),label=case.title())
    ax[0].set(ylabel=r'Median $|\Delta R_g|$ [kpc]',xlabel='Science time [Myr]')
    ax[1].set(ylabel='Change in mean $L_z$ [%]',xlabel='Science time [Myr]')
    for a in ax:a.grid(alpha=.2);a.legend(fontsize=8)
    fig.suptitle('Orbital torques and guiding-centre changes, relative to each model\'s science t=0')
    fig.savefig(out/'migration_torque_2gyr.png');plt.close(fig)

    animate_orbits(model,runs,out)
    animate_brightness(model,runs,out)


def _encode_gif(images,out,duration=110):
    if not images:raise ValueError('Empty frames')
    images[0].save(out,save_all=True,append_images=images[1:],duration=duration,loop=0,disposal=2,optimize=False)


def animate_orbits(model,runs,out):
    # 2 Gyr, 25 Myr frame interval. Display is decimated; all 5 Myr samples saved.
    frames=[];tlist=runs['barred']['t'];rate=round(25/(tlist[1]-tlist[0]))
    for i in range(0,len(tlist),rate):
        frame=Image.new('RGB',(840,450),(11,13,22));d=ImageDraw.Draw(frame)
        for case,cx in (('unbarred',210),('barred',630)):
            x=runs[case]['x'][i]
            px=np.rint(cx+x[:,0]*12.5).astype(int)
            py=np.rint(235-x[:,1]*12.5).astype(int)
            good=(px>cx-198)&(px<cx+198)&(py>34)&(py<435)
            d.rectangle((cx-200,33,cx+200,436),outline=(76,82,109),width=1)
            for xv,yv in zip(px[good],py[good]):d.point((int(xv),int(yv)),fill=(245,182,109))
        theta=model.omega_bar*tlist[i]+np.deg2rad(model.bar_angle_deg)
        dx=model.bar_length*12.5*np.cos(theta);dy=model.bar_length*12.5*np.sin(theta)
        d.line((630-dx,235+dy,630+dx,235-dy),fill=(71,225,238),width=3)
        d.text((128,15),'UNBARRED CONTROL',fill=(238,241,246))
        d.text((535,15),'EXISTING BAR',fill=(238,241,246))
        d.text((24,438),f't = {tlist[i]/1000:.3f} Gyr | 10,000 equal-weight stellar tracers | cyan: prescribed bar axis',fill=(206,209,216))
        frames.append(frame.quantize(colors=64,method=Image.Quantize.FASTOCTREE))
    _encode_gif(frames,Path(out)/'orbits_2gyr.gif',duration=110)


def animate_brightness(model,runs,out):
    # 50 Myr per frame; fixed 19.5-28.5 mag/arcsec² mapping for both cases.
    from matplotlib import colormaps
    cmap=colormaps['magma_r']
    images=[];times=runs['barred']['t'];stride=round(50/(times[1]-times[0]))
    for i in range(0,len(times),stride):
        im=Image.new('RGB',(840,440),(11,13,22));d=ImageDraw.Draw(im)
        for case,ox in (('unbarred',10),('barred',430)):
            mu=_smoothed_hist(runs[case]['x'][i],model.disc_mass,bins=90,extent=14.)
            norm=np.clip((mu-19.5)/9.0,0,1);norm=np.nan_to_num(norm,nan=1.)
            rgb=(np.asarray(cmap(norm)[:,:,:3])*255).astype('uint8')
            tile=Image.fromarray(np.flipud(rgb),'RGB').resize((400,400),Image.Resampling.BILINEAR)
            im.paste(tile,(ox,27));d.text((ox+7,9),case.upper(),fill=(242,244,247))
        d.text((20,428),f't = {times[i]/1000:.2f} Gyr  |  r-band proxy, fixed M/L=1.5 Msun/Lsun  |  common color scale',fill=(205,207,214))
        images.append(im.quantize(colors=64,method=Image.Quantize.FASTOCTREE))
    _encode_gif(images,Path(out)/'brightness_2gyr.gif',duration=130)
