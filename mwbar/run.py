"""End-to-end CLI: two 300 Myr Milky Way tracer experiments and real figures."""
import argparse,json
from dataclasses import asdict
from pathlib import Path
import numpy as np
from .model import Galaxy,KMS_TO_KPCMYR
from .integrate import Settings,integrate_pair
from .measure import diagnostics,profile
from .plots import figures

def run(n=10000,dt=.25,end=300.,warmup=500.,growth=300.,snapshot=5.,seed=2026,
        output='results',with_plots=True):
    model=Galaxy()
    config=Settings(n=n,seed=seed,timestep_myr=dt,warmup_myr=warmup,
                    growth_myr=growth,science_myr=end,snapshot_myr=snapshot)
    if not config.validate():raise ValueError('Invalid integration schedule')
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    runs=integrate_pair(model,config,progress=lambda c,m:print(f'{c}: saved {m} science frames',flush=True))
    measurement={}
    for case,r in runs.items():
        header,data,summary=diagnostics(model,r,case=='barred')
        measurement[case]=(header,data,summary)
        np.savetxt(output/(case+'_diagnostics.csv'),data,delimiter=',',
                   header=header,comments='',fmt='%.10g')
        rows=[]
        for t,x in zip(r['t'],r['x']):
            radii,sigma,counts=profile(model,x)
            rows.extend((t,ri,si,ni) for ri,si,ni in zip(radii,sigma,counts))
        np.savetxt(output/(case+'_density_profiles.csv'),rows,delimiter=',',
                   header='time_myr,radius_kpc,surface_density_Msun_pc2,tracer_count',comments='',fmt='%.9g')
        np.savez_compressed(output/(case+'_orbits.npz'),x=r['x'],v=r['v'],t=r['t'])
    if with_plots:figures(model,runs,measurement,output)
    u=measurement['unbarred'][2];b=measurement['barred'][2]
    gates={'unbarred_energy_rel_error_below_1e-4':u['max_relative_mean_integral_deviation']<1e-4,
           'unbarred_Lz_pct_below_1e-3':abs(u['mean_Lz_change_pct'])<.001,
           'barred_jacobi_rel_error_below_1e-4':b['max_relative_mean_integral_deviation']<1e-4,
           'control_R50_pct_below_10':abs(u['R50_change_pct'])<10,
           'control_zrms_pct_below_15':abs(u['zrms_change_pct'])<15}
    result={
      'status':'PASS_NUMERICAL_ORBIT_CONTROL_ONLY' if all(gates.values()) else 'FAILED_NUMERICAL_ORBIT_CONTROL',
      'not_self_gravitating':True,
      'bar_status':'imposed existing m=2 potential, constant strength/pattern speed at science t=0',
      'warmup_bar_ramp_myr':growth,'pre_science_settle_myr':warmup-growth,
      'matched_initial_particles_at_time_myr':-warmup,
      'different_phase_space_at_science_t0':True,
      'unit_system':{'length':'kpc','time':'Myr','mass':'Msun','velocity':'kpc/Myr'},
      'config':asdict(config),'galaxy_model':asdict(model),
      'derived':{'halo_mass_Msun':model.halo_mass,
                 'Vc_sun_kms':model.circular_speed([model.solar_radius])[0]/KMS_TO_KPCMYR,
                 'bar_period_myr':2*np.pi/model.omega_bar,
                 'target_max_bar_tangential_force_fraction':model.bar_force_fraction,
                 'resonance_radii_kpc':model.resonance_radii()},
      'control':u,'barred':b,
      'checks':{k:bool(v) for k,v in gates.items()},
      'scientific_limitations':[
         'Collisionless massless tracer particles; disc, bulge and halo gravitational fields are fixed and analytic',
         'The prescribed quadrupole bar is not a self-consistent 3D Ferrers bar or bar formation model',
         'No gas, star formation, live halo resonant absorption, or bar slowdown',
         'Initial epicycle DF is approximate; preparatory integration suppresses but does not prove absence of transients',
         'Tracer-weighted density is a proxy, not measured self-consistent galaxy mass redistribution',
         'Guiding radii defined relative to axisymmetric background, approximate in strongly barred regions',
         'Fitted outer Rd depends on specified radial range (4–11 kpc)',
         'Barred and unbarred t=0 states differ after matched pre-science preparation; match is at -warmup',
         'Short 300 Myr evolution is for orbital response, not many-Gyr structural change']}
    (output/'summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:result[k] for k in ('status','derived','control','barred','checks')},indent=2),flush=True)
    if result['status'].startswith('FAILED'):
        raise RuntimeError('Scientific numerical control failed; see results/summary.json')
    return result

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--n',type=int,default=10000)
    parser.add_argument('--dt',type=float,default=.25)
    parser.add_argument('--end',type=float,default=300.)
    parser.add_argument('--warmup',type=float,default=500.)
    parser.add_argument('--growth',type=float,default=300.)
    parser.add_argument('--snapshot',type=float,default=5.)
    parser.add_argument('--seed',type=int,default=2026)
    parser.add_argument('--out',default='results')
    parser.add_argument('--skip-plots',action='store_true')
    args=parser.parse_args()
    run(args.n,args.dt,args.end,args.warmup,args.growth,args.snapshot,
        args.seed,args.out,not args.skip_plots)
if __name__=='__main__':main()
