import unittest
import numpy as np
from mwbar import Galaxy,Settings,integrate_pair
from mwbar.integrate import grown_bar_amplitude
from mwbar.measure import guiding_radius,diagnostics

class TestModel(unittest.TestCase):
    def setUp(self):self.m=Galaxy()
    def test_galactic_solar_normalization(self):
        from mwbar.model import KMS_TO_KPCMYR
        self.assertAlmostEqual(self.m.circular_speed([8.2])[0]/KMS_TO_KPCMYR,232.8,places=5)
        self.assertGreater(self.m.halo_mass,0)
    def test_analytic_gravity_matches_potential_gradient(self):
        rng=np.random.default_rng(991)
        pts=rng.uniform(-10,10,size=(15,3));pts[:,2]*=.2
        eps=1e-4
        for mode in (False,True):
            for t in (0.,78.3):
                acceleration,_=self.m.force(pts,t,mode)
                for k in range(3):
                    h=np.zeros_like(pts);h[:,k]=eps
                    potential_plus=self.m.force(pts+h,t,mode)[1]
                    potential_minus=self.m.force(pts-h,t,mode)[1]
                    numeric=-(potential_plus-potential_minus)/(2*eps)
                    np.testing.assert_allclose(acceleration[:,k],numeric,rtol=1e-6,atol=1e-8)
    def test_bar_zero_azimuthal_average_and_strength(self):
        ring=np.linspace(0,2*np.pi,1000,endpoint=False)
        xyz=np.column_stack((4*np.cos(ring),4*np.sin(ring),0*ring))
        force,phi=self.m.bar(xyz,0)
        self.assertLess(abs(phi.mean()),1e-12)
        self.assertLess(abs(force[:,0].mean()),1e-12)
        self.assertGreater(self.m.bar_amplitude,0)
    def test_preparation_ramp_smooth_and_saturated(self):
        s=Settings()
        self.assertEqual(grown_bar_amplitude(-s.warmup_myr,s),0)
        self.assertEqual(grown_bar_amplitude(0.,s),1)
        self.assertEqual(grown_bar_amplitude(-s.warmup_myr+s.growth_myr,s),1)
        self.assertAlmostEqual(grown_bar_amplitude(-s.warmup_myr+s.growth_myr/2,s),.5)
    def test_resonances_ordered(self):
        r=self.m.resonance_radii()
        self.assertTrue(r['ILR'] and r['CR'] and r['OLR'])
        self.assertLess(r['ILR'][0],r['CR'][0])
        self.assertLess(r['CR'][0],r['OLR'][0])
    def test_orbit_numerics_and_cadence(self):
        s=Settings(n=600,warmup_myr=80,growth_myr=45,science_myr=50,
                   snapshot_myr=5,timestep_myr=.25)
        runs=integrate_pair(self.m,s)
        for case in ('unbarred','barred'):
            self.assertEqual(runs[case]['x'].shape,(11,600,3))
            np.testing.assert_array_equal(runs[case]['t'],np.arange(0,51,5))
            _,_,r=diagnostics(self.m,runs[case],case=='barred')
            self.assertLess(r['max_relative_mean_integral_deviation'],1e-4)
            self.assertTrue(np.isfinite(r['Rd_final_kpc']))
        self.assertGreater(runs['barred']['phase0_x'].shape[0],0)
    def test_timestep_halving_converges(self):
        # Same seed, same phase/time and particle identities in paired setups.
        a=Settings(n=300,seed=44,timestep_myr=.5,warmup_myr=60,growth_myr=40,
                   science_myr=30,snapshot_myr=5)
        b=Settings(n=300,seed=44,timestep_myr=.25,warmup_myr=60,growth_myr=40,
                   science_myr=30,snapshot_myr=5)
        ra=integrate_pair(self.m,a);rb=integrate_pair(self.m,b)
        for mode in ('unbarred','barred'):
            xa=ra[mode]['x'][-1]
            xb=rb[mode]['x'][-1]
            rms=float(np.sqrt(np.mean(np.sum((xa-xb)**2,axis=1))))
            self.assertLess(rms,.09)
if __name__=='__main__':unittest.main()
