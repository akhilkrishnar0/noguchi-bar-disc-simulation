import unittest
import numpy as np
from mwbar.model import Galaxy
from mwbar.integrate import Settings,integrate_pair
from v3.measure import mu_r_from_sigma,fourier_radial,response_extent,outer_fit,diagnostics_v3

class V3DiagnosticsTests(unittest.TestCase):
    def test_surface_brightness_consistency(self):
        # Sigma=1.5 Msun/pc2 / ML=1.5 => 1 Lsun/pc2.
        self.assertAlmostEqual(float(mu_r_from_sigma(1.5)),4.65+21.572,places=8)
        self.assertTrue(np.isnan(mu_r_from_sigma(0)))
    def test_fourier_coherent_rotating_ellipse(self):
        rng=np.random.default_rng(45);model=Galaxy();t=200
        phi=model.omega_bar*t
        a=rng.uniform(0,2*np.pi,30000)
        r=np.sqrt(rng.uniform(0.6,5.4,30000))
        # Impose coherent elongation about the bar axis.
        x=np.column_stack((r*1.4*np.cos(a),r*.70*np.sin(a),np.zeros(len(r))))
        c,s=np.cos(phi),np.sin(phi)
        x[:,:2]=np.column_stack((c*x[:,0]-s*x[:,1],s*x[:,0]+c*x[:,1]))
        center,amp,delta,count=fourier_radial(x,model,t)
        sel=(count>150)&(center>=1)&(center<4)
        self.assertGreater(np.median(amp[sel]),.2)
        self.assertLess(np.median(np.abs(delta[sel])),5)
    def test_axial_symmetry_no_false_bar(self):
        model=Galaxy(); rng=np.random.default_rng(1234);n=30000
        R=rng.gamma(shape=2,scale=2.6,size=n);ang=rng.uniform(0,2*np.pi,n)
        x=np.column_stack((R*np.cos(ang),R*np.sin(ang),np.zeros(n)))
        extent,amax,_=response_extent(x,model,0)
        self.assertTrue(np.isnan(extent))
    def test_sample_short_numerics_finite(self):
        model=Galaxy();s=Settings(n=300,seed=55,timestep_myr=.5,warmup_myr=50,growth_myr=30,
                                   science_myr=30,snapshot_myr=5)
        runs=integrate_pair(model,s)
        d,pro=diagnostics_v3(model,runs)
        for case in ('barred','unbarred'):
            h,a=d[case];self.assertEqual(len(a),7)
            self.assertTrue(np.isfinite(a[:,1]).all())
            self.assertTrue(np.isfinite(pro[case][1][:,2]).all())
if __name__=='__main__':unittest.main()
