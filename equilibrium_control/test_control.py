"""Regression tests for numerical fixed-potential control."""
import unittest
import numpy as np
from control import gravity,run
class TestControl(unittest.TestCase):
 def test_force_equals_minus_potential_gradient(self):
  xyz=np.array([[.31,.22,.035],[.08,-.31,.025],[.63,.43,-.03]])
  for forced,t in [(False,1),(True,10)]:
   acc=gravity(xyz,t,forced)
   for k in range(3):
    delta=np.zeros_like(xyz);delta[:,k]=1e-6
    plus=gravity(xyz+delta,t,forced,True)[1]
    minus=gravity(xyz-delta,t,forced,True)[1]
    np.testing.assert_allclose(acc[:,k],-(plus-minus)/2e-6,rtol=2e-8,atol=2e-8)
 def test_axisymmetric_Lz(self):
  r=run(n=1000,seed=2001,duration_myr=300,step_myr=1,output_myr=5)
  self.assertLess(abs(r['Lz_change_pct']),1e-7)
 def test_4500myr_stationarity(self):
  r=run(n=3000,seed=2001,step_myr=1,duration_myr=4500,output_myr=5)
  self.assertLess(abs(r['rms_height_change_pct']),5.)
  self.assertLess(abs(r['R50_change_pct']),2.)
  self.assertLess(abs(r['energy_change_pct']),.01)
  self.assertGreater(r['min_A2_count'],2000)
 def test_rotating_bar_jacobi(self):
  r=run(n=3000,seed=2001,step_myr=1,duration_myr=4500,output_myr=5,forced=True)
  self.assertLess(r['jacobi_post_ramp_max_fraction'],2e-4)
 def test_dt_convergence(self):
  a=run(n=3000,seed=2001,step_myr=1,duration_myr=4500,output_myr=5)
  b=run(n=3000,seed=2001,step_myr=.5,duration_myr=4500,output_myr=5)
  self.assertLess(abs(a['rms_height_change_pct']-b['rms_height_change_pct']),1)
  self.assertLess(abs(a['R50_change_pct']-b['R50_change_pct']),.5)
  self.assertLess(abs(b['energy_change_pct']),abs(a['energy_change_pct'])+.001)
if __name__=='__main__':
 unittest.main(verbosity=2)
