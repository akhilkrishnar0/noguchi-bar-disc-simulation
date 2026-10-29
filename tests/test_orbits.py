"""Small physical regression tests, independent of galaxy initialization."""
from pathlib import Path
import os,subprocess,tempfile,unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
class OrbitTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):
  cls.tmp=tempfile.TemporaryDirectory();cls.base=Path(cls.tmp.name);cls.exe=cls.base/'solver'
  subprocess.run(['g++','-O2','-fopenmp','-std=c++17',str(ROOT/'validation_run/nbody_validate.cpp'),'-o',str(cls.exe)],check=True)
 @classmethod
 def tearDownClass(cls):cls.tmp.cleanup()
 def run_case(self,frozen):
  folder=self.base/str(frozen);folder.mkdir();p=np.zeros((2,9));p[:,0]=[-.5,.5];p[:,4]=[-1,1]*np.array([.5,.5])/(1+.02**2)**.75;p[:,6]=.5;p[:,7]=.02
  initial=folder/'input.bin';p.tofile(initial)
  subprocess.run([str(self.exe),str(initial),str(folder),'1000','.001','.5',str(frozen)],check=True,stdout=subprocess.DEVNULL,env={**os.environ,'OMP_NUM_THREADS':'1'})
  d=np.genfromtxt(folder/'diagnostics.csv',delimiter=',',names=True)
  self.assertTrue(np.isfinite(d['E']).all())
  self.assertLess(np.max(abs(d['E']-d['E'][0]))/abs(d['E'][0]),1e-6)
  return d,folder
 def test_softened_circular_binary(self):
  d,folder=self.run_case(0)
  self.assertLess(np.max(abs(d['Lz']/d['Lz'][0]-1)),1e-10)
  p=np.fromfile(folder/'snapshot_1000.bin').reshape(-1,9)
  self.assertLess(abs(np.linalg.norm(p[1,:3]-p[0,:3])-1),1e-6)
 def test_frozen_external_energy(self):
  self.run_case(1)
if __name__=='__main__':unittest.main()
