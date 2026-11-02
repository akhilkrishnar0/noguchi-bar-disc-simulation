"""Fast solver regression: variable per-species softening and direct-force accuracy.

No AGAMA or long simulations needed. Run with:
python -m unittest discover -s tests -p test_agama_tree.py -v
"""
import csv,subprocess,tempfile,unittest
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
class TestBarnesHutDynamicSoftening(unittest.TestCase):
 @classmethod
 def setUpClass(cls):
  cls.temp=tempfile.TemporaryDirectory()
  cls.folder=Path(cls.temp.name)
  cls.bin=cls.folder/"nbody"
  subprocess.run(["g++","-O2","-std=c++17","-fopenmp",
      str(ROOT/"agama_fix/nbody_v2.cpp"),"-o",str(cls.bin)],check=True)
 @classmethod
 def tearDownClass(cls):cls.temp.cleanup()
 def particles(self,eps):
  rng=np.random.default_rng(4453)
  n=192
  x=rng.normal(size=(n,3))*.4;v=np.zeros_like(x)
  m=np.full(n,1/n)
  c=np.arange(n)%3
  e=np.asarray(eps)[c]
  return np.column_stack((x,v,m,e,c)).astype("<f8")
 def evaluate(self,data,label):
  infile=self.folder/f"{label}.bin";out=self.folder/label
  out.mkdir(exist_ok=True);data.tofile(infile)
  run=subprocess.run([str(self.bin),str(infile),str(out),"0","0.002",".2","1"],
       capture_output=True,text=True)
  return run,out
 def test_multiple_component_softenings(self):
  p=self.particles([.012,.025,.004])
  run,folder=self.evaluate(p,"species")
  self.assertEqual(run.returncode,0,msg=run.stderr)
  lines=list(csv.DictReader((folder/"force_validation.csv").open()))
  errs=np.array([float(v["relative_error"]) for v in lines])
  self.assertLess(np.percentile(errs,95),.03)
  self.assertTrue((folder/"snapshot_0.bin").exists())
 def test_softening_parameter_changes_gravitational_potential(self):
  q=self.particles([.02,.04,.004])
  p=self.particles([.01,.02,.002])
  a,da=self.evaluate(q,"large")
  b,db=self.evaluate(p,"small")
  self.assertEqual((a.returncode,b.returncode),(0,0))
  def energy(d):
   rows=list(csv.DictReader((d/"diagnostics.csv").open()))
   return float(rows[0]["U"])
  self.assertNotAlmostEqual(energy(da),energy(db),places=5)
 def test_reject_nonuniform_species_eps(self):
  p=self.particles([.01,.02,.003]);p[3,7]*=1.5
  run,_=self.evaluate(p,"invalid")
  self.assertNotEqual(run.returncode,0)
if __name__=="__main__":unittest.main()
