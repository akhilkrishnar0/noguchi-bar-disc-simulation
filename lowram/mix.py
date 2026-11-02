"""Controlled disc azimuthal phase mixing to remove numerical m=2 seeds after warm-up.

Rotates each disc particle's xy position and xy velocity by the same random
angle; preserves R,z,v_R,v_phi,v_z and live spheroid particles. Cannot produce
an exact distribution function or guarantee equilibrium.
"""
from pathlib import Path
import argparse,numpy as np
pa=argparse.ArgumentParser();pa.add_argument('input');pa.add_argument('output');pa.add_argument('--seed',type=int,default=2001);args=pa.parse_args()
a=np.fromfile(args.input,dtype='<f8').reshape(-1,9);d=a[:,8]==0;phi=np.random.default_rng(args.seed).uniform(0,2*np.pi,d.sum());c=np.cos(phi);s=np.sin(phi)
for i,j in [(0,1),(3,4)]:
 x=a[d,i].copy();y=a[d,j].copy();a[d,i]=x*c-y*s;a[d,j]=x*s+y*c
Path(args.output).parent.mkdir(parents=True,exist_ok=True);a.tofile(args.output)
print('axisymmetrized disc:',int(d.sum()),'particles')
