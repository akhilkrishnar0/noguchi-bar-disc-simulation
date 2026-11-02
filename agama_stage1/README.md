# Stage 1: 300 Myr live N-body equilibrium test

**Pending live validation; do not claim a successful result before the measurements.**

This experiment samples iterative self-consistent action-based distribution
functions with [AGAMA](https://github.com/GalacticDynamics-Oxford/Agama):
6,000 live disc, 2,000 live dark halo, and 2,000 live bulge particles.
It evolves 300 Myr without imposing a bar, with an internal integration
step of 0.25 Myr and snapshots every 5 Myr. The official AGAMA
example_self_consistent_model_bdh.py is the method reference.

**This is NOT an exact Noguchi (2001) Model A reproduction.**
The 0.6/0.2/0.2 component masses, disc scalelength 0.25 and bulge
Re=0.05 are targets, but the halo and truncation are independently assumed.
AGAMA uses a smooth potential for the DF, while the live N-body solver
uses softened gravity with epsilon=0.02,0.04,0.004 for the disc, halo and bulge.
This potential mismatch is a limitation that must be tested.

Outputs: exact endpoint pairwise energy, Lz drift, disc RMS thickness,
R50, radial surface density, Q profile, and m=2 structure with 0, 150
and 300 Myr figures. A machine-readable validation.json reports PASS/FAIL.
The exploratory gates are under 5% change in disc height, under 2% in R50,
under 0.1% exact energy error, under 0.05% Lz change, and sampled median Q
between 1 and 2.5 (nominal target Q=1.5).

## Run

Requires Linux, GCC, GSL, Eigen, numpy, scipy, matplotlib and AGAMA.

    sudo apt-get install -y libgsl-dev libeigen3-dev
    python -m pip install -r requirements.txt setuptools wheel
    python -m pip install agama --no-build-isolation
    bash agama_stage1/run.sh

GitHub Actions executes the actual experiment and commits only *measured*
outputs. Failures and unfinished tests are not presented as successes.

If this short screen passes, the next work is force-softening, timestep,
particle-number and seed convergence, then multi-orbit stationarity and
spontaneous bar dynamics. No conclusions about GLSB formation can be made
from a 300 Myr equilibrium screen.

Reference: Noguchi (2001), MNRAS 328, 353.
https://doi.org/10.1046/j.1365-8711.2001.04910.x
