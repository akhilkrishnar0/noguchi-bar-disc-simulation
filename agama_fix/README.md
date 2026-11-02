# Live galaxy failure diagnosis — common 10k-particle initial conditions

This diagnostic isolates the old 300 Myr AGAMA model's **20.38% thickness
increase** and **0.231% exact endpoint energy drift**. Both are measured
failures, not evidence of bar formation. This experiment retains live
disc/halo/bulge self-gravity and imposes no bar.

The same sampled AGAMA distribution-function realization is reused in all runs:

| Configuration | Epsilon scale | Tree theta | dt | Duration |
|---|---:|---:|---:|---:|
| baseline | 1.0 | 0.5 | 0.25 Myr | 300 Myr |
| refined | 1.0 | 0.35 | 0.125 Myr | 300 Myr |
| smaller_eps | 0.5 | 0.35 | 0.125 Myr | 300 Myr |

Saved snapshots have 5 Myr cadence. Integration uses the *smaller* dt
above. B changes both timestep and theta, so those two causes cannot be
isolated separately; C differs only from B in softening.

**Implementation repair:** The earlier tree solver hard-coded species
softenings in its accepted tree-node forces. The modified solver reads
species softening directly from its particle file. This is required to
perform a valid softening comparison; leaf forces already used particle
softenings.

**Force test:** AGAMA smooth-potential accelerations were saved for sampled
disc/halo/bulge particles. The new direct O(N) per sampled target calculation
compares these against the finite-N Plummer pair accelerations with three
softening choices. Force mismatch has contributions from discreteness
noise and changed smoothing; comparison is diagnostic, not a DF solution.

**Live validation:** 300 Myr tests independently measure exact-pair
endpoint energy, total Lz, disc RMS height, R50, A2 and three-epoch
images. Numerical screens: abs energy <0.1%, abs Lz <0.05%,
abs thickness change <5%, abs R50 <2%. These are exploratory criteria,
not universal scientific thresholds.

The pipeline reports physical PASS/FAIL from measurements and publishes
a validation JSON even on a physical failure. Software errors result in a
failed GitHub workflow. Binary test snapshots remain in a temporary workflow
artifact and are not checked into Git.

Reproduction from project root, on Linux with GSL, Eigen, C++17, Python
and AGAMA installed:

    bash agama_fix/run.sh

Expected output files in agama_fix/results include force_comparison.json,
fix_validation.json, three_epochs_comparison.png and
numerical_diagnostics.png.

A genuine physical repair may require generating the equilibrium DF in a
potential consistent with the finite-N softening and demonstrating
particle-number, force-opening, timestep, and multi-seed convergence.
Nothing here establishes a full Noguchi (2001) reproduction.
