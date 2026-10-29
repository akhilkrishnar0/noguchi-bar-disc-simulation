# First collisionless galaxy pilot: measured results

This is an executed, short, three-dimensional self-gravitating N-body experiment inspired by Noguchi (2001), Model A. It is **not an exact reproduction** and is **not yet a validated equilibrium model for long-term bar science**. No bar potential, gas, star formation, external torque or drag is imposed.

## What ran

59,985 particles: 30,000 disc, 10,000 halo, 19,985 bulge. G=1; component masses are 0.6, 0.2, 0.2. Disc scale length 0.25, bulge projected effective radius 0.05, component softening lengths 0.02, 0.04, 0.004. Units for displayed figures: 15 kpc, 1.5e11 solar masses, and approximately 75 Myr, following the example scaling in Noguchi (2001). The rounded time conversion is used only for display.

Two integrations use the same particle realization and tree opening angle 0.8:

- Baseline: timestep 0.005, 200 steps, t=1.
- Half-step: timestep 0.0025, 400 steps, t=1. This run supplies the displayed snapshots and animation.

Wall times printed by the integrator were approximately 86 and 155 seconds respectively, including its initial direct-force check. These are this session's timings, not portable performance guarantees. Both runs finished. No simulation continues in the background.

## Measured results

| Quantity | Initial | Half-step final |
|---|---:|---:|
| Disc half-mass cylindrical radius | 5.7626 kpc | 5.8180 kpc |
| Disc RMS vertical height | 0.2728 kpc | 0.3301 kpc |
| Global A2, 0.75 < R < 9 kpc | 0.00288 | 0.03002 |

Disc half-mass radius increases about 0.96%; RMS thickness increases about 21%. The surface-density profile remains broadly similar in the main disc, with spreading across its initially sharp edge. There is a weak growing twofold asymmetry, but no demonstrated mature bar. An A2 increase alone cannot distinguish a bar from spirals, transient adjustment or amplified sampling noise. A proper bar identification requires radial phase coherence, shape and persistence.

The short duration is 1/60 of the paper's illustrated t=60 endpoint. These results do not establish giant-LSB formation or reproduce its claimed disc-profile changes. Vertical thickening and edge spreading flag initial-equilibrium limitations; they should not be described as bar-driven redistribution.

## Numerical checks

Exact pairwise Plummer potential energies were calculated separately at the endpoints, summing every unique particle pair, without tree approximation. Kinetic energy and all momenta use all particles.

| Check | Baseline | Half-step |
|---|---:|---:|
| Exact total-energy change / absolute initial energy | +1.12865% | +0.03190% |
| Lz change / initial Lz | +0.01505% | +0.01464% |

The baseline timestep is inadequate for precision energy conservation in this reconstructed compact bulge. Halving it strongly improves this diagnostic. This does not imply the original paper was inaccurate: its particle realization, prepared bulge, force softening and integrator implementation may differ.

At t=0, tree accelerations were compared with all-particle direct summation for 192 spatial samples (particle indices spaced through all components):

| Opening angle | Median relative force error | 95th percentile | Maximum sampled |
|---|---:|---:|---:|
| 0.8 | 1.150% | 2.448% | 3.365% |
| 0.5 | 0.387% | 0.923% | 1.526% |

The 0.5 test is a force evaluation only, not a third evolutionary run. The angle-0.8 integration is adequate for this exploratory numerical check, but long-term evolution has not passed a force-convergence study. Small energy drift alone does not prove physical correctness. Octree approximations do not enforce exact pairwise force symmetry, so momentum conservation is measured rather than presumed. Energy curves in the diagnostic figure use tree-estimated potential; the endpoint table above uses the exact pair sum.

## Explicit reconstruction assumptions

1. **Halo**: a Hernquist-shaped density with scale 0.2, tapered smoothly from r=0.8 to r=1 and normalized to mass 0.2. This is an assumed halo, not a fit to Figure 1 and not a recovered published halo law.
2. **Bulge**: numerical Abel deprojection of a Sersic n=4 surface-density profile, Re=0.05, sampled on a finite radial grid to r=4 and renormalized. This follows the stated projected profile, rather than silently substituting a Hernquist bulge.
3. **Disc**: exponential radius distribution truncated at r=1, with assumed sech-squared vertical profile and z0=0.02. The sharp outer truncation is not an equilibrium distribution-function boundary and generates edge adjustment.
4. **Velocities**: Gaussian sampling of approximate Jeans moments. Spheroids use isotropic spherical moments in the combined spherical-average potential; the flattened disc is approximated by its monopole there. Disc radial dispersion targets Q=1.5 using a smooth force model, azimuthal dispersion follows epicyclic theory and mean streaming includes asymmetric drift. Vertical dispersion is a local isothermal-sheet plus external-harmonic approximation. Finite thickness in the midplane disc-force calculation is represented approximately by extra smoothing. This is not an exact distribution-function solution.
5. **Core**: negative mean-streaming-squared values in three very small-radius interpolation grid cells are clipped to zero. Epicyclic theory is not valid arbitrarily close to the origin. Gaussian spheroid velocities are not escape-truncated, so some initially unbound particles are possible. These choices require further equilibrium work before long-term use.
6. **Preparation**: the paper's staged halo/bulge relaxation was not reproduced. This pilot starts directly from approximate moment-initialized components. There is no claim that the model matches Figure 1.
7. **Force law**: symmetric Plummer pair softening, epsilon_ij²=(epsilon_i²+epsilon_j²)/2, is a documented assumption. Monopoles are accumulated separately for each component in accepted tree nodes. G=1. All components interact gravitationally.
8. **Integration**: kick-drift-kick leapfrog, rebuilt octree at every force evaluation, double precision. Random seed 2001 identifies this new realization, not Noguchi's original seed.

## Files and reproduction

Install Python packages from requirements.txt; a C++17 compiler with OpenMP is required. On Linux, run `bash run.sh`. This recreates initial conditions, both integrations, the force test, exact energy checks and figures. Set OMP_NUM_THREADS to choose CPU threads. The model uses no periodic boundaries.

- initialize.py: documented initial-condition generator.
- nbody.cpp: octree forces and leapfrog integration.
- exact_energy.cpp: independent direct-pair endpoint energy calculation.
- analyze.py: diagnostics, fixed-colour-scale snapshots, GIF and results.json.
- results.json: measured numerical results.
- initialization.json and target_profiles.npz: initialization assumptions and smooth disc targets.
- initial.bin and initial.npz: full initial particle state.
- baseline/snapshot_200.bin and half_step/snapshot_400.bin: full endpoint states.
- CSV files: conservation histories, force comparison and exact endpoint energy.
- first_run_snapshots.png, first_run_diagnostics.png, first_run.gif: visual results.

Intermediate snapshots are recreated by run.sh; they are omitted from the archive to reduce download size. Existing output directories must not contain snapshots from a different run length before rerunning analysis.

Binary particle files contain little-endian float64 rows with nine columns: x,y,z,vx,vy,vz,mass,softening,component. Components 0,1,2 are disc, halo, bulge. Read with `np.fromfile(path,dtype='<f8').reshape(-1,9)`. Snapshot numbering is timestep count, not physical time. Analysis measures positions about the total centre of mass. Maps show disc particles only; the halo and bulge are dynamically active but omitted from the rendering. Density maps are smoothed by one image pixel for display; this does not affect the simulation. Edge-on panels use an expanded vertical aspect ratio. Colours show mass surface density, not observed brightness.

## Scientific next step

Refine the equilibrium initialization and verify radial Q, density, rotation and vertical profiles before a longer run. Reproduce or explicitly replace the published staged relaxation, test force convergence over time and additional seeds, then measure coherent bar growth and angular-momentum exchange. Only after that should the t~60 surface-density evolution be compared with Noguchi's figures.

References: Noguchi (2001), MNRAS 328, 353–358, DOI 10.1046/j.1365-8711.2001.04910.x (user-provided PDF); Noguchi (1991), MNRAS 251, 360; Hernquist (1993), ApJS 86, 389. The latter papers informed interpretation of the initialization; this implementation's assumptions are stated above.
