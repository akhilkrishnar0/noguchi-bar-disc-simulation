# Validation design

Goal: test the *existing* approximate Model A realization, without changing its initial particles or claiming a reconstruction of Noguchi's equilibrium preparation.

Full particle count: 59,985. No damping, random reseeding, resizing or forced bar. The initial binary state is byte-identical to the first pilot. Pairwise Plummer softening and component masses are unchanged.

1. Live evolution to t=1 with tree opening angle 0.5 and dt=0.0025, compared with the existing theta=0.8 run. This isolates tree accuracy.
2. Live evolution to t=0.25 with theta=0.5 and dt=0.00125. Compare with test 1 at the same time. This is a shorter timestep test; no full-duration timestep convergence is claimed.
3. Frozen initial particle-source field to t=1, theta=0.5, dt=0.0025. All target particles are integrated, but source positions and masses remain fixed. Each target's original source particle is excluded, consistently in tree and direct evaluation. This is a deliberately different, externally imposed field: no back-reaction and no collective growth. It tests orbital consistency of the initial distribution. The finite-particle source field is not perfectly axisymmetric or smooth, so frozen-field evolution can include static-potential scattering and imposed-noise response.

Both live cases use KDK leapfrog and rebuilt tree forces each step; frozen-field sources use a fixed tree. An accepted monopole node is forbidden from containing the excluded source even after a target has moved. Component monopoles retain their own softening lengths. All initial positions and velocities are unchanged.

Diagnostics:

- Direct all-source acceleration comparison for the same 192 target indices at initial and final epochs.
- Exact all-pair endpoint energies. Live total energy is K + 1/2 sum(m phi); frozen-field orbital energy is K + sum(m phi_external). Frozen-field total Lz need not be conserved because its imposed particle field is not axisymmetric.
- Tensor virial support ratios sum(m v_i^2) / [-sum(m x_i a_i)] for each component; stationary equilibrium requires these ratios to be near 1, though this necessary moment condition is not sufficient to prove equilibrium. Velocities/positions are measured relative to total bulk/COM for this diagnostic.
- Disc RMS height, cylindrical half-mass radius, outer-edge population, A2, and component angular momenta. Spherical half-mass radii and inner-90%-particle inertia shape ratios for all components.
- Direct-force azimuthal-ring probes at z=0 to estimate circular velocity and epicyclic frequency. Measured radial dispersions and surface densities give a classical Q diagnostic. Finite thickness, force softening and a nonaxisymmetric live state mean this is not an exact stability threshold.

Practical screening targets for this short pilot: exact endpoint energy error <0.1%; total live Lz error <0.01%; disc RMS-height change <5% and half-mass-radius change <2%. These are chosen diagnostic tolerances, not universal physical thresholds. Component torque interpretation additionally requires numerical drift well below the inferred exchange. Passing this short screen cannot establish stability over t=60. Failing equilibrium does not invalidate the gravity law; it invalidates use of these initial conditions as a stationary starting galaxy.

The original paper's 15.6-time-unit halo relaxation and later halo-bulge preparation are NOT performed in these validation tests. No claim is made to match Figure 1. A failed test will be reported as a failed test, not relabeled as physical galaxy evolution.
