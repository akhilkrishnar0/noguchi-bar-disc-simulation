# Validation outcome

**Numerical screening improved; initial equilibrium failed.** These are completed tests of the original particle realization, not a repaired model. All simulations are stopped.

## Controlled tests

| Case | theta | dt | End t | Exact energy drift |
|---|---:|---:|---:|---:|
| Live, tighter forces | 0.5 | 0.0025 | 1.0 | +0.01112% |
| Frozen initial field | 0.5 | 0.0025 | 1.0 | +0.00673% |
| Live, smaller timestep | 0.5 | 0.00125 | 0.25 | -0.00291% |

Live Lz drift at t=1 is 0.00366%. The static particle field is not axisymmetric, so total target Lz is not a conserved quantity in the frozen-field experiment. Frozen-field energy is K+sum(m phi_external), not the live self-gravitating energy definition.

## Initial equilibrium

The initial disc tensor-virial support ratios along x,y,z are [1.012407935179993, 1.0144127631117013, 1.55841787351115]. Equilibrium requires these ratios near unity as a necessary moment condition. The vertical ratio of about 1.56 means kinetic support substantially exceeds the restoring term for the actual softened force law. The approximate initializer does not consistently solve vertical balance in that force field.

Disc RMS height rises from 0.2728 kpc to 0.3254 kpc in the live case (19.28%), and to 0.3224 kpc in the frozen field (18.17%). Similar thickening without collective response supports an initialization imbalance as a major contributor. Frozen-field graininess and approximate tree forces remain caveats; this is not a decomposition of every heating mechanism.

The final live vertical virial ratio is 0.9859, but reaching a balanced moment at one time does not establish a stationary distribution or justify relabeling the evolved model as Noguchi's initial galaxy.

## Force and timestep sensitivity

At t=1, the live force-comparison median, 95th percentile and maximum sampled relative errors are [0.3446715  0.92319115 1.31566   ] percent (192 target particles, direct all-source sum).

At t=0.25, halving dt changes disc RMS height by -0.0073%, half-mass radius by -0.000035%, and A2 by 0.0000275. This supports short-interval timestep consistency only.

At t=1, the theta=0.8 and theta=0.5 half-mass radii differ by 1.26%. The tighter run contracts slightly from 5.7626 to 5.7458 kpc, whereas the earlier run expanded slightly. Therefore the earlier ~1% expansion is not a converged physical result. A2 is much less sensitive over this interval.

## Collective response

Global A2 reaches 0.03074 in the live case and 0.00629 in the frozen case. This is consistent with self-gravitational amplification of a weak response, but not proof of a bar. Radial phase coherence, length, pattern speed and persistence were not established. Disequilibrium can itself seed collective transients.

The classical Q diagnostic starts between 1.417 and 1.551 over the sampled annuli, close to the 1.5 target. Final Q spans 1.384–1.853. This uses measured dispersions and density with an epicyclic frequency inferred from direct azimuthally averaged softened forces; it is not an exact finite-thickness/softened-gravity stability criterion.

## Decision

The tighter live run passes the chosen short-run energy (<0.1%) and total Lz (<0.01%) screening tolerances, and the half-mass-radius stationarity screen (<2%). It **fails** the vertical stationarity screen (<5%). These tolerances are practical checks, not universal standards. Component angular-momentum exchange remains more demanding than global conservation.

Do not extend this initialization to t=60 and interpret its evolution as a Noguchi reproduction. First correct vertical equilibrium, reconstruct the halo and preparation sequence, and validate density/kinematics against Figure 1. Then test several orbital periods and independent seeds before pursuing a physical bar-growth result.
