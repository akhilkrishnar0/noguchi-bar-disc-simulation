#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}"
# The pilot run creates the initial realization and comparison histories.
if [[ ! -f noguchi_pilot/half_step/snapshot_400.bin ]]; then
  bash noguchi_pilot/run.sh
fi
for d in live_tight frozen_field small_step; do
  mkdir -p "validation_run/$d"
  find "validation_run/$d" -maxdepth 1 -type f \( -name 'snapshot_*.bin' -o -name 'acceleration_*.bin' \) -delete
done
g++ -O3 -fopenmp -std=c++17 validation_run/nbody_validate.cpp -o validation_run/nbody_validate
g++ -O3 -fopenmp validation_run/probe.cpp -o validation_run/probe
g++ -O3 -fopenmp validation_run/exact_frozen_energy.cpp -o validation_run/exact_frozen_energy
g++ -O3 -fopenmp noguchi_pilot/exact_energy.cpp -o noguchi_pilot/exact_energy
validation_run/nbody_validate noguchi_pilot/initial.bin validation_run/live_tight 400 0.0025 0.5 0 > validation_run/live_tight/run.log
validation_run/nbody_validate noguchi_pilot/initial.bin validation_run/frozen_field 400 0.0025 0.5 1 > validation_run/frozen_field/run.log
validation_run/nbody_validate noguchi_pilot/initial.bin validation_run/small_step 200 0.00125 0.5 0 > validation_run/small_step/run.log
noguchi_pilot/exact_energy noguchi_pilot/initial.bin > noguchi_pilot/exact_initial.csv
noguchi_pilot/exact_energy validation_run/live_tight/snapshot_400.bin > validation_run/live_tight/exact_final.csv
noguchi_pilot/exact_energy validation_run/small_step/snapshot_200.bin > validation_run/small_step/exact_final.csv
validation_run/exact_frozen_energy validation_run/frozen_field/snapshot_400.bin noguchi_pilot/initial.bin > validation_run/frozen_field/exact_final.csv
validation_run/probe noguchi_pilot/initial.bin > validation_run/probe_initial.csv
validation_run/probe validation_run/live_tight/snapshot_400.bin > validation_run/probe_final.csv
python3 validation_run/analyze_validation.py
