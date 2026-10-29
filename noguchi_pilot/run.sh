#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}"
mkdir -p baseline half_step tight_force
for d in baseline half_step tight_force; do
  find "$d" -maxdepth 1 -type f -name "snapshot_*.bin" -delete
done
python3 initialize.py
g++ -O3 -fopenmp -std=c++17 nbody.cpp -o nbody
g++ -O3 -fopenmp exact_energy.cpp -o exact_energy
mkdir -p baseline half_step tight_force
./nbody initial.bin baseline 200 0.005 0.8
./nbody initial.bin half_step 400 0.0025 0.8
./nbody initial.bin tight_force 0 0.005 0.5
./exact_energy initial.bin > exact_initial.csv
./exact_energy baseline/snapshot_200.bin > exact_baseline.csv
./exact_energy half_step/snapshot_400.bin > exact_half_step.csv
python3 analyze.py
