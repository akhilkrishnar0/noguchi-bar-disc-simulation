#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=2
mkdir -p agama_fix/results
python3 agama_fix/initialize_fixed.py --out agama_fix/results --nd 6000 --nh 2000 --nb 2000 --iterations 4
python3 agama_fix/compare_forces.py
g++ -O3 -fopenmp -std=c++17 agama_fix/nbody_v2.cpp -o agama_fix/nbody
g++ -O3 -fopenmp -std=c++17 noguchi_pilot/exact_energy.cpp -o agama_fix/exact_energy
# Three controlled experiments, same galaxy realization, no imposed bar.
# A original epsilon / dt=0.25 Myr / theta=0.5
# B original epsilon / dt=0.125 Myr / theta=0.35
# C half epsilon     / dt=0.125 Myr / theta=0.35
run_case(){
  name="$1"; input="$2"; steps="$3"; dt="$4"; theta="$5"; cadence="$6"
  mkdir -p "agama_fix/results/$name"
  echo "Starting $name" >&2
  agama_fix/nbody "$input" "agama_fix/results/$name" "$steps" "$dt" "$theta" "$cadence" >"agama_fix/results/$name/solver.log"
  agama_fix/exact_energy "$input" >"agama_fix/results/$name/energy_initial.csv"
  agama_fix/exact_energy "agama_fix/results/$name/snapshot_$steps.bin" >"agama_fix/results/$name/energy_final.csv"
}
run_case baseline agama_fix/results/initial_eps_1p0.bin 1200 0.003333333333333333 0.5 20
run_case refined agama_fix/results/initial_eps_1p0.bin 2400 0.0016666666666666668 0.35 40
run_case smaller_eps agama_fix/results/initial_eps_0p5.bin 2400 0.0016666666666666668 0.35 40
python3 agama_fix/analyze.py
