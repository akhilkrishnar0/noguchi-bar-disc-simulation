#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=2
mkdir -p agama_stage1/results/live
python3 agama_stage1/initialize.py --out agama_stage1/results --nd 6000 --nh 2000 --nb 2000 --iterations 4
g++ -O3 -std=c++17 -fopenmp noguchi_pilot/nbody.cpp -o agama_stage1/nbody
g++ -O3 -std=c++17 -fopenmp noguchi_pilot/exact_energy.cpp -o agama_stage1/exact_energy
agama_stage1/nbody agama_stage1/results/initial.bin agama_stage1/results/live 1200 0.003333333333333333 0.5 > agama_stage1/results/live/run.log
agama_stage1/exact_energy agama_stage1/results/initial.bin > agama_stage1/results/exact_initial.csv
agama_stage1/exact_energy agama_stage1/results/live/snapshot_1200.bin > agama_stage1/results/exact_final.csv
python3 agama_stage1/analyze.py
