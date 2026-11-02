#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
export OMP_NUM_THREADS=2
export OPENBLAS_NUM_THREADS=1
python3 prepare.py --output seed.bin
g++ -O3 -std=c++17 -fopenmp nbody_long.cpp -o nbody_long
./nbody_long seed.bin stage1 halo 1170 0.013333333333333334 0.65 1170 0
./nbody_long stage1/final.bin stage2 spheroids 300 0.013333333333333334 0.65 300 0
./nbody_long stage2/final.bin warm1 live 150 0.013333333333333334 0.65 150 1
python3 mix.py warm1/final.bin mix1.bin --seed 2002
./nbody_long mix1.bin warm2 live 150 0.013333333333333334 0.65 150 1
python3 mix.py warm2/final.bin mix2.bin --seed 2003
./nbody_long mix2.bin check2 live 75 0.013333333333333334 0.65 5 1
./nbody_long mix2.bin fine live 150 0.006666666666666667 0.65 10 1
./nbody_long mix2.bin live_60 live 4500 0.013333333333333333 0.65 5 1
./nbody_long mix2.bin forced_60 forced 4500 0.013333333333333333 0.65 5 1
python3 render.py
python3 validate.py
