#!/usr/bin/bash -f

export OMP_NUM_THREADS=8
gcc -O3 -fopenmp ex.c -lm -o a.out
./a.out 10000 100 1
#gnuplot particle.plt