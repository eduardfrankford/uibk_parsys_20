#!/bin/sh

gcc -O3 -std=c99 -fopenmp monte_carlo_atomic.c -o montecarlopi_atomic
gcc -O3 -std=c99 -fopenmp monte_carlo_reduction.c -o montecarlopi_reduction
gcc -O3 -std=c99 -fopenmp monte_carlo_critical.c -o montecarlopi_critical
gcc -O3 -std=c99 -fopenmp heat_stencil_2D.c -o heat_stencil_2D
gcc -O3 -std=c99 -fopenmp heat_stencil_2D_omp.c -o heat_stencil_2D_omp