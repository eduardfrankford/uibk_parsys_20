#!/bin/bash

# Executes job in the queue "std.q" unless you have special requirements.
#$ -q std.q

# Changes to the current working directory before performing any further action
#$ -cwd

# Name of your job. Unless you use the -o and -e options, output will
# go to a unique file name.ojob_id for each job.
#$ -N Sheet4

# Redirect output stream to this file.
#$ -o sheet4.dat

# Join the error stream to the output stream.
#$ -j yes

# Use the parallel environment "openmp" with 8 job slots. Each requested
# job slot will be assigned one core on the execution host.
#$ -pe openmp 8


#$ -l h_vmem=2G

# Use gcc 8.2.0 as the default gcc
module load gcc/8.2.0

#tell OpenMP how many threads to start
export OMP_NUM_THREADS=8

for problem_size in 100000000; do
    ./montecarlopi_reduction 1 $problem_size
    ./montecarlopi_reduction 2 $problem_size
    ./montecarlopi_reduction 4 $problem_size
    ./montecarlopi_reduction 8 $problem_size
    ./montecarlopi_atomic 1 $problem_size
    ./montecarlopi_atomic 2 $problem_size
    ./montecarlopi_atomic 4 $problem_size
    ./montecarlopi_atomic 8 $problem_size
    ./montecarlopi_critical 1 $problem_size
    ./montecarlopi_critical 2 $problem_size
    ./montecarlopi_critical 4 $problem_size
    ./montecarlopi_critical 8 $problem_size
done


for problem_size in 4000 8000 12000; do
    ./heat_stencil_2D $problem_size
    ./heat_stencil_2D_omp $problem_size 1
    ./heat_stencil_2D_omp $problem_size 2
    ./heat_stencil_2D_omp $problem_size 4
    ./heat_stencil_2D_omp $problem_size 8
done