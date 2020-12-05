#!/bin/bash

# Executes job in the queue "std.q" unless you have special requirements.
#$ -q std.q

# Changes to the current working directory before performing any further action
#$ -cwd

# Name of your job. Unless you use the -o and -e options, output will
# go to a unique file name.ojob_id for each job.
#$ -N Sheet5

# Redirect output stream to this file.
#$ -o sheet5.dat

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

for num_threads in `seq 1 1 8`
do
    for problem_size in `seq 1000 1000 10000` 
    do
        ./a.out $problem_size 100 $num_threads

    done
done