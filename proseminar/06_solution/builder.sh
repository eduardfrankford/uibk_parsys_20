#!/usr/bin/bash -f

export OMP_NUM_THREADS=8
gcc -O3 -fopenmp barnes_hut_lcc2.c -lm -o a.out

for num_threads in `seq 1 1 8`
do
    for problem_size in `seq 1000 1000 10000` 
    do
        ./a.out $problem_size 100 $num_threads

    done
done

#gnuplot particle.plt