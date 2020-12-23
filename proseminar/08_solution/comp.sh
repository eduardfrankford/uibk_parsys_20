#!/bin/sh

gcc -O3 -std=c99 -fopenmp mat_mul.c -o mat_mul

./mat_mul 100 1
./mat_mul 100 2
./mat_mul 100 4
./mat_mul 100 8

./mat_mul 500 1
./mat_mul 500 2
./mat_mul 500 4
./mat_mul 500 8

./mat_mul 1000 1
./mat_mul 1000 2
./mat_mul 1000 4
./mat_mul 1000 8

./mat_mul 2000 1
./mat_mul 2000 2
./mat_mul 2000 4
./mat_mul 2000 8

./mat_mul 3000 1
./mat_mul 3000 8

gcc -O3 -std=c99 -fopenmp ex_2_sequential.c -o n_queens_seq
gcc -O3 -std=c99 -fopenmp ex_2_parallel.c -o n_queens_par

./n_queens_seq 8
./n_queens_par 8 2
./n_queens_par 8 4
./n_queens_par 8 8

./n_queens_seq 10
./n_queens_par 10 2
./n_queens_par 10 4
./n_queens_par 10 8

./n_queens_seq 12
./n_queens_par 12 2
./n_queens_par 12 4
./n_queens_par 12 8

./n_queens_seq 13
./n_queens_par 13 8