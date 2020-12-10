#!/bin/sh

mpicc -O3 -std=c99 -fopenmp heat_stencil_1D_mpi.c -o heat_stencil_1D_mpi

