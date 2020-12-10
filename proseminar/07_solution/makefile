EXECS=mpiexec --display-devel-map -n 2 heat_stencil_1D_mpi

all: barneshut ${EXECS}

barneshut: mpicc -O3 -std=c99 -fopenmp heat_stencil_1D_mpi.c -o heat_stencil_1D_mpi

clean: rm ${EXECS}