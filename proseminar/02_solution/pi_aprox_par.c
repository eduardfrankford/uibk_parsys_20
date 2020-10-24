/**
 *
 * Author     : Adrian Statescu mergesortv@gmail.com http://adrianstatescu.com
 *
 * Description: C Program to compute PI using a Monte Carlo Method.
 *
 * MIT License 
 *
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

int main(int argc, char *argv[])
{

    MPI_Init(&argc, &argv); // initialize the MPI environment

    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size); // get the number of ranks

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // get the rank of the caller

    printf("Hello world from rank %d of %d\n", rank, size);

    double startTime = (float)clock() / CLOCKS_PER_SEC;
    srand(time(NULL) + rank);
    int i, count, n;
    double x, y, z, pi;
    n = 500000000;

    int partition_n = n / size;

    count = 0;
    int global_count;

    for (i = 0; i < partition_n; ++i)
    {

        x = (double)rand() / RAND_MAX;

        y = (double)rand() / RAND_MAX;

        z = x * x + y * y;

        if (z <= 1)
            count++;
    }

    MPI_Reduce(&count, &global_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    printf("Local count: %d\n",count);
    if (rank == 0)
    {
        pi = (double)global_count / n * 4;
        printf("Global count: %d\n", global_count);
        printf("Approximate PI = %g\n", pi);
    }
    double endTime = (float)clock() / CLOCKS_PER_SEC;

    printf("time: %2.2f seconds\n", endTime - startTime);

    MPI_Finalize(); // cleanup

    return (0);
}