#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <math.h>

typedef double value_t;

#define RESOLUTION 60

// -- vector utilities --

typedef value_t *Vector;

Vector createVector(int N);

void releaseVector(Vector m);

void printTemperature(Vector m, int N);

// -- simulation code ---

int main(int argc, char **argv)
{

    MPI_Init(&argc, &argv); // initialize the MPI environment
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size); // get the number of ranks

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // get the rank of the caller

    double startTime = (float)clock() / CLOCKS_PER_SEC;
    // 'parsing' optional input parameter = problem size
    int N = 10000;
    if (argc > 1)
    {
        N = atoi(argv[1]);
    }
    int T = N * 500;

    if (N % size != 0)
    {
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    printf("Computing heat-distribution for room size N=%d for T=%d timesteps\n", N, T);

    // ---------- setup ----------

    // create a buffer for storing temperature fields
    Vector A = createVector(N);

    // set up initial conditions in A
    for (int i = 0; i < N; i++)
    {
        A[i] = 273; // temperature is 0° C everywhere (273 K)
    }

    // and there is a heat source in one corner
    int source_x = N / 4;
    A[source_x] = 273 + 60;

    printf("Initial:\t");
    printTemperature(A, N);
    printf("\n");

    // ---------- compute ----------

    // create a second buffer for the computation
    Vector B = createVector(N);
    Vector B_part = createVector(N / size);
    int number;

    // for each time step ..
    for (int t = 0; t < T; t++)
    {
        // .. we propagate the temperature
        long long lower = rank * (N / size);
        long long upper = rank * (N / size) + N / size;
        //MPI_Bcast(A, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        int j = 0;
        for (int i = lower; i < upper; i++)
        {
            // center stays constant (the heat is still on)
            if (i == source_x)
            {
                B_part[j++] = A[i];
                continue;
            }

            // get temperature at current position
            value_t tc = A[i];

            // get temperatures of adjacent cells
            value_t tl = (i != 0) ? A[i - 1] : tc;
            value_t tr = (i != N - 1) ? A[i + 1] : tc;

            if (i != 0 && i == lower)
            {
                tl = A[i - 1];
            }
            if (i != N - 1 && i == upper - 1)
            {
                tr = A[i + 1];
            }
            // compute new temperature at current position
            B_part[j++] = tc + 0.2 * (tl + tr + (-2 * tc));
        }
        MPI_Allgather(B_part, N / size, MPI_DOUBLE,B,N / size,MPI_DOUBLE,MPI_COMM_WORLD);
        //        MPI_Gather(B_part, N / size, MPI_DOUBLE,B[lower],N / size,MPI_DOUBLE,0,MPI_COMM_WORLD); -> slow

        // swap matrices (just pointers, not content)
        Vector H = A;
        A = B;
        B = H;

        // show intermediate step
        if (!(t % 100000) && rank == 0)
        {
            printf("Step t=%d:\t", t);
            printTemperature(A, N);
            printf("\n");
        }
    }

    /*for(int i = 0; i < N-1000;i++){
        if(rank == 0)
        printf("Rank = %d Heat: %f with i: %d\n",rank, A[i],i);
    }*/

    releaseVector(B);

    // ---------- check ----------

    printf("Final:\t\t");
    printTemperature(A, N);
    printf("\n");

    int success = 1;
    for (long long i = 0; i < N; i++)
    {
        value_t temp = A[i];
        if (273 <= temp && temp <= 273 + 60)
            continue;
        success = 0;
        break;
    }

    printf("Verification: %s\n", (success) ? "OK" : "FAILED");

    // ---------- cleanup ----------

    releaseVector(A);

    // done
    double endTime = (float)clock() / CLOCKS_PER_SEC;

    printf("time: %2.2f seconds\n", endTime - startTime);
    MPI_Finalize(); // cleanup
    return (success) ? EXIT_SUCCESS : EXIT_FAILURE;
}

Vector createVector(int N)
{
    // create data and index vector
    return malloc(sizeof(value_t) * N);
}

void releaseVector(Vector m) { free(m); }

void printTemperature(Vector m, int N)
{
    const char *colors = " .-:=+*^X#%@";
    const int numColors = 12;

    // boundaries for temperature (for simplicity hard-coded)
    const value_t max = 273 + 30;
    const value_t min = 273 + 0;

    // set the 'render' resolution
    int W = RESOLUTION;

    // step size in each dimension
    int sW = N / W;

    // room
    // left wall
    printf("X");
    // actual room
    for (int i = 0; i < W; i++)
    {
        // get max temperature in this tile
        value_t max_t = 0;
        for (int x = sW * i; x < sW * i + sW; x++)
        {
            max_t = (max_t < m[x]) ? m[x] : max_t;
        }
        value_t temp = max_t;

        // pick the 'color'
        int c = ((temp - min) / (max - min)) * numColors;
        c = (c >= numColors) ? numColors - 1 : ((c < 0) ? 0 : c);

        // print the average temperature
        printf("%c", colors[c]);
    }
    // right wall
    printf("X");
}