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
    int N = 1000;
    if (argc > 1)
    {
        N = atoi(argv[1]);
    }

    if (N % size != 0)
    {
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    // create a buffer for storing temperature fields
    Vector A = createVector(N);

    // set up initial conditions in A
    for (int i = 0; i < N; i++)
    {
        A[i] = i; // temperature is 0° C everywhere (273 K)
    }

    // ---------- compute ----------

    // create a second buffer for the computation
    Vector B = createVector(N);
    Vector B_part = createVector(N / size);
    int number;

    // .. we propagate the temperature
    long long lower = rank * (N / size);
    long long upper = rank * (N / size) + N / size;
    MPI_Bcast(A, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    int j = 0;
    for (int i = lower; i < upper; i++)
    {
        B_part[j++] = A[i];
    }
    MPI_Allgather(B_part, N / size, MPI_DOUBLE, B, N/size, MPI_DOUBLE, MPI_COMM_WORLD);

    if (rank == 0)
    {
        for (int i = 0; i < N; i++)
        {
            printf("%f\n", B[i]);
        }
    }
    // swap matrices (just pointers, not content)
    Vector H = A;
    A = B;
    B = H;

    /*for(int i = 0; i < N-1000;i++){
        if(rank == 0)
        printf("Rank = %d Heat: %f with i: %d\n",rank, A[i],i);
    }*/

    releaseVector(B);

    // ---------- check ----------

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