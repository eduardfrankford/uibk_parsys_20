#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <unistd.h>

typedef double value_t;

#define RESOLUTION 60

// -- vector utilities --

typedef value_t *Vector;

Vector createVector(int height, int width);

void releaseVector(Vector m);

void printTemperature(Vector m, int height, int width);

void print2DArr(Vector A, int height, int width);
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
    int N = 200;
    int height = N;
    int width = N / size;

    MPI_Datatype myColumn;
    MPI_Type_vector(height, 1, width, MPI_DOUBLE, &myColumn);
    MPI_Type_commit(&myColumn);

    if (argc > 1)
    {
        N = atoi(argv[1]);
    }
    int T = N * 500;
    printf("Computing heat-distribution for room size N=%d for T=%d timesteps\n", N, T);

    // ---------- setup ----------

    // create a buffer for storing temperature fields
    Vector A = createVector(height, width);

    // set up initial conditions in A
    for (int i = 0; i < height * width; i++)
    {
        A[i] = 273; // temperature is 0° C everywhere (273 K)
    }

    // and there is a heat source in one corner
    int source_x = width / 4;
    int source_y = height / 4;

    if (rank == 0)
        A[source_y * width + source_x] = 273 + 60;

    // //print2DArr(A, N);

    // ---------- compute ----------

    // create a second buffer for the computation
    Vector B = createVector(height, width);
    // for each time step ..
    long long global_pos;
    long subrange = N / size;
    // for each time step ..
    Vector dataRight = createVector(height, 1);
    Vector dataLeft = createVector(height, 1);
    Vector A_total = createVector(height, N);

    if (rank == 0)
    {
        printf("Initial:\t\n");
        printTemperature(A, height, width);
        printf("\n");
    }
    value_t tl, tr, tu, td;

    for (int t = 0; t < T; t++)
    {

        // .. we propagate the temperature

        // send and receieve data from neighboors
        if ((rank % 2) == 0)
        {
            // all even ranks send
            // right element
            if (rank != (size - 1))
                MPI_Send(A + width - 1, 1, myColumn, rank + 1, 0, MPI_COMM_WORLD);
            // left element
            if (rank != 0)
                MPI_Send(A, 1, myColumn, rank - 1, 0, MPI_COMM_WORLD);
            // all even ranks receieve
            //left element
            if (rank != 0)
                MPI_Recv(dataLeft, height, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            // right element
            if (rank != (size - 1))
                MPI_Recv(dataRight, height, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        else
        {
            // all odd ranks receive first
            // odd rank always has left element
            MPI_Recv(dataLeft, height, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            // right element
            if (rank != (size - 1))
                MPI_Recv(dataRight, height, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            // all odd ranks send
            // right element
            if (rank != (size - 1))
                MPI_Send(A + width - 1, 1, myColumn, rank + 1, 0, MPI_COMM_WORLD);
            // odd rank always have a left rank
            MPI_Send(A, 1, myColumn, rank - 1, 0, MPI_COMM_WORLD);
        }

        for (long long i = 0; i < height; i++)
        {
            for (long long j = 0; j < width; j++)
            {
                // center stays constant (the heat is still on)
                if (j == source_x && i == source_y && rank == 0)
                {
                    //printf("I am here %lld, %lld with %f\n",i,j,A[source_y * width + source_x] );
                    B[source_y * width + source_x] = A[source_y * width + source_x];
                    continue;
                }

                // get current temperature at (i,j)
                value_t tc = A[i * width + j];
                if (rank == 0)
                {
                    tl = (j != 0) ? A[i * width + (j - 1)] : tc;
                    tr = (j != width - 1) ? A[i * width + (j + 1)] : dataRight[i];
                    tu = (i != 0) ? A[(i - 1) * width + j] : tc;
                    td = (i != height - 1) ? A[(i + 1) * width + j] : tc;
                }
                else if (rank == size - 1)
                {
                    tl = (j != 0) ? A[i * width + (j - 1)] : dataLeft[i];
                    tr = (j != width - 1) ? A[i * width + (j + 1)] : tc;
                    tu = (i != 0) ? A[(i - 1) * width + j] : tc;
                    td = (i != height - 1) ? A[(i + 1) * width + j] : tc;
                }
                else
                {
                    tl = (j != 0) ? A[i * width + (j - 1)] : dataLeft[i];
                    tr = (j != width - 1) ? A[i * width + (j + 1)] : dataRight[i];
                    tu = (i != 0) ? A[(i - 1) * width + j] : tc;
                    td = (i != height - 1) ? A[(i + 1) * width + j] : tc;
                }
                B[i * width + j] = tc + 0.2 * (tl + tr + tu + td + (-4.0f * tc));

                // update temperature at current point
            }
        }

        // swap matrices (just pointers, not content)
        Vector H = A;
        A = B;
        B = H;

        // show intermediate step
        if (!(t % 1000))
        {
            //MPI_Gather(A, width * height, MPI_DOUBLE, A_total + rank * width * height, width * height, MPI_DOUBLE, 0, MPI_COMM_WORLD);

            // if (rank == 0)
            // {
            //     printf("Step t=%d rank = %d:\t\n", t,rank);
            //     printTemperature(A, height, width);
            //     printf("\n");
            // }

            for (int i = 0; i < height; i++)
            {
                MPI_Gather(&A[i*width], width, MPI_DOUBLE, &A_total[i*width*size + rank*width], subrange, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            }

            if (rank == 0)
            {
                printf("Step t=%d rank = %d:\t\n", t, rank);
                printTemperature(A_total, height, width*size);
                printf("\n");
            }
            // if (rank == 2)
            // {
            //     printf("Step t=%d rank = %d:\t\n", t,rank);
            //     printTemperature(A, height, width);
            //     printf("\n");
            // }
            // if (rank == 3)
            // {
            //     printf("Step t=%d rank = %d:\t\n", t,rank);
            //     printTemperature(A, height, width);
            //     printf("\n");
            // }
        }
    }

    releaseVector(B);

    // ---------- check ----------

    // for (int i = 0; i < size; i++)
    // {
    //     for (int j = 0; j < width; j++)
    //     {
    //         /* code */
    //     }

    // }

    printf("Final:\t\t\n");
    printTemperature(A, height, width);
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
    return EXIT_SUCCESS;
}

Vector createVector(int height, int width)
{
    // create data and index vector
    return malloc(sizeof(value_t) * height * width);
}

void releaseVector(Vector m) { free(m); }

void printTemperature(Vector m, int height, int width)
{
    const char *colors = " .-:=+*^X#%@";
    const int numColors = 12;

    // boundaries for temperature (for simplicity hard-coded)
    const value_t max = 273 + 30;
    const value_t min = 273 + 0;

    // set the 'render' resolution
    int W = RESOLUTION;
    int H = 20;

    // step size in each dimension
    int sH = height / H;
    int sW = width / W;

    // room
    // left wall
    // actual room
    for (int i = 0; i < H; i++)
    {
        printf("X");
        for (int j = 0; j < W; j++)
        {
            // get max temperature in this tile
            value_t max_t = 0;
            for (int x = sH * i; x < sH * i + sH; x++)
            {
                for (int y = sW * j; y < sW * j + sW; y++)
                {
                    max_t = (max_t < m[x * width + y]) ? m[x * width + y] : max_t;
                }
            }
            value_t temp = max_t;

            // pick the 'color'
            int c = ((temp - min) / (max - min)) * numColors;
            c = (c >= numColors) ? numColors - 1 : ((c < 0) ? 0 : c);

            // print the average temperature
            printf("%c", colors[c]);
        }

        // right wall
        printf("X\n");
    }
}

void print2DArr(Vector A, int height, int width)
{
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            if (A[i * width + j] > 273.01)
                printf("%f ", *(A + i * width + j));
        }
        //printf("\n");
    }
}
