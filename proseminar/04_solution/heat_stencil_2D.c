#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

typedef double value_t;

#define RESOLUTION_WIDTH 50
#define RESOLUTION_HEIGHT 50

// -- vector utilities --

typedef value_t *Vector;

Vector createVector(int N);

void releaseVector(Vector m);

void printTemperature(Vector m, int height, int width);

void print2DArr(Vector A, int N);
// -- simulation code ---

int main(int argc, char **argv)
{
    double startTime = omp_get_wtime();
    // 'parsing' optional input parameter = problem size
    int N = 4000;
    int height = N;
    int width = N;
    if (argc > 1)
    {
        N = atoi(argv[1]);
    }
    int T = 100;
    printf("Computing heat-distribution for room size N=%d for T=%d timesteps\n", N, T);

    // ---------- setup ----------

    // create a buffer for storing temperature fields
    Vector A = createVector(N);

    // set up initial conditions in A
    for (int i = 0; i < height * width; i++)
    {
        A[i] = 273; // temperature is 0° C everywhere (273 K)
    }

    // and there is a heat source in one corner
    int source_x = width / 4;
    int source_y = height / 4;
    A[source_y * width + source_x] = 273 + 60;

    // printf("Initial:\t\n");
    // printTemperature(A, height, width);
    // //print2DArr(A,N);
    // printf("\n");

    // //print2DArr(A, N);

    // ---------- compute ----------

    // create a second buffer for the computation
    Vector B = createVector(N);

    // for each time step ..
    for (int t = 0; t < T; t++)
    {
        // .. we propagate the temperature
        for (long long i = 0; i < height; i++)
        {
            for (long long j = 0; j < width; j++)
            {
                // center stays constant (the heat is still on)
                if (j == source_x && i == source_y)
                {
                    B[source_y * N + source_x] = A[source_y * N + source_x];
                    continue;
                }

                // get current temperature at (i,j)
                value_t tc = A[i * width + j];

                // get temperatures left/right and up/down
                value_t tl = (j != 0) ? A[i * width + (j - 1)] : tc;
                value_t tr = (j != width - 1) ? A[i * width + (j + 1)] : tc;
                value_t tu = (i != 0) ? A[(i - 1) * width + j] : tc;
                value_t td = (i != height - 1) ? A[(i + 1) * width + j] : tc;

                // update temperature at current point
                B[i * width + j] = tc + 0.2 * (tl + tr + tu + td + (-4.0f * tc));
            }
        }

        // swap matrices (just pointers, not content)
        Vector H = A;
        A = B;
        B = H;

        // // show intermediate step
        // if (!(t % 10))
        // {
        //     printf("Step t=%d:\t\n", t);
        //     printTemperature(A, height, width);
        //     printf("\n");
        // }
    }

    releaseVector(B);

    // ---------- check ----------

    // printf("Final:\t\t\n");
    // printTemperature(A, height, width);
    // printf("\n");

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
    double endTime = omp_get_wtime();

    printf("time: %2.2f seconds\n\n", endTime - startTime);
    return EXIT_SUCCESS;
}

Vector createVector(int N)
{
    // create data and index vector
    return malloc(sizeof(value_t) * N * N);
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
    int W = RESOLUTION_WIDTH;
    int H = RESOLUTION_HEIGHT;

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
