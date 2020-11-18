#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>

typedef double value_t;

#define RESOLUTION 120

// -- vector utilities --

typedef value_t *Vector;

Vector createVector(int N);

void printVec(Vector m, int N);

void releaseVector(Vector m);

void printTemperature(Vector m, int N);

// -- simulation code ---

int main(int argc, char **argv) {
    // 'parsing' optional input parameter = problem size
    int N = 2048;
    if (argc > 1) {
        N = atoi(argv[1]);
    }
    int T = N * 500;
    int rank, numProcs;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

    if (N % numProcs != 0) {
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    // ---------- setup ----------
    // create a buffer for storing temperature fields
    Vector A = createVector(N);

    // set up initial conditions in A
    int source_x = N / 4;
    if (rank == 0) {
        std::cout << "Computing heat-distribution for room size N=" << N << " for T= " << T << " timesteps" << std::endl;
        for (int i = 0; i < N; i++) {
            A[i] = 273;  // temperature is 0° C everywhere (273 K)
        }
        // and there is a heat source in one corner
        A[source_x] = 273 + 60;

        // std::cout << "Initial:\t";
        // printTemperature(A, N);
        // std::cout << std::endl;
    }
    MPI_Bcast(A, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // ---------- compute ----------
    // create a second buffer for the computation
    Vector B = createVector(N);

    // for each time step ..
    long long global_pos;
    long subrange = N / numProcs;
    auto start = std::chrono::high_resolution_clock::now();
    for (int t = 0; t < T; t++) {
        // send and receieve data from neighboors
        if ((rank % 2) == 0) {
            // all even ranks send
            // right element
            if (rank != (numProcs - 1))
                MPI_Send(&A[rank * subrange + subrange - 1], 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
            // left element
            if (rank != 0)
                MPI_Send(&A[rank * subrange], 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
            // all even ranks receieve
            //left element
            if (rank != 0)
                MPI_Recv(&A[rank * subrange - 1], 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            // right element
            if (rank != (numProcs - 1))
                MPI_Recv(&A[rank * subrange + subrange], 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        } else {
            // all odd ranks receive first
            // odd rank always has left element
            MPI_Recv(&A[rank * subrange - 1], 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            // right element
            if (rank != (numProcs - 1))
                MPI_Recv(&A[rank * subrange + subrange], 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            // all odd ranks send
            // right element
            if (rank != (numProcs - 1))
                MPI_Send(&A[rank * subrange + subrange - 1], 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
            // odd rank always have a left rank
            MPI_Send(&A[rank * subrange], 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
        }

        // .. we propagate the temperature
        for (long long local_pos = 0; local_pos < subrange; local_pos++) {
            global_pos = rank * subrange + local_pos;
            // center stays constant (the heat is still on)
            if (global_pos == source_x) {
                B[global_pos] = A[global_pos];
                continue;
            }

            // get temperature at current position
            value_t tc = A[global_pos];

            // get temperatures of adjacent cells
            value_t tl = (global_pos != 0) ? A[global_pos - 1] : tc;
            value_t tr = (global_pos != N - 1) ? A[global_pos + 1] : tc;

            // compute new temperature at current position
            B[global_pos] = tc + 0.2 * (tl + tr + (-2 * tc));
        }

        // swap matrices (just pointers, not content)
        Vector H = A;
        A = B;
        B = H;

        MPI_Barrier(MPI_COMM_WORLD);

        // show intermediate step
        //if (!(t % 1000) && rank == 0) {
        //     MPI_Gather(&A[rank * subrange], subrange, MPI_DOUBLE, A, subrange, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        //     if (rank == 0) {
        //         std::cout << "Step t=" << t << ":\t" << std::endl;
        //         printTemperature(A, N);
        //         std::cout << std::endl;
        //     }
        //}
    }

    MPI_Gather(&A[rank * subrange], subrange, MPI_DOUBLE, A, subrange, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    auto end = std::chrono::high_resolution_clock::now();

    releaseVector(B);

    // ---------- check ----------
    bool success = true;
    if (rank == 0) {
        // std::cout << "Final:\t\t";
        // printTemperature(A, N);
        // std::cout << std::endl;

        for (long long i = 0; i < N; i++) {
            value_t temp = A[i];
            if (273 <= temp && temp <= 273 + 60)
                continue;
            success = false;
            break;
        }

        printf("Verification: %s\n", (success) ? "OK" : "FAILED");
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        std::cout << "Time: " << duration.count() / 1e6 << "sec" << std::endl;
    }

    // ---------- cleanup ----------
    releaseVector(A);

    // done
    return (success) ? EXIT_SUCCESS : EXIT_FAILURE;
}

Vector createVector(int N) {
    // create data and index vector
    return (value_t *)malloc(sizeof(value_t) * N);
}

void releaseVector(Vector m) { free(m); }

void printVec(Vector m, int size) {
    for (long i = 0; i < size; i++)
        std::cout << m[i] << " ";
    std::cout << std::endl;
}

void printTemperature(Vector m, int N) {
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
    for (int i = 0; i < W; i++) {
        // get max temperature in this tile
        value_t max_t = 0;
        for (int x = sW * i; x < sW * i + sW; x++) {
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
