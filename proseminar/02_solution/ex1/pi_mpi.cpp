#include <mpi.h>
#include <time.h>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cout << "Expected input: ./FILE_NAME NUM_SAMPLES" << std::endl;
        return EXIT_FAILURE;
    }

    int rank, numProcs;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

    int samples = atoi(argv[1]);
    int rank_samples = samples / numProcs;
    double x, y;
    srand(time(NULL)+rank);
    int inside = 0;

    if (samples % numProcs != 0) {
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    auto start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < rank_samples; i++) {
        x = (double)(rand() % std::numeric_limits<int>::max()) / ((double)std::numeric_limits<int>::max() / 2) - 1;
        y = (double)(rand() % std::numeric_limits<int>::max()) / ((double)std::numeric_limits<int>::max() / 2) - 1;
        if (sqrt(x * x + y * y) <= 1)
            inside++;
    }

    int res = 0;
    MPI_Reduce(&inside, &res, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    auto end = std::chrono::high_resolution_clock::now();

    if (rank == 0) {
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        std::cout << "Time: " << duration.count() / 1e6 << "sec" << std::endl;
        std::cout << std::fixed;
        std::cout << "Approximation of pi: " << std::setprecision(10) << 4 * (double)res / (double)samples << std::endl;
    }

    MPI_Finalize();

    return EXIT_SUCCESS;
}
