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

    int samples = atoi(argv[1]);
    int inside = 0;
    double x, y;
    srand(time(NULL));

    auto start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < samples; i++) {
        x = (double)(rand() % std::numeric_limits<int>::max()) / ((double)std::numeric_limits<int>::max() / 2) - 1;
        y = (double)(rand() % std::numeric_limits<int>::max()) / ((double)std::numeric_limits<int>::max() / 2) - 1;
        if (sqrt(x * x + y * y) <= 1)
            inside++;
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    std::cout << "Time: " << duration.count() / 1e6 << "sec" << std::endl;
    std::cout << std::fixed;
    std::cout << "Approximation of pi: " << std::setprecision(10) << 4 * (double)inside / (double)samples << std::endl;

    return EXIT_SUCCESS;
}
