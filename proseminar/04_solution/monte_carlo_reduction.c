#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include <omp.h>

int main(int argc, char **argv)
{
    if (argc != 3)
    {
        printf("Wrong argument count.\n");
        return EXIT_FAILURE;
    }
    int num_threads = atoi(argv[1]);

    omp_set_num_threads(num_threads);

    srand(time(NULL));
    int n = atoi(argv[2]);

    double startTime = omp_get_wtime();
    unsigned int seed = num_threads + 1;
    int count = 0;

#pragma omp parallel for firstprivate(seed) reduction(+: count)
    for (int i = 0; i < n; i++)
    {
        if (seed == num_threads + 1) seed = omp_get_thread_num();
        double x = (double)rand_r(&seed) / RAND_MAX;
        double y = (double)rand_r(&seed) / RAND_MAX;
        double distance = x * x + y * y;

        if (distance > 1.0)
        {
            count++;
        }
    }

    double endTime = omp_get_wtime();

    printf("N: %d Reduction (%d threads): %2.4f s\n", n,num_threads, endTime - startTime);

    double pi = (double)(n - count) * 4 / n;
    printf("PI %f\n",pi);
    return EXIT_SUCCESS;
}