#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include <omp.h>

int main(int argc, char **argv)
{
    if (argc != 3)
    {
        fprintf(stderr, "Error: usage: %s <n_threads>\n", argv[0]);
        return EXIT_FAILURE;
    }
    int num_threads = atoi(argv[1]);
    omp_set_num_threads(num_threads);

    int n = atoi(argv[2]);
    int count = 0;

    double startTime = omp_get_wtime();
    unsigned int seed;

#pragma omp parallel shared(n, count) private(seed)
    {
        seed = omp_get_thread_num();
        int localcount = 0;
#pragma omp for
        for (int i = 0; i < n; i++)
        {
            double x = (double)rand_r(&seed) / RAND_MAX;
            double y = (double)rand_r(&seed) / RAND_MAX;
            double distance = x * x + y * y;

            if (distance > 1.0)
            {
                localcount++;
            }
        }
        
        #pragma omp flush
        #pragma omp critical
        count = localcount + count;
        #pragma omp flush

    }

    double endTime = omp_get_wtime();

    printf("N: %d Critical (%d threads): %2.4f s\n", n,num_threads, endTime - startTime);

    double pi = (double)(n - count) * 4 / n;
    printf("PI %f\n", pi);

    return EXIT_SUCCESS;
}