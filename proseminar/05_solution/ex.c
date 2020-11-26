
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <limits.h>
#include <time.h>

float distance(float x1, float y1, float x2, float y2)
{
    float x = x1 - x2;
    float y = y1 - y2;
    return (float)sqrt(x * x + y * y);
}

int main(int argc, char **argv)
{
    if (argc != 4)
    {
        fprintf(stderr, "Error: usage: %s <N> <Timesteps> <n_threads>\n", argv[0]);
        return EXIT_FAILURE;
    }

    int N = atoi(argv[1]);
    int TIMESTEPS = atoi(argv[2]);
    int num_threads = atoi(argv[3]);
    omp_set_num_threads(num_threads);

    //FILE *fp = fopen("data.dat", "w+");

    int MAX_MASS = 100;
    int DIMENSION = 100;
    float dt = (float)1 / (TIMESTEPS * N);
    //float dt = 0.05;
    float G = 1;

    float masses[N];
    float positions[N][2];
    float velocities[N][2];

    srand(time(NULL));

    for (long i = 0; i < N; ++i)
    {
        masses[i] = ((float)rand() / (float)(RAND_MAX)) * MAX_MASS;
        positions[i][0] = ((float)rand() / (float)(RAND_MAX)) * DIMENSION;
        positions[i][1] = ((float)rand() / (float)(RAND_MAX)) * DIMENSION;
        velocities[i][0] = 0;
        velocities[i][1] = 0;
    }

    double time = 0;
    double startTime = omp_get_wtime();

    for (int i = 0; i < TIMESTEPS; i++)
    {
// #pragma omp parallel for
        for (int i = 0; i < N; i++)
        {
            float forceX = 0;
            float forceY = 0;

            for (int j = 0; j < N; j++)
            {
                if (i == j)
                    continue;

                float d = distance(positions[i][0], positions[i][1], positions[j][0], positions[j][1]);

                float angle = atan2(positions[i][1] - positions[j][1], positions[i][0] - positions[j][0]);
                float force = G * masses[i] * masses[j] / d * d;
                //float force = G * (masses[i] * masses[j]) / pow(d *d + 0.0001, 3.0 / 2.0);
                // forceX -= masses[i] * masses[j] * d * cos(angle);
                // forceY -= masses[i] * masses[j] * d * sin(angle);
                forceX -= force * cos(angle);
                forceY -= force * sin(angle);
            }
            velocities[i][0] += forceX / masses[i];
            velocities[i][1] += forceY / masses[i];
        }
// #pragma omp parallel for
        for (int i = 0; i < N; i++)
        {
            positions[i][0] += dt * velocities[i][0];
            positions[i][1] += dt * velocities[i][1];
        }

        // Excluding file writing from time calculation
        //time += omp_get_wtime() - startTime;
        // for (int i = 0; i < N; i++)
        // {
        //     fprintf(fp, "%f %f\n", positions[i][0], positions[i][1]);
        // }
        // fprintf(fp, "\n\n");
        //startTime = omp_get_wtime();
    }

    //fclose(fp);
        double endTime = omp_get_wtime();

    printf("Execution time with: %2.6f s\n", endTime-startTime);

    return EXIT_SUCCESS;
}
