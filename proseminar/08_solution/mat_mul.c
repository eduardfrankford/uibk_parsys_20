
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
typedef double value_t;
typedef value_t *Vector;
Vector createVector(int N);

void releaseVector(Vector m);

void printMatrix(Vector C, int N)
{
    /* Printing the contents of third matrix. */
    printf("\n\nMatrix :");
    for (int i = 0; i < N; i++)
    {
        printf("\n\t\t\t");
        for (int j = 0; j < N; j++)
            printf("%f\t", C[i * N + j]);
    }
}

// Multiplies two matrices mat1[][] and mat2[][]
// and prints result.
// (m1) x (m2) and (n1) x (n2) are dimensions
// of given matrices.
void multiply(Vector A, Vector B, Vector C, int N)
{
    // printMatrix(A, N);

    // printf("\n");

    // printMatrix(B, N);
    /* Calculation begins for the resultant matrix. */
    #pragma omp parallel for
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            double localsum = 0.0;
            int offset = N*i;
            for (int k = 0; k < N; k++){
                
                localsum = localsum +   A[offset + k] * B[offset +k];} //remember transposed matrix B otherwise do:
                //localsum = localsum +   A[offset + k] * B[N*k + j];}
                C[offset + j] = localsum;
        }
    }
    /* Printing the contents of third matrix. */
//     printf("\n\nResultant matrix :");
    //printMatrix(C, N);
 }

// Driver code
int main(int argc, char **argv)
{
    int N = 10;
    int height = N;
    int width = N;
    int num_threads = 1;

    if (argc > 1)
    {
        N = atoi(argv[1]);
        num_threads = atoi(argv[2]);
    }
    omp_set_num_threads(num_threads);
    Vector A = createVector(N);
    Vector B = createVector(N);
    Vector C = createVector(N);

    #pragma omp parallel for
    for (int i = 0; i < N; i++)
    {
        int offset = i * N;
        for (int j = 0; j < N; j++)
        {

            A[offset + j] = i * j;
            B[N*j + i] = (i == j) ? 1 : 0; //Transposed matrix otherwise init with 
            // B[offset + j] = (i == j) ? 1 : 0;
            C[offset + j] = 0;
        }
    }
    double startTime = omp_get_wtime();
    multiply(A, B, C, N);

    // done
    double endTime = omp_get_wtime();

    printf("%d, %2.2f, %d\n ", N,endTime - startTime, num_threads);
    int success = 1;
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            value_t temp = A[i];
            if (A[i * N + j] == A[i * N + j])
                continue;
            success = 0;
            break;
        }
    }
    //printf("Verification: %s\n", (success) ? "OK" : "FAILED");

    releaseVector(A);
    releaseVector(B);
    releaseVector(C);

    return 0;
}

Vector createVector(int N)
{
    // create data and index vector
    return malloc(sizeof(value_t) * N * N);
}

void releaseVector(Vector m) { free(m); }
