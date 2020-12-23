#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include <omp.h>

int isNotAttack(int *placedQueens, int row, int col, int n)
{
    int i, j;

    // Row left
    for (i = 0; i < col; i++)
        if (placedQueens[i] == row)
            return 0;

    // Diagonal top left
    for (i = row - 1, j = col - 1; i >= 0 && j >= 0; i--, j--)
        if (placedQueens[j] == i)
            return 0;

    // Diagonal bottom left
    for (i = row + 1, j = col - 1; j >= 0 && i < n; i++, j--)
        if (placedQueens[j] == i)
            return 0;

    return 1;
}

int count = 0;

int numberOfNQueenSolutions(int *placedQueens, int col, int n)
{
    if (col >= n)
        return 1;
    int count = 0;

    for (int i = 0; i < n; i++)
    {
        {
            if (isNotAttack(placedQueens, i, col, n))
            {
                placedQueens[col] = i;
                int localCount = numberOfNQueenSolutions(placedQueens, col + 1, n);

                count += localCount;
            }
        }
    }

    #pragma omp taskwait
    return count;
}

int main(int argc, char **argv)
{
    if (argc != 2)
    {
        fprintf(stderr, "Error: usage: %s <n>\n", argv[0]);
        return EXIT_FAILURE;
    }
    int n = atoi(argv[1]);

    int placedQueens[n];
    for (int i = 0; i < n; i++)
        placedQueens[n] = 0;

    // Sequential:
    double start_time = omp_get_wtime();
    int solution = numberOfNQueenSolutions(placedQueens, 0, n);
    double end_time = omp_get_wtime();

    printf("%d, 1, %2.8f\n",n, end_time - start_time);
    // printf("Number of Solutions: %d\n", solution);

    return EXIT_SUCCESS;
}