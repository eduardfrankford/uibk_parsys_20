#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include <omp.h>

int numberOfNQueenSolutionsSeq(int *placedQueens, int col, int n);

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



int numberOfNQueenSolutions(int *placedQueens, int col, int n)
{
    if (col >= n)
    {
        return 1;
    }
    int count = 0;


    #pragma omp parallel for 
    for (int i = 0; i < n; i++)
    {
        
            if (isNotAttack(placedQueens, i, col, n))
            {
                int newPlacedQueens[n];
                for (int a = 0; a < n; a++)
                    newPlacedQueens[a] = placedQueens[a];
                newPlacedQueens[col] = i;
                int localCount = numberOfNQueenSolutionsSeq(newPlacedQueens, col + 1, n);

                // int localCount = numberOfNQueenSolutions(placedQueens, col + 1, n);
                #pragma omp atomic
                count += localCount;
            }
        
    }





#pragma omp taskwait
    return count;
}


int numberOfNQueenSolutionsSeq(int *placedQueens, int col, int n)
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
                int localCount = numberOfNQueenSolutionsSeq(placedQueens, col + 1, n);

                count += localCount;
            }
        }
    }

    #pragma omp taskwait
    return count;
}

int main(int argc, char **argv)
{
    if (argc != 3)
    {
        fprintf(stderr, "Error: usage: %s <n> <n_threads>\n", argv[0]);
        return EXIT_FAILURE;
    }
    int n = atoi(argv[1]);
    int num_threads = atoi(argv[2]);
    omp_set_num_threads(num_threads);

    int placedQueens[n];
    #pragma omp parallel for
    for (int i = 0; i < n; i++)
        placedQueens[n] = 0;

    // Parallel:
    double start_time = omp_get_wtime();


    int solution = numberOfNQueenSolutions(placedQueens, 0, n);

    double end_time = omp_get_wtime();

    printf("%d, %d, %2.8f\n", n, num_threads, end_time - start_time);
    // printf("Number of Solutions: %d\n", solution);

    return EXIT_SUCCESS;
}