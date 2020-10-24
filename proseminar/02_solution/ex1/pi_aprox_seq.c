/**
 *
 * Author     : Adrian Statescu mergesortv@gmail.com http://adrianstatescu.com
 *
 * Description: C Program to compute PI using a Monte Carlo Method.
 *
 * MIT License 
 *
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#define SEED time(NULL)

int main(int argc, char *argv[])
{
  double startTime = (float)clock() / CLOCKS_PER_SEC;
  srand(SEED);

  int i, count, n;
  double x, y, z, pi;

  n = atoi(argv[1]);

  count = 0;

  for (i = 0; i < n; ++i)
  {

    x = (double)rand() / RAND_MAX;

    y = (double)rand() / RAND_MAX;

    z = x * x + y * y;

    if (z <= 1)
      count++;
  }

  pi = (double)count / n * 4;

  printf("Sequential Version with:\n n = %d\n pi =  %g\n",n, pi);

  double endTime = (float)clock() / CLOCKS_PER_SEC;

  printf("time: %2.2f seconds\n", endTime - startTime);

  return (0);
}
