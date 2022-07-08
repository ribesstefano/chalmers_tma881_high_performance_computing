#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void row_sums(double* sums, const double** matrix, size_t nrs, size_t ncs) {
  for (size_t ix = 0; ix < nrs; ++ix) {
    double sum = 0.;
    for (size_t jx = 0; jx < ncs; ++jx) {
      sum += matrix[ix][jx];
    }
    sums[ix] = sum;
  }
}

void col_sums(double* sums, const double** matrix, size_t nrs, size_t ncs) {
  for (size_t jx = 0; jx < ncs; ++jx) {
    double sum = 0.;
    for (size_t ix = 0; ix < nrs; ++ix) {
      sum += matrix[ix][jx];
    }
    sums[jx] = sum;
  }
}

int main(int argc, char const *argv[]) {
  /**
   Write a C program, called locality, that

   * creates 1000 x 1000 doubles,
   * stores them in a matrix in row major order,
   * computes the row and the column sums.
  */
  srand(time(NULL));
  const int size = 1000;
  const int num_tests = 5000;
  double* asentries = (double*) malloc(sizeof(double) * size * size);
  double** as = (double**) malloc(sizeof(double*) * size);
  double* rsums = (double*) malloc(sizeof(double) * size);
  double* csums = (double*) malloc(sizeof(double) * size);
  for (size_t ix = 0, jx = 0; ix < size; ++ix, jx += size) {
    as[ix] = asentries + jx;
  }
  for (size_t ix = 0; ix < size; ++ix) {
    for (size_t jx = 0; jx < size; ++jx) {
      as[ix][jx] = rand();
    }
  }
  // Start computation
  // Rows
  clock_t start = clock();
  for (int t = 0; t < num_tests; ++t) {
    row_sums(rsums, (const double**)as, size, size);
  }
  clock_t end = clock();
  float tot_time = (float)(end - start) / CLOCKS_PER_SEC;
  // Printing
  printf("[ROWS] Random element: %d\n", rsums[rand() % size]);
  printf("[ROWS] Total time:     %f [s]\n", tot_time);
  printf("[ROWS] Time per iter:  %f [s]\n", tot_time / (float) (num_tests));
  // Columns
  start = clock();
  for (int t = 0; t < num_tests; ++t) {
    col_sums(csums, (const double**)as, size, size);
  }
  end = clock();
  tot_time = (float)(end - start) / CLOCKS_PER_SEC;
  // Printing
  printf("\n");
  printf("[COLS] Random element: %d\n", csums[rand() % size]);
  printf("[COLS] Total time:     %f [s]\n", tot_time);
  printf("[COLS] Time per iter:  %f [s]\n", tot_time / (float) (num_tests));
  // Free memory
  free(as);
  free(asentries);
  free(rsums);
  free(csums);
  return 0;
}