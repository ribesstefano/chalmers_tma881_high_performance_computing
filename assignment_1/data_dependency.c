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

void row_sums_unrolled2(double* sums, const double** matrix, size_t nrs,
    size_t ncs) {
  for (size_t ix = 0; ix < nrs; ++ix) {
    double sum0 = 0.;
    double sum1 = 0.;
    for (size_t jx = 0; jx < ncs; jx += 2) {
      sum0 += matrix[ix][jx];
      sum1 += matrix[ix][jx+1];
    }
    sums[ix] = sum0 + sum1;
  }
}

void row_sums_unrolled4(double* sums, const double** matrix, size_t nrs,
    size_t ncs) {
  for (size_t ix = 0; ix < nrs; ++ix) {
    double sum0 = 0.;
    double sum1 = 0.;
    double sum2 = 0.;
    double sum3 = 0.;
    for (size_t jx = 0; jx < ncs; jx += 4) {
      sum0 += matrix[ix][jx];
      sum1 += matrix[ix][jx+1];
      sum2 += matrix[ix][jx+2];
      sum3 += matrix[ix][jx+3];
    }
    sums[ix] = sum0 + sum1 + sum2 + sum3;
  }
}

void row_sums_unrolled8(double* sums, const double** matrix, size_t nrs,
    size_t ncs) {
  for (size_t ix = 0; ix < nrs; ++ix) {
    double sum0 = 0.;
    double sum1 = 0.;
    double sum2 = 0.;
    double sum3 = 0.;
    double sum4 = 0.;
    double sum5 = 0.;
    double sum6 = 0.;
    double sum7 = 0.;
    for (size_t jx = 0; jx < ncs; jx += 8) {
      sum0 += matrix[ix][jx];
      sum1 += matrix[ix][jx+1];
      sum2 += matrix[ix][jx+2];
      sum3 += matrix[ix][jx+3];
      sum4 += matrix[ix][jx+4];
      sum5 += matrix[ix][jx+5];
      sum6 += matrix[ix][jx+6];
      sum7 += matrix[ix][jx+7];
    }
    sums[ix] = sum0 + sum1 + sum2 + sum3 + sum4 + sum5 + sum6 + sum7;
  }
}

void row_sums_unrolled4_array(double* sums, const double** matrix, size_t nrs,
    size_t ncs) {
  for (size_t ix = 0; ix < nrs; ++ix) {
    double sum[4] = {0.};
    for (size_t jx = 0; jx < ncs; jx += 4) {
      sum[0] += matrix[ix][jx];
      sum[1] += matrix[ix][jx+1];
      sum[2] += matrix[ix][jx+2];
      sum[3] += matrix[ix][jx+3];
    }
    for (int i = 0; i < 4; ++i) {
      sums[ix] += sum[i];
    }
  }
}

int main(int argc, char const *argv[]) {
  /**
   Write a C program, called data_dependency, that
  
   * creates 1000 x 1000 doubles,
   * stores them in a matrix in row major order,
   * computes the row sums.
  */
  srand(time(NULL));
  const int size = 1000;
  const int num_tests = 5000;
  double* asentries = (double*) malloc(sizeof(double) * size * size);
  double** as = (double**) malloc(sizeof(double*) * size);
  double* rsums = (double*) malloc(sizeof(double) * size);
  for (size_t ix = 0, jx = 0; ix < size; ++ix, jx += size) {
    as[ix] = asentries + jx;
  }
  for (size_t ix = 0; ix < size; ++ix) {
    for (size_t jx = 0; jx < size; ++jx) {
      as[ix][jx] = rand();
    }
  }
  // Start computation
  printf("Performing sum of rows...\n");
  clock_t start = clock();
  for (int t = 0; t < num_tests; ++t) {
    row_sums(rsums, (const double**)as, size, size);
  }
  clock_t end = clock();
  float tot_time = (float)(end - start) / CLOCKS_PER_SEC;
  // Printing
  printf("[ROWS] Random element: %f\n", rsums[rand() % size]);
  printf("[ROWS] Total time:     %f [s]\n", tot_time);
  printf("[ROWS] Time per iter:  %f [s]\n", tot_time / (float) (num_tests));
  // Start computation
  printf("Performing sum of rows...\n");
  start = clock();
  for (int t = 0; t < num_tests; ++t) {
    row_sums_unrolled2(rsums, (const double**)as, size, size);
  }
  end = clock();
  tot_time = (float)(end - start) / CLOCKS_PER_SEC;
  // Printing
  printf("[ROWS - UNROLLEDx2] Random element: %f\n", rsums[rand() % size]);
  printf("[ROWS - UNROLLEDx2] Total time:     %f [s]\n", tot_time);
  printf("[ROWS - UNROLLEDx2] Time per iter:  %f [s]\n", tot_time / (float) (num_tests));
  // Start computation
  printf("Performing sum of rows...\n");
  start = clock();
  for (int t = 0; t < num_tests; ++t) {
    row_sums_unrolled4(rsums, (const double**)as, size, size);
  }
  end = clock();
  tot_time = (float)(end - start) / CLOCKS_PER_SEC;
  // Printing
  printf("[ROWS - UNROLLEDx4] Random element: %f\n", rsums[rand() % size]);
  printf("[ROWS - UNROLLEDx4] Total time:     %f [s]\n", tot_time);
  printf("[ROWS - UNROLLEDx4] Time per iter:  %f [s]\n", tot_time / (float) (num_tests));
  // Start computation
  printf("Performing sum of rows...\n");
  start = clock();
  for (int t = 0; t < num_tests; ++t) {
    row_sums_unrolled8(rsums, (const double**)as, size, size);
  }
  end = clock();
  tot_time = (float)(end - start) / CLOCKS_PER_SEC;
  // Printing
  printf("[ROWS - UNROLLEDx8] Random element: %f\n", rsums[rand() % size]);
  printf("[ROWS - UNROLLEDx8] Total time:     %f [s]\n", tot_time);
  printf("[ROWS - UNROLLEDx8] Time per iter:  %f [s]\n", tot_time / (float) (num_tests));
  // Start computation
  printf("Performing sum of rows...\n");
  start = clock();
  for (int t = 0; t < num_tests; ++t) {
    row_sums_unrolled4_array(rsums, (const double**)as, size, size);
  }
  end = clock();
  tot_time = (float)(end - start) / CLOCKS_PER_SEC;
  // Printing
  printf("[ROWS - UNROLLEDx4 (array)] Random element: %f\n", rsums[rand() % size]);
  printf("[ROWS - UNROLLEDx4 (array)] Total time:     %f [s]\n", tot_time);
  printf("[ROWS - UNROLLEDx4 (array)] Time per iter:  %f [s]\n", tot_time / (float) (num_tests));
  // Free memory
  free(as);
  free(asentries);
  free(rsums);
  return 0;
}