#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main(int argc, char const *argv[]) {
  /**
   Write a C program, called locality, that

   * creates 1000 x 1000 doubles,
   * stores them in a matrix in row major order,
   * computes the row and the column sums.
  */
  srand(time(NULL));
  const size_t size = 1000000;
  const size_t size_jump = 1000;
  const size_t num_tests = 5000;
  const double a = 3.14;
  double* x = (double*) malloc(sizeof(double) * size);
  double* y = (double*) malloc(sizeof(double) * size);
  int* p = (int*) malloc(sizeof(int) * size);
  for (size_t ix = 0; ix < size; ++ix) {
    x[ix] = rand();
    y[ix] = rand();
  }
  // Linear initialization
  for (size_t ix = 0; ix < size; ++ix) {
    p[ix] = ix;
  }
  printf("Performing linear addressing accesses...\n");
  clock_t start = clock();
  for (int t = 0; t < num_tests; ++t) {
    for ( size_t kx = 0; kx < size; ++kx ) {
      size_t jx = p[kx];
      y[jx] += a * x[jx];
    }
  }
  clock_t end = clock();
  float tot_time = (float)(end - start) / CLOCKS_PER_SEC;
  // Printing
  printf("[INDIRECT] Random element: %f\n", y[rand() % size]);
  printf("[INDIRECT] Total time:     %f [s]\n", tot_time);
  printf("[INDIRECT] Time per test:  %f [s]\n", tot_time / (float) (num_tests));
  // Indirect initialization
  for (size_t jx = 0, kx = 0; jx < size_jump; ++jx) {
    for (size_t ix = jx; ix < size; ix += size_jump, ++kx) {
      p[ix] = kx;
    }
  }
  printf("\nPerforming indirect addressing accesses...\n");
  start = clock();
  for (int t = 0; t < num_tests; ++t) {
    for ( size_t kx = 0; kx < size; ++kx ) {
      size_t jx = p[kx];
      y[jx] += a * x[jx];
    }
  }
  end = clock();
  tot_time = (float)(end - start) / CLOCKS_PER_SEC;
  // Printing
  printf("[JUMPS] Random element: %f\n", y[rand() % size]);
  printf("[JUMPS] Total time:     %f [s]\n", tot_time);
  printf("[JUMPS] Time per test:  %f [s]\n", tot_time / (float) (num_tests));
  // Linear accesses
  printf("\nPerforming linear accesses...\n");
  start = clock();
  for (int t = 0; t < num_tests; ++t) {
    for ( size_t kx = 0; kx < size; ++kx ) {
      y[kx] += a * x[kx];
    }
  }
  end = clock();
  tot_time = (float)(end - start) / CLOCKS_PER_SEC;
  // Printing
  printf("[LINEAR] Random element: %f\n", y[rand() % size]);
  printf("[LINEAR] Total time:     %f [s]\n", tot_time);
  printf("[LINEAR] Time per test:  %f [s]\n", tot_time / (float) (num_tests));
  // Free memory
  free(x);
  free(y);
  free(p);
  return 0;
}