#include <stdio.h>
#include <stdbool.h>
#include <time.h>

int main(int argc, char* const* argv) {
  const int num_tests = 1000;
  const int num_iter = 1000000;
  int sum = 0;
  clock_t start = clock();
  for (int t = 0; t < num_tests; ++t) {
    sum = 0;
    for (int i = 1; i <= num_iter; ++i) {
      sum += i;
    }
  }
  clock_t end = clock();
  float tot_time = (float)(end - start) / CLOCKS_PER_SEC;
  printf("Sum:           %d\n", sum);
  printf("Total time:    %f\n", tot_time);
  printf("Time per iter: %f\n", tot_time / (float) (num_tests));
  return 0;
}