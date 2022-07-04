#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void mul_cpx(double* a_re, double* a_im, double* b_re, double* b_im,
    double* c_re, double* c_im) {
  // NOTE: (a + ib) * (c + di) = (ac - bd) + i (ad + bc)
  *a_re = (*b_re) * (*c_re) - (*b_im) * (*c_im);
  *a_im = (*b_re) * (*c_im) + (*b_im) * (*c_re);
}

int main(int argc, char const *argv[]) {
  /**
   Write a program, called same_file, that

   * contains the declaration and definition of mul_cpx in the same file as the
     main function.
   * generates vectors as_re, as_im, bs_re, bs_im, cs_im, and cs_im of doubles,
     each of length 30,000,
   * generates entries for bs_re, bs_im, cs_re, and cs_im (any kind of entries
     will do), and then multiplies the entries of bs and cs using mul_cpx,
     saving results to as (thus call mul_\cpx 30,000 times).
   * print one random entry of as to the terminal.
  */
  srand(time(NULL));
  int num_elem = 30000;
  int num_tests = 200000;
  double* as_re = (double*) malloc(sizeof(double) * num_elem);
  double* as_im = (double*) malloc(sizeof(double) * num_elem);
  double* bs_re = (double*) malloc(sizeof(double) * num_elem);
  double* bs_im = (double*) malloc(sizeof(double) * num_elem);
  double* cs_re = (double*) malloc(sizeof(double) * num_elem);
  double* cs_im = (double*) malloc(sizeof(double) * num_elem);
  // Initialize vectors
  for (int i = 0; i < num_elem; ++i) {
    bs_re[i] = rand();
    bs_im[i] = rand();
    cs_re[i] = rand();
    cs_im[i] = rand();
  }
  // Start computation
  clock_t start = clock();
  for (int t = 0; t < num_tests; ++t) {
    for (int i = 0; i < num_elem; ++i) {
      mul_cpx(&as_re[i], &as_im[i], &bs_re[i], &bs_im[i], &cs_re[i], &cs_im[i]);
    }
  }
  clock_t end = clock();
  float tot_time = (float)(end - start) / CLOCKS_PER_SEC;
  // Printing
  int rand_idx = rand() % num_elem;
  printf("Random element: %f + %f i\n", as_re[rand_idx], as_im[rand_idx]);
  printf("Total time:     %f\n", tot_time);
  printf("Time per iter:  %f\n", tot_time / (float) (num_tests));
  free(as_re);
  free(as_im);
  free(bs_re);
  free(bs_im);
  free(cs_re);
  free(cs_im);
  return 0;
}