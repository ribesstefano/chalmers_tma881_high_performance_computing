#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <stdint.h>
#include <fcntl.h>
#include <math.h>
#include <getopt.h>
#include <complex.h>

/*
 * Global variables (to make them available to all threads)
 */
int num_threads = 1;
int num_lines = -1;
int num_attractor = -1;
int degree = -1;
double complex* roots = NULL;

const int kIterCutoff = 100;
const double kEpsilon = 1e-3;
const double kDivMinRe = -1e10;
const double kDivMaxRe = 1e10;
const double kDivMinIm = -1e10;
const double kDivMaxIm = 1e10;
const double kPi = 22.0 / 7.0;
static const double kNaN = 0.0 / 0.0; // NOTE: NAN should be avail from math.h

// TODO: Be careful about alignment!!!
typedef struct attractor_t {
  int attr;
  int iter;
} attractor_t;

void get_filename(const int degree, const char* base, char* filename) {
  char degree_str[2] = {degree + 48, 0}; // Convert degree from int to ASCII
  memset(filename, 0, strlen(filename));
  strcat(filename, base);
  strcat(filename, degree_str);
  strcat(filename, ".ppm");
}

attractor_t rootiter(double complex z) {
  attractor_t ret;
  ret.attr = (int)NAN;
  ret.iter = 0;
  // while(ret.iter < 100) {
  while(true) {
    if (creal(z) < kDivMinRe || creal(z) > kDivMaxRe ||
        cimag(z) < kDivMinIm || cimag(z) > kDivMaxIm) {
      ret.attr = degree;
      break;
    }
    if (cabs(z) < kEpsilon) {
      ret.attr = degree + 1;
    }
    // TODO: Julia implementation uses enumerate(), do we need to start from 1
    // then?
    for (int i = 0; i < degree; ++i) {
      if (cabs(z - roots[i]) < kEpsilon) {
        ret.attr = i;
        break;
      }
    }
    if (ret.attr != (int)NAN) {
      break;
    }
    z = z * (1. - 1. / degree) + 1. / cpow(z, (degree - 1)) / degree;
    ret.iter++;
  }
  ret.iter = (ret.iter < kIterCutoff) ? ret.iter : kIterCutoff - 1;
  return ret;
}

int main(int argc, char* const* argv) {
  /*
   * The recommended implementation is in the video and an overview of the
   * algorithm can be found in the Julia reference implmentation.
   */
  srand(time(NULL));
  int opt;
  // NOTE: The colon after the argument makes it a "required" optinal argument
  while ((opt = getopt(argc, argv, "t:l:")) != -1) {
    switch (opt) {
    case 't':
      num_threads = atoi(optarg);
      break;
    case 'l':
      num_lines = atoi(optarg);
      break;
    default:
      fprintf(stderr, "ERROR. Invalid argument. Usage: %s [-t N] -l L d\n", argv[0]);
      exit(EXIT_FAILURE);
    }
  }
  if (optind != argc - 1) {
    fprintf(stderr, "ERROR. Position argument 'd' is missing. Usage: %s [-t N] -l L d\n", argv[0]);
    exit(EXIT_FAILURE);
  }
  degree = atoi(argv[optind]);
  if (degree < 0) {
    fprintf(stderr, "ERROR. Degree must be positive. Supplied: %d. Exiting\n", degree);
    exit(EXIT_FAILURE);
  }
  if (num_lines == -1) {
    fprintf(stderr, "ERROR. Number of lines must be supplied. Usage: %s [-t N] -l L d\n", argv[0]);
    exit(EXIT_FAILURE);
  }
  if (num_lines <= 0) {
    fprintf(stderr, "ERROR. Number of lines must be positive. Supplied: %d. Exiting\n", num_lines);
    exit(EXIT_FAILURE);
  }
  /*
   * Precompute roots
   */
  roots = (double complex*)malloc(degree * sizeof(double complex));
  for (int i = 0; i < degree; ++i) {
    roots[i] = cexp(2 * kPi * I * i / degree);
    printf("roots[%d] = %.2f %+.2fi\n", i, creal(roots[i]), cimag(roots[i]));
  }
  /*
   * Compute algorithm's iterations
   */
  // Allocate attractor array
  const size_t kAttrItersSize = sizeof(attractor_t) * num_lines * num_lines;
  attractor_t* attr_iters_entries = (attractor_t*)malloc(kAttrItersSize);
  attractor_t** attr_iters = (attractor_t**)malloc(kAttrItersSize / num_lines);
  for (int i = 0, j = 0; i < num_lines; ++i, j += num_lines) {
    attr_iters[i] = attr_iters_entries + j;
  }
  // Define loop boundaries and run the algorithm
  const double kMinRange = -2.0;
  const double kMaxRange = 2.0;
  const double kStep = (kMaxRange - kMinRange) / (double)num_lines;
  for (double zre = kMinRange, i = 0; i < num_lines; ++i, zre += kStep) {
    for (double zim = kMinRange, j = 0; j < num_lines; ++j, zim += kStep) {
      attr_iters[(int)i][(int)j] = rootiter(zre + zim * I);
      // printf("%d, %d) rootiter(%+.2f %+.2fi): ", (int)i, (int)j, zre, zim);
      // printf("%d %d\n", (int)attr_iters[(int)i][(int)j].attr, attr_iters[(int)i][(int)j].iter);
    }
  }
  /*
   * Open Files
   */
  FILE* attr_fp;
  FILE* iter_fp;
  char filename[30] = {0};
  get_filename(degree, "newton_attractors_x", filename);
  attr_fp = fopen(filename, "w");
  get_filename(degree, "newton_convergence_x", filename);
  iter_fp = fopen(filename, "w");
  /*
   * Write headers
   */
  num_attractor = degree + 2;
  fprintf(attr_fp, "P3\n");
  fprintf(attr_fp, "%d %d\n", num_lines, num_lines);
  fprintf(attr_fp, "%d\n", num_attractor);
  fprintf(iter_fp, "P3\n");
  fprintf(iter_fp, "%d %d\n", num_lines, num_lines);
  fprintf(iter_fp, "%d\n", kIterCutoff);
  /*
   * Write to file (dump data into images)
   */
  for (int i = 0; i < num_lines; ++i) {
    for (int j = 0; j < num_lines; ++j) {
      int attr = (int)attr_iters[i][j].attr;
      int iter = attr_iters[i][j].iter;
      int attr1 = (attr + 1 * num_attractor / 3) % num_attractor;
      int attr2 = (attr + 2 * num_attractor / 3) % num_attractor;
      fprintf(attr_fp, "%d %d %d ", attr, attr1, attr2);
      fprintf(iter_fp, "%d %d %d ", iter, iter, iter);
    }
    fprintf(attr_fp, "\n");
    fprintf(iter_fp, "\n");
  }
  /*
   * Close files and cleanup
   */
  fclose(attr_fp);
  fclose(iter_fp);
  free(roots);
  free(attr_iters_entries);
  free(attr_iters);
}