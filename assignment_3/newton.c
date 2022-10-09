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
#include <stdatomic.h>

#ifdef __STDC_NO_THREADS__ 
#include "c11threads.h"
#else
#include <threads.h>
#endif
/*
 * Global variables (to make them available to all threads)
 */
int num_threads = 1;
int num_lines = -1;
int num_attractor = -1;
int degree = -1;
double complex* roots = NULL;
atomic_bool computation_done = false;
atomic_int_least16_t next_idx = 0;

const int kMaxDegree = 9;
const int kIterCutoff = 100;
const double kEpsilon = 1e-3;
const double kEpsilonSquared = 1e-6;
const double kDivMinRe = -1e10;
const double kDivMaxRe = 1e10;
const double kDivMinIm = -1e10;
const double kDivMaxIm = 1e10;
const double kPi = 22.0 / 7.0;
const int kBlockSize = 100;
const char kASCII_Char2Int = 48;
static const double kNaN = 0.0 / 0.0; // NOTE: NAN should be avail from math.h

typedef int conv_t;
typedef int attr_t;

// TODO: Be careful about alignment!!!
typedef struct attractor_t {
  int attr;
  int iter;
} attractor_t;

typedef struct {
  mtx_t mtx;
  cnd_t cnd;
  size_t cnt;
} smphr_t;

// TODO: Be careful about alignment!!!
typedef struct {
  int id;
  int start_idx;
  int step;
  double re_start_val;
  double im_start_val;
  double re_step;
  double im_step;
  mtx_t* mtx;
  cnd_t* cnd;
  attr_t** attractors;
  conv_t** convergences;
  int* computed_lines;
  smphr_t* idx_smphr;
  smphr_t* thrd_smphr;
} thrd_args_t;

/**
 * @brief      Squared complex absolute value of a complex number
 *
 * @param[in]  z  The complex number
 *
 * @return     The squared absolute value
 */
double csqabs(const double complex z) {
  return creal(z) * creal(z) + cimag(z) * cimag(z);
}

attractor_t rootiter(double complex z) {
  attractor_t ret;
  ret.attr = (int)NAN;
  ret.iter = 0;
  // while(true) {
  while(ret.iter < 50) {
    if (creal(z) < kDivMinRe || creal(z) > kDivMaxRe ||
        cimag(z) < kDivMinIm || cimag(z) > kDivMaxIm) {
      ret.attr = degree;
      break;
    }
    // TODO: The use of cabs() can be avoided. Remember cabs() it's a norm
    // operation, i.e. a distance calculation: |z| = sqrt(re^2 + im^2)
    //
    // For instance, I can avoid computing the root by comparing against a
    // squared epsilon.
    // if (cabs(z) < kEpsilon) {
    if (csqabs(z) < kEpsilonSquared) {
      ret.attr = degree + 1;
    }
    // TODO: Julia implementation uses enumerate(), do we need to start from 1
    // then?
    for (int i = 0; i < degree; ++i) {
      // if (cabs(z - roots[i]) < kEpsilon) {
      if (csqabs(z - roots[i]) < kEpsilonSquared) {
        ret.attr = i;
        break;
      }
    }
    if (ret.attr != (int)NAN) {
      break;
    }
    // TODO: The function cpow() must be removed. A switch-case on the degree
    // can help improving/optimizing the computation
    // z = z * (1. - 1. / degree) + 1. / cpow(z, (degree - 1)) / degree;
    switch (degree) {
      case 1:
        z = 1.;
        break;
      case 2:
        z = z / 2 + 1. / z / 2;
        break;
      case 3:
        z = z * (1. - 1. / 3) + 1. / (z * z) / 3;
        break;
      case 4:
        z = z * (1. - 1. / 4) + 1. / (z * z * z) / 4;
        break;
      case 5:
        z = z * (1. - 1. / 5) + 1. / (z * z * z * z) / 5;
        break;
      case 6:
        z = z * (1. - 1. / 6) + 1. / (z * z * z * z * z) / 6;
        break;
      case 7:
        z = z * (1. - 1. / 7) + 1. / (z * z * z * z * z * z) / 7;
        break;
      case 8:
        z = z * (1. - 1. / 8) + 1. / (z * z * z * z * z * z * z) / 8;
        break;
      case 9:
        z = z * (1. - 1. / 9) + 1. / (z * z * z * z * z * z * z * z) / 9;
        break;
      default:
        z = z * (1. - 1. / degree) + 1. / cpow(z, (degree - 1)) / degree;
        break;
    }
    ret.iter++;
  }
  ret.iter = (ret.iter < kIterCutoff) ? ret.iter : kIterCutoff - 1;
  return ret;
}



int smphr_init(const size_t count, smphr_t* smphr) {
  smphr->cnt = count;
  int r = mtx_init(&smphr->mtx, mtx_plain);
  if (r == thrd_success) {
    r = cnd_init(&smphr->cnd);
  }
  return r;
}

void smphr_destroy(smphr_t* smphr) {
  mtx_destroy(&smphr->mtx);
  cnd_destroy(&smphr->cnd);
}

void smphr_release(smphr_t* smphr) {
  mtx_lock(&smphr->mtx);
  ++smphr->cnt;
  cnd_signal(&smphr->cnd);
  mtx_unlock(&smphr->mtx);
}

void smphr_release_all(const size_t count, smphr_t* smphr) {
  // TODO: Double check the use of this condition. Maybe an extra one is needed
  mtx_lock(&smphr->mtx);
  while (smphr->cnt != 0) {
    cnd_wait(&smphr->cnd, &smphr->mtx);
  }
  smphr->cnt = count;
  cnd_signal(&smphr->cnd);
  mtx_unlock(&smphr->mtx);
}

int smphr_acquire(smphr_t* smphr) {
  mtx_lock(&smphr->mtx);
  while (smphr->cnt == 0) {
    cnd_wait(&smphr->cnd, &smphr->mtx);
  }
  int idx = kBlockSize - smphr->cnt;
  if (atomic_load_explicit(&computation_done, memory_order_relaxed)) {
    idx = -1;
  }
  --smphr->cnt;
  mtx_unlock(&smphr->mtx);
  return idx;
}

int compute_thread(void* args_in) {
  thrd_args_t* args = (thrd_args_t*)args_in;
  attr_t* attractor = (attr_t*)malloc(num_lines * sizeof(attr_t));
  conv_t* convergence = (conv_t*)malloc(num_lines * sizeof(conv_t));
  /*
   * Main loop
   */
  double zim = args->im_start_val;
  const double re_step = args->re_step;
  const double im_step = args->im_step;
  for (int i = args->start_idx; i < num_lines; i += args->step, zim += im_step) {
    // printf("Thread n.%d working on row %d (im: %f)\n", args->id, i, zim);
    // Start computation
    double zre = args->re_start_val;
    for (int j = 0; j < num_lines; ++j, zre += re_step) {
      const double complex x = zre + zim * I;
      attractor_t ret = rootiter(x);
      attractor[j] = ret.attr;
      convergence[j] = ret.iter;
      // zre += re_step;
    }
    // zim += im_step;
    // Update global, i.e. shared, variables
    mtx_lock(args->mtx);
    // Copy local to global
    for (int j = 0; j < num_lines; ++j) {
      args->attractors[i][j] = attractor[j];
      args->convergences[i][j] = convergence[j];
    }
    // Update thread status
    args->computed_lines[args->id] = i + args->step;
    mtx_unlock(args->mtx);
    cnd_signal(args->cnd);
  }
  free(attractor);
  free(convergence);
  return 0;
  // printf("Hello from thread n.%d\n", args->id);
  // int idx = smphr_acquire(args->idx_smphr);
  // while(!atomic_load_explicit(&computation_done, memory_order_relaxed)) {
  //   // Simulate some processing...
  //   printf("Thread n.%d acquired index n.%d\n", args->id, idx);
  //   thrd_sleep(&(struct timespec){.tv_sec=0, .tv_nsec=50000}, NULL);

  //   idx = smphr_acquire(args->idx_smphr);
  //   if (idx < 0) {
  //     printf("Thread n.%d returning.\n", args->id);
  //     return 0;
  //   }
  // }
  // printf("Thread n.%d returning.\n", args->id);
  // return 0;
}

void get_filename(const int degree, const char* base, char* filename) {
  char degree_str[2] = {degree + kASCII_Char2Int, 0}; // Convert degree from int to ASCII
  memset(filename, 0, strlen(filename));
  strcat(filename, base);
  strcat(filename, degree_str);
  strcat(filename, ".ppm");
}

int main(int argc, char* const* argv) {
  /*
   * The recommended implementation is in the video and an overview of the
   * algorithm can be found in the Julia reference implmentation.
   *
   * Example calls:
   *
   * ./newton -t2 -l1000 5
   * ./newton -t5 -l1000 7
   * ./newton -l1000 -t5 7
   */
  srand(time(NULL));
  int opt;
  // NOTE: The colon after the argument makes it a "required" optional argument
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
  if (degree <= 0) {
    fprintf(stderr, "ERROR. Degree must be positive. Supplied: %d. Exiting\n", degree);
    exit(EXIT_FAILURE);
  }
  if (degree > kMaxDegree) {
    fprintf(stderr, "ERROR. Degree must be less than %d. Supplied: %d. Exiting\n", kMaxDegree, degree);
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
  if (num_lines >= 100000) {
    fprintf(stderr, "ERROR. Number of lines must be less than 100,000. Supplied: %d. Exiting\n", num_lines);
    exit(EXIT_FAILURE);
  }
  /*
   * Precompute roots
   */
  roots = (double complex*)malloc(degree * sizeof(double complex));
  for (int i = 0; i < degree; ++i) {
    roots[i] = cexp(2 * kPi * I * i / degree);
    // printf("roots[%d] = %.2f %+.2fi\n", i, creal(roots[i]), cimag(roots[i]));
  }
#if 0
  int r;
  smphr_t idx_smphr;
  smphr_init(kBlockSize, &idx_smphr);

  thrd_sleep(&(struct timespec){.tv_sec=1, .tv_nsec=1000}, NULL);
  printf("Main thread freeing compute threads...\n");

  // printf("Releasing semaphores...\n");
  atomic_store(&computation_done, true);
  smphr_release_all(kBlockSize, &idx_smphr);
  cnd_broadcast(&(idx_smphr.cnd));

  smphr_destroy(&idx_smphr);
#endif
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
   * Allocate arrays
   */
  const size_t kAttrItersSize = sizeof(attractor_t) * num_lines * num_lines;
  attr_t* attractors_entries = (attr_t*)malloc(kAttrItersSize);
  conv_t* convergences_entries = (conv_t*)malloc(kAttrItersSize);
  attr_t** attractors = (attr_t**)malloc(kAttrItersSize / num_lines);
  conv_t** convergences = (conv_t**)malloc(kAttrItersSize / num_lines);
  for (int i = 0, j = 0; i < num_lines; ++i, j += num_lines) {
    attractors[i] = attractors_entries + j;
    convergences[i] = convergences_entries + j;
  }
  int* computed_lines = (int*) malloc(num_threads * sizeof(int));
  /*
   * Define and init mutex and conv.var
   */
  mtx_t mtx;
  cnd_t cnd;
  mtx_init(&mtx, mtx_plain);
  cnd_init(&cnd);
  /*
   * Create and launch compute threads
   */
  // Define loop boundaries and run the algorithm
  const double kMinRange = -2.0;
  const double kMaxRange = 2.0;
  const double kStep = (kMaxRange - kMinRange) / (double)num_lines;

#define USE_MULTITHREADING 1

#if USE_MULTITHREADING
  int r;
  thrd_t* compute_threads = (thrd_t*)malloc(num_threads * sizeof(thrd_t));
  thrd_args_t* thrd_args = (thrd_args_t*)malloc(num_threads * sizeof(thrd_args_t));
  // Launch threads
  for (int tx = 0; tx < num_threads; ++tx) {
    thrd_args[tx].id = tx;
    thrd_args[tx].start_idx = tx;
    thrd_args[tx].step = num_threads;
    thrd_args[tx].re_start_val = kMinRange;
    thrd_args[tx].im_start_val = kMinRange + tx * kStep;
    thrd_args[tx].re_step = kStep;
    thrd_args[tx].im_step = num_threads * kStep;
    thrd_args[tx].mtx = &mtx;
    thrd_args[tx].cnd = &cnd;
    thrd_args[tx].attractors = attractors;
    thrd_args[tx].convergences = convergences;
    thrd_args[tx].computed_lines = computed_lines;
    r = thrd_create(&compute_threads[tx], compute_thread, (void*)(&thrd_args[tx]));
    if (r != thrd_success) {
      fprintf(stderr, "ERROR. Failed to create thread %d. Exiting.\n", tx);
      exit(1);
    }
    // thrd_detach(compute_threads[tx]); // It might work as well...
  }
#else
  /*
   * Compute algorithm's iterations
   */
  for (double zre = kMinRange, i = 0; i < num_lines; ++i, zre += kStep) {
    for (double zim = kMinRange, j = 0; j < num_lines; ++j, zim += kStep) {
      attractor_t ret = rootiter(zre + zim * I);
      attractors[(int)i][(int)j] = ret.attr;
      convergences[(int)i][(int)j] = ret.iter;
    }
  }
#endif
  /*
   * Write to file (dump data into images)
   */
  char attr_str[] = "      ";
  char conv_str[] = "         ";
  const int kAttrStrLen = 6;
  const int kConvStrLen = 9;
#if USE_MULTITHREADING
  int upper_bound;
  int wr_line_idx = 0; // Written line index. In practice: lines written so far
  while (wr_line_idx < num_lines) {
    // If no new lines are available, we wait.
    mtx_lock(&mtx);
    while (1) {
      // Extract the minimum of all status variables
      upper_bound = num_lines;
      for (int tx = 0; tx < num_threads; ++tx) {
        if (upper_bound > computed_lines[tx]) {
          upper_bound = computed_lines[tx];
        }
      }
      if (upper_bound <= wr_line_idx) {
        // We rely on spurious wake-ups, which in practice happen, but are not
        // guaranteed.
        cnd_wait(&cnd, &mtx);
      } else {
        mtx_unlock(&mtx);
        break;
      }
    }
    // Write to file lines that have been computed so far
    while (wr_line_idx < upper_bound) {
      for (int j = 0; j < num_lines; ++j) {
        int iter = convergences[wr_line_idx][j];
        int attr = attractors[wr_line_idx][j];

        int attr1 = (attr + 1 * num_attractor / 3) % num_attractor;
        int attr2 = (attr + 2 * num_attractor / 3) % num_attractor;
        // Setup attractor string
        attr_str[0] = kASCII_Char2Int + attr;
        attr_str[2] = kASCII_Char2Int + attr1;
        attr_str[4] = kASCII_Char2Int + attr2;
        // Setup convergence string
        conv_str[1] = conv_str[4] = conv_str[7] = kASCII_Char2Int + iter % 10;
        conv_str[0] = conv_str[3] = conv_str[6] = kASCII_Char2Int + (iter / 10) % 10;
        // Write employing fwrite
        fwrite(attr_str, 1, kAttrStrLen, attr_fp);
        fwrite(conv_str, 1, kConvStrLen, iter_fp);
        // fprintf(attr_fp, "%d %d %d ", attr, attr1, attr2);
        // fprintf(iter_fp, "%d %d %d ", iter, iter, iter);
      }
      fwrite("\n", 1, 1, attr_fp);
      fwrite("\n", 1, 1, iter_fp);
      ++wr_line_idx;
    }
  }
  for (int tx = 0; tx < num_threads; ++tx) {
    thrd_join(compute_threads[tx], &r);
  }
  free(compute_threads);
  free(thrd_args);
#else
  for (int i = 0; i < num_lines; ++i) {
    for (int j = 0; j < num_lines; ++j) {
      int attr = attractors[i][j];
      int iter = convergences[i][j];
      int attr1 = (attr + 1 * num_attractor / 3) % num_attractor;
      int attr2 = (attr + 2 * num_attractor / 3) % num_attractor;
      // Setup attractor string
      attr_str[0] = kASCII_Char2Int + attr;
      attr_str[2] = kASCII_Char2Int + attr1;
      attr_str[4] = kASCII_Char2Int + attr2;
      // Setup convergence string
      conv_str[1] = conv_str[4] = conv_str[7] = kASCII_Char2Int + iter % 10;
      conv_str[0] = conv_str[3] = conv_str[6] = kASCII_Char2Int + (iter / 10) % 10;
      // Write employing fwrite
      fwrite(attr_str, 1, kAttrStrLen, attr_fp);
      fwrite(conv_str, 1, kConvStrLen, iter_fp);
      // fprintf(attr_fp, "%d %d %d ", attr, attr1, attr2);
      // fprintf(iter_fp, "%d %d %d ", iter, iter, iter);
    }
    fwrite("\n", 1, 1, attr_fp);
    fwrite("\n", 1, 1, iter_fp);
  }
#endif
  /*
   * Close files and cleanup (free arrays and destroy synchronization variable)
   */
  free(attractors);
  free(convergences);
  free(attractors_entries);
  free(convergences_entries);
  free(roots);
  fclose(attr_fp);
  fclose(iter_fp);
  mtx_destroy(&mtx);
  cnd_destroy(&cnd);
}