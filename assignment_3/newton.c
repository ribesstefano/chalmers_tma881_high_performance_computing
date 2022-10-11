#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <stdint.h>
#include <fcntl.h>
#include <getopt.h>
#include <math.h>
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
const double kEpsilonSquared = kEpsilon * kEpsilon;
const double kDivMinRe = -1e10;
const double kDivMaxRe = 1e10;
const double kDivMinIm = -1e10;
const double kDivMaxIm = 1e10;
const double kPi = 4.0 * atan(1.0); // 3.14159265359; // 22.0 / 7.0;
const int kBlockSize = 100; // Deprecated
const char kASCII_Char2Int = 48;
static const double kNaN = 0.0 / 0.0; // NOTE: NAN should be avail from math.h

typedef int16_t conv_t;
typedef int32_t attr_t;

// TODO: Be careful about alignment!!!
typedef struct {
  int id;
  int start_idx;
  int step;
  int* computed_lines;
  double re_start_val;
  double im_start_val;
  double re_step;
  double im_step;
  double complex* roots;
  mtx_t* mtx;
  cnd_t* cnd;
  attr_t** attractors;
  conv_t** convergences;
} thrd_args_t;

inline double absd(const double x) {
  return (x > 0.) ? x : -x;
}

inline int min(const int a, const int b) {
  return (a < b) ? a : b;
}

// bool is_below_threshold(double complex z) {
inline bool is_below_threshold(double zre, double zim) {
  // Remember: cabs() it's a norm operation, i.e. a distance calculation: |z| =
  // sqrt(re^2 + im^2)
  //
  // There are three optimizations we can exploit:
  // * triangle inequality
  // * in a rect. triangle, the hypothenuse is always less than each cathetus
  // * the square root can be avoided by comparing against a squared epsilon
  // double zre = creal(z);
  // double zim = cimag(z);
  zre = absd(zre);
  zim = absd(zim);
  if (kEpsilon < zre) {
    return false;
  }
  if (kEpsilon < zim) {
    return false;
  }
  if (zre + zim < kEpsilon) {
    return true;
  }
  return zre * zre + zim * zim < kEpsilonSquared;
}

inline void rootiter(const double complex* roots_in, const double complex z_in, attr_t* attr_final, conv_t* conv_final) {
  const attr_t kAttrDefaultVal = -1; // (attr_t)NAN;
  int attr = kAttrDefaultVal;
  int conv = 0;
  double complex z = z_in;
  while(true) {
    double zre = creal(z);
    double zim = cimag(z);
    if (zre < kDivMinRe || zre > kDivMaxRe ||
        zim < kDivMinIm || zim > kDivMaxIm) {
      attr = degree;
      break;
    }
    if (is_below_threshold(zre, zim)) {
      attr = degree + 1;
      break;
    }
    // TODO: The square norm of a complex number is the sum of two squares. When
    // computing it for a difference x - x', how can one avoid computing twice
    // the difference of the respective real and imaginary parts?
    // TODO: Julia implementation uses enumerate(), do we need to start from 1
    // then?
    for (int i = 0; i < degree; ++i) {
      double zre_diff = zre - creal(roots_in[i]);
      if (kEpsilon < absd(zre_diff)) {
        continue;
      }
      double zim_diff = zim - cimag(roots_in[i]);
      if (kEpsilon < absd(zim_diff)) {
        continue;
      }
      if (absd(zre_diff) + absd(zim_diff) < kEpsilon) {
        attr = i;
        break;
      }
      if (zre_diff * zre_diff + zim_diff * zim_diff < kEpsilonSquared) {
        attr = i;
        break;
      }
    }
    if (attr != kAttrDefaultVal) {
      break;
    }
    // TODO: The function cpow() must be removed. A switch-case on the degree
    // can help improving/optimizing the computation
    switch (degree) {
      case 1: {
        z = 1.;
        break;
      }
      case 2: {
        z = z / 2 + 1. / z / 2;
        break;
      }
      case 3: {
        z = z * (1. - 1. / 3.) + (1. / 3.) / (z * z);
        break;
      }
      case 4: {
        z = z * (1. - 1. / 4.) + (1. / 4.) / (z * z * z);
        break;
      }
      case 5: {
        double complex z2 = z * z;
        z = z * (1. - 1. / 5.) + (1. / 5.) / (z2 * z2);
        break;
      }
      case 6: {
        double complex z2 = z * z;
        z = z * (1. - 1. / 6.) + (1. / 6.) / (z2 * z2 * z);
        break;
      }
      case 7: {
        double complex z2 = z * z;
        double complex z4 = z2 * z2;
        z = z * (1. - 1. / 7.) + (1. / 7.) / (z4 * z2);
        break;
      }
      case 8: {
        double complex z2 = z * z;
        double complex z4 = z2 * z2;
        z = z * (1. - 1. / 8.) + (1. / 8.) / (z4 * z2 * z);
        break;
      }
      case 9: {
        double complex z2 = z * z;
        double complex z4 = z2 * z2;
        z = z * (1. - 1. / 9.) + (1. / 9.) / (z4 * z4);
        break;
      }
      default:
        z = z * (1. - 1. / degree) + 1. / cpow(z, (degree - 1)) / degree;
        break;
    }
    ++conv;
  }
  conv = (conv < kIterCutoff) ? conv : kIterCutoff - 1;
  *attr_final = (attr_t)attr;
  *conv_final = (conv_t)conv;
}

#if 0
typedef struct {
  mtx_t mtx;
  cnd_t cnd;
  size_t cnt;
} smphr_t;

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
#endif

int compute_thread(void* args_in) {
  double complex* loc_roots = (double complex*)malloc(degree * sizeof(double complex));
  attr_t* attractor = (attr_t*)malloc(num_lines * sizeof(attr_t));
  conv_t* convergence = (conv_t*)malloc(num_lines * sizeof(conv_t));
  // "Parse" arguments
  thrd_args_t* args = (thrd_args_t*)args_in;
  double zre;
  double zim = args->im_start_val;
  const double re_step = args->re_step;
  const double im_step = args->im_step;
  const int start_idx = args->start_idx;
  const int istep = args->step;
  const int thread_id = args->id;
  for (int i = 0; i < degree; ++i) {
    loc_roots[i] = args->roots[i];
  }
  /*
   * Main loop
   */
  for (int i = start_idx; i < num_lines; i += istep, zim += im_step) {
    // printf("Thread n.%d working on row %d (im: %f)\n", args->id, i, zim);
    // Start computation
    zre = args->re_start_val;
    for (int j = 0; j < num_lines; ++j, zre += re_step) {
      // printf("[TX %d] Working on element (%d, %d): %f + (%f i)\n", args->id, i, j, zre, zim);
      rootiter(loc_roots, zre + zim * I, &attractor[j], &convergence[j]);
    }
    // Update global, i.e. shared, variables
    mtx_lock(args->mtx);
    // Copy local to global
    for (int j = 0; j < num_lines; ++j) {
      args->attractors[i][j] = attractor[j];
      args->convergences[i][j] = convergence[j];
    }
    // Update thread status
    args->computed_lines[thread_id] = i + istep;
    mtx_unlock(args->mtx);
    cnd_signal(args->cnd);
  }
  // Free local buffers
  free(attractor);
  free(convergence);
  free(loc_roots);
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

int write_thread(void* args_in) {
//   // TODO
//   char attr_str[] = "      ";
//   char conv_str[] = "         ";
//   const int kAttrStrLen = 6;
//   const int kConvStrLen = 9;
//   int upper_bound;
//   int wr_line_idx = 0; // Written line index. In practice: lines written so far
//   while (wr_line_idx < num_lines) {
//     // If no new lines are available, we wait.
//     mtx_lock(&mtx);
//     while (1) {
//       // Extract the minimum of all status variables
//       upper_bound = num_lines;
//       for (int tx = 0; tx < num_threads; ++tx) {
//         if (upper_bound > computed_lines[tx]) {
//           upper_bound = computed_lines[tx];
//         }
//       }
//       if (upper_bound <= wr_line_idx) {
//         // We rely on spurious wake-ups, which in practice happen, but are not
//         // guaranteed.
//         cnd_wait(&cnd, &mtx);
//       } else {
//         mtx_unlock(&mtx);
//         break;
//       }
//     }
//     // Write to file lines that have been computed so far
//     while (wr_line_idx < upper_bound) {
//       for (int j = 0; j < num_lines; ++j) {
//         int iter = convergences[wr_line_idx][j];
//         int attr = attractors[wr_line_idx][j];

//         int attr1 = (attr + 1 * num_attractor / 3) % num_attractor;
//         int attr2 = (attr + 2 * num_attractor / 3) % num_attractor;
//         // Setup attractor string
//         attr_str[0] = kASCII_Char2Int + attr;
//         attr_str[2] = kASCII_Char2Int + attr1;
//         attr_str[4] = kASCII_Char2Int + attr2;
//         // Setup convergence string
//         conv_str[1] = conv_str[4] = conv_str[7] = kASCII_Char2Int + iter % 10;
//         conv_str[0] = conv_str[3] = conv_str[6] = kASCII_Char2Int + (iter / 10) % 10;
//         // Write employing fwrite
//         fwrite(attr_str, 1, kAttrStrLen, attr_fp);
//         fwrite(conv_str, 1, kConvStrLen, iter_fp);
//         // fprintf(attr_fp, "%d %d %d ", attr, attr1, attr2);
//         // fprintf(iter_fp, "%d %d %d ", iter, iter, iter);
//       }
//       fwrite("\n", 1, 1, attr_fp);
//       fwrite("\n", 1, 1, iter_fp);
//       ++wr_line_idx;
//     }
//   }
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
  // printf("kPi: %.64f\n", kPi);
  // printf("22.0 / 7.0: %.64f (diff: %.64f)\n", 22.0 / 7.0, abs(kPi - 22.0 / 7.0));
  // printf("4.0 * atan(1.0): %.64f (diff: %.64f)\n", 4.0 * atan(1.0), abs(kPi - 4.0 * atan(1.0)));
  roots = (double complex*)malloc(degree * sizeof(double complex));
  for (int dx = 0; dx < degree; ++dx) {
    roots[dx] = cexp((2. * kPi * I * dx) / ((double)degree));
    // printf("roots[%d] = %.2f %+.2fi\n", dx, creal(roots[dx]), cimag(roots[dx]));
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
  attr_t* attractors_entries = (attr_t*)malloc(sizeof(attr_t) * num_lines * num_lines);
  conv_t* convergences_entries = (conv_t*)malloc(sizeof(attr_t) * num_lines * num_lines);
  attr_t** attractors = (attr_t**)malloc(sizeof(attr_t*) * num_lines);
  conv_t** convergences = (conv_t**)malloc(sizeof(attr_t*) * num_lines);
  for (int i = 0, j = 0; i < num_lines; ++i, j += num_lines) {
    attractors[i] = attractors_entries + j;
    convergences[i] = convergences_entries + j;
  }
  int* computed_lines = (int*)malloc(num_threads * sizeof(int));
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
    thrd_args[tx].roots = roots;
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
      rootiter(roots, zre + zim * I, &attractors[(int)i][(int)j], &convergences[(int)i][(int)j]);
    }
  }
#endif
  /*
   * Write to file (dump data into images)
   */
  const int kAttrStrLen = 6;
  const int kConvStrLen = 9;
  const int kWriteBufferSize = 256;
  char* attr_wr_buffer = (char*)malloc(kAttrStrLen * kWriteBufferSize);
  char* conv_wr_buffer = (char*)malloc(kConvStrLen * kWriteBufferSize);
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
      for (int j = 0; j < num_lines; j += kWriteBufferSize) {
        for (int k = 0; j + k < min(j + kWriteBufferSize, num_lines); ++k) {
          int conv = (int)convergences[wr_line_idx][j + k];
          int attr = (int)attractors[wr_line_idx][j + k];
          int attr1 = (attr + 1 * num_attractor / 3) % num_attractor;
          int attr2 = (attr + 2 * num_attractor / 3) % num_attractor;
          // Setup attractor string
          attr_wr_buffer[k * kAttrStrLen + 0] = kASCII_Char2Int + attr;
          attr_wr_buffer[k * kAttrStrLen + 1] = ' ';
          attr_wr_buffer[k * kAttrStrLen + 2] = kASCII_Char2Int + attr1;
          attr_wr_buffer[k * kAttrStrLen + 3] = ' ';
          attr_wr_buffer[k * kAttrStrLen + 4] = kASCII_Char2Int + attr2;
          attr_wr_buffer[k * kAttrStrLen + 5] = ' ';
          // Setup convergence string
          const char conv_dec = kASCII_Char2Int + (conv / 10) % 10;
          const char conv_uni = kASCII_Char2Int + conv % 10;
          conv_wr_buffer[k * kConvStrLen + 0] = conv_dec;
          conv_wr_buffer[k * kConvStrLen + 1] = conv_uni;
          conv_wr_buffer[k * kConvStrLen + 2] = ' ';
          conv_wr_buffer[k * kConvStrLen + 3] = conv_dec;
          conv_wr_buffer[k * kConvStrLen + 4] = conv_uni;
          conv_wr_buffer[k * kConvStrLen + 5] = ' ';
          conv_wr_buffer[k * kConvStrLen + 6] = conv_dec;
          conv_wr_buffer[k * kConvStrLen + 7] = conv_uni;
          conv_wr_buffer[k * kConvStrLen + 8] = ' ';
        }
        // Get the actual size of the write buffer in order to not overshoot
        const int max_write_size = min(kWriteBufferSize, abs(num_lines - j));
        // Add carriage return if writing last buffer for the current line
        if (j + kWriteBufferSize >= num_lines) {
          attr_wr_buffer[kAttrStrLen * max_write_size - 1] = '\n';
          conv_wr_buffer[kConvStrLen * max_write_size - 1] = '\n';
        }
        // Write buffer to file using fwrite function
        fwrite(attr_wr_buffer, 1, kAttrStrLen * max_write_size, attr_fp);
        fwrite(conv_wr_buffer, 1, kConvStrLen * max_write_size, iter_fp);
      }
      ++wr_line_idx;
    }
  }
  for (int tx = 0; tx < num_threads; ++tx) {
    thrd_join(compute_threads[tx], &r);
  }
  free(compute_threads);
  free(thrd_args);
#else
  char attr_str[] = "      ";
  char conv_str[] = "         ";
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
  free(attractors_entries);
  free(convergences_entries);
  free(attractors);
  free(convergences);
  free(roots);
  free(attr_wr_buffer);
  free(conv_wr_buffer);
  fclose(attr_fp);
  fclose(iter_fp);
  mtx_destroy(&mtx);
  cnd_destroy(&cnd);
}