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
#include <threads.h>
/*
 * Global variables (to make them available to all threads)
 */
int num_threads = 1;
int num_lines = -1;
int num_attractor = -1;
int degree = -1;
double complex* roots = NULL;
/*
 * Global constants
 */
const int kMaxDegree = 9;
const int kIterCutoff = 50;
const double kEpsilon = 1e-3;
const double kEpsilonSquared = kEpsilon * kEpsilon;
const double kDivMinRe = -1e10;
const double kDivMaxRe = 1e10;
const double kDivMinIm = -1e10;
const double kDivMaxIm = 1e10;
const double kPi = 4.0 * atan(1.0); // 3.14159265359; // 22.0 / 7.0;
const char kASCII_Char2Int = 48;
/*
 * Types
 */
typedef int8_t conv_t;
typedef int8_t attr_t;

// TODO: Be careful about alignment! Or not, its just the threads arguments...
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

/**
 * @brief      Determines if complex number magnitude is below threshold.
 *
 * @param[in]  zre_in  The real part
 * @param[in]  zim_in  The imaginary part
 *
 * @return     True if magnitude below threshold, False otherwise.
 */
inline bool is_below_threshold(const double zre_in, const double zim_in) {
  // NOTE: I'm exploiting the triangle inequality to eventually skip computing
  // the actual magnitude, which is done via the (expensive) norm operation.
  //
  // In particular:
  // * |re| < |z|        => if: eps < |re| => eps < |z|
  // * |im| < |z|        => if: eps < |im| => eps < |z|
  // * |z| < |re| + |im| => if: |re| + |im| < eps => |z| < eps
  const double zre = absd(zre_in);
  if (kEpsilon < zre) {
    return false;
  }
  const double zim = absd(zim_in);
  if (kEpsilon < zim) {
    return false;
  }
  if (zre + zim < kEpsilon) {
    return true;
  }
  return zre * zre + zim * zim < kEpsilonSquared;
}

inline void rootiter(const double complex* roots_in,
    const double complex z_in, attr_t* attr_final, conv_t* conv_final) {
  const attr_t kAttrDefaultVal = -1;
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
      // NOTE: Exploiting triangle inequality like written in function
      // is_below_threshold().
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
    // NOTE: The function cpow() can be removed. A switch-case on the degree can
    // help improving/optimizing the computation by utilizing at most 5 mul
    // operations.
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
        const double complex z2 = z * z;
        z = z * (1. - 1. / 5.) + (1. / 5.) / (z2 * z2);
        break;
      }
      case 6: {
        const double complex z2 = z * z;
        z = z * (1. - 1. / 6.) + (1. / 6.) / (z2 * z2 * z);
        break;
      }
      case 7: {
        const double complex z2 = z * z;
        const double complex z4 = z2 * z2;
        z = z * (1. - 1. / 7.) + (1. / 7.) / (z4 * z2);
        break;
      }
      case 8: {
        const double complex z2 = z * z;
        const double complex z4 = z2 * z2;
        z = z * (1. - 1. / 8.) + (1. / 8.) / (z4 * z2 * z);
        break;
      }
      case 9: {
        const double complex z2 = z * z;
        const double complex z4 = z2 * z2;
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

int compute_thread(void* args_in) {
  // Allocate local arrays.
  // TODO: Check if local copy of roots is actually necessary...
  double complex* loc_roots = (double complex*)malloc(degree * sizeof(double complex));
  attr_t* attractor = (attr_t*)malloc(num_lines * sizeof(attr_t));
  conv_t* convergence = (conv_t*)malloc(num_lines * sizeof(conv_t));
  // "Parse" arguments
  thrd_args_t* args = (thrd_args_t*)args_in;
  double zre;
  double zim = args->im_start_val;
  const double re_start = args->re_start_val;
  const double re_step = args->re_step;
  const double im_step = args->im_step;
  const int start_idx = args->start_idx;
  const int istep = args->step;
  const int thread_id = args->id;
  // Copy to local roots array
  for (int i = 0; i < degree; ++i) {
    loc_roots[i] = args->roots[i];
  }
  /*
   * Main loop
   */
  for (int i = start_idx; i < num_lines; i += istep, zim += im_step) {
    zre = re_start;
    for (int j = 0; j < num_lines; ++j, zre += re_step) {
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
}

int write_thread(void* args_in) {
  // TODO
}

void get_filename(const int degree, const char* base, char* filename) {
  // Convert degree from int to ASCII
  char degree_str[2] = {degree + kASCII_Char2Int, 0};
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
  for (int dx = 0; dx < degree; ++dx) {
    roots[dx] = cexp((2. * kPi * I * dx) / ((double)degree));
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
   * Allocate arrays
   */
  attr_t* attractors_entries = (attr_t*)malloc(sizeof(attr_t) * num_lines * num_lines);
  conv_t* convergences_entries = (conv_t*)malloc(sizeof(conv_t) * num_lines * num_lines);
  attr_t** attractors = (attr_t**)malloc(sizeof(attr_t*) * num_lines);
  conv_t** convergences = (conv_t**)malloc(sizeof(conv_t*) * num_lines);
  for (int i = 0, j = 0; i < num_lines; ++i, j += num_lines) {
    attractors[i] = attractors_entries + j;
    convergences[i] = convergences_entries + j;
  }
  int* computed_lines = (int*)malloc(num_threads * sizeof(int));
  /*
   * Define and init mutex and conditional variable
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
  int r;
  thrd_t* compute_threads = (thrd_t*)malloc(num_threads * sizeof(thrd_t));
  thrd_args_t* thargs = (thrd_args_t*)malloc(num_threads * sizeof(thrd_args_t));
  // Launch threads
  for (int tx = 0; tx < num_threads; ++tx) {
    thargs[tx].id = tx;
    thargs[tx].start_idx = tx;
    thargs[tx].step = num_threads;
    thargs[tx].re_start_val = kMinRange;
    thargs[tx].im_start_val = kMinRange + tx * kStep;
    thargs[tx].re_step = kStep;
    thargs[tx].im_step = num_threads * kStep;
    thargs[tx].mtx = &mtx;
    thargs[tx].cnd = &cnd;
    thargs[tx].attractors = attractors;
    thargs[tx].convergences = convergences;
    thargs[tx].computed_lines = computed_lines;
    thargs[tx].roots = roots;
    r = thrd_create(&compute_threads[tx], compute_thread, (void*)(&thargs[tx]));
    if (r != thrd_success) {
      fprintf(stderr, "ERROR. Failed to create thread %d. Exiting.\n", tx);
      exit(1);
    }
    // thrd_detach(compute_threads[tx]); // It might work as well...
  }
  /*
   * Write to file (dump data into images)
   */
  const int kAttrStrLen = 6; // Something like: "X Y Z "
  const int kConvStrLen = 9; // Something like: "XX YY ZZ "
  const int kWriteBufferSize = 256;
  char* attr_wr_buffer = (char*)malloc(kAttrStrLen * kWriteBufferSize);
  char* conv_wr_buffer = (char*)malloc(kConvStrLen * kWriteBufferSize);
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
  free(thargs);
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