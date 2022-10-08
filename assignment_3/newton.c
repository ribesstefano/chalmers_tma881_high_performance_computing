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
#include <threads.h>
#include <stdatomic.h>

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

const int kIterCutoff = 100;
const double kEpsilon = 1e-3;
const double kEpsilonSquared = 1e-6;
const double kDivMinRe = -1e10;
const double kDivMaxRe = 1e10;
const double kDivMinIm = -1e10;
const double kDivMaxIm = 1e10;
const double kPi = 22.0 / 7.0;
const int kBlockSize = 100;
static const double kNaN = 0.0 / 0.0; // NOTE: NAN should be avail from math.h

typedef struct {
  mtx_t mtx;
  cnd_t cnd;
  size_t cnt;
} smphr_t;

typedef struct {
  int id;
  smphr_t* idx_smphr;
  smphr_t* thrd_smphr;
} thrd_args_t;

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
  printf("Hello from thread n.%d\n", args->id);
  int idx = smphr_acquire(args->idx_smphr);
  while(!atomic_load_explicit(&computation_done, memory_order_relaxed)) {

    // Simulate some processing...
    printf("Thread n.%d acquired index n.%d\n", args->id, idx);
    thrd_sleep(&(struct timespec){.tv_sec=0, .tv_nsec=50000}, NULL);

    idx = smphr_acquire(args->idx_smphr);
    if (idx < 0) {
      printf("Thread n.%d returning.\n", args->id);
      return 0;
    }
  }
  printf("Thread n.%d returning.\n", args->id);
  return 0;
}

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
  if (degree >= 10) {
    fprintf(stderr, "ERROR. Degree must be less than 10. Supplied: %d. Exiting\n", degree);
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

#if 0
  int r;
  smphr_t idx_smphr;
  smphr_init(kBlockSize, &idx_smphr);

  thrd_t* compute_threads = (thrd_t*)malloc(num_threads * sizeof(thrd_t));
  thrd_args_t* thrd_args = (thrd_args_t*)malloc(num_threads * sizeof(thrd_args_t));

  for (int tx = 0; tx < num_threads; ++tx) {
    thrd_args[tx].id = tx;
    thrd_args[tx].idx_smphr = &idx_smphr;
    int r = thrd_create(&compute_threads[tx], compute_thread, (void*)(&thrd_args[tx]));
    if (r != thrd_success) {
      fprintf(stderr, "ERROR. Failed to create thread %d. Exiting.\n", tx);
      exit(1);
    }
    // thrd_detach(compute_threads[tx]);
  }


  thrd_sleep(&(struct timespec){.tv_sec=1, .tv_nsec=1000}, NULL);
  printf("Main thread freeing compute threads...\n");


  // printf("Releasing semaphores...\n");
  atomic_store(&computation_done, true);
  smphr_release_all(kBlockSize, &idx_smphr);
  cnd_broadcast(&(idx_smphr.cnd));

  for (int tx = 0; tx < num_threads; ++tx) {
    thrd_join(compute_threads[tx], &r);
  }

  smphr_destroy(&idx_smphr);

  free(compute_threads);
  free(thrd_args);
#else
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
      // TODO: The function fprintf() is slow, prefer using fwrite().
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
  free(attr_iters_entries);
  free(attr_iters);
  free(roots);
#endif
}