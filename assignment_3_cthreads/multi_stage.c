#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "c11threads.h"

typedef struct {
  int val;
  char pad[60]; // cacheline - sizeof(int)
} int_padded;

typedef struct {
  const float **v;
  float **w;
  int ib;
  int istep;
  int sz;
  int tx;
  mtx_t *mtx;
  cnd_t *cnd;
  int_padded *status;
} thrd_info_t;

typedef struct {
  const float **v;
  float **w;
  int sz;
  int nthrds;
  mtx_t *mtx;
  cnd_t *cnd;
  int_padded *status;
} thrd_info_check_t;

int main_thrd(void *args) {
  /*
   * Get thread arguments
   */
  const thrd_info_t *thrd_info = (thrd_info_t*) args;
  const float **v = thrd_info->v;
  float **w = thrd_info->w;
  const int ib = thrd_info->ib;
  const int istep = thrd_info->istep;
  const int sz = thrd_info->sz;
  const int tx = thrd_info->tx;
  mtx_t *mtx = thrd_info->mtx;
  cnd_t *cnd = thrd_info->cnd;
  int_padded *status = thrd_info->status;
  /*
   * Main loop
   */
  for (int i = ib; i < sz; i += istep) {
    if (tx == 0) {
      printf("Thread n.%d performing iteration %d\n", tx, i);
    }
    const float *v_local = v[i];
    // We allocate the rows of the result before computing, and free them in
    // another thread.
    float *w_local = (float*) malloc(sz * sizeof(float));
    // Start computation
    for (int j = 0; j < sz; ++j) {
      w_local[j] = sqrtf(v_local[j]);
    }
    // Update global, i.e. shared, variables
    mtx_lock(mtx);
    w[i] = w_local;
    status[tx].val = i + istep;
    mtx_unlock(mtx);
    cnd_signal(cnd);
    // In order to illustrate thrd_sleep and to force more synchronization
    // points, we sleep after each line for one micro seconds.
    thrd_sleep(&(struct timespec){.tv_sec=0, .tv_nsec=1000}, NULL);
  }
  return 0;
}

int main_thrd_check(void *args) {
  /*
   * Get thread arguments
   */
  const thrd_info_check_t *thrd_info = (thrd_info_check_t*) args;
  const float **v = thrd_info->v;
  float **w = thrd_info->w;
  const int sz = thrd_info->sz;
  const int nthrds = thrd_info->nthrds;
  mtx_t *mtx = thrd_info->mtx;
  cnd_t *cnd = thrd_info->cnd;
  int_padded *status = thrd_info->status;
  // Difference threshold
  const float eps = 1e-1;
  // We do not increment ix in this loop, but in the inner one.
  // for (int ix = 0, ibound; ix < sz;) {
  int ibound;
  int i = 0;
  while (i < sz) {
    // If no new lines are available, we wait.
    // for (mtx_lock(mtx); ;) {
    mtx_lock(mtx);
    while (1) {
      // We extract the minimum of all status variables.
      ibound = sz;
      for (int tx = 0; tx < nthrds; ++tx) {
        if (ibound > status[tx].val) {
          ibound = status[tx].val;
        }
      }
      if (ibound <= i) {
        // We rely on spurious wake-ups, which in practice happen, but are not
        // guaranteed.
        cnd_wait(cnd, mtx);
      } else {
        mtx_unlock(mtx);
        break;
      }
      // Instead of employing a conditional variable, we could also invoke
      // thrd_yield or thrd_sleep in order to yield to other threads or grant a
      // specified time to the computation threads.
    }
    fprintf(stderr, "checking until %i\n", ibound);
    // We do not initialize i in this loop, but in the outer one.
    // for (; i < ibound; ++i) {
    while (i < ibound) {
      // We only check the last element in w, since we want to illustrate the
      // situation where the check thread completes earlier than the computaton
      // threads.
      int j = sz - 1;
      float diff = v[i][j] - w[i][j] * w[i][j];
      if (diff < -eps || diff > eps) {
        fprintf(stderr, "incorrect computation at %i %i: %f %f %f\n",
            i, j, diff, v[i][j], w[i][j]);
        // This exists the whole program, including all other threads.
        exit(1);
      }
      // We free the component of w, since it will never be used again.
      free(w[i]);
      ++i;
    }
  }
  return 0;
}

int main() {
  /*
   * Allocate data
   */
  const int sz = 500;
  float **v = (float**) malloc(sz*sizeof(float*));
  float **w = (float**) malloc(sz*sizeof(float*));
  float *ventries = (float*) malloc(sz*sz*sizeof(float));
  // The entries of w will be allocated in the computation threads are freed in
  // the check thread.
  for (int ix = 0, j = 0; ix < sz; ++ix, j += sz) {
    v[ix] = ventries + j;
  }
  for (int ix = 0; ix < sz*sz; ++ix) {
    ventries[ix] = ix;
  }
  /*
   * Define threads variables
   */
  const int nthrds = 8;
  thrd_t thrds[nthrds];
  thrd_t thrd_check;
  // Threads arguments
  thrd_info_t thrds_info[nthrds];
  thrd_info_check_t thrd_info_check;
  /*
   * Define and init mutex and conv.var
   */
  mtx_t mtx;
  cnd_t cnd;
  mtx_init(&mtx, mtx_plain);
  cnd_init(&cnd);
  /*
   * Create and start compute threads
   */
  int_padded status[nthrds];
  for (int tx = 0; tx < nthrds; ++tx) {
    // Setup thread arguments
    thrds_info[tx].v = (const float**) v;
    thrds_info[tx].w = w;
    thrds_info[tx].ib = tx;
    thrds_info[tx].istep = nthrds;
    thrds_info[tx].sz = sz;
    thrds_info[tx].tx = tx;
    thrds_info[tx].mtx = &mtx;
    thrds_info[tx].cnd = &cnd;
    thrds_info[tx].status = status;
    status[tx].val = -1;
    // Actual creation
    int r = thrd_create(thrds+tx, main_thrd, (void*) (thrds_info+tx));
    if (r != thrd_success) {
      fprintf(stderr, "failed to create thread\n");
      exit(1);
    }
    thrd_detach(thrds[tx]);
  }
  /*
   * Check Thread
   */
  thrd_info_check.v = (const float**) v;
  thrd_info_check.w = w;
  thrd_info_check.sz = sz;
  thrd_info_check.nthrds = nthrds;
  thrd_info_check.mtx = &mtx;
  thrd_info_check.cnd = &cnd;
  // It is important that we have initialize status in the previous for-loop,
  // since it will be consumed by the check threads.
  thrd_info_check.status = status;
  int r = thrd_create(&thrd_check, main_thrd_check, (void*) (&thrd_info_check));
  if (r != thrd_success) {
    fprintf(stderr, "failed to create thread\n");
    exit(1);
  }
  // Join check thread
  thrd_join(thrd_check, &r);
  // Free up memory
  free(ventries);
  free(v);
  free(w);
  // Destroy mutex and cond. var
  mtx_destroy(&mtx);
  cnd_destroy(&cnd);
  return 0;
}