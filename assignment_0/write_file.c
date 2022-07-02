#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>

int main(int argc, char const *argv[]) {
  const int size = 1000;
  const bool alloc_contiguosly = false;
  const bool use_bin_file = true;
  const char filename[] = "test.dat";
  FILE* fp;
  int* asentries;
  int** as = (int**) malloc(sizeof(int*) * size);
  // Allocate memory
  if (alloc_contiguosly) {
    asentries = (int*) malloc(sizeof(int) * size * size);
    for (size_t ix = 0, jx = 0; ix < size; ++ix, jx += size) {
      as[ix] = asentries + jx;
    }
  } else {
    for (size_t ix = 0; ix < size; ++ix) {
      as[ix] = (int*) malloc(sizeof(int) * size);
    }
  }
  for (size_t ix = 0; ix < size; ++ix) {
    for (size_t jx = 0; jx < size; ++jx) {
      as[ix][jx] = ix * jx;
    }
  }
  // Write to file
  fp = fopen(filename, "w+");
  if (use_bin_file) {
    if (alloc_contiguosly) {
      fwrite(asentries, sizeof(int), size * size, fp);
    } else {
      for (size_t ix = 0; ix < size; ++ix) {
        fwrite(as[ix], sizeof(int), size, fp);
      }
    }
  } else {
    for (size_t ix = 0; ix < size; ++ix) {
      for (size_t jx = 0; jx < size; ++jx) {
        fprintf(fp, "%d ", as[ix][jx]);
      }
      fprintf(fp, "\n");
    }
  }
  fclose(fp);
  // Free memory
  if (alloc_contiguosly) {
    free(as);
    free(asentries);
  } else {
    for (size_t ix = 0; ix < size; ++ix) {
      free(as[ix]);
    }
    free(as);
  }
  return 0;
}