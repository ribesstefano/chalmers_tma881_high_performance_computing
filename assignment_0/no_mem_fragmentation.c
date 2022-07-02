#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

int main(int argc, char const *argv[]) {
  // Avoiding memory fragmentation:
  int size = 10000;
  // The actual value container is the "asentries" vector, whereas as is just a
  // convenient double pointer to the entries.
  int* asentries = (int*) malloc(sizeof(int) * size * size);
  int** as = (int**) malloc(sizeof(int*) * size);
  // Instead of allocating a new vector per i index, I instead assign the i-th
  // element in "as" to a specific location (i.e. pointer) in the "asentries"
  // vector.
  for (size_t ix = 0, jx = 0; ix < size; ++ix, jx += size) {
    as[ix] = asentries + jx;
  }
  // Since "as" is just a pointer to a pointer to the entries, and the
  // "asentries" contains ALL the entries, I do not "jump" and continue to
  // access the entries "sequentially" in memory.
  for (size_t ix = 0; ix < size; ++ix) {
    for (size_t jx = 0; jx < size; ++jx) {
      as[ix][jx] = 0;
    }
  }
  printf("%d\n", as[0][0]);
  free(as);
  free(asentries);
  return 0;
}