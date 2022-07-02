#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

int main(int argc, char const *argv[]) {
  // Not avoiding memory fragmentation
  int size = 10000;
  // Here I'm allocating a NE vector at each of the as vector indexes.
  int** as = (int**) malloc(sizeof(int*) * size);
  for (size_t ix = 0; ix < size; ++ix) {
    as[ix] = (int*) malloc(sizeof(int) * size);
  }
  // Since each vector in the j-dimension is "independently allocated",
  // everytime I switch to a new index in the i-dimension, I "jump" to another
  // vector/pointer location.
  for (size_t ix = 0; ix < size; ++ix) {
    for (size_t jx = 0; jx < size; ++jx) {
      as[ix][jx] = 0;
    }
  }
  printf("%d\n", as[0][0]);
  for (size_t ix = 0; ix < size; ++ix) {
    free(as[ix]);
  }
  free(as);
  return 0;
}