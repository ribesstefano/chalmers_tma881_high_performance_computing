#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

int main(int argc, char const *argv[]) {
  long int size = 1000000000;
  int* as = (int*) malloc(sizeof(int) * size);
  for (size_t ix = 0; ix < size; ++ix) {
    as[ix] = 0;
  }
  printf("%d\n", as[0]);
  free(as);
  return 0;
}