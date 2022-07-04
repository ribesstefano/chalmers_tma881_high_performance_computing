#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

int main(int argc, char* const* argv) {
  int a, b, opt;
  bool has_a_arg = false;
  bool has_b_arg = false;
  // NOTE: The colon after the argument makes it a "required argument"
  while ((opt = getopt(argc, argv, "a:b:")) != -1) {
    switch (opt) {
    case 'a':
      has_a_arg = true;
      a = atoi(optarg);
      break;
    case 'b':
      has_b_arg = true;
      b = atoi(optarg);
      break;
    default:
      fprintf(stderr, "Usage: %s [-a A] [-b B]\n", argv[0]);
      exit(EXIT_FAILURE);
    }
  }
  if (!(has_a_arg && has_b_arg)) {
    fprintf(stderr, "Not enough arguments. Usage: %s [-a A] [-b B]\n", argv[0]);
    exit(EXIT_FAILURE);
  }
  printf("A is %d and B is %d\n", a, b);
  return 0;
}