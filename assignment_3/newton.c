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

/*
 * Global variables (to make them available to all threads)
 */
int num_threads = 1;
int num_lines = 0;

int main(int argc, char* const* argv) {
  /*
   * The recommended implementation is in the video and an overview of the
   * algorithm can be found in the Julia reference implmentation.
   */
  srand(time(NULL));
  int opt;
  // NOTE: The colon after the argument makes it a "required" optinal argument
  while ((opt = getopt(argc, argv, "t:l:")) != -1) {
    switch (opt) {
    case 't':
      num_threads = atoi(optarg);
      break;
    case 'l':
      num_lines = atoi(optarg);
      break;
    default:
      fprintf(stderr, "ERROR. Invalid optional argument. Usage: %s [-t N] [-l L] d\n", argv[0]);
      exit(EXIT_FAILURE);
    }
  }
  if (optind != argc - 1) {
    fprintf(stderr, "ERROR. Position argument 'd' is missing. Usage: %s [-t N] [-l L] d\n", argv[0]);
    exit(EXIT_FAILURE);
  }
  int degree = atoi(argv[optind]);
  printf("num_threads: %d\n", num_threads);
  printf("num_lines: %d\n", num_lines);
  printf("degree: %d\n", degree);
}