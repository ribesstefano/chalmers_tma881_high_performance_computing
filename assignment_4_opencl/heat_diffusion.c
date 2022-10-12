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

// #define CL_TARGET_OPENCL_VERSION 300
// #include <CL/cl.h>


int main(int argc, char* const* argv) {
  srand(time(NULL));
  int num_iter = -1;
  double diffusion_const = NAN;
  int opt;
  // NOTE: The colon after the argument makes it a "required" optional argument
  while ((opt = getopt(argc, argv, "n:d:")) != -1) {
    switch (opt) {
      case 'n':
        num_iter = atoi(optarg);
        break;
      case 'd':
        diffusion_const = atof(optarg);
        break;
      default:
        fprintf(stderr, "ERROR. Invalid argument. Usage: %s [-n N] [-d D]\n", argv[0]);
        exit(EXIT_FAILURE);
    }
  }
  if (num_iter == -1) {
    fprintf(stderr, "ERROR. Argument '-n' is mandatory. Usage: %s [-n N] [-d D]\n", argv[0]);
    exit(EXIT_FAILURE);
  }
  if (num_iter <= 0) {
    fprintf(stderr, "ERROR. Number of iterations must be positive. Supplied: %d. Exiting\n", num_iter);
    exit(EXIT_FAILURE);
  }
  if (isnan(diffusion_const)) {
    fprintf(stderr, "ERROR. Argument '-d' is mandatory. Usage: %s [-n N] [-d D]\n", argv[0]);
    exit(EXIT_FAILURE);
  }
  /*
   * Open file
   */
  FILE *fp;
  const char kDiffusionFilename[] = "diffusion";
  fp = fopen(kDiffusionFilename, "r");
  if (fp == NULL)  {
    fprintf(stderr, "ERROR. Error opening file. Exiting\n");
    exit(EXIT_FAILURE);
  }



  fclose(fp);
  return 0;
}