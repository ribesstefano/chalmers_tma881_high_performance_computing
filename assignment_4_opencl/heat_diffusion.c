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

#define CL_TARGET_OPENCL_VERSION 300
#include <CL/cl.h>

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
  /*
   * OpenCL Boilerplate
   */
  cl_int error;
  // Get platform (kinda the driver...)
  cl_platform_id platform_id;
  cl_uint nmb_platforms;
  if (clGetPlatformIDs(1, &platform_id, &nmb_platforms) != CL_SUCCESS) {
    fprintf(stderr, "ERROR. Cannot get platform. Exiting\n");
    return 1;
  }
  // Get device
  cl_device_id device_id;
  cl_uint nmb_devices;
  if (clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_GPU, 1, &device_id, &nmb_devices) != CL_SUCCESS) {
    fprintf(stderr, "ERROR. Cannot get device. Exiting\n");
    return 1;
  }
  // Get context
  cl_context context;
  cl_context_properties properties[] = {
    CL_CONTEXT_PLATFORM,
    (cl_context_properties) platform_id,
    0
  };
  context = clCreateContext(properties, 1, &device_id, NULL, NULL, &error);
  if (error != CL_SUCCESS) {
    fprintf(stderr, "ERROR. Cannot create context. Exiting\n");
    return 1;
  }
  // Create Command Queue
  cl_command_queue command_queue;
  command_queue = clCreateCommandQueueWithProperties(context, device_id, NULL,
    &error);
  if (error != CL_SUCCESS) {
    fprintf(stderr, "ERROR. Cannot create command queue. Exiting\n");
    return 1;
  }
  // Get kernel code from file
  // Note: We need to append the terminating string character
  FILE *clfp = fopen("./kernel.cl", "r");
  if (clfp == NULL ) {
    fprintf(stderr, "ERROR. Could not load cl source code. Exiting\n");
    return 1;
  }
  fseek(clfp, 0, SEEK_END);
  const int kClFileSize = ftell(clfp);
  char* kernel_src = (char*)malloc((kClFileSize + 1) * sizeof(char));
  fseek(clfp, 0, SEEK_SET);
  fread(kernel_src, sizeof(char), kClFileSize, clfp);
  kernel_src[kClFileSize] = 0; // Terminate string
  fclose(clfp);
  // Create OpenCL program
  const size_t src_len = strlen(kernel_src);
  cl_program program = clCreateProgramWithSource(context, 1,
    (const char**)&kernel_src, (const size_t*) &src_len, &error);
  if (error != CL_SUCCESS) {
    fprintf(stderr, "ERROR. Cannot create program. Exiting\n");
    return 1;
  }
  // Compile program
  error = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
  if (error != CL_SUCCESS) {
    fprintf(stderr, "ERROR. Cannot build program. log:\n");
    size_t log_size = 0;
    clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG, 0, NULL,
      &log_size);
    char *log = malloc(log_size*sizeof(char));
    if (log == NULL) {
      fprintf(stderr, "could not allocate memory for logging. Exiting\n");
      return 1;
    }
    clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG, log_size,
      log, NULL);
    fprintf(stderr, "%s. Exiting\n", log);
    free(log);
    return 1;
  }
  // Create required kernel
  cl_kernel kernel = clCreateKernel(program, "dot_prod_mul", &error);
  if (error != CL_SUCCESS ) {
    fprintf(stderr, "ERROR. Cannot create kernel 'dot_prod_mul'. Exiting\n");
    return 1;
  }
  // Create host memory

  // Enqueue buffers

  // Set kernel arguments

  // Enqueue work, i.e. an NDRange

  // Read buffers

  // Check finish on command queue

  // Free host memory
  free(kernel_src);

  // // Free/release buffers
  // clReleaseMemObject(input_buffer_a);
  // clReleaseMemObject(input_buffer_b);
  // clReleaseMemObject(output_buffer_c);
  // Clear OpenCL Objects
  clReleaseProgram(program);
  clReleaseKernel(kernel);
  // clReleaseCommandQueue(command_queue);
  clReleaseContext(context);
  // Close file(s)
  fclose(fp);
  return 0;
}