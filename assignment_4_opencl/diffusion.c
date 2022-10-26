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
   * Open file and parse problem dimensions
   */
  FILE *fp;
  const char kDiffusionFilename[] = "init";
  char * line = NULL;
  size_t len = 0;
  ssize_t read;
  int width = 0;
  int height = 0;
  fp = fopen(kDiffusionFilename, "r");
  if (fp == NULL)  {
    fprintf(stderr, "ERROR. Error opening file. Exiting\n");
    exit(EXIT_FAILURE);
  }
  // Read width and height from first line
  if ((read = getline(&line, &len, fp)) != -1) {
    sscanf(line, "%d %d", &width, &height);
  }
  double* box_hostbuffer = malloc(width * height * sizeof(double));
  if (box_hostbuffer == NULL) {
    fprintf(stderr, "ERROR. Unable to allocate box_hostbuffer. Exiting.\n");
    exit(EXIT_FAILURE);
  }
  while ((read = getline(&line, &len, fp)) != -1) {
    int x, y;
    double val;
    sscanf(line, "%d %d %lf", &x, &y, &val);
    box_hostbuffer[y * width + x] = val;
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
  char* kernel_src = malloc((kClFileSize + 1) * sizeof(char));
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
  // Create required kernels
  cl_kernel kernel_heat_diffusion;
  cl_kernel kernel_reduction;
  kernel_heat_diffusion = clCreateKernel(program, "heat_diffusion", &error);
  if (error != CL_SUCCESS) {
    fprintf(stderr, "ERROR. Cannot create kernel 'heat_diffusion'. Exiting\n");
    return 1;
  }
  kernel_reduction = clCreateKernel(program, "sum_reduction", &error);
  if (error != CL_SUCCESS) {
    fprintf(stderr, "ERROR. Cannot create kernel 'reduction'. Exiting\n");
    return 1;
  }
  /*
   * Define Problem constants
   */
  const int kBoxSize = width * height;
  const size_t kBoxByteSize = kBoxSize * sizeof(double);
  // NOTE: Global and local sizes must be multiples of each other
  const int kReductionGlobSize = 4096;
  const int kReductionLocSize = 1024;
  const int kReductionNumGroups = kReductionGlobSize / kReductionLocSize;
  const size_t kReductionLocByteSize = kReductionLocSize * sizeof(double);
  /*
   * Heat Diffusion Kernel
   */
  // Create host memory
  // Allocate box buffers for double-buffering
  // TODO: deprecated...
  cl_mem box_buffer = clCreateBuffer(context, CL_MEM_READ_WRITE, kBoxByteSize,
                                      NULL, &error);
  if (error != CL_SUCCESS) {
    fprintf(stderr, "ERROR. Cannot create box buffer. Exiting.\n");
    return 1;
  }
  // Enqueue writing to buffers to initialize them
  error = clEnqueueWriteBuffer(command_queue, box_buffer, CL_TRUE, 0,
                               kBoxByteSize, box_hostbuffer, 0, NULL, NULL);
  if (error != CL_SUCCESS) {
    fprintf(stderr, "ERROR. Cannot enqueue write of buffer n.0. Exiting.\n");
    return 1;
  }
  // Set heat_diffusion kernel arguments
  // TODO: What's the difference b/w "int" and "cl_int"??
  // cl_int is machine-dependent, sort of...
  clSetKernelArg(kernel_heat_diffusion, 0, sizeof(cl_mem), &box_buffer);
  clSetKernelArg(kernel_heat_diffusion, 1, sizeof(int), &num_iter);
  clSetKernelArg(kernel_heat_diffusion, 2, sizeof(int), &width);
  clSetKernelArg(kernel_heat_diffusion, 3, sizeof(int), &height);
  clSetKernelArg(kernel_heat_diffusion, 4, sizeof(double), &diffusion_const);
  // Enqueue work, i.e. an NDRange, on first kernel
  const size_t global_sz[] = {width, height};
  error = clEnqueueNDRangeKernel(command_queue, kernel_heat_diffusion, 2, NULL,
                                (const size_t*)&global_sz, NULL, 0, NULL, NULL);
  if (error != CL_SUCCESS) {
    fprintf(stderr, "ERROR. Cannot enqueue kernel. Exiting.\n");
    return 1;
  }
  /*
   * Reduction Kernel 
   */
  cl_mem output_buffer_sum;
  output_buffer_sum = clCreateBuffer(context, CL_MEM_WRITE_ONLY,
                                     kReductionNumGroups * sizeof(double),
                                     NULL, &error);
  if (error != CL_SUCCESS) {
    fprintf(stderr, "ERROR. Cannot create buffer for sum.\n");
    return 1;
  }
  char apply_abs_diff = 0;
  double null_avg = 0.0;
  cl_int box_size = (cl_int)kBoxSize;
  clSetKernelArg(kernel_reduction, 0, sizeof(cl_mem), &box_buffer);
  clSetKernelArg(kernel_reduction, 1, kReductionLocByteSize, NULL);
  clSetKernelArg(kernel_reduction, 2, sizeof(cl_int), &box_size);
  clSetKernelArg(kernel_reduction, 3, sizeof(bool), &apply_abs_diff);
  clSetKernelArg(kernel_reduction, 4, sizeof(double), &null_avg);
  clSetKernelArg(kernel_reduction, 5, sizeof(cl_mem), &output_buffer_sum);
  // Start NDRange
  size_t global_redsz_szt = (size_t)kReductionGlobSize;
  size_t local_redsz_szt = (size_t)kReductionLocSize;
  error = clEnqueueNDRangeKernel(command_queue, kernel_reduction, 1, NULL,
                                 (const size_t*)&global_redsz_szt,
                                 (const size_t*)&local_redsz_szt, 0, NULL,
                                 NULL);
  if (error != CL_SUCCESS) {
    fprintf(stderr, "ERROR. Cannot enqueue kernel reduction. Exiting.\n");
    return 1;
  }
  // Read work groups buffers
  double* groups_sums = malloc(kReductionNumGroups * sizeof(double));
  error = clEnqueueReadBuffer(command_queue, output_buffer_sum, CL_TRUE, 0,
                              kReductionNumGroups * sizeof(double),
                              groups_sums, 0, NULL, NULL);
  if (error != CL_SUCCESS) {
    fprintf(stderr, "ERROR. Cannot enqueue read of buffer c. Exiting.\n");
    return 1;
  }
  // Accumulate over the work group buffers
  double sum_total = 0.0;
  for (int i = 0; i < kReductionNumGroups; ++i) {
    sum_total += groups_sums[i];
  }
  const double kAverageTemperature = sum_total / kBoxSize;
  printf("%fl\n", kAverageTemperature);
  /*
   * Absolute difference Kernel
   */
  apply_abs_diff = 1;
  clSetKernelArg(kernel_reduction, 0, sizeof(cl_mem), &box_buffer);
  clSetKernelArg(kernel_reduction, 1, kReductionLocByteSize, NULL);
  clSetKernelArg(kernel_reduction, 2, sizeof(cl_int), &box_size);
  clSetKernelArg(kernel_reduction, 3, sizeof(bool), &apply_abs_diff);
  clSetKernelArg(kernel_reduction, 4, sizeof(double), &kAverageTemperature);
  clSetKernelArg(kernel_reduction, 5, sizeof(cl_mem), &output_buffer_sum);
  // Start NDRange
  clEnqueueNDRangeKernel(command_queue, kernel_reduction, 1, NULL,
                         (const size_t*)&global_redsz_szt,
                         (const size_t*)&local_redsz_szt, 0, NULL, NULL);
  if (error != CL_SUCCESS) {
    fprintf(stderr, "cannot enqueue kernel reduction\n");
    return 1;
  }
  // Read work groups buffers
  error = clEnqueueReadBuffer(command_queue, output_buffer_sum, CL_TRUE, 0,
                              kReductionNumGroups * sizeof(double),
                              groups_sums, 0, NULL, NULL);
  if (error != CL_SUCCESS) {
    fprintf(stderr, "ERROR. Cannot enqueue read of buffer output_buffer_sum\n");
    return 1;
  }
  // Accumulate over the work group buffers
  sum_total = 0.0;
  for (int i = 0; i < kReductionNumGroups; ++i) {
    sum_total += groups_sums[i];
  }
  const double kAbsoluteAverageTemperature = sum_total / kBoxSize;
  // Read box buffer, containing the final temperature values
  error = clEnqueueReadBuffer(command_queue, box_buffer, CL_TRUE, 0,
                              kBoxByteSize,
                              box_hostbuffer, 0, NULL, NULL);
  if (error != CL_SUCCESS) {
    fprintf(stderr, "ERROR. Cannot enqueue read of buffer box_buffer\n");
    return 1;
  }
  // Print final values
  for (int i = 0; i < height; ++i) {
    for (int j = 0; j < width; ++j) {
      printf("%fl ", box_hostbuffer[i * width + j]);
    }
    printf("\n");
  }
  printf("%fl\n", kAbsoluteAverageTemperature);
  // Check finish on command queue
  if (clFinish(command_queue) != CL_SUCCESS) {
    fprintf(stderr, "ERROR. Cannot finish command queue.\n");
    return 1;
  }
  // // Free/release buffers
  clReleaseMemObject(box_buffer);
  clReleaseMemObject(output_buffer_sum);
  // Clear OpenCL Objects
  clReleaseProgram(program);
  clReleaseKernel(kernel_heat_diffusion);
  clReleaseKernel(kernel_reduction);
  clReleaseCommandQueue(command_queue);
  clReleaseContext(context);
  // Free host memory
  free(kernel_src);
  free(box_hostbuffer);
  free(groups_sums);
  // Close file(s)
  fclose(fp);
  return 0;
}