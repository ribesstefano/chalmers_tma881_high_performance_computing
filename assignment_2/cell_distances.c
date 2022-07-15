#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <omp.h>
#include <time.h>
#include <stdint.h>
#include <fcntl.h>
#include <math.h>

#ifndef FIX_POINT_FRAC_BITS
#define FIX_POINT_FRAC_BITS 10
#endif

const float fix2float_scaler = 1.0 / (float)(1 << FIX_POINT_FRAC_BITS);

typedef struct {
  signed char x_int, y_int, z_int;
  int16_t x_dec, y_dec, z_dec;
} cell_t;

/*
 * The distance struct is encoded with the distance value and the distance
 * count. The value part is encoded in fixed point format, with 10 bits for the
 * decimal part. The count is limited to 2^16, i.e. 2 bytes, meaning that each
 * distinct distance can appear/be counted up to 2^16 times.
 */
typedef struct {
  uint16_t val;
  uint16_t count;
} dist_t;

void parse_number(const char* line_str, signed char* intp, int16_t* decp) {
  // Example line: "-04.238 -07.514 +08.942"
  const char kASCII_Offset = 48;
  *intp = (line_str[1] - kASCII_Offset) * 10 + line_str[2] - kASCII_Offset;
  *decp = (line_str[4] - kASCII_Offset) * 100 +
              (line_str[5] - kASCII_Offset) * 10 + line_str[6] - kASCII_Offset;
  if (line_str[0] == '-') {
    *intp = - (*intp);
    *decp = - (*decp);
  }
}

cell_t parse_line(const char* line_str) {
  // Example line: "-04.238 -07.514 +08.942"
  cell_t encoded;
  parse_number(&line_str[0], &encoded.x_int, &encoded.x_dec);
  parse_number(&line_str[8], &encoded.y_int, &encoded.y_dec);
  parse_number(&line_str[16], &encoded.z_int, &encoded.z_dec);
  return encoded;
}

void print_cell(const cell_t cell) {
  printf("%c%02d.%03d ", (cell.x_int > 0 ? '+' : '-'), abs(cell.x_int), abs(cell.x_dec));
  printf("%c%02d.%03d ", (cell.y_int > 0 ? '+' : '-'), abs(cell.y_int), abs(cell.y_dec));
  printf("%c%02d.%03d ", (cell.z_int > 0 ? '+' : '-'), abs(cell.z_int), abs(cell.z_dec));
}

float cell2float_x(const cell_t a) {
  float ret = a.x_int * 1.0 + a.x_dec * 0.001;
  return ret;
}

float cell2float_y(const cell_t a) {
  float ret = a.y_int * 1.0 + a.y_dec * 0.001;
  return ret;
}

float cell2float_z(const cell_t a) {
  float ret = a.z_int * 1.0 + a.z_dec * 0.001;
  return ret;
}

float cell2float(const cell_t a, const char dim) {
  float ret = 0;
  switch (dim) {
    case 'x':
      ret = a.x_int * 1.0 + a.x_dec * 0.001;
      break;
    case 'y':
      ret = a.y_int * 1.0 + a.y_dec * 0.001;
      break;
    case 'z':
      ret = a.z_int * 1.0 + a.z_dec * 0.001;
      break;
    default:
      ret = a.x_int * 1.0 + a.x_dec * 0.001;
      break;
  }
  return ret;
}

float cell_dist(const cell_t a, const cell_t b) {
  float x_diff = cell2float(a, 'x') - cell2float(b, 'x');
  float y_diff = cell2float(a, 'y') - cell2float(b, 'y');
  float z_diff = cell2float(a, 'z') - cell2float(b, 'z');
  float x_pow = x_diff * x_diff;
  float y_pow = y_diff * y_diff;
  float z_pow = z_diff * z_diff;
  return sqrt(x_pow + y_pow + z_pow);
}

uint16_t dist2fix(const float d) {
  return (uint16_t)(d * (1 << FIX_POINT_FRAC_BITS));
}

void print_dist(const dist_t a) {
  float dist = (float)a.val * fix2float_scaler;
  // TODO: The leading zeros are not printed. Why?
  if (dist < 10) {
    printf("0%02.2f %d\n", dist, a.count);
  } else {
    printf("%02.2f %d\n", dist, a.count);
  }
}

int compare_dist(const void* a, const void* b) {
  // Used in qsort: compare the fixed point value of two distance struct.
  return (((dist_t*)a)->val > ((dist_t*)b)->val) ? 1 : -1;
}

int main(int argc, char* const* argv) {
  /*
   * Performance goal:
   * 
   * Number of points   1e4   1e5   1e5   1e5
   * Number of threads   1   5   10  20
   * Maximal runtime in seconds  0.33  10.3  5.40  2.88
   * 
   * Max number of cells: 2^32 = 4,294,967,296
   * Max number of distances: 2^64 = 1.8446744e+19
   * Max number of allocated cells at the time: 582,542
   * Max number of allocated distances (4+4 bytes per element):  655,360
   * Max number of allocated distances (4 bytes per element): 1,310,720
   * 
   */
  srand(time(NULL));
  int opt;
  int num_omp_threads = 1;
  // NOTE: The colon after the argument makes it a "required argument"
  while ((opt = getopt(argc, argv, "t:")) != -1) {
    switch (opt) {
    case 't':
      num_omp_threads = atoi(optarg);
      break;
    default:
      fprintf(stderr, "Usage: %s [-t N]\n", argv[0]);
      exit(EXIT_FAILURE);
    }
  }
  printf("INFO. Number of OpenMP threads: %d\n", num_omp_threads);
  omp_set_num_threads(num_omp_threads);

  FILE *fp;
  const char kCellsFilename[] = "cells.txt";
  fp = fopen(kCellsFilename, "r");
  if (fp == NULL)  {
    fprintf(stderr, "ERROR. Error opening file. Exiting\n");
    exit(EXIT_FAILURE);
  }
  /*
   * Get number of cells in the file: since each line is in a standard format,
   * we can exploit the file size and divide by each line size. Go to the end of
   * the file, the divide by 23 + 1.
   *
   * TODO: Make the code working for files not ending with a carriage return in
   * the last row!
   */
  fseek(fp, 0, SEEK_END); // Go to end of file
  const int kBytesPerLine = 23 + 1; // 23 character plus the return char
  const long int kNumCells = ftell(fp) / kBytesPerLine;
  fseek(fp, 0, SEEK_SET); // Seek back to beginning of file
  char* line_buffer = (char*)malloc(kBytesPerLine);


  dist_t* distances = (dist_t*)malloc(sizeof(dist_t) * kNumCells * kNumCells);
  cell_t* cells = (cell_t*)malloc(sizeof(cell_t) * kNumCells);

  // Read all cells
  fread(line_buffer, kBytesPerLine, 1, fp);
  cell_t a = parse_line(line_buffer);
  print_cell(a);
  for (int i = 0; i < kNumCells; ++i) {
    fread(line_buffer, kBytesPerLine, 1, fp);
    cells[i] = parse_line(line_buffer);
  }
  printf("INFO. Number of cells: %d\n", kNumCells);

  // Naive implementation
  int num_distinct_dists = 0;
  for (int i = 0; i < kNumCells; ++i) {
    cell_t base = cells[i];
    for (int j = 0; j < kNumCells; ++j) {
      if (i == j) {
        continue;
      }
      uint16_t dist_fix = dist2fix(cell_dist(base, cells[j]));
      // printf("%f\n", cell_dist(base, cells[j]));
      bool dist_found = false;
      for (int k = 0; k < num_distinct_dists; ++k) {
        if (distances[k].val == dist_fix) {
          ++distances[k].count;
          dist_found = true;
          break;
        }
      }
      if (!dist_found) {
        distances[++num_distinct_dists].val = dist_fix;
        distances[++num_distinct_dists].count = 1;
      }
    }
  }

  qsort(distances, num_distinct_dists, sizeof(dist_t), compare_dist);
  for (int i = 0; i < num_distinct_dists; ++i) {
    print_dist(distances[i]);
  }
  fclose(fp);
  free(line_buffer);
  free(distances);
  free(cells);
  return 0;
}