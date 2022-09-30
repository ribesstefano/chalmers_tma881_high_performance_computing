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
#include <getopt.h>

// TODO: Using fixed point representation is faster, but it doesn't produce a
// nice (and expected) normal distribution of the distance counts, i.e. with
// higher counts for the distance values in the middle of the possible range [0,
// ~34.85].
#ifndef USE_FIX_POINT_CELLS
// #define USE_FIX_POINT_CELLS
#endif

// NOTE: When the precision is low, more distance values will fall
// "into the same bin", i.e. their fixed values will coincide with the same
// index. This in turn will lead to better and faster printing, avoiding many
// distance values being printed multiple times. After all, we are only printing
// the first 2 decimal places, so there is no need for a high precision.
//
// A value of 6 prevents "duplicate" values.
#ifndef FIX_POINT_FRAC_BITS
#define FIX_POINT_FRAC_BITS 6 // 10 seems to be the sweet spot
#endif

const float fix2float_scaler = 1.f / (float)(1 << FIX_POINT_FRAC_BITS);
const int kBytesPerLine = 23 + 1; // 23 character plus the return char
const int kBytesPerCell = 8; // kBytesPerLine / 3;
const char kASCII_Char2Int = 48;

typedef int16_t fix_t;
typedef uint16_t ufix_t;

typedef struct {
  fix_t x, y, z;
} cell_t;

inline int min(const int a, const int b) {
  return (a < b) ? a : b;
}

inline int max(const int a, const int b) {
  return (a > b) ? a : b;
}

inline ufix_t float2ufix(const float x) {
#ifdef USE_FIX_POINT_CELLS
  return (ufix_t)(x * (1 << FIX_POINT_FRAC_BITS));
#else
  // NOTE: For truncating the value, we add 0.5 and convert to ufix
  return (ufix_t)(x + 0.5f);
#endif
}

inline float fix2float(const fix_t x) {
  return (float)x * fix2float_scaler;
}

inline fix_t char2fix(const char* str) {
  // Expected 8 Bytes string: "-04.238 "
#ifdef USE_FIX_POINT_CELLS
  // NOTE: We are comparing against the offsetted '+' char because we offsetted
  // all the characters in `parse_lines()` function.
  const bool is_positive = str[0] + kASCII_Char2Int == '+';
  const fix_t fix = ((str[1] * 10 + str[2]) << FIX_POINT_FRAC_BITS) +
                    (str[4] >> 1) + (str[5] >> 2) + (str[6] >> 3);
  return is_positive ? fix : -fix;
#else
  fix_t intp;
  fix_t decp;
  if (str[0] + kASCII_Char2Int == '-') {
    intp = -(str[2] + str[1] * 10) * 1000;          // Integer part (scaled)
    decp = -(str[6] + str[5] * 10 + str[4] * 100);  // Decimal part
  } else {
    intp = (str[2] + str[1] * 10) * 1000;           // Integer part (scaled)
    decp = (str[6] + str[5] * 10 + str[4] * 100);   // Decimal part
  }
  return intp + decp;
#endif
}

inline void parse_lines(const int num_lines, char* line_buffer, cell_t* cells) {
  // Expected line like: "-04.238 -07.514 +08.942"
  int i;
// #pragma omp parallel for private(i) shared(line_buffer, cells)
  for (i = 0; i < num_lines * kBytesPerLine; ++i) {
    line_buffer[i] -= kASCII_Char2Int;
  }
// #pragma omp parallel for private(i) shared(line_buffer, cells)
  for (i = 0; i < num_lines; ++i) {
    cells[i].x = char2fix(&line_buffer[i * kBytesPerLine]);
    cells[i].y = char2fix(&line_buffer[i * kBytesPerLine + kBytesPerCell]);
    cells[i].z = char2fix(&line_buffer[i * kBytesPerLine + 2 * kBytesPerCell]);
  }
}

inline float cell_dist(const cell_t a, const cell_t b) {
#ifdef USE_FIX_POINT_CELLS
  float x_diff = fix2float(a.x - b.x);
  float y_diff = fix2float(a.y - b.y);
  float z_diff = fix2float(a.z - b.z);
  return sqrtf((x_diff * x_diff) + (y_diff * y_diff) + (z_diff * z_diff));
#else
  // NOTE: The multiplication result must be placed on an extended result in
  // order to not lose precision, i.e. going in underflow.
  const int32_t x2 = (a.x - b.x) * (a.x - b.x);
  const int32_t y2 = (a.y - b.y) * (a.y - b.y);
  const int32_t z2 = (a.z - b.z) * (a.z - b.z);
  return sqrtf((x2 + y2 + z2) * 0.01f);
#endif
}

int main(int argc, char* const* argv) {
  /*
   * Performance goal:
   * 
   * Number of points            1e4   1e5   1e5   1e5
   * Number of threads           1     5     10    20
   * Maximal runtime in seconds  0.33  10.3  5.40  2.88
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
  omp_set_num_threads(num_omp_threads);
  /*
   * Open file
   */
  FILE *fp;
  const char kCellsFilename[] = "cells";
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
  const long int kNumCells = ftell(fp) / kBytesPerLine;
  const int kNumDistances = kNumCells * kNumCells / 2 - kNumCells / 2;
  const int kMaxNumDistances = 3465+1;
  const int kBlockSizeX = 2048;
  const int kBlockSizeY = 2048;
  const int kBufferSize = max(kBlockSizeX, kBlockSizeY);
  /*
   * Allocate dynamic arrays
   */
  // Cell block-buffers
  cell_t* base_cells = (cell_t*)malloc(sizeof(cell_t) * kBlockSizeX);
  cell_t* block_cells = (cell_t*)malloc(sizeof(cell_t) * kBlockSizeY);
  // Line buffer for reading "large chunks" from file
  char* line_buffer = (char*)malloc(kBytesPerLine * kBufferSize);
  // Distance counts, indexed by the distance values themselves
  uint32_t* counts = (uint32_t*)calloc(sizeof(uint32_t), kMaxNumDistances);
  /*
   * Check if allocation successful
   */
  if (!base_cells) {
    fprintf(stderr, "ERROR. Error while allocating 'base_cells'. Exiting\n");
    exit(EXIT_FAILURE);
  }
  if (!block_cells) {
    fprintf(stderr, "ERROR. Error while allocating 'block_cells'. Exiting\n");
    exit(EXIT_FAILURE);
  }
  if (!line_buffer) {
    fprintf(stderr, "ERROR. Error while allocating 'line_buffer'. Exiting\n");
    exit(EXIT_FAILURE);
  }
  if (!counts) {
    fprintf(stderr, "ERROR. Error while allocating 'counts'. Exiting\n");
    exit(EXIT_FAILURE);
  }
  /*
   * Main algorithm
   */
  setvbuf(fp, NULL, _IONBF, BUFSIZ); // Using no buffering is slightly faster...
  int i, j, ii, jj;
  int block_size_x, block_size_y;
  for (i = 0; i < kNumCells; i += kBlockSizeX) {
    /*
     * Read shared base cells (the columns):
     *
     * 1. Seek to the right position in the file
     * 2. Get the block boundary (to avoid reading "outside" the file)
     * 3. Read the file into a line buffer
     * 4. Parse the lines into cell entries
     */
    block_size_x = min(i + kBlockSizeX, kNumCells) - i;
    fseek(fp, i * kBytesPerLine, SEEK_SET);
    fread(line_buffer, block_size_x * kBytesPerLine, 1, fp);
    parse_lines(block_size_x, line_buffer, base_cells);
    /*
     * Scroll over all remaining cells starting from the current base (the rows)
     *
     * NOTE: Start from the next line/cell, since the distance of a cell with
     * itself is zero. To do so, j needs to be i+1.
     */
    for (j = i + 1; j < kNumCells; j += kBlockSizeY) {
      /*
       * Similarly as above, load all the file entries block-wise:
       *
       * 1. Seek to the right position in the file
       * 2. Get the block boundary (to avoid reading "outside" the file)
       * 3. Read the file into a line buffer
       * 4. Parse the lines into cell entries
       * 5. Compute distances in the block (in a multithreaded fashion)
       * 6. Update global counts of distances (with a single thread)
       */
      // Parse file into cell block
      block_size_y = min(j + kBlockSizeY, kNumCells) - j;
      fseek(fp, j * kBytesPerLine, SEEK_SET);
      fread(line_buffer, block_size_y * kBytesPerLine, 1, fp);
      parse_lines(block_size_y, line_buffer, block_cells);
      // Compute distances in parallel: each entry in the block is independent
      //
      // TODO: Check the false-sharing problem.
      //
      // NOTE: Do the accumulation in the inner-loop already via some special
      // openmp properties over the count vector. In this way, the threads will
      // utilize a local copy that will be accumulated automatically in the end
      // In practice, it's gonna be faster having local copies on counts than
      // having looping over ~4M entries in a single loop afterwards.
      //
      // NOTE: Having local copies also saves a HUGE amount of memory. In fact,
      // the program doesn't meet the memory costraints requirements at the
      // moment.
      // 
      // NOTE: Using a reduction property is the key
#pragma omp parallel for private(ii, jj) shared(i, j, base_cells, block_cells) \
  reduction(+:counts[:kMaxNumDistances])
      for (ii = 0; ii < block_size_x; ++ii) {
        for (jj = (j == i + 1) ? ii : 0; jj < block_size_y; ++jj) {
          ++counts[float2ufix(cell_dist(base_cells[ii], block_cells[jj]))];
        }
      }
    }
  }
  /*
   * Print distance values and their counts to stdout.
   */
  unsigned char str[5] = {0};
  for (ufix_t dist = 1; dist < kMaxNumDistances; ++dist) {
    if (counts[dist] != 0) {
      // NOTE: The format on the left side of the dot specifies the total number
      // of characters to use for the printout.
      // printf("%05.2f %d\n", ufix2float(dist), counts[dist]);
      str[4] = dist % 10 + kASCII_Char2Int;
      str[3] = (dist / 10) % 10 + kASCII_Char2Int;
      str[2] = '.';
      str[1] = (dist / 100) % 10 + kASCII_Char2Int;
      str[0] = (dist / 1000) % 10 + kASCII_Char2Int;
      printf("%s %d\n", str, counts[dist]);
    }
  }
  /*
   * Free all allocated memory
   */
  free(line_buffer);
  free(base_cells);
  free(block_cells);
  free(counts);
  fclose(fp);
  return 0;
}