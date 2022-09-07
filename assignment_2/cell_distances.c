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

#ifndef DEBUG_PRINT
#define DEBUG_PRINT
#endif

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

const float fix2float_scaler = 1.0 / (float)(1 << FIX_POINT_FRAC_BITS);
const int kBytesPerLine = 23 + 1; // 23 character plus the return char
const int kBytesPerCell = 8; // kBytesPerLine / 3;
const char kASCII_Char2Int = 48;

typedef int16_t fix_t;
typedef int32_t fixd_t;
typedef uint16_t ufix_t;
typedef uint32_t ufixd_t;

typedef struct {
#ifdef USE_FIX_POINT_CELLS
  fix_t x, y, z;
#else
  signed char x_int, y_int, z_int;
  fix_t x_dec, y_dec, z_dec;
#endif
} cell_t;

inline int min(const int a, const int b) {
  return (a < b) ? a : b;
}

inline int max(const int a, const int b) {
  return (a > b) ? a : b;
}

inline ufix_t float2ufix(const float x) {
  return (ufix_t)(x * (1 << FIX_POINT_FRAC_BITS));
}

inline float fix2float(const fix_t x) {
  // return (float)x * fix2float_scaler;
  return (float)x / (float)(1 << FIX_POINT_FRAC_BITS);
}

inline float ufix2float(const ufix_t x) {
  // return (float)x * fix2float_scaler;
  return (float)x / (float)(1 << FIX_POINT_FRAC_BITS);
}

inline void parse_number(const char* line_str, signed char* intp, fix_t* decp) {
  // Example line: "-04.238 -07.514 +08.942"
  *intp = line_str[2] + line_str[1] * 10;
  *decp = line_str[6] + line_str[5] * 10 + line_str[4] * 100;
  if (line_str[0] + kASCII_Char2Int == '-') {
    *intp = - (*intp);
    *decp = - (*decp);
  }
}

inline fix_t char2fix(const unsigned char* str) {
  // Expected string example: "-04.238 ", 8 Bytes
  // NOTE: We are comparing against the offsetted '+' char because we offsetted
  // all the characters in `parse_lines()` function.
  const bool is_positive = str[0] + kASCII_Char2Int == '+';
  const fix_t fix = ((str[1] * 10 + str[2]) << FIX_POINT_FRAC_BITS) + (str[4] >> 1) + (str[5] >> 2) + (str[6] >> 3);
  return is_positive ? fix : -fix;
}

inline void parse_lines(const int num_lines, char* line_buffer, cell_t* cells) {
  // Example line: "-04.238 -07.514 +08.942"
  int i;
// #pragma omp parallel for private(i) shared(line_buffer, cells)
  for (i = 0; i < num_lines * kBytesPerLine; ++i) {
    line_buffer[i] -= kASCII_Char2Int;
  }
// #pragma omp parallel for private(i) shared(line_buffer, cells)
  for (i = 0; i < num_lines * kBytesPerLine; i += kBytesPerLine) {
    const int cx = i / kBytesPerLine;
#ifdef USE_FIX_POINT_CELLS
    cells[cx].x = char2fix(&line_buffer[i]);
    cells[cx].y = char2fix(&line_buffer[i + kBytesPerCell]);
    cells[cx].z = char2fix(&line_buffer[i + 2 * kBytesPerCell]);
#else
    parse_number(&line_buffer[i], &cells[cx].x_int, &cells[cx].x_dec);
    parse_number(&line_buffer[i + kBytesPerCell], &cells[cx].y_int, &cells[cx].y_dec);
    parse_number(&line_buffer[i + 2 * kBytesPerCell], &cells[cx].z_int, &cells[cx].z_dec);
#endif
  }
}

#ifndef USE_FIX_POINT_CELLS
inline float cell2float(const cell_t a, const char dim) {
  float ret = 0;
  switch (dim) {
    case 'x':
      ret = (float)a.x_int + (float)a.x_dec * 0.001;
      break;
    case 'y':
      ret = (float)a.y_int + (float)a.y_dec * 0.001;
      break;
    case 'z':
      ret = (float)a.z_int + (float)a.z_dec * 0.001;
      break;
    default:
      ret = (float)a.x_int + (float)a.x_dec * 0.001;
      break;
  }
  return ret;
}
#endif

inline float cell_dist(const cell_t a, const cell_t b) {
#ifdef USE_FIX_POINT_CELLS
  fix_t x_diff = a.x - b.x;
  fix_t y_diff = a.y - b.y;
  fix_t z_diff = a.z - b.z;
  ufixd_t x_pow = (fixd_t)x_diff * (fixd_t)x_diff;
  ufixd_t y_pow = (fixd_t)y_diff * (fixd_t)y_diff;
  ufixd_t z_pow = (fixd_t)z_diff * (fixd_t)z_diff;
  ufixd_t sum_fix = x_pow + y_pow + z_pow;
  float sum = (float)sum_fix * (1.0 / (float)(1 << (FIX_POINT_FRAC_BITS * 2)));
  return sqrt(sum);

  // float x_diff = fix2float(a.x) - fix2float(b.x);
  // float y_diff = fix2float(a.y) - fix2float(b.y);
  // float z_diff = fix2float(a.z) - fix2float(b.z);
  // return sqrt((x_diff * x_diff) + (y_diff * y_diff) + (z_diff * z_diff));
#else
  float x_diff = cell2float(a, 'x') - cell2float(b, 'x');
  float y_diff = cell2float(a, 'y') - cell2float(b, 'y');
  float z_diff = cell2float(a, 'z') - cell2float(b, 'z');
  return sqrt((x_diff * x_diff) + (y_diff * y_diff) + (z_diff * z_diff));
#endif
}

int main(int argc, char* const* argv) {
  /*
   * Performance goal:
   * 
   * Number of points            1e4   1e5   1e5   1e5
   * Number of threads           1     5     10    20
   * Maximal runtime in seconds  0.33  10.3  5.40  2.88
   * 
   * Max number of cells: 2^32 = 4,294,967,296
   * Max number of distances: 2^64 = 1.8446744e+19
   * Max number of allocated cells at the time: 582,542
   * Max number of allocated distances (4+4 bytes per element):  655,360
   * Max number of allocated distances (4 bytes per element): 1,310,720
   * 
   * c1  c2  c3  ... cN
   * c2  0   d23 ... d2N
   * c3  d32 0   ... d3N
   * .
   * .
   * .
   * cN  dN2
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
  // const int kMaxNumDistances = (1 << (sizeof(fix_t) * 8)) - 1; // 2^N-1, all the representable ones

  // TODO: Which block size should I use?
  const int kBlockSizeX = max(kNumCells * 0.1, 16);
  const int kBlockSizeY = max(kNumCells * 0.1, 16);
  /*
   * Allocate dynamic arrays
   */
  // Cell block-buffers
  cell_t* base_cells = (cell_t*)malloc(sizeof(cell_t) * kBlockSizeX);
  cell_t* block_cells = (cell_t*)malloc(sizeof(cell_t) * kBlockSizeY);
  // Distance block-buffer
  ufix_t* block_d_entries = (ufix_t*)malloc(sizeof(ufix_t) * kBlockSizeX * kBlockSizeY);
  ufix_t** block_distances = (ufix_t**)malloc(sizeof(ufix_t*) * kBlockSizeX);
  // Line buffer for reading "large chunks" from file
  unsigned char* line_buffer = (unsigned char*)malloc(kBytesPerLine * max(kBlockSizeX, kBlockSizeY));
  // Distance counts, indexed by the distance values themselves
  uint32_t* counts = (uint32_t*)calloc(sizeof(uint32_t), kMaxNumDistances);
  /*
   * Check if allocation successful
   */
  if (!(base_cells && block_cells && block_d_entries && block_distances && line_buffer && counts)) {
    fprintf(stderr, "ERROR. Error while allocating arrays. Exiting\n");
    exit(EXIT_FAILURE);
  }
  /*
   * "Populate" 2D array
   */
  for (int i = 0, j = 0; i < kBlockSizeX; ++i, j += kBlockSizeY) {
    block_distances[i] = block_d_entries + j;
  }
  /*
   * Main algorithm
   */
  setvbuf(fp, NULL, _IONBF, BUFSIZ); // Using no buffering is slightly faster...
#ifdef DEBUG_PRINT
  int tot_iter = 0;
  double file_rd_start = 0;
  double file_rd_end = 0;
  double dist_calc_start = 0;
  double dist_calc_end = 0;
  double update_counts_start = 0;
  double update_counts_end = 0;
  double alg_start = omp_get_wtime();
#endif
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
#ifdef DEBUG_PRINT
      file_rd_start += omp_get_wtime();
#endif
      parse_lines(block_size_y, line_buffer, block_cells);
#ifdef DEBUG_PRINT
      file_rd_end += omp_get_wtime();
      // Compute distances over the block (each entry can be computed
      // independently)
      dist_calc_start += omp_get_wtime();
#endif
#pragma omp parallel for private(ii, jj) shared(i, j, block_distances, base_cells, block_cells)
      for (ii = 0; ii < block_size_x; ++ii) {
        for (jj = (j == i + 1) ? ii : 0; jj < block_size_y; ++jj) {
          block_distances[ii][jj] =
            // float2ufix(cell_dist(base_cells[ii], block_cells[jj]));
            (ufix_t)trunc(100 * cell_dist(base_cells[ii], block_cells[jj]));
        }
      }
#ifdef DEBUG_PRINT
      dist_calc_end += omp_get_wtime();
#endif
      // Update counts (single thread)
      for (ii = 0; ii < block_size_x; ++ii) {
        for (jj = (j == i + 1) ? ii : 0; jj < block_size_y; ++jj) {
          counts[block_distances[ii][jj]]++;
        }
      }
#ifdef DEBUG_PRINT
      tot_iter++;
#endif
    }
  }
#ifdef DEBUG_PRINT
  double alg_end = omp_get_wtime();
  double alg_tot_time = alg_end - alg_start;
#endif
  /*
   * Print distance values and their counts to stdout.
   */
#ifdef DEBUG_PRINT
  double print_start = omp_get_wtime();
#endif
  for (ufix_t dist = 1; dist < kMaxNumDistances; ++dist) {
    if (counts[dist] != 0) {
      // NOTE: The format on the left side of the dot specifies the total number
      // of characters to use for the printout.
      // printf("%05.2f %d\n", ufix2float(dist), counts[dist]);

      unsigned char str[5] = {0};
      str[4] = dist % 10 + kASCII_Char2Int;
      str[3] = (dist / 10) % 10 + kASCII_Char2Int;
      str[2] = '.';
      str[1] = (dist / 10 / 10) % 10 + kASCII_Char2Int;
      str[0] = (dist / 10 / 10 / 10) % 10 + kASCII_Char2Int;
      if (strcmp(str, "00.00") != 0) {
        printf("%s %d\n", str, counts[dist]);
      }
    }
  }

#ifdef DEBUG_PRINT
  double print_end = omp_get_wtime();
  double print_tot_time = print_end - print_start;
  // Print out statistics
  // printf("INFO. Total/expected sum: %d / %d (diff: %d)\n", total_sum, kNumDistances, total_sum - kNumDistances);
  printf("INFO. Algorithm time:  %f s\n", alg_tot_time);
  printf("INFO. Dist calc time:  %f s\n", (dist_calc_end - dist_calc_start) / tot_iter);
  printf("INFO. Update counts t: %f s\n", (update_counts_end - update_counts_start) / tot_iter);
  printf("INFO. File read time:  %f s\n", (file_rd_end - file_rd_start) / tot_iter);
  printf("INFO. Printing time:   %f s\n", print_tot_time);
#endif
  /*
   * Free all allocated memory
   */
  free(line_buffer);
  free(base_cells);
  free(block_cells);
  free(block_distances);
  free(block_d_entries);
  free(counts);
  fclose(fp);
  return 0;
}