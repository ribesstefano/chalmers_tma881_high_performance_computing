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

// TODO: Using fixed point representation is faster, but it doesn't produce a
// nice (and expected) normal distribution of the distance counts, i.e. with
// higher counts for the distance values in the middle of the possible range [0,
// ~34.85].
#ifndef USE_FIX_POINT_CELLS
#define USE_FIX_POINT_CELLS
#endif

// NOTE: When the precision is low, more distance values will fall
// "into the same bin", i.e. their fixed values will coincide with the same
// index. This in turn will lead to better and faster printing, avoiding many
// distance values being printed multiple times. After all, we are only printing
// the first 2 decimal places, so there is no need for a high precision.
#ifndef FIX_POINT_FRAC_BITS
#define FIX_POINT_FRAC_BITS 6 // 10 seems to be the sweet spot
#endif

const float fix2float_scaler = 1.0 / (float)(1 << FIX_POINT_FRAC_BITS);
const int kBytesPerLine = 23 + 1; // 23 character plus the return char
const int kBytesPerCell = 8; // kBytesPerLine / 3;
const char kASCII_Char2Int = 48;

typedef struct {
#ifdef USE_FIX_POINT_CELLS
  int16_t x, y, z;
#else
  signed char x_int, y_int, z_int;
  int16_t x_dec, y_dec, z_dec;
#endif
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

inline void parse_number(const char* line_str, signed char* intp, int16_t* decp) {
  // Example line: "-04.238 -07.514 +08.942"
  *intp = line_str[2] + line_str[1] * 10;
  *decp = line_str[6] + line_str[5] * 10 + line_str[4] * 100;
  if (line_str[0] == '-' - kASCII_Char2Int) {
    *intp = - (*intp);
    *decp = - (*decp);
  }
}

inline int16_t char2fix(const char* str) {
  // Expected string example: "-04.238 ", 8 Bytes
  // NOTE: We are comparing against the offsetted '+' char because we offsetted
  // all the characters in `parse_lines()` function.
  const bool is_negative = str[0] == '+' - kASCII_Char2Int;
  // if (str[1] == '1') {
  //   return is_negative ? -(10 << FIX_POINT_FRAC_BITS) : (10 << FIX_POINT_FRAC_BITS);
  // }
  const int16_t fix = ((str[1] * 10 + str[2]) << FIX_POINT_FRAC_BITS) + (str[4] >> 1) + (str[5] >> 2) + (str[6] >> 3);
  return is_negative ? -fix : fix;
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
#ifdef USE_FIX_POINT_CELLS
    cells[i / kBytesPerLine].x = char2fix(&line_buffer[i]);
    cells[i / kBytesPerLine].y = char2fix(&line_buffer[i + kBytesPerCell]);
    cells[i / kBytesPerLine].z = char2fix(&line_buffer[i + 2 * kBytesPerCell]);
#else
    parse_number(&line_buffer[i], &cells[i / kBytesPerLine].x_int, &cells[i / kBytesPerLine].x_dec);
    parse_number(&line_buffer[i + kBytesPerCell], &cells[i / kBytesPerLine].y_int, &cells[i / kBytesPerLine].y_dec);
    parse_number(&line_buffer[i + 2 * kBytesPerCell], &cells[i / kBytesPerLine].z_int, &cells[i / kBytesPerLine].z_dec);
#endif
  }
}

#ifndef USE_FIX_POINT_CELLS
void print_cell(const cell_t cell) {
  printf("%c%02d.%03d %c%02d.%03d %c%02d.%03d ",
    (cell.x_int > 0 ? '+' : '-'), abs(cell.x_int), abs(cell.x_dec),
    (cell.y_int > 0 ? '+' : '-'), abs(cell.y_int), abs(cell.y_dec),
    (cell.z_int > 0 ? '+' : '-'), abs(cell.z_int), abs(cell.z_dec));
}

inline float cell2float(const cell_t a, const char dim) {
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
#endif

inline float cell_dist(const cell_t a, const cell_t b) {
#ifdef USE_FIX_POINT_CELLS
  int16_t x_diff = a.x - b.x;
  int16_t y_diff = a.y - b.y;
  int16_t z_diff = a.z - b.z;
  int32_t x_pow = (int32_t)x_diff * (int32_t)x_diff;
  int32_t y_pow = (int32_t)y_diff * (int32_t)y_diff;
  int32_t z_pow = (int32_t)z_diff * (int32_t)z_diff;
  int32_t sum_fix = x_pow + y_pow + z_pow;
  float sum = (float)sum_fix * (1.0 / (float)(1 << (FIX_POINT_FRAC_BITS * 2)));
  return sqrt(sum);
#else
  float x_diff = cell2float(a, 'x') - cell2float(b, 'x');
  float y_diff = cell2float(a, 'y') - cell2float(b, 'y');
  float z_diff = cell2float(a, 'z') - cell2float(b, 'z');
  float x_pow = x_diff * x_diff;
  float y_pow = y_diff * y_diff;
  float z_pow = z_diff * z_diff;
  return sqrt(x_pow + y_pow + z_pow);
#endif
}

inline uint16_t float2fix(const float d) {
  return (uint16_t)(d * (1 << FIX_POINT_FRAC_BITS));
}

inline float fix2float(const uint16_t d) {
  return (float)d * fix2float_scaler;
}

void print_dist(const dist_t a) {
  float dist = (float)a.val * fix2float_scaler;
  // TODO: The leading zeros are not printed. Why?
  if (dist < 10.0) {
    printf("0%02.2f %d\n", dist, a.count);
  } else {
    printf("%02.2f %d\n", dist, a.count);
  }
}

inline uint32_t get_hash(const uint16_t din) {
  /*
   * Check this stackoverflow question:
   * https://stackoverflow.com/questions/7666509/hash-function-for-string
   *
   * and a list of implementations in C: http://www.cse.yorku.ca/%7Eoz/hash.html
   * 
   * The following implements `djb2`, adjusted for uint16_t.
   */
  uint16_t tmp = din + 314;
  char* str = (char*)(&tmp);
  uint32_t hash = 0;
  int c;
  while (c = *str++) {
    hash = c + (hash << 6) + (hash << 16) - hash;
  }
  return hash;
}

inline dist_t* search_table(const int table_size, const uint16_t dist, dist_t** hash_table, int* hash_index) {
  // Get the hash 
  int hash_idx = get_hash(dist) % table_size;
  // Move in array until a match is found
  while(hash_table[hash_idx] != NULL) {
    if(hash_table[hash_idx]->val == dist) {
      return hash_table[hash_idx];
    }
    // Go to next cell
    ++hash_idx;
    // Wrap around the table
    hash_idx %= table_size;
  }
  *hash_index = hash_idx;
  return NULL;        
}

inline void update_dist(const int table_size, const uint16_t dist, dist_t** hash_table) {
  int hash_index;
  dist_t* found_dist = search_table(table_size, dist, hash_table, &hash_index);
  if (found_dist) {
    found_dist->count++;
  } else {
    hash_table[hash_index] = (dist_t*)malloc(sizeof(dist_t));
    hash_table[hash_index]->val = dist;
    hash_table[hash_index]->count = 1;
  }
}

void insert_dist(const int table_size, const dist_t data, dist_t** hash_table) {
  dist_t* item = (dist_t*)malloc(sizeof(dist_t));
  *item = data;
  // Get the hash
  int hash_index = data.val % table_size;
  // Move in array until an empty or deleted cell
  while(hash_table[hash_index] != NULL && hash_table[hash_index]->val != 0) {
    // Go to next cell
    ++hash_index;
    // Wrap around the table
    hash_index %= table_size;
  }
  hash_table[hash_index] = item;
}

int compare_dist(const void* a, const void* b) {
  // Used in qsort: compare the value of two distance struct.
  const dist_t* a_dist = *(const dist_t**)a;
  const dist_t* b_dist = *(const dist_t**)b;
  if (a_dist == NULL) {
    return 1;
  }
  if (b_dist == NULL) {
    return -1;
  }
  return (a_dist->val > b_dist->val) ? 1 : -1;
}

inline int min(const int a, const int b) {
  return (a < b) ? a : b;
}

inline int max(const int a, const int b) {
  return (a > b) ? a : b;
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
  const long int kNumCells = ftell(fp) / kBytesPerLine;
  const int kNumDistances = kNumCells * kNumCells / 2 - kNumCells / 2;
  // NOTE: The following +1 is for the counting sort algorithm.
  const int kMaxNumDistances = 65536 + 1; // 2^16+1, all the representable ones
  const int kBlockSize = 256;
  const int kBlockSizeX = 256;
  const int kBlockSizeY = 256;
  int block_size_x, block_size_y;
  // Cell block-buffers
  cell_t* base_cells = (cell_t*)malloc(sizeof(cell_t) * kBlockSizeX);
  cell_t* block_cells = (cell_t*)malloc(sizeof(cell_t) * kBlockSizeY);
  // Distance block-buffer
  uint16_t* block_d_entries = (uint16_t*)malloc(sizeof(uint16_t) * kBlockSizeX * kBlockSizeY);
  uint16_t** block_distances = (uint16_t**)malloc(sizeof(uint16_t*) * kBlockSizeX);
  for (int i = 0, j = 0; i < kBlockSizeX; ++i, j += kBlockSizeY) {
    block_distances[i] = block_d_entries + j;
  }
  // Line buffer for reading "large chunks" from file
  char* line_buffer = (char*)malloc(kBytesPerLine * max(kBlockSizeX, kBlockSizeY));
  /*
   * Main algorithm
   */
  setvbuf(fp, NULL, _IONBF, BUFSIZ); // Using no buffering is slightly faster...

  uint16_t* counts_global = (uint16_t*)malloc(sizeof(uint16_t) * kMaxNumDistances);
  memset(counts_global, 0, sizeof(uint16_t) * kMaxNumDistances);


  int tot_iter = 0;
  double file_rd_start = 0;
  double file_rd_end = 0;
  double dist_calc_start = 0;
  double dist_calc_end = 0;
  double update_counts_start = 0;
  double update_counts_end = 0;

  int i, j, ii, jj;
  double alg_start = omp_get_wtime();
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
      file_rd_start += omp_get_wtime();
      parse_lines(block_size_y, line_buffer, block_cells);
      file_rd_end += omp_get_wtime();
      // Compute distances over the block (each entry can be computed
      // independently)
      dist_calc_start += omp_get_wtime();
#pragma omp parallel for private(ii, jj) shared(block_distances, base_cells, block_cells)
      for (ii = 0; ii < block_size_x; ++ii) {
        for (jj = (j == i + 1) ? ii : 0; jj < block_size_y; ++jj) {
          block_distances[ii][jj] =
            float2fix(cell_dist(base_cells[ii], block_cells[jj]));
        }
      }
      dist_calc_end += omp_get_wtime();
      // Update counts (single thread)
      for (ii = 0; ii < block_size_x; ++ii) {
        for (jj = (j == i + 1) ? ii : 0; jj < block_size_y; ++jj) {
          counts_global[block_distances[ii][jj]]++;
        }
      }
      tot_iter++;
    }
  }
  double alg_end = omp_get_wtime();
  double alg_tot_time = alg_end - alg_start;
#if 0
  double alg_start = omp_get_wtime();

// #pragma omp parallel private(j) // private(block_size_y, block_cells, line_buffer)
  for (i = 0; i < kNumCells; i += kBlockSizeX) {
    /*
     * Read shared base cells (the columns):
     * 
     * 1. Seek to the right position in the file
     * 2. Get the block boundary (to avoid reading "outside" the file)
     * 3. Read the file into a line buffer
     * 4. Parse the lines into cell entries
     */
#pragma omp single
    {
#pragma omp critical
      {
        block_size_x = min(i + kBlockSizeX, kNumCells) - i;
        fseek(fp, i * kBytesPerLine, SEEK_SET);
        fread(line_buffer, block_size_x * kBytesPerLine, 1, fp);
      }
      parse_lines(block_size_x, line_buffer, base_cells);
    }
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
       * 5. Compute distances in the block
       * 6. Count local block distances
       * 7. Update global counts of distances
       */
#pragma omp critical
      {
        block_size_y = min(j + kBlockSizeY, kNumCells) - j;
        file_rd_start += omp_get_wtime();
        fseek(fp, j * kBytesPerLine, SEEK_SET);
        fread(line_buffer, block_size_y * kBytesPerLine, 1, fp);
        file_rd_end += omp_get_wtime();
      }
      parse_lines(block_size_y, line_buffer, block_cells);
      /*
       * Block computation: each element can be computed in parallel
       *
       * NOTE: If it is the left-most block, then we need to skip the first
       * elements on the left.
       *
       * NOTE: OpenMP collapse clause cannot be used because the inner loop has
       * variable loop iterations!
       */
      dist_calc_start += omp_get_wtime();
// #pragma omp parallel for private(ii, jj)
      for (ii = 0; ii < block_size_x; ++ii) {
        for (jj = (j == i + 1) ? ii : 0; jj < block_size_y; ++jj) {
          block_distances[ii][jj] =
            float2fix(cell_dist(base_cells[ii], block_cells[jj]));
        }
      }
      dist_calc_end += omp_get_wtime();

      /*
       * =======================================================================
       * TODO: IMPLEMENT A PARALLEL "COUNTING SORT" OVER THE BLOCK!
       * =======================================================================
       * Ref: https://www.journaldev.com/42355/counting-sort-algorithm
       *
       * Idea: If we encode ALL possible 16-bit distance values, each paired
       * with a 16-bit counter, then we will consume around: 2^16 * 2B / (1024 *
       * 1024) ~= 0.125 MB. In fact, we are anyway truncating the distances, so
       * storing N*N-N/2 values in the hash table might actually consume more
       * than 0.125 MB.
       */
      update_counts_start += omp_get_wtime();
      uint16_t* counts_local = (uint16_t*)malloc(sizeof(uint16_t) * kMaxNumDistances);
      memset(counts_local, 0, sizeof(uint16_t) * kMaxNumDistances);
      for (ii = 0; ii < block_size_x; ++ii) {
        for (jj = (j == i + 1) ? ii : 0; jj < block_size_y; ++jj) {
          counts_local[block_distances[ii][jj]]++;
        }
      }
#pragma omp critical
      {
        for (int k = 0; k < kMaxNumDistances; ++k) {
          counts_global[k] += counts_local[k];
        }
      }
      free(counts_local);
      update_counts_end += omp_get_wtime();
      ++tot_iter; // Just for profiling purposes...
    } // end row-scrolling
  }
  double alg_end = omp_get_wtime();
  double alg_tot_time = alg_end - alg_start;
#endif
  /*
   * Print distance values and their count to stdout.
   */
  int total_sum = 0;
  double print_start = omp_get_wtime();
  for (int dist = 0; dist < kMaxNumDistances; ++dist) {
    if (counts_global[dist] != 0) {
      // const float dist = (float)i * fix2float_scaler;
      // TODO: The leading zeros are not printed. Why?
      printf("%02.2f %d\n", fix2float((uint16_t)dist), counts_global[dist]);
      total_sum += counts_global[dist];
    }
  }
  double print_end = omp_get_wtime();
  double print_tot_time = print_end - print_start;
  // Print out statistics
  printf("INFO. Total/expected sum: %d / %d (diff: %d)\n", total_sum, kNumDistances, total_sum - kNumDistances);
  printf("INFO. Algorithm time:  %f s\n", alg_tot_time);
  printf("INFO. Dist calc time:  %f s\n", (dist_calc_end - dist_calc_start) / tot_iter);
  printf("INFO. Update counts t: %f s\n", (update_counts_end - update_counts_start) / tot_iter);
  printf("INFO. File read time:  %f s\n", (file_rd_end - file_rd_start) / tot_iter);
  printf("INFO. Printing time:   %f s\n", print_tot_time);
  // Free all allocated memory
  free(line_buffer);
  free(base_cells);
  free(block_cells);
  free(block_distances);
  free(block_d_entries);
  fclose(fp);
  return 0;
}