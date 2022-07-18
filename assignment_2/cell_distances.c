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

inline void offset_line_buffer(const int len, char* line_buffer) {
  const char kASCII_Offset = 48;
#pragma omp parallel for
  for (int i = 0; i < len; ++i) {
    line_buffer[i] -= kASCII_Offset;
  }
}

inline void parse_number(const char* line_str, signed char* intp, int16_t* decp) {
  // Example line: "-04.238 -07.514 +08.942"
  const char kASCII_Offset = 48;
  *intp = line_str[2] + line_str[1] * 10;
  *decp = line_str[6] + line_str[5] * 10 + line_str[4] * 100;
  if (line_str[0] == '-' - kASCII_Offset) {
    *intp = - (*intp);
    *decp = - (*decp);
  }
}

inline cell_t parse_line(const char* line_str) {
  // Example line: "-04.238 -07.514 +08.942"
  cell_t encoded;
  parse_number(&line_str[0], &encoded.x_int, &encoded.x_dec);
  parse_number(&line_str[8], &encoded.y_int, &encoded.y_dec);
  parse_number(&line_str[16], &encoded.z_int, &encoded.z_dec);
  return encoded;
}

inline cell_t read_cell(FILE* fp) {
  const int kBytesPerLine = 23 + 1; // 23 character plus the return char
  char line_buffer[kBytesPerLine];
  if (fread(line_buffer, kBytesPerLine, 1, fp) != kBytesPerLine) {
    // fprintf(stderr, "ERROR. Wrong number of bytes read.\n");
    // exit(EXIT_FAILURE);
  }
  return parse_line(line_buffer);
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

inline float cell_dist(const cell_t a, const cell_t b) {
  float x_diff = cell2float(a, 'x') - cell2float(b, 'x');
  float y_diff = cell2float(a, 'y') - cell2float(b, 'y');
  float z_diff = cell2float(a, 'z') - cell2float(b, 'z');
  float x_pow = x_diff * x_diff;
  float y_pow = y_diff * y_diff;
  float z_pow = z_diff * z_diff;
  return sqrt(x_pow + y_pow + z_pow);
}

inline uint16_t dist2fix(const float d) {
  return (uint16_t)(d * (1 << FIX_POINT_FRAC_BITS));
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

int compare_dist(const void* a, const void* b) {
  // Used in qsort: compare the fixed point value of two distance struct.
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

uint32_t get_hash(const uint16_t din) {
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

dist_t* search_table(const int table_size, const uint16_t dist, dist_t** hash_table, int* hash_index) {
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

void delete_table(const int table_size, dist_t** hash_table) {
  for (int i = 0; i < table_size; ++i) {
    if (hash_table[i] != NULL) {
      free(hash_table[i]);
    }
  }
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
  const int kBytesPerLine = 23 + 1; // 23 character plus the return char
  const long int kNumCells = ftell(fp) / kBytesPerLine;

  const int kBlockSize = 512;
  int m, p;

  // Cell block-buffers
  cell_t* base_cells = (cell_t*)malloc(sizeof(cell_t) * kBlockSize);
  cell_t* block_cells = (cell_t*)malloc(sizeof(cell_t) * kBlockSize);
  // Distance block-buffer
  uint16_t* block_d_entries = (uint16_t*)malloc(sizeof(uint16_t) * kBlockSize * kBlockSize);
  uint16_t** block_distances = (uint16_t**) malloc(sizeof(uint16_t*) * kBlockSize);
  for (int i = 0, j = 0; i < kBlockSize; ++i, j += kBlockSize) {
    block_distances[i] = block_d_entries + j;
  }
  // Hash table to store and count the final distances
  const int kHashTableSize = kNumCells * kNumCells / 2 - kNumCells / 2;
  dist_t** hash_table = (dist_t**) malloc(sizeof(dist_t*) * kHashTableSize);
  // Line buffer for reading "large chunks" from file
  char* line_buffer = (char*)malloc(kBytesPerLine * kBlockSize);
  /*
   * Main algorithm
   */
  clock_t alg_start = clock();

#pragma omp parallel
  for (int i = 0; i < kNumCells; i += kBlockSize) {
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
      fseek(fp, i * kBytesPerLine, SEEK_SET);
      m = min(i + kBlockSize, kNumCells) - i;
      fread(line_buffer, m * kBytesPerLine, 1, fp);
    }
    offset_line_buffer(m * kBytesPerLine, line_buffer);
#pragma omp parallel for shared(m, base_cells, line_buffer)
    for (int j = 0; j < m; ++j) {
      base_cells[j] = parse_line(&line_buffer[j * kBytesPerLine]);
    }
    /*
     * Scroll over all remaining cells starting from the current base (the rows)
     */
// #pragma omp parallel // shared(m, p) private(base_cells, block_cells, line_buffer)
    for (int j = i + 1; j < kNumCells; j += kBlockSize) {
      /*
       * Similarly as above, load all the file entries block-wise:
       * 
       * 1. Seek to the right position in the file
       * 2. Get the block boundary (to avoid reading "outside" the file)
       * 3. Read the file into a line buffer
       * 4. Parse the lines into cell entries
       */
#pragma omp single
      {
        fseek(fp, j * kBytesPerLine, SEEK_SET);
        p = min(j + kBlockSize, kNumCells) - j;
        fread(line_buffer, p * kBytesPerLine, 1, fp);
      }
      offset_line_buffer(p * kBytesPerLine, line_buffer);
#pragma omp parallel for shared(p, block_cells, line_buffer)
      for (int k = 0; k < p; ++k) {
        block_cells[k] = parse_line(&line_buffer[k * kBytesPerLine]);
      }
      /*
       * Block computation: each element can be computed in parallel
       *
       * NOTE: If it is the left-most block, then we need to skip the first
       * elements on the left.
       *
       * NOTE: OpenMP collapse clause cannot be used because the inner loop has
       * variable loop iterations!
       */
#pragma omp parallel for shared(m, p, base_cells, block_cells)
      for (int ii = 0; ii < m; ++ii) {
        for (int jj = (j == i + 1) ? ii : 0; jj < p; ++jj) {
          block_distances[ii][jj] =
            dist2fix(cell_dist(base_cells[ii], block_cells[jj]));
        }
      }
      /*
       * Update the hash table containing all distances found so far
       */
#pragma omp parallel for shared(m, p, block_distances, hash_table)
      for (int ii = 0; ii < m; ++ii) {
        for (int jj = (j == i + 1) ? ii : 0; jj < p; ++jj) {
          update_dist(kHashTableSize, block_distances[ii][jj], hash_table);
        }
      }
    } // end row-scrolling
  }
  clock_t alg_end = clock();
  float alg_tot_time = (float)(alg_end - alg_start) / CLOCKS_PER_SEC;

  // Sort hash table and finally print
  clock_t sort_start = clock();
  qsort(hash_table, kHashTableSize, sizeof(dist_t*), compare_dist);
  clock_t sort_end = clock();
  float sort_tot_time = (float)(sort_end - sort_start) / CLOCKS_PER_SEC;

  for (int i = 0; i < kHashTableSize; ++i) {
    if (hash_table[i] != NULL) {
      print_dist(*hash_table[i]);
    }
  }

  printf("INFO. Algorithm time: %f s\n", alg_tot_time);
  printf("INFO. Sorting time:   %f s\n", sort_tot_time);

  free(line_buffer);
  free(base_cells);
  free(block_cells);
  free(block_distances);
  // NOTE: Special care is needed to free the elements in the hash table
  delete_table(kHashTableSize, hash_table);
  free(hash_table);

  fclose(fp);
  return 0;
}