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

#ifndef USE_FIX_POINT_CELLS
#define USE_FIX_POINT_CELLS
#endif

#ifndef FIX_POINT_FRAC_BITS
#define FIX_POINT_FRAC_BITS 10
#endif

const float fix2float_scaler = 1.0 / (float)(1 << FIX_POINT_FRAC_BITS);
const int kBytesPerLine = 23 + 1; // 23 character plus the return char
const int kBytesPerCell = kBytesPerLine / 3;
const char kASCII_Offset = 48;

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


inline int16_t char2fix(const char* str) {
  // Expected string example: "-04.238 ", 8 Bytes
  bool is_negative = str[0] == '+' - kASCII_Offset;
  // if (str[1] == '1') {
  //   return is_negative ? -(10 << FIX_POINT_FRAC_BITS) : (10 << FIX_POINT_FRAC_BITS);
  // }
  int16_t fix = ((str[1] * 10 + str[2]) << FIX_POINT_FRAC_BITS) + (str[4] >> 1) + (str[5] >> 2) + (str[6] >> 3);
  return is_negative ? -fix : fix;
}

inline void parse_number(const char* line_str, signed char* intp, int16_t* decp) {
  // Example line: "-04.238 -07.514 +08.942"
  *intp = line_str[2] + line_str[1] * 10;
  *decp = line_str[6] + line_str[5] * 10 + line_str[4] * 100;
  if (line_str[0] == '-' - kASCII_Offset) {
    *intp = - (*intp);
    *decp = - (*decp);
  }
}

inline void parse_lines(const int num_lines, char* line_buffer, cell_t* cells) {
  // Example line: "-04.238 -07.514 +08.942"
#pragma omp parallel for
  for (int i = 0; i < num_lines * kBytesPerLine; ++i) {
    line_buffer[i] -= kASCII_Offset;
  }
#pragma omp parallel for
  for (int i = 0; i < num_lines * kBytesPerLine; i += kBytesPerLine) {
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
  printf("%c%02d.%03d ", (cell.x_int > 0 ? '+' : '-'), abs(cell.x_int), abs(cell.x_dec));
  printf("%c%02d.%03d ", (cell.y_int > 0 ? '+' : '-'), abs(cell.y_int), abs(cell.y_dec));
  printf("%c%02d.%03d ", (cell.z_int > 0 ? '+' : '-'), abs(cell.z_int), abs(cell.z_dec));
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

#ifndef TASK_SIZE
#define TASK_SIZE 200
#endif

void mergeSortAux(int *X, int n, int *tmp) {
   int i = 0;
   int j = n/2;
   int ti = 0;

   while (i<n/2 && j<n) {
      if (X[i] < X[j]) {
         tmp[ti] = X[i];
         ti++; i++;
      } else {
         tmp[ti] = X[j];
         ti++; j++;
      }
   }
   while (i<n/2) { /* finish up lower half */
      tmp[ti] = X[i];
      ti++; i++;
   }
   while (j<n) { /* finish up upper half */
      tmp[ti] = X[j];
      ti++; j++;
   }
   memcpy(X, tmp, n*sizeof(int));
}

void mergeSort(int *X, int n, int *tmp)
{
   if (n < 2) return;

   #pragma omp task shared(X) if (n > TASK_SIZE)
   mergeSort(X, n/2, tmp);

   #pragma omp task shared(X) if (n > TASK_SIZE)
   mergeSort(X+(n/2), n-(n/2), tmp + n/2);

   #pragma omp taskwait
   mergeSortAux(X, n, tmp);
}

void init(int *a, int size){
   for(int i = 0; i < size; i++)
       a[i] = 0;
}

void printArray(int *a, int size){
   for(int i = 0; i < size; i++)
       printf("%d ", a[i]);
   printf("\n");
}

int isSorted(int *a, int size){
   for(int i = 0; i < size - 1; i++)
      if(a[i] > a[i + 1])
        return 0;
   return 1;
}

void free_table_entries(const int table_size, dist_t** hash_table) {
  for (int i = 0; i < table_size; ++i) {
    if (hash_table[i] != NULL) {
      free(hash_table[i]);
    }
  }
}


// void countSort(char arr[])
// {
//     // The output character array that will have sorted arr
//     char output[strlen(arr)];

//     // Create a count array to store count of individual
//     // characters and initialize count array as 0
//     int count[RANGE + 1], i;
//     memset(count, 0, sizeof(count));

//     // Store count of each character
//     for (i = 0; arr[i]; ++i)
//         ++count[arr[i]];

//     // Change count[i] so that count[i] now contains actual
//     // position of this character in output array
//     for (i = 1; i <= RANGE; ++i)
//         count[i] += count[i - 1];

//     // Build the output character array
//     for (i = 0; arr[i]; ++i) {
//         output[count[arr[i]] - 1] = arr[i];
//         --count[arr[i]];
//     }

//     /*
//      For Stable algorithm
//      for (i = sizeof(arr)-1; i>=0; --i)
//     {
//         output[count[arr[i]]-1] = arr[i];
//         --count[arr[i]];
//     }

//     For Logic : See implementation
//     */

//     // Copy the output array to arr, so that arr now
//     // contains sorted characters
//     for (i = 0; arr[i]; ++i)
//         arr[i] = output[i];
// }


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

  const int kBlockSize = 256;
  const int kBlockSizeX = 128;
  const int kBlockSizeY = 128;
  int m, p;

  // Cell block-buffers
  cell_t* base_cells = (cell_t*)malloc(sizeof(cell_t) * kBlockSizeX);
  cell_t* block_cells = (cell_t*)malloc(sizeof(cell_t) * kBlockSizeY);
  // Distance block-buffer
  uint16_t* block_d_entries = (uint16_t*)malloc(sizeof(uint16_t) * kBlockSizeX * kBlockSizeY);
  uint16_t** block_distances = (uint16_t**) malloc(sizeof(uint16_t*) * kBlockSizeX);
  for (int i = 0, j = 0; i < kBlockSizeX; ++i, j += kBlockSizeY) {
    block_distances[i] = block_d_entries + j;
  }
  // Hash table to store and count the final distances
  const int kHashTableSize = kNumCells * kNumCells / 2 - kNumCells / 2;
  dist_t** hash_table = (dist_t**) malloc(sizeof(dist_t*) * kHashTableSize);
  // Line buffer for reading "large chunks" from file
  char* line_buffer = (char*)malloc(kBytesPerLine * max(kBlockSizeX, kBlockSizeY));
  /*
   * Main algorithm
   */
  setvbuf(fp, NULL, _IONBF, BUFSIZ);

  int tot_iter = 0;
  double file_rd_start = 0;
  double file_rd_end = 0;
  double dist_calc_start = 0;
  double dist_calc_end = 0;
  double update_table_start = 0;
  double update_table_end = 0;
  double alg_start = omp_get_wtime();

// #pragma omp parallel
  for (int i = 0; i < kNumCells; i += kBlockSizeX) {
    /*
     * Read shared base cells (the columns):
     * 
     * 1. Seek to the right position in the file
     * 2. Get the block boundary (to avoid reading "outside" the file)
     * 3. Read the file into a line buffer
     * 4. Parse the lines into cell entries
     */
#pragma omp critical
    {
      m = min(i + kBlockSizeX, kNumCells) - i;
      fseek(fp, i * kBytesPerLine, SEEK_SET);
      fread(line_buffer, m * kBytesPerLine, 1, fp);
    }
    parse_lines(m, line_buffer, base_cells);
    /*
     * Scroll over all remaining cells starting from the current base (the rows)
     *
     * NOTE: Start from the next line/cell, since the distance of a cell with
     * itself is zero. To do so, j needs to be i+1.
     */
// #pragma omp parallel // shared(m, p) private(base_cells, block_cells, line_buffer)
    for (int j = i + 1; j < kNumCells; j += kBlockSizeY) {
      /*
       * Similarly as above, load all the file entries block-wise:
       * 
       * 1. Seek to the right position in the file
       * 2. Get the block boundary (to avoid reading "outside" the file)
       * 3. Read the file into a line buffer
       * 4. Parse the lines into cell entries
       */
#pragma omp critical
      {
        p = min(j + kBlockSizeY, kNumCells) - j;
        file_rd_start += omp_get_wtime();
        fseek(fp, j * kBytesPerLine, SEEK_SET);
        fread(line_buffer, p * kBytesPerLine, 1, fp);
        file_rd_end += omp_get_wtime();
      }
      parse_lines(p, line_buffer, block_cells);
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
      int ii, jj;
#pragma omp parallel for private(ii, jj)
      for (ii = 0; ii < m; ++ii) {
        for (jj = (j == i + 1) ? ii : 0; jj < p; ++jj) {
          block_distances[ii][jj] =
            dist2fix(cell_dist(base_cells[ii], block_cells[jj]));
        }
      }
      dist_calc_end += omp_get_wtime();

      /*
       * =======================================================================
       * TODO: IMPLEMENT A PARALLEL "COUNTING SORT" OVER THE BLOCK!
       * =======================================================================
       * Ref: https://www.journaldev.com/42355/counting-sort-algorithm
       */
#pragma omp parallel
      {
        int16_t* counts_local = (int16_t*)malloc(sizeof(int16_t) * kBlockSizeX * kBlockSizeY);
        memset(counts_local, 0, sizeof(int16_t) * kBlockSizeX * kBlockSizeY);
        for (ii = 0; ii < m; ++ii) {
          for (jj = (j == i + 1) ? ii : 0; jj < p; ++jj) {
            // update_dist(kHashTableSize, block_distances[ii][jj], hash_table);
          }
        }
//         int i, b_local[10] = {0};
// #pragma omp for nowait
//         for(i = 0; i < n; i++) {
//           b_local[a[i]]++;
//         }
//         #pragma omp critical
//         for(i=0; i<10; i++) {
//           b[i] += b_local[i];
//         }
        free(counts_local);
      }



      /*
       * Update the hash table containing all distances found so far
       */
      update_table_start += omp_get_wtime();
#pragma omp critical
// #pragma omp parallel for private(ii, jj) shared(hash_table)
      for (ii = 0; ii < m; ++ii) {
        for (jj = (j == i + 1) ? ii : 0; jj < p; ++jj) {
          // update_dist(kHashTableSize, block_distances[ii][jj], hash_table);

          int hash_index;
          dist_t* found_dist = search_table(kHashTableSize, block_distances[ii][jj], hash_table, &hash_index);
          if (found_dist) {
            found_dist->count++;
          } else {
            hash_table[hash_index] = (dist_t*)malloc(sizeof(dist_t));
            hash_table[hash_index]->val = block_distances[ii][jj];
            hash_table[hash_index]->count = 1;
          }

        }
      }
      update_table_end += omp_get_wtime();

      ++tot_iter;
    } // end row-scrolling
  }
  double alg_end = omp_get_wtime();
  double alg_tot_time = alg_end - alg_start;

  // Sort hash table and finally print
  double sort_start = omp_get_wtime();
// #pragma omp parallel shared(hash_table)
  {
    qsort(hash_table, kHashTableSize, sizeof(dist_t*), compare_dist);
  }

  double sort_end = omp_get_wtime();
  double sort_tot_time = sort_end - sort_start;

  for (int i = 0; i < kHashTableSize; ++i) {
    if (hash_table[i] != NULL) {
      print_dist(*hash_table[i]);
    }
  }
  // Print out statistics
  printf("INFO. Algorithm time: %f s\n", alg_tot_time);
  printf("INFO. Dist calc time: %f s\n", (dist_calc_end - dist_calc_start) / tot_iter);
  printf("INFO. Update table t: %f s\n", (update_table_end - update_table_start) / tot_iter);
  printf("INFO. File read time: %f s\n", (file_rd_end - file_rd_start) / tot_iter);
  printf("INFO. Sorting time:   %f s\n", sort_tot_time);
  // Free all allocated memory
  free(line_buffer);
  free(base_cells);
  free(block_cells);
  free(block_distances);
  free(block_d_entries);
  // NOTE: Special care is needed to free the elements in the hash table
  free_table_entries(kHashTableSize, hash_table);
  free(hash_table);

  fclose(fp);

  int16_t x = (int16_t)(float)(9.0  * 0.5 + 9.0  * 0.25);
  int16_t max_fract = 6;
  float y = (float)max_fract / (1 << 3);

  printf("%x\n", x);
  printf("%f\n", y);

  return 0;
}