#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <omp.h>
#include <time.h>
#include <stdint.h>
#include <fctnl.h>

typedef struct {
  signed char x_int, y_int, z_int;
  signed int16_t x_dec, y_dec, z_dec;
} cell_t;

void parse_number(const char* line_str, signed char* int_part,
    signed int16_t* dec_part) {
  // Example line: "-04.238 -07.514 +08.942"
  const char kASCII_Offset = 48;
  *int_part = (line_str[1] - kASCII_Offset) * 10 + line_str[2] - kASCII_Offset;
  *dec_part = (line_str[4] - kASCII_Offset) * 100 +
              (line_str[5] - kASCII_Offset) * 10 + line_str[6] - kASCII_Offset;
  if (line_str[0] == '-') {
    *int_part = - (*int_part);
    *dec_part = - (*dec_part);
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
  printf("%c%02d.%03d ", (cell.x_int > 0 ? '+' : '-'), cell.x_int, cell.x_dec);
  printf("%c%02d.%03d ", (cell.y_int > 0 ? '+' : '-'), cell.y_int, cell.y_dec);
  printf("%c%02d.%03d\n", (cell.z_int > 0 ? '+' : '-'), cell.z_int, cell.z_dec);
}

int main(int argc, char const *argv[]) {
  srand(time(NULL));
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

  char* line = (char*)malloc(kBytesPerLine);

  // Test reading cells
  for (int i = 0; i < kNumCells; ++i) {
    fread(line, kBytesPerLine, 1, fp);
    cell_t cell = parse_line(line);
    print_cell(cell);
  }

  fclose(fp);
  free(line);
  return 0;
}