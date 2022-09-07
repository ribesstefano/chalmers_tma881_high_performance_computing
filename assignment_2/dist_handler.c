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