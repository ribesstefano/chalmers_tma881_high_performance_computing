__kernel void dot_prod_mul(__global const double *a,
                           __global const double *b,
                           __global double *c) {
  int ix = get_global_id(0);
  c[ix] = a[ix] * b[ix];
}

__kernel void mat_mul(__global const double *a,
                      __global const double *b,
                      __global double *c,
                      int width_a,
                      int width_b) {
  int ix = get_global_id(0);
  int iy = get_global_id(1);
  double value = 0;
  for (int kx = 0; kx < width_a; ++kx) {
    double a_elt = a[iy * width_a + kx];
    double b_elt = b[kx * width_b + ix];
    value += a_elt * b_elt;
  }
  c[iy * width_b + ix] = value;

  // Have to think about things like coalescence & bank conflicts

  // __local xL[256]; __local yL[256]; //declare local mem array with constant length
  // barrier(CLK_LOCAL_MEM_FENCE);  //synchronization
}

__kernel void heat_diffusion(__global double* box_glob,
                             __const int num_iter,
                             __const int width,
                             __const int height,
                             __const double diffusion_const) {
  // NOTEs and TODOs:
  // 
  // * Difference b/w "const" and "__const"?
  // * Have to think about things like coalescence & bank conflicts
  // * __local xL[256]; //declare local mem array with constant length
  const int x = get_global_id(0);
  const int y = get_global_id(1);
  for (int i = 0; i < num_iter; ++i) {
    const double n = (x > 1) ? box_glob[y * height + x - 1] : 0.0;
    const double w = (y > 1) ? box_glob[(y - 1) * height + x] : 0.0;
    const double s = (x < height - 1) ? box_glob[y * height + x + 1] : 0.0;
    const double e = (y < width - 1) ? box_glob[(y + 1) * height + x] : 0.0;
    double box_val = box_glob[y * height + x];
    box_val += diffusion_const * ((n + s + w + e) / 4 - box_val);
    barrier(CLK_GLOBAL_MEM_FENCE);
    box_glob[y * height + x] = box_val;
    barrier(CLK_GLOBAL_MEM_FENCE);
  }
}

double absval(const double x) {
  return (x > 0.0) ? x : -x;
}

__kernel void sum_reduction(__global double *buffer,
                            __local double *scratch,
                            __const int sz,
                            const char abs_difference,
                            const double avg_in,
                            __global double *result) {
  int gsz = get_global_size(0);
  int gix = get_global_id(0);
  int lsz = get_local_size(0);
  int lix = get_local_id(0);
  double acc = 0;
  for (int i = get_global_id(0); i < sz; i += gsz) {
    acc += (abs_difference) ? absval(buffer[i] - avg_in) : buffer[i];
  }
  scratch[lix] = acc;
  barrier(CLK_LOCAL_MEM_FENCE);
  for (int offset = lsz / 2; offset > 0; offset /= 2) {
    if (lix < offset) {
      scratch[lix] += scratch[lix + offset];
    }
    barrier(CLK_LOCAL_MEM_FENCE);
  }
  if (lix == 0) {
    result[get_group_id(0)] = scratch[0];
  }
}