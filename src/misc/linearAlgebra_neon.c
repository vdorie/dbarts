#include "config.h"
#include <misc/linearAlgebra.h>

#include <stdint.h> // uintptr_t

#include <misc/intrinsic.h>

#include <external/io.h>

void misc_addVectors_neon(const double* restrict x, size_t length, const double* restrict y, double* restrict z)
{
  if (length == 0) return;
  
  size_t y_offset = ((uintptr_t) y) % (8 * sizeof(double));
  size_t x_offset = ((uintptr_t) x) % (8 * sizeof(double));
  size_t z_offset = ((uintptr_t) z) % (8 * sizeof(double));
  size_t prefix = z_offset == 0 ? 0 : (8 * sizeof(double) - z_offset) / sizeof(double);

  if (prefix > length) prefix = length;
  
  size_t i = 0;
  size_t suffix = prefix + 16 * ((length - prefix) / 16);
  
  for ( ; i < prefix; ++i)
    z[i] = y[i] + x[i];
  
  if (z_offset == x_offset && z_offset == y_offset) {
    for ( ; i < suffix; i += 16) {
      float64x2x4_t x_vec = vld1q_f64_x4(x + i);
      float64x2x4_t y_vec = vld1q_f64_x4(y + i);
      y_vec.val[0] = vaddq_f64(y_vec.val[0], x_vec.val[0]);
      y_vec.val[1] = vaddq_f64(y_vec.val[1], x_vec.val[1]);
      y_vec.val[2] = vaddq_f64(y_vec.val[2], x_vec.val[2]);
      y_vec.val[3] = vaddq_f64(y_vec.val[3], x_vec.val[3]);
      vst1q_f64_x4(z + i, y_vec);
      
      x_vec = vld1q_f64_x4(x + i + 8);
      y_vec = vld1q_f64_x4(y + i + 8);
      y_vec.val[0] = vaddq_f64(y_vec.val[0], x_vec.val[0]);
      y_vec.val[1] = vaddq_f64(y_vec.val[1], x_vec.val[1]);
      y_vec.val[2] = vaddq_f64(y_vec.val[2], x_vec.val[2]);
      y_vec.val[3] = vaddq_f64(y_vec.val[3], x_vec.val[3]);
      vst1q_f64_x4(z + i + 8, y_vec);
    }
    
    suffix = prefix + 4 * ((length - prefix) / 4);
    
    for ( ; i < suffix; i += 4) {
      float64x2x2_t x_vec = vld1q_f64_x2(x + i);
      float64x2x2_t y_vec = vld1q_f64_x2(y + i);
      y_vec.val[0] = vaddq_f64(y_vec.val[0], x_vec.val[0]);
      y_vec.val[1] = vaddq_f64(y_vec.val[1], x_vec.val[1]);
      vst1q_f64_x2(z + i, y_vec);
    }
    
  } else {
    for ( ; i < suffix; i += 16) {
      vst1q_f64(z + i     , vaddq_f64(vld1q_f64(y + i     ), vld1q_f64(x + i     )));
      vst1q_f64(z + i +  2, vaddq_f64(vld1q_f64(y + i +  2), vld1q_f64(x + i +  2)));
      vst1q_f64(z + i +  4, vaddq_f64(vld1q_f64(y + i +  4), vld1q_f64(x + i +  4)));
      vst1q_f64(z + i +  6, vaddq_f64(vld1q_f64(y + i +  6), vld1q_f64(x + i +  6)));
      
      vst1q_f64(z + i +  8, vaddq_f64(vld1q_f64(y + i +  8), vld1q_f64(x + i +  8)));
      vst1q_f64(z + i + 10, vaddq_f64(vld1q_f64(y + i + 10), vld1q_f64(x + i + 10)));
      vst1q_f64(z + i + 12, vaddq_f64(vld1q_f64(y + i + 12), vld1q_f64(x + i + 12)));
      vst1q_f64(z + i + 14, vaddq_f64(vld1q_f64(y + i + 14), vld1q_f64(x + i + 14)));
    }
    
    suffix = prefix + 4 * ((length - prefix) / 4);
    
    for ( ; i < suffix; i += 4) {
      vst1q_f64(z + i     , vaddq_f64(vld1q_f64(y + i     ), vld1q_f64(x + i     )));
      vst1q_f64(z + i +  2, vaddq_f64(vld1q_f64(y + i +  2), vld1q_f64(x + i +  2)));
    }
  }
  
  for ( ; i < length; ++i)
    z[i] = y[i] + x[i];
}

void misc_subtractVectors_neon(const double* restrict x, size_t length, const double* restrict y, double* restrict z)
{
  if (length == 0) return;
  
  size_t y_offset = ((uintptr_t) y) % (8 * sizeof(double));
  size_t x_offset = ((uintptr_t) x) % (8 * sizeof(double));
  size_t z_offset = ((uintptr_t) z) % (8 * sizeof(double));
  size_t prefix = z_offset == 0 ? 0 : (8 * sizeof(double) - z_offset) / sizeof(double);

  if (prefix > length) prefix = length;
  
  size_t i = 0;
  size_t suffix = prefix + 16 * ((length - prefix) / 16);
  
  for ( ; i < prefix; ++i)
    z[i] = y[i] - x[i];
  
  if (z_offset == x_offset && z_offset == y_offset) {
    for ( ; i < suffix; i += 16) {
      float64x2x4_t x_vec = vld1q_f64_x4(x + i);
      float64x2x4_t y_vec = vld1q_f64_x4(y + i);
      y_vec.val[0] = vsubq_f64(y_vec.val[0], x_vec.val[0]);
      y_vec.val[1] = vsubq_f64(y_vec.val[1], x_vec.val[1]);
      y_vec.val[2] = vsubq_f64(y_vec.val[2], x_vec.val[2]);
      y_vec.val[3] = vsubq_f64(y_vec.val[3], x_vec.val[3]);
      vst1q_f64_x4(z + i, y_vec);
      
      x_vec = vld1q_f64_x4(x + i + 8);
      y_vec = vld1q_f64_x4(y + i + 8);
      y_vec.val[0] = vsubq_f64(y_vec.val[0], x_vec.val[0]);
      y_vec.val[1] = vsubq_f64(y_vec.val[1], x_vec.val[1]);
      y_vec.val[2] = vsubq_f64(y_vec.val[2], x_vec.val[2]);
      y_vec.val[3] = vsubq_f64(y_vec.val[3], x_vec.val[3]);
      vst1q_f64_x4(z + i + 8, y_vec);
    }
    
    suffix = prefix + 4 * ((length - prefix) / 4);
    
    for ( ; i < suffix; i += 4) {
      float64x2x2_t x_vec = vld1q_f64_x2(x + i);
      float64x2x2_t y_vec = vld1q_f64_x2(y + i);
      y_vec.val[0] = vsubq_f64(y_vec.val[0], x_vec.val[0]);
      y_vec.val[1] = vsubq_f64(y_vec.val[1], x_vec.val[1]);
      vst1q_f64_x2(z + i, y_vec);
    }
    
  } else {
    for ( ; i < suffix; i += 16) {
      vst1q_f64(z + i     , vsubq_f64(vld1q_f64(y + i     ), vld1q_f64(x + i     )));
      vst1q_f64(z + i +  2, vsubq_f64(vld1q_f64(y + i +  2), vld1q_f64(x + i +  2)));
      vst1q_f64(z + i +  4, vsubq_f64(vld1q_f64(y + i +  4), vld1q_f64(x + i +  4)));
      vst1q_f64(z + i +  6, vsubq_f64(vld1q_f64(y + i +  6), vld1q_f64(x + i +  6)));
      
      vst1q_f64(z + i +  8, vsubq_f64(vld1q_f64(y + i +  8), vld1q_f64(x + i +  8)));
      vst1q_f64(z + i + 10, vsubq_f64(vld1q_f64(y + i + 10), vld1q_f64(x + i + 10)));
      vst1q_f64(z + i + 12, vsubq_f64(vld1q_f64(y + i + 12), vld1q_f64(x + i + 12)));
      vst1q_f64(z + i + 14, vsubq_f64(vld1q_f64(y + i + 14), vld1q_f64(x + i + 14)));
    }
    
    suffix = prefix + 4 * ((length - prefix) / 4);
    
    for ( ; i < suffix; i += 4) {
      vst1q_f64(z + i     , vsubq_f64(vld1q_f64(y + i     ), vld1q_f64(x + i     )));
      vst1q_f64(z + i +  2, vsubq_f64(vld1q_f64(y + i +  2), vld1q_f64(x + i +  2)));
    }
  }
  
  for ( ; i < length; ++i)
    z[i] = y[i] - x[i];
}

void misc_addVectorsWithMultiplier_neon(const double* restrict x, size_t length, double alpha, const double* restrict y, double* restrict z)
{
  if (length == 0) return;
  
  size_t y_offset = ((uintptr_t) y) % (8 * sizeof(double));
  size_t x_offset = ((uintptr_t) x) % (8 * sizeof(double));
  size_t z_offset = ((uintptr_t) z) % (8 * sizeof(double));
  size_t prefix = z_offset == 0 ? 0 : (8 * sizeof(double) - z_offset) / sizeof(double);

  if (prefix > length) prefix = length;
  
  size_t i = 0;
  size_t suffix = prefix + 16 * ((length - prefix) / 16);

  for ( ; i < prefix; ++i)
    z[i] = y[i] + alpha * x[i];

  float64x2_t alpha_vec = vdupq_n_f64(alpha);
  if (z_offset == x_offset && z_offset == y_offset) {
    for ( ; i < suffix; i += 16) {
      float64x2x4_t x_vec = vld1q_f64_x4(x + i);
      float64x2x4_t y_vec = vld1q_f64_x4(y + i);
      y_vec.val[0] = vaddq_f64(y_vec.val[0], vmulq_f64(x_vec.val[0], alpha_vec));
      y_vec.val[1] = vaddq_f64(y_vec.val[1], vmulq_f64(x_vec.val[1], alpha_vec));
      y_vec.val[2] = vaddq_f64(y_vec.val[2], vmulq_f64(x_vec.val[2], alpha_vec));
      y_vec.val[3] = vaddq_f64(y_vec.val[3], vmulq_f64(x_vec.val[3], alpha_vec));
      vst1q_f64_x4(z + i, y_vec);
      
      x_vec = vld1q_f64_x4(x + i + 8);
      y_vec = vld1q_f64_x4(y + i + 8);
      y_vec.val[0] = vaddq_f64(y_vec.val[0], vmulq_f64(x_vec.val[0], alpha_vec));
      y_vec.val[1] = vaddq_f64(y_vec.val[1], vmulq_f64(x_vec.val[1], alpha_vec));
      y_vec.val[2] = vaddq_f64(y_vec.val[2], vmulq_f64(x_vec.val[2], alpha_vec));
      y_vec.val[3] = vaddq_f64(y_vec.val[3], vmulq_f64(x_vec.val[3], alpha_vec));
      vst1q_f64_x4(z + i + 8, y_vec);
    }
    
    suffix = prefix + 4 * ((length - prefix) / 4);
    
    for ( ; i < suffix; i += 4) {
      float64x2x2_t x_vec = vld1q_f64_x2(x + i);
      float64x2x2_t y_vec = vld1q_f64_x2(y + i);
      y_vec.val[0] = vaddq_f64(y_vec.val[0], vmulq_f64(x_vec.val[0], alpha_vec));
      y_vec.val[1] = vaddq_f64(y_vec.val[1], vmulq_f64(x_vec.val[1], alpha_vec));
      vst1q_f64_x2(z + i, y_vec);
    }
  } else {
    for ( ; i < suffix; i += 16) {
      vst1q_f64(z + i     , vaddq_f64(vld1q_f64(y + i     ), vmulq_f64(vld1q_f64(x + i     ), alpha_vec)));
      vst1q_f64(z + i +  2, vaddq_f64(vld1q_f64(y + i +  2), vmulq_f64(vld1q_f64(x + i +  2), alpha_vec)));
      vst1q_f64(z + i +  4, vaddq_f64(vld1q_f64(y + i +  4), vmulq_f64(vld1q_f64(x + i +  4), alpha_vec)));
      vst1q_f64(z + i +  6, vaddq_f64(vld1q_f64(y + i +  6), vmulq_f64(vld1q_f64(x + i +  6), alpha_vec)));
      
      vst1q_f64(z + i +  8, vaddq_f64(vld1q_f64(y + i +  8), vmulq_f64(vld1q_f64(x + i +  8), alpha_vec)));
      vst1q_f64(z + i + 10, vaddq_f64(vld1q_f64(y + i + 10), vmulq_f64(vld1q_f64(x + i + 10), alpha_vec)));
      vst1q_f64(z + i + 12, vaddq_f64(vld1q_f64(y + i + 12), vmulq_f64(vld1q_f64(x + i + 12), alpha_vec)));
      vst1q_f64(z + i + 14, vaddq_f64(vld1q_f64(y + i + 14), vmulq_f64(vld1q_f64(x + i + 14), alpha_vec)));
    }
    
    suffix = prefix + 4 * ((length - prefix) / 4);
    
    for ( ; i < suffix; i += 4) {
      vst1q_f64(z + i    , vaddq_f64(vld1q_f64(y + i    ), vmulq_f64(vld1q_f64(x + i    ), alpha_vec)));
      vst1q_f64(z + i + 2, vaddq_f64(vld1q_f64(y + i + 2), vmulq_f64(vld1q_f64(x + i + 2), alpha_vec)));
    }
  }
  
  for ( ; i < length; ++i)
    z[i] = y[i] + alpha * x[i];
}

void misc_addVectorsInPlace_neon(const double* restrict x, size_t length, double* restrict y)
{
  if (length == 0) return;
  
  size_t y_offset = ((uintptr_t) y) % (8 * sizeof(double));
  size_t x_offset = ((uintptr_t) x) % (8 * sizeof(double));
  size_t prefix = y_offset == 0 ? 0 : (8 * sizeof(double) - y_offset) / sizeof(double);

  if (prefix > length) prefix = length;
  
  size_t i = 0;
  size_t suffix = prefix + 16 * ((length - prefix) / 16);
  
  for ( ; i < prefix; ++i)
    y[i] += x[i];
  
  if (y_offset == x_offset) {
    for ( ; i < suffix; i += 16) {
      float64x2x4_t x_vec = vld1q_f64_x4(x + i);
      float64x2x4_t y_vec = vld1q_f64_x4(y + i);
      y_vec.val[0] = vaddq_f64(y_vec.val[0], x_vec.val[0]);
      y_vec.val[1] = vaddq_f64(y_vec.val[1], x_vec.val[1]);
      y_vec.val[2] = vaddq_f64(y_vec.val[2], x_vec.val[2]);
      y_vec.val[3] = vaddq_f64(y_vec.val[3], x_vec.val[3]);
      vst1q_f64_x4(y + i, y_vec);
      
      x_vec = vld1q_f64_x4(x + i + 8);
      y_vec = vld1q_f64_x4(y + i + 8);
      y_vec.val[0] = vaddq_f64(y_vec.val[0], x_vec.val[0]);
      y_vec.val[1] = vaddq_f64(y_vec.val[1], x_vec.val[1]);
      y_vec.val[2] = vaddq_f64(y_vec.val[2], x_vec.val[2]);
      y_vec.val[3] = vaddq_f64(y_vec.val[3], x_vec.val[3]);
      vst1q_f64_x4(y + i + 8, y_vec);
    }
    
    suffix = prefix + 4 * ((length - prefix) / 4);
    
    for ( ; i < suffix; i += 4) {
      float64x2x2_t x_vec = vld1q_f64_x2(x + i);
      float64x2x2_t y_vec = vld1q_f64_x2(y + i);
      y_vec.val[0] = vaddq_f64(y_vec.val[0], x_vec.val[0]);
      y_vec.val[1] = vaddq_f64(y_vec.val[1], x_vec.val[1]);
      vst1q_f64_x2(y + i, y_vec);
    }
    
  } else {
    for ( ; i < suffix; i += 16) {
      vst1q_f64(y + i     , vaddq_f64(vld1q_f64(y + i     ), vld1q_f64(x + i     )));
      vst1q_f64(y + i +  2, vaddq_f64(vld1q_f64(y + i +  2), vld1q_f64(x + i +  2)));
      vst1q_f64(y + i +  4, vaddq_f64(vld1q_f64(y + i +  4), vld1q_f64(x + i +  4)));
      vst1q_f64(y + i +  6, vaddq_f64(vld1q_f64(y + i +  6), vld1q_f64(x + i +  6)));
      
      vst1q_f64(y + i +  8, vaddq_f64(vld1q_f64(y + i +  8), vld1q_f64(x + i +  8)));
      vst1q_f64(y + i + 10, vaddq_f64(vld1q_f64(y + i + 10), vld1q_f64(x + i + 10)));
      vst1q_f64(y + i + 12, vaddq_f64(vld1q_f64(y + i + 12), vld1q_f64(x + i + 12)));
      vst1q_f64(y + i + 14, vaddq_f64(vld1q_f64(y + i + 14), vld1q_f64(x + i + 14)));
    }
    
    suffix = prefix + 4 * ((length - prefix) / 4);
    
    for ( ; i < suffix; i += 4) {
      vst1q_f64(y + i     , vaddq_f64(vld1q_f64(y + i     ), vld1q_f64(x + i     )));
      vst1q_f64(y + i +  2, vaddq_f64(vld1q_f64(y + i +  2), vld1q_f64(x + i +  2)));
    }
  }
  
  for ( ; i < length; ++i)
    y[i] += x[i];
}

void misc_subtractVectorsInPlace_neon(const double* restrict x, size_t length, double* restrict y)
{
  if (length == 0) return;
  
  size_t y_offset = ((uintptr_t) y) % (8 * sizeof(double));
  size_t x_offset = ((uintptr_t) x) % (8 * sizeof(double));
  size_t prefix = y_offset == 0 ? 0 : (8 * sizeof(double) - y_offset) / sizeof(double);

  if (prefix > length) prefix = length;
  
  size_t i = 0;
  size_t suffix = prefix + 16 * ((length - prefix) / 16);
  
  for ( ; i < prefix; ++i)
    y[i] -= x[i];
  
  if (y_offset == x_offset) {
    for ( ; i < suffix; i += 16) {
      float64x2x4_t x_vec = vld1q_f64_x4(x + i);
      float64x2x4_t y_vec = vld1q_f64_x4(y + i);
      y_vec.val[0] = vsubq_f64(y_vec.val[0], x_vec.val[0]);
      y_vec.val[1] = vsubq_f64(y_vec.val[1], x_vec.val[1]);
      y_vec.val[2] = vsubq_f64(y_vec.val[2], x_vec.val[2]);
      y_vec.val[3] = vsubq_f64(y_vec.val[3], x_vec.val[3]);
      vst1q_f64_x4(y + i, y_vec);
      
      x_vec = vld1q_f64_x4(x + i + 8);
      y_vec = vld1q_f64_x4(y + i + 8);
      y_vec.val[0] = vsubq_f64(y_vec.val[0], x_vec.val[0]);
      y_vec.val[1] = vsubq_f64(y_vec.val[1], x_vec.val[1]);
      y_vec.val[2] = vsubq_f64(y_vec.val[2], x_vec.val[2]);
      y_vec.val[3] = vsubq_f64(y_vec.val[3], x_vec.val[3]);
      vst1q_f64_x4(y + i + 8, y_vec);
    }
    
    suffix = prefix + 4 * ((length - prefix) / 4);
    
    for ( ; i < suffix; i += 4) {
      float64x2x2_t x_vec = vld1q_f64_x2(x + i);
      float64x2x2_t y_vec = vld1q_f64_x2(y + i);
      y_vec.val[0] = vsubq_f64(y_vec.val[0], x_vec.val[0]);
      y_vec.val[1] = vsubq_f64(y_vec.val[1], x_vec.val[1]);
      vst1q_f64_x2(y + i, y_vec);
    }
    
  } else {
    for ( ; i < suffix; i += 16) {
      vst1q_f64(y + i     , vsubq_f64(vld1q_f64(y + i     ), vld1q_f64(x + i     )));
      vst1q_f64(y + i +  2, vsubq_f64(vld1q_f64(y + i +  2), vld1q_f64(x + i +  2)));
      vst1q_f64(y + i +  4, vsubq_f64(vld1q_f64(y + i +  4), vld1q_f64(x + i +  4)));
      vst1q_f64(y + i +  6, vsubq_f64(vld1q_f64(y + i +  6), vld1q_f64(x + i +  6)));
      
      vst1q_f64(y + i +  8, vsubq_f64(vld1q_f64(y + i +  8), vld1q_f64(x + i +  8)));
      vst1q_f64(y + i + 10, vsubq_f64(vld1q_f64(y + i + 10), vld1q_f64(x + i + 10)));
      vst1q_f64(y + i + 12, vsubq_f64(vld1q_f64(y + i + 12), vld1q_f64(x + i + 12)));
      vst1q_f64(y + i + 14, vsubq_f64(vld1q_f64(y + i + 14), vld1q_f64(x + i + 14)));
    }
    
    suffix = prefix + 4 * ((length - prefix) / 4);
    
    for ( ; i < suffix; i += 4) {
      vst1q_f64(y + i     , vsubq_f64(vld1q_f64(y + i     ), vld1q_f64(x + i     )));
      vst1q_f64(y + i +  2, vsubq_f64(vld1q_f64(y + i +  2), vld1q_f64(x + i +  2)));
    }
  }
  
  for ( ; i < length; ++i)
    y[i] -= x[i];
}

void misc_addVectorsInPlaceWithMultiplier_neon(const double* restrict x, size_t length, double alpha, double* restrict y)
{
  if (length == 0) return;
  
  size_t y_offset = ((uintptr_t) y) % (8 * sizeof(double));
  size_t x_offset = ((uintptr_t) x) % (8 * sizeof(double));
  size_t prefix = y_offset == 0 ? 0 : (8 * sizeof(double) - y_offset) / sizeof(double);

  if (prefix > length) prefix = length;
  
  size_t i = 0;
  size_t suffix = prefix + 16 * ((length - prefix) / 16);

  for ( ; i < prefix; ++i)
    y[i] += alpha * x[i];

  float64x2_t alpha_vec = vdupq_n_f64(alpha);
  if (y_offset == x_offset) {
    for ( ; i < suffix; i += 16) {
      float64x2x4_t x_vec = vld1q_f64_x4(x + i);
      float64x2x4_t y_vec = vld1q_f64_x4(y + i);
      y_vec.val[0] = vaddq_f64(y_vec.val[0], vmulq_f64(x_vec.val[0], alpha_vec));
      y_vec.val[1] = vaddq_f64(y_vec.val[1], vmulq_f64(x_vec.val[1], alpha_vec));
      y_vec.val[2] = vaddq_f64(y_vec.val[2], vmulq_f64(x_vec.val[2], alpha_vec));
      y_vec.val[3] = vaddq_f64(y_vec.val[3], vmulq_f64(x_vec.val[3], alpha_vec));
      vst1q_f64_x4(y + i, y_vec);
      
      x_vec = vld1q_f64_x4(x + i + 8);
      y_vec = vld1q_f64_x4(y + i + 8);
      y_vec.val[0] = vaddq_f64(y_vec.val[0], vmulq_f64(x_vec.val[0], alpha_vec));
      y_vec.val[1] = vaddq_f64(y_vec.val[1], vmulq_f64(x_vec.val[1], alpha_vec));
      y_vec.val[2] = vaddq_f64(y_vec.val[2], vmulq_f64(x_vec.val[2], alpha_vec));
      y_vec.val[3] = vaddq_f64(y_vec.val[3], vmulq_f64(x_vec.val[3], alpha_vec));
      vst1q_f64_x4(y + i + 8, y_vec);
    }
    
    suffix = prefix + 4 * ((length - prefix) / 4);
    
    for ( ; i < suffix; i += 4) {
      float64x2x2_t x_vec = vld1q_f64_x2(x + i);
      float64x2x2_t y_vec = vld1q_f64_x2(y + i);
      y_vec.val[0] = vaddq_f64(y_vec.val[0], vmulq_f64(x_vec.val[0], alpha_vec));
      y_vec.val[1] = vaddq_f64(y_vec.val[1], vmulq_f64(x_vec.val[1], alpha_vec));
      vst1q_f64_x2(y + i, y_vec);
    }
  } else {
    for ( ; i < suffix; i += 16) {
      vst1q_f64(y + i     , vaddq_f64(vld1q_f64(y + i     ), vmulq_f64(vld1q_f64(x + i     ), alpha_vec)));
      vst1q_f64(y + i +  2, vaddq_f64(vld1q_f64(y + i +  2), vmulq_f64(vld1q_f64(x + i +  2), alpha_vec)));
      vst1q_f64(y + i +  4, vaddq_f64(vld1q_f64(y + i +  4), vmulq_f64(vld1q_f64(x + i +  4), alpha_vec)));
      vst1q_f64(y + i +  6, vaddq_f64(vld1q_f64(y + i +  6), vmulq_f64(vld1q_f64(x + i +  6), alpha_vec)));
      
      vst1q_f64(y + i +  8, vaddq_f64(vld1q_f64(y + i +  8), vmulq_f64(vld1q_f64(x + i +  8), alpha_vec)));
      vst1q_f64(y + i + 10, vaddq_f64(vld1q_f64(y + i + 10), vmulq_f64(vld1q_f64(x + i + 10), alpha_vec)));
      vst1q_f64(y + i + 12, vaddq_f64(vld1q_f64(y + i + 12), vmulq_f64(vld1q_f64(x + i + 12), alpha_vec)));
      vst1q_f64(y + i + 14, vaddq_f64(vld1q_f64(y + i + 14), vmulq_f64(vld1q_f64(x + i + 14), alpha_vec)));
    }
    
    suffix = prefix + 4 * ((length - prefix) / 4);
    
    for ( ; i < suffix; i += 4) {
      vst1q_f64(y + i    , vaddq_f64(vld1q_f64(y + i    ), vmulq_f64(vld1q_f64(x + i    ), alpha_vec)));
      vst1q_f64(y + i + 2, vaddq_f64(vld1q_f64(y + i + 2), vmulq_f64(vld1q_f64(x + i + 2), alpha_vec)));
    }
  }
  
  for ( ; i < length; ++i)
    y[i] += alpha * x[i];
}

void misc_addAlignedVectorsInPlace_neon(const double* restrict x, size_t length, double* restrict y)
{
  if (length == 0) return;

  size_t i = 0;
  size_t suffix = 16 * (length / 16);
  
  for ( ; i < suffix; i += 16) {
    float64x2x4_t x_vec = vld1q_f64_x4(x + i);
    float64x2x4_t y_vec = vld1q_f64_x4(y + i);
    y_vec.val[0] = vaddq_f64(y_vec.val[0], x_vec.val[0]);
    y_vec.val[1] = vaddq_f64(y_vec.val[1], x_vec.val[1]);
    y_vec.val[2] = vaddq_f64(y_vec.val[2], x_vec.val[2]);
    y_vec.val[3] = vaddq_f64(y_vec.val[3], x_vec.val[3]);
    vst1q_f64_x4(y + i, y_vec);
    
    x_vec = vld1q_f64_x4(x + i + 8);
    y_vec = vld1q_f64_x4(y + i + 8);
    y_vec.val[0] = vaddq_f64(y_vec.val[0], x_vec.val[0]);
    y_vec.val[1] = vaddq_f64(y_vec.val[1], x_vec.val[1]);
    y_vec.val[2] = vaddq_f64(y_vec.val[2], x_vec.val[2]);
    y_vec.val[3] = vaddq_f64(y_vec.val[3], x_vec.val[3]);
    vst1q_f64_x4(y + i + 8, y_vec);
  }
  
  suffix = 4 * (length / 4);
  
  for ( ; i < suffix; i += 4) {
    float64x2x2_t x_vec = vld1q_f64_x2(x + i);
    float64x2x2_t y_vec = vld1q_f64_x2(y + i);
    y_vec.val[0] = vaddq_f64(y_vec.val[0], x_vec.val[0]);
    y_vec.val[1] = vaddq_f64(y_vec.val[1], x_vec.val[1]);
    vst1q_f64_x2(y + i, y_vec);
  }
  
  for ( ; i < length; ++i)
    y[i] += x[i];
}

void misc_subtractAlignedVectorsInPlace_neon(const double* restrict x, size_t length, double* restrict y)
{
  if (length == 0) return;

  size_t i = 0;
  size_t suffix = 16 * (length / 16);
  
  for ( ; i < suffix; i += 16) {
    float64x2x4_t x_vec = vld1q_f64_x4(x + i);
    float64x2x4_t y_vec = vld1q_f64_x4(y + i);
    y_vec.val[0] = vsubq_f64(y_vec.val[0], x_vec.val[0]);
    y_vec.val[1] = vsubq_f64(y_vec.val[1], x_vec.val[1]);
    y_vec.val[2] = vsubq_f64(y_vec.val[2], x_vec.val[2]);
    y_vec.val[3] = vsubq_f64(y_vec.val[3], x_vec.val[3]);
    vst1q_f64_x4(y + i, y_vec);

    x_vec = vld1q_f64_x4(x + i + 8);
    y_vec = vld1q_f64_x4(y + i + 8);
    y_vec.val[0] = vsubq_f64(y_vec.val[0], x_vec.val[0]);
    y_vec.val[1] = vsubq_f64(y_vec.val[1], x_vec.val[1]);
    y_vec.val[2] = vsubq_f64(y_vec.val[2], x_vec.val[2]);
    y_vec.val[3] = vsubq_f64(y_vec.val[3], x_vec.val[3]);
    vst1q_f64_x4(y + i + 8, y_vec);
  }
  
  suffix = 4 * (length / 4);
  
  for ( ; i < suffix; i += 4) {
    float64x2x2_t x_vec = vld1q_f64_x2(x + i);
    float64x2x2_t y_vec = vld1q_f64_x2(y + i);
    y_vec.val[0] = vsubq_f64(y_vec.val[0], x_vec.val[0]);
    y_vec.val[1] = vsubq_f64(y_vec.val[1], x_vec.val[1]);
    vst1q_f64_x2(y + i, y_vec);
  }

  for (; i < length; ++i)
    y[i] -= x[i];
}

void misc_addAlignedVectorsInPlaceWithMultiplier_neon(const double* restrict x, size_t length, double alpha, double* restrict y)
{
  if (length == 0) return;

  size_t i = 0;
  size_t suffix = 16 * (length / 16);
  
  float64x2_t alpha_vec = vdupq_n_f64(alpha);
  for ( ; i < suffix; i += 16) {
    float64x2x4_t x_vec = vld1q_f64_x4(x + i);
    float64x2x4_t y_vec = vld1q_f64_x4(y + i);
    y_vec.val[0] = vaddq_f64(y_vec.val[0], vmulq_f64(x_vec.val[0], alpha_vec));
    y_vec.val[1] = vaddq_f64(y_vec.val[1], vmulq_f64(x_vec.val[1], alpha_vec));
    y_vec.val[2] = vaddq_f64(y_vec.val[2], vmulq_f64(x_vec.val[2], alpha_vec));
    y_vec.val[3] = vaddq_f64(y_vec.val[3], vmulq_f64(x_vec.val[3], alpha_vec));
    vst1q_f64_x4(y + i, y_vec);
    
    x_vec = vld1q_f64_x4(x + i + 8);
    y_vec = vld1q_f64_x4(y + i + 8);
    y_vec.val[0] = vaddq_f64(y_vec.val[0], vmulq_f64(x_vec.val[0], alpha_vec));
    y_vec.val[1] = vaddq_f64(y_vec.val[1], vmulq_f64(x_vec.val[1], alpha_vec));
    y_vec.val[2] = vaddq_f64(y_vec.val[2], vmulq_f64(x_vec.val[2], alpha_vec));
    y_vec.val[3] = vaddq_f64(y_vec.val[3], vmulq_f64(x_vec.val[3], alpha_vec));
    vst1q_f64_x4(y + i + 8, y_vec);
  }
  
  suffix = 4 * (length / 4);
  
  for ( ; i < suffix; i += 4) {
    float64x2x2_t x_vec = vld1q_f64_x2(x + i);
    float64x2x2_t y_vec = vld1q_f64_x2(y + i);
    y_vec.val[0] = vaddq_f64(y_vec.val[0], vmulq_f64(x_vec.val[0], alpha_vec));
    y_vec.val[1] = vaddq_f64(y_vec.val[1], vmulq_f64(x_vec.val[1], alpha_vec));
    vst1q_f64_x2(y + i, y_vec);
  }

  for (; i < length; ++i)
    y[i] += alpha * x[i];
}

void misc_addScalarToVectorInPlace_neon(double* x, size_t length, double alpha)
{
  if (length == 0) return;
  
  size_t offset = ((uintptr_t) x) % (8 * sizeof(double));
  size_t prefix = offset == 0 ? 0 : (8 * sizeof(double) - offset) / sizeof(double);
  
  if (prefix > length) prefix = length;
  
  size_t i = 0;
  for ( ; i < prefix; ++i)
    x[i] += alpha;
  
  size_t suffix = prefix + 16 * ((length - prefix) / 16);
  
  float64x2_t alpha_vec = vdupq_n_f64(alpha);
  for ( ; i < suffix; i += 16) {
    float64x2x4_t x_vec = vld1q_f64_x4(x + i);
    x_vec.val[0] = vaddq_f64(x_vec.val[0], alpha_vec);
    x_vec.val[1] = vaddq_f64(x_vec.val[1], alpha_vec);
    x_vec.val[2] = vaddq_f64(x_vec.val[2], alpha_vec);
    x_vec.val[3] = vaddq_f64(x_vec.val[3], alpha_vec);

    vst1q_f64_x4(x + i, x_vec);
    
    x_vec = vld1q_f64_x4(x + i + 8);
    x_vec.val[0] = vaddq_f64(x_vec.val[0], alpha_vec);
    x_vec.val[1] = vaddq_f64(x_vec.val[1], alpha_vec);
    x_vec.val[2] = vaddq_f64(x_vec.val[2], alpha_vec);
    x_vec.val[3] = vaddq_f64(x_vec.val[3], alpha_vec);

    vst1q_f64_x4(x + i + 8, x_vec);
  }
  
  suffix = prefix + 4 * ((length - prefix) / 4);
  
  for ( ; i < suffix; i += 4) {
    float64x2x2_t x_vec = vld1q_f64_x2(x + i);
    x_vec.val[0] = vaddq_f64(x_vec.val[0], alpha_vec);
    x_vec.val[1] = vaddq_f64(x_vec.val[1], alpha_vec);
    vst1q_f64_x2(x + i, x_vec);
  }
  
  for ( ; i < length; ++i)
    x[i] += alpha;
}

void misc_setVectorToConstant_neon(double* x, size_t length, double alpha)
{
  if (length == 0) return;
  
  size_t offset = ((uintptr_t) x) % (8 * sizeof(double));
  size_t prefix = offset == 0 ? 0 : (8 * sizeof(double) - offset) / sizeof(double);
  
  if (prefix > length) prefix = length;
  
  size_t i = 0;
  for ( ; i < prefix; ++i)
    x[i] = alpha;
  
  size_t suffix = prefix + 16 * ((length - prefix) / 16);
  
  float64x2x4_t alpha_vec_x4;
  alpha_vec_x4.val[0] = vdupq_n_f64(alpha);
  alpha_vec_x4.val[1] = vdupq_n_f64(alpha);
  alpha_vec_x4.val[2] = vdupq_n_f64(alpha);
  alpha_vec_x4.val[3] = vdupq_n_f64(alpha);

  for ( ; i < suffix; i += 16) {
    vst1q_f64_x4(x + i    , alpha_vec_x4);
    vst1q_f64_x4(x + i + 8, alpha_vec_x4);
  }
  
  suffix = prefix + 4 * ((length - prefix) / 4);
  
  float64x2x2_t alpha_vec_x2;
  alpha_vec_x2.val[0] = vdupq_n_f64(alpha);
  alpha_vec_x2.val[1] = vdupq_n_f64(alpha);
  
  for ( ; i < suffix; i += 4) {
    vst1q_f64_x2(x + i, alpha_vec_x2);
  }
  
  for ( ; i < length; ++i)
    x[i] = alpha;
}


static inline void transposeMatrixBlock(const double* restrict x, size_t ldx, double* restrict y, size_t ldy)
{
  // x:  0  4  8 12
  //     1  5  9 13
  //     2  6 10 14
  //     3  7 11 15
  // y:  0  1  2  3
  //     4  5  6  7
  //     8  9 10 11
  //    12 13 14 15
  
  float64x2_t temp0, temp1;

  temp0 = vld1q_f64(x              ); // 0, 1
  temp1 = vld1q_f64(x +         ldx); // 4, 5
  vst1q_f64(y              , vtrn1q_f64(temp0, temp1)); // 0, 4
  vst1q_f64(y +         ldy, vtrn2q_f64(temp0, temp1)); // 1 5
  
  temp0 = vld1q_f64(x + 2          ); // 2, 3
  temp1 = vld1q_f64(x + 2 +     ldx); // 6, 7
  vst1q_f64(y +     2 * ldy, vtrn1q_f64(temp0, temp1));
  vst1q_f64(y +     3 * ldy, vtrn2q_f64(temp0, temp1));
  
  temp0 = vld1q_f64(x +     2 * ldx); // 8, 9
  temp1 = vld1q_f64(x +     3 * ldx); // 12, 13
  vst1q_f64(y + 2          , vtrn1q_f64(temp0, temp1));
  vst1q_f64(y + 2 +     ldy, vtrn2q_f64(temp0, temp1));
  
  temp0 = vld1q_f64(x + 2 + 2 * ldx); // 10, 11
  temp1 = vld1q_f64(x + 2 + 3 * ldx); // 14, 15
  vst1q_f64(y + 2 + 2 * ldy, vtrn1q_f64(temp0, temp1));
  vst1q_f64(y + 2 + 3 * ldy, vtrn2q_f64(temp0, temp1));
}

void misc_transposeMatrix_neon(const double* restrict x, size_t numRows, size_t numCols, double* restrict y)
{
  if (numRows == 0 || numCols == 0) return;
  
  // We can't really ensure that loads/stores occur on 16 byte boundaries, since
  // any time there is an odd number of rows that completely screws up
  // block transposing. For fun, we start x on a boundary and transpose
  // the first row explicitly, if necessary.
  size_t x_offset = ((uintptr_t) x) % (8 * sizeof(double));
  size_t prefix = x_offset == 0 ? 0 : (8 * sizeof(double) - x_offset) / sizeof(double);
  prefix = prefix > numRows ? numRows : prefix;

  size_t row = 0;
  for ( ; row < prefix; ++row) {
    for (size_t col = 0; col < numCols; ++col) {
      y[col + row * numCols] = x[row + col * numRows];
    }
  }
  
  size_t suffix = prefix + 4 * ((numRows - prefix) / 4);
  
  if (suffix > prefix) {
    for ( ; row < suffix; row += 4) {
      size_t col = 0, colEnd = 4 * (numCols / 4);
      for ( ; col < colEnd; col += 4)
        transposeMatrixBlock(x + row + col * numRows, numRows, y + col + row * numCols, numCols);
      
      for (size_t rowInBlock = row; rowInBlock < row + 4; ++rowInBlock) {
        size_t colInBlock = col;
        for ( ; colInBlock < numCols; ++colInBlock) {
          y[colInBlock + rowInBlock * numCols] = x[rowInBlock + colInBlock * numRows];
        }
      }
    }
  }
  
  for ( ; row < numRows; ++row) {
    for (size_t col = 0; col < numCols; ++col) {
      y[col + row * numCols] = x[row + col * numRows];
    }
  }
}


