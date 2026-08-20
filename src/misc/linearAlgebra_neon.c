#include "config.h"
#include <misc/linearAlgebra.h>

#include <stdint.h> // uintptr_t

#include <misc/intrinsic.h>

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
