#include "config.h"
#include <misc/linearAlgebra.h>

#include <stdint.h> // uintptr_t

#include <misc/intrinsic.h>

#include <external/io.h>

void misc_addVectors_neon(const double* restrict x, size_t length, double alpha, const double* restrict y, double* restrict z)
{
  if (length == 0) return;
  
  size_t y_offset = ((uintptr_t) y) % (8 * sizeof(double));
  size_t x_offset = ((uintptr_t) x) % (8 * sizeof(double));
  size_t z_offset = ((uintptr_t) z) % (8 * sizeof(double));
  size_t prefix = z_offset == 0 ? 0 : (8 * sizeof(double) - z_offset) / sizeof(double);

  if (prefix > length) prefix = length;
  
  size_t i = 0;
  for ( ; i < prefix; ++i)
    z[i] = y[i] + alpha * x[i];
  

  size_t suffix = prefix + 8 * ((length - prefix) / 8);

  if (suffix > prefix) {
    float64x2_t alpha_vec = vdupq_n_f64(alpha);
    if (z_offset == x_offset && z_offset == y_offset) {
      for ( ; i < suffix; i += 8) {
        float64x2x4_t x_vec = vld1q_f64_x4(x + i);
        float64x2x4_t y_vec = vld1q_f64_x4(y + i);
        y_vec.val[0] = vaddq_f64(y_vec.val[0], vmulq_f64(x_vec.val[0], alpha_vec));
        y_vec.val[1] = vaddq_f64(y_vec.val[1], vmulq_f64(x_vec.val[1], alpha_vec));
        y_vec.val[2] = vaddq_f64(y_vec.val[2], vmulq_f64(x_vec.val[2], alpha_vec));
        y_vec.val[3] = vaddq_f64(y_vec.val[3], vmulq_f64(x_vec.val[3], alpha_vec));
        vst1q_f64_x4(z + i, y_vec);
      }
    } else {
      for ( ; i < suffix; i += 8) {
        vst1q_f64(z + i    , vaddq_f64(vld1q_f64(y + i    ), vmulq_f64(vld1q_f64(x + i    ), alpha_vec)));
        vst1q_f64(z + i + 2, vaddq_f64(vld1q_f64(y + i + 2), vmulq_f64(vld1q_f64(x + i + 2), alpha_vec)));
        vst1q_f64(z + i + 4, vaddq_f64(vld1q_f64(y + i + 4), vmulq_f64(vld1q_f64(x + i + 4), alpha_vec)));
        vst1q_f64(z + i + 6, vaddq_f64(vld1q_f64(y + i + 6), vmulq_f64(vld1q_f64(x + i + 6), alpha_vec)));
      }
    }
  }
  
  for ( ; i < length; ++i)
    z[i] = y[i] + alpha * x[i];
}

void misc_addVectorsInPlace_neon(const double* restrict x, size_t length, double alpha, double* restrict y)
{
  if (length == 0) return;

  size_t y_offset = ((uintptr_t) y) % (8 * sizeof(double));
  size_t x_offset = ((uintptr_t) x) % (8 * sizeof(double));
  size_t prefix = y_offset == 0 ? 0 : (8 * sizeof(double) - y_offset) / sizeof(double);
  
  if (prefix > length) prefix = length;
  
  size_t i = 0;
  for ( ; i < prefix; ++i)
    y[i] += alpha * x[i];
  
  size_t suffix = prefix + 8 * ((length - prefix) / 8);
  
  if (suffix > prefix) {
    float64x2_t alpha_vec = vdupq_n_f64(alpha);
    if (y_offset == x_offset) {
      for ( ; i < suffix; i += 8) {
        float64x2x4_t x_vec = vld1q_f64_x4(x + i);
        float64x2x4_t y_vec = vld1q_f64_x4(y + i);
        y_vec.val[0] = vaddq_f64(y_vec.val[0], vmulq_f64(x_vec.val[0], alpha_vec));
        y_vec.val[1] = vaddq_f64(y_vec.val[1], vmulq_f64(x_vec.val[1], alpha_vec));
        y_vec.val[2] = vaddq_f64(y_vec.val[2], vmulq_f64(x_vec.val[2], alpha_vec));
        y_vec.val[3] = vaddq_f64(y_vec.val[3], vmulq_f64(x_vec.val[3], alpha_vec));
        vst1q_f64_x4(y + i, y_vec);
      }
    } else {
      for ( ; i < suffix; i += 8) {
        vst1q_f64(y + i    , vaddq_f64(vld1q_f64(y + i    ), vmulq_f64(vld1q_f64(x + i    ), alpha_vec)));
        vst1q_f64(y + i + 2, vaddq_f64(vld1q_f64(y + i + 2), vmulq_f64(vld1q_f64(x + i + 2), alpha_vec)));
        vst1q_f64(y + i + 4, vaddq_f64(vld1q_f64(y + i + 4), vmulq_f64(vld1q_f64(x + i + 4), alpha_vec)));
        vst1q_f64(y + i + 6, vaddq_f64(vld1q_f64(y + i + 6), vmulq_f64(vld1q_f64(x + i + 6), alpha_vec)));
      }
    }
  }
  
  for ( ; i < length; ++i)
    y[i] += alpha * x[i];
}

extern void misc_addScalarToVectorInPlace_c(double* x, size_t length, double alpha);

void misc_addScalarToVectorInPlace_neon(double* x, size_t length, double alpha)
{
  if (length == 0) return;
  if (((uintptr_t) x) % sizeof(double) != 0) {
    misc_addScalarToVectorInPlace_c(x, length, alpha);
    return;
  }
  
  size_t offset = ((uintptr_t) x) % (8 * sizeof(double));
  size_t prefix = offset == 0 ? 0 : (8 * sizeof(double) - offset) / sizeof(double);
  
  if (prefix > length) prefix = length;
  
  size_t i = 0;
  for ( ; i < prefix; ++i)
    x[i] += alpha;
  
  size_t suffix = prefix + 8 * ((length - prefix) / 8);
  
  if (suffix > prefix) {
    float64x2_t alpha_vec = vdupq_n_f64(alpha);
    for ( ; i < suffix; i += 8) {
      float64x2x4_t x_vec = vld1q_f64_x4(x + i);
      x_vec.val[0] = vaddq_f64(x_vec.val[0], alpha_vec);
      x_vec.val[1] = vaddq_f64(x_vec.val[1], alpha_vec);
      x_vec.val[2] = vaddq_f64(x_vec.val[2], alpha_vec);
      x_vec.val[3] = vaddq_f64(x_vec.val[3], alpha_vec);

      vst1q_f64_x4(x + i, x_vec);
    }
  }
  
  for ( ; i < length; ++i)
    x[i] += alpha;
}

extern void misc_setVectorToConstant_c(double* x, size_t length, double alpha);

void misc_setVectorToConstant_neon(double* x, size_t length, double alpha)
{
  if (length == 0) return;
  if (((uintptr_t) x) % sizeof(double) != 0) {
    misc_setVectorToConstant_c(x, length, alpha);
    return;
  }
  
  size_t offset = ((uintptr_t) x) % (8 * sizeof(double));
  size_t prefix = offset == 0 ? 0 : (8 * sizeof(double) - offset) / sizeof(double);
  
  if (prefix > length) prefix = length;
  
  size_t i = 0;
  for ( ; i < prefix; ++i)
    x[i] = alpha;
  
  size_t suffix = prefix + 8 * ((length - prefix) / 8);
  
  if (suffix > prefix) {
    float64x2x4_t alpha_vec;
    alpha_vec.val[0] = vdupq_n_f64(alpha);
    alpha_vec.val[1] = vdupq_n_f64(alpha);
    alpha_vec.val[2] = vdupq_n_f64(alpha);
    alpha_vec.val[3] = vdupq_n_f64(alpha);

    for ( ; i < suffix; i += 8) {
      vst1q_f64_x4(x + i, alpha_vec);
    }
  }
  
  for ( ; i < length; ++i)
    x[i] = alpha;
}


