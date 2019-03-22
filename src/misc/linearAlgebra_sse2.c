#include "config.h"
#include <misc/linearAlgebra.h>

#include <stdint.h> // uintptr_t

#include <misc/intrinsic.h>

void misc_addVectors_sse2(const double* restrict x, size_t length, double alpha, const double* restrict y, double* restrict z)
{
  if (length == 0) return;
  
  size_t y_offset = ((uintptr_t) y) % (2 * sizeof(double));
  size_t x_offset = ((uintptr_t) x) % (2 * sizeof(double));
  size_t z_offset = ((uintptr_t) z) % (2 * sizeof(double));
  size_t prefix = z_offset == 0 ? 0 : (2 * sizeof(double) - z_offset) / sizeof(double);
  
  if (prefix > length) prefix = length;
  
  size_t i = 0;
  for ( ; i < prefix; ++i)
    z[i] = y[i] + alpha * x[i];
  
  size_t suffix = prefix + 6 * ((length - prefix) / 6);
  
  if (suffix > prefix) {
    __m128d alpha_vec = _mm_load1_pd(&alpha);
    if (z_offset == x_offset && z_offset == y_offset) {
      for ( ; i < suffix; i += 6) {
        _mm_stream_pd(z + i    , _mm_add_pd(_mm_load_pd(y + i    ), _mm_mul_pd(_mm_load_pd(x + i    ), alpha_vec)));
        _mm_stream_pd(z + i + 2, _mm_add_pd(_mm_load_pd(y + i + 2), _mm_mul_pd(_mm_load_pd(x + i + 2), alpha_vec)));
        _mm_stream_pd(z + i + 4, _mm_add_pd(_mm_load_pd(y + i + 4), _mm_mul_pd(_mm_load_pd(x + i + 4), alpha_vec)));
      }
    } else {
      for ( ; i < suffix; i += 6) {
        _mm_stream_pd(z + i    , _mm_add_pd(_mm_loadu_pd(y + i    ), _mm_mul_pd(_mm_loadu_pd(x + i    ), alpha_vec)));
        _mm_stream_pd(z + i + 2, _mm_add_pd(_mm_loadu_pd(y + i + 2), _mm_mul_pd(_mm_loadu_pd(x + i + 2), alpha_vec)));
        _mm_stream_pd(z + i + 4, _mm_add_pd(_mm_loadu_pd(y + i + 4), _mm_mul_pd(_mm_loadu_pd(x + i + 4), alpha_vec)));
      }
    }
  }
  
  for ( ; i < length; ++i)
    z[i] = y[i] + alpha * x[i];
}

extern void misc_setVectorToConstant_c(double* x, size_t length, double alpha);

void misc_setVectorToConstant_sse2(double* x, size_t length, double alpha)
{
  if (length == 0) return;
  if (((uintptr_t) x) % sizeof(double) != 0) {
    misc_setVectorToConstant_c(x, length, alpha);
    return;
  }
  
  size_t offset = ((uintptr_t) x) % (2 * sizeof(double));
  size_t prefix = offset == 0 ? 0 : (2 * sizeof(double) - offset) / sizeof(double);
  
  if (prefix > length) prefix = length;
  
  size_t i = 0;
  for ( ; i < prefix; ++i)
    x[i] = alpha;
  
  size_t suffix = prefix + 8 * ((length - prefix) / 8);
  
  if (suffix > prefix) {
    __m128d alpha_vec = _mm_load1_pd(&alpha);
    for ( ; i < suffix; i += 8) {
      _mm_stream_pd(x + i    , alpha_vec);
      _mm_stream_pd(x + i + 2, alpha_vec);
      _mm_stream_pd(x + i + 4, alpha_vec);
      _mm_stream_pd(x + i + 6, alpha_vec);
    }
  }
  
  for ( ; i < length; ++i)
    x[i] = alpha;
}

