#include "config.h"
#include <misc/linearAlgebra.h>

#include <stdint.h> // uintptr_t

#include <misc/alloca.h>
#include <misc/intrinsic.h>

void misc_addVectors_avx(const double* restrict x, size_t length, double alpha, const double* restrict y, double* restrict z)
{
  if (length == 0) return;
  
  size_t y_offset = ((uintptr_t) y) % (4 * sizeof(double));
  size_t x_offset = ((uintptr_t) x) % (4 * sizeof(double));
  size_t z_offset = ((uintptr_t) z) % (4 * sizeof(double));
  size_t prefix = z_offset == 0 ? 0 : (4 * sizeof(double) - z_offset) / sizeof(double);
  
  if (prefix > length) prefix = length;
    
  size_t i = 0;
  for ( ; i < prefix; ++i) z[i] = y[i] + alpha * x[i];
    
  size_t suffix = prefix + 12 * ((length - prefix) / 12);

  if (suffix > prefix) {
    __m256d alpha_vec = _mm256_set_pd(alpha, alpha, alpha, alpha);
    
    if (x_offset == y_offset && x_offset == z_offset) {
      for ( ; i < suffix; i += 12) {
        _mm256_stream_pd(z + i    , _mm256_add_pd(_mm256_load_pd(y + i    ), _mm256_mul_pd(_mm256_load_pd(x + i    ), alpha_vec)));
        _mm256_stream_pd(z + i + 4, _mm256_add_pd(_mm256_load_pd(y + i + 4), _mm256_mul_pd(_mm256_load_pd(x + i + 4), alpha_vec)));
        _mm256_stream_pd(z + i + 8, _mm256_add_pd(_mm256_load_pd(y + i + 8), _mm256_mul_pd(_mm256_load_pd(x + i + 8), alpha_vec)));
      }
    } else {
      for ( ; i < suffix; i += 12) {
        _mm256_stream_pd(z + i    , _mm256_add_pd(_mm256_loadu_pd(y + i    ), _mm256_mul_pd(_mm256_loadu_pd(x + i    ), alpha_vec)));
        _mm256_stream_pd(z + i + 4, _mm256_add_pd(_mm256_loadu_pd(y + i + 4), _mm256_mul_pd(_mm256_loadu_pd(x + i + 4), alpha_vec)));
        _mm256_stream_pd(z + i + 8, _mm256_add_pd(_mm256_loadu_pd(y + i + 8), _mm256_mul_pd(_mm256_loadu_pd(x + i + 8), alpha_vec)));
      }
    }
      
    misc_stackFree(alpha_void);
  }
  
  for ( ; i < length; ++i)
    z[i] = y[i] + alpha * x[i];
}

extern void misc_setVectorToConstant_c(double* x, size_t length, double alpha);

void misc_setVectorToConstant_avx(double* x, size_t length, double alpha)
{
  if (length == 0) return;
  if (((uintptr_t) x) % sizeof(double) != 0) {
    misc_setVectorToConstant_c(x, length, alpha);
    return;
  }
  
  size_t offset = ((uintptr_t) x) % (4 * sizeof(double));
  size_t prefix = offset == 0 ? 0 : (4 * sizeof(double) - offset) / sizeof(double);
  
  if (prefix > length) prefix = length;
  
  size_t i = 0;
  for ( ; i < prefix; ++i)
    x[i] = alpha;
  
  size_t suffix = prefix + 16 * ((length - prefix) / 16);
  
  if (suffix > prefix) {
    __m256d alpha_vec = _mm256_set_pd(alpha, alpha, alpha, alpha);
  
    for ( ; i < suffix; i += 16) {
      _mm256_stream_pd(x + i     , alpha_vec);
      _mm256_stream_pd(x + i +  4, alpha_vec);
      _mm256_stream_pd(x + i +  8, alpha_vec);
      _mm256_stream_pd(x + i + 12, alpha_vec);
    }
    
    misc_stackFree(alpha_void);
  }
  
  for ( ; i < length; ++i)
    x[i] = alpha;
}

