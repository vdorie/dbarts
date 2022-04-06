#include "config.h"
#include <misc/linearAlgebra.h>

#include <stdint.h> // uintptr_t

#include <misc/alloca.h>
#include <misc/intrinsic.h>

void misc_addVectors_avx(const double* restrict x, size_t length, const double* restrict y, double* restrict z)
{
  if (length == 0) return;
  
  size_t y_offset = ((uintptr_t) y) % (4 * sizeof(double));
  size_t x_offset = ((uintptr_t) x) % (4 * sizeof(double));
  size_t z_offset = ((uintptr_t) z) % (4 * sizeof(double));
  size_t prefix = z_offset == 0 ? 0 : (4 * sizeof(double) - z_offset) / sizeof(double);
  
  if (prefix > length) prefix = length;
    
  size_t i = 0;
  size_t suffix = prefix + 16 * ((length - prefix) / 16);
  
  for ( ; i < prefix; ++i)
    z[i] = y[i] + x[i];
    
  if (z_offset == x_offset && z_offset == y_offset) {
    for ( ; i < suffix; i += 16) {
      _mm256_stream_pd(z + i     , _mm256_add_pd(_mm256_load_pd(y + i     ), _mm256_load_pd(x + i     )));
      _mm256_stream_pd(z + i +  4, _mm256_add_pd(_mm256_load_pd(y + i +  4), _mm256_load_pd(x + i +  4)));
      _mm256_stream_pd(z + i +  8, _mm256_add_pd(_mm256_load_pd(y + i +  8), _mm256_load_pd(x + i +  8)));
      _mm256_stream_pd(z + i + 12, _mm256_add_pd(_mm256_load_pd(y + i + 12), _mm256_load_pd(x + i + 12)));
    }
    
    suffix = prefix + 4 * ((length - prefix) / 4);
    
    for ( ; i < suffix; i += 4) {
      _mm256_stream_pd(z + i     , _mm256_add_pd(_mm256_load_pd(y + i     ), _mm256_load_pd(x + i     )));
    }
  } else {
    for ( ; i < suffix; i += 16) {
      _mm256_stream_pd(z + i     , _mm256_add_pd(_mm256_loadu_pd(y + i     ), _mm256_loadu_pd(x + i     )));
      _mm256_stream_pd(z + i +  4, _mm256_add_pd(_mm256_loadu_pd(y + i +  4), _mm256_loadu_pd(x + i +  4)));
      _mm256_stream_pd(z + i +  8, _mm256_add_pd(_mm256_loadu_pd(y + i +  8), _mm256_loadu_pd(x + i +  8)));
      _mm256_stream_pd(z + i + 12, _mm256_add_pd(_mm256_loadu_pd(y + i + 12), _mm256_loadu_pd(x + i + 12)));
    }
    
    suffix = prefix + 4 * ((length - prefix) / 4);
    
    for ( ; i < suffix; i += 4) {
      _mm256_stream_pd(z + i     , _mm256_add_pd(_mm256_loadu_pd(y + i     ), _mm256_loadu_pd(x + i     )));
    }
  }
  
  for ( ; i < length; ++i)
    z[i] = y[i] + x[i];
}

void misc_subtractVectors_avx(const double* restrict x, size_t length, const double* restrict y, double* restrict z)
{
  if (length == 0) return;
  
  size_t y_offset = ((uintptr_t) y) % (4 * sizeof(double));
  size_t x_offset = ((uintptr_t) x) % (4 * sizeof(double));
  size_t z_offset = ((uintptr_t) z) % (4 * sizeof(double));
  size_t prefix = z_offset == 0 ? 0 : (4 * sizeof(double) - z_offset) / sizeof(double);
  
  if (prefix > length) prefix = length;
    
  size_t i = 0;
  size_t suffix = prefix + 16 * ((length - prefix) / 16);

  for ( ; i < prefix; ++i)
    z[i] = y[i] - x[i];
  
  if (z_offset == x_offset && z_offset == y_offset) {
    for ( ; i < suffix; i += 16) {
      _mm256_stream_pd(z + i     , _mm256_sub_pd(_mm256_load_pd(y + i     ), _mm256_load_pd(x + i     )));
      _mm256_stream_pd(z + i +  4, _mm256_sub_pd(_mm256_load_pd(y + i +  4), _mm256_load_pd(x + i +  4)));
      _mm256_stream_pd(z + i +  8, _mm256_sub_pd(_mm256_load_pd(y + i +  8), _mm256_load_pd(x + i +  8)));
      _mm256_stream_pd(z + i + 12, _mm256_sub_pd(_mm256_load_pd(y + i + 12), _mm256_load_pd(x + i + 12)));
    }
    
    suffix = prefix + 4 * ((length - prefix) / 4);
    
    for ( ; i < suffix; i += 4) {
      _mm256_stream_pd(z + i     , _mm256_sub_pd(_mm256_load_pd(y + i     ), _mm256_load_pd(x + i     )));
    }
  } else {
    for ( ; i < suffix; i += 16) {
      _mm256_stream_pd(z + i     , _mm256_sub_pd(_mm256_loadu_pd(y + i     ), _mm256_loadu_pd(x + i     )));
      _mm256_stream_pd(z + i +  4, _mm256_sub_pd(_mm256_loadu_pd(y + i +  4), _mm256_loadu_pd(x + i +  4)));
      _mm256_stream_pd(z + i +  8, _mm256_sub_pd(_mm256_loadu_pd(y + i +  8), _mm256_loadu_pd(x + i +  8)));
      _mm256_stream_pd(z + i + 12, _mm256_sub_pd(_mm256_loadu_pd(y + i + 12), _mm256_loadu_pd(x + i + 12)));
    }
    
    suffix = prefix + 4 * ((length - prefix) / 4);
    
    for ( ; i < suffix; i += 4) {
      _mm256_stream_pd(z + i     , _mm256_sub_pd(_mm256_loadu_pd(y + i     ), _mm256_loadu_pd(x + i     )));
    }
  }
  
  for ( ; i < length; ++i)
    z[i] = y[i] - x[i];
}

void misc_addVectorsWithMultiplier_avx(const double* restrict x, size_t length, double alpha, const double* restrict y, double* restrict z)
{
  if (length == 0) return;
  
  size_t y_offset = ((uintptr_t) y) % (4 * sizeof(double));
  size_t x_offset = ((uintptr_t) x) % (4 * sizeof(double));
  size_t z_offset = ((uintptr_t) z) % (4 * sizeof(double));
  size_t prefix = z_offset == 0 ? 0 : (4 * sizeof(double) - z_offset) / sizeof(double);
  
  if (prefix > length) prefix = length;
    
  size_t i = 0;
  size_t suffix = prefix + 16 * ((length - prefix) / 16);

  for ( ; i < prefix; ++i)
    z[i] = y[i] + alpha * x[i];
  
  __m256d alpha_vec = _mm256_set_pd(alpha, alpha, alpha, alpha);
  
  if (z_offset == x_offset && z_offset == y_offset) {
    for ( ; i < suffix; i += 16) {
      _mm256_stream_pd(z + i     , _mm256_add_pd(_mm256_load_pd(y + i     ), _mm256_mul_pd(_mm256_load_pd(x + i     ), alpha_vec)));
      _mm256_stream_pd(z + i +  4, _mm256_add_pd(_mm256_load_pd(y + i +  4), _mm256_mul_pd(_mm256_load_pd(x + i +  4), alpha_vec)));
      _mm256_stream_pd(z + i +  8, _mm256_add_pd(_mm256_load_pd(y + i +  8), _mm256_mul_pd(_mm256_load_pd(x + i +  8), alpha_vec)));
      _mm256_stream_pd(z + i + 12, _mm256_add_pd(_mm256_load_pd(y + i + 12), _mm256_mul_pd(_mm256_load_pd(x + i + 12), alpha_vec)));
    }
    
    suffix = prefix + 4 * ((length - prefix) / 4);
    
    for ( ; i < suffix; i += 4) {
      _mm256_stream_pd(z + i     , _mm256_add_pd(_mm256_load_pd(y + i     ), _mm256_mul_pd(_mm256_load_pd(x + i     ), alpha_vec)));
    }
  } else {
    for ( ; i < suffix; i += 16) {
      _mm256_stream_pd(z + i     , _mm256_add_pd(_mm256_loadu_pd(y + i     ), _mm256_mul_pd(_mm256_loadu_pd(x + i     ), alpha_vec)));
      _mm256_stream_pd(z + i +  4, _mm256_add_pd(_mm256_loadu_pd(y + i +  4), _mm256_mul_pd(_mm256_loadu_pd(x + i +  4), alpha_vec)));
      _mm256_stream_pd(z + i +  8, _mm256_add_pd(_mm256_loadu_pd(y + i +  8), _mm256_mul_pd(_mm256_loadu_pd(x + i +  8), alpha_vec)));
      _mm256_stream_pd(z + i + 12, _mm256_add_pd(_mm256_loadu_pd(y + i + 12), _mm256_mul_pd(_mm256_loadu_pd(x + i + 12), alpha_vec)));
    }
    
    suffix = prefix + 4 * ((length - prefix) / 4);
    
    for ( ; i < suffix; i += 4) {
      _mm256_stream_pd(z + i     , _mm256_add_pd(_mm256_loadu_pd(y + i     ), _mm256_mul_pd(_mm256_loadu_pd(x + i     ), alpha_vec)));
    }
  }
  
  for ( ; i < length; ++i)
    z[i] = y[i] + alpha * x[i];
}

void misc_addVectorsInPlace_avx(const double* restrict x, size_t length, double* restrict y)
{
  if (length == 0) return;
  
  size_t y_offset = ((uintptr_t) y) % (4 * sizeof(double));
  size_t x_offset = ((uintptr_t) x) % (4 * sizeof(double));
  size_t prefix = y_offset == 0 ? 0 : (4 * sizeof(double) - y_offset) / sizeof(double);
  
  if (prefix > length) prefix = length;
    
  size_t i = 0;
  size_t suffix = prefix + 16 * ((length - prefix) / 16);
  
  for ( ; i < prefix; ++i)
    y[i] += x[i];
    
  if (y_offset == x_offset) {
    for ( ; i < suffix; i += 16) {
      _mm256_stream_pd(y + i     , _mm256_add_pd(_mm256_load_pd(y + i     ), _mm256_load_pd(x + i     )));
      _mm256_stream_pd(y + i +  4, _mm256_add_pd(_mm256_load_pd(y + i +  4), _mm256_load_pd(x + i +  4)));
      _mm256_stream_pd(y + i +  8, _mm256_add_pd(_mm256_load_pd(y + i +  8), _mm256_load_pd(x + i +  8)));
      _mm256_stream_pd(y + i + 12, _mm256_add_pd(_mm256_load_pd(y + i + 12), _mm256_load_pd(x + i + 12)));
    }
    
    suffix = prefix + 4 * ((length - prefix) / 4);
    
    for ( ; i < suffix; i += 4) {
      _mm256_stream_pd(y + i     , _mm256_add_pd(_mm256_load_pd(y + i     ), _mm256_load_pd(x + i     )));
    }
  } else {
    for ( ; i < suffix; i += 16) {
      _mm256_stream_pd(y + i     , _mm256_add_pd(_mm256_loadu_pd(y + i     ), _mm256_loadu_pd(x + i     )));
      _mm256_stream_pd(y + i +  4, _mm256_add_pd(_mm256_loadu_pd(y + i +  4), _mm256_loadu_pd(x + i +  4)));
      _mm256_stream_pd(y + i +  8, _mm256_add_pd(_mm256_loadu_pd(y + i +  8), _mm256_loadu_pd(x + i +  8)));
      _mm256_stream_pd(y + i + 12, _mm256_add_pd(_mm256_loadu_pd(y + i + 12), _mm256_loadu_pd(x + i + 12)));
    }
    
    suffix = prefix + 4 * ((length - prefix) / 4);
    
    for ( ; i < suffix; i += 4) {
      _mm256_stream_pd(y + i     , _mm256_add_pd(_mm256_loadu_pd(y + i     ), _mm256_loadu_pd(x + i     )));
    }
  }
  
  for ( ; i < length; ++i)
    y[i] += x[i];
}

void misc_subtractVectorsInPlace_avx(const double* restrict x, size_t length, double* restrict y)
{
  if (length == 0) return;
  
  size_t y_offset = ((uintptr_t) y) % (4 * sizeof(double));
  size_t x_offset = ((uintptr_t) x) % (4 * sizeof(double));
  size_t prefix = y_offset == 0 ? 0 : (4 * sizeof(double) - y_offset) / sizeof(double);
  
  if (prefix > length) prefix = length;
    
  size_t i = 0;
  size_t suffix = prefix + 16 * ((length - prefix) / 16);

  for ( ; i < prefix; ++i)
    y[i] -= x[i];
  
  if (y_offset == x_offset) {
    for ( ; i < suffix; i += 16) {
      _mm256_stream_pd(y + i     , _mm256_sub_pd(_mm256_load_pd(y + i     ), _mm256_load_pd(x + i     )));
      _mm256_stream_pd(y + i +  4, _mm256_sub_pd(_mm256_load_pd(y + i +  4), _mm256_load_pd(x + i +  4)));
      _mm256_stream_pd(y + i +  8, _mm256_sub_pd(_mm256_load_pd(y + i +  8), _mm256_load_pd(x + i +  8)));
      _mm256_stream_pd(y + i + 12, _mm256_sub_pd(_mm256_load_pd(y + i + 12), _mm256_load_pd(x + i + 12)));
    }
    
    suffix = prefix + 4 * ((length - prefix) / 4);
    
    for ( ; i < suffix; i += 4) {
      _mm256_stream_pd(y + i     , _mm256_sub_pd(_mm256_load_pd(y + i     ), _mm256_load_pd(x + i     )));
    }
  } else {
    for ( ; i < suffix; i += 16) {
      _mm256_stream_pd(y + i     , _mm256_sub_pd(_mm256_loadu_pd(y + i     ), _mm256_loadu_pd(x + i     )));
      _mm256_stream_pd(y + i +  4, _mm256_sub_pd(_mm256_loadu_pd(y + i +  4), _mm256_loadu_pd(x + i +  4)));
      _mm256_stream_pd(y + i +  8, _mm256_sub_pd(_mm256_loadu_pd(y + i +  8), _mm256_loadu_pd(x + i +  8)));
      _mm256_stream_pd(y + i + 12, _mm256_sub_pd(_mm256_loadu_pd(y + i + 12), _mm256_loadu_pd(x + i + 12)));
    }
    
    suffix = prefix + 4 * ((length - prefix) / 4);
    
    for ( ; i < suffix; i += 4) {
      _mm256_stream_pd(y + i     , _mm256_sub_pd(_mm256_loadu_pd(y + i     ), _mm256_loadu_pd(x + i     )));
    }
  }
  
  for ( ; i < length; ++i)
    y[i] -= x[i];
}

void misc_addVectorsInPlaceWithMultiplier_avx(const double* restrict x, size_t length, double alpha, double* restrict y)
{
  if (length == 0) return;
  
  size_t y_offset = ((uintptr_t) y) % (4 * sizeof(double));
  size_t x_offset = ((uintptr_t) x) % (4 * sizeof(double));
  size_t prefix = y_offset == 0 ? 0 : (4 * sizeof(double) - y_offset) / sizeof(double);
  
  if (prefix > length) prefix = length;
    
  size_t i = 0;
  size_t suffix = prefix + 16 * ((length - prefix) / 16);

  for ( ; i < prefix; ++i)
    y[i] += alpha * x[i];
  
  __m256d alpha_vec = _mm256_set_pd(alpha, alpha, alpha, alpha);
  
  if (y_offset == x_offset) {
    for ( ; i < suffix; i += 16) {
      _mm256_stream_pd(y + i     , _mm256_add_pd(_mm256_load_pd(y + i     ), _mm256_mul_pd(_mm256_load_pd(x + i     ), alpha_vec)));
      _mm256_stream_pd(y + i +  4, _mm256_add_pd(_mm256_load_pd(y + i +  4), _mm256_mul_pd(_mm256_load_pd(x + i +  4), alpha_vec)));
      _mm256_stream_pd(y + i +  8, _mm256_add_pd(_mm256_load_pd(y + i +  8), _mm256_mul_pd(_mm256_load_pd(x + i +  8), alpha_vec)));
      _mm256_stream_pd(y + i + 12, _mm256_add_pd(_mm256_load_pd(y + i + 12), _mm256_mul_pd(_mm256_load_pd(x + i + 12), alpha_vec)));
    }
    
    suffix = prefix + 4 * ((length - prefix) / 4);
    
    for ( ; i < suffix; i += 4) {
      _mm256_stream_pd(y + i     , _mm256_add_pd(_mm256_load_pd(y + i     ), _mm256_mul_pd(_mm256_load_pd(x + i     ), alpha_vec)));
    }
  } else {
    for ( ; i < suffix; i += 16) {
      _mm256_stream_pd(y + i     , _mm256_add_pd(_mm256_loadu_pd(y + i     ), _mm256_mul_pd(_mm256_loadu_pd(x + i     ), alpha_vec)));
      _mm256_stream_pd(y + i +  4, _mm256_add_pd(_mm256_loadu_pd(y + i +  4), _mm256_mul_pd(_mm256_loadu_pd(x + i +  4), alpha_vec)));
      _mm256_stream_pd(y + i +  8, _mm256_add_pd(_mm256_loadu_pd(y + i +  8), _mm256_mul_pd(_mm256_loadu_pd(x + i +  8), alpha_vec)));
      _mm256_stream_pd(y + i + 12, _mm256_add_pd(_mm256_loadu_pd(y + i + 12), _mm256_mul_pd(_mm256_loadu_pd(x + i + 12), alpha_vec)));
    }
    
    suffix = prefix + 4 * ((length - prefix) / 4);
    
    for ( ; i < suffix; i += 4) {
      _mm256_stream_pd(y + i     , _mm256_add_pd(_mm256_loadu_pd(y + i     ), _mm256_mul_pd(_mm256_loadu_pd(x + i     ), alpha_vec)));
    }
  }
  
  for ( ; i < length; ++i)
    y[i] += alpha * x[i];
}

void misc_addAlignedVectorsInPlace_avx(const double* restrict x, size_t length, double* restrict y)
{
  if (length == 0) return;
  
  size_t i = 0;
  size_t suffix = 16 * (length / 16);

  for ( ; i < suffix; i += 16) {
    _mm256_stream_pd(y + i     , _mm256_add_pd(_mm256_load_pd(y + i     ), _mm256_load_pd(x + i     )));
    _mm256_stream_pd(y + i +  4, _mm256_add_pd(_mm256_load_pd(y + i +  4), _mm256_load_pd(x + i +  4)));
    _mm256_stream_pd(y + i +  8, _mm256_add_pd(_mm256_load_pd(y + i +  8), _mm256_load_pd(x + i +  8)));
    _mm256_stream_pd(y + i + 12, _mm256_add_pd(_mm256_load_pd(y + i + 12), _mm256_load_pd(x + i + 12)));
  }
  
  suffix = 4 * (length / 4);
  
  for ( ; i < suffix; i += 4) {
    _mm256_stream_pd(y + i     , _mm256_add_pd(_mm256_load_pd(y + i     ), _mm256_load_pd(x + i     )));
  }
  
  for ( ; i < length; ++i)
    y[i] += x[i];
}

void misc_subtractAlignedVectorsInPlace_avx(const double* restrict x, size_t length, double* restrict y)
{
  if (length == 0) return;
  
  size_t i = 0;
  size_t suffix = 16 * (length / 16);

  for ( ; i < suffix; i += 16) {
    _mm256_stream_pd(y + i     , _mm256_sub_pd(_mm256_load_pd(y + i     ), _mm256_load_pd(x + i     )));
    _mm256_stream_pd(y + i +  4, _mm256_sub_pd(_mm256_load_pd(y + i +  4), _mm256_load_pd(x + i +  4)));
    _mm256_stream_pd(y + i +  8, _mm256_sub_pd(_mm256_load_pd(y + i +  8), _mm256_load_pd(x + i +  8)));
    _mm256_stream_pd(y + i + 12, _mm256_sub_pd(_mm256_load_pd(y + i + 12), _mm256_load_pd(x + i + 12)));
  }
  
  suffix = 4 * (length / 4);
  
  for ( ; i < suffix; i += 4) {
    _mm256_stream_pd(y + i     , _mm256_sub_pd(_mm256_load_pd(y + i     ), _mm256_load_pd(x + i     )));
  }
  
  for ( ; i < length; ++i)
    y[i] -= x[i];
}

void misc_addAlignedVectorsInPlaceWithMultiplier_avx(const double* restrict x, size_t length, double alpha, double* restrict y)
{
  if (length == 0) return;
  
  size_t i = 0;
  size_t suffix = 16 * (length / 16);

  __m256d alpha_vec = _mm256_set_pd(alpha, alpha, alpha, alpha);
  for ( ; i < suffix; i += 16) {
    _mm256_stream_pd(y + i     , _mm256_sub_pd(_mm256_load_pd(y + i     ), _mm256_mul_pd(_mm256_load_pd(x + i     ), alpha_vec)));
    _mm256_stream_pd(y + i +  4, _mm256_sub_pd(_mm256_load_pd(y + i +  4), _mm256_mul_pd(_mm256_load_pd(x + i +  4), alpha_vec)));
    _mm256_stream_pd(y + i +  8, _mm256_sub_pd(_mm256_load_pd(y + i +  8), _mm256_mul_pd(_mm256_load_pd(x + i +  8), alpha_vec)));
    _mm256_stream_pd(y + i + 12, _mm256_sub_pd(_mm256_load_pd(y + i + 12), _mm256_mul_pd(_mm256_load_pd(x + i + 12), alpha_vec)));
  }
  
  suffix = 4 * (length / 4);
  
  for ( ; i < suffix; i += 4) {
    _mm256_stream_pd(y + i     , _mm256_sub_pd(_mm256_load_pd(y + i     ), _mm256_mul_pd(_mm256_load_pd(x + i     ), alpha_vec)));
  }
  
  for ( ; i < length; ++i)
    y[i] += alpha * x[i];
}

void misc_addScalarToVectorInPlace_avx(double* x, size_t length, double alpha)
{
  if (length == 0) return;
  
  size_t offset = ((uintptr_t) x) % (4 * sizeof(double));
  size_t prefix = offset == 0 ? 0 : (4 * sizeof(double) - offset) / sizeof(double);
  
  if (prefix > length) prefix = length;
  
  size_t i = 0;
  for ( ; i < prefix; ++i)
    x[i] += alpha;
  
  size_t suffix = prefix + 16 * ((length - prefix) / 16);
  
  __m256d alpha_vec = _mm256_set_pd(alpha, alpha, alpha, alpha);
  
  for ( ; i < suffix; i += 16) {
    _mm256_stream_pd(x + i     , _mm256_add_pd(_mm256_load_pd(x + i     ), alpha_vec));
    _mm256_stream_pd(x + i +  4, _mm256_add_pd(_mm256_load_pd(x + i +  4), alpha_vec));
    _mm256_stream_pd(x + i +  8, _mm256_add_pd(_mm256_load_pd(x + i +  8), alpha_vec));
    _mm256_stream_pd(x + i + 12, _mm256_add_pd(_mm256_load_pd(x + i + 12), alpha_vec));
  }
  
  suffix = prefix + 4 * ((length - prefix) / 4);
  for ( ; i < suffix; i += 4) {
    _mm256_stream_pd(x + i     , _mm256_add_pd(_mm256_load_pd(x + i     ), alpha_vec));
  }
  
  for ( ; i < length; ++i)
    x[i] += alpha;
}


void misc_setVectorToConstant_avx(double* x, size_t length, double alpha)
{
  if (length == 0) return;
  
  size_t offset = ((uintptr_t) x) % (4 * sizeof(double));
  size_t prefix = offset == 0 ? 0 : (4 * sizeof(double) - offset) / sizeof(double);
  
  if (prefix > length) prefix = length;
  
  size_t i = 0;
  for ( ; i < prefix; ++i)
    x[i] = alpha;
  
  size_t suffix = prefix + 16 * ((length - prefix) / 16);
  
  __m256d alpha_vec = _mm256_set_pd(alpha, alpha, alpha, alpha);
  
  for ( ; i < suffix; i += 16) {
    _mm256_stream_pd(x + i     , alpha_vec);
    _mm256_stream_pd(x + i +  4, alpha_vec);
    _mm256_stream_pd(x + i +  8, alpha_vec);
    _mm256_stream_pd(x + i + 12, alpha_vec);
  }
  
  suffix = prefix + 4 * ((length - prefix) / 4);
  for ( ; i < suffix; i += 4) {
    _mm256_stream_pd(x + i     , alpha_vec);
  }
  
  for ( ; i < length; ++i)
    x[i] = alpha;
}


// 8 x 8 blocks
static inline void transposeMatrixBlock(const double* restrict x, size_t ldx, double* restrict y, size_t ldy)
{
  // x:  0  8
  //     1  9
  //     2 10
  //     3 11
  //     4 12
  //     5 13
  //     6 14
  //     7 15

  
  __m256d col0, col1, col2, col3;
  __m256d temp0, temp1, temp2, temp3;

  // load upper left quadrant
  col0 = _mm256_loadu_pd(x              ); // 0, 1, 2, 3
  col1 = _mm256_loadu_pd(x +         ldx); // 8, 9, 10, 11
  col2 = _mm256_loadu_pd(x +     2 * ldx); // 16, 17, 18, 19
  col3 = _mm256_loadu_pd(x +     3 * ldx); // 24, 25, 26, 27

  // lane shuffles
  temp0 = _mm256_permute2f128_pd(col0, col2, 0x20); // 0, 1, 16, 17
  temp2 = _mm256_permute2f128_pd(col1, col3, 0x20); // 8, 9, 24, 25
  _mm256_storeu_pd(y              , _mm256_shuffle_pd(temp0, temp2, 0x00)); // 0, 8, 16, 24
  _mm256_storeu_pd(y +         ldy, _mm256_shuffle_pd(temp0, temp2, 0x0F)); // 1, 9, 17, 25

  temp1 = _mm256_permute2f128_pd(col0, col2, 0x31); // 2, 3, 18, 19
  temp3 = _mm256_permute2f128_pd(col1, col3, 0x31); // 10, 11, 26, 77
  _mm256_storeu_pd(y +     2 * ldy, _mm256_shuffle_pd(temp1, temp3, 0x00)); // 2, 10, 18, 26
  _mm256_storeu_pd(y +     3 * ldy, _mm256_shuffle_pd(temp1, temp3, 0x0F)); // 3, 11, 19, 27

  
  col0 = _mm256_loadu_pd(x + 4          ); // 4, 5, 6, 7
  col1 = _mm256_loadu_pd(x + 4 +     ldx); // 12, 13, 14, 15
  col2 = _mm256_loadu_pd(x + 4 + 2 * ldx); // 20, 21, 22, 23
  col3 = _mm256_loadu_pd(x + 4 + 3 * ldx); // 28, 29, 30, 31

  temp0 = _mm256_permute2f128_pd(col0, col2, 0x20);
  temp2 = _mm256_permute2f128_pd(col1, col3, 0x20);
  _mm256_storeu_pd(y +     4 * ldy, _mm256_shuffle_pd(temp0, temp2, 0x00));
  _mm256_storeu_pd(y +     5 * ldy, _mm256_shuffle_pd(temp0, temp2, 0x0F));

  temp1 = _mm256_permute2f128_pd(col0, col2, 0x31);
  temp3 = _mm256_permute2f128_pd(col1, col3, 0x31);
  _mm256_storeu_pd(y +     6 * ldy, _mm256_shuffle_pd(temp1, temp3, 0x00));
  _mm256_storeu_pd(y +     7 * ldy, _mm256_shuffle_pd(temp1, temp3, 0x0F));
  
  
  col0 = _mm256_loadu_pd(x +     4 * ldx);
  col1 = _mm256_loadu_pd(x +     5 * ldx);
  col2 = _mm256_loadu_pd(x +     6 * ldx);
  col3 = _mm256_loadu_pd(x +     7 * ldx);

  temp0 = _mm256_permute2f128_pd(col0, col2, 0x20);
  temp2 = _mm256_permute2f128_pd(col1, col3, 0x20);
  _mm256_storeu_pd(y + 4          , _mm256_shuffle_pd(temp0, temp2, 0x00));
  _mm256_storeu_pd(y + 4 +     ldy, _mm256_shuffle_pd(temp0, temp2, 0x0F));

  temp1 = _mm256_permute2f128_pd(col0, col2, 0x31);
  temp3 = _mm256_permute2f128_pd(col1, col3, 0x31);
  _mm256_storeu_pd(y + 4 + 2 * ldy, _mm256_shuffle_pd(temp1, temp3, 0x00));
  _mm256_storeu_pd(y + 4 + 3 * ldy, _mm256_shuffle_pd(temp1, temp3, 0x0F));


  col0 = _mm256_loadu_pd(x + 4 + 4 * ldx);
  col1 = _mm256_loadu_pd(x + 4 + 5 * ldx);
  col2 = _mm256_loadu_pd(x + 4 + 6 * ldx);
  col3 = _mm256_loadu_pd(x + 4 + 7 * ldx);

  temp0 = _mm256_permute2f128_pd(col0, col2, 0x20);
  temp2 = _mm256_permute2f128_pd(col1, col3, 0x20);
  _mm256_storeu_pd(y + 4 + 4 * ldy, _mm256_shuffle_pd(temp0, temp2, 0x00));
  _mm256_storeu_pd(y + 4 + 5 * ldy, _mm256_shuffle_pd(temp0, temp2, 0x0F));

  temp1 = _mm256_permute2f128_pd(col0, col2, 0x31);
  temp3 = _mm256_permute2f128_pd(col1, col3, 0x31);
  _mm256_storeu_pd(y + 4 + 6 * ldy, _mm256_shuffle_pd(temp1, temp3, 0x00));
  _mm256_storeu_pd(y + 4 + 7 * ldy, _mm256_shuffle_pd(temp1, temp3, 0x0F));
}

void misc_transposeMatrix_avx(const double* restrict x, size_t numRows, size_t numCols, double* restrict y)
{
  if (numRows == 0 || numCols == 0) return;
  
  // We can't really ensure that loads/stores occur on 32 byte boundaries, since
  // any time there is an odd number of rows that completely screws up
  // block transposing. For fun, we start x on a boundary and transpose
  // the first row explicitly, if necessary.
  size_t x_offset = ((uintptr_t) x) % (4 * sizeof(double));
  size_t prefix = x_offset == 0 ? 0 : (4 * sizeof(double) - x_offset) / sizeof(double);
  
  size_t row = 0;
  
  for ( ; row < prefix; ++row) {
    for (size_t col = 0; col < numCols; ++col) {
      y[col + row * numCols] = x[row + col * numRows];
    }
  }
  
  size_t suffix = prefix + 8 * ((numRows - prefix) / 8);
  
  if (suffix > prefix) {
    for ( ; row < suffix; row += 8) {
      size_t col = 0, colEnd = 8 * (numCols / 8);
      for ( ; col < colEnd; col += 8)
        transposeMatrixBlock(x + row + col * numRows, numRows, y + col + row * numCols, numCols);
      for (size_t rowInBlock = row; rowInBlock < row + 8; ++rowInBlock) {
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
