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

void misc_addVectorsInPlace_sse2(const double* restrict x, size_t length, double alpha, double* restrict y)
{
  if (length == 0) return;

  // _mm_load_pd requires that the point address be aligned on a 16 byte boundary, or
  // ((uintptr_t) y) % (2 * sizeof(double)) == 0
  // If both y and x have the same alignment, we may need to start indexing them at
  // 1 in order for them to be aligned, but after we do so we can use
  // __mm_load_pd for both. If they don't share an alignment, we have to use 
  // _mm_loadu_pd for at least x
  size_t y_offset = ((uintptr_t) y) % (2 * sizeof(double));
  size_t x_offset = ((uintptr_t) x) % (2 * sizeof(double));
  size_t prefix = y_offset == 0 ? 0 : (2 * sizeof(double) - y_offset) / sizeof(double);
  
  if (prefix > length) prefix = length;
  
  size_t i = 0;
  for ( ; i < prefix; ++i)
    y[i] += + alpha * x[i];
  
  size_t suffix = prefix + 6 * ((length - prefix) / 6);
  
  if (suffix > prefix) {
    __m128d alpha_vec = _mm_load1_pd(&alpha);
    if (y_offset == x_offset) {
      for ( ; i < suffix; i += 6) {
        _mm_stream_pd(y + i    , _mm_add_pd(_mm_load_pd(y + i    ), _mm_mul_pd(_mm_load_pd(x + i    ), alpha_vec)));
        _mm_stream_pd(y + i + 2, _mm_add_pd(_mm_load_pd(y + i + 2), _mm_mul_pd(_mm_load_pd(x + i + 2), alpha_vec)));
        _mm_stream_pd(y + i + 4, _mm_add_pd(_mm_load_pd(y + i + 4), _mm_mul_pd(_mm_load_pd(x + i + 4), alpha_vec)));
      }
    } else {
      for ( ; i < suffix; i += 6) {
        _mm_stream_pd(y + i    , _mm_add_pd(_mm_load_pd(y + i    ), _mm_mul_pd(_mm_loadu_pd(x + i    ), alpha_vec)));
        _mm_stream_pd(y + i + 2, _mm_add_pd(_mm_load_pd(y + i + 2), _mm_mul_pd(_mm_loadu_pd(x + i + 2), alpha_vec)));
        _mm_stream_pd(y + i + 4, _mm_add_pd(_mm_load_pd(y + i + 4), _mm_mul_pd(_mm_loadu_pd(x + i + 4), alpha_vec)));
      }
    }
  }
  
  for ( ; i < length; ++i)
    y[i] += alpha * x[i];
}

extern void misc_addScalarToVectorInPlace_c(double* x, size_t length, double alpha);

void misc_addScalarToVectorInPlace_sse2(double* x, size_t length, double alpha)
{
  if (length == 0) return;
  if (((uintptr_t) x) % sizeof(double) != 0) {
    misc_addScalarToVectorInPlace_sse2(x, length, alpha);
    return;
  }
  
  size_t offset = ((uintptr_t) x) % (2 * sizeof(double));
  size_t prefix = offset == 0 ? 0 : (2 * sizeof(double) - offset) / sizeof(double);
  
  if (prefix > length) prefix = length;
  
  size_t i = 0;
  for ( ; i < prefix; ++i)
    x[i] += alpha;
  
  size_t suffix = prefix + 6 * ((length - prefix) / 6);
  
  if (suffix > prefix) {
    __m128d alpha_vec = _mm_load1_pd(&alpha);
    for ( ; i < suffix; i += 6) {
      _mm_stream_pd(x + i    , _mm_add_pd(_mm_load_pd(x + i    ), alpha_vec));
      _mm_stream_pd(x + i + 2, _mm_add_pd(_mm_load_pd(x + i + 2), alpha_vec));
      _mm_stream_pd(x + i + 4, _mm_add_pd(_mm_load_pd(x + i + 4), alpha_vec));
    }
  }
  
  for ( ; i < length; ++i)
    x[i] += alpha;
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


static inline void tranposeMatrixBlock(const double* restrict x, size_t ldx, double* restrict y, size_t ldy)
{
  // x:  0  4  8 12
  //     1  5  9 13
  //     2  6 10 14
  //     3  7 11 15
  // y:  0  1  2  3
  //     4  5  6  7
  //     8  9 10 11
  //    12 13 14 15
  
  __m128d temp0, temp1;

  temp0 = _mm_loadu_pd(x              ); // 0, 1
  temp1 = _mm_loadu_pd(x +         ldx); // 4, 5
  _mm_storeu_pd(y              , _mm_shuffle_pd(temp0, temp1, 0x0)); // 0, 4
  _mm_storeu_pd(y +         ldy, _mm_shuffle_pd(temp0, temp1, 0x3)); // 1 5
  
  temp0 = _mm_loadu_pd(x + 2          ); // 2, 3
  temp1 = _mm_loadu_pd(x + 2 +     ldx); // 6, 7
  _mm_storeu_pd(y +     2 * ldy, _mm_shuffle_pd(temp0, temp1, 0x0));
  _mm_storeu_pd(y +     3 * ldy, _mm_shuffle_pd(temp0, temp1, 0x3));
  
  temp0 = _mm_loadu_pd(x +     2 * ldx); // 8, 9
  temp1 = _mm_loadu_pd(x +     3 * ldx); // 12, 13
  _mm_storeu_pd(y + 2          , _mm_shuffle_pd(temp0, temp1, 0x0));
  _mm_storeu_pd(y + 2 +     ldy, _mm_shuffle_pd(temp0, temp1, 0x3));
  
  temp0 = _mm_loadu_pd(x + 2 + 2 * ldx); // 10, 11
  temp1 = _mm_loadu_pd(x + 2 + 3 * ldx); // 14, 15
  _mm_storeu_pd(y + 2 + 2 * ldy, _mm_shuffle_pd(temp0, temp1, 0x0));
  _mm_storeu_pd(y + 2 + 3 * ldy, _mm_shuffle_pd(temp0, temp1, 0x3));
}

void misc_transposeMatrix_sse2(const double* restrict x, size_t numRows, size_t numCols, double* restrict y)
{
  if (numRows == 0 || numCols == 0) return;

  // We can't really ensure that loads/stores occur on 16 byte boundaries, since
  // any time there is an odd number of rows that completely screws up
  // block transposing. For fun, we start x on a boundary and transpose
  // the first row explicitly, if necessary.
  size_t x_offset = ((uintptr_t) x) % (2 * sizeof(double));
  size_t prefix = x_offset == 0 ? 0 : (2 * sizeof(double) - x_offset) / sizeof(double);

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
      for ( ; col < numCols; ++col) {
        y[col + row * numCols] = x[row + col * numRows];
      }
    }
  }

  for ( ; row < numRows; ++row) {
    for (size_t col = 0; col < numCols; ++col) {
      y[col + row * numCols] = x[row + col * numRows];
    }
  }
}
