#include "config.h"
#include <misc/linearAlgebra.h>

#include <stdint.h> // uintptr_t

#include <misc/intrinsic.h>

void misc_addVectorsInPlace_sse2(const double* restrict x, size_t length, double* restrict y)
{
  if (length == 0) return;
  
  size_t i = 0;
  size_t lengthMod4 = length % 4;

  for ( ; i < lengthMod4; ++i) y[i] += x[i];
  
  for ( ; i < length; i += 4) {
    y[i    ] += x[i    ];
    y[i + 1] += x[i + 1];
    y[i + 2] += x[i + 2];
    y[i + 3] += x[i + 3];
  }
}

void misc_subtractVectorsInPlace_sse2(const double* restrict x, size_t length, double* restrict y)
{
  if (length == 0) return;
  
  size_t i = 0;
  size_t lengthMod4 = length % 4;

  for ( ; i < lengthMod4; ++i) y[i] -= x[i];
  
  for ( ; i < length; i += 4) {
    y[i    ] -= x[i    ];
    y[i + 1] -= x[i + 1];
    y[i + 2] -= x[i + 2];
    y[i + 3] -= x[i + 3];
  }
}

void misc_addVectorsInPlaceWithMultiplier_sse2(const double* restrict x, size_t length, double alpha, double* restrict y)
{
  if (length == 0) return;
  
  size_t i = 0;
  size_t lengthMod4 = length % 4;

  for ( ; i < lengthMod4; ++i) y[i] += alpha * x[i];
  
  for ( ; i < length; i += 4) {
    y[i    ] += alpha * x[i    ];
    y[i + 1] += alpha * x[i + 1];
    y[i + 2] += alpha * x[i + 2];
    y[i + 3] += alpha * x[i + 3];
  }
  
}

void misc_addScalarToVectorInPlace_sse2(double* x, size_t length, double alpha)
{
  if (length == 0) return;
  
  size_t i = 0;
  size_t lengthMod4 = length % 4;

  for ( ; i < lengthMod4; ++i) x[i] += alpha;
  
  for ( ; i < length; i += 4) {
    x[i    ] += alpha;
    x[i + 1] += alpha;
    x[i + 2] += alpha;
    x[i + 3] += alpha;
  }
}

void misc_setVectorToConstant_sse2(double* x, size_t length, double alpha)
{
  if (length == 0) return;
  
  size_t i = 0;
  size_t lengthMod4 = length % 4;

  for ( ; i < lengthMod4; ++i) x[i] = alpha;
  
  for ( ; i < length; i += 4) {
    x[i    ] = alpha;
    x[i + 1] = alpha;
    x[i + 2] = alpha;
    x[i + 3] = alpha;
  }
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
