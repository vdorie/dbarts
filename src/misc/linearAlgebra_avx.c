#include "config.h"
#include <misc/linearAlgebra.h>

#include <stdint.h> // uintptr_t

#include <misc/alloca.h>
#include <misc/intrinsic.h>

void misc_addVectorsInPlace_avx(const double* restrict x, size_t length, double* restrict y)
{
  if (length == 0) return;
  
  size_t i = 0;
  size_t lengthMod8 = length % 8;

  for ( ; i < lengthMod8; ++i) y[i] += x[i];
  
  for ( ; i < length; i += 8) {
    y[i    ] += x[i    ];
    y[i + 1] += x[i + 1];
    y[i + 2] += x[i + 2];
    y[i + 3] += x[i + 3];
    y[i + 4] += x[i + 4];
    y[i + 5] += x[i + 5];
    y[i + 6] += x[i + 6];
    y[i + 7] += x[i + 7];
  }
}
  

void misc_subtractVectorsInPlace_avx(const double* restrict x, size_t length, double* restrict y)
{
  if (length == 0) return;
  
  size_t i = 0;
  size_t lengthMod8 = length % 8;

  for ( ; i < lengthMod8; ++i) y[i] -= x[i];
  
  for ( ; i < length; i += 8) {
    y[i    ] -= x[i    ];
    y[i + 1] -= x[i + 1];
    y[i + 2] -= x[i + 2];
    y[i + 3] -= x[i + 3];
    y[i + 4] -= x[i + 4];
    y[i + 5] -= x[i + 5];
    y[i + 6] -= x[i + 6];
    y[i + 7] -= x[i + 7];
  }
}
  

void misc_addVectorsInPlaceWithMultiplier_avx(const double* restrict x, size_t length, double alpha, double* restrict y)
{
  if (length == 0) return;
  
  size_t i = 0;
  size_t lengthMod8 = length % 8;

  for ( ; i < lengthMod8; ++i) y[i] += alpha * x[i];
  
  for ( ; i < length; i += 8) {
    y[i    ] += alpha * x[i    ];
    y[i + 1] += alpha * x[i + 1];
    y[i + 2] += alpha * x[i + 2];
    y[i + 3] += alpha * x[i + 3];
    y[i + 4] += alpha * x[i + 4];
    y[i + 5] += alpha * x[i + 5];
    y[i + 6] += alpha * x[i + 6];
    y[i + 7] += alpha * x[i + 7];
  }
}
  

void misc_addScalarToVectorInPlace_avx(double* x, size_t length, double alpha)
{
  if (length == 0) return;
  
  size_t i = 0;
  size_t lengthMod8 = length % 8;

  for ( ; i < lengthMod8; ++i) x[i] += alpha;
  
  for ( ; i < length; i += 8) {
    x[i    ] += alpha;
    x[i + 1] += alpha;
    x[i + 2] += alpha;
    x[i + 3] += alpha;
    x[i + 4] += alpha;
    x[i + 5] += alpha;
    x[i + 6] += alpha;
    x[i + 7] += alpha;
  }
  
}


void misc_setVectorToConstant_avx(double* x, size_t length, double alpha)
{
  if (length == 0) return;
  
  size_t i = 0;
  size_t lengthMod8 = length % 8;

  for ( ; i < lengthMod8; ++i) x[i] = alpha;
  
  for ( ; i < length; i += 8) {
    x[i    ] = alpha;
    x[i + 1] = alpha;
    x[i + 2] = alpha;
    x[i + 3] = alpha;
    x[i + 4] = alpha;
    x[i + 5] = alpha;
    x[i + 6] = alpha;
    x[i + 7] = alpha;
  }
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
  prefix = prefix > numRows ? numRows : prefix;
  
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
