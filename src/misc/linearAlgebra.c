#include "config.h"
#include <misc/linearAlgebra.h>

#include <stdbool.h>

void (*misc_addVectorsInPlace)(const double* restrict x, size_t length, double* restrict y) = 0;
void (*misc_subtractVectorsInPlace)(const double* restrict x, size_t length, double* restrict y) = 0;
void (*misc_addVectorsInPlaceWithMultiplier)(const double* restrict x, size_t length, double alpha, double* restrict y) = 0;

void (*misc_addAlignedVectorsInPlace)(const double* restrict x, size_t length, double* restrict y) = 0;
void (*misc_subtractAlignedVectorsInPlace)(const double* restrict x, size_t length, double* restrict y) = 0;

void (*misc_addScalarToVectorInPlace)(double* x, size_t length, double alpha) = 0;
void (*misc_setVectorToConstant)(double* x, size_t length, double alpha) = 0;

void (*misc_transposeMatrix)(const double* restrict x, size_t numRows, size_t numCols, double* restrict xt) = 0;


void misc_addVectors(const double* restrict x, size_t length, const double* restrict y, double* restrict z)
{
  if (length == 0) return;
  
  size_t i = 0;
  size_t lengthMod4 = length % 4;
  
  for ( ; i < lengthMod4; ++i) z[i] = y[i] + x[i];
  
  for ( ; i < length; i += 4) {
    z[i    ] = y[i    ] + x[i    ];
    z[i + 1] = y[i + 1] + x[i + 1];
    z[i + 2] = y[i + 2] + x[i + 2];
    z[i + 3] = y[i + 3] + x[i + 3];
  }
}

void misc_subtractVectors(const double* restrict x, size_t length, const double* restrict y, double* restrict z)
{
  if (length == 0) return;
  
  size_t i = 0;
  size_t lengthMod4 = length % 4;
  
  for ( ; i < lengthMod4; ++i) z[i] = y[i] - x[i];
  
  for ( ; i < length; i += 4) {
    z[i    ] = y[i    ] - x[i    ];
    z[i + 1] = y[i + 1] - x[i + 1];
    z[i + 2] = y[i + 2] - x[i + 2];
    z[i + 3] = y[i + 3] - x[i + 3];
  }
}

void misc_addVectorsWithMultiplier(const double* restrict x, size_t length, double alpha, const double* restrict y, double* restrict z)
{
  if (length == 0 || alpha == 0.0) return;
  
  size_t i = 0;
  size_t lengthMod4 = length % 4;
  
  for ( ; i < lengthMod4; ++i) z[i] = y[i] + alpha * x[i];
  
  for ( ; i < length; i += 4) {
    z[i    ] = y[i    ] + alpha * x[i    ];
    z[i + 1] = y[i + 1] + alpha * x[i + 1];
    z[i + 2] = y[i + 2] + alpha * x[i + 2];
    z[i + 3] = y[i + 3] + alpha * x[i + 3];
  }
}

void misc_addVectorsInPlace_c(const double* restrict x, size_t length, double* restrict y)
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

void misc_subtractVectorsInPlace_c(const double* restrict x, size_t length, double* restrict y)
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

void misc_addVectorsInPlaceWithMultiplier_c(const double* restrict x, size_t length, double alpha, double* restrict y)
{
  if (length == 0 || alpha == 0.0) return;

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


void misc_addScalarToVectorInPlace_c(double* x, size_t length, double alpha)
{
  if (length == 0 || alpha == 0.0) return;
  
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

void misc_setVectorToConstant_c(double* x, size_t length, double alpha)
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

bool misc_vectorIsConstant(const double* x, size_t length)
{
  if (length <= 1) return true;
  
  double x_0 = x[0];
  for (size_t i = 1; i < length; ++i) {
    if (x[i] != x_0) return false;
  }
  
  return true;
}

void misc_setIndexedVectorToConstant(double* restrict x, const size_t* restrict indices, size_t length, double alpha)
{
  if (length == 0) return;
  
  size_t i = 0;
  size_t lengthMod4 = length % 4;

  for ( ; i < lengthMod4; ++i)
    x[indices[i]] = alpha;
  
  for ( ; i < length; i += 4) {
    x[indices[i    ]] = alpha;
    x[indices[i + 1]] = alpha;
    x[indices[i + 2]] = alpha;
    x[indices[i + 3]] = alpha;
  }
}

void misc_scalarMultiplyVectorInPlace(double* x, size_t length, double alpha)
{
  if (length == 0 || alpha == 1.0) return;
  
  size_t i = 0;
  size_t lengthMod4 = length % 4;
    
  for ( ; i < lengthMod4; ++i)
    x[i] *= alpha;
  
  for ( ; i < length; i += 4) {
    x[i    ] *= alpha;
    x[i + 1] *= alpha;
    x[i + 2] *= alpha;
    x[i + 3] *= alpha;
  }
}

void misc_scalarMultiplyVector(const double* restrict x, size_t length, double alpha, double* restrict y)
{
  if (length == 0) return;
  
  size_t i = 0;
  size_t lengthMod4 = length % 4;
  
  for ( ; i < lengthMod4; ++i)
    y[i] = alpha * x[i];
  
  for ( ; i < length; i += 4) {
    y[i    ] = alpha * x[i    ];
    y[i + 1] = alpha * x[i + 1];
    y[i + 2] = alpha * x[i + 2];
    y[i + 3] = alpha * x[i + 3];
  }
}

void misc_hadamardMultiplyVectors(const double* restrict x, size_t length, const double* restrict y, double* restrict z)
{
  if (length == 0) return;
  
  size_t i = 0;
  size_t lengthMod4 = length % 4;
  
  for ( ; i < lengthMod4; ++i)
    z[i] = y[i] * x[i];
  
  for ( ; i < length; i += 4) {
    z[i    ] = y[i    ] * x[i    ];
    z[i + 1] = y[i + 1] * x[i + 1];
    z[i + 2] = y[i + 2] * x[i + 2];
    z[i + 3] = y[i + 3] * x[i + 3];
  }
}

void misc_hadamardMultiplyVectorsInPlace(double* restrict x, size_t length, const double* restrict y)
{
  if (length == 0) return;
  
  size_t i = 0;
  size_t lengthMod4 = length % 4;
  
  for ( ; i < lengthMod4; ++i)
    x[i] *= y[i];
  
  for ( ; i < length; i += 4) {
    x[i    ] *= y[i    ];
    x[i + 1] *= y[i + 1];
    x[i + 2] *= y[i + 2];
    x[i + 3] *= y[i + 3];
  }
}

double misc_sumVectorElements(const double* x, size_t length)
{
  if (length == 0) return 0.0;
  
  double result = 0.0;

  size_t i = 0;
  size_t lengthMod4 = length % 4;
  
  for ( ; i < lengthMod4; ++i)
    result += x[i];
  
  for ( ; i < length; i += 4) {
    result += x[i] + x[i + 1] + x[i + 2] + x[i + 3];
  }
  
  return result;
}

double misc_sumIndexedVectorElements(const double* x, const size_t* indices, size_t length)
{
  if (length == 0) return 0.0;
  
  double result = 0.0;

  size_t i = 0;
  size_t lengthMod4 = length % 4;
  
  for ( ; i < lengthMod4; ++i)
    result += x[indices[i]];

  for ( ; i < length; i += 4) {
    result += x[indices[i]] + x[indices[i + 1]] + x[indices[i + 2]] + x[indices[i + 3]];
  }
  
  return result;
}

void misc_transposeMatrix_c(const double* restrict x, size_t numRows, size_t numCols, double* restrict xt)
{
  if (numRows == 0 || numCols == 0) return;

  size_t colOffset = 0;
  for (size_t col = 0; col < numCols; ++col) {
    for (size_t row = 0; row < numRows; ++row) {
      xt[row * numCols + col] = x[row + colOffset];
    }
    colOffset += numRows;
  }
}

void misc_multiplyMatrixIntoVector(const double* restrict matrix, size_t numRows, size_t numCols, int useTranspose,
                                   const double* restrict vector, double* restrict result)
{
  if (!useTranspose) {
    for (size_t row = 0; row < numRows; ++row) {
      *result = 0.0;
      for (size_t col = 0; col < numCols; ++col) {
        *result += *matrix * *vector;
        matrix += numRows;
        ++vector;
      }
      ++result;
      matrix -= numRows * numCols - 1;
      vector -= numCols;
    }
  } else {
    for (size_t col = 0; col < numCols; ++col) {
      *result = 0.0;
      for (size_t row = 0; row < numRows; ++row) {
        *result += *matrix * *vector;
        ++matrix;
        ++vector;
      }
      ++result;
      vector -= numRows;
    }
  }
}


