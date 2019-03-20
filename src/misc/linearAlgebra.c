#include "config.h"
#include <misc/linearAlgebra.h>

#include <stdbool.h>

void (*misc_addVectors)(const double* restrict x, size_t length, double alpha, const double* restrict y, double* restrict z) = 0;
void (*misc_setVectorToConstant)(double* x, size_t length, double alpha) = 0;

void misc_addVectors_c(const double* restrict x, size_t length, double alpha, const double* restrict y, double* restrict z)
{
  if (length == 0) return;
  
  size_t i = 0;
  size_t lengthMod5 = length % 5;
  
  if (lengthMod5 != 0) {
    for ( ; i < lengthMod5; ++i) z[i] = y[i] + alpha * x[i];
    if (length < 5) return;
  }
  
  for ( ; i < length; i += 5) {
    z[i] = y[i] + alpha * x[i];
    z[i + 1] = y[i + 1] + alpha * x[i + 1];
    z[i + 2] = y[i + 2] + alpha * x[i + 2];
    z[i + 3] = y[i + 3] + alpha * x[i + 3];
    z[i + 4] = y[i + 4] + alpha * x[i + 4];
  }
}

void misc_setVectorToConstant_c(double* x, size_t length, double alpha)
{
  if (length == 0) return;
  
  size_t i = 0;
  size_t lengthMod5 = length % 5;
  
  if (lengthMod5 != 0) {
    for ( ; i < lengthMod5; ++i) x[i] = alpha;
    if (length < 5) return;
  }
  
  for ( ; i < length; i += 5) {
    x[i]     = alpha;
    x[i + 1] = alpha;
    x[i + 2] = alpha;
    x[i + 3] = alpha;
    x[i + 4] = alpha;
  }
}

bool misc_vectorIsConstant(const double* x, size_t length)
{
  if (length <= 1) return true;
  
  for (size_t i = 1; i < length; ++i) {
    if (x[i] != x[i - 1]) return false;
  }
  
  return true;
}

void misc_setIndexedVectorToConstant(double* restrict x, const size_t* restrict indices, size_t length, double alpha)
{
  if (length == 0) return;
  
  size_t lengthMod5 = length % 5;
  
  if (lengthMod5 != 0) {
    for (size_t i = 0; i < lengthMod5; ++i) x[indices[i]] = alpha;
    if (length < 5) return;
  }
  
  for (size_t i = lengthMod5; i < length; i += 5) {
    x[indices[i]]     = alpha;
    x[indices[i + 1]] = alpha;
    x[indices[i + 2]] = alpha;
    x[indices[i + 3]] = alpha;
    x[indices[i + 4]] = alpha;
  }
}

void misc_scalarMultiplyVectorInPlace(double* x, size_t length, double alpha)
{
  if (length == 0) return;
  
  size_t lengthMod5 = length % 5;
    
  if (lengthMod5 != 0) {
    for (size_t i = 0; i < lengthMod5; ++i) x[i] *= alpha;
    if (length < 5) return;
  }
  
  for (size_t i = lengthMod5; i < length; i += 5) {
    x[i] *= alpha;
    x[i + 1] *= alpha;
    x[i + 2] *= alpha;
    x[i + 3] *= alpha;
    x[i + 4] *= alpha;
  }
}

void misc_addScalarToVectorInPlace(double* x, size_t length, double alpha)
{
  if (length == 0) return;
  
  size_t lengthMod5 = length % 5;
  
  size_t i = 0;
  for ( /* */ ; i < lengthMod5; ++i) x[i] += alpha;
  
  for ( /* */ ; i < length; i += 5) {
    x[i]     += alpha;
    x[i + 1] += alpha;
    x[i + 2] += alpha;
    x[i + 3] += alpha;
    x[i + 4] += alpha;
  }
}

void misc_scalarMultiplyVector(const double* restrict x, size_t length, double alpha, double* restrict y)
{
  if (length == 0) return;
  
  size_t lengthMod5 = length % 5;
  
  if (lengthMod5 != 0) {
    for (size_t i = 0; i < lengthMod5; ++i) y[i] = alpha * x[i];
    if (length < 5) return;
  }
  
  for (size_t i = lengthMod5; i < length; i += 5) {
    y[i]     = alpha * x[i];
    y[i + 1] = alpha * x[i + 1];
    y[i + 2] = alpha * x[i + 2];
    y[i + 3] = alpha * x[i + 3];
    y[i + 4] = alpha * x[i + 4];
  }
}

void misc_hadamardMultiplyVectors(const double* restrict x, size_t length, const double* restrict y, double* restrict z)
{
  if (length == 0) return;
  
  size_t lengthMod5 = length % 5;
  
  if (lengthMod5 != 0) {
    for (size_t i = 0; i < lengthMod5; ++i) z[i] = y[i] * x[i];
    if (length < 5) return;
  }
  
  for (size_t i = lengthMod5; i < length; i += 5) {
    z[i] = y[i] * x[i];
    z[i + 1] = y[i + 1] * x[i + 1];
    z[i + 2] = y[i + 2] * x[i + 2];
    z[i + 3] = y[i + 3] * x[i + 3];
    z[i + 4] = y[i + 4] * x[i + 4];
  }
}

void misc_hadamardMultiplyVectorsInPlace(double* restrict x, size_t length, const double* restrict y)
{
  if (length == 0) return;
  
  size_t lengthMod5 = length % 5;
  
  if (lengthMod5 != 0) {
    for (size_t i = 0; i < lengthMod5; ++i) x[i] *= y[i];
    if (length < 5) return;
  }
  
  for (size_t i = lengthMod5; i < length; i += 5) {
    x[i] *= y[i];
    x[i + 1] *= y[i + 1];
    x[i + 2] *= y[i + 2];
    x[i + 3] *= y[i + 3];
    x[i + 4] *= y[i + 4];
  }
}

double misc_sumVectorElements(const double* x, size_t length)
{
  if (length == 0) return 0.0;
  
  double result = 0.0;
  size_t lengthMod5 = length % 5;
  
  if (lengthMod5 != 0) {
    for (size_t i = 0; i < lengthMod5; ++i) result += x[i];
    if (length < 5) return result;
  }
  
  for (size_t i = lengthMod5; i < length; i += 5) {
    result += x[i] + x[i + 1] + x[i + 2] + x[i + 3] + x[i + 4];
  }
  
  return result;
}

double misc_sumIndexedVectorElements(const double* x, const size_t* indices, size_t length)
{
  if (length == 0) return 0.0;
  
  double result = 0.0;
  size_t lengthMod5 = length % 5;
  
  size_t i;
  for (i = 0; i < lengthMod5; ++i) result += x[indices[i]];
  for (/* */; i < length; i += 5) {
    result += x[indices[i]] + x[indices[i + 1]] + x[indices[i + 2]] + x[indices[i + 3]] + x[indices[i + 4]];
  }
  
  return result;
}

void misc_transposeMatrix(const double* x, size_t numRows, size_t numCols, double* xt)
{
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

