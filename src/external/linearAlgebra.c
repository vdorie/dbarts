#include "config.h"
#include <external/linearAlgebra.h>

#include <stdbool.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>

#include <R.h>
#include <Rdefines.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <R_ext/Applic.h> // for dqrls

#include <external/alloca.h>

static const int increment = 1;

void ext_addVectorsInPlace(const double* restrict x, size_t u_length, double alpha, double* restrict y)
{
  int length = (int) u_length;
  
  F77_NAME(daxpy)(&length, &alpha, x, &increment, y, &increment);
}

void ext_addVectors(const double* restrict x, size_t length, double alpha, const double* restrict y, double* restrict z)
{
  if (length == 0) return;
  
  size_t lengthMod5 = length % 5;
  
  if (lengthMod5 != 0) {
    for (size_t i = 0; i < lengthMod5; ++i) z[i] = y[i] + alpha * x[i];
    if (length < 5) return;
  }
  
  for (size_t i = lengthMod5; i < length; i += 5) {
    z[i] = y[i] + alpha * x[i];
    z[i + 1] = y[i + 1] + alpha * x[i + 1];
    z[i + 2] = y[i + 2] + alpha * x[i + 2];
    z[i + 3] = y[i + 3] + alpha * x[i + 3];
    z[i + 4] = y[i + 4] + alpha * x[i + 4];
  }
}

void ext_setVectorToConstant(double* x, size_t length, double alpha)
{
  if (length == 0) return;
  
  size_t lengthMod5 = length % 5;
  
  if (lengthMod5 != 0) {
    for (size_t i = 0; i < lengthMod5; ++i) x[i] = alpha;
    if (length < 5) return;
  }
  
  for (size_t i = lengthMod5; i < length; i += 5) {
    x[i]     = alpha;
    x[i + 1] = alpha;
    x[i + 2] = alpha;
    x[i + 3] = alpha;
    x[i + 4] = alpha;
  }
}

void ext_setIndexedVectorToConstant(double* restrict x, const size_t* restrict indices, size_t length, double alpha)
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

void ext_scalarMultiplyVectorInPlace(double* x, size_t length, double alpha)
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

void ext_addScalarToVectorInPlace(double* x, size_t length, double alpha)
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


void ext_scalarMultiplyVector(const double* restrict x, size_t length, double alpha, double* restrict y)
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

void ext_hadamardMultiplyVectors(const double* restrict x, size_t length, const double* restrict y, double* restrict z)
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

void ext_hadamardMultiplyVectorsInPlace(double* restrict x, size_t length, const double* restrict y)
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

// b = Ax
void ext_leftMultiplyMatrixAndVector(const double* A, size_t n, size_t p, const double* x, double* b)
{
  char transpose = 'N';
  int i_n = (int) n;
  int i_p = (int) p;
  
  double d_one  = 1.0;
  double d_zero = 0.0;
  
  F77_CALL(dgemv)(&transpose,
                  &i_n, &i_p, &d_one,
                  A, &i_n,
                  x,
                  &increment, &d_zero,
                  b,
                  &increment);
}

double ext_sumVectorElements(const double* x, size_t length)
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

double ext_sumIndexedVectorElements(const double* x, const size_t* indices, size_t length)
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

double ext_sumSquaresOfVectorElements(const double* x, size_t u_length)
{
  int length = (int) u_length;
  
  double sum = 0.0;
  double scale = 1.0;
  
  F77_NAME(dlassq)(&length, x, &increment, &scale, &sum);
  
  return scale * scale * sum;
}

void ext_transposeMatrix(const double* x, size_t numRows, size_t numCols, double* xt)
{
  size_t colOffset = 0;
  for (size_t col = 0; col < numCols; ++col) {
    for (size_t row = 0; row < numRows; ++row) {
      xt[row * numCols + col] = x[row + colOffset];
    }
    colOffset += numRows;
  }
}

// LS solution to Xb = y, or b = (X'X)^-1 X'y
// destroys x and y
// returns -1 for an error, 0 for a warning, and otherwise
// the rank of the solution
int ext_findLeastSquaresFit(const double* y, size_t n, const double* x, size_t p, double* b, double tolerance, double* residuals, char** message)
{
  double* _x = malloc(n * p * sizeof(double));
  double _y[n];
  if (_x == NULL) {
    if (message != NULL) message[0] = "unable to allocate memory for least-squares fit";
    return -1;
  }
  memcpy(_x, x, n * p * sizeof(double));
  memcpy(_y, y, n * sizeof(double));
  
  int result = ext_findLeastSquaresFitInPlace(_y, n, _x, p, b, tolerance, residuals, message);
  
  free(_x);
  
  return result;
}  


int ext_findLeastSquaresFitInPlace(double* y, size_t n, double* x, size_t p, double* b, double tolerance, double* _residuals, char** message) { 
  int i_n = (int) n; int i_p = (int) p;
  int pivot[p];
  for (size_t i = 0; i < p; ++i) pivot[i] = (int) i;
  
  double* residuals = _residuals;
  if (residuals == NULL) residuals = ext_stackAllocate(n, double);
  if (residuals == NULL) {
    if (message != NULL) message[0] = "unable to allocate memory for least-squares fit";
    return - 1;
  }
  double effects[n];
  double work[2 * p];
  double qrAux[p];
  int rank;
  F77_NAME(dqrls)(x, &i_n, &i_p, y, (int*) &increment, &tolerance, b,
                  residuals, effects, &rank, pivot, qrAux, work);
  
  bool allKosher = true;
  for (size_t i = 0; i < p; ++i) {
    if (isnan(b[i]) || isinf(b[i])) {
      allKosher = false;
      break;
    }
  }
  
  if (_residuals == NULL) { ext_stackFree(residuals); }
  
  if (!allKosher) {
    message[0] = "non-finite solution to least squares fit";
    return 0;
  }
  
  double temp[p];
  for (size_t i = 0; i < p; ++i) temp[i] = b[i];
  for (size_t i = 0; i < p; ++i) b[pivot[i]] = temp[i];
  
  return (int) rank;
}

int ext_getSymmetricPositiveDefiniteTriangularFactorization(const double* restrict x,
                                                            size_t dim,
                                                            ext_triangleType triangleType,
                                                            double* restrict result)
{
  if (triangleType == EXT_TRIANGLE_TYPE_BOTH) return EINVAL;
  
  char useLowerTriangle = (triangleType == EXT_TRIANGLE_TYPE_UPPER ? 'U' : 'L');
  
  memcpy(result, x, dim * dim * sizeof(double));
  
  int i_dim = (int) dim;
  
  int lapackResult;
  
  F77_CALL(dpotrf)(&useLowerTriangle, &i_dim, result, &i_dim, &lapackResult);
  
  if (lapackResult < 0) return EINVAL;
  if (lapackResult > 0) return EDOM;
  
  return 0;
}

int ext_getSymmetricPositiveDefiniteTriangularFactorizationInPlace(double* x, size_t dim, ext_triangleType triangleType)
{
  if (triangleType == EXT_TRIANGLE_TYPE_BOTH) return EINVAL;
  
  char useLowerTriangle = (triangleType == EXT_TRIANGLE_TYPE_UPPER ? 'U' : 'L');
  
  int i_dim = (int) dim;
  
  int lapackResult;
  
  F77_CALL(dpotrf)(&useLowerTriangle, &i_dim, x, &i_dim, &lapackResult);
  
  if (lapackResult < 0) return EINVAL;
  if (lapackResult > 0) return EDOM;
  
  return 0;
}


void ext_getSingleMatrixCrossproduct(const double* restrict x, size_t numRows, size_t numCols,
                                     double* restrict result, int useTranspose, ext_triangleType triangleType)
{
  // for us, transpose is AA', for BLAS it is the reverse (A'A)
  char shouldTransposeMatrix = (useTranspose ? 'N' : 'T');
  char useLowerTriangle = (triangleType == EXT_TRIANGLE_TYPE_UPPER ? 'U' : 'L');
  
  int outputDimension  = (int) (useTranspose ? numRows : numCols);
  int inputDimension   = (int) (useTranspose ? numCols : numRows); // if user specifies transpose, he/she wants AA'
  double d_one  = 1.0;
  double d_zero = 0.0;
  
  int i_numRows = (int) numRows;
  
  F77_CALL(dsyrk)(&useLowerTriangle, &shouldTransposeMatrix,
                  &outputDimension,
                  &inputDimension, &d_one, x,
                  &i_numRows, &d_zero,
                  result, &outputDimension);
  
  if (triangleType != EXT_TRIANGLE_TYPE_BOTH) return;
  
  // copy in rest of product
  size_t u_outputDimension = (size_t) outputDimension;
  for (size_t col = 1; col < u_outputDimension; ++col) {
    for (size_t row = 0; row < col; ++row) {
      result[row + u_outputDimension * col] = result[col + u_outputDimension * row];
    }
  }
}

// solves x := A x = y, where A is triangular
int ext_solveTriangularSystemInPlace(const double* restrict lhs, size_t lhsDim, int useTranspose, ext_triangleType triangleType,
                                     double* restrict rhs, size_t rhsDim)
{
  char useLowerTriangle = (triangleType == EXT_TRIANGLE_TYPE_UPPER ? 'U' : 'L');
  char shouldTransposeMatrix = (useTranspose ? 'T' : 'N');
  char isUnitDiagonal = 'N';
  
  int i_lhsDim = (int) lhsDim;
  int i_rhsDim = (int) rhsDim;
  
  int lapackResult;
  
  F77_CALL(dtrtrs)(&useLowerTriangle, &shouldTransposeMatrix, &isUnitDiagonal,
                   &i_lhsDim, &i_rhsDim, lhs, &i_lhsDim, rhs, &i_lhsDim, &lapackResult);
  
  if (lapackResult < 0) return EINVAL;
  if (lapackResult > 0) return EDOM;
  
  return 0;
}

void ext_multiplyMatrixIntoVector(const double* restrict matrix, size_t numRows, size_t numCols, int useTranspose,
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
