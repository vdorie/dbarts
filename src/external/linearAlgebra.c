#include "config.h"
#include <external/linearAlgebra.h>

#include <errno.h>
#include <math.h> // isinf, isnan
#include <stdbool.h>
#include <stdlib.h> // malloc
#include <string.h> // memcpy

#include <misc/alloca.h>

#include <Rversion.h>

#if R_VERSION >= R_Version(3, 6, 2)
#define USE_FC_LEN_T
#endif

#if R_VERSION <= R_Version(3, 3, 1)
#  define NO_C_HEADERS
#endif

#define R_NO_REMAP
#include <R_ext/Lapack.h> // dlassq, dpotrf, dtrtrs
#include <R_ext/Applic.h> // dqrls

#undef R_NO_REMAP
#undef NO_C_HEADERS
#undef USE_FC_LEN_T

#ifndef FCONE
# define FCONE
#endif

static const int increment = 1;

void ext_addVectorsInPlace(const double* restrict x, size_t u_length, double alpha, double* restrict y)
{
  int length = (int) u_length;
  
  F77_CALL(daxpy)(&length, &alpha, x, &increment, y, &increment);
}

// b = Ax
void ext_leftMultiplyMatrixAndVector(const double* A, size_t n, size_t p, const double* x, double* b)
{
  char transpose = 'N';
  int i_n = (int) n;
  int i_p = (int) p;
  
  double d_one  = 1.0;
  double d_zero = 0.0;
  
  F77_CALL(dgemv)(&transpose, // trans
                  &i_n, &i_p, // m, n
                  &d_one, // alpha
                  A, &i_n, //a, lda
                  x, &increment, // x, incx
                  &d_zero, //beta
                  b, &increment FCONE); // y, incy, char length
}

double ext_sumSquaresOfVectorElements(const double* x, size_t u_length)
{
  int length = (int) u_length;
  
  double sum = 0.0;
  double scale = 1.0;
  
  F77_CALL(dlassq)(&length, x, &increment, &scale, &sum);
  
  return scale * scale * sum;
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
  if (residuals == NULL) residuals = misc_stackAllocate(n, double);
  if (residuals == NULL) {
    if (message != NULL) message[0] = "unable to allocate memory for least-squares fit";
    return - 1;
  }
  double effects[n];
  double work[2 * p];
  double qrAux[p];
  int rank;
  F77_CALL(dqrls)(x, &i_n, &i_p, y, (int*) &increment, &tolerance, b,
                  residuals, effects, &rank, pivot, qrAux, work);
  
  bool allKosher = true;
  for (size_t i = 0; i < p; ++i) {
    if (isnan(b[i]) || isinf(b[i])) {
      allKosher = false;
      break;
    }
  }
  
  if (_residuals == NULL) { misc_stackFree(residuals); }
  
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
  
  F77_CALL(dpotrf)(&useLowerTriangle, &i_dim, result, &i_dim, &lapackResult FCONE);
  
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
  
  F77_CALL(dpotrf)(&useLowerTriangle, &i_dim, x, &i_dim, &lapackResult FCONE);
  
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
                  result, &outputDimension FCONE FCONE);
  
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
                   &i_lhsDim, &i_rhsDim, lhs, &i_lhsDim, rhs, &i_lhsDim, &lapackResult FCONE FCONE FCONE);
  
  if (lapackResult < 0) return EINVAL;
  if (lapackResult > 0) return EDOM;
  
  return 0;
}

