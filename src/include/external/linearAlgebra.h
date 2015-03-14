#ifndef EXTERNAL_LINEAR_ALGEBRA_H
#define EXTERNAL_LINEAR_ALGEBRA_H

#include <stdbool.h>
#include "stddef.h"
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif
  
  typedef enum {
    EXT_TRIANGLE_TYPE_BOTH,
    EXT_TRIANGLE_TYPE_UPPER,
    EXT_TRIANGLE_TYPE_LOWER
  } ext_triangleType;
  
  // y := alpha * x + y
  void ext_addVectorsInPlace(const double* restrict x, ext_size_t length, double alpha, double* restrict y);
  // z = alpha * x + y; z must be distinct from x and y
  void ext_addVectors(const double* restrict x, ext_size_t length, double alpha, const double* restrict y, double* restrict z);
  
  // x: = alpha
  void ext_setVectorToConstant(double* x, ext_size_t length, double alpha);
  bool ext_vectorIsConstant(const double* d, ext_size_t length);
  void ext_setIndexedVectorToConstant(double* restrict x, const ext_size_t* restrict indices, ext_size_t length, double alpha);
  
  // x := alpha * x
  void ext_scalarMultiplyVectorInPlace(double* x, ext_size_t length, double alpha);
  // y = alpha * x; y must be distinct from x
  void ext_scalarMultiplyVector(const double* restrict x, ext_size_t length, double alpha, double* restrict y);
  
  // x := x + alpha
  void ext_addScalarToVectorInPlace(double* x, ext_size_t length, double alpha);
  
  // z = x .* y; z must be distinct from x and y
  void ext_hadamardMultiplyVectors(const double* restrict x, ext_size_t length, const double* restrict y, double* restrict z);
  // x := x. * y; x and y should be distinct
  void ext_hadamardMultiplyVectorsInPlace(double* restrict x, ext_size_t length, const double* restrict y);
  
  // b = Ax
  void ext_leftMultiplyMatrixAndVector(const double* A, ext_size_t n, ext_size_t p, const double* x, double* b);
  
  double ext_sumVectorElements(const double* x, ext_size_t length);
  double ext_sumIndexedVectorElements(const double* x, const ext_size_t* indices, ext_size_t length);
  double ext_sumSquaresOfVectorElements(const double* x, ext_size_t length);
  
  void ext_transposeMatrix(const double* x, ext_size_t numRows, ext_size_t numCols, double* xt);
  
  // least squares solution to Xb = y
  // suggested tolerance: 1.0e-7
  // if residuals or message are NULL, are not used
  int ext_findLeastSquaresFit(const double* y, ext_size_t n, const double* x, ext_size_t p, double* b,
                              double tolerance, double* residuals, char** message);
  int ext_findLeastSquaresFitInPlace(double* y, ext_size_t n, double* x, ext_size_t p, double* b,
                                     double tolerance, double* residuals, char** message);

  // returns 0 for success, EINVAL for an invalid argument, and EDOM if matrix is not symm pos def
  int ext_getSymmetricPositiveDefiniteTriangularFactorization(const double* restrict x, ext_size_t dim, ext_triangleType triangleType, double* restrict result);
  int ext_getSymmetricPositiveDefiniteTriangularFactorizationInPlace(double* x, ext_size_t dim, ext_triangleType triangleType);

  // X'X if use transpose is false, XX' otherwise
  // stores only upper, lower, or both parts of result depending
  void ext_getSingleMatrixCrossproduct(const double* restrict x, ext_size_t numRows, ext_size_t numCols,
                                       double* restrict result, int useTranspose, ext_triangleType triangleType);

  // x := solve(Ax = b) or solve(A'x = b), A triangular
  // returns 0 for success, EINVAL for an invalid argument, and EDOM if matrix is not symm pos def
  int ext_solveTriangularSystemInPlace(const double* restrict lhs, ext_size_t lhsDim, int useTranspose, ext_triangleType triangleType,
                                       double* restrict rhs, ext_size_t rhsDim);
                                       
  void ext_multiplyMatrixIntoVector(const double* restrict matrix, ext_size_t numRows, ext_size_t numCols, int useTranspose,
                                    const double* restrict vector, double* restrict result);

#ifdef __cplusplus
}
#endif
  
#endif // EXTERNAL_LINEAR_ALGEBRA_H
