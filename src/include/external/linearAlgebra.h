#ifndef EXTERNAL_LINEAR_ALGEBRA_H
#define EXTERNAL_LINEAR_ALGEBRA_H

#include <misc/stddef.h>

#ifdef __cplusplus
extern "C" {
#endif
  
  typedef enum {
    EXT_TRIANGLE_TYPE_BOTH,
    EXT_TRIANGLE_TYPE_UPPER,
    EXT_TRIANGLE_TYPE_LOWER
  } ext_triangleType;
  
  // y := alpha * x + y
  void ext_addVectorsInPlace(const double* restrict x, misc_size_t length, double alpha, double* restrict y);
   
  // x := x + alpha
  void ext_addScalarToVectorInPlace(double* x, misc_size_t length, double alpha);
  
  // b = Ax
  void ext_leftMultiplyMatrixAndVector(const double* A, misc_size_t n, misc_size_t p, const double* x, double* b);

  double ext_sumSquaresOfVectorElements(const double* x, misc_size_t length);
  
  void ext_transposeMatrix(const double* x, misc_size_t numRows, misc_size_t numCols, double* xt);
  
  // least squares solution to Xb = y
  // suggested tolerance: 1.0e-7
  // if residuals or message are NULL, are not used
  int ext_findLeastSquaresFit(const double* y, misc_size_t n, const double* x, misc_size_t p, double* b,
                              double tolerance, double* residuals, char** message);
  int ext_findLeastSquaresFitInPlace(double* y, misc_size_t n, double* x, misc_size_t p, double* b,
                                     double tolerance, double* residuals, char** message);

  // returns 0 for success, EINVAL for an invalid argument, and EDOM if matrix is not symm pos def
  int ext_getSymmetricPositiveDefiniteTriangularFactorization(const double* restrict x, misc_size_t dim, ext_triangleType triangleType, double* restrict result);
  int ext_getSymmetricPositiveDefiniteTriangularFactorizationInPlace(double* x, misc_size_t dim, ext_triangleType triangleType);

  // X'X if use transpose is false, XX' otherwise
  // stores only upper, lower, or both parts of result depending
  void ext_getSingleMatrixCrossproduct(const double* restrict x, misc_size_t numRows, misc_size_t numCols,
                                       double* restrict result, int useTranspose, ext_triangleType triangleType);

  // x := solve(Ax = b) or solve(A'x = b), A triangular
  // returns 0 for success, EINVAL for an invalid argument, and EDOM if matrix is not symm pos def
  int ext_solveTriangularSystemInPlace(const double* restrict lhs, misc_size_t lhsDim, int useTranspose, ext_triangleType triangleType,
                                       double* restrict rhs, misc_size_t rhsDim);
                                       
#ifdef __cplusplus
}
#endif
  
#endif // EXTERNAL_LINEAR_ALGEBRA_H

