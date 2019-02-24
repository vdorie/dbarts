#ifndef MISC_LINEAR_ALGEBRA_H
#define MISC_LINEAR_ALGEBRA_H

#include <stdbool.h>
#include <misc/stddef.h>

#include <external/linearAlgebra.h>

#ifdef __cplusplus
extern "C" {
#endif

// z = alpha * x + y; z must be distinct from x and y
// void misc_addVectors(const double* restrict x, misc_size_t length, double alpha, const double* restrict y, double* restrict z);
extern void (*misc_addVectors)(const double* restrict x, misc_size_t length, double alpha, const double* restrict y, double* restrict z);
 
// x: = alpha
bool misc_vectorIsConstant(const double* d, misc_size_t length);
// void misc_setVectorToConstant(double* x, misc_size_t length, double alpha);
extern void (*misc_setVectorToConstant)(double* x, misc_size_t length, double alpha);
void misc_setIndexedVectorToConstant(double* restrict x, const misc_size_t* restrict indices, misc_size_t length, double alpha);
  
// x := alpha * x
void misc_scalarMultiplyVectorInPlace(double* x, misc_size_t length, double alpha);
// y = alpha * x; y must be distinct from x
void misc_scalarMultiplyVector(const double* restrict x, misc_size_t length, double alpha, double* restrict y);
 
 
// z = x .* y; z must be distinct from x and y
void misc_hadamardMultiplyVectors(const double* restrict x, misc_size_t length, const double* restrict y, double* restrict z);
// x := x. * y; x and y should be distinct
void misc_hadamardMultiplyVectorsInPlace(double* restrict x, misc_size_t length, const double* restrict y);

 
double misc_sumVectorElements(const double* x, misc_size_t length);
double misc_sumIndexedVectorElements(const double* x, const misc_size_t* indices, misc_size_t length);

void misc_multiplyMatrixIntoVector(const double* restrict matrix, misc_size_t numRows, misc_size_t numCols, int useTranspose,
                                   const double* restrict vector, double* restrict result);

#ifdef __cplusplus
}
#endif

#define misc_addVectorsInPlace ext_addVectorsInPlace
#define misc_addScalarToVectorInPlace ext_addScalarToVectorInPlace
#define misc_sumSquaresOfVectorElements ext_sumSquaresOfVectorElements
#define misc_transposeMatrix ext_transposeMatrix

#endif // define MISC_LINEAR_ALGEBRA_H

