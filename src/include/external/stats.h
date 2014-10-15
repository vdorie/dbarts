#ifndef EXTERNAL_STATS_H
#define EXTERNAL_STATS_H

#include "stddef.h"
#include <stdint.h>

#include <Rmath.h>

#ifdef __cplusplus
extern "C" {
#endif

#define ext_quantileOfChiSquared(_P_, _NU_) qchisq((_P_), (_NU_), 1, 0)
#define ext_percentileOfChiSquared(_Q_, _NU_) pchisq((_Q_), (_NU_), 1, 0)
  
#define ext_densityOfNormal(_X_, _MU_, _SIGMA_) dnorm((_X_), (_MU_), (_SIGMA_), 0)
#define ext_cumulativeProbabilityOfNormal(_Q_, _MU_, _SIGMA_) pnorm((_Q_), (_MU_), (_SIGMA_), 1, 0)
#define ext_quantileOfNormal(_P_, _MU_, _SIGMA_) qnorm((_P_), (_MU_), (_SIGMA_), 1, 0)

double ext_computeMean(const double* x, ext_size_t length);
// weighted mean = w'x / w'1; n will be set to w'1 if not NULL
double ext_computeIndexedMean(const double* restrict x, const ext_size_t* restrict indices, ext_size_t length);
double ext_computeWeightedMean(const double* restrict x, ext_size_t length, const double* restrict w, double* restrict n);
double ext_computeIndexedWeightedMean(const double* restrict x, const ext_size_t* restrict indices, ext_size_t length, const double* restrict w, double* restrict n);
  
// variance := ssr / (n - 1)
double ext_computeVariance(const double* restrict x, ext_size_t length, double* restrict mean);
double ext_computeVarianceForKnownMean(const double* x, ext_size_t length, double mean);
double ext_computeIndexedVariance(const double* restrict x, const ext_size_t* restrict indices, ext_size_t length, double* restrict mean);
double ext_computeIndexedVarianceForKnownMean(const double* restrict x, const ext_size_t* restrict indices, ext_size_t length, double mean);
double ext_computeWeightedVarianceForKnownMean(const double* restrict x, ext_size_t length, const double* restrict w, double mean);
double ext_computeIndexedWeightedVarianceForKnownMean(const double* restrict x, const ext_size_t* restrict indices, ext_size_t length, const double* restrict w, double mean);
  
// double ext_computeAndSumSquaresOfResiduals(const double* y, ext_size_t length, double y_hat);
// double ext_computeAndSumSquaresOfResidualsForVector(const double* restrict y, ext_size_t length, const double* restrict y_hat);
double ext_computeSumOfSquaredResiduals(const double* restrict x, ext_size_t length, const double* restrict x_hat);
double ext_computeWeightedSumOfSquaredResiduals(const double* restrict x, ext_size_t length, const double* restrict w, const double* restrict x_hat);
  
#ifdef __cplusplus
}
#endif

#endif // EXTERNAL_STATS_H
