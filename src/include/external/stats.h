#ifndef EXTERNAL_STATS_H
#define EXTERNAL_STATS_H

#include "stddef.h"

#include <stdint.h>

#include <Rmath.h>

#ifdef __cplusplus
extern "C" {
#endif

#define ext_simulateChiSquared(_NU_) rchisq((_NU_))
#define ext_simulateContinuousUniform() unif_rand()
#define ext_simulateBernoulli(_P_) (unif_rand() < (_P_) ? 1u : 0u)
  
#define ext_quantileOfChiSquared(_P_, _NU_) qchisq((_P_), (_NU_), 1, 0)
#define ext_percentileOfChiSquared(_Q_, _NU_) pchisq((_Q_), (_NU_), 1, 0)
  
#define ext_densityOfNormal(_X_, _MU_, _SIGMA_) dnorm((_X_), (_MU_), (_SIGMA_), 0)
#define ext_cumulativeProbabilityOfNormal(_Q_, _MU_, _SIGMA_) pnorm((_Q_), (_MU_), (_SIGMA_), 1, 0)
#define ext_quantileOfNormal(_P_, _MU_, _SIGMA_) qnorm((_P_), (_MU_), (_SIGMA_), 1, 0)
#define ext_simulateNormal(_MU_, _SIGMA_) rnorm((_MU_), (_SIGMA_))
#define ext_simulateStandardNormal() norm_rand()
  
#define EXT_DISCRETE_DRAW_FAILURE ((ext_size_t) -1)
ext_size_t ext_drawFromDiscreteDistribution(const double* probabilities, ext_size_t length);
  
  // random in [min, min + 1, ..., max - 1, max)
int64_t ext_simulateIntegerUniformInRange(int64_t min_inclusive, int64_t max_exclusive);
uint64_t ext_simulateUnsignedIntegerUniformInRange(uint64_t min_inclusive, uint64_t max_exclusive);
  
double ext_computeMean(const double* x, ext_size_t length);
double ext_computeIndexedMean(const double* restrict x, const ext_size_t* restrict indices, size_t length);
// variance := ssr / (n - 1)
double ext_computeVariance(const double* restrict x, ext_size_t length, double* restrict mean);
double ext_computeVarianceForKnownMean(const double* x, size_t length, double mean);
double ext_computeIndexedVariance(const double* restrict x, const ext_size_t* restrict indices, size_t length, double* restrict mean);
double ext_computeIndexedVarianceForKnownMean(const double* restrict x, const ext_size_t* restrict indices, size_t length, double mean);
  
double ext_computeAndSumSquaresOfResiduals(const double* y, ext_size_t length, double y_hat);
double ext_computeAndSumSquaresOfResidualsForVector(const double* restrict y, ext_size_t length, const double* restrict y_hat);
  
#ifdef __cplusplus
}
#endif

#endif // EXTERNAL_STATS_H
