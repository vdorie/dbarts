#ifndef MISC_STATS_H
#define MISC_STATS_H

#include <misc/stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

double misc_computeMean(const double* x, misc_size_t length);
double misc_computeIndexedMean(const double* restrict x, const misc_index_t* restrict indices, misc_size_t length);
// weighted mean = w'x / w'1; n will be set to w'1 if not NULL
double misc_computeWeightedMean(const double* restrict x, misc_size_t length, const double* restrict w, double* restrict n);
double misc_computeIndexedWeightedMean(const double* restrict x, const misc_index_t* restrict indices, misc_size_t length, const double* restrict w, double* restrict n);
  
// variance := ssr / (n - 1); renormalize by (n - 1) / w'1 for a weighted estimate
double misc_computeVariance(const double* restrict x, misc_size_t length, double* restrict mean);
double misc_computeVarianceForKnownMean(const double* x, misc_size_t length, double mean);
double misc_computeIndexedVariance(const double* restrict x, const misc_index_t* restrict indices, misc_size_t length, double* restrict mean);
double misc_computeIndexedVarianceForKnownMean(const double* restrict x, const misc_index_t* restrict indices, misc_size_t length, double mean);
double misc_computeWeightedVarianceForKnownMean(const double* restrict x, misc_size_t length, const double* restrict w, double mean);
double misc_computeIndexedWeightedVarianceForKnownMean(const double* restrict x, const misc_index_t* restrict indices, misc_size_t length, const double* restrict w, double mean);
  
double misc_computeSumOfSquaredResiduals(const double* restrict x, misc_size_t length, const double* restrict x_hat);
double misc_computeWeightedSumOfSquaredResiduals(const double* restrict x, misc_size_t length, const double* restrict w, const double* restrict x_hat);

// fast paths: the unrolled accumulators the entry points above reach through
// a length-keyed online/vanilla dispatch, minus that indirection
double misc_computeMeanFast(const double* x, misc_size_t length);
double misc_computeIndexedMeanFast(const double* restrict x, const misc_index_t* restrict indices, misc_size_t length);
double misc_computeWeightedMeanFast(const double* restrict x, misc_size_t length, const double* restrict w, double* restrict n);
double misc_computeIndexedWeightedMeanFast(const double* restrict x, const misc_index_t* restrict indices, misc_size_t length, const double* restrict w, double* restrict n);

double misc_computeVarianceForKnownMeanFast(const double* x, misc_size_t length, double mean);
double misc_computeIndexedVarianceForKnownMeanFast(const double* restrict x, const misc_index_t* restrict indices, misc_size_t length, double mean);
double misc_computeWeightedVarianceForKnownMeanFast(const double* restrict x, misc_size_t length, const double* restrict w, double mean);
double misc_computeIndexedWeightedVarianceForKnownMeanFast(const double* restrict x, const misc_index_t* restrict indices, misc_size_t length, const double* restrict w, double mean);

// fused constant-leaf sufficient statistic in one order-insensitive pass:
// sumW = w'1 (= length unweighted), sumWX = w'x
void misc_computeSufficientStatisticsFast(const double* x, misc_size_t length, double* restrict sumW, double* restrict sumWX);
void misc_computeIndexedSufficientStatisticsFast(const double* restrict x, const misc_index_t* restrict indices, misc_size_t length, double* restrict sumW, double* restrict sumWX);
void misc_computeWeightedSufficientStatisticsFast(const double* restrict x, misc_size_t length, const double* restrict w, double* restrict sumW, double* restrict sumWX);
void misc_computeIndexedWeightedSufficientStatisticsFast(const double* restrict x, const misc_index_t* restrict indices, misc_size_t length, const double* restrict w, double* restrict sumW, double* restrict sumWX);

// fp32-residual (opt-in) variants: identical reduction, but the running
// residual x is stored fp32 - each element is loaded as float and promoted to
// double before accumulating, so the reduction stays fp64
void misc_computeFloatSufficientStatisticsFast(const float* x, misc_size_t length, double* restrict sumW, double* restrict sumWX);
void misc_computeIndexedFloatSufficientStatisticsFast(const float* restrict x, const misc_index_t* restrict indices, misc_size_t length, double* restrict sumW, double* restrict sumWX);
void misc_computeWeightedFloatSufficientStatisticsFast(const float* restrict x, misc_size_t length, const double* restrict w, double* restrict sumW, double* restrict sumWX);
void misc_computeIndexedWeightedFloatSufficientStatisticsFast(const float* restrict x, const misc_index_t* restrict indices, misc_size_t length, const double* restrict w, double* restrict sumW, double* restrict sumWX);

#ifdef __cplusplus
}
#endif

#endif // MISC_STATS_H

