#ifndef MISC_STATS_H
#define MISC_STATS_H

#include <misc/stddef.h>
#include <misc/thread.h>

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

// serial fast paths: the unrolled accumulators the misc_htm_* entry points
// run when given no (or a single-threaded) manager, minus the indirection
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

// multithreaded functions below

double misc_mt_computeMean               (misc_mt_manager_t restrict tm, const double* restrict x, misc_size_t length);
double misc_mt_computeIndexedMean        (misc_mt_manager_t restrict tm, const double* restrict x, const misc_index_t* restrict indices, misc_size_t length);
double misc_mt_computeWeightedMean       (misc_mt_manager_t restrict tm, const double* restrict x, misc_size_t length, const double* restrict w, double* restrict n);
double misc_mt_computeIndexedWeightedMean(misc_mt_manager_t restrict tm, const double* restrict x, const misc_index_t* restrict indices, misc_size_t length, const double* restrict w, double* restrict n);

double misc_mt_computeVariance       (misc_mt_manager_t restrict tm, const double* restrict x, misc_size_t length, double* restrict mean);
double misc_mt_computeIndexedVariance(misc_mt_manager_t restrict tm, const double* restrict x, const misc_index_t* restrict indices, misc_size_t length, double* restrict mean);

double misc_mt_computeVarianceForKnownMean               (misc_mt_manager_t restrict tm, const double* restrict x, misc_size_t length, double mean);
double misc_mt_computeIndexedVarianceForKnownMean        (misc_mt_manager_t restrict tm, const double* restrict x, const misc_index_t* restrict indices, misc_size_t length, double mean);
double misc_mt_computeWeightedVarianceForKnownMean       (misc_mt_manager_t restrict tm, const double* restrict x, misc_size_t length, const double* restrict w, double mean);
double misc_mt_computeIndexedWeightedVarianceForKnownMean(misc_mt_manager_t restrict tm, const double* restrict x, const misc_index_t* restrict indices, misc_size_t length, const double* restrict w, double mean);

double misc_mt_computeSumOfSquaredResiduals(misc_mt_manager_t restrict tm, const double* restrict x, misc_size_t length, const double* restrict x_hat);
double misc_mt_computeWeightedSumOfSquaredResiduals(misc_mt_manager_t restrict tm, const double* restrict x, misc_size_t length, const double* restrict w, const double* restrict x_hat);


// versions that use the hierarchical thread manager instead
double misc_htm_computeMean(
  misc_htm_manager_t restrict htm, misc_size_t taskId, const double* restrict x,
  misc_size_t length);

double misc_htm_computeIndexedMean(
  misc_htm_manager_t restrict htm, misc_size_t taskId, const double* restrict x,
  const misc_index_t* restrict indices, misc_size_t length);

double misc_htm_computeWeightedMean(
  misc_htm_manager_t restrict htm, misc_size_t taskId, const double* restrict x,
  misc_size_t length, const double* restrict w, double* restrict n);

double misc_htm_computeIndexedWeightedMean(
  misc_htm_manager_t restrict htm, misc_size_t taskId, const double* restrict x,
  const misc_index_t* restrict indices, misc_size_t length, const double* restrict w,
  double* restrict n);


double misc_htm_computeVarianceForKnownMean(
  misc_htm_manager_t restrict htm, misc_size_t taskId, const double* restrict x,
  misc_size_t length, double mean);

double misc_htm_computeIndexedVarianceForKnownMean(
  misc_htm_manager_t restrict htm, misc_size_t taskId, const double* restrict x,
  const misc_index_t* restrict indices, misc_size_t length, double mean);

double misc_htm_computeWeightedVarianceForKnownMean(
  misc_htm_manager_t restrict htm, misc_size_t taskId, const double* restrict x,
  misc_size_t length, const double* restrict w, double mean);

double misc_htm_computeIndexedWeightedVarianceForKnownMean(
  misc_htm_manager_t restrict htm, misc_size_t taskId, const double* restrict x,
  const misc_index_t* restrict indices, misc_size_t length, const double* restrict w,
  double mean);



double misc_htm_computeSumOfSquaredResiduals(
  misc_htm_manager_t restrict htm, misc_size_t taskId, const double* restrict x,
  misc_size_t length, const double* restrict x_hat);

double misc_htm_computeWeightedSumOfSquaredResiduals(
  misc_htm_manager_t restrict htm, misc_size_t taskId, const double* restrict x,
  misc_size_t length, const double* restrict w, const double* restrict x_hat);
  
#ifdef __cplusplus
}
#endif

#endif // MISC_STATS_MT_H

