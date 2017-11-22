#ifndef EXTERNAL_STATS_MT_H
#define EXTERNAL_STATS_MT_H

// alternative implementations for multithreaded environments

#include "stddef.h"
#include "thread.h"

#ifdef __cplusplus
extern "C" {
#endif

double ext_mt_computeMean               (ext_mt_manager_t restrict tm, const double* restrict x, ext_size_t length);
double ext_mt_computeIndexedMean        (ext_mt_manager_t restrict tm, const double* restrict x, const ext_size_t* restrict indices, ext_size_t length);
double ext_mt_computeWeightedMean       (ext_mt_manager_t restrict tm, const double* restrict x, ext_size_t length, const double* restrict w, double* restrict n);
double ext_mt_computeIndexedWeightedMean(ext_mt_manager_t restrict tm, const double* restrict x, const ext_size_t* restrict indices, ext_size_t length, const double* restrict w, double* restrict n);

double ext_mt_computeVariance       (ext_mt_manager_t restrict tm, const double* restrict x, ext_size_t length, double* restrict mean);
double ext_mt_computeIndexedVariance(ext_mt_manager_t restrict tm, const double* restrict x, const ext_size_t* restrict indices, ext_size_t length, double* restrict mean);

double ext_mt_computeVarianceForKnownMean               (ext_mt_manager_t restrict tm, const double* restrict x, ext_size_t length, double mean);
double ext_mt_computeIndexedVarianceForKnownMean        (ext_mt_manager_t restrict tm, const double* restrict x, const ext_size_t* restrict indices, ext_size_t length, double mean);
double ext_mt_computeWeightedVarianceForKnownMean       (ext_mt_manager_t restrict tm, const double* restrict x, ext_size_t length, const double* restrict w, double mean);
double ext_mt_computeIndexedWeightedVarianceForKnownMean(ext_mt_manager_t restrict tm, const double* restrict x, const ext_size_t* restrict indices, ext_size_t length, const double* restrict w, double mean);

double ext_mt_computeSumOfSquaredResiduals(ext_mt_manager_t restrict tm, const double* restrict x, ext_size_t length, const double* restrict x_hat);
double ext_mt_computeWeightedSumOfSquaredResiduals(ext_mt_manager_t restrict tm, const double* restrict x, ext_size_t length, const double* restrict w, const double* restrict x_hat);


// versions that use the hierarchical thread manager instead
double ext_htm_computeMean               (ext_htm_manager_t restrict htm, ext_size_t taskId,
                                          const double* restrict x, ext_size_t length);
double ext_htm_computeIndexedMean        (ext_htm_manager_t restrict htm, ext_size_t taskId,
                                          const double* restrict x, const ext_size_t* restrict indices, ext_size_t length);
double ext_htm_computeWeightedMean       (ext_htm_manager_t restrict htm, ext_size_t taskId,
                                          const double* restrict x, ext_size_t length, const double* restrict w, double* restrict n);
double ext_htm_computeIndexedWeightedMean(ext_htm_manager_t restrict htm, ext_size_t taskId,
                                          const double* restrict x, const ext_size_t* restrict indices, ext_size_t length,
                                          const double* restrict w, double* restrict n);

double ext_htm_computeVarianceForKnownMean               (ext_htm_manager_t restrict htm, ext_size_t taskId,
                                                          const double* restrict x, ext_size_t length, double mean);
double ext_htm_computeIndexedVarianceForKnownMean        (ext_htm_manager_t restrict htm, ext_size_t taskId,
                                                          const double* restrict x, const ext_size_t* restrict indices, ext_size_t length,
                                                          double mean);
double ext_htm_computeWeightedVarianceForKnownMean       (ext_htm_manager_t restrict htm, ext_size_t taskId,
                                                          const double* restrict x, ext_size_t length, const double* restrict w, double mean);
double ext_htm_computeIndexedWeightedVarianceForKnownMean(ext_htm_manager_t restrict htm, ext_size_t taskId,
                                                          const double* restrict x, const ext_size_t* restrict indices, ext_size_t length,
                                                          const double* restrict w, double mean);



double ext_htm_computeSumOfSquaredResiduals        (ext_htm_manager_t restrict htm, ext_size_t taskId,
                                                    const double* restrict x, ext_size_t length, const double* restrict x_hat);
double ext_htm_computeWeightedSumOfSquaredResiduals(ext_htm_manager_t restrict htm, ext_size_t taskId,
                                                    const double* restrict x, ext_size_t length, const double* restrict w,
                                                    const double* restrict x_hat);
  
#ifdef __cplusplus
}
#endif

#endif // EXTERNAL_STATS_MT_H

