#ifndef EXTERNAL_STATS_MT_H
#define EXTERNAL_STATS_MT_H

// alternative implementations for multithreaded environments

#include "stddef.h"
#include "thread.h"

#ifdef __cplusplus
extern "C" {
#endif

double ext_mt_computeMean(ext_mt_manager_t restrict threadManager, const double* restrict x, ext_size_t length);
double ext_mt_computeIndexedMean(ext_mt_manager_t restrict threadManager, const double* restrict x, const ext_size_t* restrict indices, size_t length);

double ext_mt_computeVariance(ext_mt_manager_t restrict threadManager, const double* restrict x, ext_size_t length, double* mean);
double ext_mt_computeVarianceForKnownMean(ext_mt_manager_t restrict threadManager, const double* restrict x, ext_size_t length, double mean);
  
double ext_mt_computeIndexedVariance(ext_mt_manager_t restrict threadManager, const double* restrict x, const size_t* restrict indices, ext_size_t length, double* mean);
double ext_mt_computeIndexedVarianceForKnownMean(ext_mt_manager_t restrict threadManager, const double* restrict x, const size_t* restrict indices, ext_size_t length, double mean);

double ext_mt_computeAndSumSquaresOfResiduals(ext_mt_manager_t restrict threadManager, const double* restrict x, size_t length, double x_hat);
#ifdef __cplusplus
}
#endif

#endif // EXTERNAL_STATS_MT_H
