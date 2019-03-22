#include "config.h"

#include <math.h>
#include <stdint.h>

#include <misc/stats.h>
#include <misc/thread.h>

#if (__GNUC__ > 2) || (__GNUC__ == 2 && __GNUC_MINOR__ > 4)
#  define UNUSED __attribute__ ((unused))
#else
#  define UNUSED
#endif

// This file contains mean and variance functions with a variety of implementations.
// Of those variants, we have for each a vanilla version, a vanilla version
// with unrolled loops, an "online" version that computes running averages,
// and an online version that also uses unrolled loops.
// 
// Furthermore, each is multithreaded. Depending on the success of the loop
// unrolling, often a single core is incredibly efficient hence there are
// variable points which switching to a multithreaded approach makes sense.


// around these values, multithreaded starts to beat single threaded
// #define                 MEAN_MIN_NUM_VALUES_PER_THREAD 100000
#define        UNROLLED_MEAN_MIN_NUM_VALUES_PER_THREAD 200000
// #define          ONLINE_MEAN_MIN_NUM_VALUES_PER_THREAD 7500
#define ONLINE_UNROLLED_MEAN_MIN_NUM_VALUES_PER_THREAD 25000

// #define                 VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD 75000
#define        UNROLLED_VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD 75000
// #define          ONLINE_VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD 7500
#define ONLINE_UNROLLED_VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD 25000

// #define                 VAR_MIN_NUM_VALUES_PER_THREAD 75000
#define        UNROLLED_VAR_MIN_NUM_VALUES_PER_THREAD 75000
// #define          ONLINE_VAR_MIN_NUM_VALUES_PER_THREAD 5000
#define ONLINE_UNROLLED_VAR_MIN_NUM_VALUES_PER_THREAD 12000

// #define                 WEIGHTED_MEAN_MIN_NUM_VALUES_PER_THREAD 150000
#define        UNROLLED_WEIGHTED_MEAN_MIN_NUM_VALUES_PER_THREAD 125000
// #define          ONLINE_WEIGHTED_MEAN_MIN_NUM_VALUES_PER_THREAD 7500
#define ONLINE_UNROLLED_WEIGHTED_MEAN_MIN_NUM_VALUES_PER_THREAD 35000

// #define                 WEIGHTED_VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD 150000
#define        UNROLLED_WEIGHTED_VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD 125000
// #define          ONLINE_WEIGHTED_VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD 7500
#define ONLINE_UNROLLED_WEIGHTED_VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD 35000

// #define                 INDEXED_MEAN_MIN_NUM_VALUES_PER_THREAD 100000
#define        INDEXED_UNROLLED_MEAN_MIN_NUM_VALUES_PER_THREAD 100000
// #define          INDEXED_ONLINE_MEAN_MIN_NUM_VALUES_PER_THREAD 7500
#define INDEXED_ONLINE_UNROLLED_MEAN_MIN_NUM_VALUES_PER_THREAD 25000

// #define                 INDEXED_VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD 75000
#define        INDEXED_UNROLLED_VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD 75000
// #define          INDEXED_ONLINE_VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD 7500
#define INDEXED_ONLINE_UNROLLED_VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD 25000

// #define                 INDEXED_VAR_MIN_NUM_VALUES_PER_THREAD 75000
#define        INDEXED_UNROLLED_VAR_MIN_NUM_VALUES_PER_THREAD 75000
// #define          INDEXED_ONLINE_VAR_MIN_NUM_VALUES_PER_THREAD 5000
#define INDEXED_ONLINE_UNROLLED_VAR_MIN_NUM_VALUES_PER_THREAD 12000

// #define                 INDEXED_WEIGHTED_MEAN_MIN_NUM_VALUES_PER_THREAD 60000
#define        INDEXED_UNROLLED_WEIGHTED_MEAN_MIN_NUM_VALUES_PER_THREAD 35000
// #define          INDEXED_ONLINE_WEIGHTED_MEAN_MIN_NUM_VALUES_PER_THREAD 7500
#define INDEXED_ONLINE_UNROLLED_WEIGHTED_MEAN_MIN_NUM_VALUES_PER_THREAD 30000

// #define                 INDEXED_WEIGHTED_VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD 60000
#define        INDEXED_UNROLLED_WEIGHTED_VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD 35000
// #define          INDEXED_ONLINE_WEIGHTED_VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD 7500
#define INDEXED_ONLINE_UNROLLED_WEIGHTED_VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD 25000

#define          SSR_MIN_NUM_VALUES_PER_THREAD 75000
#define WEIGHTED_SSR_MIN_NUM_VALUES_PER_THREAD 60000

// if < value, calculate straightup; otherwise use online
// algorithm to reduce round-off err
#define ONLINE_CUTOFF 10000
#define INDEXED_ONLINE_CUTOFF 10000


// various implementations
// mean
UNUSED static double computeMean                     (const double* restrict x, size_t length);
UNUSED static double computeIndexedMean              (const double* restrict x, const size_t* restrict indices, size_t length);
static double (*computeUnrolledMean)(const double* x, size_t length) = 0;
static double (*computeIndexedUnrolledMean)(const double* restrict x, const size_t* restrict indices, size_t length) = 0;
UNUSED static double computeOnlineMean               (const double* restrict x, size_t length);
UNUSED static double computeIndexedOnlineMean        (const double* restrict x, const size_t* restrict indices, size_t length);
static double (*computeOnlineUnrolledMean)(const double* restrict x, size_t length) = 0;
static double (*computeIndexedOnlineUnrolledMean)(const double* restrict x, const size_t* restrict indices, size_t length) = 0;

// var for known mean
UNUSED static double computeVarianceForKnownMean                     (const double* x, size_t length, double mean);
UNUSED static double computeIndexedVarianceForKnownMean              (const double* restrict x, const size_t* restrict indices, size_t length, double mean);
static double (*computeUnrolledVarianceForKnownMean)(const double* x, size_t length, double mean) = 0;
static double (*computeIndexedUnrolledVarianceForKnownMean)(const double* restrict x, const size_t* restrict indices, size_t length, double mean) = 0;
UNUSED static double computeOnlineVarianceForKnownMean               (const double* restrict x, size_t length, double mean);
UNUSED static double computeIndexedOnlineVarianceForKnownMean        (const double* restrict x, const size_t* restrict indices, size_t length, double mean);
static double (*computeOnlineUnrolledVarianceForKnownMean)(const double* restrict x, size_t length, double mean) = 0;
static double (*computeIndexedOnlineUnrolledVarianceForKnownMean)(const double* restrict x, const size_t* restrict indices, size_t length, double mean) = 0;

// variance and mean together
UNUSED static double computeVariance                     (const double* restrict x, size_t length, double* restrict meanPtr);
UNUSED static double computeIndexedVariance              (const double* restrict x, const size_t* restrict indices, size_t length, double* restrict meanPtr);
static double computeUnrolledVariance             (const double* restrict x, size_t length, double* restrict meanPtr);
static double computeIndexedUnrolledVariance      (const double* restrict x, const size_t* restrict indices, size_t length, double* restrict meanPtr);
UNUSED static double computeOnlineVariance               (const double* restrict x, size_t length, double* restrict meanPtr);
UNUSED static double computeIndexedOnlineVariance        (const double* restrict x, const size_t* restrict indices, size_t length, double* restrict meanPtr);
static double computeOnlineUnrolledVariance       (const double* restrict x, size_t length, double* restrict meanPtr);
static double computeIndexedOnlineUnrolledVariance(const double* restrict x, const size_t* restrict indices, size_t length, double* restrict meanPtr);

// weighted mean
UNUSED static double computeWeightedMean                     (const double* restrict x, size_t length, const double* restrict w, double* restrict nPtr);
UNUSED static double computeIndexedWeightedMean              (const double* restrict x, const size_t* restrict indices, size_t length, const double* restrict w, double* restrict nPtr);
static double (*computeUnrolledWeightedMean)(const double* restrict x, size_t length, const double* restrict w, double* restrict nPtr) = 0;
static double (*computeIndexedUnrolledWeightedMean)(const double* restrict x, const size_t* restrict indices, size_t length, const double* restrict w, double* restrict nPtr) = 0;
UNUSED static double computeOnlineWeightedMean               (const double* restrict x, size_t length, const double* restrict w, double* restrict nPtr);
UNUSED static double computeIndexedOnlineWeightedMean        (const double* restrict x, const size_t* restrict indices, size_t length, const double* restrict w, double* restrict nPtr);
static double (*computeOnlineUnrolledWeightedMean)(const double* restrict x, size_t length, const double* restrict w, double* restrict nPtr) = 0;
static double (*computeIndexedOnlineUnrolledWeightedMean)(const double* restrict x, const size_t* restrict indices, size_t length, const double* restrict w, double* restrict nPtr) = 0;

// weighted variance for known mean
UNUSED static double computeWeightedVarianceForKnownMean                     (const double* restrict x, size_t length, const double* restrict w, double mean);
UNUSED static double computeIndexedWeightedVarianceForKnownMean              (const double* restrict x, const size_t* restrict indices, size_t length, const double* restrict w, double mean);
static double (*computeUnrolledWeightedVarianceForKnownMean)(const double* restrict x, size_t length, const double* restrict w, double mean) = 0;
static double (*computeIndexedUnrolledWeightedVarianceForKnownMean)(const double* restrict x, const size_t* restrict indices, size_t length, const double* restrict w, double mean) = 0;
UNUSED static double computeOnlineWeightedVarianceForKnownMean               (const double* restrict x, size_t length, const double* restrict w, double mean);
UNUSED static double computeIndexedOnlineWeightedVarianceForKnownMean        (const double* restrict x, const size_t* restrict indices, size_t length, const double* restrict w, double mean);
static double (*computeOnlineUnrolledWeightedVarianceForKnownMean)(const double* restrict x, size_t length, const double* restrict w, double mean) = 0;
static double (*computeIndexedOnlineUnrolledWeightedVarianceForKnownMean)(const double* restrict x, const size_t* restrict indices, size_t length, const double* restrict w, double mean) = 0;

// static double mt_computeMean(misc_mt_manager_t restrict threadManager, const double* restrict x, size_t length);
// static double mt_computeIndexedMean(misc_mt_manager_t restrict threadManager, const double* restrict x, const size_t* indices, size_t length);
static double mt_computeUnrolledMean(misc_mt_manager_t restrict threadManager, const double* restrict x, size_t length);
static double mt_computeIndexedUnrolledMean(misc_mt_manager_t restrict threadManager, const double* restrict x, const size_t* indices, size_t length);
// static double mt_computeOnlineMean(misc_mt_manager_t restrict threadManager, const double* restrict x, size_t length);
// static double mt_computeIndexedOnlineMean(misc_mt_manager_t restrict threadManager, const double* restrict x, const size_t* indices, size_t length);
static double mt_computeOnlineUnrolledMean(misc_mt_manager_t restrict threadManager, const double* restrict x, size_t length);
static double mt_computeIndexedOnlineUnrolledMean(misc_mt_manager_t restrict threadManager, const double* restrict x, const size_t* indices, size_t length);

// static double mt_computeVarianceForKnownMean(misc_mt_manager_t restrict threadManager, const double* restrict x, size_t length, double mean);
// static double mt_computeIndexedVarianceForKnownMean(misc_mt_manager_t restrict threadManager, const double* restrict x, const size_t* restrict indices, size_t length, double mean);
static double mt_computeUnrolledVarianceForKnownMean(misc_mt_manager_t restrict threadManager, const double* restrict x, size_t length, double mean);
static double mt_computeIndexedUnrolledVarianceForKnownMean(misc_mt_manager_t restrict threadManager, const double* restrict x, const size_t* restrict indices, size_t length, double mean);
// static double mt_computeOnlineVarianceForKnownMean(misc_mt_manager_t restrict threadManager, const double* restrict x, size_t length, double mean);
// static double mt_computeIndexedOnlineVarianceForKnownMean(misc_mt_manager_t restrict threadManager, const double* restrict x, const size_t* restrict indices, size_t length, double mean);
static double mt_computeOnlineUnrolledVarianceForKnownMean(misc_mt_manager_t restrict threadManager, const double* restrict x, size_t length, double mean);
static double mt_computeIndexedOnlineUnrolledVarianceForKnownMean(misc_mt_manager_t restrict threadManager, const double* restrict x, const size_t* restrict indices, size_t length, double mean);

// static double mt_computeVariance(misc_mt_manager_t restrict threadManager, const double* restrict x, size_t length, double* restrict meanPtr);
// static double mt_computeIndexedVariance(misc_mt_manager_t restrict threadManager, const double* restrict x, const size_t* restrict indices, size_t length, double* restrict mean);
static double mt_computeUnrolledVariance(misc_mt_manager_t restrict threadManager, const double* restrict x, size_t length, double* restrict meanPtr);
static double mt_computeIndexedUnrolledVariance(misc_mt_manager_t restrict threadManager, const double* restrict x, const size_t* restrict indices, size_t length, double* restrict mean);
// static double mt_computeOnlineVariance(misc_mt_manager_t restrict threadManager, const double* restrict x, size_t length, double* restrict meanPtr);
// static double mt_computeIndexedOnlineVariance(misc_mt_manager_t restrict threadManager, const double* restrict x, const size_t* restrict indices, size_t length, double* restrict mean);
static double mt_computeOnlineUnrolledVariance(misc_mt_manager_t restrict threadManager, const double* restrict x, size_t length, double* restrict meanPtr);
static double mt_computeIndexedOnlineUnrolledVariance(misc_mt_manager_t restrict threadManager, const double* restrict x, const size_t* restrict indices, size_t length, double* restrict mean);

// static double mt_computeWeightedMean(misc_mt_manager_t restrict threadManager, const double* restrict x, size_t length, const double* restrict w, double* restrict nPtr);
// static double mt_computeIndexedWeightedMean(misc_mt_manager_t restrict threadManager, const double* restrict x, const size_t* restrict indices, size_t length, const double* restrict w, double* restrict nPtr);
static double mt_computeUnrolledWeightedMean(misc_mt_manager_t restrict threadManager, const double* restrict x, size_t length, const double* restrict w, double* restrict nPtr);
static double mt_computeIndexedUnrolledWeightedMean(misc_mt_manager_t restrict threadManager, const double* restrict x, const size_t* restrict indices, size_t length, const double* restrict w, double* restrict nPtr);
// static double mt_computeOnlineWeightedMean(misc_mt_manager_t restrict threadManager, const double* restrict x, size_t length, const double* restrict w, double* restrict nPtr);
// static double mt_computeIndexedOnlineWeightedMean(misc_mt_manager_t restrict threadManager, const double* restrict x, const size_t* restrict indices, size_t length, const double* restrict w, double* restrict nPtr);
static double mt_computeOnlineUnrolledWeightedMean(misc_mt_manager_t restrict threadManager, const double* restrict x, size_t length, const double* restrict w, double* restrict nPtr);
static double mt_computeIndexedOnlineUnrolledWeightedMean(misc_mt_manager_t restrict threadManager, const double* restrict x, const size_t* restrict indices, size_t length, const double* restrict w, double* restrict nPtr);

// static double mt_computeWeightedVarianceForKnownMean(misc_mt_manager_t restrict threadManager, const double* restrict x, size_t length, const double* restrict w, double mean);
// static double mt_computeIndexedWeightedVarianceForKnownMean(misc_mt_manager_t restrict threadManager, const double* restrict x, const size_t* restrict indices, size_t length, const double* restrict w, double mean);
static double mt_computeUnrolledWeightedVarianceForKnownMean(misc_mt_manager_t restrict threadManager, const double* restrict x, size_t length, const double* restrict w, double mean);
static double mt_computeIndexedUnrolledWeightedVarianceForKnownMean(misc_mt_manager_t restrict threadManager, const double* restrict x, const size_t* restrict indices, size_t length, const double* restrict w, double mean);
// static double mt_computeOnlineWeightedVarianceForKnownMean(misc_mt_manager_t restrict threadManager, const double* restrict x, size_t length, const double* restrict w, double mean);
// static double mt_computeIndexedOnlineWeightedVarianceForKnownMean(misc_mt_manager_t restrict threadManager, const double* restrict x, const size_t* restrict indices, size_t length, const double* restrict w, double mean);
static double mt_computeOnlineUnrolledWeightedVarianceForKnownMean(misc_mt_manager_t restrict threadManager, const double* restrict x, size_t length, const double* restrict w, double mean);
static double mt_computeIndexedOnlineUnrolledWeightedVarianceForKnownMean(misc_mt_manager_t restrict threadManager, const double* restrict x, const size_t* restrict indices, size_t length, const double* restrict w, double mean);

// interface functions that dispatch to the workers
double misc_computeMean(const double* x, size_t length)
{
  if (length > ONLINE_CUTOFF) return computeOnlineUnrolledMean(x, length);
  return computeUnrolledMean(x, length);
}

double misc_computeIndexedMean(const double* restrict x, const size_t* restrict indices, size_t length)
{
  if (length > ONLINE_CUTOFF) return computeIndexedOnlineUnrolledMean(x, indices, length);
  return computeIndexedUnrolledMean(x, indices, length);
}

double misc_computeVarianceForKnownMean(const double* x, size_t length, double mean)
{
  if (length > ONLINE_CUTOFF) return computeOnlineUnrolledVarianceForKnownMean(x, length, mean);
  return computeUnrolledVarianceForKnownMean(x, length, mean);
}

double misc_computeIndexedVarianceForKnownMean(const double* restrict x, const misc_size_t* restrict indices, size_t length, double mean)
{
  if (length > ONLINE_CUTOFF) computeIndexedUnrolledVarianceForKnownMean(x, indices, length, mean);
  return computeIndexedOnlineUnrolledVarianceForKnownMean(x, indices, length, mean);
}

// The two-pass is faster when using online algorithms, one-pass when not
double misc_computeVariance(const double* restrict x, size_t length, double* restrict meanPtr)
{
  if (length > ONLINE_CUTOFF) {
    double mean = computeOnlineUnrolledMean(x, length);
    if (meanPtr != NULL) *meanPtr = mean;
    return computeOnlineUnrolledVarianceForKnownMean(x, length, mean);
  }
  
  return computeUnrolledVariance(x, length, meanPtr);
}

double misc_computeIndexedVariance(const double* restrict x, const size_t* restrict indices, size_t length, double* restrict meanPtr)
{
  if (length > ONLINE_CUTOFF) return computeIndexedOnlineUnrolledVariance(x, indices, length, meanPtr);
  return computeIndexedUnrolledVariance(x, indices, length, meanPtr);
}

double misc_computeWeightedMean(const double* restrict x, size_t length, const double* restrict w, double* restrict n)
{
  if (length > ONLINE_CUTOFF) return computeOnlineUnrolledWeightedMean(x, length, w, n);
  return computeUnrolledWeightedMean(x, length, w, n);
}

double misc_computeIndexedWeightedMean(const double* restrict x, const size_t* restrict indices, size_t length, const double* restrict w, double* restrict n)
{
  if (length > ONLINE_CUTOFF) return computeIndexedOnlineUnrolledWeightedMean(x, indices, length, w, n);
  return computeIndexedUnrolledWeightedMean(x, indices, length, w, n);
}

double misc_computeWeightedVarianceForKnownMean(const double* restrict x, size_t length, const double* restrict w, double mean)
{
  if (length > ONLINE_CUTOFF) return computeOnlineUnrolledWeightedVarianceForKnownMean(x, length, w, mean);
  return computeUnrolledWeightedVarianceForKnownMean(x, length, w, mean);
}

double misc_computeIndexedWeightedVarianceForKnownMean(const double* restrict x, const size_t* restrict indices, size_t length, const double* restrict w, double mean)
{
  if (length > ONLINE_CUTOFF) return computeIndexedOnlineUnrolledWeightedVarianceForKnownMean(x, indices, length, w, mean);
  return computeIndexedUnrolledWeightedVarianceForKnownMean(x, indices, length, w, mean);
}

#define minimum(_A_, _B_) ((_A_) < (_B_) ? (_A_) : (_B_))

// if the data for any thread would, by itself, trigger a fall-back to single threaded
// and that single-threaded function equiv would prefer the non-online version, do that instead
double misc_mt_computeMean(misc_mt_manager_t restrict threadManager, const double* restrict x, size_t length)
{
  size_t numThreads = misc_mt_getNumThreads(threadManager);
  size_t onlineCutoff = minimum(ONLINE_UNROLLED_MEAN_MIN_NUM_VALUES_PER_THREAD, ONLINE_CUTOFF);
  
  if (length / numThreads >= onlineCutoff) return mt_computeOnlineUnrolledMean(threadManager, x, length);
  return mt_computeUnrolledMean(threadManager, x, length);
}

double misc_mt_computeIndexedMean(misc_mt_manager_t restrict threadManager, const double* restrict x, const size_t* restrict indices, size_t length)
{
  size_t numThreads = misc_mt_getNumThreads(threadManager);
  size_t onlineCutoff = minimum(INDEXED_UNROLLED_MEAN_MIN_NUM_VALUES_PER_THREAD, INDEXED_ONLINE_CUTOFF);
  
  if (length / numThreads >= onlineCutoff) return mt_computeIndexedOnlineUnrolledMean(threadManager, x, indices, length);
  return mt_computeIndexedUnrolledMean(threadManager, x, indices, length);
}

double misc_mt_computeVarianceForKnownMean(misc_mt_manager_t restrict threadManager, const double* restrict x, size_t length, double mean)
{
  size_t numThreads = misc_mt_getNumThreads(threadManager);
  size_t onlineCutoff = minimum(UNROLLED_VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD, ONLINE_CUTOFF);
  
  if (length / numThreads >= onlineCutoff) return mt_computeOnlineUnrolledVarianceForKnownMean(threadManager, x, length, mean);
  return mt_computeUnrolledVarianceForKnownMean(threadManager, x, length, mean);
}

double misc_mt_computeIndexedVarianceForKnownMean(misc_mt_manager_t restrict threadManager, const double* restrict x, const size_t* restrict indices, size_t length, double mean)
{
  size_t numThreads = misc_mt_getNumThreads(threadManager);
  size_t onlineCutoff = minimum(INDEXED_UNROLLED_VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD, INDEXED_ONLINE_CUTOFF);
  
  if (length / numThreads >= onlineCutoff) return mt_computeIndexedOnlineUnrolledVarianceForKnownMean(threadManager, x, indices, length, mean);
  return mt_computeIndexedUnrolledVarianceForKnownMean(threadManager, x, indices, length, mean);
}

double misc_mt_computeVariance(misc_mt_manager_t restrict threadManager, const double* restrict x, size_t length, double* meanPtr)
{
  size_t numThreads = misc_mt_getNumThreads(threadManager);
  size_t onlineCutoff = minimum(UNROLLED_VAR_MIN_NUM_VALUES_PER_THREAD, ONLINE_CUTOFF);
  
  if (length / numThreads >= onlineCutoff) return mt_computeOnlineUnrolledVariance(threadManager, x, length, meanPtr);
  return mt_computeUnrolledVariance(threadManager, x, length, meanPtr);
}

double misc_mt_computeIndexedVariance(misc_mt_manager_t restrict threadManager, const double* restrict x, const size_t* restrict indices, size_t length, double* meanPtr)
{
  size_t numThreads = misc_mt_getNumThreads(threadManager);
  size_t onlineCutoff = minimum(INDEXED_UNROLLED_VAR_MIN_NUM_VALUES_PER_THREAD, INDEXED_ONLINE_CUTOFF);
  
  if (length / numThreads >= onlineCutoff) return mt_computeIndexedOnlineUnrolledVariance(threadManager, x, indices, length, meanPtr);
  return mt_computeIndexedUnrolledVariance(threadManager, x, indices, length, meanPtr);
}

double misc_mt_computeWeightedMean(misc_mt_manager_t restrict threadManager, const double* restrict x, size_t length, const double* restrict w, double* restrict nPtr)
{
  size_t numThreads = misc_mt_getNumThreads(threadManager);
  size_t onlineCutoff = minimum(UNROLLED_WEIGHTED_MEAN_MIN_NUM_VALUES_PER_THREAD, ONLINE_CUTOFF);
  
  if (length / numThreads >= onlineCutoff) return mt_computeOnlineUnrolledWeightedMean(threadManager, x, length, w, nPtr);
  return mt_computeUnrolledWeightedMean(threadManager, x, length, w, nPtr);
}

double misc_mt_computeIndexedWeightedMean(misc_mt_manager_t restrict threadManager, const double* restrict x, const size_t* restrict indices, size_t length, const double* restrict w, double* restrict nPtr)
{
  size_t numThreads = misc_mt_getNumThreads(threadManager);
  size_t onlineCutoff = minimum(INDEXED_UNROLLED_WEIGHTED_MEAN_MIN_NUM_VALUES_PER_THREAD, ONLINE_CUTOFF);
  
  if (length / numThreads >= onlineCutoff) return mt_computeIndexedOnlineUnrolledWeightedMean(threadManager, x, indices, length, w, nPtr);
  return mt_computeIndexedUnrolledWeightedMean(threadManager, x, indices, length, w, nPtr);
}

double misc_mt_computeWeightedVarianceForKnownMean(misc_mt_manager_t restrict threadManager, const double* restrict x, misc_size_t length, const double* restrict w, double mean)
{
  size_t numThreads = misc_mt_getNumThreads(threadManager);
  size_t onlineCutoff = minimum(UNROLLED_WEIGHTED_VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD, ONLINE_CUTOFF);
  
  if (length / numThreads >= onlineCutoff) return mt_computeOnlineUnrolledWeightedVarianceForKnownMean(threadManager, x, length, w, mean);
  return mt_computeUnrolledWeightedVarianceForKnownMean(threadManager, x, length, w, mean);
}

double misc_mt_computeIndexedWeightedVarianceForKnownMean(misc_mt_manager_t restrict threadManager, const double* restrict x, const size_t* restrict indices, size_t length, const double* restrict w, double mean)
{
  size_t numThreads = misc_mt_getNumThreads(threadManager);
  size_t onlineCutoff = minimum(INDEXED_UNROLLED_WEIGHTED_VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD, ONLINE_CUTOFF);
  
  if (length / numThreads >= onlineCutoff) return mt_computeIndexedOnlineUnrolledWeightedVarianceForKnownMean(threadManager, x, indices, length, w, mean);
  return mt_computeIndexedUnrolledWeightedVarianceForKnownMean(threadManager, x, indices, length, w, mean);
}

// work-horse implementations; indexed versions follow regular b/c implementations are highly similar

static double computeMean(const double* x, size_t length)
{
  if (length == 0.0) return 0.0;
  
  double result = 0.0;
  for (size_t i = 0; i < length; ++i) result += x[i];
  return result / (double) length;
}

static double computeIndexedMean(const double* restrict x, const size_t* restrict indices, size_t length)
{
  if (length == 0.0) return 0.0;
  
  double result = 0.0;
  for (size_t i = 0; i < length; ++i) result += x[indices[i]];
  return result / (double) length;
}

static double computeUnrolledMean_c(const double* x, size_t length)
{
  if (length == 0) return 0.0;
  
  size_t i = 0;
  size_t lengthMod5 = length % 5;
  
  double result = 0.0;
  if (lengthMod5 != 0) {
    for ( ; i < lengthMod5; ++i) result += x[i];
    if (length < 5) return result / (double) length;
  }
  
  for ( ; i < length; i += 5) {
    result += x[i] + x[i + 1] + x[i + 2] + x[i + 3] + x[i + 4];
  }
  
  return result / (double) length;
}

static double computeIndexedUnrolledMean_c(const double* restrict x, const size_t* restrict indices, size_t length)
{
  if (length == 0) return 0.0;
  
  size_t i = 0;
  size_t lengthMod5 = length % 5;
  
  double result = 0.0;
  if (lengthMod5 != 0) {
    for ( ; i < lengthMod5; ++i) result += x[indices[i]];
    if (length < 5) return result / (double) length;
  }
  
  for ( ; i < length; i += 5) {
    result += x[indices[i]] + x[indices[i + 1]] + x[indices[i + 2]] + x[indices[i + 3]] + x[indices[i + 4]];
  }
  
  return result / (double) length;
}

static double computeOnlineMean(const double* x, size_t length)
{
  if (length == 0) return 0.0;
  
  double result = x[0];
  for (size_t i = 0; i < length; ++i) result += (x[i] - result) / (double) (i + 1);
  return result;
}

static double computeIndexedOnlineMean(const double* restrict x, const size_t* restrict indices, size_t length)
{
  if (length == 0) return 0.0;
  
  double result = x[indices[0]];
  for (size_t i = 0; i < length; ++i) result += (x[indices[i]] - result) / (double) (i + 1);
  return result;
}

static double computeOnlineUnrolledMean_c(const double* x, size_t length)
{
  if (length == 0) return 0.0;
  
  size_t i = 1;
  size_t lengthMod5 = (length - 1) % 5;
  
  double result = x[0];
  if (lengthMod5++ != 0) {
    for ( ; i < lengthMod5; ++i) result += (x[i] - result) / (double) (i + 1);
    if (length < 6) return result;
  }
  
  for ( ; i < length; i += 5) {
    result += (x[i] + x[i + 1] + x[i + 2] + x[i + 3] + x[i + 4] - 5.0 * result) / (double) (i + 5);
  }
  
  return result;
}

static double computeIndexedOnlineUnrolledMean_c(const double* restrict x, const size_t* restrict indices, size_t length)
{
  if (length == 0) return 0.0;
  
  size_t i = 1;
  size_t lengthMod5 = (length - 1) % 5;
  
  double result = x[indices[0]];
  if (lengthMod5++ != 0) {
    for ( ; i < lengthMod5; ++i) result += (x[indices[i]] - result) / (double) (i + 1);
    if (length < 6) return result;
  }
  
  for ( ; i < length; i += 5) {
    result += (x[indices[i]] + x[indices[i + 1]] + x[indices[i + 2]] + x[indices[i + 3]] + x[indices[i + 4]] - 5.0 * result) / (double) (i + 5);
  }
  
  return result;
}

static double computeVarianceForKnownMean(const double* x, size_t length, double mean)
{
  if (length == 0 || isnan(mean)) return nan("");
  if (length == 1) return 0.0;
  
  double result = 0.0;
  for (size_t i = 0; i < length; ++i) result += (x[i] - mean) * (x[i] - mean);
  return result / (double) (length - 1);
}

static double computeIndexedVarianceForKnownMean(const double* restrict x, const size_t* restrict indices, size_t length, double mean)
{
  if (length == 0 || isnan(mean)) return nan("");
  if (length == 1) return 0.0;
  
  double result = 0.0;
  for (size_t i = 0; i < length; ++i) result += (x[indices[i]] - mean) * (x[indices[i]] - mean);
  return result / (double) (length - 1);
}

static double computeUnrolledVarianceForKnownMean_c(const double* x, size_t length, double mean)
{
  if (length == 0 || isnan(mean)) return nan("");
  if (length == 1) return 0.0;
  
  size_t i = 0;
  size_t lengthMod5 = length % 5;
  
  double result = 0.0;
  if (lengthMod5 != 0) {
    for ( ; i < lengthMod5; ++i) result += (x[i] - mean) * (x[i] - mean);
    if (length < 5) return result / (double) (length - 1);
  }
  
  for ( ; i < length; i += 5) {
    result += (x[i    ] - mean) * (x[i    ] - mean) +
              (x[i + 1] - mean) * (x[i + 1] - mean) +
              (x[i + 2] - mean) * (x[i + 2] - mean) +
              (x[i + 3] - mean) * (x[i + 3] - mean) +
              (x[i + 4] - mean) * (x[i + 4] - mean);
  }
  
  return result / (double) (length - 1);
}

static double computeIndexedUnrolledVarianceForKnownMean_c(const double* restrict x, const size_t* restrict indices, size_t length, double mean)
{
  if (length == 0 || isnan(mean)) return nan("");
  if (length == 1) return 0.0;
  
  size_t i = 0;
  size_t lengthMod5 = length % 5;
  
  double result = 0.0;
  if (lengthMod5 != 0) {
    for ( ; i < lengthMod5; ++i) result += (x[indices[i]] - mean) * (x[indices[i]] - mean);
    if (length < 5) return result / (double) (length - 1);
  }
  
  for ( ; i < length; i += 5) {
    result += (x[indices[i]] - mean) * (x[indices[i]] - mean) +
    (x[indices[i + 1]] - mean) * (x[indices[i + 1]] - mean) +
    (x[indices[i + 2]] - mean) * (x[indices[i + 2]] - mean) +
    (x[indices[i + 3]] - mean) * (x[indices[i + 3]] - mean) +
    (x[indices[i + 4]] - mean) * (x[indices[i + 4]] - mean);
  }
  
  return result / (double) (length - 1);
}

static double computeOnlineVarianceForKnownMean(const double* x, size_t length, double mean)
{
  if (length == 0 || isnan(mean)) return nan("");
  if (length == 1) return 0.0;
  
  double result = (x[0] - mean) * (x[0] - mean) + (x[1] - mean) * (x[1] - mean);
  for (size_t i = 2; i < length; ++i) result += ((x[i] - mean) * (x[i] - mean) - result) / (double) i;
  
  return result;
}

static double computeIndexedOnlineVarianceForKnownMean(const double* restrict x, const size_t* restrict indices, size_t length, double mean)
{
  if (length == 0 || isnan(mean)) return nan("");
  if (length == 1) return 0.0;
  
  double result = (x[indices[0]] - mean) * (x[indices[0]] - mean) + (x[indices[1]] - mean) * (x[indices[1]] - mean);
  for (size_t i = 2; i < length; ++i) result += ((x[indices[i]] - mean) * (x[indices[i]] - mean) - result) / (double) i;
  
  return result;
}

static double computeOnlineUnrolledVarianceForKnownMean_c(const double* x, size_t length, double mean)
{
  if (length == 0 || isnan(mean)) return nan("");
  if (length == 1) return 0.0;
  
  size_t i = 2;
  size_t lengthMod5 = (length - 2) % 5;
  
  double result = (x[0] - mean) * (x[0] - mean) + (x[1] - mean) * (x[1] - mean);
  if (lengthMod5 != 0) {
    for ( ; i < lengthMod5 + 2; ++i) result += ((x[i] - mean) * (x[i] - mean) - result) / (double) i;
    if (length < 7) return result;
  }
  
  for ( ; i < length; i += 5) {
    result += ((x[i] - mean) * (x[i] - mean) +
               (x[i + 1] - mean) * (x[i + 1] - mean) +
               (x[i + 2] - mean) * (x[i + 2] - mean) +
               (x[i + 3] - mean) * (x[i + 3] - mean) +
               (x[i + 4] - mean) * (x[i + 4] - mean) - 5.0 * result) / (double) (i + 4);
  }
  
  return result;
}

static double computeIndexedOnlineUnrolledVarianceForKnownMean_c(const double* restrict x, const size_t* restrict indices, size_t length, double mean)
{
  if (length == 0 || isnan(mean)) return nan("");
  if (length == 1) return 0.0;
  
  size_t i = 2;
  size_t lengthMod5 = (length - 2) % 5;
  
  double result = (x[indices[0]] - mean) * (x[indices[0]] - mean) + (x[indices[1]] - mean) * (x[indices[1]] - mean);
  if (lengthMod5 != 0) {
    for ( ; i < lengthMod5 + 2; ++i) result += ((x[indices[i]] - mean) * (x[indices[i]] - mean) - result) / (double) i;
    if (length < 7) return result;
  }
  
  for ( ; i < length; i += 5) {
    result += ((x[indices[i]] - mean) * (x[indices[i]] - mean) +
               (x[indices[i + 1]] - mean) * (x[indices[i + 1]] - mean) +
               (x[indices[i + 2]] - mean) * (x[indices[i + 2]] - mean) +
               (x[indices[i + 3]] - mean) * (x[indices[i + 3]] - mean) +
               (x[indices[i + 4]] - mean) * (x[indices[i + 4]] - mean) - 5.0 * result) / (double) (i + 4);
  }
  
  return result;
}

static double computeVariance(const double* restrict x, size_t length, double* restrict meanPtr)
{
  if (length == 0) { if (meanPtr != NULL) *meanPtr = 0.0; return nan(""); }
  if (length == 1) { if (meanPtr != NULL) *meanPtr = x[0]; return 0.0; }
  
  double mean = 0.0;
  double s_sq = 0.0;
  
  for (size_t i = 0; i < length; ++i) {
    mean += x[i];
    s_sq += x[i] * x[i];
  }
  mean /= (double) length;
  
  if (meanPtr != NULL) *meanPtr = mean;
  return (s_sq - mean * mean * (double) length) / (double) (length - 1);
}

static double computeIndexedVariance(const double* restrict x, const size_t* restrict indices, size_t length, double* restrict meanPtr)
{
  if (length == 0) { if (meanPtr != NULL) *meanPtr = 0.0; return nan(""); }
  if (length == 1) { if (meanPtr != NULL) *meanPtr = x[indices[0]]; return 0.0; }
  
  double mean = 0.0;
  double s_sq = 0.0;
  
  for (size_t i = 0; i < length; ++i) {
    mean += x[indices[i]];
    s_sq += x[indices[i]] * x[indices[i]];
  }
  mean /= (double) length;
  
  if (meanPtr != NULL)  *meanPtr = mean;
  return (s_sq - mean * mean * (double) length) / (double) (length - 1);
}

static double computeUnrolledVariance(const double* restrict x, size_t length, double* restrict meanPtr)
{
  if (length == 0) { if (meanPtr != NULL) *meanPtr = 0.0; return nan(""); }
  if (length == 1) { if (meanPtr != NULL) *meanPtr = x[0]; return 0.0; }
  
  double mean = 0.0;
  double x_sq = 0.0;
  
  size_t i = 0;
  size_t lengthMod5 = length % 5;
  
  if (lengthMod5 != 0) {
    for ( ; i < lengthMod5; ++i) {
      mean += x[i];
      x_sq += x[i] * x[i];
    }
    if (length < 5) { mean /= (double) length; if (meanPtr != NULL) *meanPtr = mean; return (x_sq - mean * mean * (double) length) / (double) (length - 1); }
  }
  
  for ( ; i < length; i += 5) {
    mean += x[i] + x[i + 1] + x[i + 2] + x[i + 3] + x[i + 4];
    x_sq += x[i] * x[i] +
    x[i + 1] * x[i + 1] +
    x[i + 2] * x[i + 2] +
    x[i + 3] * x[i + 3] +
    x[i + 4] * x[i + 4];
  }
  mean /= (double) length;
  
  if (meanPtr != NULL) *meanPtr = mean;
  return (x_sq - mean * mean * (double) length) / (double) (length - 1);
}

static double computeIndexedUnrolledVariance(const double* restrict x, const size_t* restrict indices, size_t length, double* restrict meanPtr)
{
  if (length == 0) { if (meanPtr != NULL) *meanPtr = 0.0; return nan(""); }
  if (length == 1) { if (meanPtr != NULL) *meanPtr = x[indices[0]]; return 0.0; }
  
  double mean = 0.0;
  double x_sq = 0.0;
  
  size_t i = 0;
  size_t lengthMod5 = length % 5;
  
  if (lengthMod5 != 0) {
    for ( ; i < lengthMod5; ++i) {
      mean += x[indices[i]];
      x_sq += x[indices[i]] * x[indices[i]];
    }
    if (length < 5) { mean /= (double) length; if (meanPtr != NULL) *meanPtr = mean; return (x_sq - mean * mean * (double) length) / (double) (length - 1); }
  }
  
  for ( ; i < length; i += 5) {
    mean += x[indices[i]] + x[indices[i + 1]] + x[indices[i + 2]] + x[indices[i + 3]] + x[indices[i + 4]];
    x_sq += (x[indices[i]] * x[indices[i]] +
             x[indices[i + 1]] * x[indices[i + 1]] +
             x[indices[i + 2]] * x[indices[i + 2]] +
             x[indices[i + 3]] * x[indices[i + 3]] +
             x[indices[i + 4]] * x[indices[i + 4]]);
  }
  mean /= (double) length;
  
  if (meanPtr != NULL) *meanPtr = mean;
  return (x_sq - mean * mean * (double) length) / (double) (length - 1);
}

static double computeOnlineVariance(const double* restrict x, size_t length, double* restrict meanPtr)
{
  if (length == 0) { if (meanPtr != NULL) *meanPtr = 0.0; return nan(""); }
  if (length == 1) { if (meanPtr != NULL) *meanPtr = x[0]; return 0.0; }
  
  double mean = x[0];
  double var  = 0.0;
  
  for (size_t i = 1; i < length; ++i) {
    double dev = x[i] - mean;
    mean += dev / (double) (i + 1);
    var  += (dev * (x[i] - mean) - var) / (double) i;
  }
  
  if (meanPtr != NULL) *meanPtr = mean;
  return var;
}

static double computeIndexedOnlineVariance(const double* restrict x, const size_t* restrict indices, size_t length, double* restrict meanPtr)
{
  if (length == 0) { if (meanPtr != NULL) *meanPtr = 0.0; return nan(""); }
  if (length == 1) { if (meanPtr != NULL) *meanPtr = x[indices[0]]; return 0.0; }
  
  double mean = x[indices[0]];
  double var  = 0.0;
  
  for (size_t i = 1; i < length; ++i) {
    double dev = x[indices[i]] - mean;
    mean += dev / (double) (i + 1);
    var  += (dev * (x[indices[i]] - mean) - var) / (double) i;
  }
  
  if (meanPtr != NULL) *meanPtr = mean;
  return var;
}

// for this and this alone, "variance" is divided by n
// we rescale before returning, however
static double computeOnlineUnrolledVariance(const double* restrict x, size_t length, double* restrict meanPtr)
{
  if (length == 0) { if (meanPtr != NULL) *meanPtr = 0.0; return nan(""); }
  if (length == 1) { if (meanPtr != NULL) *meanPtr = x[0]; return 0.0; }
  
  double mean = x[0];
  double var  = 0.0;
  double nScale = (double) length / (double) (length - 1);
  
  size_t i = 1;
  size_t lengthMod5 = (length - 1) % 5;
  
  if (lengthMod5++ != 0) {
    for ( ; i < lengthMod5; ++i) {
      double dev = x[i] - mean;
      mean += dev / (double) (i + 1);
      var  += (dev * (x[i] - mean) - var) / (double) (i + 1);
    }
    if (length < 6) { if (meanPtr != NULL) *meanPtr = mean; return nScale * var; }
  }
  
  for ( ; i < length; i += 5) {
#define n1 ((double) i)
#define n2 5.0
#define n  ((double) (i + 5))
    
    double meanNext5 = (x[i] + x[i + 1] + x[i + 2] + x[i + 3] + x[i + 4]) / n2;
    double varNext5  = ((x[i] - meanNext5) * (x[i] - meanNext5) + (x[i + 1] - meanNext5) * (x[i + 1] - meanNext5) +
                        (x[i + 2] - meanNext5) * (x[i + 2] - meanNext5) + (x[i + 3] - meanNext5) * (x[i + 3] - meanNext5) + (x[i + 4] - meanNext5) * (x[i + 4] - meanNext5)) / n2;
    
    var  += (n2 / n) * (varNext5 - var) + ((meanNext5 - mean) * (n1 / n)) * ((meanNext5 - mean) * (n2 / n));
    mean += (n2 / n) * (meanNext5 - mean);
#undef n1
#undef n2
#undef n
  }
  
  if (meanPtr != NULL) *meanPtr = mean;
  return nScale * var;
}

static double computeIndexedOnlineUnrolledVariance(const double* restrict x, const size_t* restrict indices, size_t length, double* restrict meanPtr)
{
  if (length == 0) { if (meanPtr != NULL) *meanPtr = 0.0; return nan(""); }
  if (length == 1) { if (meanPtr != NULL) *meanPtr = x[indices[0]]; return 0.0; }
  
  double mean = x[indices[0]];
  double var  = 0.0;
  double nScale = (double) length / (double) (length - 1);
  
  size_t i = 1;
  size_t lengthMod5 = (length - 1) % 5;
  
  if (lengthMod5++ != 0) {
    for ( ; i < lengthMod5; ++i) {
      double dev = x[indices[i]] - mean;
      mean += dev / (double) (i + 1);
      var  += (dev * (x[indices[i]] - mean) - var) / (double) (i + 1);
    }
    if (length < 6) { if (meanPtr != NULL) *meanPtr = mean; return nScale * var; }
  }
  
  for ( ; i < length; i += 5) {
#define n1 ((double) i)
#define n2 5.0
#define n  ((double) (i + 5))
    
    double meanNext5 = (x[indices[i]] + x[indices[i + 1]] + x[indices[i + 2]] + x[indices[i + 3]] + x[indices[i + 4]]) / n2;
    double varNext5  = ((x[indices[i]] - meanNext5) * (x[indices[i]] - meanNext5) +
                        (x[indices[i + 1]] - meanNext5) * (x[indices[i + 1]] - meanNext5) +
                        (x[indices[i + 2]] - meanNext5) * (x[indices[i + 2]] - meanNext5) +
                        (x[indices[i + 3]] - meanNext5) * (x[indices[i + 3]] - meanNext5) +
                        (x[indices[i + 4]] - meanNext5) * (x[indices[i + 4]] - meanNext5)) / n2;
    
    var  += (n2 / n) * (varNext5 - var) + ((meanNext5 - mean) * (n1 / n)) * ((meanNext5 - mean) * (n2 / n));
    mean += (n2 / n) * (meanNext5 - mean);
#undef n1
#undef n2
#undef n
  }
  
  if (meanPtr != NULL) *meanPtr = mean;
  return nScale * var;
}

static double computeWeightedMean(const double* restrict x, size_t length, const double* restrict w, double* restrict nPtr)
{
  double n = 0.0;
  double result = 0.0;
  for (size_t i = 0; i < length; ++i) {
    result += x[i] * w[i];
    n += w[i];
  }
  if (nPtr != NULL) *nPtr = n;
  return result / n;
}

static double computeIndexedWeightedMean(const double* restrict x, const size_t* restrict indices, size_t length, const double* restrict w, double* restrict nPtr)
{
  double n = 0.0;
  double result = 0.0;
  for (size_t i = 0; i < length; ++i) {
    result += x[indices[i]] * w[indices[i]];
    n += w[indices[i]];
  }
  if (nPtr != NULL) *nPtr = n;
  return result / n;
}

static double computeUnrolledWeightedMean_c(const double* restrict x, size_t length, const double* restrict w, double* restrict nPtr)
{
  if (length == 0) { if (nPtr != NULL) *nPtr = 0.0; return 0.0; }
  
  size_t i = 0;
  size_t lengthMod5 = length % 5;
  
  double result = 0.0;
  double n = 0.0;
  if (lengthMod5 != 0) {
    for ( ; i < lengthMod5; ++i) { result += x[i] * w[i]; n += w[i]; }
    if (length < 5) { if (nPtr != NULL) *nPtr = n; return result / n; }
  }
  
  for ( ; i < length; i += 5) {
    result += x[i] * w[i] + x[i + 1] * w[i + 1] + x[i + 2] * w[i + 2] + x[i + 3] * w[i + 3] + x[i + 4] * w[i + 4];
    n += w[i] + w[i + 1] + w[i + 2] + w[i + 3] + w[i + 4];
  }
  
  if (nPtr != NULL) *nPtr = n;
  return result / n;
}

static double computeIndexedUnrolledWeightedMean_c(const double* restrict x, const size_t* restrict indices, size_t length, const double* restrict w, double* restrict nPtr)
{
  if (length == 0) { if (nPtr != NULL) *nPtr = 0.0; return 0.0; }
  
  size_t i = 0;
  size_t lengthMod5 = length % 5;
  
  double result = 0.0;
  double n = 0.0;
  if (lengthMod5 != 0) {
    for ( ; i < lengthMod5; ++i) { result += x[indices[i]] * w[indices[i]]; n += w[indices[i]]; }
    if (length < 5) { if (nPtr != NULL) *nPtr = n; return result / n; }
  }
  
  for ( ; i < length; i += 5) {
    result += x[indices[i]] * w[indices[i]] + x[indices[i + 1]] * w[indices[i + 1]] + x[indices[i + 2]] * w[indices[i + 2]] +
              x[indices[i + 3]] * w[indices[i + 3]] + x[indices[i + 4]] * w[indices[i + 4]];
    n += w[indices[i]] + w[indices[i + 1]] + w[indices[i + 2]] + w[indices[i + 3]] + w[indices[i + 4]];
  }
  
  if (nPtr != NULL) *nPtr = n;
  return result / n;
}

static double computeOnlineWeightedMean(const double* restrict x, size_t length, const double* restrict w, double* restrict nPtr)
{
  if (length == 0) { if (nPtr != NULL) *nPtr = 0.0; return 0.0; }
  
  double n = w[0];
  double result = x[0];
  for (size_t i = 1; i < length; ++i) {
    n += w[i];
    result += (x[i] * w[i] - w[i] * result) / n;
  }
  
  if (nPtr != NULL) *nPtr = n;
  return result;
}

static double computeIndexedOnlineWeightedMean(const double* restrict x, const size_t* restrict indices, size_t length, const double* restrict w, double* restrict nPtr)
{
  if (length == 0) { if (nPtr != NULL) *nPtr = 0.0; return 0.0; }
  
  double n = w[indices[0]];
  double result = x[indices[0]];
  for (size_t i = 1; i < length; ++i) {
    n += w[indices[i]];
    result += (x[indices[i]] * w[indices[i]] - w[indices[i]] * result) / n;
  }
  
  if (nPtr != NULL) *nPtr = n;
  return result;
}

static double computeOnlineUnrolledWeightedMean_c(const double* restrict x, size_t length, const double* restrict w, double* restrict nPtr)
{
  if (length == 0) { if (nPtr != NULL) *nPtr = 0.0; return 0.0; }
  
  size_t i = 1;
  size_t lengthMod5 = (length - 1) % 5;
  
  double n = w[0];
  double result = x[0];
  if (lengthMod5++ != 0) {
    for ( ; i < lengthMod5; ++i) {
      n += w[i];
      result += (x[i] - result) * (w[i] / n);
    }
    if (length < 6) {
      if (nPtr != NULL) *nPtr = n;
      return result;
    }
  }
  
  for ( ; i < length; i += 5) {
    double delta_n = w[i] + w[i + 1] + w[i + 2] + w[i + 3] + w[i + 4];
    n += delta_n;
    result += (x[i] * w[i] + x[i + 1] * w[i + 1] + x[i + 2] * w[i + 2] + x[i + 3] * w[i + 3] + x[i + 4] * w[i + 4] -
               delta_n * result) / n;
  }
  
  if (nPtr != NULL) *nPtr = n;
  return result;
}

static double computeIndexedOnlineUnrolledWeightedMean_c(const double* restrict x, const size_t* restrict indices, size_t length, const double* restrict w, double* restrict nPtr)
{
  if (length == 0) { if (nPtr != NULL) *nPtr = 0.0; return 0.0; }
  
  size_t i = 1;
  size_t lengthMod5 = (length - 1) % 5;
  
  double n = w[indices[0]];
  double result = x[indices[0]];
  if (lengthMod5++ != 0) {
    for ( ; i < lengthMod5; ++i) {
      n += w[indices[i]];
      result += (x[indices[i]] - result) * (w[indices[i]] / n);
    }
    if (length < 6) {
      if (nPtr != NULL) *nPtr = n;
      return result;
    }
  }
  
  for ( ; i < length; i += 5) {
    double delta_n = w[indices[i]] + w[indices[i + 1]] + w[indices[i + 2]] + w[indices[i + 3]] + w[indices[i + 4]];
    n += delta_n;
    result += (x[indices[i]] * w[indices[i]] + x[indices[i + 1]] * w[indices[i + 1]] +
               x[indices[i + 2]] * w[indices[i + 2]] + x[indices[i + 3]] * w[indices[i + 3]] +
               x[indices[i + 4]] * w[indices[i + 4]] - delta_n * result) / n;
  }
  
  if (nPtr != NULL) *nPtr = n;
  return result;
}

static double computeWeightedVarianceForKnownMean(const double* restrict x, size_t length, const double* restrict w, double mean)
{
  if (length == 0 || isnan(mean)) return nan("");
  if (length == 1) return 0.0;
  
  double result = 0.0;
  for (size_t i = 0; i < length; ++i) result += w[i] * (x[i] - mean) * (x[i] - mean);
  return result / (double) (length - 1);
}

static double computeIndexedWeightedVarianceForKnownMean(const double* restrict x, const size_t* restrict indices, size_t length, const double* restrict w, double mean)
{
  if (length == 0 || isnan(mean)) return nan("");
  if (length == 1) return 0.0;
  
  double result = 0.0;
  for (size_t i = 0; i < length; ++i) result += w[indices[i]] * (x[indices[i]] - mean) * (x[indices[i]] - mean);
  return result / (double) (length - 1);
}

static double computeUnrolledWeightedVarianceForKnownMean_c(const double* restrict x, size_t length, const double* restrict w, double mean)
{
  if (length == 0 || isnan(mean)) return nan("");
  if (length == 1) return 0.0;
  
  size_t i = 0;
  size_t lengthMod5 = length % 5;
  
  double result = 0.0;
  if (lengthMod5 != 0) {
    for ( ; i < lengthMod5; ++i) result += w[i] * (x[i] - mean) * (x[i] - mean);
    if (length < 5) return result / (double) (length - 1);
  }
  
  for ( ; i < length; i += 5) {
    result += w[i] * (x[i] - mean) * (x[i] - mean) +
              w[i + 1] * (x[i + 1] - mean) * (x[i + 1] - mean) +
              w[i + 2] * (x[i + 2] - mean) * (x[i + 2] - mean) +
              w[i + 3] * (x[i + 3] - mean) * (x[i + 3] - mean) +
              w[i + 4] * (x[i + 4] - mean) * (x[i + 4] - mean);
  }
  
  return result / (double) (length - 1);
}

static double computeIndexedUnrolledWeightedVarianceForKnownMean_c(const double* restrict x, const size_t* restrict indices, size_t length, const double* restrict w, double mean)
{
  if (length == 0 || isnan(mean)) return nan("");
  if (length == 1) return 0.0;
  
  size_t i = 0;
  size_t lengthMod5 = length % 5;
  
  double result = 0.0;
  if (lengthMod5 != 0) {
    for ( ; i < lengthMod5; ++i) result += w[indices[i]] * (x[indices[i]] - mean) * (x[indices[i]] - mean);
    if (length < 5) return result / (double) (length - 1);
  }
  
  for ( ; i < length; i += 5) {
    result += w[indices[i]] * (x[indices[i]] - mean) * (x[indices[i]] - mean) +
              w[indices[i + 1]] * (x[indices[i + 1]] - mean) * (x[indices[i + 1]] - mean) +
              w[indices[i + 2]] * (x[indices[i + 2]] - mean) * (x[indices[i + 2]] - mean) +
              w[indices[i + 3]] * (x[indices[i + 3]] - mean) * (x[indices[i + 3]] - mean) +
              w[indices[i + 4]] * (x[indices[i + 4]] - mean) * (x[indices[i + 4]] - mean);
  }
  
  return result / (double) (length - 1);
}

static double computeOnlineWeightedVarianceForKnownMean(const double* restrict x, size_t length, const double* restrict w, double mean)
{
  if (length == 0 || isnan(mean)) return nan("");
  if (length == 1) return 0.0;
  
  double result = w[0] * (x[0] - mean) * (x[0] - mean) + w[1] * (x[1] - mean) * (x[1] - mean);
  for (size_t i = 2; i < length; ++i) {
    result += (w[i] * (x[i] - mean) * (x[i] - mean) - result) / (double) i;
  }
  
  return result;
}

static double computeIndexedOnlineWeightedVarianceForKnownMean(const double* restrict x, const size_t* restrict indices, size_t length, const double* restrict w, double mean)
{
  if (length == 0 || isnan(mean)) return nan("");
  if (length == 1) return 0.0;
  
  double result = w[indices[0]] * (x[indices[0]] - mean) * (x[indices[0]] - mean) + w[indices[1]] * (x[indices[1]] - mean) * (x[indices[1]] - mean);
  for (size_t i = 2; i < length; ++i) {
    result += (w[indices[i]] * (x[indices[i]] - mean) * (x[indices[i]] - mean) - result) / (double) i;
  }
  
  return result;
}

static double computeOnlineUnrolledWeightedVarianceForKnownMean_c(const double* restrict x, size_t length, const double* restrict w, double mean)
{
  if (length == 0 || isnan(mean)) return nan("");
  if (length == 1) return 0.0;
  
  size_t i = 2;
  size_t lengthMod5 = (length - 2) % 5;
  
  double result = w[0] * (x[0] - mean) * (x[0] - mean) + w[1] * (x[1] - mean) * (x[1] - mean);
  if (lengthMod5 != 0) {
    for ( ; i < lengthMod5 + 2; ++i) result += (w[i] * (x[i] - mean) * (x[i] - mean) - result) / (double) i;
    if (length < 7) return result;
  }
  
  for ( ; i < length; i += 5) {
    result += (w[i] * (x[i] - mean) * (x[i] - mean) +
               w[i + 1] * (x[i + 1] - mean) * (x[i + 1] - mean) +
               w[i + 2] * (x[i + 2] - mean) * (x[i + 2] - mean) +
               w[i + 3] * (x[i + 3] - mean) * (x[i + 3] - mean) +
               w[i + 4] * (x[i + 4] - mean) * (x[i + 4] - mean) - 5.0 * result) / (double) (i + 4);
  }
  
  return result;
}

static double computeIndexedOnlineUnrolledWeightedVarianceForKnownMean_c(const double* restrict x, const size_t* restrict indices, size_t length, const double* restrict w, double mean)
{
  if (length == 0 || isnan(mean)) return nan("");
  if (length == 1) return 0.0;
  
  size_t i = 2;
  size_t lengthMod5 = (length - 2) % 5;
  
  double result = w[indices[0]] * (x[indices[0]] - mean) * (x[indices[0]] - mean) + w[indices[1]] * (x[indices[1]] - mean) * (x[indices[1]] - mean);
  if (lengthMod5 != 0) {
    for ( ; i < lengthMod5 + 2; ++i) result += (w[indices[i]] * (x[indices[i]] - mean) * (x[indices[i]] - mean) - result) / (double) i;
    if (length < 7) return result;
  }
  
  for ( ; i < length; i += 5) {
    result += (w[indices[i]] * (x[indices[i]] - mean) * (x[indices[i]] - mean) +
               w[indices[i + 1]] * (x[indices[i + 1]] - mean) * (x[indices[i + 1]] - mean) +
               w[indices[i + 2]] * (x[indices[i + 2]] - mean) * (x[indices[i + 2]] - mean) +
               w[indices[i + 3]] * (x[indices[i + 3]] - mean) * (x[indices[i + 3]] - mean) +
               w[indices[i + 4]] * (x[indices[i + 4]] - mean) * (x[indices[i + 4]] - mean) - 5.0 * result) / (double) (i + 4);
  }
  
  return result;
}

#include <misc/simd.h>

#ifdef HAVE_SSE2
extern double misc_computeUnrolledMean_sse2(const double* x, size_t length);
extern double misc_computeIndexedUnrolledMean_sse2(const double* restrict x, const size_t* restrict indices, size_t length);
extern double misc_computeOnlineUnrolledMean_sse2(const double* x, size_t length);
extern double misc_computeIndexedOnlineUnrolledMean_sse2(const double* restrict x, const size_t* restrict indices, size_t length);
extern double misc_computeUnrolledVarianceForKnownMean_sse2(const double* x, size_t length, double mean);
extern double misc_computeIndexedUnrolledVarianceForKnownMean_sse2(const double* restrict x, const size_t* restrict indices, size_t length, double mean);
extern double misc_computeOnlineUnrolledVarianceForKnownMean_sse2(const double* x, size_t length, double mean);
extern double misc_computeIndexedOnlineUnrolledVarianceForKnownMean_sse2(const double* restrict x, const size_t* restrict indices, size_t length, double mean);
extern double misc_computeUnrolledWeightedMean_sse2(const double* restrict x, size_t length, const double* restrict w, double* restrict nPtr);
extern double misc_computeIndexedUnrolledWeightedMean_sse2(const double* restrict x, const size_t* restrict indices, size_t length, const double* restrict w, double* restrict nPtr);
extern double misc_computeOnlineUnrolledWeightedMean_sse2(const double* restrict x, size_t length, const double* restrict w, double* restrict nPtr);
extern double misc_computeIndexedOnlineUnrolledWeightedMean_sse2(const double* restrict x, const size_t* restrict indices, size_t length, const double* restrict w, double* restrict nPtr);
extern double misc_computeUnrolledWeightedVarianceForKnownMean_sse2(const double* restrict x, size_t length, const double* restrict w, double mean);
extern double misc_computeIndexedUnrolledWeightedVarianceForKnownMean_sse2(const double* restrict x, const size_t* restrict indices, size_t length, const double* restrict w, double mean);
extern double misc_computeOnlineUnrolledWeightedVarianceForKnownMean_sse2(const double* restrict x, size_t length, const double* restrict w, double mean);
extern double misc_computeIndexedOnlineUnrolledWeightedVarianceForKnownMean_sse2(const double* restrict x, const size_t* restrict indices, size_t length, const double* restrict w, double mean);
#endif

void misc_stat_setSIMDInstructionSet(misc_simd_instructionSet i)
{
#ifdef HAVE_SSE2
  if (i >= MISC_INST_SSE2) {
    computeUnrolledMean = &misc_computeUnrolledMean_sse2;
    computeOnlineUnrolledMean = &misc_computeOnlineUnrolledMean_sse2;
    computeIndexedUnrolledMean = &misc_computeIndexedUnrolledMean_sse2;
    computeIndexedOnlineUnrolledMean = &misc_computeIndexedOnlineUnrolledMean_sse2;
    computeUnrolledWeightedMean = &misc_computeUnrolledWeightedMean_sse2;
    computeIndexedUnrolledWeightedMean = &misc_computeIndexedUnrolledWeightedMean_sse2;
    computeOnlineUnrolledWeightedMean = &misc_computeOnlineUnrolledWeightedMean_sse2;
    computeIndexedOnlineUnrolledWeightedMean = &misc_computeIndexedOnlineUnrolledWeightedMean_sse2;
    
    computeUnrolledVarianceForKnownMean = &misc_computeUnrolledVarianceForKnownMean_sse2;
    computeIndexedUnrolledVarianceForKnownMean = &misc_computeIndexedUnrolledVarianceForKnownMean_sse2;
    computeOnlineUnrolledVarianceForKnownMean = &misc_computeOnlineUnrolledVarianceForKnownMean_sse2;
    computeIndexedOnlineUnrolledVarianceForKnownMean = &misc_computeIndexedOnlineUnrolledVarianceForKnownMean_sse2;
    computeUnrolledWeightedVarianceForKnownMean = &misc_computeUnrolledWeightedVarianceForKnownMean_sse2;
    computeIndexedUnrolledWeightedVarianceForKnownMean = &misc_computeIndexedUnrolledWeightedVarianceForKnownMean_sse2;
    computeOnlineUnrolledWeightedVarianceForKnownMean = &misc_computeOnlineUnrolledWeightedVarianceForKnownMean_sse2;
    computeIndexedOnlineUnrolledWeightedVarianceForKnownMean = &misc_computeIndexedOnlineUnrolledWeightedVarianceForKnownMean_sse2;
  } else
#endif
  {
    computeUnrolledMean = &computeUnrolledMean_c;
    computeOnlineUnrolledMean = &computeOnlineUnrolledMean_c;
    computeIndexedUnrolledMean = &computeIndexedUnrolledMean_c;
    computeIndexedOnlineUnrolledMean = &computeIndexedOnlineUnrolledMean_c;
    computeUnrolledWeightedMean = &computeUnrolledWeightedMean_c;
    computeIndexedUnrolledWeightedMean = &computeIndexedUnrolledWeightedMean_c;
    computeOnlineUnrolledWeightedMean = &computeOnlineUnrolledWeightedMean_c;
    computeIndexedOnlineUnrolledWeightedMean = &computeIndexedOnlineUnrolledWeightedMean_c;
    
    computeUnrolledVarianceForKnownMean = &computeUnrolledVarianceForKnownMean_c;
    computeIndexedUnrolledVarianceForKnownMean = &computeIndexedUnrolledVarianceForKnownMean_c;
    computeOnlineUnrolledVarianceForKnownMean = &computeOnlineUnrolledVarianceForKnownMean_c;
    computeIndexedOnlineUnrolledVarianceForKnownMean = &computeIndexedOnlineUnrolledVarianceForKnownMean_c;
    computeUnrolledWeightedVarianceForKnownMean = &computeUnrolledWeightedVarianceForKnownMean_c;
    computeIndexedUnrolledWeightedVarianceForKnownMean = &computeIndexedUnrolledWeightedVarianceForKnownMean_c;
    computeOnlineUnrolledWeightedVarianceForKnownMean = &computeOnlineUnrolledWeightedVarianceForKnownMean_c;
    computeIndexedOnlineUnrolledWeightedVarianceForKnownMean = &computeIndexedOnlineUnrolledWeightedVarianceForKnownMean_c;
  }

}

// below this, multithreaded madness

typedef double (*meanFunction)(const double* restrict x, size_t length);
typedef double (*indexedMeanFunction)(const double* restrict x, const size_t* restrict indices, size_t length);

typedef struct {
  const double* x;
  size_t length;
  double result;
  meanFunction function;
} MeanData;

typedef struct {
  const double* x;
  const size_t* indices;
  size_t length;
  double result;
  indexedMeanFunction function;
} IndexedMeanData;

static void computeMeanTask(void* v_data)
{
  MeanData* data = (MeanData*) v_data;
  data->result = data->function(data->x, data->length);
}

static void computeIndexedMeanTask(void* v_data)
{
  IndexedMeanData* data = (IndexedMeanData*) v_data;
  data->result = data->function(data->x, data->indices, data->length);
}

static void setupMeanData(MeanData* restrict threadData, size_t numThreads, const double* restrict x,
                          size_t numValuesPerThread, size_t offByOneIndex, meanFunction function)
{
  size_t i = 0;
  for ( ; i < offByOneIndex; ++i) {
    threadData[i].x = x + i * numValuesPerThread;
    threadData[i].length = numValuesPerThread;
    threadData[i].function = function;
  }
  for ( ; i < numThreads; ++i) {
    threadData[i].x = x + offByOneIndex * numValuesPerThread + (i - offByOneIndex) * (numValuesPerThread - 1);
    threadData[i].length = numValuesPerThread - 1;
    threadData[i].function = function;
  }
}

static void setupIndexedMeanData(IndexedMeanData* restrict threadData, size_t numThreads, const double* restrict x,
                                 const size_t* restrict indices, size_t numValuesPerThread, size_t offByOneIndex,
                                 indexedMeanFunction function)
{
  size_t i = 0;
  for ( ; i < offByOneIndex; ++i) {
    threadData[i].x = x;
    threadData[i].indices = indices + i * numValuesPerThread;
    threadData[i].length = numValuesPerThread;
    threadData[i].function = function;
  }
  for ( ; i < numThreads; ++i) {
    threadData[i].x = x;
    threadData[i].indices = indices + offByOneIndex * numValuesPerThread + (i - offByOneIndex) * (numValuesPerThread - 1);
    threadData[i].length = numValuesPerThread - 1;
    threadData[i].function = function;
  }
}

static double aggregateMeanResults(const MeanData* threadData, size_t numThreads)
{
  double result = threadData[0].result;
  size_t setSize = threadData[0].length;
  for (size_t i = 1; i < numThreads; ++i) {
    setSize += threadData[i].length;
    result += ((double) threadData[i].length / (double) setSize) * (threadData[i].result - result);
  }
  return result;
}

static double aggregateIndexedMeanResults(const IndexedMeanData* threadData, size_t numThreads)
{
  double result = threadData[0].result;
  size_t setSize = threadData[0].length;
  for (size_t i = 1; i < numThreads; ++i) {
    setSize += threadData[i].length;
    result += ((double) threadData[i].length / (double) setSize) * (threadData[i].result - result);
  }
  return result;
}

/* static double mt_computeMean(misc_mt_manager_t restrict threadManager, const double* restrict x, size_t length)
{
  size_t numThreads, numValuesPerThread, offByOneIndex;
  misc_mt_getNumThreadsForJob(threadManager, length, MEAN_MIN_NUM_VALUES_PER_THREAD,
                             &numThreads, &numValuesPerThread, &offByOneIndex);
  
  if (numThreads <= 1) return computeMean(x, length);
  
  MeanData threadData[numThreads];
  setupMeanData(threadData, numThreads, x, numValuesPerThread, offByOneIndex, computeMean);
  
  void* threadDataPtrs[numThreads];
  for (size_t i = 0; i < numThreads; ++i) threadDataPtrs[i] = (void*) &threadData[i];
  
  
  
  misc_mt_runTasks(threadManager, &computeMeanTask, threadDataPtrs, numThreads);
  
  
  
  return aggregateMeanResults(threadData, numThreads);
}

static double mt_computeIndexedMean(misc_mt_manager_t restrict threadManager, const double* restrict x, const size_t* restrict indices, size_t length)
{
  size_t numThreads, numValuesPerThread, offByOneIndex;
  misc_mt_getNumThreadsForJob(threadManager, length, INDEXED_MEAN_MIN_NUM_VALUES_PER_THREAD,
                             &numThreads, &numValuesPerThread, &offByOneIndex);
  
  if (numThreads <= 1) return computeIndexedMean(x, indices, length);
  
  IndexedMeanData threadData[numThreads];
  setupIndexedMeanData(threadData, numThreads, x, indices, numValuesPerThread, offByOneIndex, computeIndexedMean);
  
  void* threadDataPtrs[numThreads];
  for (size_t i = 0; i < numThreads; ++i) threadDataPtrs[i] = (void*) &threadData[i];
  
  
  
  misc_mt_runTasks(threadManager, &computeIndexedMeanTask, threadDataPtrs, numThreads);
  
  
  
  return aggregateIndexedMeanResults(threadData, numThreads);
} */

static double mt_computeUnrolledMean(misc_mt_manager_t restrict threadManager, const double* restrict x, size_t length)
{
  size_t numThreads, numValuesPerThread, offByOneIndex;
  misc_mt_getNumThreadsForJob(threadManager, length, UNROLLED_MEAN_MIN_NUM_VALUES_PER_THREAD,
                           &numThreads, &numValuesPerThread, &offByOneIndex);
  
  if (numThreads <= 1) return computeUnrolledMean(x, length);
  
  MeanData threadData[numThreads];
  setupMeanData(threadData, numThreads, x, numValuesPerThread, offByOneIndex, computeUnrolledMean);
  
  void* threadDataPtrs[numThreads];
  for (size_t i = 0; i < numThreads; ++i) threadDataPtrs[i] = (void*) &threadData[i];
  
  
  
  misc_mt_runTasks(threadManager, &computeMeanTask, threadDataPtrs, numThreads);
  
  
  
  return aggregateMeanResults(threadData, numThreads);
}

static double mt_computeIndexedUnrolledMean(misc_mt_manager_t restrict threadManager, const double* restrict x, const size_t* restrict indices, size_t length)
{
  size_t numThreads, numValuesPerThread, offByOneIndex;
  misc_mt_getNumThreadsForJob(threadManager, length, INDEXED_UNROLLED_MEAN_MIN_NUM_VALUES_PER_THREAD,
                             &numThreads, &numValuesPerThread, &offByOneIndex);
  
  if (numThreads <= 1) return computeIndexedUnrolledMean(x, indices, length);
  
  IndexedMeanData threadData[numThreads];
  setupIndexedMeanData(threadData, numThreads, x, indices, numValuesPerThread, offByOneIndex, computeIndexedUnrolledMean);
  
  void* threadDataPtrs[numThreads];
  for (size_t i = 0; i < numThreads; ++i) threadDataPtrs[i] = (void*) &threadData[i];
  
  
  
  misc_mt_runTasks(threadManager, &computeIndexedMeanTask, threadDataPtrs, numThreads);
  
  
  
  return aggregateIndexedMeanResults(threadData, numThreads);
}

/* static double mt_computeOnlineMean(misc_mt_manager_t restrict threadManager, const double* restrict x, size_t length)
{
  size_t numThreads, numValuesPerThread, offByOneIndex;
  misc_mt_getNumThreadsForJob(threadManager, length, ONLINE_MEAN_MIN_NUM_VALUES_PER_THREAD,
                           &numThreads, &numValuesPerThread, &offByOneIndex);
  
  if (numThreads <= 1) return computeOnlineMean(x, length);
  
  MeanData threadData[numThreads];
  setupMeanData(threadData, numThreads, x, numValuesPerThread, offByOneIndex, computeOnlineMean);
  
  void* threadDataPtrs[numThreads];
  for (size_t i = 0; i < numThreads; ++i) threadDataPtrs[i] = (void*) &threadData[i];
  
  
  
  misc_mt_runTasks(threadManager, &computeMeanTask, threadDataPtrs, numThreads);
  
  
  
  return aggregateMeanResults(threadData, numThreads);
}

static double mt_computeIndexedOnlineMean(misc_mt_manager_t restrict threadManager, const double* restrict x, const size_t* restrict indices, size_t length)
{
  size_t numThreads, numValuesPerThread, offByOneIndex;
  misc_mt_getNumThreadsForJob(threadManager, length, INDEXED_ONLINE_MEAN_MIN_NUM_VALUES_PER_THREAD,
                             &numThreads, &numValuesPerThread, &offByOneIndex);
  
  if (numThreads <= 1) return computeIndexedOnlineMean(x, indices, length);
  
  IndexedMeanData threadData[numThreads];
  setupIndexedMeanData(threadData, numThreads, x, indices, numValuesPerThread, offByOneIndex, computeIndexedOnlineMean);
  
  void* threadDataPtrs[numThreads];
  for (size_t i = 0; i < numThreads; ++i) threadDataPtrs[i] = (void*) &threadData[i];
  
  
  
  misc_mt_runTasks(threadManager, &computeIndexedMeanTask, threadDataPtrs, numThreads);
  
  
  
  return aggregateIndexedMeanResults(threadData, numThreads);
} */

static double mt_computeOnlineUnrolledMean(misc_mt_manager_t restrict threadManager, const double* restrict x, size_t length)
{
  size_t numThreads, numValuesPerThread, offByOneIndex;
  misc_mt_getNumThreadsForJob(threadManager, length, ONLINE_UNROLLED_MEAN_MIN_NUM_VALUES_PER_THREAD,
                             &numThreads, &numValuesPerThread, &offByOneIndex);
  
  if (numThreads <= 1) return computeOnlineUnrolledMean(x, length);
  
  MeanData threadData[numThreads];
  setupMeanData(threadData, numThreads, x, numValuesPerThread, offByOneIndex, computeOnlineUnrolledMean);
  
  void* threadDataPtrs[numThreads];
  for (size_t i = 0; i < numThreads; ++i) threadDataPtrs[i] = (void*) &threadData[i];
  
  
  
  misc_mt_runTasks(threadManager, &computeMeanTask, threadDataPtrs, numThreads);
  
  
  
  return aggregateMeanResults(threadData, numThreads);
}

static double mt_computeIndexedOnlineUnrolledMean(misc_mt_manager_t restrict threadManager, const double* restrict x, const size_t* restrict indices, size_t length)
{
  size_t numThreads, numValuesPerThread, offByOneIndex;
  misc_mt_getNumThreadsForJob(threadManager, length, INDEXED_ONLINE_UNROLLED_MEAN_MIN_NUM_VALUES_PER_THREAD,
                             &numThreads, &numValuesPerThread, &offByOneIndex);
  
  if (numThreads <= 1) return computeIndexedOnlineUnrolledMean(x, indices, length);
  
  IndexedMeanData threadData[numThreads];
  setupIndexedMeanData(threadData, numThreads, x, indices, numValuesPerThread, offByOneIndex, computeIndexedOnlineUnrolledMean);
  
  void* threadDataPtrs[numThreads];
  for (size_t i = 0; i < numThreads; ++i) threadDataPtrs[i] = (void*) &threadData[i];
  
  
  
  misc_mt_runTasks(threadManager, &computeIndexedMeanTask, threadDataPtrs, numThreads);
  
  
  
  return aggregateIndexedMeanResults(threadData, numThreads);
}


typedef double (*varianceForKnownMeanFunction)(const double* x, size_t length, double mean);
typedef double (*indexedVarianceForKnownMeanFunction)(const double* restrict x, const size_t* restrict indices, size_t length, double mean);

typedef struct {
  const double* x;
  size_t length;
  double mean;
  double result;
  varianceForKnownMeanFunction function;
} VarianceForKnownMeanData;

typedef struct {
  const double* x;
  const size_t* indices;
  size_t length;
  double mean;
  double result;
  indexedVarianceForKnownMeanFunction function;
} IndexedVarianceForKnownMeanData;


static void computeVarianceForKnownMeanTask(void* v_data)
{
  VarianceForKnownMeanData* data = (VarianceForKnownMeanData*) v_data;
  data->result = data->function(data->x, data->length, data->mean);
}

static void computeIndexedVarianceForKnownMeanTask(void* v_data)
{
  IndexedVarianceForKnownMeanData* data = (IndexedVarianceForKnownMeanData*) v_data;
  data->result = data->function(data->x, data->indices, data->length, data->mean);
}

static void setupVarianceForKnownMeanData(VarianceForKnownMeanData* restrict threadData, size_t numThreads, const double* restrict x,
                                         size_t numValuesPerThread, size_t offByOneIndex, varianceForKnownMeanFunction function, double mean)
{
  size_t i = 0;
  for ( ; i < offByOneIndex; ++i) {
    threadData[i].x = x + i * numValuesPerThread;
    threadData[i].length = numValuesPerThread;
    threadData[i].function = function;
    threadData[i].mean = mean;
  }
  for ( ; i < numThreads; ++i) {
    threadData[i].x = x + offByOneIndex * numValuesPerThread + (i - offByOneIndex) * (numValuesPerThread - 1);
    threadData[i].length = numValuesPerThread - 1;
    threadData[i].function = function;
    threadData[i].mean = mean;
  }
}

static void setupIndexedVarianceForKnownMeanData(IndexedVarianceForKnownMeanData* restrict threadData, size_t numThreads, const double* restrict x,
                                                 const size_t* restrict indices, size_t numValuesPerThread, size_t offByOneIndex, indexedVarianceForKnownMeanFunction function, double mean)
{
  size_t i = 0;
  for ( ; i < offByOneIndex; ++i) {
    threadData[i].x = x;
    threadData[i].indices = indices + i * numValuesPerThread;
    threadData[i].length = numValuesPerThread;
    threadData[i].function = function;
    threadData[i].mean = mean;
  }
  for ( ; i < numThreads; ++i) {
    threadData[i].x = x;
    threadData[i].indices = indices + offByOneIndex * numValuesPerThread + (i - offByOneIndex) * (numValuesPerThread - 1);
    threadData[i].length = numValuesPerThread - 1;
    threadData[i].function = function;
    threadData[i].mean = mean;
  }
}

static double aggregateVarianceForKnownMeanData(const VarianceForKnownMeanData* threadData, size_t numThreads)
{
  double result = threadData[0].result;
  size_t setSize = threadData[0].length;
  for (size_t i = 1; i < numThreads; ++i) {
    double n2 = (double) threadData[i].length;
    setSize += threadData[i].length;
    double n =  (double) setSize;
        
    result += ((n2 - 1.0) / (n - 1.0)) * (threadData[i].result - result) - result / (n - 1.0);
  }
  return result;
}

static double aggregateIndexedVarianceForKnownMeanData(const IndexedVarianceForKnownMeanData* threadData, size_t numThreads)
{
  double result = threadData[0].result;
  size_t setSize = threadData[0].length;
  for (size_t i = 1; i < numThreads; ++i) {
    double n2 = (double) threadData[i].length;
    setSize += threadData[i].length;
    double n =  (double) setSize;
        
    result += ((n2 - 1.0) / (n - 1.0)) * (threadData[i].result - result) - result / (n - 1.0);
  }
  return result;
}

/*
static double mt_computeVarianceForKnownMean(misc_mt_manager_t restrict threadManager, const double* restrict x,
                                      size_t length, double mean)
{
  size_t numThreads, numValuesPerThread, offByOneIndex;
  misc_mt_getNumThreadsForJob(threadManager, length, VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD,
                           &numThreads, &numValuesPerThread, &offByOneIndex);
  
  if (numThreads <= 1) return computeVarianceForKnownMean(x, length, mean);
  
  VarianceForKnownMeanData threadData[numThreads];
  setupVarianceForKnownMeanData(threadData, numThreads, x, numValuesPerThread, offByOneIndex, &computeVarianceForKnownMean, mean);
  
  void* threadDataPtrs[numThreads];
  for (size_t i = 0; i < numThreads; ++i) threadDataPtrs[i] = (void*) &threadData[i];
  
    
  misc_mt_runTasks(threadManager, &computeVarianceForKnownMeanTask, threadDataPtrs, numThreads);
  
  return aggregateVarianceForKnownMeanData(threadData, numThreads);
}

static double mt_computeIndexedVarianceForKnownMean(misc_mt_manager_t restrict threadManager, const double* restrict x,
                                                    const size_t* restrict indices, size_t length, double mean)
{
  size_t numThreads, numValuesPerThread, offByOneIndex;
  misc_mt_getNumThreadsForJob(threadManager, length, INDEXED_VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD,
                             &numThreads, &numValuesPerThread, &offByOneIndex);
  
  if (numThreads <= 1) return computeIndexedVarianceForKnownMean(x, indices, length, mean);
  
  IndexedVarianceForKnownMeanData threadData[numThreads];
  setupIndexedVarianceForKnownMeanData(threadData, numThreads, x, indices, numValuesPerThread, offByOneIndex, &computeIndexedVarianceForKnownMean, mean);
  
  void* threadDataPtrs[numThreads];
  for (size_t i = 0; i < numThreads; ++i) threadDataPtrs[i] = (void*) &threadData[i];
  
  
  misc_mt_runTasks(threadManager, &computeIndexedVarianceForKnownMeanTask, threadDataPtrs, numThreads);
  
  return aggregateIndexedVarianceForKnownMeanData(threadData, numThreads);
}*/

static double mt_computeUnrolledVarianceForKnownMean(misc_mt_manager_t restrict threadManager, const double* restrict x,
                                              size_t length, double mean)
{
  size_t numThreads, numValuesPerThread, offByOneIndex;
  misc_mt_getNumThreadsForJob(threadManager, length, UNROLLED_VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD,
                           &numThreads, &numValuesPerThread, &offByOneIndex);
  
  if (numThreads <= 1) return computeUnrolledVarianceForKnownMean(x, length, mean);
  
  VarianceForKnownMeanData threadData[numThreads];
  setupVarianceForKnownMeanData(threadData, numThreads, x, numValuesPerThread, offByOneIndex, computeUnrolledVarianceForKnownMean, mean);
  
  void* threadDataPtrs[numThreads];
  for (size_t i = 0; i < numThreads; ++i) threadDataPtrs[i] = (void*) &threadData[i];
  
  
  misc_mt_runTasks(threadManager, &computeVarianceForKnownMeanTask, threadDataPtrs, numThreads);
  
  return aggregateVarianceForKnownMeanData(threadData, numThreads);
}

static double mt_computeIndexedUnrolledVarianceForKnownMean(misc_mt_manager_t restrict threadManager, const double* restrict x,
                                                            const size_t* restrict indices, size_t length, double mean)
{
  size_t numThreads, numValuesPerThread, offByOneIndex;
  misc_mt_getNumThreadsForJob(threadManager, length, INDEXED_UNROLLED_VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD,
                             &numThreads, &numValuesPerThread, &offByOneIndex);
  
  if (numThreads <= 1) return computeIndexedUnrolledVarianceForKnownMean(x, indices, length, mean);
  
  IndexedVarianceForKnownMeanData threadData[numThreads];
  setupIndexedVarianceForKnownMeanData(threadData, numThreads, x, indices, numValuesPerThread, offByOneIndex, computeIndexedUnrolledVarianceForKnownMean, mean);
  
  void* threadDataPtrs[numThreads];
  for (size_t i = 0; i < numThreads; ++i) threadDataPtrs[i] = (void*) &threadData[i];
  
  
  misc_mt_runTasks(threadManager, &computeIndexedVarianceForKnownMeanTask, threadDataPtrs, numThreads);
  
  return aggregateIndexedVarianceForKnownMeanData(threadData, numThreads);
}

/*
static double mt_computeOnlineVarianceForKnownMean(misc_mt_manager_t restrict threadManager, const double* restrict x,
                                                 size_t length, double mean)
{
  size_t numThreads, numValuesPerThread, offByOneIndex;
  misc_mt_getNumThreadsForJob(threadManager, length, ONLINE_VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD,
                           &numThreads, &numValuesPerThread, &offByOneIndex);
  
  if (numThreads <= 1) return computeOnlineVarianceForKnownMean(x, length, mean);
  
  VarianceForKnownMeanData threadData[numThreads];
  setupVarianceForKnownMeanData(threadData, numThreads, x, numValuesPerThread, offByOneIndex, &computeOnlineVarianceForKnownMean, mean);
  
  void* threadDataPtrs[numThreads];
  for (size_t i = 0; i < numThreads; ++i) threadDataPtrs[i] = (void*) &threadData[i];
  
  
  misc_mt_runTasks(threadManager, &computeVarianceForKnownMeanTask, threadDataPtrs, numThreads);
  
  return aggregateVarianceForKnownMeanData(threadData, numThreads);
}

static double mt_computeIndexedOnlineVarianceForKnownMean(misc_mt_manager_t restrict threadManager, const double* restrict x,
                                                          const size_t* restrict indices, size_t length, double mean)
{
  size_t numThreads, numValuesPerThread, offByOneIndex;
  misc_mt_getNumThreadsForJob(threadManager, length, INDEXED_ONLINE_VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD,
                             &numThreads, &numValuesPerThread, &offByOneIndex);
  
  if (numThreads <= 1) return computeIndexedOnlineVarianceForKnownMean(x, indices, length, mean);
  
  IndexedVarianceForKnownMeanData threadData[numThreads];
  setupIndexedVarianceForKnownMeanData(threadData, numThreads, x, indices, numValuesPerThread, offByOneIndex, &computeIndexedOnlineVarianceForKnownMean, mean);
  
  void* threadDataPtrs[numThreads];
  for (size_t i = 0; i < numThreads; ++i) threadDataPtrs[i] = (void*) &threadData[i];
  
  
  misc_mt_runTasks(threadManager, &computeIndexedVarianceForKnownMeanTask, threadDataPtrs, numThreads);
  
  return aggregateIndexedVarianceForKnownMeanData(threadData, numThreads);
} */

static double mt_computeOnlineUnrolledVarianceForKnownMean(misc_mt_manager_t restrict threadManager, const double* restrict x,
                                                    size_t length, double mean)
{
  size_t numThreads, numValuesPerThread, offByOneIndex;
  misc_mt_getNumThreadsForJob(threadManager, length, ONLINE_UNROLLED_VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD,
                             &numThreads, &numValuesPerThread, &offByOneIndex);
  
  if (numThreads <= 1) return computeOnlineUnrolledVarianceForKnownMean(x, length, mean);
  
  VarianceForKnownMeanData threadData[numThreads];
  setupVarianceForKnownMeanData(threadData, numThreads, x, numValuesPerThread, offByOneIndex, computeOnlineUnrolledVarianceForKnownMean, mean);
  
  void* threadDataPtrs[numThreads];
  for (size_t i = 0; i < numThreads; ++i) threadDataPtrs[i] = (void*) &threadData[i];
  
  
  misc_mt_runTasks(threadManager, &computeVarianceForKnownMeanTask, threadDataPtrs, numThreads);
  
  return aggregateVarianceForKnownMeanData(threadData, numThreads);
}

static double mt_computeIndexedOnlineUnrolledVarianceForKnownMean(misc_mt_manager_t restrict threadManager, const double* restrict x,
                                                           const size_t* restrict indices, size_t length, double mean)
{
  size_t numThreads, numValuesPerThread, offByOneIndex;
  misc_mt_getNumThreadsForJob(threadManager, length, INDEXED_ONLINE_UNROLLED_VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD,
                             &numThreads, &numValuesPerThread, &offByOneIndex);
  
  if (numThreads <= 1) return computeIndexedOnlineUnrolledVarianceForKnownMean(x, indices, length, mean);
  
  IndexedVarianceForKnownMeanData threadData[numThreads];
  setupIndexedVarianceForKnownMeanData(threadData, numThreads, x, indices, numValuesPerThread, offByOneIndex, computeIndexedOnlineUnrolledVarianceForKnownMean, mean);
  
  void* threadDataPtrs[numThreads];
  for (size_t i = 0; i < numThreads; ++i) threadDataPtrs[i] = (void*) &threadData[i];
  
  
  misc_mt_runTasks(threadManager, &computeIndexedVarianceForKnownMeanTask, threadDataPtrs, numThreads);
  
  return aggregateIndexedVarianceForKnownMeanData(threadData, numThreads);
}

typedef double (*varianceFunction)(const double* restrict x, size_t length, double* restrict meanPtr);
typedef double (*indexedVarianceFunction)(const double* restrict x, const size_t* restrict indices, size_t length, double* restrict meanPtr);

typedef struct {
  const double* x;
  size_t length;
  double mean;
  double var;
  varianceFunction function;
} VarianceData;

typedef struct {
  const double* x;
  const size_t* indices;
  size_t length;
  double mean;
  double var;
  indexedVarianceFunction function;
} IndexedVarianceData;

static void computeVarianceTask(void* v_data) {
  VarianceData* data = (VarianceData*) v_data;
  data->var = data->function(data->x, data->length, &data->mean);
}

static void computeIndexedVarianceTask(void* v_data) {
  IndexedVarianceData* data = (IndexedVarianceData*) v_data;
  data->var = data->function(data->x, data->indices, data->length, &data->mean);
}

static void setupVarianceData(VarianceData* restrict threadData, size_t numThreads, const double* restrict x,
                              size_t numValuesPerThread, size_t offByOneIndex, varianceFunction function)
{
  size_t i = 0;
  for ( ; i < offByOneIndex; ++i) {
    threadData[i].x = x + i * numValuesPerThread;
    threadData[i].length = numValuesPerThread;
    threadData[i].function = function;
  }
  for ( ; i < numThreads; ++i) {
    threadData[i].x = x + offByOneIndex * numValuesPerThread + (i - offByOneIndex) * (numValuesPerThread - 1);
    threadData[i].length = numValuesPerThread - 1;
    threadData[i].function = function;
  }
}

static void setupIndexedVarianceData(IndexedVarianceData* restrict threadData, size_t numThreads, const double* restrict x,
                                     const size_t* restrict indices,
                                     size_t numValuesPerThread, size_t offByOneIndex, indexedVarianceFunction function)
{
  size_t i = 0;
  for ( ; i < offByOneIndex; ++i) {
    threadData[i].x = x;
    threadData[i].indices = indices + i * numValuesPerThread;
    threadData[i].length = numValuesPerThread;
    threadData[i].function = function;
  }
  for ( ; i < numThreads; ++i) {
    threadData[i].x = x;
    threadData[i].indices = indices + offByOneIndex * numValuesPerThread + (i - offByOneIndex) * (numValuesPerThread - 1);
    threadData[i].length = numValuesPerThread - 1;
    threadData[i].function = function;
  }
}

static double aggregateVarianceData(const VarianceData* restrict threadData, size_t numThreads, double* restrict meanPtr)
{
  double var  = threadData[0].var;
  double mean = threadData[0].mean;
  size_t setSize = threadData[0].length;
  for (size_t i = 1; i < numThreads; ++i) {
    double n1 = (double) setSize;
    double n2 = (double) threadData[i].length;
    double n =  (double) (threadData[i].length + setSize);
    setSize += threadData[i].length;
    
    var  += (threadData[i].var - var) * (n2 / (n - 1.0)) - threadData[i].var / (n - 1.0) + 
            ((threadData[i].mean - mean) * (n1 / (n - 1.0))) * ((threadData[i].mean - mean) * (n2 / n));
    mean += (n2 / n) * (threadData[i].mean - mean);
  }
  if (meanPtr != NULL) *meanPtr = mean;
  return var;
}

static double aggregateIndexedVarianceData(const IndexedVarianceData* restrict threadData, size_t numThreads, double* restrict meanPtr)
{
  double var  = threadData[0].var;
  double mean = threadData[0].mean;
  size_t setSize = threadData[0].length;
  for (size_t i = 1; i < numThreads; ++i) {
    double n1 = (double) setSize;
    double n2 = (double) threadData[i].length;
    double n =  (double) (threadData[i].length + setSize);
    setSize += threadData[i].length;
    
    var  += (threadData[i].var - var) * (n2 / (n - 1.0)) - threadData[i].var / (n - 1.0) +
    ((threadData[i].mean - mean) * (n1 / (n - 1.0))) * ((threadData[i].mean - mean) * (n2 / n));
    mean += (n2 / n) * (threadData[i].mean - mean);
  }
  if (meanPtr != NULL) *meanPtr = mean;
  return var;
}

/*
static double mt_computeVariance(misc_mt_manager_t restrict threadManager, const double* restrict x,
                          size_t length, double* restrict meanPtr)
{
  size_t numThreads, numValuesPerThread, offByOneIndex;
  misc_mt_getNumThreadsForJob(threadManager, length, VAR_MIN_NUM_VALUES_PER_THREAD,
                           &numThreads, &numValuesPerThread, &offByOneIndex);
  
  if (numThreads <= 1) return computeVariance(x, length, meanPtr);
  
  VarianceData threadData[numThreads];
  setupVarianceData(threadData, numThreads, x, numValuesPerThread, offByOneIndex, &computeVariance);
  
  void* threadDataPtrs[numThreads];
  for (size_t i = 0; i < numThreads; ++i) threadDataPtrs[i] = (void*) &threadData[i];
  
  
  misc_mt_runTasks(threadManager, &computeVarianceTask, threadDataPtrs, numThreads);
  
  
  return aggregateVarianceData(threadData, numThreads, meanPtr);
}

static double mt_computeIndexedVariance(misc_mt_manager_t restrict threadManager, const double* restrict x,
                                        const size_t* restrict indices, size_t length, double* restrict meanPtr)
{
  size_t numThreads, numValuesPerThread, offByOneIndex;
  misc_mt_getNumThreadsForJob(threadManager, length, INDEXED_VAR_MIN_NUM_VALUES_PER_THREAD,
                             &numThreads, &numValuesPerThread, &offByOneIndex);
  
  if (numThreads <= 1) return computeIndexedVariance(x, indices, length, meanPtr);
  
  IndexedVarianceData threadData[numThreads];
  setupIndexedVarianceData(threadData, numThreads, x, indices, numValuesPerThread, offByOneIndex, &computeIndexedVariance);
  
  void* threadDataPtrs[numThreads];
  for (size_t i = 0; i < numThreads; ++i) threadDataPtrs[i] = (void*) &threadData[i];
  
  
  misc_mt_runTasks(threadManager, &computeIndexedVarianceTask, threadDataPtrs, numThreads);
  
  
  return aggregateIndexedVarianceData(threadData, numThreads, meanPtr);
} */


static double mt_computeUnrolledVariance(misc_mt_manager_t restrict threadManager, const double* restrict x,
                                  size_t length, double* restrict meanPtr)
{
  size_t numThreads, numValuesPerThread, offByOneIndex;
  misc_mt_getNumThreadsForJob(threadManager, length, UNROLLED_VAR_MIN_NUM_VALUES_PER_THREAD,
                           &numThreads, &numValuesPerThread, &offByOneIndex);
  
  if (numThreads <= 1) return computeUnrolledVariance(x, length, meanPtr);
  
  VarianceData threadData[numThreads];
  setupVarianceData(threadData, numThreads, x, numValuesPerThread, offByOneIndex, &computeUnrolledVariance);
  
  void* threadDataPtrs[numThreads];
  for (size_t i = 0; i < numThreads; ++i) threadDataPtrs[i] = (void*) &threadData[i];
  
  
  misc_mt_runTasks(threadManager, &computeVarianceTask, threadDataPtrs, numThreads);
  
  
  return aggregateVarianceData(threadData, numThreads, meanPtr);
}

static double mt_computeIndexedUnrolledVariance(misc_mt_manager_t restrict threadManager, const double* restrict x,
                                                const size_t* restrict indices, size_t length, double* restrict meanPtr)
{
  size_t numThreads, numValuesPerThread, offByOneIndex;
  misc_mt_getNumThreadsForJob(threadManager, length, INDEXED_UNROLLED_VAR_MIN_NUM_VALUES_PER_THREAD,
                             &numThreads, &numValuesPerThread, &offByOneIndex);
  
  if (numThreads <= 1) return computeIndexedUnrolledVariance(x,  indices, length, meanPtr);
  
  IndexedVarianceData threadData[numThreads];
  setupIndexedVarianceData(threadData, numThreads, x, indices, numValuesPerThread, offByOneIndex, &computeIndexedUnrolledVariance);
  
  void* threadDataPtrs[numThreads];
  for (size_t i = 0; i < numThreads; ++i) threadDataPtrs[i] = (void*) &threadData[i];
  
  
  misc_mt_runTasks(threadManager, &computeIndexedVarianceTask, threadDataPtrs, numThreads);
  
  
  return aggregateIndexedVarianceData(threadData, numThreads, meanPtr);
}

/*
static double mt_computeOnlineVariance(misc_mt_manager_t restrict threadManager, const double* restrict x,
                                size_t length, double* restrict meanPtr)
{
  size_t numThreads, numValuesPerThread, offByOneIndex;
  misc_mt_getNumThreadsForJob(threadManager, length, ONLINE_VAR_MIN_NUM_VALUES_PER_THREAD,
                           &numThreads, &numValuesPerThread, &offByOneIndex);
  
  if (numThreads <= 1) return computeOnlineVariance(x, length, meanPtr);
  
  VarianceData threadData[numThreads];
  setupVarianceData(threadData, numThreads, x, numValuesPerThread, offByOneIndex, &computeOnlineVariance);
  
  void* threadDataPtrs[numThreads];
  for (size_t i = 0; i < numThreads; ++i) threadDataPtrs[i] = (void*) &threadData[i];
  
  
  misc_mt_runTasks(threadManager, &computeVarianceTask, threadDataPtrs, numThreads);
  
  
  return aggregateVarianceData(threadData, numThreads, meanPtr);
}

static double mt_computeIndexedOnlineVariance(misc_mt_manager_t restrict threadManager, const double* restrict x,
                                              const size_t* restrict indices, size_t length, double* restrict meanPtr)
{
  size_t numThreads, numValuesPerThread, offByOneIndex;
  misc_mt_getNumThreadsForJob(threadManager, length, INDEXED_ONLINE_VAR_MIN_NUM_VALUES_PER_THREAD,
                             &numThreads, &numValuesPerThread, &offByOneIndex);
  
  if (numThreads <= 1) return computeIndexedOnlineVariance(x, indices, length, meanPtr);
  
  IndexedVarianceData threadData[numThreads];
  setupIndexedVarianceData(threadData, numThreads, x, indices, numValuesPerThread, offByOneIndex, &computeIndexedOnlineVariance);
  
  void* threadDataPtrs[numThreads];
  for (size_t i = 0; i < numThreads; ++i) threadDataPtrs[i] = (void*) &threadData[i];
  
  
  misc_mt_runTasks(threadManager, &computeIndexedVarianceTask, threadDataPtrs, numThreads);
  
  
  return aggregateIndexedVarianceData(threadData, numThreads, meanPtr);
} */

static double mt_computeOnlineUnrolledVariance(misc_mt_manager_t restrict threadManager, const double* restrict x,
                                        size_t length, double* restrict meanPtr)
{
  size_t numThreads, numValuesPerThread, offByOneIndex;
  misc_mt_getNumThreadsForJob(threadManager, length, ONLINE_UNROLLED_VAR_MIN_NUM_VALUES_PER_THREAD,
                           &numThreads, &numValuesPerThread, &offByOneIndex);
  
  if (numThreads <= 1) return computeOnlineUnrolledVariance(x, length, meanPtr);
  
  VarianceData threadData[numThreads];
  setupVarianceData(threadData, numThreads, x, numValuesPerThread, offByOneIndex, &computeOnlineUnrolledVariance);
  
  void* threadDataPtrs[numThreads];
  for (size_t i = 0; i < numThreads; ++i) threadDataPtrs[i] = (void*) &threadData[i];
  
  
  misc_mt_runTasks(threadManager, &computeVarianceTask, threadDataPtrs, numThreads);
  
  
  return aggregateVarianceData(threadData, numThreads, meanPtr);
}

static double mt_computeIndexedOnlineUnrolledVariance(misc_mt_manager_t restrict threadManager, const double* restrict x,
                                               const size_t* restrict indices, size_t length, double* restrict meanPtr)
{
  size_t numThreads, numValuesPerThread, offByOneIndex;
  misc_mt_getNumThreadsForJob(threadManager, length, INDEXED_ONLINE_UNROLLED_VAR_MIN_NUM_VALUES_PER_THREAD,
                           &numThreads, &numValuesPerThread, &offByOneIndex);
  
  if (numThreads <= 1) return computeIndexedOnlineUnrolledVariance(x, indices, length, meanPtr);
  
  IndexedVarianceData threadData[numThreads];
  setupIndexedVarianceData(threadData, numThreads, x, indices, numValuesPerThread, offByOneIndex, &computeIndexedOnlineUnrolledVariance);
  
  void* threadDataPtrs[numThreads];
  for (size_t i = 0; i < numThreads; ++i) threadDataPtrs[i] = (void*) &threadData[i];
  
  
  misc_mt_runTasks(threadManager, &computeIndexedVarianceTask, threadDataPtrs, numThreads);
  
  
  return aggregateIndexedVarianceData(threadData, numThreads, meanPtr);
}

typedef double (*weightedMeanFunction)(const double* restrict x, size_t length, const double* restrict w, double* restrict n);
typedef double (*indexedWeightedMeanFunction)(const double* restrict x, const size_t* restrict indices, size_t length, const double* restrict w, double* restrict n);

typedef struct {
  const double* x;
  size_t length;
  const double* w;
  double n;
  double result;
  weightedMeanFunction function;
} WeightedMeanData;

typedef struct {
  const double* x;
  const size_t* indices;
  size_t length;
  const double* w;
  double n;
  double result;
  indexedWeightedMeanFunction function;
} IndexedWeightedMeanData;

static void computeWeightedMeanTask(void* v_data)
{
  WeightedMeanData* data = (WeightedMeanData*) v_data;
  data->result = data->function(data->x, data->length, data->w, &data->n);
}

static void computeIndexedWeightedMeanTask(void* v_data)
{
  IndexedWeightedMeanData* data = (IndexedWeightedMeanData*) v_data;
  data->result = data->function(data->x, data->indices, data->length, data->w, &data->n);
}

static void setupWeightedMeanData(WeightedMeanData* restrict threadData, size_t numThreads, const double* restrict x, const double* restrict w,
                                  size_t numValuesPerThread, size_t offByOneIndex, weightedMeanFunction function)
{
  size_t i = 0;
  for ( ; i < offByOneIndex; ++i) {
    threadData[i].x = x + i * numValuesPerThread;
    threadData[i].length = numValuesPerThread;
    threadData[i].w = w + i * numValuesPerThread;
    threadData[i].function = function;
  }
  for ( ; i < numThreads; ++i) {
    threadData[i].x = x + offByOneIndex * numValuesPerThread + (i - offByOneIndex) * (numValuesPerThread - 1);
    threadData[i].length = numValuesPerThread - 1;
    threadData[i].w = w + offByOneIndex * numValuesPerThread + (i - offByOneIndex) * (numValuesPerThread - 1);
    threadData[i].function = function;
  }
}

static void setupIndexedWeightedMeanData(IndexedWeightedMeanData* restrict threadData, size_t numThreads, const double* restrict x,
                                         const size_t* restrict indices, const double* restrict w,
                                         size_t numValuesPerThread, size_t offByOneIndex, indexedWeightedMeanFunction function)
{
  size_t i = 0;
  for ( ; i < offByOneIndex; ++i) {
    threadData[i].x = x;
    threadData[i].indices = indices + i * numValuesPerThread;
    threadData[i].length = numValuesPerThread;
    threadData[i].w = w;
    threadData[i].function = function;
  }
  for ( ; i < numThreads; ++i) {
    threadData[i].x = x;
    threadData[i].indices = indices + offByOneIndex * numValuesPerThread + (i - offByOneIndex) * (numValuesPerThread - 1);
    threadData[i].length = numValuesPerThread - 1;
    threadData[i].w = w;
    threadData[i].function = function;
  }
}

static double aggregateWeightedMeanResults(const WeightedMeanData* restrict threadData, size_t numThreads, double* restrict nPtr)
{
  double result = threadData[0].result;
  double n      = threadData[0].n;
  for (size_t i = 1; i < numThreads; ++i) {
    n += threadData[i].n;
    result += (threadData[i].n / n) * (threadData[i].result - result);
  }
  if (nPtr != NULL) *nPtr = n;
  return result;
}

static double aggregateIndexedWeightedMeanResults(const IndexedWeightedMeanData* restrict threadData, size_t numThreads, double* restrict nPtr)
{
  double result = threadData[0].result;
  double n      = threadData[0].n;
  for (size_t i = 1; i < numThreads; ++i) {
    n += threadData[i].n;
    result += (threadData[i].n / n) * (threadData[i].result - result);
  }
  if (nPtr != NULL) *nPtr = n;
  return result;
}

/* 
static double mt_computeWeightedMean(misc_mt_manager_t restrict threadManager, const double* restrict x, size_t length, const double* restrict w, double* restrict nPtr)
{
  size_t numThreads, numValuesPerThread, offByOneIndex;
  misc_mt_getNumThreadsForJob(threadManager, length, WEIGHTED_MEAN_MIN_NUM_VALUES_PER_THREAD,
                             &numThreads, &numValuesPerThread, &offByOneIndex);
  
  if (numThreads <= 1) return computeWeightedMean(x, length, w, nPtr);
  
  WeightedMeanData threadData[numThreads];
  setupWeightedMeanData(threadData, numThreads, x, w, numValuesPerThread, offByOneIndex, computeWeightedMean);
  
  void* threadDataPtrs[numThreads];
  for (size_t i = 0; i < numThreads; ++i) threadDataPtrs[i] = (void*) &threadData[i];
  
  
  
  misc_mt_runTasks(threadManager, &computeWeightedMeanTask, threadDataPtrs, numThreads);
  
  
  
  return aggregateWeightedMeanResults(threadData, numThreads, nPtr);
}

static double mt_computeIndexedWeightedMean(misc_mt_manager_t restrict threadManager, const double* restrict x, const size_t* restrict indices, size_t length,
                                            const double* restrict w, double* restrict nPtr)
{
  size_t numThreads, numValuesPerThread, offByOneIndex;
  misc_mt_getNumThreadsForJob(threadManager, length, INDEXED_WEIGHTED_MEAN_MIN_NUM_VALUES_PER_THREAD,
                             &numThreads, &numValuesPerThread, &offByOneIndex);
  
  if (numThreads <= 1) return computeIndexedWeightedMean(x, indices, length, w, nPtr);
  
  IndexedWeightedMeanData threadData[numThreads];
  setupIndexedWeightedMeanData(threadData, numThreads, x, indices, w, numValuesPerThread, offByOneIndex, computeIndexedWeightedMean);
  
  void* threadDataPtrs[numThreads];
  for (size_t i = 0; i < numThreads; ++i) threadDataPtrs[i] = (void*) &threadData[i];
  
  
  
  misc_mt_runTasks(threadManager, &computeIndexedWeightedMeanTask, threadDataPtrs, numThreads);
  
  
  
  return aggregateIndexedWeightedMeanResults(threadData, numThreads, nPtr);
} */

static double mt_computeUnrolledWeightedMean(misc_mt_manager_t restrict threadManager, const double* restrict x, size_t length, const double* restrict w, double* restrict nPtr)
{
  size_t numThreads, numValuesPerThread, offByOneIndex;
  misc_mt_getNumThreadsForJob(threadManager, length, UNROLLED_WEIGHTED_MEAN_MIN_NUM_VALUES_PER_THREAD,
                             &numThreads, &numValuesPerThread, &offByOneIndex);
  
  if (numThreads <= 1) return computeUnrolledWeightedMean(x, length, w, nPtr);
  
  WeightedMeanData threadData[numThreads];
  setupWeightedMeanData(threadData, numThreads, x, w, numValuesPerThread, offByOneIndex, computeUnrolledWeightedMean);
  
  void* threadDataPtrs[numThreads];
  for (size_t i = 0; i < numThreads; ++i) threadDataPtrs[i] = (void*) &threadData[i];
  
  
  
  misc_mt_runTasks(threadManager, &computeWeightedMeanTask, threadDataPtrs, numThreads);
  
  
  
  return aggregateWeightedMeanResults(threadData, numThreads, nPtr);
}

static double mt_computeIndexedUnrolledWeightedMean(misc_mt_manager_t restrict threadManager, const double* restrict x, const size_t* restrict indices, size_t length,
                                                    const double* restrict w, double* restrict nPtr)
{
  size_t numThreads, numValuesPerThread, offByOneIndex;
  misc_mt_getNumThreadsForJob(threadManager, length, INDEXED_UNROLLED_WEIGHTED_MEAN_MIN_NUM_VALUES_PER_THREAD,
                             &numThreads, &numValuesPerThread, &offByOneIndex);
  
  if (numThreads <= 1) return computeIndexedUnrolledWeightedMean(x, indices, length, w, nPtr);
  
  IndexedWeightedMeanData threadData[numThreads];
  setupIndexedWeightedMeanData(threadData, numThreads, x, indices, w, numValuesPerThread, offByOneIndex, computeIndexedUnrolledWeightedMean);
  
  void* threadDataPtrs[numThreads];
  for (size_t i = 0; i < numThreads; ++i) threadDataPtrs[i] = (void*) &threadData[i];
  
  
  
  misc_mt_runTasks(threadManager, &computeIndexedWeightedMeanTask, threadDataPtrs, numThreads);
  
  
  
  return aggregateIndexedWeightedMeanResults(threadData, numThreads, nPtr);
}

/* 
static double mt_computeOnlineWeightedMean(misc_mt_manager_t restrict threadManager, const double* restrict x, size_t length, const double* restrict w, double* restrict nPtr)
{
  size_t numThreads, numValuesPerThread, offByOneIndex;
  misc_mt_getNumThreadsForJob(threadManager, length, ONLINE_WEIGHTED_MEAN_MIN_NUM_VALUES_PER_THREAD,
                             &numThreads, &numValuesPerThread, &offByOneIndex);
  
  if (numThreads <= 1) return computeOnlineWeightedMean(x, length, w, nPtr);
  
  WeightedMeanData threadData[numThreads];
  setupWeightedMeanData(threadData, numThreads, x, w, numValuesPerThread, offByOneIndex, computeOnlineWeightedMean);
  
  void* threadDataPtrs[numThreads];
  for (size_t i = 0; i < numThreads; ++i) threadDataPtrs[i] = (void*) &threadData[i];
  
  
  
  misc_mt_runTasks(threadManager, &computeWeightedMeanTask, threadDataPtrs, numThreads);
  
  
  
  return aggregateWeightedMeanResults(threadData, numThreads, nPtr);
}

static double mt_computeIndexedOnlineWeightedMean(misc_mt_manager_t restrict threadManager, const double* restrict x, const size_t* restrict indices, size_t length,
                                                  const double* restrict w, double* restrict nPtr)
{
  size_t numThreads, numValuesPerThread, offByOneIndex;
  misc_mt_getNumThreadsForJob(threadManager, length, INDEXED_ONLINE_WEIGHTED_MEAN_MIN_NUM_VALUES_PER_THREAD,
                             &numThreads, &numValuesPerThread, &offByOneIndex);
  
  if (numThreads <= 1) return computeIndexedOnlineWeightedMean(x, indices, length, w, nPtr);
  
  IndexedWeightedMeanData threadData[numThreads];
  setupIndexedWeightedMeanData(threadData, numThreads, x, indices, w, numValuesPerThread, offByOneIndex, computeIndexedOnlineWeightedMean);
  
  void* threadDataPtrs[numThreads];
  for (size_t i = 0; i < numThreads; ++i) threadDataPtrs[i] = (void*) &threadData[i];
  
  
  
  misc_mt_runTasks(threadManager, &computeIndexedWeightedMeanTask, threadDataPtrs, numThreads);
  
  
  
  return aggregateIndexedWeightedMeanResults(threadData, numThreads, nPtr);
} */

static double mt_computeOnlineUnrolledWeightedMean(misc_mt_manager_t restrict threadManager, const double* restrict x, size_t length, const double* restrict w, double* restrict nPtr)
{
  size_t numThreads, numValuesPerThread, offByOneIndex;
  misc_mt_getNumThreadsForJob(threadManager, length, ONLINE_UNROLLED_WEIGHTED_MEAN_MIN_NUM_VALUES_PER_THREAD,
                             &numThreads, &numValuesPerThread, &offByOneIndex);
  
  if (numThreads <= 1) return computeOnlineUnrolledWeightedMean(x, length, w, nPtr);
  
  WeightedMeanData threadData[numThreads];
  setupWeightedMeanData(threadData, numThreads, x, w, numValuesPerThread, offByOneIndex, computeOnlineUnrolledWeightedMean);
  
  void* threadDataPtrs[numThreads];
  for (size_t i = 0; i < numThreads; ++i) threadDataPtrs[i] = (void*) &threadData[i];
  
  
  
  misc_mt_runTasks(threadManager, &computeWeightedMeanTask, threadDataPtrs, numThreads);
  
  
  
  return aggregateWeightedMeanResults(threadData, numThreads, nPtr);
}

static double mt_computeIndexedOnlineUnrolledWeightedMean(misc_mt_manager_t restrict threadManager, const double* restrict x, const size_t* restrict indices, size_t length,
                                                          const double* restrict w, double* restrict nPtr)
{
  size_t numThreads, numValuesPerThread, offByOneIndex;
  misc_mt_getNumThreadsForJob(threadManager, length, INDEXED_ONLINE_UNROLLED_WEIGHTED_MEAN_MIN_NUM_VALUES_PER_THREAD,
                             &numThreads, &numValuesPerThread, &offByOneIndex);
  
  if (numThreads <= 1) return computeIndexedOnlineUnrolledWeightedMean(x, indices, length, w, nPtr);
  
  IndexedWeightedMeanData threadData[numThreads];
  setupIndexedWeightedMeanData(threadData, numThreads, x, indices, w, numValuesPerThread, offByOneIndex, computeIndexedOnlineUnrolledWeightedMean);
  
  void* threadDataPtrs[numThreads];
  for (size_t i = 0; i < numThreads; ++i) threadDataPtrs[i] = (void*) &threadData[i];
  
  
  
  misc_mt_runTasks(threadManager, &computeIndexedWeightedMeanTask, threadDataPtrs, numThreads);
  
  
  
  return aggregateIndexedWeightedMeanResults(threadData, numThreads, nPtr);
}

typedef double (*weightedVarianceForKnownMeanFunction)(const double* restrict x, size_t length, const double* restrict w, double mean);
typedef double (*indexedWeightedVarianceForKnownMeanFunction)(const double* restrict x, const size_t* restrict indices, size_t length,  const double* restrict w, double mean);

typedef struct {
  const double* x;
  size_t length;
  const double* w;
  double mean;
  double result;
  weightedVarianceForKnownMeanFunction function;
} WeightedVarianceForKnownMeanData;

typedef struct {
  const double* x;
  const size_t* indices;
  size_t length;
  const double* w;
  double mean;
  double result;
  indexedWeightedVarianceForKnownMeanFunction function;
} IndexedWeightedVarianceForKnownMeanData;


static void computeWeightedVarianceForKnownMeanTask(void* v_data)
{
  WeightedVarianceForKnownMeanData* data = (WeightedVarianceForKnownMeanData*) v_data;
  data->result = data->function(data->x, data->length, data->w, data->mean);
}

static void computeIndexedWeightedVarianceForKnownMeanTask(void* v_data)
{
  IndexedWeightedVarianceForKnownMeanData* data = (IndexedWeightedVarianceForKnownMeanData*) v_data;
  data->result = data->function(data->x, data->indices, data->length, data->w, data->mean);
}

static void setupWeightedVarianceForKnownMeanData(WeightedVarianceForKnownMeanData* restrict threadData, size_t numThreads,
                                                  const double* restrict x, const double* restrict w,
                                                  size_t numValuesPerThread, size_t offByOneIndex, weightedVarianceForKnownMeanFunction function, double mean)
{
  size_t i = 0;
  for ( ; i < offByOneIndex; ++i) {
    threadData[i].x = x + i * numValuesPerThread;
    threadData[i].length = numValuesPerThread;
    threadData[i].w = w + i * numValuesPerThread;
    threadData[i].function = function;
    threadData[i].mean = mean;
  }
  for ( ; i < numThreads; ++i) {
    threadData[i].x = x + offByOneIndex * numValuesPerThread + (i - offByOneIndex) * (numValuesPerThread - 1);
    threadData[i].length = numValuesPerThread - 1;
    threadData[i].w = w + offByOneIndex * numValuesPerThread + (i - offByOneIndex) * (numValuesPerThread - 1);
    threadData[i].function = function;
    threadData[i].mean = mean;
  }
}

static void setupIndexedWeightedVarianceForKnownMeanData(IndexedWeightedVarianceForKnownMeanData* restrict threadData, size_t numThreads, const double* restrict x,
                                                         const size_t* restrict indices, const double* restrict w,
                                                         size_t numValuesPerThread, size_t offByOneIndex, indexedWeightedVarianceForKnownMeanFunction function, double mean)
{
  size_t i = 0;
  for ( ; i < offByOneIndex; ++i) {
    threadData[i].x = x;
    threadData[i].indices = indices + i * numValuesPerThread;
    threadData[i].length = numValuesPerThread;
    threadData[i].w = w;
    threadData[i].function = function;
    threadData[i].mean = mean;
  }
  for (; i < numThreads; ++i) {
    threadData[i].x = x;
    threadData[i].indices = indices + offByOneIndex * numValuesPerThread + (i - offByOneIndex) * (numValuesPerThread - 1);
    threadData[i].length = numValuesPerThread - 1;
    threadData[i].w = w;
    threadData[i].function = function;
    threadData[i].mean = mean;
  }
}

static double aggregateWeightedVarianceForKnownMeanData(const WeightedVarianceForKnownMeanData* threadData, size_t numThreads)
{
  double result = threadData[0].result;
  size_t setSize = threadData[0].length;
  for (size_t i = 1; i < numThreads; ++i) {
    double n2 = (double) threadData[i].length;
    setSize += threadData[i].length;
    double n =  (double) setSize;
        
    result += ((n2 - 1.0) / (n - 1.0)) * (threadData[i].result - result) - result / (n - 1.0);
  }
  return result;
}

static double aggregateIndexedWeightedVarianceForKnownMeanData(const IndexedWeightedVarianceForKnownMeanData* threadData, size_t numThreads)
{
  double result = threadData[0].result;
  size_t setSize = threadData[0].length;
  for (size_t i = 1; i < numThreads; ++i) {
    double n2 = (double) threadData[i].length;
    setSize += threadData[i].length;
    double n =  (double) setSize;
        
    result += ((n2 - 1.0) / (n - 1.0)) * (threadData[i].result - result) - result / (n - 1.0);
  }
  return result;
}

/* 
static double mt_computeWeightedVarianceForKnownMean(misc_mt_manager_t restrict threadManager, const double* restrict x,
                                              size_t length, const double* restrict w, double mean)
{
  size_t numThreads, numValuesPerThread, offByOneIndex;
  misc_mt_getNumThreadsForJob(threadManager, length, WEIGHTED_VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD,
                             &numThreads, &numValuesPerThread, &offByOneIndex);
  
  if (numThreads <= 1) return computeWeightedVarianceForKnownMean(x, length, w, mean);
  
  WeightedVarianceForKnownMeanData threadData[numThreads];
  setupWeightedVarianceForKnownMeanData(threadData, numThreads, x, w, numValuesPerThread, offByOneIndex, &computeWeightedVarianceForKnownMean, mean);
  
  void* threadDataPtrs[numThreads];
  for (size_t i = 0; i < numThreads; ++i) threadDataPtrs[i] = (void*) &threadData[i];
  
  
  misc_mt_runTasks(threadManager, &computeWeightedVarianceForKnownMeanTask, threadDataPtrs, numThreads);
  
  return aggregateWeightedVarianceForKnownMeanData(threadData, numThreads);
}

static double mt_computeIndexedWeightedVarianceForKnownMean(misc_mt_manager_t restrict threadManager, const double* restrict x,
                                                     const size_t* restrict indices, size_t length, const double* restrict w, double mean)
{
  size_t numThreads, numValuesPerThread, offByOneIndex;
  misc_mt_getNumThreadsForJob(threadManager, length, INDEXED_WEIGHTED_VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD,
                             &numThreads, &numValuesPerThread, &offByOneIndex);
  
  if (numThreads <= 1) return computeIndexedWeightedVarianceForKnownMean(x, indices, length, w, mean);
  
  IndexedWeightedVarianceForKnownMeanData threadData[numThreads];
  setupIndexedWeightedVarianceForKnownMeanData(threadData, numThreads, x, indices, w, numValuesPerThread, offByOneIndex, &computeIndexedWeightedVarianceForKnownMean, mean);
  
  void* threadDataPtrs[numThreads];
  for (size_t i = 0; i < numThreads; ++i) threadDataPtrs[i] = (void*) &threadData[i];
  
  
  misc_mt_runTasks(threadManager, &computeIndexedWeightedVarianceForKnownMeanTask, threadDataPtrs, numThreads);
  
  return aggregateIndexedWeightedVarianceForKnownMeanData(threadData, numThreads);
} */

static double mt_computeUnrolledWeightedVarianceForKnownMean(misc_mt_manager_t restrict threadManager, const double* restrict x,
                                              size_t length, const double* restrict w, double mean)
{
  size_t numThreads, numValuesPerThread, offByOneIndex;
  misc_mt_getNumThreadsForJob(threadManager, length, UNROLLED_WEIGHTED_VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD,
                             &numThreads, &numValuesPerThread, &offByOneIndex);
  
  if (numThreads <= 1) return computeUnrolledWeightedVarianceForKnownMean(x, length, w, mean);
  
  WeightedVarianceForKnownMeanData threadData[numThreads];
  setupWeightedVarianceForKnownMeanData(threadData, numThreads, x, w, numValuesPerThread, offByOneIndex, computeUnrolledWeightedVarianceForKnownMean, mean);
  
  void* threadDataPtrs[numThreads];
  for (size_t i = 0; i < numThreads; ++i) threadDataPtrs[i] = (void*) &threadData[i];
  
  
  misc_mt_runTasks(threadManager, &computeWeightedVarianceForKnownMeanTask, threadDataPtrs, numThreads);
  
  return aggregateWeightedVarianceForKnownMeanData(threadData, numThreads);
}

static double mt_computeIndexedUnrolledWeightedVarianceForKnownMean(misc_mt_manager_t restrict threadManager, const double* restrict x,
                                                             const size_t* restrict indices, size_t length, const double* restrict w, double mean)
{
  size_t numThreads, numValuesPerThread, offByOneIndex;
  misc_mt_getNumThreadsForJob(threadManager, length, INDEXED_UNROLLED_WEIGHTED_VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD,
                             &numThreads, &numValuesPerThread, &offByOneIndex);
  
  if (numThreads <= 1) return computeIndexedUnrolledWeightedVarianceForKnownMean(x, indices, length, w, mean);
  
  IndexedWeightedVarianceForKnownMeanData threadData[numThreads];
  setupIndexedWeightedVarianceForKnownMeanData(threadData, numThreads, x, indices, w, numValuesPerThread, offByOneIndex, computeIndexedUnrolledWeightedVarianceForKnownMean, mean);
  
  void* threadDataPtrs[numThreads];
  for (size_t i = 0; i < numThreads; ++i) threadDataPtrs[i] = (void*) &threadData[i];
  
  
  misc_mt_runTasks(threadManager, &computeIndexedWeightedVarianceForKnownMeanTask, threadDataPtrs, numThreads);
  
  return aggregateIndexedWeightedVarianceForKnownMeanData(threadData, numThreads);
}

/*
static double mt_computeOnlineWeightedVarianceForKnownMean(misc_mt_manager_t restrict threadManager, const double* restrict x,
                                                      size_t length, const double* restrict w, double mean)
{
  size_t numThreads, numValuesPerThread, offByOneIndex;
  misc_mt_getNumThreadsForJob(threadManager, length, ONLINE_WEIGHTED_VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD,
                             &numThreads, &numValuesPerThread, &offByOneIndex);
  
  if (numThreads <= 1) return computeOnlineWeightedVarianceForKnownMean(x, length, w, mean);
  
  WeightedVarianceForKnownMeanData threadData[numThreads];
  setupWeightedVarianceForKnownMeanData(threadData, numThreads, x, w, numValuesPerThread, offByOneIndex, &computeOnlineWeightedVarianceForKnownMean, mean);
  
  void* threadDataPtrs[numThreads];
  for (size_t i = 0; i < numThreads; ++i) threadDataPtrs[i] = (void*) &threadData[i];
  
  
  misc_mt_runTasks(threadManager, &computeWeightedVarianceForKnownMeanTask, threadDataPtrs, numThreads);
  
  return aggregateWeightedVarianceForKnownMeanData(threadData, numThreads);
}

static double mt_computeIndexedOnlineWeightedVarianceForKnownMean(misc_mt_manager_t restrict threadManager, const double* restrict x,
                                                           const size_t* restrict indices, size_t length, const double* restrict w, double mean)
{
  size_t numThreads, numValuesPerThread, offByOneIndex;
  misc_mt_getNumThreadsForJob(threadManager, length, INDEXED_ONLINE_WEIGHTED_VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD,
                             &numThreads, &numValuesPerThread, &offByOneIndex);
  
  if (numThreads <= 1) return computeIndexedOnlineWeightedVarianceForKnownMean(x, indices, length, w, mean);
  
  IndexedWeightedVarianceForKnownMeanData threadData[numThreads];
  setupIndexedWeightedVarianceForKnownMeanData(threadData, numThreads, x, indices, w, numValuesPerThread, offByOneIndex, &computeIndexedOnlineWeightedVarianceForKnownMean, mean);
  
  void* threadDataPtrs[numThreads];
  for (size_t i = 0; i < numThreads; ++i) threadDataPtrs[i] = (void*) &threadData[i];
  
  
  misc_mt_runTasks(threadManager, &computeIndexedWeightedVarianceForKnownMeanTask, threadDataPtrs, numThreads);
  
  return aggregateIndexedWeightedVarianceForKnownMeanData(threadData, numThreads);
} */

static double mt_computeOnlineUnrolledWeightedVarianceForKnownMean(misc_mt_manager_t restrict threadManager, const double* restrict x,
                                                                   size_t length, const double* restrict w, double mean)
{
  size_t numThreads, numValuesPerThread, offByOneIndex;
  misc_mt_getNumThreadsForJob(threadManager, length, ONLINE_UNROLLED_WEIGHTED_VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD,
                             &numThreads, &numValuesPerThread, &offByOneIndex);
  
  if (numThreads <= 1) return computeOnlineUnrolledWeightedVarianceForKnownMean(x, length, w, mean);
  
  WeightedVarianceForKnownMeanData threadData[numThreads];
  setupWeightedVarianceForKnownMeanData(threadData, numThreads, x, w, numValuesPerThread, offByOneIndex, computeOnlineUnrolledWeightedVarianceForKnownMean, mean);
  
  void* threadDataPtrs[numThreads];
  for (size_t i = 0; i < numThreads; ++i) threadDataPtrs[i] = (void*) &threadData[i];
  
  
  misc_mt_runTasks(threadManager, &computeWeightedVarianceForKnownMeanTask, threadDataPtrs, numThreads);
  
  return aggregateWeightedVarianceForKnownMeanData(threadData, numThreads);
}

static double mt_computeIndexedOnlineUnrolledWeightedVarianceForKnownMean(misc_mt_manager_t restrict threadManager, const double* restrict x,
                                                                         const size_t* restrict indices, size_t length, const double* restrict w, double mean)
{
  size_t numThreads, numValuesPerThread, offByOneIndex;
  misc_mt_getNumThreadsForJob(threadManager, length, INDEXED_ONLINE_UNROLLED_WEIGHTED_VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD,
                             &numThreads, &numValuesPerThread, &offByOneIndex);
  
  if (numThreads <= 1) return computeIndexedOnlineUnrolledWeightedVarianceForKnownMean(x, indices, length, w, mean);
  
  IndexedWeightedVarianceForKnownMeanData threadData[numThreads];
  setupIndexedWeightedVarianceForKnownMeanData(threadData, numThreads, x, indices, w, numValuesPerThread, offByOneIndex, computeIndexedOnlineUnrolledWeightedVarianceForKnownMean, mean);
  
  void* threadDataPtrs[numThreads];
  for (size_t i = 0; i < numThreads; ++i) threadDataPtrs[i] = (void*) &threadData[i];
  
  
  misc_mt_runTasks(threadManager, &computeIndexedWeightedVarianceForKnownMeanTask, threadDataPtrs, numThreads);
  
  return aggregateIndexedWeightedVarianceForKnownMeanData(threadData, numThreads);
}




double misc_computeSumOfSquaredResiduals(const double* restrict x, size_t length, const double* restrict x_hat)
{
  if (length == 0) return 0.0;
  
  double result = 0.0;
  size_t lengthMod5 = length % 5;
  
  size_t i = 0;
  for ( /* */ ; i < lengthMod5; ++i) result += (x[i] - x_hat[i]) * (x[i] - x_hat[i]);
  
  for ( /* */ ; i < length; i += 5) {
    result += (x[i] - x_hat[i]) * (x[i] - x_hat[i]) + 
              (x[i + 1] - x_hat[i + 1]) * (x[i + 1] - x_hat[i + 1]) + 
              (x[i + 2] - x_hat[i + 2]) * (x[i + 2] - x_hat[i + 2]) +
              (x[i + 3] - x_hat[i + 3]) * (x[i + 3] - x_hat[i + 3]) +
              (x[i + 4] - x_hat[i + 4]) * (x[i + 4] - x_hat[i + 4]);
  }
  
  return result;
}

double misc_computeWeightedSumOfSquaredResiduals(const double* restrict x, size_t length, const double* restrict w, const double* restrict x_hat)
{
  if (length == 0) return 0.0;
  
  double result = 0.0;
  size_t lengthMod5 = length % 5;
  
  size_t i = 0;
  for ( /* */ ; i < lengthMod5; ++i) result += w[i] * (x[i] - x_hat[i]) * (x[i] - x_hat[i]);
  
  for ( /* */ ; i < length; i += 5) {
    result += w[i] * (x[i] - x_hat[i]) * (x[i] - x_hat[i]) + 
              w[i + 1] * (x[i + 1] - x_hat[i + 1]) * (x[i + 1] - x_hat[i + 1]) + 
              w[i + 2] * (x[i + 2] - x_hat[i + 2]) * (x[i + 2] - x_hat[i + 2]) +
              w[i + 3] * (x[i + 3] - x_hat[i + 3]) * (x[i + 3] - x_hat[i + 3]) +
              w[i + 4] * (x[i + 4] - x_hat[i + 4]) * (x[i + 4] - x_hat[i + 4]);
  }
  
  return result;
}

typedef double (*sumOfSquaredResidualsFunction)(const double* restrict x, size_t length, const double* restrict x_hat);
typedef double (*weightedSumOfSquaredResidualsFunction)(const double* restrict x, size_t length, const double* restrict w, const double* restrict x_hat);

typedef struct {
  const double* x;
  size_t length;
  const double* x_hat;
  double result;
  sumOfSquaredResidualsFunction function;
} SumOfSquaredResidualsData;

typedef struct {
  const double* x;
  size_t length;
  const double* w;
  const double* x_hat;
  double result;
  weightedSumOfSquaredResidualsFunction function;
} WeightedSumOfSquaredResidualsData;


static void computeSumOfSquaredResidualsTask(void* v_data)
{
  SumOfSquaredResidualsData* data = (SumOfSquaredResidualsData*) v_data;
  data->result = data->function(data->x, data->length, data->x_hat);
}

static void computeWeightedSumOfSquaredResidualsTask(void* v_data)
{
  WeightedSumOfSquaredResidualsData* data = (WeightedSumOfSquaredResidualsData*) v_data;
  data->result = data->function(data->x, data->length, data->w, data->x_hat);
}

static void setupSumOfSquaredResidualsData(SumOfSquaredResidualsData* restrict threadData, size_t numThreads,
                                           const double* restrict x, const double* restrict x_hat,
                                           size_t numValuesPerThread, size_t offByOneIndex, sumOfSquaredResidualsFunction function)
{
  size_t i = 0;
  for ( ; i < offByOneIndex; ++i) {
    threadData[i].x = x + i * numValuesPerThread;
    threadData[i].length = numValuesPerThread;
    threadData[i].x_hat = x_hat + i * numValuesPerThread;
    threadData[i].function = function;
  }
  for ( ; i < numThreads; ++i) {
    threadData[i].x = x + offByOneIndex * numValuesPerThread + (i - offByOneIndex) * (numValuesPerThread - 1);
    threadData[i].length = numValuesPerThread - 1;
    threadData[i].x_hat = x_hat + offByOneIndex * numValuesPerThread + (i - offByOneIndex) * (numValuesPerThread - 1);
    threadData[i].function = function;
  }
}

static void setupWeightedSumOfSquaredResidualsData(WeightedSumOfSquaredResidualsData* restrict threadData, size_t numThreads,
                                                   const double* restrict x,
                                                   const double* restrict w, const double* restrict x_hat,
                                                   size_t numValuesPerThread, size_t offByOneIndex, weightedSumOfSquaredResidualsFunction function)
{
  size_t i = 0;
  for ( ; i < offByOneIndex; ++i) {
    threadData[i].x = x + i * numValuesPerThread;
    threadData[i].length = numValuesPerThread;
    threadData[i].w = w + i * numValuesPerThread;
    threadData[i].x_hat = x_hat + i * numValuesPerThread;
    threadData[i].function = function;
  }
  for ( ; i < numThreads; ++i) {
    threadData[i].x = x + offByOneIndex * numValuesPerThread + (i - offByOneIndex) * (numValuesPerThread - 1);
    threadData[i].length = numValuesPerThread - 1;
    threadData[i].w = w + offByOneIndex * numValuesPerThread + (i - offByOneIndex) * (numValuesPerThread - 1);
    threadData[i].x_hat = x_hat + offByOneIndex * numValuesPerThread + (i - offByOneIndex) * (numValuesPerThread - 1);
    threadData[i].function = function;
  }
}

static double aggregateSumOfSquaredResidualsData(const SumOfSquaredResidualsData* threadData, size_t numThreads)
{
  double result = threadData[0].result;
  for (size_t i = 1; i < numThreads; ++i) result += threadData[i].result;
  return result;
}

static double aggregateWeightedSumOfSquaredResidualsData(const WeightedSumOfSquaredResidualsData* threadData, size_t numThreads)
{
  double result = threadData[0].result;
  for (size_t i = 1; i < numThreads; ++i) result += threadData[i].result;
  return result;
}


double misc_mt_computeSumOfSquaredResiduals(misc_mt_manager_t restrict threadManager, const double* restrict x, size_t length,
                                           const double* restrict x_hat)
{
  size_t numThreads, numValuesPerThread, offByOneIndex;
  misc_mt_getNumThreadsForJob(threadManager, length, SSR_MIN_NUM_VALUES_PER_THREAD, &numThreads, &numValuesPerThread, &offByOneIndex);
  
  if (numThreads <= 1) return misc_computeSumOfSquaredResiduals(x, length, x_hat);
  
  SumOfSquaredResidualsData threadData[numThreads];
  setupSumOfSquaredResidualsData(threadData, numThreads, x, x_hat, numValuesPerThread, offByOneIndex, &misc_computeSumOfSquaredResiduals);
  
  void* threadDataPtrs[numThreads];
  for (size_t i = 0; i < numThreads; ++i) threadDataPtrs[i] = (void*) &threadData[i];
  
  
  misc_mt_runTasks(threadManager, &computeSumOfSquaredResidualsTask, threadDataPtrs, numThreads);
  
  return aggregateSumOfSquaredResidualsData(threadData, numThreads);
}

double misc_mt_computeWeightedSumOfSquaredResiduals(misc_mt_manager_t restrict threadManager, const double* restrict x, size_t length,
                                                   const double* restrict w, const double* restrict x_hat)
{
  size_t numThreads, numValuesPerThread, offByOneIndex;
  misc_mt_getNumThreadsForJob(threadManager, length, WEIGHTED_SSR_MIN_NUM_VALUES_PER_THREAD, &numThreads, &numValuesPerThread, &offByOneIndex);
  
  if (numThreads <= 1) return misc_computeWeightedSumOfSquaredResiduals(x, length, w, x_hat);
  
  WeightedSumOfSquaredResidualsData threadData[numThreads];
  setupWeightedSumOfSquaredResidualsData(threadData, numThreads, x, w, x_hat, numValuesPerThread, offByOneIndex, &misc_computeWeightedSumOfSquaredResiduals);
  
  void* threadDataPtrs[numThreads];
  for (size_t i = 0; i < numThreads; ++i) threadDataPtrs[i] = (void*) &threadData[i];
  
  
  misc_mt_runTasks(threadManager, &computeWeightedSumOfSquaredResidualsTask, threadDataPtrs, numThreads);
  
  return aggregateWeightedSumOfSquaredResidualsData(threadData, numThreads);
}

double misc_htm_computeSumOfSquaredResiduals(misc_htm_manager_t restrict threadManager, size_t taskId, const double* restrict x, size_t length,
                                            const double* restrict x_hat)
{
  size_t numPieces, numValuesPerPiece, offByOneIndex;
  
  misc_htm_getNumPiecesForSubTask(threadManager, taskId, length, SSR_MIN_NUM_VALUES_PER_THREAD,
                                 &numPieces, &numValuesPerPiece, &offByOneIndex);
  
  if (numPieces <= 1) return misc_computeSumOfSquaredResiduals(x, length, x_hat);
  
  SumOfSquaredResidualsData data[numPieces];
  setupSumOfSquaredResidualsData(data, numPieces, x, x_hat, numValuesPerPiece, offByOneIndex, &misc_computeSumOfSquaredResiduals);
  
  void* dataPtrs[numPieces];
  for (size_t i = 0; i < numPieces; ++i) dataPtrs[i] = (void*) &data[i];
  
  misc_htm_runSubTask(threadManager, taskId, &computeSumOfSquaredResidualsTask, dataPtrs, numPieces);
  
  return aggregateSumOfSquaredResidualsData(data, numPieces);
}

double misc_htm_computeWeightedSumOfSquaredResiduals(misc_htm_manager_t restrict threadManager, size_t taskId,
                                                    const double* restrict x, size_t length,
                                                    const double* restrict w, const double* restrict x_hat)
{
  size_t numPieces, numValuesPerPiece, offByOneIndex;
  misc_htm_getNumPiecesForSubTask(threadManager, taskId, length, WEIGHTED_SSR_MIN_NUM_VALUES_PER_THREAD,
                                 &numPieces, &numValuesPerPiece, &offByOneIndex);
  
  if (numPieces <= 1) return misc_computeWeightedSumOfSquaredResiduals(x, length, w, x_hat);
  
  WeightedSumOfSquaredResidualsData data[numPieces];
  setupWeightedSumOfSquaredResidualsData(data, numPieces, x, w, x_hat, numValuesPerPiece, offByOneIndex, &misc_computeWeightedSumOfSquaredResiduals);
  
  void* dataPtrs[numPieces];
  for (size_t i = 0; i < numPieces; ++i) dataPtrs[i] = (void*) &data[i];
  
  misc_htm_runSubTask(threadManager, taskId,  &computeWeightedSumOfSquaredResidualsTask, dataPtrs, numPieces);
  
  return aggregateWeightedSumOfSquaredResidualsData(data, numPieces);
}

// forward declarations used in computeMean and computeIndexedMean
static double htm_computeUnrolledMean(misc_htm_manager_t restrict htm, size_t taskId, const double* restrict x, size_t length);
static double htm_computeIndexedUnrolledMean(misc_htm_manager_t restrict htm, size_t taskId, const double* restrict x, const size_t* restrict indices, size_t length);
static double htm_computeOnlineUnrolledMean(misc_htm_manager_t restrict threadManager, size_t taskId, const double* restrict x, size_t length);
static double htm_computeIndexedOnlineUnrolledMean(misc_htm_manager_t restrict threadManager, size_t taskId, const double* restrict x, const size_t* restrict indices, size_t length);


double misc_htm_computeMean(misc_htm_manager_t restrict threadManager, size_t taskId, const double* restrict x, size_t length)
{
  size_t numThreads = misc_htm_getNumThreadsForTopLevelTask(threadManager, taskId);
  size_t onlineCutoff = minimum(ONLINE_UNROLLED_MEAN_MIN_NUM_VALUES_PER_THREAD, ONLINE_CUTOFF);
  
  if (length / numThreads >= onlineCutoff) return htm_computeOnlineUnrolledMean(threadManager, taskId, x, length);
  return htm_computeUnrolledMean(threadManager, taskId, x, length);
}

double misc_htm_computeIndexedMean(misc_htm_manager_t restrict threadManager, size_t taskId, const double* restrict x, const size_t* restrict indices, size_t length)
{
  size_t numThreads = misc_htm_getNumThreadsForTopLevelTask(threadManager, taskId);
  size_t onlineCutoff = minimum(INDEXED_UNROLLED_MEAN_MIN_NUM_VALUES_PER_THREAD, INDEXED_ONLINE_CUTOFF);
  
  if (length / numThreads >= onlineCutoff) return htm_computeIndexedOnlineUnrolledMean(threadManager, taskId, x, indices, length);
  return htm_computeIndexedUnrolledMean(threadManager, taskId, x, indices, length);
}

static double htm_computeUnrolledMean(misc_htm_manager_t restrict threadManager, size_t taskId, const double* restrict x, size_t length)
{
  size_t numPieces, numValuesPerPiece, offByOneIndex;
  misc_htm_getNumPiecesForSubTask(threadManager, taskId, length, UNROLLED_MEAN_MIN_NUM_VALUES_PER_THREAD,
                                 &numPieces, &numValuesPerPiece, &offByOneIndex);
  
  if (numPieces <= 1) return computeUnrolledMean(x, length);
  
  MeanData data[numPieces];
  setupMeanData(data, numPieces, x, numValuesPerPiece, offByOneIndex, computeUnrolledMean);
  
  void* dataPtrs[numPieces];
  for (size_t i = 0; i < numPieces; ++i) dataPtrs[i] = (void*) &data[i];
  
  misc_htm_runSubTask(threadManager, taskId, &computeMeanTask, dataPtrs, numPieces);
  
  return aggregateMeanResults(data, numPieces);
}

static double htm_computeIndexedUnrolledMean(misc_htm_manager_t restrict threadManager, size_t taskId, const double* restrict x, const size_t* restrict indices, size_t length)
{
  size_t numPieces, numValuesPerPiece, offByOneIndex;
  misc_htm_getNumPiecesForSubTask(threadManager, taskId, length, INDEXED_UNROLLED_MEAN_MIN_NUM_VALUES_PER_THREAD,
                                 &numPieces, &numValuesPerPiece, &offByOneIndex);
  
  if (numPieces <= 1) return computeIndexedUnrolledMean(x, indices, length);
  
  IndexedMeanData data[numPieces];
  setupIndexedMeanData(data, numPieces, x, indices, numValuesPerPiece, offByOneIndex, computeIndexedUnrolledMean);
  
  void* dataPtrs[numPieces];
  for (size_t i = 0; i < numPieces; ++i) dataPtrs[i] = (void*) &data[i];
  
  misc_htm_runSubTask(threadManager, taskId, &computeIndexedMeanTask, dataPtrs, numPieces);
  
  return aggregateIndexedMeanResults(data, numPieces);
}

static double htm_computeOnlineUnrolledMean(misc_htm_manager_t restrict threadManager, size_t taskId, const double* restrict x, size_t length)
{
  size_t numPieces, numValuesPerPiece, offByOneIndex;
  misc_htm_getNumPiecesForSubTask(threadManager, taskId, length, ONLINE_UNROLLED_MEAN_MIN_NUM_VALUES_PER_THREAD,
                                 &numPieces, &numValuesPerPiece, &offByOneIndex);
  
  if (numPieces <= 1) return computeOnlineUnrolledMean(x, length);
  
  MeanData data[numPieces];
  setupMeanData(data, numPieces, x, numValuesPerPiece, offByOneIndex, computeOnlineUnrolledMean);
  
  void* dataPtrs[numPieces];
  for (size_t i = 0; i < numPieces; ++i) dataPtrs[i] = (void*) &data[i];
  
  misc_htm_runSubTask(threadManager, taskId, &computeMeanTask, dataPtrs, numPieces);
  
  return aggregateMeanResults(data, numPieces);
}

static double htm_computeIndexedOnlineUnrolledMean(misc_htm_manager_t restrict threadManager, size_t taskId, const double* restrict x, const size_t* restrict indices, size_t length)
{
  size_t numPieces, numValuesPerPiece, offByOneIndex;
  misc_htm_getNumPiecesForSubTask(threadManager, taskId, length, INDEXED_ONLINE_UNROLLED_MEAN_MIN_NUM_VALUES_PER_THREAD,
                                 &numPieces, &numValuesPerPiece, &offByOneIndex);
  
  if (numPieces <= 1) return computeIndexedOnlineUnrolledMean(x, indices, length);
  
  IndexedMeanData data[numPieces];
  setupIndexedMeanData(data, numPieces, x, indices, numValuesPerPiece, offByOneIndex, computeIndexedOnlineUnrolledMean);
  
  void* dataPtrs[numPieces];
  for (size_t i = 0; i < numPieces; ++i) dataPtrs[i] = (void*) &data[i];
  
  misc_htm_runSubTask(threadManager, taskId, &computeIndexedMeanTask, dataPtrs, numPieces);
  
  return aggregateIndexedMeanResults(data, numPieces);
}

// forward declarations used in computeWeightedMean and computeIndexedWeightedMean
static double htm_computeUnrolledWeightedMean(misc_htm_manager_t restrict htm, size_t taskId, const double* restrict x, size_t length, const double* restrict w, double* restrict nPtr);
static double htm_computeIndexedUnrolledWeightedMean(misc_htm_manager_t restrict htm, size_t taskId, const double* restrict x, const size_t* restrict indices, size_t length, const double* restrict w, double* restrict nPtr);
static double htm_computeOnlineUnrolledWeightedMean(misc_htm_manager_t restrict htm, size_t taskId, const double* restrict x, size_t length, const double* restrict w, double* restrict nPtr);
static double htm_computeIndexedOnlineUnrolledWeightedMean(misc_htm_manager_t restrict htm, size_t taskId, const double* restrict x, const size_t* restrict indices, size_t length, const double* restrict w, double* restrict nPtr);

double misc_htm_computeWeightedMean(misc_htm_manager_t restrict threadManager, size_t taskId, const double* restrict x,
                                   size_t length, const double* restrict w, double* restrict nPtr)
{
  size_t numThreads = misc_htm_getNumThreadsForTopLevelTask(threadManager, taskId);
  size_t onlineCutoff = minimum(UNROLLED_WEIGHTED_MEAN_MIN_NUM_VALUES_PER_THREAD, ONLINE_CUTOFF);
  
  if (length / numThreads >= onlineCutoff) return htm_computeOnlineUnrolledWeightedMean(threadManager, taskId, x, length, w, nPtr);
  return htm_computeUnrolledWeightedMean(threadManager, taskId, x, length, w, nPtr);
}

double misc_htm_computeIndexedWeightedMean(misc_htm_manager_t restrict threadManager, size_t taskId, const double* restrict x, 
                                          const size_t* restrict indices, size_t length, const double* restrict w, double* restrict nPtr)
{
  size_t numThreads = misc_htm_getNumThreadsForTopLevelTask(threadManager, taskId);
  size_t onlineCutoff = minimum(INDEXED_UNROLLED_WEIGHTED_MEAN_MIN_NUM_VALUES_PER_THREAD, ONLINE_CUTOFF);
  
  if (length / numThreads >= onlineCutoff) return htm_computeIndexedOnlineUnrolledWeightedMean(threadManager, taskId, x, indices, length, w, nPtr);
  return htm_computeIndexedUnrolledWeightedMean(threadManager, taskId, x, indices, length, w, nPtr);
}

static double htm_computeUnrolledWeightedMean(misc_htm_manager_t restrict threadManager, size_t taskId, const double* restrict x,
                                              size_t length, const double* restrict w, double* restrict nPtr)
{
  size_t numPieces, numValuesPerPiece, offByOneIndex;
  misc_htm_getNumPiecesForSubTask(threadManager, taskId, length, UNROLLED_WEIGHTED_MEAN_MIN_NUM_VALUES_PER_THREAD,
                                 &numPieces, &numValuesPerPiece, &offByOneIndex);
  
  if (numPieces <= 1) return computeUnrolledWeightedMean(x, length, w, nPtr);
  
  WeightedMeanData data[numPieces];
  setupWeightedMeanData(data, numPieces, x, w, numValuesPerPiece, offByOneIndex, computeUnrolledWeightedMean);
  
  void* dataPtrs[numPieces];
  for (size_t i = 0; i < numPieces; ++i) dataPtrs[i] = (void*) &data[i];
  
  misc_htm_runSubTask(threadManager, taskId, &computeWeightedMeanTask, dataPtrs, numPieces);
  
  return aggregateWeightedMeanResults(data, numPieces, nPtr);
}

static double htm_computeIndexedUnrolledWeightedMean(misc_htm_manager_t restrict threadManager, size_t taskId, const double* restrict x,
                                                     const size_t* restrict indices, size_t length,
                                                     const double* restrict w, double* restrict nPtr)
{
  size_t numPieces, numValuesPerPiece, offByOneIndex;
  misc_htm_getNumPiecesForSubTask(threadManager, taskId, length, INDEXED_UNROLLED_WEIGHTED_MEAN_MIN_NUM_VALUES_PER_THREAD,
                                 &numPieces, &numValuesPerPiece, &offByOneIndex);
  
  if (numPieces <= 1) return computeIndexedUnrolledWeightedMean(x, indices, length, w, nPtr);
  
  IndexedWeightedMeanData data[numPieces];
  setupIndexedWeightedMeanData(data, numPieces, x, indices, w, numValuesPerPiece, offByOneIndex, computeIndexedUnrolledWeightedMean);
  
  void* dataPtrs[numPieces];
  for (size_t i = 0; i < numPieces; ++i) dataPtrs[i] = (void*) &data[i];
  
  misc_htm_runSubTask(threadManager, taskId, &computeIndexedWeightedMeanTask, dataPtrs, numPieces);
  
  return aggregateIndexedWeightedMeanResults(data, numPieces, nPtr);
}

static double htm_computeOnlineUnrolledWeightedMean(misc_htm_manager_t restrict threadManager, size_t taskId, const double* restrict x,
                                                   size_t length, const double* restrict w, double* restrict nPtr)
{
  size_t numPieces, numValuesPerPiece, offByOneIndex;
  misc_htm_getNumPiecesForSubTask(threadManager, taskId, length, ONLINE_UNROLLED_WEIGHTED_MEAN_MIN_NUM_VALUES_PER_THREAD,
                                 &numPieces, &numValuesPerPiece, &offByOneIndex);
  
  if (numPieces <= 1) return computeOnlineUnrolledWeightedMean(x, length, w, nPtr);
  
  WeightedMeanData data[numPieces];
  setupWeightedMeanData(data, numPieces, x, w, numValuesPerPiece, offByOneIndex, computeOnlineUnrolledWeightedMean);
  
  void* dataPtrs[numPieces];
  for (size_t i = 0; i < numPieces; ++i) dataPtrs[i] = (void*) &data[i];
  
  misc_htm_runSubTask(threadManager, taskId, &computeWeightedMeanTask, dataPtrs, numPieces);
  
  return aggregateWeightedMeanResults(data, numPieces, nPtr);
}

static double htm_computeIndexedOnlineUnrolledWeightedMean(misc_htm_manager_t restrict threadManager, size_t taskId, const double* restrict x,
                                                           const size_t* restrict indices, size_t length,
                                                           const double* restrict w, double* restrict nPtr)
{
  size_t numPieces, numValuesPerPiece, offByOneIndex;
  misc_htm_getNumPiecesForSubTask(threadManager, taskId, length, INDEXED_ONLINE_UNROLLED_WEIGHTED_MEAN_MIN_NUM_VALUES_PER_THREAD,
                                 &numPieces, &numValuesPerPiece, &offByOneIndex);
  
  if (numPieces <= 1) return computeIndexedOnlineUnrolledWeightedMean(x, indices, length, w, nPtr);
  
  IndexedWeightedMeanData data[numPieces];
  setupIndexedWeightedMeanData(data, numPieces, x, indices, w, numValuesPerPiece, offByOneIndex, computeIndexedOnlineUnrolledWeightedMean);
  
  void* dataPtrs[numPieces];
  for (size_t i = 0; i < numPieces; ++i) dataPtrs[i] = (void*) &data[i];
  
  misc_htm_runSubTask(threadManager, taskId, &computeIndexedWeightedMeanTask, dataPtrs, numPieces);
  
  return aggregateIndexedWeightedMeanResults(data, numPieces, nPtr);
}

// forward declarations for variance/indexed var with known mean
static double htm_computeUnrolledVarianceForKnownMean(misc_htm_manager_t restrict htm, size_t taskId, const double* restrict x, size_t length, double mean);
static double htm_computeIndexedUnrolledVarianceForKnownMean(misc_htm_manager_t restrict htm, size_t taskId, const double* restrict x, const size_t* restrict indices, size_t length, double mean);
static double htm_computeOnlineUnrolledVarianceForKnownMean(misc_htm_manager_t restrict htm, size_t taskId, const double* restrict x, size_t length, double mean);
static double htm_computeIndexedOnlineUnrolledVarianceForKnownMean(misc_htm_manager_t restrict htm, size_t taskId, const double* restrict x, const size_t* restrict indices, size_t length, double mean);

double misc_htm_computeVarianceForKnownMean(misc_htm_manager_t restrict threadManager, size_t taskId, const double* restrict x,
                                           size_t length, double mean)
{
  size_t numThreads = misc_htm_getNumThreadsForTopLevelTask(threadManager, taskId);
  size_t onlineCutoff = minimum(UNROLLED_VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD, ONLINE_CUTOFF);
  
  if (length / numThreads >= onlineCutoff) return htm_computeOnlineUnrolledVarianceForKnownMean(threadManager, taskId, x, length, mean);
  return htm_computeUnrolledVarianceForKnownMean(threadManager, taskId, x, length, mean);
}

double misc_htm_computeIndexedVarianceForKnownMean(misc_htm_manager_t restrict threadManager, size_t taskId, const double* restrict x,
                                                  const size_t* restrict indices, size_t length, double mean)
{
  size_t numThreads = misc_htm_getNumThreadsForTopLevelTask(threadManager, taskId);
  size_t onlineCutoff = minimum(INDEXED_UNROLLED_VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD, INDEXED_ONLINE_CUTOFF);
  
  if (length / numThreads >= onlineCutoff) return htm_computeIndexedOnlineUnrolledVarianceForKnownMean(threadManager, taskId, x, indices, length, mean);
  return htm_computeIndexedUnrolledVarianceForKnownMean(threadManager, taskId, x, indices, length, mean);
}

static double htm_computeUnrolledVarianceForKnownMean(misc_htm_manager_t restrict threadManager, size_t taskId, const double* restrict x,
                                                      size_t length, double mean)
{
  size_t numPieces, numValuesPerPiece, offByOneIndex;
  misc_htm_getNumPiecesForSubTask(threadManager, taskId, length, UNROLLED_VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD,
                                 &numPieces, &numValuesPerPiece, &offByOneIndex);
  
  if (numPieces <= 1) return computeUnrolledVarianceForKnownMean(x, length, mean);
  
  VarianceForKnownMeanData data[numPieces];
  setupVarianceForKnownMeanData(data, numPieces, x, numValuesPerPiece, offByOneIndex, computeUnrolledVarianceForKnownMean, mean);
  
  void* dataPtrs[numPieces];
  for (size_t i = 0; i < numPieces; ++i) dataPtrs[i] = (void*) &data[i];
  
  misc_htm_runSubTask(threadManager, taskId, &computeVarianceForKnownMeanTask, dataPtrs, numPieces);
  
  return aggregateVarianceForKnownMeanData(data, numPieces);
}


static double htm_computeOnlineUnrolledVarianceForKnownMean(misc_htm_manager_t restrict threadManager, size_t taskId, const double* restrict x,
                                                            size_t length, double mean)
{
  size_t numPieces, numValuesPerPiece, offByOneIndex;
  misc_htm_getNumPiecesForSubTask(threadManager, taskId, length, ONLINE_UNROLLED_VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD,
                                 &numPieces, &numValuesPerPiece, &offByOneIndex);
  
  if (numPieces <= 1) return computeOnlineUnrolledVarianceForKnownMean(x, length, mean);
  
  VarianceForKnownMeanData data[numPieces];
  setupVarianceForKnownMeanData(data, numPieces, x, numValuesPerPiece, offByOneIndex, computeOnlineUnrolledVarianceForKnownMean, mean);
  
  void* dataPtrs[numPieces];
  for (size_t i = 0; i < numPieces; ++i) dataPtrs[i] = (void*) &data[i];
  
  misc_htm_runSubTask(threadManager, taskId, &computeVarianceForKnownMeanTask, dataPtrs, numPieces);
  
  return aggregateVarianceForKnownMeanData(data, numPieces);
}

static double htm_computeIndexedUnrolledVarianceForKnownMean(misc_htm_manager_t restrict threadManager, size_t taskId, const double* restrict x,
                                                             const size_t* restrict indices, size_t length, double mean)
{
  size_t numPieces, numValuesPerPiece, offByOneIndex;
  misc_htm_getNumPiecesForSubTask(threadManager, taskId, length, INDEXED_UNROLLED_VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD,
                                 &numPieces, &numValuesPerPiece, &offByOneIndex);
  
  if (numPieces <= 1) return computeIndexedUnrolledVarianceForKnownMean(x, indices, length, mean);
  
  IndexedVarianceForKnownMeanData data[numPieces];
  setupIndexedVarianceForKnownMeanData(data, numPieces, x, indices, numValuesPerPiece, offByOneIndex, computeIndexedUnrolledVarianceForKnownMean, mean);
  
  void* dataPtrs[numPieces];
  for (size_t i = 0; i < numPieces; ++i) dataPtrs[i] = (void*) &data[i];
  
  misc_htm_runSubTask(threadManager, taskId, &computeIndexedVarianceForKnownMeanTask, dataPtrs, numPieces);
  
  return aggregateIndexedVarianceForKnownMeanData(data, numPieces);
}

static double htm_computeIndexedOnlineUnrolledVarianceForKnownMean(misc_htm_manager_t restrict threadManager, size_t taskId, const double* restrict x,
                                                                   const size_t* restrict indices, size_t length, double mean)
{
  size_t numPieces, numValuesPerPiece, offByOneIndex;
  misc_htm_getNumPiecesForSubTask(threadManager, taskId, length, INDEXED_ONLINE_UNROLLED_VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD,
                                 &numPieces, &numValuesPerPiece, &offByOneIndex);
  
  if (numPieces <= 1) return computeIndexedOnlineUnrolledVarianceForKnownMean(x, indices, length, mean);
  
  IndexedVarianceForKnownMeanData data[numPieces];
  setupIndexedVarianceForKnownMeanData(data, numPieces, x, indices, numValuesPerPiece, offByOneIndex, computeIndexedOnlineUnrolledVarianceForKnownMean, mean);
  
  void* dataPtrs[numPieces];
  for (size_t i = 0; i < numPieces; ++i) dataPtrs[i] = (void*) &data[i];
  
  misc_htm_runSubTask(threadManager, taskId, &computeIndexedVarianceForKnownMeanTask, dataPtrs, numPieces);
  
  return aggregateIndexedVarianceForKnownMeanData(data, numPieces);
}

// forward decs
static double htm_computeUnrolledWeightedVarianceForKnownMean(
  misc_htm_manager_t restrict htm, size_t taskId, const double* restrict x, size_t length,
  const double* restrict w, double mean);
static double htm_computeOnlineUnrolledWeightedVarianceForKnownMean(
  misc_htm_manager_t restrict htm, size_t taskId, const double* restrict x, size_t length,
  const double* restrict w, double mean);
static double htm_computeIndexedUnrolledWeightedVarianceForKnownMean(
  misc_htm_manager_t restrict htm, size_t taskId, const double* restrict x,
  const size_t* restrict indices, size_t length,
  const double* restrict w, double mean);
static double htm_computeIndexedOnlineUnrolledWeightedVarianceForKnownMean(
  misc_htm_manager_t restrict htm, size_t taskId, const double* restrict x,
  const size_t* restrict indices, size_t length,
  const double* restrict w, double mean);


double misc_htm_computeWeightedVarianceForKnownMean(
  misc_htm_manager_t restrict threadManager, size_t taskId, const double* restrict x,
  size_t length, const double* restrict w, double mean)
{
  size_t numThreads = misc_htm_getNumThreadsForTopLevelTask(threadManager, taskId);
  size_t onlineCutoff = minimum(UNROLLED_WEIGHTED_VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD, ONLINE_CUTOFF);
  
  if (length / numThreads >= onlineCutoff) return htm_computeOnlineUnrolledWeightedVarianceForKnownMean(threadManager, taskId, x, length, w, mean);
  return htm_computeUnrolledWeightedVarianceForKnownMean(threadManager, taskId, x, length, w, mean);
}

double misc_htm_computeIndexedWeightedVarianceForKnownMean(
  misc_htm_manager_t restrict threadManager, size_t taskId, const double* restrict x,
  const size_t* restrict indices, size_t length,
  const double* restrict w, double mean)
{
  size_t numThreads = misc_htm_getNumThreadsForTopLevelTask(threadManager, taskId);
  size_t onlineCutoff = minimum(INDEXED_UNROLLED_WEIGHTED_VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD, ONLINE_CUTOFF);
  
  if (length / numThreads >= onlineCutoff) return htm_computeIndexedOnlineUnrolledWeightedVarianceForKnownMean(threadManager, taskId, x, indices, length, w, mean);
  return htm_computeIndexedUnrolledWeightedVarianceForKnownMean(threadManager, taskId, x, indices, length, w, mean);
}

static double htm_computeUnrolledWeightedVarianceForKnownMean(
  misc_htm_manager_t restrict threadManager, size_t taskId, const double* restrict x, size_t length,
  const double* restrict w, double mean)
{
  size_t numPieces, numValuesPerPiece, offByOneIndex;
  misc_htm_getNumPiecesForSubTask(threadManager, taskId, length, UNROLLED_WEIGHTED_VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD,
                                 &numPieces, &numValuesPerPiece, &offByOneIndex);
  
  if (numPieces <= 1) return computeUnrolledWeightedVarianceForKnownMean(x, length, w, mean);
  
  WeightedVarianceForKnownMeanData data[numPieces];
  setupWeightedVarianceForKnownMeanData(data, numPieces, x, w, numValuesPerPiece, offByOneIndex, computeUnrolledWeightedVarianceForKnownMean, mean);
  
  void* dataPtrs[numPieces];
  for (size_t i = 0; i < numPieces; ++i) dataPtrs[i] = (void*) &data[i];
  
  misc_htm_runSubTask(threadManager, taskId, &computeWeightedVarianceForKnownMeanTask, dataPtrs, numPieces);
  
  return aggregateWeightedVarianceForKnownMeanData(data, numPieces);
}

static double htm_computeOnlineUnrolledWeightedVarianceForKnownMean(
  misc_htm_manager_t restrict threadManager, size_t taskId, const double* restrict x, size_t length,
  const double* restrict w, double mean)
{
  size_t numPieces, numValuesPerPiece, offByOneIndex;
  misc_htm_getNumPiecesForSubTask(threadManager, taskId, length, ONLINE_UNROLLED_WEIGHTED_VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD,
                                 &numPieces, &numValuesPerPiece, &offByOneIndex);
  
  if (numPieces <= 1) return computeOnlineUnrolledWeightedVarianceForKnownMean(x, length, w, mean);
  
  WeightedVarianceForKnownMeanData data[numPieces];
  setupWeightedVarianceForKnownMeanData(data, numPieces, x, w, numValuesPerPiece, offByOneIndex, computeOnlineUnrolledWeightedVarianceForKnownMean, mean);
  
  void* dataPtrs[numPieces];
  for (size_t i = 0; i < numPieces; ++i) dataPtrs[i] = (void*) &data[i];
  
  misc_htm_runSubTask(threadManager, taskId, &computeWeightedVarianceForKnownMeanTask, dataPtrs, numPieces);
  
  return aggregateWeightedVarianceForKnownMeanData(data, numPieces);
}


static double htm_computeIndexedUnrolledWeightedVarianceForKnownMean(
  misc_htm_manager_t restrict threadManager, size_t taskId, const double* restrict x,
  const size_t* restrict indices, size_t length,
  const double* restrict w, double mean)
{
  size_t numPieces, numValuesPerPiece, offByOneIndex;
  misc_htm_getNumPiecesForSubTask(threadManager, taskId, length, INDEXED_UNROLLED_WEIGHTED_VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD,
                                 &numPieces, &numValuesPerPiece, &offByOneIndex);
  
  if (numPieces <= 1) return computeIndexedUnrolledWeightedVarianceForKnownMean(x, indices, length, w, mean);
  
  IndexedWeightedVarianceForKnownMeanData data[numPieces];
  setupIndexedWeightedVarianceForKnownMeanData(data, numPieces, x, indices, w, numValuesPerPiece, offByOneIndex, computeIndexedUnrolledWeightedVarianceForKnownMean, mean);
  
  void* dataPtrs[numPieces];
  for (size_t i = 0; i < numPieces; ++i) dataPtrs[i] = (void*) &data[i];
  
  misc_htm_runSubTask(threadManager, taskId, &computeIndexedWeightedVarianceForKnownMeanTask, dataPtrs, numPieces);
  
  return aggregateIndexedWeightedVarianceForKnownMeanData(data, numPieces);
}

static double htm_computeIndexedOnlineUnrolledWeightedVarianceForKnownMean(
  misc_htm_manager_t restrict threadManager, size_t taskId, const double* restrict x,
  const size_t* restrict indices, size_t length,
  const double* restrict w, double mean)
{
  size_t numPieces, numValuesPerPiece, offByOneIndex;
  misc_htm_getNumPiecesForSubTask(threadManager, taskId, length, INDEXED_ONLINE_UNROLLED_WEIGHTED_VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD,
                                 &numPieces, &numValuesPerPiece, &offByOneIndex);
  
  if (numPieces <= 1) return computeIndexedOnlineUnrolledWeightedVarianceForKnownMean(x, indices, length, w, mean);
  
  IndexedWeightedVarianceForKnownMeanData data[numPieces];
  setupIndexedWeightedVarianceForKnownMeanData(data, numPieces, x, indices, w, numValuesPerPiece, offByOneIndex, computeIndexedOnlineUnrolledWeightedVarianceForKnownMean, mean);
  
  void* dataPtrs[numPieces];
  for (size_t i = 0; i < numPieces; ++i) dataPtrs[i] = (void*) &data[i];
  
  misc_htm_runSubTask(threadManager, taskId, &computeIndexedWeightedVarianceForKnownMeanTask, dataPtrs, numPieces);
  
  return aggregateIndexedWeightedVarianceForKnownMeanData(data, numPieces);
}

