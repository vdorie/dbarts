#include <external/stats.h>
#include <external/stats_mt.h>
#include <external/thread.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
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
#define                 MEAN_MIN_NUM_VALUES_PER_THREAD 100000
#define        UNROLLED_MEAN_MIN_NUM_VALUES_PER_THREAD 200000
#define          ONLINE_MEAN_MIN_NUM_VALUES_PER_THREAD 7500
#define ONLINE_UNROLLED_MEAN_MIN_NUM_VALUES_PER_THREAD 25000

#define                 VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD 75000
#define        UNROLLED_VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD 75000
#define          ONLINE_VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD 7500
#define ONLINE_UNROLLED_VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD 25000

#define                 VAR_MIN_NUM_VALUES_PER_THREAD 75000
#define        UNROLLED_VAR_MIN_NUM_VALUES_PER_THREAD 75000
#define          ONLINE_VAR_MIN_NUM_VALUES_PER_THREAD 5000
#define ONLINE_UNROLLED_VAR_MIN_NUM_VALUES_PER_THREAD 12000

#define                 INDEXED_MEAN_MIN_NUM_VALUES_PER_THREAD 100000
#define        INDEXED_UNROLLED_MEAN_MIN_NUM_VALUES_PER_THREAD 100000
#define          INDEXED_ONLINE_MEAN_MIN_NUM_VALUES_PER_THREAD 7500
#define INDEXED_ONLINE_UNROLLED_MEAN_MIN_NUM_VALUES_PER_THREAD 25000

#define                 INDEXED_VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD 75000
#define        INDEXED_UNROLLED_VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD 75000
#define          INDEXED_ONLINE_VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD 7500
#define INDEXED_ONLINE_UNROLLED_VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD 25000

#define                 INDEXED_VAR_MIN_NUM_VALUES_PER_THREAD 75000
#define        INDEXED_UNROLLED_VAR_MIN_NUM_VALUES_PER_THREAD 75000
#define          INDEXED_ONLINE_VAR_MIN_NUM_VALUES_PER_THREAD 5000
#define INDEXED_ONLINE_UNROLLED_VAR_MIN_NUM_VALUES_PER_THREAD 12000

// if < value, calculate straightup; otherwise use online
// algorithm to reduce round-off err
#define ONLINE_CUTOFF 10000
#define INDEXED_ONLINE_CUTOFF 10000

// various implementations
static double computeMean(const double* x, size_t length);
static double computeUnrolledMean(const double* x, size_t length);
static double computeOnlineMean(const double* x, size_t length);
static double computeOnlineUnrolledMean(const double* x, size_t length);

static double computeVarianceForKnownMean(const double* x, size_t length, double mean);
static double computeUnrolledVarianceForKnownMean(const double* x, size_t length, double mean);
static double computeOnlineVarianceForKnownMean(const double* x, size_t length, double mean);
static double computeOnlineUnrolledVarianceForKnownMean(const double* x, size_t length, double mean);

static double computeVariance(const double* restrict x, size_t length, double* restrict meanPtr);
static double computeUnrolledVariance(const double* restrict x, size_t length, double* restrict meanPtr);
static double computeOnlineVariance(const double* restrict x, size_t length, double* restrict meanPtr);
static double computeOnlineUnrolledVariance(const double* restrict x, size_t length, double* restrict meanPtr);

static double computeIndexedMean(const double* restrict x, const size_t* restrict indices, size_t length);
static double computeIndexedUnrolledMean(const double* restrict x, const size_t* restrict indices, size_t length);
static double computeIndexedOnlineMean(const double* restrict x, const size_t* restrict indices, size_t length);
static double computeIndexedOnlineUnrolledMean(const double* restrict x, const size_t* restrict indices, size_t length);

static double computeIndexedVarianceForKnownMean(const double* restrict x, const size_t* restrict indices, size_t length, double mean);
static double computeIndexedUnrolledVarianceForKnownMean(const double* restrict x, const size_t* restrict indices, size_t length, double mean);
static double computeIndexedOnlineVarianceForKnownMean(const double* restrict x, const size_t* restrict indices, size_t length, double mean);
static double computeIndexedOnlineUnrolledVarianceForKnownMean(const double* restrict x, const size_t* restrict indices, size_t length, double mean);

static double computeIndexedVariance(const double* restrict x, const size_t* restrict indices, size_t length, double* restrict meanPtr);
static double computeIndexedUnrolledVariance(const double* restrict x, const size_t* restrict indices, size_t length, double* restrict meanPtr);
static double computeIndexedOnlineVariance(const double* restrict x, const size_t* restrict indices, size_t length, double* restrict meanPtr);
static double computeIndexedOnlineUnrolledVariance(const double* restrict x,const size_t* restrict indices, size_t length, double* restrict meanPtr);

UNUSED static double mt_computeMean(ext_mt_manager_t restrict threadManager, const double* restrict x, size_t length);
static double mt_computeUnrolledMean(ext_mt_manager_t restrict threadManager, const double* restrict x, size_t length);
UNUSED static double mt_computeOnlineMean(ext_mt_manager_t restrict threadManager, const double* restrict x, size_t length);
static double mt_computeOnlineUnrolledMean(ext_mt_manager_t restrict threadManager, const double* restrict x, size_t length);

UNUSED static double mt_computeVarianceForKnownMean(ext_mt_manager_t restrict threadManager, const double* restrict x, size_t length, double mean);
static double mt_computeUnrolledVarianceForKnownMean(ext_mt_manager_t restrict threadManager, const double* restrict x, size_t length, double mean);
UNUSED static double mt_computeOnlineVarianceForKnownMean(ext_mt_manager_t restrict threadManager, const double* restrict x, size_t length, double mean);
static double mt_computeOnlineUnrolledVarianceForKnownMean(ext_mt_manager_t restrict threadManager, const double* restrict x, size_t length, double mean);

UNUSED static double mt_computeVariance(ext_mt_manager_t restrict threadManager, const double* restrict x, size_t length, double* restrict meanPtr);
static double mt_computeUnrolledVariance(ext_mt_manager_t restrict threadManager, const double* restrict x, size_t length, double* restrict meanPtr);
UNUSED static double mt_computeOnlineVariance(ext_mt_manager_t restrict threadManager, const double* restrict x, size_t length, double* restrict meanPtr);
static double mt_computeOnlineUnrolledVariance(ext_mt_manager_t restrict threadManager, const double* restrict x, size_t length, double* restrict meanPtr);

UNUSED static double mt_computeIndexedMean(ext_mt_manager_t restrict threadManager, const double* restrict x, const size_t* indices, size_t length);
static double mt_computeIndexedUnrolledMean(ext_mt_manager_t restrict threadManager, const double* restrict x, const size_t* indices, size_t length);
UNUSED static double mt_computeIndexedOnlineMean(ext_mt_manager_t restrict threadManager, const double* restrict x, const size_t* indices, size_t length);
static double mt_computeIndexedOnlineUnrolledMean(ext_mt_manager_t restrict threadManager, const double* restrict x, const size_t* indices, size_t length);

UNUSED static double mt_computeIndexedVarianceForKnownMean(ext_mt_manager_t restrict threadManager, const double* restrict x, const size_t* restrict indices, size_t length, double mean);
static double mt_computeIndexedUnrolledVarianceForKnownMean(ext_mt_manager_t restrict threadManager, const double* restrict x, const size_t* restrict indices, size_t length, double mean);
UNUSED static double mt_computeIndexedOnlineVarianceForKnownMean(ext_mt_manager_t restrict threadManager, const double* restrict x, const size_t* restrict indices, size_t length, double mean);
static double mt_computeIndexedOnlineUnrolledVarianceForKnownMean(ext_mt_manager_t restrict threadManager, const double* restrict x, const size_t* restrict indices, size_t length, double mean);

UNUSED static double mt_computeIndexedVariance(ext_mt_manager_t restrict threadManager, const double* restrict x, const size_t* restrict indices, size_t length, double* restrict mean);
static double mt_computeIndexedUnrolledVariance(ext_mt_manager_t restrict threadManager, const double* restrict x, const size_t* restrict indices, size_t length, double* restrict mean);
UNUSED static double mt_computeIndexedOnlineVariance(ext_mt_manager_t restrict threadManager, const double* restrict x, const size_t* restrict indices, size_t length, double* restrict mean);
static double mt_computeIndexedOnlineUnrolledVariance(ext_mt_manager_t restrict threadManager, const double* restrict x, const size_t* restrict indices, size_t length, double* restrict mean);

// interface functions that dispatch to the workers
double ext_computeMean(const double* x, size_t length)
{
  if (length > ONLINE_CUTOFF) return computeOnlineUnrolledMean(x, length);
  return computeUnrolledMean(x, length);
}

double ext_computeIndexedMean(const double* restrict x, const size_t* restrict indices, size_t length)
{
    if (length > ONLINE_CUTOFF) return computeIndexedOnlineUnrolledMean(x, indices, length);
  return computeIndexedUnrolledMean(x, indices, length);
}

// The two-pass is faster when using online algorithms, one-pass when not
double ext_computeVariance(const double* restrict x, size_t length, double* restrict meanPtr)
{
  if (length > ONLINE_CUTOFF) {
    double mean = computeOnlineUnrolledMean(x, length);
    if (meanPtr != NULL) *meanPtr = mean;
    return computeOnlineUnrolledVarianceForKnownMean(x, length, mean);
  }
  
  return computeUnrolledVariance(x, length, meanPtr);
}

// one pass pretty much always faster when indexed
double ext_computeIndexedVariance(const double* restrict x, const size_t* restrict indices, size_t length, double* restrict meanPtr)
{
  if (length > ONLINE_CUTOFF) return computeIndexedOnlineUnrolledVariance(x, indices, length, meanPtr);
  
  return computeIndexedUnrolledVariance(x, indices, length, meanPtr);
}

double ext_computeVarianceForKnownMean(const double* x, size_t length, double mean)
{
  if (length > ONLINE_CUTOFF) return computeOnlineUnrolledVarianceForKnownMean(x, length, mean);
  return computeUnrolledVarianceForKnownMean(x, length, mean);
}

double ext_computeIndexedVarianceForKnownMean(const double* restrict x, const ext_size_t* restrict indices, size_t length, double mean)
{
  if (length > ONLINE_CUTOFF) computeIndexedUnrolledVarianceForKnownMean(x, indices, length, mean);
  return computeIndexedOnlineUnrolledVarianceForKnownMean(x, indices, length, mean);
}



#define minimum(_A_, _B_) ((_A_) < (_B_) ? (_A_) : (_B_))

// if the data for any thread would, by itself, trigger a fall-back to single threaded
// and that single-threaded function equiv would prefer the non-online version, do that instead
double ext_mt_computeMean(ext_mt_manager_t restrict threadManager, const double* restrict x, size_t length)
{
  size_t numThreads = ext_mt_getNumThreads(threadManager);
  size_t onlineCutoff = minimum(ONLINE_UNROLLED_MEAN_MIN_NUM_VALUES_PER_THREAD, ONLINE_CUTOFF);
  
  if (length / numThreads >= onlineCutoff) return mt_computeOnlineUnrolledMean(threadManager, x, length);
  return mt_computeUnrolledMean(threadManager, x, length);
}

double ext_mt_computeIndexedMean(ext_mt_manager_t restrict threadManager, const double* restrict x, const size_t* restrict indices, size_t length)
{
  size_t numThreads = ext_mt_getNumThreads(threadManager);
  size_t onlineCutoff = minimum(INDEXED_UNROLLED_MEAN_MIN_NUM_VALUES_PER_THREAD, INDEXED_ONLINE_CUTOFF);
  
  if (length / numThreads >= onlineCutoff) return mt_computeIndexedOnlineUnrolledMean(threadManager, x, indices, length);
  return mt_computeIndexedUnrolledMean(threadManager, x, indices, length);
}

double ext_mt_computeVariance(ext_mt_manager_t restrict threadManager, const double* restrict x, size_t length, double* meanPtr)
{
  size_t numThreads = ext_mt_getNumThreads(threadManager);
  size_t onlineCutoff = minimum(UNROLLED_VAR_MIN_NUM_VALUES_PER_THREAD, ONLINE_CUTOFF);
  
  if (length / numThreads >= onlineCutoff) return mt_computeOnlineUnrolledVariance(threadManager, x, length, meanPtr);
  return mt_computeUnrolledVariance(threadManager, x, length, meanPtr);
}

double ext_mt_computeIndexedVariance(ext_mt_manager_t restrict threadManager, const double* restrict x, const size_t* restrict indices, size_t length, double* meanPtr)
{
  size_t numThreads = ext_mt_getNumThreads(threadManager);
  size_t onlineCutoff = minimum(INDEXED_UNROLLED_VAR_MIN_NUM_VALUES_PER_THREAD, INDEXED_ONLINE_CUTOFF);
  
  if (length / numThreads >= onlineCutoff) return mt_computeIndexedOnlineUnrolledVariance(threadManager, x, indices, length, meanPtr);
  return mt_computeIndexedUnrolledVariance(threadManager, x, indices, length, meanPtr);
}

double ext_mt_computeVarianceForKnownMean(ext_mt_manager_t restrict threadManager, const double* restrict x, size_t length, double mean)
{
  size_t numThreads = ext_mt_getNumThreads(threadManager);
  size_t onlineCutoff = minimum(UNROLLED_VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD, ONLINE_CUTOFF);
  
  if (length / numThreads >= onlineCutoff) return mt_computeOnlineUnrolledVarianceForKnownMean(threadManager, x, length, mean);
  return mt_computeUnrolledVarianceForKnownMean(threadManager, x, length, mean);
}

double ext_mt_computeIndexedVarianceForKnownMean(ext_mt_manager_t restrict threadManager, const double* restrict x, const size_t* restrict indices, size_t length, double mean)
{
  size_t numThreads = ext_mt_getNumThreads(threadManager);
  size_t onlineCutoff = minimum(INDEXED_UNROLLED_VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD, INDEXED_ONLINE_CUTOFF);
  
  if (length / numThreads >= onlineCutoff) return mt_computeIndexedOnlineUnrolledVarianceForKnownMean(threadManager, x, indices, length, mean);
  return mt_computeIndexedUnrolledVarianceForKnownMean(threadManager, x, indices, length, mean);
}


// work-horse implementations

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

static double computeUnrolledMean(const double* x, size_t length)
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

static double computeIndexedUnrolledMean(const double* restrict x, const size_t* restrict indices, size_t length)
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

static double computeOnlineUnrolledMean(const double* x, size_t length)
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

static double computeIndexedOnlineUnrolledMean(const double* restrict x, const size_t* restrict indices, size_t length)
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

static double computeUnrolledVarianceForKnownMean(const double* x, size_t length, double mean)
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
    result += (x[i] - mean) * (x[i] - mean) +
              (x[i + 1] - mean) * (x[i + 1] - mean) +
              (x[i + 2] - mean) * (x[i + 2] - mean) +
              (x[i + 3] - mean) * (x[i + 3] - mean) +
              (x[i + 4] - mean) * (x[i + 4] - mean);
  }
  
  return result / (double) (length - 1);
}

static double computeIndexedUnrolledVarianceForKnownMean(const double* restrict x, const size_t* restrict indices, size_t length, double mean)
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

static double computeOnlineUnrolledVarianceForKnownMean(const double* x, size_t length, double mean)
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

static double computeIndexedOnlineUnrolledVarianceForKnownMean(const double* restrict x, const size_t* restrict indices, size_t length, double mean)
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
  
  /* size_t i = 1;
  size_t lengthMod4 = (length - 1) % 4;
  
  if (lengthMod4++ != 0) {
    for ( ; i < lengthMod4; ++i) {
      double dev = x[indices[i]] - mean;
      mean += dev / (double) (i + 1);
      var  += (dev * (x[indices[i]] - mean) - var) / (double) (i + 1);
    }
    if (length < 5) { if (meanPtr != NULL) *meanPtr = mean; return nScale * var; }
  }
  
  for ( ; i < length; i += 4) {
#define n1 ((double) i)
#define n2 4.0
#define n  ((double) (i + 4))
    
    double meanNext4 = (x[indices[i]] + x[indices[i + 1]] + x[indices[i + 2]] + x[indices[i + 3]]) / n2;
    double varNext4  = ((x[indices[i]] - meanNext4) * (x[indices[i]] - meanNext4) +
                        (x[indices[i + 1]] - meanNext4) * (x[indices[i + 1]] - meanNext4) +
                        (x[indices[i + 2]] - meanNext4) * (x[indices[i + 2]] - meanNext4) +
                        (x[indices[i + 3]] - meanNext4) * (x[indices[i + 3]] - meanNext4)) / n2;
    
    var  += (n2 / n) * (varNext4 - var) + ((meanNext4 - mean) * (n1 / n)) * ((meanNext4 - mean) * (n2 / n));
    mean += (n2 / n) * (meanNext4 - mean);
#undef n1
#undef n2
#undef n
  } */
  
  if (meanPtr != NULL) *meanPtr = mean;
  return nScale * var;
}





// below this, multithreaded madness

typedef double (*meanFunction)(const double* x, size_t length);

typedef struct {
  const double* x;
  size_t length;
  double result;
  meanFunction function;
} MeanData;

static void computeMeanTask(void* v_data)
{
  MeanData* data = (MeanData*) v_data;
  data->result = data->function(data->x, data->length);
}

static void setupMeanData(MeanData* restrict threadData, size_t numThreads, const double* restrict x, size_t length,
                          size_t numValuesPerThread, meanFunction function)
{
  size_t i = 0;
  for ( ; i < numThreads - 1; ++i) {
    threadData[i].x = x + i * numValuesPerThread;
    threadData[i].length = numValuesPerThread;
    threadData[i].function = function;
  }
  threadData[i].x = x + i * numValuesPerThread;
  threadData[i].length = length - i * numValuesPerThread;
  threadData[i].function = function;
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

static double mt_computeMean(ext_mt_manager_t restrict threadManager, const double* restrict x, size_t length)
{
  size_t numThreads, numValuesPerThread;
  ext_mt_getNumThreadsForJob(threadManager, length, MEAN_MIN_NUM_VALUES_PER_THREAD,
                             &numThreads, &numValuesPerThread);
  
  if (numThreads <= 1) return computeMean(x, length);
  
  MeanData threadData[numThreads];
  setupMeanData(threadData, numThreads, x, length, numValuesPerThread, computeMean);
  
  void* threadDataPtrs[numThreads];
  for (size_t i = 0; i < numThreads; ++i) threadDataPtrs[i] = (void*) &threadData[i];
  
  
  
  ext_mt_runTasks(threadManager, &computeMeanTask, threadDataPtrs, numThreads);
  
  
  
  return aggregateMeanResults(threadData, numThreads);
}

static double mt_computeUnrolledMean(ext_mt_manager_t restrict threadManager, const double* restrict x, size_t length)
{
  size_t numThreads, numValuesPerThread;
  ext_mt_getNumThreadsForJob(threadManager, length, UNROLLED_MEAN_MIN_NUM_VALUES_PER_THREAD,
                           &numThreads, &numValuesPerThread);
  
  if (numThreads <= 1) return computeUnrolledMean(x, length);
  
  MeanData threadData[numThreads];
  setupMeanData(threadData, numThreads, x, length, numValuesPerThread, computeUnrolledMean);
  
  void* threadDataPtrs[numThreads];
  for (size_t i = 0; i < numThreads; ++i) threadDataPtrs[i] = (void*) &threadData[i];
  
  
  
  ext_mt_runTasks(threadManager, &computeMeanTask, threadDataPtrs, numThreads);
  
  
  
  return aggregateMeanResults(threadData, numThreads);
}

static double mt_computeOnlineMean(ext_mt_manager_t restrict threadManager, const double* restrict x, size_t length)
{
  size_t numThreads, numValuesPerThread;
  ext_mt_getNumThreadsForJob(threadManager, length, ONLINE_MEAN_MIN_NUM_VALUES_PER_THREAD,
                           &numThreads, &numValuesPerThread);
  
  if (numThreads <= 1) return computeOnlineMean(x, length);
  
  MeanData threadData[numThreads];
  setupMeanData(threadData, numThreads, x, length, numValuesPerThread, computeOnlineMean);
  
  void* threadDataPtrs[numThreads];
  for (size_t i = 0; i < numThreads; ++i) threadDataPtrs[i] = (void*) &threadData[i];
  
  
  
  ext_mt_runTasks(threadManager, &computeMeanTask, threadDataPtrs, numThreads);
  
  
  
  return aggregateMeanResults(threadData, numThreads);
}

static double mt_computeOnlineUnrolledMean(ext_mt_manager_t restrict threadManager, const double* restrict x, size_t length)
{
  size_t numThreads, numValuesPerThread;
  ext_mt_getNumThreadsForJob(threadManager, length, ONLINE_UNROLLED_MEAN_MIN_NUM_VALUES_PER_THREAD,
                           &numThreads, &numValuesPerThread);
  
  if (numThreads <= 1) return computeOnlineUnrolledMean(x, length);
  
  MeanData threadData[numThreads];
  setupMeanData(threadData, numThreads, x, length, numValuesPerThread, computeOnlineUnrolledMean);
  
  void* threadDataPtrs[numThreads];
  for (size_t i = 0; i < numThreads; ++i) threadDataPtrs[i] = (void*) &threadData[i];
  
  
  
  ext_mt_runTasks(threadManager, &computeMeanTask, threadDataPtrs, numThreads);
  
  
  
  return aggregateMeanResults(threadData, numThreads);
}

typedef double (*indexedMeanFunction)(const double* restrict x, const size_t* restrict indices, size_t length);

typedef struct {
  const double* x;
  const size_t* indices;
  size_t length;
  double result;
  indexedMeanFunction function;
} IndexedMeanData;

static void computeIndexedMeanTask(void* v_data)
{
  IndexedMeanData* data = (IndexedMeanData*) v_data;
  data->result = data->function(data->x, data->indices, data->length);
}

static void setupIndexedMeanData(IndexedMeanData* restrict threadData, size_t numThreads, const double* restrict x,
                                 const size_t* restrict indices, size_t length, size_t numValuesPerThread, indexedMeanFunction function)
{
  size_t i = 0;
  for ( ; i < numThreads - 1; ++i) {
    threadData[i].x = x;
    threadData[i].indices = indices + i * numValuesPerThread;
    threadData[i].length = numValuesPerThread;
    threadData[i].function = function;
  }
  threadData[i].x = x;
  threadData[i].indices = indices + i * numValuesPerThread;
  threadData[i].length = length - i * numValuesPerThread;
  threadData[i].function = function;
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

static double mt_computeIndexedMean(ext_mt_manager_t restrict threadManager, const double* restrict x, const size_t* restrict indices, size_t length)
{
  size_t numThreads, numValuesPerThread;
  ext_mt_getNumThreadsForJob(threadManager, length, INDEXED_MEAN_MIN_NUM_VALUES_PER_THREAD,
                           &numThreads, &numValuesPerThread);
  
  if (numThreads <= 1) return computeIndexedMean(x, indices, length);
  
  IndexedMeanData threadData[numThreads];
  setupIndexedMeanData(threadData, numThreads, x, indices, length, numValuesPerThread, computeIndexedMean);
  
  void* threadDataPtrs[numThreads];
  for (size_t i = 0; i < numThreads; ++i) threadDataPtrs[i] = (void*) &threadData[i];
  
  
  
  ext_mt_runTasks(threadManager, &computeIndexedMeanTask, threadDataPtrs, numThreads);
  
  
  
  return aggregateIndexedMeanResults(threadData, numThreads);
}

static double mt_computeIndexedUnrolledMean(ext_mt_manager_t restrict threadManager, const double* restrict x, const size_t* restrict indices, size_t length)
{
  size_t numThreads, numValuesPerThread;
  ext_mt_getNumThreadsForJob(threadManager, length, INDEXED_UNROLLED_MEAN_MIN_NUM_VALUES_PER_THREAD,
                           &numThreads, &numValuesPerThread);
  
  if (numThreads <= 1) return computeIndexedUnrolledMean(x, indices, length);
  
  IndexedMeanData threadData[numThreads];
  setupIndexedMeanData(threadData, numThreads, x, indices, length, numValuesPerThread, computeIndexedUnrolledMean);
  
  void* threadDataPtrs[numThreads];
  for (size_t i = 0; i < numThreads; ++i) threadDataPtrs[i] = (void*) &threadData[i];
  
  
  
  ext_mt_runTasks(threadManager, &computeIndexedMeanTask, threadDataPtrs, numThreads);
  
  
  
  return aggregateIndexedMeanResults(threadData, numThreads);
}

static double mt_computeIndexedOnlineMean(ext_mt_manager_t restrict threadManager, const double* restrict x, const size_t* restrict indices, size_t length)
{
  size_t numThreads, numValuesPerThread;
  ext_mt_getNumThreadsForJob(threadManager, length, ONLINE_MEAN_MIN_NUM_VALUES_PER_THREAD,
                           &numThreads, &numValuesPerThread);
  
  if (numThreads <= 1) return computeIndexedOnlineMean(x, indices, length);
  
  IndexedMeanData threadData[numThreads];
  setupIndexedMeanData(threadData, numThreads, x, indices, length, numValuesPerThread, computeIndexedOnlineMean);
  
  void* threadDataPtrs[numThreads];
  for (size_t i = 0; i < numThreads; ++i) threadDataPtrs[i] = (void*) &threadData[i];
  
  
  
  ext_mt_runTasks(threadManager, &computeIndexedMeanTask, threadDataPtrs, numThreads);
  
  
  
  return aggregateIndexedMeanResults(threadData, numThreads);
}

static double mt_computeIndexedOnlineUnrolledMean(ext_mt_manager_t restrict threadManager, const double* restrict x, const size_t* restrict indices, size_t length)
{
  size_t numThreads, numValuesPerThread;
  ext_mt_getNumThreadsForJob(threadManager, length, ONLINE_UNROLLED_MEAN_MIN_NUM_VALUES_PER_THREAD,
                           &numThreads, &numValuesPerThread);
  
  if (numThreads <= 1) return computeIndexedOnlineUnrolledMean(x, indices, length);
  
  IndexedMeanData threadData[numThreads];
  setupIndexedMeanData(threadData, numThreads, x, indices, length, numValuesPerThread, computeIndexedOnlineUnrolledMean);
  
  void* threadDataPtrs[numThreads];
  for (size_t i = 0; i < numThreads; ++i) threadDataPtrs[i] = (void*) &threadData[i];
  
  
  
  ext_mt_runTasks(threadManager, &computeIndexedMeanTask, threadDataPtrs, numThreads);
  
  
  
  return aggregateIndexedMeanResults(threadData, numThreads);
}


typedef double (*varianceForKnownMeanFunction)(const double* x, size_t length, double mean);

typedef struct {
  const double* x;
  size_t length;
  double mean;
  double result;
  varianceForKnownMeanFunction function;
} VarianceForKnownMeanData;


static void computeVarianceForKnownMeanTask(void* v_data)
{
  VarianceForKnownMeanData* data = (VarianceForKnownMeanData*) v_data;
  data->result = data->function(data->x, data->length, data->mean);
}

static void setupVarianceForKnownMeanData(const double* restrict x, size_t length, VarianceForKnownMeanData* restrict threadData, size_t numThreads,
                                         size_t numValuesPerThread, varianceForKnownMeanFunction function, double mean)
{
  size_t i = 0;
  for ( ; i < numThreads - 1; ++i) {
    threadData[i].x = x + i * numValuesPerThread;
    threadData[i].length = numValuesPerThread;
    threadData[i].function = function;
    threadData[i].mean = mean;
  }
  threadData[i].x = x + i * numValuesPerThread;
  threadData[i].length = length - i * numValuesPerThread;
  threadData[i].function = function;
  threadData[i].mean = mean;
}

static double aggregateVarianceForKnownMeanData(const VarianceForKnownMeanData* threadData, size_t numThreads)
{
  double result = threadData[0].result;
  size_t setSize = threadData[0].length - 1;
  for (size_t i = 1; i < numThreads; ++i) {
    setSize += threadData[i].length;
    result += ((double) (threadData[i].length - 1) / (double) setSize) * (threadData[i].result - result);
  }
  return result;
}


static double mt_computeVarianceForKnownMean(ext_mt_manager_t restrict threadManager, const double* restrict x,
                                      size_t length, double mean)
{
  size_t numThreads, numValuesPerThread;
  ext_mt_getNumThreadsForJob(threadManager, length, VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD,
                           &numThreads, &numValuesPerThread);
  
  if (numThreads <= 1) return computeVarianceForKnownMean(x, length, mean);
  
  VarianceForKnownMeanData threadData[numThreads];
  setupVarianceForKnownMeanData(x, length, threadData, numThreads, numValuesPerThread, &computeVarianceForKnownMean, mean);
  
  void* threadDataPtrs[numThreads];
  for (size_t i = 0; i < numThreads; ++i) threadDataPtrs[i] = (void*) &threadData[i];
  
    
  ext_mt_runTasks(threadManager, &computeVarianceForKnownMeanTask, threadDataPtrs, numThreads);
  
  return aggregateVarianceForKnownMeanData(threadData, numThreads);
}

static double mt_computeUnrolledVarianceForKnownMean(ext_mt_manager_t restrict threadManager, const double* restrict x,
                                              size_t length, double mean)
{
  size_t numThreads, numValuesPerThread;
  ext_mt_getNumThreadsForJob(threadManager, length, UNROLLED_VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD,
                           &numThreads, &numValuesPerThread);
  
  if (numThreads <= 1) return computeUnrolledVarianceForKnownMean(x, length, mean);
  
  VarianceForKnownMeanData threadData[numThreads];
  setupVarianceForKnownMeanData(x, length, threadData, numThreads, numValuesPerThread, &computeUnrolledVarianceForKnownMean, mean);
  
  void* threadDataPtrs[numThreads];
  for (size_t i = 0; i < numThreads; ++i) threadDataPtrs[i] = (void*) &threadData[i];
  
  
  ext_mt_runTasks(threadManager, &computeVarianceForKnownMeanTask, threadDataPtrs, numThreads);
  
  return aggregateVarianceForKnownMeanData(threadData, numThreads);
}

static double mt_computeOnlineVarianceForKnownMean(ext_mt_manager_t restrict threadManager, const double* restrict x,
                                                 size_t length, double mean)
{
  size_t numThreads, numValuesPerThread;
  ext_mt_getNumThreadsForJob(threadManager, length, ONLINE_VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD,
                           &numThreads, &numValuesPerThread);
  
  if (numThreads <= 1) return computeOnlineVarianceForKnownMean(x, length, mean);
  
  VarianceForKnownMeanData threadData[numThreads];
  setupVarianceForKnownMeanData(x, length, threadData, numThreads, numValuesPerThread, &computeOnlineVarianceForKnownMean, mean);
  
  void* threadDataPtrs[numThreads];
  for (size_t i = 0; i < numThreads; ++i) threadDataPtrs[i] = (void*) &threadData[i];
  
  
  ext_mt_runTasks(threadManager, &computeVarianceForKnownMeanTask, threadDataPtrs, numThreads);
  
  return aggregateVarianceForKnownMeanData(threadData, numThreads);
}

static double mt_computeOnlineUnrolledVarianceForKnownMean(ext_mt_manager_t restrict threadManager, const double* restrict x,
                                                    size_t length, double mean)
{
  size_t numThreads, numValuesPerThread;
  ext_mt_getNumThreadsForJob(threadManager, length, ONLINE_UNROLLED_VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD,
                           &numThreads, &numValuesPerThread);
  
  if (numThreads <= 1) return computeOnlineUnrolledVarianceForKnownMean(x, length, mean);
  
  VarianceForKnownMeanData threadData[numThreads];
  setupVarianceForKnownMeanData(x, length, threadData, numThreads, numValuesPerThread, &computeOnlineUnrolledVarianceForKnownMean, mean);
  
  void* threadDataPtrs[numThreads];
  for (size_t i = 0; i < numThreads; ++i) threadDataPtrs[i] = (void*) &threadData[i];
  
  
  ext_mt_runTasks(threadManager, &computeVarianceForKnownMeanTask, threadDataPtrs, numThreads);
  
  return aggregateVarianceForKnownMeanData(threadData, numThreads);
}


typedef double (*indexedVarianceForKnownMeanFunction)(const double* restrict x, const size_t* restrict indices, size_t length, double mean);

typedef struct {
  const double* x;
  const size_t* indices;
  size_t length;
  double mean;
  double result;
  indexedVarianceForKnownMeanFunction function;
} IndexedVarianceForKnownMeanData;


static void computeIndexedVarianceForKnownMeanTask(void* v_data)
{
  IndexedVarianceForKnownMeanData* data = (IndexedVarianceForKnownMeanData*) v_data;
  data->result = data->function(data->x, data->indices, data->length, data->mean);
}

static void setupIndexedVarianceForKnownMeanData(IndexedVarianceForKnownMeanData* restrict threadData, size_t numThreads, const double* restrict x,
                                                 const size_t* restrict indices, size_t length, size_t numValuesPerThread, indexedVarianceForKnownMeanFunction function, double mean)
{
  size_t i = 0;
  for ( ; i < numThreads - 1; ++i) {
    threadData[i].x = x;
    threadData[i].indices = indices + i * numValuesPerThread;
    threadData[i].length = numValuesPerThread;
    threadData[i].function = function;
    threadData[i].mean = mean;
  }
  threadData[i].x = x;
  threadData[i].indices = indices + i * numValuesPerThread;
  threadData[i].length = length - i * numValuesPerThread;
  threadData[i].function = function;
  threadData[i].mean = mean;
}

static double aggregateIndexedVarianceForKnownMeanData(const IndexedVarianceForKnownMeanData* threadData, size_t numThreads)
{
  double result = threadData[0].result;
  size_t setSize = threadData[0].length - 1;
  for (size_t i = 1; i < numThreads; ++i) {
    setSize += threadData[i].length;
    result += ((double) (threadData[i].length - 1) / (double) setSize) * (threadData[i].result - result);
  }
  return result;
}


static double mt_computeIndexedVarianceForKnownMean(ext_mt_manager_t restrict threadManager, const double* restrict x,
                                             const size_t* restrict indices, size_t length, double mean)
{
  size_t numThreads, numValuesPerThread;
  ext_mt_getNumThreadsForJob(threadManager, length, INDEXED_VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD,
                           &numThreads, &numValuesPerThread);
  
  if (numThreads <= 1) return computeIndexedVarianceForKnownMean(x, indices, length, mean);
  
  IndexedVarianceForKnownMeanData threadData[numThreads];
  setupIndexedVarianceForKnownMeanData(threadData, numThreads, x, indices, length, numValuesPerThread, &computeIndexedVarianceForKnownMean, mean);
  
  void* threadDataPtrs[numThreads];
  for (size_t i = 0; i < numThreads; ++i) threadDataPtrs[i] = (void*) &threadData[i];
  
  
  ext_mt_runTasks(threadManager, &computeIndexedVarianceForKnownMeanTask, threadDataPtrs, numThreads);
  
  return aggregateIndexedVarianceForKnownMeanData(threadData, numThreads);
}

static double mt_computeIndexedUnrolledVarianceForKnownMean(ext_mt_manager_t restrict threadManager, const double* restrict x,
                                                     const size_t* restrict indices, size_t length, double mean)
{
  size_t numThreads, numValuesPerThread;
  ext_mt_getNumThreadsForJob(threadManager, length, INDEXED_UNROLLED_VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD,
                           &numThreads, &numValuesPerThread);
  
  if (numThreads <= 1) return computeIndexedUnrolledVarianceForKnownMean(x, indices, length, mean);
  
  IndexedVarianceForKnownMeanData threadData[numThreads];
  setupIndexedVarianceForKnownMeanData(threadData, numThreads, x, indices, length, numValuesPerThread, &computeIndexedUnrolledVarianceForKnownMean, mean);
  
  void* threadDataPtrs[numThreads];
  for (size_t i = 0; i < numThreads; ++i) threadDataPtrs[i] = (void*) &threadData[i];
  
  
  ext_mt_runTasks(threadManager, &computeIndexedVarianceForKnownMeanTask, threadDataPtrs, numThreads);
  
  return aggregateIndexedVarianceForKnownMeanData(threadData, numThreads);
}

static double mt_computeIndexedOnlineVarianceForKnownMean(ext_mt_manager_t restrict threadManager, const double* restrict x,
                                                   const size_t* restrict indices, size_t length, double mean)
{
  size_t numThreads, numValuesPerThread;
  ext_mt_getNumThreadsForJob(threadManager, length, INDEXED_ONLINE_VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD,
                           &numThreads, &numValuesPerThread);
  
  if (numThreads <= 1) return computeIndexedOnlineVarianceForKnownMean(x, indices, length, mean);
  
  IndexedVarianceForKnownMeanData threadData[numThreads];
  setupIndexedVarianceForKnownMeanData(threadData, numThreads, x, indices, length, numValuesPerThread, &computeIndexedOnlineVarianceForKnownMean, mean);
  
  void* threadDataPtrs[numThreads];
  for (size_t i = 0; i < numThreads; ++i) threadDataPtrs[i] = (void*) &threadData[i];
  
  
  ext_mt_runTasks(threadManager, &computeIndexedVarianceForKnownMeanTask, threadDataPtrs, numThreads);
  
  return aggregateIndexedVarianceForKnownMeanData(threadData, numThreads);
}

static double mt_computeIndexedOnlineUnrolledVarianceForKnownMean(ext_mt_manager_t restrict threadManager, const double* restrict x,
                                                           const size_t* restrict indices, size_t length, double mean)
{
  size_t numThreads, numValuesPerThread;
  ext_mt_getNumThreadsForJob(threadManager, length, INDEXED_ONLINE_UNROLLED_VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD,
                           &numThreads, &numValuesPerThread);
  
  if (numThreads <= 1) return computeIndexedOnlineUnrolledVarianceForKnownMean(x, indices, length, mean);
  
  IndexedVarianceForKnownMeanData threadData[numThreads];
  setupIndexedVarianceForKnownMeanData(threadData, numThreads, x, indices, length, numValuesPerThread, &computeIndexedOnlineUnrolledVarianceForKnownMean, mean);
  
  void* threadDataPtrs[numThreads];
  for (size_t i = 0; i < numThreads; ++i) threadDataPtrs[i] = (void*) &threadData[i];
  
  
  ext_mt_runTasks(threadManager, &computeIndexedVarianceForKnownMeanTask, threadDataPtrs, numThreads);
  
  return aggregateIndexedVarianceForKnownMeanData(threadData, numThreads);
}

typedef double (*varianceFunction)(const double* restrict x, size_t length, double* restrict meanPtr);

typedef struct {
  const double* x;
  size_t length;
  double mean;
  double var;
  varianceFunction function;
} VarianceData;

static void computeVarianceTask(void* v_data) {
  VarianceData* data = (VarianceData*) v_data;
  data->var = data->function(data->x, data->length, &data->mean);
}


static void setupVarianceData(VarianceData* restrict threadData, size_t numThreads, const double* restrict x, size_t length,
                              size_t numValuesPerThread, varianceFunction function)
{
  size_t i = 0;
  for ( ; i < numThreads - 1; ++i) {
    threadData[i].x = x + i * numValuesPerThread;
    threadData[i].length = numValuesPerThread;
    threadData[i].function = function;
  }
  threadData[i].x = x + i * numValuesPerThread;
  threadData[i].length = length - i * numValuesPerThread;
  threadData[i].function = function;
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

static double mt_computeVariance(ext_mt_manager_t restrict threadManager, const double* restrict x,
                          size_t length, double* restrict meanPtr)
{
  size_t numThreads, numValuesPerThread;
  ext_mt_getNumThreadsForJob(threadManager, length, VAR_MIN_NUM_VALUES_PER_THREAD,
                           &numThreads, &numValuesPerThread);
  
  if (numThreads <= 1) return computeVariance(x, length, meanPtr);
  
  VarianceData threadData[numThreads];
  setupVarianceData(threadData, numThreads, x, length, numValuesPerThread, &computeVariance);
  
  void* threadDataPtrs[numThreads];
  for (size_t i = 0; i < numThreads; ++i) threadDataPtrs[i] = (void*) &threadData[i];
  
  
  ext_mt_runTasks(threadManager, &computeVarianceTask, threadDataPtrs, numThreads);
  
  
  return aggregateVarianceData(threadData, numThreads, meanPtr);
}

static double mt_computeUnrolledVariance(ext_mt_manager_t restrict threadManager, const double* restrict x,
                                  size_t length, double* restrict meanPtr)
{
  size_t numThreads, numValuesPerThread;
  ext_mt_getNumThreadsForJob(threadManager, length, UNROLLED_VAR_MIN_NUM_VALUES_PER_THREAD,
                           &numThreads, &numValuesPerThread);
  
  if (numThreads <= 1) return computeUnrolledVariance(x, length, meanPtr);
  
  VarianceData threadData[numThreads];
  setupVarianceData(threadData, numThreads, x, length, numValuesPerThread, &computeUnrolledVariance);
  
  void* threadDataPtrs[numThreads];
  for (size_t i = 0; i < numThreads; ++i) threadDataPtrs[i] = (void*) &threadData[i];
  
  
  ext_mt_runTasks(threadManager, &computeVarianceTask, threadDataPtrs, numThreads);
  
  
  return aggregateVarianceData(threadData, numThreads, meanPtr);
}

static double mt_computeOnlineVariance(ext_mt_manager_t restrict threadManager, const double* restrict x,
                                size_t length, double* restrict meanPtr)
{
  size_t numThreads, numValuesPerThread;
  ext_mt_getNumThreadsForJob(threadManager, length, ONLINE_VAR_MIN_NUM_VALUES_PER_THREAD,
                           &numThreads, &numValuesPerThread);
  
  if (numThreads <= 1) return computeOnlineVariance(x, length, meanPtr);
  
  VarianceData threadData[numThreads];
  setupVarianceData(threadData, numThreads, x, length, numValuesPerThread, &computeOnlineVariance);
  
  void* threadDataPtrs[numThreads];
  for (size_t i = 0; i < numThreads; ++i) threadDataPtrs[i] = (void*) &threadData[i];
  
  
  ext_mt_runTasks(threadManager, &computeVarianceTask, threadDataPtrs, numThreads);
  
  
  return aggregateVarianceData(threadData, numThreads, meanPtr);
}

static double mt_computeOnlineUnrolledVariance(ext_mt_manager_t restrict threadManager, const double* restrict x,
                                        size_t length, double* restrict meanPtr)
{
  size_t numThreads, numValuesPerThread;
  ext_mt_getNumThreadsForJob(threadManager, length, ONLINE_UNROLLED_VAR_MIN_NUM_VALUES_PER_THREAD,
                           &numThreads, &numValuesPerThread);
  
  if (numThreads <= 1) return computeOnlineUnrolledVariance(x, length, meanPtr);
  
  VarianceData threadData[numThreads];
  setupVarianceData(threadData, numThreads, x, length, numValuesPerThread, &computeOnlineUnrolledVariance);
  
  void* threadDataPtrs[numThreads];
  for (size_t i = 0; i < numThreads; ++i) threadDataPtrs[i] = (void*) &threadData[i];
  
  
  ext_mt_runTasks(threadManager, &computeVarianceTask, threadDataPtrs, numThreads);
  
  
  return aggregateVarianceData(threadData, numThreads, meanPtr);
}


typedef double (*indexedVarianceFunction)(const double* restrict x, const size_t* restrict indices, size_t length, double* restrict meanPtr);

typedef struct {
  const double* x;
  const size_t* indices;
  size_t length;
  double mean;
  double var;
  indexedVarianceFunction function;
} IndexedVarianceData;

static void computeIndexedVarianceTask(void* v_data) {
  IndexedVarianceData* data = (IndexedVarianceData*) v_data;
  data->var = data->function(data->x, data->indices, data->length, &data->mean);
}


static void setupIndexedVarianceData(IndexedVarianceData* restrict threadData, size_t numThreads, const double* restrict x,
                                     const size_t* restrict indices, size_t length,
                                     size_t numValuesPerThread, indexedVarianceFunction function)
{
  size_t i = 0;
  for ( ; i < numThreads - 1; ++i) {
    threadData[i].x = x;
    threadData[i].indices = indices + i * numValuesPerThread;
    threadData[i].length = numValuesPerThread;
    threadData[i].function = function;
  }
  threadData[i].x = x;
  threadData[i].indices = indices + i * numValuesPerThread;
  threadData[i].length = length - i * numValuesPerThread;
  threadData[i].function = function;
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

static double mt_computeIndexedVariance(ext_mt_manager_t restrict threadManager, const double* restrict x,
                                 const size_t* restrict indices, size_t length, double* restrict meanPtr)
{
  size_t numThreads, numValuesPerThread;
  ext_mt_getNumThreadsForJob(threadManager, length, INDEXED_VAR_MIN_NUM_VALUES_PER_THREAD,
                           &numThreads, &numValuesPerThread);
  
  if (numThreads <= 1) return computeIndexedVariance(x, indices, length, meanPtr);
  
  IndexedVarianceData threadData[numThreads];
  setupIndexedVarianceData(threadData, numThreads, x, indices, length, numValuesPerThread, &computeIndexedVariance);
  
  void* threadDataPtrs[numThreads];
  for (size_t i = 0; i < numThreads; ++i) threadDataPtrs[i] = (void*) &threadData[i];
  
  
  ext_mt_runTasks(threadManager, &computeIndexedVarianceTask, threadDataPtrs, numThreads);
  
  
  return aggregateIndexedVarianceData(threadData, numThreads, meanPtr);
}

static double mt_computeIndexedUnrolledVariance(ext_mt_manager_t restrict threadManager, const double* restrict x,
                                         const size_t* restrict indices, size_t length, double* restrict meanPtr)
{
  size_t numThreads, numValuesPerThread;
  ext_mt_getNumThreadsForJob(threadManager, length, INDEXED_UNROLLED_VAR_MIN_NUM_VALUES_PER_THREAD,
                           &numThreads, &numValuesPerThread);
  
  if (numThreads <= 1) return computeIndexedUnrolledVariance(x,  indices, length, meanPtr);
  
  IndexedVarianceData threadData[numThreads];
  setupIndexedVarianceData(threadData, numThreads, x, indices, length, numValuesPerThread, &computeIndexedUnrolledVariance);
  
  void* threadDataPtrs[numThreads];
  for (size_t i = 0; i < numThreads; ++i) threadDataPtrs[i] = (void*) &threadData[i];
  
  
  ext_mt_runTasks(threadManager, &computeIndexedVarianceTask, threadDataPtrs, numThreads);
  
  
  return aggregateIndexedVarianceData(threadData, numThreads, meanPtr);
}

static double mt_computeIndexedOnlineVariance(ext_mt_manager_t restrict threadManager, const double* restrict x,
                                       const size_t* restrict indices, size_t length, double* restrict meanPtr)
{
  size_t numThreads, numValuesPerThread;
  ext_mt_getNumThreadsForJob(threadManager, length, INDEXED_ONLINE_VAR_MIN_NUM_VALUES_PER_THREAD,
                           &numThreads, &numValuesPerThread);
  
  if (numThreads <= 1) return computeIndexedOnlineVariance(x, indices, length, meanPtr);
  
  IndexedVarianceData threadData[numThreads];
  setupIndexedVarianceData(threadData, numThreads, x, indices, length, numValuesPerThread, &computeIndexedOnlineVariance);
  
  void* threadDataPtrs[numThreads];
  for (size_t i = 0; i < numThreads; ++i) threadDataPtrs[i] = (void*) &threadData[i];
  
  
  ext_mt_runTasks(threadManager, &computeIndexedVarianceTask, threadDataPtrs, numThreads);
  
  
  return aggregateIndexedVarianceData(threadData, numThreads, meanPtr);
}

static double mt_computeIndexedOnlineUnrolledVariance(ext_mt_manager_t restrict threadManager, const double* restrict x,
                                               const size_t* restrict indices, size_t length, double* restrict meanPtr)
{
  size_t numThreads, numValuesPerThread;
  ext_mt_getNumThreadsForJob(threadManager, length, INDEXED_ONLINE_UNROLLED_VAR_MIN_NUM_VALUES_PER_THREAD,
                           &numThreads, &numValuesPerThread);
  
  if (numThreads <= 1) return computeIndexedOnlineUnrolledVariance(x, indices, length, meanPtr);
  
  IndexedVarianceData threadData[numThreads];
  setupIndexedVarianceData(threadData, numThreads, x, indices, length, numValuesPerThread, &computeIndexedOnlineUnrolledVariance);
  
  void* threadDataPtrs[numThreads];
  for (size_t i = 0; i < numThreads; ++i) threadDataPtrs[i] = (void*) &threadData[i];
  
  
  ext_mt_runTasks(threadManager, &computeIndexedVarianceTask, threadDataPtrs, numThreads);
  
  
  return aggregateIndexedVarianceData(threadData, numThreads, meanPtr);
}
