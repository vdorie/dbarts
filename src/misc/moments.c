#include "config.h"

#include <math.h>
#include <stdint.h>

#include <misc/stats.h>

// Mean, variance and sum-of-squared-residual reductions. Each has a vanilla
// unrolled form and, where the accumulated round-off would matter, an "online"
// running-average form; the unrolled kernels themselves are selected at runtime
// by misc_stat_setSIMDInstructionSet. Every entry point here is serial.


// the lengths past which the *Fast entry points take the online form, held per
// kernel because the unrolling pays off differently in each
#define ONLINE_UNROLLED_MEAN_MIN_NUM_VALUES_PER_THREAD 25000

#define UNROLLED_VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD 75000

#define UNROLLED_WEIGHTED_MEAN_MIN_NUM_VALUES_PER_THREAD 125000

#define UNROLLED_WEIGHTED_VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD 125000

#define INDEXED_UNROLLED_MEAN_MIN_NUM_VALUES_PER_THREAD 100000

#define INDEXED_UNROLLED_VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD 75000

#define INDEXED_UNROLLED_WEIGHTED_MEAN_MIN_NUM_VALUES_PER_THREAD 35000

#define INDEXED_UNROLLED_WEIGHTED_VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD 35000

// if < value, calculate straightup; otherwise use online
// algorithm to reduce round-off err
#define ONLINE_CUTOFF 10000
#define INDEXED_ONLINE_CUTOFF 10000


// various implementations
// mean
static double (*computeUnrolledMean)(const double* x, size_t length) = 0;
static double (*computeIndexedUnrolledMean)(const double* restrict x, const misc_index_t* restrict indices, size_t length) = 0;
static double (*computeOnlineUnrolledMean)(const double* restrict x, size_t length) = 0;
static double (*computeIndexedOnlineUnrolledMean)(const double* restrict x, const misc_index_t* restrict indices, size_t length) = 0;

// var for known mean
static double (*computeUnrolledVarianceForKnownMean)(const double* x, size_t length, double mean) = 0;
static double (*computeIndexedUnrolledVarianceForKnownMean)(const double* restrict x, const misc_index_t* restrict indices, size_t length, double mean) = 0;
static double (*computeOnlineUnrolledVarianceForKnownMean)(const double* restrict x, size_t length, double mean) = 0;
static double (*computeIndexedOnlineUnrolledVarianceForKnownMean)(const double* restrict x, const misc_index_t* restrict indices, size_t length, double mean) = 0;

// variance and mean together
static double computeUnrolledVariance             (const double* restrict x, size_t length, double* restrict meanPtr);
static double computeIndexedUnrolledVariance      (const double* restrict x, const misc_index_t* restrict indices, size_t length, double* restrict meanPtr);
static double computeIndexedOnlineUnrolledVariance(const double* restrict x, const misc_index_t* restrict indices, size_t length, double* restrict meanPtr);

// weighted mean
static double (*computeUnrolledWeightedMean)(const double* restrict x, size_t length, const double* restrict w, double* restrict nPtr) = 0;
static double (*computeIndexedUnrolledWeightedMean)(const double* restrict x, const misc_index_t* restrict indices, size_t length, const double* restrict w, double* restrict nPtr) = 0;
static double (*computeOnlineUnrolledWeightedMean)(const double* restrict x, size_t length, const double* restrict w, double* restrict nPtr) = 0;
static double (*computeIndexedOnlineUnrolledWeightedMean)(const double* restrict x, const misc_index_t* restrict indices, size_t length, const double* restrict w, double* restrict nPtr) = 0;

// weighted variance for known mean
static double (*computeUnrolledWeightedVarianceForKnownMean)(const double* restrict x, size_t length, const double* restrict w, double mean) = 0;
static double (*computeIndexedUnrolledWeightedVarianceForKnownMean)(const double* restrict x, const misc_index_t* restrict indices, size_t length, const double* restrict w, double mean) = 0;
static double (*computeOnlineUnrolledWeightedVarianceForKnownMean)(const double* restrict x, size_t length, const double* restrict w, double mean) = 0;
static double (*computeIndexedOnlineUnrolledWeightedVarianceForKnownMean)(const double* restrict x, const misc_index_t* restrict indices, size_t length, const double* restrict w, double mean) = 0;

// interface functions that dispatch to the workers
double misc_computeMean(const double* x, size_t length)
{
  if (length > ONLINE_CUTOFF) return computeOnlineUnrolledMean(x, length);
  return computeUnrolledMean(x, length);
}

double misc_computeIndexedMean(const double* restrict x, const misc_index_t* restrict indices, size_t length)
{
  if (length > ONLINE_CUTOFF) return computeIndexedOnlineUnrolledMean(x, indices, length);
  return computeIndexedUnrolledMean(x, indices, length);
}

double misc_computeVarianceForKnownMean(const double* x, size_t length, double mean)
{
  if (length > ONLINE_CUTOFF) return computeOnlineUnrolledVarianceForKnownMean(x, length, mean);
  return computeUnrolledVarianceForKnownMean(x, length, mean);
}

double misc_computeIndexedVarianceForKnownMean(const double* restrict x, const misc_index_t* restrict indices, size_t length, double mean)
{
  if (length > ONLINE_CUTOFF) return computeIndexedOnlineUnrolledVarianceForKnownMean(x, indices, length, mean);
  return computeIndexedUnrolledVarianceForKnownMean(x, indices, length, mean);
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

double misc_computeIndexedVariance(const double* restrict x, const misc_index_t* restrict indices, size_t length, double* restrict meanPtr)
{
  if (length > ONLINE_CUTOFF) return computeIndexedOnlineUnrolledVariance(x, indices, length, meanPtr);
  return computeIndexedUnrolledVariance(x, indices, length, meanPtr);
}

double misc_computeWeightedMean(const double* restrict x, size_t length, const double* restrict w, double* restrict n)
{
  if (length > ONLINE_CUTOFF) return computeOnlineUnrolledWeightedMean(x, length, w, n);
  return computeUnrolledWeightedMean(x, length, w, n);
}

double misc_computeIndexedWeightedMean(const double* restrict x, const misc_index_t* restrict indices, size_t length, const double* restrict w, double* restrict n)
{
  if (length > ONLINE_CUTOFF) return computeIndexedOnlineUnrolledWeightedMean(x, indices, length, w, n);
  return computeIndexedUnrolledWeightedMean(x, indices, length, w, n);
}

double misc_computeWeightedVarianceForKnownMean(const double* restrict x, size_t length, const double* restrict w, double mean)
{
  if (length > ONLINE_CUTOFF) return computeOnlineUnrolledWeightedVarianceForKnownMean(x, length, w, mean);
  return computeUnrolledWeightedVarianceForKnownMean(x, length, w, mean);
}

double misc_computeIndexedWeightedVarianceForKnownMean(const double* restrict x, const misc_index_t* restrict indices, size_t length, const double* restrict w, double mean)
{
  if (length > ONLINE_CUTOFF) return computeIndexedOnlineUnrolledWeightedVarianceForKnownMean(x, indices, length, w, mean);
  return computeIndexedUnrolledWeightedVarianceForKnownMean(x, indices, length, w, mean);
}

#define minimum(_A_, _B_) ((_A_) < (_B_) ? (_A_) : (_B_))

// fast paths: the unrolled accumulators the entry points above reach through
// a length-keyed online/vanilla dispatch, minus that indirection
double misc_computeMeanFast(const double* x, size_t length)
{
  if (length >= minimum(ONLINE_UNROLLED_MEAN_MIN_NUM_VALUES_PER_THREAD, ONLINE_CUTOFF))
    return computeOnlineUnrolledMean(x, length);
  return computeUnrolledMean(x, length);
}

double misc_computeIndexedMeanFast(const double* restrict x, const misc_index_t* restrict indices, size_t length)
{
  if (length >= minimum(INDEXED_UNROLLED_MEAN_MIN_NUM_VALUES_PER_THREAD, INDEXED_ONLINE_CUTOFF))
    return computeIndexedOnlineUnrolledMean(x, indices, length);
  return computeIndexedUnrolledMean(x, indices, length);
}

double misc_computeWeightedMeanFast(const double* restrict x, size_t length, const double* restrict w, double* restrict n)
{
  if (length >= minimum(UNROLLED_WEIGHTED_MEAN_MIN_NUM_VALUES_PER_THREAD, ONLINE_CUTOFF))
    return computeOnlineUnrolledWeightedMean(x, length, w, n);
  return computeUnrolledWeightedMean(x, length, w, n);
}

double misc_computeIndexedWeightedMeanFast(const double* restrict x, const misc_index_t* restrict indices, size_t length, const double* restrict w, double* restrict n)
{
  if (length >= minimum(INDEXED_UNROLLED_WEIGHTED_MEAN_MIN_NUM_VALUES_PER_THREAD, ONLINE_CUTOFF))
    return computeIndexedOnlineUnrolledWeightedMean(x, indices, length, w, n);
  return computeIndexedUnrolledWeightedMean(x, indices, length, w, n);
}

double misc_computeVarianceForKnownMeanFast(const double* x, size_t length, double mean)
{
  if (length >= minimum(UNROLLED_VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD, ONLINE_CUTOFF))
    return computeOnlineUnrolledVarianceForKnownMean(x, length, mean);
  return computeUnrolledVarianceForKnownMean(x, length, mean);
}

double misc_computeIndexedVarianceForKnownMeanFast(const double* restrict x, const misc_index_t* restrict indices, size_t length, double mean)
{
  if (length >= minimum(INDEXED_UNROLLED_VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD, INDEXED_ONLINE_CUTOFF))
    return computeIndexedOnlineUnrolledVarianceForKnownMean(x, indices, length, mean);
  return computeIndexedUnrolledVarianceForKnownMean(x, indices, length, mean);
}

double misc_computeWeightedVarianceForKnownMeanFast(const double* restrict x, size_t length, const double* restrict w, double mean)
{
  if (length >= minimum(UNROLLED_WEIGHTED_VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD, ONLINE_CUTOFF))
    return computeOnlineUnrolledWeightedVarianceForKnownMean(x, length, w, mean);
  return computeUnrolledWeightedVarianceForKnownMean(x, length, w, mean);
}

double misc_computeIndexedWeightedVarianceForKnownMeanFast(const double* restrict x, const misc_index_t* restrict indices, size_t length, const double* restrict w, double mean)
{
  if (length >= minimum(INDEXED_UNROLLED_WEIGHTED_VAR_FOR_MEAN_MIN_NUM_VALUES_PER_THREAD, ONLINE_CUTOFF))
    return computeIndexedOnlineUnrolledWeightedVarianceForKnownMean(x, indices, length, w, mean);
  return computeIndexedUnrolledWeightedVarianceForKnownMean(x, indices, length, w, mean);
}

// Fused (sum w, sum wx) suffstat, unrolled like the mean kernels above; a
// SIMD specialization would slot in behind the same signature if profiling
// asked for one. The raw sums are order-insensitive, unlike the
// mean-then-centered-variance pair these replace.
void misc_computeSufficientStatisticsFast(const double* x, size_t length, double* restrict sumW, double* restrict sumWX)
{
  *sumW = (double) length;

  size_t i = 0;
  size_t lengthMod5 = length % 5;

  double swx = 0.0;
  if (lengthMod5 != 0) {
    for ( ; i < lengthMod5; ++i) { swx += x[i]; }
    if (length < 5) { *sumWX = swx; return; }
  }

  for ( ; i < length; i += 5) {
    swx += x[i] + x[i + 1] + x[i + 2] + x[i + 3] + x[i + 4];
  }

  *sumWX = swx;
}

void misc_computeIndexedSufficientStatisticsFast(const double* restrict x, const misc_index_t* restrict indices, size_t length, double* restrict sumW, double* restrict sumWX)
{
  *sumW = (double) length;

  size_t i = 0;
  size_t lengthMod5 = length % 5;

  double swx = 0.0;
  if (lengthMod5 != 0) {
    for ( ; i < lengthMod5; ++i) { swx += x[indices[i]]; }
    if (length < 5) { *sumWX = swx; return; }
  }

  for ( ; i < length; i += 5) {
    double v0 = x[indices[i]], v1 = x[indices[i + 1]], v2 = x[indices[i + 2]],
           v3 = x[indices[i + 3]], v4 = x[indices[i + 4]];
    swx += v0 + v1 + v2 + v3 + v4;
  }

  *sumWX = swx;
}

void misc_computeWeightedSufficientStatisticsFast(const double* restrict x, size_t length, const double* restrict w, double* restrict sumW, double* restrict sumWX)
{
  size_t i = 0;
  size_t lengthMod5 = length % 5;

  double sw = 0.0, swx = 0.0;
  if (lengthMod5 != 0) {
    for ( ; i < lengthMod5; ++i) {
      double wi = w[i], v = x[i];
      sw += wi; swx += wi * v;
    }
    if (length < 5) { *sumW = sw; *sumWX = swx; return; }
  }

  for ( ; i < length; i += 5) {
    double wv0 = w[i] * x[i], wv1 = w[i + 1] * x[i + 1],
           wv2 = w[i + 2] * x[i + 2], wv3 = w[i + 3] * x[i + 3],
           wv4 = w[i + 4] * x[i + 4];
    sw += w[i] + w[i + 1] + w[i + 2] + w[i + 3] + w[i + 4];
    swx += wv0 + wv1 + wv2 + wv3 + wv4;
  }

  *sumW = sw;
  *sumWX = swx;
}

void misc_computeIndexedWeightedSufficientStatisticsFast(const double* restrict x, const misc_index_t* restrict indices, size_t length, const double* restrict w, double* restrict sumW, double* restrict sumWX)
{
  size_t i = 0;
  size_t lengthMod5 = length % 5;

  double sw = 0.0, swx = 0.0;
  if (lengthMod5 != 0) {
    for ( ; i < lengthMod5; ++i) {
      size_t j = indices[i];
      double wi = w[j], v = x[j];
      sw += wi; swx += wi * v;
    }
    if (length < 5) { *sumW = sw; *sumWX = swx; return; }
  }

  for ( ; i < length; i += 5) {
    size_t j0 = indices[i], j1 = indices[i + 1], j2 = indices[i + 2],
           j3 = indices[i + 3], j4 = indices[i + 4];
    double wv0 = w[j0] * x[j0], wv1 = w[j1] * x[j1], wv2 = w[j2] * x[j2],
           wv3 = w[j3] * x[j3], wv4 = w[j4] * x[j4];
    sw += w[j0] + w[j1] + w[j2] + w[j3] + w[j4];
    swx += wv0 + wv1 + wv2 + wv3 + wv4;
  }

  *sumW = sw;
  *sumWX = swx;
}

// fp32-residual variants of the four suffstat kernels above: the running
// residual x is stored fp32, so these load float and PROMOTE each element
// to double before
// accumulating - the reduction stays fp64 and, because float->double promotion
// is exact and the summation order mirrors the double kernels byte-for-byte,
// the sums equal round-to-fp32-then-sum-in-fp64 exactly. Only the opt-in
// gaussian constant-leaf path (tree.hpp computeLeafStats) selects them.
void misc_computeFloatSufficientStatisticsFast(const float* x, size_t length, double* restrict sumW, double* restrict sumWX)
{
  *sumW = (double) length;

  size_t i = 0;
  size_t lengthMod5 = length % 5;

  double swx = 0.0;
  if (lengthMod5 != 0) {
    for ( ; i < lengthMod5; ++i) { swx += x[i]; }
    if (length < 5) { *sumWX = swx; return; }
  }

  for ( ; i < length; i += 5) {
    swx += x[i] + x[i + 1] + x[i + 2] + x[i + 3] + x[i + 4];
  }

  *sumWX = swx;
}

void misc_computeIndexedFloatSufficientStatisticsFast(const float* restrict x, const misc_index_t* restrict indices, size_t length, double* restrict sumW, double* restrict sumWX)
{
  *sumW = (double) length;

  size_t i = 0;
  size_t lengthMod5 = length % 5;

  double swx = 0.0;
  if (lengthMod5 != 0) {
    for ( ; i < lengthMod5; ++i) { swx += x[indices[i]]; }
    if (length < 5) { *sumWX = swx; return; }
  }

  for ( ; i < length; i += 5) {
    double v0 = x[indices[i]], v1 = x[indices[i + 1]], v2 = x[indices[i + 2]],
           v3 = x[indices[i + 3]], v4 = x[indices[i + 4]];
    swx += v0 + v1 + v2 + v3 + v4;
  }

  *sumWX = swx;
}

void misc_computeWeightedFloatSufficientStatisticsFast(const float* restrict x, size_t length, const double* restrict w, double* restrict sumW, double* restrict sumWX)
{
  size_t i = 0;
  size_t lengthMod5 = length % 5;

  double sw = 0.0, swx = 0.0;
  if (lengthMod5 != 0) {
    for ( ; i < lengthMod5; ++i) {
      double wi = w[i], v = x[i];
      sw += wi; swx += wi * v;
    }
    if (length < 5) { *sumW = sw; *sumWX = swx; return; }
  }

  for ( ; i < length; i += 5) {
    double wv0 = w[i] * x[i], wv1 = w[i + 1] * x[i + 1],
           wv2 = w[i + 2] * x[i + 2], wv3 = w[i + 3] * x[i + 3],
           wv4 = w[i + 4] * x[i + 4];
    sw += w[i] + w[i + 1] + w[i + 2] + w[i + 3] + w[i + 4];
    swx += wv0 + wv1 + wv2 + wv3 + wv4;
  }

  *sumW = sw;
  *sumWX = swx;
}

void misc_computeIndexedWeightedFloatSufficientStatisticsFast(const float* restrict x, const misc_index_t* restrict indices, size_t length, const double* restrict w, double* restrict sumW, double* restrict sumWX)
{
  size_t i = 0;
  size_t lengthMod5 = length % 5;

  double sw = 0.0, swx = 0.0;
  if (lengthMod5 != 0) {
    for ( ; i < lengthMod5; ++i) {
      size_t j = indices[i];
      double wi = w[j], v = x[j];
      sw += wi; swx += wi * v;
    }
    if (length < 5) { *sumW = sw; *sumWX = swx; return; }
  }

  for ( ; i < length; i += 5) {
    size_t j0 = indices[i], j1 = indices[i + 1], j2 = indices[i + 2],
           j3 = indices[i + 3], j4 = indices[i + 4];
    double wv0 = w[j0] * x[j0], wv1 = w[j1] * x[j1], wv2 = w[j2] * x[j2],
           wv3 = w[j3] * x[j3], wv4 = w[j4] * x[j4];
    sw += w[j0] + w[j1] + w[j2] + w[j3] + w[j4];
    swx += wv0 + wv1 + wv2 + wv3 + wv4;
  }

  *sumW = sw;
  *sumWX = swx;
}

// work-horse implementations; indexed versions follow regular b/c implementations are highly similar

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

static double computeIndexedUnrolledMean_c(const double* restrict x, const misc_index_t* restrict indices, size_t length)
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

static double computeIndexedOnlineUnrolledMean_c(const double* restrict x, const misc_index_t* restrict indices, size_t length)
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

static double computeIndexedUnrolledVarianceForKnownMean_c(const double* restrict x, const misc_index_t* restrict indices, size_t length, double mean)
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

static double computeIndexedOnlineUnrolledVarianceForKnownMean_c(const double* restrict x, const misc_index_t* restrict indices, size_t length, double mean)
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

static double computeIndexedUnrolledVariance(const double* restrict x, const misc_index_t* restrict indices, size_t length, double* restrict meanPtr)
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

static double computeIndexedOnlineUnrolledVariance(const double* restrict x, const misc_index_t* restrict indices, size_t length, double* restrict meanPtr)
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
    const double n1 = (double) i;
    const double n2 = 5.0;
    const double n  = (double) (i + 5);
    
    double meanNext5 = (x[indices[i]] + x[indices[i + 1]] + x[indices[i + 2]] + x[indices[i + 3]] + x[indices[i + 4]]) / n2;
    double varNext5  = ((x[indices[i]] - meanNext5) * (x[indices[i]] - meanNext5) +
                        (x[indices[i + 1]] - meanNext5) * (x[indices[i + 1]] - meanNext5) +
                        (x[indices[i + 2]] - meanNext5) * (x[indices[i + 2]] - meanNext5) +
                        (x[indices[i + 3]] - meanNext5) * (x[indices[i + 3]] - meanNext5) +
                        (x[indices[i + 4]] - meanNext5) * (x[indices[i + 4]] - meanNext5)) / n2;
    
    var  += (n2 / n) * (varNext5 - var) + ((meanNext5 - mean) * (n1 / n)) * ((meanNext5 - mean) * (n2 / n));
    mean += (n2 / n) * (meanNext5 - mean);
  }
  
  if (meanPtr != NULL) *meanPtr = mean;
  return nScale * var;
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

static double computeIndexedUnrolledWeightedMean_c(const double* restrict x, const misc_index_t* restrict indices, size_t length, const double* restrict w, double* restrict nPtr)
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

static double computeIndexedOnlineUnrolledWeightedMean_c(const double* restrict x, const misc_index_t* restrict indices, size_t length, const double* restrict w, double* restrict nPtr)
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

static double computeIndexedUnrolledWeightedVarianceForKnownMean_c(const double* restrict x, const misc_index_t* restrict indices, size_t length, const double* restrict w, double mean)
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

static double computeIndexedOnlineUnrolledWeightedVarianceForKnownMean_c(const double* restrict x, const misc_index_t* restrict indices, size_t length, const double* restrict w, double mean)
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

#ifdef COMPILER_SUPPORTS_SSE2
extern double misc_computeUnrolledMean_sse2(const double* x, size_t length);
extern double misc_computeIndexedUnrolledMean_sse2(const double* restrict x, const misc_index_t* restrict indices, size_t length);
extern double misc_computeOnlineUnrolledMean_sse2(const double* x, size_t length);
extern double misc_computeIndexedOnlineUnrolledMean_sse2(const double* restrict x, const misc_index_t* restrict indices, size_t length);
extern double misc_computeUnrolledVarianceForKnownMean_sse2(const double* x, size_t length, double mean);
extern double misc_computeIndexedUnrolledVarianceForKnownMean_sse2(const double* restrict x, const misc_index_t* restrict indices, size_t length, double mean);
extern double misc_computeOnlineUnrolledVarianceForKnownMean_sse2(const double* x, size_t length, double mean);
extern double misc_computeIndexedOnlineUnrolledVarianceForKnownMean_sse2(const double* restrict x, const misc_index_t* restrict indices, size_t length, double mean);
extern double misc_computeUnrolledWeightedMean_sse2(const double* restrict x, size_t length, const double* restrict w, double* restrict nPtr);
extern double misc_computeIndexedUnrolledWeightedMean_sse2(const double* restrict x, const misc_index_t* restrict indices, size_t length, const double* restrict w, double* restrict nPtr);
extern double misc_computeOnlineUnrolledWeightedMean_sse2(const double* restrict x, size_t length, const double* restrict w, double* restrict nPtr);
extern double misc_computeIndexedOnlineUnrolledWeightedMean_sse2(const double* restrict x, const misc_index_t* restrict indices, size_t length, const double* restrict w, double* restrict nPtr);
extern double misc_computeUnrolledWeightedVarianceForKnownMean_sse2(const double* restrict x, size_t length, const double* restrict w, double mean);
extern double misc_computeIndexedUnrolledWeightedVarianceForKnownMean_sse2(const double* restrict x, const misc_index_t* restrict indices, size_t length, const double* restrict w, double mean);
extern double misc_computeOnlineUnrolledWeightedVarianceForKnownMean_sse2(const double* restrict x, size_t length, const double* restrict w, double mean);
extern double misc_computeIndexedOnlineUnrolledWeightedVarianceForKnownMean_sse2(const double* restrict x, const misc_index_t* restrict indices, size_t length, const double* restrict w, double mean);
#endif

void misc_stat_setSIMDInstructionSet(misc_simd_instructionSet i)
{
  (void) i; // unused when no x86 SIMD moments kernels are compiled in
#ifdef COMPILER_SUPPORTS_SSE2
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

double misc_computeSumOfSquaredResiduals(const double* restrict x, size_t length, const double* restrict x_hat)
{
  if (length == 0) return 0.0;
  
  double result = 0.0;
  size_t lengthMod5 = length % 5;
  
  size_t i = 0;
  for ( ; i < lengthMod5; ++i) result += (x[i] - x_hat[i]) * (x[i] - x_hat[i]);
  
  for ( ; i < length; i += 5) {
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
  for ( ; i < lengthMod5; ++i) result += w[i] * (x[i] - x_hat[i]) * (x[i] - x_hat[i]);
  
  for ( ; i < length; i += 5) {
    result += w[i] * (x[i] - x_hat[i]) * (x[i] - x_hat[i]) + 
              w[i + 1] * (x[i + 1] - x_hat[i + 1]) * (x[i + 1] - x_hat[i + 1]) + 
              w[i + 2] * (x[i + 2] - x_hat[i + 2]) * (x[i + 2] - x_hat[i + 2]) +
              w[i + 3] * (x[i + 3] - x_hat[i + 3]) * (x[i + 3] - x_hat[i + 3]) +
              w[i + 4] * (x[i + 4] - x_hat[i + 4]) * (x[i + 4] - x_hat[i + 4]);
  }
  
  return result;
}
