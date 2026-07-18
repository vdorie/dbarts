#ifndef EXTERNAL_RANDOM_H
#define EXTERNAL_RANDOM_H

// create an RNG and then seed it, or create it using a state constructed from
// somewhere else

// this is largely a thread-safe duplication of R's random number generator code

#include <misc/stddef.h>
#include <stdbool.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct ext_rng ext_rng;

typedef enum {
  EXT_RNG_ALGORITHM_MERSENNE_TWISTER = 0,
  EXT_RNG_ALGORITHM_USER_UNIFORM,
  EXT_RNG_ALGORITHM_INVALID // must be last
} ext_rng_algorithm_t;

typedef enum {
  EXT_RNG_STANDARD_NORMAL_USER_NORM = 0,
  EXT_RNG_STANDARD_NORMAL_INVERSION,
  EXT_RNG_STANDARD_NORMAL_INVALID // must be last
} ext_rng_standardNormal_t;

// state can be null; size is determined by algorithm; see below for a few state
// structs
ext_rng* ext_rng_create(ext_rng_algorithm_t algorithm, const void* state);
void ext_rng_destroy(ext_rng* generator);

void ext_rng_setState(ext_rng* generator, const void* state);

// createDefault seeds the result so it is ready to use
// useNative attempts to use the rng in the embedded environment; generally not
// thread safe
ext_rng* ext_rng_createDefault(bool useNative);

int ext_rng_createAndSeed(
  ext_rng** result,
  ext_rng_algorithm_t algorithm,
  ext_rng_standardNormal_t standardNormalAlgorithm
);

// returns what will be used in createDefault, unless useNative is specified;
ext_rng_algorithm_t ext_rng_getDefaultAlgorithmType(void);
ext_rng_standardNormal_t ext_rng_getDefaultStandardNormalType(void);
const char* ext_rng_getAlgorithmName(ext_rng_algorithm_t algorithm);

bool ext_rng_seedsAreEqual(const ext_rng* rng1, const ext_rng* rng2);
unsigned int ext_rng_getState0(const ext_rng* rng);

// state can be null; for USER_NORM, it should be a userFunction outlined below
int ext_rng_setStandardNormalAlgorithm(
  ext_rng* generator,
  ext_rng_standardNormal_t standardNormalAlgorithm,
  const void* state
);
int ext_rng_setSeed(ext_rng* generator, uint_least32_t seed);
int ext_rng_setSeedFromClock(ext_rng* generator);

// not the same as "state" above, since it also includes the status of the
// standard normal algorithm and gamma simulation;
// guarantees that the result is aligned to sizeof(int), however the length is
// in characters user is responsible for seralizing user functions in algorithm
// or standardNormalAlgorithm
misc_size_t ext_rng_getSerializedStateLength(const ext_rng* generator);
void ext_rng_writeSerializedState(const ext_rng* generator, void* state);
void ext_rng_readSerializedState(ext_rng* generator, const void* state);

bool ext_rng_usesNativeRNG(const ext_rng* generator);

double ext_rng_simulateContinuousUniform(ext_rng* generator); // randomBase.c
double ext_rng_simulateStandardNormal(ext_rng* generator);    // randomNorm.c

// standard normal truncated below at lowerBound, using Robert (1995)
double ext_rng_simulateLowerTruncatedStandardNormal(
  ext_rng* generator,
  double lowerBound
);
// use the previous to generate truncated normals with sd 1 and nonzero mean
double ext_rng_simulateLowerTruncatedNormalScale1(
  ext_rng* generator,
  double mean,
  double bound
);
double ext_rng_simulateUpperTruncatedNormalScale1(
  ext_rng* generator,
  double mean,
  double bound
);
double ext_rng_simulateLowerTruncatedNormal(
  ext_rng* generator,
  double mean,
  double sd,
  double bound
);
double ext_rng_simulateUpperTruncatedNormal(
  ext_rng* generator,
  double mean,
  double sd,
  double bound
);
// standard normal (sd 1) with mean, truncated to the interval (lower, upper]:
// inverse-CDF in the bulk, Robert (1995) rejection when the tail probability
// gap underflows
double ext_rng_simulateTruncatedNormalScale1(
  ext_rng* generator,
  double mean,
  double lower,
  double upper
);

// subsequent in random.c
double ext_rng_simulateExponential(ext_rng* generator, double scale);
double ext_rng_simulateGamma(ext_rng* generator, double shape, double scale);

// Generalized inverse Gaussian GIG(p, a, b): density proportional to
// x^(p - 1) exp(-(a x + b / x) / 2), x > 0. b == 0 is the Gamma(p, rate a/2)
// limit (needs p > 0), a == 0 the inverse-gamma limit (needs p < 0); general
// a, b > 0 uses one ratio-of-uniforms region for all p (see random.c).
double ext_rng_simulateGeneralizedInverseGaussian(
  ext_rng* generator,
  double p,
  double a,
  double b
);

// Polya-Gamma PG(1, psi), using the exact alternating-series method of
// Devroye (2009) as adapted by Polson, Scott, and Windle (2013)
double ext_rng_simulatePolyaGamma(ext_rng* generator, double psi);

#define ext_rng_simulateChiSquared(_GENERATOR_, _DF_)                          \
  ext_rng_simulateGamma(_GENERATOR_, (_DF_) / 2.0, 2.0)
#define ext_rng_simulateBernoulli(_GENERATOR_, _P_)                            \
  (ext_rng_simulateContinuousUniform(_GENERATOR_) < (_P_) ? 1u : 0u)
#define ext_rng_simulateNormal(_GENERATOR_, _MU_, _SIGMA_)                     \
  (ext_rng_simulateStandardNormal(_GENERATOR_) * (_SIGMA_) + (_MU_))

#define EXT_DISCRETE_DRAW_FAILURE ((misc_size_t) - 1)
misc_size_t ext_rng_drawFromDiscreteDistribution(
  ext_rng* generator,
  const double* probabilities,
  misc_size_t length
);

// random in [min, min + 1, ..., max - 1, max)
int64_t ext_rng_simulateIntegerUniformInRange(
  ext_rng* generator,
  int64_t min_inclusive,
  int64_t max_exclusive
);
uint64_t ext_rng_simulateUnsignedIntegerUniformInRange(
  ext_rng* generator,
  uint64_t min_inclusive,
  uint64_t max_exclusive
);

void ext_rng_drawPermutation(
  ext_rng* generator,
  misc_size_t* x,
  misc_size_t length
);

#define EXT_RNG_MERSENNE_TWISTER_NUM_RANDOM 624

typedef struct {
  int_least32_t info;
  uint_least32_t state[EXT_RNG_MERSENNE_TWISTER_NUM_RANDOM];
} ext_rng_mersenneTwisterState;

// used for EXT_RNG_ALGORITHM_USER_UNIFORM and EXT_RNG_STANDARD_NORMAL_USER_NORM
typedef struct {
  union {
    double (*stateless)(
      void
    ); // used if state is NULL to avoid a pointless function call
    double (*stateful)(void*);
  } f;

  void* state;
} ext_rng_userFunction;

#ifdef __cplusplus
}
#endif

#endif // EXTERNAL_RANDOM_H
