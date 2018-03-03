#ifndef EXTERNAL_RANDOM_H
#define EXTERNAL_RANDOM_H

// create an RNG and then seed it, or create it using a state constructed from somewhere else

#include <stdbool.h>
#include "stddef.h"
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct ext_rng ext_rng;

typedef enum {
    EXT_RNG_ALGORITHM_WICHMANN_HILL = 0,
    EXT_RNG_ALGORITHM_MARSAGLIA_MULTICARRY,
    EXT_RNG_ALGORITHM_SUPER_DUPER,
    EXT_RNG_ALGORITHM_MERSENNE_TWISTER,
    EXT_RNG_ALGORITHM_KNUTH_TAOCP,
    EXT_RNG_ALGORITHM_USER_UNIFORM,
    EXT_RNG_ALGORITHM_KNUTH_TAOCP2,
    EXT_RNG_ALGORITHM_LECUYER_CMRG,
    EXT_RNG_ALGORITHM_INVALID // must be last
} ext_rng_algorithm_t;

typedef enum {
    EXT_RNG_STANDARD_NORMAL_BUGGY_KINDERMAN_RAMAGE = 0,
    EXT_RNG_STANDARD_NORMAL_AHRENS_DIETER,
    EXT_RNG_STANDARD_NORMAL_BOX_MULLER,
    EXT_RNG_STANDARD_NORMAL_USER_NORM,
    EXT_RNG_STANDARD_NORMAL_INVERSION,
    EXT_RNG_STANDARD_NORMAL_KINDERMAN_RAMAGE,
    EXT_RNG_STANDARD_NORMAL_INVALID // must be last
} ext_rng_standardNormal_t;

// state can be null; size is determined by algorithm; see below for a few state structs
ext_rng* ext_rng_create(ext_rng_algorithm_t algorithm, const void* state);
void ext_rng_destroy(ext_rng* generator);

void ext_rng_setState(ext_rng* generator, const void* state);

// createDefault seeds the result so it is ready to use
// useNative attempts to use the rng in the embedded environment; generally not thread safe
ext_rng* ext_rng_createDefault(bool useNative);

int ext_rng_createAndSeed(ext_rng** result, ext_rng_algorithm_t algorithm, ext_rng_standardNormal_t standardNormalAlgorithm);

// returns what will be used in createDefault, unless useNative is specified;
ext_rng_algorithm_t ext_rng_getDefaultAlgorithmType();
ext_rng_standardNormal_t ext_rng_getDefaultStandardNormalType();

// state can be null; for BOX_MULLER, it should point to a double that is the next number, or 0.0 if that isn't set yet
// for USER_NORM, it should be a userFunction outlined below
int ext_rng_setStandardNormalAlgorithm(ext_rng* generator, ext_rng_standardNormal_t standardNormalAlgorithm, const void* state);
int ext_rng_setSeed(ext_rng* generator, uint_least32_t seed);
int ext_rng_setSeedFromClock(ext_rng* generator);

// not the same as "state" above, since it also includes the status of the
// standard normal algorithm and gamma simulation;
// guarantees that the result is aligned to sizeof(int), however the length is in characters
// user is responsible for seralizing user functions in algorithm or standardNormalAlgorithm
ext_size_t ext_rng_getSerializedStateLength(const ext_rng* generator);
void ext_rng_writeSerializedState(const ext_rng* generator, void* state);
void ext_rng_readSerializedState(ext_rng* generator, const void* state);

double ext_rng_simulateContinuousUniform(ext_rng* generator); // randomBase.c
double ext_rng_simulateStandardNormal(ext_rng* generator);    // randomNorm.c

// standard normal truncated below at lowerBound, using Robert (1995)
double ext_rng_simulateLowerTruncatedStandardNormal(ext_rng* generator, double lowerBound); 
// use the previous to generate truncated normals with sd 1 and nonzero mean
double ext_rng_simulateLowerTruncatedNormalScale1(ext_rng* generator, double mean, double bound);
double ext_rng_simulateUpperTruncatedNormalScale1(ext_rng* generator, double mean, double bound);

// subsequent in random.c
double ext_rng_simulateExponential(ext_rng* generator, double scale);
double ext_rng_simulateGamma(ext_rng* generator, double shape, double scale);

#define ext_rng_simulateChiSquared(_GENERATOR_, _DF_) ext_rng_simulateGamma(_GENERATOR_, (_DF_) / 2.0, 2.0)
#define ext_rng_simulateBernoulli(_GENERATOR_, _P_) (ext_rng_simulateContinuousUniform(_GENERATOR_) < (_P_) ? 1u : 0u)
#define ext_rng_simulateNormal(_GENERATOR_, _MU_, _SIGMA_) (ext_rng_simulateStandardNormal(_GENERATOR_) * (_SIGMA_) + (_MU_))


#define EXT_DISCRETE_DRAW_FAILURE ((ext_size_t) -1)
ext_size_t ext_rng_drawFromDiscreteDistribution(ext_rng* generator, const double* probabilities, ext_size_t length);
  
// random in [min, min + 1, ..., max - 1, max)
int64_t ext_rng_simulateIntegerUniformInRange(ext_rng* generator, int64_t min_inclusive, int64_t max_exclusive);
uint64_t ext_rng_simulateUnsignedIntegerUniformInRange(ext_rng* generator, uint64_t min_inclusive, uint64_t max_exclusive);



#define EXT_RNG_KNUTH_NUM_RANDOM 100
#define EXT_RNG_KNUTH_QUALITY 1009

typedef struct {
  uint_least32_t state1[EXT_RNG_KNUTH_NUM_RANDOM];
  int_least32_t info;
  uint_least32_t state2[EXT_RNG_KNUTH_QUALITY];
} ext_rng_knuthState;

#define EXT_RNG_MERSENNE_TWISTER_NUM_RANDOM 624

typedef struct {
  int_least32_t info;
  uint_least32_t state[EXT_RNG_MERSENNE_TWISTER_NUM_RANDOM];
} ext_rng_mersenneTwisterState;

// used for EXT_RNG_ALGORITHM_USER_UNIFORM and EXT_RNG_STANDARD_NORMAL_USER_NORM
typedef struct {
  union {
    double (*stateless)(void); // used if state is NULL to avoid a pointless function call
    double (*stateful)(void*);
  } f;

  void* state;
} ext_rng_userFunction;

#ifdef __cplusplus
}
#endif

#endif // EXTERNAL_RANDOM_H

