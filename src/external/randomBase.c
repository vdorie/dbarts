/*
 * Much of this is based off of R's RNG.c, available at r-project.org and published under
 * the GNU General Public License (http://www.r-project.org/Licenses/).
 */

#include <external/random.h>
#include "config.h"

#include <errno.h>
#include <stdbool.h>
#include <stddef.h> // size_t, malloc
#include <stdint.h> // uint32_t
#include <stdlib.h>
#include <string.h>

#ifdef HAVE_UNISTD_H
#include <unistd.h> // getpid
#elif defined(_WIN32)
#include <process.h>
#endif

// clock_gettime + CLOCK_REALTIME are in time.h, gettimeofday is in sys/time.h; plain time() is in time.h too
#if (!defined(HAVE_CLOCK_GETTIME) || !defined(CLOCK_REALTIME)) && defined(HAVE_GETTIMEOFDAY)
#  include <sys/time.h>
#else
#  include <time.h>
#endif

#include <external/alloca.h>
#include <external/io.h>

#define STANDARD_NORMAL_DEFAULT EXT_RNG_STANDARD_NORMAL_INVERSION

#if defined(__GNUC__) && (\
  (!defined(__clang__) && (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 6))) || \
  ( defined(__clang__) && (__GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 7))))
#  define SUPPRESS_DIAGNOSTIC 1
#endif

#include <external/Rinternals.h> // SEXP

// should match enum order
static const char* const rngNames[] = {
  "Wichmann-Hill",
  "Marsaglia-Multicarry",
  "Super-Duper",
  "Mersenne-Twister",
  "Knuth-TAOCP",
  "user-supplied",
  "Knuth-TAOCP-2002",
  "L'Ecuyer-CMRG",
  "invalid"
};

/* static const char* const standardNormalNames[] = {
  "Buggy Kinderman-Ramage",
  "Ahrens-Dieter",
  "Box-Muller",
  "user-supplied",
  "Inversion",
  "Kinderman-Ramage",
  "invalid"
}; */


typedef ext_rng_mersenneTwisterState MersenneTwisterState;
typedef ext_rng_knuthState KnuthState;
typedef ext_rng_userFunction UserFunction;

static const size_t stateLengths[] = {
  3 * sizeof(uint_least32_t),
  2 * sizeof(uint_least32_t),
  2 * sizeof(uint_least32_t),
  sizeof(MersenneTwisterState),
  sizeof(KnuthState),
  sizeof(UserFunction),
  sizeof(KnuthState),
  6 * sizeof(uint_least32_t)
};

// this is duplicated in randomBase.c, randomNorm.c, and random.c
struct ext_rng {
  ext_rng_algorithm_t algorithm;
  ext_rng_standardNormal_t standardNormalAlgorithm;
  void* state;
  
  union {
    double nextNormal; // used in BOX_MULLER
    ext_rng_userFunction simulateNormal;
  } normalState;
  double gammaState[9];
};

ext_rng* ext_rng_create(ext_rng_algorithm_t algorithm, const void* v_state)
{
  if (algorithm >= EXT_RNG_ALGORITHM_INVALID) {
    errno = EINVAL;
    return NULL;
  }
  
  ext_rng* result = (ext_rng*) malloc(sizeof(ext_rng));
  if (result == NULL) return NULL;
  
  result->algorithm = algorithm;
  if ((errno = ext_rng_setStandardNormalAlgorithm(result, STANDARD_NORMAL_DEFAULT, NULL)) != 0) {
    free(result);
    return NULL;
  }
  
  size_t stateLength = stateLengths[algorithm];
  
  result->state = malloc(stateLength);
  if (result->state == NULL) {
    free(result);
    return NULL;
  }
  
  if (v_state != NULL) {
    memcpy(result->state, v_state, stateLength);
  } else {
    if (algorithm == EXT_RNG_ALGORITHM_MERSENNE_TWISTER) {
      MersenneTwisterState* state = (MersenneTwisterState*) result->state;
      state->info = EXT_RNG_MERSENNE_TWISTER_NUM_RANDOM + 1;
    }
  }
  
  for (size_t i = 0; i < 9; ++i) result->gammaState[i] = 0.0;
  
  return result;
}

ext_rng* ext_rng_createDefault(bool useNative)
{
  ext_rng* result;
  
  if (useNative) {
    ext_rng_userFunction uniformFunction;
    uniformFunction.f.stateless = &unif_rand;
    uniformFunction.state = NULL;
    result = ext_rng_create(EXT_RNG_ALGORITHM_USER_UNIFORM, &uniformFunction);
    if (result == NULL) return NULL;
  
    ext_rng_userFunction normalFunction;
    normalFunction.f.stateless = &norm_rand;
    normalFunction.state = NULL;
    ext_rng_setStandardNormalAlgorithm(result, EXT_RNG_STANDARD_NORMAL_USER_NORM, &normalFunction);
    return result;
  }
  
  // if not useNative, we at least seed from native and match its type
  SEXP seedsExpr = Rf_findVarInFrame(R_GlobalEnv, R_SeedsSymbol);
  if (seedsExpr == R_UnboundValue) {
    GetRNGstate();
    PutRNGstate();
    seedsExpr = Rf_findVarInFrame(R_GlobalEnv, R_SeedsSymbol);
  }
  if (TYPEOF(seedsExpr) == PROMSXP) seedsExpr = Rf_eval(R_SeedsSymbol, R_GlobalEnv);
  
  bool seedFound = true;
  if (seedsExpr == R_UnboundValue) {
    seedFound = false;
    ext_issueWarning("seeds still unbound after calling GetRNGstate/PutRNGstate");
  } else if (!Rf_isInteger(seedsExpr)) {
    seedFound = false;
    if (seedsExpr == R_MissingArg)
      ext_issueWarning("'.Random.seed' is a missing argument with no default");
    else 
      ext_issueWarning("'.Random.seed' is not an integer vector but of type '%s', so ignored", Rf_type2char(TYPEOF(seedsExpr)));
  }
  
  if (!seedFound) {
    // use defaults
    result = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
    if (result != NULL) ext_rng_setSeedFromClock(result);
    return result;
  }  
  
  uint_least32_t seed0 = (uint_least32_t) INTEGER(seedsExpr)[0];
  
  ext_rng_algorithm_t algorithmType      = (ext_rng_algorithm_t) (seed0 % 100);
  ext_rng_standardNormal_t stdNormalType = (ext_rng_standardNormal_t) (seed0 / 100);
  
  void* state = (void*) (1 + INTEGER(seedsExpr));
  
#ifdef SUPPRESS_DIAGNOSTIC
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wswitch-enum"
#endif
  switch (algorithmType) {
    case EXT_RNG_ALGORITHM_KNUTH_TAOCP:
    case EXT_RNG_ALGORITHM_KNUTH_TAOCP2:
    {
      ext_rng_knuthState* kt = (ext_rng_knuthState*) malloc(sizeof(ext_rng_knuthState));
      if (kt == NULL) return NULL;
      memcpy(kt->state1, state, EXT_RNG_KNUTH_NUM_RANDOM * sizeof(uint_least32_t));
      kt->info = EXT_RNG_KNUTH_NUM_RANDOM; // this is a static var which we cannot access
      for (size_t i = 0; i < EXT_RNG_KNUTH_QUALITY; ++i) kt->state2[i] = 0; // also static
      state = kt;
    }
    break;
    case EXT_RNG_ALGORITHM_USER_UNIFORM:
    {
      ext_rng_userFunction* uniformFunction = (ext_rng_userFunction*) malloc(sizeof(ext_rng_userFunction));
      if (uniformFunction == NULL) return NULL;
      uniformFunction->f.stateless = &unif_rand;
      uniformFunction->state = NULL;
      state = uniformFunction;
    }
    break;
    default:
    break;
  }    
  result = ext_rng_create(algorithmType, state);
  if (algorithmType == EXT_RNG_ALGORITHM_KNUTH_TAOCP  ||
      algorithmType == EXT_RNG_ALGORITHM_KNUTH_TAOCP2 ||
      algorithmType == EXT_RNG_ALGORITHM_USER_UNIFORM) free(state);
  if (result == NULL) return NULL; 
  
  void* normalState = NULL;
  switch (stdNormalType) {
    case EXT_RNG_STANDARD_NORMAL_BOX_MULLER:
    normalState = malloc(sizeof(double));
    if (normalState == NULL) { ext_rng_destroy(result); return NULL; }
    *((double*) normalState) = 0.0; // static var, again
    break;
    case EXT_RNG_STANDARD_NORMAL_USER_NORM:
    {
      ext_rng_userFunction* normalFunction = (ext_rng_userFunction*) malloc(sizeof(ext_rng_userFunction));
      if (normalFunction == NULL) { ext_rng_destroy(result); return NULL; }
      normalFunction->f.stateless = &norm_rand;
      normalFunction->state = NULL;
      normalState = normalFunction;
    }
    break;
    default:
    break;
  }
  
#ifdef SUPPRESS_DIAGNOSTIC
#  pragma GCC diagnostic pop
#endif
  
  int errorCode = ext_rng_setStandardNormalAlgorithm(result, stdNormalType, normalState);
  if (stdNormalType == EXT_RNG_STANDARD_NORMAL_BOX_MULLER || stdNormalType == EXT_RNG_STANDARD_NORMAL_USER_NORM) free(normalState);
  
  if (errorCode != 0) {
    ext_rng_destroy(result);
    return NULL;
  }
  
  return result;
}


void ext_rng_destroy(ext_rng* generator)
{
  if (generator == NULL) return;
  
  if (generator->state != NULL) free(generator->state);
  free(generator);
}

int ext_rng_createAndSeed(ext_rng** result, ext_rng_algorithm_t algorithm, ext_rng_standardNormal_t standardNormalAlgorithm)
{
  if ((*result = ext_rng_create(algorithm, NULL)) == NULL)
    return ENOMEM;
  
  if (ext_rng_setStandardNormalAlgorithm(*result, standardNormalAlgorithm, NULL) != 0) {
    ext_rng_destroy(*result);
    *result = NULL;
    return EINVAL;
  }
  
  int errorCode = ext_rng_setSeedFromClock(*result);
  if (errorCode != 0) {
    ext_rng_destroy(*result);
    *result = NULL;
  }
  return errorCode;
}

ext_rng_algorithm_t ext_rng_getDefaultAlgorithmType()
{
  SEXP seedsExpr = Rf_findVarInFrame(R_GlobalEnv, R_SeedsSymbol);
  if (seedsExpr == R_UnboundValue) {
    GetRNGstate();
    PutRNGstate();
    seedsExpr = Rf_findVarInFrame(R_GlobalEnv, R_SeedsSymbol);
  }
  if (TYPEOF(seedsExpr) == PROMSXP) seedsExpr = Rf_eval(R_SeedsSymbol, R_GlobalEnv);
  
  if (seedsExpr == R_UnboundValue || !Rf_isInteger(seedsExpr))
    return EXT_RNG_ALGORITHM_MERSENNE_TWISTER;
  
  uint_least32_t seed0 = (uint_least32_t) INTEGER(seedsExpr)[0];
  
  return (ext_rng_algorithm_t) (seed0 % 100);
}

ext_rng_standardNormal_t ext_rng_getDefaultStandardNormalType()
{
  SEXP seedsExpr = Rf_findVarInFrame(R_GlobalEnv, R_SeedsSymbol);
  if (seedsExpr == R_UnboundValue) {
    GetRNGstate();
    PutRNGstate();
    seedsExpr = Rf_findVarInFrame(R_GlobalEnv, R_SeedsSymbol);
  }
  if (TYPEOF(seedsExpr) == PROMSXP) seedsExpr = Rf_eval(R_SeedsSymbol, R_GlobalEnv);
  
  if (seedsExpr == R_UnboundValue || !Rf_isInteger(seedsExpr))
    return EXT_RNG_STANDARD_NORMAL_INVERSION;
  
  uint_least32_t seed0 = (uint_least32_t) INTEGER(seedsExpr)[0];
  
  return (ext_rng_standardNormal_t) (seed0 / 100);
}

int ext_rng_setStandardNormalAlgorithm(ext_rng* generator, ext_rng_standardNormal_t standardNormalAlgorithm, const void* state)
{
  if (generator == NULL) return EFAULT;
  
  if (standardNormalAlgorithm >= EXT_RNG_STANDARD_NORMAL_INVALID) return EINVAL;
  
  generator->standardNormalAlgorithm = standardNormalAlgorithm;
  
  if (standardNormalAlgorithm == EXT_RNG_STANDARD_NORMAL_BOX_MULLER) {
    generator->normalState.nextNormal = (state != NULL ? *((double*) state) : 0.0);
  } else if (standardNormalAlgorithm == EXT_RNG_STANDARD_NORMAL_USER_NORM) {
    if (state != NULL) {
      memcpy(&generator->normalState, state, sizeof(ext_rng_userFunction));
    } else {
      return EINVAL;
    }
  }
  
  return 0;
}


#define LECUYER_M1 4294967087
#define LECUYER_M2 4294944443
static void validateSeed(ext_rng* generator, bool isFirstRun);
static void knuth_setSeed(KnuthState* kt, uint_least32_t seed);
static void knuth2_setSeed(KnuthState* kt, uint_least32_t seed);

int ext_rng_setSeed(ext_rng* generator, uint_least32_t seed)
{
  if (generator == NULL) return EFAULT;
  
  size_t stateLength = stateLengths[generator->algorithm];
  uint_least32_t* state = (uint_least32_t*) generator->state;
  if (generator->standardNormalAlgorithm == EXT_RNG_STANDARD_NORMAL_USER_NORM) generator->normalState.nextNormal = 0.0;

  // initial scrambling
  for (size_t j = 0; j < 50; ++j) seed = (69069 * seed + 1);
  
  switch (generator->algorithm) {
    case EXT_RNG_ALGORITHM_WICHMANN_HILL:
    case EXT_RNG_ALGORITHM_MARSAGLIA_MULTICARRY:
    case EXT_RNG_ALGORITHM_SUPER_DUPER:
    case EXT_RNG_ALGORITHM_MERSENNE_TWISTER:
    {
      // this is broken for MT, as it treats the info slot as any-old-int, but is necessary
      // for backwards compatibility
      size_t numInts = stateLength / sizeof(uint_least32_t);
      for (size_t j = 0; j < numInts; ++j) {
        seed = (69069 * seed + 1);
        state[j] = seed;
      }
      validateSeed(generator, true);
    }
    break;
    case EXT_RNG_ALGORITHM_KNUTH_TAOCP:
    knuth_setSeed((KnuthState*) generator->state, seed);
    break;
    case EXT_RNG_ALGORITHM_KNUTH_TAOCP2:
    knuth2_setSeed((KnuthState*) generator->state, seed);
    break;
    case EXT_RNG_ALGORITHM_LECUYER_CMRG:
    {
      size_t numInts = stateLength / sizeof(uint_least32_t);
      for (size_t j = 0; j < numInts; ++j) {
        seed = (69069 * seed + 1);
        while (seed >= LECUYER_M2) seed = (69069 * seed + 1);
        state[j] = seed;
      }
    }
    break;
    case EXT_RNG_ALGORITHM_USER_UNIFORM:
    case EXT_RNG_ALGORITHM_INVALID:
    return EINVAL;
  }
  
  return 0;
}

int ext_rng_setSeedFromClock(ext_rng* generator)
{
  uint_least32_t seed;
  unsigned int pid = (unsigned int) getpid();
  
#if defined(HAVE_CLOCK_GETTIME) && defined(CLOCK_REALTIME)
  struct timespec tp;
  clock_gettime(CLOCK_REALTIME, &tp);
  seed = (uint_least32_t) (((uint_least64_t) tp.tv_nsec << 16) ^ tp.tv_sec);
#elif defined(HAVE_GETTIMEOFDAY)
  struct timeval tv;
  gettimeofday(&tv, NULL);
  seed = (uint_least32_t) (((uint_least64_t) tv.tv_usec << 16) ^ tv.tv_sec);
#else
  // C89, so must work
  seed = (uint_least32_t) time(NULL);
#endif
  seed ^= (pid << 16);
  
  return ext_rng_setSeed(generator, seed);
}

size_t ext_rng_getSerializedStateLength(const ext_rng* generator)
{
  if (generator == NULL) return(0);
  
  size_t stateLength = 0;
  
  // can't serialize a function pointer and a state whose length we don't know
  if (generator->algorithm != EXT_RNG_ALGORITHM_USER_UNIFORM) {
    stateLength = stateLengths[generator->algorithm];
    stateLength += stateLength % sizeof(int);
  }
  
  if (generator->standardNormalAlgorithm == EXT_RNG_STANDARD_NORMAL_BOX_MULLER) {
    stateLength += sizeof(double);
    stateLength += stateLength % sizeof(int);
  }
  
  stateLength += sizeof(((ext_rng*) 0)->gammaState);
  stateLength += stateLength % sizeof(int);
  
  return(stateLength);
}

void ext_rng_writeSerializedState(const ext_rng* generator, void* state)
{
  if (generator->algorithm != EXT_RNG_ALGORITHM_USER_UNIFORM) {
    size_t stateLength = stateLengths[generator->algorithm];
    
    memcpy(state, generator->state, stateLength);
    
    state = (void*) ((char*) state + stateLength);
    
    // fill out any 0s if necessary
    for (size_t i = 0; i < stateLength % sizeof(int); ++i) {
      *((char*) state) = 0;
      state = (void*) ((char*) state + 1);
    }
  }
  if (generator->standardNormalAlgorithm == EXT_RNG_STANDARD_NORMAL_BOX_MULLER) {
    memcpy(state, &generator->normalState, sizeof(double));
    
    state = (void*) ((char*) state + sizeof(double));
    
#if defined(__clang__) && (__GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 5))
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Wtautological-compare"
#endif
    for (size_t i = 0; i < sizeof(double) % sizeof(int); ++i) {
      *((char*) state) = 0;
      state = (void*) ((char*) state + 1);
    }
  }

  memcpy(state, &generator->gammaState, sizeof(((ext_rng*) 0)->gammaState));
  state = (void*) ((char*) state + sizeof(((ext_rng*) 0)->gammaState));
  
  for (size_t i = 0; i < sizeof(((ext_rng*) 0)->gammaState) % sizeof(int); ++i) {
     *((char*) state) = 0;
     state = (void*) ((char*) state + 1);
  }
#if defined(__clang__) && (__GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 5))
#  pragma clang diagnostic pop
#endif
}

void ext_rng_readSerializedState(ext_rng* generator, const void* state)
{
  if (generator->algorithm != EXT_RNG_ALGORITHM_USER_UNIFORM) {
    size_t stateLength = stateLengths[generator->algorithm];
    
    memcpy(generator->state, state, stateLength);
    
    state = (const void*) ((const char*) state + stateLength + stateLength % sizeof(int));
  }
  if (generator->standardNormalAlgorithm == EXT_RNG_STANDARD_NORMAL_BOX_MULLER) {
    memcpy(&generator->normalState, state, sizeof(double));
    
    state = (void*) ((char*) state + sizeof(double) + sizeof(double) % sizeof(int));
  }

  memcpy(&generator->gammaState, state, sizeof(((ext_rng*) 0)->gammaState));
}

static void validateSeed(ext_rng* generator, bool isFirstRun)
{
  bool notAllAreZero = false;

  uint_least32_t* state = (uint_least32_t*) generator->state;

  switch (generator->algorithm) {
    case EXT_RNG_ALGORITHM_WICHMANN_HILL:
    state[0] = state[0] % 30269;
    state[1] = state[1] % 30307;
    state[2] = state[2] % 30323;
    
  	// map values equal to 0 mod modulus to 1
  	if (state[0] == 0) state[0] = 1;
  	if (state[1] == 0) state[1] = 1;
  	if (state[2] == 0) state[2] = 1;
    
    break;
    
    case EXT_RNG_ALGORITHM_SUPER_DUPER:
  	if (state[0] == 0) state[0] = 1;
  	// state[1] = Congruential: must be ODD
  	state[1] |= 1;
    break;
    


    case EXT_RNG_ALGORITHM_MARSAGLIA_MULTICARRY:
    if (state[0] == 0) state[0] = 1;
    if (state[1] == 0) state[1] = 1;
    break;
    
    case EXT_RNG_ALGORITHM_MERSENNE_TWISTER:
    {
      MersenneTwisterState* mt = (MersenneTwisterState*) generator->state;
      
      if (isFirstRun || mt->info <= 0) mt->info = EXT_RNG_MERSENNE_TWISTER_NUM_RANDOM;
    
      for (size_t j = 0; j < EXT_RNG_MERSENNE_TWISTER_NUM_RANDOM; ++j) {
	      if (mt->state[j] != 0) {
		      notAllAreZero = true; break;
        }
      }
      if (!notAllAreZero) ext_rng_setSeedFromClock(generator);
    }
    break;
    
    case EXT_RNG_ALGORITHM_KNUTH_TAOCP:
    case EXT_RNG_ALGORITHM_KNUTH_TAOCP2:
    {
      KnuthState* kt = (KnuthState*) generator->state;
    
      if (kt->info <= 0) kt->info = EXT_RNG_KNUTH_NUM_RANDOM;
	    for (size_t j = 0; j < EXT_RNG_KNUTH_NUM_RANDOM; ++j) {
        if (state[j] != 0) {
          notAllAreZero = true; break;
        }
      }
      if (!notAllAreZero) ext_rng_setSeedFromClock(generator);
    }
    break;
    case EXT_RNG_ALGORITHM_LECUYER_CMRG:
    // first set: not all zero, in [0, m1)
    // second set: not all zero, in [0, m2)
    {
      uint_least32_t temp;
      bool allAreValid = true;
      for (size_t j = 0; j < 3; ++j) {
        temp = state[j];
        if (temp != 0) notAllAreZero = true;
        if (temp >= LECUYER_M1) allAreValid = false;
	    }
      if (!notAllAreZero || !allAreValid) ext_rng_setSeedFromClock(generator);
      for (size_t j = 3; j < 6; ++j) {
        temp = state[j];
        if (temp != 0) notAllAreZero = true;
        if (temp >= LECUYER_M2) allAreValid = false;
      }
      if (!notAllAreZero || !allAreValid) ext_rng_setSeedFromClock(generator);
    }
    break;
    case EXT_RNG_ALGORITHM_USER_UNIFORM:
    case EXT_RNG_ALGORITHM_INVALID:
    break;
  }
}

// #define THIRTY_TWO_BIT_MAX     4294967296.0
#define THIRTY_TWO_BIT_INVERSE 2.328306437080797e-10 /* = 1/(2^32 - 1) */
#define KNUTH_CONSTANT         9.31322574615479e-10

#ifdef SUPPRESS_DIAGNOSTIC
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wuninitialized"
#  if defined(__GNUC__) && !defined(__clang__) && (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7))
#    pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#  endif
#endif
// guarantees results in (0, 1)
inline static double truncateToUnitInterval(double x) {
  if (x <= 0.0) return 0.5 * THIRTY_TWO_BIT_INVERSE;
  if ((1.0 - x) <= 0.0) return 1.0 - 0.5 * THIRTY_TWO_BIT_INVERSE;
  return x;
}

static double mersenneTwister_getNext(MersenneTwisterState* mt);
static uint_least32_t knuth_getNext(KnuthState* kt);

double ext_rng_simulateContinuousUniform(ext_rng* generator)
{
  double result; // for some reason GCC thinks this might get used unitialized, but the switch is exhaustive
  
  uint_least32_t* state = (uint_least32_t*) generator->state;
  
  switch (generator->algorithm) {
    case EXT_RNG_ALGORITHM_WICHMANN_HILL:
    {
      state[0] = state[0] * 171 % 30269;
      state[1] = state[1] * 172 % 30307;
    	state[2] = state[2] * 170 % 30323;
      result = state[0] / 30269.0 + state[1] / 30307.0 + state[2] / 30323.0;
      result -= (int) result;
    }
    break;
    case EXT_RNG_ALGORITHM_MARSAGLIA_MULTICARRY: // 0177777(octal) == 65535(decimal)
    {
	    state[0] = 36969 * (state[0] & 0177777) + (state[0] >> 16);
      state[1] = 18000 * (state[1] & 0177777) + (state[1] >> 16);
      result = ((double) ((state[0] << 16) ^ (state[1] & 0177777))) * THIRTY_TWO_BIT_INVERSE;
    }
    break;
    case EXT_RNG_ALGORITHM_SUPER_DUPER:
	    /* This is Reeds et al (1984) implementation;
	     * modified using __unsigned__ states instead of signed ones
	      */
    {
      state[0] ^= ((state[0] >> 15) & 0377777); /* Tausworthe */
    	state[0] ^= state[0] << 17;
    	state[1] *= 69069;  /* Congruential */
      result = (state[0] ^ state[1]) * THIRTY_TWO_BIT_INVERSE;
    }
    break;
    case EXT_RNG_ALGORITHM_MERSENNE_TWISTER:
    {
      result = mersenneTwister_getNext((MersenneTwisterState*) generator->state);
    }
    break;
    case EXT_RNG_ALGORITHM_KNUTH_TAOCP:
    case EXT_RNG_ALGORITHM_KNUTH_TAOCP2:
    {
      result = knuth_getNext((KnuthState*) generator->state) * KNUTH_CONSTANT;
    }
    break;
    case EXT_RNG_ALGORITHM_USER_UNIFORM:
    {
      UserFunction* function = (UserFunction*) generator->state;
      result = (function->state == NULL ? function->f.stateless() : function->f.stateful(function->state));
    }
    break;
    case EXT_RNG_ALGORITHM_LECUYER_CMRG:
    {
      /* Based loosely on the GPL-ed version of
         http://www.iro.umontreal.ca/~lecuyer/myftp/streams00/c2010/RngStream.c
         but using int_least64_t, which C99 guarantees.
      */
      
      int_least32_t k;
      int_least64_t p1, p2;

#define normc  2.328306549295727688e-10
#define a12     (int_least64_t) 1403580
#define a13n    (int_least64_t) 810728
#define a21     (int_least64_t) 527612
#define a23n    (int_least64_t) 1370589

	    p1 = a12 * state[1] - a13n * state[0];
	    /* p1 % m1 would surely do */
	    k = (int_least32_t) (p1 / LECUYER_M1);
	    p1 -= k * LECUYER_M1;
	    if (p1 < 0) p1 += LECUYER_M1;
      state[0] = state[1]; state[1] = state[2]; state[2] = (uint_least32_t) p1;
      
      p2 = a21 * state[5] - a23n * state[3];
      k = (int_least32_t) (p2 / LECUYER_M2);
      p2 -= k * LECUYER_M2;
      if (p2 < 0) p2 += LECUYER_M2;
      state[3] = state[4]; state[4] = state[5]; state[5] = (uint_least32_t) p2;

      result = (double) ((p1 > p2) ? (p1 - p2) : (p1 - p2 + LECUYER_M1)) * normc;
#undef normc
#undef a12
#undef a13n
#undef a21
#undef a23n
    }
    break;
    case EXT_RNG_ALGORITHM_INVALID:
	  ext_throwError("ext_rng_simulateContinuousUniform: unimplemented rng kind %s", rngNames[generator->algorithm]);
  }

  return truncateToUnitInterval(result);
}

#ifdef SUPPRESS_DIAGNOSTIC
#  pragma GCC diagnostic pop
#endif

/* ===================  Mersenne Twister ========================== */
/* From http://www.math.keio.ac.jp/~matumoto/emt.html */

/* A C-program for MT19937: Real number version([0,1)-interval)
   (1999/10/28)
     genrand() generates one pseudorandom real number (double)
   which is uniformly distributed on [0,1)-interval, for each
   call. sgenrand(seed) sets initial values to the working area
   of 624 words. Before genrand(), sgenrand(seed) must be
   called once. (seed is any 32-bit integer.)
   Integer generator is obtained by modifying two lines.
     Coded by Takuji Nishimura, considering the suggestions by
   Topher Cooper and Marc Rieffel in July-Aug. 1997.

   Copyright (C) 1997, 1999 Makoto Matsumoto and Takuji Nishimura.
   When you use this, send an email to: matumoto@math.keio.ac.jp
   with an appropriate reference to your work.

   REFERENCE
   M. Matsumoto and T. Nishimura,
   "Mersenne Twister: A 623-Dimensionally Equidistributed Uniform
   Pseudo-Random Number Generator",
   ACM Transactions on Modeling and Computer Simulation,
   Vol. 8, No. 1, January 1998, pp 3--30.
*/

/* Period parameters */
#define N EXT_RNG_MERSENNE_TWISTER_NUM_RANDOM
#define M 397
#define MmN (397 - EXT_RNG_MERSENNE_TWISTER_NUM_RANDOM)
#define MATRIX_A 0x9908b0df   /* constant vector a */
#define UPPER_MASK 0x80000000 /* most significant w-r bits */
#define LOWER_MASK 0x7fffffff /* least significant r bits */

/* Tempering parameters */
#define TEMPERING_MASK_B 0x9d2c5680
#define TEMPERING_MASK_C 0xefc60000
#define TEMPERING_SHIFT_U(y)  (y >> 11)
#define TEMPERING_SHIFT_S(y)  (y << 7)
#define TEMPERING_SHIFT_T(y)  (y << 15)
#define TEMPERING_SHIFT_L(y)  (y >> 18)


/* Initializing the array with a seed */
static void
mersenneTwister_generateStateFromSeed(MersenneTwisterState* mt, uint_least32_t seed)
{
  for (uint_least32_t i = 0; i < N; ++i) {
    mt->state[i] = seed & 0xffff0000;
    seed = 69069 * seed + 1;
    mt->state[i] |= (seed & 0xffff0000) >> 16;
    seed = 69069 * seed + 1;
  }
  mt->info = N;
}

/* Initialization by "sgenrand()" is an example. Theoretically,
   there are 2^19937-1 possible states as an intial state.
   Essential bits in "seed_array[]" is following 19937 bits:
    (seed_array[0]&UPPER_MASK), seed_array[1], ..., seed_array[N-1].
   (seed_array[0]&LOWER_MASK) is discarded.
   Theoretically,
    (seed_array[0]&UPPER_MASK), seed_array[1], ..., seed_array[N-1]
   can take any values except all zeros. */

static const uint_least32_t mag01[2] = { 0x0, MATRIX_A };
// mag01[x] = x * MATRIX_A  for x = 0,1
  
static double mersenneTwister_getNext(MersenneTwisterState* mt)
{
  uint_least32_t y;
  
  uint_least32_t* state = mt->state;

  if (mt->info >= N) { /* generate N words at one time */
	  size_t kk;
    
    if (mt->info == N + 1)  /* if sgenrand() has not been called, */
	    mersenneTwister_generateStateFromSeed(mt, 4357); /* a default initial seed is used */
    
    for (kk = 0; kk < N - M; kk++) {
	    y = (state[kk] & UPPER_MASK) | (state[kk + 1] & LOWER_MASK);
	    state[kk] = state[kk + M] ^ (y >> 1) ^ mag01[y & 0x1];
	  }
    for (; kk < N - 1; kk++) {
	    y = (state[kk] & UPPER_MASK) | (state[kk + 1] & LOWER_MASK);
	    state[kk] = state[(int_least32_t) kk + MmN] ^ (y >> 1) ^ mag01[y & 0x1];
    }
    y = (state[N - 1] & UPPER_MASK) | (state[0] & LOWER_MASK);
    state[N - 1] = state[M - 1] ^ (y >> 1) ^ mag01[y & 0x1];

    mt->info = 0;
  }
  
  y = state[mt->info++];
  y ^= TEMPERING_SHIFT_U(y);
  y ^= TEMPERING_SHIFT_S(y) & TEMPERING_MASK_B;
  y ^= TEMPERING_SHIFT_T(y) & TEMPERING_MASK_C;
  y ^= TEMPERING_SHIFT_L(y);

  return ((double) y * 2.3283064365386963e-10); /* reals: [0,1)-interval */
}

#undef N
#undef M
#undef MATRIX_A
#undef UPPER_MASK
#undef LOWER_MASK

#undef TEMPERING_MASK_B
#undef TEMPERING_MASK_C
#undef TEMPERING_SHIFT_U
#undef TEMPERING_SHIFT_S
#undef TEMPERING_SHIFT_T
#undef TEMPERING_SHIFT_L


/*
   The following code was taken from earlier versions of
   http://www-cs-faculty.stanford.edu/~knuth/programs/rng.c-old
   http://www-cs-faculty.stanford.edu/~knuth/programs/rng.c
*/

/* ===================  Knuth TAOCP  2002 ========================== */

/*    This program by D E Knuth is in the public domain and freely copyable.
 *    It is explained in Seminumerical Algorithms, 3rd edition, Section 3.6
 *    (or in the errata to the 2nd edition --- see
 *        http://www-cs-faculty.stanford.edu/~knuth/taocp.html
 *    in the changes to Volume 2 on pages 171 and following).              */

/*    N.B. The MODIFICATIONS introduced in the 9th printing (2002) are
      included here; there's no backwards compatibility with the original. */

#define longLag EXT_RNG_KNUTH_NUM_RANDOM
#define shortLag 37
#define modulus 1073741824 // 2^30
#define KKK (longLag + longLag - 1) // someone really needs better variable names than KKK
#define KKL (longLag - shortLag)
#define differenceModulo(_X_, _Y_) (((_X_) - (_Y_)) & (modulus - 1))
#define isOdd(_X_) ((_X_) & 1)

#define TT 70 // guaranteed separation between streams

static void knuth_setSeed(KnuthState* kt, uint_least32_t seed)
{
  seed %= 1073741821;
  
  uint_least32_t ss = seed - (seed % 2) + 2;
  
  uint_least32_t temp[KKK];
    
  for (size_t j = 0; j < longLag; ++j) {
    temp[j] = ss;
    ss <<= 1;
    if (ss >= modulus) ss -= modulus - 2;
  }
  temp[1]++;
  
  for (size_t j = longLag; j < KKK; ++j) temp[j] = 0;
  
  ss = seed;
  size_t t = TT - 1;
  while (t > 0) {
    for (size_t j = longLag - 1; j > 0; --j) temp[j + j] = temp[j];
    
    for (size_t j = KKK - 1; j >= KKL + 1; j -= 2) temp[KKK - j] = temp[j] - (temp[j] % 2);
    
    for (size_t j = KKK - 1; j >= longLag; --j) if (temp[j] % 2 == 1) {
      temp[j - KKL] = differenceModulo(temp[j - KKL], temp[j]);
      temp[j - longLag] = differenceModulo(temp[j - longLag], temp[j]);
    }
    if (isOdd(ss)) {
      for (size_t j = longLag; j > 0; --j) temp[j] = temp[j - 1];
      temp[0] = temp[longLag];
      if (temp[longLag] % 2 == 1)
        temp[shortLag] =  differenceModulo(temp[shortLag], temp[longLag]);
    }
    if (ss != 0) ss /= 2;
    else t--;
  }
  
  uint_least32_t* state = kt->state1;
  
  memcpy(state, temp + shortLag, KKL * sizeof(uint_least32_t));
  memcpy(state + KKL, temp, shortLag * sizeof(uint_least32_t));
  
  kt->info = EXT_RNG_KNUTH_NUM_RANDOM;
}

static void knuth_randomizeArray(KnuthState* kt, uint_least32_t* array, size_t length);

static void knuth2_setSeed(KnuthState* kt, uint_least32_t seed)
{
  seed %= 1073741821;
  
  uint_least32_t* temp = ext_stackAllocate(longLag + longLag - 1, uint_least32_t);
  
  uint_least32_t ss = (seed + 2) & (modulus - 2);
  
  for (size_t j = 0; j < longLag; ++j) {
    temp[j] = ss;
    ss <<= 1;
    if (ss >= modulus) ss -= modulus - 2;
  }
  temp[1]++;
  
  ss = seed & (modulus - 1);
  size_t t = TT - 1;
  while (t > 0) {
    for (size_t j = longLag - 1; j > 0; --j) { temp[j + j] = temp[j]; temp[j + j - 1] = 0; }
    for (size_t j = KKK - 1; j >= longLag; --j) {
      temp[j - KKL] = differenceModulo(temp[j - KKL], temp[j]);
      temp[j - longLag] = differenceModulo(temp[j - longLag], temp[j]);
    }
    if (isOdd(ss)) {
      for (size_t j = longLag; j > 0; --j) temp[j] = temp[j - 1];
      temp[0] = temp[longLag];
      temp[shortLag] = differenceModulo(temp[shortLag], temp[longLag]);
    }
    if (ss != 0) ss /= 2;
    else t--;
  }
  
  uint_least32_t* state = kt->state1;
    
  memcpy(state, temp + shortLag, KKL * sizeof(uint_least32_t));
  memcpy(state + KKL, temp, shortLag * sizeof(uint_least32_t));
  
  for (size_t j = 0; j < 10; ++j) knuth_randomizeArray(kt, temp, KKK);
  
  ext_stackFree(temp);
  
  kt->info = EXT_RNG_KNUTH_NUM_RANDOM;
}

static void knuth_randomizeArray(KnuthState* kt, uint_least32_t* array, size_t length)
{
  size_t i, j;
  
  for (j = 0; j < longLag; ++j) array[j] = kt->state1[j];
  for ( ; j < length; ++j) array[j] = differenceModulo(array[j - longLag], array[j - shortLag]);
  
  for (i = 0; i < shortLag; ++i, ++j) kt->state1[i] = differenceModulo(array[j - longLag], array[j - shortLag]);
  for ( ; i < longLag; ++i, ++j) kt->state1[i] = differenceModulo(array[j - longLag], kt->state1[i - shortLag]);
}

static void knuth_cycleArray(KnuthState* kt);

static uint_least32_t knuth_getNext(KnuthState* kt)
{
  if (kt->info >= EXT_RNG_KNUTH_NUM_RANDOM) {
    knuth_cycleArray(kt);
    kt->info = 0;
  }
  
  return kt->state1[kt->info++];
}

static void knuth_cycleArray(KnuthState* kt)
{
  knuth_randomizeArray(kt, kt->state2, EXT_RNG_KNUTH_QUALITY);
  kt->state2[longLag] = (uint_least32_t) -1;
}


#undef longLag
#undef shortLag
#undef modulus
#undef KKK
#undef KKL
#undef differenceModulo
#undef isOdd
