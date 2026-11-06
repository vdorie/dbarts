#ifndef EXTERNAL_RNG_INTERNAL_H
#define EXTERNAL_RNG_INTERNAL_H

// Private layout of the opaque ext_rng type, shared by the RNG translation
// units (random.c, randomBase.c, randomNorm.c). It is deliberately not exposed
// in the public <external/random.h>. The gammaState[9] / normalState layout is
// load-bearing: it is serialized verbatim by ext_rng_writeSerializedState.
#include <external/random.h>

struct ext_rng {
  ext_rng_algorithm_t algorithm;
  ext_rng_standardNormal_t standardNormalAlgorithm;
  void* state;

  ext_rng_userFunction normalState;
  double gammaState[9];
};

/// Recomputes the constants of \p memo, a gammaState, from the shapes it
/// records them for, leaving a shape the gamma sampler never matches and its
/// constants as they are. The memo is a cache of those shapes, so a memo the
/// sampler wrote comes back bit for bit and one read from elsewhere is made
/// whole (random.c).
void ext_rng_rederiveGammaMemo(double* memo);

#endif // EXTERNAL_RNG_INTERNAL_H
