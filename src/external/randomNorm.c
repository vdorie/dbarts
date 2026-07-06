/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 *  Copyright (C) 1997--2019  The R Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 */

#include <external/random.h>

#include <math.h> // NAN
#include <stdbool.h> // true

#include <external/io.h>
#include <external/stats.h> // qnorm

// this is duplicated in randomBase.c, randomNorm.c, and random.c
struct ext_rng {
  ext_rng_algorithm_t algorithm;
  ext_rng_standardNormal_t standardNormalAlgorithm;
  void* state;

  ext_rng_userFunction normalState;
  double gammaState[9];
};


double ext_rng_simulateStandardNormal(ext_rng* generator)
{
  switch (generator->standardNormalAlgorithm) {
    case EXT_RNG_STANDARD_NORMAL_INVERSION:
    {
#define BIG 134217728 /* 2^27 */
      /* unif_rand() alone is not of high enough precision */
      double u1 = ext_rng_simulateContinuousUniform(generator);
      u1 = (double) ((int_least32_t) (BIG * u1)) + ext_rng_simulateContinuousUniform(generator);
      return ext_quantileOfNormal(u1 / BIG, 0.0, 1.0);
    }
#undef BIG
    case EXT_RNG_STANDARD_NORMAL_USER_NORM:
    return (generator->normalState.state == NULL ? generator->normalState.f.stateless() :
              generator->normalState.f.stateful(generator->normalState.state));
    case EXT_RNG_STANDARD_NORMAL_INVALID:
    ext_throwError("unsupported standard normal generator type");
  }
  
  return NAN;
}
