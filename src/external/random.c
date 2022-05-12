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

/*
 * drawFromDiscreteDistribution and below are original to this project
 * and inherit its copyright, also GPL-2.
 */

#include <external/random.h>

#include <math.h> // exp, log, expm1, fabs, nan
#include <stdbool.h>

#if (defined(__clang__) && (__clang_major__ > 3 || (__clang_major__ == 3 && __clang_minor__ >= 7))) || \
    (defined(__GNUC__) && (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 6)))
#  define SUPPRESS_DIAGNOSTIC 1
#endif

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

static double simulateStandardExponential(ext_rng* generator);

double ext_rng_simulateExponential(ext_rng* generator, double scale)
{
  if (!isfinite(scale) || scale <= 0.0) {
    return (scale == 0.0) ? 0.0 : NAN;
  }
  
  return simulateStandardExponential(generator) * scale;
}

double ext_rng_simulateGamma(ext_rng* generator, double shape, double scale)
{
  static const double sqrt32 = 5.656854;
  static const double exp_m1 = 0.36787944117144232159; /* exp(-1) = 1/e */

  /* Coefficients q[k] - for q0 = sum(q[k]*a^(-k))
   * Coefficients a[k] - for q = q0+(t*t/2)*sum(a[k]*v^k)
   * Coefficients e[k] - for exp(q)-1 = sum(e[k]*q^k)
   */
  static const double q1 = 0.04166669;
  static const double q2 = 0.02083148;
  static const double q3 = 0.00801191;
  static const double q4 = 0.00144121;
  static const double q5 = -7.388e-5;
  static const double q6 = 2.4511e-4;
  static const double q7 = 2.424e-4;

  static const double a1 = 0.3333333;
  static const double a2 = -0.250003;
  static const double a3 = 0.2000062;
  static const double a4 = -0.1662921;
  static const double a5 = 0.1423657;
  static const double a6 = -0.1367177;
  static const double a7 = 0.1233795;

#define oldShape (generator->gammaState[0])
#define oldShape2 (generator->gammaState[1])
#define s   (generator->gammaState[2]) // no. 1 (step 1)
#define s2  (generator->gammaState[3]) //
#define d   (generator->gammaState[4]) //
#define q0  (generator->gammaState[5]) // no. 2 (step 4) */
#define b   (generator->gammaState[6]) //
#define si  (generator->gammaState[7]) //
#define c   (generator->gammaState[8]) //
  
  double e, p, q, r, t, u, v, w, x, ret_val;

  if (!isfinite(shape) || !isfinite(scale) || shape < 0.0 || scale <= 0.0)
    return (scale == 0.0) ? 0.0 : NAN;
  
  if (shape < 1.0) { /* GS algorithm for parameters a < 1 */
    if (shape == 0.0) return 0.0;
    e = 1.0 + exp_m1 * shape;
    do {
      p = e * ext_rng_simulateContinuousUniform(generator);
      if (p >= 1.0) {
        x = -log((e - p) / shape);
        if (simulateStandardExponential(generator) >= (1.0 - shape) * log(x))
          break;
      } else {
        x = exp(log(p) / shape);
        if (simulateStandardExponential(generator) >= x)
          break;
      }
    } while (true);
    
    return scale * x;
  }

  /* --- a >= 1 : GD algorithm --- */
  
  /* Step 1: Recalculations of s2, s, d if a has changed */
#ifdef SUPPRESS_DIAGNOSTIC
#  ifdef __clang__
#    pragma clang diagnostic push
#    pragma clang diagnostic ignored "-Wfloat-equal"
#  else
#    pragma GCC diagnostic push
#    pragma GCC diagnostic ignored "-Wfloat-equal"
#  endif
#endif
  if (shape != oldShape) {
#ifdef SUPPRESS_DIAGNOSTIC
#  ifdef __clang__
#    pragma clang diagnostic pop
#  else
#    pragma GCC diagnostic pop
#  endif
#endif
    oldShape = shape;
    s2 = shape - 0.5;
    s = sqrt(s2);
    d = sqrt32 - s * 12.0;
  }
  /* Step 2: t = standard normal deviate,
             x = (s,1/2) -normal deviate. */

  /* immediate acceptance (i) */
  t = ext_rng_simulateStandardNormal(generator);
  x = s + 0.5 * t;
  ret_val = x * x;
  if (t >= 0.0) return scale * ret_val;

  /* Step 3: u = 0,1 - uniform sample. squeeze acceptance (s) */
  u = ext_rng_simulateContinuousUniform(generator);
  if (d * u <= t * t * t) return scale * ret_val;

  /* Step 4: recalculations of q0, b, si, c if necessary */
#ifdef SUPPRESS_DIAGNOSTIC
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wfloat-equal"
#endif
  if (shape != oldShape2) {
#ifdef SUPPRESS_DIAGNOSTIC
#  pragma GCC diagnostic pop
#endif
    oldShape2 = shape;
    r = 1.0 / shape;
    q0 = ((((((q7 * r + q6) * r + q5) * r + q4) * r + q3) * r
	         + q2) * r + q1) * r;

    /* Approximation depending on size of parameter a */
    /* The constants in the expressions for b, si and c */
    /* were established by numerical experiments */
    
    if (shape <= 3.686) {
      b = 0.463 + s + 0.178 * s2;
      si = 1.235;
      c = 0.195 / s - 0.079 + 0.16 * s;
    } else if (shape <= 13.022) {
      b = 1.654 + 0.0076 * s2;
      si = 1.68 / s + 0.275;
      c = 0.062 / s + 0.024;
    } else {
      b = 1.77;
      si = 0.75;
      c = 0.1515 / s;
    }
  }
  
  /* Step 5: no quotient test if x not positive */
  if (x > 0.0) {
    
    /* Step 6: calculation of v and quotient q */
    v = t / (s + s);
    if (fabs(v) <= 0.25)
      q = q0 + 0.5 * t * t * ((((((a7 * v + a6) * v + a5) * v + a4) * v
                                + a3) * v + a2) * v + a1) * v;
    else
      q = q0 - s * t + 0.25 * t * t + (s2 + s2) * log(1.0 + v);
    
    
    /* Step 7: quotient acceptance (q) */
    if (log(1.0 - u) <= q) return scale * ret_val;
  }

  do {
    /* Step 8: e = standard exponential deviate
     *	u =  0,1 -uniform deviate
     *	t = (b,si)-double exponential (laplace) sample */
    e = simulateStandardExponential(generator);
    u = ext_rng_simulateContinuousUniform(generator);
    u = u + u - 1.0;
    if (u < 0.0) t = b - si * e;
    else t = b + si * e;
    /* Step	 9:  rejection if t < tau(1) = -0.71874483771719 */
    if (t >= -0.71874483771719) {
      /* Step 10:	 calculation of v and quotient q */
      v = t / (s + s);
      if (fabs(v) <= 0.25)
        q = q0 + 0.5 * t * t * ((((((a7 * v + a6) * v + a5) * v + a4) * v + a3) * v + a2) * v + a1) * v;
      else
        q = q0 - s * t + 0.25 * t * t + (s2 + s2) * log(1.0 + v);
      /* Step 11:	 hat acceptance (h) */
      /* (if q not positive go to step 8) */
      if (q > 0.0) {
        w = expm1(q);
        /*  ^^^^^ original code had approximation with rel.err < 2e-7 */
        /* if t is rejected sample again at step 8 */
        if (c * fabs(u) <= w * exp(e - 0.5 * t * t))
          break;
      }
    }
  } while (true);
  x = s + 0.5 * t;
  return scale * x * x;
}

#undef oldShape
#undef oldShape2
#undef s
#undef s2
#undef d
#undef q0
#undef b
#undef si
#undef c


static double simulateStandardExponential(ext_rng* generator)
{ 
  /*  REFERENCE
   *
   *    Ahrens, J.H. and Dieter, U. (1972).
   *    Computer methods for sampling from the exponential and
   *    normal distributions.
   *    Comm. ACM, 15, 873-882.
   */
  
    /* q[k-1] = sum(log(2)^k / k!)  k=1,..,n, */
    /* The highest n (here 16) is determined by q[n-1] = 1.0 */
    /* within standard precision */
  static const double q[] = {
    0.6931471805599453, 0.9333736875190459, 0.9888777961838675, 0.9984959252914960,
    0.9998292811061389, 0.9999833164100727, 0.9999985691438767, 0.9999998906925558,
    0.9999999924734159, 0.9999999995283275, 0.9999999999728814, 0.9999999999985598,
    0.9999999999999289, 0.9999999999999968, 0.9999999999999999, 1.0000000000000000 };

  double a = 0.;
  double u = ext_rng_simulateContinuousUniform(generator); /* precaution if u = 0 is ever returned */
  while (u <= 0.0 || u >= 1.0) u = ext_rng_simulateContinuousUniform(generator);
  do {
    u += u;
    if (u > 1.0) break;
    a += q[0];
  } while (true);
  u -= 1.0;

  if (u <= q[0]) return a + u;

  int_least32_t i = 0;
  double ustar = ext_rng_simulateContinuousUniform(generator);
  double umin = ustar;
  do {
    ustar = ext_rng_simulateContinuousUniform(generator);
    if (umin > ustar) umin = ustar;
    ++i;
  } while (u > q[i]);
    
  return a + umin * q[0];
}

size_t ext_rng_drawFromDiscreteDistribution(ext_rng* generator, const double* probabilities, size_t length)
{
  if (length == 0) return EXT_DISCRETE_DRAW_FAILURE;
  
  double u = ext_rng_simulateContinuousUniform(generator);
  
  size_t result = 0;
  double sum = probabilities[0];
  
  while (sum < u && result < length - 1) {
    sum += probabilities[++result];
  }
  
  if (result == length - 1 && sum < u) return EXT_DISCRETE_DRAW_FAILURE;
  return result;
}

int64_t ext_rng_simulateIntegerUniformInRange(ext_rng* generator, int64_t min_inclusive, int64_t max_exclusive)
{
  double range = fabs((double) (max_exclusive - min_inclusive));
  int64_t actualMin = (min_inclusive < max_exclusive ? min_inclusive : max_exclusive);
  
  double u = ext_rng_simulateContinuousUniform(generator);
  
  return actualMin + (int64_t) (u * range);
}

uint64_t ext_rng_simulateUnsignedIntegerUniformInRange(ext_rng* generator, uint64_t min_inclusive, uint64_t max_exclusive)
{
  uint64_t actualMin, actualMax;
  if (min_inclusive < max_exclusive) {
    actualMin = min_inclusive;
    actualMax = max_exclusive;
  } else {
    actualMin = max_exclusive;
    actualMax = min_inclusive;
  }
  
  double range = (double) (actualMax - actualMin);
  
  double u = ext_rng_simulateContinuousUniform(generator);
  
  return actualMin + (uint64_t) (u * range);
}

#define MAX_ITER 1000
double ext_rng_simulateLowerTruncatedStandardNormal(ext_rng* generator, double lowerBound)
{
  double x;
  int iter = 0;
  if (lowerBound < 0.0) {
    x = ext_rng_simulateStandardNormal(generator);
    while (x < lowerBound && iter++ < MAX_ITER) x = ext_rng_simulateStandardNormal(generator);
    if (iter == MAX_ITER && x >= lowerBound) return nan("");
  } else {
    double a = 0.5 * (lowerBound + sqrt(lowerBound * lowerBound + 4.0));
    double u, r;
    
    do {
      x = ext_rng_simulateExponential(generator, 1.0 / a) + lowerBound;
      u = ext_rng_simulateContinuousUniform(generator);
      double diff = x - a;
      r = exp(-0.5 * diff * diff);
    } while (u > r && iter++ < MAX_ITER);
    if (iter == MAX_ITER && u <= r) return nan("");
  }
  return x;
}

double ext_rng_simulateLowerTruncatedNormalScale1(ext_rng* generator, double mean, double bound) {
  return mean + ext_rng_simulateLowerTruncatedStandardNormal(generator, bound - mean);
}

double ext_rng_simulateUpperTruncatedNormalScale1(ext_rng* generator, double mean, double bound) {
  return mean - ext_rng_simulateLowerTruncatedStandardNormal(generator, mean - bound);
}
