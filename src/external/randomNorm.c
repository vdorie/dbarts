/*
 * Much of this is based off of R's RNG.c, available at r-project.org and published under
 * the GNU General Public License (http://www.r-project.org/Licenses/).
 */

#include <external/random.h>

#include <math.h> // log, sqrt, exp, M_PI, sin, cos
#include <float.h> // DBL_MIN
#include <stdbool.h> // true

#include <external/io.h>
#include <external/stats.h> // qnorm

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


/*
 *  REFERENCE
 *
 *    Ahrens, J.H. and Dieter, U.
 *    Extensions of Forsythe's method for random sampling from
 *    the normal distribution.
 *    Math. Comput. 27, 927-937.
 *
 *    The definitions of the constants a[k], d[k], t[k] and
 *    h[k] are according to the abovementioned article
 */
const static double a[32] = {
  0.0000000, 0.03917609, 0.07841241, 0.1177699,
  0.1573107, 0.19709910, 0.23720210, 0.2776904,
  0.3186394, 0.36012990, 0.40225010, 0.4450965,
  0.4887764, 0.53340970, 0.57913220, 0.6260990,
  0.6744898, 0.72451440, 0.77642180, 0.8305109,
  0.8871466, 0.94678180, 1.00999000, 1.0775160,
  1.1503490, 1.22985900, 1.31801100, 1.4177970,
  1.5341210, 1.67594000, 1.86273200, 2.1538750
};

const static double d[31] = {
  0.0000000, 0.0000000, 0.0000000, 0.0000000,
  0.0000000, 0.2636843, 0.2425085, 0.2255674,
  0.2116342, 0.1999243, 0.1899108, 0.1812252,
  0.1736014, 0.1668419, 0.1607967, 0.1553497,
  0.1504094, 0.1459026, 0.1417700, 0.1379632,
  0.1344418, 0.1311722, 0.1281260, 0.1252791,
  0.1226109, 0.1201036, 0.1177417, 0.1155119,
  0.1134023, 0.1114027, 0.1095039
};

const static double t[31] = {
  7.673828e-4, 0.002306870, 0.003860618, 0.005438454,
  0.007050699, 0.008708396, 0.010423570, 0.012209530,
  0.014081250, 0.016055790, 0.018152900, 0.020395730,
  0.022811770, 0.025434070, 0.028302960, 0.031468220,
  0.034992330, 0.038954830, 0.043458780, 0.048640350,
  0.054683340, 0.061842220, 0.070479830, 0.081131950,
  0.094624440, 0.112300100, 0.136498000, 0.171688600,
  0.227624100, 0.330498000, 0.584703100
};

const static double h[31] = {
  0.03920617, 0.03932705, 0.03950999, 0.03975703,
  0.04007093, 0.04045533, 0.04091481, 0.04145507,
  0.04208311, 0.04280748, 0.04363863, 0.04458932,
  0.04567523, 0.04691571, 0.04833487, 0.04996298,
  0.05183859, 0.05401138, 0.05654656, 0.05953130,
  0.06308489, 0.06737503, 0.07264544, 0.07926471,
  0.08781922, 0.09930398, 0.11555990, 0.14043440,
  0.18361420, 0.27900160, 0.70104740
};

double ext_rng_simulateStandardNormal(ext_rng* generator)
{
  switch (generator->standardNormalAlgorithm) {
    case EXT_RNG_STANDARD_NORMAL_AHRENS_DIETER:
    {
      double u1, u2;
      double s, w;
      double aa, tt;
      int_least32_t i;
      
      u1 = ext_rng_simulateContinuousUniform(generator);
      s = 0.0;
      s = u1 > 0.5 ? 1.0 : 0.0;
      u1 = u1 + u1 - s;
      u1 *= 32.0;
      i = (int_least32_t) u1;
      if (i == 32) i = 31;
      if (i != 0) {
        u2 = u1 - (double) i;
        aa = a[i - 1];
        while (u2 <= t[i - 1]) {
          u1 = ext_rng_simulateContinuousUniform(generator);
          w = u1 * (a[i] - aa);
          tt = (w * 0.5 + aa) * w;
          do {
            if (u2 > tt) goto deliver;
            u1 = ext_rng_simulateContinuousUniform(generator);
            if (u2 < u1) break;
            tt = u1;
            u2 = ext_rng_simulateContinuousUniform(generator);
          } while (true);
          u2 = ext_rng_simulateContinuousUniform(generator);
        }
        w = (u2 - t[i - 1]) * h[i - 1];
      } else {
        i = 6;
        aa = a[31];
        do {
          u1 = u1 + u1;
          if (u1 >= 1.0) break;
          aa = aa + d[i - 1];
          i = i + 1;
        } while (true);
        u1 = u1 - 1.0;
        do {
          w = u1 * d[i - 1];
          tt = (w * 0.5 + aa) * w;
          do {
            u2 = ext_rng_simulateContinuousUniform(generator);
            if (u2 > tt) goto jump;
            u1 = ext_rng_simulateContinuousUniform(generator);
            if (u2 < u1) break;
            tt = u1;
          } while (true);
          u1 = ext_rng_simulateContinuousUniform(generator);
        } while (true);
jump:   /* */ ;
    	}

deliver:
      return (s == 1.0) ? -(aa + w) : (aa + w);
    }
    break;
    
    case EXT_RNG_STANDARD_NORMAL_BUGGY_KINDERMAN_RAMAGE:
    /* see Reference above */
    /* note: this has problems, but is retained for
     * reproducibility of older codes, with the same
     * numeric code */
#define C1		0.398942280401433
#define C2		0.180025191068563
#define g(x)  (C1 * exp(-x * x / 2.0) - C2 * (A - x))
#define A     2.216035867166471
    {
      double u1, u2, u3;
      double tt;
      
      u1 = ext_rng_simulateContinuousUniform(generator);
      if (u1 < 0.884070402298758) {
        u2 = ext_rng_simulateContinuousUniform(generator);
        return A * (1.13113163544180 * u1 + u2 - 1.0);
    	}
      
      if (u1 >= 0.973310954173898) { /* tail: */
        do {
          u2 = ext_rng_simulateContinuousUniform(generator);
          u3 = ext_rng_simulateContinuousUniform(generator);
          tt = (A * A - 2.0 * log(u3));
          if (u2 * u2 < (A * A) / tt)
            return (u1 < 0.986655477086949) ? sqrt(tt) : -sqrt(tt);
        } while (true);
      }
      
      if (u1 >= 0.958720824790463) { /* region3: */
        do {
          u2 = ext_rng_simulateContinuousUniform(generator);
          u3 = ext_rng_simulateContinuousUniform(generator);
    		  tt = A - 0.630834801921960 * fmin(u2, u3);
          if (fmax(u2, u3) <= 0.755591531667601)
            return (u2 < u3) ? tt : -tt;
          if (0.034240503750111 * fabs(u2 - u3) <= g(tt))
            return (u2<u3) ? tt : -tt;
        } while (true);
      }
      
      if (u1 >= 0.911312780288703) { /* region2: */
        do {
          u2 = ext_rng_simulateContinuousUniform(generator);
          u3 = ext_rng_simulateContinuousUniform(generator);
          tt = 0.479727404222441 + 1.105473661022070 * fmin(u2, u3);
          if (fmax2(u2, u3) <= 0.872834976671790)
            return (u2 < u3) ? tt : -tt;
          if (0.049264496373128 * fabs(u2 - u3) <= g(tt))
            return (u2 < u3) ? tt : -tt;
        } while (true);
      }
      
      /* ELSE	 region1: */
      do {
        u2 = ext_rng_simulateContinuousUniform(generator);
        u3 = ext_rng_simulateContinuousUniform(generator);
        tt = 0.479727404222441 - 0.595507138015940 * fmin(u2, u3);
        if (fmax(u2, u3) <= 0.805577924423817)
          return (u2 < u3) ? tt : -tt;
      } while (true);
    }
    break;
    
    
    case EXT_RNG_STANDARD_NORMAL_KINDERMAN_RAMAGE: // see Reference above
    // corrected version from Josef Leydold
    {
      double u1, u2, u3;
      double tt;
      
      u1 = ext_rng_simulateContinuousUniform(generator);
      if (u1 < 0.884070402298758) {
        u2 = ext_rng_simulateContinuousUniform(generator);
        return A * (1.131131635444180 * u1 + u2 - 1.0);
      }

      if (u1 >= 0.973310954173898) { /* tail: */
        do {
          u2 = ext_rng_simulateContinuousUniform(generator);
          u3 = ext_rng_simulateContinuousUniform(generator);
          tt = (A * A - 2 * log(u3));
          if (u2 * u2 < (A * A) / tt)
            return (u1 < 0.986655477086949) ? sqrt(tt) : -sqrt(tt);
        } while (true);
      }
      
      if (u1 >= 0.958720824790463) { /* region3: */
        do {
          u2 = ext_rng_simulateContinuousUniform(generator);
          u3 = ext_rng_simulateContinuousUniform(generator);
          tt = A - 0.630834801921960 * fmin(u2, u3);
          if (fmax(u2, u3) <= 0.755591531667601)
            return (u2 < u3) ? tt : -tt;
          if (0.034240503750111 * fabs(u2 - u3) <= g(tt))
            return (u2 < u3) ? tt : -tt;
        } while (true);
      }

      if (u1 >= 0.911312780288703) { /* region2: */
        do {
          u2 = ext_rng_simulateContinuousUniform(generator);
          u3 = ext_rng_simulateContinuousUniform(generator);
          tt = 0.479727404222441+1.105473661022070*fmin2(u2,u3);
          if (fmax(u2, u3) <= 0.872834976671790)
            return (u2 < u3) ? tt : -tt;
          if (0.049264496373128 * fabs(u2 - u3) <= g(tt))
            return (u2 < u3) ? tt : -tt;
        } while (true);
      }
      
      /* ELSE	 region1: */
      do {
        u2 = ext_rng_simulateContinuousUniform(generator);
        u3 = ext_rng_simulateContinuousUniform(generator);
        tt = 0.479727404222441 - 0.595507138015940 * fmin(u2, u3);
        if (tt < 0.0) continue;
        if (fmax(u2, u3) <= 0.805577924423817)
          return (u2 < u3) ? tt : -tt;
        if (0.053377549506886 * fabs(u2 - u3) <= g(tt))
          return (u2 < u3) ? tt : -tt;
      } while (true);
    }
    break;
#undef A
#undef g
#undef C2
#undef C1
    
    case EXT_RNG_STANDARD_NORMAL_BOX_MULLER:
    {
	    if (generator->normalState.nextNormal != 0.0) { /* An exact test is intentional */
        double s = generator->normalState.nextNormal;
        generator->normalState.nextNormal = 0.0;
        return s;
      } else {
	      double theta = 2.0 * M_PI * ext_rng_simulateContinuousUniform(generator);
        double R = sqrt(-2.0 * log(ext_rng_simulateContinuousUniform(generator))) + 10.0 * DBL_MIN; /* ensure non-zero */
        generator->normalState.nextNormal = R * sin(theta);
	      return R * cos(theta);
      }
    }
    break;
    
    case EXT_RNG_STANDARD_NORMAL_INVERSION:
    {
#define BIG 134217728 /* 2^27 */
      /* unif_rand() alone is not of high enough precision */
      double u1 = ext_rng_simulateContinuousUniform(generator);
      u1 = (double) ((int_least32_t) (BIG * u1)) + ext_rng_simulateContinuousUniform(generator);
      return ext_quantileOfNormal(u1 / BIG, 0.0, 1.0);
    }
    break;
#undef BIG
    case EXT_RNG_STANDARD_NORMAL_USER_NORM:
    return (generator->normalState.simulateNormal.state == NULL ? generator->normalState.simulateNormal.f.stateless() :
              generator->normalState.simulateNormal.f.stateful(generator->normalState.simulateNormal.state));
    break;
    default:
    ext_throwError("unsupported standard normal generator type");
    break;
  }
  
  return NAN;
}
