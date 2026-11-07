/* The standalone-math target's driver: the support libraries (external, misc)
 * and the copy of R's math compiled with DBARTS_HOST_STANDALONE, linked with no
 * libR. It calls each of the seven math names through the support library's own
 * header, draws from the generators that reach them (the inversion normal, the
 * truncated normals), prints through the standalone output functions, checks
 * values that hold exactly in any correct build, and checks that each name
 * hands every argument and flag to its renamed routine. Exit code 0 on
 * success. The copy's agreement with R's own is the libR suite's nmath grid. */

#include <external/io.h>
#include <external/random.h>
#include <external/stats.h>
#include <misc/simd.h>

#include <math.h>
#include <stdio.h>

static int failures = 0;

static void check(int condition, const char* what) {
  if (!condition) {
    ++failures;
    printf("FAIL: %s\n", what);
  }
}

int main(void) {
  misc_simd_init();

  check(ext_pnorm(0.0, 0.0, 1.0, 1, 0) == 0.5, "pnorm(0) is 1/2");
  check(ext_pnorm(1.0, 1.0, 2.0, 0, 1) == log(0.5), "pnorm's upper tail on the log scale");
  check(ext_qnorm(0.5, 0.0, 1.0, 1, 0) == 0.0, "qnorm(1/2) is 0");
  check(ext_qnorm(0.975, 0.0, 1.0, 1, 0) > 1.9599 && ext_qnorm(0.975, 0.0, 1.0, 1, 0) < 1.9600,
        "qnorm(0.975)");
  check(ext_dnorm(0.0, 0.0, 1.0, 0) == 0.398942280401432677939946059934, "dnorm(0) is 1 / sqrt(2 pi)");
  check(fabs(ext_dt(0.0, 1.0, 0) * M_PI - 1.0) < 1e-15, "dt(0, 1) is 1 / pi");
  check(ext_qt(0.5, 3.0, 1, 0) == 0.0, "qt(1/2) is 0");
  check(fabs(ext_qchisq(0.5, 2.0, 1, 0) - 2.0 * log(2.0)) < 1e-14, "qchisq(1/2, 2) is 2 log 2");
  check(fabs(ext_quantileOfStudentT(0.75, 1.0) - 1.0) < 1e-15, "the t prior's quartile at 1 df");

  /* Each name passes every argument and flag through to its routine: at a
     point where no tail or scale is symmetric, a dropped or fixed flag, or two
     swapped arguments, gives a different value. */
  for (int lowerTail = 0; lowerTail <= 1; ++lowerTail) {
    for (int logP = 0; logP <= 1; ++logP) {
      double p = logP ? log(0.3) : 0.3;
      char what[64];
      snprintf(what, sizeof(what), "ext_pnorm at lower.tail %d, log.p %d", lowerTail, logP);
      check(ext_pnorm(0.7, 0.2, 1.3, lowerTail, logP) ==
              dbarts_nmath_pnorm5(0.7, 0.2, 1.3, lowerTail, logP), what);
      snprintf(what, sizeof(what), "ext_qnorm at lower.tail %d, log.p %d", lowerTail, logP);
      check(ext_qnorm(p, 0.2, 1.3, lowerTail, logP) ==
              dbarts_nmath_qnorm5(p, 0.2, 1.3, lowerTail, logP), what);
      snprintf(what, sizeof(what), "ext_qt at lower.tail %d, log.p %d", lowerTail, logP);
      check(ext_qt(p, 3.5, lowerTail, logP) == dbarts_nmath_qt(p, 3.5, lowerTail, logP), what);
      snprintf(what, sizeof(what), "ext_qchisq at lower.tail %d, log.p %d", lowerTail, logP);
      check(ext_qchisq(p, 3.5, lowerTail, logP) == dbarts_nmath_qchisq(p, 3.5, lowerTail, logP),
            what);
    }
  }
  for (int logP = 0; logP <= 1; ++logP) {
    char what[64];
    snprintf(what, sizeof(what), "ext_dnorm at log %d", logP);
    check(ext_dnorm(0.7, 0.2, 1.3, logP) == dbarts_nmath_dnorm4(0.7, 0.2, 1.3, logP), what);
    snprintf(what, sizeof(what), "ext_dt at log %d", logP);
    check(ext_dt(0.7, 3.5, logP) == dbarts_nmath_dt(0.7, 3.5, logP), what);
  }
  check(ext_quantileOfChiSquared(0.3, 3.5) == dbarts_nmath_qchisq(0.3, 3.5, 1, 0),
        "ext_quantileOfChiSquared");
  check(ext_densityOfNormal(0.7, 0.2, 1.3) == dbarts_nmath_dnorm4(0.7, 0.2, 1.3, 0),
        "ext_densityOfNormal");
  check(ext_cumulativeProbabilityOfNormal(0.7, 0.2, 1.3) ==
          dbarts_nmath_pnorm5(0.7, 0.2, 1.3, 1, 0), "ext_cumulativeProbabilityOfNormal");
  check(ext_quantileOfNormal(0.3, 0.2, 1.3) == dbarts_nmath_qnorm5(0.3, 0.2, 1.3, 1, 0),
        "ext_quantileOfNormal");
  check(ext_quantileOfStudentT(0.3, 3.5) == dbarts_nmath_qt(0.3, 3.5, 1, 0),
        "ext_quantileOfStudentT");
  /* lchoose, the binomial row's constant: exact where the routine's small
     branches answer, and its two arguments passed in order at a point where
     swapping them gives -Inf */
  check(ext_lchoose(7.0, 0.0) == 0.0 && ext_lchoose(7.0, 1.0) == log(7.0) &&
          ext_lchoose(7.0, 6.0) == log(7.0), "lchoose at k = 0, 1 and n - 1");
  check(ext_lchoose(1e6, 1000.0) == dbarts_nmath_lchoose(1e6, 1000.0) &&
          ext_lchoose(1000.0, 1e6) == -INFINITY, "ext_lchoose");

  ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, NULL);
  check(rng != NULL && ext_rng_setSeed(rng, 42) == 0, "a seeded Mersenne Twister");
  if (rng == NULL) return 1;

  double sum = 0.0, sumSq = 0.0;
  int inBounds = 1;
  const int numDraws = 20000;
  for (int i = 0; i < numDraws; ++i) {
    double z = ext_rng_simulateStandardNormal(rng);
    sum += z;
    sumSq += z * z;
    double t = ext_rng_simulateTruncatedNormalScale1(rng, 0.3, -0.5, 1.25);
    if (!(t > -0.5 && t <= 1.25)) inBounds = 0;
    double tail = ext_rng_simulateLowerTruncatedNormalScale1(rng, 0.0, 4.0);
    if (!(tail > 4.0)) inBounds = 0;
  }
  double mean = sum / numDraws, variance = sumSq / numDraws - mean * mean;
  check(fabs(mean) < 0.05 && fabs(variance - 1.0) < 0.05, "the inversion normal's moments");
  check(inBounds, "every truncated draw inside its bounds");
  ext_rng_destroy(rng);

  ext_printMessage("standalone-math: %d draws, mean %.4f, variance %.4f", numDraws, mean,
                   variance);
  if (failures > 0) {
    printf("%d failure(s)\n", failures);
    return 1;
  }
  printf("standalone-math: all checks passed\n");
  return 0;
}
