/* The standalone-math target's driver: the support libraries (external, misc)
 * and the copy of R's math compiled with DBARTS_HOST_STANDALONE, linked with no
 * libR. It calls each of the six math names through the support library's own
 * header, draws from the generators that reach them (the inversion normal, the
 * truncated normals), prints through the standalone output functions, and
 * checks values that hold exactly in any correct build. Exit code 0 on
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
