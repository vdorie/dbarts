/* Time today's exact sum: ext_rng_simulatePolyaGamma (one PG(1, psi) draw) from the package's sources,
 * compiled standalone. A PG(b, psi) draw by the sum is b such draws. Usage: ./pgbench */
#include <stdio.h>
#include <time.h>
#include <external/random.h>
static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + 1e-9 * t.tv_nsec; }
int main(void) {
  ext_rng* rng = NULL;
  if (ext_rng_createAndSeed(&rng, EXT_RNG_ALGORITHM_MERSENNE_TWISTER, EXT_RNG_STANDARD_NORMAL_INVERSION) != 0) return 1;
  const double psis[] = {0.0, 1.0, 8.0, 50.0, 150.0};
  for (int j = 0; j < 5; ++j) {
    const int n = 5000000; double s = 0.0, t = now();
    for (int i = 0; i < n; ++i) s += ext_rng_simulatePolyaGamma(rng, psis[j]);
    printf("psi %g: %.1f ns per PG(1) draw (mean %.4g)\n", psis[j], 1e9 * (now() - t) / n, s / n);
  }
  ext_rng_destroy(rng);
  return 0;
}
