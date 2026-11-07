// Times the certified accept-or-keep Polya-Gamma step on one row
// (src/bartcore/pgdensity.hpp) against today's exact sum of b PG(1, psi)
// draws (ext_rng_simulatePolyaGamma), in one process on one core.
//
// Per |psi| in {0, 8, 50, 150} and b in {0.2, 1, 2.5, 16, 64, 1e3, 1e4, 1e5}:
// - c_R, the four-evaluation cost of a step from a fixed point x (the K + 1
//   gammas, the uniform, the saddle solves and the four densities through the
//   levels the decision needs), the maximum over x = E + z sd, z in
//   {-3, 0, 3} (E / 20 where that is not positive);
// - the three-evaluation figure, c_R less the cost of log f at x, as a
//   cache of log f at the current point would save it (the cache's own
//   double-double update is not built and not charged);
// - the proposal's build (K by the rule, the rates, the tail sums), reported
//   beside c_R and not in it;
// - c_S, the sum's cost at whole b.
// Then b*(|psi|): the smallest whole b with c_R <= c_S / 2, log-interpolated
// between whole-b grid points and rounded up (infinite if never reached), for
// the four- and the three-evaluation costs, and the rule's max(b*, 64).
//
// Usage: ./pg_step [seconds per cell, default 0.2]
// Build: make pg_step (after an in-place package build).

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include <external/random.h>

#include <bartcore/pgdensity.hpp>

namespace pg = bartcore::pgdensity;

namespace {

double now() {
  return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

double pgMean(double b, double psi) {
  return psi == 0.0 ? b / 4.0 : b * std::tanh(psi / 2.0) / (2.0 * psi);
}

double pgSd(double b, double psi) {
  if (psi == 0.0) return std::sqrt(b / 24.0);
  double c = std::cosh(psi / 2.0);
  return std::sqrt(b * (std::sinh(psi) - psi) / (4.0 * psi * psi * psi * c * c));
}

volatile double sink = 0.0;

/// Mean seconds per call of fn, repeated until `budget` seconds pass (at
/// least 3 calls), the best of three such runs.
template <class Fn>
double timeIt(Fn fn, double budget) {
  double best = 1e300;
  for (int rep = 0; rep < 3; ++rep) {
    int n = 0;
    double t0 = now(), t = t0;
    while (n < 3 || t - t0 < budget / 3.0) {
      fn();
      ++n;
      if ((n & 7) == 0 || n < 8) t = now();
    }
    t = now();
    best = std::min(best, (t - t0) / n);
  }
  return best;
}

/// log density at x through the levels until finite, as the step evaluates
/// it when that density's enclosure alone is wanted.
template <class F>
void evaluateOnce(F& fam, double x) {
  fam.prepare(x);
  pg::Enclosure e;
  for (int level = 0; level < pg::kNumLevels && !e.finite(); ++level) e = pg::logDensity(fam, level);
  sink = sink + e.mid;
}

/// b*: the smallest whole b where cR <= cS / 2, log-interpolated in b
/// between whole-b grid points.
double bStar(const std::vector<double>& bs, const std::vector<double>& cR, const std::vector<double>& cS) {
  double prevB = 0.0, prevG = 0.0;
  bool have = false;
  for (std::size_t i = 0; i < bs.size(); ++i) {
    if (bs[i] != std::floor(bs[i])) continue;
    double g = std::log(cR[i]) - std::log(0.5 * cS[i]);
    if (g <= 0.0) {
      if (!have) return bs[i];
      double lb = std::log(prevB) + (std::log(bs[i]) - std::log(prevB)) * prevG / (prevG - g);
      return std::ceil(std::exp(lb));
    }
    prevB = bs[i];
    prevG = g;
    have = true;
  }
  return INFINITY;
}

}  // namespace

int main(int argc, char** argv) {
  const double budget = argc > 1 ? std::atof(argv[1]) : 0.2;
  ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, nullptr);
  ext_rng_setSeed(rng, 20261107u);
  pg::ExtRngWordSource words{rng};
  pg::StepWorkspace ws;

  const double psis[] = {0.0, 8.0, 50.0, 150.0};
  const std::vector<double> bs = {0.2, 1.0, 2.5, 16.0, 64.0, 1e3, 1e4, 1e5};
  const double zs[] = {-3.0, 0.0, 3.0};

  printf("# seconds per cell %g; times in us per row step\n", budget);
  printf("# psi b K | c_R (z at max) | 3-eval | build | c_S | c_R / c_S | close calls\n");
  for (double psi : psis) {
    std::vector<double> c4(bs.size()), c3(bs.size()), cs(bs.size(), INFINITY);
    for (std::size_t i = 0; i < bs.size(); ++i) {
      const double b = bs[i];
      pg::Proposal prop(b, psi);
      double worst4 = 0.0, worst3 = 0.0, zAt = 0.0;
      long closeCalls = 0;
      for (double z : zs) {
        double x = pgMean(b, psi) + z * pgSd(b, psi);
        if (!(x > 0.0)) x = pgMean(b, psi) / 20.0;
        double t4 = timeIt(
          [&] {
            pg::StepResult r = pg::stepRow(prop, x, rng, words, ws);
            closeCalls += r.closeCall;
            sink = sink + r.omega;
          },
          budget);
        pg::TargetDensity fx(b, psi);
        double tfx = timeIt([&] { evaluateOnce(fx, x); }, budget);
        if (t4 > worst4) {
          worst4 = t4;
          zAt = z;
        }
        worst3 = std::max(worst3, t4 - tfx);
      }
      double tb = timeIt(
        [&] {
          pg::Proposal p(b, psi);
          sink = sink + p.beta();
        },
        budget);
      c4[i] = worst4;
      c3[i] = worst3;
      if (b == std::floor(b)) {
        const long n = static_cast<long>(b);
        cs[i] = timeIt(
          [&] {
            double s = 0.0;
            for (long k = 0; k < n; ++k) s += ext_rng_simulatePolyaGamma(rng, psi);
            sink = sink + s;
          },
          budget);
      }
      printf("%g %g %d | %.2f (%g) | %.2f | %.2f | ", psi, b, prop.numTerms(), 1e6 * worst4, zAt,
             1e6 * worst3, 1e6 * tb);
      if (std::isfinite(cs[i])) printf("%.2f | %.3f", 1e6 * cs[i], worst4 / cs[i]);
      else printf("- | -");
      printf(" | %ld%s\n", closeCalls, b != std::floor(b) && 1e6 * worst4 > 70.0 ? " | above 70 us" : "");
      fflush(stdout);
    }
    double s4 = bStar(bs, c4, cs), s3 = bStar(bs, c3, cs);
    printf("# psi %g: timed b* %g (four evaluations), %g (three); the rule's b* %g; three-evaluation "
           "b* %.1f%% below\n",
           psi, s4, s3, std::max(s4, 64.0), std::isfinite(s4) ? 100.0 * (1.0 - s3 / s4) : 0.0);
  }
  ext_rng_destroy(rng);
  return 0;
}
