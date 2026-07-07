// Sweep-level profile of the linear-leaf hot paths: what fraction of a
// linear-leaves MCMC sweep goes to the crossproduct accumulation and to
// chol(V), the two candidate reuse optimizations of
// docs/plans/linear-leaf-reuse.md.
//
// The denominator is the wall time of uninstrumented sweeps. The leaf-method
// time is measured directly - a shadowing leaf clocks the real leaf bodies,
// delegating all math and rng to the base, so no cost model enters the headline
// fractions. Standalone kernels attribute that time: choleskyDecompose prices
// chol(V), and a crossproduct pass timed with vs. without its U'WU inner loop
// prices the part a persistent stats cache could actually skip (U'Wz and z'Wz
// still need the observation pass every sweep, since the residual moves).
//
// Build the package first (R CMD INSTALL . leaves src/*.a), then
// make linear_leaf && ./linear_leaf.

#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <algorithm>
#include <vector>

#include <misc/simd.h>
#include <external/random.h>

#include <bartcore/bartcore.hpp>

using namespace bartcore;
using std::size_t;

namespace {

constexpr size_t kMaxParams = LinearGaussianLeaf::maxNumCovariates + 1;

uint64_t nowNs() {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return static_cast<uint64_t>(ts.tv_sec) * 1000000000ull +
         static_cast<uint64_t>(ts.tv_nsec);
}

uint64_t rngState = 0x9E3779B97F4A7C15ull;
double nextUniform() {
  rngState ^= rngState << 13;
  rngState ^= rngState >> 7;
  rngState ^= rngState << 17;
  return static_cast<double>(rngState >> 11) * 0x1.0p-53;
}

double median(std::vector<double>& v) {
  std::sort(v.begin(), v.end());
  size_t m = v.size() / 2;
  return v.size() % 2 ? v[m] : 0.5 * (v[m - 1] + v[m]);
}

// Per-sweep counts and in-method wall time, harvested by shadowing the two hot
// leaf methods; the base does the arithmetic, so evolution and rng are identical.
struct Counters {
  uint64_t scoreCalls = 0, scoreObs = 0, scoreNs = 0;
  uint64_t drawCalls = 0, drawObs = 0, drawNs = 0;
};
Counters g;

struct TimedLinearLeaf : LinearGaussianLeaf {
  double logIntegratedLikelihoodForNode(const Tree& tree, const double* y,
                                        const double* w, double k, double rv,
                                        int32_t node) const {
    size_t n = tree.at(node).numObservations();
    uint64_t t0 = nowNs();
    double r = LinearGaussianLeaf::logIntegratedLikelihoodForNode(tree, y, w, k,
                                                                  rv, node);
    if (n > 0) {
      g.scoreNs += nowNs() - t0;
      ++g.scoreCalls;
      g.scoreObs += n;
    }
    return r;
  }
  void drawFromPosteriorForNode(ext_rng* rng, const Tree& tree, const double* y,
                                const double* w, double k, double rv,
                                int32_t node, double* out) const {
    size_t n = tree.at(node).numObservations();
    uint64_t t0 = nowNs();
    LinearGaussianLeaf::drawFromPosteriorForNode(rng, tree, y, w, k, rv, node,
                                                 out);
    if (n > 0) {
      g.drawNs += nowNs() - t0;
      ++g.drawCalls;
      g.drawObs += n;
    }
  }
};

// One node-statistics pass, copied from
// LinearGaussianLeaf::accumulateNodeStatistics. withCrossproduct=false drops
// the U'WU inner loop, leaving the U'Wz/z'Wz that a cache cannot avoid; the two
// timings' difference is the cacheable work. u is column-major n x q.
double timeAccumulate(size_t q, const double* u, size_t n, const size_t* seg,
                      size_t segLen, const double* y, size_t reps,
                      bool withCrossproduct) {
  size_t p = q + 1;
  double crossproduct[kMaxParams * kMaxParams], projection[kMaxParams];
  double row[kMaxParams];
  volatile double sink = 0.0;
  uint64_t t0 = nowNs();
  for (size_t r = 0; r < reps; ++r) {
    for (size_t a = 0; a < p * p; ++a) crossproduct[a] = 0.0;
    for (size_t a = 0; a < p; ++a) projection[a] = 0.0;
    double responseSumOfSquares = 0.0;
    row[0] = 1.0;
    for (size_t m = 0; m < segLen; ++m) {
      size_t i = seg[m];
      double z = y[i];
      for (size_t j = 0; j < q; ++j) row[j + 1] = u[i + j * n];
      for (size_t a = 0; a < p; ++a) {
        double scaled = row[a];
        projection[a] += scaled * z;
        if (withCrossproduct)
          for (size_t b = a; b < p; ++b) crossproduct[a * p + b] += scaled * row[b];
      }
      responseSumOfSquares += z * z;
    }
    sink += crossproduct[0] + projection[0] + responseSumOfSquares;
  }
  (void) sink;
  return static_cast<double>(nowNs() - t0) / static_cast<double>(reps);
}

// In-place lower Cholesky, copied from LinearGaussianLeaf::choleskyDecompose;
// priced net of the per-rep SPD refresh (chol is destructive). Ns per call.
void choleskyInPlace(double* m, size_t p) {
  for (size_t j = 0; j < p; ++j) {
    double diagonal = m[j * p + j];
    for (size_t a = 0; a < j; ++a) diagonal -= m[j * p + a] * m[j * p + a];
    diagonal = std::sqrt(diagonal);
    m[j * p + j] = diagonal;
    for (size_t i = j + 1; i < p; ++i) {
      double value = m[i * p + j];
      for (size_t a = 0; a < j; ++a) value -= m[i * p + a] * m[j * p + a];
      m[i * p + j] = value / diagonal;
    }
  }
}

double timeChol(size_t p, size_t reps) {
  double base[kMaxParams * kMaxParams], m[kMaxParams * kMaxParams];
  for (size_t a = 0; a < p; ++a)
    for (size_t b = 0; b < p; ++b)
      base[a * p + b] = a == b ? static_cast<double>(p + 1) : 0.3;
  volatile double sink = 0.0;
  uint64_t t0 = nowNs();
  for (size_t r = 0; r < reps; ++r) {
    for (size_t k = 0; k < p * p; ++k) m[k] = base[k];
    choleskyInPlace(m, p);
    sink += m[0];
  }
  uint64_t both = nowNs() - t0;
  t0 = nowNs();
  for (size_t r = 0; r < reps; ++r) {
    for (size_t k = 0; k < p * p; ++k) m[k] = base[k];
    sink += m[0];
  }
  double copyOnly = static_cast<double>(nowNs() - t0);
  (void) sink;
  double d = static_cast<double>(both) - copyOnly;
  return (d > 0.0 ? d : 0.0) / static_cast<double>(reps);
}

struct Config {
  size_t n;
  size_t q;
  int signal;  // 0: linear (shallow trees); 1: step (deep trees)
};

void measure(const Config& cfg, size_t numTrees, size_t warmup, size_t sweeps,
             size_t reps) {
  size_t n = cfg.n, q = cfg.q, p = q;

  std::vector<double> x(n * p), y(n);
  for (size_t j = 0; j < p; ++j)
    for (size_t i = 0; i < n; ++i) x[i + j * n] = nextUniform();
  for (size_t i = 0; i < n; ++i) {
    double f = 0.0;
    if (cfg.signal == 0) {
      f = 2.0 * x[i] + 1.5 * x[i + n];
      if (p > 2) f += (x[i + 2 * n] > 0.5 ? 1.0 : -1.0);
    } else {
      for (size_t j = 0; j < p; ++j) f += (x[i + j * n] > 0.5 ? 0.8 : -0.8);
    }
    y[i] = f + 0.5 * (nextUniform() - 0.5);
  }
  double mean = 0.0;
  for (double v : y) mean += v;
  mean /= static_cast<double>(n);
  double sd = 0.0;
  for (double v : y) sd += (v - mean) * (v - mean);
  sd = std::sqrt(sd / static_cast<double>(n - 1));

  std::vector<size_t> columns(q);
  for (size_t j = 0; j < q; ++j) columns[j] = j;
  SamplerOptions options;
  options.numTrees = numTrees;
  options.leafCovariateColumns = columns.data();
  options.numLeafCovariates = q;
  Results empty;

  // Denominator: uninstrumented sweeps.
  ext_rng* rngT = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, nullptr);
  ext_rng_setSeed(rngT, 20250706u);
  SamplerFacade<LinearGaussianLeaf> plain(x.data(), y.data(), n, p, nullptr,
                                          nullptr, ResponseFamily::gaussian, sd,
                                          3.0, 0.9, options, &rngT);
  plain.run(warmup, 0, empty);
  std::vector<double> sweepNs;
  for (size_t r = 0; r < reps; ++r) {
    uint64_t t0 = nowNs();
    plain.run(sweeps, 0, empty);
    sweepNs.push_back(static_cast<double>(nowNs() - t0) /
                      static_cast<double>(sweeps));
  }
  double tSweep = median(sweepNs);
  ext_rng_destroy(rngT);

  // Leaf-method time and counts, same seed and data.
  ext_rng* rngC = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, nullptr);
  ext_rng_setSeed(rngC, 20250706u);
  SamplerFacade<TimedLinearLeaf> timed(x.data(), y.data(), n, p, nullptr,
                                       nullptr, ResponseFamily::gaussian, sd,
                                       3.0, 0.9, options, &rngC);
  timed.run(warmup, 0, empty);
  g = Counters{};
  timed.run(sweeps, 0, empty);
  ext_rng_destroy(rngC);
  double s = static_cast<double>(sweeps);
  double scoreCalls = g.scoreCalls / s, scoreObs = g.scoreObs / s;
  double drawCalls = g.drawCalls / s, drawObs = g.drawObs / s;
  double scoreNs = g.scoreNs / s, drawNs = g.drawNs / s;
  double drawMean = drawCalls > 0 ? drawObs / drawCalls : 1.0;

  // Standalone kernel costs at the observed leaf scale.
  std::vector<size_t> perm(n);
  for (size_t i = 0; i < n; ++i) perm[i] = i;
  for (size_t i = n - 1; i > 0; --i)
    std::swap(perm[i], perm[static_cast<size_t>(nextUniform() *
                                                static_cast<double>(i + 1))]);
  size_t seg = std::max<size_t>(64, static_cast<size_t>(drawMean));
  if (seg > n) seg = n;
  size_t repsA = std::max<size_t>(20, 30000000 / seg);
  double full = timeAccumulate(q, x.data(), n, perm.data(), seg, y.data(),
                               repsA, true);
  double reduced = timeAccumulate(q, x.data(), n, perm.data(), seg, y.data(),
                                  repsA, false);
  double cacheableFrac = full > 0.0 ? (full - reduced) / full : 0.0;
  double tChol = timeChol(q + 1, 2000000);

  // Crossproduct = leaf-method time net of chol (draw also carries small
  // solve+rng terms); only the U'WU inner loop of it is cacheable.
  double cholScore = scoreCalls * tChol, cholDraw = drawCalls * tChol;
  double xpScore = scoreNs - cholScore, xpDraw = drawNs - cholDraw;
  double cacheable = cacheableFrac * (xpScore + xpDraw);

  auto pct = [&](double x) { return 100.0 * x / tSweep; };
  const char* sig = cfg.signal == 0 ? "linear" : "step  ";
  printf("n=%-7zu q=%zu %s | sweep=%.0f us leaves/sweep=%.1f obs/leaf=%.0f\n", n,
         q, sig, tSweep / 1000.0, drawCalls, drawMean);
  printf("  crossproduct  score %6.2f%%  draw %6.2f%%  total %6.2f%%\n",
         pct(xpScore), pct(xpDraw), pct(xpScore + xpDraw));
  printf("  chol(V)       score %6.3f%%  draw %6.3f%%  total %6.3f%%\n",
         pct(cholScore), pct(cholDraw), pct(cholScore + cholDraw));
  printf("  ceiling(a+b: score-xprod + draw-chol) %6.2f%%   "
         "cacheable U'WU share of xprod %5.1f%% -> realistic xprod ceiling "
         "%6.2f%%\n\n",
         pct(xpScore + cholDraw), 100.0 * cacheableFrac, pct(cacheable));
  fflush(stdout);
}

}  // namespace

int main(int argc, char** argv) {
  misc_simd_init();
  setvbuf(stdout, nullptr, _IOLBF, 0);
  size_t numTrees = 50, warmup = 30, sweeps = 40, reps = 5;
  if (argc > 1) sweeps = static_cast<size_t>(std::atol(argv[1]));

  printf("# linear-leaf sweep profile: numTrees=%zu warmup=%zu sweeps=%zu "
         "reps=%zu, single chain/thread\n",
         numTrees, warmup, sweeps, reps);
  printf("# leaf basis capped at LinearGaussianLeaf::maxNumCovariates=%zu; "
         "step signal drives deeper trees\n\n",
         LinearGaussianLeaf::maxNumCovariates);

  Config configs[] = {{10000, 4, 0},  {10000, 8, 0},  {100000, 4, 0},
                      {100000, 8, 0}, {10000, 8, 1},  {100000, 8, 1}};
  for (const Config& c : configs) measure(c, numTrees, warmup, sweeps, reps);
  return 0;
}
