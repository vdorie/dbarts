// The certified Polya-Gamma evaluator (pgdensity.hpp) and its double-double
// functions (ddouble.hpp) against reference tables generated in mpmath and
// Arb (pgdensity_refs.inc):
//
// - the pins: at every pin and every level the enclosure of log f, or of
//   log q at the rule's K from the sampler's own rates, contains a reference
//   independent of the evaluator, and the engine's radius is reported next to
//   Arb's;
// - the floor: per (b, psi) of the f pins, the step floor the engine's
//   largest radii over z imply, at most 2.3e-11;
// - the K rule's K at every cell of the K table, and r_1 below d_1;
// - the double-double gate: each of the five functions within its stated
//   relative error bound (all at most 2^-100) over its range;
// - the reach: log f and log q decided far into both tails past the pins;
// - the one-row step: decided steps and a mean that matches PG(b, psi)'s.
//
// Draws come from the suite's own generator, so it neither shifts nor is
// shifted by any other suite's.

#include "assert.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <vector>

#include <external/random.h>

#include <bartcore/pgdensity.hpp>

namespace {

namespace pg = bartcore::pgdensity;
namespace dd = bartcore::ddouble;

constexpr double kInf = std::numeric_limits<double>::infinity();

struct FPin {
  double b, psi, z, x, refHi, refLo, arbMid, arbRad;
};
struct QPin {
  double b, psi, z;
  int K;
  double x, r1, rK, alpha, beta, refHi, refLo, arbMid, arbRad;
  bool moschopoulos;
};
struct KCell {
  double b, psi;
  int K;
  double bDBelow, bDAt;
};
struct DdCase {
  double xHi, xLo, yHi, yLo;
};

#include "pgdensity_refs.inc"

template <class T, std::size_t N>
constexpr std::size_t countOf(const T (&)[N]) {
  return N;
}

/// |ref - mid| as a double, ref a double-double: exact to rounding.
double distance(double refHi, double refLo, double mid) {
  return std::fabs((refHi - mid) + refLo);
}

/// Checks every level's enclosure (and their intersection) against the
/// reference; returns the intersection through the last level.
template <class F>
pg::Enclosure checkLevels(const F& fam, double refHi, double refLo, const char* tag, int& failed,
                          int& undecided, double& worstRatio) {
  pg::Enclosure e;
  for (int level = 0; level < pg::kNumLevels; ++level) {
    pg::Enclosure el = pg::logDensity(fam, level);
    if (!el.finite()) continue;
    double d = distance(refHi, refLo, el.mid);
    worstRatio = std::max(worstRatio, d / el.rad);
    if (!(d <= el.rad)) {
      if (failed++ < 10)
        printf("FAIL: pgdensity %s level %d: |ref - mid| %.3g > rad %.3g\n", tag, level, d, el.rad);
    }
    e = pg::intersect(e, el);
  }
  if (!e.finite()) ++undecided;
  return e;
}

void testFPins(std::vector<double>& radF) {
  int failed = 0, undecided = 0;
  double worstRatio = 0.0, worstVsArb = 0.0, maxRad = 0.0;
  const FPin* worstPin = nullptr;
  radF.assign(countOf(kFPins), kInf);
  for (std::size_t i = 0; i < countOf(kFPins); ++i) {
    const FPin& p = kFPins[i];
    pg::TargetDensity f(p.b, p.psi);
    f.prepare(p.x);
    char tag[96];
    snprintf(tag, sizeof tag, "f b=%g psi=%g z=%g", p.b, p.psi, p.z);
    pg::Enclosure e = checkLevels(f, p.refHi, p.refLo, tag, failed, undecided, worstRatio);
    // Arb's enclosure must overlap the engine's: both contain the value
    if (e.finite() && !(std::fabs(e.mid - p.arbMid) <= e.rad + p.arbRad)) {
      if (failed++ < 10) printf("FAIL: pgdensity %s: disjoint from Arb's enclosure\n", tag);
    }
    radF[i] = e.rad;
    maxRad = std::max(maxRad, e.rad);
    if (e.finite() && e.rad / p.arbRad > worstVsArb) {
      worstVsArb = e.rad / p.arbRad;
      worstPin = &p;
    }
  }
  printf("pgdensity f pins: %zu, outside %d, undecided %d, max |ref - mid| / rad %.3f, max rad %.3g, "
         "max rad / Arb's %.2f (b=%g psi=%g z=%g)\n",
         countOf(kFPins), failed, undecided, worstRatio, maxRad, worstVsArb,
         worstPin ? worstPin->b : 0.0, worstPin ? worstPin->psi : 0.0, worstPin ? worstPin->z : 0.0);
  check(failed == 0, "pgdensity: every f pin's enclosure contains its reference at every level");
  check(undecided == 0, "pgdensity: every f pin is decided by the last level");
}

void testQPins() {
  int failed = 0, undecided = 0, rates = 0;
  double worstRatio = 0.0, worstVsArb = 0.0, maxRad = 0.0, worstAB = 0.0;
  for (const QPin& p : kQPins) {
    pg::Proposal prop(p.b, p.psi);
    char tag[96];
    snprintf(tag, sizeof tag, "q b=%g psi=%g z=%g K=%d", p.b, p.psi, p.z, p.K);
    if (prop.numTerms() != p.K || prop.rate(0) != p.r1 || prop.rate(p.K - 1) != p.rK) {
      if (rates++ < 10)
        printf("FAIL: pgdensity %s: the proposal's K or rates differ from the table's (K %d, r1 %a, "
               "rK %a)\n",
               tag, prop.numTerms(), prop.rate(0), prop.rate(prop.numTerms() - 1));
      continue;
    }
    worstAB = std::max({worstAB, std::fabs(prop.alpha() / p.alpha - 1.0),
                        std::fabs(prop.beta() / p.beta - 1.0)});
    // the reference is at the table's (alpha, beta), formed from 40-digit
    // tail sums; the engine's own agree to about 1e-13
    std::vector<double> shapes(prop.shapes(), prop.shapes() + p.K + 1);
    std::vector<double> rs(prop.rates(), prop.rates() + p.K + 1);
    shapes[p.K] = p.alpha;
    rs[p.K] = p.beta;
    pg::GammaSumDensity q;
    q.setParameters(shapes.data(), rs.data(), shapes.size());
    q.prepare(p.x);
    pg::Enclosure e = checkLevels(q, p.refHi, p.refLo, tag, failed, undecided, worstRatio);
    if (e.finite() && !(std::fabs(e.mid - p.arbMid) <= e.rad + p.arbRad)) {
      if (failed++ < 10) printf("FAIL: pgdensity %s: disjoint from Arb's enclosure\n", tag);
    }
    maxRad = std::max(maxRad, e.rad);
    if (e.finite()) worstVsArb = std::max(worstVsArb, e.rad / p.arbRad);
  }
  printf("pgdensity q pins: %zu, outside %d, undecided %d, max |ref - mid| / rad %.3f, max rad %.3g, "
         "max rad / Arb's %.2f, alpha and beta within %.2g of the table's\n",
         countOf(kQPins), failed, undecided, worstRatio, maxRad, worstVsArb, worstAB);
  check(rates == 0, "pgdensity: the proposal's K and rates are the table's, bit for bit");
  check(worstAB <= 1e-11, "pgdensity: the proposal's alpha and beta match the table's");
  check(failed == 0, "pgdensity: every q pin's enclosure contains its reference at every level");
  check(undecided == 0, "pgdensity: every q pin is decided by the last level");
}

/// Per (b, psi) of the f pins: the largest radius of log f and of log q
/// (at the rule's K, the engine's own proposal) over z, and the step floor.
// c0's radius charges (17 + 3n) u^2 (sum + x): the quotients and products,
// one double-double addition a term and the final subtraction.
void testResidual() {
  int failed = 0;
  std::size_t largestN = 0;
  for (const QPin& p : kQPins) {
    pg::Proposal prop(p.b, p.psi);
    pg::GammaSumDensity q;
    prop.density(q);
    q.prepare(p.x);
    const std::size_t n = q.numTerms();
    double sigma = 0.0;
    for (std::size_t j = 0; j < n; ++j) sigma += prop.shapes()[j] / (prop.rates()[j] + q.saddle());
    const double need = (17.0 + 3.0 * n) * 0x1p-106 * (sigma + p.x) * (1.0 - 0x1p-40);
    if (!(q.residual().r >= need)) {
      if (failed++ < 5)
        printf("FAIL: pgdensity c0 at b=%g psi=%g z=%g: radius %.3g below %.3g\n", p.b, p.psi, p.z,
               q.residual().r, need);
    }
    largestN = std::max(largestN, n);
  }
  printf("pgdensity c0: %zu q pins, up to %zu terms, %d charged short\n", countOf(kQPins), largestN, failed);
  check(failed == 0, "pgdensity: c0's radius is at least (17 + 3n) u^2 (sum + x) at every q pin");
}

void testFloor(const std::vector<double>& radF) {
  double worst = 0.0;
  int undecided = 0;
  printf("pgdensity floor: b psi K | max rad f (Arb's) | max rad q | floor\n");
  for (std::size_t i = 0; i < countOf(kFPins);) {
    std::size_t j = i;
    double rf = 0.0, arb = 0.0, rq = 0.0;
    pg::Proposal prop(kFPins[i].b, kFPins[i].psi);
    pg::GammaSumDensity q;
    prop.density(q);
    for (; j < countOf(kFPins) && kFPins[j].b == kFPins[i].b && kFPins[j].psi == kFPins[i].psi; ++j) {
      rf = std::max(rf, radF[j]);
      arb = std::max(arb, kFPins[j].arbRad);
      q.prepare(kFPins[j].x);
      pg::Enclosure e = pg::logDensityThrough(q, pg::kNumLevels - 1);
      if (!e.finite()) ++undecided;
      rq = std::max(rq, e.rad);
    }
    double fl = pg::stepFloor(rf, rq);
    worst = std::max(worst, fl);
    printf("  %g %g %d | %.3g (%.3g) | %.3g | %.3g\n", kFPins[i].b, kFPins[i].psi, prop.numTerms(), rf,
           arb, rq, fl);
    i = j;
  }
  printf("pgdensity floor: max %.3g against 2.3e-11\n", worst);
  check(undecided == 0, "pgdensity: log q decided at every f pin's point");
  check(worst <= 2.3e-11, "pgdensity: the step floor is at most 2.3e-11 at every (b, psi)");
  // the floor's formula: twice the sum of four radii plus the uniform's grain
  check(pg::stepFloor(0.0, 0.0) == 0x1p-52, "pgdensity: the floor charges the uniform's grain");
  check(pg::stepFloor(0x1p-40, 0x1p-41) == 2.0 * (2.0 * 0x1p-40 + 2.0 * 0x1p-41) + 0x1p-52,
        "pgdensity: the floor is 2 (2 rad_f + 2 rad_q) + 2^(1 - beta)");
}

void testKTable() {
  int failed = 0;
  double worstD = 0.0;
  for (const KCell& c : kKTable) {
    int K = pg::ruleK(c.b, c.psi);
    if (K != c.K) {
      if (failed++ < 10)
        printf("FAIL: pgdensity K rule at b=%g psi=%g: %d, the table's %d (b D at K-1 %.6g, at K %.6g)\n",
               c.b, c.psi, K, c.K, c.bDBelow, c.bDAt);
    }
    worstD = std::max(worstD, std::fabs(c.b * pg::ruleD(c.psi, c.K) / c.bDAt - 1.0));
  }
  printf("pgdensity K table: %zu cells, %d differ, b D in double within %.2g of 40 digits\n",
         countOf(kKTable), failed, worstD);
  check(failed == 0, "pgdensity: the K rule's K at every cell of the K table");
}

/// r_1 below d_1 = pi^2 / 2 + psi^2 / 2, compared in double-double.
void testFirstRate() {
  const dd::DD pi{0x1.921fb54442d18p+1, 0x1.1a62633145c07p-53};
  int failed = 0;
  for (int i = 0; i <= 3000; ++i) {
    double psi = i < 1500 ? i * 0.01 : 15.0 + (i - 1500) * 0.09;
    pg::Proposal prop(1.0, psi, 1);
    dd::DD d1 = dd::add(dd::mul(dd::mul(pi, pi), 0.5), dd::mul(dd::twoProd(psi, psi), 0.5));
    dd::DD gap = dd::sub(d1, dd::fromDouble(prop.rate(0)));
    if (!(gap.hi > 0.0)) {
      if (failed++ < 5) printf("FAIL: pgdensity r_1 %a not below d_1 at psi=%g\n", prop.rate(0), psi);
    }
  }
  check(failed == 0, "pgdensity: r_1 is below d_1 at every psi tested");
}

template <class Fn>
void ddGate(const char* name, Fn fn, const DdCase* cases, std::size_t n, const double* range,
            double stated) {
  double worst = 0.0;
  int outside = 0;
  for (std::size_t i = 0; i < n; ++i) {
    const DdCase& c = cases[i];
    double x = c.xHi + c.xLo;
    if (!(x >= range[0] && x <= range[1])) ++outside;
    dd::DD y = fn(dd::DD{c.xHi, c.xLo});
    worst = std::max(worst, dd::relativeDifference(y, dd::DD{c.yHi, c.yLo}));
  }
  // the reference's own rounding to a double-double is at most 2^-106
  printf("pgdensity dd %s: %zu cases over [%a, %a], largest relative error 2^%.2f, stated 2^%.2f, "
         "assumed 2^-100\n",
         name, n, range[0], range[1], std::log2(worst), std::log2(stated));
  char what[96];
  snprintf(what, sizeof what, "pgdensity: dd %s within its stated bound", name);
  check(outside == 0 && worst <= stated + 0x1p-106, what);
  snprintf(what, sizeof what, "pgdensity: dd %s's stated bound is at most 2^-100", name);
  check(stated <= dd::kFunctionRelErrMax && dd::kFunctionRelErrMax <= 0x1p-100, what);
}

void testDoubleDouble() {
  ddGate("exp", [](dd::DD x) { return dd::exp(x); }, kDdExp, countOf(kDdExp), kDdExpRange,
         dd::kExpRelErr);
  ddGate("log1p", [](dd::DD x) { return dd::log1p(x); }, kDdLog1p, countOf(kDdLog1p), kDdLog1pRange,
         dd::kLog1pRelErr);
  ddGate("log", [](dd::DD x) { return dd::log(x); }, kDdLog, countOf(kDdLog), kDdLogRange,
         dd::kLogRelErr);
  ddGate("sin", [](dd::DD x) { return dd::sin(x); }, kDdSin, countOf(kDdSin), kDdSinRange,
         dd::kSinRelErr);
  ddGate("sqrt", [](dd::DD x) { return dd::sqrt(x); }, kDdSqrt, countOf(kDdSqrt), kDdSqrtRange,
         dd::kSqrtRelErr);
}

/// Past the pins: log f and log q decided far into both tails, where the
/// proposal's draws at small b reach (no reference; finite and decided).
void testReach() {
  int undecided = 0, n = 0;
  for (double b : {0.2, 1.0, 16.0, 1e3, 1e5}) {
    for (double psi : {0.0, 8.0, 50.0, 150.0}) {
      pg::Proposal prop(b, psi);
      pg::TargetDensity f(b, psi);
      pg::GammaSumDensity q;
      prop.density(q);
      const double mean = psi == 0.0 ? b / 4.0 : b * std::tanh(psi / 2.0) / (2.0 * psi);
      const double sd = psi == 0.0 ? std::sqrt(b / 24.0)
                                   : std::sqrt(b * (std::sinh(psi) - psi) /
                                               (4.0 * psi * psi * psi * std::cosh(psi / 2.0) *
                                                std::cosh(psi / 2.0)));
      for (double x : {mean * 1e-6, mean * 1e-4, mean * 1e-2, mean + 20.0 * sd, mean + 40.0 * sd}) {
        f.prepare(x);
        q.prepare(x);
        pg::Enclosure ef = pg::logDensityThrough(f, pg::kNumLevels - 1);
        pg::Enclosure eq = pg::logDensityThrough(q, pg::kNumLevels - 1);
        n += 2;
        if (!ef.finite() || !eq.finite()) {
          if (undecided++ < 10)
            printf("FAIL: pgdensity undecided at b=%g psi=%g x=%g: f %g +- %g, q %g +- %g\n", b, psi, x,
                   ef.mid, ef.rad, eq.mid, eq.rad);
        }
      }
    }
  }
  printf("pgdensity reach: %d evaluations, %d undecided\n", n, undecided);
  check(undecided == 0, "pgdensity: log f and log q decided far into both tails");
}

/// The step from the proposal's own draws: every step decided at moderate
/// (b, psi), and the chain's mean within six standard errors (inflated for
/// the chain's correlation) of PG(b, psi)'s.
void testStep() {
  ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, nullptr);
  ext_rng_setSeed(rng, 20261107u);
  pg::ExtRngWordSource words{rng};
  pg::StepWorkspace ws;
  const double cases[][2] = {{2.5, 1.0}, {0.3, 8.0}, {200.0, 0.0}, {0.2, 0.0}, {0.2, 50.0}};
  for (const auto& c : cases) {
    const double b = c[0], psi = c[1];
    pg::Proposal prop(b, psi);
    const double mean = psi == 0.0 ? b / 4.0 : b * std::tanh(psi / 2.0) / (2.0 * psi);
    const double var = psi == 0.0 ? b / 24.0
                                  : b * (std::sinh(psi) - psi) /
                                      (4.0 * psi * psi * psi * std::cosh(psi / 2.0) * std::cosh(psi / 2.0));
    double x = prop.draw(rng);
    const int n = 4000;
    int accepted = 0, close = 0, beyond0 = 0;
    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
      pg::StepResult r = pg::stepRow(prop, x, rng, words, ws);
      x = r.omega;
      accepted += r.accepted;
      close += r.closeCall;
      beyond0 += r.level > 0;
      sum += x;
    }
    const double z = (sum / n - mean) / std::sqrt(4.0 * var / n);
    printf("pgdensity step b=%g psi=%g K=%d: %d steps, accepted %d, past level 0 %d, close calls %d, "
           "mean z %.2f\n",
           b, psi, prop.numTerms(), n, accepted, beyond0, close, z);
    char what[96];
    snprintf(what, sizeof what, "pgdensity: the step at b=%g psi=%g has no close call", b, psi);
    check(close == 0, what);
    snprintf(what, sizeof what, "pgdensity: the step's mean at b=%g psi=%g matches PG's", b, psi);
    check(std::fabs(z) <= 6.0, what);
  }
  ext_rng_destroy(rng);

  // the interval's 53 bits are the first word's 32 and the second's top 21
  std::uint32_t seq[2] = {0xFFFFFFFFu, 0xFFFFF800u};
  int at = 0;
  auto src = [&]() { return seq[at++]; };
  pg::UniformInterval u = pg::drawUniformInterval(src);
  check(u.k == (std::uint64_t{1} << 53) - 1 && u.upper() == 1.0,
        "pgdensity: the uniform interval takes 53 bits from two words");
}

}  // namespace

void runPgdensityTests() {
  std::vector<double> radF;
  testDoubleDouble();
  testKTable();
  testFirstRate();
  testFPins(radF);
  testQPins();
  testResidual();
  testFloor(radF);
  testReach();
  testStep();
}
