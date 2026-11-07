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
// - the one-row step: decided steps and a mean that matches PG(b, psi)'s;
// - the decision on made-up bounds, the uniform's thresholds under the worst
//   library, the non-finite rules, the K rule's guard, and the pieces'
//   bounds where the pins cannot reach (small operands, point balls, the
//   boundary images, the far tail, the words, sumsT near its guard).
//
// Draws come from the suite's own generator, so it neither shifts nor is
// shifted by any other suite's.

#include "assert.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <random>
#include <utility>
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

void testFPins(std::vector<double>& radF, std::vector<double>& midF) {
  int failed = 0, undecided = 0;
  double worstRatio = 0.0, worstVsArb = 0.0, maxRad = 0.0;
  const FPin* worstPin = nullptr;
  radF.assign(countOf(kFPins), kInf);
  midF.assign(countOf(kFPins), kInf);
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
    midF[i] = std::fabs(e.mid);
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

void testFloor(const std::vector<double>& radF, const std::vector<double>& midF) {
  double worst = 0.0;
  int undecided = 0;
  printf("pgdensity floor: b psi K | max rad f (Arb's) | max rad q | floor\n");
  for (std::size_t i = 0; i < countOf(kFPins);) {
    std::size_t j = i;
    double rf = 0.0, arb = 0.0, rq = 0.0, mf = 0.0, mq = 0.0;
    pg::Proposal prop(kFPins[i].b, kFPins[i].psi);
    pg::GammaSumDensity q;
    prop.density(q);
    for (; j < countOf(kFPins) && kFPins[j].b == kFPins[i].b && kFPins[j].psi == kFPins[i].psi; ++j) {
      rf = std::max(rf, radF[j]);
      mf = std::max(mf, midF[j]);
      arb = std::max(arb, kFPins[j].arbRad);
      q.prepare(kFPins[j].x);
      pg::Enclosure e = pg::logDensityThrough(q, pg::kNumLevels - 1);
      if (!e.finite()) ++undecided;
      rq = std::max(rq, e.rad);
      mq = std::max(mq, std::fabs(e.mid));
    }
    double fl = pg::stepFloor(rf, rq, mf, mq);
    worst = std::max(worst, fl);
    printf("  %g %g %d | %.3g (%.3g) | %.3g | %.3g\n", kFPins[i].b, kFPins[i].psi, prop.numTerms(), rf,
           arb, rq, fl);
    i = j;
  }
  printf("pgdensity floor: max %.3g against 2.3e-11\n", worst);
  check(undecided == 0, "pgdensity: log q decided at every f pin's point");
  check(worst <= 2.3e-11, "pgdensity: the step floor is at most 2.3e-11 at every (b, psi)");
  // the floor's formula: twice the half-width, the thresholds' widening and
  // the uniform's grain
  check(pg::stepFloor(0.0, 0.0, 0.0, 0.0) >= 0x1p-52 + 0x1p-50,
        "pgdensity: the floor charges the uniform's grain and the thresholds' widening");
  const double hw = 2.0 * 0x1p-40 + 2.0 * 0x1p-41 + 4.0 * 0x1p-53 * (2.0 * 30.0 + 2.0 * 20.0);
  check(pg::stepFloor(0x1p-40, 0x1p-41, 30.0, 20.0) >= 2.0 * hw + 0x1p-52 + 0x1p-50 &&
            pg::stepFloor(0x1p-40, 0x1p-41, 30.0, 20.0) <= (2.0 * hw + 0x1p-52 + 0x1p-50) * (1.0 + 0x1p-40),
        "pgdensity: the floor is 2 (2 rad_f + 2 rad_q + 8 u (|mid_f| + |mid_q|)) + 2 u_L + 2^(1 - beta)");
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
  // distinct words, so taking the first twice shows
  std::uint32_t seq2[2] = {0x80000000u, 0x00000800u};
  at = 0;
  auto src2 = [&]() { return seq2[at++]; };
  u = pg::drawUniformInterval(src2);
  check(u.k == ((std::uint64_t{1} << 52) | 1u) && at == 2,
        "pgdensity: the uniform interval takes both words, the first word high");
}

pg::Enclosure enc(double mid, double rad) {
  pg::Enclosure e;
  e.mid = mid;
  e.rad = rad;
  return e;
}

/// The accept-or-keep decision on made-up bounds: one level's verdict, the
/// step's thresholds from the uniform, the close calls, and an invalid draw.
void testDecision() {
  const pg::Enclosure z = enc(0.0, 0.0), none;
  using V = pg::Verdict;
  // log r = mid of efy here; accept needs its lower end at accTh
  check(pg::decideLevel(enc(-1.0, 1e-3), z, z, z, -1.5, -3.0) == V::accept &&
            pg::decideLevel(enc(-1.0, 0.1), z, z, z, -1.0, -3.0) == V::undecided &&
            pg::decideLevel(enc(-1.0, 0.1), z, z, z, -0.95, -3.0) == V::undecided,
        "pgdensity: a level accepts only when log r's lower end reaches the threshold");
  check(pg::decideLevel(enc(-2.0, 1e-3), z, z, z, -1.0, -1.5) == V::keep &&
            pg::decideLevel(enc(-2.0, 0.1), z, z, z, -1.0, -2.0) == V::undecided &&
            pg::decideLevel(enc(-2.0, 0.1), z, z, z, -1.0, -2.05) == V::undecided,
        "pgdensity: a level keeps only when log r's upper end is at most the threshold");
  // log r = efy + eqx - efx - eqy: each enclosure enters with its sign
  check(pg::decideLevel(z, enc(-2.0, 1e-3), z, z, -1.0, -1.5) == V::keep &&
            pg::decideLevel(z, z, enc(2.0, 1e-3), z, -1.0, -1.5) == V::keep &&
            pg::decideLevel(z, z, z, enc(-2.0, 1e-3), 1.5, -3.0) == V::accept,
        "pgdensity: log r takes f(y) and q(x) with a plus sign, f(x) and q(y) with a minus");
  // the radius charges 4 u of every midpoint: with 2^31 in a plus slot and
  // in a minus slot, mid is 0 and the charge 8 u 2^31 = 2^-19 keeps it from
  // accepting at -1e-6; a charge that misses either slot is 2^-20 and accepts
  const double M = 0x1p31;
  const int pairs[4][2] = {{0, 2}, {0, 3}, {1, 2}, {1, 3}};
  bool charged = true;
  for (const auto& p : pairs) {
    pg::Enclosure e[4] = {z, z, z, z};
    e[p[0]] = e[p[1]] = enc(M, 0.0);
    charged &= pg::decideLevel(e[0], e[1], e[2], e[3], -1e-6, -1.0) == V::undecided;
  }
  check(charged, "pgdensity: a level's radius charges each midpoint's rounding");
  bool allUndecided = true;
  for (int i = 0; i < 4; ++i) {
    pg::Enclosure e[4] = {z, z, z, z};
    e[i] = none;
    allUndecided &= pg::decideLevel(e[0], e[1], e[2], e[3], -1.0, 1.0) == V::undecided;
    e[i] = enc(kInf, 0.0);
    allUndecided &= pg::decideLevel(e[0], e[1], e[2], e[3], -1.0, 1.0) == V::undecided;
  }
  check(allUndecided, "pgdensity: a level with a non-finite enclosure is undecided");

  // the step: scripted levels, made-up thresholds from k
  struct Script {
    pg::Enclosure lr[pg::kNumLevels];  // log f(y) at each level; the rest 0
    int prepared = 0, calls = 0;
  };
  auto run = [](Script& sc, double x, double y, std::uint64_t k) {
    pg::UniformInterval U{k};
    return pg::decideStep(
        x, y, U, [&]() { ++sc.prepared; },
        [&](int l, pg::Enclosure* e) {
          ++sc.calls;
          e[0] = sc.lr[l];
          e[1] = e[2] = e[3] = enc(0.0, 0.0);
        });
  };
  const double x = 0.25, y = 0.75;
  // k = 1: U in [2^-53, 2^-52); log r inside (log U-, log U+) never decides
  Script s1;
  for (auto& e : s1.lr) e = enc(std::log(1.5 * 0x1p-53), 1e-6);
  pg::StepResult r = run(s1, x, y, 1);
  check(r.closeCall && !r.accepted && r.omega == x && r.level == pg::kNumLevels - 1 &&
            s1.prepared == 1 && s1.calls == pg::kNumLevels,
        "pgdensity: log r between log U- and log U+ at every level is a close call that keeps x");
  Script s2;
  for (auto& e : s2.lr) e = enc(std::log(2.0 * 0x1p-52), 1e-6);
  r = run(s2, x, y, 1);
  check(r.accepted && !r.closeCall && r.omega == y && r.level == 0,
        "pgdensity: log r above log U+ accepts y at level 0");
  Script s3;
  for (auto& e : s3.lr) e = enc(std::log(0.5 * 0x1p-53), 1e-6);
  r = run(s3, x, y, 1);
  check(!r.accepted && !r.closeCall && r.omega == x && r.level == 0,
        "pgdensity: log r below log U- keeps x at level 0");
  // undecided levels go on: non-finite, then a straddle, then accept
  Script s4;
  s4.lr[0] = enc(-40.0, kInf);
  for (int l = 1; l < pg::kNumLevels; ++l) s4.lr[l] = enc(std::log(2.0 * 0x1p-52), 1e-6);
  r = run(s4, x, y, 1);
  check(r.accepted && r.omega == y && r.level == 1,
        "pgdensity: a non-finite level is undecided and the step goes on to the next");
  Script s5;
  for (auto& e : s5.lr) e = enc(std::log(0.5 * 0x1p-53), 1e-6);
  r = run(s5, x, y, 0);
  check(r.closeCall && r.omega == x,
        "pgdensity: at k = 0 nothing keeps (log U- is -inf): a close call");
  // an invalid draw: a close call, no level evaluated
  bool bad = true;
  for (double yb : {std::nan(""), kInf, -kInf, 0.0, -1.0}) {
    Script s6;
    for (auto& e : s6.lr) e = enc(0.0, 0.0);  // would accept
    r = run(s6, x, yb, 12345);
    bad &= r.closeCall && !r.accepted && r.omega == x && s6.prepared == 0 && s6.calls == 0;
  }
  check(bad, "pgdensity: a draw that is not positive and finite is a close call that keeps x");

  // an invalid proposal: the rule finds no K (psi 1e5 at b 1e5, psi 1e6) or
  // psi is not finite; the step is a close call and draws nothing
  int noK = 0, invalid = 0, stepped = 0, taken = 0;
  ext_rng* rng = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, nullptr);
  ext_rng_setSeed(rng, 7u);
  auto words = [&]() {
    ++taken;
    return 0x80000000u;
  };
  pg::StepWorkspace ws;
  const double cases[][2] = {{1e5, 1e5}, {0.2, 1e6}, {1.0, 1e6}, {1e5, 1e6}, {1.0, kInf},
                             {1.0, std::nan("")}, {kInf, 1.0}, {0.0, 1.0}, {-1.0, 1.0}};
  for (const auto& c : cases) {
    noK += pg::ruleK(c[0], c[1]) == -1;
    pg::Proposal prop(c[0], c[1]);
    invalid += !prop.valid() && std::isnan(prop.draw(rng));
    pg::StepResult sr = pg::stepRow(prop, 0.5, rng, words, ws);
    stepped += sr.closeCall && sr.omega == 0.5 && !sr.accepted;
  }
  const int n = static_cast<int>(countOf(cases));
  printf("pgdensity: no K at %d of %d inputs, invalid proposals %d, close-call steps %d, "
         "words taken %d\n",
         noK, n, invalid, stepped, taken);
  check(noK == n && invalid == n && stepped == n,
        "pgdensity: where the K rule finds no K the proposal is invalid and the step a close call");
  check(taken == 0, "pgdensity: a step on an invalid proposal takes no uniform words");
  // a stated K does not make a bad shape or psi valid
  check(!pg::Proposal(0.0, 1.0, 3).valid() && !pg::Proposal(1.0, kInf, 3).valid() &&
            !pg::Proposal(-1.0, 1.0, 3).valid() && !pg::Proposal(1.0, std::nan(""), 3).valid(),
        "pgdensity: a proposal with a stated K and a bad b or psi is invalid");
  check(pg::Proposal(2.5, 1.0).valid(), "pgdensity: a proposal the rule builds is valid");
  ext_rng_destroy(rng);
}

/// The uniform's thresholds under the worst library (A2s): for l anywhere
/// within u_L |log v| of log v the widened bounds are on their side of
/// log v, and their slack is at most (2 u_L + 2 u) |log v| + 2^-999.
void testThresholds() {
  std::mt19937_64 g(5);
  const double uL = pg::kLibRelErr;
  int outside = 0, wide = 0, n = 0, ownOutside = 0;
  double worstSlack = 0.0;
  for (int i = 0; i < 200000; ++i) {
    std::uint64_t k = (g() >> 11) >> (g() % 53);
    if (k == 0) continue;
    pg::UniformInterval U{k};
    for (int side = 0; side < 2; ++side) {
      const double v = side == 0 ? U.upper() : U.lower();
      if (v == 1.0) continue;
      const dd::DD L = dd::log(dd::DD{v, 0.0});
      const double tolRef = 2.0 * dd::kLogRelErr * std::fabs(L.hi);
      const double lim = uL * std::fabs(L.hi) * (1.0 - 0x1p-50);
      auto dist = [&](double d) { return (d - L.hi) - L.lo; };
      // the doubles farthest below and above log v within u_L |log v|
      double lo = L.hi - lim, hi = L.hi + lim;
      while (std::fabs(dist(lo)) > lim) lo = std::nextafter(lo, L.hi);
      while (std::fabs(dist(hi)) > lim) hi = std::nextafter(hi, L.hi);
      while (std::fabs(dist(std::nextafter(lo, -kInf))) <= lim) lo = std::nextafter(lo, -kInf);
      while (std::fabs(dist(std::nextafter(hi, kInf))) <= lim) hi = std::nextafter(hi, kInf);
      ++n;
      if (side == 0) {
        outside += dist(pg::logWidenUp(lo)) < -tolRef;
        double s = dist(pg::logWidenUp(hi));
        worstSlack = std::max(worstSlack, s / std::fabs(L.hi));
        wide += s > (2.0 * uL + 2.0 * pg::kU) * (1.0 + 0x1p-40) * std::fabs(L.hi) + 0x1p-999;
        ownOutside += dist(U.logUpperUp()) < -tolRef;
      } else {
        outside += dist(pg::logWidenDown(hi)) > tolRef;
        double s = -dist(pg::logWidenDown(lo));
        worstSlack = std::max(worstSlack, s / std::fabs(L.hi));
        wide += s > (2.0 * uL + 2.0 * pg::kU) * (1.0 + 0x1p-40) * std::fabs(L.hi) + 0x1p-999;
        ownOutside += dist(U.logLowerDown()) > tolRef;
      }
    }
  }
  printf("pgdensity thresholds: %d bounds under the worst library, %d on the wrong side, %d too wide, "
         "worst slack %.4f (2 u_L + 2 u) |log v|; with this library %d on the wrong side\n",
         n, outside, wide, worstSlack / (2.0 * uL + 2.0 * pg::kU), ownOutside);
  check(outside == 0 && ownOutside == 0,
        "pgdensity: the uniform's log bounds hold under any library within u_L");
  check(wide == 0, "pgdensity: the uniform's log bounds are within the floor's slack");

  // the widening's sufficient condition, exact in doubles and with no
  // reference: d = widenUp(l) - l and l - widenDown(l) are exact (Sterbenz),
  // and d (1 - u_L) >= u_L |l| is decided by a twoSum, at l on both sides of
  // powers of two and at random l
  auto enough = [&](double d, double l) {
    double b = d * uL, c = uL * std::fabs(l);
    double s = b + c, bb = s - b, e = (b - (s - bb)) + (c - bb);
    return d > s || (d == s && e <= 0.0);
  };
  int nExact = 0, short_ = 0;
  auto exact = [&](double l) {
    nExact += 2;
    short_ += !enough(pg::logWidenUp(l) - l, l);
    short_ += !enough(l - pg::logWidenDown(l), l);
  };
  for (int k = -60; k < 7; ++k) {
    double a = -std::ldexp(1.0, k), c = a;
    for (int j = 0; j < 300; ++j) {
      exact(a);
      a = std::nextafter(a, -kInf);
      c = std::nextafter(c, 0.0);
      exact(c);
    }
  }
  std::mt19937_64 g2(3);
  std::uniform_real_distribution<double> U12(1.0, 2.0);
  for (int i = 0; i < 100000; ++i) exact(-std::ldexp(U12(g2), static_cast<int>(g2() % 66) - 60));
  printf("pgdensity thresholds: the widening checked exactly at %d bounds, %d short of u_L |l| / (1 - u_L)\n",
         nExact, short_);
  check(short_ == 0, "pgdensity: the widening covers u_L |l| / (1 - u_L), exactly in doubles");
  check(pg::UniformInterval{0}.logLowerDown() == -kInf &&
            pg::UniformInterval{(std::uint64_t{1} << 53) - 1}.logUpperUp() > 0.0,
        "pgdensity: log U- is -inf at k = 0 and log U+ above 0 at U+ = 1");
}

/// logDensity's non-finite rule, with a family whose truncation has no
/// finite bound: the far tail is infinite and the nodes past the main sum
/// are not finite, so the enclosure must be undecided.
struct TailBlindTarget : pg::TargetDensity {
  using pg::TargetDensity::TargetDensity;
  mutable bool pastSum = false;
  double farTail(double, double) const {
    pastSum = true;
    return kInf;
  }
  pg::CBall nodeExponent(double y) const {
    if (pastSum) return pg::CBall{0.0, 0.0, kInf};
    return pg::TargetDensity::nodeExponent(y);
  }
};

void testNonFinite() {
  int finiteReal = 0, finiteBlind = 0;
  for (double b : {0.5, 2.5, 30.0}) {
    for (double xm : {0.3, 1.0, 3.0}) {
      const double x = b / 4.0 * xm;
      pg::TargetDensity f(b, 1.0);
      f.prepare(x);
      TailBlindTarget t(b, 1.0);
      t.prepare(x);
      for (int l = 0; l < pg::kNumLevels; ++l) {
        finiteReal += pg::logDensity(f, l).finite();
        t.pastSum = false;
        finiteBlind += pg::logDensity(t, l).finite();
      }
    }
  }
  printf("pgdensity non-finite truncation: %d of %d real evaluations finite, %d of the blind ones\n",
         finiteReal, 9 * pg::kNumLevels, finiteBlind);
  check(finiteReal == 9 * pg::kNumLevels && finiteBlind == 0,
        "pgdensity: a truncation without a finite bound leaves the level undecided");

  pg::LogSumExpUpper a, b, c, d, e;
  a.add(1.0);
  a.add(std::nan(""));
  a.add(2.0);
  b.add(1.0);
  b.add(2.0);
  c.add(1.0);
  c.add(-kInf);
  c.add(2.0);
  d.add(1.0);
  d.add(kInf);
  check(a.value() == kInf && c.value() == b.value() && d.value() == kInf && e.value() == -kInf &&
            b.value() >= 2.0 + std::log1p(std::exp(-1.0)),
        "pgdensity: the log-sum bound is infinite with a NaN or +inf term and skips -inf");
}

/// cmul with an operand below 2^-500 (down to subnormal): the error term
/// 2.25 u |a| |b| must not underflow before the product.
void testSmallOperands() {
  std::mt19937_64 g(31);
  std::uniform_real_distribution<double> U(0.0, 1.0);
  int outside = 0, n = 0;
  for (int i = 0; i < 4000; ++i) {
    const int ea = -975 - static_cast<int>(95.0 * U(g));  // |a| near 2^-975 to 2^-1070
    const int ep = -600 - static_cast<int>(300.0 * U(g));  // |a b| normal
    const double s = U(g) < 0.5 ? -1.0 : 1.0;
    pg::CBall a{std::ldexp(0.5 + U(g), ea), s * std::ldexp(0.5 + U(g), ea - static_cast<int>(4.0 * U(g))),
                0.0};
    pg::CBall b{std::ldexp(0.5 + U(g), ep - ea), -s * std::ldexp(0.5 + U(g), ep - ea), 0.0};
    if (i % 2) std::swap(a, b);
    pg::CBall p = pg::cmul(a, b);
    // the exact product's parts as double-doubles (the products' low parts are normal)
    dd::DD re = dd::sub(dd::twoProd(a.re, b.re), dd::twoProd(a.im, b.im));
    dd::DD im = dd::add(dd::twoProd(a.re, b.im), dd::twoProd(a.im, b.re));
    double dre = (re.hi - p.re) + re.lo, dim = (im.hi - p.im) + im.lo;
    ++n;
    if (!(std::hypot(dre, dim) <= p.r)) ++outside;
  }
  printf("pgdensity small operands: cmul %d products, %d outside\n", n, outside);
  check(outside == 0, "pgdensity: cmul's bound holds with an operand below 2^-500");

  // cdiv with a divisor below 2^-400 or above 2^400, the quotient normal:
  // the reference from operands scaled exactly to near 1
  int qOutside = 0, qn = 0;
  for (int i = 0; i < 4000; ++i) {
    const int eb = (i % 2 ? 1 : -1) * (420 + static_cast<int>(480.0 * U(g)));
    pg::CBall a{0.5 + U(g), U(g) - 0.5, 0.0};
    pg::CBall b{std::ldexp(0.5 + U(g), eb), std::ldexp(U(g) - 0.5, eb), 0.0};
    pg::CBall q = pg::cdiv(a, b);
    if (!std::isfinite(q.r)) continue;
    const double br = std::ldexp(b.re, -eb), bi = std::ldexp(b.im, -eb);
    dd::DD n2 = dd::add(dd::twoProd(br, br), dd::twoProd(bi, bi));
    dd::DD re = dd::div(dd::add(dd::twoProd(a.re, br), dd::twoProd(a.im, bi)), n2);
    dd::DD im = dd::div(dd::sub(dd::twoProd(a.im, br), dd::twoProd(a.re, bi)), n2);
    const double dre = (re.hi - std::ldexp(q.re, eb)) + re.lo, dim = (im.hi - std::ldexp(q.im, eb)) + im.lo;
    ++qn;
    if (!(std::hypot(dre, dim) <= std::ldexp(q.r, eb))) ++qOutside;
  }
  printf("pgdensity small operands: cdiv %d quotients, %d outside\n", qn, qOutside);
  check(qn > 0 && qOutside == 0, "pgdensity: cdiv's bound holds with a divisor outside [2^-400, 2^400]");
}

/// Point balls (radius 0) against double-double references: each result
/// ball holds the exact value, so the midpoint's rounding charge suffices.
void testPointBalls() {
  std::mt19937_64 g(41);
  std::uniform_real_distribution<double> U(-1.0, 1.0);
  int divOut = 0, sqrtOut = 0, n = 0;
  for (int i = 0; i < 20000; ++i) {
    const double s = std::ldexp(1.0, static_cast<int>(20.0 * U(g)));
    pg::CBall a{s * U(g), s * U(g), 0.0}, b{U(g), U(g), 0.0};
    if (!(std::fabs(b.re) + std::fabs(b.im) > 1e-3)) continue;
    ++n;
    pg::CBall q = pg::cdiv(a, b);
    dd::DD n2 = dd::add(dd::twoProd(b.re, b.re), dd::twoProd(b.im, b.im));
    dd::DD re = dd::div(dd::add(dd::twoProd(a.re, b.re), dd::twoProd(a.im, b.im)), n2);
    dd::DD im = dd::div(dd::sub(dd::twoProd(a.im, b.re), dd::twoProd(a.re, b.im)), n2);
    divOut += !(std::hypot((re.hi - q.re) + re.lo, (im.hi - q.im) + im.lo) <= q.r);
    // the square root in the closed upper half plane
    pg::CBall w{a.re, std::fabs(a.im), 0.0};
    pg::CBall r = pg::csqrtUpper(w);
    dd::DD mod = dd::sqrt(dd::add(dd::twoProd(w.re, w.re), dd::twoProd(w.im, w.im)));
    dd::DD p, t;
    if (w.re >= 0.0) {
      p = dd::sqrt(dd::mul(dd::add(mod, w.re), 0.5));
      t = dd::div(dd::DD{w.im, 0.0}, dd::mul(p, 2.0));
    } else {
      t = dd::sqrt(dd::mul(dd::sub(mod, w.re), 0.5));
      p = dd::div(dd::DD{w.im, 0.0}, dd::mul(t, 2.0));
    }
    sqrtOut += !(std::hypot((p.hi - r.re) + p.lo, (t.hi - r.im) + t.lo) <= r.r);
  }
  printf("pgdensity point balls: %d cases, cdiv %d outside, csqrt %d outside\n", n, divOut, sqrtOut);
  check(divOut == 0 && sqrtOut == 0, "pgdensity: cdiv's and csqrt's midpoint charges hold at point balls");
}

/// Past the saddle cap: log f decided far to the right, x 400 b and 1000 b,
/// where the saddle's v is capped below the closed forms' guard (log f is
/// undecided from about 1e4 in x at b <= 2, 1,422 b at b 16 and 403 b at
/// b 1e3, log f there below -4.7e4: no draw at a fixed b gets there).
void testFarRight() {
  int undecided = 0, n = 0;
  for (double b : {0.2, 1.0, 16.0}) {
    for (double psi : {0.0, 8.0}) {
      pg::TargetDensity f(b, psi);
      for (double m : {400.0, 1000.0}) {
        f.prepare(m * b);
        ++n;
        undecided += !pg::logDensityThrough(f, pg::kNumLevels - 1).finite();
      }
    }
  }
  printf("pgdensity far right: %d evaluations, %d undecided\n", n, undecided);
  check(undecided == 0, "pgdensity: log f decided past the saddle cap");
}

/// The input's radius carried through: images of points on the input disk's
/// boundary lie in the output ball.
void testBoundary() {
  std::mt19937_64 g(23);
  std::uniform_real_distribution<double> U(0.0, 1.0);
  auto probe = [&](const char* name, auto f, double lo, double hi) {
    int miss = 0;
    for (int i = 0; i < 2000; ++i) {
      double a = lo + (hi - lo) * U(g), b = lo + (hi - lo) * U(g);
      double r = 1e-6 * (std::fabs(a) + std::fabs(b));
      pg::CBall o = f(pg::CBall{a, b, r});
      if (!std::isfinite(o.r)) continue;
      for (int k = 0; k < 16; ++k) {
        double t = 2.0 * pg::kPi * k / 16.0;
        pg::CBall q = f(pg::CBall{a + r * std::cos(t), b + r * std::sin(t), 0.0});
        if (std::hypot(q.re - o.re, q.im - o.im) > o.r + 1e-13 * std::hypot(o.re, o.im)) {
          ++miss;
          break;
        }
      }
    }
    printf("pgdensity boundary: %-6s %d of 2000 balls miss a boundary image\n", name, miss);
    char what[96];
    snprintf(what, sizeof what, "pgdensity: %s carries its input's radius", name);
    check(miss == 0, what);
  };
  probe("l1pm", [](pg::CBall z) { return pg::clog1pMinus(z); }, -0.17, 0.17);
  probe("l1pm+", [](pg::CBall z) { return pg::clog1pMinus(z); }, -0.6, 0.6);
  probe("expm1", [](pg::CBall z) { return pg::cexpm1(z); }, -0.3, 0.3);
  probe("exp", [](pg::CBall z) { return pg::cexp(z); }, -5.0, 5.0);
  probe("log1p", [](pg::CBall z) { return pg::clog1p(z); }, -0.5, 2.0);
  probe("sqrt", [](pg::CBall z) { return pg::csqrtUpper(pg::CBall{z.re, z.im + 10.0, z.r}); }, -3.0, 3.0);
  probe("mul", [](pg::CBall z) { return pg::cmul(z, pg::CBall{0.7, -1.3, 1e-7}); }, -3.0, 3.0);
  probe("div", [](pg::CBall z) { return pg::cdiv(pg::CBall{1.1, 0.4, 1e-7}, z); }, 0.5, 3.0);
}

/// q's dyadic far tail at least the explicit node sum sum_{k h >= y0} |G(k h)|.
void testFarTail() {
  double worst = kInf;
  int n = 0, below = 0;
  for (double b : {0.2, 2.5, 1e3}) {
    for (double psi : {8.0, 50.0, 150.0}) {
      for (double zs : {-3.0, 0.5, 12.0}) {
        pg::Proposal prop(b, psi);
        pg::GammaSumDensity q;
        prop.density(q);
        const double A = 0.5 * psi, mean = b / (2.0 * psi) * std::tanh(A);
        const double sd = std::sqrt(b / (4.0 * psi * psi * psi) * (std::sinh(psi) - psi) /
                                    (std::cosh(A) * std::cosh(A)));
        double x = mean + zs * sd;
        if (!(x > 0.0)) x = mean / 20.0;
        q.prepare(x);
        pg::EvalCounts c;
        pg::logDensity(q, 0, &c);
        const double h = c.h, rho = q.rho();
        if (!(h > 0.0)) continue;
        const double yW = 0.5 * (-rho + std::sqrt(rho * rho + 1.0 / x));
        for (double m : {1.0, 2.0, 8.0}) {
          const double y0 = std::max(yW * m, h);
          const long k0 = static_cast<long>(std::ceil(y0 / h));
          double s = 0.0;
          for (long k = k0; k < k0 + 4000000; ++k) {
            const double y = k * h;
            pg::CBall E = q.nodeExponent(y);
            double t = std::exp(E.re) * std::hypot(rho, y);
            s += t;
            if (t < 1e-30 * s && k > k0 + 100) break;
            if (!(s > 0.0) && E.re < -800.0) break;
          }
          if (!(s > 0.0)) continue;
          const double ft = q.farTail(k0 * h, h);
          ++n;
          worst = std::min(worst, ft / s);
          below += !(ft >= s);
        }
      }
    }
  }
  printf("pgdensity far tail: %d cases, %d below the node sum, smallest far tail / node sum %.4g\n", n,
         below, worst);
  check(n > 0 && below == 0, "pgdensity: q's far tail is at least the explicit node sum");
}

/// The word source against the Mersenne Twister's uniforms: each word w is
/// the uniform's 32 bits, u = w 2^-32 exactly.
void testWordSource() {
  ext_rng* r1 = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, nullptr);
  ext_rng* r2 = ext_rng_create(EXT_RNG_ALGORITHM_MERSENNE_TWISTER, nullptr);
  ext_rng_setSeed(r1, 99u);
  ext_rng_setSeed(r2, 99u);
  pg::ExtRngWordSource words{r1};
  int differ = 0;
  for (int i = 0; i < 100000; ++i) {
    std::uint32_t w = words();
    double u = ext_rng_simulateContinuousUniform(r2);
    differ += !(w == 0 ? u < 0x1p-32 : std::ldexp(static_cast<double>(w), -32) == u);
  }
  ext_rng_destroy(r1);
  ext_rng_destroy(r2);
  check(differ == 0, "pgdensity: the word source's words are the Mersenne Twister's uniforms' bits");
}

/// sumsT near its guard against a double-double reference: accepted where
/// pi/2 - sqrt(-W) > 2^-11 and not below, within kTRelErr where accepted.
void testSumsT() {
  std::mt19937_64 g(11);
  std::uniform_real_distribution<double> U(0.0, 1.0);
  int n = 0, accepted = 0, wrongGuard = 0, outside = 0;
  double worst = 0.0;
  for (int i = 0; i < 20000; ++i) {
    const double gap = std::ldexp(1.0, -21) * std::pow(2.0, 20.0 * U(g));  // 2^-21 to 2^-1
    const double v0 = 0.5 * pg::kPi - gap;
    const double W = -(v0 * v0);
    double t1, t2;
    const bool ok = pg::sumsT(W, t1, t2);
    const dd::DD v = dd::sqrt(dd::DD{-W, 0.0});
    const dd::DD gp = dd::sub(dd::kPiOver2, v);
    if (!(gp.hi > 0.0)) continue;
    ++n;
    if (gp.hi > 0x1p-11 * (1.0 + 0x1p-30) && !ok) ++wrongGuard;
    if (gp.hi < 0x1p-11 * (1.0 - 0x1p-30) && ok) ++wrongGuard;
    if (!ok) continue;
    ++accepted;
    const dd::DD sn = dd::sin(v), cs = dd::sin(gp);
    const dd::DD tn = dd::div(sn, cs);
    const dd::DD T1 = dd::div(tn, dd::mul(v, 2.0));
    const dd::DD T2 = dd::div(dd::sub(dd::div(v, dd::mul(cs, cs)), tn), dd::mul(dd::mul(v, 4.0), -W));
    const double e = std::max(std::fabs(dd::toDouble(dd::sub(dd::DD{t1, 0.0}, T1)) / T1.hi),
                              std::fabs(dd::toDouble(dd::sub(dd::DD{t2, 0.0}, T2)) / T2.hi));
    worst = std::max(worst, e);
    outside += !(e <= pg::kTRelErr);
  }
  printf("pgdensity sumsT near the guard: %d W, %d accepted, %d on the wrong side of the guard, "
         "%d outside, worst %.3g of kTRelErr\n",
         n, accepted, wrongGuard, outside, worst / pg::kTRelErr);
  check(wrongGuard == 0 && outside == 0 && accepted > 0,
        "pgdensity: sumsT is within kTRelErr where its guard accepts W");
}

}  // namespace

void runPgdensityTests() {
  std::vector<double> radF, midF;
  testDoubleDouble();
  testKTable();
  testFirstRate();
  testFPins(radF, midF);
  testQPins();
  testResidual();
  testFloor(radF, midF);
  testReach();
  testStep();
  testDecision();
  testThresholds();
  testNonFinite();
  testSmallOperands();
  testPointBalls();
  testFarRight();
  testBoundary();
  testFarTail();
  testWordSource();
  testSumsT();
}
