// The standalone build's copy of R's math routines (standalone/nmath) against
// R's own (libR's Rmath), on a grid over the seven routines the engine and the
// support library call and their flag combinations.
//
// This binary links both: the copy's symbols carry a dbarts_nmath_ prefix, so
// neither shadows the other. On the reference host (arm64 macOS, R 4.6.x, the
// copy compiled with R's own C flags) every value must be bitwise equal, which
// is what makes a standalone build's draws R's draws; elsewhere the gate is a
// relative difference of 1e-13. Moving R's minor version moves the pin: the
// copy is re-taken from that R and this grid rerun.
//
// The grid's inputs come from the suite's own generator, so it neither shifts
// nor is shifted by any other suite's draws.

#include "assert.hpp"

#include <cinttypes>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <vector>

#include <Rversion.h>

#include <external/nmath.h>  // the R build's mapping: ext_ names -> Rf_ names

extern "C" {
double dbarts_nmath_dnorm4(double, double, double, int);
double dbarts_nmath_pnorm5(double, double, double, int, int);
double dbarts_nmath_qnorm5(double, double, double, int, int);
double dbarts_nmath_dt(double, double, int);
double dbarts_nmath_qt(double, double, int, int);
double dbarts_nmath_qchisq(double, double, int, int);
double dbarts_nmath_lchoose(double, double);
}

#if defined(__APPLE__) && defined(__aarch64__) && R_VERSION >= R_Version(4, 6, 0) && \
  R_VERSION < R_Version(4, 7, 0)
#  define NMATH_GRID_BITWISE 1
#endif

namespace {

enum Routine { PNORM, QNORM, DNORM, DT, QT, QCHISQ, LCHOOSE, NUM_ROUTINES };
const char* const routineNames[NUM_ROUTINES] = {
  "pnorm", "qnorm", "dnorm", "dt", "qt", "qchisq", "lchoose"};

struct Point {
  Routine routine;
  double a, b, c;  // x or p; then mu and sd, or df; lchoose's n and k
  int lowerTail, logP;
};

double evalR(const Point& p) {
  switch (p.routine) {
    case PNORM: return ext_pnorm(p.a, p.b, p.c, p.lowerTail, p.logP);
    case QNORM: return ext_qnorm(p.a, p.b, p.c, p.lowerTail, p.logP);
    case DNORM: return ext_dnorm(p.a, p.b, p.c, p.logP);
    case DT: return ext_dt(p.a, p.b, p.logP);
    case QT: return ext_qt(p.a, p.b, p.lowerTail, p.logP);
    case QCHISQ: return ext_qchisq(p.a, p.b, p.lowerTail, p.logP);
    case LCHOOSE: return ext_lchoose(p.a, p.b);
    case NUM_ROUTINES: break;
  }
  return std::numeric_limits<double>::quiet_NaN();
}

double evalCopy(const Point& p) {
  switch (p.routine) {
    case PNORM: return dbarts_nmath_pnorm5(p.a, p.b, p.c, p.lowerTail, p.logP);
    case QNORM: return dbarts_nmath_qnorm5(p.a, p.b, p.c, p.lowerTail, p.logP);
    case DNORM: return dbarts_nmath_dnorm4(p.a, p.b, p.c, p.logP);
    case DT: return dbarts_nmath_dt(p.a, p.b, p.logP);
    case QT: return dbarts_nmath_qt(p.a, p.b, p.lowerTail, p.logP);
    case QCHISQ: return dbarts_nmath_qchisq(p.a, p.b, p.lowerTail, p.logP);
    case LCHOOSE: return dbarts_nmath_lchoose(p.a, p.b);
    case NUM_ROUTINES: break;
  }
  return std::numeric_limits<double>::quiet_NaN();
}

// splitmix64: the grid is the same on every host and compiler
struct GridRng {
  uint64_t state = 0x6A09E667F3BCC909ull;
  double unif() {
    uint64_t z = (state += 0x9E3779B97F4A7C15ull);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
    z ^= z >> 31;
    return ((double) (z >> 11) + 0.5) * 0x1.0p-53;  // in (0, 1)
  }
  double normal(double sd) {  // Box-Muller; exactness is not needed here
    double u = unif(), v = unif();
    return sd * std::sqrt(-2.0 * std::log(u)) * std::cos(6.283185307179586 * v);
  }
  double logUniform(double lo, double hi) {
    return std::exp(std::log(lo) + (std::log(hi) - std::log(lo)) * unif());
  }
};

void addTails(std::vector<Point>& grid, Routine r, double a, double b, double c,
              bool logInput) {
  for (int lt = 0; lt <= 1; ++lt)
    for (int lg = 0; lg <= 1; ++lg)
      grid.push_back({r, lg && logInput ? std::log(a) : a, b, c, lt, lg});
}

std::vector<Point> buildGrid() {
  GridRng g;
  std::vector<Point> grid;

  // pnorm: the probe's grid, x on [-40, 40] and N(0, 5), every flag pair;
  // then a location and scale, as the censored likelihood calls it
  for (int i = 0; i < 1601; ++i)
    addTails(grid, PNORM, -40.0 + 80.0 * i / 1600.0, 0.0, 1.0, false);
  for (int i = 0; i < 400; ++i) addTails(grid, PNORM, g.normal(5.0), 0.0, 1.0, false);
  for (int i = 0; i < 500; ++i) {
    double x = g.normal(5.0), mu = g.normal(2.0), sd = g.logUniform(0.1, 10.0);
    addTails(grid, PNORM, x, mu, sd, false);
  }

  // qnorm: uniform p, p down to 1e-300, 1 - p down to 1e-15, and the
  // inversion generator's u1 / BIG-like values; every flag pair, the log
  // ones on log(p)
  std::vector<double> probabilities;
  for (int i = 0; i < 1500; ++i) probabilities.push_back(g.unif());
  for (int i = 0; i < 300; ++i)
    probabilities.push_back(std::pow(10.0, -(1.0 + 299.0 * i / 299.0)));
  for (int i = 0; i < 100; ++i)
    probabilities.push_back(1.0 - std::pow(10.0, -(1.0 + 14.0 * i / 99.0)));
  for (int i = 1; i <= 1000; ++i) probabilities.push_back(i / 1001.0);
  for (double p : probabilities) addTails(grid, QNORM, p, 0.0, 1.0, true);
  for (int i = 0; i < 500; ++i)
    addTails(grid, QNORM, g.unif(), g.normal(2.0), g.logUniform(0.1, 10.0), true);

  // dnorm: N(0, 10), and with a location and scale; both log flags
  for (int i = 0; i < 1000; ++i) {
    double x = g.normal(10.0);
    for (int lg = 0; lg <= 1; ++lg) grid.push_back({DNORM, x, 0.0, 1.0, 1, lg});
  }
  for (int i = 0; i < 500; ++i) {
    double x = g.normal(10.0), mu = g.normal(2.0), sd = g.logUniform(0.1, 10.0);
    for (int lg = 0; lg <= 1; ++lg) grid.push_back({DNORM, x, mu, sd, 1, lg});
  }

  // dt: N(0, 5) at df log-uniform on [0.2, 1e5]; both log flags
  for (int i = 0; i < 1000; ++i) {
    double x = g.normal(5.0), df = g.logUniform(0.2, 1e5);
    for (int lg = 0; lg <= 1; ++lg) grid.push_back({DT, x, df, 0.0, 1, lg});
  }

  // qt and qchisq: the 0.75 and 0.5 points the engine reads at df from 0.05
  // to 1e7 and 1 to 200; then uniform p at df log-uniform on [0.5, 1000],
  // every flag pair
  std::vector<double> dfs;
  for (int i = 0; i < 1500; ++i)
    dfs.push_back(std::exp(std::log(0.05) + (std::log(1e7) - std::log(0.05)) * i / 1499.0));
  for (int i = 1; i <= 200; ++i) dfs.push_back(i);
  for (double df : dfs) grid.push_back({QT, 0.75, df, 0.0, 1, 0});
  for (double df : dfs) grid.push_back({QCHISQ, 0.5, df, 0.0, 1, 0});
  for (int i = 0; i < 500; ++i) {
    double p = g.unif(), df = g.logUniform(0.5, 1000.0);
    addTails(grid, QT, p, df, 0.0, true);
  }
  for (int i = 0; i < 500; ++i) {
    double p = g.unif(), df = g.logUniform(0.5, 1000.0);
    addTails(grid, QCHISQ, p, df, 0.0, true);
  }

  // edges each routine answers without a warning: infinities, NaN, the ends
  // of the unit interval and past them, a zero or negative scale or df
  const double inf = std::numeric_limits<double>::infinity();
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double xs[] = {-inf, inf, nan, 0.0, -0.0};
  for (double x : xs) {
    addTails(grid, PNORM, x, 0.0, 1.0, false);
    for (int lg = 0; lg <= 1; ++lg) {
      grid.push_back({DNORM, x, 0.0, 1.0, 1, lg});
      grid.push_back({DT, x, 3.0, 0.0, 1, lg});
    }
  }
  const double ps[] = {0.0, 1.0, -0.5, 1.5, nan};
  for (double p : ps) {
    for (int lt = 0; lt <= 1; ++lt) {
      grid.push_back({QNORM, p, 0.0, 1.0, lt, 0});
      grid.push_back({QT, p, 4.0, 0.0, lt, 0});
      grid.push_back({QCHISQ, p, 4.0, 0.0, lt, 0});
    }
  }
  for (int lt = 0; lt <= 1; ++lt) {
    grid.push_back({PNORM, 1.0, 0.0, 0.0, lt, 0});
    grid.push_back({PNORM, 1.0, 0.0, -1.0, lt, 0});
    grid.push_back({QNORM, 0.3, 0.0, 0.0, lt, 0});
    grid.push_back({QNORM, 0.3, 0.0, -1.0, lt, 0});
    grid.push_back({QT, 0.3, -1.0, 0.0, lt, 0});
    grid.push_back({QT, 0.3, inf, 0.0, lt, 0});
    grid.push_back({QCHISQ, 0.3, -1.0, 0.0, lt, 0});
    grid.push_back({QCHISQ, 0.3, 0.0, 0.0, lt, 0});
  }
  grid.push_back({DNORM, 1.0, 0.0, -1.0, 1, 0});
  grid.push_back({DT, 1.0, -1.0, 0.0, 1, 0});
  grid.push_back({DT, 1.0, inf, 0.0, 1, 0});

  // lchoose at the binomial row's log choose constants: whole n up to the
  // trial cap of one million and whole k in [0, n], log-uniform in n; every
  // k at small n, where the routine's branches (k < 2, n - k < 2) lie; and
  // the cap's k = n / 2 and n / 1000
  for (int n = 0; n <= 12; ++n)
    for (int k = 0; k <= n; ++k) grid.push_back({LCHOOSE, (double) n, (double) k, 0.0, 1, 0});
  for (int i = 0; i < 2000; ++i) {
    double n = std::floor(g.logUniform(1.0, 1e6));
    double k = std::floor((n + 1.0) * g.unif());
    if (k > n) k = n;
    grid.push_back({LCHOOSE, n, k, 0.0, 1, 0});
    grid.push_back({LCHOOSE, n, n - k, 0.0, 1, 0});
  }
  const double caps[] = {1e6, 1e6 - 1.0, 999983.0};
  for (double n : caps) {
    const double ks[] = {0.0, 1.0, 2.0, std::floor(n / 1000.0), std::floor(n / 2.0), n - 2.0,
                         n - 1.0, n};
    for (double k : ks) grid.push_back({LCHOOSE, n, k, 0.0, 1, 0});
  }

  return grid;
}

uint64_t ulpDistance(double a, double b) {
  int64_t ia, ib;
  std::memcpy(&ia, &a, sizeof ia);
  std::memcpy(&ib, &b, sizeof ib);
  if (ia < 0) ia = INT64_MIN - ia;
  if (ib < 0) ib = INT64_MIN - ib;
  return ia > ib ? (uint64_t) ia - (uint64_t) ib : (uint64_t) ib - (uint64_t) ia;
}

// NaN payloads are not part of the contract (R's NaN is a global R sets at
// startup, the copy's is 0.0 / 0.0): two NaNs agree.
bool agrees(double r, double copy) {
  if (std::isnan(r) || std::isnan(copy)) return std::isnan(r) && std::isnan(copy);
#ifdef NMATH_GRID_BITWISE
  return std::memcmp(&r, &copy, sizeof r) == 0;
#else
  if (r == copy) return true;
  if (std::isinf(r) || std::isinf(copy)) return false;
  return std::fabs(copy - r) <= 1e-13 * std::fabs(r);
#endif
}

}  // namespace

void runNmathTests() {
  std::vector<Point> grid = buildGrid();

  long n[NUM_ROUTINES] = {0}, exact[NUM_ROUTINES] = {0}, failed[NUM_ROUTINES] = {0};
  uint64_t worst[NUM_ROUTINES] = {0};
  int printed = 0;
  for (const Point& p : grid) {
    double r = evalR(p), copy = evalCopy(p);
    ++n[p.routine];
    bool bothNan = std::isnan(r) && std::isnan(copy);
    if (bothNan || std::memcmp(&r, &copy, sizeof r) == 0) ++exact[p.routine];
    else if (!std::isnan(r) && !std::isnan(copy)) {
      uint64_t u = ulpDistance(r, copy);
      if (u > worst[p.routine]) worst[p.routine] = u;
    }
    if (!agrees(r, copy)) {
      ++failed[p.routine];
      if (printed++ < 20)
        printf("FAIL: nmath %s(%a, %a, %a, lower %d, log %d): R %a, copy %a\n",
               routineNames[p.routine], p.a, p.b, p.c, p.lowerTail, p.logP, r, copy);
    }
  }

#ifdef NMATH_GRID_BITWISE
  const char* mode = "bitwise (reference host)";
#else
  const char* mode = "relative 1e-13";
#endif
  long total = 0, totalFailed = 0;
  for (int k = 0; k < NUM_ROUTINES; ++k) {
    printf("nmath %-6s n=%5ld exact=%5ld worst_ulp=%" PRIu64 " failed=%ld\n", routineNames[k],
           n[k], exact[k], worst[k], failed[k]);
    total += n[k];
    totalFailed += failed[k];
  }
  printf("nmath grid: %ld points, gate %s, %ld outside it\n", total, mode, totalFailed);
  for (int k = 0; k < NUM_ROUTINES; ++k) {
    char what[96];
    snprintf(what, sizeof what, "nmath %s: the copy agrees with R's own (%s)", routineNames[k],
             mode);
    check(n[k] > 0 && failed[k] == 0, what);
  }
}
