#ifndef BARTCORE_PGDENSITY_HPP
#define BARTCORE_PGDENSITY_HPP

// Certified log densities of a Polya-Gamma variable PG(b, psi) (the target f)
// and of the proposal q that approximates it (K gammas Ga(b, r_k) plus one
// Ga(alpha, beta) for the remainder), and one row's independence
// Metropolis-Hastings step between them with a decision taken only when the
// enclosures certify it.
//
// Every density is the trapezoid sum on a parabolic contour through the
// saddle point, z = rho + i y in the variable whose square moves the Laplace
// argument, t = t0 + 2 (z^2 - rho^2):
//
//   p(x) = (2 / pi) exp(x t0) L(t0) S,   S = int Re G(y) dy,
//   G(y) = exp(2 x Delta + log(L(t) / L(t0))) (rho + i y),
//   Delta = z^2 - rho^2 = -y^2 + 2 i rho y.
//
// The sum's error is bounded a priori: discretization from two edge lines
// Re z = rho -+ a (their integrals of |G| bounded by cells), truncation by a
// monotone tail. Rounding is carried in midpoint-radius ball arithmetic in
// binary64: every operation's rounding is charged, library functions are
// assumed to return results within kLibRelErr relative error, and every
// radius is computed with its own rounding covered (up()). The lead
// x t0 + log L(t0) is formed in double-double (ddouble.hpp). A log density is
// returned as an Enclosure, the true value within [mid - rad, mid + rad]; a
// non-finite midpoint or radius means undecided.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

#include <external/random.h>

#include "ddouble.hpp"

namespace bartcore {
namespace pgdensity {

namespace dd = ddouble;

/// Unit roundoff of binary64.
inline constexpr double kU = 0x1p-53;
/// Relative error the math library's exp, expm1, log, log1p, sin, cos, tan,
/// sinh, cosh, tanh and atan2 are assumed to stay within (two ulps).
inline constexpr double kLibRelErr = 0x1p-51;
inline constexpr double kPi = 0x1.921fb54442d18p+1;
inline constexpr double kPiSquared = 0x1.3bd3cc9be45dep+3;
inline constexpr double kLog2 = 0x1.62e42fefa39efp-1;
/// log(2 / pi) rounded to a double (error below 2^-55).
inline constexpr double kLogTwoOverPi = -0x1.ce6bb25aa1316p-2;
inline constexpr double kInf = std::numeric_limits<double>::infinity();

/// A computed radius r grown to cover the rounding of the up to about 120
/// operations that formed it and any underflow.
inline double up(double r) { return r * (1.0 + 0x1p-46) + 0x1p-1020; }

/// An enclosure of a real value: within [mid - rad, mid + rad].
struct Enclosure {
  double mid = 0.0;
  double rad = kInf;

  bool finite() const { return std::isfinite(mid) && std::isfinite(rad); }
  double lower() const { return mid - rad; }
  double upper() const { return mid + rad; }
};

/// The intersection of two enclosures of the same value. Neither being
/// finite, or the two not meeting (only possible if one is invalid), returns
/// the other or an undecided enclosure.
inline Enclosure intersect(Enclosure a, Enclosure b) {
  if (!a.finite()) return b;
  if (!b.finite()) return a;
  double lo = std::max(a.lower(), b.lower());
  double hi = std::min(a.upper(), b.upper());
  if (!(lo <= hi)) return Enclosure{};
  double mid = 0.5 * (lo + hi);
  return Enclosure{mid, up(std::max(hi - mid, mid - lo) + kU * std::fabs(mid))};
}

// ---------------------------------------------------------------------------
// Ball arithmetic. A complex ball holds the midpoint (re, im) and a disk
// radius; a real ball a midpoint and a radius.

struct CBall {
  double re;
  double im;
  double r;
};

struct RBall {
  double m;
  double r;
};

inline double l1(double re, double im) { return std::fabs(re) + std::fabs(im); }
/// An upper bound on |re + i im| (the rounding of the sum of squares and the
/// root is at most 2.5 u relative).
inline double modUp(double re, double im) {
  return std::sqrt(re * re + im * im) * (1.0 + 4.0 * kU);
}
inline double modDown(double re, double im) {
  return std::sqrt(re * re + im * im) * (1.0 - 4.0 * kU);
}

inline CBall cadd(CBall a, CBall b) {
  double re = a.re + b.re;
  double im = a.im + b.im;
  return CBall{re, im, up(a.r + b.r + kU * l1(re, im))};
}

inline CBall csub(CBall a, CBall b) {
  double re = a.re - b.re;
  double im = a.im - b.im;
  return CBall{re, im, up(a.r + b.r + kU * l1(re, im))};
}

/// Product of two balls: |ab - m_a m_b| <= |m_a| r_b + |m_b| r_a + r_a r_b,
/// plus at most sqrt(5) u |m_a| |m_b| for the midpoint's rounding (Brent,
/// Percival and Zimmermann 2007; 2 u when the compiler contracts to fused
/// multiply-adds).
inline CBall cmul(CBall a, CBall b) {
  double re = a.re * b.re - a.im * b.im;
  double im = a.re * b.im + a.im * b.re;
  double ma = modUp(a.re, a.im);
  double mb = modUp(b.re, b.im);
  double r = ma * b.r + mb * a.r + a.r * b.r + 2.25 * kU * ma * mb;
  return CBall{re, im, up(r)};
}

/// Product of a ball by a real ball.
inline CBall cmulReal(CBall a, RBall s) {
  double re = a.re * s.m;
  double im = a.im * s.m;
  double r = modUp(a.re, a.im) * s.r + std::fabs(s.m) * a.r + a.r * s.r +
             kU * l1(re, im);
  return CBall{re, im, up(r)};
}

/// Quotient a / b for r_b < |m_b|; undecided otherwise.
inline CBall cdiv(CBall a, CBall b) {
  double n = b.re * b.re + b.im * b.im;
  double re = (a.re * b.re + a.im * b.im) / n;
  double im = (a.im * b.re - a.re * b.im) / n;
  double mb = modDown(b.re, b.im);
  if (!(b.r < mb)) return CBall{re, im, kInf};
  double ma = modUp(a.re, a.im);
  double prop = (ma * b.r + modUp(b.re, b.im) * a.r) / (mb * (mb - b.r));
  // the numerator's products (sqrt(5) u |a| |b|), |b|^2 (2 u) and the two
  // divisions (u)
  double mid = 5.5 * kU * ma / mb;
  return CBall{re, im, up(prop + mid)};
}

/// exp(z): |e^z - e^m| <= |e^m| expm1(r); the midpoint's error is at most
/// (2 u_L + u) |e^m| (exp, then cos or sin, then the product).
inline CBall cexp(CBall z) {
  double ex = std::exp(z.re);
  double c = std::cos(z.im);
  double s = std::sin(z.im);
  double re = ex * c;
  double im = ex * s;
  double exu = ex * (1.0 + kLibRelErr);
  double r = exu * (std::expm1(z.r) * (1.0 + kLibRelErr) + 2.0 * kLibRelErr + kU);
  return CBall{re, im, up(r)};
}

/// expm1(z) for |z| below 1/2 by expm1(a) cos b - 2 sin^2(b/2) +
/// i e^a sin b (relative accuracy at small |z|), else e^z - 1.
inline CBall cexpm1(CBall z) {
  if (l1(z.re, z.im) >= 0.5) {
    CBall e = cexp(z);
    double re = e.re - 1.0;
    return CBall{re, e.im, up(e.r + kU * std::fabs(re))};
  }
  double em = std::expm1(z.re);
  double c = std::cos(z.im);
  double s = std::sin(z.im);
  double sh = std::sin(0.5 * z.im);
  double t1 = em * c;
  double t2 = 2.0 * sh * sh;
  double ex = 1.0 + em;
  double re = t1 - t2;
  double im = ex * s;
  double mid = (2.0 * kLibRelErr + 2.0 * kU) * (std::fabs(t1) + t2 + std::fabs(im)) +
               kU * (std::fabs(re) + std::fabs(im));
  double r = ex * (1.0 + kLibRelErr) * std::expm1(z.r) * (1.0 + kLibRelErr) + mid;
  return CBall{re, im, up(r)};
}

/// log1p(z), principal branch, when the disk of 1 + z misses the closed
/// negative real axis; undecided otherwise. |log1p(z) - log1p(m)| <=
/// -log1p(-r / |1 + m|) for r < |1 + m|.
inline CBall clog1p(CBall z) {
  double opre = 1.0 + z.re;
  if (!(opre > z.r || std::fabs(z.im) > z.r)) return CBall{0.0, 0.0, kInf};
  double m2 = z.re * (2.0 + z.re) + z.im * z.im;  // |1 + z|^2 - 1
  double em2 = 3.0 * kU * (std::fabs(z.re) * std::fabs(2.0 + z.re) + z.im * z.im);
  double l = std::log1p(m2);
  double re = 0.5 * l;
  double im = std::atan2(z.im, opre);
  double onep = 1.0 + m2 - em2;  // a lower bound on |1 + z|^2
  if (!(onep > 0.0)) return CBall{re, im, kInf};
  double mod = std::sqrt(onep) * (1.0 - 2.0 * kU);
  if (!(z.r < mod)) return CBall{re, im, kInf};
  // atan2's argument 1 + re carries u |1 + re|, moving the angle by at most
  // u |im z| / |1 + z|
  double mid = 0.5 * (em2 / onep + kLibRelErr * std::fabs(l)) + kU * std::fabs(re) +
               kLibRelErr * std::fabs(im) + 2.0 * kU * std::fabs(z.im) / mod;
  // -log1p(-e) <= e / (1 - e)
  double e = z.r / mod;
  double prop = e / (1.0 - e) * (1.0 + 4.0 * kU);
  return CBall{re, im, up(mid + prop)};
}

/// Principal square root of a ball whose midpoint and true value are known to
/// lie in the closed upper half plane (so the branch is continuous there):
/// |sqrt(W) - sqrt(m)| <= r / (sqrt(2) sqrt(|m| - r)) for r < |m|.
inline CBall csqrtUpper(CBall w) {
  double a = w.re;
  double c = w.im;
  double mod = std::sqrt(a * a + c * c);
  double p, q;
  if (a >= 0.0) {
    p = std::sqrt(0.5 * (mod + a));
    q = (p > 0.0) ? c / (2.0 * p) : 0.0;
  } else {
    q = std::sqrt(0.5 * (mod - a));
    if (c < 0.0) q = -q;
    p = std::fabs(c) / (2.0 * std::fabs(q));
  }
  double modLo = mod * (1.0 - 4.0 * kU);
  if (!(w.r < modLo)) return CBall{p, q, kInf};
  // both roots in the closed first quadrant: |sqrt(W) + sqrt(m)| >=
  // sqrt(2) sqrt(|m| - r)
  double prop = w.r / (0x1.6a09e667f3bccp+0 * std::sqrt(modLo - w.r) * (1.0 - 2.0 * kU));
  // |a + i c| to 1.5 u, p to 2.25 u, q = c / (2 p) to 3.25 u (or the
  // reverse for a < 0)
  double mid = 3.5 * kU * (std::fabs(p) + std::fabs(q));
  return CBall{p, q, up(prop + mid)};
}

/// log1p(z) - z. For |z| <= 1/4 by the series -z^2 sum_k (-z)^k / (k + 2),
/// with n terms where |z| <= 2^-e gives |z|^n <= 2^-54: the truncation is
/// below 2^-54 |z|^2 / 1.5 and Horner's rounding below 16 u |z|^2; the
/// input's radius passes through |d/dz| = |z| / |1 + z| <= m / (1 - m) on the
/// disk, m = |z| + r. Beyond, log1p(z) - z directly.
inline CBall clog1pMinus(CBall z) {
  const double m = modUp(z.re, z.im) + z.r;
  if (m <= 0.25) {
    if (m == 0.0) return CBall{0.0, 0.0, 0.0};
    int e;
    std::frexp(m, &e);  // m < 2^e, e <= -2
    const int n = (54 + (-e) - 1) / (-e);
    double pr = 1.0 / (n + 1), pi = 0.0;
    for (int k = n - 1; k >= 0; --k) {
      // p <- 1 / (k + 2) - z p
      double tr = z.re * pr - z.im * pi;
      double ti = z.re * pi + z.im * pr;
      pr = 1.0 / (k + 2) - tr;
      pi = -ti;
    }
    double z2r = z.re * z.re - z.im * z.im;
    double z2i = 2.0 * z.re * z.im;
    double re = -(z2r * pr - z2i * pi);
    double im = -(z2r * pi + z2i * pr);
    double m2 = m * m;
    double r = m2 * (16.0 * kU + 0x1p-54) + z.r * m / (1.0 - m) * (1.0 + 4.0 * kU);
    return CBall{re, im, up(r)};
  }
  CBall L = clog1p(z);
  double re = L.re - z.re;
  double im = L.im - z.im;
  return CBall{re, im, up(L.r + z.r + kU * (l1(L.re, L.im) + l1(re, im)))};
}

/// An upper bound on the log of a sum of exp(v_i), each v_i an upper bound.
class LogSumExpUpper {
 public:
  void add(double v) {
    if (!(v > -kInf)) return;
    if (std::isnan(v) || v == kInf) {
      bad_ = true;
      return;
    }
    if (v > max_) {
      sum_ = sum_ * std::exp(max_ - v) * (1.0 + kLibRelErr) + 1.0;
      max_ = v;
    } else {
      sum_ += std::exp(v - max_) * (1.0 + kLibRelErr);
    }
    ++n_;
  }
  double value() const {
    if (bad_) return kInf;
    if (n_ == 0) return -kInf;
    double s = sum_ * (1.0 + (2.0 * n_ + 4.0) * kU);
    double l = std::log(s);
    return max_ + l + kLibRelErr * std::fabs(l) + kU * std::fabs(max_) + 0x1p-1000;
  }

 private:
  double max_ = -kInf;
  double sum_ = 0.0;
  std::size_t n_ = 0;
  bool bad_ = false;
};

/// log(a) + log(b) style accumulation in the log domain with the sum rule:
/// returns an upper bound of v given a computed value and the sum of
/// magnitudes of the terms that formed it.
inline double logUpper(double value, double magnitude) {
  return value + 4.0 * kLibRelErr * magnitude + 0x1p-1000;
}

/// Rounds rho down to 24 significant bits: then rho y and y^2 are exact for
/// the nodes' y = k h (h with 4 significant bits, k < 2^20).
inline double roundDown24(double v) {
  int e;
  double f = std::frexp(v, &e);
  return std::ldexp(std::floor(std::ldexp(f, 24)), e - 24);
}

/// Rounds a positive step down to the ladder m 2^-j, m in 8..15, so every
/// node k h is exact.
inline double roundDownStep(double h) {
  int e;
  double f = std::frexp(h, &e);  // h = f 2^e, f in [0.5, 1)
  double m = std::floor(f * 16.0);  // 8..15
  return std::ldexp(m, e - 4);
}

// ---------------------------------------------------------------------------
// The settings a refinement level evaluates with, fixed before any
// evaluation: the edge lines (a_+ as fractions of rho, a_- as absolute
// distances, negative values meaning that multiple of rho), the cell counts
// per line, the depth of splitting, and the target of the discretization
// bound as a power of two relative to the estimated sum.

struct LevelSettings {
  double fracPlus[4];
  int numPlus;
  double aMinus[6];
  int numMinus;
  int cellsMin;
  int cellsMax;
  int maxDepth;
  int targetBits;
};

/// Level 0, one line each side at rho/2, 10 to 40 cells, no splitting;
/// level 1, two lines each side, 32 cells, two levels of splitting; level 2,
/// the full default. The q family's default uses fewer lines and cells.
inline constexpr int kNumLevels = 3;
inline constexpr LevelSettings kLevelsF[kNumLevels] = {
  {{0.5}, 1, {-0.5}, 1, 10, 40, 0, 56},
  {{0.5, 0.2}, 2, {-0.5, 1.0}, 2, 32, 32, 2, 60},
  {{0.5, 0.2, 0.05, 0.01}, 4, {-0.5, -1.0, -0.1, 0.25, 1.0, 4.0}, 6, 200, 3000, 3, 64},
};
inline constexpr LevelSettings kLevelsQ[kNumLevels] = {
  {{0.5}, 1, {-0.5}, 1, 10, 40, 0, 56},
  {{0.5, 0.2}, 2, {-0.5, 1.0}, 2, 32, 32, 2, 60},
  {{0.5, 0.2, 0.05}, 3, {-0.5, -0.1, 0.25, 1.0}, 4, 100, 600, 3, 64},
};

/// Counts of the work an evaluation did, for the timing harness and the
/// cost record.
struct EvalCounts {
  std::size_t nodes = 0;
  std::size_t truncationNodes = 0;
  std::size_t cells = 0;
  int tailPath = 0;  // 1 monotone tail, 2 node by node
  double h = 0.0;
};

/// The aggregated majorant's quantities on one line Re z = u (Lemma
/// "Aggregated majorant on a line", with J0 = {1..n0}, every j with
/// A_j(u) <= 2 u^2, so each J1 term has k_j > 0): kappa = max over J1 of k_j,
/// s = sum over J1 of a_j k_j (both upper bounds), and c an upper bound of
/// sum_j a_j log(A_j(rho) / A_j(u)).
struct LineAgg {
  double u = 0.0;
  int n0 = 0;
  double kappa = kInf;
  double s = kInf;
  double c = kInf;
  bool valid = false;
};

// ---------------------------------------------------------------------------
// T1(W) = sum_j 1 / (W + c'_j) = tanh(r) / (2 r) and T2(W) = sum_j
// 1 / (W + c'_j)^2 = (tanh r - r sech^2 r) / (4 r^3), r = sqrt(W),
// c'_j = pi^2 (j - 1/2)^2, continued to -pi^2/4 < W < 0; series in W for
// |W| <= 0.1. Returned values carry relative error below kTRelErr.

inline constexpr double kTRelErr = 0x1p-40;

namespace detail {
// P_m = sum_j c'_j^-m, m = 1..15.
inline constexpr double kPowerSum[15] = {
  0x1.0000000000000p-1, 0x1.5555555555555p-3, 0x1.1111111111111p-4,
  0x1.ba1ba1ba1ba1cp-6, 0x1.664f4882c10fap-7, 0x1.226e355e6c23dp-8,
  0x1.d6d3d0e157de0p-10, 0x1.7da36452b75e3p-11, 0x1.3558248036744p-12,
  0x1.f57d7734d1664p-14, 0x1.967e18afcafadp-15, 0x1.497d8eea25259p-16,
  0x1.0b132d39a6050p-17, 0x1.b0f72d3ee24e9p-19, 0x1.5ef2da474e5b7p-20,
};
}  // namespace detail

inline bool sumsT(double W, double& t1, double& t2) {
  if (std::fabs(W) <= 0.1) {
    // T1 = sum_n (-W)^n P_{n+1}, T2 = sum_n (n + 1) (-W)^n P_{n+2}; the
    // terms decrease by at least |W| / c'_1 <= 0.0406, so 13 terms leave a
    // relative error below 2^-60.
    double a = 0.0, b = 0.0;
    for (int n = 12; n >= 0; --n) {
      a = detail::kPowerSum[n] - W * a;
      b = (n + 1) * detail::kPowerSum[n + 1] - W * b;
    }
    t1 = a;
    t2 = b;
    return true;
  }
  if (W > 0.0) {
    double r = std::sqrt(W);
    double th = std::tanh(r);
    double ch = std::cosh(r);
    t1 = th / (2.0 * r);
    t2 = (th - r / (ch * ch)) / (4.0 * r * W);
    return std::isfinite(t1) && std::isfinite(t2);
  }
  double v = std::sqrt(-W);
  // tan's relative sensitivity to v's rounding is v / (sin v cos v)
  if (!(0.5 * kPi - v > 0x1p-20)) return false;
  double tn = std::tan(v);
  double cs = std::cos(v);
  t1 = tn / (2.0 * v);
  t2 = (v / (cs * cs) - tn) / (4.0 * v * (-W));
  return true;
}

/// log cosh sqrt(W) for W >= 0 (lambda(sqrt W)) and log cos sqrt(-W) for
/// -pi^2/4 < W < 0, with an absolute error bound returned in err.
inline double logCoshSqrt(double W, double& err) {
  if (W >= 0.0) {
    double w = std::sqrt(W);
    double e = std::exp(-2.0 * w);
    double l = std::log1p(e);
    double v = w - kLog2 + l;
    err = 8.0 * kLibRelErr * (w + kLog2 + l) + 0x1p-1000;
    return v;
  }
  double v = std::sqrt(-W);
  double c = std::cos(v);
  double l = std::log(c);
  // cos at v near pi/2 loses only relative accuracy kLibRelErr plus the
  // root's rounding times tan(v) v
  err = 4.0 * kLibRelErr * (std::fabs(l) + 1.0 + std::tan(v) * v) + 0x1p-1000;
  return l;
}

// ---------------------------------------------------------------------------
// The target f = PG(b, psi): its Laplace transform is
// L(t) = cosh^b(A) / cosh^b sqrt(A^2 + t / 2), A = |psi| / 2. In the frame of
// the contour, W = A^2 + t / 2 = W0 + Delta, and a node's log ratio is
// -b D with D = log cosh w - log cosh w0 in the cancellation-free form
// log1p(expm1(delta) (1 - e^{-(w + w0)}) / (1 + e^{-2 w0})),
// delta = Delta / (w + w0).

class TargetDensity {
 public:
  TargetDensity() = default;
  TargetDensity(double b, double psi) { setParameters(b, psi); }

  void setParameters(double b, double psi) {
    b_ = b;
    a_ = 0.5 * std::fabs(psi);
    prepared_ = false;
  }

  double shape() const { return b_; }
  double halfPsi() const { return a_; }

  /// Chooses the saddle and contour for point x and forms the lead.
  void prepare(double x);

  double point() const { return x_; }
  double rho() const { return rho_; }
  double w0Squared() const { return W0_; }
  double tiltedVariance() const { return vt_; }
  Enclosure lead() const { return lead_; }

  /// The node exponent 2 x Delta + log(L(t) / L(t0)) at y > 0.
  CBall nodeExponent(double y) const;
  /// The direct form's D, for the branch and as a fallback.
  CBall directD(CBall w) const;

  void lineSetup(double u, LineAgg& agg) const;
  double cellLogMajor(double u, double ya, double yb, const LineAgg& agg) const;
  double tailLogMajor(double u, double Y) const;
  /// An upper bound on the log of the J0 terms of the aggregated majorant on
  /// the contour line at y^2 = Y (each term at its value, not its minimum).
  double j0TermsAt(double Y, const LineAgg& agg) const;
  /// An upper bound on sum_{k h >= y0} |G(k h)| on the contour line, for
  /// y0 >= (-rho + sqrt(rho^2 + 1 / x)) / 2: the monotone tail.
  double farTail(double y0, double h) const;
  /// log |G| at u + i y on its midpoint only (a splitting heuristic).
  double pointLog(double u, double y) const;

 private:
  double cellLogMajorDirect(double u, double ya, double yb) const;
  double sPrime() const { return rho2_ - W0_; }

  double b_ = 1.0;
  double a_ = 0.0;
  double x_ = 0.0;
  double W0_ = 0.0;   // exact: +w0^2 or -v^2
  double w0_ = 0.0;   // sqrt(W0) when W0 >= 0
  double v_ = 0.0;    // sqrt(-W0) when W0 < 0
  double rho_ = 0.0;  // 24 bits
  double rho2_ = 0.0; // exact
  double vt_ = 0.0;
  CBall invDen_{1.0, 0.0, 0.0};  // 1 / (1 + e^{-2 w0})
  RBall sigma0_{0.0, 0.0};       // 1 / (1 + e^{2 w0}), for w0 >= 1
  RBall c1_{0.0, 0.0};           // 4 x w0 - b, for w0 >= 1
  bool lean_ = false;
  CBall lambda0_{0.0, 0.0, 0.0}; // log cosh w0, real
  double logCosh0_ = 0.0;        // an upper bound of log cosh w0
  double logCosh0Err_ = 0.0;
  Enclosure lead_;
  bool prepared_ = false;
};

// ---------------------------------------------------------------------------
// A sum of independent gammas, X = sum_j Ga(a_j, mu_j) with rates ascending
// (the proposal: a_j = b at the K rates r_k and alpha at beta). Its Laplace
// transform is prod_j (1 + t / mu_j)^-a_j; a node's log ratio is
// -sum_j a_j log1p(2 Delta / (mu_j + t0)).

class GammaSumDensity {
 public:
  GammaSumDensity() = default;

  /// The shapes and rates are copied; rates must be positive and ascending.
  void setParameters(const double* shapes, const double* rates, std::size_t n) {
    shape_.assign(shapes, shapes + n);
    rate_.assign(rates, rates + n);
    inv_.resize(n);
    invRad_.resize(n);
    logA_.resize(n);
    tau2_.resize(n);
    tau2Err_.resize(n);
    prepared_ = false;
  }

  std::size_t numTerms() const { return shape_.size(); }

  void prepare(double x);

  double point() const { return x_; }
  double rho() const { return rho_; }
  double saddle() const { return t0_; }
  double tiltedVariance() const { return vt_; }
  Enclosure lead() const { return lead_; }

  CBall nodeExponent(double y) const;

  void lineSetup(double u, LineAgg& agg) const;
  double cellLogMajor(double u, double ya, double yb, const LineAgg& agg) const;
  double tailLogMajor(double u, double Y) const;
  double j0TermsAt(double Y, const LineAgg& agg) const;
  double pointLog(double u, double y) const;
  /// An upper bound on sum_{k h >= y0} |G(k h)| on the contour line, for
  /// y0 >= (-rho + sqrt(rho^2 + 1 / x)) / 2, by blocks [Y, 2Y] in y^2: on
  /// each, a bound linear in Y, -2 x Y plus, per factor, (4/7) a_j k_j Y while
  /// k_j Y_hi <= 1/8 (-log(1 - c) / 2 <= c / (2 (1 - c))), else its minimum
  /// over the block, taken at the block's worse end and times its node count
  /// and largest |z|; past
  /// Y >= 2 max_j (A_j - 2 rho^2) every factor is at most 1 and the monotone
  /// tail closes the sum.
  double farTail(double y0, double h) const;

  /// A lower bound on the minimum over Y' in [ylo, yhi] (yhi infinite
  /// allowed) of (A - Y')^2 + 4 u^2 Y', A = u^2 + tau2, given tau2 within
  /// tau2Err: the convex function's ends and, when its vertex
  /// V = tau2 - u^2 may lie inside, its vertex value 4 u^2 tau2.
  static double termMin(double u2, double ylo, double yhi, double tau2, double tau2Err);

 private:
  /// An upper bound on sum_j a_j [log A_j(rho) - log(min |z^2 + tau_j'^2|) ]
  /// over Re z = u, y^2 in [ylo, yhi], by each term's minimum.
  double perTermLogMajor(double u2, double ylo, double yhi) const;

  std::vector<double> shape_;
  std::vector<double> rate_;
  std::vector<double> inv_;      // 1 / (mu_j + t0)
  std::vector<double> invRad_;   // its absolute radius
  std::vector<double> logA_;     // an upper bound of log((mu_j + t0) / 2)
  std::vector<double> tau2_;     // (mu_j + t0) / 2 - rho^2
  std::vector<double> tau2Err_;  // its absolute error
  std::size_t runEnd_ = 0;       // terms [0, runEnd_) share shape_[0]
  RBall c0_{0.0, 0.0};           // x - sum_j a_j / (mu_j + t0)
  double x_ = 0.0;
  double t0_ = 0.0;
  double rho_ = 0.0;
  double rho2_ = 0.0;
  double vt_ = 0.0;
  Enclosure lead_;
  bool prepared_ = false;
};

// ---------------------------------------------------------------------------
// TargetDensity implementation.

inline void TargetDensity::prepare(double x) {
  x_ = x;
  // Saddle: x = (b / 4) g(W0), g(W) = tanh(sqrt W) / sqrt W, decreasing from
  // infinity at -pi^2/4 to 0; any W0 is valid (double precision suffices),
  // and it is then rounded to 26 bits of sqrt so W0 = +-w0^2 exactly.
  double c = 4.0 * x / b_;
  if (c < 1.0) {
    // tanh(w) / w = c on w in (0, 1/c]
    double lo = 0.0, hi = 1.0 / c;
    double w = (c > 0.99) ? std::sqrt(3.0 * (1.0 - c)) : std::min(hi, 1.0 / c);
    for (int it = 0; it < 100; ++it) {
      double th = std::tanh(w);
      double f = th - c * w;  // decreasing through the root
      if (f > 0.0) lo = w; else hi = w;
      double ch = std::cosh(w);
      double fp = 1.0 / (ch * ch) - c;
      double wn = w - f / fp;
      if (!(wn > lo && wn < hi)) wn = 0.5 * (lo + hi);
      if (std::fabs(wn - w) <= 1e-15 * w) { w = wn; break; }
      w = wn;
    }
    int e;
    double f = std::frexp(w, &e);
    w0_ = std::ldexp(std::nearbyint(std::ldexp(f, 26)), e - 26);
    v_ = 0.0;
    W0_ = w0_ * w0_;
  } else if (c > 1.0) {
    // tan(v) / v = c on v in (0, pi/2)
    double lo = 0.0, hi = 0.5 * kPi;
    double v = (c < 1.01) ? std::sqrt(3.0 * (c - 1.0)) : 0.5 * kPi * (1.0 - 1.0 / c);
    for (int it = 0; it < 100; ++it) {
      double f = std::tan(v) - c * v;  // increasing through the root
      if (f < 0.0) lo = v; else hi = v;
      double cs = std::cos(v);
      double fp = 1.0 / (cs * cs) - c;
      double vn = v - f / fp;
      if (!(vn > lo && vn < hi)) vn = 0.5 * (lo + hi);
      if (std::fabs(vn - v) <= 1e-15 * v) { v = vn; break; }
      v = vn;
    }
    int e;
    double f = std::frexp(v, &e);
    v_ = std::ldexp(std::floor(std::ldexp(f, 26)), e - 26);
    if (!(v_ < 0x1.921fb5p+0)) v_ = 0x1.921fb5p+0;
    w0_ = 0.0;
    W0_ = -v_ * v_;
  } else {
    w0_ = v_ = W0_ = 0.0;
  }

  // Contour: T1 = W0 + pi^2/4 (> 0); rho^2 = max(W0, T1 / 2) < T1.
  double T1 = W0_ + 0.25 * kPiSquared;
  double target = std::max(W0_, 0.5 * T1);
  rho_ = roundDown24(std::sqrt(target) * (1.0 - 4.0 * kU));
  rho2_ = rho_ * rho_;

  double t1, t2;
  sumsT(W0_, t1, t2);
  vt_ = 0.25 * b_ * t2;

  // 1 / (1 + e^{-2 w0}) and log cosh w0
  if (W0_ >= 0.0) {
    double e = std::exp(-2.0 * w0_);
    double d = 1.0 / (1.0 + e);
    invDen_ = CBall{d, 0.0, up(d * (2.0 * kU + kLibRelErr * e))};
    lean_ = w0_ >= 1.0;
    if (lean_) {
      // e / (1 + e): exp, sum and quotient; an underflowed e leaves 2^-1074
      double sg = e * d;
      sigma0_ = RBall{sg, up(sg * (kLibRelErr + 3.0 * kU) + 0x1p-1074)};
      double c1 = std::fma(4.0 * x, w0_, -b_);
      c1_ = RBall{c1, up(kU * std::fabs(c1))};
    }
  } else {
    double t = std::tan(v_);
    invDen_ = CBall{0.5, 0.5 * t, up(0.5 * std::fabs(t) * (kLibRelErr + kU))};
    lean_ = false;
  }
  double lerr;
  double lc = logCoshSqrt(W0_, lerr);
  lambda0_ = CBall{lc, 0.0, lerr};
  logCosh0_ = lc;
  logCosh0Err_ = lerr;

  // The lead x t0 - b [lambda(w0) - lambda(A)] in double-double, t0 =
  // 2 (W0 - A^2) exact as a double-double.
  dd::DD A2 = dd::twoProd(a_, a_);
  dd::DD t0 = dd::ldexp(dd::sub(dd::fromDouble(W0_), A2), 1);
  dd::DD xt0 = dd::mul(t0, x);
  auto lambdaDD = [](double w, double& mag) {
    // w - log 2 + log1p(e^{-2w}); past w = 325 the last term is below e^-650
    dd::DD v = dd::sub(dd::fromDouble(w), dd::kLn2);
    double l = 0.0;
    if (w <= 325.0) {
      dd::DD e = dd::exp(dd::fromDouble(-2.0 * w));
      dd::DD lp = dd::log1p(e);
      v = dd::add(v, lp);
      l = lp.hi;
    }
    mag = w + kLog2 + l;
    return v;
  };
  double magA, mag0;
  dd::DD lamA = lambdaDD(a_, magA);
  dd::DD lam0;
  double extra = 0.0;
  if (W0_ >= 0.0) {
    lam0 = lambdaDD(w0_, mag0);
  } else {
    // log cos v = log sin(pi/2 - v), pi/2 in double-double; the argument's
    // error (pi/2's rounding, the subtraction) passes through cot(arg)
    dd::DD arg = dd::sub(dd::kPiOver2, dd::fromDouble(v_));
    dd::DD s = dd::sin(arg);
    lam0 = dd::log(s);
    mag0 = std::fabs(lam0.hi) + 1.0;
    extra = (0x1p-104 + dd::kU2 * arg.hi) / arg.hi;
  }
  dd::DD diff = dd::sub(lam0, lamA);
  dd::DD bdiff = dd::mul(diff, b_);
  dd::DD lead = dd::sub(xt0, bdiff);
  double ulf = dd::kLogRelErr;  // the largest of the functions' bounds
  double mag = std::fabs(xt0.hi) + b_ * (mag0 + magA + 2.0);
  double rad = (16.0 * dd::kU2 + 3.0 * ulf) * mag + b_ * extra * 2.0 +
               (W0_ >= 0.0 ? 0.0 : 0.0);
  // The lead plus log(2 / pi): the double midpoint's rounding is charged.
  double mid = lead.hi + (lead.lo + kLogTwoOverPi);
  rad += kU * std::fabs(mid) + 0x1p-55 + kU * std::fabs(lead.lo + kLogTwoOverPi);
  lead_ = Enclosure{mid, up(rad)};
  prepared_ = true;
}

inline CBall TargetDensity::directD(CBall w) const {
  // lambda(w) - lambda(w0), lambda(w) = w - log 2 + log1p(e^{-2w})
  CBall m2w{-2.0 * w.re, -2.0 * w.im, 2.0 * w.r};
  CBall e = cexp(m2w);
  CBall l = clog1p(e);
  CBall lw{w.re - kLog2 + l.re, w.im + l.im,
           up(w.r + l.r + kU * (std::fabs(w.re) + kLog2 + std::fabs(l.re)) * 2.0 +
              kU * std::fabs(w.im + l.im) + 0x1p-54)};
  return csub(lw, lambda0_);
}

inline CBall TargetDensity::nodeExponent(double y) const {
  // Delta = -y^2 + 2 i rho y, exact; W = W0 + Delta, real part rounded.
  double dre = -y * y;
  double dim = 2.0 * rho_ * y;
  double wre = W0_ + dre;
  CBall W{wre, dim, kU * std::fabs(wre)};
  CBall w = csqrtUpper(W);
  CBall s = (W0_ >= 0.0) ? CBall{w.re + w0_, w.im, 0.0} : CBall{w.re, w.im + v_, 0.0};
  s.r = up(w.r + kU * l1(s.re, s.im));
  CBall delta = cdiv(CBall{dre, dim, 0.0}, s);
  if (lean_) {
    // D = delta + log1p(sigma0 expm1(-2 delta)), delta = w - w0, and
    // 2 x Delta - b delta = delta (2 x delta + 4 x w0 - b), so
    // E = 2 x delta^2 + c1 delta - b log1p(sigma0 expm1(-2 delta)); the
    // principal log1p is the continuous branch (1 + e^{-2w} has positive real
    // part on Re w > 0).
    CBall em = cexpm1(CBall{-2.0 * delta.re, -2.0 * delta.im, 2.0 * delta.r});
    CBall L = clog1p(cmulReal(em, sigma0_));
    CBall d2 = cmul(delta, delta);
    CBall t1 = cmulReal(d2, RBall{2.0 * x_, 0.0});
    CBall t2 = cmulReal(delta, c1_);
    CBall t3 = cmulReal(L, RBall{b_, 0.0});
    return csub(cadd(t1, t2), t3);
  }
  CBall em = cexpm1(delta);
  CBall ems = cexpm1(CBall{-s.re, -s.im, s.r});
  CBall om{-ems.re, -ems.im, ems.r};  // 1 - e^{-s}
  CBall arg = cmul(cmul(em, om), invDen_);
  CBall D = clog1p(arg);
  bool useDirect = !(std::isfinite(D.r) && D.r < 0.25);
  if (!useDirect) {
    // The principal log1p is the continuous branch while Im w <= pi/2 (then
    // Im D lies in (Im w - pi/2, Im w + pi/2)); otherwise fix the branch
    // integer from the direct form.
    if (w.im + w.r >= 0.5 * kPi - 0.01) {
      CBall Dd = directD(w);
      if (!std::isfinite(Dd.r) || Dd.r > 1.0) {
        D.r = kInf;
      } else {
        double k = std::nearbyint((Dd.im - D.im) / (2.0 * kPi));
        double res = Dd.im - D.im - k * 2.0 * kPi;
        if (std::fabs(res) > 1.0) {
          D.r = kInf;
        } else if (k != 0.0) {
          double shift = k * 2.0 * kPi;
          D.im += shift;
          D.r = up(D.r + std::fabs(k) * 2.0 * 0x1p-52 + kU * std::fabs(D.im));
        }
      }
    }
  } else {
    D = directD(w);
  }
  // E = 2 x Delta - b D
  double ere2 = 2.0 * x_ * dre;
  double eim2 = 2.0 * x_ * dim;
  double bre = b_ * D.re;
  double bim = b_ * D.im;
  double ere = ere2 - bre;
  double eim = eim2 - bim;
  double r = b_ * D.r + kU * (l1(ere2, eim2) + l1(bre, bim) + l1(ere, eim));
  return CBall{ere, eim, up(r)};
}

inline void TargetDensity::lineSetup(double u, LineAgg& agg) const {
  agg = LineAgg{};
  agg.u = u;
  double u2 = u * u;
  // W_u = W0 + u^2 - rho^2 and A_j(u) = c'_j + W_u
  double Wu = W0_ + (u2 - rho2_);
  double eWu = 4.0 * kU * (std::fabs(W0_) + u2 + rho2_);
  if (!(Wu - eWu > -0.25 * kPiSquared)) return;
  auto Aj = [&](int j) { double h = j - 0.5; return kPiSquared * h * h + Wu; };
  auto kj = [&](double A) { return 2.0 * (A - 2.0 * u2) / (A * A); };
  double eA = eWu + 8.0 * kU * (std::fabs(Wu) + kPiSquared);
  int n0 = 0;
  double sumJ0 = 0.0, magJ0 = 0.0, errJ0 = 0.0;
  while (Aj(n0 + 1) - eA * (n0 + 1) * (n0 + 1) <= 2.0 * u2 * (1.0 + 4.0 * kU)) {
    double A = Aj(n0 + 1);
    double eAj = eA * (n0 + 1) * (n0 + 1);
    double Alo = A - eAj;
    if (!(Alo > 0.0)) return;
    double k = kj(A);
    sumJ0 += b_ * k;
    magJ0 += b_ * std::fabs(k);
    // |dk / dA| = |8 u^2 - 2 A| / A^3
    errJ0 += b_ * eAj * (8.0 * u2 + 2.0 * A + 2.0 * eAj) / (Alo * Alo * Alo);
    ++n0;
    if (n0 > 100000) return;
  }
  // kappa: g(A) = 2 (A - 2u^2) / A^2 rises on (2u^2, 4u^2) and falls past
  // 4u^2; take the max over J1's terms up to the first past 4u^2.
  double kappa = 0.0;
  for (int j = n0 + 1;; ++j) {
    double A = Aj(j);
    double Alo = A - eA * j * j;
    double Ahi = A + eA * j * j;
    // an upper bound of g on [Alo, Ahi]
    double g = 2.0 * (Ahi - 2.0 * u2) / (Alo * Alo) * (1.0 + 8.0 * kU);
    if (Alo <= 4.0 * u2 && Ahi >= 4.0 * u2) g = std::max(g, 0.25 / u2 * (1.0 + 4.0 * kU));
    kappa = std::max(kappa, g);
    if (Alo > 4.0 * u2) break;
    if (j > n0 + 100000) return;
  }
  double t1, t2;
  if (!sumsT(Wu, t1, t2)) return;
  // sum_J b k_j = 2b [T1 - 2 u^2 T2]; S = that minus J0's
  double tot = 2.0 * b_ * (t1 - 2.0 * u2 * t2);
  double totErr = 2.0 * b_ * (kTRelErr * (t1 + 2.0 * u2 * t2) +
                              (eWu + eA) * (t2 + 4.0 * u2 * t2 / (Wu + 0.25 * kPiSquared)));
  double S = tot - sumJ0;
  double SErr = totErr + 4.0 * kU * (std::fabs(tot) + magJ0) + errJ0;
  // C = b [log cosh sqrt(W0) - log cosh sqrt(W_u)]
  double eu;
  double lu = logCoshSqrt(Wu, eu);
  // W_u's own error moves log cosh sqrt by at most eWu T1-like: bound the
  // derivative tanh(r) / (2r) <= T1(W) <= T1 at Wu - eWu
  double d1, d2;
  if (!sumsT(Wu - eWu, d1, d2)) return;
  double C = b_ * (logCosh0_ - lu);
  double CErr = b_ * (logCosh0Err_ + eu + eWu * d1 * (1.0 + kTRelErr)) +
                4.0 * kU * b_ * (std::fabs(logCosh0_) + std::fabs(lu));
  agg.n0 = n0;
  agg.kappa = kappa;
  agg.s = std::max(0.0, S + SErr);
  agg.c = C + CErr;
  agg.valid = std::isfinite(agg.s) && std::isfinite(agg.c) && kappa > 0.0;
}

namespace detail {

/// |cosh w|^2 lower bound pieces for the target's line majorant: w =
/// sqrt(W(y)), W(y) = W0 + u^2 - rho^2 - y^2 + 2 i u y. Returns p - r and
/// q -+ r with r the root's ball radius.
struct PQ {
  double plo, qlo, qhi, reW, reWErr;
  bool ok;
};

/// log sinh(p) for p > 0, from below up to the final log's own rounding,
/// which the caller charges; past p = 20 as p - log 2 + log1p(-e^{-2p}), so
/// it holds beyond sinh's range of a double.
inline double logSinhLower(double p) {
  if (p <= 20.0) return std::log(std::sinh(p) * (1.0 - 2.0 * kLibRelErr));
  double v = (p - kLog2) + std::log1p(-std::exp(-2.0 * p));
  return v - 4.0 * kU * (p + 1.0);
}

inline PQ pqAt(double W0, double rho2, double u, double y) {
  double u2 = u * u;
  double reW = W0 + (u2 - rho2) - y * y;
  double reWErr = 4.0 * kU * (std::fabs(W0) + u2 + rho2 + y * y);
  double imW = 2.0 * u * y;
  if (y == 0.0 && reW - reWErr < 0.0) {
    // the root on the cut's side: p = 0, q = sqrt(-reW)
    double q = std::sqrt(std::max(0.0, -reW));
    double qe = std::sqrt(reWErr) + 4.0 * kU * q;
    return PQ{0.0, std::max(0.0, q - qe), q + qe, reW, reWErr, true};
  }
  CBall w = csqrtUpper(CBall{reW, imW, reWErr + kU * imW});
  if (!std::isfinite(w.r)) return PQ{0, 0, 0, reW, reWErr, false};
  return PQ{std::max(0.0, w.re - w.r), std::max(0.0, w.im - w.r), w.im + w.r, reW, reWErr, true};
}

}  // namespace detail

inline double TargetDensity::cellLogMajorDirect(double u, double ya, double yb) const {
  detail::PQ a = detail::pqAt(W0_, rho2_, u, ya);
  detail::PQ c = detail::pqAt(W0_, rho2_, u, yb);
  if (!a.ok || !c.ok) return kInf;
  double plo = std::min(a.plo, c.plo);
  double qlo = a.qlo;
  double qhi = c.qhi;
  double c2 = 0.0;
  if (qhi - qlo < kPi) {
    // a zero of cos, (j + 1/2) pi, inside [qlo, qhi] puts the minimum at 0
    double j = std::ceil(qlo / kPi - 0.5 - 1e-12);
    bool hit = (j + 0.5) * kPi <= qhi * (1.0 + 1e-12) + 1e-12;
    if (!hit) {
      double ca = std::cos(qlo), cb = std::cos(qhi);
      c2 = std::min(ca * ca, cb * cb) * (1.0 - 4.0 * kLibRelErr) - 1e-300;
      if (c2 < 0.0) c2 = 0.0;
    }
  }
  double m2c = 1.0 + c.reW - c.reWErr;
  double lm;
  if (plo > 20.0) {
    // |cosh w|^2 >= sinh^2 p, in logs: sinh^2 p leaves a double's range
    lm = 2.0 * detail::logSinhLower(plo);
    if (m2c > 0.0) lm = std::max(lm, std::log(m2c * (1.0 - 2.0 * kU)));
  } else {
    double sp = std::sinh(plo);
    double m2 = sp * sp * (1.0 - 4.0 * kLibRelErr) + c2;
    m2 = std::max(m2, m2c * (1.0 - 2.0 * kU));
    if (!(m2 > 0.0)) return kInf;
    lm = std::log(m2);
  }
  double v = b_ * (logCosh0_ - 0.5 * lm);
  return v + b_ * (logCosh0Err_ + kLibRelErr * std::fabs(lm)) +
         4.0 * kU * b_ * (std::fabs(logCosh0_) + 0.5 * std::fabs(lm)) + 0x1p-1000;
}

inline double TargetDensity::cellLogMajor(double u, double ya, double yb,
                                          const LineAgg& agg) const {
  double direct = cellLogMajorDirect(u, ya, yb);
  if (!agg.valid) return direct;
  double Yb = yb * yb;
  if (!(agg.kappa * Yb <= 0.5)) return direct;
  // Aggregated: C - sum_J0 (b/2) log(|z^2 + tau'^2|^2 / A_j^2) + (S/2) Y / (1 - kappa Y)
  double u2 = u * u;
  double Wu = W0_ + (u2 - rho2_);
  double tauBase = W0_ - rho2_;  // tau_j'^2 = c'_j + W0 - rho^2
  double eTau = 4.0 * kU * (std::fabs(W0_) + rho2_);
  double v = agg.c + 0.5 * agg.s * Yb / (1.0 - agg.kappa * Yb) * (1.0 + 8.0 * kU);
  double Ya = ya * ya;
  for (int j = 1; j <= agg.n0; ++j) {
    double h = j - 0.5;
    double cj = kPiSquared * h * h;
    double tau2 = cj + tauBase;
    double A = cj + Wu;
    double e = eTau + 8.0 * kU * (cj + std::fabs(tauBase));
    double m = GammaSumDensity::termMin(u2, Ya, Yb, tau2, e);
    if (!(m > 0.0)) return direct;
    double Ahi = A + e + 4.0 * kU * std::fabs(A);
    double lv = 2.0 * std::log(Ahi) - std::log(m);
    v += 0.5 * b_ * lv + b_ * kLibRelErr * (2.0 * std::fabs(std::log(Ahi)) + std::fabs(std::log(m)));
  }
  return std::min(direct, v + 4.0 * kU * std::fabs(v) + 0x1p-1000);
}

inline double TargetDensity::tailLogMajor(double u, double Y) const {
  // y >= Y on Re z = u: |cosh w| >= sinh p, p >= u Y / sqrt(Y^2 + s') when
  // s' > 0, p >= u when s' <= 0
  double sp = sPrime();
  double spUp = sp + 4.0 * kU * (rho2_ + std::fabs(W0_));
  double plo;
  if (spUp <= 0.0) {
    plo = u;
  } else {
    plo = u * Y / std::sqrt(Y * Y + spUp) * (1.0 - 4.0 * kU);
  }
  if (!(plo > 0.0)) return kInf;
  double ls = detail::logSinhLower(plo);
  double v = b_ * (logCosh0_ - ls);
  return v + b_ * (logCosh0Err_ + 2.0 * kLibRelErr * std::fabs(ls)) +
         4.0 * kU * b_ * (std::fabs(logCosh0_) + std::fabs(ls)) + 0x1p-1000;
}

inline double TargetDensity::j0TermsAt(double Y, const LineAgg& agg) const {
  // On the contour line every J0 term has A_j <= 2 rho^2 and so is at most 0;
  // its value -(b/2) log(((A - Y)^2 + 4 rho^2 Y) / A^2) bounded above.
  double tauBase = W0_ - rho2_;
  double eTau = 4.0 * kU * (std::fabs(W0_) + rho2_);
  double v = 0.0;
  for (int j = 1; j <= agg.n0; ++j) {
    double h = j - 0.5;
    double cj = kPiSquared * h * h;
    double tau2 = cj + tauBase;
    double A = cj + W0_;
    double e = eTau + 8.0 * kU * (cj + std::fabs(tauBase));
    double m = GammaSumDensity::termMin(rho2_, Y, Y, tau2, e);
    if (!(m > 0.0)) return 0.0;
    double Ahi = A + e + 4.0 * kU * std::fabs(A);
    double lv = 2.0 * std::log(Ahi) - std::log(m);
    double t = 0.5 * b_ * lv + b_ * kLibRelErr * (2.0 * std::fabs(std::log(Ahi)) + std::fabs(std::log(m)));
    v += std::min(0.0, t);
  }
  return v;
}

inline double TargetDensity::pointLog(double u, double y) const {
  // midpoint only: Re[2 x Delta_u - b D] + log |z|, Delta_u = u^2 - rho^2 - y^2 + 2 i u y
  double dre = u * u - rho2_ - y * y;
  double dim = 2.0 * u * y;
  CBall W{W0_ + dre, dim, 0.0};
  CBall w = csqrtUpper(W);
  CBall s = (W0_ >= 0.0) ? CBall{w.re + w0_, w.im, 0.0} : CBall{w.re, w.im + v_, 0.0};
  CBall delta = cdiv(CBall{dre, dim, 0.0}, CBall{s.re, s.im, 0.0});
  CBall em = cexpm1(CBall{delta.re, delta.im, 0.0});
  CBall ems = cexpm1(CBall{-s.re, -s.im, 0.0});
  CBall arg = cmul(cmul(em, CBall{-ems.re, -ems.im, 0.0}), CBall{invDen_.re, invDen_.im, 0.0});
  CBall D = clog1p(CBall{arg.re, arg.im, 0.0});
  if (!std::isfinite(D.r)) D = directD(CBall{w.re, w.im, 0.0});
  return 2.0 * x_ * dre - b_ * D.re + 0.5 * std::log(u * u + y * y);
}


inline double TargetDensity::farTail(double y0, double h) const {
  double lt = tailLogMajor(rho_, y0);
  if (!std::isfinite(lt)) return kInf;
  double g = -2.0 * x_ * y0 * y0;
  double e = lt + g + 4.0 * kLibRelErr * (std::fabs(lt) + std::fabs(g)) + 0x1p-1000;
  double f = (rho_ + y0) + (rho_ / y0 + 1.0) / (4.0 * x_ * h);
  return std::exp(e) * f * (1.0 + 2.0 * kLibRelErr + 8.0 * kU);
}

// ---------------------------------------------------------------------------
// GammaSumDensity implementation.

inline double GammaSumDensity::farTail(double y0, double h) const {
  const std::size_t n = shape_.size();
  double Amax = 0.0;
  for (std::size_t j = 0; j < n; ++j) Amax = std::max(Amax, tau2_[j] + tau2Err_[j] + rho2_);
  const double Yend = 2.0 * (Amax - 2.0 * rho2_) * (1.0 + 8.0 * kU);
  double Ya = y0 * y0 * (1.0 - 2.0 * kU);
  LogSumExpUpper acc;
  for (int blk = 0; blk < 2100; ++blk) {
    if (Ya >= Yend) {
      double ya = std::sqrt(Ya) * (1.0 - 2.0 * kU);
      double g = -2.0 * x_ * Ya;
      double lf = std::log((rho_ + ya) + (rho_ / ya + 1.0) / (4.0 * x_ * h));
      acc.add(g + lf + 4.0 * kLibRelErr * (std::fabs(g) + std::fabs(lf)) + 0x1p-1000);
      double v = acc.value();
      return std::exp(v) * (1.0 + 2.0 * kLibRelErr);
    }
    double Yb = std::min(2.0 * Ya, Yend);
    if (!(Yb > Ya)) Yb = Yend;
    // B(Y) <= Csw + slope Y on the block: switched factors (k_j Y_hi > 1/8)
    // at their minimum over the block, the others by
    // -log(1 - c) / 2 <= c / (2 (1 - 1/8)), c = k_j Y
    double Csw = 0.0, slope = 0.0, mag = 0.0;
    for (std::size_t j = 0; j < n; ++j) {
      double A = tau2_[j] + rho2_;
      double eA = tau2Err_[j] + kU * A;
      if (A + eA <= 2.0 * rho2_) continue;  // k_j <= 0: the factor is at most 1
      double Alo = A - eA;
      double kUp = 2.0 * (A + eA - 2.0 * rho2_) / (Alo * Alo) * (1.0 + 8.0 * kU);
      if (kUp * Yb <= 0.125) {
        slope += shape_[j] * kUp * (4.0 / 7.0) * (1.0 + 4.0 * kU);
      } else {
        double m = termMin(rho2_, Ya, Yb, tau2_[j], tau2Err_[j]);
        if (!(m > 0.0)) return kInf;
        double lm = std::log(m);
        double t = shape_[j] * (logA_[j] - 0.5 * lm + 0.5 * kLibRelErr * std::fabs(lm));
        Csw += t;
        mag += std::fabs(t);
      }
    }
    slope = slope * (1.0 + 2.0 * n * kU) - 2.0 * x_ * (1.0 - 2.0 * kU);
    double B = Csw + slope * (slope < 0.0 ? Ya : Yb);
    mag += std::fabs(slope) * Yb;
    double cnt = (std::sqrt(Yb) - std::sqrt(Ya)) / h + 2.0;
    double lc = std::log(cnt * (rho_ + std::sqrt(Yb)) * (1.0 + 8.0 * kU));
    acc.add(B + lc + 4.0 * kU * mag + kLibRelErr * std::fabs(lc) + 0x1p-1000);
    Ya = Yb;
  }
  return kInf;
}

inline double GammaSumDensity::termMin(double u2, double ylo, double yhi, double tau2,
                                       double tau2Err) {
  double A = u2 + tau2;
  double eA = tau2Err + 2.0 * kU * std::fabs(A);
  auto fLow = [&](double Y) {
    double d = std::fabs(A - Y);
    double dlo = std::max(0.0, d - eA - 2.0 * kU * (std::fabs(A) + Y));
    return (dlo * dlo + 4.0 * u2 * Y) * (1.0 - 8.0 * kU);
  };
  double V = tau2 - u2;
  double eV = tau2Err + 2.0 * kU * (std::fabs(tau2) + u2);
  bool finiteHi = std::isfinite(yhi);
  if (V + eV <= ylo) return fLow(ylo);
  if (finiteHi && V - eV >= yhi) return fLow(yhi);
  double m = fLow(ylo);
  if (finiteHi) m = std::min(m, fLow(yhi));
  return std::min(m, 4.0 * u2 * std::max(0.0, tau2 - tau2Err) * (1.0 - 4.0 * kU));
}

inline void GammaSumDensity::prepare(double x) {
  x_ = x;
  prepared_ = false;
  lead_ = Enclosure{};
  const std::size_t n = shape_.size();
  if (n == 0 || !(x > 0.0)) return;
  runEnd_ = 1;
  while (runEnd_ < n && shape_[runEnd_] == shape_[0]) ++runEnd_;

  // Saddle: g(t) = sum_j a_j / (mu_j + t) - x is convex and decreasing on
  // (-mu_1, inf); Newton from a point with g >= 0 rises monotonically to the
  // root. Any t0 > -mu_1 gives a valid contour.
  const double mu1 = rate_[0];
  double t = -mu1 + shape_[0] / x;
  if (!(mu1 + t > 0.0)) t = -mu1 * (1.0 - 0x1p-50);
  for (int it = 0; it < 100; ++it) {
    double g = -x, gp = 0.0;
    for (std::size_t j = 0; j < n; ++j) {
      double iv = 1.0 / (rate_[j] + t);
      g += shape_[j] * iv;
      gp -= shape_[j] * iv * iv;
    }
    double tn = t - g / gp;
    if (!(tn >= t)) break;
    bool done = tn - t <= 1e-13 * (mu1 + t);
    t = tn;
    if (done) break;
  }
  if (!(mu1 + t > 0.0)) return;
  t0_ = t;

  // Contour: T1 = (mu_1 + t0) / 2 = A_1(rho); rho^2 below it.
  double T1 = 0.5 * (mu1 + t0_);
  double target = std::max(T1 - 0.25 * kPiSquared, 0.5 * T1);
  rho_ = roundDown24(std::sqrt(target) * (1.0 - 4.0 * kU));
  rho2_ = rho_ * rho_;
  if (!(rho_ > 0.0)) return;

  vt_ = 0.0;
  for (std::size_t j = 0; j < n; ++j) {
    double s = rate_[j] + t0_;  // relative error at most u
    double iv = 1.0 / s;
    inv_[j] = iv;
    invRad_[j] = iv * (2.0 * kU + 0x1p-100);
    double la = std::log(0.5 * s);
    logA_[j] = la + kLibRelErr * std::fabs(la) + 1.01 * kU + 0x1p-1000;
    double tau2 = 0.5 * s - rho2_;
    tau2_[j] = tau2;
    tau2Err_[j] = kU * (0.5 * s) + kU * std::fabs(tau2) + 0x1p-1060;
    vt_ += shape_[j] * iv * iv;
  }

  // c0 = x - sum_j a_j / (mu_j + t0) in double-double: mu_j + t0 exact, each
  // quotient and product 15 u^2 + 2 u^2, the sum 3 u^2 per term
  {
    dd::DD acc{0.0, 0.0};
    double mag = 0.0;
    for (std::size_t j = 0; j < n; ++j) {
      dd::DD sj = dd::twoSum(rate_[j], t0_);
      dd::DD q = dd::mul(dd::div(dd::DD{1.0, 0.0}, sj), shape_[j]);
      acc = dd::add(acc, q);
      mag += q.hi;
    }
    dd::DD c = dd::sub(dd::DD{x, 0.0}, acc);
    double c0 = c.hi + c.lo;
    c0_ = RBall{c0, up(kU * std::fabs(c0) + 24.0 * 0x1p-106 * (mag + x) * (1.0 + n * 0x1p-50))};
  }

  // The lead x t0 - sum_j a_j log(1 + t0 / mu_j), each run of equal shapes as
  // a log of a double-double product rescaled by powers of two.
  dd::DD lead = dd::twoProd(x, t0_);
  double mag = std::fabs(lead.hi);
  double rad = 0.0;
  for (std::size_t from = 0; from < n;) {
    std::size_t to = from + 1;
    while (to < n && shape_[to] == shape_[from]) ++to;
    dd::DD p{1.0, 0.0};
    int e = 0;
    for (std::size_t j = from; j < to; ++j) {
      dd::DD num = dd::twoSum(rate_[j], t0_);
      p = dd::mul(p, dd::div(num, rate_[j]));
      int k;
      std::frexp(p.hi, &k);
      if (k > 256 || k < -256) {
        p = dd::ldexp(p, -k);
        e += k;
      }
    }
    dd::DD l = dd::log(p);
    if (e != 0) l = dd::add(l, dd::mul(dd::kLn2, static_cast<double>(e)));
    double a = shape_[from];
    lead = dd::sub(lead, dd::mul(l, a));
    double cnt = static_cast<double>(to - from);
    // each factor's division and product, 7 u^2; the log's bound; ln 2
    rad += a * (cnt * 8.0 * 0x1p-106 + dd::kLogRelErr * std::fabs(l.hi) +
                std::fabs(static_cast<double>(e)) * 0x1p-105);
    mag += a * std::fabs(l.hi);
    from = to;
  }
  rad += 16.0 * dd::kU2 * mag;
  double mid = lead.hi + (lead.lo + kLogTwoOverPi);
  rad += kU * std::fabs(mid) + 0x1p-55 + kU * std::fabs(lead.lo + kLogTwoOverPi);
  lead_ = Enclosure{mid, up(rad)};
  prepared_ = true;
}

inline CBall GammaSumDensity::nodeExponent(double y) const {
  // zeta_j = 2 Delta / (mu_j + t0), Delta exact; E = 2 x Delta - sum_j a_j
  // log1p(zeta_j) = 2 Delta c0 - sum_j a_j (log1p(zeta_j) - zeta_j), the
  // first-order terms summed exactly into c0.
  const double dre = -y * y;
  const double dim = 2.0 * rho_ * y;
  const double tre = 2.0 * dre;
  const double tim = 2.0 * dim;
  const double tmod = modUp(tre, tim);
  double ere = tre * c0_.m;
  double eim = tim * c0_.m;
  double r = tmod * c0_.r + kU * l1(ere, eim);
  const std::size_t n = shape_.size();
  for (std::size_t from = 0; from < n;) {
    std::size_t to = (from == 0) ? runEnd_ : from + 1;
    double lre = 0.0, lim = 0.0, lr = 0.0;
    for (std::size_t j = from; j < to; ++j) {
      double iv = inv_[j];
      CBall z{tre * iv, tim * iv, 0.0};
      z.r = tmod * invRad_[j] + kU * l1(z.re, z.im);
      CBall L = clog1pMinus(z);
      lre += L.re;
      lim += L.im;
      lr += L.r + kU * l1(lre, lim);
    }
    double a = shape_[from];
    double are = a * lre;
    double aim = a * lim;
    ere -= are;
    eim -= aim;
    r += a * lr + kU * (l1(are, aim) + l1(ere, eim));
    from = to;
  }
  return CBall{ere, eim, up(r)};
}

inline void GammaSumDensity::lineSetup(double u, LineAgg& agg) const {
  agg = LineAgg{};
  agg.u = u;
  const double u2 = u * u;
  const double dlt = rho2_ - u2;  // A_j(rho) - A_j(u)
  const double edlt = kU * (std::fabs(dlt) + u2);
  const std::size_t n = shape_.size();
  int n0 = 0;
  double kappa = 0.0, S = 0.0, C = 0.0;
  for (std::size_t j = 0; j < n; ++j) {
    double a = shape_[j];
    double A = tau2_[j] + u2;
    double eA = tau2Err_[j] + kU * (std::fabs(A) + u2);
    double Alo = A - eA;
    if (!(Alo > 0.0)) return;
    double Ahi = A + eA;
    if (Alo <= 2.0 * u2 * (1.0 + 4.0 * kU)) {
      if (static_cast<int>(j) != n0) return;  // J0 must be a prefix
      ++n0;
    } else {
      double g = 2.0 * (Ahi - 2.0 * u2) / (Alo * Alo) * (1.0 + 8.0 * kU);
      if (Alo <= 4.0 * u2 && Ahi >= 4.0 * u2) g = std::max(g, 0.25 / u2 * (1.0 + 4.0 * kU));
      kappa = std::max(kappa, g);
      S += a * g * (1.0 + 2.0 * kU);
    }
    // log(A_j(rho) / A_j(u)) = log1p(dlt / A_j(u)), bounded above
    double arg = dlt / A;
    double argUp = (dlt + edlt) / Alo;
    if (argUp < 0.0) argUp = (dlt + edlt) / Ahi;
    argUp = argUp + kU * std::fabs(argUp);
    if (!(argUp > -1.0)) return;
    double lp = std::log1p(argUp);
    C += a * (lp + kLibRelErr * std::fabs(lp)) + kU * (std::fabs(C) + a * std::fabs(lp));
    (void) arg;
  }
  agg.n0 = n0;
  agg.kappa = kappa;
  agg.s = S * (1.0 + 2.0 * n * kU);
  agg.c = C + 0x1p-1000;
  agg.valid = std::isfinite(agg.s) && std::isfinite(agg.c) && kappa > 0.0;
}

inline double GammaSumDensity::perTermLogMajor(double u2, double ylo, double yhi) const {
  const std::size_t n = shape_.size();
  double v = 0.0, mag = 0.0;
  for (std::size_t j = 0; j < n; ++j) {
    double m = termMin(u2, ylo, yhi, tau2_[j], tau2Err_[j]);
    if (!(m > 0.0)) return kInf;
    double lm = std::log(m);
    double t = shape_[j] * (logA_[j] - 0.5 * lm + 0.5 * kLibRelErr * std::fabs(lm));
    v += t;
    mag += std::fabs(t);
  }
  return v + 4.0 * kU * mag + 0x1p-1000;
}

inline double GammaSumDensity::cellLogMajor(double u, double ya, double yb,
                                            const LineAgg& agg) const {
  const double u2 = u * u;
  const double Ya = ya * ya * (1.0 - 2.0 * kU);
  const double Yb = yb * yb * (1.0 + 2.0 * kU);
  if (agg.valid && agg.kappa * Yb <= 0.5) {
    double v = agg.c + 0.5 * agg.s * Yb / (1.0 - agg.kappa * Yb) * (1.0 + 8.0 * kU);
    double mag = std::fabs(v);
    for (int j = 0; j < agg.n0; ++j) {
      double m = termMin(u2, Ya, Yb, tau2_[j], tau2Err_[j]);
      if (!(m > 0.0)) return kInf;
      double A = tau2_[j] + u2;
      double Ahi = A + tau2Err_[j] + 2.0 * kU * (std::fabs(A) + u2);
      double la = std::log(Ahi);
      double lm = std::log(m);
      double t = shape_[j] * (la - 0.5 * lm + kLibRelErr * (std::fabs(la) + 0.5 * std::fabs(lm)));
      v += t;
      mag += std::fabs(t);
    }
    return v + 4.0 * kU * mag + 0x1p-1000;
  }
  return perTermLogMajor(u2, Ya, Yb);
}

inline double GammaSumDensity::tailLogMajor(double u, double Y) const {
  return perTermLogMajor(u * u, Y * Y * (1.0 - 2.0 * kU), kInf);
}

inline double GammaSumDensity::j0TermsAt(double Y, const LineAgg& agg) const {
  const double Ylo = Y * (1.0 - 2.0 * kU);
  const double Yhi = Y * (1.0 + 2.0 * kU);
  double v = 0.0;
  for (int j = 0; j < agg.n0; ++j) {
    double m = termMin(rho2_, Ylo, Yhi, tau2_[j], tau2Err_[j]);
    if (!(m > 0.0)) return kInf;
    double lm = std::log(m);
    double t = shape_[j] * (logA_[j] - 0.5 * lm + 0.5 * kLibRelErr * std::fabs(lm));
    t += 4.0 * kU * std::fabs(t);
    v += std::min(0.0, t);
  }
  return v;
}

inline double GammaSumDensity::pointLog(double u, double y) const {
  double dre = u * u - rho2_ - y * y;
  double dim = 2.0 * u * y;
  double v = 2.0 * x_ * dre;
  for (std::size_t j = 0; j < shape_.size(); ++j) {
    double zr = 2.0 * dre * inv_[j];
    double zi = 2.0 * dim * inv_[j];
    v -= 0.5 * shape_[j] * std::log1p(zr * (2.0 + zr) + zi * zi);
  }
  return v + 0.5 * std::log(u * u + y * y);
}

// ---------------------------------------------------------------------------
// The evaluator.

inline const LevelSettings& levelSettings(const TargetDensity&, int level) {
  return kLevelsF[level];
}
inline const LevelSettings& levelSettings(const GammaSumDensity&, int level) {
  return kLevelsQ[level];
}

namespace detail {

/// log of an upper bound on 2 int_0^inf |G(u + i y)| dy on the line Re z = u
/// (the edge line's integral of the absolute integrand, both halves):
/// cells on [0, Y] each bounded by its majorant at the cell's lower end of
/// e^{-2 x y^2} and upper end of |z|, then the tail past Y.
template <class F>
double lineLogMajor(const F& fam, const LevelSettings& s, double u, double Yinit,
                    EvalCounts* counts) {
  const double x = fam.point();
  const double rho = fam.rho();
  const double u2 = u * u;
  const double rho2 = rho * rho;
  const double pref = 2.0 * x * (u2 - rho2);
  const double prefErr = 6.0 * kU * 2.0 * x * (u2 + rho2);
  LineAgg agg;
  fam.lineSetup(u, agg);
  double Y = Yinit;
  double lt0 = fam.tailLogMajor(u, Y);
  if (!std::isfinite(lt0)) return kInf;
  double need = (pref + lt0 + 62.0) / (2.0 * x);
  if (need > Y * Y) Y = std::sqrt(need);
  int n0 = static_cast<int>(std::min<double>(s.cellsMax, std::max<double>(s.cellsMin, 10.0 * Y * std::sqrt(x))));
  const double dy = Y / n0;
  const double Yend = n0 * dy;

  LogSumExpUpper acc;
  auto cellTerm = [&](double ya, double yb) {
    double c = fam.cellLogMajor(u, ya, yb, agg);
    double g = -2.0 * x * ya * ya;
    double lz = 0.5 * std::log(u2 + yb * yb);
    double lw = std::log(yb - ya);
    double v = g + c + lz + lw;
    return v + 4.0 * kLibRelErr * (std::fabs(g) + std::fabs(c) + std::fabs(lz) + std::fabs(lw)) + 0x1p-1000;
  };
  if (s.maxDepth == 0) {
    for (int j = 0; j < n0; ++j) {
      double ya = j * dy;
      double yb = (j + 1 == n0) ? Yend : (j + 1) * dy;
      acc.add(cellTerm(ya, yb));
    }
    if (counts) counts->cells += n0;
  } else {
    // split a cell where its majorant passes the integrand's value at its
    // ends by more than e^2 and it is not negligible (a choice only)
    std::vector<double> pts(n0 + 1);
    double est = -kInf;
    for (int j = 0; j <= n0; ++j) {
      pts[j] = fam.pointLog(u, j * dy);
      if (std::isfinite(pts[j])) est = std::max(est, pts[j]);
    }
    est += std::log(Y);
    struct Cell { double ya, yb, pv; int depth; };
    std::vector<Cell> stack;
    stack.reserve(n0 + 4 * s.maxDepth);
    for (int j = n0 - 1; j >= 0; --j) {
      double yb = (j + 1 == n0) ? Yend : (j + 1) * dy;
      stack.push_back(Cell{j * dy, yb, std::max(pts[j], pts[j + 1]), 0});
    }
    while (!stack.empty()) {
      Cell c = stack.back();
      stack.pop_back();
      double m = cellTerm(c.ya, c.yb);
      if (counts) ++counts->cells;
      double w = std::log(c.yb - c.ya);
      if (c.depth < s.maxDepth && m - w > est - 40.0 - w && m - w - c.pv > 2.0) {
        double ym = 0.5 * (c.ya + c.yb);
        double pm = fam.pointLog(u, ym);
        double pv = std::isfinite(pm) ? std::max(c.pv, pm) : c.pv;
        stack.push_back(Cell{ym, c.yb, pv, c.depth + 1});
        stack.push_back(Cell{c.ya, ym, pv, c.depth + 1});
        continue;
      }
      acc.add(m);
    }
  }
  double lt = fam.tailLogMajor(u, Yend);
  double gt = -2.0 * x * Yend * Yend;
  double lg = std::log((u / Yend + 1.0) / (4.0 * x));
  acc.add(lt + gt + lg + 4.0 * kLibRelErr * (std::fabs(lt) + std::fabs(gt) + std::fabs(lg)) + 0x1p-1000);
  double v = acc.value();
  return kLog2 + pref + v + prefErr + 2.0 * kU * (std::fabs(pref) + std::fabs(v) + 1.0);
}

}  // namespace detail

/// The enclosure of log p(x) at one refinement level for a family prepared
/// at x; undecided (infinite radius) when the level's certificate does not
/// close.
template <class F>
Enclosure logDensity(const F& fam, int level, EvalCounts* counts = nullptr) {
  const LevelSettings& s = levelSettings(fam, level);
  const Enclosure lead = fam.lead();
  if (!lead.finite()) return Enclosure{};
  const double x = fam.point();
  const double rho = fam.rho();
  const double vt = fam.tiltedVariance();
  if (!(vt > 0.0) || !(rho > 0.0)) return Enclosure{};
  const double cEst = 8.0 * vt * rho * rho;
  const double sEst = std::sqrt(kPi / (8.0 * vt));
  const double yInit = std::max(std::sqrt(44.0 / cEst), 8.0 / std::sqrt(x));
  const double logEps = std::log(sEst) - s.targetBits * kLog2;

  // The edge lines and the step from their bounds: on a line at distance a
  // the discretization error is at most M / (e^{2 pi a / h} - 1), so
  // h = 2 pi a / log(1 + 2 M / eps) on each side's best line. A line's
  // distance is capped at 4 / sqrt(x): log M grows about as 2 x a^2, so
  // farther lines only shrink h (a choice; any distance inside the strip is
  // valid).
  const double aCap = 4.0 / std::sqrt(x);
  double logMp[4], aP[4], logMm[6], aM[6];
  double hPlus = 0.0, hMinus = 0.0;
  for (int i = 0; i < s.numPlus; ++i) {
    double u = rho * s.fracPlus[i] <= aCap ? rho * (1.0 - s.fracPlus[i]) : rho - aCap;
    aP[i] = rho - u;  // exact: u >= rho / 2
    logMp[i] = detail::lineLogMajor(fam, s, u, yInit, counts);
    if (std::isfinite(logMp[i])) {
      double v = logMp[i] - logEps;
      double d = v > 30.0 ? v + kLog2 : std::log1p(2.0 * std::exp(v));
      hPlus = std::max(hPlus, 2.0 * kPi * aP[i] / d);
    }
  }
  for (int i = 0; i < s.numMinus; ++i) {
    double a0 = std::min(s.aMinus[i] < 0.0 ? -s.aMinus[i] * rho : s.aMinus[i], aCap);
    double u = rho + a0;
    aM[i] = (u - rho) - 2.0 * kU * u;
    logMm[i] = detail::lineLogMajor(fam, s, u, yInit, counts);
    if (std::isfinite(logMm[i]) && aM[i] > 0.0) {
      double v = logMm[i] - logEps;
      double d = v > 30.0 ? v + kLog2 : std::log1p(2.0 * std::exp(v));
      hMinus = std::max(hMinus, 2.0 * kPi * aM[i] / d);
    }
  }
  if (!(hPlus > 0.0) || !(hMinus > 0.0)) return Enclosure{};
  const double h = roundDownStep(std::min(hPlus, hMinus) * (1.0 - 0x1p-40));
  if (counts) counts->h = h;

  auto sideBound = [&](const double* logM, const double* a, int n) {
    double best = kInf;
    for (int i = 0; i < n; ++i) {
      if (!std::isfinite(logM[i]) || !(a[i] > 0.0)) continue;
      double arg = 2.0 * kPi * a[i] / h * (1.0 - 4.0 * kU);
      // log(e^arg - 1) from below
      double le = arg > 1.0 ? arg + std::log1p(-std::exp(-arg)) : std::log(std::expm1(arg));
      le -= kLibRelErr * (std::fabs(le) + 2.0) + 4.0 * kU * arg;
      double v = logM[i] - le;
      best = std::min(best, std::exp(v + kLibRelErr * std::fabs(v)) * (1.0 + 2.0 * kLibRelErr));
    }
    return best;
  };
  const double Ed = sideBound(logMp, aP, s.numPlus) + sideBound(logMm, aM, s.numMinus);
  if (!std::isfinite(Ed)) return Enclosure{};

  // The truncation's pieces on the contour line (the monotone tail of the
  // aggregated majorant) and the depth they need the main sum to reach.
  LineAgg agg0;
  fam.lineSetup(rho, agg0);
  bool corOk = false;
  double Y1 = 0.0, delta = 0.0, yReq = 0.0;
  if (agg0.valid && agg0.s < 4.0 * x * (1.0 - 16.0 * kU)) {
    const double S = agg0.s;
    const double kappa = agg0.kappa;
    if (S < 2.0 * x * (1.0 - 16.0 * kU)) {
      Y1 = (1.0 - std::sqrt(S / (2.0 * x))) / kappa;
    } else {
      Y1 = 0.5 * (1.0 - std::sqrt(S / (4.0 * x))) / kappa;
    }
    Y1 *= 1.0 - 0x1p-30;
    double om = (1.0 - kappa * Y1) * (1.0 - 4.0 * kU);
    delta = 2.0 * x * (1.0 - kU) - 0.5 * S / (om * om) * (1.0 + 6.0 * kU);
    if (om > 0.0 && delta > 0.0) {
      corOk = true;
      yReq = 0.5 * (-rho + std::sqrt(rho * rho + 2.0 / delta)) * (1.0 + 0x1p-40) + 0x1p-1000;
    }
  }
  const double yW = 0.5 * (-rho + std::sqrt(rho * rho + 1.0 / x)) * (1.0 + 0x1p-40);
  const double yStop = std::max(1.0 / std::sqrt(x), corOk ? yReq : 0.0);

  // The main sum: h [G(0) + 2 sum_k Re G(k h)], G(0) = rho.
  dd::DD sum{rho, 0.0};
  double rad = 0.0, ddErr = 0.0, absSum = rho;
  long k = 0;
  for (;;) {
    ++k;
    if (k >= (1L << 20)) return Enclosure{};
    const double y = k * h;
    CBall E = fam.nodeExponent(y);
    if (!std::isfinite(E.r) || !std::isfinite(E.re)) return Enclosure{};
    CBall G = cmul(cexp(E), CBall{rho, y, 0.0});
    sum = dd::add(sum, 2.0 * G.re);
    ddErr += std::fabs(sum.hi);
    rad += 2.0 * G.r;
    double gabs = modUp(G.re, G.im) + G.r;
    absSum += 2.0 * gabs;
    if (k > 4 && gabs < 0x1p-61 * absSum && y > yStop && (k + 1) * h >= yStop) break;
  }
  const long N = k;
  if (counts) counts->nodes += N;
  dd::DD sh = dd::mul(sum, h);
  const double sHat = sh.hi + sh.lo;
  double rRound = up(h * rad * (1.0 + 4.0 * kU) + h * ddErr * dd::kU2 + 4.0 * dd::kU2 * h * std::fabs(sum.hi) +
                     kU * std::fabs(sHat));
  if (!(sHat > 0.0)) return Enclosure{};
  const double epsT = std::ldexp(sHat, -s.targetBits);

  // Truncation.
  double Et = kInf;
  if (corOk) {
    const double y1 = (N + 1) * h;
    const double Y = y1 * y1;
    if (Y < Y1) {
      const double S = agg0.s;
      const double kappa = agg0.kappa;
      double j0 = fam.j0TermsAt(Y, agg0);
      double phi = -2.0 * x * Y * (1.0 - 2.0 * kU) + agg0.c + j0 +
                   0.5 * S * Y / ((1.0 - kappa * Y) * (1.0 - 4.0 * kU)) * (1.0 + 6.0 * kU);
      phi += 4.0 * kU * (2.0 * x * Y + std::fabs(agg0.c) + std::fabs(j0) + S * Y) + 0x1p-1000;
      double b1 = std::exp(phi) * (1.0 + 2.0 * kLibRelErr) *
                  ((rho + y1) + (rho / y1 + 1.0) / (2.0 * delta * h)) * (1.0 + 8.0 * kU);
      double yp = std::sqrt(Y1) * (1.0 - 4.0 * kU);
      double b2 = fam.farTail(std::max(yp, y1), h);
      Et = 2.0 * h * (b1 + b2) * (1.0 + 4.0 * kU);
    } else {
      Et = 2.0 * h * fam.farTail(y1, h) * (1.0 + 4.0 * kU);
    }
    if (counts) counts->tailPath = 1;
  }
  if (!(Et <= epsT)) {
    // node by node, with the monotone tail checked every 16 nodes
    double tot = 0.0;
    double T = kInf;
    long kk = N;
    for (;;) {
      ++kk;
      if (kk - N > 400000 || kk >= (1L << 20)) break;
      const double y = kk * h;
      if ((kk - N) % 16 == 1 && y >= yW) {
        double t = fam.farTail(y, h);
        if (t <= 0.0625 * epsT) {
          T = t;
          break;
        }
      }
      CBall E = fam.nodeExponent(y);
      if (!std::isfinite(E.r)) break;
      double g = std::exp(E.re + E.r) * (1.0 + 2.0 * kLibRelErr) * modUp(rho, y);
      tot += g * (1.0 + 2.0 * kU);
    }
    if (counts) {
      counts->truncationNodes += static_cast<std::size_t>(kk - N);
      counts->tailPath = 2;
    }
    double e2 = 2.0 * h * (tot + T) * (1.0 + 4.0 * kU);
    Et = std::min(Et, e2);
  }
  if (!std::isfinite(Et)) return Enclosure{};

  const double R = up(rRound + Ed + Et);
  if (!(sHat > R)) return Enclosure{};
  const double ls = std::log(sHat);
  const double e = R / sHat * (1.0 + 2.0 * kU);
  const double radS = e / (1.0 - e) * (1.0 + 4.0 * kU) + kLibRelErr * std::fabs(ls);
  const double mid = lead.mid + ls;
  return Enclosure{mid, up(lead.rad + radS + kU * std::fabs(mid))};
}

/// The enclosure intersected over levels 0 through level.
template <class F>
Enclosure logDensityThrough(const F& fam, int level, EvalCounts* counts = nullptr) {
  Enclosure e;
  for (int l = 0; l <= level; ++l) e = intersect(e, logDensity(fam, l, counts));
  return e;
}

// ---------------------------------------------------------------------------
// The proposal: q = sum_{k <= K} Ga(b, r_k) + Ga(alpha, beta), r_k the
// target's rates d_k = 2 pi^2 (k - 1/2)^2 + psi^2 / 2 in binary64 with r_1
// rounded down (so r_1 < d_1), and (alpha, beta) matching the mean and
// variance of the remainder sum_{k > K} Ga(b, d_k).

namespace detail {

/// d_k in binary64 in a fixed order of operations.
inline double targetRate(int k, double halfPsiSquared) {
  double h = k - 0.5;
  // unfused, so the rates are the same doubles on every build
  return dd::detail::opaque(((2.0 * kPi) * kPi) * (h * h)) + halfPsiSquared;
}

/// sum_{k > K} d_k^-n, n = 1 or 2: direct to K + 128, then the midpoint
/// Euler-Maclaurin formula with the g' and g''' terms (relative error below
/// 1e-13).
inline double tailPowerSum(int K, double sigma, int n) {
  const int M = K + 128;
  const double a = (2.0 * kPi) * kPi;
  double s = 0.0;
  for (int k = M; k > K; --k) {
    double d = targetRate(k, sigma);
    s += (n == 1) ? 1.0 / d : 1.0 / (d * d);
  }
  // g(s') = (sigma + a s'^2)^-n at s' = t - 1/2; the remainder from t = M + 1/2
  const double s0 = M;
  const double G = sigma + a * s0 * s0;
  double integral;
  if (sigma == 0.0) {
    integral = (n == 1) ? 1.0 / (a * s0) : 1.0 / (3.0 * a * a * s0 * s0 * s0);
  } else {
    double eps = sigma / (a * s0 * s0);
    double at = std::atan(std::sqrt(sigma / a) / s0) / std::sqrt(a * sigma);
    if (n == 1) {
      integral = at;
    } else if (eps >= 0.1) {
      integral = (at - s0 / G) / (2.0 * sigma);
    } else {
      // (1 / (a^2 s0^3)) sum_m (-1)^m (m + 1) eps^m / (2m + 3)
      double t = 0.0;
      for (int m = 20; m >= 0; --m) t = (m + 1.0) / (2.0 * m + 3.0) - eps * t;
      integral = t / (a * a * s0 * s0 * s0);
    }
  }
  const double nn = n;
  const double g1 = -2.0 * nn * a * s0 * std::pow(G, -nn - 1.0);
  const double g3 = 12.0 * nn * (nn + 1.0) * a * a * s0 * std::pow(G, -nn - 2.0) -
                    8.0 * nn * (nn + 1.0) * (nn + 2.0) * a * a * a * s0 * s0 * s0 * std::pow(G, -nn - 3.0);
  return s + integral + g1 / 24.0 - 7.0 / 5760.0 * g3;
}

/// log cosh(a) for a >= 0.
inline double logCosh(double a) { return a - kLog2 + std::log1p(std::exp(-2.0 * a)); }

}  // namespace detail

/// D(psi, K), the K rule's closed form in binary64:
/// -log(d_1 / (2 pi cosh(psi / 2))) + sum_{k=2}^K log1p(-d_1 / d_k)
/// + (m^2 / v) log1p(-d_1 / beta), m and v the remainder's tail sums.
inline double ruleD(double psi, int K) {
  const double sigma = (0.5 * psi) * psi;
  const double d1 = detail::targetRate(1, sigma);
  double s = 0x1.d67f1c864beb5p+0 + detail::logCosh(0.5 * std::fabs(psi)) - std::log(d1);
  for (int k = 2; k <= K; ++k) s += std::log1p(-d1 / detail::targetRate(k, sigma));
  const double m = detail::tailPowerSum(K, sigma, 1);
  const double v = detail::tailPowerSum(K, sigma, 2);
  const double beta = m / v;
  return s + (m * m / v) * std::log1p(-d1 / beta);
}

/// The smallest K with b D(psi, K) <= tol: a guess from D's leading order,
/// then an exponential search and bisection (D decreasing in K).
inline int ruleK(double b, double psi, double tol = 0.01) {
  const double sigma = (0.5 * psi) * psi;
  const double d1 = detail::targetRate(1, sigma);
  const double pi6 = kPiSquared * kPiSquared * kPiSquared;
  int guess = static_cast<int>(std::ceil(std::pow(b * d1 * d1 * d1 * 0.0037 / (tol * pi6), 0.2)));
  guess = std::max(1, std::min(guess, 1 << 20));
  auto passes = [&](int K) { return b * ruleD(psi, K) <= tol; };
  int lo, hi;  // lo fails (or 0), hi passes
  if (passes(guess)) {
    hi = guess;
    int step = 1;
    lo = hi - step;
    while (lo >= 1 && passes(lo)) {
      hi = lo;
      step *= 2;
      lo = hi - step;
    }
    if (lo < 1) lo = 0;
  } else {
    lo = guess;
    int step = 1;
    hi = lo + step;
    while (!passes(hi)) {
      lo = hi;
      step *= 2;
      hi = lo + step;
      if (hi > (1 << 20)) return -1;
    }
  }
  while (hi - lo > 1) {
    int mid = lo + (hi - lo) / 2;
    if (mid >= 1 && passes(mid)) hi = mid; else lo = mid;
  }
  return hi;
}

/// One row's proposal for PG(b, psi).
class Proposal {
 public:
  Proposal() = default;
  Proposal(double b, double psi, int K = 0) { build(b, psi, K); }

  /// Builds the proposal at the rule's K, or at K when positive.
  void build(double b, double psi, int K = 0) {
    b_ = b;
    psi_ = std::fabs(psi);
    K_ = K > 0 ? K : ruleK(b, psi_);
    const double sigma = (0.5 * psi_) * psi_;
    rates_.resize(K_ + 1);
    shapes_.assign(K_ + 1, b);
    for (int k = 1; k <= K_; ++k) rates_[k - 1] = detail::targetRate(k, sigma);
    // r_1 below d_1: the rate's relative rounding error is at most about
    // 3.7 u, so a reduction by 2^-50 = 8 u clears it
    rates_[0] *= 1.0 - 0x1p-50;
    const double m = detail::tailPowerSum(K_, sigma, 1);
    const double v = detail::tailPowerSum(K_, sigma, 2);
    alpha_ = b * m * m / v;
    beta_ = m / v;
    shapes_[K_] = alpha_;
    rates_[K_] = beta_;
  }

  double shape() const { return b_; }
  double psi() const { return psi_; }
  int numTerms() const { return K_; }
  double rate(int k) const { return rates_[k]; }
  double alpha() const { return alpha_; }
  double beta() const { return beta_; }
  const double* shapes() const { return shapes_.data(); }
  const double* rates() const { return rates_.data(); }

  /// A draw: sum_k Ga(b, 1) / r_k + Ga(alpha, 1) / beta, K + 1 gammas in order.
  double draw(ext_rng* rng) const {
    double s = 0.0;
    for (int k = 0; k < K_; ++k) s += ext_rng_simulateGamma(rng, b_, 1.0) / rates_[k];
    return s + ext_rng_simulateGamma(rng, alpha_, 1.0) / beta_;
  }

  /// The proposal's density family.
  void density(GammaSumDensity& q) const { q.setParameters(shapes_.data(), rates_.data(), K_ + 1); }

 private:
  double b_ = 1.0;
  double psi_ = 0.0;
  int K_ = 1;
  double alpha_ = 0.0;
  double beta_ = 0.0;
  std::vector<double> rates_;
  std::vector<double> shapes_;
};

// ---------------------------------------------------------------------------
// The one-row step.

/// The step floor: twice the half-width of log r (the sum of four radii,
/// each density's largest radius taken twice) plus the uniform's grain.
inline double stepFloor(double radF, double radQ, int beta = 53) {
  return 2.0 * (2.0 * radF + 2.0 * radQ) + std::ldexp(1.0, 1 - beta);
}

/// 32-bit words from an ext_rng: its uniform is word / (2^32 - 1), so
/// floor(u 2^32) recovers the word.
struct ExtRngWordSource {
  ext_rng* rng;
  std::uint32_t operator()() {
    return static_cast<std::uint32_t>(ext_rng_simulateContinuousUniform(rng) * 4294967296.0);
  }
};

/// A uniform known to lie in [lower, upper) = [k 2^-53, (k + 1) 2^-53).
struct UniformInterval {
  std::uint64_t k = 0;
  double lower() const { return std::ldexp(static_cast<double>(k), -53); }
  double upper() const { return std::ldexp(static_cast<double>(k + 1), -53); }
  /// An upper bound on log(upper()).
  double logUpperUp() const {
    double l = std::log(upper());
    return l + kLibRelErr * std::fabs(l) + 0x1p-1000;
  }
  /// A lower bound on log(lower()); -inf at k = 0.
  double logLowerDown() const {
    if (k == 0) return -kInf;
    double l = std::log(lower());
    return l - kLibRelErr * std::fabs(l) - 0x1p-1000;
  }
};

template <class WordSource>
UniformInterval drawUniformInterval(WordSource& words) {
  std::uint64_t w1 = words();
  std::uint64_t w2 = words();
  return UniformInterval{(w1 << 21) | (w2 >> 11)};
}

struct StepWorkspace {
  TargetDensity fx, fy;
  GammaSumDensity qx, qy;
};

struct StepResult {
  double omega = 0.0;
  bool accepted = false;
  bool closeCall = false;
  int level = 0;  // the last level evaluated
};

/// One independence Metropolis-Hastings step for PG(b, psi) from the
/// current point x: a proposal draw y, a uniform interval, then levels until
/// log r = log f(y) + log q(x) - log f(x) - log q(y) is certified above the
/// uniform's log (accept) or below it (keep); undecided at the last level
/// keeps x as a close call.
template <class WordSource>
StepResult stepRow(const Proposal& prop, double x, ext_rng* rng, WordSource& words,
                   StepWorkspace& ws) {
  StepResult res;
  res.omega = x;
  const double y = prop.draw(rng);
  const UniformInterval U = drawUniformInterval(words);
  if (!(y > 0.0) || !std::isfinite(y)) {
    res.closeCall = true;
    return res;
  }
  ws.fx.setParameters(prop.shape(), prop.psi());
  ws.fy.setParameters(prop.shape(), prop.psi());
  ws.fx.prepare(x);
  ws.fy.prepare(y);
  prop.density(ws.qx);
  prop.density(ws.qy);
  ws.qx.prepare(x);
  ws.qy.prepare(y);
  const double accTh = std::min(0.0, U.logUpperUp());
  const double keepTh = U.logLowerDown();
  Enclosure efy, efx, eqx, eqy;
  for (int level = 0; level < kNumLevels; ++level) {
    res.level = level;
    efy = intersect(efy, logDensity(ws.fy, level));
    eqx = intersect(eqx, logDensity(ws.qx, level));
    efx = intersect(efx, logDensity(ws.fx, level));
    eqy = intersect(eqy, logDensity(ws.qy, level));
    if (!(efy.finite() && eqx.finite() && efx.finite() && eqy.finite())) continue;
    double mid = (efy.mid + eqx.mid) - (efx.mid + eqy.mid);
    double rad = up(efy.rad + eqx.rad + efx.rad + eqy.rad +
                    4.0 * kU * (std::fabs(efy.mid) + std::fabs(eqx.mid) + std::fabs(efx.mid) +
                                std::fabs(eqy.mid)));
    if (mid - rad >= accTh) {
      res.omega = y;
      res.accepted = true;
      return res;
    }
    if (mid + rad <= keepTh) return res;
  }
  res.closeCall = true;
  return res;
}

}  // namespace pgdensity
}  // namespace bartcore

#endif  // BARTCORE_PGDENSITY_HPP
