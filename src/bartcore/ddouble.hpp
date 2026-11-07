#ifndef BARTCORE_DDOUBLE_HPP
#define BARTCORE_DDOUBLE_HPP

// Double-double arithmetic: a value is the unevaluated sum hi + lo of two
// binary64 doubles with |lo| <= ulp(hi) / 2. Each arithmetic operation below
// is one of the algorithms of Joldes, Muller and Popescu (2017), whose
// relative error bounds (stated per function, all at most 16 u^2 = 2^-102
// with u = 2^-53) hold in round-to-nearest binary64 with a correctly rounded
// fused multiply-add and no overflow or underflow. The five functions (exp,
// log, log1p, sin, sqrt) each state a relative error bound, kExpRelErr and
// the others, derived in the comment at each, all at most 2^-100; a test
// measures each against 200-digit references over the ranges the density
// evaluator reaches.

#include <cmath>
#include <cstdint>

namespace bartcore {
namespace ddouble {

/// Unit roundoff of binary64.
inline constexpr double kU = 0x1p-53;
/// Relative error bound shared by every arithmetic operation (16 u^2).
inline constexpr double kU2 = 0x1p-102;

namespace detail {

/// v, hidden from the optimizer: a product passed through it cannot be fused
/// into a later sum. The error-free transformations need each sum and product
/// rounded as written; GCC's default in GNU modes contracts a * b + c into a
/// fused multiply-add even across statements once inlined, which would break
/// them. Costs no instruction.
inline double opaque(double v) {
#if defined(__GNUC__) && (defined(__x86_64__) || defined(__i386__)) && defined(__SSE2__)
  __asm__("" : "+x"(v));
#elif defined(__GNUC__) && defined(__aarch64__)
  __asm__("" : "+w"(v));
#elif defined(__GNUC__)
  __asm__("" : "+m"(v));
#endif
  return v;
}

}  // namespace detail

/// A double-double value hi + lo.
struct DD {
  double hi;
  double lo;
};

inline constexpr DD fromDouble(double a) { return DD{a, 0.0}; }
inline double toDouble(DD a) { return a.hi + a.lo; }

/// Exact sum: s + e == a + b (Knuth's 2Sum).
inline DD twoSum(double a, double b) {
  a = detail::opaque(a);
  b = detail::opaque(b);
  double s = a + b;
  double bb = s - a;
  double e = (a - (s - bb)) + (b - bb);
  return DD{s, e};
}

/// Exact sum when a == 0 or |a| >= |b| (Dekker's Fast2Sum).
inline DD fastTwoSum(double a, double b) {
  a = detail::opaque(a);
  b = detail::opaque(b);
  double s = a + b;
  double e = b - (s - a);
  return DD{s, e};
}

/// Exact product: p + e == a * b, by a fused multiply-add.
inline DD twoProd(double a, double b) {
  double p = detail::opaque(a * b);
  double e = std::fma(a, b, -p);
  return DD{p, e};
}

/// DWPlusFP (JMP Algorithm 4): relative error at most 2 u^2.
inline DD add(DD x, double y) {
  DD s = twoSum(x.hi, y);
  double v = x.lo + s.lo;
  return fastTwoSum(s.hi, v);
}

/// AccurateDWPlusDW (JMP Algorithm 6): relative error at most 3 u^2 + 13 u^3.
inline DD add(DD x, DD y) {
  DD s = twoSum(x.hi, y.hi);
  DD t = twoSum(x.lo, y.lo);
  double c = s.lo + t.hi;
  DD v = fastTwoSum(s.hi, c);
  double w = t.lo + v.lo;
  return fastTwoSum(v.hi, w);
}

inline DD neg(DD x) { return DD{-x.hi, -x.lo}; }
inline DD sub(DD x, DD y) { return add(x, neg(y)); }
inline DD sub(DD x, double y) { return add(x, -y); }

/// DWTimesFP3 (JMP Algorithm 9): relative error at most 2 u^2.
inline DD mul(DD x, double y) {
  DD c = twoProd(x.hi, y);
  double cl3 = std::fma(x.lo, y, c.lo);
  return fastTwoSum(c.hi, cl3);
}

/// DWTimesDW3 (JMP Algorithm 12): relative error at most 4 u^2.
inline DD mul(DD x, DD y) {
  DD c = twoProd(x.hi, y.hi);
  double tl0 = x.lo * y.lo;
  double tl1 = std::fma(x.hi, y.lo, tl0);
  double cl2 = std::fma(x.lo, y.hi, tl1);
  double cl3 = c.lo + cl2;
  return fastTwoSum(c.hi, cl3);
}

/// DWDivFP3 (JMP Algorithm 15): relative error at most 3 u^2.
inline DD div(DD x, double y) {
  double th = x.hi / y;
  DD p = twoProd(th, y);
  double dh = x.hi - p.hi;
  double dt = dh - p.lo;
  double d = dt + x.lo;
  double tl = d / y;
  return fastTwoSum(th, tl);
}

/// DWTimesFP1 (JMP Algorithm 7): relative error at most 1.5 u^2 + 4 u^3.
inline DD mulAccurateFP1(DD x, double y) {
  DD c = twoProd(x.hi, y);
  double cl2 = x.lo * y;
  DD t = fastTwoSum(c.hi, cl2);
  double tl2 = t.lo + c.lo;
  return fastTwoSum(t.hi, tl2);
}

/// DWDivDW2 (JMP Algorithm 17): relative error at most 15 u^2 + 56 u^3.
inline DD div(DD x, DD y) {
  double th = x.hi / y.hi;
  DD r = mulAccurateFP1(y, th);
  DD p = twoSum(x.hi, -r.hi);
  double dl = x.lo - r.lo;
  double d = p.lo + dl;
  d = p.hi + d;
  double tl = d / y.hi;
  return fastTwoSum(th, tl);
}

/// x * 2^e, exact while neither part leaves the normal range.
inline DD ldexp(DD x, int e) { return DD{std::ldexp(x.hi, e), std::ldexp(x.lo, e)}; }

/// |x - y| / |y| as a double, for tests.
inline double relativeDifference(DD x, DD y) {
  double d = (x.hi - y.hi) + (x.lo - y.lo);
  return std::fabs(d) / std::fabs(y.hi);
}

// ln 2 = kLn2 (hi + lo, a double-double rounded to 106 bits) and as three
// parts for the exponential's reduction: kLn2Hi30 holds 30 bits, so k times
// it is exact for |k| < 2^23.
inline constexpr DD kLn2{0x1.62e42fefa39efp-1, 0x1.abc9e3b39803fp-56};
inline constexpr double kLn2Hi30 = 0x1.62e42fe800000p-1;
inline constexpr double kLn2Mid = 0x1.e8e7bcd5e4f1ep-31;
inline constexpr double kLn2Lo = -0x1.8cff81a12a17ep-85;
/// pi / 2 rounded to a double-double (relative error at most 2^-107).
inline constexpr DD kPiOver2{0x1.921fb54442d18p+0, 0x1.1a62633145c07p-54};

namespace detail {

// 1 / (i + 1)! for i = 0..15, rounded to double-double.
inline constexpr DD kInvFactorial[16] = {
  {0x1.0000000000000p+0, 0x0.0p+0},
  {0x1.0000000000000p-1, 0x0.0p+0},
  {0x1.5555555555555p-3, 0x1.5555555555555p-57},
  {0x1.5555555555555p-5, 0x1.5555555555555p-59},
  {0x1.1111111111111p-7, 0x1.1111111111111p-63},
  {0x1.6c16c16c16c17p-10, -0x1.f49f49f49f49fp-65},
  {0x1.a01a01a01a01ap-13, 0x1.a01a01a01a01ap-73},
  {0x1.a01a01a01a01ap-16, 0x1.a01a01a01a01ap-76},
  {0x1.71de3a556c734p-19, -0x1.c154f8ddc6c00p-73},
  {0x1.27e4fb7789f5cp-22, 0x1.cbbc05b4fa99ap-76},
  {0x1.ae64567f544e4p-26, -0x1.c062e06d1f209p-80},
  {0x1.1eed8eff8d898p-29, -0x1.2aec959e14c06p-83},
  {0x1.6124613a86d09p-33, 0x1.f28e0cc748ebep-87},
  {0x1.93974a8c07c9dp-37, 0x1.05d6f8a2efd1fp-92},
  {0x1.ae7f3e733b81fp-41, 0x1.1d8656b0ee8cbp-97},
  {0x1.ae7f3e733b81fp-45, 0x1.1d8656b0ee8cbp-101},
};

// 1 / (2 i + 1) for i = 0..20, rounded to double-double.
inline constexpr DD kInvOdd[21] = {
  {0x1.0000000000000p+0, 0x0.0p+0},
  {0x1.5555555555555p-2, 0x1.5555555555555p-56},
  {0x1.999999999999ap-3, -0x1.999999999999ap-57},
  {0x1.2492492492492p-3, 0x1.2492492492492p-57},
  {0x1.c71c71c71c71cp-4, 0x1.c71c71c71c71cp-58},
  {0x1.745d1745d1746p-4, -0x1.745d1745d1746p-59},
  {0x1.3b13b13b13b14p-4, -0x1.3b13b13b13b14p-58},
  {0x1.1111111111111p-4, 0x1.1111111111111p-60},
  {0x1.e1e1e1e1e1e1ep-5, 0x1.e1e1e1e1e1e1ep-61},
  {0x1.af286bca1af28p-5, 0x1.af286bca1af28p-59},
  {0x1.8618618618618p-5, 0x1.8618618618618p-59},
  {0x1.642c8590b2164p-5, 0x1.642c8590b2164p-60},
  {0x1.47ae147ae147bp-5, -0x1.eb851eb851eb8p-61},
  {0x1.2f684bda12f68p-5, 0x1.2f684bda12f68p-59},
  {0x1.1a7b9611a7b96p-5, 0x1.1a7b9611a7b96p-61},
  {0x1.0842108421084p-5, 0x1.0842108421084p-60},
  {0x1.f07c1f07c1f08p-6, -0x1.f07c1f07c1f08p-61},
  {0x1.d41d41d41d41dp-6, 0x1.0750750750750p-60},
  {0x1.bacf914c1bad0p-6, -0x1.bacf914c1bad0p-60},
  {0x1.a41a41a41a41ap-6, 0x1.0690690690690p-60},
  {0x1.8f9c18f9c18fap-6, -0x1.f3831f3831f38p-61},
};

// (-1)^i / (2 i + 1)! for i = 0..17, rounded to double-double.
inline constexpr DD kSinCoef[18] = {
  {0x1.0000000000000p+0, 0x0.0p+0},
  {-0x1.5555555555555p-3, -0x1.5555555555555p-57},
  {0x1.1111111111111p-7, 0x1.1111111111111p-63},
  {-0x1.a01a01a01a01ap-13, -0x1.a01a01a01a01ap-73},
  {0x1.71de3a556c734p-19, -0x1.c154f8ddc6c00p-73},
  {-0x1.ae64567f544e4p-26, 0x1.c062e06d1f209p-80},
  {0x1.6124613a86d09p-33, 0x1.f28e0cc748ebep-87},
  {-0x1.ae7f3e733b81fp-41, -0x1.1d8656b0ee8cbp-97},
  {0x1.952c77030ad4ap-49, 0x1.ac981465ddc6cp-103},
  {-0x1.2f49b46814157p-57, -0x1.2650f61dbdcb4p-112},
  {0x1.71b8ef6dcf572p-66, -0x1.d043ae40c4647p-120},
  {-0x1.761b41316381ap-75, 0x1.3423c7d91404fp-130},
  {0x1.3f3ccdd165fa9p-84, -0x1.58ddadf344487p-139},
  {-0x1.d1ab1c2dccea3p-94, -0x1.054d0c78aea14p-149},
  {0x1.259f98b4358adp-103, 0x1.eaf8c39dd9bc5p-157},
  {-0x1.434d2e783f5bcp-113, -0x1.0b87b91be9affp-167},
  {0x1.3981254dd0d52p-123, -0x1.2b1f4c8015a2fp-177},
  {-0x1.0dc59c716d91fp-133, -0x1.419e3fad3f031p-188},
};

/// expm1(r) for |r| <= 0.35: Taylor of degree 15 at s = r / 8, then three
/// steps q <- q (2 + q). Relative error: the reduced Taylor sum s P(s),
/// |s| <= 0.0438, has truncation below 2^-111 and Horner error at most
/// 4.5 u^2 (each step's error damped by |s| / (i + 1)), 9 u^2 after the
/// product by s; a step multiplies the relative error by 1 + |q| / |2 + q|
/// and adds 6 u^2 (one DWPlusFP, 2 u^2 of 2 + q, and one DWTimesDW3, 4 u^2).
/// Entering step i (i = 1, 2, 3), |q| <= e^{0.35 / 2^(4 - i)} - 1, at most
/// 0.0448, 0.0915 and 0.1912, so the factors are at most 1.0230, 1.0479 and
/// 1.1057 and three steps give at most ((9 1.0230 + 6) 1.0479 + 6) 1.1057
/// + 6 <= 30.3 u^2.
inline DD expm1Reduced(DD r) {
  DD s = ldexp(r, -3);
  DD p = kInvFactorial[15];
  for (int i = 14; i >= 0; --i) p = add(kInvFactorial[i], mul(s, p));
  DD q = mul(s, p);
  for (int i = 0; i < 3; ++i) q = mul(q, add(q, 2.0));
  return q;
}

/// 2 atanh(s) for |s| <= 0.1716 (s^2 <= 0.02944), the series of 21 terms in
/// z = s^2 by Horner: truncation below 2^-112 relative; Horner error at most
/// 4.6 u^2 (one AccurateDWPlusDW and a coefficient rounding per step, each
/// later step damped by z / 3); the product by s adds 4 u^2. The input's own
/// relative error passes through with factor
/// s / ((1 - s^2) atanh s) <= 1.0202 at |s| = 0.1716.
inline DD twoAtanh(DD s) {
  DD z = mul(s, s);
  DD p = kInvOdd[20];
  for (int i = 19; i >= 0; --i) p = add(kInvOdd[i], mul(z, p));
  return ldexp(mul(s, p), 1);
}

/// log of y = 2^k m, m in [sqrt(2)/2, sqrt(2)), given y > 0 finite and
/// normal. When k == 0 and ym1 (y - 1, exact) is given, s = (y - 1)/(y + 1)
/// is formed from it directly.
inline DD logCore(DD y, const DD* ym1) {
  int k;
  double mh = std::frexp(y.hi, &k);  // y.hi = mh 2^k, mh in [0.5, 1)
  if (mh < 0x1.6a09e667f3bcdp-1) {
    mh *= 2.0;
    --k;
  }
  DD m = ldexp(y, -k);
  DD mm1;
  if (k == 0 && ym1 != nullptr) {
    mm1 = *ym1;
  } else {
    // m.hi in [0.70, 1.42): m.hi - 1 is exact (Sterbenz), then the exact
    // sum with m.lo
    mm1 = twoSum(m.hi - 1.0, m.lo);
  }
  DD mp1 = add(m, 1.0);
  DD s = div(mm1, mp1);
  DD a = twoAtanh(s);
  if (k == 0) return a;
  return add(mul(kLn2, static_cast<double>(k)), a);
}

}  // namespace detail

/// exp(x) for x in [-650, 700]. Reduction: k = nearest(x / ln 2),
/// r = x - k ln 2 with ln 2 in three parts (k * kLn2Hi30 exact, x.hi minus it
/// exact by Sterbenz), absolute error at most 5 u^2 |r| + 2^-125; then
/// 2^k (1 + expm1(r)). The sum 1 + q, |q| / (1 + q) <= 0.42, carries 0.42 of
/// expm1's 30.3 u^2 (12.8 u^2), r's error (at most 1.8 u^2 relative at
/// |r| <= 0.3466) and its own 2 u^2: at most 16.6 u^2, stated as 24 u^2.
inline constexpr double kExpRelErr = 0x1.8p-102;  // 24 u^2
inline DD exp(DD x) {
  double k = std::nearbyint(x.hi * 0x1.71547652b82fep+0);
  double a = x.hi - k * kLn2Hi30;
  DD r = twoSum(a, x.lo);
  r = sub(r, twoProd(k, kLn2Mid));
  r = sub(r, k * kLn2Lo);
  DD q = detail::expm1Reduced(r);
  DD e = add(q, 1.0);
  return ldexp(e, static_cast<int>(k));
}

/// log1p(t) for t = 0 or t > -1 finite with |t| >= 2^-960 (below that the
/// low words of t / (2 + t) and its square go subnormal and the bound
/// fails; the callers pass e^{-2w} >= e^-650 > 2^-938). With y = 1 + t (one
/// DWPlusFP, 2 u^2) written 2^k m: for k == 0 (t in [-0.293, 0.414))
/// s = t / (2 + t) from t directly, relative error at most 17 u^2, and
/// 2 atanh(s) carries 1.0202 (17 u^2) + 8.6 u^2 <= 26 u^2 in all; for
/// k != 0, |k ln 2| >= 0.693 >= 2 |2 atanh s|, so the sum's relative error
/// is at most 2 (3 u^2) + 26 u^2 + 2 u^2 / 0.3466 + 3 u^2, about 41 u^2.
inline constexpr double kLog1pRelErr = 0x1.8p-101;  // 48 u^2
inline DD log1p(DD t) {
  DD y = add(t, 1.0);
  int k;
  double mh = std::frexp(y.hi, &k);
  if (mh < 0x1.6a09e667f3bcdp-1) --k;
  if (k == 0) {
    DD s = div(t, add(t, 2.0));
    return detail::twoAtanh(s);
  }
  return detail::logCore(y, nullptr);
}

/// log(x) for x in [2^-960, 2^1023] (the low words stay normal through the
/// reduction): the same reduction with y = x and,
/// for k == 0, x - 1 exact; relative error as log1p's.
inline constexpr double kLogRelErr = 0x1.8p-101;  // 48 u^2
inline DD log(DD x) {
  if (x.hi > 0.5 && x.hi < 2.0) {
    DD xm1 = twoSum(x.hi - 1.0, x.lo);
    return detail::logCore(x, &xm1);
  }
  return detail::logCore(x, nullptr);
}

/// sin(x) for x in [2^-400, pi/2]: x P(z), z = x^2, P the Taylor series of
/// 18 terms (truncation below 2^-111 relative at x = pi/2). z = x * x
/// carries 4 u^2 relative, which reaches P with factor
/// |z P'(z) / P(z)| = (1 - x cot x) / 2 <= 1/2: 2 u^2. The alternating Horner
/// sum adds at most 7.4 u^2 (the first step's |z P_1 / 6| is at most 0.57 of
/// P_0), the coefficients' rounding (2^-107 relative each) at most 1.1 u^2,
/// and the product by x 4 u^2: at most 14.5 u^2.
inline constexpr double kSinRelErr = 0x1p-102;  // 16 u^2
inline DD sin(DD x) {
  DD z = mul(x, x);
  DD p = detail::kSinCoef[17];
  for (int i = 16; i >= 0; --i) p = add(detail::kSinCoef[i], mul(z, p));
  return mul(x, p);
}

/// sqrt(x) for x in [2^-960, 2^1000] (below that s^2's low word goes
/// subnormal): one Newton step from the double square root s, x - s^2 formed
/// exactly to the rounding of two additions. s is within 1.5 u of sqrt(x)
/// (its own rounding and x.lo's u / 2), so the step's quadratic term is at
/// most (1.5 u)^2 / 2 <= 1.13 u^2 relative, and the correction's rounding adds
/// 3.2 u^2.
inline constexpr double kSqrtRelErr = 0x1p-103;  // 8 u^2
inline DD sqrt(DD x) {
  double s = std::sqrt(x.hi);
  DD p = twoProd(s, s);
  double r = ((x.hi - p.hi) - p.lo) + x.lo;
  double c = r / (2.0 * s);
  return fastTwoSum(s, c);
}

/// The largest relative error any of the five functions may carry (the
/// density evaluator's assumption): every stated bound must be at most this.
inline constexpr double kFunctionRelErrMax = 0x1p-100;

}  // namespace ddouble
}  // namespace bartcore

#endif  // BARTCORE_DDOUBLE_HPP
