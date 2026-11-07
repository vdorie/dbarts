"""Double-precision reference implementation of the certified density rules
of proof.tex (Section 5, Appendices A and B), used to check the bounds against
mpmath references. Not the engine; a specification by example.

A family is the distribution of X = sum_j Ga(a_j, 2 tau_j^2) (finitely or
infinitely many terms), Laplace transform L(t) = prod (1 + t/(2 tau_j^2))^-a_j.
  F0(b):   PG(b, 0), tau_k = pi (k - 1/2), a_k = b, L(t) = cosh^-b sqrt(t/2).
  QT(...): finitely many terms (the untilted proposal).
Contour (Lemma 5.1): with a shift s in [0, tau_min^2) and rho > 0,
  p(x) = (2/pi) e^{2 x r2} L(2 r2) * S,   r2 = rho^2 - s,
  S = int_R Re[e^{2x(z^2 - rho^2)} Lt(z)/Lt(rho) z] dy,  z = rho + i y,
  Lt(z) = L(2 z^2 - 2 s), singular only at z = +-i tau_j', tau_j'^2 = tau_j^2 - s.
Every evaluation returns log p_hat and a bound on |log p_hat - log p|.

Floating-point model (Appendix B): +,-,*,/ and real sqrt are correctly
rounded (relative error <= U = 2^-53, IEEE 754 binary64); the library
functions exp, log, log1p, expm1, sin, cos, tan, sinh, cosh, tanh, atan,
atan2, erfc and the complex sqrt return results with relative error at most
ULIB = 2U (absolute error ULIB for sin and cos). Majorants (M, tails) are
computed in floating point and then inflated by INFL, which dominates the
relative rounding of those short computations.
"""
import math, cmath

U = 2.0 ** -53
PI = math.pi
ULIB = 2.0 * U
INFL = 1.0 + 2.0 ** -30          # inflation of computed majorants
DEFL = 1.0 - 2.0 ** -30          # deflation of computed minorants
TAU1SQ_LO = 2.4674011002723395   # (pi/2)^2 rounded down, for the shift check


LOG2 = math.log(2.0)


def lcosh(p):
    """log cosh p, p >= 0, without overflow."""
    return p + math.log1p(math.exp(-2 * p)) - LOG2


def lsinh(p):
    """log sinh p, p > 0, without overflow."""
    if p > 20:
        return p - LOG2 + math.log1p(-math.exp(-2 * p))
    return math.log(math.sinh(p))


def logcosh_pq(p, q):
    """log cosh(p + i q), p >= 0, on the branch continuous from the real axis:
    re = 0.5 log(sinh^2 p + cos^2 q), im = q + atan(N/D) with
    N = (tanh p - 1) sin q cos q, D = cos^2 q + tanh p sin^2 q > 0.
    Returns (value, bound on |computed - exact at (p, q)|, m2 = |cosh w|^2).
    Appendix B, Lemma B.2."""
    if p > 20:
        # log cosh w = w - log 2 + log(1 + e^{-2w}), |e^{-2w}| = e^{-2p} < 5e-18
        e = math.exp(-2 * p)
        c2, s2 = math.cos(2 * q), math.sin(2 * q)
        re = p - LOG2 + 0.5 * math.log1p(2 * e * c2 + e * e)
        im = q + math.atan(-e * s2 / (1 + e * c2))
        err = U * (abs(re) + p) + ULIB + 6 * ULIB * e + U * abs(im) + ULIB * 4 * e + ULIB
        return complex(re, im), err * INFL, re
    sp, cq, sq = math.sinh(p), math.cos(q), math.sin(q)
    m2 = sp * sp + cq * cq
    re = 0.5 * math.log(m2)
    tp = math.tanh(p)
    tm1 = -2.0 / (math.exp(2 * p) + 1.0)                 # tanh p - 1
    N = tm1 * sq * cq
    D = cq * cq + tp * sq * sq
    r = N / D
    im = q + math.atan(r)
    dm2 = 3 * ULIB * (sp * sp + abs(cq) + m2)
    err_re = (0.5 * dm2 / (m2 - dm2) if m2 > dm2 else math.inf) + ULIB * abs(re)
    dN = abs(tm1) * 2.01 * ULIB + abs(N) * 4 * ULIB
    dD = 2.01 * ULIB * (abs(cq) + tp * abs(sq)) + D * 3 * ULIB
    dr = ((dN + abs(r) * dD) / (D - dD) if D > dD else math.inf) + U * abs(r)
    err_im = dr + ULIB * abs(math.atan(r)) + U * abs(im)
    return complex(re, im), (err_re + err_im) * INFL, re


class F0:
    """PG(b, 0); contour (rho, s)."""
    def __init__(self, b):
        self.b = float(b)
        self.mean = self.b / 4.0
        self.tmin = PI * PI / 4

    def set_contour(self, rho, s):
        assert rho > 0 and 0 <= s < TAU1SQ_LO * DEFL
        self.rho, self.s = rho, s
        r2 = rho * rho - s
        self.r2 = r2
        self.dr2 = 2.2 * U * (rho * rho + abs(r2))
        if r2 >= 0:
            v = math.sqrt(r2)
            lc0 = v + math.log1p(math.exp(-2 * v)) - math.log(2.0)
            dlc0 = 6 * ULIB * (abs(lc0) + v + 1.0) + self.dr2 * 0.5
        else:
            v = math.sqrt(-r2)
            lc0 = math.log(math.cos(v))
            dlc0 = 6 * ULIB * (abs(lc0) + 1.0) / math.cos(v) + self.dr2 * 0.51 * math.tan(v) / v
        self.lc0, self.dlc0 = lc0, dlc0 * INFL

    def shape_total(self):
        return self.b

    def log_L_vertex(self):
        """log L(2 r2) = -b log cosh sqrt(r2), with an absolute error bound."""
        return -self.b * self.lc0, U * self.b * abs(self.lc0)

    def _w(self, u, y):
        W = complex((u * u - self.s) - y * y, 2 * u * y)
        dW = 3.3 * U * (u * u + self.s + y * y + u * abs(y))
        return W, dW, cmath.sqrt(W)

    def log_ratio(self, y):
        """log(Lt(rho + i y)/Lt(rho)), y >= 0, with an absolute error bound."""
        W, dW, w = self._w(self.rho, y)
        p, q = abs(w.real), w.imag
        lc, elc, lm = logcosh_pq(p, q)          # lm = log |cosh w| (computed)
        lm_lo = lm - elc - 1e-15
        absW = abs(W)
        # (i) rounding of W: Phi(W) = log cosh sqrt W is analytic in W,
        # |Phi'(W')| = |tanh w'/(2 w')| over the disk |W' - W| <= dW.
        R = math.sqrt(absW + dW)
        if R <= 1.0:
            dphi = math.tan(R) / (2 * R) if R > 1e-8 else 0.5000001
        else:
            # a branch of sqrt is analytic on the disk (it misses 0), so
            # |w' - w| <= dw = dW / (2 sqrt(|W| - dW)) for the branch near w;
            # |tanh w'| <= cosh(p + dw) / (|cosh w| - dw cosh(p + dw))
            rw = math.sqrt(absW - dW)
            dw = dW / (2 * rw)
            rat = math.exp(lm_lo - lcosh(p + dw))
            dphi = 1.0 / (2 * rw * (rat - dw)) if rat > dw else math.inf
        # (ii) the complex sqrt: |w_hat - sqrt W_hat| <= 4 ULIB |w|
        dl = 4 * ULIB * abs(w)
        rat = math.exp(lm_lo - lcosh(p + dl))
        dth = 1.0 / (rat - dl) if rat > dl else math.inf
        err = self.b * (dW * dphi + dl * dth + elc) * INFL
        val = -self.b * (lc - self.lc0)
        return val, err + 2 * U * abs(val)

    def _pq(self, u, y):
        """p = Re w, q = Im w at (u, y) and an error allowance."""
        W, dW, w = self._w(u, y)
        return abs(w.real), w.imag, 1e-12 * abs(w) + math.sqrt(dW)

    def abs_L_cell_log(self, u, ya, yb):
        """upper bound on max_{ya<=y<=yb} log|Lt(u+iy)| - log Lt(rho), from
        |cosh w|^2 = sinh^2 p + cos^2 q with p, q nondecreasing in y."""
        pa, qa, ea = self._pq(u, ya)
        pb, qb, eb = self._pq(u, yb)
        plo = max(pa - ea, 0.0)
        qlo, qhi = qa - ea, qb + eb
        j0 = math.floor(qlo / PI - 0.5)
        if qhi - qlo >= PI or any(qlo <= (j + 0.5) * PI <= qhi for j in range(j0, j0 + 3)):
            c2 = 0.0
        else:
            c2 = min(math.cos(qlo) ** 2, math.cos(qhi) ** 2) * DEFL
        # log(sinh^2 p + c2), overflow-safe
        if plo > 0:
            ls2 = 2 * lsinh(plo) + math.log(DEFL)
            lm2 = ls2 + math.log1p(c2 * math.exp(-ls2)) if ls2 < 700 else ls2
        else:
            lm2 = math.log(c2) if c2 > 0 else -math.inf
        # coupled bound: |cosh w|^2 = cosh^2 p - sin^2 q >= 1 + p^2 - q^2
        # = 1 + Re W = 1 + u^2 - s - y^2 (cosh^2 p >= 1 + p^2, sin^2 q <= q^2)
        m2c = 1.0 + (u * u - self.s) - yb * yb - 8 * U * (1 + u * u + self.s + yb * yb)
        if m2c > 0:
            lm2 = max(lm2, math.log(m2c))
        if lm2 == -math.inf:
            return math.inf
        return -0.5 * self.b * lm2 + self.b * self.lc0 + 1e-12

    def abs_L_tail_log(self, u, Y):
        """upper bound on sup_{y >= Y} log|Lt(u+iy)| - log Lt(rho), from
        |cosh w| >= sinh p and p nondecreasing in y."""
        pa, qa, ea = self._pq(u, Y)
        plo = max(pa - ea, 0.0)
        if plo <= 0:
            return math.inf
        return -self.b * (lsinh(plo) + math.log(DEFL)) + self.b * self.lc0 + 1e-12

    def saddle_r2(self, x):
        """r2 in (-pi^2/4, inf) with x = b tanh(v)/(4 v), v = sqrt(r2)."""
        target = 4 * x / self.b
        def g(r2):
            if r2 > 1e-12:
                v = math.sqrt(r2); return math.tanh(v) / v
            if r2 < -1e-12:
                v = math.sqrt(-r2); return math.tan(v) / v
            return 1.0
        lo, hi = -self.tmin * (1 - 1e-15), 1.0
        while g(hi) > target:
            hi *= 4
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if g(mid) > target:
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi)


class QT:
    """sum_j Ga(a_j, 2 tau2_j), tau2_j known to relative error rel_j; the a_j
    are the exact doubles the sampler uses. Contour (rho, s)."""
    def __init__(self, tau2, a, rel=None):
        self.tau2 = [float(t) for t in tau2]
        self.a = [float(v) for v in a]
        self.rel = rel if rel is not None else [2 * U] * len(self.tau2)
        self.mean = sum(av / (2 * t) for av, t in zip(self.a, self.tau2))
        self.tmin = min(self.tau2)

    @classmethod
    def proposal(cls, b, K, alpha, beta_t):
        """K exact terms Ga(b, c_k), c_k = 2 pi^2 (k-1/2)^2, and Ga(alpha, beta_t)."""
        tau2 = [PI * PI * (k - 0.5) ** 2 for k in range(1, K + 1)] + [beta_t / 2.0]
        return cls(tau2, [b] * K + [alpha], [4 * U] * K + [2 * U])

    def set_contour(self, rho, s):
        assert rho > 0 and 0 <= s < self.tmin * DEFL
        self.rho, self.s = rho, s
        self.r2 = rho * rho - s
        self.dr2 = 2.2 * U * (rho * rho + abs(self.r2))
        self.tp = [t - s for t in self.tau2]
        self.dtp = [(r + U) * (t + s) * 1.01 for r, t in zip(self.rel, self.tau2)]
        self.den = [rho * rho + t for t in self.tp]

    def shape_total(self):
        # effective curvature weight (cell counts only; any value is valid)
        return sum(av * min(1.0, 2.5 / den) for av, den in zip(self.a, self.den))

    def log_L_vertex(self):
        """log L(2 r2) = -sum a log1p(r2/tau2), with an absolute error bound."""
        val, err = 0.0, 0.0
        for av, t, rl in zip(self.a, self.tau2, self.rel):
            v = self.r2 / t
            l = math.log1p(v)
            val -= av * l
            dv = (abs(v) * (2 * U + rl) + self.dr2 / t) * 1.01
            err += av * (ULIB * abs(l) + dv / (1 + v - dv)) + 2 * U * abs(val)
        return val, err * INFL

    def log_ratio(self, y):
        rho = self.rho
        tot = 0j
        err = 0.0
        for av, tp, dtp, den in zip(self.a, self.tp, self.dtp, self.den):
            v = complex(-y * y, 2 * rho * y) / den
            re = 0.5 * math.log1p(2 * v.real + v.real * v.real + v.imag * v.imag)
            im = math.atan2(v.imag, 1.0 + v.real)
            tot += -av * complex(re, im)
            m = abs(1 + v)           # |z^2 + tau'^2| / (rho^2 + tau'^2)
            e_eval = (8 * U * abs(v) * (1 + abs(v)) / (m * m) + ULIB * abs(re)
                      + ULIB * abs(im) + 4 * U * abs(v) / m + 2 * U)
            e_tau = dtp * (1 / (m * den) + 1 / den) + 2 * U
            err += (av * (e_eval + e_tau) + U * av * (abs(re) + abs(im))) * INFL + 2 * U * abs(tot)
        return tot, err

    def abs_L_cell_log(self, u, ya, yb):
        out = 0.0
        for av, tp, den in zip(self.a, self.tp, self.den):
            # |z^2 + tau'^2|^2 = (A - Y)^2 + 4 u^2 Y, A = u^2 + tau'^2, Y = y^2:
            # convex in Y, least at Y = tau'^2 - u^2; clamp to the cell
            A = u * u + tp
            Yc = min(max(tp - u * u, ya * ya), yb * yb)
            m = ((A - Yc) ** 2 + 4 * u * u * Yc) * DEFL
            out += -0.5 * av * math.log(m) + av * math.log(den * INFL)
        return out * (1 + 1e-14) + 1e-12

    def abs_L_tail_log(self, u, Y):
        out = 0.0
        for av, tp, den in zip(self.a, self.tp, self.den):
            t = math.sqrt(tp)
            if Y >= t:
                m = (u * u + (Y - t) ** 2) * (u * u + (Y + t) ** 2)
            else:
                m = 4 * tp * u * u if u * u < tp else (u * u + tp) ** 2
            out += -0.5 * av * math.log(m * DEFL) + av * math.log(den * INFL)
        return out * (1 + 1e-14) + 1e-12

    def saddle_r2(self, x):
        """r2 in (-tmin, inf) with 2x = sum a_j/(r2 + tau2_j)."""
        g = lambda r2: sum(av / (r2 + t) for av, t in zip(self.a, self.tau2))
        lo, hi = -self.tmin * (1 - 1e-15), 1.0
        while g(hi) > 2 * x:
            hi *= 4
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if g(mid) > 2 * x:
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi)


def choose_contour(fam, x, c=0.5):
    """Vertex at the real saddle t* = 2 r2; rho^2 = max(r2, c (r2 + tmin)),
    s = rho^2 - r2. Any (rho, s) is valid; this choice is about cost only."""
    r2 = fam.saddle_r2(x)
    D = r2 + fam.tmin
    if r2 >= c * D:
        return math.sqrt(r2), 0.0
    rho = math.sqrt(c * D)
    s = max(rho * rho - r2, 0.0)
    s = min(s, fam.tmin * (1 - 1e-6))
    return rho, s


def _ratio(logM, sexp):
    """exp(logM) / (exp(sexp) - 1) without overflow."""
    v = logM - sexp - math.log1p(-math.exp(-sexp))
    return math.exp(v) if v < 700 else math.inf


def sexp(v):
    """exp that saturates to inf instead of raising."""
    return math.exp(v) if v < 709 else math.inf


def logsumexp(vals):
    m = max(vals)
    if m == -math.inf or m == math.inf:
        return m
    return m + math.log(sum(math.exp(v - m) for v in vals))


def gauss_tail_log(u, x, Y):
    """log of an upper bound on int_Y^inf e^{-2 x y^2} (u + y) dy. For Y > 0,
    erfc(a) <= e^{-a^2}/(a sqrt(pi)) gives e^{-2xY^2} (u/Y + 1)/(4x)."""
    if Y > 0:
        return -2 * x * Y * Y + math.log((u / Y + 1.0) / (4 * x)) + 1e-12
    return math.log(u * math.sqrt(PI / (8 * x)) + 1.0 / (4 * x)) + 1e-12


def line_l1_log(fam, x, u, Y, ncell):
    """log of an upper bound on int_R e^{2x(u^2 - y^2 - rho^2)} |Lt(u+iy)|/Lt(rho)
    |u + i y| dy on the single line Re z = u (Lemma 5.3 needs only the two
    edge lines of the strip): Darboux cells on [0, Y], a closed-form tail."""
    rho = fam.rho
    pref = 2 * x * (u * u - rho * rho)
    terms = []
    # extend Y until the Gaussian beats the global sup of |Lt| on the line
    # (|Lt| need not decay in y: |cosh(u+iy)| is periodic in y when s = 0)
    lsup0 = fam.abs_L_tail_log(u, Y)        # nonincreasing in Y
    if 0 < lsup0 < math.inf:
        Y = max(Y, math.sqrt((lsup0 + 60.0) / (2 * x)))
    ncell = int(min(40000, max(ncell, 40 * Y * math.sqrt(x + fam.shape_total()))))
    dy = Y / ncell
    for j in range(ncell):
        ya, yb = j * dy, (j + 1) * dy
        terms.append(-2 * x * ya * ya + fam.abs_L_cell_log(u, ya, yb)
                     + 0.5 * math.log(u * u + yb * yb) + math.log(dy))
    terms.append(fam.abs_L_tail_log(u, Y) + gauss_tail_log(u, x, Y))
    return (math.log(2.0) + pref + logsumexp(terms)) * (1 + 1e-14) + 1e-12


def zline(fam, x, rho, s, h, theta=0.5, verbose=False, cache=None):
    """Certified evaluation of log p(x) on the contour (rho, s), step h.
    h has at most 2 significant bits, so the nodes k h are exact."""
    fam.set_contour(rho, s)

    def g(y):
        lr, elr = fam.log_ratio(y)
        E = complex(-2 * x * y * y, 4 * x * rho * y) + lr
        val = cmath.exp(E) * complex(rho, y)
        dE = 6 * U * 2 * x * (y * y + 2 * rho * y) + elr + 2 * U * abs(E)
        return val, dE

    S, absum, errround = 0.0, 0.0, 0.0
    k = 0
    while True:
        y = k * h
        v, dE = g(y)
        wgt = 1.0 if k == 0 else 2.0
        S += wgt * v.real
        absum += wgt * abs(v)
        errround += wgt * abs(v) * ((math.expm1(dE) if dE < 1 else math.inf) + 8 * ULIB) * 1.01
        k += 1
        if k > 4 and abs(v) < 1e-18 * absum and y > 1.0 / math.sqrt(x):
            break
        if k > 2_000_000:
            return dict(logp=None, err=math.inf, N=k, rho=rho, s=s, h=h)
    N = k - 1
    S *= h
    absum *= h
    errround = errround * h + (N + 4) * U * absum
    # truncation: cell majorants on [N h, KD h], then the Gaussian tail
    KD = 4 * N + 50
    lo, YD = N * h, KD * h
    nc = 4 * N + 50
    dy = (YD - lo) / nc
    tl = [-2 * x * (lo + j * dy) ** 2 + fam.abs_L_cell_log(rho, lo + j * dy, lo + (j + 1) * dy)
          + 0.5 * math.log(rho * rho + (lo + (j + 1) * dy) ** 2) for j in range(nc)]
    tail = 2 * h * (dy / h + 1) * sexp(logsumexp(tl)) * INFL
    yW = 0.5 * (-rho + math.sqrt(rho * rho + 1.0 / x))
    YD = max(YD, yW)
    lsup0 = fam.abs_L_tail_log(rho, YD)
    if 0 < lsup0 < math.inf and 2 * x * YD * YD < lsup0 + 60.0:
        # cells up to where the Gaussian beats the global sup
        Y2 = math.sqrt((lsup0 + 60.0) / (2 * x))
        nc2 = int(min(40000, 40 * (Y2 - YD) * math.sqrt(x + fam.shape_total())) + 10)
        dy2 = (Y2 - YD) / nc2
        tl2 = [-2 * x * (YD + j * dy2) ** 2 + fam.abs_L_cell_log(rho, YD + j * dy2, YD + (j + 1) * dy2)
               + 0.5 * math.log(rho * rho + (YD + (j + 1) * dy2) ** 2) for j in range(nc2)]
        tail += 2 * h * (dy2 / h + 1) * sexp(logsumexp(tl2)) * INFL
        YD = Y2
    lsup = fam.abs_L_tail_log(rho, YD)
    lt = lsup + logsumexp([gauss_tail_log(rho, x, YD),
                           math.log(h * (rho + YD)) - 2 * x * YD * YD])
    tail += 2 * math.exp(lt) * INFL if lt < 700 else math.inf
    # discretization (TW Thm 5.1, unsymmetric strip)
    Ydb = max(N * h, 8.0 / math.sqrt(x))
    ncell = int(min(20000, max(400, 40 * Ydb * math.sqrt(x + fam.shape_total()))))
    # Lemma 5.3 (trapezoid error from the two edge lines of the strip)
    # holds for every a+ in (0, rho) and a- > 0: take the best of a few.
    fr = [0.5, 0.25, 0.1, 0.04, 0.015, 0.005]
    if cache is None:
        cache = {}
    def L1(u):
        if u not in cache:
            cache[u] = line_l1_log(fam, x, u, Ydb, ncell)
        return cache[u]
    bp = min(_ratio(L1(rho - f * rho), 2 * PI * f * rho / h) for f in fr if f <= theta)
    bm = min(_ratio(L1(rho + am), 2 * PI * am / h)
             for am in [f * rho for f in [2.0, 1.0] + fr] + [0.25, 0.5, 1.0, 2.0])
    disc = bp + bm
    E = disc + tail + errround
    if not S > E:
        return dict(logp=None, err=math.inf, N=N, rho=rho, s=s, h=h, S=S, E=E)
    lLv, elLv = fam.log_L_vertex()
    r2 = fam.r2
    lead = math.log(2 / PI) + 2 * x * r2 + lLv
    dlead = (ULIB * 0.5 + 2 * x * fam.dr2 + 4 * U * 2 * x * abs(r2) + elLv
             + 3 * U * abs(lead)) * INFL
    logp = lead + math.log(S)
    err = -math.log1p(-E / S) + dlead + ULIB * abs(math.log(S)) + 2 * U * abs(logp)
    out = dict(logp=logp, err=err, N=N, rho=rho, s=s, h=h, disc=disc / S,
               tail=tail / S, rnd=errround / S, cond=absum / S)
    if verbose:
        print(out)
    return out


H_SEQ = [m * 2.0 ** -j for j in range(0, 40) for m in (1.0, 0.75)]


def auto(fam, x, target=1e-13, c=0.5):
    """Walk h down H_SEQ until the bound meets the target; stop once the
    bound stops improving. Returns the best result."""
    rho, s = choose_contour(fam, x, c)
    best, worse = None, 0
    cache = {}
    for h in H_SEQ:
        r = zline(fam, x, rho, s, h, cache=cache)
        if best is None or r['err'] < 0.9 * best['err']:
            best, worse = r, 0
        elif best['err'] < 1e-6:
            worse += 1
            if worse >= 2:
                break
        if r['logp'] is not None and r['err'] <= target:
            return r
    return best


def left_series(b, x, target=1e-15):
    """Left alternating series for PG(b,0) at x. Returns (log f_hat, err, n)
    or None if the start condition (every term ratio <= 1/2) fails."""
    b = float(b)
    Cb = (b + 2.0) if b >= 1 else (b + 2.0) / b
    if Cb * math.exp(-(b + 1) / (2 * x)) > 0.5 * DEFL:
        return None
    s, t, absum, n = 1.0, 1.0, 1.0, 0
    while True:
        r = (n + b) / (n + 1) * (2 * n + b + 2) / (2 * n + b) * math.exp(-(2 * n + b + 1) / (2 * x))
        t *= r
        n += 1
        if t <= target * 0.25 * s:
            trunc = t * 1.01
            break
        s += (-1) ** n * t
        absum += t
    loga0 = ((b - 1) * math.log(2.0) + math.log(b) - 1.5 * math.log(x) - b * b / (8 * x)
             - 0.5 * math.log(2 * PI))
    rel_term = n * (8 * U + ULIB * (1 + (2 * n + b + 1) / (2 * x)))
    rnd = (rel_term + (n + 2) * U) * absum * 1.01
    E = trunc + rnd
    da0 = 4 * ULIB * (abs((b - 1) * math.log(2.0)) + abs(math.log(b)) + abs(1.5 * math.log(x))
                      + b * b / (8 * x) + 2.0)
    logp = loga0 + math.log(s)
    return logp, -math.log1p(-E / s) + da0 + ULIB * abs(math.log(s)) + 2 * U * abs(logp), n
