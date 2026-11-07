"""Certified log densities for the accept-or-keep Polya-Gamma step, in ball
arithmetic (python-flint arb). Every quantity is an arb/acb ball; arb's
contract (each operation returns a ball containing the exact result for every
point of its input balls, elementary functions included) makes every
enclosure rigorous with no assumption about a math library. The mathematical
bounds (discretization, truncation) are those of proof.tex, Section 4.

Families: X = sum_j Ga(a_j, lam_j), Laplace transform L(t) = prod (1+t/lam_j)^-a_j.
  PGFam(b, psi): PG(b, psi), lam_k = 2 pi^2 (k-1/2)^2 + psi^2/2,
                 L(t) = (cosh A / cosh sqrt(A^2 + t/2))^b, A = |psi|/2.
                 PGFam(b, 0) is f0 (the untilted target of the old note).
  GamFam(a, lam): a finite sum of gammas (the proposal q, tilted: the
                 sampler's own rates; or the old untilted q-tilde).
Contour (Lemma 4.1 in the tilted frame): vertex t0 > -lam_1, rho > 0 with
2 rho^2 < lam_1 + t0; z = rho + i y, t(y) = t0 + 2 (z^2 - rho^2),
  p(x) = (2/pi) e^{x t0} L(t0) S,  S = int_R Re[e^{x (t - t0)} L(t)/L(t0) z] dy.
Formulations:
  'new':  the node exponent x (t - t0) + log(L(t)/L(t0)) is formed from
          Dh = z^2 - rho^2 = -y^2 + 2 i rho y (exact up to rounding of y^2, rho y)
          and, for PG, D = log cosh w - log cosh w0 by the cancellation-free
          D = log1p(expm1(delta) (1 - e^{-(w+w0)}) / (1 + e^{-2 w0})),
          delta = Dh / (w + w0), the branch fixed by the direct difference.
  'old':  D as the difference of two separately rounded log cosh values
          (the formula of the first note), for comparing floors.
The density is evaluated at the double x; inputs b, psi, a_j, lam_j, x, t0,
rho, h are doubles and enter exactly.
"""
import math
from flint import arb, acb, ctx


def setprec(p):
    ctx.prec = p


def _lc(w):
    """log cosh w on the branch analytic on Re w > 0 (and at w = i v, |v| < pi/2)."""
    return w - arb(2).log() + (-2 * w).exp().log1p()


class PGFam:
    def __init__(self, b, psi, mode='new'):
        self.b = float(b)
        self.psi = abs(float(psi))
        self.mode = mode
        self.tag = 'PG'

    # ---- double-precision helpers (choices only; any choice is valid) ----
    def lam1_d(self):
        return math.pi ** 2 / 2 + self.psi ** 2 / 2

    def mean_d(self):
        a = self.psi / 2
        return self.b * (math.tanh(a) / (4 * a) if a > 1e-8 else 0.25)

    def tilt_mean_d(self, t):
        """E X under the exponential tilt by t: b tanh(w)/(4w), w^2 = A^2 + t/2."""
        W = (self.psi / 2) ** 2 + t / 2
        if W > 1e-14:
            w = math.sqrt(W)
            return self.b * math.tanh(w) / (4 * w)
        if W < -1e-14:
            v = math.sqrt(-W)
            return self.b * math.tan(v) / (4 * v)
        return self.b / 4

    def shape_d(self):
        return self.b

    # ---- arb pieces ----
    def lam1(self):
        return arb.pi() ** 2 / 2 + arb(self.psi) ** 2 / 2

    def _w0(self, W0):
        if W0 > 0:
            return acb(W0.sqrt())
        return acb(0, (-W0).sqrt())

    def _coshw0(self, W0):
        """cosh sqrt(W0) > 0 (W0 > -pi^2/4)."""
        if W0 > 0:
            return W0.sqrt().cosh()
        return (-W0).sqrt().cos()

    def prepare(self, x, t0, rho):
        A = arb(self.psi) / 2
        self.A2 = A * A
        self.t0 = arb(t0)
        self.rho = arb(rho)
        self.W0 = self.A2 + self.t0 / 2
        if self.W0.contains(0):
            raise ValueError('W0 straddles 0; nudge t0')
        self.w0 = self._w0(self.W0)
        self.cw0 = self._coshw0(self.W0)
        self.e2w0 = (-2 * self.w0).exp()
        self.lc0 = _lc(self.w0)

    def D(self, Dh):
        """log cosh sqrt(W0 + Dh) - log cosh sqrt(W0), the branch continuous
        along the contour from Dh = 0 (Im Dh > 0)."""
        W = self.W0 + Dh
        w = W.sqrt()
        Ddir = _lc(w) - self.lc0
        if self.mode == 'old':
            return Ddir
        s = w + self.w0
        delta = Dh / s
        v = delta.expm1() * (1 - (-s).exp()) / (1 + self.e2w0)
        Dcf = v.log1p()
        if not Dcf.is_finite():
            return Ddir
        twopi = 2 * arb.pi()
        k = round(float((Ddir - Dcf).imag.mid()) / float(twopi.mid()))
        Dfix = Dcf + acb(0, twopi * k)
        if not (Ddir - Dfix).abs_upper() < 3:
            raise ValueError('branch not resolved')
        return Dfix if Dfix.rad() < Ddir.rad() else Ddir

    def log_ratio(self, Dh):
        return -self.b * self.D(Dh)

    def log_L_t0(self):
        """log L(t0) = -b [log cosh w0 - log cosh A], real path."""
        if self.t0 == 0:
            return arb(0)
        A = arb(self.psi) / 2
        if self.mode == 'old':
            return -self.b * (self.lc0.real - _lc(acb(A)).real)
        s = self.w0 + A
        delta = acb(self.t0 / 2) / s
        v = delta.expm1() * (1 - (-s).exp()) / (1 + (-2 * A).exp())
        Dcf = v.log1p().real
        Ddir = self.lc0.real - _lc(acb(A)).real
        return -self.b * (Dcf if Dcf.rad() < Ddir.rad() else Ddir)

    # majorants: |L(t)/L(t0)| on the line Re z = u
    def _pq(self, u, y):
        """w = p + i q = sqrt(z^2 - s'), z = u + i y, s' = rho^2 - W0 (principal)."""
        W = acb(self.W0 + u * u - self.rho * self.rho - y * y, 2 * u * y)
        if y == 0 and W.real < 0:
            return arb(0), (-W.real).sqrt(), W
        w = W.sqrt()
        return w.real, w.imag, W

    def cell_log_major(self, u, ya, yb):
        """Lemma 4.3(a) extended: on Re z = u > 0, q = Im w is nondecreasing in
        y, and p = Re w is monotone (nondecreasing if s' >= 0, nonincreasing if
        s' <= 0), so on [ya, yb]: p >= min(p(ya), p(yb)), q in [q(ya), q(yb)];
        |cosh w|^2 = sinh^2 p + cos^2 q and |cosh w|^2 >= 1 + Re W."""
        u = arb(u)
        pa, qa, Wa = self._pq(u, arb(ya))
        pb, qb, Wb = self._pq(u, arb(yb))
        plo = pa.min(pb)
        qlo, qhi = qa.lower(), qb.upper()
        pi = arb.pi()
        j0 = int(math.floor(float(qlo) / math.pi - 0.5)) - 1
        hit = any(((j + 0.5) * pi).overlaps(arb((qlo + qhi) / 2, (qhi - qlo) / 2) if True else 0)
                  or (qlo <= (j + 0.5) * pi <= qhi) for j in range(j0, j0 + 4))
        if float(qhi - qlo) >= math.pi:
            hit = True
        c2 = arb(0) if hit else (qa.cos() ** 2).min(qb.cos() ** 2).lower()
        if c2 < 0:
            c2 = arb(0)
        m2 = plo.sinh() ** 2 + c2 if plo > 0 else c2
        m2c = 1 + Wb.real
        if m2c > 0:
            m2 = m2.max(m2c)
        m2 = m2.lower()
        if not m2 > 0:
            return arb('inf')
        return self.b * (self.cw0.log() - m2.log() / 2)

    def tail_log_major(self, u, Y):
        """sup over y >= Y on the line Re z = u: |cosh w| >= sinh p, with
        p >= u if s' <= 0, p >= u Y / sqrt(Y^2 + s') if s' > 0, s' = rho^2 - W0."""
        u = arb(u)
        sp = self.rho * self.rho - self.W0
        if sp.upper() <= 0:
            plo = u
        else:
            Y = arb(Y)
            plo = u * Y / (Y * Y + sp.upper()).sqrt()
        if not plo > 0:
            return arb('inf')
        return self.b * (self.cw0.log() - plo.sinh().log())


class GamFam:
    def __init__(self, a, lam, mode='new'):
        self.a = [float(v) for v in a]
        self.lam = [v if isinstance(v, arb) else arb(float(v)) for v in lam]
        self.lamd = [float(v.mid()) for v in self.lam]
        self.mode = mode
        self.tag = 'Gam'

    def lam1_d(self):
        return min(self.lamd)

    def mean_d(self):
        return sum(a / l for a, l in zip(self.a, self.lamd))

    def tilt_mean_d(self, t):
        return sum(a / (l + t) for a, l in zip(self.a, self.lamd))

    def shape_d(self):
        return sum(self.a)

    def lam1(self):
        m = self.lam[0]
        for l in self.lam[1:]:
            m = m.min(l)
        return m

    def prepare(self, x, t0, rho):
        self.t0 = arb(t0)
        self.rho = arb(rho)
        self.lt = [l + self.t0 for l in self.lam]
        self.tau2 = [lt / 2 - self.rho * self.rho for lt in self.lt]
        self.loglt = [lt.log() for lt in self.lt]

    def log_ratio(self, Dh):
        tot = acb(0)
        two = 2 * Dh
        for a, lt in zip(self.a, self.lt):
            tot -= a * (two / lt).log1p()
        return tot

    def log_L_t0(self):
        tot = arb(0)
        for a, l in zip(self.a, self.lam):
            tot -= a * (self.t0 / l).log1p()
        return tot

    def _fmin(self, u2, Ylo, Yhi, tau2):
        # min over Y' in [Ylo, Yhi] of (A - Y')^2 + 4 u^2 Y', A = u^2 + tau2
        # (convex; vertex V = tau2 - u^2, value there 4 u^2 tau2)
        A = u2 + tau2
        V = tau2 - u2
        f = lambda Yp: (A - Yp) ** 2 + 4 * u2 * Yp
        cands = []
        if Yhi is None:
            if V.upper() <= Ylo:
                return f(Ylo)
            cands = [f(Ylo), 4 * u2 * tau2]
        else:
            if V.upper() <= Ylo:
                return f(Ylo)
            if V.lower() >= Yhi:
                return f(Yhi)
            cands = [f(Ylo), f(Yhi), 4 * u2 * tau2]
        m = cands[0]
        for c in cands[1:]:
            m = m.min(c)
        return m

    def cell_log_major(self, u, ya, yb):
        u = arb(u)
        u2 = u * u
        Ylo, Yhi = arb(ya) * arb(ya), arb(yb) * arb(yb)
        tot = arb(0)
        for a, llt, tau2 in zip(self.a, self.loglt, self.tau2):
            m = self._fmin(u2, Ylo, Yhi, tau2)
            if not m > 0:
                return arb('inf')
            tot += a * (llt - (4 * m).log() / 2)
        return tot

    def tail_log_major(self, u, Y):
        u = arb(u)
        u2 = u * u
        Ylo = arb(Y) * arb(Y)
        tot = arb(0)
        for a, llt, tau2 in zip(self.a, self.loglt, self.tau2):
            m = self._fmin(u2, Ylo, None, tau2)
            if not m > 0:
                return arb('inf')
            tot += a * (llt - (4 * m).log() / 2)
        return tot

    def tail_log_combined(self, x, u, Y, h):
        """log of a bound on sum_{k h >= Y} |G(kh)| for the finite family (Lemma maj(b'),
        the note). With A_j = (lam_j + t0)/2, s = y^2 and c_j = k_j s, k_j = 2(A_j - 2u^2)/A_j^2,
        each factor obeys log|...|_j <= (a_j/2) c_j/(1-c_j) while c_j <= 1/2, and its clamped
        bound C_j (over all s >= 1/(2k_j)) beyond. Between consecutive switch points the bound
        B(s) = -2xs + sum of the active pieces is convex in s, so its sup on each interval is at
        an end; past the last switch point B decreases and the Gaussian tail is integrated."""
        u = arb(u); u2 = u * u; xa = arb(x); s0 = float(Y) * float(Y)
        terms = []
        for a, lt, llt, tau2 in zip(self.a, self.lt, self.loglt, self.tau2):
            Aj = lt / 2
            kj = 2 * (Aj - 2 * u2) / (Aj * Aj)
            ku = float(kj.upper())
            sj = 0.5 / ku if ku > 0 else math.inf
            if sj <= s0:
                sj = s0
            if sj < math.inf:
                m = self._fmin(u2, arb(sj), None, tau2)
                if not m > 0:
                    return arb('inf')
                Cj = a * (llt - (4 * m).log() / 2)
            else:
                Cj = None
            terms.append((a, kj, sj, Cj))
        brk = sorted(set([s0] + [t[2] for t in terms if t[2] < math.inf]))
        def B(sv, s_mid):
            sa = arb(sv); tot = -2 * xa * sa
            for a, kj, sj, Cj in terms:
                if s_mid < sj:
                    c = kj * sa
                    tot += a * c / (1 - c) / 2
                else:
                    tot += Cj
            return tot
        parts = []
        for lo, hi in zip(brk, brk[1:]):
            mid = 0.5 * (lo + hi)
            bm = B(lo, mid).max(B(hi, mid))
            ylo, yhi = math.sqrt(lo), math.sqrt(hi)
            nodes = arb((yhi - ylo) / h + 1)
            parts.append(bm + (nodes * (u + arb(yhi))).log())
        sl = brk[-1]; Yl = arb(math.sqrt(sl))
        # beyond sl: terms with k_j > 0 are switched (constants C_j); those with k_j <= 0 give
        # (a/2) c/(1-c) <= 0, nonincreasing in s, so their value at sl bounds them
        Bl = B(sl, sl + 1) + 2 * xa * arb(sl)
        tail = Bl - 2 * xa * Yl * Yl + ((u + Yl) + (u / Yl + 1) / (4 * xa) / h).log()
        parts.append(tail)
        m = parts[0]
        for q in parts[1:]:
            m = m.max(q)
        tot = arb(0)
        for q in parts:
            tot += (q - m).exp()
        return m + tot.log()

# ---------------------------------------------------------------------------
def saddle_t0(fam, x):
    """t0 in (-lam_1, inf) with E_t X = x (double bisection; a choice only)."""
    lam1 = fam.lam1_d()
    lo, hi = -lam1 * (1 - 1e-15), 1.0
    while fam.tilt_mean_d(hi) > x:
        hi = hi * 4 if hi > 0 else 1.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if fam.tilt_mean_d(mid) > x:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def choose_contour(fam, x, c=0.5):
    t0 = saddle_t0(fam, x)
    if isinstance(fam, PGFam):
        W0 = (fam.psi / 2) ** 2 + t0 / 2
        if abs(W0) < 1e-9 * (1 + (fam.psi / 2) ** 2):
            t0 += 4e-9 * (1 + (fam.psi / 2) ** 2)
    # T1 = (t0 + lam_1)/2 is tau_1'^2 + rho^2, frame-invariant; rho^2 = T1 - pi^2/4
    # is the first note's saddle vertex in the untilted frame (s' = 0 for PG)
    T1 = (t0 + fam.lam1_d()) / 2
    rho2 = max(T1 - math.pi ** 2 / 4, c * T1)
    return t0, math.sqrt(rho2)


def _lse(vals):
    """log sum exp of arb values (rigorous: arb exp/log)."""
    tot = arb(0)
    for v in vals:
        tot += v.exp()
    return tot.log()


def _point_log(fam, xa, u, y):
    """log |e^{x(t - t0)} L(t)/L(t0) z| at z = u + i y (a split heuristic only)."""
    ua, ya = arb(u), arb(y)
    Dh = acb(ua * ua - fam.rho * fam.rho - ya * ya, 2 * ua * ya)
    try:
        lr = fam.log_ratio(Dh).real
    except ValueError:
        return float('inf')
    return float((2 * xa * Dh.real + lr + (ua * ua + ya * ya).log() / 2).mid())


def line_major_log(fam, x, u, Y0, ncell, maxdepth=3):
    """log of an upper bound on int_R |G| on the line Re z = u: rigorous cell
    majorants on [0, Y] (cells split adaptively where the majorant exceeds the
    integrand's value at the cell's ends by more than e^2), then the tail."""
    rho = fam.rho
    xa = arb(x)
    pref = 2 * xa * (arb(u) * arb(u) - rho * rho)
    lt0 = fam.tail_log_major(u, Y0)
    Y = Y0
    if lt0.is_finite():
        need = (float((pref + lt0).upper()) + 0.7 * ctx.prec + 25.0) / (2 * x)
        if need > Y * Y:
            Y = math.sqrt(need)
    n0 = int(min(ncell[1], max(ncell[0], 10 * Y * math.sqrt(x))))
    dy = Y / n0
    pts = [_point_log(fam, xa, u, j * dy) for j in range(n0 + 1)]
    est = max(pts) + math.log(Y)
    stack = [(j * dy, (j + 1) * dy, 0, max(pts[j], pts[j + 1])) for j in range(n0)]
    terms = []
    while stack:
        ya, yb, dep, pv = stack.pop()
        yaa = arb(ya)
        m = (-2 * xa * yaa * yaa + fam.cell_log_major(u, ya, yb)
             + (arb(u) * arb(u) + arb(yb) * arb(yb)).log() / 2)
        mu = float(m.upper())
        if dep < maxdepth and mu > est - 40 - math.log(yb - ya) and mu - pv > 2:
            ym = 0.5 * (ya + yb)
            pm = _point_log(fam, xa, u, ym)
            stack.append((ya, ym, dep + 1, max(pv, pm)))
            stack.append((ym, yb, dep + 1, max(pv, pm)))
            continue
        terms.append(m + arb(yb - ya).log())
    lt = fam.tail_log_major(u, Y)
    Ya = arb(Y)
    gt = -2 * xa * Ya * Ya + ((arb(u) / Ya + 1) / (4 * xa)).log()
    terms.append(lt + gt)
    return arb(2).log() + pref + _lse(terms)


def evaluate(fam, x, h, t0, rho, cache, ncell=(200, 3000), full=True, lead_prec=None):
    """One rung: the trapezoid sum at step h and its bounds. Returns dict with
    'lead','S' (arb balls), 'Ed','Et' (arb upper bounds), 'N'."""
    fam.prepare(x, t0, rho)
    xa = arb(x)
    rhoa = fam.rho
    # node sum
    S = rhoa  # k = 0 node: G(0) = rho
    absS = float(rho)
    k = 0
    while True:
        k += 1
        y = k * h
        ya = arb(y)
        Dh = acb(-ya * ya, 2 * rhoa * ya)
        E = 2 * xa * Dh + fam.log_ratio(Dh)
        G = E.exp() * acb(rhoa, ya)
        S += 2 * G.real
        gabs = float(G.abs_upper())
        absS += 2 * gabs
        if k > 4 and gabs < 2.0 ** (-ctx.prec - 8) * absS and y > 1.0 / math.sqrt(x):
            break
        if k > 20000:
            raise RuntimeError('node count')
    N = k
    S = S * h
    out = dict(S=S, N=N, h=h, cond=absS * h / abs(float(S.mid())))
    if lead_prec:
        p0 = ctx.prec
        ctx.prec = lead_prec
        fam.prepare(x, t0, rho)
        lead = arb(2 / arb.pi()).log() + arb(x) * fam.t0 + fam.log_L_t0()
        ctx.prec = p0
        fam.prepare(x, t0, rho)
    else:
        lead = arb(2 / arb.pi()).log() + xa * fam.t0 + fam.log_L_t0()
    out['lead'] = lead
    if not full:
        return out
    # truncation: nodes k > N one by one (rigorous |G| upper bounds), until the
    # monotone tail (Lemma 4.3(c)) past Y2 with the majorant sup_{y >= Y2} is
    # negligible
    yW = 0.5 * (-rho + math.sqrt(rho * rho + 1.0 / x))
    eps_log = math.log(abs(float(S.mid()))) - (ctx.prec + 8) * math.log(2.0)
    tot = arb(0)
    k = N
    while True:
        k += 1
        y = k * h
        if (k - N) % 16 == 1 and y >= yW:
            Ya = arb(y)
            lt = fam.tail_log_major(rho, y)
            tail_nodes = (rhoa + Ya) * (-2 * xa * Ya * Ya).exp() \
                + (-2 * xa * Ya * Ya).exp() * (rhoa / Ya + 1) / (4 * xa) / h
            logT = lt + tail_nodes.log()
            if fam.tag == 'Gam' and not (logT.is_finite() and float(logT.upper()) < eps_log - math.log(h) - 3):
                logT = fam.tail_log_combined(x, rho, y, h)
            if logT.is_finite() and float(logT.upper()) < eps_log - math.log(h) - 3:
                T = logT.exp()
                break
        ya = arb(y)
        Dh = acb(-ya * ya, 2 * rhoa * ya)
        G = (2 * xa * Dh + fam.log_ratio(Dh)).exp() * acb(rhoa, ya)
        tot += G.abs_upper()
        if k - N > 400000:
            raise RuntimeError('truncation nodes')
    Et = 2 * h * (tot + T)
    out['Et'] = Et
    # discretization: two edge lines, best of a few (Lemma 4.2)
    Ydb = max(N * h, 8.0 / math.sqrt(x))
    if LEAN is not None:
        fr_p = cache.setdefault('fr_p', LEAN['fr_p'])
        am = cache.setdefault('am', [-v * rho if v < 0 else v for v in LEAN['am']])
        ncell = LEAN['ncell']
    elif fam.tag == 'PG':
        fr_p = cache.setdefault('fr_p', [0.5, 0.2, 0.05, 0.01])
        am = cache.setdefault('am', [0.5 * rho, rho, 0.1 * rho, 0.25, 1.0, 4.0])
    else:
        fr_p = cache.setdefault('fr_p', [0.5, 0.2, 0.05])
        am = cache.setdefault('am', [0.5 * rho, 0.1 * rho, 0.25, 1.0])
        ncell = (100, 600)
    best_p = None
    for f in fr_p:
        key = ('p', f)
        if key not in cache:
            cache[key] = line_major_log(fam, x, rho * (1 - f), Ydb, ncell, maxdepth=LEAN['maxdepth'] if LEAN else 3)
        a = arb(rho * f)
        val = (cache[key] - (2 * arb.pi() * a / h).expm1().log()).exp()
        best_p = val if best_p is None else best_p.min(val)
    best_m = None
    for a_ in am:
        key = ('m', a_)
        if key not in cache:
            cache[key] = line_major_log(fam, x, rho + a_, Ydb, ncell, maxdepth=LEAN['maxdepth'] if LEAN else 3)
        a = arb(a_)
        val = (cache[key] - (2 * arb.pi() * a / h).expm1().log()).exp()
        best_m = val if best_m is None else best_m.min(val)
    out['Ed'] = best_p + best_m
    return out


# A lean discretization bound (one edge line each side, few cells, no splitting), for the cost
# study: fr_p fractions of rho for a+, am values for a- (negative = that multiple of rho).
LEAN = None

H_SEQ = [m * 2.0 ** -j for j in range(0, 40) for m in (1.0, 0.75)]


def logdensity(fam, x, prec=53, c=0.5, ncell=(200, 3000), maxsteps=40, verbose=False, lead_prec=None):
    """Walk h down until the mathematical bound is below a hundredth of the
    rounding radius (or stops improving). Returns a dict with the enclosure
    of log p(x) and its parts: rad_round (rounding, from arb), Emath/S."""
    old = ctx.prec
    ctx.prec = prec
    try:
        t0, rho = choose_contour(fam, x, c)
        cache = {}
        best = None
        for h in H_SEQ[:maxsteps]:
            try:
                r = evaluate(fam, x, h, t0, rho, cache, ncell, lead_prec=lead_prec)
            except RuntimeError:
                break
            S = r['S']
            Em = (r['Ed'] + r['Et']).upper()
            if not (S.lower() > Em):
                continue
            Sb = S + arb(0, float(Em.upper()) * (1 + 1e-15))
            logp = r['lead'] + Sb.log()
            logp_round = r['lead'] + S.log()
            res = dict(logp=logp, mid=float(logp.mid()), rad=float(logp.rad()),
                       rad_round=float(logp_round.rad()), Emath=float(Em) / float(S.mid()),
                       N=r['N'], h=h, t0=t0, rho=rho, cond=r['cond'],
                       rad_lead=float(r['lead'].rad()))
            if verbose:
                print(res)
            if best is None or res['rad'] < best['rad']:
                best = res
            if res['Emath'] < 0.01 * res['rad_round']:
                break
        return best
    finally:
        ctx.prec = old
