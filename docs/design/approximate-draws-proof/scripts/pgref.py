"""High-precision reference densities for PG(b, psi), its truncated series
remainder R_K, the moment-matched proposal q, and helpers.

Not a proof: an independent numerical reference (mpmath, 50 digits) used to
sanity-check the proven statements and to exhibit counterexamples.

Notation (proof.tex): c_k = 2 pi^2 (k - 1/2)^2, d_k = c_k + psi^2/2,
PG(b, psi) = sum_k g_k / d_k, g_k ~ Ga(b, 1).
"""
import mpmath as mp

mp.mp.dps = 50
PI = mp.pi


def c(k):
    return 2 * PI**2 * (k - mp.mpf(1) / 2) ** 2


def d(k, psi):
    return c(k) + mp.mpf(psi) ** 2 / 2


def pg1_mean(psi):
    psi = mp.mpf(psi)
    if psi == 0:
        return mp.mpf(1) / 4
    return mp.tanh(psi / 2) / (2 * psi)


def pg1_var(psi):
    psi = mp.mpf(psi)
    if abs(psi) < mp.mpf('1e-8'):
        return mp.mpf(1) / 24
    return (mp.sinh(psi) - psi) / (4 * psi**3 * mp.cosh(psi / 2) ** 2)


def tail_moments(psi, K):
    """sum_{k>K} 1/d_k and sum_{k>K} 1/d_k^2 (per unit shape)."""
    m = pg1_mean(psi) - mp.fsum(1 / d(k, psi) for k in range(1, K + 1))
    v = pg1_var(psi) - mp.fsum(1 / d(k, psi) ** 2 for k in range(1, K + 1))
    return m, v


def proposal_params(b, psi, K):
    """alpha, beta of the moment-matched gamma remainder."""
    m, v = tail_moments(psi, K)
    b = mp.mpf(b)
    return b * m * m / v, m / v


def D(psi, K):
    """Right-tail log ratio per unit b: lim f/q = exp(b D)."""
    m, v = tail_moments(psi, K)
    d1 = d(1, psi)
    beta = m / v
    alpha_over_b = m * m / v
    s = mp.nsum(lambda k: -mp.log(1 - d1 / d(int(k), psi)), [K + 1, mp.inf])
    return s + alpha_over_b * mp.log(1 - d1 / beta)


def log_cosh(z):
    """Branch of log cosh z continuous on Re z > 0, real on the real axis."""
    return z + mp.log(1 + mp.exp(-2 * z)) - mp.log(2)


def logL(t, b, psi, terms=(), tail_from=None, gam=None):
    """log Laplace transform E exp(-t X) of
    X = sum_{(a, lam) in terms} Ga(a, lam) + [sum_{k>tail_from} g_k/d_k] + Ga(gam).
    tail_from = 0 gives the full PG(b, psi)."""
    out = mp.mpf(0)
    for a, lam in terms:
        out -= a * mp.log(1 + t / lam)
    if tail_from is not None:
        s = t + mp.mpf(psi) ** 2 / 2
        z = mp.sqrt(s / 2)
        out += -b * log_cosh(z) + b * log_cosh(mp.mpf(psi) / 2)
        for k in range(1, tail_from + 1):
            out += b * mp.log(1 + t / d(k, psi))
    if gam is not None:
        a, lam = gam
        out -= a * mp.log(1 + t / lam)
    return out


def density_zline(x, b, psi, rho=None, **kw):
    """(4/pi) Re int_0^inf exp(t x) L(t) z dy on z = rho + i y,
    t = 2 z^2 - psi^2/2. Valid when every singularity of L, written in z,
    lies on the imaginary axis (true for all families here)."""
    x = mp.mpf(x)
    psi = mp.mpf(psi)
    if rho is None:
        rho = 1 / (2 * mp.sqrt(x)) + mp.mpf('0.05')

    def g(y):
        z = rho + 1j * y
        t = 2 * z * z - psi**2 / 2
        return mp.re(mp.exp(t * x + logL(t, b, psi, **kw)) * z)

    w = 1 / mp.sqrt(x)
    pts = [0] + [w * j for j in (0.5, 1, 2, 4, 8, 16, 32)] + [mp.inf]
    return 4 / PI * mp.quad(g, pts)


def pg0_left_series(x, b, nterms=None):
    """Left alternating series for PG(b, 0) at x (PG scale: WPS eq. (11) on
    the J* scale with x_J = 4 x)."""
    x = mp.mpf(x)
    b = mp.mpf(b)
    tot = mp.mpf(0)
    n = 0
    lead = None
    while True:
        a = (mp.exp(mp.loggamma(n + b) - mp.loggamma(n + 1)) * (2 * n + b)
             / mp.sqrt(2 * PI * x**3) * mp.exp(-(2 * n + b) ** 2 / (8 * x)))
        tot += (-1) ** n * a
        if lead is None:
            lead = a
        n += 1
        if nterms is not None and n >= nterms:
            break
        if nterms is None and a < lead * mp.mpf(10) ** (-45) and n > 5:
            break
    return 2**(b - 1) / mp.gamma(b) * tot


def pg_density(x, b, psi):
    """PG(b, psi) density by exact tilt of PG(b, 0)."""
    psi = mp.mpf(psi)
    return (mp.cosh(psi / 2) ** b * mp.exp(-psi**2 * mp.mpf(x) / 2)
            * density_zline(x, b, 0, tail_from=0))


def q_density(x, b, psi, K):
    al, be = proposal_params(b, psi, K)
    terms = [(mp.mpf(b), d(k, psi)) for k in range(1, K + 1)]
    return density_zline(x, b, psi, terms=terms, gam=(al, be))


def f_density(x, b, psi):
    return density_zline(x, b, psi, tail_from=0)


def density_shifted(x, logLt, rho, s, scale=None):
    """p(x) = e^{-2 s x} (2/pi) int_R Re[e^{2x z^2} Lt(z) z] dy, z = rho + i y,
    Lt(z) = L(2 z^2 - 2 s), given log Lt as a function of z (Lemma 5.1)."""
    x, rho, s = mp.mpf(x), mp.mpf(rho), mp.mpf(s)
    l0 = logLt(rho)

    def g(y):
        z = rho + 1j * y
        return mp.re(mp.exp(2 * x * (z * z - rho * rho) + logLt(z) - l0) * z)

    w = (scale if scale is not None else 1 / mp.sqrt(x))
    pts = [0] + [w * j for j in (0.25, 0.5, 1, 2, 3, 4, 6, 8, 12, 16, 32, 64)] + [mp.inf]
    I = 2 * mp.quad(g, pts)
    return mp.re(mp.log(2 / PI) - 2 * s * x + 2 * x * rho * rho + l0 + mp.log(I))


def logf0_shifted(x, b, rho, s):
    """log PG(b, 0) density at x on the contour (rho, s)."""
    b = mp.mpf(b)
    s = mp.mpf(s)
    return density_shifted(x, lambda z: -b * log_cosh(mp.sqrt(z * z - s)), rho, s,
                           scale=1 / mp.sqrt(mp.mpf(x) + b))


def logqt_shifted(x, terms, rho, s):
    """log density of sum Ga(a, lam) at x, terms = [(a, lam)], contour (rho, s)."""
    s = mp.mpf(s)
    tot = mp.fsum(mp.mpf(a) for a, lam in terms)
    return density_shifted(
        x, lambda z: -mp.fsum(mp.mpf(a) * mp.log(1 + (2 * z * z - 2 * s) / mp.mpf(lam))
                              for a, lam in terms), rho, s, scale=1 / mp.sqrt(mp.mpf(x) + tot))
