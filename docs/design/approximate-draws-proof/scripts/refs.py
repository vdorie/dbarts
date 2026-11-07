"""Independent high-precision references (mpmath), returning mpf.
f = PG(b, psi): log f = log f0 + b log cosh(psi/2) - psi^2 x / 2, with log f0 by
  (i) the left alternating series in adaptive precision (formula of
      the second referee's check, here returning mpf), and
  (ii) the parabolic contour of pgref.py (mpmath quad).
q = sum_j Ga(a_j, lam_j) (the sampler's own rates): pgref.logqt_shifted."""
import sys, math
import mpmath as mp
import pgref

def logf0_left(b, x, est):
    b = mp.mpf(b); x = mp.mpf(x)
    nmax = max(0, int((math.sqrt(8 * float(x) * max(1.0, float(b))) - float(b)) / 2))
    lt = lambda n: (mp.loggamma(n + b) - mp.loggamma(n + 1) - mp.loggamma(b) + mp.log(2 * n + b)
                    - (2 * n + b) ** 2 / (8 * x))
    with mp.workdps(30):
        lmax = max(lt(n) for n in set([0, nmax // 2, nmax, 2 * nmax + 1]))
        lpref = (b - 1) * mp.log(2) - mp.log(2 * mp.pi * x ** 3) / 2
    need = int((float(lmax + lpref) - est) / 2.30 + 45 + math.log10(abs(est) + 10))
    with mp.workdps(max(45, need)):
        b = mp.mpf(b); x = mp.mpf(x)
        s = mp.mpf(0); n = 0
        while True:
            t = mp.exp(lt(n))
            s += t if n % 2 == 0 else -t
            if n > nmax and t < mp.mpf(10) ** (-max(45, need) + 10) * abs(s):
                break
            n += 1
        return +(mp.log(s) + (b - 1) * mp.log(2) - mp.log(2 * mp.pi * x ** 3) / 2)

def tilt(b, psi, x):
    psi = mp.mpf(psi)
    return mp.mpf(b) * mp.log(mp.cosh(psi / 2)) - psi ** 2 * mp.mpf(x) / 2

def logf_left(b, psi, x, est_f):
    with mp.workdps(60):
        T = tilt(b, psi, x)
        return +(logf0_left(b, x, est_f - float(T)) + T)

def logf_contour(b, psi, x, rho, s, dps=50):
    with mp.workdps(dps):
        T = tilt(b, psi, x)
        return +(pgref.logf0_shifted(mp.mpf(x), mp.mpf(b), rho, s) + T)

def logq_contour(a, lam, x, rho, s, dps=40):
    with mp.workdps(dps):
        return +pgref.logqt_shifted(mp.mpf(x), [(mp.mpf(ai), mp.mpf(li)) for ai, li in zip(a, lam)], rho, s)

def logq_bromwich(a, lam, x, dps=30):
    """log q(x) by the Bromwich integral on the vertical line through the saddle t0 (independent of
    the parabolic contour of Lemma 4.1): q(x) = (1/pi) int_0^inf Re[e^{x(t0+i eta)} L(t0+i eta)] d eta,
    L(t) = prod_j (1 + t/lam_j)^{-a_j}; mpmath quadrature split at multiples of the tilted sd."""
    with mp.workdps(dps):
        al = [mp.mpf(v) for v in a]; ll = [mp.mpf(v) for v in lam]; x = mp.mpf(x)
        g = lambda t: mp.fsum(av / (l + t) for av, l in zip(al, ll)) - x
        lo = -min(ll) * (1 - mp.mpf(10) ** -12); hi = mp.mpf(1)
        while g(hi) > 0:
            hi *= 4
        for _ in range(int(3.4 * dps) + 20):
            mid = (lo + hi) / 2
            if g(mid) > 0: lo = mid
            else: hi = mid
        t0 = (lo + hi) / 2
        logL = lambda t: -mp.fsum(av * mp.log1p(t / l) for av, l in zip(al, ll))
        L0 = logL(t0)
        sd = 1 / mp.sqrt(mp.fsum(av / (l + t0) ** 2 for av, l in zip(al, ll)))
        f = lambda eta: mp.re(mp.exp(1j * eta * x + logL(t0 + 1j * eta) - L0))
        pts = [0] + [sd * k for k in (0.25, 0.5, 1, 2, 3, 4, 6, 8, 12, 16, 32, 64)] + [mp.inf]
        I = mp.quad(f, pts)
        return +(t0 * x + L0 + mp.log(I / mp.pi))


def logq_mosch(a, lam, x, dps=40, maxterms=4000):
    """log q(x) by Moschopoulos's series (Ann. Inst. Statist. Math. 37, 1985): the density of
    sum_j Ga(a_j, lam_j) is C sum_k delta_k g(rho+k, lam_max)(x), with C = prod (lam_j/lam_max)^a_j,
    rho = sum a_j, gamma_k = sum_j a_j (1 - lam_j/lam_max)^k / k, delta_0 = 1,
    delta_{k+1} = (1/(k+1)) sum_{i=1}^{k+1} i gamma_i delta_{k+1-i}. All terms are positive, so no
    cancellation; independent of any contour. Returns None when lam_max * x is too large to sum."""
    with mp.workdps(dps):
        al = [mp.mpf(v) for v in a]; ll = [mp.mpf(v) for v in lam]; x = mp.mpf(x)
        lm = max(ll)
        if lm * x > maxterms / 2.5:
            return None
        rho = mp.fsum(al)
        logC = mp.fsum(av * mp.log(l / lm) for av, l in zip(al, ll))
        om = [1 - l / lm for l in ll]
        gam = [None]; dl = [mp.mpf(1)]
        lx = mp.log(x); llm = mp.log(lm)
        def lterm(k):
            return (rho + k) * llm + (rho + k - 1) * lx - lm * x - mp.loggamma(rho + k)
        total = dl[0] * mp.exp(lterm(0)); small = 0
        for k in range(1, maxterms):
            gam.append(mp.fsum(av * w ** k for av, w in zip(al, om)) / k)
            dl.append(mp.fsum(i * gam[i] * dl[k - i] for i in range(1, k + 1)) / k)
            t = dl[k] * mp.exp(lterm(k))
            total += t
            if k > lm * x + 10 and t < total * mp.mpf(10) ** (-dps + 2):
                small += 1
                if small > 20:
                    return +(logC + mp.log(total))
            else:
                small = 0
        return None
