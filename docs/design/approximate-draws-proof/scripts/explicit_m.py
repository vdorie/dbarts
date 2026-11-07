"""Proposition B.1 (proof.tex, Appendix B): the explicit bound max(M_A, M_B, M_C)
on sup f/q, evaluated in 40-digit arithmetic at a few (b, psi, K) as a sanity
check of the formula against the certified lower bounds of kd.out (the bound
must exceed them). Exact rates (r_k = d_k). theta and X are optimized on a grid.
Not a certified computation (no interval arithmetic); a consistency check."""
import mpmath as mp
from pgref import c, d, tail_moments
mp.mp.dps = 40
PI = mp.pi


def lg(a, mu, x):
    """log of the Ga(a, mu) density at x"""
    return a * mp.log(mu) + (a - 1) * mp.log(x) - mu * x - mp.loggamma(a)


def bound(b, psi, K):
    b = mp.mpf(b); psi = mp.mpf(psi); sig = psi ** 2 / 2
    m, v = tail_moments(psi, K)
    al, be = b * m * m / v, m / v
    rt = [c(k) for k in range(1, K + 1)]          # untilted rates r_k - sigma
    bt = be - sig
    A = K * b + al
    J = int(mp.ceil(1 / b))
    logkappa = b * mp.log(mp.cosh(psi / 2)) + sum(b * mp.log(c(k) / d(k, psi)) for k in range(1, K + 1)) \
        + al * mp.log(bt / be)
    logct = sum(b * mp.log(r / bt) for r in rt)
    c1 = c(1)
    logPiJ = sum(b * mp.log(c(k) / c1) for k in range(2, J + 1))
    # E_J = prod_{k>J} (1 - c1/c_k)^-b, with prod_{k>=2} = pi/4
    logEJ = -b * (mp.log(PI / 4) - sum(mp.log(1 - c1 / c(k)) for k in range(2, J + 1)))
    Cb = (b + 2) if b >= 1 else (b + 2) / b
    xs = (b + 1) / (2 * mp.log(2 * Cb))
    mu_t = sum(b / (r - rt[0]) for r in rt[1:]) + al / (bt - rt[0])
    logEq = -sum(b * mp.log(1 - rt[0] / r) for r in rt[1:]) - al * mp.log(1 - rt[0] / bt)
    # M_A
    lMA = (logkappa + (b - 1) * mp.log(2) + mp.log(b) + mp.loggamma(A) - logct - 0.5 * mp.log(2 * PI)
           - A * mp.log(bt) + bt * xs + (A + 0.5) * mp.log(8 * (A + 0.5) / b ** 2) - (A + 0.5))
    if b < 1:
        # E' = prod_{k >= J+2} (1 - c2/c_k)^-b: partial sum plus a tail allowance
        # (-log(1 - t) <= t/(1 - t) and sum_{k>n} 1/c_k <= 1/(2 pi^2 (n - 1/2)))
        c2 = c(2); n = 2000
        s = mp.fsum(mp.log(1 - c2 / c(k)) for k in range(J + 2, n + 1))
        s -= (c2 / (2 * PI ** 2 * (n - 0.5))) / (1 - c2 / c(n + 1))
        logEp = -b * s
    best = None
    for th in ([mp.mpf(1)] if b <= 1 else []) + [mp.mpf(i) / 40 for i in range(1, 40)]:
        for fx in [1.01, 1.1, 1.5, 2, 3, 5, 10, 30, 100]:
            X = mu_t / th * fx
            # M_B over [xs, X]
            if xs < X:
                lr = lambda x: lg(J * b, c1, x) - lg(A, bt, x)
                lMB = logkappa + logPiJ + logEJ - logct + max(lr(xs), lr(X))
            else:
                lMB = -mp.inf
            # M_C
            if b >= 1:
                lB = -b * mp.log(PI / 4)
            else:
                c2 = c(2)
                logPip = sum(b * mp.log(c(k) / c2) for k in range(3, J + 2))
                y = max(X / 2, J * b / (c2 - c1))
                lsup = mp.log(2 * y) + c1 * y + lg(J * b, c2, y)
                lB = mp.log(2 ** (1 - b) * (PI / 4) ** (-b) + mp.exp(logPip + logEp + lsup) / (b * 2 ** b))
            lMC = (logkappa + b * mp.log(c1 / rt[0]) + lB - logEq - (max(b - 1, 0) * mp.log(1 - th) if b > 1 else 0)
                   - mp.log(1 - mu_t / (th * X))) if th < 1 or b <= 1 else mp.inf
            tot = max(lMA, lMB, lMC)
            if best is None or tot < best[0]:
                best = (tot, lMA, lMB, lMC, th, X)
    # limit e^{bD} for reference
    return best, xs, mu_t


if __name__ == "__main__":
    for b, psi, K, lower in [(1, 0, 1, 1.0700), (1, 8, 1, 1.0909), (0.3, 8, 1, 1.2632), (2, 0, 1, 1.0004),
                             (16, 8, 3, 1.0048), (0.2, 150, 1, 35.61)]:
        (tot, a_, b_, c_, th, X), xs, mu = bound(b, psi, K)
        print(f"b={b} psi={psi} K={K}: log10 bound = {float(tot / mp.log(10)):.3g} "
              f"(log10 M_A {float(a_ / mp.log(10)):.3g}, M_B {float(b_ / mp.log(10)) if b_ > -mp.inf else '-'}, "
              f"M_C {float(c_ / mp.log(10)):.3g}; theta {float(th):.3g}, X {float(X):.3g}, x_s {float(xs):.3g}) "
              f"certified lower bound {lower}: {'consistent' if tot >= mp.log(lower) else 'VIOLATED'}")
