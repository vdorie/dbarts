"""The proposal's parameters as the sampler forms them (the rule proof.tex,
Section 12, specifies): c_k = 2 pi^2 (k - 1/2)^2 and sigma = psi^2/2 in double,
r_k = fl(c_k + sigma), r_1 then multiplied by (1 - 2^-50) so that r_1 <= d_1
exactly; m, v the remainder's tail sums (here in 40 digits, then rounded),
alpha = fl(b m^2 / v), beta = fl(m / v). Any positive values give a valid
proposal; only r_1 <= d_1 matters for uniform ergodicity."""
import math, sys
import mpmath as mp
from pgref import d, pg1_mean, pg1_var

def tail_mv(psi, K):
    with mp.workdps(40):
        m = pg1_mean(psi) - mp.fsum(1 / d(k, psi) for k in range(1, K + 1))
        v = pg1_var(psi) - mp.fsum(1 / d(k, psi) ** 2 for k in range(1, K + 1))
        return m, v

def proposal(b, psi, K):
    b = float(b); psi = float(psi)
    sig = 0.5 * psi * psi
    r = [2.0 * math.pi * math.pi * (k - 0.5) ** 2 + sig for k in range(1, K + 1)]
    r[0] = r[0] * (1.0 - 2.0 ** -50)
    m, v = tail_mv(psi, K)
    with mp.workdps(40):
        alpha = float(mp.mpf(b) * m * m / v); beta = float(m / v)
    with mp.workdps(60):
        d1 = mp.pi ** 2 / 2 + mp.mpf(psi) ** 2 / 2
        assert mp.mpf(r[0]) <= d1
    return [b] * K + [alpha], r + [beta]

def untilted(b, psi, K):
    """the first note's q-tilde: rates fl(r_k - sigma), fl(beta - sigma)."""
    a, lam = proposal(b, psi, K)
    sig = 0.5 * float(psi) ** 2
    return a, [l - sig for l in lam]
