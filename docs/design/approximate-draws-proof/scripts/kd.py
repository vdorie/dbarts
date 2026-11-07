"""The K rule (b D(psi, K) <= 0.01) and the proposal weight w = f/q.

D from its closed form (proof.tex, Lemma 3.2):
  D = -log(d1 / (2 pi cosh(psi/2))) + sum_{k=2}^K log(1 - d1/d_k) + (m^2/v) log(1 - d1/beta)
with the remainder's moments m, v as tail sums. Checked against direct
summation. Then, for a few (b, psi, K), the grid maximum of log w with the
certified density bounds: a certified LOWER bound on sup w (the sup over the
grid can only be smaller than the sup over (0, inf)), against exp(b D).
Not a certified upper bound on sup w (Section 9, Not proven)."""
import mpmath as mp, math, sys
from pgref import c, d, pg1_mean, pg1_var, tail_moments
from certified import F0, QT, auto, left_series
mp.mp.dps = 40

def D_closed(psi, K):
    psi = mp.mpf(psi)
    d1 = d(1, psi)
    m, v = tail_moments(psi, K)
    beta = m / v
    s = -mp.log(d1 / (2 * mp.pi * mp.cosh(psi / 2)))
    s += mp.fsum(mp.log(1 - d1 / d(k, psi)) for k in range(2, K + 1))
    return s + (m * m / v) * mp.log(1 - d1 / beta)

def D_direct(psi, K):
    """direct summation (mpmath nsum); reliable for small K only (used K <= 10)."""
    psi = mp.mpf(psi)
    d1 = d(1, psi)
    m, v = tail_moments(psi, K)
    s = mp.nsum(lambda k: -mp.log(1 - d1 / d(int(k), psi)), [K + 1, mp.inf])
    return s + (m * m / v) * mp.log(1 - d1 / (m / v))

def K_rule(b, psi, tol=0.01, Kmax=1000):
    for K in range(1, Kmax + 1):
        if b * D_closed(psi, K) <= tol:
            return K
    return None

def logw_grid(b, psi, K, zs):
    """certified log w = log f - log q at x = mean + z sd of PG(b, psi)."""
    m, v = tail_moments(psi, K)
    al = float(b * m * m / v); be = m / v
    sig = mp.mpf(psi) ** 2 / 2
    bt = float(be - sig)
    f0, qt = F0(b), QT.proposal(b, K, al, bt)
    # constants: log cosh^b(psi/2) - sum b log(d_k/c_k) - alpha log(beta/beta_t)
    const = float(b * mp.log(mp.cosh(mp.mpf(psi) / 2))
                  - mp.fsum(b * mp.log(d(k, psi) / c(k)) for k in range(1, K + 1))
                  - al * mp.log(be / mp.mpf(bt)))
    mean = float(b * pg1_mean(psi)); sd = float(mp.sqrt(b * pg1_var(psi)))
    out = []
    for z in zs:
        x = mean + z * sd
        if x <= 0:
            continue
        lf = auto(f0, x); lq = auto(qt, x)
        ls = left_series(b, x)
        lfv, lfe = (ls[0], ls[1]) if ls and ls[1] < lf['err'] else (lf['logp'], lf['err'])
        out.append((z, x, lfv - lq['logp'] + const, lfe + lq['err']))
    return out

if __name__ == '__main__':
    print('# D closed form vs direct summation')
    for psi in [0, 8, 64, 150]:
        for K in [1, 3, 10]:
            a, bb = D_closed(psi, K), D_direct(psi, K)
            print(f'psi={psi} K={K} D={mp.nstr(a, 12)} |closed-direct|={mp.nstr(abs(a - bb), 3)}')
    print('# K rule: smallest K >= 1 with b D(psi, K) <= 0.01')
    for psi in [0, 8, 50, 64, 100, 150]:
        row = []
        for b in [0.2, 1, 16, 100, 1e3, 1e4, 1e5]:
            row.append(f'{K_rule(b, psi)}')
        print(f'psi={psi}: ' + ' '.join(f'b={b:g}:K={k}' for b, k in zip([0.2, 1, 16, 100, 1e3, 1e4, 1e5], row)), flush=True)
    print('# b at which b D(psi, 48) = 0.01 (the K = 48 cap fails above it)')
    for psi in [64, 100, 150]:
        print(f'psi={psi}: b* = {float(0.01 / D_closed(psi, 48)):.4g}')
    print('# K rule at b = 1e5, psi = 150:', K_rule(1e5, 150))
    print('# grid max of w (certified lower bound on sup w) vs exp(b D)')
    zs = [-2, -1.5, -1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 1, 1.5, 2, 3, 4, 6, 8, 12, 16]
    for b, psi, K in [(0.3, 8, 1), (0.3, 8, 2), (1, 0, 1), (1, 8, 1), (2, 0, 1), (16, 8, 3), (0.2, 150, 1), (4, 100, 2)]:
        g = logw_grid(b, psi, K, zs)
        zmax, xmax, lw, e = max(g, key=lambda t: t[2])
        print(f'b={b} psi={psi} K={K}: exp(bD)={float(mp.exp(b * D_closed(psi, K))):.6f} '
              f'max_grid w={math.exp(lw):.6f} (+-{e:.1e} in log) at z={zmax}, '
              f'w at z=16: {math.exp(g[-1][2]):.6f}', flush=True)
