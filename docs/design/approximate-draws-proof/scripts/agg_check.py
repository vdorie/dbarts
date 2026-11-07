"""Numerical check of Lemma 5.9 (aggregated majorant) for f: closed forms T1, T2 against direct sums, and the
bound against log|Lt(z)/Lt(rho)| on lines u, at y with kappa Y < 1. mpmath, 30 digits."""
import mpmath as mp
mp.mp.dps = 30
def T1(W):
    r = mp.sqrt(mp.mpc(W)); return mp.re(mp.tanh(r) / (2 * r)) if W != 0 else mp.mpf(1) / 2
def T2(W):
    r = mp.sqrt(mp.mpc(W)); return mp.re((mp.tanh(r) - r * mp.sech(r) ** 2) / (4 * r ** 3)) if W != 0 else mp.mpf(1) / 6
c = lambda j: mp.pi ** 2 * (j - mp.mpf(1) / 2) ** 2
for W in [-2.0, 0.0, 0.3, 5.0, 400.0]:
    s1 = mp.nsum(lambda j: 1 / (W + c(j)), [1, mp.inf]); s2 = mp.nsum(lambda j: 1 / (W + c(j)) ** 2, [1, mp.inf])
    print(f'W={W}: T1 err {float(abs(T1(W) - s1)):.1e}, T2 err {float(abs(T2(W) - s2)):.1e}')
worst = -mp.inf
for b, sp, rho, u in [(16, 0.5, 1.2, 0.6), (16, 0.5, 1.2, 1.8), (1e5, 2.0, 1.6, 0.8), (0.2, 0.0, 0.9, 2.5), (64, 2.4, 0.4, 0.2)]:
    b = mp.mpf(b); Wr = rho ** 2 - sp; Wu = u ** 2 - sp
    A = lambda j: Wu + c(j)
    J0 = [j for j in range(1, 200) if A(j) < 4 * u ** 2]
    j1 = (max(J0) + 1) if J0 else 1
    kap = 2 * (A(j1) - 2 * u ** 2) / A(j1) ** 2
    Stot = 2 * b * (T1(Wu) - 2 * u ** 2 * T2(Wu))
    S = Stot - sum(b * 2 * (A(j) - 2 * u ** 2) / A(j) ** 2 for j in J0)
    lc = lambda W: mp.log(mp.cosh(mp.sqrt(W))) if W >= 0 else mp.log(mp.cos(mp.sqrt(-W)))
    const = b * (lc(Wr) - lc(Wu))
    for frac in [0.1, 0.5, 0.9, 0.99]:
        Y = frac / kap; y = mp.sqrt(Y); z = mp.mpc(u, y)
        w = mp.sqrt(z * z - sp); wr = mp.sqrt(mp.mpc(Wr))
        lhs = b * (mp.log(abs(mp.cosh(wr))) - mp.log(abs(mp.cosh(w))))
        j0t = sum(-(b / 2) * mp.log(abs(z * z + c(j) - sp) ** 2 / A(j) ** 2) for j in J0)
        rhs = const + j0t + S / 2 * Y / (1 - kap * Y)
        worst = max(worst, lhs - rhs)
        print(f'b={float(b):g} s\'={sp} rho={rho} u={u} |J0|={len(J0)} kY={frac}: lhs {float(lhs):.6g} rhs {float(rhs):.6g} ok {lhs <= rhs}')
print('max lhs - rhs', float(worst))
