"""Toy 6 (Proposition 6.4): the fixed decision "keep x" on an undecided band, in exact rational arithmetic.
The three-state chain is the first referee's (round 2, chk2/keep.py), rebuilt here on the routines of
exact_toys.py; the two-state chain is the second referee's argument made concrete."""
import os, sys
from fractions import Fraction as Fr
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'exact_toys.py')).read().split('# Toy 1')[0]
exec(src)

def keep_kernel(pi, q, band):
    n = len(pi)
    P = [[Fr(0)] * n for _ in range(n)]
    for x in range(n):
        for y in range(n):
            if y == x:
                P[x][x] += q[y]; continue
            a = min(Fr(1), pi[y] * q[x] / (pi[x] * q[y]))
            # a < 1: the move is taken only when U is certainly below a (U < a - d); U in [a - d, a) keeps x.
            # a = 1: decided once the band is below log(1/r), so the move is always taken.
            acc = max(Fr(0), a - band(x, y, a)) if a < 1 else Fr(1)
            P[x][y] += q[y] * acc; P[x][x] += q[y] * (1 - acc)
    return P

for name, pi, q, band in [
        ('two states, pi=(2/3,1/3), q uniform, band 1/100', [Fr(2, 3), Fr(1, 3)], [Fr(1, 2)] * 2, lambda x, y, a: Fr(1, 100)),
        ('three states, relative band a/100', [Fr(1, 2), Fr(1, 3), Fr(1, 6)], [Fr(1, 3)] * 3, lambda x, y, a: a / 100),
        ('three states, band (1+x+2y)/100', [Fr(1, 2), Fr(1, 3), Fr(1, 6)], [Fr(1, 3)] * 3, lambda x, y, a: Fr(1 + x + 2 * y, 100))]:
    P = keep_kernel(pi, q, band)
    st = stationary(P)
    print(f'{name}: max|piP-pi| = {residual(pi, P)} = {float(residual(pi, P)):.4g}; stationary {st}; TV = {float(tv(st, pi)):.4g}')
