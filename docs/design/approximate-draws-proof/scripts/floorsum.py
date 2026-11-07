"""Summaries of the 53-bit floor from pins_f.out and pins_q.out: max radius by z band, and the
averaged undecided probability per step, 2 * (E rf(x) + E rq(x) + E rf(y) + E rq(y)), with x and y
taken as N(0,1) in z (both near-normal at the shapes where the radius matters) and the radius
interpolated linearly in z between pins (held constant beyond the last pin).
Usage: floorsum.py pins_f.out pins_q.out"""
import sys, math
from collections import defaultdict

def load(path):
    d = defaultdict(dict)
    for line in open(path):
        if 'rad=' not in line or line.startswith('SUMMARY'):
            continue
        kv = dict(t.split('=', 1) for t in line.split() if '=' in t)
        d[(float(kv['b']), float(kv['psi']))][float(kv['z'])] = float(kv['rad'])
    return d

def interp(tab, z):
    zs = sorted(tab)
    if z <= zs[0]: return tab[zs[0]]
    if z >= zs[-1]: return tab[zs[-1]]
    for a, b in zip(zs, zs[1:]):
        if a <= z <= b:
            return tab[a] + (tab[b] - tab[a]) * (z - a) / (b - a)

def mean_rad(tab):
    n = 4000; tot = 0.0
    for j in range(n):
        z = -8 + 16 * (j + 0.5) / n
        tot += interp(tab, z) * math.exp(-z * z / 2) / math.sqrt(2 * math.pi) * 16 / n
    return tot

F, Q = load(sys.argv[1]), load(sys.argv[2])
for name, D in (('f', F), ('q', Q)):
    for lo, hi in ((0, 16), (16.1, 1000), (1000.1, 1e5)):
        sel = [(k, t) for k, t in D.items() if lo <= k[0] <= hi]
        if not sel: continue
        m0 = max((t.get(0.0, t.get(0.5, 0)) for k, t in sel))
        m3 = max(max(v for z, v in t.items() if abs(z) <= 3) for k, t in sel)
        ma = max(max(t.values()) for k, t in sel)
        print(f'{name} b in [{lo:g},{hi:g}]: max rad z~0 {m0:.2e}, |z|<=3 {m3:.2e}, all {ma:.2e}')
worst = 0
for k in sorted(set(F) & set(Q)):
    p = 4 * (mean_rad(F[k]) + mean_rad(Q[k]))
    worst = max(worst, p)
    print(f'b={k[0]:g} psi={k[1]:g}: averaged undecided prob <= {p:.1e}')
print(f'worst averaged over common shapes: {worst:.1e}')
