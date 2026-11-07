"""Certified lower bounds on M* = sup f/q from grid points: log w(x) = log f(x) - log q(x),
both enclosed by arbpg (53-bit), q from the sampler's rates (proposal.py).
Usage: wgrid.py b psi K z_lo z_hi n"""
import sys, math
import mpmath as mp
from pgref import pg1_mean, pg1_var
from kd import D_closed
from arbpg import PGFam, GamFam, logdensity
from proposal import proposal
b, psi, K = float(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3])
zl, zh, n = float(sys.argv[4]), float(sys.argv[5]), int(sys.argv[6])
mean = float(b * pg1_mean(psi)); sd = float(mp.sqrt(b * pg1_var(psi)))
a, lam = proposal(b, psi, K)
F, Q = PGFam(b, psi), GamFam(a, lam)
best = None
for j in range(n + 1):
    z = zl + (zh - zl) * j / n
    x = mean + z * sd
    if x <= 0: continue
    lw = logdensity(F, x)['logp'] - logdensity(Q, x)['logp']
    lo = float(lw.lower())
    if best is None or lo > best[0]: best = (lo, z, float(lw.rad()))
print(f"b={b:g} psi={psi:g} K={K} exp(bD)={float(mp.exp(b * D_closed(psi, K))):.5f} M*>={math.exp(best[0]):.5f} at z={best[1]:.3f} (rad {best[2]:.1e})")
