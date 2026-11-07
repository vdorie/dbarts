"""Radius of log q at b = 1e5, |psi| = 150 (K = 560 by the rule), x = mean + z sd, lead at 106 bits.
No reference (out of reach). Usage: q560.py z1,z2"""
import sys, time
import mpmath as mp
from pgref import pg1_mean, pg1_var
from kd import K_rule
from arbpg import GamFam, logdensity
from proposal import proposal
b, psi = 1e5, 150.0
K = K_rule(b, psi); a, lam = proposal(b, psi, K); fam = GamFam(a, lam)
mean = float(b * pg1_mean(psi)); sd = float(mp.sqrt(b * pg1_var(psi)))
for z in [float(v) for v in sys.argv[1].split(',')]:
    t = time.time(); r = logdensity(fam, mean + z * sd, lead_prec=106)
    print(f"q560 K={K} z={z:g} rad={r['rad']:.2e} lead={r['rad_lead']:.2e} N={r['N']} h={r['h']:.3g} t={time.time()-t:.0f}", flush=True)
