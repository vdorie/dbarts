"""Radius of the enclosure of log f (and log q with 'q') at 106 bits away from the mean.
Usage: hi106.py f|q b1,b2 psi1,psi2 z1,z2"""
import sys, time
import mpmath as mp
from pgref import pg1_mean, pg1_var
from kd import K_rule
from arbpg import PGFam, GamFam, logdensity
from proposal import proposal
which = sys.argv[1]
for b in [float(v) for v in sys.argv[2].split(',')]:
    for psi in [float(v) for v in sys.argv[3].split(',')]:
        mean = float(b * pg1_mean(psi)); sd = float(mp.sqrt(b * pg1_var(psi)))
        if which == 'q':
            K = K_rule(b, psi); a, lam = proposal(b, psi, K); fam = GamFam(a, lam)
        else:
            K = ''; fam = PGFam(b, psi)
        for z in [float(v) for v in sys.argv[4].split(',')]:
            x = mean + z * sd
            if x <= 0: x = mean / 20
            t = time.time(); r = logdensity(fam, x, prec=106)
            print(f"{which}106 b={b:g} psi={psi:g} K={K} z={z:g} x={x:.6g} rad={r['rad']:.2e} t={time.time()-t:.1f}",
                  flush=True)
