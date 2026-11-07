"""Rounding floors (arb enclosure radii) of log f and log q over the reachable
range, old untilted formulation against the new tilted one.
x = mean + z sd of PG(b, psi). Columns: radii of the certified enclosure of
  f_old:  log f0 (untilted, first note's node formula), 53-bit
  f53:    log f (tilted, cancellation-free), 53-bit
  f53L:   same, lead in 106-bit (double-double analog), nodes 53-bit
  f106:   everything 106-bit
  q_old:  log q-tilde (untilted rates fl(r_k - sigma)), 53-bit
  q53:    log q (the sampler's rates), 53-bit
  q53L:   same, lead 106-bit
and LR = 2 (f + q), the log-ratio floor with x and y both at this point.
Usage: floor.py b1,b2 psi1,psi2 z [skipq]"""
import sys, math, time
import mpmath as mp
from pgref import pg1_mean, pg1_var
from kd import K_rule
from arbpg import PGFam, GamFam, logdensity
from proposal import proposal, untilted

def main(bs, psis, z, skipq=False):
    for b in bs:
        for psi in psis:
            mean = float(b * pg1_mean(psi)); sd = float(mp.sqrt(b * pg1_var(psi)))
            x = mean + z * sd
            if x <= 0:
                x = mean / 20
            t = time.time()
            fo = logdensity(PGFam(b, 0.0, 'old'), x)
            fn = logdensity(PGFam(b, psi), x)
            fl = logdensity(PGFam(b, psi), x, lead_prec=106)
            f6 = logdensity(PGFam(b, psi), x, prec=106)
            line = (f"b={b:g} psi={psi:g} z={z:g} x={x:.6g} f_old={fo['rad']:.1e} f53={fn['rad']:.1e} "
                    f"f53L={fl['rad']:.1e} f106={f6['rad']:.1e} Nf={fn['N']}")
            if not skipq:
                K = K_rule(b, psi)
                a, lam = proposal(b, psi, K)
                qo = logdensity(GamFam(*untilted(b, psi, K), mode='old'), x)
                qn = logdensity(GamFam(a, lam), x)
                ql = logdensity(GamFam(a, lam), x, lead_prec=106)
                line += (f" K={K} q_old={qo['rad']:.1e} q53={qn['rad']:.1e} q53L={ql['rad']:.1e} "
                         f"LR_old={2 * (fo['rad'] + qo['rad']):.1e} LR53={2 * (fn['rad'] + qn['rad']):.1e} "
                         f"LR53L={2 * (fl['rad'] + ql['rad']):.1e}")
            print(line + f" t={time.time() - t:.0f}", flush=True)

if __name__ == '__main__':
    main([float(v) for v in sys.argv[1].split(',')], [float(v) for v in sys.argv[2].split(',')],
         float(sys.argv[3]), len(sys.argv) > 4)
