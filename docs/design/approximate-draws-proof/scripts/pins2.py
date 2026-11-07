"""Pins for the tilted evaluator (arbpg.py, 53-bit): the enclosure of log f
(f = PG(b, psi)) or log q (the proposal from the sampler's own rates, K by the
rule) must contain an independent mpmath reference. x = mean + z sd of
PG(b, psi) (x = mean/20 if that is not positive). References: pgref's
parabolic contour (mpmath quad) on two contours (c = 0.3 and 0.6), and for f
with b <= 16 also the left alternating series. 'act' is |mid - ref| (first
reference), 'rad' the enclosure radius; ok means act <= rad.
Usage: pins2.py f|q b1,b2 psi1,psi2 z1,z2"""
import sys, math, time
import mpmath as mp
from pgref import pg1_mean, pg1_var
from kd import K_rule
from arbpg import PGFam, GamFam, logdensity, choose_contour
from proposal import proposal
import refs, os
LEAD = int(os.environ['LEAD']) if os.environ.get('LEAD') else None  # lead at this precision (106)

def main(which, bs, psis, zs):
    n = nfail = 0; worst = 0.0; maxrad = 0.0
    for b in bs:
        for psi in psis:
            mean = float(b * pg1_mean(psi)); sd = float(mp.sqrt(b * pg1_var(psi)))
            if which == 'q':
                K = K_rule(b, psi)
                a, lam = proposal(b, psi, K)
                fam = GamFam(a, lam); tag = f'b={b:g} psi={psi:g} K={K}'
            else:
                fam = PGFam(b, psi); tag = f'b={b:g} psi={psi:g}'
            for z in zs:
                x = mean + z * sd
                if x <= 0:
                    x = mean / 20
                t = time.time()
                r = logdensity(fam, x, lead_prec=LEAD)
                t1 = time.time() - t
                if r is None:
                    print(f'{which} {tag} z={z:g} x={x:.6g} NO ENCLOSURE', flush=True); nfail += 1; n += 1
                    continue
                refl = []
                for c in (0.3, 0.6):
                    t0, rho = choose_contour(fam, x, c)
                    if which == 'q':
                        refl.append(refs.logq_contour(a, lam, x, rho, rho * rho - t0 / 2, dps=34))
                    else:
                        W0 = (psi / 2) ** 2 + t0 / 2
                        refl.append(refs.logf_contour(b, psi, x, rho, rho * rho - W0, dps=50))
                if which == 'f' and b <= 16:
                    refl.append(refs.logf_left(b, psi, x, r['mid']))
                act = float(abs(refl[0] - mp.mpf(r['mid'])))
                dis = max(float(abs(v - refl[0])) for v in refl[1:])
                ok = act <= r['rad']
                n += 1; nfail += (not ok); worst = max(worst, act / r['rad']); maxrad = max(maxrad, r['rad'])
                print(f"{which} {tag} z={z:g} x={x:.6g} N={r['N']} h={r['h']:.3g} rad={r['rad']:.2e} "
                      f"round={r['rad_round']:.2e} act={act:.2e} ok={ok} refdis={dis:.0e} t={t1:.1f}", flush=True)
    print(f'SUMMARY {which}: pins={n} failures={nfail} max act/rad={worst:.3f} max rad={maxrad:.2e}', flush=True)

if __name__ == '__main__':
    main(sys.argv[1], [float(v) for v in sys.argv[2].split(',')], [float(v) for v in sys.argv[3].split(',')],
         [float(v) for v in sys.argv[4].split(',')])
