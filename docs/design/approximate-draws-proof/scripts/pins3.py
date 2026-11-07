"""pins3: as pins2, but the two reference contours differ from the enclosure's
(rho^2 = max(T1 - pi^2/8, 0.75 T1) and 0.3 T1, against the enclosure's max(T1 - pi^2/4, T1/2)).
Pins for the tilted evaluator (arbpg.py, 53-bit): the enclosure of log f
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
                t0, _ = choose_contour(fam, x, 0.5)
                T1 = (t0 + fam.lam1_d()) / 2
                for rho2 in (max(T1 - math.pi ** 2 / 8, 0.75 * T1), 0.3 * T1):
                    rho = math.sqrt(rho2)
                    if which == 'q':
                        refl.append(refs.logq_contour(a, lam, x, rho, rho * rho - t0 / 2, dps=34))
                    else:
                        W0 = (psi / 2) ** 2 + t0 / 2
                        refl.append(refs.logf_contour(b, psi, x, rho, rho * rho - W0, dps=50))
                if which == 'q' and os.environ.get('BROM'):
                    rm = refs.logq_mosch(a, lam, x); r3 = 'mosch' if rm is not None else 'brom'
                    refl.append(rm if rm is not None else refs.logq_bromwich(a, lam, x, dps=30))
                if which == 'f' and b <= 16:
                    refl.append(refs.logf_left(b, psi, x, r['mid']))
                act = max(float(abs(v - mp.mpf(r['mid']))) for v in refl)
                dis = max(float(abs(v - refl[0])) for v in refl[1:])
                ok = act <= r['rad']
                n += 1; nfail += (not ok); worst = max(worst, act / r['rad']); maxrad = max(maxrad, r['rad'])
                print(f"{which} {tag} z={z:g} x={x:.6g} N={r['N']} h={r['h']:.3g} rad={r['rad']:.2e} "
                      f"round={r['rad_round']:.2e} act={act:.2e} ok={ok} refdis={dis:.0e} t={t1:.1f}" + (f" ref3={r3}" if which == 'q' and os.environ.get('BROM') else ''), flush=True)
    print(f'SUMMARY {which}: pins={n} failures={nfail} max act/rad={worst:.3f} max rad={maxrad:.2e}', flush=True)

if __name__ == '__main__':
    main(sys.argv[1], [float(v) for v in sys.argv[2].split(',')], [float(v) for v in sys.argv[3].split(',')],
         [float(v) for v in sys.argv[4].split(',')])
