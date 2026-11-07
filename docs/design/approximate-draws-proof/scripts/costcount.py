"""Operation counts for one certified evaluation at the final step h that logdensity selects:
main nodes N, truncation nodes evaluated one by one, and line-majorant cells for the
discretization bound (each a per-term bound: K+1 terms for q, O(1) for f). Lead at 106 bits.
Usage: costcount.py b1,b2 psi1,psi2 z"""
import sys, time
import mpmath as mp
from flint import ctx
from pgref import pg1_mean, pg1_var
from kd import K_rule
import arbpg
from arbpg import PGFam, GamFam, logdensity, evaluate, choose_contour
from proposal import proposal
z = float(sys.argv[3])
import os
if os.environ.get('LEAN'):
    arbpg.LEAN = eval(os.environ['LEAN'])
for b in [float(v) for v in sys.argv[1].split(',')]:
    for psi in [float(v) for v in sys.argv[2].split(',')]:
        mean = float(b * pg1_mean(psi)); sd = float(mp.sqrt(b * pg1_var(psi)))
        x = mean + z * sd
        if x <= 0: x = mean / 20
        K = K_rule(b, psi); a, lam = proposal(b, psi, K)
        out = []
        for fam in (PGFam(b, psi), GamFam(a, lam)):
            if os.environ.get('TARGET'):
                # coarse level: the first h whose relative mathematical bound is below TARGET
                tgt = float(os.environ['TARGET']); ctx.prec = 53
                t0_, rho_ = choose_contour(fam, x); cache_ = {}; r = None
                for h_ in arbpg.H_SEQ:
                    try:
                        e_ = evaluate(fam, x, h_, t0_, rho_, cache_, lead_prec=106)
                    except RuntimeError:
                        break
                    Em_ = (e_['Ed'] + e_['Et']).upper()
                    if e_['S'].lower() > Em_ and float(Em_) / float(e_['S'].mid()) < tgt:
                        r = dict(h=h_, rad=float(Em_) / float(e_['S'].mid())); break
            else:
                r = logdensity(fam, x, lead_prec=106)
            if r is None:
                out.append(f'{fam.tag}: no enclosure'); continue
            cnt = {'lr': 0, 'cell': 0}
            olr, oc = fam.log_ratio, fam.cell_log_major
            def lr(D, _o=olr): cnt['lr'] += 1; return _o(D)
            def cl(*a_, _o=oc): cnt['cell'] += 1; return _o(*a_)
            fam.log_ratio, fam.cell_log_major = lr, cl
            ctx.prec = 53
            t0, rho = choose_contour(fam, x)
            t = time.time(); e = evaluate(fam, x, r['h'], t0, rho, {}, lead_prec=106); dt = time.time() - t
            ctx.prec = 53
            out.append(f"{fam.tag}: N={e['N']} nodes_total={cnt['lr']} cells={cnt['cell']} rad={r['rad']:.1e} t={dt:.2f}s")
        print(f"b={b:g} psi={psi:g} K={K} z={z:g} | " + " | ".join(out), flush=True)
