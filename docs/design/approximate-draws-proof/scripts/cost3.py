"""Round-3 cost study: operation counts of one certified evaluation of log f and log q, under cell
configurations FIXED HERE BEFORE ANY EVALUATION (CONFIGS below; none is chosen per point), at
(b, psi, z) from the command line. For each configuration the full h-walk of logdensity runs under that
configuration, and the script prints that run's own radius, its rounding part and its mathematical part,
whether it reached the floor (mathematical part below the rounding part), the counts at the final h
(main nodes N, truncation nodes evaluated one by one, line-majorant cells) and the counts summed over
the whole h-walk. Term evaluations: an f node or cell is one term; a q node or cell is K+1 terms.
DEF (slow) runs only when no lean configuration reached the floor. The a-priori cascade RULE tries the configurations in the order given and stops at the first that
reaches the floor; its cost is the sum over the configurations tried.
ONLY=name runs that configuration alone (always, DEF included).
Usage: cost3.py b1,b2 psi1,psi2 z1,z2"""
import sys, time
import mpmath as mp
from flint import ctx
from pgref import pg1_mean, pg1_var
from kd import K_rule
import arbpg
from arbpg import PGFam, GamFam, logdensity
from proposal import proposal

CONFIGS = [
    ('L10', dict(fr_p=[0.5], am=[-0.5], ncell=(10, 40), maxdepth=0)),      # second referee's round-3 setting
    ('L32', dict(fr_p=[0.5, 0.2], am=[-0.5, 1.0], ncell=(32, 32), maxdepth=2)),
    ('L8', dict(fr_p=[0.3], am=[-0.5], ncell=(8, 8), maxdepth=1)),
    ('DEF', None),                                                          # arbpg's default
]
RULE = ['L10', 'L32', 'DEF']

def run(fam, x, cfg):
    arbpg.LEAN = cfg
    cnt = {'lr': 0, 'cell': 0}
    olr, oc = fam.log_ratio, fam.cell_log_major
    fam.log_ratio = lambda D, _o=olr: (cnt.__setitem__('lr', cnt['lr'] + 1), _o(D))[1]
    fam.cell_log_major = lambda *a_, _o=oc: (cnt.__setitem__('cell', cnt['cell'] + 1), _o(*a_))[1]
    t = time.time()
    r = logdensity(fam, x, lead_prec=106)
    walk = dict(cnt); dt = time.time() - t
    out = None
    if r is not None:
        # recount at the final h alone, fresh cache (what a build that picks h directly would pay)
        cnt['lr'] = cnt['cell'] = 0
        ctx.prec = 53
        arbpg.evaluate(fam, x, r['h'], r['t0'], r['rho'], {}, lead_prec=106)
        out = dict(rad=r['rad'], rround=r['rad_round'], emath=r['Emath'] * 1.0, N=r['N'], h=r['h'],
                   nodes=cnt['lr'], cells=cnt['cell'], wnodes=walk['lr'], wcells=walk['cell'], t=dt,
                   floor=(r['Emath'] <= r['rad_round']))
    fam.log_ratio, fam.cell_log_major = olr, oc
    arbpg.LEAN = None
    return out

bs = [float(v) for v in sys.argv[1].split(',')]
psis = [float(v) for v in sys.argv[2].split(',')]
zs = [float(v) for v in sys.argv[3].split(',')]
import os
if os.environ.get('ONLY'):
    CONFIGS = [c for c in CONFIGS if c[0] == os.environ['ONLY']]   # e.g. ONLY=DEF: the coarsest certified h
print('CONFIGS', CONFIGS, 'RULE', RULE, flush=True)
for b in bs:
    for psi in psis:
        K = K_rule(b, psi); a, lam = proposal(b, psi, K)
        for z in zs:
            mean = float(b * pg1_mean(psi)); sd = float(mp.sqrt(b * pg1_var(psi)))
            x = mean + z * sd
            if x <= 0: x = mean / 20
            for fam, w in ((PGFam(b, psi), 1), (GamFam(a, lam), K + 1)):
                ok_any = False
                for name, cfg in CONFIGS:
                    if name == 'DEF' and ok_any:
                        print(f'{fam.tag} b={b:g} psi={psi:g} K={K} z={z:g} cfg=DEF skipped (a lean configuration reached the floor)', flush=True)
                        continue
                    o = run(fam, x, cfg)
                    ok_any = ok_any or (o is not None and o['floor'])
                    if o is None:
                        print(f'{fam.tag} b={b:g} psi={psi:g} K={K} z={z:g} cfg={name} NO ENCLOSURE', flush=True)
                        continue
                    tot = o['nodes'] + o['cells']
                    print(f"{fam.tag} b={b:g} psi={psi:g} K={K} z={z:g} cfg={name} rad={o['rad']:.2e} "
                          f"round={o['rround']:.2e} Emath={o['emath']:.1e} floor={o['floor']} N={o['N']} h={o['h']:.3g} "
                          f"nodes={o['nodes']} cells={o['cells']} total={tot} terms={tot * w} "
                          f"walk_nodes={o['wnodes']} walk_cells={o['wcells']} t={o['t']:.1f}", flush=True)
