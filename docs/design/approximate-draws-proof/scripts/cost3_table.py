"""Aggregate cost3_*.out into the per-row-step table of Sec. 16 item 8.
Per evaluation the a-priori cascade RULE of cost3.py is charged: the configurations tried in order until one
reaches the floor, their final-h counts summed. A row step is f(y) + q(y) + q(x) (three fresh evaluations;
the cached f(x) is moved by eq. (4)). Term evaluations: an f node or cell is one term, a q node or cell K+1
terms. Kernel times (termbench.c, this machine, plain double, no radius bookkeeping): f node 116 ns (f cells
priced the same, an upper figure), q node term 12.6 ns, q cell term 3.1 ns. Today's sum: pgbench.c on this
machine, 136 to 160 ns per PG(1) draw (157 at psi = 0, 136 at 150), times b. Gammas: the design's 41 ns.
The 'main' column is the main sum alone (N nodes, N(K+1) for q) at the coarsest h any configuration
certified: what any certificate costs at least, at that h.
Usage: cost3_table.py files..."""
import re, sys
NSF, NSQ, NSQC = 116.0, 12.6, 3.1
RULE = ['L10', 'L32', 'DEF']
rows = {}
pat = re.compile(r'^(PG|Gam) b=(\S+) psi=(\S+) K=(\d+) z=(\S+) cfg=(\S+) rad=(\S+) round=(\S+) Emath=(\S+) floor=(\S+) N=(\d+) h=(\S+) nodes=(\d+) cells=(\d+)')
for fn in sys.argv[1:]:
    for l in open(fn):
        m = pat.match(l)
        if not m: continue
        tag, b, psi, K, z, cfg, rad, rr, em, fl, N, h, nd, cl = m.groups()
        rows.setdefault((float(b), float(psi), int(K), float(z), tag), {})[cfg] = dict(
            rad=float(rad), floor=fl == 'True', N=int(N), h=float(h), nodes=int(nd), cells=int(cl))
def rule(d, w, nsn, nsc):
    terms = ns = 0; used = []
    for cfg in RULE:
        if cfg not in d: return None
        e = d[cfg]; terms += (e['nodes'] + e['cells']) * w; ns += (e['nodes'] * nsn + e['cells'] * nsc) * w
        used.append(cfg)
        if e['floor']: return terms, ns, '+'.join(used), e['rad']
    return None
print('| b | psi | K | z | f terms | q terms | step terms | step us | main-sum terms | main-sum us | gammas us | sum us | rad f, q |')
print('|---|---|---|---|---|---|---|---|---|---|---|---|---|')
for (b, psi, K, z) in sorted({k[:4] for k in rows}):
    df, dq = rows.get((b, psi, K, z, 'PG'), {}), rows.get((b, psi, K, z, 'Gam'), {})
    rf, rq = rule(df, 1, NSF, NSF), rule(dq, K + 1, NSQ, NSQC)
    if rf is None or rq is None:
        print(f'| {b:g} | {psi:g} | {K} | {z:g} | rule failed |'); continue
    step = rf[0] + 2 * rq[0]; us = (rf[1] + 2 * rq[1]) / 1000
    nf = min(e['N'] for e in df.values() if e['floor']); nq = min(e['N'] for e in dq.values() if e['floor'])
    main = nf + 2 * nq * (K + 1); mus = (nf * NSF + 2 * nq * (K + 1) * NSQ) / 1000
    gam = (K + 1) * 0.041
    s = f'{(0.157 if psi == 0 else 0.136) * b:.3g}' if abs(b - round(b)) < 1e-9 else 'none'
    print(f'| {b:g} | {psi:g} | {K} | {z:g} | {rf[0]} ({rf[2]}) | {rq[0]} ({rq[2]}) | {step:.3g} | {us:.3g} | '
          f'{main:.3g} | {mus:.3g} | {gam:.3g} | {s} | {rf[3]:.1e}, {rq[3]:.1e} |')
