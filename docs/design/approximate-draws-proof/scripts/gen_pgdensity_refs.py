"""Write the engine evaluator's reference tables, tests/cpp/pgdensity_refs.inc: C initializers of hex-float
doubles, included by tests/cpp/test_pgdensity.cpp after it declares the row types.

- kFPins: log f, f = PG(b, psi), at the 500 pins of pins2.py (b in {0.2, 0.5, 1, 2.5, 8, 16, 64, 1e3, 1e4, 1e5},
  |psi| in {0, 8, 50, 100, 150}, z in {-4, -2, -1, 0, 0.5, 1.5, 3, 5, 8, 12}): x as pins2.py forms it, the
  reference (refs.logf_contour, 50 digits) rounded to a
  double-double, and arbpg's 53-bit enclosure with the lead in 106 bits (midpoint, radius). The reference is
  quadrature on a parabolic contour other than arbpg's and the engine's (pins3.py's rho^2 = max(T1 - pi^2/8,
  0.75 T1)); the contour rho^2 = 0.3 T1 and, for b <= 16, the left alternating series check it, their largest
  disagreement printed.
- kQPins: log q at the rule's K from the sampler's own rates (proposal.py) at b in {0.2, 2.5, 64, 1e3, 1e4},
  |psi| in {8, 150}, z in {-4, 0.5, 12}: K, r_1, r_K, alpha, beta, the reference independent of the parabolic
  contour (Moschopoulos's series where it converges, else the vertical line through the saddle), and arbpg's
  enclosure. Two parabolic contours other than arbpg's check the reference; their disagreement is printed.
- kKTable: the K rule's K (kd.K_rule, D in 40 digits) at |psi| in {0, 8, 50, 64, 100, 150} by b in {0.2, 1, 16,
  100, 1e3, 1e4, 1e5}, with b D at K - 1 and at K.
- kDd*: arguments for the double-double gate (ddouble.hpp's exp, log1p, log, sin, sqrt) and their values at 200
  digits rounded to double-doubles, over each function's range (below; each covers the reach the lead's formula
  gives at the pins, printed and written as comments).

Usage, from this folder, in the environment of ../requirements.txt:
    python3 -B gen_pgdensity_refs.py ../../../../tests/cpp/pgdensity_refs.inc [NPROC]
NPROC (default 2) worker processes."""
import math, os, random, sys, time
from multiprocessing import Pool
import mpmath as mp
from pgref import pg1_mean, pg1_var
from kd import K_rule, D_closed
from arbpg import PGFam, GamFam, logdensity, choose_contour
from proposal import proposal
import refs

F_B = [0.2, 0.5, 1, 2.5, 8, 16, 64, 1e3, 1e4, 1e5]
F_PSI = [0, 8, 50, 100, 150]
F_Z = [-4, -2, -1, 0, 0.5, 1.5, 3, 5, 8, 12]
Q_B = [0.2, 2.5, 64, 1e3, 1e4]
Q_PSI = [8, 150]
Q_Z = [-4, 0.5, 12]
K_PSI = [0, 8, 50, 64, 100, 150]
K_B = [0.2, 1, 16, 100, 1e3, 1e4, 1e5]
N_DD = 512
LEAD = 106


def pin_x(b, psi, z):
    mean = float(b * pg1_mean(psi)); sd = float(mp.sqrt(b * pg1_var(psi)))
    x = mean + z * sd
    return x if x > 0 else mean / 20


def dd_of(v):
    hi = float(v)
    lo = float(mp.mpf(v) - hi)
    return hi, lo


def f_pin(args):
    b, psi, z = args
    t = time.time()
    x = pin_x(b, psi, z)
    fam = PGFam(b, psi)
    r = logdensity(fam, x, lead_prec=LEAD)
    refl = []
    t0, _ = choose_contour(fam, x, 0.5)
    T1 = (t0 + fam.lam1_d()) / 2
    W0 = (psi / 2) ** 2 + t0 / 2
    for rho2 in (max(T1 - math.pi ** 2 / 8, 0.75 * T1), 0.3 * T1):
        rho = math.sqrt(rho2)
        refl.append(refs.logf_contour(b, psi, x, rho, rho * rho - W0, dps=50))
    if b <= 16:
        refl.append(refs.logf_left(b, psi, x, r['mid']))
    dis = max(float(abs(v - refl[0])) for v in refl[1:])
    return dict(b=b, psi=psi, z=z, x=x, ref=dd_of(refl[0]), mid=r['mid'], rad=r['rad'], dis=dis,
                t0=t0, secs=time.time() - t)


def q_pin(args):
    b, psi, z = args
    t = time.time()
    x = pin_x(b, psi, z)
    K = K_rule(b, psi)
    a, lam = proposal(b, psi, K)
    fam = GamFam(a, lam)
    r = logdensity(fam, x, lead_prec=LEAD)
    rm = refs.logq_mosch(a, lam, x)
    kind = 'mosch' if rm is not None else 'brom'
    ref = rm if rm is not None else refs.logq_bromwich(a, lam, x, dps=30)
    t0, _ = choose_contour(fam, x, 0.5)
    T1 = (t0 + fam.lam1_d()) / 2
    chk = []
    for rho2 in (max(T1 - math.pi ** 2 / 8, 0.75 * T1), 0.3 * T1):
        rho = math.sqrt(rho2)
        chk.append(refs.logq_contour(a, lam, x, rho, rho * rho - t0 / 2, dps=34))
    dis = max(float(abs(v - ref)) for v in chk)
    # the engine's lead takes log of products of (r + t0) / r over each run of equal shapes, rescaled by
    # powers of two whenever the exponent leaves [-256, 256]; record the products it would see
    prods = []
    runs = [(0, K), (K, K + 1)] if a[K] != b else [(0, K + 1)]
    for lo, hi in runs:
        p = 1.0
        for j in range(lo, hi):
            p *= (lam[j] + t0) / lam[j]
            m, e = math.frexp(p)
            if e > 256 or e < -256:
                p = m
        prods.append(p)
    return dict(b=b, psi=psi, z=z, x=x, K=K, r1=lam[0], rK=lam[K - 1], alpha=a[K], beta=lam[K],
                ref=dd_of(ref), kind=kind, mid=r['mid'], rad=r['rad'], dis=dis, prods=prods,
                secs=time.time() - t)


def k_cell(args):
    psi, b = args
    K = K_rule(b, psi)
    with mp.workdps(40):
        at = float(b * D_closed(psi, K))
        below = float(b * D_closed(psi, K - 1)) if K > 1 else float('inf')
    return dict(b=b, psi=psi, K=K, at=at, below=below)


def hx(v):
    if v == float('inf'):
        return 'kInf'
    return float(v).hex()


def dd_cases(fn, lo, hi, rng, near=None, logscale=True):
    """N_DD double-double arguments in [lo, hi]: half log-uniform in magnitude (or uniform), half uniform,
    the ends, and points near `near`; half the arguments with a nonzero low part."""
    xs = [lo, hi]
    if near is not None:
        for k in range(1, 60, 2):
            for s in (-1, 1):
                v = near + s * 2.0 ** -k
                if lo <= v <= hi:
                    xs.append(v)
    while len(xs) < N_DD:
        if logscale and len(xs) % 2 == 0:
            sgn = -1.0 if hi <= 0 else 1.0
            a, b = sorted((abs(lo) if lo != 0 else abs(hi) * 2.0 ** -60, abs(hi) if hi != 0 else abs(lo) * 2.0 ** -60))
            v = sgn * math.exp(rng.uniform(math.log(a), math.log(b)))
        else:
            v = rng.uniform(lo, hi)
        if lo <= v <= hi:
            xs.append(v)
    out = []
    with mp.workdps(200):
        for i, h in enumerate(xs):
            l = 0.0
            if i % 2 == 1 and h != 0.0:
                u = math.ulp(h)
                l = (rng.random() - 0.5) * u * 0.99
                if not (lo <= h + l <= hi):
                    l = 0.0
            xv = mp.mpf(h) + mp.mpf(l)
            out.append((h, l) + dd_of(fn(xv)))
    return out


def main(out, nproc):
    t = time.time()
    fargs = [(b, psi, z) for b in F_B for psi in F_PSI for z in F_Z]
    qargs = [(b, psi, z) for b in Q_B for psi in Q_PSI for z in Q_Z]
    kargs = [(psi, b) for psi in K_PSI for b in K_B]
    with Pool(nproc) as pool:
        qp = pool.map(q_pin, qargs, chunksize=1)
        print(f'q pins done, {time.time() - t:.0f} s', flush=True)
        kt = pool.map(k_cell, kargs, chunksize=1)
        print(f'K table done, {time.time() - t:.0f} s', flush=True)
        fp = pool.map(f_pin, fargs, chunksize=4)
        print(f'f pins done, {time.time() - t:.0f} s', flush=True)

    # the reach of the lead's arguments at the pins
    expArgs, sinArgs, sqrtArgs, logArgs = [], [], [], []
    for p in fp:
        A = abs(p['psi']) / 2
        W0 = A * A + p['t0'] / 2
        ws = [A]
        if W0 >= 0:
            ws.append(math.sqrt(W0)); sqrtArgs.append(W0)
        else:
            s = math.pi / 2 - math.sqrt(-W0); sinArgs.append(s); logArgs.append(math.sin(s))
        expArgs += [-2 * w for w in ws if w <= 325]
    for p in qp:
        logArgs += p['prods']
    reach = dict(exp=(min(expArgs), max(expArgs)), sin=(min(sinArgs), max(sinArgs)),
                 sqrt=(min(sqrtArgs), max(sqrtArgs)), log=(min(logArgs), max(logArgs)))
    # the ranges tested: the hard bounds of each call (exp's argument is -2 w with w <= 325, so log1p's
    # is in [e^-650, 1]; sin's is pi/2 - v, v in (0, pi/2); log's is sin's value or a product rescaled into
    # [2^-257, 2^257]); sin's and sqrt's lower ends and sqrt's upper end widened from the reach
    ranges = dict(exp=(-650.0, 0.0), log1p=(math.exp(-650.0), 1.0),
                  sin=(min(2.0 ** -40, reach['sin'][0] / 16), math.pi / 2),
                  log=(min(2.0 ** -257, reach['log'][0] / 16), 2.0 ** 257),
                  sqrt=(min(2.0 ** -60, reach['sqrt'][0] / 16), max(2.0 ** 20, 16 * reach['sqrt'][1])))
    rng = random.Random(20261107)
    cases = dict(
        exp=dd_cases(mp.exp, *ranges['exp'], rng, near=0.0, logscale=True),
        log1p=dd_cases(mp.log1p, *ranges['log1p'], rng, logscale=True),
        log=dd_cases(mp.log, *ranges['log'], rng, near=1.0, logscale=True),
        sin=dd_cases(mp.sin, *ranges['sin'], rng, logscale=True),
        sqrt=dd_cases(mp.sqrt, *ranges['sqrt'], rng, logscale=True))

    fdis = max(p['dis'] for p in fp); qdis = max(p['dis'] for p in qp)
    with open(out, 'w') as o:
        w = o.write
        w('// Generated by docs/design/approximate-draws-proof/scripts/gen_pgdensity_refs.py; do not edit.\n')
        w('// Regenerate from that folder: python3 -B gen_pgdensity_refs.py ../../../../tests/cpp/pgdensity_refs.inc\n')
        w(f'// mpmath {mp.__version__}. Largest reference disagreement: f {fdis:.1e}, q {qdis:.1e}.\n')
        w('// Reach at the pins: ' + ', '.join(f'{k} [{v[0]:.3e}, {v[1]:.3e}]' for k, v in reach.items()) + '.\n\n')
        w('// {b, psi, z, x, ref.hi, ref.lo, arbMid, arbRad}\n')
        w('static const FPin kFPins[] = {\n')
        for p in fp:
            w('  {' + ', '.join(hx(v) for v in (p['b'], p['psi'], p['z'], p['x'], *p['ref'], p['mid'], p['rad'])) + '},\n')
        w('};\n\n')
        w('// {b, psi, z, K, x, r1, rK, alpha, beta, ref.hi, ref.lo, arbMid, arbRad, moschopoulos}\n')
        w('static const QPin kQPins[] = {\n')
        for p in qp:
            w('  {' + ', '.join(hx(v) for v in (p['b'], p['psi'], p['z'])) + f", {p['K']}, " +
              ', '.join(hx(v) for v in (p['x'], p['r1'], p['rK'], p['alpha'], p['beta'], *p['ref'], p['mid'], p['rad'])) +
              f", {'true' if p['kind'] == 'mosch' else 'false'}" + '},\n')
        w('};\n\n')
        w('// {b, psi, K, b D(psi, K - 1), b D(psi, K)}, D in 40 digits\n')
        w('static const KCell kKTable[] = {\n')
        for c in kt:
            w(f"  {{{hx(c['b'])}, {hx(c['psi'])}, {c['K']}, {hx(c['below'])}, {hx(c['at'])}}},\n")
        w('};\n\n')
        for name, key in (('Exp', 'exp'), ('Log1p', 'log1p'), ('Log', 'log'), ('Sin', 'sin'), ('Sqrt', 'sqrt')):
            lo, hi = ranges[key]
            w(f'// {key} over [{hx(lo)}, {hx(hi)}]: {{x.hi, x.lo, y.hi, y.lo}}\n')
            w(f'static const double kDd{name}Range[2] = {{{hx(lo)}, {hx(hi)}}};\n')
            w(f'static const DdCase kDd{name}[] = {{\n')
            for c in cases[key]:
                w('  {' + ', '.join(hx(v) for v in c) + '},\n')
            w('};\n\n')
    worst = sorted(fp + qp, key=lambda p: -p['secs'])[:3]
    print(f'wrote {out}: {len(fp)} f pins, {len(qp)} q pins, {len(kt)} K cells, {N_DD} dd cases each; '
          f'reference disagreement f {fdis:.1e}, q {qdis:.1e}; reach {reach}; slowest '
          + ', '.join(f"(b={p['b']:g}, psi={p['psi']:g}, z={p['z']:g}) {p['secs']:.0f} s" for p in worst)
          + f'; {time.time() - t:.0f} s', flush=True)


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 2)
