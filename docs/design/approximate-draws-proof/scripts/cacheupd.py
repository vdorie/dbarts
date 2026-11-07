"""The cached-f(x) route: log f_{psi'}(x) = log f_psi(x) + b[phi(A'^2) - phi(A^2)] - x(sigma' - sigma),
phi(W) = log cosh sqrt(W), A = |psi|/2, sigma = psi^2/2 (an exact identity, the tilt). The bracket is
formed cancellation-free as b*D with D from eq. (cf) on the real path (PGFam's lead machinery), the
x-term from psi' - psi exactly. Prints the update's radius at 53 bits, and checks the updated
enclosure against a fresh evaluation at psi'. Usage: cacheupd.py b psi z dpsi1,dpsi2"""
import sys
import mpmath as mp
from flint import arb, ctx
from pgref import pg1_mean, pg1_var
from arbpg import PGFam, logdensity

def logcosh_diff(a1, a2):
    """log cosh(a2) - log cosh(a1), real a1, a2 >= 0, cancellation-free:
    (a2 - a1) + log1p(e^{-2 a2}) - log1p(e^{-2 a1}) for the large part, i.e. lambda(a2) - lambda(a1);
    with lambda(w) = w - log 2 + log1p(e^{-2w}) (Prop. cf's direct form; here real, no branch)."""
    A1, A2 = arb(a1), arb(a2)
    return (A2 - A1) + (-2 * A2).exp().log1p() - (-2 * A1).exp().log1p()

import os
UPREC = int(os.environ.get('UPREC', '53'))  # precision of the update scalar (106: double-double)
b, psi, z = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])
mean = float(b * pg1_mean(psi)); sd = float(mp.sqrt(b * pg1_var(psi))); x = mean + z * sd
ctx.prec = 53
r0 = logdensity(PGFam(b, psi), x, lead_prec=106)
print(f"b={b:g} psi={psi:g} z={z:g} x={x:.6g}: fresh rad {r0['rad']:.2e}")
for dpsi in [float(v) for v in sys.argv[4].split(',')]:
    ctx.prec = 53
    p2 = psi + dpsi
    ctx.prec = UPREC
    upd = arb(b) * logcosh_diff(abs(psi) / 2, abs(p2) / 2) - arb(x) * (arb(p2) - arb(psi)) * (arb(p2) + arb(psi)) / 2
    lf = r0['logp'] + upd          # the addition at UPREC too (double-double in the engine: R5, round 3)
    ctx.prec = 53
    r1 = logdensity(PGFam(b, p2), x, lead_prec=106)
    ok = lf.overlaps(r1['logp'])
    print(f"  dpsi={dpsi:g}: update rad {float(upd.rad()):.2e}, updated rad {float(lf.rad()):.2e}, fresh rad {r1['rad']:.2e}, "
          f"|upd mid - fresh mid| {abs(float(lf.mid()) - r1['mid']):.1e}, overlap {ok}")
