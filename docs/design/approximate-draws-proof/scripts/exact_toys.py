"""Exact finite-state checks (rational arithmetic, no rounding) of the chain
statements of proof.tex: Theorems 4-7 and Propositions 2.4, 7.3, 9.1.

Each toy builds the exact transition matrix of a sweep and reports
  max |pi P - pi|  (0 means pi is invariant, exactly)
and, where a scheme is not invariant, the total variation distance between
its stationary distribution and the target. Arbitrary positive rational
distributions are used, so a zero is structural, not a coincidence.

Run: python exact_toys.py
"""
from fractions import Fraction as Fr
import itertools, random

random.seed(20261105)


def rpos(n):
    """random positive rational distribution on n points"""
    w = [Fr(random.randint(1, 9)) for _ in range(n)]
    s = sum(w)
    return [x / s for x in w]


def imh(target, prop):
    """exact independence MH kernel on a finite set"""
    n = len(target)
    P = [[Fr(0)] * n for _ in range(n)]
    for x in range(n):
        for y in range(n):
            if y == x:
                continue
            a = min(Fr(1), target[y] * prop[x] / (target[x] * prop[y]))
            P[x][y] = prop[y] * a
        P[x][x] = 1 - sum(P[x])
    return P


def stationary(P):
    """solve pi P = pi, sum pi = 1, exactly"""
    n = len(P)
    A = [[P[j][i] - (1 if i == j else 0) for j in range(n)] for i in range(n)]
    A[-1] = [Fr(1)] * n
    b = [Fr(0)] * (n - 1) + [Fr(1)]
    M = [row[:] + [bb] for row, bb in zip(A, b)]
    for c in range(n):
        p = next(r for r in range(c, n) if M[r][c] != 0)
        M[c], M[p] = M[p], M[c]
        piv = M[c][c]
        M[c] = [v / piv for v in M[c]]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c]
                M[r] = [a - f * bb for a, bb in zip(M[r], M[c])]
    return [M[i][n] for i in range(n)]


def residual(pi, P):
    n = len(pi)
    return max(abs(sum(pi[i] * P[i][j] for i in range(n)) - pi[j]) for j in range(n))


def tv(p, q):
    return sum(abs(a - b) for a, b in zip(p, q)) / 2


def matmul(A, B):
    n = len(A)
    return [[sum(A[i][k] * B[k][j] for k in range(n)) for j in range(n)] for i in range(n)]


# --------------------------------------------------------------------------
# Toy 1 (Theorem 5): estimated shape r drawn collapsed over omega, then omega
# fresh where r moved and by IMH where r stayed; then the fits f given (r, omega).
def toy_nb():
    F, R, W = 3, 3, 3
    joint = rpos(F * R * W)
    idx = lambda f, r, w: (f * R + r) * W + w
    pi = joint
    pfr = lambda f, r: sum(joint[idx(f, r, w)] for w in range(W))
    pf = lambda f: sum(pfr(f, r) for r in range(R))
    pw = lambda f, r: [joint[idx(f, r, w)] / pfr(f, r) for w in range(W)]
    prop = {(f, r): rpos(W) for f in range(F) for r in range(R)}
    out = {}
    for rule in ("designed", "always_imh"):
        n = F * R * W
        P1 = [[Fr(0)] * n for _ in range(n)]
        for f, r, w in itertools.product(range(F), range(R), range(W)):
            for r2 in range(R):
                pr = pfr(f, r2) / pf(f)
                if r2 != r and rule == "designed":
                    for w2 in range(W):
                        P1[idx(f, r, w)][idx(f, r2, w2)] += pr * pw(f, r2)[w2]
                else:
                    K = imh(pw(f, r2), prop[(f, r2)])
                    for w2 in range(W):
                        P1[idx(f, r, w)][idx(f, r2, w2)] += pr * K[w][w2]
        P2 = [[Fr(0)] * n for _ in range(n)]
        for f, r, w in itertools.product(range(F), range(R), range(W)):
            tot = sum(joint[idx(f2, r, w)] for f2 in range(F))
            for f2 in range(F):
                P2[idx(f, r, w)][idx(f2, r, w)] = joint[idx(f2, r, w)] / tot
        P = matmul(P1, P2)
        st = stationary(P)
        out[rule] = (residual(pi, P), tv(st, pi))
    return out


# --------------------------------------------------------------------------
# Toy 2 (Proposition 2.4): an unresolved decision replaced by a fresh exact
# draw. Decision unresolved iff U in (a - d, a + d), d depending on (x, y).
def toy_fresh_on_unresolved():
    n = 3
    pi = [Fr(1, 2), Fr(1, 3), Fr(1, 6)]
    q = [Fr(1, 3)] * 3
    d = {(x, y): Fr(1, 100) * (1 + x + 2 * y) for x in range(n) for y in range(n)}
    P = [[Fr(0)] * n for _ in range(n)]
    for x in range(n):
        for y in range(n):
            if y == x:
                P[x][x] += q[y]
                continue
            a = min(Fr(1), pi[y] * q[x] / (pi[x] * q[y]))
            lo, hi = max(Fr(0), a - d[(x, y)]), min(Fr(1), a + d[(x, y)])
            unres = hi - lo if a < 1 else Fr(0)
            acc = lo if a < 1 else Fr(1)
            P[x][y] += q[y] * acc
            P[x][x] += q[y] * (1 - acc - unres)
            for z in range(n):
                P[x][z] += q[y] * unres * pi[z]
    st = stationary(P)
    return residual(pi, P), tv(st, pi)


# --------------------------------------------------------------------------
# Toy 3 (Proposition 9.1): the interleaved two-category augmentation. Target
# p(t1, t2); category k's latent w_k ~ c_k(w | t1, t2); forest k draws
# t_k ~ p(t_k | t_-k, w_k) from p(t) c_k(w_k | t). Exact: w_k fresh before
# forest k. Route: w_k stored and moved by IMH against c_k(. | t).
def toy_multinomial():
    T, W = 3, 3
    pt = rpos(T * T)
    c = {k: {(t1, t2): rpos(W) for t1 in range(T) for t2 in range(T)} for k in (0, 1)}
    prop = {k: {(t1, t2): rpos(W) for t1 in range(T) for t2 in range(T)} for k in (0, 1)}
    res = {}
    # exact: chain on t alone
    def step_exact(k):
        n = T * T
        P = [[Fr(0)] * n for _ in range(n)]
        for t1, t2 in itertools.product(range(T), range(T)):
            for w in range(W):
                pw = c[k][(t1, t2)][w]
                if k == 0:
                    tot = sum(pt[s * T + t2] * c[0][(s, t2)][w] for s in range(T))
                    for s in range(T):
                        P[t1 * T + t2][s * T + t2] += pw * pt[s * T + t2] * c[0][(s, t2)][w] / tot
                else:
                    tot = sum(pt[t1 * T + s] * c[1][(t1, s)][w] for s in range(T))
                    for s in range(T):
                        P[t1 * T + t2][t1 * T + s] += pw * pt[t1 * T + s] * c[1][(t1, s)][w] / tot
        return P
    P = matmul(step_exact(0), step_exact(1))
    res["fresh"] = (residual(pt, P), tv(stationary(P), pt))
    # route: state (t1, t2, w1, w2)
    n = T * T * W * W
    idx = lambda t1, t2, w1, w2: ((t1 * T + t2) * W + w1) * W + w2
    def step_route(k):
        P = [[Fr(0)] * n for _ in range(n)]
        for t1, t2, w1, w2 in itertools.product(range(T), range(T), range(W), range(W)):
            wk = (w1, w2)[k]
            K = imh(c[k][(t1, t2)], prop[k][(t1, t2)])
            for wn in range(W):
                pw = K[wk][wn]
                if pw == 0:
                    continue
                if k == 0:
                    tot = sum(pt[s * T + t2] * c[0][(s, t2)][wn] for s in range(T))
                    for s in range(T):
                        P[idx(t1, t2, w1, w2)][idx(s, t2, wn, w2)] += pw * pt[s * T + t2] * c[0][(s, t2)][wn] / tot
                else:
                    tot = sum(pt[t1 * T + s] * c[1][(t1, s)][wn] for s in range(T))
                    for s in range(T):
                        P[idx(t1, t2, w1, w2)][idx(t1, s, w1, wn)] += pw * pt[t1 * T + s] * c[1][(t1, s)][wn] / tot
        return P
    P = matmul(step_route(0), step_route(1))
    st = stationary(P)
    marg = [sum(st[idx(t1, t2, a, b)] for a in range(W) for b in range(W))
            for t1 in range(T) for t2 in range(T)]
    res["route"] = (None, tv(marg, pt))
    return res


# --------------------------------------------------------------------------
# Toy 4 (Theorem 7): an outer parameter th (offset), the fits f, one row's
# latent w. Target pi(th, f, w). Sweep: outer step on th, then the latent
# step, then f ~ pi(f | th, w).
#   contract I: th ~ pi(th | f)   (w integrated out; w is stale afterwards)
#   contract C: th ~ pi(th | f, w)
#   latent: fresh w ~ pi(w | th, f), or IMH against pi(w | th, f).
def toy_setters():
    A, F, W = 2, 3, 3
    joint = rpos(A * F * W)
    idx = lambda a, f, w: (a * F + f) * W + w
    prop = {(a, f): rpos(W) for a in range(A) for f in range(F)}
    n = A * F * W
    def cond_w(a, f):
        s = sum(joint[idx(a, f, w)] for w in range(W))
        return [joint[idx(a, f, w)] / s for w in range(W)]
    out = {}
    for contract in ("I", "C"):
        for latent in ("fresh", "route"):
            P1 = [[Fr(0)] * n for _ in range(n)]
            for a, f, w in itertools.product(range(A), range(F), range(W)):
                if contract == "I":
                    pa = [sum(joint[idx(a2, f, ww)] for ww in range(W)) for a2 in range(A)]
                else:
                    pa = [joint[idx(a2, f, w)] for a2 in range(A)]
                s = sum(pa)
                for a2 in range(A):
                    P1[idx(a, f, w)][idx(a2, f, w)] += pa[a2] / s
            P2 = [[Fr(0)] * n for _ in range(n)]
            for a, f, w in itertools.product(range(A), range(F), range(W)):
                cw = cond_w(a, f)
                K = imh(cw, prop[(a, f)]) if latent == "route" else [cw] * W
                for w2 in range(W):
                    P2[idx(a, f, w)][idx(a, f, w2)] += K[w][w2]
            P3 = [[Fr(0)] * n for _ in range(n)]
            for a, f, w in itertools.product(range(A), range(F), range(W)):
                s = sum(joint[idx(a, f2, w)] for f2 in range(F))
                for f2 in range(F):
                    P3[idx(a, f, w)][idx(a, f2, w)] = joint[idx(a, f2, w)] / s
            P = matmul(matmul(P1, P2), P3)
            st = stationary(P)
            out[(contract, latent)] = (residual(joint, P), tv(st, joint))
    return out


# --------------------------------------------------------------------------
# Toy 5 (Theorem 6, Proposition 6.2): a per-row choice between a fresh draw
# and IMH made from the conditioning variables keeps pi(w | f); a choice
# made from w itself does not.
def toy_choice():
    W = 4
    target = rpos(W)
    prop = rpos(W)
    K = imh(target, prop)
    fresh = [target] * W
    bad = [[(fresh[x][y] if x >= 2 else K[x][y]) for y in range(W)] for x in range(W)]
    return residual(target, K), residual(target, fresh), residual(target, bad), tv(stationary(bad), target)


if __name__ == "__main__":
    nb = toy_nb()
    print("Toy 1, estimated r: designed rule  |piP - pi| =", float(nb["designed"][0]),
          " TV(stationary, pi) =", float(nb["designed"][1]))
    print("                    IMH after any r |piP - pi| =", float(nb["always_imh"][0]),
          " TV(stationary, pi) =", float(nb["always_imh"][1]))
    r, t = toy_fresh_on_unresolved()
    print("Toy 2, fresh draw on unresolved decisions: |piP - pi| =", float(r), " TV =", float(t))
    m = toy_multinomial()
    print("Toy 3, interleaved augmentation: fresh  |piP - pi| =", float(m["fresh"][0]),
          " TV =", float(m["fresh"][1]))
    print("                                 route  TV(theta-marginal, p) =", float(m["route"][1]))
    s = toy_setters()
    for key, (r, t) in sorted(s.items()):
        print("Toy 4, contract %s, latent %-5s |piP - pi| = %.3g  TV = %.3g" % (key[0], key[1], float(r), float(t)))
    a, b, c, d = toy_choice()
    print("Toy 5, IMH %g, fresh %g, choice read from w: |piP - pi| = %.3g, TV = %.3g" % (float(a), float(b), float(c), float(d)))
