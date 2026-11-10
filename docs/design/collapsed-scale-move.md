# The collapsed rescaling move under logistic

Status: BUILT ON A BRANCH for pricing, 2026-11-10 (collapsed-move-rebuild; plan
[collapsed-scale-move.md](../plans/collapsed-scale-move.md), step 0 of k-move-pricing). This note holds the
move's exactness derivation, checked by the engine review and by one independent referee before landing;
the rest of the note (schedule, measurements, cost) is written at landing.

## Exactness

The construction is ours: a new combination of published steps, each used as published. These are the
generalized Gibbs step on a group of transformations (Liu and Sabatti 2000, Theorems 1 and 2), the Polya-Gamma
identity (Polson, Scott and Windle 2013, Theorem 1) and the stepping-out slice update (Neal 2003; figure and
section numbers below are those of its preprint, arXiv physics/0009028). The move
is [`Chain::drawCollapsedRescale`](../../src/bartcore/chain.hpp).

**The collapsed target.** The move acts on x = (T, mu, k): the trees T, the M occupied leaves mu of the
orbit's forests (one plain constant-leaf forest, or every forest that shares one prior total) and their one
k. A leaf is occupied when at least one row reaches it, kept by the mask or not; the k draw counts the same
leaves, and an empty leaf carries no prior or likelihood term. Leaf j has prior N(0, (c_j / k)^2)
with c_j > 0 fixed by the model: the node scale, its tree's factor under blocks() and, on the sharers,
sqrt(1 / K). The prior on k is s chi(nu), density proportional to k^(nu-1) exp(-k^2 / (2 s^2)), with
nu > 0, which chi() enforces. Row i in the
likelihood (kept by the mask, of positive count w_i and of at least one trial n_i, n_i = 1 for a binary row)
has log-odds eta_i = f_i(mu) + u_i + o_i. Here f_i is the orbit forests' summed fit, linear in mu; o_i is the
offset; and u_i is every other forest's contribution, a function of parameters the move does not touch.
Theorem 1 of Polson, Scott and Windle, with b = w_i n_i and a = w_i y_i, integrates each latent omega_i out of
the augmented joint and leaves the binomial likelihood exp(l_i(eta_i)), l_i = w_i log dbinom(y_i, n_i,
plogis(eta_i)) up to the constant choose(n_i, y_i)^w_i 2^(-w_i n_i). The same theorem states that
omega_i | eta_i ~ PG(w_i n_i, eta_i). The R layer enforces whole counts w_i and whole trials n_i, so the code's
sum of w_i n_i draws from PG(1, eta_i) is a draw from PG(w_i n_i, eta_i). So the move's target, the marginal
of x, is proportional to

    p(T) prod_j N(mu_j; 0, (c_j / k)^2) k^(nu-1) exp(-k^2 / (2 s^2)) prod_i exp(l_i(f_i(mu) + u_i + o_i)).

**The group and the step.** Gamma = {r > 0} acts by g_r(T, mu, k) = (T, r mu, k / r), with g_r g_q = g_(rq).
This is the scale group of Liu and Sabatti's section 2: locally compact, unimodular, Haar measure
L(dr) = dr / r. The Jacobian of (mu, k) -> (r mu, k / r) is r^M r^-1 = r^(M-1). Their Theorem 1 draws r from
p_x(r) dr proportional to pi(g_r x) r^(M-1) dr / r. The state space changes dimension with the trees, but
g_r fixes the trees and every parameter outside (mu, k), so it maps each tree configuration's component onto
itself, with Jacobian r^(M-1) there. The theorems apply on each component, under counting measure on the trees
times Lebesgue measure on (mu, k), and the step is their mixture. At g_r x:
- each leaf's density is r^-1 times its value at x: N(r mu_j; 0, (c_j r / k)^2) = r^-1 N(mu_j; 0, (c_j / k)^2);
- the k prior is (k / r)^(nu-1) exp(-(k / (r s))^2 / 2), that is, r^-(nu-1) times k^(nu-1) and the
  exponential at k / r;
- the likelihood is exp(sum_i l_i(r f_i + u_i + o_i)), u_i and o_i unmoved.

The powers multiply to r^(M-1) r^-M r^-(nu-1) r^-1 dr = r^-(nu+1) dr = r^-nu dv in v = log r, leaving the log
density

    log p_x(v) = -nu v - (k e^-v / s)^2 / 2 + sum_i l_i(e^v f_i + u_i + o_i) + const,

in which M cancels. It is proper for a finite s: as v -> -inf the k term goes to -inf like -e^(-2v), and as
v -> +inf the term -nu v goes to -inf while each l_i <= 0. With every row out it is k's prior along the orbit.
At s = inf it is improper (as v -> -inf the sum stays bounded and -nu v grows), so the move declines there.
At nu = inf, which chi() admits, log p_x(0) is not finite and the move declines without drawing.

**The slice step and condition (1).** The move replaces Theorem 1's exact draw by one transition
T_x(v, v') of [`sliceFromZero`](../../src/bartcore/chain.hpp) from v = 0, the identity. That transition is
Neal's stepping out (his Figure 3) and shrinkage (his Figure 5) on the log scale, the level
log p_x(0) - Exp(1) standing for his vertical draw from (0, p_x(0)). His section 4.3 shows the update
leaves p_x invariant, the step limit split at random as his Figure 3 does. Liu and Sabatti's Theorem 2 lets such a
transition stand in for the exact draw if it satisfies their condition (1), T_x(g, g') = T_(g0^-1 x)(g g0,
g' g0). In v, with c = log g0, that reads T_x(a, b) = T_(g0^-1 x)(a + c, b + c). The group is abelian, so
g_(r g0)(g0^-1 x) = g_r x, and log p_(g0^-1 x)(a + c) = log p_x(a) + const. The condition therefore asks
that the slice transition started at a + c on the target shifted by c be the transition started at a,
shifted by c. That holds:
- the level is the log density at the current point less an Exp(1) draw;
- the initial interval [a - wU, a - wU + w] and the split of the step limit into J and m - 1 - J are placed
  relative to the current point;
- shrinkage moves an end to each rejected point;
- the target enters only through comparisons with the level.

Every one of these commutes with a shift. So by Theorem 2, g_(e^v')(x) is distributed as the marginal of x
whenever x is.

**The latents.** The move reads no latent and writes none. The sweep's augmentation draw follows it before
anything reads a latent and draws each omega_i from PG(w_i n_i, eta_i) at the moved state, its exact full
conditional. If (x, omega) is distributed as the joint, x is distributed as the marginal; after the move x'
is too; and omega' drawn from omega | x' makes (x', omega') distributed as the joint.

**The sharers.** Each sharer's leaves have sd c_j / k with c_j including its sqrt(1 / K), so the group's one k
scales all of them alike. A sharer's coefficient is held at 1 on an all-ones basis, so eta_i is linear in the
sharers' leaves with slope f_i, as the density above takes it. u_i, the forests outside the group, holds
their leaves and coefficients, which g_r does not touch and whose priors do not involve the group's k. The
derivation above therefore holds with M the pooled count of occupied leaves and f_i the sharers' summed fit.

**Why the standardized leaves are invariant.** A leaf's standardized value is mu_j k / c_j. Under g_r it
becomes (r mu_j)(k / r) / c_j = mu_j k / c_j. That holds exactly because every leaf's scale is a fixed
multiple c_j of 1 / k, with c_j free of mu and k: a per-tree factor under blocks(), a sharer's sqrt(1 / K)
or 1. So the leaf prior contributes a pure r^-M and nothing else, and the move needs to read no leaf scale.
A leaf scale that did not move with 1 / k would leave a factor of r in the leaf prior that the density
above omits.

**The schedule and the declines.** A sweep takes the move on a coin of 1/4, drawn independently of the
state. The move declines with no draw where the family, the forest shape, whether k is drawn or an
infinite s says so; each of these is constant on an orbit. It also declines on a stale leaf map, which is
history rather than state: a stale map comes only from a prior draw the caller asks for, so the sweeps on
which the move declines for it are fixed by the call sequence and independent of the state. Last, it
declines where k is not finite and positive or where log p_x(0) is not finite, which happens only on a null
set. A mixture of the move and the
identity, weighted by a function constant on orbits, keeps the marginal.
