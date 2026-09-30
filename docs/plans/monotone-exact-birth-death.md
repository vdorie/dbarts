# monotone-exact-birth-death: the monotone chain targets the chosen prior exactly

Status: PLANNED 2026-09-29 (dec-B144), revised after a blind critique and a few-tree measurement. The maintainer
then ruled that the package offers both monotone priors, "leaf" and "joint", chosen by monotone(prior = )
(dec-B145 to dec-B148); the default stays open for the feel study below. Every move is then counted on the finer
tree's side (Counting: algorithm), so a death proposal never counts the component its merge creates; a later
move touching an accepted merge does. Under "leaf" a slow count warns and never stops the run (dec-B149): there
is no count limit. A blind critique of the counting then fixed its gaps, and a whole-plan critique found it not
ready: the plan now stages the work in seven reviewed commits (Staging), with a checkpoint on the corrected
engine's counts after the fourth, and records three orchestrator calls (dec-A128). The monotone() signature is
monotone(directions, prior = ) (dec-B150). Derivation and gate verified on an R prototype; not implemented.

agent: opus (engine numerics: move seam, order counting, exact pair redraw, gate)
rng: posterior-changing for every fit with an active monotone constraint (all of its draws move, prior draws
included); unconstrained fits byte-identical, since every engine change lives in the monotone instantiation
window: before 1.0-0 (TODO monotone-exact-birth-death)
budget: ~2,000 lines (engine ~700, bridge and R ~330, tests/cpp ~470, tinytest ~200, gate and SBC wiring ~150,
docs ~150; the gate script and the feel study's script excluded). The second prior adds ~150-250 of that, the
monotone() constructor with its vocabulary and caller sweep ~160, counting on the finer tree's side ~100, the
slow-count warning with the interrupt, allocation and rebuild paths ~150, the SBC arms ~120, and the lazy
cache ~80. Plan estimates have run 1.5-2x low: expect up to ~4,000.

## Goal

The package offers two monotone priors, chosen by monotone(prior = ) (dec-B145, dec-B148), and the chain
targets the chosen one exactly:

- "leaf": the CGM tree prior and, given the tree, iid (c-inflated) normal leaves restricted to the cone C(T),
  normalized per tree by Z_T = P(unconstrained leaves lie in C(T)). The constraint restricts the leaf-value
  prior's support; the tree prior is unchanged.
- "joint": p(T, M) proportional to p_CGM(T) prod phi 1{M in C(T)}, the tree's structure and leaf values
  conditioned on the cone together, so the tree marginal is p_CGM(T) Z_T.

Which is the default is open (Default: feel study). Exactness covers four pieces under both:

- the birth/death acceptance;
- the redraw of a birth's two children;
- the treatment of leaves left empty;
- the prior leaf draw.

No path reaches a structural move other than birth/death. Every state the sampler holds lies in the cone. An
exact gate over multi-split trees, N-shaped leaf orders included, fails the current engine decisively and
guards the fix.

## Context

- The defect. [`MonotoneConstantGaussianLeaf::oneLeafLogMarginal`](../../src/bartcore/model.hpp) and
  [`MonotoneConstantGaussianLeaf::twoLeafCoupledLogMarginal`](../../src/bartcore/model.hpp) divide the touched
  leaves' constrained marginal by d, their prior cone mass given the frozen neighbours. d agrees with the
  whole-tree normalizer only when no touched leaf has a frozen constrained neighbour. That covers every move of
  part (a) of the existing gate ([9. Gates](../design/monotone.md#9-gates)), which is why it passed.
- Derivation. pi(T, M) is proportional to p(T) prod_k phi_k(mu_k) 1{M in C(T)} lik / Z_T, with
  phi_k = N(0, s_k^2) and s_k = c scale / k for a leaf with a constrained neighbour
  ([`MonotoneConstantGaussianLeaf::priorSd`](../../src/bartcore/model.hpp)). Integrating the touched leaves
  given the rest ("same") gives pi(T, same) = p(T) / Z_T * prod_same phi lik * I_T(same), where I_T is today's
  score numerator with the division removed.
  - Birth, exactly:
    log alpha = log[p(T*) / p(T0)] + log[q(T* -> T0) / q(T0 -> T*)] + log I_T*(same) - log I_T0(same)
    + log Z_T0 - log Z_T*. Death is the reverse.
  - Against today's value this adds + log d* - log d0 + log Z_T0 - log Z_T*.
  - The frozen phi_k cancel because a birth or death never changes whether a frozen leaf has a constrained
    neighbour: 0 of 101 births in a 3x2 enumeration.
  - The proposal terms in [`birthOrDeathMove`](../../src/bartcore/moves.hpp) are unchanged.
  - An empty merged cone keeps its -HUGE_VAL sentinel, which is correct: there pi(T*, same) = 0.
- Constrained pairs. [`monotoneNeighborBounds`](../../src/bartcore/model.hpp) relates j < k when, along a
  constrained axis, j's code box ends one code below k's start (direction -1 flips it), and the boxes share a
  code on every other path axis.
  - Edge or corner contact is not a pair; its order follows by transitivity through the leaves between.
  - The order is acyclic: at a split on axis v the only cross pairs run along v, one way. A split on a free
    axis creates no cross pairs.
  - Related leaves share the sd c scale / k, so Z_T = e(P_T) / L!, with e the number of linear extensions.
  - The brute-force oracle agrees (zcheck mode below): 32 trees, 1-3 constrained axes, mixed directions, worst
    |z| 1.64.
- No closed form. A tree splitting on one constrained axis only is a chain (e = 1), but one constrained plus
  one free axis already yields N-shaped orders (e = 5 of 4!). The exact count is a DP over down-sets,
  O(down-sets x L).
- A merged leaf's relations are the union of its children's, so a death maps the new tree's down-sets
  injectively into the old tree's, and the whole-tree down-set count never grows under a death (0 of 6,575
  deaths across three enumerations).
- Measured sizes (current sampler, n 1000, p 5, 75 trees, 1/2/3 constrained predictors, 7,500 trees each):
  at most 9 leaves and 43 down-sets. The C++ count averages 0.14 us, against ~120 us per monotone tree step
  today (2.4 ms per 20-tree sweep, quadrature-bound).
- Stress (random guillotine partitions): two constrained axes, 116 us at 40 leaves. One constrained plus one
  free axis, 52 ms at 32 leaves and 0.1 s at 40 (850k down-sets).
- Redraw defect (blocking). [`MonotoneConstantGaussianLeaf::redrawAfterBirth`](../../src/bartcore/model.hpp)
  gives up after 100 tries and keeps mu[upper]. That slot is stale: freed node slots are recycled, and the
  resize zero-fills only new ones. It then draws the lower leaf on [aL, min(bL, stale)], which may be empty,
  and the clamp in [`MonotoneConstantGaussianLeaf::drawTruncatedNormal`](../../src/bartcore/model.hpp) then
  runs with lo > hi (undefined behaviour), leaving an infeasible state. The critique measured the fallback on
  about 0.5% of accepted births on flat data, 41% under a mildly decreasing partial residual, and 97% on a
  steeper one.
- Empty leaves. [`dbartsSampler$setPredictor`](../../man/dbartsSampler-class.Rd) with forceUpdate,
  updatePredictor, [`dbartsSampler$setData`](../../man/dbartsSampler-class.Rd) and
  [`dbartsSampler$setCutPoints`](../../man/dbartsSampler-class.Rd) can strand a member-empty leaf.
  [`MonotoneConstantGaussianLeaf::drawOneLeaf`](../../src/bartcore/model.hpp) pins it at 0 with no bound
  check, so the cone breaks and Z_T is no longer e / L!. This recurs in embedded use.
- Other paths. At creation R forces birth/death ([`resolveSamplerSpec`](../../R/spec.R)).
  [`xbart()`](../../R/xbart.R) refuses `monotone`, and [`ruleGibbsMove`](../../src/bartcore/moves.hpp)
  compiles out. Everything below was probed.
  - [`dbartsSampler$setControl`](../../man/dbartsSampler-class.Rd) installs a change/swap mixture.
    [`dbartsSampler$setModel`](../../man/dbartsSampler-class.Rd) with a model that lacks the monotone
    attribute does the same.
  - [`dbartsSampler$installTrees`](../../man/dbartsSampler-class.Rd) and
    [`dbartsSampler$setState`](../../man/dbartsSampler-class.Rd) accept an unconstrained donor's leaves. The
    fit then decreases by up to 1.4 along the constrained axis, and still by 0.48 after one sweep.
  - [`Chain::growForestFromRoot`](../../src/bartcore/chain.hpp) draws against the old tree's mu at reused
    node ids. It can also grow, as can installTrees and setState install, a tree of any size.
- Prior draw. [`MonotoneConstantGaussianLeaf::drawFromPriorForTree`](../../src/bartcore/model.hpp) samples the
  right law by rejection, but its 1e6-attempt cap fails ~6% of the time on a 9-leaf chain (1/9! = 2.8e-6).
- mBART (Chipman, George, McCulloch and Shively, arXiv 1612.01619v3). Section 3 defines only the "leaf"
  prior: eq. 3.3 incorporates the constraint "by constraining the CGM10 BART independence form p(Mj | Tj) ...
  to have support only over C", and the tree prior of section 3.1 "is the same form used for unconstrained
  BART" (section 4, after eq. 4.9). Section 4.3 then sets the move's normalizing constants d~* and d~0
  (eqs. 4.13 and 4.19) to one ("we reduce the computational burden"), which samples the "joint" prior, and
  compensates with base .25 and power .8 in place of .95 and 2: "we get tree sizes comparable to those obtained
  in unconstrained BART". All its examples use those values. The software (remcc/mBART_shlib, bd.cpp,
  coninteg1 and coninteg2) matches section 4.3: it accumulates the prior mass sumpr and never uses it. mBART
  users have therefore fit the "joint" prior with that retuned tree prior, on a grid.
- Unaffected: the leaf Gibbs sweep and the level-fibre shift (conditionals given T), and prediction.
  Monotone is new in 1.0-0, so NEWS gets no new entry; its existing 1.0-0 item, which shows
  monotone = c(x1 = "+", x2 = "-"), is edited to the new vocabulary (step 14).
- Claims to reword: [4. Decision - marginal likelihood for the structure moves](../design/monotone.md#4-decision---marginal-likelihood-for-the-structure-moves),
  [11. Costs, risks, and confidence](../design/monotone.md#11-costs-risks-and-confidence),
  [Plan-vs-code note](../design/monotone.md#plan-vs-code-note), and dec-B16 in [decisions.md](../decisions.md).

## Decision

This section records the options considered and the maintainer's rulings. The question was which prior, and
under the normalized ("leaf") prior, whether the exact count runs under a limit, and if so, one that changes the
model or one that stops the run.

A limit that gives trees zero prior mass must be closed under deaths. Then every allowed tree reaches the root
through allowed trees, the restricted chain is irreducible, and it targets the stated restricted prior. A limit
that is not closed can trap trees away from the root, and the chain's target then depends on where it starts.
A limit that stops the run changes no model and needs no closure.

Measured with [monotone-order-size.R](../../benchmarks/R/monotone-order-size.R) (fits, closure checks) and
[monotone_count.cpp](../../benchmarks/kernels/monotone_count.cpp) (the layered count in C++): the current
engine, n 5000 (one fit 20000), 1-2 constrained and 1-2 free axes, 1-3 seeds, 200 kept sweeps per fit. The C++
count runs at ~6 ns per down-set x leaf on these fits' components, and a move's components are within 10% of
the state metric the table uses. Sweep time is per kept sweep.

| trees | fits | sweeps where a 2^22 whole-tree budget binds | largest component: state / created by a move | sweep time with the count: per-fit medians (worst sweep) |
|---|---|---|---|---|
| 1 | 7 | all, in 1 fit (60 leaves, whole-tree up to 2.2e9) | 1.1e6 / 1.9e6 down-sets | 1.4x-260x, median 6x (710x) |
| 5 | 5 | all, in 4 fits; 4% in the fifth | 6.1e5 / 3.6e8 | 1.06x-14x, median 1.5x (41x) |
| 10 | 3 | none (up to 1.8e6) | 1.7e4 / 2.4e4 | 1.01x-1.22x (1.7x) |
| 20 | 2 | none (up to 1.5e4) | 54 / 86 | about 1.001x |

- Births create components up to 1.9x the largest state component. A death merges components: in the 5-tree,
  2-constrained + 1-free fit, 10% of trees have a death whose merged component passes 2^24 (up to 3.6e8), and
  such a death is proposed about once in 30 sweeps. Counting: algorithm counts every move on the finer tree's
  side, so the proposal does not count the merged component; if the death is accepted, the next move touching
  that component does. On the current engine's states, the moves' counted orders range up to 1.9e6 down-sets
  in these fits and 3.19e6 in a 1-tree fit with 3 free axes (Counting: algorithm).
- At ~2^24 down-sets the count is slower per unit: a 25-leaf star takes 15 s (36 ns per unit), a 60-leaf order
  of eight chains 17 s and 424 MB keeping two layers. Keeping every layer, as step 6's draw and the position
  laws do, took 23 s and 1.09 GB peak on a 1.68e7-down-set star plus a singleton.
- The current engine targets the d-normalized law, so the corrected one may grow different trees: on monotone
  data its 1/Z_T factor favours more constrained orders. That moves how much a budget truncates, and under
  option 5 only the cost.

1. Whole-tree down-set budget, B = 2^22 (the first recommendation). Closed: a death maps the new tree's
   down-sets injectively into the old tree's. Cost at most ~0.5 s and ~200 MB per count. What a fit sees:
   every sweep of 4 of 5 five-tree fits and of 1 of 7 one-tree fits has a tree past B, so births are refused
   and the posterior is not the documented model's, with no signal to the user. Its premise held only at 75
   trees.
2. Whole-tree budget at a larger B. Closed. Clearing the measured fits needs B > 2.4e9 (~2^31), and more data
   or more free axes go past any fixed B, since the product grows with the number of components. At such a B
   the product no longer bounds the work: one component near B costs minutes and ~100 GB per count. Slow counts
   still need option 5's handling, and then its truncation buys nothing.
3. Per-component budget: the largest or the touched component's count, or its width. It tracks the work but
   is not closed: a free-axis death merges two components, and the merged count can approach their product.
   4 of 1,164 deaths of random trees raise the largest component's count. The tree that splits x1
   (constrained) and then x2 (free) on both sides has two 3-down-set components, and both of its deaths make
   one of 5, so under a budget of 4 it cannot reach the root. The smallest closed bound is the largest
   component count over every pruning of the tree; no cheap way to compute it below the whole-tree count is
   known. Whole-tree width is closed (the same injection maps antichains), but it sums over components as the
   count multiplies, so a cap on it binds on many-component trees whose counting is cheap, as option 1 does.
   Not recommended.
4. Leaf-count cap. Closed. It does not bound the work (one component of L leaves can have 2^(L-1) + 1
   down-sets), so slow counts still need option 5's handling, and a cap low enough to matter truncates 1-tree
   fits (60 leaves seen). Not recommended.
5. No budget (taken for "leaf"). The documented prior holds verbatim for every fit. Two ways to handle a slow
   count were weighed.
   - A work guard (proposed by the agents, not taken): a dbartsControl limit, default 2^24 down-sets in one
     component, past which the run stops with an error, keeping the tree it held so a caller can raise the
     limit and continue. It changes no model, but exact then holds only for runs that complete, and rerunning
     with new seeds until one completes selects smaller trees.
   - Warn and continue (taken, dec-B149): no limit; the run finishes on the exact model however slow a count
     is, and R warns once after the run (step 15). Base-R-style fitters mostly cap up front or warn: glm's
     maxit warns, rstan's max_treedepth caps and warns, rpart's maxdepth and ranger's max.depth cap the tree.
     Here a cap is options 1-4; the ruling warns as glm does, but finishes the count.
   - What a fit sees: at 20 or more trees nothing, at 10 up to 1.2x slower sweeps, at 1-5 trees the table's
     slowdowns.
   - A cheap exact shortcut cuts counts (steps 1 and 2). Series-parallel and twin reductions were measured
     (Counting: algorithm): the first barely splits these orders, the second cuts down-sets 2-10x; neither is
     planned.
6. The unnormalized prior (dec-B144's rejected alternative, reopened by the maintainer; taken as "joint"):
   p(T, M) proportional to p_CGM(T) prod phi 1{M in C(T)}, BART's prior conditioned on every tree being
   monotone.
   - The move is exact with no count: today's score with d dropped. No budget and no slow counts. Under
     "joint" the move seam of step 2 and step 15's tally are inactive and sampleTreesFromPrior draws by
     count-free joint rejection (acceptance is the prior mean of Z_T); step 1's counter runs only for step 6's
     given-T prior draw (sampleNodeParametersFromPrior), off the MCMC path; steps 3-5, 7, 9-14 and 16 apply as
     to "leaf", and the gate retargets by dropping log Z_T from its weight.
   - What a fit sees: the tree marginal becomes p_CGM(T) Z_T, so each constrained split costs weight. Measured
     with monotone-order-size.R's logz mode: at 200 trees (n 5000, 1 constrained + 1 free), 20% of trees
     carry a constrained split, at 0.85 nats (~2.3x) each; at 20 trees 1.0 nats. In the exact one-tree
     enumerations ([monotone-exact-enumeration.R](../../benchmarks/R/monotone-exact-enumeration.R)
     unnormalized mode) constrained splits per tree fall from 1.47 to 1.15 (c1) and 0.87 to 0.55 (c3), and the
     posteriors differ by total variation 0.11-0.28. The effect on fitted functions at 200 trees is
     unmeasured.
   - The normalized prior instead keeps CGM's tree marginal: the constraint changes leaf values, not which
     trees are likely.
   - mBART's examples and software sample this prior, with base .25 and power .8 (Context).

The agents recommended the normalized prior alone with option 5: no budget, the guard as a control argument,
and the free-bound shortcut, since it costs nothing at the default tree count and keeps the documented model
for every fit that completes; they stated the unnormalized prior's case as no count, no guard, what mBART
samples, and a plain reading as BART conditioned on monotone trees, at the price of fewer splits on the
constrained predictors.

Rulings, 2026-09-29:

- Both priors (dec-B145). Asked "Why not both?", then: "Yes, proceed on both. I can already tell you for 1.
  that I want to know how users will 'feel' the default, beyond just the run time." The normalized prior
  takes option 5's no budget and the free-bound shortcut; the unnormalized prior needs no count. Open: the
  default, for the feel study below.
- Placement (dec-B146): a monotone() constructor in the existing monotone = argument, as interactions = and
  blocks = take interactions() and blocks(), carrying the per-predictor directions and prior = . Its
  signature was later ruled below (dec-B150). The plain
  vector monotone = c(x1 = "increasing") stays as shorthand for the default prior. cgm() and normal() do not
  change. Not chosen: a separate formal monotone.prior on dbarts(), bart() and dbartsSpec() (the LightGBM
  monotone_constraints_method precedent), an argument of cgm() or normal(), and a dbartsControl option. "Use
  your recommendation for the vocab, and I guess for the monotone() constructor too."
- Vocabulary (dec-B147): "increasing" and "decreasing", and 1, -1 and 0, with 0 unconstrained in the unnamed
  full-length positional form; "+" and "-" are dropped, and matching is case-sensitive, as base R's match.arg
  is. Not chosen: keeping the current vocabulary, which adds "+" and "-" and matches case-insensitively.
  Precedent: gbm, xgboost and LightGBM take 1, -1 and 0; mboost's bmono takes words. Monotone is new in 1.0-0,
  so there is no NEWS entry or deprecation.
- Values (dec-B148): prior = "leaf" (the normalized prior: the constraint restricts the leaf-value prior's
  support, the tree prior unchanged) or "joint" (the unnormalized prior: the tree's structure and leaf values
  conditioned together). "Go with \"leaf\" and \"joint\"." Not chosen: "per.tree"/"joint", "leaf"/"tree",
  "normalized"/"conditional". Open: their order, which is the default.
- Slow counts (dec-B149): no count limit under "leaf". The run always continues on the exact model; the engine
  records counts that took more than about a second, and after the run R warns once that some counts were
  slow, why, and the remedies (more trees, or prior = "joint"). A long count polls for a user interrupt, and
  an allocation failure in it becomes an ordinary R error, not a crash. Not chosen: the guard above, the hybrid
  Barker move ([monotone-barker-hybrid.md](../design/monotone-barker-hybrid.md), proposed, not adopted; the
  upgrade if slow counts show in practice), and a budget. The maintainer chose the no-limit option with a
  post-run warning, which the presentation numbered option 2 and which is this section's option 5, warn and
  continue: "OK, let's do option 2 now." What a cancel leaves behind is an orchestrator call (dec-A128, step 15).
- Signature (dec-B150): monotone(directions, prior = ), the directions a first-argument vector exactly as the
  plain-vector shorthand takes them, following blocks(groups, trees.per.group = NULL):
  monotone(c(x1 = "increasing", x2 = "decreasing"), prior = "leaf"). Not chosen: directions as named
  arguments with prior = reserved, which collides with a predictor named "prior", and accepting both forms.
  "Use option 1."

## Default: feel study

The maintainer wants the default chosen on how users will feel it, beyond run time. This study runs on the real
engine once both priors are built, in about an hour, and is not run before then. Its script lands beside the
gate in benchmarks/R.

- Arms: "leaf" and "joint" at cgm() defaults (base .95, power 2), "joint" at cgm(power = 0.8, base = 0.25), the
  mBART paper's setting, and an unconstrained fit as the reference a user compares against. Every monotone arm
  names its prior explicitly. Tree counts 200, 50 and 5.
- Prior predictive: 500 draws of f per arm and tree count, one constrained and one free predictor on a grid,
  through samplePriorPredictive. Along the constrained axis at fixed free values: the number of distinct
  levels (steps), the largest jump as a share of the total rise, the share of flat grid intervals, and the
  total rise; and the share of splits on the constrained predictor. This is what a user sees drawing from
  the prior.
- Fits against truth: x1 constrained increasing and x2 free on [0, 1]^2 plus one noise predictor, noise sd a
  third of the truth's sd. Truths: a smooth ramp (x1); a steep step (1{x1 > 0.5}); flat then rising
  (max(0, x1 - 0.6) / 0.4); and a strong free-axis interaction (x1 (1 + 2 x2) + sin(2 pi x2)), monotone in x1 for
  every x2. n 200 and 2000, 4 replicates, 500 burn-in and 500 kept sweeps, run in parallel over cores. Per
  fit:
  - RMSE against the truth on a test grid, and 95% interval coverage and mean width there;
  - the partial dependence along x1 (averaged over x2): its error, and its number of visible steps;
  - varcount shares on constrained, free and noise predictors, the variable importance a user reads;
  - held-out log predictive score on fresh data;
  - time per sweep; for "leaf" the largest count time per fit and whether the slow-count warning fired.
- Sensitivity to n.trees is read across 200, 50 and 5 in every measure above.
- What favours "leaf" as default: RMSE, coverage and held-out score at least as good as "joint" at either
  tree prior on the steep-step and interaction truths; varcount shares close to an unconstrained fit's, where
  "joint" under-reports the constrained predictor; prior draws along the constrained axis that look like
  unconstrained BART's with sorted leaves; sweep time within 10% of "joint" at 50 and 200 trees, and no
  slow-count warning there.
- What favours "joint": "leaf" over-splitting the constrained axis (prior draws with many small steps, wider or
  under-covering intervals, worse held-out score) or its counts turning slow at 50 trees or more; "joint"
  matching "leaf" on fit quality. If only the mBART-tuned "joint" matches, a "joint" default also needs a
  different tree-prior default under monotone, which users would meet as a surprise; the study reports that
  case separately.
- Slow counts in practice under "leaf", at whatever default, make the hybrid Barker move
  ([monotone-barker-hybrid.md](../design/monotone-barker-hybrid.md)) the upgrade; the largest count time per
  fit is the evidence.

## Counting: algorithm

Every move's ratio is counted on the side of the finer tree T*, the tree that holds the move's pair as two
leaves c1 and c2: for a birth the proposal, for a death the current tree. c1 is the child lower in the ORDER,
not in code: on an increasing axis the lower-code child, on a decreasing axis the higher-code child (taking the
lower-code one there gives theta = 0); on a free split either child, fixed as the left. A one-component death
counts its merged component C0, which costs no more than C*; a two-component death does not count the merged
component, which can be far larger, but if accepted it leaves that component in the state uncounted, and the
next move touching it pays its count (step 1).

- Identity. The linear extensions of T0 correspond one to one with those of T* in which c2 immediately follows
  c1: replace the merged leaf by c1 c2, or merge them back. Both directions hold because the merged leaf's
  relations are the union of its children's (Context). Let U be the union of T*'s components holding c1 or c2,
  m = |U| (the merged leaf's component size plus one), and theta the probability that c2 immediately follows c1
  in a uniform linear extension of U. Then Z_T0 / Z_T* = m theta, and theta <= 1 gives step 2's free bounds.
  A blind check confirmed the identity on 4,219 random guillotine births (1-3 axes, increasing, decreasing and
  free; 316 with the pair in two components), and theta against permutation enumeration on 2,861 moves.
- The pair in one component C* of T*: theta = e(C0) / e(C*), with C0 the order of C* with c1 and c2 merged
  (the union of their relations, closed). The layered DP counts both. C0's down-sets map injectively into C*'s,
  so C0 costs no more than C*, and when the current tree's count is kept, a move counts one of the two.
- The pair in two components C1 (holding c1, size a) and C2 (holding c2, size b), which a free split can
  produce: with P1(i) the probability that c1 is i-th in a uniform extension of C1, and P2(j) that c2 is j-th in
  C2, theta = sum over i and j of P1(i) P2(j) C(i+j-2, i-1) C(a+b-i-j, a-i) / C(a+b, a), the share of
  interleavings that put c1 just before c2 (C(n, k) binomials, held as logs). P1(i) sums
  f(D) g(D + c1) / e(C1) over the down-sets D of size i-1 that c1 can extend, where f(D) counts the orderings of
  D and g(D) the extensions of the rest: both come from the layered DP keeping every layer, as step 6 already
  does. C0, whose down-sets can approach the product of C1's and C2's (it exceeded both in 104 of the 316
  checked two-component births, by up to 15.5x), is not counted by this move.
- Scale: e(C) passes the double range (1.8e308) near 1,000 leaves, for example two ~515-leaf chains on a shared
  minimum, only ~2.7e5 down-sets; C(a + b, a) does so near a + b = 1,030. Each DP layer is held scaled by its
  largest value with the log scale carried beside it, and binomials as logs. Every term is positive, so nothing
  cancels.
- A move's own work is bounded by T*'s components, symmetrically in the move pair, since both directions count
  T*, plus any uncounted component of the current tree it touches (step 1). No state is guaranteed cheap.
- Measured with the pairs mode of [monotone-order-size.R](../../benchmarks/R/monotone-order-size.R) and
  [monotone_ratio.cpp](../../benchmarks/kernels/monotone_ratio.cpp): eight of the Decision's fits (1 tree: four,
  5 trees: three, 10 trees: one), all deaths and 10 random births per tree on every fifth kept sweep, 13,383
  moves, 1,395 of them with the pair in two components.
  - Exact: on the 13,361 moves whose C0 the direct count reaches, theta matches e(C0) / e(U) to relative
    2.3e-14.
  - The death that merges two components into 3.6e8 down-sets (5 trees, 2 constrained + 1 free) takes 34 ms,
    from lattices of 4.6e3 and 1.0e5 down-sets. Counting C0 directly takes 1.9e10 down-set x leaf units: minutes
    at the Decision's 6-36 ns per unit, and about 9 GB at its scaling.
  - On these states the counted orders range up to 1.9e6 down-sets (a birth in the 1-tree, 60-leaf fit): 0.8 s
    with both sides counted, and about half with the current tree's count kept. In the 5-tree, 2-constrained
    fit up to 1.2e6 and 0.5-0.6 s.
  - The same mode on a 1-tree fit with 1 constrained and 3 free axes (p 4, seed 1, beyond the Decision's 1-2
    free axes): up to 3.19e6 down-sets, 1.3-3.4 s per move with both sides counted, and 209 of 505 moves past
    2^20 down-sets.
  - Accepted merges: in the 5-tree, 2-constrained fit, tree 2 proposes deaths whose merged component reaches up
    to 3.6e8 down-sets on 24 of 40 sampled sweeps. These figures are proposals on the current engine's states;
    whether the corrected move accepts such deaths, and so what later moves pay, is unmeasured (checkpoint).
- What does not help, measured on the same orders:
  - Series-parallel decomposition. The merged 3.6e8 component is prime (53 of 53 elements); in the 5-tree,
    2-constrained fit it trims 118 of 662 single components of T* by at most 5 elements.
  - Twins, incomparable elements with the same relations (the antichain case of modular decomposition, where
    e = k! e(with the twins chained)): 2-10x fewer down-sets (3.6e8 to 9.5e7, 6.1e5 to 6.1e4). A later
    optimization at most.
  - Width reaches 22 in the merged component, so the treewidth of its incomparability graph is at least 21. This
    defeats both the DP's O(n^w) bound and the algorithm that is fixed-parameter tractable in that treewidth
    (Eiben, Ganian, Kangas and Ordyniak, ESA 2016, Theorem 17). A DP over the guillotine tree must carry how each
    cut face's leaves interleave, which is the down-set lattice again.
- Literature, checked in the papers named: counting linear extensions is #P-complete (Brightwell and Winkler 1991, as
  Dittmer and Pak and Eiben et al. state it); it stays #P-complete for height two, for dimension two, and for
  incidence posets of graphs (Dittmer and Pak, Electron. J. Combin. 27(4), 2020); and it has no algorithm
  fixed-parameter tractable in the treewidth of the cover graph unless FPT = W[1] (Eiben et al., Theorem 7). It
  is polynomial for series-parallel orders and for orders whose cover graph is a polytree (as Eiben et al. cite).
  Whether guillotine leaf orders form a tractable class is not known. The ratio is no easier than the count: a
  chain of births from the root multiplies ratios into e(T).
- Count-free alternative, not planned now. The identity makes theta a coin: draw a uniform linear extension of U
  exactly, and call it heads when c2 immediately follows c1.
  - Barker's acceptance R / (1 + R) is then exact by the two-coin algorithm (Goncalves, Latuszynski and
    Roberts, Braz. J. Probab. Stat. 31, 2017), with r1 the rest of the ratio. A birth has R = c theta with
    c = r1 m: with probability c / (1 + c) draw the coin and accept on heads, repeating on tails; otherwise
    reject. A death has R = 1 / (c theta) with c = m / r1, so the coin sits on the reject side: with
    probability c / (1 + c) draw the coin and reject on heads, repeating on tails; otherwise accept. Neither
    counts. Each takes c / (1 + c theta) <= min(c, 1 / theta) draws on average, and Barker's acceptance is at
    least half of Metropolis-Hastings'.
  - Exact draws come from Huber's bounding-chain coupling from the past (Discrete Math. 306, 2006), expected
    O(n^3 log n) steps. For the 3.6e8 death it takes 4.4 ms per draw of T*'s two components (24 and 30
    elements), and theta is 0.050 there, so at most about 20 draws (90 ms). Its draws matched the exact theta
    over 692 moves (mean z 0.007, mean z^2 0.94).
  - Used only when a move's T* components (or an uncounted current component it touches) pass a size
    threshold, a symmetric rule so each move pair keeps one exact kernel, it would bound the time of slow counts
    with an exact move that mixes somewhat slower, and give step 6 a draw that keeps no layers. It is dec-B149's
    upgrade if slow counts show in practice; [monotone-barker-hybrid.md](../design/monotone-barker-hybrid.md)
    proposes it, and the checkpoint (Staging) is the evidence for whether it goes in before release. This plan
    does not adopt it.
  - Simpler exact schemes lose acceptance. Accepting a birth with min(1, r1 m) times the coin, and a death with
    min(1, 1 / (r1 m)) and no draw, matches Metropolis-Hastings when r1 m <= 1 but otherwise scales acceptance
    by theta (0.02-0.2 on the large moves). The exchange algorithm (Murray, Ghahramani and MacKay, UAI 2006),
    with an exact leaf draw from T*'s cone, likewise scales acceptance by the chance of an indicator. The
    pseudo-marginal route needs an unbiased estimate of 1 / Z_T. One exists: adding the elements in turn, each
    one's chance p_k of landing above its predecessors has a coin (one exact draw of the elements before it), and
    the product of the geometric counts of trials to heads is unbiased. Its relative variance is the product of
    (2 - p_k) less one, exponential in the component size, so the chain would stick.

## Constraints

- Exact for each stated prior (dec-B14). Counts are held as scaled doubles with a log scale (Counting:
  algorithm): exact to 2^53, relative error under 1e-12 beyond, and no overflow.
- Unconstrained samplers byte-identical: the new seam, the prior flag and the lazy cache compile out, like the
  three existing monotone seams. No dbarts.h change. Not compiled out, and so checked for identical draws: a
  new facade virtual for the slow-count tally (a --preclean rebuild), the try/catch in run()'s worker bodies,
  and the cancel function passed down from Chain::run to the moves.
- The default prior lives in one constant, read by monotone() and by the plain-vector shorthand; until the
  default is ruled its value is "leaf", provisionally, matching today's target. Every harness or test that
  encodes the "leaf" target names prior = "leaf" explicitly (monotone-reference.R part (a) and
  test-calibration-prior-draws.R's monotone block use the shorthand today), so a default ruling changes no
  gate.
- Out of scope: change moves under the constraint, quadrature speed (TODO monotone-leaf-quadrature), and
  reconciling a chi k hyperprior with the truncated law.

## Staging

Reviewed commits, in this order. Landing is two-phase: commits 1-5 land with the default provisional ("leaf",
Constraints) and every harness pinned, and step 12's docs wait for the default ruling.

1. Redraw fix, empty leaves and reachability (steps 4, 5, and of step 7: the proposal-mix refusal, setControl,
   growForestFromRoot's reseed, and reseed-then-validate with setState's up-front check), with those bullets'
   tests from step 9.
2. The monotone() constructor, vocabulary and caller sweep (step 14). prior = "joint" is refused with a message
   saying it arrives with the corrected move, until commit 4; the default constant is "leaf", provisionally.
3. The counter and the linear-extension draw with their tests/cpp checks, not wired in (steps 1, 6 and 8's
   counter checks), counting both sides of every move.
4. The seam, the dropped d terms, the prior flag and the gate under both priors (steps 2, 3, 13, 10), step 6's
   switch of `drawFromPriorForTree` to the exact draw, and step 7's sampleTreesFromPrior bullet for both
   priors; prior = "joint" is accepted from here. Per-move count timing is recorded through the move census
   (BARTCORE_MOVE_CENSUS; step 2). Then the checkpoint below, and stop.
4b. The lazy cache (step 1), after the checkpoint: the checkpoint sees the uncached worst case, since the cache
   changes only cost, and its fits are re-timed with the cache.
5. Slow-count warning, interrupt, allocation and rebuild (step 15), and step 7's setModel refusal (dec-A128).
   Commit 4 alone cannot be interrupted mid-count, and an allocation failure on a worker thread there calls
   std::terminate, so commits 4 and 5 are pushed together; before commit 5, commit 4 runs only on the
   implementer's machine, under the checkpoint's caps.
6. Equivalence scenarios (step 11), the SBC arms (step 16), and docs (step 12, after the default ruling).

Checkpoint (stop and report), after commit 4, before the feel study and before step 12:

- Fits: [monotone-order-size.R](../../benchmarks/R/monotone-order-size.R)'s fits on the corrected engine under
  prior = "leaf", built with the move census: 1, 5, 10 and 20 trees, 1-3 free axes, the p 4 case (1
  constrained, 3 free) included, 1000 burn-in and 200 kept sweeps.
- Caps (macOS has no `timeout`, and `ulimit -v` fails there): a driver script starts each fit as a child
  Rscript with system2(wait = FALSE), the child writing Sys.getpid() to a file first. The driver polls every
  second with `ps -o rss= -p <pid>`, and kills the child with tools::pskill when it passes 60 min of wall time
  (a 1-tree, 3-free-axis fit measured 1.3-3.4 s per count on 41% of moves, ~30 min for 1,200 sweeps) or 8 GB
  of resident memory (a direct count of the 3.6e8-down-set merge measured ~9 GB). A killed fit is recorded as
  capped, with the cap it hit.
- Record per fit, from the census: the largest per-move count time and memory, the share of sweeps with a
  count over 1 s, and how often accepted deaths leave a merged component past 2^20 down-sets; and the hybrid
  design's quantities, the share of moves it would switch at B = 2^22 and a_B / a_MH on those moves.
- Trigger: the hybrid Barker move ([monotone-barker-hybrid.md](../design/monotone-barker-hybrid.md)) goes to the
  maintainer for a before-release call if any fit at 10 or more trees has a count over 1 s (the measured
  10-tree fits stay under 2.4e4 down-sets, milliseconds, so this is a change in kind at the tree counts users
  run); or any 1-5 tree fit has a count over 1 s in more than 10% of its sweeps (at 0.1-1 ms per sweep
  otherwise, such a fit runs over 100x slower on those sweeps); or any count takes over 60 s or 4 GB (a stall
  or a memory risk on a laptop); or any fit hits a cap. Otherwise the checkpoint's numbers go to the
  maintainer as a report, and the plan continues. Either way the implementer stops here and reports.

## Steps

1. Order counter (engine, beside the monotone geometry), used by the "leaf" prior's move and by both priors'
   given-T leaf draw (step 6). Down-set keys are multi-word bitsets.
   - Relation: build the tree's relation matrix with the adjacency test of `monotoneNeighborBounds`, and split
     it into components.
   - Count: run the layered down-set DP per component, keeping two layers (every layer for step 6's draw and
     for position laws), and return the whole-tree log Z_T = sum over components of log e(C) - lgamma(|C| + 1).
     There is no limit (dec-B149); step 15 times each count, lets it be interrupted, and handles its
     allocation failure.
   - Ratio (Counting: algorithm): a move changes only the component of the merged leaf in T0 and the
     component(s) holding the two children in T*, and is counted on T*'s side. With the children in one
     component C*, count C* and C0; the current tree's component counts are kept, so a move counts one of the
     two. With the children in two components, take theta from the position laws of c1 and c2, and do not count
     C0. c1 is the child lower in the order (Counting: algorithm).
   - First version (commit 3): count both sides of every move; nothing is cached.
   - Lazy cache (commit 4b): per-component counts of the current trees, owned by the monotone leaf
     model's per-chain scratch, not by the shared Tree, and compiled out elsewhere. An accepted two-component
     death leaves its merged component uncounted, marked so; the next move that needs that count (a
     one-component birth or death inside it) counts it then, at the merged component's full cost, and keeps
     it. The cache is invalidated wholesale by setState, installForest, growForestFromRoot,
     sampleTreesFromPrior, setData, setPredictor, setCutPoints and copy. No claim is made that a state is cheap.
   - Scale: layers held scaled with a log scale, binomials as logs (Counting: algorithm).
2. Seam, active under "leaf" only (step 13). Add an optional leaf concept, `logNormalizerRatio`, declared by
   the monotone leaf. [`birthOrDeathMove`](../../src/bartcore/moves.hpp) evaluates it once per move while T*
   is in place, after `tree.birth` or before `orphanChildren`, given the pair, and it returns
   log Z_T0 - log Z_T* counted on T*'s side (Counting: algorithm). A birth adds it to the log prior ratio and a
   death subtracts it (logged in the census's prior column). Under BARTCORE_MOVE_CENSUS the census also logs
   each count's wall time, peak down-sets and component size, whether B = 2^22 would switch the move, and
   a_B / a_MH from the same ratio; that census build is how R reads per-move timing for the checkpoint. Step
   15's tally times counts in every build.
   - Free bounds: Z_T0 / Z_T* = m theta with theta <= 1 (Counting: algorithm), m the size of the union of T*'s
     components holding the pair, known without counting. So a birth's ratio is at most m and a death's,
     Z_T* / Z_T0, at least 1 / m. This is tighter than the whole-tree bounds L0 + 1 and 1 / L0 (a birth never
     lowers e; 0 of 1,164 random births do, monotone-order-size.R closure mode). Draw u first; with r1 the rest
     of the ratio, reject a birth without counting when u > r1 m, and accept a death without counting when
     u < r1 / m. This is plain Metropolis-Hastings with no loss of acceptance; the savings are
     unmeasured. (A two-stage acceptance, min(1, a) min(1, b), is also exact but lowers acceptance.)
   - A slow count never rejects or stops the move; an interrupt or allocation failure during it cancels the
     move and the run (step 15).
3. Drop the d terms: `priorMass` in `oneLeafLogMarginal`, `denom` in `twoLeafCoupledLogMarginal`. Under "joint"
   this alone is the exact move.
4. Exact pair redraw in `redrawAfterBirth`:
   - Keep the capped rejection for the upper leaf; it is exact whenever it accepts.
   - When the cap is reached or acceptMax underflows, draw the upper leaf by inverting the CDF of its marginal,
     phi_U(u) [Phi_L(min(bL, u)) - Phi_L(aL)] on [max(aL, aU), bU]. Use the `coneProbability` quadrature's
     cumulative with safeguarded Newton, evaluate the density exactly, and work on the log scale in the tail.
   - Then draw the lower leaf on [aL, min(bL, mu_U)], which is never empty.
   - Mixing a capped exact rejection with an exact fallback is exact. Never read mu[upper] before writing it,
     and never call the clamp with lo > hi.
5. Empty leaves are leaves with no data. `drawOneLeaf`, both redraws and the prior draw give an empty leaf its
   prior truncated to its neighbour bounds instead of pinning it at 0, and
   [`monotoneTreeIsFeasible`](../../src/bartcore/model.hpp) checks empty leaves too.
   - The cone and Z_T = e / L! then hold in every state. Trees with an empty leaf still carry zero posterior
     mass through the veto, and the veto ranks move the chain out.
   - The cost is that a test row routed to an empty leaf predicts that leaf's draw instead of 0.
   - An empty constrained sibling pair (both children of a birth empty) is drawn by step 4's coupled pair draw
     with no data, not leaf by leaf.
   - The alternative, refusing mutations that strand a leaf, breaks the embedded use the sampler exists for.
6. Exact prior leaf draw. Per component, draw a uniform linear extension by backward sampling on the DP
   counts, every layer kept: remove a maximal element x with probability e(D - x) / e(D). Draw |C| iid
   N(0, c scale / k), sort them, and assign them in extension order; isolated leaves draw alone. This replaces
   the rejection loop and its 1e6 cap in `drawFromPriorForTree` (switched in by commit 4).
   - Out of a run (sampleNodeParametersFromPrior, sampleTreesFromPrior): the count runs before any leaf of the
     tree is written, so an allocation failure leaves the tree's leaves as they were, and the bridge's
     existing conversion of C++ exceptions makes it an ordinary R error.
7. Reachability:
   Each bullet names its commit (Staging).
   - Commit 1. A bridge refusal after [`parseProposalProbs`](../../src/R_interface_bartcore.cpp), at creation and in
     [`bartcore_setModel`](../../src/R_interface_bartcore.cpp). It is keyed on the engine's leaf kind
     ([`LeafModelKind`](../../src/bartcore/model.hpp)), not on the incoming model, and refuses any nonzero swap,
     change, perturb or rule_gibbs share. The all-zero frozen mixture stays allowed.
   - Commit 1. setControl mirrors creation: a defaulted mixture is rewritten to birth/death silently, and an explicit
     non-birth/death one is refused.
   - Commit 5. setModel refuses a change of monotone directions or of the monotone prior against the sampler's
     own, as setControl mirrors creation, refuses a model lacking the monotone attribute on a monotone sampler,
     and refuses one adding the attribute to an unconstrained sampler (dec-A128). Today a model lacking it is
     accepted, and copy() then fits unconstrained: a 0.34 drop along x1 in the critique's probe.
   - Commit 1. growForestFromRoot reseeds mu to the all-zero feasible seed before its draw.
   - Commit 1. Reseed, then validate (dec-A128). [`Chain::installForest`](../../src/bartcore/chain.hpp) reaches
     the trees through [`Chain::rebuildLiveForest`](../../src/bartcore/chain.hpp) or
     `rebuildLiveForestRemapped`. On the installTrees path, every installed tree that fails
     `monotoneTreeIsFeasible` has all its leaves set to 0, the all-zero feasible seed growForestFromRoot uses
     (equal values satisfy every constraint, and no RNG is drawn); feasible trees keep the donor's values. The
     next sweep's leaf Gibbs step redraws them.
   - Commit 1. setState, copy and reload are refused up front: Sampler::setState checks a new
     monotoneStateFeasible predicate on every chain's live trees before any chain is touched, beside
     interactionStateFeasible and columnMaskStateFeasible, so a refusal leaves the sampler exactly as it was, as
     setState promises. With steps 4 and 5 every state the sampler produces passes it.
   - Commit 4. [`Chain::sampleTreesFromPrior`](../../src/bartcore/chain.hpp) under "leaf" needs no new predicate: the
     prior is not restricted. Under "joint" it draws each tree jointly by rejection, a CGM tree and iid
     unconstrained leaves kept only when the leaves lie in its cone (acceptance is the prior mean of Z_T), with
     no count, and then discards the leaves, as its contract returns trees without leaf values.
8. tests/cpp:
   - The count against brute-force permutation counts on hand-built trees (1-3 axes, mixed directions, N, a
     star, a chain over 64 leaves) and 200 random trees.
   - The ratio on T*'s side against e(C0) / e(U) from direct counts, over every birth and death of those random
     trees, with the pair in one component and in two, and the position law against brute force. Hand cases,
     x1 constrained and x2 free: x1 cut twice into a < b < c, then b cut on x2 puts the pair in one component
     (theta 1/2, ratio 2); x1 cut once, then both halves cut on x2 at the same code, puts a sibling pair in two
     components (theta 1/3, ratio 4/3). Both again with x1 decreasing, giving the same thetas; and a single x1
     cut with x1 decreasing, a constrained pair whose c1 is the higher-code child: theta 1 (ratio 2), where
     the lower-code child would give 0.
   - Lazy counts: an accepted two-component death marks the merged component uncounted, and the next move
     inside it counts it and matches a direct count.
   - Scale: two 515-leaf chains on a shared minimum give log e and a two-component theta equal to their closed
     forms, with no overflow.
   - Log Z against the enumeration's e / L!.
   - RNG-free ratio tests on A < B -> A < B1 < B2 (constrained split) and A < B -> {A < B1, B2} (free split)
     with pinned mu_A. Compare the move's log alpha with the closed form of the corrected statement in Context
     (it must differ from the current code's value).
   - The pair redraw on a design with P_post(lower <= upper) ~ 1e-4, where the old loop exhausts almost every
     time: every draw is feasible, and the upper leaf's draws match the quadrature CDF (KS).
   - The linear-extension draw is uniform over extensions (chi-square on a 5-leaf N-plus-chain).
   - Slow counts (step 15): with the threshold lowered, a count enters the tally with its time and size, and
     the tally resets at the next run. A cancel function that turns true after k polls aborts a long count,
     inline and on a worker chain; an allocation failure injected into the count reaches the caller as an
     exception, inline and rethrown from a worker after the join. After either, the touched tree is T0, the
     chain's derived state equals a from-scratch rebuild of the same trees, and the next run is valid.
   - Free bounds: over random trees and moves, a move decided without counting gets the same decision as with
     the count, for the same u.
   - Under "joint" the move's log alpha equals the closed form without log Z_T0 - log Z_T*, the count is never
     called, and the joint prior draw's tree marginal matches p_CGM(T) Z_T on small enumerated trees.
   - [`testMonotoneMarginal`](../../tests/cpp/test_model.cpp) loses its d_* = 1/2 check.
9. tinytest ([test-monotone.R](../../inst/tinytest/test-monotone.R)):
   - a setControl change mix errors, and a defaulted one is rewritten;
   - setModel refuses a model without the monotone attribute, or with other directions or another prior, and
     refuses adding a constraint to an unconstrained sampler;
   - installTrees from an unconstrained donor zeroes its infeasible trees' leaves and leaves the fit monotone at
     once; setState of the unconstrained donor's state is refused, and the sampler is unchanged (its state and
     predictions equal those before the call);
   - growFromRoot plus one sweep is monotone;
   - after 2,000 sweeps on data decreasing along the constrained axis, copy() and setState round-trip;
   - setPredictor(forceUpdate = TRUE) stranding a leaf keeps the fit monotone;
   - with the slow-count threshold lowered through its test hook, a run warns once (dbartsSlowCountWarning)
     naming the remedies, and a run with the default threshold does not;
   - an allocation failure injected through a test hook is an ordinary R error, and the sampler runs on
     afterwards with a valid state; an interrupt is covered in tests/cpp only;
   - every check above runs under both priors, and a "joint" fit's moves never call the count (the given-T
     prior draw of step 6 counts under both);
   - step 14's constructor and vocabulary checks.
   A statistical check does not fit: the most sensitive cheap functional sat at |z| 0.6-0.8 against the current
   move at 50k-100k draws.
10. Gate: [`buildDesign`](../../benchmarks/R/monotone-exact-enumeration.R), [`prototypeKeys`](../../benchmarks/R/monotone-exact-enumeration.R) (details under Verification).
    - The script enumerates every one-tree structure on cell grids with the CGM prior, weights each by
      exp(sum base) P_post(C) / Z_T, and computes P_post(C) by the same down-set recursion carrying
      Richardson-extrapolated trapezoid quadrature, self-checked to 2e-12 against the closed form.
    - It compares draws per root rule, which only the root-only tree can change, against each group's
      conditional law, using a Hotelling T^2 on batch-mean cell frequencies with its F reference. A group fails
      at p < 1e-4.
    - Designs, 10 rows per cell: c1, x1 constrained, 4 cells; c2, x1 and x2 constrained, 3x2; c3, x1
      constrained and x2 free, 3x2; cN, x1 constrained and x2 free, 2x3, where the two halves change level at
      different cuts, giving 26% N mass in the tested x1-rooted group.
    - The script replaces the planned part (c) of monotone-reference.R. It joins exact-gates.yaml's list in the
      fix commit, since it fails the current engine by design.
    - The general DP beyond these sizes rests on step 8's brute-force checks.
    - Both targets: engine and prototype modes take the prior, weighting each tree with log Z_T in its weight
      ("leaf") or without it ("joint"), and every group must pass under each. Its unnormalized mode already
      prints the two exact laws side by side. Every run names its prior.
    - A mirrored design with x1 decreasing (cN's cells reflected along x1), so a c1 taken by code instead of
      by order fails the gate.
    - Runtime: quick mode measured 8 min 15 s per prior over four designs, so with the mirrored design about
      10.5 min per prior and 21 min for both. The monotone gate runs as its own CI job with a 40-min timeout,
      and exact-gates.yaml keeps its job and timeout.
11. Two monotone scenarios in benchmarks/R/equivalence.R (x1 and x2 constrained, 20 trees, one per prior): a
    55-scenario re-record, the other 53 bitwise, and MANIFEST rows naming the enumeration gate as their ORACLE
    (P17).
12. Docs:
    - monotone.md section 2 drops "+" and "-" for the new vocabulary (dec-B147), and sections 4, 9 and 11 and
      the Plan-vs-code note state both priors, B' with the whole-tree
      normalizer under "leaf", and what mBART's paper and software sample (Context).
    - The help states both priors and names the chosen default with the feel study's reason; under "leaf" it
      states the few-tree cost and the slow-count warning with its remedies, and that a count's memory grows
      with it (1.09 GB measured for 1.68e7 down-sets with every layer kept).
    - [Monotone arm: design](sbc-family-tiers.md#monotone-arm-design) states one arm per prior (step 16).
    - dec-B16 is marked superseded in part by dec-B144.
    - Status lines and INDEX at landing, and the TODO item removed. These docs land in phase two, after the
      default ruling (Staging).
13. The prior switch, end to end. R resolves monotone(prior = ) (step 14); the bridge passes it with the
    directions to the engine as a flag on the monotone leaf. The flag rides the model's monotone attribute with
    the directions, not the saved state, so stateFormatVersion does not change; copy and reload rebuild from the
    model, and setModel refuses a change (step 7). Under "leaf" the Z seam (step 2)
    and step 15 are active; under "joint" `logNormalizerRatio` is not evaluated, no
    count runs, and sampleTreesFromPrior draws jointly (step 7). Steps 3-6 apply to both. The fit object
    records the prior, and print and summary show it.
14. The monotone() constructor and vocabulary (dec-B146, dec-B147).
    - Signature (dec-B150): monotone(directions, prior = ), for example
      monotone(c(x1 = "increasing", x2 = "decreasing"), prior = "leaf"), following blocks(groups,
      trees.per.group = NULL). directions is exactly the plain-vector shorthand, named or positional, so a
      predictor named "prior" needs nothing special. prior = is one of "leaf" and "joint" (their order, the
      default, open), checked with match.arg, and the default comes from the one constant (Constraints). Class
      dbartsMonotone, in R/model.R beside interactions() and blocks().
    - Like interactions() it is not exported itself: it joins dbartsForests, its exported face, and resolves
      by bare name inside monotone = on dbarts(), bart() and dbartsSpec() through resolveForestArguments
      (FOREST_ARGUMENT_VOCABULARIES), with a bare name the caller bound to a value being that value.
      interactions() has no print or format method, so neither does monotone().
    - `resolveMonotone` takes a dbartsMonotone or the plain vector, which is shorthand for the default prior,
      and returns the directions and the prior.
    - `parseMonotoneSign` takes "increasing", "decreasing", 1, -1 and 0, matched case-sensitively; "+", "-" and
      other cases are errors naming the vocabulary. 0 means unconstrained in both forms: dec-B147 names it for
      the positional form, and c(x1 = 0) stays accepted as unconstrained in the named form (test-monotone.R
      relies on it), since 0 is the same value there.
    - bart.R's multinomial branch routes monotone through resolveForestArguments, as it does variance, so a
      bare monotone() resolves there too.
    - Rd: a monotone topic beside interactions and blocks, with a _pkgdown.yml entry; the monotone items of
      dbarts.Rd, bart.Rd (both drop "case-insensitive") and dbartsSpec.Rd, and dbartsForests.Rd, updated.
    - Caller sweep for the vocabulary (every "+", "-" or case-variant direction): test-monotone.R,
      test-blocks.R, test-argument-surface.R, test-proposal-probs.R, inst/NEWS.Rd's 1.0-0 item (edited, not a
      new entry), and the SBC design text in sbc-family-tiers.md. benchmarks/R/binary-hyperprior.R and
      surfaces-common.R carry no monotone call and are not touched.
    - tinytest: bare-name resolution and a caller-bound shadow, the vector shorthand equal to the default
      prior's monotone(), an unknown prior value, "+" and "Increasing" refused, 0 accepted in the positional
      form, a named 0 accepted as unconstrained, a predictor named prior constrained through directions, the
      default print of a monotone() object showing its directions and prior (no print or format method), and
      each prior fitting monotone.
15. Slow counts (dec-B149).
    - Tally: each chain times every count and records those over a threshold, default one second: how many,
      the slowest, and its component's size and down-sets. It resets at the start of each run, as the GP
      fallback tally does, and the threshold is settable only through a test hook.
    - Warning: the bridge attaches the tally to run's result as it attaches the GP fallback census, and R warns
      once after the run, in warnOnGPFallback's pattern (class c("dbartsSlowCountWarning", "dbartsWarning")):
      some counts took more than about a second, because under prior = "leaf" a large tree's leaf order is
      costly to count in time and memory, and the remedies are more trees or prior = "joint". bart() and the
      sampler's run both warn.
    - Interrupt: the DP polls the chain's cancel function, the one Chain::run checks between sweeps, every
      2^16 down-sets. On inline chains it calls the host's throttled pollInterrupt on the main thread; on
      worker chains it reads the atomic cancel flag the main thread sets, so no worker calls into R, and SIGINT
      stays blocked on workers as run() already arranges.
    - Mechanics. A cancel inside the count throws bartcore::CountCancelled (a new std::exception type); an
      allocation failure throws std::bad_alloc. birthOrDeathMove catches both, undoes the birth (the rollback
      its reject path already takes) or leaves the death unapplied (the count runs before `orphanChildren`),
      and rethrows. Chain::run catches both, rebuilds (below), and returns cancelled for CountCancelled or
      rethrows bad_alloc. In Sampler::run each worker's try/catch sits inside its per-chain body, around
      chains_[c]->run, so numChainsRunning is still decremented and the join cannot deadlock; it stores the first
      exception_ptr and sets the cancel flag, and run() rethrows it after the join on the caller's thread, as
      fanOutPredictSlabs does. The bridge turns it into an ordinary R error.
    - After a cancel (dec-A128): the chain restores the touched tree to T0 (undoing the birth, or before the
      death's orphaning), then rebuilds its derived state from the trees: the running residual, the total
      fits, and every cached fit the sweep had advanced. The cost is one fit pass over every tree, O(n x trees),
      as setState's rebuild. Every state after a cancel is then consistent and valid, but the sweep may be
      partly applied: earlier trees of that sweep hold their new draws. The run's results are discarded as on
      any interrupt, and the sampler can run again.
16. SBC arms (benchmarks/R/sbc.R has none today; ~120 lines), as
    [Monotone arm: design](sbc-family-tiers.md#monotone-arm-design) lays them out: the burn-monotone run, the
    monotone-1 and 20-tree arms once per prior, each naming its prior, and the unconstrained monotone-bd twin
    once. monotone-1 at 0.4 s per replicate, the 20-tree arm ~85 min at R 200.

## Verification

- `cd tests/cpp && make && ./test_bartcore`: the count, ratio, lazy-count, scale, redraw, extension-draw,
  slow-count, interrupt, allocation and "joint" checks pass.
- The checkpoint (Staging) is reported before commits 5 and 6.
- `R CMD INSTALL --preclean -l <lib> .`, then `R_LIBS=<lib> Rscript -e 'tinytest::test_package("dbarts")'`.
- `R_LIBS=<lib> Rscript benchmarks/R/monotone-exact-enumeration.R quick`, under each prior, the mirrored
  decreasing design included: every group passes. This takes about 8 min per prior; full mode runs 900k draws.
  The mutation run restores the d divisions and drops the Z term (then `touch` the header and reinstall), and
  must fail like the current engine. A second mutation takes c1 by code instead of by order, and must fail the
  mirrored decreasing design.
  - Current engine, quick: c1 p 1.6e-24 and 6.5e-35 (two root rules); c2 3.8e-10 and 4.3e-5; c3 1.4e-15;
    cN 1.9e-29.
  - `prototype` mode (the corrected move in R, with an exact pair redraw): every group p >= 0.07.
  - `prototype-old` reproduces the current engine (c1 T^2 196 and 246, against the engine's 163 and 253).
  - Against an e(N) miscount of 10%, cN's noncentrality is ~80, so power is about 1.
  - `zcheck` mode verifies Z_T = e / L!.
- `R_LIBS=<lib> Rscript benchmarks/R/monotone-reference.R quick`, with prior = "leaf" named: parts (a) and (b)
  still pass. Then run the
  whole exact-gates.yaml list with `quick`.
- `Rscript benchmarks/R/equivalence.R compare <current>`: the 53 unconstrained scenarios "identical draws (same
  RNG stream)" and no "max |z|", after the 55-scenario re-record. BCF and multinomial compare identical. The
  snapshot files carry no monotone fit.
- Release level: step 16's SBC arms, once per prior. monotone-1 (0.4 s per replicate) flags the current move and
  must pass; then the 20-tree arm (~85 min at R 200) must pass before admission.
- Speed: on a quiet machine, under each prior, monotone sweep time at 20 trees, 1 and 2 constrained predictors,
  within 5% of today; at 1 and 5 trees under "leaf", the slowdown is recorded against the Decision's estimates.
  bench-sampler compare unchanged on the unconstrained paths.
- The feel study (Default: feel study) runs once both priors pass the gates above, before the default is
  ruled.
- `Rscript tools/check-doc-freshness.R .` passes.
