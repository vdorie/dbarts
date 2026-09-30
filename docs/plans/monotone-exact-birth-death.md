# monotone-exact-birth-death: the monotone birth/death move targets the documented prior

Status: PLANNED 2026-09-29 (dec-B144), revised after a blind critique; the Decision was reworked after a
few-tree measurement and its critique, adds the unnormalized prior as an option, and remains the maintainer's
call. Derivation and gate verified on an R prototype; not
implemented.

agent: opus (engine numerics: move seam, order counting, exact pair redraw, gate)
rng: posterior-changing for every fit with an active monotone constraint (all of its draws move, prior draws
included); unconstrained fits byte-identical, since every engine change lives in the monotone instantiation
window: before 1.0-0 (TODO monotone-exact-birth-death)
budget: ~1,050 lines (engine ~450, bridge and R ~110, tests/cpp ~320, tinytest ~90, gate wiring ~10, docs ~70;
the gate script itself is already in the tree). Plan estimates have run 1.5-2x low: expect up to ~2,000.

## Goal

The monotone chain targets the documented prior exactly. That prior is the CGM tree prior and, given the tree,
iid (c-inflated) normal leaves restricted to the cone C(T), normalized per tree by
Z_T = P(unconstrained leaves lie in C(T)). This covers four pieces:

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
- Unaffected: the leaf Gibbs sweep and the level-fibre shift (conditionals given T), and prediction.
  Monotone is new in 1.0-0, so NEWS gets no entry.
- Claims to reword: [4. Decision - marginal likelihood for the structure moves](../design/monotone.md#4-decision---marginal-likelihood-for-the-structure-moves),
  [11. Costs, risks, and confidence](../design/monotone.md#11-costs-risks-and-confidence),
  [Plan-vs-code note](../design/monotone.md#plan-vs-code-note), and dec-B16 in [decisions.md](../decisions.md).

## Decision

Question: which prior, and under the documented (normalized) prior, should the exact count run under a limit,
and if so, one that changes the model or one that stops the run?

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
  such a death is proposed about once in 30 sweeps.
- At ~2^24 down-sets the count is slower per unit: a 25-leaf star takes 15 s (36 ns per unit), a 60-leaf order
  of eight chains 17 s and 424 MB keeping two layers. Keeping every layer, as step 6's draw does, costs about
  1 GB.
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
   the product no longer bounds the work: one component near B costs minutes and ~100 GB per count. It needs
   option 5's guard anyway, and then its truncation buys nothing.
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
   down-sets), so it needs option 5's guard too, and a cap low enough to matter truncates 1-tree fits (60
   leaves seen). Not recommended.
5. No budget, with a work guard. The documented prior holds verbatim for every fit.
   - The guard is a dbartsControl argument, default G = 2^24 down-sets in one component. A count past it stops
     the run with an error naming the order's size and suggesting more trees or a larger G. It never rejects
     a move, so it changes no model and needs no closure. The sampler keeps the tree it held before the move,
     so an embedded caller can catch the error, raise G and continue.
   - What a fit sees: at 20 or more trees nothing, at 10 up to 1.2x slower sweeps, at 1-5 trees the table's
     slowdowns. The margin over proposals is ~9x in the 1-tree fits, and the 5-tree, 2-constrained fit above
     would stop within its first few dozen sweeps at any G that can be counted.
   - Exact holds for runs that complete. Rerunning with new seeds until one completes selects smaller trees;
     the help page says so.
   - Base-R-style fitters mostly cap up front or warn instead: glm's maxit warns, rstan's max_treedepth caps
     and warns, rpart's maxdepth and ranger's max.depth cap the tree. Here a cap is options 1-4, and a warning
     cannot continue without the count.
   - A cheap exact shortcut cuts counts (steps 1 and 2). Splitting a component by series-parallel or modular
     decomposition could shrink the counts further; unmeasured, future work.
6. The unnormalized prior (dec-B144's rejected alternative, reopened as an option by the maintainer):
   p(T, M) proportional to p_CGM(T) prod phi 1{M in C(T)}, BART's prior conditioned on every tree being
   monotone.
   - The move is exact with no count: today's score with d dropped. No budget or guard; steps 1 and 2 and the
     guard pruning in step 7 go. sampleTreesFromPrior becomes a count-free joint rejection (acceptance is the
     prior mean of Z_T). Step 6's counter survives only for the given-T prior draw
     (sampleNodeParametersFromPrior), off the MCMC path. Steps 3-5 and 9-12 stay, and the gate retargets by
     dropping log Z_T from its weight.
   - What a fit sees: the tree marginal becomes p_CGM(T) Z_T, so each constrained split costs weight. Measured
     with monotone-order-size.R's logz mode: at 200 trees (n 5000, 1 constrained + 1 free), 20% of trees
     carry a constrained split, at 0.85 nats (~2.3x) each; at 20 trees 1.0 nats. In the exact one-tree
     enumerations ([monotone-exact-enumeration.R](../../benchmarks/R/monotone-exact-enumeration.R)
     unnormalized mode) constrained splits per tree fall from 1.47 to 1.15 (c1) and 0.87 to 0.55 (c3), and the
     posteriors differ by total variation 0.11-0.28. The effect on fitted functions at 200 trees is
     unmeasured.
   - The normalized prior instead keeps CGM's tree marginal: the constraint changes leaf values, not which
     trees are likely.
   - mBART: the paper (arXiv 1612.01619v3, eqs. 3.1 and 3.3) states the per-tree normalized prior, and its move
     (eqs. 4.11 and 4.13) divides by the local d*. The software (remcc/mBART_shlib, bd.cpp, coninteg1 and
     coninteg2) accumulates the prior mass sumpr and never uses it, so it samples the unnormalized target on a
     grid.

Recommendation: keep the normalized prior with option 5: no budget, the guard as a control argument, and the
free-bound shortcut. It costs nothing at the default tree count and keeps the documented model for every fit
that completes. Its price is speed in 1-5 tree fits and a stop in some of them. The unnormalized prior's case
is real: no count, no guard, what mBART's software samples, and a plain reading as BART conditioned on
monotone trees; its price is fewer splits on the constrained predictors. Evidence that would change the call: a
fit at 10 or more trees whose counting outweighs its sweep or that reaches the guard.

## Constraints

- Exact for the stated prior (dec-B14). Counts are held in double: exact to 2^53, with relative error under
  1e-12 beyond.
- Unconstrained samplers byte-identical: the new seam compiles out, like the three existing monotone seams.
  No dbarts.h change.
- Out of scope: change moves under the constraint, quadrature speed (TODO monotone-leaf-quadrature), and
  reconciling a chi k hyperprior with the truncated law.

## Steps

1. Order counter (engine, beside the monotone geometry). Down-set keys are multi-word bitsets.
   - Relation: build the tree's relation matrix with the adjacency test of `monotoneNeighborBounds`, and split
     it into components.
   - Count: run the layered down-set DP per component, keeping two layers (every layer only for step 6's
     draw), and return the whole-tree log Z_T = sum over components of log e(C) - lgamma(|C| + 1). A component
     whose down-sets pass the guard G (a dbartsControl argument, default 2^24; Decision) stops the count at
     once and fails.
   - A move changes only the component of the touched leaf in T0 and the component(s) holding the two children
     in T*. After a free-axis birth the children may share a component: take the distinct components. Keep the
     current tree's component counts, so a move counts only the components it creates.
2. Seam. Add an optional leaf concept, `logTreeNormalizer`, declared by the monotone leaf.
   [`birthOrDeathMove`](../../src/bartcore/moves.hpp) evaluates it in each state: before and after
   `tree.birth`, before and after `orphanChildren`. It adds log Z_T0 - log Z_T* to the log prior ratio (logged
   in the census's prior column).
   - Free bounds: a birth never lowers e (a linear extension of T0 with the split leaf replaced by its two
     children in order is one of T*, and distinct extensions stay distinct), so a birth's Z_T0 / Z_T* =
     (L0 + 1) e(T0) / e(T*) is at most L0 + 1 and a death's is at least 1 / L0, L0 the current leaf count (0
     of 1,164 random births lower e; monotone-order-size.R closure mode). Draw u first; with r1 the rest of
     the ratio, reject a birth without counting when u > r1 (L0 + 1), and accept a death without counting
     when u < r1 / L0. This is plain Metropolis-Hastings with no loss of acceptance; the savings are
     unmeasured. (A two-stage acceptance, min(1, a) min(1, b), is also exact but lowers acceptance.)
   - A count over the guard stops the run with an error the bridge raises; it never rejects, and the chain
     keeps T0, so a caller can catch the error, raise G and continue.
3. Drop the d terms: `priorMass` in `oneLeafLogMarginal`, `denom` in `twoLeafCoupledLogMarginal`.
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
   - The alternative, refusing mutations that strand a leaf, breaks the embedded use the sampler exists for.
6. Exact prior leaf draw. Per component, draw a uniform linear extension by backward sampling on the DP
   counts, every layer kept: remove a maximal element x with probability e(D - x) / e(D). Draw |C| iid
   N(0, c scale / k), sort them, and assign them in extension order; isolated leaves draw alone. This replaces
   the rejection loop and its 1e6 cap in `drawFromPriorForTree`.
7. Reachability:
   - A bridge refusal after [`parseProposalProbs`](../../src/R_interface_bartcore.cpp), at creation and in
     [`bartcore_setModel`](../../src/R_interface_bartcore.cpp). It is keyed on the engine's leaf kind
     ([`LeafModelKind`](../../src/bartcore/model.hpp)), not on the incoming model, and refuses any nonzero swap,
     change, perturb or rule_gibbs share. The all-zero frozen mixture stays allowed.
   - setControl mirrors creation: a defaulted mixture is rewritten to birth/death silently, and an explicit
     non-birth/death one is refused.
   - growForestFromRoot reseeds mu to the all-zero feasible seed before its draw.
   - [`Chain::installForest`](../../src/bartcore/chain.hpp) reseeds infeasible leaves.
   - growForestFromRoot and installForest both prune a tree with a component over the guard by deaths of its
     deepest nodes until every component is within it. This ends at the root at worst, and a start state is
     not part of the target.
   - [`Chain::rebuildLiveForest`](../../src/bartcore/chain.hpp) (setState, copy, reload) refuses an infeasible
     tree or one over the guard. With steps 4 and 5 every state the sampler produces passes this check.
   - [`Chain::sampleTreesFromPrior`](../../src/bartcore/chain.hpp) needs no new predicate: the prior is not
     restricted, and passing the guard takes 25 or more leaves in one component, which carry less than 1.5e-18
     under CGM(0.95, 2).
8. tests/cpp:
   - The count against brute-force permutation counts on hand-built trees (1-3 axes, mixed directions, N, a
     star, a chain over 64 leaves) and 200 random trees.
   - Log Z against the enumeration's e / L!.
   - RNG-free ratio tests on A < B -> A < B1 < B2 (constrained split) and A < B -> {A < B1, B2} (free split)
     with pinned mu_A. Compare the move's log alpha with the closed form of the corrected statement in Context
     (it must differ from the current code's value).
   - The pair redraw on a design with P_post(lower <= upper) ~ 1e-4, where the old loop exhausts almost every
     time: every draw is feasible, and the upper leaf's draws match the quadrature CDF (KS).
   - The linear-extension draw is uniform over extensions (chi-square on a 5-leaf N-plus-chain), and the guard
     fires.
   - Free bounds: over random trees and moves, a move decided without counting gets the same decision as with
     the count, for the same u.
   - [`testMonotoneMarginal`](../../tests/cpp/test_model.cpp) loses its d_* = 1/2 check.
9. tinytest ([test-monotone.R](../../inst/tinytest/test-monotone.R)):
   - a setControl change mix errors, and a defaulted one is rewritten;
   - setModel without the monotone attribute cannot install a change mix;
   - installTrees from an unconstrained donor leaves the fit monotone at once, and setState of that state is
     refused;
   - growFromRoot plus one sweep is monotone;
   - after 2,000 sweeps on data decreasing along the constrained axis, copy() and setState round-trip;
   - setPredictor(forceUpdate = TRUE) stranding a leaf keeps the fit monotone;
   - a guard set low errors, leaves the state unchanged, and the run continues after the guard is raised.
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
11. A monotone scenario in benchmarks/R/equivalence.R (x1 and x2 constrained, 20 trees), with the other 53
    scenarios bitwise. Its MANIFEST row names the enumeration gate as the ORACLE (P17).
12. Docs:
    - monotone.md sections 4, 9 and 11 and the Plan-vs-code note restate B' with the whole-tree normalizer, and
      record that mBART's d-normalized eq. 4.11 targets neither this prior nor the software's.
    - The prior statement stays verbatim; the `monotone` argument's help names the guard, its error, the
      few-tree cost, and that rerunning with new seeds until a run completes selects smaller trees.
    - dec-B16 is marked superseded in part by dec-B144.
    - Status lines and INDEX at landing, and the TODO item removed.

## Verification

- `cd tests/cpp && make && ./test_bartcore`: the count, ratio, redraw, extension-draw and guard checks pass.
- `R CMD INSTALL --preclean -l <lib> .`, then `R_LIBS=<lib> Rscript -e 'tinytest::test_package("dbarts")'`.
- `R_LIBS=<lib> Rscript benchmarks/R/monotone-exact-enumeration.R quick`: every group passes. This takes about
  6 min; full mode runs 900k draws. The mutation run restores the d divisions and drops the Z term (then
  `touch` the header and reinstall), and must fail like the current engine.
  - Current engine, quick: c1 p 1.6e-24 and 6.5e-35 (two root rules); c2 3.8e-10 and 4.3e-5; c3 1.4e-15;
    cN 1.9e-29.
  - `prototype` mode (the corrected move in R, with an exact pair redraw): every group p >= 0.07.
  - `prototype-old` reproduces the current engine (c1 T^2 196 and 246, against the engine's 163 and 253).
  - Against an e(N) miscount of 10%, cN's noncentrality is ~80, so power is about 1.
  - `zcheck` mode verifies Z_T = e / L!.
- `R_LIBS=<lib> Rscript benchmarks/R/monotone-reference.R quick`: parts (a) and (b) still pass. Then run the
  whole exact-gates.yaml list with `quick`.
- `Rscript benchmarks/R/equivalence.R compare <current>`: 53 scenarios "identical draws (same RNG stream)" and
  no "max |z|". BCF and multinomial compare identical. The snapshot files carry no monotone fit.
- Release level: the SBC arm in [Monotone arm: design](sbc-family-tiers.md#monotone-arm-design). monotone-1
  (0.4 s per replicate) flags the current move and must pass; then the 20-tree arm (~85 min at R 200) must pass
  before admission.
- Speed: on a quiet machine, monotone sweep time at 20 trees, 1 and 2 constrained predictors, within 5% of
  today; at 1 and 5 trees, the slowdown is recorded against the Decision's estimates. bench-sampler compare
  unchanged on the unconstrained paths.
- `Rscript tools/check-doc-freshness.R .` passes.
