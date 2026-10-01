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
- Empty leaves. [`dbartsSampler$installTrees`](../../man/dbartsSampler-class.Rd) strands a member-empty
  leaf: [`Chain::rebuildLiveForest`](../../src/bartcore/chain.hpp) repartitions a donor's trees over this
  sampler's rows and collapses nothing. [`MonotoneConstantGaussianLeaf::drawOneLeaf`](../../src/bartcore/model.hpp)
  pins such a leaf at 0 with no bound check, so the cone breaks and Z_T is no longer e / L!. This recurs in
  embedded use.
  - Corrected in implementation: [`dbartsSampler$setPredictor`](../../man/dbartsSampler-class.Rd) with
    forceUpdate, [`dbartsSampler$setData`](../../man/dbartsSampler-class.Rd) and
    [`dbartsSampler$setCutPoints`](../../man/dbartsSampler-class.Rd) do not strand one: they collapse emptied
    subtrees, and setData also remaps splits onto a new grid (an unforced update that would empty a leaf
    rolls back). A collapse merges leaves at their weighted mean, and the merged leaf can border neighbours
    it did not border before (a split on a free axis whose one side empties), so a tree can leave the cone
    with no empty leaf in it; a remap can relate leaves the same way. Step 7's reseed covers these paths.
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
- Default (dec-B151, 2026-09-30): "joint". "I think we'll want to use \"joint\" since it seems
  indistinguishable." Not chosen: "leaf". Its tree prior, cgm()'s defaults or mBART's base 0.25 and power 0.8
  whenever a constraint is present, is measured first ("Measure first.") by the second study in "Tree prior
  under joint" below; the maintainer makes the final ruling on its results: "Bring the results to me and let
  me make the final ruling."
- Tree prior under "joint" (dec-B152, 2026-10-01): cgm()'s defaults, with mBART's values documented in the
  monotone() help as an opt-in. "Yes, use option 1." Not chosen: mBART's values whenever a constraint is
  present, and a search for an intermediate setting.

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

### Results, 2026-09-30

Run by [monotone-feel-study.R](../../benchmarks/R/monotone-feel-study.R) on a library built at d494eb41, arm64
macOS, eight fits at once on a machine shared with other jobs (1-minute load average 80-160 on 10 cores).
Every arm, tree count and truth ran; nothing was cut. Fits ran 8 replicates, twice the spec's 4, since the
leaf-joint gap on the step truth at n 200 was a question more replicates could settle; about 1.1 CPU-hours.
Per fit and per prior arm in
[monotone-feel-d494eb41.csv](../../benchmarks/baselines/monotone-feel-d494eb41.csv). "joint, mBART" is "joint"
under cgm(power = 0.8, base = 0.25); "free" is the unconstrained fit. A paired difference is an arm minus
"leaf" on the same data set; its standard error is over the 8 replicates.

Prior predictive, along x1 at five fixed x2 values, 500 draws; the response spans 1, so a rise of 1 is the
whole observed range of y. Levels are distinct values on 201 points; largest jump is a share of the rise; flat
intervals a share of the 200. The free arm's curves are read sorted (the draw with its values rearranged into
increasing order); unsorted, 18%, 7% and 1% of its intervals fall at 200, 50 and 5 trees.

| trees | arm | levels | largest jump | flat intervals | rise | curves with no rise | sd of f | x1 split share | leaves per tree |
|---|---|---|---|---|---|---|---|---|---|
| 200 | leaf | 74.3 | 0.05 | 0.63 | 2.94 | 0.00 | 0.28 | 0.50 | 2.48 |
| 200 | joint | 51.1 | 0.07 | 0.75 | 1.64 | 0.00 | 0.26 | 0.31 | 2.24 |
| 200 | joint, mBART | 13.1 | 0.22 | 0.94 | 0.31 | 0.00 | 0.26 | 0.32 | 1.22 |
| 200 | free | 74.3 | 0.10 | 0.63 | 0.40 | 0.00 | 0.26 | 0.50 | 2.48 |
| 50 | leaf | 28.8 | 0.12 | 0.86 | 1.45 | 0.00 | 0.28 | 0.50 | 2.48 |
| 50 | joint | 16.9 | 0.18 | 0.92 | 0.82 | 0.00 | 0.28 | 0.31 | 2.24 |
| 50 | joint, mBART | 4.2 | 0.59 | 0.98 | 0.16 | 0.04 | 0.27 | 0.32 | 1.22 |
| 50 | free | 29.1 | 0.16 | 0.86 | 0.38 | 0.00 | 0.25 | 0.50 | 2.49 |
| 5 | leaf | 4.3 | 0.60 | 0.98 | 0.45 | 0.02 | 0.28 | 0.52 | 2.48 |
| 5 | joint | 2.8 | 0.78 | 0.99 | 0.27 | 0.13 | 0.26 | 0.32 | 2.25 |
| 5 | joint, mBART | 1.4 | 0.95 | 1.00 | 0.06 | 0.70 | 0.25 | 0.32 | 1.23 |
| 5 | free | 4.4 | 0.60 | 0.98 | 0.28 | 0.02 | 0.25 | 0.50 | 2.46 |

Fit quality on the test grid: RMSE of the posterior mean against the truth and 95% interval coverage.

| truth | n | trees | RMSE leaf | joint | joint, mBART | free | coverage leaf | joint | joint, mBART | free |
|---|---|---|---|---|---|---|---|---|---|---|
| ramp | 2000 | 200 | 0.022 | 0.023 | 0.018 | 0.024 | 0.99 | 0.98 | 0.94 | 0.98 |
| ramp | 2000 | 50 | 0.019 | 0.021 | 0.019 | 0.019 | 0.94 | 0.91 | 0.78 | 0.95 |
| ramp | 2000 | 5 | 0.020 | 0.022 | 0.022 | 0.021 | 0.74 | 0.66 | 0.67 | 0.71 |
| ramp | 200 | 200 | 0.044 | 0.043 | 0.032 | 0.045 | 0.99 | 0.99 | 0.98 | 0.99 |
| ramp | 200 | 50 | 0.040 | 0.042 | 0.035 | 0.041 | 0.99 | 0.98 | 0.93 | 0.99 |
| ramp | 200 | 5 | 0.043 | 0.044 | 0.045 | 0.043 | 0.84 | 0.79 | 0.76 | 0.79 |
| step | 2000 | 200 | 0.056 | 0.055 | 0.046 | 0.058 | 0.96 | 0.96 | 0.98 | 0.98 |
| step | 2000 | 50 | 0.050 | 0.049 | 0.044 | 0.051 | 0.97 | 0.98 | 0.99 | 0.98 |
| step | 2000 | 5 | 0.044 | 0.044 | 0.051 | 0.049 | 0.96 | 0.97 | 0.97 | 0.98 |
| step | 200 | 200 | 0.117 | 0.104 | 0.075 | 0.091 | 0.91 | 0.94 | 0.95 | 0.97 |
| step | 200 | 50 | 0.086 | 0.084 | 0.063 | 0.085 | 0.96 | 0.97 | 0.98 | 0.98 |
| step | 200 | 5 | 0.070 | 0.077 | 0.044 | 0.061 | 0.96 | 0.98 | 0.98 | 0.99 |
| hinge | 2000 | 200 | 0.028 | 0.027 | 0.018 | 0.028 | 0.96 | 0.97 | 0.96 | 0.98 |
| hinge | 2000 | 50 | 0.020 | 0.022 | 0.018 | 0.023 | 0.96 | 0.93 | 0.90 | 0.93 |
| hinge | 2000 | 5 | 0.020 | 0.021 | 0.020 | 0.019 | 0.84 | 0.82 | 0.71 | 0.80 |
| hinge | 200 | 200 | 0.048 | 0.045 | 0.034 | 0.047 | 0.98 | 0.99 | 0.97 | 1.00 |
| hinge | 200 | 50 | 0.042 | 0.045 | 0.036 | 0.045 | 0.97 | 0.98 | 0.92 | 0.99 |
| hinge | 200 | 5 | 0.037 | 0.042 | 0.044 | 0.037 | 0.86 | 0.84 | 0.81 | 0.81 |
| interaction | 2000 | 200 | 0.085 | 0.090 | 0.078 | 0.084 | 0.94 | 0.92 | 0.87 | 0.96 |
| interaction | 2000 | 50 | 0.082 | 0.081 | 0.088 | 0.081 | 0.88 | 0.85 | 0.76 | 0.88 |
| interaction | 2000 | 5 | 0.106 | 0.119 | 0.119 | 0.115 | 0.68 | 0.62 | 0.58 | 0.65 |
| interaction | 200 | 200 | 0.148 | 0.146 | 0.130 | 0.141 | 0.97 | 0.98 | 0.96 | 0.99 |
| interaction | 200 | 50 | 0.144 | 0.159 | 0.161 | 0.140 | 0.96 | 0.95 | 0.88 | 0.97 |
| interaction | 200 | 5 | 0.215 | 0.207 | 0.237 | 0.219 | 0.74 | 0.74 | 0.63 | 0.68 |

Paired, "joint" minus "leaf": the RMSE difference is within two standard errors of zero, or under 0.003, in
all but five cells. "joint" is worse in four (interaction n 200 at 50 trees +0.014, se 0.006; interaction n
2000 at 5 trees +0.013, se 0.004; hinge n 200 at 5 trees +0.005 and at 50 trees +0.003, se 0.002 and 0.001)
and better in one, the step at n 200 and 200 trees (-0.012, se 0.002). Mean interval width, averaged over the
truths, as a share of "leaf"'s:

| n | trees | joint | joint, mBART | free |
|---|---|---|---|---|
| 2000 | 200 | 1.00 | 0.62 | 1.11 |
| 2000 | 50 | 1.00 | 0.62 | 1.09 |
| 2000 | 5 | 0.99 | 0.85 | 1.03 |
| 200 | 200 | 0.99 | 0.70 | 1.07 |
| 200 | 50 | 1.00 | 0.65 | 1.07 |
| 200 | 5 | 0.99 | 0.81 | 1.00 |

Partial dependence along x1 (the posterior mean averaged over the grid's rows), averaged over the truths: its
RMSE against the truth's; the number of steps a single draw's curve takes over the 100 grid intervals; and the
visible steps of the posterior-mean curve, increments over 5% of the truth's rise. Every constrained curve is
nondecreasing; the free arm's posterior-mean curve falls, at its largest dip, by 2-5% of the truth's rise on
average at 50 and 200 trees, and under 1% at 5.

| n | trees | PD RMSE leaf | joint | joint, mBART | free | steps per draw leaf | joint | joint, mBART | free | visible steps leaf | joint | joint, mBART | free |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 2000 | 200 | 0.027 | 0.027 | 0.029 | 0.032 | 34 | 27 | 17 | 60 | 1.8 | 2.6 | 4.9 | 3.9 |
| 2000 | 50 | 0.029 | 0.030 | 0.032 | 0.031 | 19 | 16 | 12 | 24 | 4.7 | 5.7 | 6.9 | 5.1 |
| 2000 | 5 | 0.034 | 0.036 | 0.035 | 0.035 | 13 | 13 | 12 | 11 | 7.4 | 7.2 | 7.5 | 6.9 |
| 200 | 200 | 0.054 | 0.049 | 0.043 | 0.050 | 39 | 29 | 16 | 62 | 1.7 | 1.5 | 2.3 | 2.1 |
| 200 | 50 | 0.047 | 0.051 | 0.053 | 0.051 | 17 | 13 | 8 | 24 | 3.1 | 4.0 | 4.8 | 3.8 |
| 200 | 5 | 0.064 | 0.067 | 0.065 | 0.061 | 7 | 6 | 5 | 6 | 5.4 | 4.8 | 4.6 | 4.9 |

The step truth at n 200 and 200 trees drives "leaf"'s larger PD RMSE there (0.104, against "joint" 0.091,
"joint, mBART" 0.071, free 0.077).

Variable importance: varcount shares, averaged over the truths, and splits per tree.

| n | trees | x1 leaf | joint | joint, mBART | free | x3 (noise) leaf | joint | joint, mBART | free | splits per tree leaf | joint | joint, mBART | free |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 2000 | 200 | 0.19 | 0.14 | 0.39 | 0.36 | 0.39 | 0.42 | 0.25 | 0.31 | 1.27 | 1.21 | 0.25 | 1.37 |
| 2000 | 50 | 0.36 | 0.31 | 0.67 | 0.46 | 0.27 | 0.29 | 0.10 | 0.22 | 1.21 | 1.14 | 0.43 | 1.26 |
| 2000 | 5 | 0.73 | 0.74 | 0.70 | 0.74 | 0.07 | 0.06 | 0.07 | 0.06 | 4.56 | 4.14 | 4.00 | 3.58 |
| 200 | 200 | 0.20 | 0.14 | 0.31 | 0.35 | 0.40 | 0.42 | 0.32 | 0.32 | 1.36 | 1.29 | 0.30 | 1.46 |
| 200 | 50 | 0.28 | 0.22 | 0.49 | 0.39 | 0.34 | 0.37 | 0.20 | 0.29 | 1.36 | 1.27 | 0.35 | 1.43 |
| 200 | 5 | 0.67 | 0.64 | 0.78 | 0.69 | 0.11 | 0.12 | 0.05 | 0.09 | 2.23 | 1.91 | 1.51 | 1.83 |

Held-out log predictive score, mean per point on 1000 fresh points (higher is better):

| truth | n | trees | leaf | joint | joint, mBART | free |
|---|---|---|---|---|---|---|
| ramp | 2000 | 200 | 0.904 | 0.901 | 0.915 | 0.899 |
| ramp | 2000 | 50 | 0.912 | 0.909 | 0.911 | 0.912 |
| ramp | 2000 | 5 | 0.914 | 0.906 | 0.907 | 0.905 |
| ramp | 200 | 200 | 0.832 | 0.833 | 0.869 | 0.823 |
| ramp | 200 | 50 | 0.839 | 0.829 | 0.852 | 0.837 |
| ramp | 200 | 5 | 0.823 | 0.824 | 0.818 | 0.817 |
| step | 2000 | 200 | 0.317 | 0.317 | 0.332 | 0.314 |
| step | 2000 | 50 | 0.328 | 0.329 | 0.334 | 0.328 |
| step | 2000 | 5 | 0.334 | 0.337 | 0.331 | 0.332 |
| step | 200 | 200 | 0.179 | 0.213 | 0.285 | 0.243 |
| step | 200 | 50 | 0.259 | 0.260 | 0.303 | 0.267 |
| step | 200 | 5 | 0.306 | 0.298 | 0.308 | 0.306 |
| hinge | 2000 | 200 | 0.846 | 0.846 | 0.861 | 0.844 |
| hinge | 2000 | 50 | 0.858 | 0.855 | 0.861 | 0.856 |
| hinge | 2000 | 5 | 0.856 | 0.858 | 0.859 | 0.860 |
| hinge | 200 | 200 | 0.758 | 0.760 | 0.809 | 0.757 |
| hinge | 200 | 50 | 0.773 | 0.763 | 0.802 | 0.767 |
| hinge | 200 | 5 | 0.801 | 0.787 | 0.773 | 0.797 |
| interaction | 2000 | 200 | -0.126 | -0.128 | -0.126 | -0.129 |
| interaction | 2000 | 50 | -0.126 | -0.127 | -0.138 | -0.130 |
| interaction | 2000 | 5 | -0.172 | -0.175 | -0.163 | -0.165 |
| interaction | 200 | 200 | -0.209 | -0.209 | -0.200 | -0.205 |
| interaction | 200 | 50 | -0.225 | -0.237 | -0.251 | -0.220 |
| interaction | 200 | 5 | -0.361 | -0.343 | -0.368 | -0.359 |

Paired, the only score difference between "joint" and "leaf" both over 0.015 and beyond two standard errors is
the step at n 200 and 200 trees (+0.034, se 0.003).

Time and counts. CPU per sweep over "joint"'s in the same job, median over its 32 jobs (the four arms of a job
ran back to back; the machine's load inflates absolute times, so only these within-job ratios compare); and
"joint"'s and the free arm's CPU ms per sweep at n 2000 for scale.

| n | trees | leaf | joint, mBART | free | joint ms (n 2000) | free ms (n 2000) |
|---|---|---|---|---|---|---|
| 2000 | 200 | 1.04 | 1.29 | 0.18 | 12.8 | 2.3 |
| 2000 | 50 | 1.05 | 1.16 | 0.16 | 3.8 | 0.59 |
| 2000 | 5 | 1.01 | 0.99 | 0.22 | 0.59 | 0.13 |
| 200 | 200 | 1.04 | 1.27 | 0.13 | | |
| 200 | 50 | 1.06 | 1.19 | 0.12 | | |
| 200 | 5 | 0.95 | 0.97 | 0.13 | | |

Under "leaf" a fit ran 0.3-0.6 order counts per tree per sweep (2,800 at 5 trees, 26,600-28,800 at 50,
110,000-120,000 at 200). The largest count in any fit held 16 leaves and 598 down-sets (5 trees, n 2000); at
50 trees 9 leaves and 126 down-sets, at 200 trees 4 and 14. The slow-count warning fired in no fit. The
largest count's measured time was 0.08 s, on a 2-leaf count; counts this size take microseconds on an idle
thread, so these times are the load descheduling the thread, not counting.

Reading, against the criteria above. For "leaf":

- Fit quality against "joint" at cgm() defaults: RMSE, coverage and held-out score match in most cells, and
  where they differ "leaf" is ahead more often (RMSE on the interaction at n 200 and 50 trees and n 2000 and 5
  trees, and on the hinge at n 200 and 5 or 50 trees; coverage on the ramp at n 2000 and 50 or 5 trees, 0.94
  and 0.74 against 0.91 and 0.66). The exception is the steep step at n 200 and 200 trees, where "leaf" is
  behind "joint" and the free fit alike (RMSE 0.117 against 0.104 and 0.091, coverage 0.91 against 0.94 and
  0.97, score 0.179 against 0.213 and 0.243). Against the mBART-tuned "joint" the criterion fails on the step:
  that arm has the lowest RMSE and best score there at every n and at 50 and 200 trees.
- Variable importance: at 50 and 200 trees both priors report the constrained predictor below the free fit,
  and the noise predictor above it. "leaf" is closer (x1 0.19 against free 0.36 at n 2000 and 200 trees;
  "joint" 0.14), but neither is close. The mBART-tuned "joint" is closest at 200 trees and over-reports x1 at
  50 (0.67 against 0.46).
- Prior draws: "leaf" matches the free prior in levels, flat intervals, split share on x1 and tree size,
  exactly as a per-tree sort of BART's leaves would. It differs in the rise: each tree's sorted leaves add up,
  so the prior curve climbs 0.45, 1.45 and 2.94 response ranges at 5, 50 and 200 trees, with a pointwise sd of
  0.28, where the free draws rearranged climb 0.3-0.4. "joint" has fewer levels (51 against 74 at 200 trees)
  and about half as many splits on x1 (split share 0.31 against 0.50, in smaller trees), and climbs 1.64 at
  200 trees.
- Time: within 4-6% of "joint" at 50 and 200 trees, and no slow-count warning anywhere.

For "joint": "leaf" does not over-split the constrained axis (posterior splits per tree 1.27 against 1.21,
free 1.37), its intervals are no wider (width within 1% of "joint"'s), its held-out score is worse only on the
one step cell, and its counts stay tiny at 50 trees and more. "joint" at cgm() defaults matches "leaf" on fit
quality everywhere except that it is slightly behind in the cells listed above and ahead on that step cell.

The mBART-tuned case, reported separately: that tree prior gives a different model, not a closer match. Trees
hold 0.2-0.4 splits at 50 and 200 trees against 1.2-1.4; the prior curve takes 13 levels at 200 trees and is
flat in 70% of draws at 5; intervals are 30-38% narrower at 50 and 200 trees. On RMSE and score it is best or
tied on every truth at 200 trees, and on the step, hinge and ramp at 50, but under-covers the smooth truths at
50 trees (ramp 0.78 and interaction 0.76 at n 2000), fits the interaction worse at 50 and 5 trees, and runs
16-29% slower per sweep than "joint" at 50 and 200 trees.

How a user would feel each default. Under either prior a constrained fit runs 4.5-8.5 times the CPU per sweep
of the same fit without the constraint; that cost is the constraint's, not the prior's. With "leaf" a user
sees fits, intervals and scores like "joint"'s and close to the free fit's; a smoother partial dependence (the
fewest visible steps at n 2000 and 200 trees, and no dips, which the free fit shows); the constrained
predictor's importance about half the free fit's at 200 trees and three quarters at 50; and, if they draw from the prior, curves that rise steeply
with the tree count. With "joint" at cgm() defaults a user sees nearly the same fits with a lower importance
on the constrained predictor and a gentler prior rise. A "joint" default tuned as mBART does would feel
different from both: sharper, fewer steps, tighter intervals that under-cover smooth truths at 50 trees, and a
tree prior that differs from the one users get without a constraint.

### Tree prior under "joint"

The feel study's one design (a constrained, a free and a noise predictor) left the tree prior under "joint"
open: at 200 trees mBART's cgm(power = 0.8, base = 0.25) had the lowest error and best score on every truth,
with intervals 30-38% narrower, but at 50 trees it under-covered the ramp (0.78) and the interaction (0.76)
and fit the interaction worse. This second study widens the designs. A rule fixed before it runs gives a
verdict; the maintainer rules on the results.

- Arms: "joint" under cgm() (base .95, power 2); "joint" under cgm(power = 0.8, base = 0.25); the
  unconstrained fit under cgm() as reference. Tree counts 200 and 50.
- Designs, predictors uniform on [0, 1]^p:
  - friedman: p 10, f = 10 sin(pi x1 x2 / 2) + 20 (x3 - 0.5)^2 + 10 x4 + 5 x5, increasing in x1, x2, x4 and
    x5 (constrained), x3 free, x6-x10 noise;
  - additive: p 10, f = 2 x1 + 1{x2 > 0.5} + exp(2 x3) / e^2 - sin(2 pi x4), increasing in x1, x2 and x3
    (constrained), x4 free, x5-x10 noise;
  - interaction: p 5, f = x1 (1 + 4 x2) + x1 x3 + sin(2 pi x2), increasing in x1 (constrained), x2 and x3
    free, x4 and x5 noise.
- Noise sd a third of the truth's sd, and equal to it; n 200 and 2000; 8 replicates; one chain of 500 burn-in
  and 500 kept sweeps. Measures as in the feel study's fits: RMSE of the posterior mean against the truth and
  95% interval coverage on 2000 held-out points from the design's law, held-out log predictive score, interval
  width, varcount shares, time per sweep.
- Rule (advisory): mBART's values are favoured as the tree prior under every monotone fit if (a) at 200 trees, in every design, n
  and noise cell, their RMSE and score are no worse than cgm()'s defaults beyond two paired standard errors,
  and (b) their 95% coverage is at least 0.90 in every cell at 200 and 50 trees. Otherwise cgm()'s defaults
  are.

### Results, 2026-09-30

Run by [monotone-treeprior-study.R](../../benchmarks/R/monotone-treeprior-study.R) on a library built at
7593a0d4, arm64 macOS, six fits at once; all 192 jobs, 8 replicates; per fit in
[monotone-treeprior-7593a0d4.csv](../../benchmarks/baselines/monotone-treeprior-7593a0d4.csv). "Low" noise is
a third of the truth's sd, "high" equal to it.

- At 200 trees mBART's values have lower RMSE (4-25%) and a better held-out score than cgm()'s defaults in all
  12 cells, mostly beyond two paired standard errors. Their intervals are 23-34% narrower and cover 0.84, 0.88
  and 0.86 (friedman, additive, interaction) at n 2000 and low noise, where the defaults cover 0.93-0.96 and
  the free fit 0.96; elsewhere 0.93-0.97.
- At 50 trees mBART's values are worse at n 2000 and low noise in all three designs beyond two standard
  errors, better at n 200, and cover 0.70-0.94, under 0.90 in 9 of 12 cells; the defaults miss 0.90 in 5 and
  the free fit in 3.
- The constrained predictors' varcount share at 200 trees: defaults 0.15-0.17, mBART 0.29-0.38, free
  0.31-0.34. Splits per tree 1.2-1.3 against 0.3-0.66. mBART costs 4-21% more CPU per sweep. No monotone draw
  fell along any line.
- The advisory rule: clause (a) passes in every cell, clause (b) fails in 12 of 24, so it favours cgm()'s
  defaults. The maintainer ruled the same (dec-B152).

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
3b. Leaf geometry on factor and missing-value axes, before the counter is wired in. monotoneLeafBox reads
   every axis through Tree::splitInterval, which on a subset split (an unordered factor) returns the low
   bits of the level mask with no cuts, and threshold boxes ignore where rows with a missing value go
   (missingGoesRight). So monotoneNeighborBounds, monotoneTreeIsFeasible and buildMonotoneLeafOrder both
   miss pairs the constraint requires and add pairs it does not: found in stage 3's review, a one-tree fit
   with a free 4-level factor decreased along x1 for one level in 1.5% of draws, and with 40% missing in a
   free predictor the fit at the missing value decreased along x1 in every draw under two of three seeds.
   Ordered factors and cut remaps are unaffected apart from missing values; pooled factors are affected.
   The fix: one helper gives each leaf, per axis, either a code interval plus whether it reaches missing
   values (threshold axes) or its reachable level set (subset axes, pooled ones included); two leaves share
   a free axis when their intervals overlap or both reach missing values, or their level sets intersect;
   adjacency stays interval-touching on constrained axes. Leaf regions stay products, so the order
   argument and the union-of-children merge hold. The bounds, the feasibility check, the builder and the
   tests all use the helper. Oracles: a point-based geometry oracle in tests/cpp with factor and missing
   columns (today's oracles share the bug), an exact-gate design with a free factor, and a tinytest that a
   fitted surface is monotone per factor level and at the missing value. Draws change only for fits that
   split on a factor or carry missing values in a free predictor. What a missing value in a CONSTRAINED
   predictor means for the order is settled in this stage's design note, not assumed.
4. The seam, the dropped d terms, the prior flag and the gate under both priors (steps 2, 3, 13, 10), step 6's
   switch of `drawFromPriorForTree` to the exact draw, and step 7's sampleTreesFromPrior bullet for both
   priors; prior = "joint" is accepted from here. Per-move count timing is recorded through the move census
   (BARTCORE_MOVE_CENSUS; step 2). Then the checkpoint below, and stop. Until commit 5's setModel refusal,
   setModel silently accepts a changed monotone prior (the engine keeps the old one); commits 4 and 5 are
   pushed together, so no pushed tip carries it.
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

### Stage 3b: leaf geometry

A point is one value per predictor: a code on a threshold axis (numeric or ordered factor), a level on a
subset axis (unordered factor, pooled or not), or missing where the training column has missing values.
Every rule tests one predictor, so the
points prediction routes to a leaf form a product over predictors, and one helper
(`MonotoneLeafGeometry`) records it per leaf and split variable:

- threshold axis: the code interval [lo, hi], and whether a missing value reaches the leaf (every ancestor
  rule on the axis sends missing values to the leaf's side);
- subset axis: the reachable level set, with the missing position when the column has missing values, each
  ancestor mask filtering it, as `Tree::reachableCategories` and its pooled analogue do.

Calls:

- Two leaves share an axis when some value reaches both: the intervals overlap or both reach missing values,
  or the level sets intersect. j is below k along a constrained axis when j's interval ends one code below
  k's start (direction -1 flips it) and they share every other axis. Only a threshold axis can be
  constrained: R refuses a direction on an unordered factor, and the facade drops the whole constraint when
  any direction sits on a subset axis (`monotoneConstraintIsActive`, pre-existing).
- A missing value in the constrained predictor x1 has no position along x1, so x1's own missing flag takes
  no part in adjacency along x1. The order claims that the fit is monotone in x1 along every line of points
  with x1 observed, the other predictors at any values, missing included, and claims nothing at x1 missing.
  A missing value in another constrained predictor is shared like any value.
- Missing is a value of an axis only when its training column has missing values, as the rules' reachable
  sets already hold it; predict refuses a missing test value in a column that had none. Treating it as a
  value everywhere was tried and rejected in review: where no rule was drawn for it, a missing value goes
  left at every split, and which side of a level split is called left is an arbitrary label, so the order
  (and log Z) of one partition would depend on its labelling, and ordinary factor fits would carry
  constraints at points no one can predict at.
- A setPredictor, a row-by-row update or setData that brings a column its first missing values gives its
  axis that value, which can relate leaves the order did not. The accepted update paths reseed any tree
  that leaves the cone to all-zero, as stage 1's collapse and remap paths do; the unforced and row-by-row
  paths gained that reseed here.

Why the order argument stands:

- The box test is the point test: regions are products, so two leaves are related along x1 exactly when a
  point of one and a point of the other differ only in x1, by one code. Transitivity then orders every two
  points on a line with x1 observed.
- Acyclic: a rule sends each value, missing included, to one side, so the two sides of a split share no value
  on its axis and relate only along it, one way. Context's cut argument applies unchanged, and a pair relates
  along one axis only.
- A merged leaf's region is the union of its children's, so its point pairs, and relations, are theirs.

What a user sees: the constraint holds as documented, on every factor level and wherever another predictor
is missing. Draws change for monotone fits whose trees split on an unordered factor, or route a missing value
right on a predictor other than the only constrained one (missing values in a free predictor, or in one of
two or more constrained ones); every other fit draws as before. The help's wording for a missing
constrained value waits for step 12.

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
   - As implemented (commit 3), in model.hpp beside `monotoneNeighborBounds`: `buildMonotoneLeafOrder`,
     `monotoneLogExtensions`, `monotoneLogNormalizer`, `monotoneLogNormalizerRatio`, `monotoneDrawExtension`
     and step 6's `monotoneDrawPriorLeaves`, none called yet. The relation holds adjacencies only, not their
     closure: the DP adds an element once its direct predecessors are in, which yields the same down-sets.
     Each layer is scaled by a power of two, so counts stay exact to 2^53 (a 12-leaf star's 12! is exact);
     a position law sums its f g terms relative to the largest term's binary exponent, since the two scaled
     factors can underflow as a product where the term still counts. Down-sets live in an open-addressing
     table of 32-bit indices, so a layer holds fewer than 2^32 of them.
   - c1's choice is inert in this ratio (deviation from Counting: algorithm, found in implementation): with
     the pair in one component theta is e(C0) / e(C*), and the merge is symmetric in the pair; the pair sits
     in two components only on a free split, where c1 is the left child under either rule. Taking c1 by code
     leaves every tests/cpp check passing. The order matters only for a theta taken as an adjacency (the
     brute-force oracle, a coin), where by code gives 0 on a decreasing split, and tests/cpp shows that.
   - Measured (arm64 macOS, same machine): the synthetic 25-leaf star 4.6 s against the prototype kernel's
     10.4 s, the 60-leaf eight-chain order 4.9 s against 10.0 s (351 MB against 404 MB peak); the star with
     every layer kept 4.6 s and 524 MB. A 5-tree, 2-constrained + 1-free fit's 2,437 moves (102 with the pair
     in two components, up to 7.8e4 down-sets): 3.5 s in all, worst 14 ms, the prototype 4.4 s, theta equal to
     print precision. The 1-tree, 3-free-axis fit (p 4, seed 1) on the stage 1 engine's draws puts 45-54
     leaves in every move's U: the first has 8.2e8 down-sets, and neither the counter nor the prototype
     finished it in 10 minutes (checkpoint).
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
   - As implemented (commit 4): the seam is `NormalizedLeafModel`, split in two. `prepareLogNormalizerRatio`
     builds T*'s order while T* is in place (after `tree.birth`, before `orphanChildren`) and returns log m;
     `logNormalizerRatio` counts off that order alone, never the tree, and runs only when the bound leaves the
     decision open (`decideNormalizedMove`). A death's count therefore runs after `orphanChildren`, since r1
     needs the merged leaf's score: a cancel there restores the node as the reject path does (step 15). The
     count is held to log m, so rounding never lets the bound and the count disagree for one u. A death whose
     merged cone is empty (r1 = 0) rejects without counting, and a pair of two isolated leaves gives 0 without
     a count. The census writes a 'z' record per normalized move (m, the pair's component count, whether the
     count was needed and ran, log Z_T0 - log Z_T*, time, down-sets, W(U), peak bytes, whether B = 2^22 would
     switch it, a_B / a_MH, and on an accepted death the merged component's down-sets counted up to 2^20);
     BARTCORE_MOVE_CENSUS_COUNT_ALL makes it count moves the bound decided too, the decision unchanged. The
     'p' record's logPrior carries the Z term, NA when no count ran.
3. Drop the d terms: `priorMass` in `oneLeafLogMarginal`, `denom` in `twoLeafCoupledLogMarginal`. Under "joint"
   this alone is the exact move.
4. Exact pair redraw in `redrawAfterBirth`:
   - Keep the capped rejection for the upper leaf; it is exact whenever it accepts.
   - When the cap is reached or acceptMax underflows, draw the upper leaf by inverting the CDF of its marginal,
     phi_U(u) [Phi_L(min(bL, u)) - Phi_L(aL)] on [max(aL, aU), bU]. Use the `coneProbability` quadrature's
     cumulative with safeguarded Newton, evaluate the density exactly, and work on the log scale in the tail.
   - Then draw the lower leaf on [aL, min(bL, mu_U)], which is never empty.
   - As implemented: the marginal is log-concave, so the fallback finds its mode, cuts its support where the
     log density is 50 below the peak, and inverts a per-side composite adaptive Simpson cumulative with
     safeguarded Newton, all relative to the peak (`monotoneInvertLogConcave`), rather than reusing
     `coneProbability`'s fixed window. About 70-110 us per fallback draw. `normalMass` takes upper tails above
     the mean, so the rejection's admitted mass does not cancel there.
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
   - A forced predictor update, setCutPoints and setData leave no empty leaf but can leave a tree outside
     the cone through a collapse or a remap (Context); step 7's reseed, run after each, puts it back.
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
     As implemented: the bridge compares against the engine's own directions and prior, which SamplerShape now
     carries, so a model all-zero in directions counts as unconstrained.
   - Commit 1. growForestFromRoot reseeds mu to the all-zero feasible seed before its draw.
   - Commit 1. Reseed, then validate (dec-A128). [`Chain::installForest`](../../src/bartcore/chain.hpp) reaches
     the trees through [`Chain::rebuildLiveForest`](../../src/bartcore/chain.hpp) or
     `rebuildLiveForestRemapped`. On the installTrees path, every installed tree that fails
     `monotoneTreeIsFeasible` has all its leaves set to 0, the all-zero feasible seed growForestFromRoot uses
     (equal values satisfy every constraint, and no RNG is drawn); feasible trees keep the donor's values. The
     next sweep's leaf Gibbs step redraws them.
   - Commit 1, found in implementation (orchestrator call). The same reseed runs after a forced predictor
     update, setCutPoints and setData ([`Chain::forceRefreshTrees`](../../src/bartcore/chain.hpp),
     [`Chain::applyNewData`](../../src/bartcore/chain.hpp)): a collapse or a remap can relate merged leaves
     to new neighbours, not only strand empty ones. Any feasible state is a valid continuation point, and the
     tree is reseeded to all-zero with no RNG draw.
   - Commit 1. setState, copy and reload are refused up front: Sampler::setState checks a new
     monotoneStateFeasible predicate on every chain's live trees before any chain is touched, beside
     interactionStateFeasible and columnMaskStateFeasible, so a refusal leaves the sampler exactly as it was, as
     setState promises. The refusal is named: the state's leaf values violate the monotone constraint. With
     steps 4, 5 and the reseeds every state the sampler produces passes this predicate, but not setState
     overall: installTrees on the same grid can leave member-empty leaves, which setState's validity check
     refuses, so a copy or reload fails until the moves clear them. That is pre-existing and holds for
     unconstrained samplers too; it is queued separately.
   - Commit 4. [`Chain::sampleTreesFromPrior`](../../src/bartcore/chain.hpp) under "leaf" needs no new predicate: the
     prior is not restricted. Under "joint" it draws each tree jointly by rejection, a CGM tree and iid
     unconstrained leaves kept only when the leaves lie in its cone (acceptance is the prior mean of Z_T), with
     no count, and then discards the leaves, as its contract returns trees without leaf values.
     As implemented: the cone test joins the empty-leaf test in the one rejection loop and shares its attempt
     cap, which still bounds acceptance below by 1 - base, since the bare root always lands in its cone.
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
      different cuts, giving 26% N mass in the tested x1-rooted group; cF (stage 3b), c3's cells with x2 a
      free two-level factor with no missing values, x1 rising at one level and falling at the other. A
      two-level factor's one partition is the enumeration's one cut, and with no missing value its labelling
      does not change the order, so the gate reads the engine's subset rule without enumerating level sets;
      wider factors and missing values rest on tests/cpp's point oracle and the tinytest.
    - The script replaces the planned part (c) of monotone-reference.R. It joins exact-gates.yaml's list in the
      fix commit, since it fails the current engine by design.
    - The general DP beyond these sizes rests on step 8's brute-force checks.
    - Both targets: engine and prototype modes take the prior, weighting each tree with log Z_T in its weight
      ("leaf") or without it ("joint"), and every group must pass under each. Its unnormalized mode already
      prints the two exact laws side by side. Every run names its prior.
    - A mirrored design with x1 decreasing (cN's cells reflected along x1), guarding the move's and the
      redraws' handling of a decreasing axis. It cannot see the counter's: a c1 taken by code is inert in
      the ratio (step 1), and an order and its reverse have the same count, so one constrained axis's
      direction does not change e. tests/cpp's mixed-direction checks cover the counter.
    - Runtime: quick mode measured 8 min 15 s per prior over four designs, so with the mirrored design about
      10.5 min per prior and 21 min for both. The monotone gate runs as its own CI job with a 40-min timeout,
      and exact-gates.yaml keeps its job and timeout.
    - As implemented (commit 4): the mirrored design is cM. With no prior named the script runs both, each
      run naming its prior; `leaf` or `joint` on the command line runs one. The CI job is a second job in
      exact-gates.yaml and always runs quick: full mode triples the draws and would pass the timeout. Quick,
      both priors, six designs: about 21 minutes on arm64 macOS.
11. Two monotone scenarios in benchmarks/R/equivalence.R (x1 and x2 constrained, 20 trees, one per prior): a
    55-scenario re-record, the other 53 bitwise, and MANIFEST rows naming the enumeration gate as their ORACLE
    (P17).
    - As implemented (stage 6): monotoneleaf and monotonejoint share one design, n 400 and p 5, x1 increasing and
      x2 decreasing (the surface decreases in x2, so the decreasing axis is exercised), through the sampler API.
      The new baseline is named after the newest package-code commit, as its predecessors are; one MANIFEST row
      covers both scenarios. The 53 others reproduce the previous baseline bitwise, and the row states that the
      enumeration gate reaches one tree on small grids only, given step 16's flag.
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
    - As implemented (commit 4): the prior rides a second model attribute, monotone.prior ("leaf" or "joint"),
      beside the monotone directions attribute rather than inside it, so the directions stay a plain integer
      vector for the callers that read them. The bridge requires it wherever directions are present and holds
      no default of its own. copy() and reload rebuild the sampler from the model, so both carry the prior;
      setState installs trees and leaves only, so the sampler keeps its own. The fit carries monotone.prior,
      absent without a constraint; print shows "monotone prior:" and summary "Monotone prior:".
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
    - As implemented (commit 2): the one constant is MONOTONE_PRIORS, the prior values with the default first,
      installed as monotone()'s prior formal, so match.arg and the shorthand both read it. The direction
      vocabulary also takes "1", "-1" and "0" as strings, the values c() makes of the codes when words share
      the vector ("0" was already accepted); a number must be exactly -1, 0 or 1 (0.6 no longer rounds to 1),
      and a logical is refused. prior is matched at construction; "joint" is refused at fit time, in
      resolveMonotone, whose prior goes no further than R until commit 4.
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
      its reject path already takes) or restores the dying node (the count reads T*'s order, built before
      `orphanChildren`, but runs after it; step 2), and rethrows. Chain::run catches both, rebuilds (below), and returns cancelled for CountCancelled or
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
    - As implemented (commit 5): the poll sits in the forward DP, the backward counts and the position laws, one
      per down-set, every 2^16 by default. Only totalFits is rebuilt, re-summed from the trees in tree order: the
      running residual needs nothing, since each sweep's first roll rewrites it whole from totalFits, and the test
      totals are rebuilt by every recorded sweep before they are read. A move's allocation failure is rethrown as
      bartcore::CountOutOfMemory, still a bad_alloc, whose message names the remedies. sampleNodeParametersFromPrior
      rebuilds totalFits the same way before its allocation failure goes on. The tally and the threshold live on the
      monotone leaf; the threshold, the poll interval and a one-shot injected allocation failure are process-wide
      test knobs (monotoneCountHooks), R reaching the first and third through an unexported bridge entry. The
      bridge attaches the tally as a "slow.count" attribute; a burn-only run, which returns NULL, returns an empty
      list carrying it when a count was slow. bart()'s burn-in and kept runs merge their tallies and warn once.
      The burn-only and per-sweep-callback run entries now capture engine exceptions as the main one does.
16. SBC arms (benchmarks/R/sbc.R has none today; ~120 lines), as
    [Monotone arm: design](sbc-family-tiers.md#monotone-arm-design) lays them out: the burn-monotone run, the
    monotone-1 and 20-tree arms once per prior, each naming its prior, and the unconstrained monotone-bd twin
    once. monotone-1 at 0.4 s per replicate, the 20-tree arm ~85 min at R 200.
    - As implemented (stage 6): arms monotone-leaf, monotone-joint (20 trees, n 150, p 3), monotone-1-leaf,
      monotone-1-joint and monotone-bd (1 tree, p 1, the design's mini SBC; n 100, then n 20 below), and the ladders
      burn-monotone-leaf and burn-monotone-joint. Each runs the gaussian replication as a family-spec arm, so
      every rank takes sbcDiscreteRank's tie-break (mono.local has an atom at 0) and the ladder applies. The
      20-tree arms read their band at 0.05 / (57 + 20), both arms joining the matrix, not 0.05 / (57 + 10).
    - The design's one-tree setting (thin 10, 1000 burn sweeps) is too short: there the twin flags sigma
      (ecdfDiff 0.101 against a band of 0.066) and f.star3 (0.075). At thin 50 and 5000 burn sweeps it passes
      all nine functionals, so the one-tree arms run at thin 50.
    - First result, R 400, L 100, thin 50, n 100, on arm64 macOS: monotone-1 flagged under both priors ("leaf" eight
      of nine functionals, mono.wide 0.181 against a band of 0.066; "joint" four), ranks piling at 0: the posterior
      contrast too steep. The twin passed.
    - Finding (dec-A133): no engine defect; a one-tree birth/death chain does not mix its structure on an
      informative design, the twin included.
      - The kernel is exact. A successive-conditional check (Geweke 2004: draw theta0 from the prior, simulate y,
        run K sweeps from theta0, compare paired functionals) shows no drift at K 1 (100,000 replications), 20
        (20,000) and 200 (21,000) under both priors and for the twin, every |z| under 2. The SBC with the tree
        held at the generating one (the all-zero mixture) passes every functional under both priors: the leaf
        Gibbs sweep and the sigma step are exact given T.
      - The chain is trapped. On one dataset, a chain started at the truth stays near the true leaf count while
        one started from a prior draw, or from the root, holds 4 to 9.5 leaves against 2 to 4.4 for 20,000
        sweeps. The twin traps the same way: its leaf count fails SBC (rank z -13) while its f functionals pass,
        since extra leaves at one level do not bias an unconstrained fit. Under the constraint they do: the
        ordered, truncated leaves of a flat stretch spread apart, so f steepens and sigma grows, the flagged
        pattern under both priors, which is why "joint", which never counts, flags as well.
      - Larger enumerations measure mixing, not exactness: a six-cell one-predictor design with a step and
        strong data fails the enumeration gate under "leaf" (T2 up to 1.2e6), and so does its unconstrained
        birth/death twin (T2 up to 9.7e5), each chain holding an over-split tree whose exact conditional mass is
        0.06 at 0.88. A flat five-cell design passes.
    - So monotone-1 is a mixing diagnostic, not a pass requirement (dec-A133), at n 20, where the likelihood is
      weak enough to let the chain move, with the leaf count among its functionals. Result, R 400, L 100, thin
      50: "leaf" passes the f functionals, sigma and mono.wide (0.060) and flags mono.local (0.107) and the leaf
      count (0.268); "joint" flags mono.wide (0.093) and the leaf count (0.230); the twin flags only the leaf
      count (0.193).
    - Exactness at depth is benchmarks/R/monotone-successive-conditional.R's: one tree, n 100, tree prior power
      0.5 so most moves touch a leaf with a frozen constrained neighbor, K 20, under both priors and the twin.
      It needs no mixing, since the chain starts at a posterior draw. Quick mode (6,000 replications, about a
      minute) passes at worst |z| 2.1 and fails the old ratio (the d divisions restored, the Z term dropped) at
      |z| 9.0 under "leaf" and 31 under "joint" on the leaf count; at the default tree prior the old ratio read
      only |z| 2.3. It runs in the monotone CI job, per prior. Full mode (20,000 replications, about 2.5 minutes) passes at worst |z| 1.9.
    - The 20-tree arms: the ladders (40000 sweeps x 3 datasets) put the transient in the first two or three
      4000-sweep blocks and the slowest functionals past ACF 0.1 at lag ~200 ("leaf") and ~140 ("joint"), so
      the arms run at 12000 burn sweeps and thin 100. Result, R 200, L 150, band 0.137 at 0.05 / (57 + 20), on arm64 macOS:
      both PASS every functional. "leaf" worst ecdfDiff 0.075 (mono.local), mono.wide 0.073, sigma 0.042;
      "joint" worst 0.082 (f.star1), mono.wide 0.036, sigma 0.045. About 31 s per replicate of CPU at one
      thread, under unrelated load.

## Verification

- `cd tests/cpp && make && ./test_bartcore`: the count, ratio, lazy-count, scale, redraw, extension-draw,
  slow-count, interrupt, allocation and "joint" checks pass.
- Stage 3b: tests/cpp's point oracle routes every combination of the split variables' values, each level
  and the missing value where the column has one, through random trees over a numeric store and one with
  a 4-level and a pooled 70-level factor and missing values; the order builder's relation must equal the
  pairs of points
  one code apart along a constrained predictor, and the bounds and the feasibility check must read them.
  The count and ratio brute-force checks also run on that store. A hand-built tree relates two leaves
  through the pooled factor's missing value alone, and one level partition labelled two ways must give one
  order. Restoring the cut-interval reading of a subset axis, ignoring where missing values go, or moving
  the pooled missing position off bit K fails them. An unforced and a row-by-row predictor update that
  bring a factor its first missing value must reseed a tree the new value takes out of the cone. test-monotone.R checks that one-tree
  fits are monotone along x1 at every level of a free factor and at a free predictor's missing value; the
  engine before 3b fails both.
- The checkpoint (Staging) is reported before commits 5 and 6.
- `R CMD INSTALL --preclean -l <lib> .`, then `R_LIBS=<lib> Rscript -e 'tinytest::test_package("dbarts")'`.
- `R_LIBS=<lib> Rscript benchmarks/R/monotone-exact-enumeration.R quick`, under each prior, the mirrored
  decreasing design included: every group passes. This takes about 8 min per prior; full mode runs 900k draws.
  The mutation run restores the d divisions and drops the Z term (then `touch` the header and reinstall), and
  must fail like the current engine. (The planned second mutation, c1 by code instead of by order, is inert in
  the engine's ratio, step 1.)
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
- `R_LIBS=<lib> Rscript benchmarks/R/monotone-successive-conditional.R quick`: every arm passes; with the d
  divisions restored and the Z term dropped, "leaf" and "joint" fail.
- Release level: step 16's 20-tree SBC arms, once per prior, must pass before admission. monotone-1 is a
  mixing diagnostic (dec-A133): its leaf count flags for the twin too, and it gates nothing.
- Speed: on a quiet machine, under each prior, monotone sweep time at 20 trees, 1 and 2 constrained predictors,
  within 5% of today; at 1 and 5 trees under "leaf", the slowdown is recorded against the Decision's estimates.
  bench-sampler compare unchanged on the unconstrained paths.
- The feel study (Default: feel study) runs once both priors pass the gates above, before the default is
  ruled.
- `Rscript tools/check-doc-freshness.R .` passes.

## Landing

Stage 1, 2026-09-30: the redraw fix, empty leaves and reachability (a0e0100a), the setControl fix that
stored a control before the engine accepted it (565dcd2b), and the truncated-normal primitive reporting
a stall as NaN with an exact narrow-tail proposal (4f7a5f02). The monotone draw reflects intervals above
the mean. Reviewed twice by an independent reader; tests/cpp with ASan/UBSan clean, tinytest 10,756/0,
the lint chain clean, all 25 exact gates quick, monotone-reference.R quick; on the reference build the
three equivalence compares bitwise 53/15/11 and the snapshots unchanged, so no recorded draw moved;
ordinal-exact quick byte-identical; R CMD check --as-cran one NOTE (Date); stan4bart 491/491, bartCause
0 failures. The enumeration gate still fails as the plan expects until stage 4. Mutations of the
fallback, the empty-leaf pin, the up-front setState check, the growForestFromRoot reseed, the collapse
reseed and the reflection each fail the new tests. Calls made while implementing: dec-A129.

Stage 2, 2026-09-30: the monotone() constructor, the direction vocabulary and the caller sweep
(bcab58b2), and the review follow-ups: a partly named or doubly named vector and a non-string prior
refused, bart.Rd's bare-name note (dc15a41b). prior = "joint" is refused until stage 4; the default
constant MONOTONE_PRIORS holds "leaf" first, provisionally. Reviewed by an independent reader; tinytest
10,795/0, the lint chain clean, pkgdown clean, R CMD check --as-cran one NOTE (Date); on the reference
build the three equivalence compares bitwise 53/15/11 and the snapshots unchanged; stan4bart 566/566,
bartCause 0 failures. Mutations restoring case-insensitive matching and the sign glyphs fail the tests.
Calls made while implementing: dec-A130.

Stage 3, 2026-09-30: the leaf-order builder, the down-set count, log Z_T, the move ratio on the finer
tree's side, the position laws, the linear-extension draw and the exact prior leaf draw, none called yet
(e33a79cb). Reviewed by an independent reader: count, ratio, draws and scaling verified; mutations of the
position law and the extension weights fail the tests. tests/cpp with ASan/UBSan clean, tinytest
10,795/0, the lint chain clean, R CMD check --as-cran one NOTE (Date); on the reference build the three
equivalence compares bitwise 53/15/11. The engine counts the 25-leaf star (1.68e7 down-sets) in 4.6 s,
about twice the prototype's speed. The review found the geometry defect that stage 3b now fixes, and the
1-tree p 4 fit re-run on stage 1's draws needs 8.2e8 down-sets on its first move, which neither the engine
nor the prototype counted in 10 minutes: the checkpoint's trigger is likely to fire.

Stage 3b, 2026-09-30: the leaf geometry reads level sets and missing-value routing (08613432 design
note, e3d11f53, d9d338db). A factor split and a free predictor's missing values now relate leaves as the
constraint requires: the fitted surface is monotone at every factor level and at the missing value,
where the stage-3 engine decreased in up to every draw. An axis has a missing position only when its
training column has missing values, as tree routing does; the first design gave every factor axis one,
which made the order depend on how a level split was labelled. A tree that an update's new missing
values make infeasible is reseeded, on the unforced update paths too. Reviewed twice by an independent
reader: the point oracle, which routes points through the tree's own rules, agrees with the relation on
800 random trees (229 pairs through missing values); numeric fits, and fits whose only missing values are
in one constrained predictor, draw bitwise as before; tests/cpp with ASan/UBSan clean, tinytest 10,809/0,
the lint chain clean, R CMD check --as-cran one NOTE (Date); reference-build compares bitwise 53/15/11;
the exact gates quick; stan4bart 566/0, bartCause 1145/0. Calls made: dec-A131.


Checkpoint, 2026-09-30: stage 4 (08d978a1) built with the move census, prior = "leaf", seed 1, n 5000,
1000 burn-in and 200 kept sweeps, at most three fits at once on arm64 macOS; per fit in
[monotone-checkpoint-08d978a1.csv](../../benchmarks/baselines/monotone-checkpoint-08d978a1.csv). Count
time and memory are the largest per-move count the sampler needed; "sweeps > 1 s" is the share of sweeps
with such a count over 1 s; "big merges" is accepted deaths leaving a merged component past 2^20
down-sets, out of all accepted deaths.

| fit | trees | free axes | sweeps completed | max count | max memory | sweeps > 1 s | big merges | capped |
|---|---|---|---|---|---|---|---|---|
| 1 constrained + 3 free (p 4) | 1 | 3 | 882 of 1200 | 70.1 s | 361 MB | 20% (178) | 13 of 32 | stopped at 36 min |
| 1 + 2 | 1 | 2 | 1200 | 0.06 s | 1.0 MB | 0 | 0 of 20 | no (4 s) |
| 1 + 1 | 1 | 1 | 1200 | 23.3 s | 143 MB | 9.8% (118) | 8 of 46 | no (11 min) |
| 2 + 1 | 1 | 1 | 1200 | 4.4 s | 41 MB | 7.3% (88) | 5 of 36 | no (3.3 min) |
| all four designs | 5 | 1-3 | 1200 | 0.09 s | 2.1 MB | 0 | 0 | no (2-16 s) |
| all four designs | 10 | 1-3 | 1200 | 0.5 ms | 0.03 MB | 0 | 0 | no (2 s) |
| all four designs | 20 | 1-3 | 1200 | 0.7 ms | 0.02 MB | 0 | 0 | no (3 s) |

- The p 4 fit was stopped by hand once settled: its 178 slow sweeps already exceed 10% of 1200, six
  counts took over 60 s (66-70 s), and its last ten minutes ran about 4 sweeps a minute with 318 sweeps
  left in 24 minutes, so it would have hit the 60-minute cap. A first run, stopped at 863 sweeps, drew
  the same moves. Its slowest counts held 1.5e8 down-sets; no count came near 4 GB, and the process
  peaked at 618 MB resident. In the 1 + 1 fit every slow count falls after sweep 945.
- Hybrid (the same fits with every move counted, 15 minutes each): the share of counted moves the
  hybrid would switch at B = 2^22, and a_B / a_MH on those moves (mean, median, minimum; moves the free
  bound decided carry no ratio and are left out). One tree: p 4 32.5% (794 sweeps; 0.960, 0.996,
  0.559), 1 + 2 10.0% (0.976, 0.996, 0.596), 1 + 1 45.6% (1176 sweeps; 0.966, 0.995, 0.513), 2 + 1 28.7%
  (0.964, 0.997, 0.508). Five trees: 1 + 1 0.3% (0.967, 0.989, 0.736), 2 + 1 3.2% (0.953, 0.988, 0.544),
  the others none. Ten and twenty trees: none. Every fit not named with a sweep count completed its 1200
  sweeps. Switched moves would keep 95-98% of Metropolis-Hastings' acceptance on average, and the
  median move nearly all of it.
- Trigger: fires. Conditions that fired: a 1-5 tree fit with a count over 1 s in more than 10% of its
  sweeps (1-tree p 4, at least 15% of 1200); a count over 60 s (the same fit); and a fit on course for
  the time cap (the same fit, stopped before it). Not fired: no fit at 10 or more trees has a count over
  1 s (the largest is 0.7 ms), no count reached 4 GB, and the 1-tree 1 + 1 fit sits just under 10%
  (9.8%). The hybrid goes to the maintainer for a before-release call.
- Run times: timing 10:59-11:23 and 12:10-12:46 (the p 4 fit re-run after the first was stopped early),
  hybrid 11:15-11:22 and 12:10-12:25; runs that spanned a machine sleep were discarded and re-run.

Stages 4 and 5, 2026-09-30, pushed together: the exact move under both priors, the prior flag from
monotone(prior = ) to the engine, the exact prior leaf draw and the joint prior's structure draw
(720bf389, 08d978a1), the checkpoint driver's fixes (0591b67c) and results (15613958, above); then the
order count polling the chain's cancel function, a cancel or allocation failure restoring the move's
tree and rebuilding the fits, the slow-count warning, allocation failures as R errors on every run
entry, and setModel refusing a changed constraint or prior (3f7177ea, 90074738, 28dad8bb, 7a72aeb7,
a50f6afb). Each stage reviewed by an independent reader. The enumeration gate passes every design under
both priors (p 0.24-0.99) and fails with the old ratio restored (p <= 1e-10); tests/cpp clean under
ASan/UBSan and TSan; tinytest 10,852/0; the lint chain clean; R CMD check --as-cran one NOTE (Date);
reference-build compares bitwise 53/15/11 and snapshots unchanged; the 25 exact gates quick; stan4bart
491/0, bartCause 0 failures. A real interrupt mid-count leaves a state equal to a rebuild of the same
trees and takes effect in 33-47 ms. Calls made: dec-A132.

Stage 6 (all but step 12's docs) and the feel study, 2026-09-30: the two monotone equivalence scenarios,
one per prior, and equivalence-90074738 (55 scenarios, the other 53 bitwise as before; d80f27c7,
4ff51229); the SBC arms (4d3c4041, f354dc6f, 90de0ba3); the successive-conditional check, which runs in
the monotone CI job and fails the old move at |z| 9 under "leaf" and 31 under "joint" (d9807526,
61bcc779); the one-tree SBC arm demoted to a mixing diagnostic after its flag was traced to birth/death
mixing shared with the unconstrained twin, not to the kernel (dec-A133; 2ae44594, 2819d47a); and the feel
study (7178b075, 922035df, 95eb4113). Both 20-tree SBC arms pass every functional (worst ecdf distance
0.075 "leaf", 0.082 "joint", band 0.137). Lint chain clean; the equivalence file reproduces 55/55 bitwise
on the reference build. Step 12's docs wait for the maintainer's default ruling.

