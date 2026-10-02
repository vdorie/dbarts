# Review 3 - verification of multiforest.md

Tree .claude/worktrees/review3 at 01dee4b4, library r3-lib. Own probes: scratchpad r3-verify-mf-01.R ..
r3-verify-mf-05.R (new data and seeds, not the finder's scripts). Checked against TODO, docs/decisions.md
(through dec-B175), feature-matrix.md, multinomial.md, interaction-constraints.md, linear-leaves.md,
gp-leaves.md, heteroscedastic.md. None of the five is filed or ruled.

| id | verdict | severity | fix kind | public surface | recorded draws |
|---|---|---|---|---|---|
| 01 | CONFIRMED | BLOCKER (silently wrong model) | refusal or implementation (maintainer's choice) | yes, either way | no |
| 02 | CONFIRMED | BLOCKER by the rubric; narrow reach | implementation (state block) | state object gains a block | no |
| 03 | CONFIRMED | MAJOR | up-front gate plus bridge-03's atomicity | no (message only) | no |
| 04 | CONFIRMED, plus plot() | MINOR | docs, or drop the element (maintainer) | yes if dropped | no |
| 05 | CONFIRMED (all four bullets) | MINOR | refusals and messages | refusals only | no |

None is a regression against 0.9-34: multinomial, interactions()/blocks(), linear/gp leaves, variance
forests, installTrees and forest bases are all new in 1.0-0.

---

## multiforest-01 - CONFIRMED, BLOCKER

Code. Chain::buildMultinomialForest (src/bartcore/chain.hpp) sets tree counts, the proposal mix, k and the leaf
scale, and nothing else. MultinomialForestSpec (combiner.hpp) has no interaction or block fields. Compare
ForestStructureSpec and buildSpecifiedForest, which install both. buildMultinomialSampler
(R_interface_bartcore.cpp) copies only base, power and the proposal probabilities from the parsed model.
model.interactionMaxOrder, the forbidden pairs, blockOfColumn and blockTreeCounts are parsed, then dropped.
R's `unsupportedMultinomial` (R/spec.R) names DART, split.probs, monotone, linear/gp, a k hyperprior, a
named sd and single storage. It omits these two.

Intent. multinomial.md: "Unsupported surface is refused by name, never silently reshaped". interactions.Rd:
"hard (an availability ban), applied per forest". No document says multinomial is meant to ignore them. By
git history, interactions()/blocks() landed 2026-08-20 and the public multinomial surface on 2026-08-24.
The refusal list was written from what MultinomialForestSpec lacks, and these two were missed. This was an
oversight, not a design choice.

Own probe (r3-verify-mf-01.R): n = 300, 4 predictors, 3 categories, 6 trees, 100 sweeps, seed 1.
```
multinomial max.order=1: trees with >1 distinct var: 6 of 18
gaussian control max.order=1: trees with >1 distinct var: 0 of 6
multinomial blocks: trees mixing blocks: 2 of 18
bart() multinomial forbid accepted: TRUE
```

The two options:

1. Apply the constraint to every category forest. MultinomialForestSpec gains the same six fields
   ForestStructureSpec carries (or embeds a ForestStructureSpec). The bridge copies them from the parsed
   model. buildMultinomialForest installs them through a helper factored out of buildSpecifiedForest's
   interaction and installBlockMasks block. About 40 lines of C++. The setState and installTrees containment
   checks (interactionStateFeasible, columnMaskStateFeasible) already loop over every forest, so the
   multinomial gets them for free. The softmax is symmetric in its K forests, so applying one constraint to
   all K raises no identifiability question.
   Costs: interactions.Rd and blocks.Rd gain a sentence. Trust in the result needs a constrained arm in the
   multinomial exactness gate (multinomial-exact.R), or a Geweke arm like the finder's BCF ones. That is
   the larger part of the work.
   Draws: an unconstrained fit is byte-identical, since the fields default off. The multinomial-equivalence
   baseline does not move.
2. Refuse. Add "interactions()" and "blocks()" to `unsupportedMultinomial`, with a backstop in
   buildMultinomialSampler that errors when any of the four parsed fields is set. About 10 lines plus two
   expect_error tests. One sentence in interactions.Rd and blocks.Rd, and in multinomial.md's refusal list.
   A user loses a combination that never worked.

Recommendation: option 1, unless the exactness arm cannot be built before release, in which case refuse
now (option 2) and file option 1 in TODO. The engine work is plumbing over code that is already generic per
forest. The finder's Geweke checks of the identical per-forest install on BCF passed. What a user wants from
"interactions = max.order 1" on a multinomial fit is clear: each category's log-odds additive.

Tests (either option): a multinomial fit with max.order = 1, forbid and blocks. For option 1, assert that no
tree in any category forest violates the constraint after run and growFromRoot, and that setState and
installTrees of a violating donor are refused. For option 2, expect_error naming the argument, on both
dbarts() and bart().

Why gates missed: test-interactions.R and test-blocks.R never use family = "multinomial", and the
multinomial tests never pass a constraint.

## multiforest-02 - CONFIRMED, BLOCKER by the rubric (silently wrong predictions), narrow reach

Code. LinearGaussianLeaf and GPGaussianLeaf (src/bartcore/model.hpp) hold means_ and sds_. GP also holds
lengthscales_, which come from a median-pairwise-distance heuristic over the standardized training values
when not supplied. regatherTrainingCovariates keeps all of them on setPredictor, by design (linear-leaves.md
"Mutation semantics"; gp-leaves.md). reinitialize recomputes them from the data. ChainStateData and
ForestStateData carry none of them. Re-creation (copy(), readRDS, a dead pointer) builds a fresh sampler
over data@x, which setPredictor has already updated, and then setState. The saved slopes and GP function
values are therefore read under different constants. The cut grid is the same kind of calibration, and it
does ride the state (Sampler::setState installs state.cutPoints), so the omission of these constants is
inconsistent.

Own probe (r3-verify-mf-02.R): n = 100, linear(c("x1","x2")) or gp(..., max.leaf.size = 100), 5 trees, keepTrees,
run(40, 4). Then setPredictor(x2^3 with ends pinned to 0 and 1, column 2), which returned TRUE, then run(0, 2). sd(y) = 0.95.
```
linear: live pred-vs-recorded 2.2e-15   copy 1.4    readRDS 1.4    next sweep live-vs-copy 0.207
gp:     live pred-vs-recorded 0.0026    copy 2.28   readRDS 2.28   next sweep live-vs-copy 0.214
constant-leaf control, same mutation: copy pred-vs-recorded 1.8e-15
```
Side note: a setPredictor whose values would leave a leaf empty returns FALSE and changes nothing. A probe
that does not check the return value sees the old column still in data@x.

Reach: dbarts() sampler users who call setPredictor on a designated leaf-covariate column, then copy, save or
lose the pointer. bart() fits never mutate. The gp continuation gap without mutation is documented and is
separate.

Fix sketch:
- Engine. ForestStateData gains `leafCovariateCenters`, `leafCovariateScales` and, for gp,
  `leafLengthscales`. These are per forest because the leaf is per forest, and empty means absent. The
  capture step fills them from forest.leaf. Chain::setState, when they are present, assigns them and then
  regathers the training and test covariates under them. GP also clears its kernel caches, as
  regatherTrainingCovariates does. stateIsValid requires length q and finite positive scales, and finite
  positive lengthscales.
- installForests must not copy them, since a warm start reinterprets the donor on the recipient's data, as
  setData does. It builds dst field by field, so the work is a comment and a test.
- Bridge. Append forest slots ("leaf.covariate.center", "leaf.covariate.scale", "leaf.lengthscales") to
  forestSlotNames. The registry is append-only and read by name, so there is no formatVersion bump and no
  change to minReadableStateFormatVersion.
- Versions. A state without the block (anything saved earlier in 1.0 development; none were released) loads
  as today and recomputes from its data, which is correct for every sampler that never mutated a designated
  column. Requiring the block for linear/gp states is possible, since nothing released carries one, but it
  buys nothing.
- Rejected alternative: re-deriving the constants at setPredictor, as setData does. It reverses the
  documented sticky calibration. It would move the prior and reinterpret every live slope on each in-place
  update, which is exactly the imputation-inside-Gibbs use the sampler exists for. It also moves draws for
  anyone using that path today.

Draws: no recorded baseline moves. The default path is untouched, and states without a mutation carry the
same constants the recipient would recompute.
Public surface: the state object (storeState/$state) gains named blocks, additively.
Tests: extend test-mutate-then-serialize.R with the copy and readRDS route, the one where data@x already holds
the mutated column. Cover linear and gp, and assert that the copy's predict(xn) matches the live one (bitwise
for linear, to the documented nugget tolerance for gp), and for linear that the next run is identical given
the carried RNG state. Add a hand-edited state with a zero scale, which must be refused.
Why gates missed: as the finder says, the existing test restores into a cold sampler over the original data and
replays the mutation.

## multiforest-03 - CONFIRMED, MAJOR

Code. Sampler::installForests' per-chain gate checks `src.varianceTrees.empty() == hasVarianceForest()`,
which is presence only. The slot path checks the saved buffer against the donor's own count (nvt), never
against this sampler's. The count check `trees.size() != vf.numTrees` sits inside Chain::installVarianceForest,
which runs after the mean-forest commit loop.

Own probe (r3-verify-mf-03.R): donor variance n.trees 4, keepTrees, run(50, 3). Target n.trees 5, run(30, 1).
installTrees(D, samples = 2).
```
chains 1 : warm-start donor's variance trees cannot be installed on this sampler's data (a rebuilt variance tree leaves a ...
  mean fits changed: 0.181  variance changed: 0   usable after: TRUE
chains 2 : same message; mean fits changed: 0.181  variance changed: 0
mean-tree-count mismatch (control): "not shape-compatible ..."  changed: 0
```

Does bridge-03 cover it? Only half of it. bridge-verify.md already notes that the variance loop shares the
partial-commit defect. A scratch-then-commit fix that includes installVarianceForest makes this refusal
atomic. It would still give the wrong reason: the message blames empty leaves or non-positive scales for a
plain count mismatch.
The cross-grid arm of installVarianceForest is also non-atomic within the variance forest: it initializes
tree j before buildFromFlat can fail. The same snapshot/commit covers that too.

Fix: in the per-chain gate, `if (hasVarianceForest() && src.varianceTrees.size() !=
chains_[c]->numVarianceTrees()) return shapeMismatch`. A dedicated result would be better, so the message can
name "variance trees". This also bounds nvt on the slot path. Land it with bridge-03's atomic install.
Test: in test-heteroscedastic-warm-start.R, install a donor with 4 variance trees into a sampler with 5.
expect_error naming the variance tree count, and getForestFits()/getVariance() identical after the refusal.
Draws: none move.

## multiforest-04 - CONFIRMED, MINOR, one more symptom

Own probe (r3-verify-mf-04.R): y = 10 x1 + noise with sd 0.2 + x2. fit$sigma and fit$first.sigma are the single
value 10.7705, which is diff(range(y)). run()$sigma and getSigmas() give the same. mean(s.train) is 0.62.
extract(fit, "sigma") correctly refuses. summary() omits sigma.
New: plot(fit) draws the sigma trace panel for a heteroscedastic fit. fitHasResidual is TRUE for gaussian, and
plotSigmaTrace is called with a constant series at range(y), a meaningless flat line.

Intent. bart.Rd's `type` paragraph states the exception ("a fixed unit residual times the range of the
response and is not its residual scale"). bartBT.Rd's Value `sigma` and dbartsSampler-class.Rd's getSigmas
do not. The statement exists, but in the wrong place.

Fix options (a maintainer call, since it touches what a fit carries):
(a) Drop sigma and first.sigma on a heteroscedastic fit, so they are absent (NULL = absent). Have getSigmas
return NULL, as getDispersion does off nbinom. Guard plot's sigma panel with !fitIsHeteroscedastic.
fitFamily's legacy `is.null(sigma)` fallback is reached only when $family is missing, which a 1.0 fit never
is. Run the exact gates' fit-carried checks.
(b) Docs only: one clause in bartBT.Rd's Value sigma and in getSigmas, and skip the plot panel.
Recommend (a). The element is a constant that cannot be read as what its name says.

## multiforest-05 - CONFIRMED, MINOR (each bullet)

Own probe (r3-verify-mf-05.R):
```
fit-time basis factor(z, levels = 0:2): ACCEPTED; empty-level amplitude row of glue has sd 0.605 (prior only)
explicit cbind(z, 0): REFUSED "a 'basis' column of all zeros contributes nothing"
setForestBasis(rep(1e300)) and fit-time basis 1e300: ACCEPTED; next run train finite FALSE
variance = TRUE + linear("x1"): REFUSED "invalid sampler specification: either the leaf covariate designation ... or a variance forest ..."
NA y, basis pre-dropped to 49 rows, no subset: "matching 'subset' (49) but not the full data (50 rows)"
setState of a max.order-violating donor: "state is not consistent with this sampler"
installTrees of the same donor: names the interaction constraint
```
Cause of the last one: Sampler::setState separates only columnMaskStateFeasible and monotoneStateFeasible.
The interaction check sits inside Chain::stateIsValid. Related wording issue: setState's blocks refusal says
"warm-start donor" although no warm start is involved.
Fixes: refuse an empty factor level (the zero-column rule) and a non-finite basis row norm, at fit time and
in setForestBasis. Add the leaf model to R's variance-forest refusal. Name the NA-response row drop instead
of 'subset' when no subset was given. Separate an interactionRefused flag in Sampler::setState, as for the
column mask, and drop "warm-start donor" from the setState wording. None moves draws. Each takes one
expect_error test.
