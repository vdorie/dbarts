# Review 3, wave 2 - lens: multiforest

Pinned tree 01dee4b4, library r3-lib. Probes are under the scratchpad as r3-multiforest-*.R.

Covered:
- Successive-conditional (Geweke) checks I wrote for this review. Each draws theta from the sampler's
  own prior, simulates y, runs 10 sweeps, and z-tests paired differences (|z| bound about 3 over
  10-20 functionals). Arms:
  - BCF with glue fixed (update.amplitude = FALSE): multi-tree, tau blocks(), interactions()
    (max.order, forbid), vars, three forests (3-level factor plus numeric basis), numeric 2-column
    basis, weights, missing predictors, a categorical factor.
  - BCF with mutations between draws: setForestBasis (same width and numeric), setPredictor
    (whole matrix and one column), setOffset, setWeights, setActiveRows, setLeafPrior(forests =).
  - BCF probit and logistic.
  - BCF with prior-drawn amplitudes (half-Cauchy mixture variance included), installed through the
    state's glue block: K = 2, K = 3 with wide and numeric bases, all-basis, probit, K = 3 probit,
    offset swap, weight swap.
  - Heteroscedastic (sampleVarianceForestFromPrior): 1-3 variance trees, vars restriction,
    weights, missing predictors, active rows.
  - gp() and linear() leaves: gaussian, weighted, probit, logistic, student(df = 4), active rows,
    gp fallback regime.
  Each run used 20k-80k replications. Arms near |z| 2.7-2.9 were rerun on a new seed, and every
  one cleared.
- The exact gates' blind spots, listed so you know what my Geweke arms filled in: multi-tree
  forests, K > 2, non-indicator bases, blocks/interactions on a forest, mutation between draws,
  gp leaves (no gate at all), and leaf models under latent families.
- Per-forest bookkeeping: copy-equivalence after each mutation, i.e. the live sampler versus
  copy() of its stored state, run on bitwise draws. Covered BCF gaussian and probit (11
  mutators), hetero (9), multinomial (setCounts, setCategoryOffset, setPredictor,
  setActiveRows), gp and linear leaves (11 mutators). Also checked the documented
  fit-reconstruction identity and forest/contribution shapes and names, combined and
  uncombined. The per-forest varcount agrees with getTrees.
- Constraint enforcement: whether a blocks() or interactions() violation can be reached through
  growFromRoot, run, installTrees, setState (single forest and per BCF forest), setControl,
  vars plus blocks, trees.per.group against per-forest tree counts, the multinomial, ordinal,
  nbinom and hurdle families, and xbart (which refuses these arguments).
- Refusal matrix for BCF, gp/linear leaves and the variance forest. The forest.Rd sd claims
  (1.484, 0.699 -> 0.989) checked against getLeafPrior. setForestBasis edge inputs. Basis row
  alignment under subset and NA rows.

Not covered: drawn k under gp/linear (no prior setter), ordinal/nbinom/AFT with a variance
forest, multinomial exactness beyond its own gate, flat C API (bridge lens).

## multiforest-01 - BLOCKER - multinomial fits silently ignore interactions() and blocks()

Location: R/spec.R resolveSamplerSpec (the `unsupportedMultinomial` list);
src/R_interface_bartcore.cpp (the multinomial MultinomialSpec build);
src/bartcore/chain.hpp Chain::buildMultinomialForest ("no split restriction").

Claim: `family = "multinomial"` (dbarts() and bart()) accepts `interactions =` and `blocks =`
without complaint. The multinomial forest builder then installs neither, so every category forest
splits freely. interactions.Rd describes the constraint as "hard (an availability ban)" and
"applied per forest". The R refusal list next to it exists precisely so that a setting the
multinomial factory drops is "named instead of dropped in silence"; it names DART, split.probs,
monotone, linear/gp, a k hyperprior, a named sd and single storage, but not these two.

Probe (r3-multiforest-con5.R, con6.R): n = 300, 4 predictors, 3 categories, 6 trees, 100 sweeps.

```
dbarts(x, yc, family = "multinomial", interactions = interactions(max.order = 1L))
  -> multinomial max.order1 violations: 8 trees   (forests 1, 2, 3)
     forest 1 tree 3: root splits x3, its child splits x4  (two predictors on one path)
dbarts(..., blocks = blocks(list(c("x1","x3"), c("x2","x4"))))  -> 12 trees violate
bart(x, yc, family = "multinomial", interactions = interactions(forbid = c("x1","x2")))
  -> 21 trees put x1 and x2 on one path
same probes on ordinal, nbinom, hurdle.lognormal: 0 violations
```

Why gates missed: the tinytest constraint tests and the multinomial gate never combine the two,
and the R refusal list was built from what MultinomialForestSpec lacks, without checking it
against every constraint argument resolveSamplerSpec resolves.

Fix: add "interactions()" and "blocks()" to `unsupportedMultinomial`, with a bridge backstop.
Alternatively, carry the interaction and column-mask fields into MultinomialForestSpec and
install them per category forest, as buildSpecifiedForest does.

## multiforest-02 - BLOCKER - after setPredictor on a gp or linear leaf column, a copied or reloaded sampler predicts wrongly

Location: src/bartcore (the leaf covariate standardization constants, "sticky" under in-place
predictor changes per docs/design/linear-leaves.md "Mutation semantics" and gp-leaves.md),
ChainStateData (does not carry them); R/dbarts.R copy and getPointer re-creation (rebuild from
data@x, which setPredictor already updated).

Claim: setPredictor on a designated linear() or gp() column keeps the creation-time
standardization constants in the live engine, by design. Those constants are not in the stored
state. Any re-creation recomputes them from the mutated data@x and reinterprets the saved
coefficients and kernels under different constants. Re-creation covers $copy(), saveRDS/readRDS,
and a dead pointer. The reloaded sampler's current fit jumps, and predict() replays every saved
draw wrongly, silently.

Probe (r3-multiforest-lin2.R, gp2.R, lin.R): n = 80, linear(c("x1","x2")) or
gp(c("x1","x2")), keepTrees, run(30, 3), setPredictor(2 * runif(n), column = 2L), run(0, 3).
sd(y) = 0.95.

```
linear: live predict(new x) vs recorded train          2.2e-15
        copy()  predict(new x) vs live recorded train  1.18
        saveRDS/readRDS predict vs recorded train      1.18
        live vs copy current fits (no run)             0.047
gp:     live predict vs recorded train                 0.0013 (the documented nugget)
        copy() predict vs live recorded train          2.32
setData route (constants re-derived on both sides): live and copy both 8.9e-16
```

Why gates missed: test-mutate-then-serialize.R restores into a cold sampler built over the
ORIGINAL data and replays the mutation, so both sides hold the same sticky constants. It
compares states, not predictions. The copy() and readRDS route, where data@x already holds the
mutated column, is never exercised. The first wave's copy-continuation sweep did not mutate a
designated column.

Fix: serialize the per-column standardization constants in the chain state and restore them on
setState, so re-creation reproduces the live engine. Alternatively, make setPredictor on a
designated column re-derive them, as setData does.

## multiforest-03 - MAJOR - a variance-tree-count mismatch in installTrees is refused after the mean forests are replaced

Location: src/bartcore/sampler.hpp Sampler::installForests (the variance half,
`installVarianceForest` after the per-chain `installForest` loop); src/R_interface_bartcore.cpp
installForests (the varianceMismatch message).

Claim: installing a heteroscedastic donor whose variance forest has a different tree count is
not caught by the up-front shape gate. That gate pairs only presence with a non-empty block.
The mismatch fails in installVarianceForest after every chain's mean forest has been rebuilt from
the donor. The refusal then says the donor's variance trees "leave a leaf empty, a scale leaf is
not positive, or a flat tree failed to rebuild", when the cause is a tree count. This is distinct
from bridge-03: one chain suffices, no grid change is involved, and the trigger is a plain shape
mismatch that should be refused before any commit.

Probe (r3-multiforest-inst.R, inst2.R): donor varianceForest(n.trees = 4) with keepTrees, run(50, 3).
Target varianceForest(n.trees = 5) with n.chains = 1, run(30, 1).

```
Error in S2$installTrees(D): warm-start donor's variance trees cannot be installed on this
  sampler's data (a rebuilt variance tree leaves a leaf empty, ...)
mean fits changed by refused install: 0.18   (2 chains: 0.31)
variance changed by refused install:  0
```

Why gates missed: warm-start refusal tests use mean-tree-count, DART and grid mismatches, which
are checked up front. No test installs a donor that differs only in its variance tree count, and
none checks the sampler after a variance refusal.

Fix: compare the donor's and destination's variance tree counts in the shape gate (returning
shapeMismatch) before the commit loop. bridge-03's scratch-then-commit fix should include the
variance half.

## multiforest-04 - MINOR - a heteroscedastic fit's sigma channel carries range(y), contrary to its Value entry

Location: R/bart.R packaging of `sigma`/`first.sigma`; man/bartBT.Rd Value `sigma`;
man/dbartsSampler-class.Rd getSigmas.

Claim: under `variance =` the engine pins sigma at 1 on the internal scale. fit$sigma,
fit$first.sigma, run()$sigma and $getSigmas() therefore all report a constant equal to the
response range. bartBT.Rd's Value entry calls `sigma` "posterior samples of sigma, the
residual/error standard deviation". The getSigmas doc says "each chain's current residual
standard deviation". Only the `type` argument paragraph says it "is not its residual scale".
extract(fit, "sigma") correctly refuses.

Probe (r3-multiforest-hs2.R): y = 10 x1 + noise with sd 0.2 or 1.
`fit$sigma` = 11.31 11.31 ... = diff(range(y)); mean sqrt(s.train) = 0.72.

Why gates missed: tests check extract's refusal, not the raw element.

Fix: return NULL (or NA) for sigma/first.sigma on a heteroscedastic fit and from getSigmas, or
state the exception in both Value entries.

## multiforest-05 - MINOR - edge inputs and messages

- setForestBasis and fit-time bases accept a factor with an empty level. That leaves an all-zero
  column, which the "column of all zeros contributes nothing" refusal exists to stop, so an
  amplitude is left that only its prior moves. A basis of 1e300 is accepted and poisons the
  sampler: getLeafPrior shows anchor 0 and basis.row.norm Inf, and the next run returns NaN glue
  and non-finite train. Probe: r3-multiforest-basis.R. Fix: refuse empty factor levels, and
  refuse a non-finite row norm.
- variance = TRUE with linear() or gp() leaves is refused only by the bridge's disjunctive
  "invalid sampler specification: either the leaf covariate designation ... or a variance forest
  is combined with ...". The R variance check names Student-t and monotone but not leaf
  covariates, though varianceForest.Rd says "constant leaves only". Probe: r3-multiforest-comp3.R.
  Fix: add the leaf model to the R variance-forest refusal.
- A `forests =` basis already shortened by an NA response, with no `subset` given, is refused
  with "matching 'subset' (19)". Probe: r3-multiforest-sub.R. Fix: name the dropped rows rather
  than subset.
- setState of a donor that violates an interactions() constraint (single forest, or either BCF
  forest) is refused with the generic "state is not consistent with this sampler". installTrees
  and the blocks/vars refusals name the constraint. Probe: r3-multiforest-con1.R, ss.R.

## Checked and found correct

- Every Geweke arm listed above: BCF fixed and drawn glue, K = 2/3, factor/numeric/wide bases,
  blocks, interactions, vars, weights, offset, active rows, missing, categorical, probit,
  logistic. Hetero with weights, active rows, missing. gp and linear under gaussian, weighted,
  probit, logistic, student and active rows.
- Copy-equivalence after every BCF, hetero and multinomial mutator, and after gp and linear
  mutators other than the designated-column setPredictor. The gp continuation divergence with no
  mutation is documented (gp-leaves.md stage 4 addendum). The gp predict-at-training-rows nugget
  gap is documented.
- Constraints hold through run, growFromRoot (BCF per forest included), installTrees refusals,
  setState refusals, vars plus blocks by name and index, trees.per.group per forest; ordinal,
  nbinom and hurdle honour interactions/blocks.
- bart() forest() term fits: extract/predict type = "forest", contribution, combineChains both
  ways, names forest1/forest2, glue "forest" attribute, the reconstruction identity
  yhat = shift + sum (B glue) forestFits, predict(newdata) = yhat.train.
- Basis rows align with subset, NA-response and na.omit row drops. forest.Rd's induced-prior
  claims hold to 4 digits. setForestBasis width remap matches its doc.
- BCF refusals: variance, monotone, Student-t, non-default k, named sd, split.probs, single
  storage, test predictors, setData, setResponse(updateScale = TRUE).
- Variance-forest installTrees from a matching donor installs the saved sample's surface exactly.
