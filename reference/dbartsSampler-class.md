# Class "dbartsSampler" of Discrete Bayesian Additive Regression Trees Sampler

A reference class object that contains a Bayesian Additive Regression
Trees sampler in such a way that it can be modified, stopped, and
started all while maintaining its own state.

## Usage

``` r
# S4 method for class 'dbartsSampler'
run(
  numBurnIn, numSamples, updateState = NA, n.threads = control@n.threads
)
# S4 method for class 'dbartsSampler'
sampleTreesFromPrior(updateState = NA)
# S4 method for class 'dbartsSampler'
sampleNodeParametersFromPrior(updateState = NA)
# S4 method for class 'dbartsSampler'
growFromRoot(n.sweeps = 2L, updateState = NA)
# S4 method for class 'dbartsSampler'
copy(shallow = FALSE)
# S4 method for class 'dbartsSampler'
show()
# S4 method for class 'dbartsSampler'
predict(x.test, offset.test, n.threads = control@n.threads)
# S4 method for class 'dbartsSampler'
setControl(control)
# S4 method for class 'dbartsSampler'
setModel(model)
# S4 method for class 'dbartsSampler'
setData(data)
# S4 method for class 'dbartsSampler'
setResponse(y, updateScale = FALSE, updateState = NA)
# S4 method for class 'dbartsSampler'
setOffset(offset, updateScale = FALSE, updateState = NA)
# S4 method for class 'dbartsSampler'
setWeights(weights, updateState = NA)
# S4 method for class 'dbartsSampler'
setActiveRows(active, updateState = NA)
# S4 method for class 'dbartsSampler'
setForestWeights(forest, weights, updateState = NA)
# S4 method for class 'dbartsSampler'
setForestBasis(forest, basis, updateState = NA)
# S4 method for class 'dbartsSampler'
setSigma(sigma, updateState = NA)
# S4 method for class 'dbartsSampler'
setPredictor(x, column, forceUpdate, updateCutPoints = FALSE, updateState = NA)
# S4 method for class 'dbartsSampler'
setCutPoints(cuts, column, updateState = NA)
# S4 method for class 'dbartsSampler'
setTestPredictor(x.test, column)
# S4 method for class 'dbartsSampler'
setTestPredictorAndOffset(x.test, offset.test)
# S4 method for class 'dbartsSampler'
setTestOffset(offset.test)
# S4 method for class 'dbartsSampler'
printTrees(treeNums, chainNums, sampleNums)
# S4 method for class 'dbartsSampler'
getTrees(
  treeNums, chainNums, sampleNums, current = FALSE, newdata = NULL
)
# S4 method for class 'dbartsSampler'
getSigmas(result)
# S4 method for class 'dbartsSampler'
getLatents(result)
# S4 method for class 'dbartsSampler'
getSumsOfSquaredResiduals(result)
# S4 method for class 'dbartsSampler'
getForestFits(forest)
# S4 method for class 'dbartsSampler'
getForestAmplitudes(forest = NULL)
# S4 method for class 'dbartsSampler'
getForestVariableCounts(forest)
# S4 method for class 'dbartsSampler'
getCalibration(forest = 1L)
# S4 method for class 'dbartsSampler'
setCalibration(
  prior.scale, prior.sd, prior.mean, forest = 1L, updateState = NA
)
# S4 method for class 'dbartsSampler'
installTrees(donor, samples = NULL)
# S4 method for class 'dbartsSampler'
storeState()
# S4 method for class 'dbartsSampler'
setState(newState)
# S4 method for class 'dbartsSampler'
plotTree(
  treeNum, chainNum, sampleNum, treePlotPars = c(
    nodeHeight = 12, nodeWidth = 40, nodeGap = 8),
  ...
)
```

## Note

`dbartsSampler` is a reference class: its methods are not called as free
functions but as `$`-dispatched calls on a sampler instance, e.g.
`sampler$run(100, 100)` or `sampler$setResponse(newY)`. The S4-method
notation in ‘Usage’ below is an artifact of how reference-class methods
are documented and does not reflect the calling syntax; see ‘Examples’.

## Arguments

- numBurnIn:

  A non-negative integer determining how many iterations the sampler
  should skip before storing results. If missing or `NA`, the default is
  filled in from the sampler's
  [`control`](https://vdorie.github.io/dbarts/reference/dbartsControl.md)
  object.

- numSamples:

  A positive integer determining how many posterior samples should be
  returned. If missing or `NA`, the default is also filled in from the
  control object.

- updateState:

  A logical determining if the local cache of the sampler's state should
  be updated after the call completes. Two conventions apply, by method:
  for `run`, `sampleTreesFromPrior`, and
  `sampleNodeParametersFromPrior`, `NA` (the default) fills in the
  sampler's
  [`control`](https://vdorie.github.io/dbarts/reference/dbartsControl.md)
  object's `updateState`, and explicit `TRUE`/`FALSE` override it. For
  the mutators - `setData`, `setResponse`, `setOffset`, `setWeights`,
  `setForestBasis`, `setSigma`, `setCalibration`, `setPredictor`, and
  `setCutPoints` - the state is stored only on explicit `TRUE`; `NA`
  (the default) and `FALSE` both store nothing, regardless of
  `control@updateState`. These are typically called once per sweep
  inside a larger Gibbs/MH loop (as `dbartsSampler` is designed for),
  where storing state on every mutation would be wasted work whenever
  the loop only reads `state` occasionally (or never); an unforced
  `state` promise materializes the sampler's *current* state on first
  access regardless, so a mutate-then-first-read sequence needs no
  explicit store. Pass `TRUE` explicitly when `state` was already forced
  (read or saved) earlier and a later mutation must be reflected in the
  next save.

- shallow:

  A logical determining if the copy should retain the underlying data of
  the sampler (`TRUE`) or have its own copies (`FALSE`).

- control:

  An object inheriting from
  [`dbartsControl`](https://vdorie.github.io/dbarts/reference/dbartsControl.md).
  When passed to `setControl`, it cannot change `n.trees`, `n.chains`,
  `useQuantiles`, or `rngSeed` from the values the sampler was created
  with, and cannot set `keepTrees = TRUE` without also giving
  `n.samples`; either is an error.

- model:

  An object inheriting from `dbartsModel`. When passed to `setModel`, it
  cannot switch a DART tree prior on or off relative to the sampler's
  creation-time model; recreate the sampler instead.

- data:

  An object inheriting from
  [`dbartsData`](https://vdorie.github.io/dbarts/reference/dbartsData.md).

- y:

  A numeric response vector of length equal to that with which the
  sampler was created, and lying in the response family's support: 0/1
  for `probit` and `logistic`, an integer category index in \\\[1, K\]\\
  for `ordinal`, and a finite non-negative integer count for `nbinom`.
  Values off the support are refused, as they are at creation;
  `gaussian` and `aft` (log survival times) constrain nothing.

- x:

  A numeric predictor vector of length equal to that with which the
  sampler was created. Can be of a distinct number of rows for
  `setTestPredictor`.

- x.test:

  A new matrix, data frame, or sparse-bearing test set (the column types
  accepted by
  [`dbartsData`](https://vdorie.github.io/dbarts/reference/dbartsData.md)),
  of the number of columns equal to that in the current model. A
  data-frame/sparse `x.test` is coded against the training levels;
  `setTestPredictor` and `setTestPredictorAndOffset` then install it as
  a resident, whole-object replacement of the test set (see `column`),
  while `predict` materializes it to a numeric matrix before evaluating.

- offset:

  A numeric vector of length equal to that with which the sampler was
  created, or `NULL`. If `offset.test` was set from `offset`, will
  attempt to update that as well.

- updateScale:

  Logical indicating whether BART's internal scale should re-anchor to
  the new offset (`setOffset`) or response (`setResponse`). Defaults to
  `FALSE`, locking the scale set at creation; should only be `TRUE`
  during burn-in, as re-anchoring mid-run makes the fits across
  iterations no longer comparable. `TRUE` is refused by samplers whose
  leaf calibrations are stated against the transform fixed at creation
  and are never restated: a multi-forest (`forests`) sampler - refused
  whatever its response family, the calibration map being pinned at
  creation whether or not there is a transform to re-anchor - and a
  heteroscedastic (`variance`) one, whose variance forest would
  otherwise keep reporting \\s^2(x)\\ on the abandoned scale.

- offset.test:

  A numeric vector of length equal to that of the test matrix, or
  `NULL`. Can be missing for `setTestPredictorAndOffset`.

- n.threads:

  Currently has no effect: `run` and `predict` both execute serially
  regardless of the value passed here. The sampler's own thread count is
  fixed when it is created, from the `n.threads` of its
  [`control`](https://vdorie.github.io/dbarts/reference/dbartsControl.md)
  object; this argument is reserved for a future per-call override.

- sigma:

  A single positive residual standard deviation, on the original
  response scale; it is assigned to every chain. `setSigma` only applies
  to a sampler whose residual standard deviation is a free parameter: it
  is refused on a probit-, logistic-, multinomial-, ordinal-, or
  count-family sampler, which fix it at 1 by the model definition, and
  on a heteroscedastic (`variance`) sampler, whose variance forest owns
  the residual scale. A gaussian sampler with `resid.prior = fixed()` is
  not fixed in this sense and is still accepted: suppressing the
  internal draw so an outer sampler owns `sigma` is what that prior is
  for.

- weights:

  A numeric vector of non-negative case weights, of length equal to that
  with which the sampler was created. Weights of zero exclude an
  observation from the likelihood while keeping its fitted values.
  `setWeights` only applies to a gaussian-family sampler; calling it on
  a probit- or logistic-family sampler is refused, since a weighted
  probit has no tractable latent-variable form and logistic weights
  (observation counts) are fixed when the sampler is created. See
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md) for
  the family-specific weight rules that apply at creation time.

  For `setForestWeights`, a distinct per-FOREST case weight on a
  Bayesian causal forest (built with `forests = `, see
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md)),
  refused off one: it composes with `weights` and `active` as \\(w_i
  a_i) m_f^2 s_i\\, where \\m_f\\ is the named forest's own multiplier,
  rather than widening either channel, so `s_i = 0` excludes row \\i\\
  from *that forest's own leaf conditionals only* - its occupancy, its
  place in the combination, and the residual degrees of freedom are
  unaffected, and the reported location still carries the forest's full
  contribution even when every row is excluded. The third leg of a
  three-way degenerate-value contrast: `weights` installs and is
  measurably distinct from carrying none, `active` at all-ones installs
  nothing and clears any mask in force, and `setForestWeights` at
  all-ones INSTALLS - a round trip reports it - but is bitwise IDENTICAL
  to carrying no per-forest weight at all, the null gate its
  multiplicative composition (\\m_f \times 1 = m_f\\) is built to
  guarantee. The weight does not ride the sampler's saved `state` - a
  sampler rebuilt with `setState` from a stored state silently drops it
  while `statesAgree` still reports agreement - so it is additionally
  mirrored on an R5 field that `getPointer`, `setState`, and `copy` all
  reinstall on every re-creation; a caller never reinstalls it by hand.
  `setData` needs no clearing rule here the way it does for `active`: it
  is refused outright on any multi-forest sampler (“a multi-forest
  sampler fixes its data at creation”), so the two channels never
  interact.

- active:

  A numeric vector of per-observation 0/1 membership indicators - "row
  \\i\\ is not in the data set for this sampler this sweep" - of length
  equal to that with which the sampler was created, or `NULL`. Unlike
  `weights`, `setActiveRows` reaches families that refuse case weights
  entirely: gaussian, Student-t, `probit`, `ordinal`, `logistic`,
  `nbinom`, and `aft` samplers all accept it, including a multi-forest
  (`forests`) sampler, whose forests share one response of whichever
  family it was built under. It is absolute and independent of
  `weights`, composing with it in either call order (\\w_i \times a_i\\
  in effect). `NULL` clears the mask (every row active). Length and `NA`
  are validated as for `weights`; the only legal values are exactly `0`
  and `1`, and a fractional element refuses the WHOLE call and installs
  nothing. An all-ones vector reports success but installs NOTHING,
  clearing any previously installed mask - the OPPOSITE of `weights`,
  where an all-ones vector installs and is measurably distinct from
  carrying no weights at all. One channel states membership, the other
  precision, and the two deliberately carry opposite degenerate-value
  policies.

  An inactive row (`active[i] == 0`) contributes nothing to any leaf
  sufficient statistic, branch log-likelihood, birth-scan weight total,
  leaf parameter draw, or family-level parameter update that sums over
  rows (residual degrees of freedom, a dispersion statistic, a group's
  per-group sums, and so on). It still occupies its leaf for every
  count-based accounting - `numObservations`, the empty-leaf veto, leaf
  collapsing - and it still receives a fitted value: `run()$train`,
  `getForestFits`, and `predict` report \\f(x_i)\\ at an inactive row
  exactly as at an active one, which is what makes this channel worth
  more than physically dropping the row.

  An inactive row's own latent draw is SKIPPED - no random numbers are
  consumed for it - for every family except Student-t, whose per-row
  \\\lambda_i\\ is drawn unconditionally at every row regardless of the
  mask; there the mask instead ANNIHILATES its contribution through the
  composed weight (\\w_i \lambda_i a_i\\). For a skipping family,
  `getLatents` at an inactive row returns its LAST DRAWN value, which is
  stale; the correct read at an inactive row is its fit, not its latent.
  Reactivating a row (a mask change from inactive to active) is a
  one-sweep MODEL hazard, not only a read hazard: that sweep's tree
  moves and leaf draws run against the row's stale latent (or, for
  Student-t, a \\\lambda_i\\ drawn while the row was out) before that
  sweep's own latent refresh updates it; this does not disturb the
  posterior while a mask is held fixed, but a caller that moves the mask
  every sweep should expect it.

  An inactive row's pointwise log-likelihood - the value the flat C
  API's `logLikelihood` results field reports (see ‘Mutation cost’ above
  for the C interface) - is `NaN` rather than a finite value or `-Inf`:
  a row that is not in the model has no likelihood to report.

  Three quantities do NOT follow the mask, by design, so a caller
  wanting them recomputed on the active subset must re-create the
  sampler instead: the split-candidate cut grid stays the FULL-data
  grid, holding the tree prior fixed while the row set moves; the
  response transform and residual prior scale stay the FULL-data
  calibration (a masked gaussian - and so `aft`, and any `variance` or
  grouped decorator - keeps the range and residual-prior scale of every
  row); and an ordinal sampler keeps its FULL-data number of categories,
  free cutpoints, and log-gap prior even when a mask empties a boundary
  category - only `setData` changes them. `setData` also CLEARS an
  installed mask - the mask is a per-row vector and `setData` may change
  the number of rows - so a caller replacing the data must reinstall any
  mask it wants kept.

  The mask does not ride a sampler's saved `state`: a sampler rebuilt
  with `setState` from a stored state silently drops the mask and
  computes different draws while `statesAgree` still reports agreement.
  Conversely, the latents DO ride the state, and under an installed mask
  they are stale at every inactive row for a skipping family, so two
  samplers at the same posterior state but different mask histories
  carry different stored latents and `statesAgree` correctly reports
  DISAGREEMENT; a Student-t sampler is affected the same way, since an
  inactive row's stored \\\lambda_i\\ is drawn from a conditional its
  own data no longer informs.

  On a sampler with grouped random effects, a group all of whose rows
  are inactive draws its effect from the PRIOR through the group-effect
  update's own formula - the coherent answer, and NOT what a compacted
  sampler with that group physically deleted would do.

- basis:

  The data the named forest's amplitudes multiply, with as many rows as
  the sampler was created with and one amplitude per column: a factor
  (or a one-sided formula naming one) expands to its level indicators
  with no reference level dropped - a two-level factor giving the pair
  whose amplitudes are \\(b_0, b_1)\\ - and a numeric vector or matrix
  is already those columns, which is the same expansion
  [`forest`](https://vdorie.github.io/dbarts/reference/forest.md)
  applies at creation. Any forest takes a basis of any width.
  `setForestBasis` is the *sole* route by which a basis changes after
  creation, and it applies only to a sampler whose forests carry
  amplitudes, built with `forests = ` (see
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md)). The
  amplitudes are preserved and remapped: a width-preserving install
  leaves every one of them bitwise, and a width change carries each
  forest's block to its new offset and enters the added coordinates at
  `1`. A basis that is not one of the canonical shapes (a dense all-ones
  column, or a complementary 0/1 pair) moves that forest onto the
  general amplitude conditional, which is a model fact and not
  revertible by a later data swap. The matrix is written to both the
  engine and `data@bases`, so `getPointer`'s transparent re-creation
  after save/load carries the current assignment rather than the one the
  sampler was created with; because the bases ride creation, a widening
  applied after a state restore preserves the *restored* amplitudes.

- forest:

  A single positive integer indexing the forests from `1`.
  `getCalibration` and `setCalibration` default to `1`, the only forest
  of an ordinary sampler, and `getForestAmplitudes` defaults to `NULL`,
  every forest's block stacked; `getForestFits`,
  `getForestVariableCounts`, `setForestWeights`, and `setForestBasis`
  have no default. A Bayesian causal forest's prognostic forest is `1`
  and its basis forest `2`; `setForestWeights` and `setForestBasis` are
  refused on a sampler whose forests carry no amplitudes, and
  `setForestBasis` accepts any forest of one that does. `getForestFits`
  and `getForestVariableCounts` accept `forest = 1` on any sampler - it
  selects the only forest - and refuse only an out-of-range index.
  `getCalibration` is likewise served on every forest of a multi-forest
  sampler, and its calibration-map columns are that forest's own.

- prior.scale:

  For `setCalibration`, the prior standard deviation of the forest total
  at `k = 1`, in response units (the family's latent units where the
  response is not rescaled). This is the identified quantity: only the
  ratio of the leaf scale to `k` enters a draw law, so under a fixed `k`
  the pair has one degree of freedom, and `prior.scale` is the half of
  it a hyperprior on `k` leaves alone. Must be a single positive finite
  number; exactly one of `prior.scale` and `prior.sd` is given.

- prior.sd:

  For `setCalibration`, the same statement at the `k` currently in
  force, so `prior.scale = prior.sd * k`. Refused when `k` is drawn from
  a hyperprior - it would name a different prior every sweep - and
  refused when the chains' `k` have diverged, since one number would
  then mean a different scale on each; `prior.scale` serves in both
  cases. Note the binary families default to a sampled `k`.

- prior.mean:

  Not writable, and present only to say so with its remedy: the leaf
  values it would shift are already drawn and stored. The prior mean of
  the forest total is the response transform's shift, and the lever that
  moves the modelled quantity is the offset channel,
  `setOffset(rep_len(-getCalibration()[1, "prior.mean"], n))`.

- cuts:

  Vector of cut points for use with `setCutPoints`.

- column:

  An integer or character string vector specifying which column/columns
  of the predictor matrix is to be replaced or for which cuts should be
  used. When updating predictors, it can be missing and the entire
  matrix will be substituted; a sparse or mixed dense/sparse design
  accepts that too, and keeps its storage (the replacement is spliced
  into the container rather than replacing it with a plain matrix). Such
  a design may equally be mutated one column at a time by naming
  `column` (dense-backed and sparse-backed columns alike, and a single
  call may name both kinds), which is the cheaper spelling when only a
  few columns move. Note that replacing a sparse-backed column densifies
  its storage permanently: every row then differs from the column's
  implicit value, so it stores `n` entries from that point on. For
  `setTestPredictor`, a per-column update is refused when the current
  test set is a data frame or sparse container rather than a plain
  matrix; replace the whole test set instead.

- forceUpdate:

  For `setPredictor`, controls how the trees respond to a new predictor
  that would leave a leaf empty in any tree of any forest of any chain -
  the ordinary mean forest, and, on a sampler with more than one forest,
  all of them together: a Bayesian causal forest's second (treatment)
  forest, a multinomial sampler's per-category forests, and a
  heteroscedastic (`variance`) sampler's variance forest are checked on
  exactly the same terms as the mean forest, not independently (see
  ‘Multi-forest and heteroscedastic predictor mutation’ below). `FALSE`
  rolls the whole update back and leaves the sampler unchanged
  (rejection sampling); `TRUE` forces the change, collapsing any empty
  nodes; the string `"partial"` installs each observation's new value
  individually and rolls back only those observations whose value would
  empty a leaf in any of those trees, returning a per-observation
  logical of what was installed. `"partial"` requires a single `column`
  and cannot be combined with `updateCutPoints`; it also requires a
  dense-backed column, since a sparse column's storage fixes its nonzero
  pattern per cell - replace such a column whole instead. When missing,
  defaults to `TRUE` when the whole predictor matrix is being replaced
  (`column` is missing) and `FALSE` when a single column is being
  replaced.

- updateCutPoints:

  For `setPredictor`, a logical; when `TRUE` the cut points (split
  candidate locations) for the replaced column(s) are recomputed from
  the new values, otherwise the existing cut points are kept. Defaults
  to `FALSE` and cannot be combined with `forceUpdate = "partial"`.

- treeNums:

  An integer vector listing the indices of the trees to print or return.

- chainNums:

  An integer vector listing the chains to return from `getTrees`. When
  missing, all chains are returned.

- sampleNums:

  An integer vector listing the saved samples to return from `getTrees`.
  Applies only when `keepTrees` is `TRUE` and `current` is `FALSE`;
  otherwise the live working trees have no sample dimension and it is
  ignored.

- current:

  For `getTrees`, a logical; when `TRUE` the live working trees (the
  current sampler position) are returned even if `keepTrees` is `TRUE`,
  rather than the saved samples.

- newdata:

  For `getTrees`, an optional matrix of predictors (in the form accepted
  by `predict`). When supplied, its observations are routed through each
  tree so the `n` column counts them instead of the training data; the
  tree structure, splits, and leaf values are unchanged.

- result:

  For `getLatents`, an optional pre-allocated numeric vector (or, with
  multiple chains, matrix) that it fills in place and returns instead of
  allocating a new one; its length must equal the number of observations
  times the number of chains. Omit it to let `getLatents` allocate.
  Accepted but not used by `getSigmas` and `getSumsOfSquaredResiduals`.

- treeNum:

  An integer listing the indices of the tree to plot.

- chainNum:

  For `plotTree`, an integer giving the chain to plot from. Required
  when the sampler has more than one chain; defaults to `1` when there
  is only one.

- sampleNum:

  For `plotTree`, an integer giving the saved sample to plot. Applies
  only when `keepTrees` is `TRUE`, in which case it defaults to the most
  recently drawn sample; ignored (with the live working trees plotted
  instead) when `keepTrees` is `FALSE`.

- treePlotPars:

  A named numeric vector containing the quantities `nodeHeight`,
  `nodeWidth`, and `nodeGap`, all of which control aspects of the
  resulting plot.

- donor:

  For `installTrees`, the fit whose forests seed this sampler: another
  `dbartsSampler`, or a `bart` object fit with `keepSampler = TRUE`. The
  donor must share this sampler's predictors, tree count, and DART
  setting; a different cut grid is allowed, its splits remapped onto
  this sampler's grid.

- samples:

  For `installTrees`, an optional integer vector with one entry per
  chain, each a 1-based index into the donor's pool of samples (its
  saved trees when it kept them, else its final trees per chain). `NULL`
  spreads the chains evenly across the pool.

- n.sweeps:

  For `growFromRoot`, a single positive integer giving the number of
  grow-from-root sweeps run in place (default `2L`). Each sweep rebuilds
  every tree in the forest against the current residual, so a small
  handful reaches a good fit.

- newState:

  For `setState`, a state object previously produced by this sampler
  (its `state` field, or the return of `storeState`) over the same
  model. Must inherit from `bartcoreState`.

- ...:

  Extra arguments to
  [`plot`](https://rdrr.io/r/graphics/plot.default.html).

## Fields

Read with `$`, as the methods are called. Treat them as read-only:
assigning to a field changes only the R-side copy and leaves the
underlying engine untouched, so the sampler's behavior does not follow.
Route changes through the `set*` methods instead.

- `control`:

  The sampler's
  [`dbartsControl`](https://vdorie.github.io/dbarts/reference/dbartsControl.md),
  carrying `n.chains`, `n.threads`, `n.trees`, `keepTrees`,
  `updateState`, and the rest. Replace with `setControl`.

- `model`:

  A `dbartsModel` holding the parsed tree, node, and residual priors,
  together with any monotonicity, interaction, and block constraints.
  Replace with `setModel`.

- `data`:

  A
  [`dbartsData`](https://vdorie.github.io/dbarts/reference/dbartsData.md)
  holding the current predictors, response, test data, weights, and
  offset. Its `x` slot is stored columnar for data-frame input and so
  need not be a plain matrix; use
  [`extract`](https://vdorie.github.io/dbarts/reference/extract.dbartsSampler.md)`(sampler, "predictors")`
  for the numeric predictor matrix, which is what tree replay and
  `getTrees` count against. Replace with `setData`, or a piece at a time
  with the other `set*` methods.

- `state`:

  The cached, serializable engine state, or `NULL` until one is
  materialized. Reading it forces the sampler's *current* state;
  `storeState` refreshes it on demand and `updateState` governs when the
  methods do so themselves. It is the only field
  [`save`](https://rdrr.io/r/base/save.html) needs, and restoring one
  requires `setState` - see ‘Saving’.

- `pointer`:

  The external pointer to the C++ engine. Internal; never manipulate it
  directly.

## Details

A `dbartsSampler` is a mutable object which contains information
pertaining to fitting a Bayesian additive regression tree model. The
sampler is first created and then, in a separate instruction, run or
modified. In this way, MCMC samplers can be constructed with BART
components filling arbitrary roles.

### Saving

[`save`](https://rdrr.io/r/base/save.html)-ing and
[`load`](https://rdrr.io/r/base/load.html)ing a `dbarts` sampler for
future use requires that R's serialization mechanism be able to access
the state of the sampler which, for memory purposes, is only made
available to R on request. To make it available, call `storeState()` on
the sampler before saving, e.g. for the object `sampler`, execute
`sampler$storeState()`; this captures the sampler's current state into
the serializable `state` field regardless of the `updateState` setting.

The state object is opaque and engine-specific, and carries a format
version and the package version that wrote it. A state loads only within
the format version of the `dbarts` that wrote it; loading it with an
incompatible version refuses cleanly, naming both versions, rather than
risk a silent misread. There is no cross-version migration: re-fit the
model, or restore the state with the `dbarts` release that wrote it.

To restore a saved state into a sampler, call `setState(newState)`: it
validates that `newState` inherits from `bartcoreState`, re-creates the
underlying engine if needed, pushes the state into it, and caches it on
the `state` field. Assigning the field directly
(`sampler$state <- newState`) does *not* restore the sampler - it only
overwrites the R-side cache, leaving the engine untouched, so the next
run continues from the engine's own state rather than the assigned one.
Always route a restore through `setState`.

### Mutation cost

Each accepted predictor mutation does two things: the engine updates its
internal representation, and the R layer collects the accepted values
into the sampler's data object so that replay, saving, and
re-quantization see current data. The collection is by reference for
data-frame-built samplers (only the affected column changes hands) but
requires copying the full predictor matrix when the sampler was built
from a matrix, so tight loops that mutate predictors every iteration are
better served by data-frame input. The remaining per-call overhead is a
few tens of microseconds of R method dispatch and bookkeeping - the
price of R-level consistency. Clients for which that matters can drive
the sampler through the C interface (see the `dbarts.h` header installed
with the package), which invokes the engine directly and performs no
R-side collection; such clients supply their own current predictor
matrix when replaying saved trees.

### Multi-forest and heteroscedastic predictor mutation

`setPredictor` - whole-matrix, single- or multi-column (via `column`),
and per-observation (`forceUpdate = "partial"`) - and
[`updatePredictorPerObservationJointly`](https://vdorie.github.io/dbarts/reference/updatePredictorPerObservationJointly.md)
accept a Bayesian causal forest, a multinomial sampler, and a
heteroscedastic (`variance`) sampler on the same terms as an ordinary
single-forest sampler. The acceptance criterion is a single conjunction
over every tree of every forest of every chain of the sampler: a
treatment forest, a multinomial sampler's per-category forests, and a
variance forest are all checked together with the mean forest, never
independently, so a change that a mean-forest tree would tolerate can
still be rejected because it would empty a leaf in the treatment or
variance forest, and vice versa. A rejected transactional whole-matrix
or column update rolls *every* forest back to its pre-call state; under
`forceUpdate = "partial"` the per-observation install mask reflects that
same joint criterion, so an observation installs only where it empties
no leaf anywhere in the sampler, and a rejected observation is rolled
back in every forest, leaving the rest of the run untouched. The one
exemption is a forest whose trees never split on the column being
changed - a treatment or moderator forest restricted to a column subset,
or a `variance` forest declared over its own predictor subset - which
structurally cannot veto a change to a column outside its own reach;
nothing needs to be done to obtain that exemption; it falls out of which
columns the forest's trees actually use.
`updatePredictorPerObservationJointly` applies the identical per-sampler
conjunction across every sampler passed to it as well: an observation
installs only where it would remain valid in the union of every
sampler's forests, and it is declined in all of them otherwise.

Column names after a whole-matrix replacement: a transactional or forced
whole-matrix `setPredictor(x, ...)` (`column` missing) installs `x` into
`sampler$data@x` as a plain matrix, which does not carry `x`'s column
names even when the sampler was originally built with named columns.
Once that has happened, every later *by-name* column reference on that
sampler fails: `setPredictor`'s own `column` argument, given as a
character name, resolves it through `colnames(sampler$data@x)` (now
`NULL`) and errors with “column names not specified at initialization”,
and so does
[`updatePredictorPerObservationJointly`](https://vdorie.github.io/dbarts/reference/updatePredictorPerObservationJointly.md)'s
shared-column match, whether its own `column` is given by name or by
integer (it still resolves a shared name across every referenced
sampler). An *integer* `column` to `setPredictor` itself is unaffected,
since it skips the name lookup entirely. Two reliable workarounds: pass
integer column indices to any later `setPredictor` call on the affected
sampler; and, for `updatePredictorPerObservationJointly` specifically
(integer `column` does not help there), either replace predictors a
column at a time (`setPredictor(x, column = ...)`), which merges into
the existing matrix in place and so keeps its dimnames rather than
replacing the whole matrix, or call
`updatePredictorPerObservationJointly` before the first whole-matrix
replacement in a script that needs both.

### Warm starts

`installTrees` seeds the sampler's forests from a `donor` instead of
drawing trees from the prior, for scaling to more chains or embedding a
fit in a larger sampler. Only the donor's trees, `sigma`, and `k`
transfer; each chain keeps its own random-number stream and redraws
everything else, so several chains seeded from one donor stay
overdispersed. A donor with a different tree count or DART setting is
refused rather than silently reshaped; a donor fit on a different cut
grid is instead remapped onto this sampler's grid, collapsing any splits
the grid starves, the same way a data replacement remaps existing
splits. A warm start biases the early draws toward the donor, so it
shortens burn-in rather than removing it; keep drawing a non-zero number
of burn-in samples before treating the chain as converged.

`growFromRoot` instead builds the sampler's initial forest by
XBART-style root-down stochastic tree construction (He, Yalov and Hahn
2019): each of `n.sweeps` sweeps rebuilds every tree from the root,
sampling each split from the integrated-likelihood weight of every
candidate cut under the tree prior. This reaches a good fit in far fewer
sweeps than the exact sampler, so it is a fast starting point rather
than a posterior sampler - the exact MCMC sweeps own stationarity once
`run` begins, and the posterior is unchanged. It is available for the
constant-leaf model only; calling it on a `linear` or `gp` node prior is
an error, not a silent fall-back, so initialize those forests with
`sampleTreesFromPrior` instead. As with `installTrees`, the grown forest
biases the early draws toward its fit, so shorten burn-in rather than
skipping it. Each chain grows on its own random-number stream, so the
result is independent of `n.threads`. A composable cross-sampler
workflow follows for free: `donor$growFromRoot(k)` then
`target$installTrees(donor)`.

## Value

For `run`, a named-list with contents `sigma`, `train`, `test`, and
`varcount` (plus `k`, `varprobs`, `tau`, and `ranef` when applicable).
`train` is an array of dimension n.obs x n.samples x n.chains, and
likewise `test` (or `NULL` if the sampler has no test data) and
`varcount` are n.predictors x n.samples x n.chains; `sigma` is n.samples
x n.chains. When `n.chains` is `1` the trailing chain dimension is
dropped, so `train` is a plain n.obs x n.samples matrix and `sigma` a
plain vector of length n.samples. A Bayesian causal forest adds two
more: `forestFits`, an n.obs x n.forests x n.samples x n.chains array of
each forest's fitted values on the internal scale (the prognostic
\\\mu\\ first, the treatment \\\tau\\ second), and `glue`, a sum(q) x
n.samples x n.chains array of the amplitudes each draw combines them
through, stacked forest-major and as wide as each forest's own basis (a
Bayesian causal forest's \\(a, b_0, b_1)\\, three rows), so both
surfaces and their recombination come from a single call; `train`
carries the combination \\a \mu(x_i) + b\_{z_i} \tau(x_i)\\ on the
response scale, or on the latent scale when the sampler was built under
`"probit"` or `"logistic"`, and `test` is filled with `NaN` (there is no
test treatment vector to combine off-sample). No other model reports
either element. A run can be interrupted with `Ctrl-C`: it stops between
iterations - joining any worker threads first - and signals an error,
returning no samples from the interrupted run. The sampler's chains are
left at the iteration they reached, which is a valid state to run again
from.

For `setPredictor`, `TRUE`/`FALSE` depending on whether or not the
operation was successful. The operation can fail if the new predictor
results in an empty leaf-node in any tree of any forest of the sampler
(see ‘Multi-forest and heteroscedastic predictor mutation’). If only
single columns were replaced, the update is rolled back so that the
sampler remains in a valid state. When `forceUpdate` is `"partial"`,
instead returns a logical vector of length equal to the number of
observations, `TRUE` where that observation's new value was installed
and `FALSE` where it was rolled back to its previous value to keep every
tree of every forest valid.

`predict` keeps the current test matrix in place and uses the current
set of tree splits. This function has two use cases. The first is when
`keepTrees` of
[`dbartsControl`](https://vdorie.github.io/dbarts/reference/dbartsControl.md)
is `TRUE`, in which case the sampler should be run to completion and the
function can be used to interrogate the existing fit. When `keepTrees`
is `FALSE`, the function can be used to obtain the likelihood as part of
a proposed new set of covariates in a Metropolis-Hastings step in a
full-Bayes sampler. This would typically be followed by a call to
`setPredictor` if the step is accepted.

For `getTrees`, a `data.frame` with one row per tree node in
depth-first, left-hand-side pre-order, with columns `chain`, `sample`
(present only for saved samples), `tree`, `n` (the number of
observations in the node), `var` (the splitting variable, or -1 at a
leaf), and `value` (the split value, or the leaf prediction). An ordinal
rule's value is its cut point and observations with values less than or
equal to it go left; a categorical rule carries no data value (its
`value` is `NA`) and its split is reported in the `directions` column
instead. When the sampler has any categorical predictors the result
gains a `directions` column decoding each categorical rule into one
`"L"`/`"R"` character per level, in level order (level `k` goes right
when its character is `"R"`); ordinal rules and leaves are `NA`. When
any predictor contains missing values the result gains a `missing`
column giving the branch (`"L"`/`"R"`) each rule sends missing values
down; rules on complete columns and leaves are `NA`. Under a `linear`
node prior each leaf's `value` is its intercept and the result gains one
`beta.<column>` column per designated covariate holding that leaf's
slope on the internal standardized scale; internal nodes are `NA`.

For `getSigmas`, a numeric vector of length equal to the number of
chains, giving each chain's current residual standard deviation on the
original response scale.

For `getSumsOfSquaredResiduals`, a numeric vector of length equal to the
number of chains, giving each chain's residual sum of squares \\\sum
(y - \hat{y})^2\\ on the original response scale; a binary-response
sampler reports on the latent scale instead. On a multi-forest
(`forests`) sampler \\\hat{y}\\ is the COMBINED location \\\sum_f
m_f(x_i) f_f(x_i)\\ - the same quantity `run()$train` records - and not
any one forest's own total. Two qualifications, both keyed to the
response family rather than to the forest count. A `logistic` sampler's
working response is \\w_i(y_i - 1/2)/\omega_i - o_i\\ with \\\omega_i\\
a Polya-Gamma variate redrawn every sweep, so what it reports is a
function of that sweep's auxiliary variables and is not a residual sum
of squares in any scale. A multinomial fit reports `NaN`: its coupling
reports a category-by-observation probability slab rather than a
location, so no residual sum of squares is defined - the `NaN` is
deliberate, in place of the sum of squared category probabilities a bare
substitution would have given.

For `getLatents`, `NULL` when the model has no latent-variable
representation (e.g. a gaussian response); otherwise the sampler's
current latent values - a plain vector of length equal to the number of
observations when there is a single chain, or an observations-by-chains
matrix otherwise - written into `result` when one was supplied.

For `getForestFits`, a multi-forest sampler's requested forest's current
internal-scale fitted values, an n.observations x n.chains matrix. For
`getForestAmplitudes`, the named forest's amplitudes - one per basis
column - as a q x n.chains matrix, or, at the default `forest = NULL`,
every forest's stacked forest-major into a sum(q) x n.chains matrix. The
vector is *ragged*: a Bayesian causal forest's forest `1` carries the
single \\a\\ on its implicit intercept and its forest `2` the pair
\\(b_0, b_1)\\, so the stacked read is the \\(a, b_0, b_1)\\ of \\y = a
\mu(x) + b_z \tau(x) + \epsilon\\. For `getForestVariableCounts`, the
requested forest's current per-predictor split counts, an n.predictors x
n.chains integer matrix.

For `getCalibration`, the leaf-prior calibration a forest currently runs
under, as a numeric matrix with one row per chain and the columns
`prior.scale` (the prior standard deviation of the forest total at
`k = 1`, in response units), `prior.sd` (`prior.scale / k`),
`prior.mean`, `k`, `k.has.hyperprior` (1 when this forest's `k` is drawn
every sweep, in which case `prior.sd` moves every sweep while
`prior.scale` does not), `response.scale`, and `response.shift`. The
leaf model rides as a `"leaf.model"` attribute, one of `"constant"`,
`"monotone"`, `"linear"`, or `"gp"`, and qualifies what `prior.sd`
means.

Five further columns report the multi-forest CALIBRATION MAP that fixed
`prior.scale` on a sampler built with `forests =` or
`dbartsData(bases = )` (see
[`forest`](https://vdorie.github.io/dbarts/reference/forest.md)), and
are `NaN` on every forest whose scale that map does not own - any
single-forest sampler, and a multinomial one, whose scale is not
map-derived. `amplitude.prior.variance` and `amplitude.prior.scale` are
EXCLUSIVE per forest: a forest whose amplitudes carry a fixed prior
variance reports that variance and a `NaN` scale, and one whose
amplitude carries the half-Cauchy scale mixture (a forest declaring no
basis) reports its median and a `NaN` variance. Each is a prior the
caller may set; neither moves with the scale mixture's own variance
auxiliary, which is a drawn quantity rather than a prior.
`node.scale.factor` and `node.scale.divisor` are the map's two factors
and `basis.row.norm` the median nonzero row norm of the forest's basis
IN FORCE, which `setForestBasis` re-derives.

Together they decompose the reported scale as
`prior.scale = node.scale.factor * s / (node.scale.divisor * basis.row.norm)`,
so the family's own latent anchor \\s\\ - the only quantity of the map
with no column of its own, and data-dependent under a gaussian
response - is recovered as
`prior.scale * node.scale.divisor * basis.row.norm / node.scale.factor`.
That identity holds whenever `node.scale.factor` is not `NaN`. It
becomes `NaN`, on both `node.scale` columns, when `setState` or
`installTrees` installs a leaf scale differing from the one in force:
the donor's trees arrive with the donor's calibration, and the stored
factor and divisor no longer decompose it. The other three columns are
unaffected - the amplitude prior FOLLOWS the installed state, and the
bases are not state - and a `setForestBasis` on that forest re-imposes
the map and restores both columns. Restoring a sampler's own state
installs a bitwise-identical scale and so changes nothing.

`prior.scale` and `prior.sd` describe the LEAF-PARAMETER scale of the
forest total. They equal the prior standard deviation of \\f(x)\\ at
every \\x\\ for the constant leaf only; for the other three the prior of
\\f(x)\\ is x-dependent and `prior.sd` bounds it in a leaf-specific
direction. Under `"linear"` it is a LOWER bound attained at the
standardized covariate origin, with \\sd(f(x)) = \\ `prior.sd`
\\\sqrt{1 + \\z(x)\\^2}\\ in the internally standardized leaf covariates
(a missing value maps to \\z_j = 0\\); `prior.mean` is exact. Under
`"gp"` it is an UPPER bound over \\x\\, attained at rows that reproduce
a leaf member and on over-cap leaves, and elsewhere decaying to zero as
\\x\\ leaves the leaf's data cloud, where every prior draw equals
`prior.mean` exactly. Under `"monotone"` it is a LOWER bound in the
interior - the realized standard deviation runs a few per cent to about
twenty per cent above it - and `prior.mean` is NOT the prior mean of
\\f(x)\\ under an active constraint: that marginal is skew, with an
x-dependent mean tracking the constraint's direction across several
`prior.sd` (see `monotone` in
[`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md)).

This is the authoritative reader of the calibration in force. A model's
`prior.scale` slot records the named intent and is never rewritten by
the engine, so a channel that re-anchors the response transform -
`setResponse` or `setOffset` at `updateScale = TRUE`, or `setData` -
moves what is in force while leaving the intent alone; `getCalibration`
shows the move, and `setCalibration` or `setModel(sampler$model)`
re-issues the intent.

For `setCalibration`, `NULL` invisibly. The write lands on every chain
and takes effect on the next sweep, reinterpreting no leaf value already
drawn; a write that reproduces what is already in force is skipped
bitwise, so a read followed by a write cannot perturb a draw. It is
total over the four leaf models, each of which carries the one scale it
writes. It is refused on a Bayesian causal forest and on a multinomial
sampler, whose per-forest leaf scales come from their own calibration
maps, and on the host shell of a fit whose model lives elsewhere; a
value that is not a single positive finite number is an error. Nothing
else moves - not `k`, not the response transform, not `sigma`, not the
tree prior, not a DART split prior - which is the difference from
`setModel`, and the reason a DART sampler is served here and refused
there. A heteroscedastic sampler's variance forest is a separate leaf
model and is not addressable.

For `storeState`, `NULL` invisibly; it is called for its side effect of
capturing the sampler's current engine state into the serializable
`state` field, which [`save`](https://rdrr.io/r/base/save.html) then
writes out. See ‘Saving’.

## See also

[`updatePredictorPerObservationJointly`](https://vdorie.github.io/dbarts/reference/updatePredictorPerObservationJointly.md)
for applying a shared per-observation predictor update across several
samplers at once.

[`samplePriorPredictive`](https://vdorie.github.io/dbarts/reference/samplePriorPredictive.md)
for repeated
`sampleTreesFromPrior`/`sampleNodeParametersFromPrior`/`predict` draws
on a private sampler, for calibrating priors before fitting.

## Examples

``` r
# the embedded-Gibbs pattern: BART as a conditional model inside a larger
# sampler, alternating its own draws with updates to its response between
# calls to run()
n <- 100
x <- matrix(runif(n * 2), n, 2)
y <- x[, 1] - x[, 2] + rnorm(n, 0, 0.1)

sampler <- dbarts(
  y ~ x,
  control = dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.burn = 0L,
    n.samples = 1L,
    updateState = TRUE
  )
)

## first draw: response as given at creation
samples <- sampler$run()
str(samples$train) # a plain n.obs x n.samples matrix (n.chains == 1)
#>  num [1:100, 1] -0.8207 0.8956 0.2326 -0.3679 0.0244 ...

## an outer Gibbs step changes the response (e.g. after updating some
## other part of a joint model); the sampler picks the new target up on
## the next run() without being recreated
newY <- y + rnorm(n, 0, 0.05)
sampler$setResponse(newY)
samples <- sampler$run()
str(samples$train)
#>  num [1:100, 1] -0.6833 0.6015 0.3523 -0.4685 -0.0563 ...
```
