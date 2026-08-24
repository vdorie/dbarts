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
predictForests(x.test, offset.test)
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
setCounts(counts, updateState = NA)
# S4 method for class 'dbartsSampler'
setCategoryOffset(offset, updateState = NA)
# S4 method for class 'dbartsSampler'
setCategoryTestOffset(offset.test, updateState = NA)
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
getDispersion()
# S4 method for class 'dbartsSampler'
getLatents(result)
# S4 method for class 'dbartsSampler'
getSumsOfSquaredResiduals(result)
# S4 method for class 'dbartsSampler'
getFitsWithoutOffset()
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
  `useQuantiles`, or `seed` from the values the sampler was created
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
  for `ordinal`, and a finite non-negative integer count no larger than
  \\10^6\\ for `nbinom` (the dispersion grid's count histogram is sized
  from the largest count, so a larger one allocates without bound).
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
  while `predict` and `predictForests` materialize it to a numeric
  matrix before evaluating, retaining nothing.

- offset:

  A numeric vector of length equal to that with which the sampler was
  created, or `NULL`. If `offset.test` was set from `offset`, will
  attempt to update that as well.

  For `setCategoryOffset`, a multinomial sampler's per-category shift
  instead: an \\n \times K\\ numeric matrix, or `NULL` to clear one. The
  latent becomes \\f\_{ik} + o\_{ik}\\, so the shift enters the
  log-sum-exp margins, every category's working response and the
  reported probabilities, and never a leaf value. This is the
  response-side counterpart of `setCounts` rather than of `setOffset`,
  whose flat shift is added AFTER the categories are blended - the wrong
  side of the nonlinearity - and is the softmax's own null direction
  besides, which is why `setOffset` itself is refused on such a sampler.
  Only the row-centred part is identified: adding a constant to a whole
  row leaves every reported probability unchanged, and the entrance
  leaves the matrix as given rather than re-centring it. It shifts the
  TRAINING latent only. Mirrored into `data@offset.category`, so a
  re-created sampler carries it.

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
  otherwise keep reporting \\s^2(x)\\ on the abandoned scale. A sampler
  carrying grouped random effects (see
  [`rbart_vi`](https://vdorie.github.io/dbarts/reference/rbart.md)) is
  refused on the same grounds, but only where there is a data-derived
  transform to abandon: under a gaussian, Student-t
  (`resid.dist = student`) or `"aft"` response the random intercepts
  \\b\\ and their scale \\\tau\\ are held against the transform fixed at
  creation and nothing converts them, so `TRUE` would silently restate
  both in response units while `sigma` moved with the scale; under
  `"probit"` or `"logistic"` the transform is the link's own and `TRUE`
  is accepted as the no-op it already is.

- offset.test:

  A numeric vector of length equal to that of the test matrix, or
  `NULL`. Can be missing for `setTestPredictorAndOffset`. Refused, by
  name, for `predictForests`: an offset shifts the combination of the
  forests, exactly as the response transform's shift does, and neither
  is any one forest's own total.

  For `setCategoryTestOffset`, a multinomial sampler's
  \\n\_{\mathrm{test}} \times K\\ per-category test shift, or `NULL` to
  clear one: the recorded test channel becomes
  \\\mathrm{softmax}(f\_{\mathrm{test}} + o\_{\mathrm{test}})\\, formed
  where the training blend forms its own. The test fits enter no
  likelihood, so this moves the reported test probabilities and nothing
  else. Its rows are the CURRENT test rows, so replacing those rows
  while it is installed is refused rather than silently reinterpreted -
  clear it first - and out-of-sample `predict` does not read it at all,
  taking its own matrix for the rows it is given. Mirrored into
  `data@offset.category.test`. For `predict` on a multinomial sampler,
  likewise a per-category matrix, one row per PREDICTED row: a flat
  vector is refused there for the same reason, and a sampler holding
  either resident category offset refuses an unnamed call rather than
  reporting the offset-free surface, which an explicit all-zero matrix
  asks for on purpose.

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

- counts:

  For `setCounts`, the replacement response of a multinomial (softmax)
  sampler: an \\n \times K\\ matrix of non-negative integer counts,
  column \\k\\ holding category \\k\\'s successes, with trials \\n_i =
  \sum_k\\ `counts[i, ]` at least 1. Both \\n\\ and \\K\\ are fixed at
  creation - every combiner buffer is sized by \\n\\, and \\K\\ is the
  forest count - so only the values may change; a matrix of the wrong
  shape is refused naming the count in force. The trees carry over,
  fitted to the previous counts exactly as `setResponse` leaves a
  single-forest sampler's, and the next `run` forms every category's
  working response against the new matrix. The matrix is written to both
  the engine and `data@counts` (its row sums to `data@y`), so
  `getPointer`'s transparent re-creation after save and load carries the
  current response rather than the one the sampler was created with.
  Cost, not a defect: the sweep draws \\n_i\\ Polya-Gamma variates per
  observation per category, so replacing single-trial labels with
  grouped counts multiplies sweep cost by `mean(n_i)`. Refused, naming
  the reason, on any sampler that carries no count response.

- weights:

  A numeric vector of case weights, of length equal to that with which
  the sampler was created. What a weight MEANS is family-specific: for a
  gaussian-family sampler it is a precision, non-negative, and a weight
  of zero excludes an observation from the likelihood while keeping its
  fitted values; for a `logistic`-family sampler it is an observation
  COUNT and must be a positive integer, a zero count being a dropped row
  rather than a down-weighted one (`active` is the supported way to take
  a row out of the data set for a sweep). Installed BETWEEN samples on a
  forest that has already grown, a vector of zeros can leave leaves that
  no positive-weight row reaches - the trees were drawn before it
  existed, and weights are not part of the saved `state` - which is a
  legal, transient state rather than an error: such a tree keeps moving
  under the tree prior at a constant likelihood and returns to leaves
  the likelihood reaches. `setWeights` applies to gaussian- and
  `logistic`-family samplers. A logistic swap is a model change with a
  defined meaning rather than a reweighting: the counts are the shape of
  the Polya-Gamma latents, which are redrawn against the new counts
  before the call returns - from the sampler's own generators, not R's
  stream - so an outer sampler can vary exposure between runs. `probit`,
  `ordinal`, `aft` and `nbinom` refuse it by identification: a weighted
  probit has no tractable latent-variable form, `ordinal` inherits that,
  `aft` fixes its censoring structure at creation, and `nbinom`'s
  Polya-Gamma shape is \\y_i + r\\ with no weight slot. `setData`
  carries the same rule on the whole-data conduit, and redraws the
  latents against the counts the replacement data carries; replacement
  data given without weights is single-trial, as at creation, so a
  logistic sampler built with counts and handed weightless data becomes
  an unweighted one. Under an installed mask a swap redraws only the
  ACTIVE rows - an inactive row consumes no random numbers and returns
  to its deterministic cold start against the new count. The weights
  themselves are not part of the saved `state`, but a digest of the ones
  in force when it was stored is: `setState` compares it against the
  destination's own weights and, where the two disagree, re-derives the
  weight-dependent latents against the DESTINATION's - so a restore
  lands where the same `setWeights` call would rather than pairing one
  vector's latents with another's counts, silently and off the restored
  generators. A matched round trip re-derives nothing and installs the
  stored state unchanged. Only a family whose augmentation is stated
  against the weights moves under it (`logistic`); for gaussian,
  Student-t and every weight-refusing family it is a no-op. See
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
  equal to that with which the sampler was created, or `NULL`.
  `setActiveRows` reaches every family, the four that refuse case
  weights entirely included: gaussian, Student-t, `probit`, `ordinal`,
  `logistic`, `nbinom`, and `aft` samplers all accept it, including a
  multi-forest (`forests`) sampler, whose forests share one response of
  whichever family it was built under. It is absolute and independent of
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
  per-group sums, and so on). It still occupies its leaf for COUNT-based
  accounting - `numObservations`, the birth scan's own member count, and
  leaf collapsing, which triggers on member count regardless of
  weights - and it still receives a fitted value: `run()$train`,
  `getForestFits`, `getFitsWithoutOffset`, and `predict` report
  \\f(x_i)\\ at an inactive row exactly as at an active one, which is
  what makes this channel worth more than physically dropping the row.
  The empty-leaf VETO is the one exception, and it is conditional: it
  counts POSITIVE-WEIGHT members rather than members, so it degenerates
  to the member count exactly when no weight vector is installed - but a
  mask IS a weight vector, so with one installed a leaf all of whose
  rows are inactive is vetoed rather than counted as occupied. A mask
  installed mid-run on a grown forest can therefore strand a leaf the
  trees were drawn without it: the affected tree keeps moving under the
  tree prior at a constant likelihood and is absorbed back into the
  vetoed-free set as those leaves clear. With every row inactive that is
  the whole forest, which is the exact sense in which it “sits at its
  prior” - a distribution over structures, not a frozen one.

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
  It is DROPPED rather than reconciled, which is where it parts from the
  case weights: a destination that never had a mask has no other mask
  for the restore to re-derive against, only an absent channel, so a
  caller that wants one reinstalls it by hand. Conversely, the latents
  DO ride the state, and under an installed mask they are stale at every
  inactive row for a skipping family, so two samplers at the same
  posterior state but different mask histories carry different stored
  latents and `statesAgree` correctly reports DISAGREEMENT; a Student-t
  sampler is affected the same way, since an inactive row's stored
  \\\lambda_i\\ is drawn from a conditional its own data no longer
  informs.

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

  An integer vector listing the saved samples to return from `getTrees`
  or print from `printTrees`. Applies only when `keepTrees` is `TRUE`
  and `current` is `FALSE`; otherwise the live working trees have no
  sample dimension and it is ignored. Sample numbers address recorded
  DRAWS on the oldest-first axis `predict` reports - `1` to the number
  of draws recorded, at most `control@n.samples` - and not the store's
  internal slots.

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
  materialized. The saved-tree store's write position and the number of
  draws it has recorded both ride it, so `predict` after a `setState`
  reports the same draws, in the same order, as before the store.
  Reading it forces the sampler's *current* state; `storeState`
  refreshes it on demand and `updateState` governs when the methods do
  so themselves. It is the only field
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

### Reading the fit

Six methods report a fitted quantity and no two report the same one.
They differ along three axes - the SCALE the values live on, whether the
installed OFFSET is in them, and whether they answer from the CURRENT
state or from stored draws:

|  |  |  |  |
|----|----|----|----|
| **method** | **scale** | **offset** | **current or stored** |
| `run()$train` | response (latent under a binary family) | included | one slab per recorded draw |
| `$getFitsWithoutOffset()` | response (latent under a binary family) | excluded | current |
| `$getForestFits(f)` | internal, and one forest only | excluded | current |
| `$predict(x)` | response (latent under a binary family) | whatever you pass it | the saved samples under `keepTrees`, else the current trees |
| `$predictForests(x)` | internal, one channel per forest | excluded, and refused if passed | the saved samples under `keepTrees`, else the current trees |
| `$getLatents()` | the family's own augmentation variable, which is a location for some families and a precision for others | not a fit; see ‘Value’ | current |

A multinomial (softmax) sampler reads and writes a strict subset of this
surface, by the model rather than by omission. Its response channels are
`$setCounts`, `$setCategoryOffset` and `$setCategoryTestOffset`; the
predictor family (`$setPredictor` in all four shapes, `$setCutPoints`,
`$setTestPredictor`) and the global `$setActiveRows` are open, and
`$predict` answers with the \\K\\ category probabilities and takes the
matching per-category offset. Everything else on the mutation surface is
refused, each by a message naming the capability and, where one exists,
the channel that serves the caller instead: `$setResponse` and
`$setOffset` (a flat vector cannot express the matrix, and a flat shift
is the softmax's null direction), `$setWeights` (an integer case weight
is already row-wise replication in the count response, and a non-integer
one has no exact augmentation sampler), `$setSigma` (no residual scale),
`$setData` and `$setModel` (the K category forests fix their data and
their calibration at creation), `$setCalibration` (the softmax map owns
every category forest's leaf scale), `$setForestWeights` and
`$setForestBasis` (its forests are its categories, which carry no
amplitudes), and `$getFitsWithoutOffset` (its reported channels are
probabilities, not one additive location). Per-forest masking is refused
permanently rather than pending: a category's margin is a log-sum-exp
over the other \\K-1\\, so a mask on one forest alone has no conditional
to be a mask of.

The composition rule follows from the table: an outer block conditions
on \\f(x_i)\\, so it reads `$getFitsWithoutOffset()` and adds back
whatever offset it installed, rather than differencing `$getLatents()`
against `run()$train`, which mixes an offset-bearing quantity with one
that does not carry the offset.

Where the boundary runs. These are the methods of the mutable sampler,
and their audience is a model writer composing BART into a larger
scheme: they return engine primitives on the engine's own terms. The R
modelling conventions -
[`fitted`](https://rdrr.io/r/stats/fitted.values.html),
[`predict`](https://rdrr.io/r/stats/predict.html),
[`extract`](https://vdorie.github.io/dbarts/reference/bart.md) - apply
to FIT objects instead, the results of
[`bart`](https://vdorie.github.io/dbarts/reference/bart.md),
[`bart2`](https://vdorie.github.io/dbarts/reference/bart2.md) and
[`rbart_vi`](https://vdorie.github.io/dbarts/reference/rbart.md). Those
accessors read the stored `yhat.train`/`yhat.test` channels the engine
already wrote the offset into, and no `type` arm removes it:
`type = "bart"` returns those draws as they stand, `"ev"` maps them
through the response transform, `"ppd"` samples from them, and
`"loglik"` evaluates against them - all offset-inclusive. (`extract` on
a `dbartsSampler` itself takes `type = "predictors"` only and returns
the design matrix; it is not a fitted-value accessor.) So an offset-free
fit is available from the sampler surface and nowhere on the fit-object
surface, which is the division rather than an omission.

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
the `state` field. Validation covers every forest the state carries: a
state whose trees split outside the recipient forest's allowed columns -
a `blocks`-constrained or moderator-restricted mean forest, or a
`variance = ~ x1 + x2` variance forest - is refused with the message
`installTrees` gives the same donor, the two entries sharing one rule so
neither admits what the other refuses. Every check runs before any live
state is touched, so a refused restore leaves the sampler exactly as it
was. Assigning the field directly (`sampler$state <- newState`) does
*not* restore the sampler - it only overwrites the R-side cache, leaving
the engine untouched, so the next run continues from the engine's own
state rather than the assigned one. Always route a restore through
`setState`. The case weights are not in the state, so a restore is
reconciled against the DESTINATION's own rather than governed by the
source's: where they differ from the weights the state was stored under,
the weight-dependent latents are re-derived against the destination's
before `setState` returns (see `weights` above).

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
whole-matrix `setPredictor(x, ...)` (`column` missing) carries column
names onto the installed `sampler$data@x`. When `x` itself has
`dimnames`, those names are installed; when it does not (a bare matrix,
or a vector reshaped through `dim`), the sampler's existing column names
are kept instead. Either way, a later *by-name* column reference on that
sampler - `setPredictor`'s own `column` argument given as a character
string, or
[`updatePredictorPerObservationJointly`](https://vdorie.github.io/dbarts/reference/updatePredictorPerObservationJointly.md)'s
shared-column match - continues to resolve against
`colnames(sampler$data@x)` as it did before the replacement.

### Grouped random effects

A sampler carrying grouped random intercepts - the one
[`rbart_vi`](https://vdorie.github.io/dbarts/reference/rbart.md) returns
in `$fit` when its prior is built in - is mutable on the response side.
`setResponse` and `setOffset` accept a same-length replacement at the
pinned scale (`updateScale = FALSE`, the default): the group effects and
their scale are a Gibbs block the swap deliberately carries across,
exactly as they are carried across a tree sweep, and the group indices
are per-observation and unchanged, so every row stays in its group.
`updateScale = TRUE` is refused under a re-anchoring family (see
`updateScale` above). `setData` stays refused outright: it may change
the number of rows, and the grouping is fixed at creation, so there is
no coherent reading of the new rows' group membership - re-create the
sampler instead.

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
plain vector of length n.samples. On a multi-forest (`forests`) sampler
`varcount` gains a forest axis between the predictors and the samples -
n.predictors x n.forests x n.samples x n.chains, forest-major within a
draw and the prognostic forest first - so each forest's own per-draw
split counts arrive from one call rather than only the reported
forest's; a single-forest sampler's array keeps exactly its n.predictors
x n.samples x n.chains shape, and the same widening is what a
`family = "multinomial"` sampler's per-category counts ride. This is the
RAW run shape; the packaged fit
[`bart2`](https://vdorie.github.io/dbarts/reference/bart.md) builds
reshapes it draws-first with the forest names on the trailing margin, as
it does for the multinomial channel. `$getForestVariableCounts` reads
the same quantity for the CURRENT state, one forest at a time. A
Bayesian causal forest adds two more: `forestFits`, an n.obs x n.forests
x n.samples x n.chains array of each forest's fitted values on the
internal scale (the prognostic \\\mu\\ first, the treatment \\\tau\\
second), and `glue`, a sum(q) x n.samples x n.chains array of the
amplitudes each draw combines them through, stacked forest-major and as
wide as each forest's own basis (a Bayesian causal forest's \\(a, b_0,
b_1)\\, three rows), so both surfaces and their recombination come from
a single call; `train` carries the combination \\a \mu(x_i) + b\_{z_i}
\tau(x_i)\\ on the response scale, or on the latent scale when the
sampler was built under `"probit"` or `"logistic"`, and `test` is filled
with `NaN` (there is no test treatment vector to combine off-sample). No
other model reports either element. A `"nbinom"` sampler adds one of its
own: `dispersion`, the negative-binomial \\r\\ each draw is conditioned
on, shaped exactly as `sigma` (a length n.samples vector at one chain,
an n.samples x n.chains matrix otherwise) because it is the count analog
of it - fixed at the value the sampler was created with under a fixed
`dispersion`, and that sweep's grid draw otherwise. It is written from
the same state `storeState` serializes and consumes no random numbers,
so reading it costs a run nothing. No other family carries the element
at all: it is absent from the list, not `NULL` within it, so
`run()$dispersion` is `NULL` on every non-`"nbinom"` sampler and a test
of the channel must be `!is.null(...)` rather than a comparison, which
`NULL` would satisfy vacuously. A sampler built with
`resid.dist = student()` adds `resid.df` on exactly the same terms - the
degrees of freedom \\\nu\\ each draw is conditioned on, shaped as
`sigma`, written from settled state and consuming no random numbers,
absent from the list under any other error law - fixed at the value
supplied to `student(df = )`, and that sweep's grid draw when the
degrees of freedom are estimated. A run can be interrupted with
`Ctrl-C`: it stops between iterations - joining any worker threads
first - and signals an error, returning no samples from the interrupted
run. The sampler's chains are left at the iteration they reached, which
is a valid state to run again from.

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

Under `keepTrees` the draws `predict`, `predictForests`, `getTrees` and
`printTrees` report come out OLDEST FIRST: the store keeps the most
recent `control@n.samples` recorded draws, so a sampler driven through
several recorded `run` calls replays them in the order they were drawn
and lines up draw for draw with the `train`/`sigma` channels those runs
returned. A sampler that has recorded fewer draws than that reports only
the draws it holds, and one that has recorded none - `keepTrees` set but
nothing run since creation, a `setControl` that resized the store, or an
`installTrees` warm start - is refused rather than answered from slots
nothing has written.

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
depth-first, left-hand-side pre-order, with columns `chain` (present
only when `n.chains > 1`), `sample` (present only for saved samples, and
reporting the draw number asked for), `tree`, `n` (the number of
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

For `getDispersion`, the negative-binomial dispersion \\r\\ currently in
force, a numeric vector of length equal to the number of chains, and
`NULL` on every other family - the count analog of `getSigmas`. It is
the same scalar `run()$dispersion` records once per kept draw, read
mid-sweep and without serializing state, so a host driving the sampler
one sweep at a time reads it here rather than through `storeState()` and
`state[[chain]]$dispersion`. Because the refusal for a family carrying
no dispersion is a `NULL` and not an error, a caller distinguishing “no
dispersion” from a value must test `!is.null(...)`; comparing to an
expected number would pass vacuously.

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

For `getLatents`, `NULL` when the model augments nothing - a plain
gaussian response - and, deliberately, on a multinomial one, which
augments but reports nothing: its \\K \times n\\ Polya-Gamma \\\omega\\
matrix is drawn every sweep, but each \\\omega\_{ik}\\ is an interleaved
ONE-VS-REST augmentation variable refreshed against a margin that the
other \\K-1\\ forests move within the same sweep, so a host reading it
between sweeps reads a quantity whose conditioning it cannot reproduce,
and the composition recipes reach the same posterior through the fits
and the category offset without touching it. That decline is scoped to
the augmentation families whose latent is a PRECISION refreshed against
a moving margin; the location-valued latents below - probit's and
ordinal's \\z_i\\, and aft's imputed log survival time - are reported as
documented. Otherwise the sampler's current draw of the augmentation
variable, a plain vector of length equal to the number of observations
when there is a single chain, or an observations-by-chains matrix
otherwise, written into `result` when one was supplied. **What that
variable IS depends on the family and is not uniform.** It is a
LOCATION, on the sampler's own latent scale, for `"probit"` (the
truncated normal \\z_i\\), `"ordinal"` (the same \\z_i\\ under the cut
points) and `"aft"` (the imputed log survival time), all of which a host
regresses on directly. It is a PRECISION, one per observation, for
`"logistic"` and `"nbinom"` (the Polya-Gamma \\\omega_i\\) and for a
Student-t residual distribution (the scale-mixing \\\lambda_i\\); these
WEIGHT a working response and are not on the response scale at all, so
differencing them against a fit is meaningless. Note the last case: a
sampler whose `family` is `"gaussian"` but whose `resid.dist` is
`student()` (see
[`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md)) does
report latents, and they are precisions.

For `getFitsWithoutOffset`, the sampler's combined per-observation
location on the response scale and *without* the installed offset, an
n.observations x n.chains matrix. `run()$train` reports the same
quantity with the offset folded in, so this value plus the installed
offset is that one; there is deliberately no with-offset twin, since the
offset is the caller's own. It is the read a host driving the sampler
one sweep at a time needs: an incremental scheme updates its outer block
against \\f(x_i)\\, and taking `getLatents()` minus `run()$train`
instead is BIASED, because the training channel carries the offset the
latent scale does not. On a multi-forest (`forests`) sampler this is the
only route from R to the combined fit at all - `getForestFits` gives one
forest's internal-scale totals and `predict` is refused. It is refused,
naming the reason, on a multinomial sampler, whose reported channels are
per-category softmax probabilities rather than one additive location;
`predict(x)` serves that read, but reports the SAVED samples rather than
the current state when the sampler was built with `keepTrees`.

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
n.chains integer matrix; its rows are named by the training matrix's
column names when that matrix carries them, and left unnamed otherwise.

For `predictForests`, the off-sample twin of `getForestFits`: every
forest replayed separately at `x.test`, as an n.test x n.forests x
n.samples (x n.chains) array on the forests' internal scale,
forest-major within a draw. It reports the saved samples when the
sampler was built with `keepTrees`, and the current trees otherwise,
exactly as `predict` does. Nothing but the forests' own totals is in it:
no amplitude glue, no response transform and no offset, the last refused
rather than ignored. That is the point of the channel rather than an
omission - an amplitude coupling's location is `response.shift` +
\\\sum_f (B_f a_f) \times (\mathrm{response.scale} \times f_f)\\, and
off the training rows only the caller knows the bases \\B_f\\, so the
whole recombination is theirs, with `$getForestAmplitudes()` and
`$getCalibration()` supplying the rest. At the fit level it is packaged:
`predict` on a
[`bart`](https://vdorie.github.io/dbarts/reference/bart.md) object with
`type = "ev"`, `"ppd"` or `"bart"` performs the same recombination off
this replay, taking the bases at the predicted rows through its own
`bases` argument. Only a sampler that composes its forests through
scalar amplitude glue reports per-forest fits; every other one is
refused by name, a multinomial sampler included, whose raw per-category
totals are perfectly well defined but whose reported quantity is a
softmax probability that `predict` serves. The refusal is off-sample
only: `getForestFits` is not gated the same way, so those same
per-category totals stay readable AT the training rows, for the current
trees and one draw, and it is the off-sample replay - which would have
to hand back a level nothing identifies at rows the sampler never saw -
that has no counterpart.

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
maps; a value that is not a single positive finite number is an error.
Nothing else moves - not `k`, not the response transform, not `sigma`,
not the tree prior, not a DART split prior - which is the difference
from `setModel`, and the reason a DART sampler is served here and
refused there. A heteroscedastic sampler's variance forest is a separate
leaf model and is not addressable.

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

[`dbarts-embedding`](https://vdorie.github.io/dbarts/reference/dbarts-embedding.md)
for which channel an outer block writes, which fit it reads back, and
whose prior is in force after a swap, and
[`vignette("dbarts-as-a-component", package = "dbarts")`](https://vdorie.github.io/dbarts/articles/dbarts-as-a-component.md)
for six worked compositions built from these methods.

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
