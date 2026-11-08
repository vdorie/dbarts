# Class "dbartsSampler" of Discrete Bayesian Additive Regression Trees Sampler

A reference class object that contains a Bayesian Additive Regression
Trees sampler in such a way that it can be modified, stopped, and
started all while maintaining its own state.

## Usage

``` r
# S4 method for class 'dbartsSampler'
run(
  numBurnIn = NULL, numSamples = NULL, updateState = NULL, ..., callback = NULL
)
# S4 method for class 'dbartsSampler'
sampleTreesFromPrior(updateState = NULL)
# S4 method for class 'dbartsSampler'
sampleLeafParametersFromPrior(updateState = NULL)
# S4 method for class 'dbartsSampler'
sampleVarianceForestFromPrior(updateState = NULL)
# S4 method for class 'dbartsSampler'
growFromRoot(n.sweeps = 2L, updateState = NULL)
# S4 method for class 'dbartsSampler'
copy(shallow = FALSE)
# S4 method for class 'dbartsSampler'
show()
# S4 method for class 'dbartsSampler'
predict(x.test, offset.test, n.threads = control@n.threads)
# S4 method for class 'dbartsSampler'
predictForests(
  x.test, offset.test, n.threads = control@n.threads
)
# S4 method for class 'dbartsSampler'
setControl(newControl)
# S4 method for class 'dbartsSampler'
setModel(newModel)
# S4 method for class 'dbartsSampler'
setData(
  newData,
  updateState = NULL,
  updateCutPoints = NULL,
  updateScale = NULL,
  updateStandardization = NULL
)
# S4 method for class 'dbartsSampler'
setResponse(
  y, updateScale = FALSE, updateState = NULL, status = NULL
)
# S4 method for class 'dbartsSampler'
setOffset(offset, updateScale = FALSE, updateState = NULL)
# S4 method for class 'dbartsSampler'
setWeights(
  weights, updateState = NULL, updateScale = FALSE,
  updateStandardization = FALSE
)
# S4 method for class 'dbartsSampler'
setCategoryOffset(offset, updateState = NULL)
# S4 method for class 'dbartsSampler'
setCategoryTestOffset(offset.test)
# S4 method for class 'dbartsSampler'
setActiveRows(
  active, updateState = NULL, updateScale = FALSE,
  updateStandardization = FALSE
)
# S4 method for class 'dbartsSampler'
setForestWeights(forest, weights, updateState = NULL)
# S4 method for class 'dbartsSampler'
setForestBasis(forest, basis, updateState = NULL,
    updateBasisScale = FALSE)
# S4 method for class 'dbartsSampler'
setSigma(sigma, updateState = NULL)
# S4 method for class 'dbartsSampler'
setPredictor(
  x,
  column,
  forceUpdate,
  updateCutPoints = "none",
  updateState = NULL,
  updateStandardization = FALSE
)
# S4 method for class 'dbartsSampler'
setCutPoints(cuts, column, updateState = NULL, splits = "position")
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
  treeNums, chainNums, sampleNums, current = FALSE, newdata = NULL,
  forest = NULL
)
# S4 method for class 'dbartsSampler'
getSigmas(result)
# S4 method for class 'dbartsSampler'
getShape()
# S4 method for class 'dbartsSampler'
getLatents(result)
# S4 method for class 'dbartsSampler'
getSumsOfSquaredResiduals(result)
# S4 method for class 'dbartsSampler'
getFitsWithoutOffset()
# S4 method for class 'dbartsSampler'
getVariance(test = FALSE)
# S4 method for class 'dbartsSampler'
getForestFits(forest = NULL)
# S4 method for class 'dbartsSampler'
getForestAmplitudes(forest = NULL)
# S4 method for class 'dbartsSampler'
getForestVariableCounts(forest = NULL)
# S4 method for class 'dbartsSampler'
getLeafPrior(forest = NULL)
# S4 method for class 'dbartsSampler'
getK(forest = NULL)
# S4 method for class 'dbartsSampler'
setLeafPrior(leaf.prior, updateState = NULL)
# S4 method for class 'dbartsSampler'
getForest(forest = NULL)
# S4 method for class 'dbartsSampler'
setForest(forest = NULL, sd = NULL, power = NULL, base = NULL,
    amplitude = NULL, proposal.probs = NULL, n.perturb.cuts = NULL, treeShift = NULL,
    updateState = NULL, n.trees = NULL, vars = NULL, interactions = NULL,
    blocks = NULL, basis = NULL, ...)
# S4 method for class 'dbartsSampler'
setTreePrior(power = NULL, base = NULL, split.probs = NULL)
# S4 method for class 'dbartsSampler'
getVarianceForest()
# S4 method for class 'dbartsSampler'
setVarianceForest(power = NULL, base = NULL,
    proposal.probs = NULL, n.perturb.cuts = NULL, updateState = NULL, n.trees = NULL,
    vars = NULL, ...)
# S4 method for class 'dbartsSampler'
installTrees(donor, samples = NULL, forceUpdate = FALSE)
# S4 method for class 'dbartsSampler'
sampleNodeParametersFromPrior(updateState = NULL)
# S4 method for class 'dbartsSampler'
startThreads(n.threads = control@n.threads)
# S4 method for class 'dbartsSampler'
stopThreads()
# S4 method for class 'dbartsSampler'
storeState(ptr = getPointer())
# S4 method for class 'dbartsSampler'
setState(newState, forceUpdate = FALSE)
# S4 method for class 'dbartsSampler'
plotTree(
  treeNum, chainNum, sampleNum, forest = NULL, treePlotPars = c(
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
  should skip before storing results. If `NULL` (the default), it is
  filled in from the sampler's
  [`control`](https://vdorie.github.io/dbarts/reference/dbartsControl.md)
  object; an explicit `NA` reads as `NULL` for one release, with a
  once-per-session warning.

- numSamples:

  A positive integer determining how many posterior samples should be
  returned. If `NULL` (the default), it is also filled in from the
  control object; an explicit `NA` reads as `NULL` for one release, with
  a once-per-session warning.

- callback:

  `NULL` (the default) or a list with two elements, `fn` and `context`,
  each an R external pointer (`externalptr`): `fn` the address of a
  compiled `dbarts_draw_callback` (`inst/include/dbarts/dbarts.h`) and
  `context` an opaque pointer handed back to it unchanged, or `NULL`. It
  fires once per SAVED draw of this call, per chain, on the thread that
  owns that chain, over const pointers into this run's own result
  buffers or, for a channel this run does not keep, into a per-chain
  one-draw scratch buffer the chain's next draw overwrites: the
  per-observation channels under `control@keepFits == FALSE`, and always
  the draw's log-likelihood, which this run never returns and computes
  only for the callback. A nonzero return aborts the run. A draw's sigma
  is `NaN` on a sampler with a variance forest, which reports the
  variance channels instead. Only that both list elements are external
  pointers is checked - the address is dereferenced exactly as handed,
  so a malformed pair or a stale pointer crashes the session with no
  condition to catch. Registering one is a compiled task (Rcpp or plain
  C), not an R closure: see
  [`vignette("dbarts-as-a-component")`](https://vdorie.github.io/dbarts/articles/dbarts-as-a-component.md).
  **Three warnings.** A callback that crashes takes the whole session
  down with it. An interrupt (`Ctrl-C`) cannot land while a call is
  running, so a callback that blocks or loops hangs the session. And
  under `control@keepFits == FALSE` the per-observation channels this
  call returns are `NULL` - present in the returned list, holding
  nothing - so whatever the callback itself accumulated is the only
  record of the run's per-observation fits. On a
  `family = "hazard"`/`"hazard.probit"`/`"hazard.logistic"` sampler the
  draw struct's `numObservations` counts the person-period-EXPANDED
  rows, one per subject per at-risk period, not the original subjects;
  the row-to-(subject, period) mapping is built in R and is not carried
  across to the callback. A hook set through the flat C entry
  `dbarts_sampler_setDrawCallback` does NOT fire on this route:
  `bartcore_run` (which this method calls) builds its own per-run hook
  from THIS argument, so the two channels are independent - a C caller
  mixing both must register the same callback on each to see it fire on
  both.

- updateState:

  A logical determining if the local cache of the sampler's state should
  be updated after the call completes. One convention applies to every
  method that takes this argument - `run`, `sampleTreesFromPrior`,
  `sampleLeafParametersFromPrior`, `sampleVarianceForestFromPrior`,
  `growFromRoot`, and every mutator (`setData`, `setResponse`,
  `setOffset`, `setWeights`, `setCategoryOffset`, `setActiveRows`,
  `setForestWeights`, `setForestBasis`, `setSigma`, `setLeafPrior`,
  `setPredictor`, and `setCutPoints`): `NULL` (the default) resolves to
  the sampler's
  [`control`](https://vdorie.github.io/dbarts/reference/dbartsControl.md)
  object's `updateState`, and an explicit `TRUE` or `FALSE`, or a single
  number read as `as.logical` reads it (0 is `FALSE`), overrides it.
  Anything else is an error, raised before the sampler changes; an
  explicit `NA`, the default of 0.9-x, reads as `NULL` for one release,
  with a once-per-session warning. The mutators are typically called
  once per sweep inside a larger Gibbs/MH loop (as `dbartsSampler` is
  designed for); a control with `updateState = FALSE` keeps such a loop
  free of the store's cost on every call, while the default `TRUE` keeps
  `state` current after each one. Under `control@updateState = TRUE` the
  first read of a `state` no store has yet filled stores the sampler's
  *current* state, so a mutate-then-first-read sequence needs no
  explicit store; under `FALSE` that first read is `NULL`, and
  `$storeState()` makes the store. Pass `TRUE` explicitly to force a
  store under a `control@updateState = FALSE` sampler, or `FALSE` to
  skip one under `control@updateState = TRUE`.

- shallow:

  A logical determining if the copy should retain the underlying data of
  the sampler (`TRUE`) or have its own copies (`FALSE`). Either way, the
  copy keeps the original's model and its data object, on whose record
  of the cut points, their weights and the response range it is built
  (see
  [`dbartsData`](https://vdorie.github.io/dbarts/reference/dbartsData.md)).
  A copy of a sampler whose engine exists in this R session takes that
  engine's current chain, as stored, converting nothing, whatever the
  original's `state` field holds, and the copy's `state` field is that
  chain; the original's own field and draws are untouched. A copy of one
  whose engine is gone, as after a
  [`load`](https://rdrr.io/r/base/load.html), installs its stored state
  as a reload does, and without one is refused, as a reload is (see
  ‘Saving’). Either way it continues the same generator streams as the
  original - not an independent chain - typically matching an
  uninterrupted run to the last few digits, not bitwise; see the
  Reproducibility section of
  [bart](https://vdorie.github.io/dbarts/reference/bart.md). When
  `n.samples` was changed after the chain was stored, the copy's
  saved-tree store takes the chain's capacity, as a reload does, and
  keeps the saved trees it holds rather than truncating or resetting
  them.

- newControl:

  An object inheriting from
  [`dbartsControl`](https://vdorie.github.io/dbarts/reference/dbartsControl.md).
  When passed to `setControl`, it cannot change `n.chains` or `seed`
  from the values the sampler was created with, a `n.trees`,
  `useQuantiles` or `n.cuts` its caller wrote cannot differ from the
  sampler's, and it cannot set `keepTrees = TRUE` without also giving
  `n.samples`; each is an error. A control taken from another fit that
  holds another of those three is taken with a warning, and a fresh
  control silently, the sampler keeping its own. Under
  `treeShift = "auto"` the step follows the mixture, which `setControl`
  may freeze, on every forest that states no mixture of its own.

- newModel:

  An object inheriting from `dbartsModel`. When passed to `setModel`, it
  changes the numbers inside the sampler's model and is refused, by name
  and before anything is stored, where it differs in structure: the
  response family (a Student-t error against a Gaussian one included),
  the kind of leaf and its covariate columns, a Gaussian-process leaf's
  `max.leaf.size`, the interaction limits and blocks, DART against the
  ordinary tree prior, the columns a forest may split on (its
  [`forest`](https://vdorie.github.io/dbarts/reference/forest.md)
  `vars`), the tree count of a forest or the variance forest, and the
  monotone directions; make a new sampler instead. A count the model
  does not state, as in `dbartsSpec(data, power = 3)$model`, leaves the
  sampler's. It takes `power` and `base`, the split probabilities, the
  leaf scale as `k` or `sd`, the residual prior, the Student-t degrees
  of freedom, a negative-binomial shape, a Gaussian-process lengthscale
  and DART's own settings. A degrees of freedom, shape or DART `alpha`
  handed from held to drawn starts at the value in force, not where a
  new sampler would start it; a new lengthscale is refused on a sampler
  holding kept draws, which replay only under the kernel they were drawn
  with. A model that states no family, degrees of freedom or lengthscale
  keeps the sampler's. A split probability of zero on a column some
  chain's trees split on is refused, naming the column. A residual prior
  other than the sampler's own is refused on a family with no residual
  scale (probit, logistic, ordinal, negative binomial, hazard), a
  `chisq()` residual prior of `df` 2 or below on a sampler whose
  variance forest has two or more trees, each tree's prior being
  calibrated to the prior mean of \\\sigma^2\\, which is infinite there,
  and so is a model whose leaf prior's `k` or `sd` does not resolve to
  the leaf hyperprior and prior scale it carries: `setLeafPrior`
  restates a leaf prior. The sampler's own model, as after `setResponse`
  or `setOffset`, is always accepted. On a sampler of several forests or
  a multinomial one it takes the fitting function's level, as
  `setTreePrior` and `setLeafPrior` write it: a changed `power` or
  `base` reaches the forests with no basis that state none of their own,
  and, where one forest has no basis, or on a multinomial sampler's
  category forests, those forests whatever they state; the leaf prior
  goes where `setLeafPrior` sends it; and the residual prior is
  installed as on one forest. A forest's own value that differs, its
  coefficient's prior included, or a change in which values it states,
  is refused, naming `setForest`, which changes one forest; so is a leaf
  scale other than the sampler's, and split probabilities other than
  equal ones. A variance forest stays as it is, whatever the model's
  record of it says; `setVarianceForest` changes it. A new fixed
  negative-binomial shape keeps the latents in force. This draws
  nothing; the next sweep draws the augmentation variables against the
  change.

- newData:

  An object inheriting from
  [`dbartsData`](https://vdorie.github.io/dbarts/reference/dbartsData.md).
  A sampler created with one forest and none refuses data that carry
  `bases`, having no forest to multiply them; make a new sampler from
  such data with
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md).
  Columns are matched by name: the new data carry the sampler's names,
  in any order, and names that are not the sampler's are refused, as
  `predict` refuses them, leaving the sampler as it was. A matrix with
  no column names is read by position, with no message, as is any data
  handed to a sampler that has none, whose record then takes the data's
  names, which `predict` needs; a sampler whose names repeat takes the
  same names in the same places only. A factor column is read by label
  against the levels the sampler was created with, which it keeps: fewer
  levels, the levels in another order and declared levels no row has are
  taken, each row at its own label; a row, training or test, at a label
  the sampler never had is refused by name, as `predict` refuses one,
  and so is a factor column given as numbers with no level table. A
  factor whose labels were changed on purpose is relabeled to the
  sampler's labels first, with `levels<-`, or given to a new sampler. A
  column that could hold missing values before the call still can, and
  one the new data first give a missing value becomes able to, as
  ‘Missing values in predictors’ in ‘Details’ describes. The cut points,
  the response range and a `linear` or `gp` leaf's covariate
  standardization follow `updateCutPoints`, `updateScale` and
  `updateStandardization`; with each left at `NULL`, a value `newData`
  records (its `cutPoints`, `response.range` and `standardization`
  slots, see
  [`dbartsData`](https://vdorie.github.io/dbarts/reference/dbartsData.md))
  is used and one it does not record is derived from it. A cut count or
  quantile rule `newData` states wins over the sampler's (see `n.cuts`
  and `useQuantiles` on
  [`dbartsData`](https://vdorie.github.io/dbarts/reference/dbartsData.md)),
  and the sampler derives by what it takes from then on, in this call
  and in any later refresh; an object carrying `cutPoints` brings the
  count and rule that made them, an unstated rule reading as even; an
  object stating neither keeps the sampler's, a count stated to
  [`dbartsData()`](https://vdorie.github.io/dbarts/reference/dbartsData.md)
  is read over the predictors the sampler holds, a hazard sampler's
  `period` column among them, which takes a count only by its own name
  and is otherwise 100, and `updateCutPoints = "none"` keeps the
  sampler's whatever the object states. The sampler writes its values
  into its own data object, so `sampler$setData(sampler$data)`, the
  object handed back whole, keeps the cut points, the range and the
  standardization the sampler held, where 0.9-x derived them again; to
  change a column and keep them,
  `sampler$setPredictor(newValues, 1, forceUpdate = TRUE)` takes the
  values and keeps the cut points and the record; and a data object made
  by
  [`dbartsData`](https://vdorie.github.io/dbarts/reference/dbartsData.md)
  records none and derives every one, as before. Draws the sampler has
  saved (`keepTrees`) stay the functions they were, so `predict` returns
  the same values for them before and after the call. To that end a
  `linear` leaf's coefficients, saved and current, are rewritten into
  the center and scale the call derives or takes for its leaf
  covariates, and every saved draw into the response scale it derives or
  takes; a value kept, or taken equal to the one in force, rewrites
  nothing. A row whose leaf covariate is missing is read at that
  covariate's center, so what a linear leaf predicts for it follows the
  center to the new data's. A standardization the call derives is read
  over the new data's rows in the likelihood, those of positive weight;
  where fewer than two rows are in (subjects, for a hazard fit, whose
  data hold each subject's rows together, at periods 1, 2, 3 and on in
  order, as
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md)
  expands them, data laid out otherwise being refused), or the call
  reads the response's scale from every row because those rows hold
  fewer than two distinct values of it, the leaf covariates are read
  from every row, under one warning; where the rows read hold one value
  of a leaf covariate, or observe none, the call is refused before
  anything changes, unless a `gp` leaf states its `lengthscale`, and
  `updateStandardization = FALSE` takes the data under the sampler's
  own. Under such a family the sampler's own data object carries no
  weights, so `setData(sampler$data)` drops the mask in force: hand the
  0/1 weights back on the data object, or call `setActiveRows` after. A
  leaf covariate that holds a single value is centered at it and divided
  by 1, so a leaf reads zero for it on every row and its slope adds
  nothing to the fit. Where a covariate is such a column on one side of
  the call only, the current coefficients are still rewritten by the
  formula, with 1 for that side's scale, so each leaf keeps its function
  of the covariate: onto such a column it takes the value its function
  has at the new constant, and off one the slope no row informed is
  carried onto the new data's spread, so the fit may move until the
  sampler has run on them. Where it is such a column on both sides the
  current coefficients are left as they are. A column that
  `setPredictor` has since given other values is read as what it then
  holds: one moved to another single value no longer sits at its center
  and is rewritten as any other column is, by the formula, which keeps
  the function of the covariate. The data fixed that function at one
  point only, so new data that spread the covariate widely can carry the
  current fit far from where it was until the sampler has run on them.
  Saved draws are rewritten by the formula in every case, with 1 for the
  scale. The saved draws of a Gaussian-process (`gp`) leaf cannot be
  rewritten: a sampler with such leaves that holds saved draws refuses a
  call that would derive the standardization or the range, as one given
  a data object that records neither does, or that would take a recorded
  one other than its own; it takes one that keeps both, as
  `setData(sampler$data)` or
  `updateStandardization = FALSE, updateScale = FALSE` do, whatever cut
  points it is given, a saved draw holding its thresholds by value. Pass
  those two, or make a new sampler. A family with no response range
  (`probit`, `ordinal`, `logistic`) has none to change. One that holds
  no saved draw takes any call as any other sampler does. The ordinal
  cutpoints, the Student-t degrees of freedom and the negative-binomial
  shape the sampler has drawn are kept, as are its trees, and the
  latents start as creation starts them, the number of rows possibly
  differing (an ordinal or negative-binomial sampler from the kept
  cutpoints or shape); the next sweep draws them from the kept values.
  This draws nothing; the next sweep draws the augmentation variables
  against the change.

- y:

  A numeric response vector, or on a count sampler the matrix below, of
  length equal to that with which the sampler was created, and lying in
  the response family's support: 0/1 for `probit` and `logistic`, an
  integer category index in \\\[1, K\]\\ for `ordinal`, and a finite
  non-negative integer count no larger than \\10^6\\ for `nbinom` (the
  shape grid's count histogram is sized from the largest count, so a
  larger one allocates without bound). Values off the support are
  refused, as they are at creation; `gaussian` and `aft` (log survival
  times) constrain nothing. A missing value is refused in every family:
  to leave a row out of the likelihood mid-chain, mark it inactive with
  `setActiveRows` and give it any value in the support, which enters no
  draw while the row stays inactive unless the call re-anchors the scale
  (`updateScale = TRUE`), which reads the rows in the likelihood. On a
  multinomial (softmax) sampler `y` is its response as its formula
  writes it: an \\n \times K\\ matrix of non-negative integer counts,
  column \\k\\ holding category \\k\\'s successes, with trials \\n_i =
  \sum_k\\ `y[i, ]` at least 0 and no larger than \\10^6\\: a row with
  no trial enters no likelihood and still receives fitted probabilities,
  and the first such row in a session warns. Both \\n\\ and \\K\\ are
  fixed at creation - every combiner buffer is sized by \\n\\, and \\K\\
  is the forest count - so only the values may change; a matrix of the
  wrong shape is refused naming the count in force. The trees carry
  over, fitted to the previous counts exactly as `setResponse` leaves
  any sampler's, and the next `run` forms every category's working
  response against the new matrix. The matrix is written to both the
  engine and `data@counts` (its row sums to `data@y`), so `getPointer`'s
  transparent re-creation after save and load carries the current
  response rather than the one the sampler was created with. Cost, not a
  defect: the sweep draws \\n_i\\ Polya-Gamma variates per observation
  per category, so replacing single-trial labels with grouped counts
  multiplies sweep cost by `mean(n_i)`. A missing count is refused, as a
  missing response is; to leave a row out of the likelihood mid-chain,
  mark it inactive with `setActiveRows`, as in every family, or write
  zeros into it, which is the same as far as the likelihood goes. A
  vector is refused there, since with each row's trials held it cannot
  say how they split among the categories. On a binomial sampler, a
  `logistic` one whose data carry the response as counts of successes
  and failures, `y` is the \\n \times 2\\ matrix of successes and
  failures, or a vector of successes with each row's trials held, a
  count above its row's trials refused; `setData` takes new data whose
  counts give both, the number of rows free. A count sampler has no
  response scale to re-anchor and ignores `updateScale`. This draws
  nothing; the next sweep draws the augmentation variables against the
  change.

- status:

  For `setResponse` on an `aft` (survival) sampler, a replacement
  per-observation censoring status: a numeric vector of length equal to
  that with which the sampler was created, 1 an uncensored event and 0 a
  right-censored observation whose response value is its truncation
  lower bound. `NULL` (the default) leaves the censoring structure in
  force and is the call every other family makes; a non-null `status` on
  any other family is refused by name. The response and the status
  change in ONE call, and the status is installed first, so a row newly
  censored is bounded at its own NEW observed log time and held there
  until the next sweep draws it above that bound against the fit then in
  force, as `setResponse` leaves every censored row at its new bound,
  while a row that becomes an event is restored to that time - data,
  rather than a draw. A status-only change is
  `setResponse(sampler$data@y, status = s)`. The number of observations
  is fixed, the status being stated over it; `setData`, whose
  replacement may change it, refuses an `aft` sampler and names this
  channel instead. The vector is validated - real, of length \\n\\,
  every element exactly 0 or 1, `NA` refused - before anything installs,
  so a refused call leaves the sampler, `data@y` and the recorded status
  exactly as they were; it is written to the sampler's data object only
  after the engine accepts it, so `getPointer`'s re-creation after save
  and load carries the current censoring rather than the one the sampler
  was created with. Under an installed mask (`active`) the response
  install is mask-blind but the sweep's draw is not, so an INACTIVE
  censored row's latent stays at its bound until the mask clears; that
  is what `setResponse` leaves such a row at with or without a status.

- x:

  For `setPredictor` with no `column`, the whole of the predictors: a
  data frame, whose columns are matched to the sampler's by name, in any
  order, a name that is not the sampler's being refused (a sampler with
  no names reads them by position), and whose factor columns are coded
  by label as a column update's are, or a numeric matrix, taken only
  where every predictor column is numeric - on a sampler with a factor
  column a matrix is refused by name, since its codes could not say
  which values are labels. A matrix that carries column names is matched
  by name as a frame is, and one with none is read by position, with no
  message. Otherwise, a numeric predictor vector of length equal to that
  with which the sampler was created. Can be of a distinct number of
  rows for `setTestPredictor`. A column update (`column` given) on a
  column coded from a factor takes that column's labels instead, for
  `setPredictor` and `setTestPredictor` alike: a factor, character
  vector or
  [`sparseFactor`](https://vdorie.github.io/dbarts/reference/sparseFactor.md),
  matched by label against the training levels, as a whole data-frame
  update is coded; several columns take a data frame. A label the column
  does not declare is refused by name, as is a number, which could only
  be read as an internal code. A missing label is a missing value:
  `setPredictor` takes it whether or not the column has held one, and
  `setTestPredictor` refuses a missing value, a label or a number, in a
  column that has never held one in training, as the whole test set is
  refused. What a column's first missing value does is under ‘Missing
  values in predictors’ in ‘Details’. On a hazard fit, `setPredictor`
  takes the period column, whole, as a column or row by row
  (`"partial"`), only unchanged, refusing any other values before
  anything moves: the column is the expansion of each subject's time at
  risk, set with the response, and `setData` replaces it. `setPredictor`
  draws nothing: the next sweep draws the augmentation variables against
  the accepted predictors.

- x.test:

  A new matrix, data frame, or sparse-bearing test set (the column types
  accepted by
  [`dbartsData`](https://vdorie.github.io/dbarts/reference/dbartsData.md)),
  of the number of columns equal to that in the current model; a numeric
  matrix is taken only where every predictor column is numeric, and is
  refused by name on a sampler with a factor column, whose test rows
  come as a data frame (the same holds for `predict` on a fit). A
  data-frame/sparse `x.test` is coded against the training levels;
  `setTestPredictor` and `setTestPredictorAndOffset` then install it as
  a resident, whole-object replacement of the test set (see `column`),
  while `predict` and `predictForests` route its rows through the trees
  off that same storage, retaining nothing.

- offset:

  A numeric vector of length equal to that with which the sampler was
  created, or `NULL`. If `offset.test` was set from `offset`, will
  attempt to update that as well.

  For `setCategoryOffset`, a multinomial sampler's per-category shift
  instead: an \\n \times K\\ numeric matrix, or `NULL` to clear one. The
  latent becomes \\f\_{ik} + o\_{ik}\\, so the shift enters the
  log-sum-exp margins, every category's working response and the
  reported probabilities, and never a leaf value. This is the
  response-side counterpart of `setResponse`'s count matrix rather than
  of `setOffset`, whose flat shift is added AFTER the categories are
  blended - the wrong side of the nonlinearity - and is the softmax's
  own null direction besides, which is why `setOffset` itself is refused
  on such a sampler. Only the row-centered part is identified: adding a
  constant to a whole row leaves every reported probability unchanged,
  and the entrance leaves the matrix as given rather than re-centering
  it. It shifts the TRAINING latent only. Mirrored into
  `data@offset.category`, so a re-created sampler carries it. Column
  names, when present, are matched to the sampler's category names (the
  column names of its counts; `"1"`, ..., `"K"` when it has none): a
  permutation is reordered by name and other names are refused, except
  that a sampler without count column names reorders a permutation of
  its indices and ignores other names. Unnamed columns stay positional.
  `setOffset` draws nothing: the next sweep draws the augmentation
  variables against the new offset, on a logistic sampler as on every
  other family.

- updateScale:

  Logical indicating whether BART's internal scale should re-anchor to
  the new offset (`setOffset`), response (`setResponse`), weights
  (`setWeights`) or active rows (`setActiveRows`). On the last two it is
  the argument after `updateState`, so the positional
  `setWeights(w, TRUE)` of earlier releases still names `updateState`.
  With `TRUE` they install the weights or the mask first and then take
  the response's scale (a count sampler's center) again from the rows
  then in the likelihood, as creation on those rows would, whether or
  not the call moved a row in or out; the response and the offset are
  not changed, and neither is a test offset, a cut point, a covariate
  standardization or the residual prior. Nothing is taken again by a
  one-forest `probit`, `ordinal`, `logistic` or `hazard` sampler, which
  accepts the argument and ignores it, also with `gp` leaves holding
  saved draws, none of which it moves. A call refused for any reason
  installs nothing; a call that keeps the scale (below) installs the
  weights or the mask. A stored state does not undo a re-derivation, and
  with `updateState` resolving to `TRUE` the state is stored after it.
  Defaults to `FALSE`, locking the scale set at creation; should only be
  `TRUE` during burn-in, as re-anchoring mid-run makes the fits across
  iterations no longer comparable. Draws already saved (`keepTrees`) are
  rewritten into the new scale, a heteroscedastic sampler's variance
  draws included, so `predict` returns the values it returned before the
  call; the current fit follows the new scale. `TRUE` is refused by a
  sampler with Gaussian-process (`gp`) leaves that holds saved draws,
  whatever the new response, offset, weights or mask, where its family
  has a scale to re-derive (`gaussian`, `student`, `aft` and `nbinom`):
  make a new sampler, or leave `updateScale` at `FALSE`. A `probit`,
  `logistic` or `ordinal` sampler has none, so it takes `TRUE` as it
  takes it in `setData`, and its saved draws predict what they did. A
  draw is saved by every sweep run as a sample while trees are kept,
  warm-up sweeps included, so such a sampler takes a re-anchor during
  burn-in only as long as its sweeps are run as burn-in, which saves
  nothing: `run(1, 0)` between re-anchors, not `run(0, 1)`. For a
  survival (`aft`) sampler the scale is always that of the observed
  times, a censored row's censoring time among them, and never of a time
  the sampler has drawn for a censored row. `TRUE` is refused by a
  multi-forest (`forests`) sampler, whose calibration map is pinned at
  creation - whatever its response family, and whether or not there is a
  transform to re-anchor - and is never restated; the response and a
  basis scale may span different rows, which is why.
  `setActiveRows(updateScale = TRUE)` is refused by a multinomial
  sampler as its other response channels are. A heteroscedastic
  (`variance`) sampler takes it: the variance forest's leaf prior is
  re-derived from the residual prior on the new scale, exactly as
  creation on the new response would have stated it, and the drawn
  \\s^2(x)\\ is carried into the same units, so the surface keeps its
  meaning in response units across the swap. `setData` carries that
  restatement too wherever the scale moves. On `setData` the argument
  defaults to `NULL`, which takes the `response.range` the new data
  object records and otherwise derives the scale; `TRUE` derives it from
  the new data, over its rows in the likelihood as creation would, and
  `FALSE` keeps the sampler's, the new response going in under it with
  no draw rewritten. A sampler made by hand from another's parts takes
  the record in the same way. Where rows are out and the rows in hold
  fewer than two distinct values of the response, a `setData` that
  derives the scale takes every row, under a warning of class
  `dbartsScaleFallbackWarning`; `setResponse`, `setOffset`, `setWeights`
  and `setActiveRows` with `updateScale = TRUE` instead keep the scale
  and the center in force, as with `updateScale = FALSE`, and warn with
  the same class once for a sampler (a copy or a reload is a new
  sampler, `setData` is not). On `setResponse`, supplying this argument
  positionally (`setResponse(y, TRUE)`) rather than by name warns once
  per session - it is the second argument, where a caller porting code
  written before this order had `updateState` there instead.

- offset.test:

  A numeric vector of length equal to that of the test matrix, or `NULL`
  for no offset; a missing value is refused. Can be missing for
  `setTestPredictorAndOffset`. Refused, by name, for `predictForests`:
  an offset shifts the combination of the forests, exactly as the
  response transform's shift does, and neither is any one forest's own
  total.

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
  asks for on purpose. Column names, when present, are matched to the
  sampler's category names (the column names of its counts; `"1"`, ...,
  `"K"` when it has none): a permutation is reordered by name and other
  names are refused, except that a sampler without count column names
  reorders a permutation of its indices and ignores other names. Unnamed
  columns stay positional.

- n.threads:

  For `predict` and `predictForests`, a per-call worker count that does
  not persist, defaulting to the sampler's own (the `n.threads` of its
  [`control`](https://vdorie.github.io/dbarts/reference/dbartsControl.md)
  object, fixed at creation). The replay is partitioned by (chain, saved
  draw), each partition writing its own rows and nothing being reduced
  across workers, so the result is identical bit for bit at every value
  and only the time taken changes. A replay small enough that the
  partition could not pay for itself runs inline whatever is passed.

  These are native threads, so `_R_CHECK_LIMIT_CORES_` does not govern
  them: that variable is read only by `parallel:::.check_ncores`, which
  counts simultaneously spawned processes and never sees a thread.

  `run` takes no thread count: a run uses the sampler's own, which
  `setControl` changes, and its draws are the same at every count.

- sigma:

  A single positive residual standard deviation, on the original
  response scale; it is assigned to every chain. `setSigma` only applies
  to a sampler whose residual standard deviation is a free parameter: it
  is refused on a probit-, logistic-, multinomial-, ordinal-, or
  count-family sampler, which fix it at 1 by the model definition, and
  on a heteroscedastic (`variance`) sampler, whose variance forest owns
  the residual scale. A gaussian sampler with
  `family = gaussian(sigma = fixed())` is not fixed in this sense and is
  still accepted: suppressing the internal draw so an outer sampler owns
  `sigma` is what that prior is for. There the write rewrites the
  model's fixed value, recorded on the `model` field, so a `copy` or a
  save and load keeps it; a write through the C API reaches the engine
  alone, and a sampler re-created after one holds the model's value
  until the next write. This draws nothing; the next sweep draws the
  augmentation variables against the change.

- weights:

  A numeric vector of case weights, of length equal to that with which
  the sampler was created. What a weight MEANS is family-specific: for a
  gaussian-family sampler it is a precision, non-negative, or under
  `weights.type = "frequency"` (fixed at creation, see
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md)) a
  count of identical observations, real and non-negative, and a weight
  of zero takes an observation out of the likelihood and out of the
  residual degrees of freedom while keeping it in the design, its place
  in a leaf and its fitted values (see below); for a `logistic`-,
  `nbinom`- or `multinomial`-family sampler it is an observation COUNT
  and must be a non-negative whole number, no larger than \\10^6\\ times
  the row's trials, \\y + r\\ or row total, a count of 0 taking the row
  out of the likelihood: it consumes no random number, and the row draws
  as one masked by `active` does. A zero-weight row is IN THE DESIGN BUT
  NOT IN THE LIKELIHOOD, and takes no part in the response scale, the
  count center or a forest's default sd, which are taken over the rows
  of positive weight when the sampler is created (a row is in when its
  weight is positive and no mask switches it off) and held when the
  weights or the mask change, unless the call carries
  `updateScale = TRUE`, so the three ways to leave a row out - at
  creation, by `setWeights` or `setActiveRows`, and by `subset` - differ
  by when the row left. It counts as a member of its leaf, so the rule
  that rejects a tree move leaving a leaf empty counts it too: the trees
  a sampler may hold do not depend on the weights, and a leaf may hold
  zero-weight rows only, in which case no likelihood term reaches it and
  its value is a draw from the leaf prior. A fit with zero weights is
  therefore close to the fit on the remaining rows alone but is not that
  fit, since the trees can spend splits on regions of zero-weight rows.
  In one measured case, 300 rows with 50 trees and the cut points held
  equal and the scale that of the remaining rows: with a random half of
  the rows at weight zero the posterior mean fit at the remaining rows
  moved by at most 0.02 of a posterior standard deviation, and with a
  whole region at weight zero by at most 0.05; the residual standard
  deviation did not move, and the fit at the zero-weight rows was 1 to 2
  percent more variable.
  [`dbartsData`](https://vdorie.github.io/dbarts/reference/dbartsData.md)'s
  `subset` drops rows outright, and with them their part in the cut
  points and the response scale. In exchange, an outer sampler that
  redraws 0/1 weights or `active` every sweep from its own conditional
  samples the joint model it assumes, which it would not if the
  admissible trees moved with the weights. A Student-t sampler holds a
  scale \\\lambda_i\\ per row, which `setWeights` leaves as it is: the
  next sweep draws every scale from its conditional given the fit and
  the weights then in force, so an outer sampler that draws the weight
  given the scales and one that draws it with the scales integrated out
  are both exact. `setWeights` sets a weighted likelihood on gaussian-,
  `logistic`- and `nbinom`-family samplers. A logistic or nbinom swap is
  a model change with a defined meaning rather than a reweighting: the
  counts multiply the shape of the Polya-Gamma latents, which the next
  sweep draws against the new counts, so an outer sampler can vary
  exposure between runs. `probit`, `ordinal` and `aft` refuse a weighted
  likelihood by identification: a weighted probit has no tractable
  latent-variable form, `ordinal` inherits that, and `aft` fixes its
  censoring structure at creation. On a `probit`, `ordinal` or `nbinom`
  sampler a vector of 0s and 1s is accepted all the same: it names the
  rows in the data set rather than a precision, is held in the data
  object's `weights` slot, all ones as no weights, and composes with
  `active` below, a row being in when its weight is 1 and the mask keeps
  it; any other value keeps the refusal on `probit` and `ordinal`, and
  is a count on `nbinom`. `setData` carries the same rule on the
  whole-data conduit, the next sweep drawing the latents against the
  counts the replacement data carries; replacement data given without
  weights is single-trial, as at creation, so a logistic sampler built
  with counts and handed weightless data becomes an unweighted one. The
  weights themselves are not part of the saved `state`, and the state
  records nothing of them: the latents it carries are the values its
  trees were drawn from, which `setState` installs without drawing, and
  the next sweep draws them against the weights then in force. A weight
  change is undone by `setWeights` with the old weights and `setState`
  with the stored state, in either order, the sampler then continuing
  bit for bit as one that installed the stored state under the old
  weights would. This draws nothing; the next sweep draws the
  augmentation variables against the change. See
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
  unaffected, a leaf of that forest may hold such rows only (its value
  is then a draw from the leaf prior), and the reported location still
  carries the forest's full contribution even when every row is
  excluded. The third leg of a three-way degenerate-value contrast:
  `weights` at all-ones installs and, on gaussian and logistic, is
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
  `setActiveRows` reaches every family, the four that refuse a weighted
  likelihood included: gaussian, Student-t, `probit`, `ordinal`,
  `logistic`, `nbinom`, and `aft` samplers all accept it, including a
  multi-forest (`forests`) sampler, whose forests share one response of
  whichever family it was built under. It is absolute and independent of
  `weights`, composing with it in either call order (\\w_i \times a_i\\
  in effect). `NULL` clears the mask (every row active). Length and `NA`
  are validated as for `weights`; the only legal values are exactly `0`
  and `1`, and a fractional element refuses the WHOLE call and installs
  nothing. An all-ones vector reports success but installs NOTHING,
  clearing any previously installed mask - the OPPOSITE of `weights` on
  gaussian and logistic, where an all-ones vector installs and is
  measurably distinct from carrying no weights at all. One channel
  states membership, the other precision, and the two deliberately carry
  opposite degenerate-value policies. The response scale, the count
  center and a forest's default sd are held across a mask as across
  weights, unless the call carries `updateScale = TRUE`, and a `linear`
  or `gp` leaf's covariate standardization unless it carries
  `updateStandardization = TRUE`. On a `probit`, `ordinal` or `nbinom`
  sampler, whose 0/1 weights state membership too, `activeRows` holds
  only what `setActiveRows` installed, and a row either channel brings
  back in draws nothing at the call, the next sweep drawing its latent
  first.

  An inactive row (`active[i] == 0`) contributes nothing to any leaf
  sufficient statistic, branch log-likelihood, birth-scan weight total,
  leaf parameter draw, or family-level parameter update that sums over
  rows (residual degrees of freedom, a shape statistic, a group's
  per-group sums, and so on). It still occupies its leaf for COUNT-based
  accounting - `numObservations`, the birth scan's own member count, and
  leaf collapsing, which triggers on member count regardless of
  weights - and it still receives a fitted value: `run()$train`,
  `getForestFits`, `getFitsWithoutOffset`, and `predict` report
  \\f(x_i)\\ at an inactive row exactly as at an active one, which is
  what makes this channel worth more than physically dropping the row.
  The rule that rejects a tree move leaving a leaf empty is count-based
  too: an inactive row is in the design but not in the likelihood, so a
  leaf all of whose rows are inactive is legal, takes no likelihood term
  and has its value drawn from the leaf prior, and the trees a sampler
  may hold do not depend on the mask. That is what lets a larger sampler
  that redraws the mask every sweep, from each row's conditional given
  the fit, sample the joint model it assumes. Its cost is that a fit
  under a fixed mask is close to the fit on the active rows alone but is
  not that fit, by the amounts given under `weights`. With every row
  inactive every leaf is such a leaf, which is the exact sense in which
  the forest “sits at its prior” - a distribution over structures, not a
  frozen one.

  An inactive row's own latent draw is SKIPPED - no random numbers are
  consumed for it - for every family except Student-t, whose per-row
  \\\lambda_i\\ is drawn unconditionally at every row regardless of the
  mask; there the mask instead ANNIHILATES its contribution through the
  composed weight (\\w_i \lambda_i a_i\\). For a skipping family,
  `getLatents` at an inactive row returns its LAST DRAWN value, which is
  stale; the correct read at an inactive row is its fit, not its latent.
  A stale latent does not come back into the likelihood: each sweep
  draws the latent of every row the mask holds in, for Student-t its
  \\\lambda_i\\ of every row, from the row's conditional given the fit
  and before its trees, so a row switched back in is drawn afresh
  against the current fit with the rest. `setActiveRows` itself draws
  nothing and consumes no random numbers, whatever it switches, so a
  caller that draws the mask from the fit with the latents integrated
  out (a probit row's membership from its Bernoulli likelihood, say)
  need not call `setResponse` to refresh them. A `multinomial` sampler
  holds no latent between sweeps. A state installed over a different
  mask keeps, at rows that were out when it was stored, the latents they
  left with until the next sweep draws the rows the mask now holds in.

  An inactive row's pointwise log-likelihood - the value the flat C
  API's `logLikelihood` results field reports (see ‘Mutation cost’ above
  for the C interface) - is `NaN` rather than a finite value or `-Inf`:
  a row that is not in the model has no likelihood to report.

  Three quantities do NOT follow the mask, by design, so a caller
  wanting them recomputed on the active subset must re-create the
  sampler instead: the split-candidate cut grid stays the FULL-data
  grid, holding the tree prior fixed while the row set moves; the
  response transform and residual prior scale stay creation's
  calibration (a masked gaussian - and so `aft`, and any `variance`
  decoration - keeps the range and residual-prior scale creation set,
  which were taken over the rows in the likelihood then, every row when
  none was out); and an ordinal sampler keeps its FULL-data number of
  categories, free thresholds, and log-gap prior even when a mask
  empties a boundary category - only `setData` changes them. `setData`
  also CLEARS an installed mask - the mask is a per-row vector and
  `setData` may change the number of rows - so a caller replacing the
  data must reinstall any mask it wants kept.

  The mask does not ride a sampler's saved `state` - the state carries
  derived quantities and no raw conditioning vector - so, exactly as
  `setForestWeights`'s per-forest weight is, it is mirrored on an R5
  field that `getPointer`, `setState` and `copy` all reinstall on every
  re-creation: a masked sampler saved and reloaded goes on masked, and a
  caller never reinstalls it by hand. What the mirror cannot do is
  travel with the state itself. Installing a stored state into a
  DIFFERENT sampler object leaves that object's own mask in force,
  whatever the donor's was - a destination that never had one stays
  unmasked and computes different draws while `statesAgree` still
  reports agreement - since the state holds no mask for the restore to
  re-derive against. Conversely, the latents DO ride the state, and
  under an installed mask they are stale at every inactive row for a
  skipping family, so two samplers at the same posterior state but
  different mask histories carry different stored latents and
  `statesAgree` correctly reports DISAGREEMENT; a Student-t sampler is
  affected the same way, since an inactive row's stored \\\lambda_i\\ is
  drawn from a conditional its own data no longer informs.

- basis:

  The data the named forest's amplitudes multiply, with as many rows as
  the sampler was created with and one amplitude per column. It is a
  value, or a one-sided formula, which is read as the `basis` of a
  [`forest`](https://vdorie.github.io/dbarts/reference/forest.md) is,
  its names found where the formula was written. It must be of the class
  the forest was created with, which decides what the forest multiplies
  and never changes. A forest of levels, created with a factor, a
  character vector or a logical vector, takes one of those three,
  matched to the forest's levels by name and not by position: the
  replacement may list its levels in another order or leave one out, a
  level of the forest that no row has is a column of zeros, and a row at
  a level the forest does not have, or at `NA`, is refused. Coefficient
  \\j\\ multiplies level \\j\\'s rows after any swap, whatever the
  replacement declares. Numbers are refused for such a forest, the level
  indicators the sampler stores among them: to hand a basis back, hand
  back the factor. A forest of numbers takes one numeric column, a
  vector or a matrix of one column, a column of zeros among them, and
  refuses several columns, a factor, a character vector and a logical
  vector; [`as.numeric()`](https://rdrr.io/r/base/numeric.html) of a
  logical vector is its 0/1 column. A forest created with no basis takes
  none. The sampler keeps the forest's block, one column for each level
  of a factor with no reference level dropped, in `data@bases[[f]]`, and
  its levels in `data@basis.levels[[f]]`, `f` the forest's position
  (neither list has names, so a label does not index them), and not the
  factor it was handed: a caller that will hand a basis back, to restore
  a proposal it rejected, keeps the factor it handed over. No forest
  changes width, held or drawn, so the names its columns were created
  with stay, whatever names the replacement has. The forest's label does
  not change. Unlike creation, a swap may leave a level of the forest
  with no row, and a numeric basis may hold a column of all zeros, the
  amplitude moving under its prior alone until a row takes it again, so
  a swap inside a larger sampler can leave a level or an arm momentarily
  unobserved; creation drops a level no row has and refuses a column of
  zeros. `setForestBasis` is the *sole* route by which a basis changes
  after creation, and it applies only to a sampler whose forests carry
  amplitudes, built with `forests = ` (see
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md)). The
  amplitudes are preserved: an install leaves every one of them bitwise.
  A swap moves the basis and nothing else: no forest's prior moves, so a
  numeric forest keeps the `sd` per unit its creation gave it, a default
  one the scale of the column it was created on (`data@basis.scale`),
  unless `updateBasisScale = TRUE` takes the default of the new column.
  A basis that is not one of the canonical shapes (a dense all-ones
  column, or a complementary 0/1 pair) moves that forest onto the
  general amplitude conditional, which is a model fact and not
  revertible by a later data swap. The matrix is written to both the
  engine and `data@bases`, so `getPointer`'s transparent re-creation
  after save/load carries the current assignment rather than the one the
  sampler was created with; because the bases ride creation, a swap
  applied after a state restore preserves the *restored* amplitudes. A
  state records the levels of each forest a factor multiplies, and
  `setState` puts each of its coefficients on this sampler's level of
  the same name, refusing a state whose names differ (see ‘Saving’).
  This draws nothing; the next sweep draws the augmentation variables
  against the change.

- updateBasisScale:

  For `setForestBasis`, `TRUE` or `FALSE`, the default. `TRUE`, meant
  for burn-in, also takes the default of a numeric forest that states no
  `sd` from the column installed: the column's standard deviation is
  read over the observations in the likelihood at the call (a positive
  weight, kept by the active-row mask; a forest's own weights from
  `setForestWeights` are not read), and the forest's `sd` becomes one
  unit of the response's scale for one standard deviation of the new
  column, from the next sweep, as for a sampler created on that column
  over those observations (see `sd` in
  [`forest`](https://vdorie.github.io/dbarts/reference/forest.md)). The
  trees, the coefficients and every other forest are left as they are,
  and no draw is taken. The rows read are those in at the call, which
  may not be the rows the response's scale was read from at creation,
  and re-installing the column in force with `TRUE` is how to take its
  scale again after the weights or the mask change. The new scale
  replaces the forest's entry of `data@basis.scale`, so a copy, a saved
  and reloaded sampler and a sampler made again from its own control,
  model and data run the new prior. Where fewer than two observations
  are in the likelihood, the column takes one value over them, or the
  weights are counts that total at most 1 over them, the column is
  installed and the scale in force kept, with a warning the first time
  the sampler keeps a scale; a sampler created on such a column with no
  `sd` is refused. `TRUE` is refused, before anything is installed, on a
  forest that states its `sd`, at creation or through `setLeafPrior`,
  and on a forest multiplied by a factor, neither having a default taken
  from a column.

- forest:

  A single positive integer indexing the forests from `1`, or a string,
  the forest's label (see
  [`forest`](https://vdorie.github.io/dbarts/reference/forest.md);
  `forest<i>` names position *i*, a sampler of one forest has no labels,
  and a multinomial one's are the response's levels), or - for
  `getForestFits`, `getForestAmplitudes`, `getForestVariableCounts`,
  `getLeafPrior`, and `getK` - `NULL` (the default on all five) for
  every forest, stacked (see ‘Value’ for each reader's stacked shape). A
  single-forest sampler's `NULL` read is bitwise its `forest = 1` read,
  so the default costs nothing on an ordinary sampler.
  `setForestWeights` and `setForestBasis` have no default: a writer
  names one target rather than reading every one. `setLeafPrior` takes
  no `forest`: a multinomial sampler's `normal(k = )` states every
  category forest, and a sampler whose forests carry amplitudes states
  one forest's size through `setForest`. `setForest` and `getForest`
  read their `forest` as the others do, left out on a sampler of one
  forest and on a multinomial one, whose `setForest` refuses a position
  or a level: its forests take one statement. Declared in that order, a
  Bayesian causal forest's prognostic forest is `1` and its basis forest
  `2`; `setForestWeights` and `setForestBasis` are both refused on a
  sampler whose forests carry no amplitudes, but not with the same
  message - `setForestWeights` names the missing capability, while
  `setForestBasis` raises `"forest index out of range"` - and
  `setForestBasis` accepts any forest of one that does. `getForestFits`
  and `getForestVariableCounts` accept `forest = 1` on any sampler - it
  selects the only forest - and refuse only an out-of-range index.
  `getLeafPrior` and `getK` are likewise served on every forest of a
  multi-forest sampler, and the calibration-map entries are that
  forest's own. `getTrees` takes the same `forest = NULL` default,
  stacking every forest forest-major with a leading `forest` column on a
  sampler with several forests or one declared with `forests =`, and no
  such column on any other (see ‘Value’); it additionally accepts a
  vector of forest indices or of labels, each validated as a lone
  `forest` argument is here and stacked in the order given. `plotTree`
  takes a single `forest`, defaulting to the sampler's only one and
  required when there is more than one.

- leaf.prior:

  For `setLeafPrior`, a leaf prior specification in the vocabulary a
  fitting function's `leaf.prior` takes, resolved inside the call as it
  is there: `normal(sd = )`, `normal(k = )`,
  `normal(sd = invchi(df, scale))`, or `linear(sd = )` and `gp(sd = )`
  on those leaf models (see
  [`dbartsPriors`](https://vdorie.github.io/dbarts/reference/dbartsPriors.md)).
  It must name the sampler's own leaf model - `normal` for the constant
  and monotone leaves - and may omit leaf-model details such as a linear
  leaf's `columns`; any it states must match. Changing the leaf model,
  or anything else about the model, is `setModel`'s. A multinomial
  sampler takes `normal(k = )` with a fixed `k`, `Inf` included, as its
  creation does. On any sampler `normal(k = Inf)` holds every tree at
  zero, so the fit is the offset, plus on a continuous response the
  center of the range the response less its offset spans. On a sampler
  whose forests carry amplitudes it is the fitting function's leaf
  prior; see ‘Value’.

- sd, power, base, amplitude, proposal.probs, n.perturb.cuts, treeShift:

  For `setForest`, and of them `power`, `base`, `proposal.probs` and
  `n.perturb.cuts` for `setVarianceForest`, a value in the terms
  [`forest`](https://vdorie.github.io/dbarts/reference/forest.md) takes
  it, a coefficient's prior resolved inside the call as it is there. For
  `setTreePrior`, `power` and `base` in the terms a fitting function
  ([`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md))
  takes them. `NULL`, the default, is the same as leaving the name out.
  See ‘Value’.

- split.probs:

  For `setTreePrior`, the split probabilities in the terms a fitting
  function takes them: a vector, `1` for equal probabilities, or on a
  DART sampler `dart(...)`. `NULL`, the default, leaves them.

- cuts:

  For `setCutPoints`, a numeric vector of cut points for the column
  named, or a list of such vectors, one for each column named; a data
  frame is taken as the list of its columns. A vector holds no `NA` or
  `NaN` and is sorted if it is out of order. A cut point must be finite;
  `-.Machine$double.xmax` separates values at `-Inf`. A point may appear
  only once, so a grid with a repeated point is refused; for more splits
  near a value, give a denser grid around it. No grid a sampler holds
  repeats a point, so the grids it reports can be handed back to it
  whole. A grid given here is unweighted, its cut points equally likely,
  except the grid the column holds handed back bit for bit, which keeps
  the weights the default rule gives a column with fewer distinct values
  than `n.cuts`. With `column` missing, `cuts` is a list with one entry
  for every column of the predictors, and the entries of factor columns
  are not read, whatever they are. Cut points apply to numeric
  predictors only: naming a factor column of either kind is refused,
  naming the kind, its grid following its level table rather than any
  externally chosen one. A later `setState` does not undo the change: a
  state goes onto the cut points the sampler holds (see ‘Saving’). The
  sampler records a grid set this way in its data object, so a later
  `setData(sampler$data)` keeps it; a later `setData` that derives the
  cut points, as one given a data object made by
  [`dbartsData`](https://vdorie.github.io/dbarts/reference/dbartsData.md)
  does, derives every numeric column's cut points again whatever grid
  was set, at most as many as the `n.cuts` the sampler was created with:
  with `useQuantiles`, a column with no more distinct values than
  `n.cuts` gets one cut point fewer than it has distinct values, and a
  constant column gets one.

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
  few columns move. A replaced sparse-backed column stays sparse: it
  stores the replacement's entries that differ from its implicit value
  (zero, or a sparse factor's reference level), in the storage layout
  chosen at creation, and a sparse column of the replacement is read as
  those entries rather than expanded, so for sparse columns the call's
  memory follows their nonzeros; a dense column of a mixed replacement
  is copied, as in any dense replacement. A replacement of any other
  sparse `Matrix` class is read as the `dgCMatrix` it holds, as is such
  a test set. For `setTestPredictor`, a per-column update reaches a
  data-frame or sparse test container too: the replacement is spliced
  into its storage rather than densifying the whole test set.

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
  logical of what was installed. Under a
  [`monotone`](https://vdorie.github.io/dbarts/reference/monotone.md)
  constraint a new predictor can also leave a tree's leaf values out of
  order along a constrained predictor (that page says what this is).
  `TRUE` completes such an update, setting every leaf value of that tree
  to zero without a message; the next iteration draws them again.
  `FALSE` refuses it and rolls it back, as for an empty leaf, and
  `"partial"` declines each row that would do it. The one unforced
  change that can do this is an unordered factor column's first missing
  value, under the sides drawn for the column's rules when it arrives
  (see
  [`monotone`](https://vdorie.github.io/dbarts/reference/monotone.md)).
  `"partial"` requires a single `column` and cannot be combined with
  `updateCutPoints`; it also requires a dense-backed column, since a
  sparse column's storage fixes its nonzero pattern per cell - replace
  such a column whole instead. When missing, defaults to `TRUE` when the
  whole predictor matrix is being replaced (`column` is missing) and
  `FALSE` when a single column is being replaced.

  For `setState`, a single `TRUE` or `FALSE`. `FALSE`, the default,
  installs a state only where no tree of it has to be changed to fit the
  sampler, and otherwise leaves the sampler as it was; `TRUE` installs
  it, repairing the trees that do not fit (see ‘Saving’). `"partial"` is
  not taken.

  For `installTrees`, `FALSE` (the default) takes the donor only if
  every tree it gives can stand in this sampler as the donor holds it,
  and `TRUE` takes it changing what must be changed (see ‘Warm starts’).

- updateCutPoints:

  For `setPredictor`, one of `"none"`, `"position"` or `"value"`, or a
  unique abbreviation of one. `"none"`, the default, keeps the cut
  points (split candidate locations) of the replaced column(s) and the
  weights they carry. The other two derive them again from the new
  values, as a new sampler would for those values: up to the `n.cuts`
  the sampler holds, whatever number the column held, and fewer where
  the values supply fewer, a constant column getting one, and under the
  default rule one per gap, weighted by width, on a column with fewer
  distinct values than `n.cuts`. They differ in where the splits already
  on such a column go. Under `"position"` each keeps its position on the
  grid, which is what `TRUE` did in 0.9-x; when the number of cut points
  changes, the position is rescaled to the new number. Under `"value"`
  each moves to the new cut point nearest its old threshold. In either
  case a split stays between the splits above it in its tree; one left
  with no cut point there is merged when the update is forced, and
  otherwise the update is refused and rolled back, trees and cut points
  included. A later `setState` does not bring back the cut points a
  state was stored on (see ‘Saving’). Anything but `"none"` cannot be
  combined with `forceUpdate = "partial"`. For `setData` the same three
  words, or `NULL`, the default: `"none"` keeps the sampler's cut points
  and their weights, the new values read on them, `"position"` and
  `"value"` derive them from the new data as above, and `NULL` takes the
  `cutPoints` the new data object records, a split on a column whose cut
  points that changes moving to the nearest of them, and otherwise
  derives them as `"value"` does, as `setData` always did. For
  `setPredictor`, a logical is taken until dbarts 1.1-0, `TRUE` as
  `"position"` and `FALSE` as `"none"`, with a warning once per session
  ([`dbarts-deprecated`](https://vdorie.github.io/dbarts/reference/dbarts-deprecated.md)).

- updateStandardization:

  For `setPredictor`, `FALSE`, the default, keeps a `linear` or `gp`
  leaf's covariate standardization, each covariate's center and scale
  and a `gp` leaf's heuristic lengthscale, whatever values the update
  installs; `TRUE` derives it again for the covariates among the columns
  the update replaces, every column where `column` is missing, over
  their new values at the rows in the likelihood under the weights and
  mask in force, as creation would over them, once the update has gone
  in, keeps the other covariates', and rewrites a `linear` leaf's
  coefficients, current and saved, into it, so that each leaf is the
  function of the covariate it was and `predict` returns what it did. A
  replaced covariate those rows hold at one value, or observe at none,
  keeps its numbers, and the call warns, with class
  `dbartsScaleFallbackWarning`, on the sampler's first keep, of this or
  of the response's scale (`updateScale`); on a hazard fit where fewer
  than two subjects are in, every replaced covariate keeps, the period
  among them, the warning counting subjects. For `setData`, `NULL`, the
  default, takes the `standardization` the new data object records and
  otherwise derives it from the new data's rows in, `TRUE` derives it,
  both reading every row where fewer than two rows are in or the
  response's scale is read from every row, and refused where two or more
  rows in hold one value of a covariate unless a `gp` leaf states its
  `lengthscale`, which takes the placeholder (see `newData`), and
  `FALSE` keeps the sampler's, the new rows read under it. A leaf that
  takes no covariate has nothing to derive. A `gp` sampler that holds
  saved draws refuses `TRUE` on `setPredictor`, and on `setData` any
  call that would change the standardization (see `newData`). It cannot
  be combined with `forceUpdate = "partial"`.

  For `setWeights` and `setActiveRows`, the argument after
  `updateScale`: `FALSE`, the default, keeps the standardization
  whatever rows the call leaves in; `TRUE` installs the weights or the
  mask first and then derives every leaf covariate's center and scale,
  and a `gp` leaf's heuristic lengthscale, again over the rows then in
  the likelihood, as `setPredictor` with `TRUE` would on every leaf
  covariate's values in force, rewriting a `linear` leaf's coefficients
  in the same way. A covariate those rows hold at one value, or observe
  at none, keeps its numbers, the others re-deriving, and the call warns
  as `setPredictor` does, once per sampler. On a hazard fit, where fewer
  than two subjects are in after the call, every covariate keeps its
  numbers, the period among them, which varies within a subject, under
  one warning that counts subjects. Where the call also passes
  `updateScale = TRUE` and keeps the response's scale, those rows
  holding fewer than two distinct values of the response, the prior in
  force is kept whole: every covariate keeps its numbers too, whatever
  its spread over those rows, and one warning names the response and the
  covariates, as creation on those rows would read the response's scale
  and the covariates from every row together. A `gp` sampler holding
  saved draws refuses it before anything is installed. A sampler whose
  leaves take no covariate, one of several forests among them, takes it
  and has nothing to derive. A `probit`, `ordinal` or `nbinom` sampler's
  `setWeights` passes it to `setActiveRows` with the 0/1 weights.

  A loop whose rows in the likelihood change, such as a mixture
  component or an imputation step, can let the response's scale and the
  leaf covariates' center and scale follow the rows in during burn-in
  and then hold them, so every kept draw comes from one prior: pass
  `updateScale` and `updateStandardization` as `TRUE` on the call that
  changes the rows while burning in, and `FALSE` after. A `gp` leaf that
  has kept a draw refuses both, so burn in with sweeps that keep
  nothing:

      sampler <- dbarts(y ~ x1 + x2, data, weights = w,
                        leaf.prior = linear("x1"),
                        control = dbartsControl(n.chains = 1L, keepTrees = TRUE))
      for (i in seq_len(n.burn + n.keep)) {
        w <- drawMembership(sampler)            # the outer sampler's step
        burning <- i <= n.burn
        sampler$setWeights(w, updateScale = burning,
                           updateStandardization = burning)
        draw <- if (burning) sampler$run(1L, 0L) else sampler$run(0L, 1L)
      }

  Creation and a `setData` that derives read every row where fewer than
  two rows are in, or where the response's scale is read from every row,
  and refuse two or more rows in that hold one value of a leaf
  covariate, a `gp` leaf that states its `lengthscale` taking the
  placeholder (see `newData`); these setters keep instead, one row or
  one subject in included.

- splits:

  For `setCutPoints`, `"position"` (the default) or `"value"`, or a
  unique abbreviation: where the splits already on a column go when its
  cut points are replaced, as for `updateCutPoints`. The change is
  always forced, so a split left with no cut point between the splits
  above it is merged.

- treeNums:

  An integer vector listing the indices of the trees to print or return.

- chainNums:

  An integer vector listing the chains to return from `getTrees`. When
  missing, all chains are returned.

- sampleNums:

  An integer vector listing the saved samples to return from `getTrees`
  or print from `printTrees`. Applies only when `keepTrees` is `TRUE`
  and `current` is `FALSE`; otherwise the live working trees have no
  sample dimension, and supplying it explicitly is ignored, with a
  warning. Sample numbers address recorded DRAWS on the oldest-first
  axis `predict` reports - `1` to the number of draws recorded, at most
  `control@n.samples` - and not the store's internal slots.

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
  Refused by name by `getSigmas` and `getSumsOfSquaredResiduals`: both
  readers allocate their own vector, and filling a caller's buffer in
  place is `getLatents`'s contract alone.

- test:

  For `getVariance`, a logical; `FALSE`, the default, reads the variance
  surface at the training rows and `TRUE` at the installed test rows.

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
  donor must share this sampler's predictors; a donor of more trees
  gives some of them, and a donor on other cut points has its splits
  moved onto this sampler's under `forceUpdate = TRUE`, and needs no
  move where every split it holds is one of this sampler's cut points,
  as on a finer grid that holds the donor's; and a different center and
  scale of a `linear` leaf's covariates is allowed (see ‘Warm starts’).

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

- ptr:

  For `storeState`, the sampler's engine pointer to capture from;
  defaults to `getPointer()` and is rarely passed explicitly.

- newState:

  For `setState`, a state object previously produced by this sampler or
  another over the same data: its `state` field, calling `storeState()`
  first to refresh it if a later mutation may have gone unstored
  (`storeState` itself returns `NULL` invisibly, so its return value
  cannot be passed here). A sampler's own `state` field is never refused
  for the cut points it records: where a reload is made on other cut
  points than its stored state records, as after a `setCutPoints` with
  no `storeState()` since (see ‘Saving’), the field holds the state of
  its engine in place of the stored one. A reload's field is brought to
  its engine's at first use, and a state read from it before that is the
  stored one. Must inherit from `bartcoreState`. `setState` installs the
  latents the state carries, the values its trees were drawn from, and
  draws none.

- n.trees, vars, interactions, blocks:

  For `setForest` (and `n.trees` and `vars` for `setVarianceForest`),
  the names
  [`forest`](https://vdorie.github.io/dbarts/reference/forest.md) takes
  that say what the model is. They are refused by name whatever their
  value, being part of the model, and are arguments so that a short name
  matches as it does at creation.

- ...:

  For `setForest` and `setVarianceForest`, nothing is taken: any name
  not listed is refused as an unused argument. For `plotTree`, extra
  arguments to [`plot`](https://rdrr.io/r/graphics/plot.default.html).
  For `run`, only 0.9-x's per-call thread count, as `n.threads`,
  `numThreads` or the one unnamed argument after `updateState`, which is
  ignored with a once-per-session warning until dbarts 1.1-0
  ([`dbarts-deprecated`](https://vdorie.github.io/dbarts/reference/dbarts-deprecated.md));
  anything else is refused. A `callback` to `run` must therefore be
  named.

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

  A `dbartsModel` holding the parsed tree, leaf, and residual priors,
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
  materialized. It holds the chains and not the model: trees, leaf
  values, the quantities the sampler draws, the cut points the trees'
  split thresholds were stored on (the attribute `cutPoints`, a list
  with one entry per predictor column, `NULL` for a factor column that
  splits by subsets of its levels; with their weights, the attribute
  `cutMass`, where a column's cut points are not equally likely) and the
  leaf-covariate standardization they are read through, a record of the
  response mapping they were read under, which predictor columns could
  hold a missing value when it was stored, and the levels of each forest
  a factor multiplies (the attribute `basisLevels`) - no leaf prior and
  no value the sampler holds fixed (see ‘Saving’). It also carries each
  chain's generator state, so restoring it continues the same streams
  the sampler was drawing from, typically to the last few digits, not
  bitwise. The saved-tree store's write position and the number of draws
  it has recorded both ride it, so `predict` after a `setState` reports
  the same draws, in the same order, as before the store. Reading it
  forces the sampler's *current* state only the first time, before any
  value has been materialized; once set, it is a cached snapshot that a
  later mutation does not refresh automatically - call `storeState`
  again, or pass `updateState = TRUE` to the mutating call (see
  `updateState` above), to bring it forward. After a `setState` that put
  a state onto cut points, or weights of them, other than the ones the
  state records, the field holds the sampler's state as `storeState()`
  gives it and not the state as it was given, so that a `copy` or a
  reload has the sampler's cut points and trees. It is the only field
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

Besides [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md),
a sampler is created by
`new("dbartsSampler", control, model, data, active = NULL)` from the
list
[`dbartsSpec`](https://vdorie.github.io/dbarts/reference/dbartsSpec.md)
returns. The list's data object carries the fit's weights, a `"probit"`,
`"ordinal"` or `"nbinom"` fit's 0/1 weights among them, so a sampler
created from it is the sampler `dbarts` builds
([`dbartsSpec`](https://vdorie.github.io/dbarts/reference/dbartsSpec.md),
`data`, gives the by-hand route). `active` is a 0/1 mask installed once
the sampler is built, as `setActiveRows(active)` installs it, composing
with the weights; a sampler created with it reads over the rows in,
those the weights and the mask both keep, the center of an `"nbinom"`
count, a `linear` or `gp` leaf's covariate center and scale, and under
`"probit"` and `"logistic"` the default `sd` of a forest on a numeric
column that states none.

A sampler made this way takes its tree count from the model and its
quantile rule and cut count from the data object, where they state them.
A control that is fresh, or copied from a fit, says nothing of them, and
a copied count or rule that differs draws a warning saying it is not
used. A count or rule written on the control is used where the model or
the data object states none and is refused, by name, where it states
another, or, for a count, where every forest has a `basis`. Whether the
family is binary is the model's family alone.

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
`$setResponse`, which takes the \\n \times K\\ count matrix,
`$setCategoryOffset` and `$setCategoryTestOffset`; the predictor family
(`$setPredictor` in all four shapes, `$setCutPoints`,
`$setTestPredictor`) and the global `$setActiveRows` are open, and
`$predict` answers with the \\K\\ category probabilities and takes the
matching per-category offset. Everything else on the mutation surface is
refused, each by a message naming the capability and, where one exists,
the channel that serves the caller instead: `$setResponse` given a
vector and `$setOffset` (a flat vector cannot express the matrix, and a
flat shift is the softmax's null direction), `$setActiveRows` at
`updateScale = TRUE` (the categories keep their scale at creation),
`$setWeights` (an integer case weight is already row-wise replication in
the count response, and a non-integer one has no exact augmentation
sampler), `$setSigma` (no residual scale), `$setData` and `$setModel`
(the K category forests fix their data and their calibration at
creation), a named sd or a `k` hyperprior on `$setLeafPrior`, which
takes `normal(k = )` with a fixed `k` (the softmax map owns every
category forest's leaf scale), `$setForestWeights` and `$setForestBasis`
(its forests are its categories, which carry no amplitudes), and
`$getFitsWithoutOffset` (its reported channels are probabilities, not
one additive location). Per-forest masking is refused permanently rather
than pending: a category's margin is a log-sum-exp over the other
\\K-1\\, so a mask on one forest alone has no conditional to be a mask
of.

The composition rule follows from the table: an outer block conditions
on \\f(x_i)\\, so it reads `$getFitsWithoutOffset()` and adds back
whatever offset it installed, rather than differencing `$getLatents()`
against `run()$train`: a probit, ordinal or aft latent is on the
offset-inclusive scale, \\z = f(x) + o + e\\, as `run()$train` is, so
the difference is the residual \\e\\ with none of the offset's signal. A
block that owns the offset regresses `z - $getFitsWithoutOffset()`; a
precision latent carries no location and no offset, and
[`dbartsWorkingResponse`](https://vdorie.github.io/dbarts/reference/dbartsAugmentation.md)
gives the offset-free working response.

Where the boundary runs. These are the methods of the mutable sampler,
and their audience is a model writer composing BART into a larger
scheme: they return engine primitives on the engine's own terms. The R
modeling conventions -
[`fitted`](https://rdrr.io/r/stats/fitted.values.html),
[`predict`](https://rdrr.io/r/stats/predict.html),
[`extract`](https://vdorie.github.io/dbarts/reference/bartBT.md) - apply
to FIT objects instead, the results of
[`bartBT`](https://vdorie.github.io/dbarts/reference/bartBT.md) and
[`bart`](https://vdorie.github.io/dbarts/reference/bart.md). Those
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
model, or restore the state with the `dbarts` release that wrote it. A
sampler saved by `dbarts` 0.9-x refuses every method, naming that
version and saying to refit it with this one; its `data`, `control` and
`model` fields still read.

To restore a saved state into a sampler, call `setState(newState)` and
check what it returns, as in
`if (!sampler$setState(state)) stop("the state no longer fits")`.
`setState` validates that `newState` inherits from `bartcoreState` and
re-creates the underlying engine if needed: a sampler whose engine is
gone, as after a [`load`](https://rdrr.io/r/base/load.html), is first
rebuilt on the cut points, weights, response range and leaf
standardization its data object records and from its own saved state, as
any other call rebuilds it, so that `setState` gives what it gives a
sampler that was never saved; one with no saved state, or one whose
saved state is refused, is created anew on the same record, and
`newState` goes onto its cut points as onto a live sampler's. Then,
where every tree of the state fits the sampler as it stands, it installs
the state, caches it on the `state` field and returns `TRUE`; where one
does not, it installs nothing, leaves the sampler, its engine and its
`state` field exactly as they were, and returns `FALSE`. `TRUE` does not
mean bit for bit. It means the trees and leaf values are the stored
ones, and the quantities the sampler draws are the stored ones where the
state holds them. Of the kept draws (`keepTrees`) the newest that fit
the sampler's own store are taken, whatever number the state holds. The
latent variables it carries go in as stored and none is drawn, the next
sweep drawing them against the weights and censoring status then in
force, and the chain continues the state's own generator streams from
where they left off, typically to the last few digits, not bitwise (see
the Reproducibility section of
[bart](https://vdorie.github.io/dbarts/reference/bart.md)).

A tree does not fit, and `setState` returns `FALSE`, when it would have
to be changed before the sampler could hold it. The trees are judged
against the predictors and the cut points the sampler holds at the call
and against its own constraints: a split whose threshold is not one of
the sampler's cut points for its column; a leaf that no row of those
predictors reaches; a split the state holds outside the range the splits
above it leave; a rule that records a side for missing values on a
column that has never held one in this sampler (see ‘Missing values in
predictors’); a split on a column the tree's forest may not use - a
forest's `vars`, a `blocks` group, a moderator forest's columns or a
`variance = ~ x1 + x2` variance forest's; a split that breaks the
forest's `interactions` limit; and, under a
[`monotone`](https://vdorie.github.io/dbarts/reference/monotone.md)
constraint, leaf values that are out of order as stored.
`setState(newState, forceUpdate = TRUE)` installs such a state all the
same, as a forced `setPredictor` leaves the trees it meets, and returns
`NULL` invisibly: a split whose threshold is not a cut point moves to
the nearest cut point inside the range the splits above it leave, to the
upper one of two equally near, and one with no cut point left there
becomes a single leaf, along with everything beneath it, at the mean of
the leaves it replaces, weighted by the rows each held under the stored
thresholds; a split with a side no row reaches becomes a single leaf,
along with everything beneath it, at the mean of the leaves it replaces,
weighted by their rows; so does a split outside its range; a tree is cut
back in the same way at the first split from its root that the forest's
columns or its interaction limit do not allow, the splits above it kept;
a rule stops recording a side its column cannot use; and a monotone tree
whose values are out of order has every leaf value set to zero, to be
drawn again at the next iteration. After a forced `setState` that
repaired a tree of a state stored on the sampler's own cut points and
weights the `state` field holds the state as it was given, not as
repaired: installing the field again without `forceUpdate` returns
`FALSE`, and a `copy` or a reload repairs it the same way. A state that
went onto other cut points or weights than it records, forced or not, is
not kept: the field then holds what the sampler holds after the install.
Likewise after an install that drew a side for missing values (see
‘Missing values in predictors’) the field holds the state as given,
without the sides, and installing it again, as a `copy` or a reload
does, draws the same ones. `storeState()` stores what the sampler now
holds. Forcing checks nothing first, so a loop that restores its own
state onto rows it has already put back, where nothing can need repair,
can pass `forceUpdate = TRUE` and skip the check, which costs about as
much as one iteration of the sampler. To undo a predictor change
exactly, put the old predictor back first, with
`setPredictor(..., forceUpdate = TRUE)`, and then call `setState`. An
unforced `setPredictor` can be refused, returning `FALSE` and leaving
the changed predictor in place; `setState` then judges the state against
the changed predictor, as it does when called before the predictor is
put back.

Some states are refused, with an error, whether forced or not, and the
sampler is left exactly as it was: one that is not a `bartcoreState` or
is of an older encoding, one with another number of chains, forests or
trees or another leaf model, one whose recorded levels are not the
sampler's (below), and one with a malformed block or tree. A state in
which a numeric column's cut points repeat a value is refused, naming
the column and the value: no sampler of this version writes one, and
with every point held once a stored split's value names one position. A
Student-t, logistic or negative-binomial sampler also refuses a state
whose latent variables are not all positive and finite. So is one
holding a value no chain could run from, which no sampler writes: a
`sigma` that is not positive and finite; a `k` that is negative, or is 0
or infinite under a prior whose scale could not draw that value (a draw
reaches 0 only under a scale below about 3e-145 and infinity only under
an infinite scale or one above about 1e146); a response transform that
is not finite; a random-number generator state out of range or all zero;
DART split weights that are negative or not finite, or 0 on a column the
forest may split on, or a DART concentration that is not positive and
finite; and ordinal cutpoints whose first is not 0 or that do not
increase. Beyond that a state is not checked against the response family
of the sampler that stored it: one whose contents fit this sampler is
installed, and what the sampler then holds of another family's state is
not promised. To start a sampler from a fit of another family, use
`installTrees`, which takes the donor's trees and, where this sampler
draws them, its `sigma`, `k` and DART state, and none of its latent
variables. A state's `k` is installed as recorded, whatever leaf prior
this sampler holds; the leaf sd is then `k.scale / k` against this
sampler's `k.scale`.

Assigning the field directly (`sampler$state <- newState`) does *not*
restore the sampler - it only overwrites the R-side cache, leaving the
engine untouched, so the next run continues from the engine's own state
rather than the assigned one. A reload, and a `copy` of a sampler whose
engine is gone, install the field as `setState` does with
`forceUpdate = TRUE`, having no other state to fall back on: a tree that
no longer fits is repaired, and nothing is reported. Such an object
whose field holds a state that is refused cannot be copied, and after a
reload raises that refusal at each use until a state of its own is
assigned; one whose field holds no state is refused the same ways. A
`copy` of a sampler whose engine exists takes its current chain and not
the field. Always route a restore through `setState`.

A state carries no data, a forest's basis among them, but it records the
level names of each forest a factor multiplies, as the sampler's data
object holds them (the attribute `basisLevels`, a list with one entry
per forest, `NULL` for a forest of no factor; absent where no forest
multiplies a factor). `setState`, a `copy` and a reload put each of the
state's coefficients on the recipient's level of the same name, so a
state installed in a sampler whose factor has the same levels in another
order, as after `relevel`, gives each level its own coefficient; the
recipient's levels and block are its own before and after. A state whose
names differ from the recipient's, forest by forest - another name,
another number of levels, or a factor where the recipient's forest
multiplies numbers or nothing - is refused with an error, whatever
`forceUpdate` says, and the sampler is left as it was. A state that
records no levels goes in by position.

A state never changes the sampler's model. The leaf prior, its `k` where
fixed, a fixed `sigma`, Student-t df, negative-binomial shape or DART
concentration, a fixed amplitude prior or amplitude, and a supplied
Gaussian-process lengthscale are the sampler's own, before and after; a
value the sampler draws is installed where the state holds one and kept
where it does not, so a state written before these blocks were dropped
still installs. The response transform a `k`-named prior is anchored to,
and which every leaf value is stored relative to, is not in the state
either: it is derived from the data when the sampler is created and when
`setResponse` or `setOffset` at `updateScale = TRUE`, `setWeights` or
`setActiveRows` at `updateScale = TRUE`, or a `setData` that derives it,
re-anchors it, and recorded as the `response.range` of the sampler's
data object (see
[`dbartsData`](https://vdorie.github.io/dbarts/reference/dbartsData.md)),
so a `copy` or a reload is re-created in it. A state goes in as stored:
its trees, leaf values, linear slopes, `k`, `sigma` and kept draws are
numbers on the sampler's internal scale, and an install - `setState`,
`copy`, a reload - reads them against this sampler's response mapping
and leaf standardization and never converts them. A state stored under
the mapping in force, as a sampler's own, a `copy`'s or a reload's is
unless a re-anchor came between, puts the chain back where it was
stored, value for value. A state stored before a re-anchor, or by a
sampler on another response, is legal and has no special meaning: the
same internal numbers read on this sampler's scale are a position the
chain could hold, and the next draws move the parameters where the data
say. A response that was only rescaled is then no change at all on the
internal scale, the fit, the leaf prior's spread and `sigma` all
rescaling with it; one whose center moved as well is fitted wrongly
until the draws catch up. A kept draw restored across mappings predicts
on the new scale, where a re-anchor converts the kept draws so that they
predict what they did. No leaf model refuses such a state; one whose
saved Gaussian-process draws were made under another leaf
standardization or other lengthscales than this sampler holds is
refused, a saved Gaussian-process draw replaying only under its own. The
state records the mapping it was read under, which no install reads.
Because a re-anchor is not in the state, restoring a state stored before
one does not undo it; to roll a re-anchoring
`setResponse(y2, updateScale = TRUE)` back, re-anchor to the old
response with `setResponse(y, updateScale = TRUE)` and then `setState`
the saved state, so that the response is restored with the state. The
cut points and their weights are the sampler's, and `setState` writes
neither: a state records the ones its thresholds were stored on (its
`cutPoints` and `cutMass` attributes) and goes onto the sampler's. Where
every threshold of the state is one of the sampler's cut points, as on a
grid that only gained points, the trees are the stored ones and
`setState` returns `TRUE` if they fit; where one is not, the state is
declined or, forced, has that split moved (above). Other weights on
equal cut points change no tree: the state goes in, `TRUE`, and the
chain runs on under the sampler's weights, which is not the chain the
state was stored from. So restoring a state does not undo a change of
cut points made by `setCutPoints` or by
`setPredictor(updateCutPoints = )`; 0.9-x's `setState` did. To restore
the cut points too, hold the state in a variable, since `setCutPoints`
replaces the `state` field under `updateState`, and set them first, the
whole set in one call: `sampler$setCutPoints(attr(st, "cutPoints"))` and
then `sampler$setState(st)`. A column whose cut points had been changed
and for which `attr(st, "cutMass")` holds an entry then has its points
back and not its weights: `setState` still returns `TRUE`, and the chain
runs with that column's cut points equally likely.
`setPredictor(x[, j], j, updateCutPoints = "value")` gives back both
while the column holds the values its cut points were derived from. A
`copy` and a reload, which create a new engine, take the cut points,
weights and leaf standardization the sampler's data object records,
which are the ones in force (see
[`dbartsData`](https://vdorie.github.io/dbarts/reference/dbartsData.md)),
and none from a state. A reload of a sampler saved after a change of cut
points with no `storeState()` since has the stored state's splits moved
onto the recorded cut points by value, as a forced `setState` moves
them, and its `state` field then holds what its engine holds; call
`storeState()` before saving after a change of cut points to have it
reload as it stands. A `copy` of a sampler whose engine exists takes its
current chain, which is on those cut points already. On a column the
sampler holds as an ordered factor, whose cut points its number of
levels fixes, a state recording other cut points is refused, by
`setState`, a `copy` and a reload alike. On a sampler made from a data
frame, `setData` reads the new data's factors by label and records the
sampler's own levels, so a `copy` and a reload hold the levels the
sampler was created with and take the state it stored. A linear or
Gaussian-process leaf's covariate standardization, each covariate's
center and scale and a Gaussian-process leaf's heuristic lengthscale, is
data too: derived when the sampler is created, by a `setData` that
derives it and by `setPredictor`, `setWeights` or `setActiveRows` with
`updateStandardization = TRUE`, and recorded as the `standardization` of
the sampler's data object, which a `copy` and a reload take. A state
records the standardization it was stored under, which no install reads:
the sampler keeps its own, so kept linear draws restored across a change
of standardization predict under the sampler's until the draws replace
them, as kept draws restored across a change of response mapping predict
on the new scale. A leaf covariate that held a single value when its
standardization was derived has no scale: the state and the record hold
`NA` for it, and the covariate is divided by 1.

A state brings back none of these values, so a sampler that changed one
and is to go back to a stored state has it put back first. The data
object held in a variable before the change, `d <- sampler$data`, is the
general answer: it records every one, the weights of the cut points
included, so `setData(d)` puts back the observations and every derived
value with them.

|  |  |  |
|----|----|----|
| value | changed by | put back by |
| cut points | `setCutPoints`; `setPredictor(updateCutPoints = )`; a `setData` that derives them | `setCutPoints(d@cutPoints)`, the whole set; or `setData(d)` |
| their weights | the same | `setPredictor(x[, j], j, updateCutPoints = "value")` while column `j` holds the values they came from; or `setData(d)` |
| leaf center, scale, lengthscale | `updateStandardization = TRUE` on `setPredictor`, `setWeights` or `setActiveRows`; a `setData` that derives them | `setPredictor` of the old values with `updateStandardization = TRUE`, under the weights and mask they were read over; or `setData(d)` |
| response range | `updateScale = TRUE` on `setResponse`, `setOffset`, `setWeights` or `setActiveRows`; a `setData` that derives it | the same setter on the old response, offset, weights or mask with `updateScale = TRUE`; or `setData(d)` |
| forest anchor and unit | nothing | \- |

Then `setState(st)` of a state `st` stored before the change installs it
and returns `TRUE`, and the chain is the one stored.

The case weights and an `aft` sampler's censoring status are not in the
state, and the state records nothing of either: the latents it carries
are the values its trees were drawn from, `setState` installs them and
draws none, and the next sweep draws them against the weights and
censoring status then in force. A weight change is undone by
`setWeights` with the old weights and `setState` with the stored state,
in either order, the sampler then continuing bit for bit as one that
installed the stored state under the old weights would, and `setState`
returns `TRUE` in both. A sampler saved after its weights changed
without a store is re-created from the last stored state under the
weights it then holds; it is a continuation of the stored chain and not
of the live sampler, carrying no sweep made since the store. Call
`storeState()` after the change to have it reload as it stands; a `copy`
of a sampler whose engine exists takes its current chain.

### Missing values in predictors

A predictor column can hold missing values from the first time the
sampler is given one in it - at creation, or later by `setPredictor`,
`setData` or
[`updatePredictorPerObservationJointly`](https://vdorie.github.io/dbarts/reference/updatePredictorPerObservationJointly.md) -
and from then on for the life of the sampler, whatever values the column
holds afterwards. Each split rule on such a column sends missing values
down one side, and the sampler draws that side as part of the model.
When a column first becomes able to hold missing values, every rule
already on it, in the current trees and in every saved draw
(`keepTrees`), has its side drawn at random, left or right with
probability one half each, from each chain's own generator; a change
that is refused undoes these draws with the rest of it. A rule on an
unordered factor column has its side drawn only where a missing value
can reach it: beneath a rule on the same column that sends missing
values down the other side none can arrive, and the rule holds no drawn
side, exactly as on a sampler whose column held missing values at
creation. A saved draw was made when no value in the column was missing,
so nothing it was drawn from bears on that side, and one half is that
draw's posterior for it. `predict` then routes a missing value in the
column through every saved draw, those made before its first missing
value and those made after. A column that has never held a missing value
has no side to route by, and `predict` refuses a missing value in it.
`copy` and a reload keep which columns can hold missing values.

A stored state records which columns could hold a missing value when it
was stored. When `setState`, a `copy` or a reload installs it into a
sampler where a further column can - a state stored before that column's
first missing value, or by a sampler whose column never held one - every
rule on that column, in the state's current trees and in its saved
draws, has its side drawn in the same way before the state is judged,
from the generators the state carries; a state stored just before a
single `setPredictor` call gave the column its first missing value, and
installed just after it, draws the sides that call drew. `setState`
returns `TRUE` for such an install unless a drawn side leaves a leaf
with no row or a
[`monotone`](https://vdorie.github.io/dbarts/reference/monotone.md) tree
out of order, where it returns `FALSE` as for any state that does not
fit (see ‘Saving’). A state stored before this record was kept draws
nothing, and its rules on such a column send missing values left. A warm
start (`installTrees`, or `warm.start` of
[`bart`](https://vdorie.github.io/dbarts/reference/bart.md)) from a fit
whose data could not hold a missing value in a column draws, from each
receiving chain's own generator, the side of every rule on that column
in the trees it installs.

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
price of R-level consistency. Predictor mutation has no counterpart in
the flat C interface (`dbarts.h`), which reaches the per-sweep
conditioning channels and leaves every whole-data channel to these
methods, so a client that mutates predictors every iteration is served
by data-frame input rather than by dropping below R.

### Infrastructure methods

`getPointer()` returns the sampler's underlying engine pointer, checking
it for validity first. If the engine the cached pointer refers to is no
longer live (typically after a
[`load`](https://rdrr.io/r/base/load.html)), it transparently re-creates
one from `control`, `model`, `data` and the stored `state`, reinstalls
every recorded forest weight and any recorded active-row mask (see
`reapplyForestWeights` below, since neither rides the saved state) so
the replacement matches a freshly built engine, and only then replaces
the cached pointer; a re-creation that fails leaves the sampler exactly
as it was. `reapplyForestWeights(ptr)` re-installs every per-forest
weight mirrored on the sampler onto a freshly (re-)created `ptr`, and
`reapplyActiveRows(ptr)` the mirrored active-row mask; they are called
only from `getPointer`, `setState` and `copy`, never directly.

`getPointer()` is the only R-level primitive for exchanging engine
pointers with code that drives the engine directly through `dbarts.h`:
it is where a compiled consumer obtains the handle that header's entry
points take, since that header declares no creation entry of its own.
There is no other R-level path to that boundary: the low-level `.Call`
wrappers this class's methods delegate to are internal to the package
and unexported, so a consumer below this reference class has `dbarts.h`
and nothing else.

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
names onto the installed `sampler$data@x`. When `x` itself has column
names they are matched to the sampler's, which it keeps, a name that is
not the sampler's being refused, and a sampler with no column names
takes `x`'s; when it does not (a bare matrix, or a vector reshaped
through `dim`), the sampler's existing column names are kept instead.
Either way, a later *by-name* column reference on that sampler -
`setPredictor`'s own `column` argument given as a character string, or
[`updatePredictorPerObservationJointly`](https://vdorie.github.io/dbarts/reference/updatePredictorPerObservationJointly.md)'s
shared-column match - continues to resolve against
`colnames(sampler$data@x)` as it did before the replacement.

### Prior draws

`sampleTreesFromPrior` and `sampleLeafParametersFromPrior` are
MEAN-forest entries by contract: each leaves the variance forest of a
heteroscedastic (`variance`) sampler exactly as it finds it, along with
`sigma`, `k` and the family's latent block.
`sampleVarianceForestFromPrior` is the variance forest's own, drawing
its tree structures from the tree prior conditioned on carrying no empty
leaf - the distribution the variance moves themselves price - and each
leaf's factor from the calibrated inverse-chi-squared scale prior, then
rebuilding \\s^2(x)\\ from the drawn factors. A caller wanting a whole
heteroscedastic chain at its prior calls all three; on a homoscedastic
sampler this one draws nothing and returns, so a loop over families
needs no guard.

### Warm starts

`installTrees` seeds the sampler's forests from a `donor` instead of
drawing trees from the prior, for scaling to more chains or embedding a
fit in a larger sampler. Only the donor's trees, and its `sigma`, `k`
and DART state where this sampler draws them, transfer: the sampler
keeps its own model, as `setState` does. The donor's trees and leaf
values, and its `k` and `sigma` where this sampler draws them, go in as
stored on the internal scale and are not converted to this sampler's
response: a donor fitted to a response on another scale seeds the same
internal numbers, which the next draws move where this sampler's data
say, and a donor whose response differs only by a rescaling seeds its
own fit, rescaled. The predictors' frame goes the same way. A `linear`
leaf's coefficients go in as stored, numbers per unit of this sampler's
own centered and scaled covariates, and this sampler keeps its own
center and scale, as a Gaussian-process leaf keeps its own center, scale
and lengthscale, whichever cut grid the donor is on; the donor's are not
read. Where the two standardize a covariate differently the seeded
forest is then not the donor's function of it, and the next draws move
it where this sampler's data say. A covariate that holds a single value
is read at the center it holds, so between two single values the seeded
fit is the donor's on the training rows. Each chain keeps its own
random-number stream and redraws everything else, so several chains
seeded from one donor stay overdispersed. DART weights go only between
two DART fits and a variance forest only where both have one: a donor's
that this sampler has no place for is left behind, and this sampler's
own, where the donor has none or where the donor's hold 0 on a column
this sampler's forest may split on, starts as a new fit's once the donor
goes in, DART drawing its weights from the trees taken at the first
sweep; neither is judged, and a donor that is declined leaves them as
they were. A donor holding a value no chain could run from, as
`setState` refuses one, is an error in both forms, and nothing is
installed. A donor of constant leaves into linear ones gives each leaf
its value as intercept and slopes of zero, the donor as it stands. A
donor of linear leaves into constant ones, or into linear ones on
another number of covariates, leaves its slopes behind, each leaf
starting at the intercept the donor stores. A donor of Gaussian-process
leaves where this sampler's are of another kind, or the reverse, or one
fit to another number of rows, is an error naming the two, as is a donor
fit to another number of predictors or one whose column is an unordered
factor where this sampler's is not, or the reverse. A donor of fewer
trees than this sampler is taken whole, its trees in its order and the
rest starting as single leaves at zero, or at one in a variance forest.
Under `blocks`, whatever the donor's tree count, each block in turn
takes, in the donor's order and up to its count, the donor's trees not
yet taken that it holds as they stand, the rest of its trees starting as
single leaves, and where that leaves a tree out, the same is tried once
more with the donor's single leaves offered last; a donor all of whose
trees go in so is taken whole, and one fitted with the same blocks and
counts goes in tree for tree. A donor of more trees, or under `blocks`
one whose trees do not all go in so, gives those whose fitted values
vary most over this sampler's rows, each row weighted by its case
weight, in the donor's order, the mean and the variance forest choosing
apart and each chain among the trees of the donor sample it takes; under
`blocks` each block in turn takes up to its count of the trees that vary
most within it as it holds them, never one it would hold as a single
leaf when the donor holds it as more, and starts the rest as single
leaves.

`installTrees(donor)` asks first. A donor every tree of which fits this
sampler as stored is installed and the call returns `TRUE`. One that
would have to be changed is not installed: the call returns `FALSE` and
the sampler is as it was, its trees, drawn values, generators, kept
draws and `state` field. A tree has to be changed when a split's
threshold is not one of this sampler's cut points, as for a donor fitted
to other rows; when a leaf holds none of this sampler's rows; when a
split is on a column, or in an interaction, its forest does not allow;
when it records a side for missing values on a column that has never
held one here; or when a monotone forest's leaf values are out of order.
A donor is also declined when it holds more trees than this sampler's
forest or, under `blocks`, trees the blocks cannot take as they stand,
when it holds linear slopes this sampler's leaves do not take, constant
leaves or linear ones on another number of covariates, and when its
Gaussian-process fits were recorded on other cut points.
`installTrees(donor, forceUpdate = TRUE)` installs such a donor and
returns `NULL` invisibly, whether or not anything was changed: a split
off the cut points moves to the nearest one inside the range the splits
above it leave or, with none left, becomes one leaf with everything
beneath it; a split with a side no row reaches, and the first split from
the root that its forest does not allow, likewise become one leaf, at
the mean of the leaves it replaces weighted by the rows each held; a
side is dropped; a monotone tree whose leaf values are out of order is
cut back where they are; Gaussian-process fits start at zero; slopes the
sampler does not take are left behind; and of a donor holding more
trees, or under `blocks` one whose trees the blocks cannot all take as
they stand, the trees whose fitted values vary most go in. The sampler's
cut points and their weights are the same before and after. A side drawn
for missing values is how a donor is read, not a change to it. On a
column this sampler holds as an ordered factor a donor recorded on other
cut points, a numeric column or a factor of another number of levels, is
an error in both forms. `bart`'s `warm.start` always forces.

A warm start biases the early draws toward the donor, so it shortens
burn-in rather than removing it; keep drawing a non-zero number of
burn-in samples before treating the chain as converged. Single-forest
samplers only: a sampler carrying two or more forests (a Bayesian causal
forest or any other amplitude model, and a multinomial sampler's K
category forests) refuses the call by name, reporting its forest count,
because the install is not tested at more than one forest - it
reassembles the trees from a saved sample but takes the amplitudes from
the donor's current state. Use `growFromRoot` to initialize such a
sampler instead. A variance forest is not one of the counted forests, so
a heteroscedastic sampler takes a warm start as any single-forest one
does.

`growFromRoot` instead builds the sampler's initial forest by
XBART-style root-down stochastic tree construction (He, Yalov and Hahn
2019): each of `n.sweeps` sweeps rebuilds every tree from the root,
sampling each split from the integrated-likelihood weight of every
candidate cut under the tree prior. This reaches a good fit in far fewer
sweeps than the exact sampler, so it is a fast starting point rather
than a posterior sampler - the exact MCMC sweeps own stationarity once
`run` begins, and the posterior is unchanged. It is available for the
constant-leaf model only; calling it on a `linear` or `gp` leaf prior is
an error, not a silent fall-back, so initialize those forests with
`sampleTreesFromPrior` instead. As with `installTrees`, the grown forest
biases the early draws toward its fit, so shorten burn-in rather than
skipping it. Each chain grows on its own random-number stream, so the
result is independent of `n.threads`. Unlike `installTrees` it runs on a
multi-forest sampler, each forest grown against its own residual inside
the sweep and the forests' amplitudes drawn with them. A composable
cross-sampler workflow follows for free: `donor$growFromRoot(k)` then
`target$installTrees(donor)`.

## Value

For `run`, a named-list with contents `sigma` (absent on a
heteroscedastic sampler, below), `train`, `test`, and `varcount` (plus
`k`, when the reported forest draws `k` in that run and absent when it
is fixed, and `varprobs` when applicable). `train` is an array of
dimension n.obs x n.samples x n.chains, and likewise `test` (or `NULL`
if the sampler has no test data) and `varcount` are n.predictors x
n.samples x n.chains; `sigma` is n.samples x n.chains. When `n.chains`
is `1` the trailing chain dimension is dropped, so `train` is a plain
n.obs x n.samples matrix and `sigma` a plain vector of length n.samples.
On a multi-forest (`forests`) sampler `varcount` gains a forest axis
between the predictors and the samples - n.predictors x n.forests x
n.samples x n.chains, forest-major within a draw and the forests in the
model's order - so each forest's own per-draw split counts arrive from
one call rather than only the reported forest's; a single-forest
sampler's array keeps exactly its n.predictors x n.samples x n.chains
shape, and the same widening is what a `family = "multinomial"`
sampler's per-category counts ride. This is the RAW run shape; the
packaged fit [`bart`](https://vdorie.github.io/dbarts/reference/bart.md)
builds reshapes it draws-first with the forest names on the trailing
margin, as it does for the multinomial channel.
`$getForestVariableCounts` reads the same quantity for the CURRENT
state, one forest at a time. A Bayesian causal forest adds two more:
`forestFits`, an n.obs x n.forests x n.samples x n.chains array of each
forest's fitted values on the internal scale (the forests in the model's
order), and `glue`, a sum(q) x n.samples x n.chains array of the
amplitudes each draw combines them through, stacked forest-major and as
wide as each forest's own basis (a Bayesian causal forest's \\(a, b_0,
b_1)\\, three rows), so both surfaces and their recombination come from
a single call; `train` carries the combination \\a \mu(x_i) + b\_{z_i}
\tau(x_i)\\ on the response scale, or on the latent scale when the
sampler was built under `"probit"` or `"logistic"`, and `test` is filled
with `NaN` (there is no test treatment vector to combine off-sample). No
other model reports either element. A `"nbinom"` sampler adds one of its
own: `shape`, the negative-binomial \\r\\ each draw is conditioned on,
shaped exactly as `sigma` (a length n.samples vector at one chain, an
n.samples x n.chains matrix otherwise) because it is the count analog of
it - fixed at the value the sampler was created with under a fixed
`shape`, and that sweep's grid draw otherwise. It is written from the
same state `storeState` serializes and consumes no random numbers, so
reading it costs a run nothing. No other family carries the element at
all: it is absent from the list, not `NULL` within it, so `run()$shape`
is `NULL` on every non-`"nbinom"` sampler and a test of the channel must
be `!is.null(...)` rather than a comparison, which `NULL` would satisfy
vacuously. A sampler built with `family = student()` adds `resid.df` on
exactly the same terms - the degrees of freedom \\\nu\\ each draw is
conditioned on, shaped as `sigma`, written from settled state and
consuming no random numbers, absent from the list under any other error
distribution - fixed at the value supplied to `student(df = )`, and that
sweep's grid draw when the degrees of freedom are estimated. A
heteroscedastic (`variance` forest) sampler carries no `sigma` element -
absent, as `shape` is off `"nbinom"` - since its engine holds the scalar
fixed; it adds `variance`, the per-observation variance surface
\\s^2(x_i)\\ on the original response scale, shaped as `train`, whose
square root is the residual scale, and `varianceTest`, the same at the
test rows (`NULL` without test data). A run can be interrupted with
`Ctrl-C`: it stops between iterations, or inside an iteration's
Polya-Gamma latent draws or a monotone leaf-order count, whose time
grows with a logistic, negative-binomial or multinomial fit's counts or
with a large tree, joining any worker threads first. It then raises R's
interrupt condition, class `c("interrupt", "condition")`, as base R's
own long computations do: nothing is printed as an error,
[`try()`](https://rdrr.io/r/base/try.html) and error handlers do not
catch it, `tryCatch(interrupt = )` does, and calling handlers, `on.exit`
and `finally` run. A loop of fits inside
[`try()`](https://rdrr.io/r/base/try.html) stops at the interrupted one.
The call returns no samples. The sampler's chains are left where they
stopped and it can be run again; where the stop falls inside a pass over
the trees (a slow monotone count, a latent refresh) that pass may be
partly applied, each latent at this iteration's draw or the previous
one, which is a valid state to run on from. The saved trees the stopped
call wrote are dropped: `predict` and `getTrees` serve the earlier
calls' draws it did not overwrite, oldest first. Unlike a base R
interrupt, it offers no `"resume"` restart: the run is over. Under
`control@keepFits == FALSE`, `train`, `test`, and (when the model
carries them) the variance and forest channels come back `NULL` rather
than an array - present in the list, holding nothing, unlike
`shape`/`resid.df` above, which are absent outright when inapplicable;
see `callback` above and
[`dbartsControl`](https://vdorie.github.io/dbarts/reference/dbartsControl.md)'s
`keepFits`. A `callback` that returns nonzero ABORTS the run the same
way an interrupt does: the results and the saved trees the call wrote
are discarded and `run` signals an error, naming that a callback stopped
it; that is an error, where an interrupt is not.

For `setCutPoints`, `NULL` invisibly. The change is always made, as in
0.9-x: there is no unforced form and no `forceUpdate` argument, and a
split left with no cut point between the splits above it is merged, so
the trees can change. It sets a grid, typically once before sampling.
Inside a larger sampler whose draws move a column's values, let the grid
follow them with `setPredictor(updateCutPoints = )`, which derives the
grid a creation would for the new values, up to `n.cuts`, each time, and
whose unforced form is refused, with the sampler as it was, when the
trees cannot hold the new grid.

For `setPredictor`, `TRUE` if the new predictor was installed and
`FALSE` if it was refused. An unforced update is refused, with no error
and no warning, when it would leave a leaf empty in any tree of any
forest of the sampler (see ‘Multi-forest and heteroscedastic predictor
mutation’) or, under a
[`monotone`](https://vdorie.github.io/dbarts/reference/monotone.md)
constraint, leave a tree's leaf values out of order. A refused update is
rolled back, whether it named columns or replaced the whole matrix: the
sampler's trees, fits, predictors and random number generator are as
they were before the call. A forced update always installs and returns
`NULL` invisibly: the value says whether a valid update was made, which
has no meaning when it is forced. When `forceUpdate` is `"partial"`,
instead returns a logical vector of length equal to the number of
observations, `TRUE` where that observation's new value was installed
and `FALSE` where it was rolled back to its previous value to keep every
tree of every forest valid.

For `installTrees`, `TRUE` if the donor was installed and `FALSE` if it
was not, the sampler then being exactly as it was: it is `FALSE` when a
tree a chain would take, or a leaf value, would have to be changed (see
‘Warm starts’). No warning is given. With `forceUpdate = TRUE` the donor
is installed with the changes and the value is `NULL`, invisibly,
whether or not anything was changed.

For `setState`, `TRUE` if the state was installed and `FALSE` if it was
not, the sampler then being exactly as it was. It is `FALSE` when a tree
of any chain would have to be changed to fit the sampler: a split whose
threshold is not one of the sampler's cut points, a leaf that no row of
the current predictors reaches, a split that the state holds outside the
range the splits above it leave, a rule with a side for missing values
on a column that has never held one, a split its forest's columns or
interaction limit do not allow, or monotone leaf values out of order
(see ‘Saving’). No warning is given. `TRUE` means the trees and leaf
values are the stored ones. It does not promise the cut points or their
weights, which stay the sampler's whatever the state records, so that a
chain restored onto other weights or a finer grid is not the stored
chain continued; the case weights and an `aft` censoring status, which
stay the sampler's; the generator, which is left as the sampler's own
when the state's is of another kind; or every kept draw, of which the
newest that fit the sampler's store are taken. Nor does it cover a value
the sampler holds fixed, which a state never changes. A state stored
under another response mapping is installed as stored and reports
`TRUE`. With `forceUpdate = TRUE` the state is installed, repaired where
it has to be, and the value is `NULL`, invisibly, whether or not
anything was repaired: a caller who needs to know calls without it
first. `copy` and a reload install a state that way and report nothing.

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

Its return shape mirrors `run()`'s own `train`/`test` channel for a
single-forest sampler, with two exceptions. An amplitude-coupled
(multi-forest, including a Bayesian causal forest) sampler refuses
`predict` outright rather than answering here - see
`getFitsWithoutOffset` below, which already records the same refusal -
so everything that follows describes a single-forest or multinomial
sampler only. For a single-forest sampler: with saved trees to replay
(`keepTrees` `TRUE` and at least one recorded draw), the result is an
n.test x n.samples x n.chains array, the trailing chain dimension
dropped to a plain n.test x n.samples matrix when `n.chains` is `1`;
without saved trees to replay, the CURRENT (live) trees are replayed
once per chain instead, with no samples axis at all - a plain vector of
length n.test at one chain, an n.test x n.chains matrix otherwise. A
`family = "multinomial"` sampler's surface reports the K category
probabilities instead: an n.test x K x n.samples (x n.chains) array, the
K axis inserted between the rows and the samples. Unlike the
single-forest case, this samples axis is never dropped: it is sized `1`
rather than omitted when there are no saved trees to replay. A
heteroscedastic (variance forest) single-forest sampler with saved trees
to replay returns not a bare array but a named list
`list(mean = ..., variance = ...)`, both members the single-forest shape
above; without saved trees there is nothing for the variance forest to
replay, so the bare mean array returns alone, exactly as a
non-heteroscedastic sampler's would.

For `getTrees`, a `data.frame` with one row per tree node in
depth-first, left-hand-side pre-order, with columns `forest` (indexed
from 1, present only on a sampler with more than one forest - a
multinomial one, or one declared with `forests =` naming several - and
then leading every other column, forest-major, at the default
`forest = NULL` or a vector; any other sampler returns no `forest`
column and accepts `forest = 1`), `chain` (present only when
`n.chains > 1`), `sample` (present only for saved samples, and reporting
the draw number asked for), `tree`, `n` (the number of observations in
the node), `var` (the splitting variable, or -1 at a leaf), and `value`
(the split value, or the leaf prediction). An ordinal rule's value is
its cut point and observations with values less than or equal to it go
left; on an ordered-factor column that cut point is an exact
half-integer at a level boundary, so `levels(x)[floor(value) + 1]` is
the last level going left (the split value is a 0-based code, R's
`levels` is 1-based); a categorical rule carries no data value (its
`value` is `NA`) and its split is reported in the `directions` column
instead. When the sampler has any categorical predictors the result
gains a `directions` column decoding each categorical rule into one
`"L"`/`"R"` character per level, in level order (level `k` goes right
when its character is `"R"`); ordinal rules and leaves are `NA`. When
any predictor contains missing values the result gains a `missing`
column giving the branch (`"L"`/`"R"`) each rule sends missing values
down; rules on complete columns and leaves are `NA`. Under a `linear`
leaf prior each leaf's `value` is its intercept and the result gains one
`beta.<column>` column per designated covariate holding that leaf's
slope on the internal standardized scale; internal nodes are `NA`.

For `getSigmas`, a numeric vector of length equal to the number of
chains, giving each chain's current residual standard deviation on the
original response scale, or `NULL` on a heteroscedastic sampler, whose
residual scale is the surface `getVariance` reports.

For `getK`, each chain's current `k`, the value `run()$k` records per
draw when `k` is drawn in that run (a run records no `k` while it is
fixed), read without running, as `getSigmas` reports `sigma`: after a
run it is bitwise the last draw. A fixed `k` is reported once per chain,
and a forest whose scale a calibration map sets reports 1; a forest with
no basis stating neither `sd` nor a coefficient's prior reports the leaf
prior's `k`, and those that share one prior total (see
[`forest`](https://vdorie.github.io/dbarts/reference/forest.md)) each
the group's one. It is `k` whatever terms the prior was named in,
relative to `getLeafPrior()$k.scale`: under an `sd`-named prior, that
scale over the spread in force. A numeric vector of length equal to the
number of chains at one forest, or, at the default `forest = NULL` on a
multi-forest sampler, an n.forests x n.chains matrix; a single-forest
sampler's `NULL` read is bitwise its `forest = 1` read.

For `getShape`, the shape parameter of the sampler's family currently in
force, a numeric vector of length equal to the number of chains, and
`NULL` on a family with none - the count analog of `getSigmas`. Only the
negative-binomial family has one today, its \\r\\: the variance is
\\\mu + \mu^2 / r\\, a larger shape is closer to Poisson, and \\r\\ is
`size` in `rnbinom` and `theta` in MASS and mgcv. It is the same scalar
`run()$shape` records once per kept draw, read mid-sweep and without
serializing state, so an outer sampler driving this one a sweep at a
time reads it here; `state[[chain]]$shape` holds it only where it is
drawn. Because the refusal for a family carrying no shape is a `NULL`
and not an error, a caller distinguishing “no shape” from a value must
test `!is.null(...)`; comparing to an expected number would pass
vacuously.

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
other \\K-1\\ forests move within the same sweep, so the code driving
the sampler, reading it between sweeps, reads a quantity whose
conditioning it cannot reproduce, and the composition recipes reach the
same posterior through the fits and the category offset without touching
it. That decline is scoped to the augmentation families whose latent is
a PRECISION refreshed against a moving margin; the location-valued
latents below - probit's and ordinal's \\z_i\\, and aft's imputed log
survival time - are reported as documented. Otherwise the values the
last sweep's trees were drawn from. Each sweep draws them first, after
any change made between sweeps and before its trees, so between a setter
and the next run they are the values from before the setter (`setData`
starts them as creation does, an ordinal or negative-binomial sampler
from the kept cutpoints or shape; an `aft` `setResponse` leaves a
censored row at its new bound) and the next sweep draws them afresh.
They are a plain vector of length equal to the number of observations
when there is a single chain, or an observations-by-chains matrix
otherwise, written into `result` when one was supplied. **What that
variable IS depends on the family and is not uniform.** It is a
LOCATION, on the sampler's own latent scale, for `"probit"` (the
truncated normal \\z_i\\), `"ordinal"` (the same \\z_i\\ under the
ordinal thresholds) and `"aft"` (the imputed log survival time), all of
which the code driving the sampler regresses on directly. It is a
PRECISION, one per observation, for `"logistic"` and `"nbinom"` (the
Polya-Gamma \\\omega_i\\) and for a Student-t residual distribution (the
scale-mixing \\\lambda_i\\); these WEIGHT a working response and are not
on the response scale at all, so differencing them against a fit is
meaningless. Note the last case: a sampler whose `family` is
`"gaussian"` but built with `family = student()` (see
[`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md)) does
report latents, and they are precisions.

For `getFitsWithoutOffset`, the sampler's combined per-observation
location on the response scale and *without* the installed offset, an
n.observations x n.chains matrix. `run()$train` reports the same
quantity with the offset folded in, so this value plus the installed
offset is that one; there is deliberately no with-offset twin, since the
offset is the caller's own. It is the read an outer sampler driving this
one a sweep at a time needs: an incremental scheme updates its outer
block against \\f(x_i)\\, and taking `getLatents()` minus `run()$train`
instead is BIASED: a location latent is on the offset-inclusive scale,
as the training channel is, so that difference is the residual alone and
carries none of the offset's signal. A block that owns the offset
regresses `getLatents() - getFitsWithoutOffset()` on its covariate. On a
multi-forest (`forests`) sampler this is the only route, short of
running the sampler, from R to the combined fit - `getForestFits` gives
one forest's internal-scale totals and `predict` is refused. It is
refused, naming the reason, on a multinomial sampler, whose reported
channels are per-category softmax probabilities rather than one additive
location; `predict(x)` serves that read, but reports the SAVED samples
rather than the current state when the sampler was built with
`keepTrees`.

For `getVariance`, a heteroscedastic (`variance` forest) sampler's
current variance surface \\s^2(x_i)\\ on the original response scale, an
n.observations x n.chains matrix at the default `test = FALSE` and an
n.test x n.chains matrix at `test = TRUE`. It is the mid-sweep read of
the channels `run()` records as `variance` and `varianceTest`: at the
state a recorded sweep left, the accessor and the channel agree exactly.
What it reports is a VARIANCE, so a residual scale is its square root -
it is the surface analog of `getSigmas`, which reports the scalar
\\\sigma\\ a homoscedastic sampler carries. `NULL` where those channels
report nothing: on a homoscedastic sampler, and at `test = TRUE` with no
test rows installed. Unlike `predict`, which addresses the saved
samples, it reads the trees currently in force, so it answers after
`$sampleVarianceForestFromPrior()` and needs no `keepTrees`.

For `getForestFits`, a multi-forest sampler's requested forest's current
internal-scale fitted values, an n.observations x n.chains matrix, or,
at the default `forest = NULL`, every forest stacked with the forest
margin between the observations and the chains, n.observations x
n.forests x n.chains (bitwise the `forest = 1` read on a single-forest
sampler). For `getForestAmplitudes`, the named forest's amplitudes - one
per basis column - as a q x n.chains matrix, or, at the default
`forest = NULL`, every forest's stacked forest-major into a sum(q) x
n.chains matrix. A drawn one has prior median size 1, on a factor the
difference between two levels' amplitudes, a level's own having 0.71.
The vector is *ragged*: a Bayesian causal forest's forest `1` carries
the single \\a\\ on its implicit intercept and its forest `2` the pair
\\(b_0, b_1)\\, so the stacked read is the \\(a, b_0, b_1)\\ of \\y = a
\mu(x) + b_z \tau(x) + \epsilon\\. For `getForestVariableCounts`, the
requested forest's current per-predictor split counts, an n.predictors x
n.chains integer matrix, or, at the default `forest = NULL`, every
forest stacked the same way as `getForestFits`, n.predictors x n.forests
x n.chains; its rows are named by the training matrix's column names
when that matrix carries them (on margin 1 in both shapes), and left
unnamed otherwise.

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
`$getLeafPrior()` supplying the rest. At the fit level it is packaged:
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

For `getLeafPrior`, the leaf prior a forest runs under, alone, as a
named list: `leaf.prior`, the specification in the terms it was named
in - a `normal()`, `linear()` or `gp()` object (see
[`dbartsPriors`](https://vdorie.github.io/dbarts/reference/dbartsPriors.md))
carrying one of `k`, a number or a `chi()` prior, or `sd`, a number or
an `invchi()` prior, the family default when none was named - which goes
back into `setLeafPrior` or a fitting function's `leaf.prior` as is;
`leaf.model`, one of `"constant"`, `"monotone"`, `"linear"`, or `"gp"`;
`prior.sd.of`, what the sd is the sd of: `"forest"`, `"coefficient"`
(linear, per standardized covariate) or `"amplitude"` (gp);
`prior.mean`; `k.scale`, the value `k` is relative to, so the spread in
force on each chain is `k.scale / getK()` - the data's scale, which the
family and the response transform fix and no spelling of the prior
moves, so that under an `sd`-named prior `k` is `k.scale` over the sd;
and `response.scale` and `response.shift`. A forest the calibration map
does not set also carries, after `leaf.prior`, `sd`, its size in force
in the response's units - `k.scale` over `k` at a fixed `k`, the number
as stated where an `sd` was, and `NA` where `k` is drawn or the `sd` is
under `invchi()`, which have no one number - and `sd.stated`, whether
the leaf prior is named by an `sd`, at creation or through
`setLeafPrior`. A fixed `k` is the one in force, read off the engine; an
`sd`, a number or under `invchi()`, is reported as it was named. Every
chain runs under the sampler's one prior and one response transform,
which no state install moves, so every value is shared by the chains. A
drawn `k` is chain state, read by `getK`. At the default `forest = NULL`
a multi-forest sampler returns a list of one prior per forest, named by
the forests' labels, a multinomial sampler's by the response's levels; a
single-forest sampler's `NULL` read is bitwise its `forest = 1` read. An
`NA` spread is refused by `setLeafPrior`, `setModel` and the fitting
functions. On a forest whose scale a multi-forest calibration map sets,
`k` is pinned at 1 and `leaf.prior` is the statement of the forest's
size:
[`forest`](https://vdorie.github.io/dbarts/reference/forest.md)`(sd = v)`,
`v` the number as it was written, where an `sd` was stated, and
[`forest()`](https://vdorie.github.io/dbarts/reference/forest.md) where
none was. Its `sd` is the one `setForest` states. The list then begins
with `leaf.prior`, `sd`, the forest's `sd` in force in the response's
units - `v` itself, bit for bit, or the default's multiple of the
response's scale times the unit the sampler derived at creation -
`sd.stated`, `TRUE` or `FALSE`, `multiplier`, what the forest
multiplies, `"none"`, `"factor"` or `"numeric"`, and `amplitude`,
`fixed()` where the forest's coefficient is held, `student(df)` where it
is drawn under the Student-t prior its creation stated, `normal()` where
a forest with no basis is drawn under the normal prior its creation
stated, and absent where it is drawn under its default prior.
`prior.sd.of` says what `sd` is the size of: `"forest"` for a forest
with no basis, `"level difference"` for a factor and
`"forest, per unit of basis"` for a numeric column, the size exactly
where the coefficient is held, with `" (at a coefficient of 1)"` where
it is drawn. A numeric forest's `sd` is per unit of its column and
carries the column's name where the basis has one, and the list adds
`basis.scale`, where the data object records one, the standard deviation
of the column the default was taken from at creation, named as `sd` is.
`k.scale` is the forest's own prior sd before its coefficient, in the
response's units: the forest's `sd`, on every forest but a sharer. A
forest with no basis stating neither `sd` nor a coefficient's prior,
alone or sharing one prior total with others (see
[`forest`](https://vdorie.github.io/dbarts/reference/forest.md)), reads
as one forest does, under the fitting function's leaf prior:
`leaf.prior` is that prior as named, `sd` the forest's part of the total
in force, `NA` where `k` is drawn, `sd.stated` `FALSE`, the sd being the
fitting function's and not its own (a one-forest sampler reads `TRUE`
when its leaf prior is named by an sd), `multiplier` `"none"`,
`amplitude` `fixed()`, `prior.sd.of` `"forest"` alone and
`"forest, shared"` in a group, one left of a group included, and
`k.scale` over the `k` its part.

On a forest whose scale the multi-forest CALIBRATION MAP sets, on a
sampler built with `forests =` or `dbartsData(bases = )` (see
[`forest`](https://vdorie.github.io/dbarts/reference/forest.md)), `k` is
pinned at 1, `k.scale` is the map's leaf scale, and the entries above
that the map adds are absent (so `NULL`) on every other forest, a
multinomial one included. A drawn coefficient is a multiplier with prior
median size 1 - on a factor the difference between two levels'
coefficients, a level's own having 0.71 - and on a numeric column
through bcf's 0.674.

The spread (`k.scale / k`) describes the LEAF-PARAMETER scale of the
forest total. It equals the prior standard deviation of \\f(x)\\ at
every \\x\\ for the constant leaf only; for the other three the prior of
\\f(x)\\ is x-dependent and the spread bounds it in a leaf-specific
direction. Under `"linear"` it is a LOWER bound attained at the
standardized covariate origin, with \\sd(f(x)) = \\ the spread
\\\sqrt{1 + \\z(x)\\^2}\\ in the internally standardized leaf covariates
(a missing value maps to \\z_j = 0\\); `prior.mean` is exact. Under
`"gp"` it is an UPPER bound over \\x\\, attained at rows that reproduce
a leaf member and on over-cap leaves, and elsewhere decaying to zero as
\\x\\ leaves the leaf's data cloud, where every prior draw equals
`prior.mean` exactly. Under `"monotone"` it is a LOWER bound in the
interior - the realized standard deviation runs a few per cent to about
twenty per cent above it - and `prior.mean` is NOT the prior mean of
\\f(x)\\ under an active constraint: that marginal is skew, with an
x-dependent mean tracking the constraint's direction across several the
spread (see `monotone` in
[`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md)).

This is the authoritative reader of the prior in force. A named `sd` is
absolute: the sampler holds it as `k = k.scale / sd`, or as the `chi`
prior on `k` at that scale, and restates it whenever the response
transform is re-anchored - `setResponse`, `setOffset`, `setWeights` or
`setActiveRows` at `updateScale = TRUE`, `setData`, and `setResponse`
and `setOffset` made by compiled code through `dbarts.h`, which has no
entry for the other two - so a fixed sd, and the prior on a drawn one,
keep their meaning in response units. A drawn sd's current value is
`k`'s: at a re-anchor it stays where it was on the internal scale,
stretching with the leaves, and lags until the next draw, as a drawn `k`
does. A `k` is relative to the response transform recorded at creation
or at the last re-anchor, and moves only with it. A state install
(`setState`, `installTrees`, `copy`, a reload) leaves the prior as it
was.

For `setLeafPrior`, `NULL` invisibly. It changes the leaf prior's spread
or the hyperprior it is drawn under, on every chain, and nothing else:
not the response transform, not `sigma` (which `setModel` re-pins on a
fixed-sigma gaussian sampler), not the tree prior or a DART split prior.
Into a drawn prior - a `chi` or `invchi` prior - each chain keeps the
spread in force at the call, whatever the prior before it, fixed or
drawn and in either the `k` or the `sd` form: `k.scale` is the data's
under every prior, so `k` is not touched, `getK` does not move, and the
new prior acts from the next draw of `k`. A fixed `k` or `sd` stated
sets the spread as written, a fixed `sd` as `k = k.scale / sd`. Because
the write keeps the chains' spread, on a sampler that has not yet run it
is not the start that creation under the new prior would make. The write
takes effect on the next sweep, reinterpreting no leaf value already
drawn; a write that reproduces what is already in force is skipped
bitwise, so a read followed by a write cannot perturb a draw. It is
total over the four leaf models, and a DART sampler is served. A
multinomial sampler takes `normal(k = )` and applies the `k` to every
category forest, leaving the softmax map's leaf scale; a named sd or a
`k` hyperprior is refused there, as at creation. On a sampler whose
forests carry amplitudes the `leaf.prior` is the fitting function's,
written on the forests it reaches as creation under it gives them. Where
one forest has no basis it is written there whatever that forest states,
a stated `sd` going: a `k` or a named `sd` as at creation, and on a
forest stating a coefficient's prior, which has no `k`, the size the
prior gives at the scale creation set - `normal(sd = v)` the size `v`,
`normal(k = v)` that scale over `v`, a `chi()` prior on `k` that scale
over its median, once, with no `k` drawn there - leaving the
coefficient's prior and every leaf and coefficient drawn. Where several
forests have no basis, those sharing one prior total take it, each
keeping the part of it its creation gave it, so after a forest has left
the group (see `setForest`) the sizes no longer add up to it, and those
stating their own `sd` keep it. Where none shares, as where every forest
has a basis or every forest with no basis states its `sd`, the default
restated changes nothing and any other `leaf.prior` is refused, naming
`setForest(forest, sd = )`. Into a drawn `k` a forest starts at the
spread in force. A write changes the model. The write is recorded where
re-creation reads it. The map's k.scale is recorded at creation too, so
a re-creation after `setResponse` or `setOffset` at
`updateScale = FALSE` states every forest against the same k.scale as
the live sampler. A value that is not a single positive finite number is
an error. A heteroscedastic sampler's variance forest is a separate leaf
model, which `setLeafPrior` does not reach; `setVarianceForest` changes
its settings. `setModel` changes everything else.

For `getForest`, every value in force on one forest as a plain named
list, no class and no attribute, in `setForest`'s order and with the
names its kind takes: `sd`, the size `getLeafPrior` reports with no
column's name on it, `NULL` where `k` is drawn or no one positive number
is in force; `power` and `base`; `amplitude`, the coefficient's prior as
`fixed()`, `normal()` or `student(df)`, a forest with a basis stating
none reading `normal()`, and a forest with no coefficient `fixed()`; and
`proposal.probs`, `n.perturb.cuts` and `treeShift`, a forest stating
none reading the control's. A value the forest takes from above reads
like any other. `do.call(s$setForest, c(list(f), s$getForest(f)))`
changes nothing, whether the forest states a value included. On a
multinomial sampler the list has no `sd` or `amplitude` and covers every
category forest, whichever is named. Tree counts and columns are not in
it: `print` and `getTrees` show them.

For `setForest`, `NULL` invisibly. It changes the forest named, and on
it only the names given, each meaning what it means in
[`forest`](https://vdorie.github.io/dbarts/reference/forest.md) on that
forest at creation. It acts from the next sweep and leaves the trees,
leaf values and coefficients already drawn as they are; a write changes
the model. A name left out or `NULL` changes nothing, and so does a
value equal to the one in force, which leaves the forest following what
it followed: a forest taking its depth from the fitting function goes on
taking it, a forest sharing one prior total stays in the group, and a
numeric forest's default goes on being read from its column. A forest
that states its `sd` or its depth is changed by stating another; nothing
returns it to following the fitting function, short of a new sampler.

`sd` on a model of one forest is its size, written as
`setLeafPrior(normal(sd = ))` writes it, in the terms of its leaf model.
On a model of several it sets the forest's trees' scale: its size where
the coefficient is held and its size at a coefficient of 1 where it is
drawn, on a factor the size of the difference between two levels and on
a numeric column the size per unit. A forest with no basis alone under
the fitting function's leaf prior leaves that prior at an `sd` stated
for it, its `k` then pinned at 1. A forest sharing one prior total with
others leaves the group the same way, and those still sharing keep their
parts of the total, set once at creation, and their sizes; nothing puts
a forest back in the group but a new sampler. A multinomial sampler
refuses `sd` as its creation does. `power` and `base` are the forest's
depth, every category's on a multinomial sampler. `amplitude` changes
the coefficient's prior between held, `fixed()`, and drawn, `normal()`
or `student(df)`, and between those two, on every chain. A hold sets the
coefficient to its held value at the call - 1, or 0 at a factor's first
level and 1 at its second - so `getForestAmplitudes` reads it at once; a
release draws it from that value from the next sweep; a Student-t
prior's mixing variance starts at the normal prior's and is kept where
only the degrees of freedom change. No leaf and no trees' scale moves,
except on a forest whose size the fitting function's leaf prior gives: a
forest with no basis holding the whole of that prior, released, takes
the size creation gives the statement, its `k` pinned at 1, and held
again follows the leaf prior again, a drawn `k` starting at the spread
in force; the last forest left of a group, released, states its part of
the total in force as its `sd`, which `getForest` then reads. `normal()`
writes a forest with a basis back to its default and `fixed()` one with
no basis. The prior in force changes nothing, and `fixed()` on a forest
with no coefficient, as a model of one forest or a sharer has none.
`student(df)` or `normal()` is refused as at creation there, a sharer
taking it with an `sd` in the same call, which takes the forest out of
the group, and where the leaf prior would give the size and gives none,
as a `chi()` prior on `k` with an infinite scale. A multinomial sampler
takes `fixed()`, which changes nothing, as its creation does, and
refuses `student(df)` and `normal()` as its creation does.
`proposal.probs`, `n.perturb.cuts` and `treeShift` are the forest's own,
and a later `setControl` leaves them.

`vars`, `n.trees`, `interactions` and `blocks` are refused by name
whatever their value, being part of the model, and `basis` is refused
naming `setForestBasis`, a basis being data. What creation refuses of a
forest's statement is refused here for the same reason. Every check
comes before any write, so a refused call changes nothing, its good
names included. A multinomial model's forests take one statement:
`s$setForest(power = 4)` reaches every category, with `forest` left out,
and a position or a level is refused. Every change is recorded on the
model, as creation with the statement records it, so a copy and a reload
run what the sampler runs.

Where `setForest` and `setTreePrior`, `setLeafPrior` or `setModel` reach
one setting of one forest, as on a model of one forest, the later write
wins, whatever the forest stated. A heteroscedastic sampler's variance
forest takes the mean forest's depth at creation and is on its own after
that: `setForest`, `setTreePrior` and `setModel` leave it, and
`setVarianceForest` changes it.

For `setTreePrior`, `NULL` invisibly. It changes the fitting function's
tree prior, as `power`, `base` and `split.probs` mean at creation,
through `setModel`'s route, so the two cannot differ. A name left out or
`NULL` changes nothing, and so does a value equal to the one in force.
`power` and `base` reach the forests with no basis that state none of
their own and, where one forest has no basis, or on a multinomial
sampler's category forests, those forests whatever they state, the later
write winning; they make no forest state them, and nothing returns a
forest that states its own to following, short of a new sampler. A name
that reaches no forest is refused, naming the forests that state it and
`setForest`. A `split.probs` vector changes the values and `1` restores
equal ones; on a DART sampler `dart(...)` changes DART's settings, a
setting left out taking creation's default; `dart()` on an ordinary
sampler and a vector on a DART one change the model's structure and are
refused, and a sampler of several forests takes equal probabilities
only. A variance forest stays as it is.

For `getVarianceForest`, every value in force on the variance forest as
a plain named list, `power`, `base`, `proposal.probs` and
`n.perturb.cuts`, or `NULL` on a sampler with none. For
`setVarianceForest`, `NULL` invisibly. It changes those settings by
`setForest`'s rules, each meaning what it means in
[`varianceForest`](https://vdorie.github.io/dbarts/reference/varianceForest.md).
`vars` and `n.trees` are refused by name, and a sampler with no variance
forest refuses the call: whether there is one is part of the model.

For `storeState`, `NULL` invisibly; it is called for its side effect of
capturing the sampler's current engine state into the serializable
`state` field, which [`save`](https://rdrr.io/r/base/save.html) then
writes out. See ‘Saving’. What it writes is one draw's trees, about 13
bytes per node per tree per chain, which is 3.5 pct of a single draw's
training predictions at 20,000 observations and 75 trees and falls like
\\1/n\\; it exceeds those predictions only below a few hundred
observations, and it is a small fraction of a whole run's `yhat.train`,
not a multiple of it.

## See also

[`updatePredictorPerObservationJointly`](https://vdorie.github.io/dbarts/reference/updatePredictorPerObservationJointly.md)
for applying a shared per-observation predictor update across several
samplers at once.

[`samplePriorPredictive`](https://vdorie.github.io/dbarts/reference/samplePriorPredictive.md)
for repeated prior draws of the hyperparameters, `sampleTreesFromPrior`,
`sampleLeafParametersFromPrior` and `predict` on a private sampler, for
calibrating priors before fitting.

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
#>  num [1:100, 1] -0.709 0.304 0.291 -0.472 -0.126 ...

## an outer Gibbs step changes the response (e.g. after updating some
## other part of a joint model); the sampler picks the new target up on
## the next run() without being recreated
newY <- y + rnorm(n, 0, 0.05)
sampler$setResponse(newY)
samples <- sampler$run()
str(samples$train)
#>  num [1:100, 1] -0.828054 0.520308 0.337472 -0.489366 0.000512 ...
```
