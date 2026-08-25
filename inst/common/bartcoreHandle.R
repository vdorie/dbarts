# Low-level bartcore handle wrappers, moved out of the package namespace:
# each is a thin .Call onto one C_dbarts_bartcore_* entry point, exercising
# the handle layer directly rather than through a dbartsSampler. They exist
# for tests and benchmarks only - the package itself drives the handle
# through the R5 sampler's own methods - so keeping them out of the
# namespace means there is no dbarts::: path to them; a caller below the R5
# class has only the shipped C header. Reached with dbarts:::C_dbarts_* for
# the .Call targets and dbarts:::<name> for the three package-internal
# validators the bodies still use (asCountMatrix, validateCategoryOffset,
# validateLiveScale). bartcoreRun, bartcorePredict and bartcoreSetModel have
# live in-package callers and stay defined in R/bartcore.R; the aliases
# below just bind those definitions under the same names so a caller here
# cannot see two divergent copies.

# Replaces a multinomial sampler's response: an n x K matrix of nonnegative
# integer counts in the same layout the count creation entry takes (column k is
# category k), with trials n_i = sum_k counts[i, k] >= 1. n and K are fixed at
# creation - every combiner buffer is sized by n, and K is the forest count - so
# only the values may change. The trees carry over, fitted to the previous
# counts, exactly as a single-forest setResponse leaves them; the next run forms
# every category's working response against the new counts.
#
# The sweep draws n_i Polya-Gamma variates per observation per category, so
# replacing single-trial labels with grouped counts multiplies sweep cost by
# mean(n_i).
#
# Validated R-side (safe over fast); the engine re-derives the trials and
# re-checks the invariants, and refuses everything before it installs anything,
# so a rejected matrix leaves the sampler untouched.
bartcoreSetCounts <- function(bcSampler, counts) {
  counts <- as.matrix(counts)
  if (!is.numeric(counts)) {
    stop("multinomial counts must be a numeric matrix of non-negative integers")
  }
  if (anyNA(counts)) {
    stop("multinomial counts must not contain missing values")
  }
  if (any(counts < 0)) {
    stop("multinomial counts must be non-negative")
  }
  if (any(counts != round(counts))) {
    stop("multinomial counts must be whole numbers")
  }
  if (any(rowSums(counts) < 1)) {
    stop("every multinomial count row must have at least one trial (n_i >= 1)")
  }
  counts <- dbarts:::asCountMatrix(counts)
  # counts is column-major (category-major), the layout the combiner reads
  invisible(.Call(dbarts:::C_dbarts_bartcore_setCounts, bcSampler$ptr, counts))
}

# Installs (or clears, at NULL) a multinomial sampler's n x K category offset:
# the latent becomes f_ik + o_ik, so the offset enters the log-sum-exp margins,
# each category's working response and the reported softmax probabilities, and
# never a leaf value. n and K are fixed at creation, as they are for the counts.
#
# This is the response-side counterpart of bartcoreSetCounts, not
# bartcoreSetOffset: the response model's offset is added to every reported
# channel AFTER the K forests are blended, which for a softmax is the wrong side
# of the nonlinearity, and a flat per-observation offset is the softmax's own
# null direction in any case. Only the row-centred part is identified - adding a
# constant to a whole row of the matrix leaves every reported probability
# unchanged - and the entrance leaves the input as given rather than re-centring
# it.
#
# This one shifts the TRAIN latent only. The test rows are other rows, so they
# carry their own offset (bartcoreSetCategoryTestOffset), and predict takes its
# own matrix per call; neither is derived from this one.
bartcoreSetCategoryOffset <- function(bcSampler, offset) {
  # the shape check needs the handle's own K, which only a multinomial handle
  # carries; off one, the entry's capability probe is the refusal, and it names
  # the family situation rather than a shape
  if (!is.null(bcSampler$K)) {
    offset <- dbarts:::validateCategoryOffset(
      offset,
      nrow(bcSampler$x),
      bcSampler$K
    )
  }
  invisible(.Call(
    dbarts:::C_dbarts_bartcore_setCategoryOffset,
    bcSampler$ptr,
    offset
  ))
}

# Installs (or clears, at NULL) a multinomial sampler's nTest x K category test
# offset: the recorded test channel becomes softmax(f_test + o_test), formed
# where the train blend forms softmax(f + o). The test fits enter no
# likelihood, so this moves the reported test probabilities and nothing else -
# no draw, no working response, no train channel.
#
# Its rows are the CURRENT test rows. Replacing those rows while it is
# installed is refused rather than silently reinterpreted (clear it first);
# out-of-sample predict does not read it at all, taking its own matrix for the
# rows it is given.
bartcoreSetCategoryTestOffset <- function(bcSampler, offset.test) {
  # as for the train offset, the shape check needs the handle's own K; off a
  # multinomial handle the entry's capability probe is the refusal. The ROW
  # count belongs to the sampler's test store, which the handle does not track,
  # so the engine is what pins it - here only K and finiteness are checked
  if (!is.null(bcSampler$K) && !is.null(offset.test)) {
    offset.test <- as.matrix(offset.test)
    offset.test <- dbarts:::validateCategoryOffset(
      offset.test,
      nrow(offset.test),
      bcSampler$K,
      "category test offset"
    )
  }
  invisible(.Call(
    dbarts:::C_dbarts_bartcore_setCategoryTestOffset,
    bcSampler$ptr,
    offset.test
  ))
}

# Forest `forest`'s (0-based) amplitude basis, an n x q numeric matrix; re-forms
# every multiplier and both residuals on the next run. The sole basis-mutation
# route, at any forest and any width: the amplitudes are preserved and remapped,
# and a width-preserving install is the bitwise identity on all of them.
bartcoreSetForestBasis <- function(bcSampler, forest, basis) {
  basis <- as.matrix(basis)
  storage.mode(basis) <- "double"
  invisible(.Call(
    dbarts:::C_dbarts_bartcore_setForestBasis,
    bcSampler$ptr,
    as.integer(forest),
    basis
  ))
}

# A per-forest, per-observation weight s on forest `forest` (0 prognostic, 1
# treatment): a multiplicative precision factor on that forest's own leaf
# conditionals, composing with the case weights and the active-row mask so
# the forest's draws see (w_i * a_i) * m_f^2 * s_i. It does not remove the
# row from occupancy, from the empty-leaf veto, from the combination (the
# row still receives m_f f_f(x_i)), or from the residual sigma degrees of
# freedom; s_i = 0 says only that row i carries no information about forest
# f, and that forest's leaves over such rows stay well-defined prior draws.
#
# Two edges, or a consumer is misled. At s_i = 0 with a nonzero multiplier only
# the WEIGHT is zeroed - the response stays the reparameterized residual - so
# the reported-fit exactness an exactly zero multiplier buys does not follow
# this channel. And the weight lives on the sampler, not in its state, so a
# pipeline that REBUILDS a sampler and restores a stored state silently drops
# the weight and fits a different model while the states still agree.
bartcoreSetForestWeights <- function(bcSampler, forest, weights) {
  weights <- as.double(weights)
  if (!all(is.finite(weights)) || any(weights < 0)) {
    stop("forest weights must be finite and non-negative")
  }
  invisible(.Call(
    dbarts:::C_dbarts_bartcore_setForestWeights,
    bcSampler$ptr,
    as.integer(forest),
    weights
  ))
}

# Forest `forest`'s (0-based) amplitudes on the combining response, a
# q_forest x n.chains matrix, or - at forest = NULL - the whole vector stacked
# forest-major, sum q_f x n.chains, which is the shape the run's own `glue`
# channel carries. Ragged forest by forest, which is why a forest is named:
# bcf's forest 0 carries a, its forest 1 the pair (b0, b1).
bartcoreForestAmplitudes <- function(bcSampler, forest = NULL) {
  .Call(
    dbarts:::C_dbarts_bartcore_getForestAmplitudes,
    bcSampler$ptr,
    if (is.null(forest)) NULL else as.integer(forest)
  )
}

# A forest's internal-scale function values (0 prognostic, 1 treatment),
# n.observations x n.chains. For a multinomial handle (forest = category,
# 0-based) this is the OFFSET-FREE totalFits: with a category offset
# installed, softmax_k(bartcoreForestFits(bc, k) + offset[, k]) reproduces
# the reported train channel, not softmax_k(bartcoreForestFits(bc, k)) alone.
bartcoreForestFits <- function(bcSampler, forest) {
  .Call(
    dbarts:::C_dbarts_bartcore_getForestFits,
    bcSampler$ptr,
    as.integer(forest)
  )
}

# The combined per-observation location on the RESPONSE scale and without the
# offset, n.observations x n.chains - the quantity the run's train channel
# carries with the offset folded in. Refused on a multinomial handle, whose
# reported channels are per-category softmax probabilities and not one additive
# location; bartcorePredict(bc, x) serves that read instead, and reports the
# saved samples rather than the current state under keepTrees.
bartcoreFitsWithoutOffset <- function(bcSampler) {
  .Call(dbarts:::C_dbarts_bartcore_getFitsWithoutOffset, bcSampler$ptr)
}

# A forest's leaf-prior calibration in RESPONSE units, one row per chain
# (columns prior.scale, prior.sd, prior.mean, k, k.has.hyperprior,
# response.scale, response.shift, then the calibration map's own
# amplitude.prior.variance, amplitude.prior.scale, node.scale.factor,
# node.scale.divisor and basis.row.norm, NaN off the map) with the leaf model
# on a "leaf.model" attribute. Total over forests: a combiner's forests report
# the calibration its own map fixed. Forest is 0-based here, as everywhere on
# this layer; the R-level $getCalibration indexes from 1.
bartcoreForestCalibration <- function(bcSampler, forest) {
  .Call(
    dbarts:::C_dbarts_bartcore_getCalibration,
    bcSampler$ptr,
    as.integer(forest)
  )
}

# Restates a forest's leaf prior on every chain so the forest total's prior sd
# at k = 1 is prior.scale, response units. Refused when a combiner owns the
# calibration; nothing else moves, and a write reproducing what is in force is
# skipped bitwise.
bartcoreSetForestPriorScale <- function(bcSampler, forest, prior.scale) {
  invisible(.Call(
    dbarts:::C_dbarts_bartcore_setCalibration,
    bcSampler$ptr,
    as.integer(forest),
    dbarts:::validateLiveScale(prior.scale, "prior.scale")
  ))
}

# A forest's current per-predictor split counts (0 prognostic, 1 treatment),
# n.predictors x n.chains; the per-forest analog of the run's varcount channel.
bartcoreForestVariableCounts <- function(bcSampler, forest) {
  .Call(
    dbarts:::C_dbarts_bartcore_getForestVariableCounts,
    bcSampler$ptr,
    as.integer(forest)
  )
}

bartcoreSetOffset <- function(bcSampler, offset, updateScale = FALSE) {
  if (!is.null(offset)) {
    offset <- as.double(offset)
  }
  invisible(.Call(
    dbarts:::C_dbarts_bartcore_setOffset,
    bcSampler$ptr,
    offset,
    as.logical(updateScale)
  ))
}

bartcoreSetResponse <- function(bcSampler, y, updateScale = FALSE) {
  invisible(.Call(
    dbarts:::C_dbarts_bartcore_setResponse,
    bcSampler$ptr,
    as.double(y),
    as.logical(updateScale)
  ))
}

# A per-observation 0/1 mask saying which rows are in the data set for this
# sampler this sweep. An inactive row leaves every sufficient statistic, every
# family-level parameter update and its own latent draw, but keeps its leaf
# occupancy and still receives a fitted value. Absolute and independent of the
# case weights: the family serves w * a in either call order.
#
# NULL clears. An all-ones mask reports success and installs nothing - the
# opposite of setWeights, where all-ones installs - because one channel is
# membership and the other precision. An all-zeros mask is accepted and runs,
# with every forest at its prior. Values other than 0 and 1 refuse the whole
# call and install nothing.
bartcoreSetActiveRows <- function(bcSampler, active) {
  if (!is.null(active)) {
    active <- as.double(active)
  }
  invisible(.Call(
    dbarts:::C_dbarts_bartcore_setActiveRows,
    bcSampler$ptr,
    active
  ))
}

bartcoreSetWeights <- function(bcSampler, weights) {
  invisible(.Call(
    dbarts:::C_dbarts_bartcore_setWeights,
    bcSampler$ptr,
    as.double(weights)
  ))
}

bartcoreSetTestOffset <- function(bcSampler, offset.test) {
  if (!is.null(offset.test)) {
    offset.test <- as.double(offset.test)
  }
  invisible(.Call(
    dbarts:::C_dbarts_bartcore_setTestOffset,
    bcSampler$ptr,
    offset.test
  ))
}

bartcoreSetData <- function(bcSampler, data) {
  invisible(.Call(dbarts:::C_dbarts_bartcore_setData, bcSampler$ptr, data))
}

bartcoreSetTestPredictor <- function(bcSampler, x.test) {
  x.test <- as.matrix(x.test)
  storage.mode(x.test) <- "double"
  invisible(.Call(
    dbarts:::C_dbarts_bartcore_setTestPredictor,
    bcSampler$ptr,
    x.test
  ))
}

bartcoreSetPredictor <- function(
  bcSampler,
  x,
  forceUpdate = FALSE,
  updateCutPoints = FALSE
) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  result <- .Call(
    dbarts:::C_dbarts_bartcore_setPredictor,
    bcSampler$ptr,
    x,
    as.logical(forceUpdate),
    as.logical(updateCutPoints)
  )
  # track the current predictors R-side (the engine keeps none) when installed
  if (isTRUE(result)) {
    bcSampler$x <- x
  }
  result
}

# Overwrite columns of the current predictors; the engine keeps no matrix, so
# the wrapper maintains the R-side copy the re-quantize entry points read.
bartcoreUpdatePredictor <- function(
  bcSampler,
  x,
  columns,
  forceUpdate = FALSE,
  updateCutPoints = FALSE
) {
  columns <- as.integer(columns)
  result <- .Call(
    dbarts:::C_dbarts_bartcore_updatePredictor,
    bcSampler$ptr,
    as.double(x),
    columns,
    as.logical(forceUpdate),
    as.logical(updateCutPoints)
  )
  if (isTRUE(result) && !is.null(bcSampler$x)) {
    bcSampler$x[, columns] <- matrix(as.double(x), nrow(bcSampler$x))
  }
  result
}

bartcoreUpdatePredictorPerObservation <- function(bcSampler, x, column) {
  x <- as.double(x)
  column <- as.integer(column)
  installed <- .Call(
    dbarts:::C_dbarts_bartcore_updatePredictorPerObservation,
    bcSampler$ptr,
    x,
    column
  )
  if (!is.null(bcSampler$x)) {
    bcSampler$x[installed, column] <- x[installed]
  }
  installed
}

bartcoreUpdatePredictorPerObservationJointly <- function(
  bcSamplers,
  x,
  columns
) {
  x <- as.double(x)
  columns <- as.integer(columns)
  installed <- .Call(
    dbarts:::C_dbarts_bartcore_updatePredictorPerObservationJointly,
    lapply(bcSamplers, function(s) s$ptr),
    x,
    columns
  )
  for (k in seq_along(bcSamplers)) {
    if (!is.null(bcSamplers[[k]]$x)) {
      bcSamplers[[k]]$x[installed, columns[k]] <- x[installed]
    }
  }
  installed
}

# cutPoints is a list of strictly increasing numeric vectors, one per column;
# trees are refreshed unconditionally, collapsing splits the new cuts orphan
bartcoreSetCutPoints <- function(bcSampler, cutPoints, columns) {
  if (!is.list(cutPoints)) {
    cutPoints <- list(cutPoints)
  }
  cutPoints <- lapply(cutPoints, as.double)
  invisible(.Call(
    dbarts:::C_dbarts_bartcore_setCutPoints,
    bcSampler$ptr,
    cutPoints,
    as.integer(columns),
    bcSampler$x
  ))
}

bartcoreGetLatents <- function(bcSampler) {
  .Call(dbarts:::C_dbarts_bartcore_getLatents, bcSampler$ptr, NULL)
}

# Each forest's own RAW fits at new rows: an n.new x n.forests x n.samples
# (x n.chains) array on the forests' INTERNAL scale, the off-sample twin of
# bartcoreForestFits. No amplitude glue, no response transform and no offset are
# folded in - an amplitude coupling's location is
# shift + sum_f (basis_f %*% glue_f) * (scale * f_f), and off the training rows
# only the caller knows the bases, so the whole recombination is theirs.
# offset.test exists to be refused by name, since a shift belongs to that
# recombination rather than to any one forest's total. Refused on a handle whose
# coupling reports no per-forest fits, which is every single-forest sampler and
# a multinomial one (its raw f_k are defined, but bartcorePredict reports the
# softmax probabilities that handle's surface is stated in).
bartcorePredictPerForest <- function(bcSampler, x.test, offset.test = NULL) {
  x.test <- as.matrix(x.test)
  storage.mode(x.test) <- "double"
  .Call(
    dbarts:::C_dbarts_bartcore_predictPerForest,
    bcSampler$ptr,
    x.test,
    offset.test
  )
}

bartcoreGetTrees <- function(
  bcSampler,
  chainNums,
  sampleNums = NULL,
  treeNums,
  current = FALSE,
  newdata = NULL,
  forest = 0L
) {
  if (!is.null(newdata)) {
    newdata <- as.matrix(newdata)
    storage.mode(newdata) <- "double"
  }
  # forest is 0-based, as for bartcoreForestFits (0 prognostic, 1 treatment)
  .Call(
    dbarts:::C_dbarts_bartcore_getTrees,
    bcSampler$ptr,
    as.integer(chainNums),
    if (is.null(sampleNums)) NULL else as.integer(sampleNums),
    as.integer(treeNums),
    as.logical(current),
    newdata,
    bcSampler$x,
    as.integer(forest)
  )
}

bartcoreStoreState <- function(bcSampler) {
  .Call(dbarts:::C_dbarts_bartcore_storeState, bcSampler$ptr)
}

bartcoreSetState <- function(bcSampler, state) {
  invisible(.Call(
    dbarts:::C_dbarts_bartcore_setState,
    bcSampler$ptr,
    state,
    bcSampler$x
  ))
}

bartcoreSetModel <- dbarts:::bartcoreSetModel
bartcoreRun <- dbarts:::bartcoreRun
bartcorePredict <- dbarts:::bartcorePredict
