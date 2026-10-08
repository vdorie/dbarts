# The bartcore engine behind dbartsSampler. The dbartsSampler methods
# delegate to the bartcoreSamplerRun and bartcoreSamplerSet* functions below;
# the C side borrows vectors and pins them in the external pointer's
# protection slot.
#
# Not supported (methods error): weights with binary responses (weighted
# probit has no coherent latent-variable form), and setControl changes to
# anything fixed at creation (chain/tree counts, generators, and the cut
# grid).
#
# keepTrees/getTrees/predict use test-pinned formats. State serialization
# (storeState/setState, or runs with updateState) produces an engine-specific
# opaque object; restoring it into a sampler over the same data - including
# the transparent re-creation getPointer performs after save/load - continues
# the chains bitwise identically.

# A public dbartsSampler created through the forests = spec branch carries
# its forests' amplitude bases on data@bases, mirroring data@weights; this
# is the R-level capability probe, cheaper than
# a round trip through the bridge's own (totalAmplitudes-based) one. It is a
# CAPABILITY test - "carries amplitudes" - and deliberately not a forest count:
# a K-forest multinomial carries several forests and no amplitudes at all, so a
# numForests probe would misfire on it, as the bridge and the flat C entry each
# record independently.
samplerCarriesAmplitudes <- function(sampler) {
  !is.null(sampler$data@bases)
}

# Capability-specific wording for a mutation the bridge refuses through a guard
# shared with every multi-forest model (refuseMultiForestMutation and its
# siblings in R_interface_bartcore.cpp also cover the multinomial creation
# route, so their own message cannot name the amplitudes by itself). Raised
# R-side, before the .Call, so an amplitude-carrying sampler never reaches the
# bridge's generic "multi-forest" phrasing. The message names the CAPABILITY
# rather than the argument that declared it: dbartsData(bases = ) reaches the
# same samplers without a forests = declaration.
refuseAmplitudeMutation <- function(sampler, what, ...) {
  if (samplerCarriesAmplitudes(sampler)) {
    stop(
      what,
      " does not support a sampler that carries forest amplitudes: ",
      ...
    )
  }
}

# The multinomial capability probe on a SAMPLER, the counts analog of
# samplerCarriesAmplitudes: a K-forest softmax sampler is exactly the one whose
# data object carries the n x K count response.
# A CAPABILITY test, deliberately not a forest count, for the reason the
# amplitude probe is one - the two multi-forest models are indistinguishable by
# numForests. The read goes through the migration guard, so a sampler restored
# from a fit saved before the slot existed answers FALSE rather than raising.
samplerCarriesCounts <- function(sampler) {
  !is.null(dataCounts(sampler$data))
}

# Capability-specific wording for a channel the K-forest softmax gives no
# meaning to. Raised R-side, before the .Call, so a multinomial sampler never
# reaches the bridge's generic multi-forest phrasing, and so the R-canonical
# message can name the R5 method that DOES serve the caller. The message names
# the CAPABILITY, never a C entry point: those are not callable from R.
refuseCountsMutation <- function(sampler, what, ...) {
  if (samplerCarriesCounts(sampler)) {
    stop(what, " is not available on a multinomial sampler: ", ...)
  }
}

# The inverse probe: a channel only the K-forest softmax has, named on a
# sampler that carries no count response for it to write.
requireCountsCapability <- function(sampler, what) {
  if (!samplerCarriesCounts(sampler)) {
    stop(
      what,
      " is not available on a sampler that carries no count response: only a ",
      "multinomial (softmax) sampler has one"
    )
  }
}

# The three multinomial response channels, R5-side. Each mirrors its argument
# into the data object as setWeights mirrors weights: those slots are what
# CREATION reads, so getPointer's re-creation branch, setState's, and a
# save/load round trip all carry the current value with no reapply step of
# their own. Validation is R-side (safe over fast) and total before the .Call,
# which itself refuses before it installs anything, so the mirror runs only on
# a write that took.
bartcoreSamplerSetCounts <- function(sampler, counts) {
  current <- dataCounts(sampler$data)
  # K is the forest count, so a wrong category count is a different sampler,
  # not a malformed matrix; ask that before any per-row invariant, which a
  # truncated matrix would fail first and misattribute
  if (NCOL(counts) != ncol(current)) {
    stop("'counts' must have ", ncol(current), " categories")
  }
  counts <- validateMultinomialCounts(
    counts,
    nrow(current),
    seq_len(nrow(current))
  )
  ptr <- sampler$getPointer()
  .Call(C_dbarts_bartcore_setCounts, ptr, counts)
  sampler$data@counts <- counts
  # y is the trials the counts imply, never an independent quantity
  sampler$data@y <- as.double(rowSums(counts))
  # after the mirror, so a warning promoted to an error leaves them in step
  warnZeroTrials(counts)
  invisible(ptr)
}

# The sampler's own category names are its counts' column names, or "1".."K"
# (synthesized) when the counts carry none; the matching rule is
# alignCategoryColumns's.
alignSamplerCategoryColumns <- function(offset, counts, argument) {
  nms <- colnames(counts)
  synthesized <- is.null(nms)
  if (synthesized) {
    nms <- as.character(seq_len(ncol(counts)))
  }
  alignCategoryColumns(offset, nms, argument, synthesized)
}

bartcoreSamplerSetCategoryOffset <- function(sampler, offset) {
  current <- dataCounts(sampler$data)
  offset <- validateDataCategoryOffset(
    offset,
    nrow(current),
    seq_len(nrow(current)),
    "offset"
  )
  if (!is.null(offset) && ncol(offset) != ncol(current)) {
    stop("'offset' must have ", ncol(current), " categories")
  }
  offset <- alignSamplerCategoryColumns(offset, current, "offset")
  ptr <- sampler$getPointer()
  .Call(C_dbarts_bartcore_setCategoryOffset, ptr, offset)
  sampler$data@offset.category <- offset
  invisible(ptr)
}

bartcoreSamplerSetCategoryTestOffset <- function(sampler, offset.test) {
  current <- dataCounts(sampler$data)
  offset.test <- validateDataCategoryOffset(
    offset.test,
    NULL,
    NULL,
    "offset.test"
  )
  if (!is.null(offset.test) && ncol(offset.test) != ncol(current)) {
    stop("'offset.test' must have ", ncol(current), " categories")
  }
  offset.test <- alignSamplerCategoryColumns(
    offset.test,
    current,
    "offset.test"
  )
  ptr <- sampler$getPointer()
  .Call(C_dbarts_bartcore_setCategoryTestOffset, ptr, offset.test)
  sampler$data@offset.category.test <- offset.test
  invisible(ptr)
}

# Validates a 'callback' argument's shape - NULL, or a list carrying 'fn'
# and 'context' elements, 'fn' an external pointer and 'context' an external
# pointer or NULL - and returns list(fn = , context = ) with both NULL when
# 'callback' itself is NULL, the shape .Call's two separate arguments take.
# Nothing else is checked: the address is dereferenced exactly as handed, on
# the chain's own worker thread, so a callback pointing at the wrong kind of
# function crashes the session with no condition to catch.
validateCallback <- function(callback) {
  if (is.null(callback)) {
    return(list(fn = NULL, context = NULL))
  }
  if (!is.list(callback) || !all(c("fn", "context") %in% names(callback))) {
    stop(
      "'callback' must be NULL or a list with 'fn' and 'context' elements"
    )
  }
  fn <- callback$fn
  context <- callback$context
  if (is.null(fn) || typeof(fn) != "externalptr") {
    stop("'callback$fn' must be an external pointer")
  }
  if (!is.null(context) && typeof(context) != "externalptr") {
    stop("'callback$context' must be an external pointer or NULL")
  }
  list(fn = fn, context = context)
}

# Drives a dbartsSampler (the R-level sampler layer), reading its control
# defaults and delegating through its external pointer; cf. bartcoreRun,
# which drives a low-level bartcore handle directly.
bartcoreSamplerRun <- function(
  sampler,
  numBurnIn,
  numSamples,
  callback = NULL
) {
  control <- sampler$control
  numBurnIn <- coerceOrError(numBurnIn, "integer")
  numSamples <- coerceOrError(numSamples, "integer")
  if (length(numBurnIn) != 1L) {
    stop("'numBurnIn' must be a single integer", call. = FALSE)
  }
  if (length(numSamples) != 1L) {
    stop("'numSamples' must be a single integer", call. = FALSE)
  }
  if (is.na(numBurnIn)) {
    numBurnIn <- control@n.burn
  }
  if (is.na(numSamples)) {
    numSamples <- control@n.samples
  }
  if (is.na(numSamples)) {
    stop("bartcore engine samplers require 'numSamples' to be specified")
  }
  # as 0.9-x refused them: a negative count would otherwise reach the engine
  # as a wrapped size and report draws it never recorded
  if (numBurnIn < 0L) {
    stop(
      "number of burn-in steps must be greater than or equal to 0",
      call. = FALSE
    )
  }
  if (numSamples < 0L) {
    stop("number of samples must be greater than or equal to 0", call. = FALSE)
  }
  if (numBurnIn == 0L && numSamples == 0L) {
    stop("either number of burn-in or samples must be positive", call. = FALSE)
  }

  resolved <- validateCallback(callback)

  # keepFits is a per-RUN argument at the .Call boundary, but this layer's
  # own callers all read it off the control slot rather than passing an
  # independent value: the sampler's run method takes 'callback' but not
  # its own 'keepFits', so the control's own value is always what a run
  # honours
  result <- .Call(
    C_dbarts_bartcore_run,
    sampler$getPointer(),
    numBurnIn,
    numSamples,
    resolved$fn,
    resolved$context,
    control@keepFits
  )
  # a burn-only run returns NULL, or an empty list carrying only the
  # slow-count tally
  if (length(result) == 0L) {
    warnOnSlowCount(result)
    return(invisible(NULL))
  }
  warnOnGPFallback(result)
  warnOnSlowCount(result)
  result
}

# A Gaussian-process leaf larger than max.leaf.size is scored and drawn as a
# CONSTANT leaf: the fit stays coherent, but over most of the data it is not
# the model that was asked for, and nothing else says so. The engine counts
# every GP leaf evaluation and every one that took that path, and attaches the
# pair to the run; above a quarter of evaluations this says so, with the share.
# The counts ride the result either way, so a caller who wants a different
# threshold reads them directly. One site, so both front doors and a bare
# sampler are covered; it therefore fires once per run() call, which is once
# per fit for bart() and once per step for a sampler driven in a Gibbs loop.
warnOnGPFallback <- function(result) {
  tally <- attr(result, "gp.fallback")
  if (is.null(tally) || tally[["evaluations"]] <= 0) {
    return(invisible(NULL))
  }
  share <- tally[["fallbacks"]] / tally[["evaluations"]]
  if (share <= 0.25) {
    return(invisible(NULL))
  }
  warning(warningCondition(
    sprintf(
      paste0(
        "%.1f%% of Gaussian-process leaf evaluations fell back to a ",
        "constant leaf because the leaf held more than 'max.leaf.size' ",
        "observations, so much of this fit is not a Gaussian process; ",
        "raise the cap with gp(max.leaf.size = ), which costs roughly ten ",
        "times per doubling, or use FEWER trees, which grows deeper trees ",
        "and so smaller leaves"
      ),
      100.0 * share
    ),
    class = c("dbartsGPFallbackWarning", "dbartsWarning")
  ))
  invisible(NULL)
}

# Under monotone(prior = "leaf") every structure move counts its tree's leaf
# order, with no limit, so a large tree's count can take seconds and a lot of
# memory. The engine records the counts over about a second and attaches them
# to the run; this warns once per run() call, and runWithBurnIn merges its two
# runs' tallies so one fit warns once. The tally rides the condition.
warnOnSlowCount <- function(result) {
  tally <- attr(result, "slow.count")
  if (is.null(tally) || tally[["counts"]] <= 0) {
    return(invisible(NULL))
  }
  warning(warningCondition(
    sprintf(
      paste0(
        "%d monotone leaf-order count%s took more than about a second ",
        "(the slowest %.1f seconds, over %d leaves): under prior = \"leaf\" ",
        "a large tree's leaf order is costly to count, in time and in ",
        "memory; use more trees, which keeps trees small, or ",
        "monotone(prior = \"joint\"), which counts nothing"
      ),
      as.integer(tally[["counts"]]),
      if (tally[["counts"]] == 1) "" else "s",
      tally[["slowest.seconds"]],
      as.integer(tally[["slowest.leaves"]])
    ),
    class = c("dbartsSlowCountWarning", "dbartsWarning"),
    tally = tally
  ))
  invisible(NULL)
}

# What the bridge evaluates once a cancelled run has joined its workers and
# left nothing live: the interrupt condition, as base R signals one, handed to
# any calling or exiting handler for "interrupt". Only an interrupt nothing
# handled gets what follows, as R's own would: options("interrupt"), the blank
# line, options("error") (unless an interrupt function took its place), then
# the first of the "browser", "tryRestart" and "abort" restarts. (The poll's
# own check runs under an exiting handler, so R has done none of that already.)
# It is not an error, so try() and error handlers leave it alone. Unlike a
# real one it offers no "resume" restart: the run is over.
signalInterrupt <- function() {
  signalCondition(structure(class = c("interrupt", "condition"), list()))
  handler <- getOption("interrupt")
  if (is.function(handler)) {
    handler()
  }
  cat("\n", file = stderr())
  # as R's own, an options("interrupt") function stands in for options("error")
  handler <- getOption("error")
  if (!is.function(getOption("interrupt")) && !is.null(handler)) {
    eval(handler, globalenv())
  }
  for (restart in computeRestarts()) {
    if (restart[[1L]] %in% c("browser", "tryRestart", "abort")) {
      invokeRestart(restart)
    }
  }
}

# Sums slow-count tallies: counts add, the slowest count wins.
mergeSlowCountTallies <- function(a, b) {
  if (is.null(a)) {
    return(b)
  }
  if (is.null(b)) {
    return(a)
  }
  merged <- if (b[["slowest.seconds"]] > a[["slowest.seconds"]]) b else a
  merged[["counts"]] <- a[["counts"]] + b[["counts"]]
  merged
}

# Resolves a character 'column' against source's colnames into a 1-based
# integer index (or indices); NULL or an already-numeric 'column' passes
# through unchanged. 'what' names source for the not-found message. A missing
# index is refused by name here, ahead of the range checks it would otherwise
# reach as a bare missing condition.
resolveColumnIndex <- function(source, column, what) {
  if (anyNA(column)) {
    stop("'column' contains missing values", call. = FALSE)
  }
  if (is.null(column) || !is.character(column)) {
    return(column)
  }
  if (is.null(colnames(source))) {
    stop(
      "column names not specified at initialization, so cannot be ",
      "replaced by name"
    )
  }
  column <- match(column, colnames(source))
  if (anyNA(column)) {
    stop("column name not found in names of ", what)
  }
  column
}

# A column update addressing a column the training design coded from a factor
# takes that column's labels: a factor, character vector or sparseFactor,
# matched by label against the training levels and installed as the engine's
# 0-based codes, as a whole-frame update installs them. A label the column
# does not declare, and a missing value where the column has none (no route
# was learned for one), are refused by name; so is a number, which could only
# be read as a code, and a coded matrix or container. A data frame addresses
# several columns; other columns pass through as numbers.
codeCategoricalColumnUpdate <- function(x.train, x, column) {
  factorLevels <- attr(x.train, "factor.levels")
  if (is.null(factorLevels) || !is.numeric(column) || anyNA(column)) {
    return(x)
  }
  coded <- vapply(
    column,
    function(j) {
      j >= 1L && j <= length(factorLevels) && !is.null(factorLevels[[j]])
    },
    FALSE
  )
  if (!any(coded)) {
    return(x)
  }
  columnNames <- colnames(x.train)
  label <- function(j) {
    if (is.null(columnNames)) {
      paste0("column ", j)
    } else {
      paste0("column '", columnNames[j], "'")
    }
  }
  values <- if (is.data.frame(x)) {
    as.list(x)
  } else if (length(column) == 1L && is.null(dim(x))) {
    list(x)
  } else {
    NULL
  }
  if (is.null(values) || length(values) != length(column)) {
    stop(
      label(column[which(coded)[1L]]),
      " is categorical; give its values as a factor or character vector of ",
      "its labels, or a data frame for several columns"
    )
  }
  result <- matrix(NA_real_, NROW(values[[1L]]), length(column))
  for (k in seq_along(column)) {
    value <- values[[k]]
    j <- column[k]
    if (!coded[k]) {
      result[, k] <- as.double(value)
      next
    }
    if (
      !is.factor(value) &&
        !is.character(value) &&
        !methods::is(value, "sparseFactor")
    ) {
      stop(
        label(j),
        " is categorical; give its values as a factor or character vector ",
        "of its labels, not numbers"
      )
    }
    labels <- as.character(value)
    codes <- match(labels, factorLevels[[j]]) - 1L
    unknown <- unique(labels[!is.na(labels) & is.na(codes)])
    if (length(unknown) > 0L) {
      stop(
        label(j),
        " has ",
        if (length(unknown) > 1L) "labels " else "label ",
        quotedNameList(unknown),
        " not among its training levels"
      )
    }
    if (
      anyNA(labels) &&
        !sourceColumnHasNA(x.train, j, ncol(x.train), nrow(x.train))
    ) {
      stop(label(j), " has missing values, which its training values do not")
    }
    result[, k] <- as.double(codes)
  }
  if (length(column) == 1L) result[, 1L] else result
}

# The joint row-by-row update's values for its one shared column, as the one
# vector every sampler installs. One vector is one level in every sampler only
# if they hold the column with the same level table, so samplers that hold it
# coded from a factor with different levels are refused, for numbers as for
# labels. Then labels - a factor, character vector or sparseFactor - are
# matched to the levels by codeCategoricalColumnUpdate, in its words, with a
# refusal raised by a later sampler naming it; a number for a categorical
# column is refused in the words setPredictor uses, as is anything else that
# is not labels, where as.double would read a logical as two codes. Labels for
# a column only some samplers hold as a factor are refused. A column no
# sampler holds as a factor takes numbers, and a factor, a sparseFactor or
# text that is not numerals is refused for it.
codeJointColumnUpdate <- function(samplers, x, columnIndices, columnName) {
  levelTables <- lapply(seq_along(samplers), function(i) {
    factorLevels <- attr(samplers[[i]]$data@x, "factor.levels")
    if (columnIndices[i] <= length(factorLevels)) {
      factorLevels[[columnIndices[i]]]
    }
  })
  categorical <- !vapply(levelTables, is.null, FALSE)
  isLabels <- is.factor(x) ||
    is.character(x) ||
    methods::is(x, "sparseFactor")
  if (!any(categorical)) {
    if (
      is.factor(x) ||
        methods::is(x, "sparseFactor") ||
        (is.character(x) &&
          any(is.na(suppressWarnings(as.double(x))) & !is.na(x)))
    ) {
      stop("column '", columnName, "' is numeric and cannot take labels")
    }
    return(x)
  }
  first <- which(categorical)[1L]
  for (i in which(categorical)) {
    if (!identical(levelTables[[i]], levelTables[[first]])) {
      stop(
        "the samplers hold column '",
        columnName,
        "' with different levels (sampler ",
        i,
        " differs from sampler ",
        first,
        "), so one value would be a different level in each; update them in ",
        "separate calls, or create them with the same levels in the same order"
      )
    }
  }
  if (is.numeric(x) && is.null(dim(x)) && !all(categorical)) {
    stop(
      "column '",
      columnName,
      "' is categorical in sampler ",
      first,
      " and not in sampler ",
      which(!categorical)[1L],
      ", so numbers cannot be installed in both; update them in separate calls"
    )
  }
  if (!isLabels || !is.null(dim(x))) {
    stop(
      "column '",
      columnName,
      "' is categorical; give its values as a factor or character vector of ",
      "its labels",
      if (is.numeric(x) && is.null(dim(x))) ", not numbers"
    )
  }
  if (!all(categorical)) {
    stop(
      "column '",
      columnName,
      "' is categorical in sampler ",
      first,
      " and not in sampler ",
      which(!categorical)[1L],
      ", so its labels cannot be installed in both; update them in separate ",
      "calls"
    )
  }
  codes <- NULL
  for (i in seq_along(samplers)) {
    codesHere <- tryCatch(
      codeCategoricalColumnUpdate(
        samplers[[i]]$data@x,
        x,
        columnIndices[i]
      ),
      error = function(e) {
        stop(
          conditionMessage(e),
          if (length(samplers) > 1L) paste0(" (sampler ", i, ")"),
          call. = FALSE
        )
      }
    )
    if (is.null(codes)) {
      codes <- codesHere
    }
  }
  codes
}

## What a change of a column's cut grid does to the grid and to the splits
## on it, in the words setPredictor's updateCutPoints and setCutPoints' splits
## take and the codes the bridge reads: "none" keeps the grid, "position"
## keeps each split's position on the new one, rescaled when its count
## changes, and "value" moves each split to the new point nearest its
## threshold.
cutPointRuleCodes <- c(none = 0L, position = 1L, value = 2L)

## One word of choices, matched as match.arg matches; anything else is
## refused by the argument's name.
matchCutPointRule <- function(x, choices, name) {
  if (is.character(x) && length(x) == 1L && !is.na(x)) {
    matched <- pmatch(x, choices)
    if (!is.na(matched)) {
      return(choices[[matched]])
    }
  }
  stop(
    "'",
    name,
    "' must be one of ",
    paste0('"', choices, '"', collapse = ", "),
    call. = FALSE
  )
}

## setPredictor's updateCutPoints as one of its three words. 0.9-x took a
## logical, TRUE re-deriving the grid with every split left on its position:
## one is still taken, with a warning once in a session.
resolveUpdateCutPoints <- function(updateCutPoints) {
  if (
    is.logical(updateCutPoints) &&
      length(updateCutPoints) == 1L &&
      !is.na(updateCutPoints)
  ) {
    warnOnce(
      "tombstone.updateCutPoints.logical",
      "'updateCutPoints' is now one of \"none\", \"position\" or \"value\"; ",
      "TRUE was taken as \"position\" and FALSE as \"none\". A logical is no ",
      "longer taken in dbarts ",
      tombstoneExpiry,
      ".",
      class = "dbartsDeprecatedWarning"
    )
    return(if (updateCutPoints) "position" else "none")
  }
  matchCutPointRule(
    updateCutPoints,
    names(cutPointRuleCodes),
    "updateCutPoints"
  )
}

bartcoreSamplerSetPredictor <- function(
  sampler,
  x,
  column,
  forceUpdate,
  updateCutPoints
) {
  updateCutPoints <- resolveUpdateCutPoints(updateCutPoints)

  # read once: each sampler$data is a typed reference-class field's active
  # binding, a few microseconds a read, against a rejected update's whole
  # cost of under a hundred. Nothing below writes data@x before its last read.
  currentX <- sampler$data@x

  # A sparse design - a pure dgCMatrix or a mixed dense/sparse container -
  # accepts column-granular and whole-matrix mutation, maintained R-side by
  # installPredictorColumns rather than by the dense branch's pointer swap;
  # only per-observation replacement of a sparse-backed column stays fixed at
  # creation. Read before data@x is swapped.
  sparseSource <- predictorSourceIsSparse(currentX)

  # no BCF pre-check on the partial path either: the session's cell guard
  # caches every forest, pruned to the trees the column can move, so a row
  # installs only if it empties no leaf anywhere and a two-forest sampler
  # takes it
  partialUpdate <- !is.null(forceUpdate) &&
    is.character(forceUpdate) &&
    length(forceUpdate) == 1L &&
    !is.na(forceUpdate) &&
    forceUpdate == "partial"

  column <- resolveColumnIndex(currentX, column, "current X")
  if (!is.null(column)) {
    x <- codeCategoricalColumnUpdate(currentX, x, column)
  }

  # a triplet, row-compressed, symmetric, triangular, logical or pattern
  # sparse argument becomes the dgCMatrix the sparse path takes, rather than
  # densifying below
  x <- asDgCMatrix(x)

  ptr <- sampler$getPointer()

  if (partialUpdate) {
    if (is.null(column)) {
      stop("partial updates require a single 'column' to be specified")
    }
    if (length(column) != 1L) {
      stop("partial updates can only be applied to a single column")
    }
    column <- coerceOrError(column, "integer")
    # a CSC-backed column's rank storage cannot take a cell-at-a-time write
    # without an O(nnz) shift per cell; a DENSE-backed column of a mixed design
    # can, and is the motivating IRT latent case, so the refusal is per
    # column rather than per design
    if (predictorColumnIsSparseBacked(currentX, column)) {
      stop(
        "per-observation updates require a dense-backed column; replace a ",
        "sparse column wholesale with a non-partial update"
      )
    }
    if (updateCutPoints != "none") {
      stop("partial updates cannot also update cut points")
    }

    x <- as.double(x)
    installed <- .Call(
      C_dbarts_bartcore_updatePredictorPerObservation,
      ptr,
      x,
      as.integer(column)
    )
    # the engine keeps no predictor matrix, so maintain data@x R-side for the
    # observations the scan installed; install by reference - the merge
    # starts from the old column and overwrites only the installed rows
    sampler$data@x <- installPredictorColumns(
      currentX,
      installed,
      column,
      x[installed]
    )
    return(installed)
  }

  forceUpdate <- if (is.null(forceUpdate)) {
    is.null(column)
  } else {
    coerceOrError(forceUpdate, "logical")
  }
  if (length(forceUpdate) != 1L || is.na(forceUpdate)) {
    stop("'forceUpdate' must be TRUE, FALSE or \"partial\"")
  }
  updateCutPoints <- cutPointRuleCodes[[updateCutPoints]]

  # no BCF pre-check here either: a transactional whole-matrix or column
  # update revalidates every forest and rolls the whole change back if any
  # leaf of any tree of any forest would empty, so a two-forest sampler takes
  # it - the same as the per-observation session above.

  # dim(), not is.matrix(): the latter is FALSE for every Matrix class, so a
  # transposed dgCMatrix argument (same total length, wrong shape) fell
  # through to the length-only check below and was silently reinterpreted
  # column-major
  xDim <- dim(x)
  if (is.null(column)) {
    if (!is.null(xDim)) {
      if (xDim[2L] != ncol(currentX)) {
        stop("dimension of x must be equal to ", ncol(currentX))
      }
      if (xDim[1L] != nrow(currentX)) {
        stop("dimension of x must be equal to ", nrow(currentX))
      }
    } else if (length(x) != prod(dim(currentX))) {
      stop("'x' must have length ", prod(dim(currentX)))
    }
    # a sparse-valued argument onto a sparse-backed design rides to the bridge
    # as supplied: the bridge hands its sparse columns to the engine as stored
    # entries, under the store's own implicit rule. Every other argument - a
    # plain vector or a sparseVector - keeps the as.double path, as does a
    # plain-matrix design.
    if (!(sparseSource && predictorSourceIsSparse(x))) {
      # matrix(as.double(x), ...) strips every attribute, so the incoming
      # dimnames would otherwise vanish; carry them onto the replacement,
      # falling back to the sampler's current names when x supplies none -
      # the shapes already agree by the checks above
      xDimnames <- dimnames(x)
      x <- if (!is.null(xDim)) {
        matrix(as.double(x), xDim[1L])
      } else {
        matrix(as.double(x), nrow(currentX))
      }
      dimnames(x) <- if (!is.null(xDimnames)) {
        xDimnames
      } else {
        dimnames(currentX)
      }
      # the design's builder attributes describe its columns, which a
      # replacement of the rows does not change. A re-creation reads a
      # factor's declared levels from them, and counting levels from the
      # codes instead comes up short when the new rows miss the top one,
      # leaving the sampler's own state unbuildable on its copy or reload.
      for (name in hazardDesignAttrs) {
        if (is.null(attr(x, name))) {
          attr(x, name) <- attr(currentX, name)
        }
      }
    }
    if (!sparseSource) {
      # a pointer swap: the engine borrows data@x, so install there first and
      # revert if the transaction rolls back
      oldX <- currentX
      sampler$data@x <- x
      tryResult <- tryCatch(
        updateSuccessful <- .Call(
          C_dbarts_bartcore_setPredictor,
          ptr,
          sampler$data@x,
          forceUpdate,
          updateCutPoints
        ),
        error = function(e) {
          sampler$data@x <- oldX
          e
        }
      )
      if (inherits(tryResult, "error")) {
        stop(tryResult)
      }
      if (!forceUpdate && !updateSuccessful) sampler$data@x <- oldX
    } else {
      # a sparse-bearing source: the engine borrows the argument matrix rather
      # than data@x, so nothing needs installing until it accepts. Splice the
      # replacement columns into the container BEFORE the call, so a throw
      # there cannot leave data@x describing the old design (sampler
      # re-creation after save/load reads it). A replaced sparse column stays
      # sparse: the engine and the splice both keep only the entries that
      # differ from its implicit value.
      newX <- installPredictorColumns(
        currentX,
        NULL,
        seq_len(ncol(currentX)),
        x
      )
      updateSuccessful <- .Call(
        C_dbarts_bartcore_setPredictor,
        ptr,
        x,
        forceUpdate,
        updateCutPoints
      )
      if (isTRUE(updateSuccessful)) sampler$data@x <- newX
    }
  } else {
    column <- coerceOrError(column, "integer")
    if (any(column < 1L | column > ncol(currentX))) {
      stop(
        "column '",
        column[which(column < 1L | column > ncol(currentX))[1L]],
        "' is out of range"
      )
    }
    # length() counts a container's fields rather than its cells, so the shape
    # check reads dim() wherever the argument carries one; only a dimensionless
    # argument falls back to the total-length check
    if (!is.null(xDim)) {
      if (xDim[2L] != length(column)) {
        stop("'x' must have ", length(column), " column(s)")
      }
      if (xDim[1L] != nrow(currentX)) {
        stop("'x' must have ", nrow(currentX), " row(s)")
      }
    } else if (length(x) != nrow(currentX) * length(column)) {
      stop("'x' must have length ", nrow(currentX) * length(column))
    }
    if (!(sparseSource && predictorSourceIsSparse(x))) {
      x <- as.double(x)
    }
    # the engine keeps no predictor matrix, so maintain data@x R-side when the
    # update is applied (forceUpdate, or a non-rolled-back transaction);
    # install by reference - only the addressed columns change, the rest of the
    # container is shared. Build the replacement BEFORE the engine commits: a
    # CSC-backed column's install rewrites the container's sparse slots, and a
    # throw there after acceptance would leave data@x describing the old design
    # (sampler re-creation after save/load reads it).
    newX <- installPredictorColumns(currentX, NULL, column, x)
    updateSuccessful <- .Call(
      C_dbarts_bartcore_updatePredictor,
      ptr,
      x,
      column,
      forceUpdate,
      updateCutPoints
    )
    if (isTRUE(updateSuccessful)) {
      sampler$data@x <- newX
    }
  }

  # a forced update always installs, and what it returns says nothing of
  # validity (dec-B310): NULL, invisibly, as 0.9-34's did
  if (!forceUpdate) updateSuccessful else invisible(NULL)
}

# The response conduits' updateScale: a single TRUE or FALSE. NA or 1 would
# otherwise skip the isTRUE pre-checks below and reach the engine as a
# different answer than the R side acted on.
checkUpdateScale <- function(updateScale) {
  if (
    !is.logical(updateScale) ||
      length(updateScale) != 1L ||
      is.na(updateScale)
  ) {
    stop("'updateScale' must be TRUE or FALSE", call. = FALSE)
  }
  updateScale
}

bartcoreSamplerSetResponse <- function(
  sampler,
  y,
  updateScale = FALSE,
  status = NULL
) {
  updateScale <- checkUpdateScale(updateScale)
  y <- as.double(y)
  if (anyNA(y)) {
    stop("response contains missing values")
  }
  # as creation refuses one: a single infinite value leaves sigma and every
  # fit NaN, even after the response is put back
  if (any(is.infinite(y))) {
    stop("response contains non-finite values")
  }
  if (!is.null(status)) {
    status <- as.double(status)
    if (anyNA(status)) {
      stop("survival status cannot be NA")
    }
  }
  if (isTRUE(updateScale)) {
    refuseAmplitudeMutation(
      sampler,
      "setResponse(updateScale = TRUE)",
      "every forest keeps its leaf calibration stated against the scale ",
      "fixed at creation; use updateScale = FALSE instead"
    )
  }
  # validate (the C length, family and support checks) before installing, so a
  # rejected y or status never leaves data@y or the control attribute holding
  # the bad replacement
  .Call(
    C_dbarts_bartcore_setResponse,
    sampler$getPointer(),
    y,
    updateScale,
    status
  )
  sampler$data@y <- y
  if (!is.null(status)) {
    # re-creation after save and load rebuilds the sampler from the control,
    # model and data it holds, so the status it reads must be the current one -
    # data@y's rule, on the channel the status travels
    newControl <- sampler$control
    attr(newControl, "bartcore.survival") <- status
    sampler$control <- newControl
  }
  invisible(NULL)
}

bartcoreSamplerSetOffset <- function(sampler, offset, updateScale) {
  updateScale <- checkUpdateScale(updateScale)
  if (updateScale) {
    refuseAmplitudeMutation(
      sampler,
      "setOffset(updateScale = TRUE)",
      "every forest keeps its leaf calibration stated against the scale ",
      "fixed at creation; use updateScale = FALSE instead"
    )
  }
  # a synced test offset follows the regular one;
  # NA marks "leave the test offset alone"
  offset.test <- NA
  if (is.null(offset)) {
    if (identical(sampler$data@testUsesRegularOffset, TRUE)) {
      offset.test <- NULL
    }
  } else {
    offset <- as.double(offset)
    if (anyNA(offset)) {
      stop("'offset' contains missing values")
    }
    if (any(is.infinite(offset))) {
      stop("'offset' contains non-finite values")
    }
    if (length(offset) == 1L) {
      if (identical(sampler$data@testUsesRegularOffset, TRUE)) {
        offset.test <- if (!is.null(sampler$data@x.test)) {
          rep_len(offset, nrow(sampler$data@x.test))
        } else {
          NULL
        }
      }
      offset <- rep_len(offset, length(sampler$data@y))
    } else {
      if (length(offset) != length(sampler$data@y)) {
        stop(
          "length of replacement offset is not equal to number of observations"
        )
      }
      if (identical(sampler$data@testUsesRegularOffset, TRUE)) {
        offset.test <- if (
          !is.null(sampler$data@x.test) &&
            length(offset) == nrow(sampler$data@x.test)
        ) {
          offset
        } else {
          NULL
        }
      }
    }
  }

  ptr <- sampler$getPointer()

  # the engine copies the offset, so the mirror is installed only once the
  # engine has taken it: a refused swap leaves data@offset what the live
  # sampler holds, which a save and load re-creates from
  .Call(C_dbarts_bartcore_setOffset, ptr, offset, updateScale)
  sampler$data@offset <- offset

  if (!identical(offset.test, NA)) {
    oldOffset.test <- sampler$data@offset.test
    sampler$data@offset.test <- offset.test
    tryResult <- tryCatch(
      .Call(C_dbarts_bartcore_setTestOffset, ptr, sampler$data@offset.test),
      error = function(e) {
        sampler$data@offset.test <- oldOffset.test
        e
      }
    )
    if (inherits(tryResult, "error")) stop(tryResult)
  }

  invisible(NULL)
}

bartcoreSamplerSetData <- function(sampler, newData) {
  if (!inherits(newData, "dbartsData")) {
    stop("'data' must inherit from dbartsData")
  }
  refuseAmplitudeMutation(
    sampler,
    "setData",
    "every forest is calibrated against the data at creation; make a new ",
    "sampler instead"
  )
  if (ncol(newData@x) != ncol(sampler$data@x)) {
    stop("$setData: requires the same predictors")
  }

  newData@n.cuts <- sampler$data@n.cuts
  newData@sigma <- sampler$data@sigma

  # a probit or ordinal sampler carries no weight channel for the whole-data
  # conduit to fill, and 0/1 weights there name the rows in the data set: pull
  # them off the object and install them as the active-row mask instead, sized
  # by the replacement's own n. The mask goes in AFTER the swap, which clears
  # whatever was in force (n may change with the data); all-ones resolves to
  # no mask at all, as it does at creation.
  active <- NULL
  if (isMaskedWeightFamily(sampler$model@family) && !is.null(newData@weights)) {
    w <- newData@weights
    if (anyNA(w) || any(w != 0 & w != 1)) {
      stop(
        sampler$model@family,
        " models do not support case weights other than 0 and 1, which mark ",
        "rows in and out of the likelihood: such a vector installs as the ",
        "active-row mask, and a weighted truncated-normal latent likelihood ",
        "is not a coherent model"
      )
    }
    if (!all(w == 1)) {
      active <- w
    }
    newData@weights <- NULL
  }

  ptr <- sampler$getPointer()

  oldData <- sampler$data
  sampler$data <- newData
  tryResult <- tryCatch(
    .Call(C_dbarts_bartcore_setData, ptr, sampler$data),
    error = function(e) {
      sampler$data <- oldData
      e
    }
  )
  if (inherits(tryResult, "error")) {
    stop(tryResult)
  }
  # the swap cleared whatever mask was in force, so the mirror that would
  # otherwise re-apply one at re-creation goes with it; setActiveRows records
  # the replacement's own below
  sampler$activeRows <- NULL
  # the caller's own updateState decides the store, once, after the swap
  if (!is.null(active)) {
    sampler$setActiveRows(active, updateState = FALSE)
  }

  invisible(NULL)
}

bartcoreSamplerSetCutPoints <- function(sampler, cuts, column, splits) {
  splits <- matchCutPointRule(splits, c("position", "value"), "splits")
  # a missing column stays NULL: the bridge then takes one entry per column
  # and skips those of factor columns
  column <- resolveColumnIndex(sampler$data@x, column, "current X")
  if (!is.null(column)) {
    column <- coerceOrError(column, "integer")
  }

  # a data frame is a list of columns and is taken as one
  cuts <- if (is.list(cuts)) as.list(cuts) else list(cuts)
  # a list of another length is the bridge's to refuse, and goes to it with
  # no entry read. Of one of the right length, the entries the bridge skips
  # are dropped unread, so nothing in a factor column's place warns or fails,
  # and every other entry is a vector of numbers: as.double would take a
  # factor's codes, a Date or a logical as a grid
  varTypes <- sampler$data@varTypes
  numEntries <- if (is.null(column)) length(varTypes) else length(column)
  if (length(cuts) == numEntries) {
    isRead <- rep_len(TRUE, numEntries)
    if (is.null(column)) {
      isRead <- varTypes == ORDINAL_VARIABLE
    }
    if (!all(vapply(cuts[isRead], is.numeric, NA))) {
      stop("$setCutPoints: 'cuts' must be numeric", call. = FALSE)
    }
    cuts[!isRead] <- list(NULL)
    # a grid is a set of thresholds, so one out of order is sorted; one with a
    # NaN is left for the bridge to refuse, which sort() would drop it from
    cuts[isRead] <- lapply(cuts[isRead], function(grid) {
      grid <- as.double(grid)
      if (!anyNA(grid) && is.unsorted(grid)) sort(grid) else grid
    })
  }

  # the engine re-quantizes the transient borrow of the current predictors
  .Call(
    C_dbarts_bartcore_setCutPoints,
    sampler$getPointer(),
    cuts,
    column,
    rawPredictorMatrix(sampler$data@x),
    cutPointRuleCodes[[splits]]
  )
  invisible(NULL)
}

bartcoreSamplerSetTestPredictor <- function(sampler, x.test, column) {
  column <- resolveColumnIndex(
    sampler$data@x.test,
    column,
    "current test predictor matrix"
  )

  # a column update keeps the rows, and so their names
  testRowNames <- dataRowNames(sampler$data, "test")
  if (is.null(column)) {
    # NULL removes the test data; a frame/sparse input becomes a container the
    # bridge codes against the training cuts. The bridge clears any test offset
    # with a NULL removal.
    testRowNames <- observationRowNames(x.test)
    x.test <- validateXTest(x.test, sampler$data@x)
  } else {
    column <- coerceOrError(column, "integer")
    if (any(column < 1L | column > ncol(sampler$data@x.test))) {
      stop(
        "column '",
        column[which(column < 1L | column > ncol(sampler$data@x.test))[1L]],
        "' is out of range"
      )
    }
    x.test <- codeCategoricalColumnUpdate(sampler$data@x, x.test, column)
    xTestDim <- dim(x.test)
    if (!is.null(xTestDim) && xTestDim[2L] != length(column)) {
      stop("'x.test' must have ", length(column), " column(s)")
    }
    if (length(x.test) != nrow(sampler$data@x.test) * length(column)) {
      stop(
        "'x.test' must have length ",
        nrow(sampler$data@x.test) * length(column)
      )
    }
    if (inherits(sampler$data@x.test, "dbartsMixedMatrix")) {
      # a container's per-column storage decision (dense vs CSC-backed) is
      # preserved; installPredictorColumns splices the replacement in place,
      # canonicalizing a CSC-backed target column against its implicit
      x.test <- installPredictorColumns(
        sampler$data@x.test,
        NULL,
        column,
        x.test
      )
    } else {
      # the engine replaces the whole matrix; column updates copy-modify it
      new.x.test <- sampler$data@x.test
      new.x.test[, column] <- as.double(x.test)
      x.test <- new.x.test
    }
  }

  # install the new test set R-side, then roll back if the bridge refuses it
  # (a container whose leaf-covariate column is CSC-backed), keeping the
  # R-level object and the engine's prior test store consistent
  oldX.test <- sampler$data@x.test
  oldOffset.test <- sampler$data@offset.test
  sampler$data@x.test <- x.test
  if (is.null(x.test)) {
    sampler$data@offset.test <- NULL
  }
  tryResult <- tryCatch(
    .Call(
      C_dbarts_bartcore_setTestPredictor,
      sampler$getPointer(),
      sampler$data@x.test
    ),
    error = function(e) {
      sampler$data@x.test <- oldX.test
      sampler$data@offset.test <- oldOffset.test
      e
    }
  )
  if (inherits(tryResult, "error")) {
    stop(tryResult)
  }
  sampler$data <- setDataRowNames(sampler$data, "test", testRowNames)
  invisible(NULL)
}

# A built predictor store (cuts + codes) shared across row-subset samplers;
# internal and unserializable. control contributes useQuantiles; data
# contributes x, the column types, and
# n.cuts. leafCovariateColumns names (1-based) the columns whose raw values
# a view's leaf model will read; the handle owns raw only for those, so a
# constant-leaf caller passes none. A view designating an undeclared column
# is refused when it is built.
bartcoreDataHandle <- function(control, data, leafCovariateColumns = NULL) {
  result <- new.env(parent = emptyenv())
  result$ptr <- .Call(
    C_dbarts_bartcore_createDataHandle,
    control,
    data,
    if (!is.null(leafCovariateColumns)) as.integer(leafCovariateColumns)
  )
  result
}

# A sampler over a row subset of a handle: it copies the handle's cut grid
# and gathers its rows' codes, so folds bin identically to the full data.
# data is the full data object the handle was built from; trainRows and
# testRows index its rows, y/weights/offset are sliced by trainRows, and a
# test offset comes from offset[testRows] (xbart's fold semantics). columns,
# when given, are 1-based indices restricting the view to a column subset (NULL
# spans every column). The result refuses raw-predictor mutation (setPredictor
# and friends, setData, setCutPoints, setState); family overrides the
# response model ("" keeps the bridge's own dispatch, "logistic" selects it
# for a binary response).
bartcoreSamplerFromHandle <- function(
  handle,
  control,
  model,
  data,
  trainRows,
  testRows = NULL,
  family = "",
  columns = NULL
) {
  result <- new.env(parent = emptyenv())
  result$ptr <- .Call(
    C_dbarts_bartcore_createFromHandle,
    control,
    model,
    data,
    handle$ptr,
    as.integer(trainRows),
    if (!is.null(testRows)) as.integer(testRows),
    as.character(family),
    if (!is.null(columns)) as.integer(columns)
  )
  # a view refuses the raw-predictor and re-quantize surface, and its live-tree
  # getTrees needs no training replay, so it carries no predictor matrix
  result$x <- NULL
  result
}

# Integer coercion for an already-validated count matrix. The RANGE check comes
# BEFORE the coercion: storage.mode() turns a value past .Machine$integer.max
# into NA with a warning, which the engine would then report as a negative
# count - a true refusal naming the wrong reason. A missing row (na.action's
# na.pass) has nothing to range-check, so it is ignored here rather than
# tripping an ambiguous "missing value where TRUE/FALSE needed".
asCountMatrix <- function(counts) {
  if (any(counts > .Machine$integer.max, na.rm = TRUE)) {
    stop("multinomial counts must be representable as integers")
  }
  storage.mode(counts) <- "integer"
  counts
}

# The n x K category offset a multinomial predict() validates its own offset
# argument against (generics.R, predict.bartMultinomial): NULL means none,
# anything else must be a numeric n x K matrix with every entry finite (an
# infinite entry propagates through the log-sum-exp margin into a NaN for
# every category of that row). Only the row-centred part of the offset is
# identified - adding a constant to every entry of a row leaves the softmax
# unchanged - so the input is passed through as given rather than silently
# re-centred. `what` names the offset in the refusals, since the train rows
# and the test rows each carry their own.
validateCategoryOffset <- function(offset, n, K, what = "category offset") {
  if (is.null(offset)) {
    return(NULL)
  }
  offset <- as.matrix(offset)
  if (!is.numeric(offset)) {
    stop(sprintf("multinomial %s must be a numeric matrix", what))
  }
  if (nrow(offset) != n || ncol(offset) != K) {
    stop(sprintf(
      "multinomial %s must be a %d x %d matrix",
      what,
      n,
      K
    ))
  }
  if (!all(is.finite(offset))) {
    stop(sprintf("multinomial %s must be finite", what))
  }
  storage.mode(offset) <- "double"
  offset
}

# A string as the code it spells, so that "I(dose / 30)" finds the label
# "I(dose/30)"; NA where it is not one expression of R.
forestLabelCode <- function(text) {
  parsed <- tryCatch(
    parse(text = text, keep.source = FALSE),
    error = function(e) NULL
  )
  if (length(parsed) != 1L) {
    return(NA_character_)
  }
  paste(deparse(parsed[[1L]], width.cutoff = 500L), collapse = " ")
}

# The refusal of a string that is the label of forest `labelAt` and also the
# name forest<i> of position `position`.
forestLabelConflict <- function(text, labelAt, position) {
  paste0(
    "'forest' (",
    encodeString(text, quote = "\""),
    ") is the label of forest ",
    labelAt,
    " and the name of position ",
    position,
    "; select by position, as forest = ",
    labelAt
  )
}

# The refusal of a name in a list given forest by forest that is another
# forest's label where it is also this position's forest<i>; NULL otherwise.
# `labels` is every forest's label, "" or NULL where none is recorded.
# `argument` is the list's own name, which the remedy speaks of: the list is
# read by position and has no forest argument to select by.
forestNameTaken <- function(name, position, labels, argument) {
  other <- setdiff(which(labels == name), position)
  if (length(other) == 0L) {
    return(NULL)
  }
  paste0(
    "'",
    argument,
    "' names forest ",
    position,
    " ",
    encodeString(name, quote = "\""),
    ", which is the label of forest ",
    other[[1L]],
    "; name each entry by its own forest's label, or leave the names off"
  )
}

# The one reader of a 'forest' argument that is not NULL, given the labels
# the model's forests carry (NULL where it records none: a sampler of one
# forest, a multinomial one) and the number of forests. A number is a
# position, as `[[` reads one, and is returned untouched for the caller's own
# checks. A string is never a position: it is the forest with that label;
# failing that, the forest whose label is the same code; and "forest<i>" is
# forest i, the name position i has on every per-forest margin. A string that
# is one forest's label and another position's name, or is the same code as
# two labels, is refused naming both. NA, an empty string, a factor, a
# logical and a list are refused by name. With `several` a vector of strings
# returns a position for each, in the order given.
selectForest <- function(forest, labels, numForests, several = FALSE) {
  if (is.numeric(forest)) {
    return(forest)
  }
  if (is.logical(forest) && length(forest) == 1L && is.na(forest)) {
    stop("'forest' must not be NA or an empty string")
  }
  if (!is.character(forest)) {
    stop(
      "'forest' must be a number, the forest's position, or a string, its ",
      "label; not ",
      if (is.factor(forest)) {
        "a factor"
      } else if (is.logical(forest)) {
        "a logical"
      } else if (is.list(forest)) {
        "a list"
      } else {
        paste0("a ", class(forest)[[1L]])
      }
    )
  }
  if (anyNA(forest) || any(!nzchar(forest))) {
    stop("'forest' must not be NA or an empty string")
  }
  if (!several && length(forest) != 1L) {
    stop("'forest' must be a single number or a single label")
  }
  # labels and the forest count come from one record, so they agree; this
  # fires only if a record were written with the two out of step (a corrupt
  # restore), when no label can be trusted and forest<i> alone is taken
  if (!is.character(labels) || length(labels) != numForests) {
    labels <- NULL
  }
  codes <- NULL
  quote <- function(text) encodeString(text, quote = "\"")
  listForests <- function(index) {
    if (length(index) == 2L) {
      paste(index, collapse = " and ")
    } else {
      paste0(
        paste(index[-length(index)], collapse = ", "),
        " and ",
        index[length(index)]
      )
    }
  }
  one <- function(text) {
    named <- suppressWarnings(
      if (grepl("^forest[0-9]+$", text)) {
        as.integer(substring(text, 7L))
      } else {
        NA_integer_
      }
    )
    if (!is.na(named) && (named < 1L || named > numForests)) {
      named <- NA_integer_
    }
    exact <- if (!is.null(labels)) which(labels == text)
    if (length(exact) == 1L) {
      if (!is.na(named) && named != exact) {
        stop(forestLabelConflict(text, exact, named))
      }
      return(exact)
    }
    if (!is.na(named)) {
      return(named)
    }
    if (!is.null(labels)) {
      if (is.null(codes)) {
        codes <<- vapply(labels, forestLabelCode, "", USE.NAMES = FALSE)
      }
      target <- forestLabelCode(text)
      hits <- if (!is.na(target)) which(codes == target)
      if (length(hits) == 1L) {
        return(hits)
      }
      if (length(hits) > 1L) {
        stop(
          "'forest' (",
          quote(text),
          ") is the label of forests ",
          listForests(hits),
          " (",
          paste(quote(labels[hits]), collapse = ", "),
          "); give one exactly, or select by position"
        )
      }
    }
    stop(
      "'forest' names no forest of this model: ",
      quote(text),
      if (is.null(labels)) {
        "; this sampler's forests have no labels, so select one by position"
      } else {
        paste0(
          "; its forests are ",
          paste(quote(labels), collapse = ", "),
          if (grepl("^[0-9]+$", text)) {
            paste0("; a position is given as a number, forest = ", text)
          }
        )
      }
    )
  }
  vapply(forest, one, 0L, USE.NAMES = FALSE)
}

# The labels a sampler's forests carry, from the forests' record on its
# control; NULL for a sampler that records none.
samplerForestLabels <- function(control) {
  attr(control, "bartcore.forests", exact = TRUE)$labels
}

# A sampler's reading of its 'forest' argument: its own labels, and its forest
# count from the engine, asked only when a string needs it. `ptr` is not
# evaluated for a number.
samplerForestIndex <- function(forest, control, ptr, several = FALSE) {
  labels <- samplerForestLabels(control)
  numForests <- if (is.character(forest)) {
    bartcoreNumForests(ptr)
  } else {
    NA_integer_
  }
  if (several) {
    return(selectForest(forest, labels, numForests, TRUE) - 1L)
  }
  resolveForestIndex(forest, labels, numForests)
}

# The R-level calibration surface indexes forests from 1, as R indexes; the
# bridge and the flat C entries count from 0, as the engine does. A string is
# a label, read by selectForest; `labels` and `numForests` are the sampler's
# and are read only for one.
resolveForestIndex <- function(
  forest,
  labels = NULL,
  numForests = NA_integer_
) {
  forest <- selectForest(forest, labels, numForests)
  forest <- coerceOrError(forest, "integer")
  if (length(forest) != 1L || is.na(forest) || forest < 1L) {
    stop("'forest' must be a single positive integer (1 selects the first)")
  }
  forest - 1L
}

# The sampler's forest count. A COUNT, not a capability probe:
# samplerCarriesAmplitudes and samplerCarriesCounts each answer only for their
# own model, and neither sees a plain single-forest sampler.
bartcoreNumForests <- function(ptr) .Call(C_dbarts_bartcore_numForests, ptr)

# One 0-based forest's own tree count, the engine's authoritative answer
# (getTrees defaults and validates treeNums against it, forest by forest):
# forest 1's count need not equal control@n.trees once a forest() term or a
# forests = entry gives a later forest its own n.trees.
bartcoreForestTreeCount <- function(ptr, forest) {
  .Call(C_dbarts_bartcore_numTreesInForest, ptr, forest)
}

# The forest-count refusal on a DONOR warm start. At more than one forest the
# install would answer rather than raise - the trees arrive from a saved slot
# and the amplitudes from the donor's live state - leaving a legal-looking fit
# with a miscalibrated start that nothing covers. Grow-from-root is NOT held
# to this: it composes through the combiner every sweep and is covered at two
# forests. A COUNT, not a capability probe: the amplitude coupling and the
# K-forest softmax are both uncovered here, and the count is what the message
# can name. Raised R-side so a caller reads the surface it wrote rather than
# the bridge's own entry point name; the bridge keeps the same refusal as a
# backstop for the routes that skip this layer.
refuseMultiForestWarmStart <- function(ptr, what) {
  numForests <- bartcoreNumForests(ptr)
  if (numForests >= 2L) {
    stop(
      what,
      " does not support a multi-forest sampler: this one carries ",
      numForests,
      " forests, which have no tested warm start from a donor; draw from the ",
      "prior, or grow from the root, instead"
    )
  }
  invisible(NULL)
}

# The control travels with the model: the tree-move mixture is a control slot
# the prior install reads.
bartcoreSetModel <- function(bcSampler, model, data, control) {
  invisible(.Call(
    C_dbarts_bartcore_setModel,
    bcSampler$ptr,
    model,
    data,
    control
  ))
}

# Drives a low-level bartcore handle (a bcSampler env holding $ptr) directly,
# not a dbartsSampler; cf. bartcoreSamplerRun, the R-level sampler-layer entry.
bartcoreRun <- function(bcSampler, numBurnIn = 0L, numSamples = 1L) {
  .Call(
    C_dbarts_bartcore_run,
    bcSampler$ptr,
    as.integer(numBurnIn),
    as.integer(numSamples),
    NULL,
    NULL,
    TRUE
  )
}
