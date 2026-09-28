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
  invisible(ptr)
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
  if (is.na(numBurnIn)) {
    numBurnIn <- control@n.burn
  }
  if (is.na(numSamples)) {
    numSamples <- control@n.samples
  }
  if (is.na(numSamples)) {
    stop("bartcore engine samplers require 'numSamples' to be specified")
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
  if (is.null(result)) {
    return(invisible(NULL))
  }
  warnOnGPFallback(result)
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
        "observations, so most of this fit is not a Gaussian process; ",
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

# Resolves a character 'column' against source's colnames into a 1-based
# integer index (or indices); NULL or an already-numeric 'column' passes
# through unchanged. 'what' names source for the not-found message.
resolveColumnIndex <- function(source, column, what) {
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

bartcoreSamplerSetPredictor <- function(
  sampler,
  x,
  column,
  forceUpdate,
  updateCutPoints
) {
  # A sparse design - a pure dgCMatrix or a mixed dense/sparse container -
  # accepts column-granular and whole-matrix mutation, maintained R-side by
  # installPredictorColumns rather than by the dense branch's pointer swap;
  # only per-observation replacement of a sparse-backed column stays fixed at
  # creation. Read before data@x is swapped.
  sparseSource <- predictorSourceIsSparse(sampler$data@x)

  # no BCF pre-check on the partial path either: the session's cell guard
  # caches every forest, pruned to the trees the column can move, so a row
  # installs only if it empties no leaf anywhere and a two-forest sampler
  # takes it
  partialUpdate <- !is.null(forceUpdate) &&
    is.character(forceUpdate) &&
    length(forceUpdate) == 1L &&
    !is.na(forceUpdate) &&
    forceUpdate == "partial"

  column <- resolveColumnIndex(sampler$data@x, column, "current X")

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
    if (predictorColumnIsSparseBacked(sampler$data@x, column)) {
      stop(
        "per-observation updates require a dense-backed column; replace a ",
        "sparse column wholesale with a non-partial update"
      )
    }
    if (isTRUE(coerceOrError(updateCutPoints, "logical"))) {
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
      sampler$data@x,
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
  updateCutPoints <- coerceOrError(updateCutPoints, "logical")

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
      if (xDim[2L] != ncol(sampler$data@x)) {
        stop("dimension of x must be equal to ", ncol(sampler$data@x))
      }
      if (xDim[1L] != nrow(sampler$data@x)) {
        stop("dimension of x must be equal to ", nrow(sampler$data@x))
      }
    } else if (length(x) != prod(dim(sampler$data@x))) {
      stop("'x' must have length ", prod(dim(sampler$data@x)))
    }
    # a sparse-valued argument onto a sparse-backed design rides to the bridge
    # as supplied: it materializes there, under the store's own implicit rule,
    # rather than being densified here. Every other argument - a plain vector,
    # a sparseVector, any Matrix class the bridge does not ingest - keeps the
    # as.double path, as does a plain-matrix design.
    if (!(sparseSource && predictorSourceIsSparse(x))) {
      # matrix(as.double(x), ...) strips every attribute, so the incoming
      # dimnames would otherwise vanish; carry them onto the replacement,
      # falling back to the sampler's current names when x supplies none -
      # the shapes already agree by the checks above
      xDimnames <- dimnames(x)
      x <- if (!is.null(xDim)) {
        matrix(as.double(x), xDim[1L])
      } else {
        matrix(as.double(x), nrow(sampler$data@x))
      }
      dimnames(x) <- if (!is.null(xDimnames)) {
        xDimnames
      } else {
        dimnames(sampler$data@x)
      }
    }
    if (!sparseSource) {
      # a pointer swap: the engine borrows data@x, so install there first and
      # revert if the transaction rolls back
      oldX <- sampler$data@x
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
      # re-creation after save/load reads it). Replacing a sparse column
      # densifies its storage - every row now differs from the implicit value.
      newX <- installPredictorColumns(
        sampler$data@x,
        NULL,
        seq_len(ncol(sampler$data@x)),
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
    if (any(column < 1L | column > ncol(sampler$data@x))) {
      stop(
        "column '",
        column[which(column < 1L | column > ncol(sampler$data@x))[1L]],
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
      if (xDim[1L] != nrow(sampler$data@x)) {
        stop("'x' must have ", nrow(sampler$data@x), " row(s)")
      }
    } else if (length(x) != nrow(sampler$data@x) * length(column)) {
      stop("'x' must have length ", nrow(sampler$data@x) * length(column))
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
    newX <- installPredictorColumns(sampler$data@x, NULL, column, x)
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

  if (!forceUpdate) updateSuccessful else invisible(NULL)
}

bartcoreSamplerSetResponse <- function(
  sampler,
  y,
  updateScale = FALSE,
  status = NULL
) {
  y <- as.double(y)
  if (anyNA(y)) {
    stop("response contains missing values")
  }
  if (!is.null(status)) {
    status <- as.double(status)
  }
  if (isTRUE(updateScale)) {
    refuseAmplitudeMutation(
      sampler,
      "setResponse(updateScale = TRUE)",
      "every forest keeps its leaf calibration stated against the anchor ",
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
  if (isTRUE(updateScale)) {
    refuseAmplitudeMutation(
      sampler,
      "setOffset(updateScale = TRUE)",
      "every forest keeps its leaf calibration stated against the anchor ",
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

  sampler$data@offset <- offset
  .Call(
    C_dbarts_bartcore_setOffset,
    ptr,
    sampler$data@offset,
    as.logical(updateScale)
  )

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
    stop("bartcore setData requires the same predictors")
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

bartcoreSamplerSetCutPoints <- function(sampler, cuts, column) {
  column <- resolveColumnIndex(sampler$data@x, column, "current X")
  if (is.null(column)) {
    column <- seq_len(ncol(sampler$data@x))
  }

  if (!is.list(cuts)) {
    cuts <- list(cuts)
  }
  cuts <- lapply(cuts, as.double)

  # the engine re-quantizes the transient borrow of the current predictors
  column <- coerceOrError(column, "integer")
  .Call(
    C_dbarts_bartcore_setCutPoints,
    sampler$getPointer(),
    cuts,
    column,
    rawPredictorMatrix(sampler$data@x)
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

# The R-level calibration surface indexes forests from 1, as R indexes; the
# bridge and the flat C entries count from 0, as the engine does.
resolveForestIndex <- function(forest) {
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
