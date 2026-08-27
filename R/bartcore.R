# The bartcore engine behind dbartsSampler. The dbartsSampler methods
# delegate to the bartcoreSampler* functions below; the C side borrows
# vectors and pins them in the external pointer's protection slot.
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

# Drives a dbartsSampler (the R-level sampler layer), reading its control
# defaults and delegating through its external pointer; cf. bartcoreRun,
# which drives a low-level bartcore handle directly.
bartcoreSamplerRun <- function(sampler, numBurnIn, numSamples) {
  control <- sampler$control
  if (is.na(numBurnIn)) {
    numBurnIn <- control@n.burn
  }
  if (is.na(numSamples)) {
    numSamples <- control@n.samples
  }
  if (is.na(numSamples)) {
    stop("bartcore engine samplers require 'numSamples' to be specified")
  }

  result <- .Call(
    C_dbarts_bartcore_run,
    sampler$getPointer(),
    as.integer(numBurnIn),
    as.integer(numSamples)
  )
  if (is.null(result)) {
    return(invisible(NULL))
  }
  result
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
    column <- as.integer(column)
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
    column <- as.integer(column)
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

bartcoreSamplerSetResponse <- function(sampler, y, updateScale = FALSE) {
  y <- as.double(y)
  if (anyNA(y)) {
    stop("response contains missing values")
  }
  if (isTRUE(updateScale)) {
    refuseAmplitudeMutation(
      sampler,
      "setResponse(updateScale = TRUE)",
      "every forest keeps its leaf calibration stated against the anchor ",
      "fixed at creation; use updateScale = FALSE instead"
    )
  }
  # validate (the C length check) before installing, so a rejected y never
  # leaves data@y holding the bad replacement
  .Call(
    C_dbarts_bartcore_setResponse,
    sampler$getPointer(),
    y,
    updateScale
  )
  sampler$data@y <- y
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
  .Call(
    C_dbarts_bartcore_setCutPoints,
    sampler$getPointer(),
    cuts,
    as.integer(column),
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

  if (is.null(column)) {
    # NULL removes the test data; a frame/sparse input becomes a container the
    # bridge codes against the training cuts. The bridge clears any test offset
    # with a NULL removal.
    x.test <- validateXTest(x.test, sampler$data@x)
  } else {
    column <- as.integer(column)
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
  invisible(NULL)
}

# Internal interface used by the tests and the equivalence harness; a
# bartcore handle is constructed from the validated control/model/data of an
# existing dbartsSampler regardless of that sampler's own engine. family
# selects the response model for binary responses ("probit", the default,
# or "logistic"), the same choice bart2's own family argument exposes
# publicly. Columns marked CATEGORICAL (1L) in data@varTypes split by
# category subset rather than by threshold; factors = "categorical", the
# default on the public data constructor, sets this directly, so callers
# never flip the type by hand. Values must be integer codes 0..K-1, K <= 32.

bartcoreSampler <- function(sampler, family = "") {
  result <- new.env(parent = emptyenv())
  result$ptr <- .Call(
    C_dbarts_bartcore_create,
    sampler$control,
    sampler$model,
    sampler$data,
    as.character(family)
  )
  # the engine keeps no predictor matrix, so the low-level wrappers track the
  # current predictors R-side to feed the re-quantize entry points
  # (setCutPoints, setState, saved-tree getTrees); null for sparse designs
  result$x <- rawPredictorMatrix(sampler$data@x)
  result
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
# and friends, setData, setCutPoints, setState); family is as bartcoreSampler's.
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

# A BCF two-forest sampler, internal:
# the model spec is the prognostic forest mu(x, pihat) - the caller supplies
# pihat as an ordinary predictor column - the arguments below the treatment
# forest tau(x) and the glue. z is the 0/1 treatment. sd.control and
# sd.moderate are the prognostic and effect magnitudes in units of the
# family's own latent scale (sd(y) under gaussian, 1 under probit); the
# calibration map converts them to the per-forest leaf scales at creation,
# overriding the host model's node prior and k for the mu forest. A NULL
# sd.control takes the family's own default median, 2 under gaussian and 1
# under the latent families, exactly as the public route does; sd.moderate's 1
# is the K-aware default's own K = 2 value, sqrt(2/2), so this two-forest
# spelling states it as the literal it has always been. update.a /
# update.b hold the matching glue block fixed when FALSE. moderators restricts
# the treatment forest to a subset of the shared design's columns (a character
# vector matched to colnames(data@x), or a 1-based numeric index vector; NULL =
# all columns); the prognostic forest always reads the full store. In the BCF
# of Hahn, Murray and Carvalho (2020) the estimated propensity score enters
# mu's design but not tau's - with a shared design that means carrying pihat as
# a data column and leaving it out of moderators.
#
# family names the response family the two forests are combined under
# (gaussian, probit or logistic), and is written into the model this call
# hands the bridge, which is the ONE place the family is read from. Supplying
# it here rather than only on the model keeps the family visible at the call
# site, as the other internal creators do, while leaving one source of truth.
bartcoreBCFSampler <- function(
  sampler,
  z,
  family = NULL,
  n.trees.treatment = 50L,
  treatment.base = 0.25,
  treatment.power = 3,
  sd.control = NULL,
  sd.moderate = 1,
  b.prior.variance = 0.5,
  update.a = TRUE,
  update.b = TRUE,
  moderators = NULL,
  mu.interactions = NULL,
  tau.interactions = NULL,
  mu.blocks = NULL,
  tau.blocks = NULL
) {
  moderators <- resolveModerators(moderators, sampler$data)
  # the bridge derives the family from this model's own slot, so a supplied name
  # is written into the copy this call hands it; an unknown one is refused
  # there, against the response shape, as every other route's is. Resolved
  # BEFORE the transport below rather than beside the .Call, because
  # sd.control's default is the family's: a NULL left inside as.double(c(...))
  # collapses that vector to length 7 and surfaces as the bridge's length-8
  # refusal, naming a shape rather than the family.
  model <- sampler$model
  if (!is.null(family)) {
    model@family <- family
  }
  if (is.null(sd.control)) {
    # the same helper the public route takes (R/model.R), so the two routes
    # stay on one set of defaults, and a family this route cannot build still
    # meets the bridge's own refusal rather than one from here
    sd.control <- defaultAmplitudePriorScale(model@family)
  }
  # the K-length per-forest transport, at K = 2: the prognostic forest carries
  # a plain scalar amplitude under the half-Cauchy scale mixture sd.control is
  # the median of, and the treatment forest a fixed-variance amplitude pair
  # whose magnitude rides its node scale through the half-normal median
  bcfParams <- list(
    as.double(c(0, 0, 0, 1, 1, 1, sd.control, update.a)),
    as.double(c(
      n.trees.treatment,
      treatment.base,
      treatment.power,
      sd.moderate,
      0.674,
      b.prior.variance,
      0,
      update.b
    ))
  )
  # per-forest interaction constraints: mu and tau each resolve their own
  # interactions() prior against the shared
  # design, so the treatment forest can be capped additive-or-low-order while
  # the prognostic forest stays free (the calibrated-additivity causal use)
  muInteractions <- resolveInteractions(mu.interactions, sampler$data)
  tauInteractions <- resolveInteractions(tau.interactions, sampler$data)

  # per-forest block-additive constraints (each whole tree confined to one
  # declared column group): mu partitions the full
  # design over its own tree count; tau partitions its available columns (the
  # moderator subset if restricted, else the full design) over the treatment
  # tree count, and the engine intersects tau's block rows with the moderator
  # mask at install. The deterministic capacity consumes no rng.
  muBlocks <- resolveBlocks(mu.blocks, sampler$data, sampler$control@n.trees)
  tauBlocks <- resolveBlocks(
    tau.blocks,
    sampler$data,
    as.integer(n.trees.treatment),
    availableColumns = moderators
  )

  # the treatment forest's basis is the two-level factor basis of z, whose
  # amplitudes are exactly (b0, b1); the prognostic forest declares none, so
  # the engine gives it the implicit intercept its single amplitude scales
  z <- as.double(z != 0)
  bases <- list(NULL, cbind(1 - z, z))

  result <- new.env(parent = emptyenv())
  result$ptr <- .Call(
    C_dbarts_bartcore_createBCF,
    sampler$control,
    model,
    sampler$data,
    bases,
    bcfParams,
    # resolved 1-based column indices per forest, or NULL for unrestricted
    list(NULL, moderators),
    # resolved interactions() lists (max.order + 0-based forbidden pairs), or
    # NULL; forest 1 is prognostic, forest 2 the treatment forest
    list(muInteractions, tauInteractions),
    # resolved blocks() lists (0-based per-column group + per-group tree
    # capacity), or NULL; per-forest as with the interactions above
    list(muBlocks, tauBlocks)
  )
  # BCF requires dense predictors; track them R-side for the re-quantize surface
  result$x <- rawPredictorMatrix(sampler$data@x)
  result
}

# Integer coercion for an already-validated count matrix. The RANGE check comes
# BEFORE the coercion: storage.mode() turns a value past .Machine$integer.max
# into NA with a warning, which the engine would then report as a negative
# count - a true refusal naming the wrong reason.
asCountMatrix <- function(counts) {
  if (any(counts > .Machine$integer.max)) {
    stop("multinomial counts must be representable as integers")
  }
  storage.mode(counts) <- "integer"
  counts
}

# The n x K category offset shared by the two multinomial creators and the
# setter: NULL means none, anything else must be a numeric n x K matrix with
# every entry finite (an infinite entry propagates through the log-sum-exp
# margin into a NaN for every category of that row). Only the row-centred part
# of the offset is identified - adding a constant to every entry of a row
# leaves the softmax unchanged - so the input is passed through as given rather
# than silently re-centred. `what` names the offset in the refusals, since the
# train rows and the test rows each carry their own.
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

# The creation-time test twin, shaped against the host data object's test rows:
# without them there is nothing for a per-test-row offset to describe, and
# accepting one would leave it silently unread.
validateCategoryTestOffset <- function(offset.test, sampler, K) {
  if (is.null(offset.test)) {
    return(NULL)
  }
  if (is.null(sampler$data@x.test)) {
    stop(
      "multinomial category test offset requires test data: it carries one ",
      "row per test row"
    )
  }
  validateCategoryOffset(
    offset.test,
    nrow(sampler$data@x.test),
    K,
    "category test offset"
  )
}

# A K-forest multinomial (softmax) sampler, internal and single-trial: the
# K symmetric category forests couple through a
# softmax likelihood with an interleaved one-vs-rest Polya-Gamma augmentation
# and a level-centering move. labels are the 0-based category codes (0..K-1),
# one per observation - the response, so the host sampler's own response is
# ignored. K defaults to one past the largest code. The per-forest leaf scale
# follows the multinomial calibration (pi*sqrt(3)/sqrt(2)), not the host node
# prior. The run's train channel carries the K softmax probabilities
# (n.observations x K x n.samples); per-category fits and split counts read
# through bartcoreForestFits / bartcoreForestVariableCounts (0-based forest =
# category). When the data object carries x.test, the run's test channel
# reports the same K softmax probabilities on the held-out rows (the K forests'
# totalTestFits blended by softmax). Under keepTrees, out-of-sample predict()
# replays all K forests' saved trees and softmaxes them (bartcorePredict).
# offset, when given, is the n x K category offset bartcoreSetCategoryOffset
# installs; offset.test the nTest x K one bartcoreSetCategoryTestOffset
# installs over the host data object's x.test. See those two for the semantics
# and for what they refuse.
#
# A thin shim: one-hot expands labels into counts and routes through
# bartcoreMultinomialDataSampler, the same public bartcore_create dispatch a
# direct bart2(family = "multinomial") fit reaches, rather than a dedicated
# creation entry.
bartcoreMultinomialSampler <- function(
  sampler,
  labels,
  K = NULL,
  offset = NULL,
  offset.test = NULL
) {
  labels <- as.integer(labels)
  if (anyNA(labels)) {
    stop("multinomial labels must be integer category codes 0..K-1")
  }
  if (any(labels < 0L)) {
    stop("multinomial labels must be non-negative category codes 0..K-1")
  }
  if (is.null(K)) {
    K <- max(labels) + 1L
  }
  K <- as.integer(K)
  if (any(labels >= K)) {
    stop("multinomial label out of range 0..K-1")
  }
  # one-hot, category-major: the n_i = 1 special case of the grouped-count
  # combiner, byte-identical to it
  counts <- matrix(0L, length(labels), K)
  counts[cbind(seq_along(labels), labels + 1L)] <- 1L
  bartcoreMultinomialDataSampler(sampler, counts, K, offset, offset.test)
}

# The grouped-count analog of bartcoreMultinomialSampler: the response is
# an n x K matrix of nonnegative integer
# counts, column k holding category k's success counts per observation, with
# trials n_i = sum_k counts[i, k] (>= 1). K defaults to the column count. Same
# K-forest softmax engine; the single-trial label path is the special case of a
# one-hot matrix with every trial 1. Validated R-side (safe over fast); the
# engine re-derives the trials and re-checks the invariants. A thin shim;
# see bartcoreMultinomialSampler.
bartcoreMultinomialCountSampler <- function(
  sampler,
  counts,
  K = NULL,
  offset = NULL,
  offset.test = NULL
) {
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
  if (is.null(K)) {
    K <- ncol(counts)
  }
  if (ncol(counts) != K) {
    stop("multinomial counts must have K columns")
  }
  if (K < 2L) {
    stop("multinomial requires at least two categories")
  }
  if (any(rowSums(counts) < 1)) {
    stop("every multinomial count row must have at least one trial (n_i >= 1)")
  }
  bartcoreMultinomialDataSampler(
    sampler,
    asCountMatrix(counts),
    as.integer(K),
    offset,
    offset.test
  )
}

# Shared factory for the two shims above: builds a counts-carrying copy of
# sampler's own data object and creates the K-forest engine through
# bartcore_create's public multinomial dispatch
# (src/R_interface_bartcore.cpp createMultinomialDataHolder), the same entry
# dbarts(family = "multinomial") and getPointer's re-creation branch reach.
# counts is already validated, category-major (n x K), integer-typed.
# offset/offset.test are validated here, against the same creation-time
# checks a dedicated entry point would apply.
bartcoreMultinomialDataSampler <- function(
  sampler,
  counts,
  K,
  offset,
  offset.test
) {
  offset <- validateCategoryOffset(offset, nrow(counts), K)
  offset.test <- validateCategoryTestOffset(offset.test, sampler, K)
  newData <- sampler$data
  newData@counts <- counts
  newData@y <- as.double(rowSums(counts))
  newData@offset.category <- offset
  newData@offset.category.test <- offset.test
  result <- new.env(parent = emptyenv())
  result$ptr <- .Call(
    C_dbarts_bartcore_create,
    sampler$control,
    sampler$model,
    newData,
    "multinomial"
  )
  # the engine keeps no predictor matrix; track it R-side for the re-quantize
  # surface, as the other dense-predictor wrappers do
  result$x <- rawPredictorMatrix(sampler$data@x)
  # K, so the response-side entries can state the expected shape R-side; the
  # engine reads its own (numReportedLocations) and re-checks
  result$K <- K
  result
}

# The R-level calibration surface indexes forests from 1, as R indexes; the
# bridge and the flat C entries count from 0, as the engine does.
resolveForestIndex <- function(forest) {
  forest <- as.integer(forest)
  if (length(forest) != 1L || is.na(forest) || forest < 1L) {
    stop("'forest' must be a single positive integer (1 selects the first)")
  }
  forest - 1L
}

# The sampler's forest count. A COUNT, not a capability probe:
# samplerCarriesAmplitudes and samplerCarriesCounts each answer only for their
# own model, and neither sees a plain single-forest sampler.
bartcoreNumForests <- function(ptr) .Call(C_dbarts_bartcore_numForests, ptr)

bartcoreSetModel <- function(bcSampler, model, data) {
  invisible(.Call(C_dbarts_bartcore_setModel, bcSampler$ptr, model, data))
}

# Drives a low-level bartcore handle (a bcSampler env holding $ptr) directly,
# not a dbartsSampler; cf. bartcoreSamplerRun, the R-level sampler-layer entry.
bartcoreRun <- function(bcSampler, numBurnIn = 0L, numSamples = 1L) {
  .Call(
    C_dbarts_bartcore_run,
    bcSampler$ptr,
    as.integer(numBurnIn),
    as.integer(numSamples)
  )
}

# offset.test takes the shape of the surface: a flat per-row vector for every
# additive family, and for a multinomial handle an nNew x K matrix, one row per
# predicted row, entering the raw fits before the softmax. It is never taken
# from the sampler - these rows are the caller's - so a handle carrying a
# category offset must be given one here, an all-zero matrix being how the
# offset-free surface is asked for.
bartcorePredict <- function(
  bcSampler,
  x.test,
  offset.test = NULL,
  n.threads = 1L
) {
  x.test <- as.matrix(x.test)
  storage.mode(x.test) <- "double"
  if (!is.null(offset.test)) {
    if (!is.null(bcSampler$K)) {
      offset.test <- validateCategoryOffset(
        offset.test,
        nrow(x.test),
        bcSampler$K,
        "predict category offset"
      )
    } else {
      offset.test <- as.double(offset.test)
    }
  }
  .Call(
    C_dbarts_bartcore_predict,
    bcSampler$ptr,
    x.test,
    offset.test,
    n.threads
  )
}
