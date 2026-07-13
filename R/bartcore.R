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

bartcoreSamplerSetPredictor <- function(
  sampler,
  x,
  column,
  forceUpdate,
  updateCutPoints
) {
  # guard before data@x is swapped: the C entry point refuses too, but the
  # R5 object must not be left holding a half-installed dense matrix
  if (predictorSourceIsSparse(sampler$data@x)) {
    stop(
      "sparse predictors fix the design at creation; make a new sampler instead"
    )
  }

  partialUpdate <- !is.null(forceUpdate) &&
    is.character(forceUpdate) &&
    length(forceUpdate) == 1L &&
    !is.na(forceUpdate) &&
    forceUpdate == "partial"

  if (!is.null(column) && is.character(column)) {
    if (is.null(colnames(sampler$data@x))) {
      stop(
        "column names not specified at initialization, so cannot be ",
        "replaced by name"
      )
    }
    column <- match(column, colnames(sampler$data@x))
    if (anyNA(column)) stop("column name not found in names of current X")
  }

  ptr <- sampler$getPointer()

  if (partialUpdate) {
    if (is.null(column)) {
      stop("partial updates require a single 'column' to be specified")
    }
    if (length(column) != 1L) {
      stop("partial updates can only be applied to a single column")
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
    # observations the scan installed (the interim of design plans 1-2)
    sampler$data@x <- assignIntoPredictorSource(
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

  if (is.null(column)) {
    if (is.matrix(x)) {
      if (ncol(x) != ncol(sampler$data@x)) {
        stop("dimension of x must be equal to ", ncol(sampler$data@x))
      }
      if (nrow(x) != nrow(sampler$data@x)) {
        stop("dimension of x must be equal to ", nrow(sampler$data@x))
      }
    } else if (length(x) != prod(dim(sampler$data@x))) {
      stop("length of new x does not match old")
    }
    # a pointer swap: the engine borrows data@x, so install there first and
    # revert if the transaction rolls back
    x <- if (is.matrix(x)) {
      matrix(as.double(x), nrow(x))
    } else {
      matrix(as.double(x), nrow(sampler$data@x))
    }
    x.old <- sampler$data@x
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
        sampler$data@x <- x.old
        e
      }
    )
    if (inherits(tryResult, "error")) {
      stop(tryResult)
    }
    if (!forceUpdate && !updateSuccessful) sampler$data@x <- x.old
  } else {
    column <- as.integer(column)
    if (any(column < 1L | column > ncol(sampler$data@x))) {
      stop(
        "column '",
        column[which(column < 1L | column > ncol(sampler$data@x))[1L]],
        "' is out of range"
      )
    }
    if (is.matrix(x) && ncol(x) != length(column)) {
      stop(
        "number of columns of new x does not match length of columns to replace"
      )
    }
    if (length(x) != nrow(sampler$data@x) * length(column)) {
      stop("length of new x does not match y")
    }
    updateSuccessful <- .Call(
      C_dbarts_bartcore_updatePredictor,
      ptr,
      as.double(x),
      column,
      forceUpdate,
      updateCutPoints
    )
    # the engine keeps no predictor matrix, so maintain data@x R-side when the
    # update is applied (forceUpdate, or a non-rolled-back transaction); the
    # interim of design plans 1-2, until plan 3 rewires mutation
    if (isTRUE(updateSuccessful)) {
      sampler$data@x <- assignIntoPredictorSource(
        sampler$data@x,
        NULL,
        column,
        as.double(x)
      )
    }
  }

  if (!forceUpdate) updateSuccessful else invisible(NULL)
}

bartcoreSamplerSetResponse <- function(sampler, y, updateScale = FALSE) {
  sampler$data@y <- as.double(y)
  .Call(
    C_dbarts_bartcore_setResponse,
    sampler$getPointer(),
    sampler$data@y,
    updateScale
  )
  invisible(NULL)
}

bartcoreSamplerSetOffset <- function(sampler, offset, updateScale) {
  # a synced test offset follows the regular one;
  # NA marks "leave the test offset alone"
  offset.test <- NA
  if (is.null(offset)) {
    if (identical(sampler$data@testUsesRegularOffset, TRUE)) {
      offset.test <- NULL
    }
  } else {
    offset <- as.double(offset)
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
  if (!is.null(column) && is.character(column)) {
    if (is.null(colnames(sampler$data@x))) {
      stop(
        "column names not specified at initialization, so cannot be ",
        "replaced by name"
      )
    }
    column <- match(column, colnames(sampler$data@x))
    if (anyNA(column)) stop("column name not found in names of current X")
  }
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
  if (!is.null(column) && is.character(column)) {
    if (is.null(colnames(sampler$data@x.test))) {
      stop(
        "column names not specified at initialization, so cannot be ",
        "replaced by name"
      )
    }
    column <- match(column, colnames(sampler$data@x.test))
    if (anyNA(column)) {
      stop("column name not found in names of current test predictor matrix")
    }
  }

  if (is.null(column)) {
    # NULL removes the test data; the bridge clears any test offset with it
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
    if (length(x.test) != nrow(sampler$data@x.test) * length(column)) {
      stop("length of new x does not match old x.test")
    }
    # the engine replaces the whole matrix; column updates copy-modify it
    new.x.test <- sampler$data@x.test
    new.x.test[, column] <- as.double(x.test)
    x.test <- new.x.test
  }

  sampler$data@x.test <- x.test
  if (is.null(x.test)) {
    sampler$data@offset.test <- NULL
  }
  .Call(
    C_dbarts_bartcore_setTestPredictor,
    sampler$getPointer(),
    sampler$data@x.test
  )
  invisible(NULL)
}

# Internal interface used by the tests and the equivalence harness; a
# bartcore handle is constructed from the validated control/model/data of an
# existing dbartsSampler regardless of that sampler's own engine. family
# selects the response model for binary responses ("probit", the default,
# or "logistic"); it is the only access to the logistic sampler until the
# new public surface lands. Likewise, columns marked CATEGORICAL (1L) in
# data@varTypes split by category subset rather than by threshold - a
# bartcore-only capability with no public ingestion path yet, so callers
# flip the type by hand; values must be integer codes 0..K-1, K <= 32.

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
# public-surface.md section 5, internal and unserializable. control
# contributes useQuantiles; data contributes x, the column types, and
# n.cuts.
bartcoreDataHandle <- function(control, data) {
  result <- new.env(parent = emptyenv())
  result$ptr <- .Call(C_dbarts_bartcore_createDataHandle, control, data)
  result
}

# A sampler over a row subset of a handle: it copies the handle's cut grid
# and gathers its rows' codes, so folds bin identically to the full data.
# data is the full data object the handle was built from; trainRows and
# testRows index its rows, y/weights/offset are sliced by trainRows, and a
# test offset comes from offset[testRows] (xbart's fold semantics). The
# result refuses raw-predictor mutation (setPredictor and friends, setData,
# setCutPoints, setState); family is as bartcoreSampler's.
bartcoreSamplerFromHandle <- function(
  handle,
  control,
  model,
  data,
  trainRows,
  testRows = NULL,
  family = ""
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
    as.character(family)
  )
  # a view refuses the raw-predictor and re-quantize surface, and its live-tree
  # getTrees needs no training replay, so it carries no predictor matrix
  result$x <- NULL
  result
}

# A BCF two-forest sampler (docs/design/bcf.md), internal and gaussian only:
# the model spec is the prognostic forest mu(x, pihat) - the caller supplies
# pihat as an ordinary predictor column - the arguments below the treatment
# forest tau(x) and the glue. z is the 0/1 treatment. sd.control and
# sd.moderate are the prognostic and effect magnitudes in sd(y) units; the
# calibration map converts them to the per-forest leaf scales at creation,
# overriding the host model's node prior and k for the mu forest. update.a /
# update.b hold the matching glue block fixed when FALSE.
bartcoreBCFSampler <- function(
  sampler,
  z,
  n.trees.treatment = 50L,
  treatment.base = 0.25,
  treatment.power = 3,
  sd.control = 2,
  sd.moderate = 1,
  b.prior.variance = 0.5,
  update.a = TRUE,
  update.b = TRUE
) {
  bcfParams <- as.double(c(
    n.trees.treatment,
    treatment.base,
    treatment.power,
    sd.control,
    sd.moderate,
    b.prior.variance,
    update.a,
    update.b
  ))
  result <- new.env(parent = emptyenv())
  result$ptr <- .Call(
    C_dbarts_bartcore_createBCF,
    sampler$control,
    sampler$model,
    sampler$data,
    as.double(z),
    bcfParams
  )
  # BCF requires dense predictors; track them R-side for the re-quantize surface
  result$x <- rawPredictorMatrix(sampler$data@x)
  result
}

# The 0/1 treatment the treatment forest contrasts on; re-forms b_{z_i} and
# both residuals on the next run.
bartcoreSetTreatment <- function(bcSampler, z) {
  invisible(.Call(C_dbarts_bartcore_setTreatment, bcSampler$ptr, as.double(z)))
}

# The glue on the combining response, a 3 x n.chains matrix of (a, b0, b1).
bartcoreBCFGlue <- function(bcSampler) {
  .Call(C_dbarts_bartcore_getBCFGlue, bcSampler$ptr)
}

# A forest's internal-scale function values (0 prognostic, 1 treatment),
# n.observations x n.chains.
bartcoreForestFits <- function(bcSampler, forest) {
  .Call(C_dbarts_bartcore_getForestFits, bcSampler$ptr, as.integer(forest))
}

bartcoreSetModel <- function(bcSampler, model, control, data) {
  invisible(.Call(
    C_dbarts_bartcore_setModel,
    bcSampler$ptr,
    model,
    control,
    data
  ))
}

bartcoreRun <- function(bcSampler, numBurnIn = 0L, numSamples = 1L) {
  .Call(
    C_dbarts_bartcore_run,
    bcSampler$ptr,
    as.integer(numBurnIn),
    as.integer(numSamples)
  )
}

bartcoreSetOffset <- function(bcSampler, offset, updateScale = FALSE) {
  if (!is.null(offset)) {
    offset <- as.double(offset)
  }
  invisible(.Call(
    C_dbarts_bartcore_setOffset,
    bcSampler$ptr,
    offset,
    as.logical(updateScale)
  ))
}

bartcoreSetResponse <- function(bcSampler, y, updateScale = FALSE) {
  invisible(.Call(
    C_dbarts_bartcore_setResponse,
    bcSampler$ptr,
    as.double(y),
    as.logical(updateScale)
  ))
}

bartcoreSetWeights <- function(bcSampler, weights) {
  invisible(.Call(
    C_dbarts_bartcore_setWeights,
    bcSampler$ptr,
    as.double(weights)
  ))
}

bartcoreSetTestOffset <- function(bcSampler, offset.test) {
  if (!is.null(offset.test)) {
    offset.test <- as.double(offset.test)
  }
  invisible(.Call(C_dbarts_bartcore_setTestOffset, bcSampler$ptr, offset.test))
}

bartcoreSetSigma <- function(bcSampler, sigma) {
  invisible(.Call(C_dbarts_bartcore_setSigma, bcSampler$ptr, as.double(sigma)))
}

bartcoreSetData <- function(bcSampler, data) {
  invisible(.Call(C_dbarts_bartcore_setData, bcSampler$ptr, data))
}

bartcoreSetTestPredictor <- function(bcSampler, x.test) {
  x.test <- as.matrix(x.test)
  storage.mode(x.test) <- "double"
  invisible(.Call(C_dbarts_bartcore_setTestPredictor, bcSampler$ptr, x.test))
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
    C_dbarts_bartcore_setPredictor,
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
    C_dbarts_bartcore_updatePredictor,
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
    C_dbarts_bartcore_updatePredictorPerObservation,
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
    C_dbarts_bartcore_updatePredictorPerObservationJointly,
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
    C_dbarts_bartcore_setCutPoints,
    bcSampler$ptr,
    cutPoints,
    as.integer(columns),
    bcSampler$x
  ))
}

bartcoreGetLatents <- function(bcSampler) {
  .Call(C_dbarts_bartcore_getLatents, bcSampler$ptr, NULL)
}

bartcorePredict <- function(bcSampler, x.test, offset.test = NULL) {
  x.test <- as.matrix(x.test)
  storage.mode(x.test) <- "double"
  if (!is.null(offset.test)) {
    offset.test <- as.double(offset.test)
  }
  .Call(C_dbarts_bartcore_predict, bcSampler$ptr, x.test, offset.test)
}

bartcoreGetTrees <- function(
  bcSampler,
  chainNums,
  sampleNums = NULL,
  treeNums,
  current = FALSE,
  newdata = NULL
) {
  if (!is.null(newdata)) {
    newdata <- as.matrix(newdata)
    storage.mode(newdata) <- "double"
  }
  .Call(
    C_dbarts_bartcore_getTrees,
    bcSampler$ptr,
    as.integer(chainNums),
    if (is.null(sampleNums)) NULL else as.integer(sampleNums),
    as.integer(treeNums),
    as.logical(current),
    newdata,
    bcSampler$x
  )
}

bartcoreStoreState <- function(bcSampler) {
  .Call(C_dbarts_bartcore_storeState, bcSampler$ptr)
}

bartcoreSetState <- function(bcSampler, state) {
  invisible(.Call(
    C_dbarts_bartcore_setState,
    bcSampler$ptr,
    state,
    bcSampler$x
  ))
}

bartcoreSampleTreesFromPrior <- function(bcSampler) {
  invisible(.Call(C_dbarts_bartcore_sampleTreesFromPrior, bcSampler$ptr))
}

bartcoreSampleNodeParametersFromPrior <- function(bcSampler) {
  invisible(.Call(
    C_dbarts_bartcore_sampleNodeParametersFromPrior,
    bcSampler$ptr
  ))
}

bartcorePrintTrees <- function(
  bcSampler,
  chainNums = NULL,
  sampleNums = NULL,
  treeNums = NULL
) {
  invisible(.Call(
    C_dbarts_bartcore_printTrees,
    bcSampler$ptr,
    if (is.null(chainNums)) NULL else as.integer(chainNums),
    if (is.null(sampleNums)) NULL else as.integer(sampleNums),
    if (is.null(treeNums)) NULL else as.integer(treeNums)
  ))
}
