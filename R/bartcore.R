# The bartcore engine behind dbartsSampler(control = dbartsControl(engine =
# "bartcore")). The dbartsSampler methods delegate to the bartcoreSampler*
# functions below, which mirror the classic methods' semantics; the C side
# borrows vectors and pins them in the external pointer's protection slot.
#
# Not yet supported (methods error): setData, setWeights, test offsets,
# keepTrees/predict/getTrees/plotTree/printTrees, state serialization
# (storeState/setState; samplers cannot be restored after save/load), and
# sampling from the prior.

assertClassicEngine <- function(control, methodName) {
  if (control@engine != "classic") {
    stop("'", methodName, "' is not supported by the bartcore engine")
  }
  invisible(NULL)
}

bartcoreSamplerRun <- function(sampler, numBurnIn, numSamples) {
  control <- sampler$control
  if (is.na(numBurnIn)) numBurnIn <- control@n.burn
  if (is.na(numSamples)) numSamples <- control@n.samples
  if (is.na(numSamples)) {
    stop("bartcore engine samplers require 'numSamples' to be specified")
  }

  result <- .Call(C_dbarts_bartcore_run, sampler$getPointer(),
                  as.integer(numBurnIn), as.integer(numSamples))
  if (is.null(result)) return(invisible(NULL))
  result
}

bartcoreSamplerSetPredictor <- function(sampler, x, column, forceUpdate,
                                        updateCutPoints) {
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

    return(.Call(C_dbarts_bartcore_updatePredictorPerObservation, ptr,
                 as.double(x), as.integer(column)))
  }

  forceUpdate <- if (is.null(forceUpdate)) {
    is.null(column)
  } else {
    coerceOrError(forceUpdate, "logical")
  }
  updateCutPoints <- coerceOrError(updateCutPoints, "logical")

  if (is.null(column)) {
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
      updateSuccessful <- .Call(C_dbarts_bartcore_setPredictor, ptr,
                                sampler$data@x, forceUpdate, updateCutPoints),
      error = function(e) {
        sampler$data@x <- x.old
        e
      }
    )
    if (inherits(tryResult, "error")) stop(tryResult)
    if (!forceUpdate && !updateSuccessful) sampler$data@x <- x.old
  } else {
    # written in place through the matrix the engine borrows, so data@x
    # already reflects the change
    updateSuccessful <- .Call(C_dbarts_bartcore_updatePredictor, ptr,
                              as.double(x), as.integer(column), forceUpdate,
                              updateCutPoints)
  }

  if (!forceUpdate) updateSuccessful else invisible(NULL)
}

bartcoreSamplerSetResponse <- function(sampler, y) {
  sampler$data@y <- as.double(y)
  .Call(C_dbarts_bartcore_setResponse, sampler$getPointer(), sampler$data@y)
  invisible(NULL)
}

bartcoreSamplerSetOffset <- function(sampler, offset, updateScale) {
  if (!is.null(offset)) {
    offset <- as.double(offset)
    if (length(offset) == 1L) {
      offset <- rep_len(offset, length(sampler$data@y))
    } else if (length(offset) != length(sampler$data@y)) {
      stop("length of replacement offset is not equal to number of observations")
    }
    if (identical(sampler$data@testUsesRegularOffset, TRUE) &&
        !is.null(sampler$data@x.test)) {
      stop("test offsets are not supported by the bartcore engine")
    }
  }

  sampler$data@offset <- offset
  .Call(C_dbarts_bartcore_setOffset, sampler$getPointer(),
        sampler$data@offset, as.logical(updateScale))
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
  if (is.null(column)) column <- seq_len(ncol(sampler$data@x))

  if (!is.list(cuts)) cuts <- list(cuts)
  cuts <- lapply(cuts, as.double)

  .Call(C_dbarts_bartcore_setCutPoints, sampler$getPointer(), cuts,
        as.integer(column))
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
    x.test <- validateXTest(x.test, attr(sampler$data@x, "term.labels"),
                            ncol(sampler$data@x), colnames(sampler$data@x),
                            attr(sampler$data@x, "drop"))
    if (is.null(x.test)) {
      stop("removing test data is not supported by the bartcore engine")
    }
  } else {
    # the engine replaces the whole matrix; column updates copy-modify it
    new.x.test <- sampler$data@x.test
    new.x.test[, column] <- as.double(x.test)
    x.test <- new.x.test
  }

  sampler$data@x.test <- x.test
  .Call(C_dbarts_bartcore_setTestPredictor, sampler$getPointer(),
        sampler$data@x.test)
  invisible(NULL)
}

# Internal interface used by the tests and the equivalence harness; a
# bartcore handle is constructed from the validated control/model/data of an
# existing dbartsSampler regardless of that sampler's own engine.

bartcoreSampler <- function(sampler) {
  result <- new.env(parent = emptyenv())
  result$ptr <- .Call(C_dbarts_bartcore_create, sampler$control, sampler$model,
                      sampler$data)
  result
}

bartcoreRun <- function(bcSampler, numBurnIn = 0L, numSamples = 1L)
  .Call(C_dbarts_bartcore_run, bcSampler$ptr, as.integer(numBurnIn),
        as.integer(numSamples))

bartcoreSetOffset <- function(bcSampler, offset, updateScale = FALSE) {
  if (!is.null(offset)) offset <- as.double(offset)
  invisible(.Call(C_dbarts_bartcore_setOffset, bcSampler$ptr, offset,
                  as.logical(updateScale)))
}

bartcoreSetResponse <- function(bcSampler, y)
  invisible(.Call(C_dbarts_bartcore_setResponse, bcSampler$ptr, as.double(y)))

bartcoreSetSigma <- function(bcSampler, sigma)
  invisible(.Call(C_dbarts_bartcore_setSigma, bcSampler$ptr, as.double(sigma)))

bartcoreSetTestPredictor <- function(bcSampler, x.test) {
  x.test <- as.matrix(x.test)
  storage.mode(x.test) <- "double"
  invisible(.Call(C_dbarts_bartcore_setTestPredictor, bcSampler$ptr, x.test))
}

bartcoreSetPredictor <- function(bcSampler, x, forceUpdate = FALSE,
                                 updateCutPoints = FALSE) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  .Call(C_dbarts_bartcore_setPredictor, bcSampler$ptr, x,
        as.logical(forceUpdate), as.logical(updateCutPoints))
}

# In-place overwrite of columns in the matrix the sampler borrows, visible
# through the originating data object, exactly like the classic engine.
bartcoreUpdatePredictor <- function(bcSampler, x, columns, forceUpdate = FALSE,
                                    updateCutPoints = FALSE)
  .Call(C_dbarts_bartcore_updatePredictor, bcSampler$ptr, as.double(x),
        as.integer(columns), as.logical(forceUpdate),
        as.logical(updateCutPoints))

bartcoreUpdatePredictorPerObservation <- function(bcSampler, x, column)
  .Call(C_dbarts_bartcore_updatePredictorPerObservation, bcSampler$ptr,
        as.double(x), as.integer(column))

bartcoreUpdatePredictorPerObservationJointly <- function(bcSamplers, x, columns)
  .Call(C_dbarts_bartcore_updatePredictorPerObservationJointly,
        lapply(bcSamplers, function(s) s$ptr), as.double(x),
        as.integer(columns))

# cutPoints is a list of strictly increasing numeric vectors, one per column;
# trees are refreshed unconditionally, collapsing splits the new cuts orphan
bartcoreSetCutPoints <- function(bcSampler, cutPoints, columns) {
  if (!is.list(cutPoints)) cutPoints <- list(cutPoints)
  cutPoints <- lapply(cutPoints, as.double)
  invisible(.Call(C_dbarts_bartcore_setCutPoints, bcSampler$ptr, cutPoints,
                  as.integer(columns)))
}

bartcoreGetLatents <- function(bcSampler)
  .Call(C_dbarts_bartcore_getLatents, bcSampler$ptr)
