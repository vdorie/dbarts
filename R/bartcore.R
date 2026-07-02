# Internal interface to the generalized core (src/bartcore/); phase-2
# scaffolding, not exported. A bartcore sampler is constructed from the
# validated control/model/data of an existing dbartsSampler; the flag-based
# opt-in replaces this once the mutation surface reaches parity.
#
# The C side borrows vectors, so anything passed to a setter is retained in
# the handle's environment to pin it against gc.

bartcoreSampler <- function(sampler) {
  result <- new.env(parent = emptyenv())
  result$ptr <- .Call(C_dbarts_bartcore_create, sampler$control, sampler$model,
                      sampler$data)
  result$retained <- list(data = sampler$data)
  result
}

bartcoreRun <- function(bcSampler, numBurnIn = 0L, numSamples = 1L)
  .Call(C_dbarts_bartcore_run, bcSampler$ptr, as.integer(numBurnIn),
        as.integer(numSamples))

bartcoreSetOffset <- function(bcSampler, offset, updateScale = FALSE) {
  if (!is.null(offset)) offset <- as.double(offset)
  bcSampler$retained$offset <- offset
  invisible(.Call(C_dbarts_bartcore_setOffset, bcSampler$ptr, offset,
                  as.logical(updateScale)))
}

bartcoreSetResponse <- function(bcSampler, y) {
  y <- as.double(y)
  bcSampler$retained$y <- y
  invisible(.Call(C_dbarts_bartcore_setResponse, bcSampler$ptr, y))
}

bartcoreSetSigma <- function(bcSampler, sigma)
  invisible(.Call(C_dbarts_bartcore_setSigma, bcSampler$ptr, as.double(sigma)))

bartcoreSetTestPredictor <- function(bcSampler, x.test) {
  x.test <- as.matrix(x.test)
  storage.mode(x.test) <- "double"
  bcSampler$retained$x.test <- x.test
  invisible(.Call(C_dbarts_bartcore_setTestPredictor, bcSampler$ptr, x.test))
}

bartcoreSetPredictor <- function(bcSampler, x, forceUpdate = FALSE,
                                 updateCutPoints = FALSE) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  updated <- .Call(C_dbarts_bartcore_setPredictor, bcSampler$ptr, x,
                   as.logical(forceUpdate), as.logical(updateCutPoints))
  # a pointer swap: pin the new matrix only if it was installed
  if (updated) bcSampler$retained$x <- x
  updated
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
