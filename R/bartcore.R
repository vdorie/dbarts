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

bartcoreGetLatents <- function(bcSampler)
  .Call(C_dbarts_bartcore_getLatents, bcSampler$ptr)
