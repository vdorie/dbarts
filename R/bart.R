packageBartResults <- function(fit, samples, burnInSigma = NULL)
{
  responseIsBinary <- fit$control@binary
  
  yhat.train <- NULL
  yhat.train.mean <- NULL
  if (fit$control@keepTrainingFits) {
    yhat.train <- if (fit$control@n.chains <= 1L) t(samples$train) else aperm(samples$train, c(3L, 2L, 1L))
    if (!responseIsBinary) yhat.train.mean <- apply(yhat.train, length(dim(yhat.train)), mean)
  }

  yhat.test <- NULL
  yhat.test.mean <- NULL
  if (NROW(fit$data@x.test) > 0) {
    yhat.test <- if (fit$control@n.chains <= 1L) t(samples$test) else aperm(samples$test, c(3L, 2L, 1L))
    if (!responseIsBinary) yhat.test.mean <- apply(yhat.test, 2, mean)
  }

  if (!responseIsBinary) sigma <- if (fit$control@n.chains <= 1L) samples$sigma else t(samples$sigma)
    
  varcount <- if (fit$control@n.chains <= 1L) t(samples$varcount) else aperm(samples$varcount, c(3L, 2L, 1L))
  
  if (!is.null(burnInSigma) && fit$control@n.chains > 1L) burnInSigma <- t(burnInSigma)
  
  if (responseIsBinary) {
    result <- list(
      call = fit$control@call,
      yhat.train = yhat.train,
      yhat.test = yhat.test,
      varcount = varcount,
      binaryOffset = fit$data@offset)
  } else {
    result <- list(
      call = fit$control@call,
      first.sigma = burnInSigma,
      sigma = sigma,
      sigest = fit$data@sigma,
      yhat.train = yhat.train,
      yhat.train.mean = yhat.train.mean,
      yhat.test = yhat.test,
      yhat.test.mean = yhat.test.mean,
      varcount = varcount,
      y = fit$data@y)
  }
  class(result) <- 'bart'
  invisible(result)
}

bart2 <- function(
  formula, data, test, subset, weights, offset, offset.test = offset,
  sigest = NA_real_, sigdf = 3.0, sigquant = 0.90,
  k = 2.0,
  power = 2.0, base = 0.95,
  n.trees = 75L,
  n.samples = 500L, n.burn = 500L,
  n.chains = 4L, n.threads = min(guessNumCores(), n.chains),
  n.cuts = 100L, useQuantiles = FALSE,
  n.thin = 1L, keepTrainingFits = TRUE,
  printEvery = 100L, printCutoffs = 0L,
  verbose = TRUE, keepCall = TRUE, ...
)
{
  matchedCall <- match.call()
  callingEnv <- parent.frame()
  
  argNames <- names(matchedCall)[-1L]
  unknownArgs <- argNames %not_in% names(formals(dbarts::bart2)) & argNames %not_in% names(formals(dbarts::dbartsControl))
  if (any(unknownArgs))
    stop("unknown arguments: '", paste0(argNames[unknownArgs], collapse = "', '"), "'")
  
  controlCall <- redirectCall(matchedCall, dbarts::dbartsControl)
  control <- eval(controlCall, envir = callingEnv)
  
  control@call <- if (keepCall) matchedCall else call("NULL")
  control@n.burn     <- control@n.burn     %/% control@n.thin
  control@n.samples  <- control@n.samples  %/% control@n.thin
  control@printEvery <- control@printEvery %/% control@n.thin
  
  tree.prior <- quote(cgm(power, base))
  tree.prior[[2L]] <- power; tree.prior[[3L]] <- base

  node.prior <- quote(normal(k))
  node.prior[[2L]] <- k

  resid.prior <- quote(chisq(sigdf, sigquant))
  resid.prior[[2L]] <- sigdf; resid.prior[[3L]] <- sigquant
  
  samplerCall <- redirectCall(matchedCall, dbarts::dbarts)
  samplerCall$control <- control
  samplerCall$sigma <- as.numeric(sigest)
  
  sampler <- eval(samplerCall, envir = callingEnv)
  
  control <- sampler$control
  
  sampler$sampleTreesFromPrior()
  
  burnInSigma <- NULL
  if (n.burn > 0L) {
    oldX.test <- sampler$data@x.test
    oldOffset.test <- sampler$data@offset.test
    
    oldKeepTrainingFits <- control@keepTrainingFits
    oldVerbose <- control@verbose

    if (length(oldX.test) > 0) sampler$setTestPredictorAndOffset(NULL, NULL, updateState = FALSE)
    control@keepTrainingFits <- FALSE
    control@verbose <- FALSE
    sampler$setControl(control)

    burnInSigma <- sampler$run(0L, control@n.burn, FALSE)$sigma
    
    if (length(oldX.test) > 0) sampler$setTestPredictorAndOffset(oldX.test, oldOffset.test, updateState = FALSE)
    control@keepTrainingFits <- oldKeepTrainingFits
    control@verbose <- oldVerbose
    sampler$setControl(control)

    samples <- sampler$run(0L, control@n.samples)
  } else {
    samples <- sampler$run()
  }

  result <- packageBartResults(sampler, samples, burnInSigma)
  rm(sampler, samples)
  
  result
}

bart <- function(
  x.train, y.train, x.test = matrix(0.0, 0L, 0L),
  sigest = NA_real_, sigdf = 3.0, sigquant = 0.90, 
  k = 2.0,
  power = 2.0, base = 0.95,
  binaryOffset = 0.0, weights = NULL,
  ntree = 200L,
  ndpost = 1000L, nskip = 100L,
  printevery = 100L, keepevery = 1L, keeptrainfits = TRUE,
  usequants = FALSE, numcut = 100L, printcutoffs = 0L,
  verbose = TRUE, nchain = 1L, nthread = 1L, keepcall = TRUE
)
{
  control <- dbartsControl(keepTrainingFits = as.logical(keeptrainfits), useQuantiles = as.logical(usequants),
                           n.burn = as.integer(nskip), n.trees = as.integer(ntree), n.chains = as.integer(nchain),
                           n.threads = as.integer(nthread), n.thin = as.integer(keepevery),
                           printEvery = as.integer(printevery), printCutoffs = as.integer(printcutoffs),
                           n.cuts = numcut)
  matchedCall <- if (keepcall) match.call() else call("NULL")
  control@call <- matchedCall
  control@n.burn <- control@n.burn %/% control@n.thin
  control@printEvery <- control@printEvery %/% control@n.thin
  ndpost <- as.integer(ndpost) %/% control@n.thin

  tree.prior <- quote(cgm(power, base))
  tree.prior[[2L]] <- power; tree.prior[[3L]] <- base

  node.prior <- quote(normal(k))
  node.prior[[2L]] <- k

  resid.prior <- quote(chisq(sigdf, sigquant))
  resid.prior[[2L]] <- sigdf; resid.prior[[3L]] <- sigquant
  
  args <- list(formula = x.train, data = y.train, test = x.test, subset = NULL, weights = weights,
               offset = binaryOffset, verbose = as.logical(verbose), n.samples = as.integer(ndpost),
               tree.prior = tree.prior, node.prior = node.prior,
               resid.prior = resid.prior, control = control, sigma = as.numeric(sigest))
  sampler <- do.call(dbarts::dbarts, args, envir = parent.frame(1L))

  control <- sampler$control
  
  burnInSigma <- NULL
  if (nskip > 0L) {
    oldX.test <- sampler$data@x.test
    oldOffset.test <- sampler$data@offset.test
    
    oldKeepTrainingFits <- control@keepTrainingFits
    oldVerbose <- control@verbose

    if (length(x.test) > 0) sampler$setTestPredictorAndOffset(NULL, NULL, updateState = FALSE)
    control@keepTrainingFits <- FALSE
    control@verbose <- FALSE
    sampler$setControl(control)

    burnInSigma <- sampler$run(0L, control@n.burn, FALSE)$sigma
    
    if (length(x.test) > 0) sampler$setTestPredictorAndOffset(oldX.test, oldOffset.test, updateState = FALSE)
    control@keepTrainingFits <- oldKeepTrainingFits
    control@verbose <- oldVerbose
    sampler$setControl(control)

    samples <- sampler$run(0L, control@n.samples)
  } else {
    samples <- sampler$run()
  }

  result <- packageBartResults(sampler, samples, burnInSigma)
  rm(sampler, samples)
  
  result
}

makeind <- function(x, all = TRUE)
{
  ignored <- all ## for R check
  makeModelMatrixFromDataFrame(x, TRUE)
}
