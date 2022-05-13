# expects n.pars x n.samples x n.chains if n.chains > 1
# returns n.chains x n.samples x n.pars or (n.chains x n.samples) x n.pars
#
# for quantities such as yhat.train, n.pars is n.obs; for sigma, it is 1
# and the dimension is dropped
#
# preserves the per-parameter dimnames if they exist (for ranef)
convertSamplesFromDbartsToBart <-
  function(samples, n.chains = dim(samples)[length(dim(samples))], combineChains = FALSE)
{
  if (!combineChains) {
    ifelse_3(is.null(dim(samples)), length(dim(samples)) == 2L,
             samples, t(samples), aperm(samples, c(3L, 2L, 1L)))
  } else {
    ifelse_3(is.null(dim(samples)), length(dim(samples)) == 2L,
             samples,
             if (n.chains <= 1L) t(samples) else as.vector(t(samples)),
             {
               x <- NULL
               res <- evalx(dim(samples), t(matrix(samples, x[1L], prod(x[-1L]))))
               if (!is.null(dimnames(samples))) colnames(res) <- dimnames(samples)[[1L]]
               res
             })
  }
}

# expects n.chains x n.samples x n.pars or (n.chains x n.samples) x n.pars
# returns n.pars x n.samples x n.chains or n.pars x (n.samples x n.chains)
convertSamplesFromBartsToDbarts <- function(samples, n.chains, uncombineChains = FALSE)
{
  if (!uncombineChains) {
    ifelse_3(is.null(dim(samples)), is.matrix(samples),
             samples, t(samples), aperm(samples, c(3L, 2L, 1L)))
  } else {
    x <- NULL
    if (is.null(dim(samples))) {
      res <- if (n.chains == 1L) samples else matrix(samples, length(samples) %/% n.chains, n.chains)
      evalx(dimnames(samples), if (!is.null(x)) dimnames(res) <- list(x[[length(x)]], NULL))
      res
    } else {
      res <- if (n.chains == 1L) samples else array(t(samples), c(ncol(samples), nrow(samples) %/% n.chains, n.chains))
      evalx(dimnames(samples), if (!is.null(x)) dimnames(res) <- list(x[[length(x)]], NULL, NULL))
      res
    }
  } 
}

# input n.samples x n.chains x n.pars, or n.samples x n.pars when n.chains = 1
# output (n.samples * n.chains) x n.pars
combineChains <- function(samples) {
  ifelse_3(is.null(dim(samples)), length(dim(samples)) == 2L,
           samples,
           as.vector(samples),
           {
             x <- NULL
             res <- evalx(dim(samples), matrix(aperm(samples, c(2L, 1L, 3L)), prod(x[1L:2L]), x[3L]))
             if (!is.null(dimnames(samples)))
               dimnames(res) <- evalx(dimnames(samples), list(NULL, x[[length(x)]]))
             res
           })
}

uncombineChains <- function(samples, n.chains) {
  if (is.null(dim(samples))) {
    if (n.chains == 1L) samples else matrix(samples, n.chains, length(samples) %/% n.chains)
  } else {
    res <- if (n.chains == 1L) samples else aperm(array(samples, c(dim(samples)[1L] %/% n.chains, n.chains, ncol(samples))), c(2L, 1L, 3L))
    if (!is.null(dimnames(samples))) dimnames(res) <- list(NULL, NULL, dimnames(samples)[[2L]])
    res
  }
}

packageBartResults <- function(fit, samples, burnInSigma, burnInK, combineChains, keepSampler)
{
  responseIsBinary <- fit$control@binary
  n.chains <- fit$control@n.chains
  
  yhat.train <- NULL
  yhat.train.mean <- NULL
  if (fit$control@keepTrainingFits) {
    yhat.train <- convertSamplesFromDbartsToBart(samples$train, n.chains, combineChains)
    if (!responseIsBinary) yhat.train.mean <- apply(yhat.train, length(dim(yhat.train)), mean)
  }

  yhat.test <- NULL
  yhat.test.mean <- NULL
  if (NROW(fit$data@x.test) > 0) {
    yhat.test <- convertSamplesFromDbartsToBart(samples$test, n.chains, combineChains)
    if (!responseIsBinary) yhat.test.mean <- apply(yhat.test, length(dim(yhat.test)), mean)
  }

  if (!responseIsBinary) sigma <- convertSamplesFromDbartsToBart(samples$sigma, n.chains, combineChains)
    
  varcount <- convertSamplesFromDbartsToBart(samples$varcount, n.chains, combineChains)
  if (!is.null(colnames(fit$data@x)) && !is.null(dim(varcount)))
   dimnames(varcount) <- if (length(dim(varcount)) > 2L) list(NULL, NULL, colnames(fit$data@x)) else list(NULL, colnames(fit$data@x))
  
  if (!is.null(burnInSigma)) burnInSigma <- convertSamplesFromDbartsToBart(burnInSigma, n.chains, combineChains)
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
  
  if (keepSampler)
    result$fit <- fit
  else
    result$n.chains <- n.chains
  if (!is.null(samples[["k"]])) {
    result[["k"]] <- convertSamplesFromDbartsToBart(samples[["k"]], n.chains, combineChains)
    result[["first.k"]] <- convertSamplesFromDbartsToBart(burnInK, n.chains, combineChains)
  }
  
  class(result) <- 'bart'
  invisible(result)
}

.kDefault <- quote(if (control@binary) quote(chi(1.25, Inf)) else 2)

bart2 <- function(
  formula, data, test, subset, weights, offset, offset.test = offset,
  sigest = NA_real_, sigdf = 3.0, sigquant = 0.90,
  k = NULL,
  power = 2.0, base = 0.95, split.probs = 1 / num.vars,
  n.trees = 75L,
  n.samples = 500L, n.burn = 500L,
  n.chains = 4L, n.threads = min(guessNumCores(), n.chains), combineChains = FALSE,
  n.cuts = 100L, useQuantiles = FALSE,
  n.thin = 1L, keepTrainingFits = TRUE,
  printEvery = 100L, printCutoffs = 0L,
  verbose = TRUE, keepTrees = FALSE, keepCall = TRUE,
  samplerOnly = FALSE,
  seed = NA_integer_,
  proposal.probs = NULL,
  keepSampler = keepTrees,
  ...
)
{
  matchedCall <- match.call()
  callingEnv <- parent.frame()
  
  argNames <- names(matchedCall)[-1L]
  unknownArgs <- argNames %not_in% names(formals(dbarts::bart2)) & argNames %not_in% names(formals(dbarts::dbartsControl))
  if (any(unknownArgs))
    stop("unknown arguments: '", paste0(argNames[unknownArgs], collapse = "', '"), "'")
  
  controlCall <- redirectCall(matchedCall, dbarts::dbartsControl)
  missingDefaultArgs <- names(formals(bart2))[names(formals(bart2)) %in% names(formals(dbarts::dbartsControl)) &
                                              names(formals(bart2)) %not_in% names(matchedCall)]
  if (length(missingDefaultArgs) > 0L) {
    currentEnv <- sys.frame(sys.nframe())
    controlCall[missingDefaultArgs] <- lapply(formals(bart2)[missingDefaultArgs], eval, envir = currentEnv)
  }
  if (!is.na(seed)) controlCall[["rngSeed"]] <- seed
  control <- eval(controlCall, envir = callingEnv)
  
  control@call <- if (keepCall) matchedCall else call("NULL")
  control@n.burn     <- control@n.burn     %/% control@n.thin
  control@n.samples  <- control@n.samples  %/% control@n.thin
  control@printEvery <- control@printEvery %/% control@n.thin
  keepSampler <- keepSampler || control@keepTrees

  if (control@n.burn == 0L && keepTrees == TRUE) control@keepTrees <- TRUE
  if (control@n.burn > 0L) control@keepTrees <- FALSE
  
  tree.prior <- quote(cgm(power, base, split.probs))
  tree.prior[[2L]] <- power; tree.prior[[3L]] <- base
  if ("split.probs" %in% names(matchedCall))
    tree.prior[[4L]] <- matchedCall$split.probs
  else
    tree.prior[[4L]] <- formals(dbarts::bart2)[["split.probs"]]
  
  if (!is.null(matchedCall[["k"]])) {
    node.prior <- quote(normal(k))
    node.prior[[2L]] <- matchedCall[["k"]]
  } else {
    node.prior <- NULL
  }
  
  resid.prior <- quote(chisq(sigdf, sigquant))
  resid.prior[[2L]] <- sigdf; resid.prior[[3L]] <- sigquant
  
  samplerCall <- redirectCall(matchedCall, dbarts::dbarts)
  samplerCall$control <- control
  samplerCall$n.samples <- NULL
  samplerCall$tree.prior <- tree.prior
  samplerCall$node.prior <- node.prior
  samplerCall$resid.prior <- resid.prior
  samplerCall$sigma <- as.numeric(sigest)
  
  sampler <- eval(samplerCall, envir = callingEnv)
  if (samplerOnly == TRUE) return(sampler)
  
  control <- sampler$control
  
  sampler$sampleTreesFromPrior(updateState = FALSE)
  
  burnInSigma <- NULL
  burnInK     <- NULL
  if (n.burn > 0L) {
    oldX.test <- sampler$data@x.test
    oldOffset.test <- sampler$data@offset.test
    
    oldKeepTrainingFits <- control@keepTrainingFits
    oldVerbose <- control@verbose

    if (length(oldX.test) > 0L)
      sampler$setTestPredictorAndOffset(NULL, NULL)
    control@keepTrainingFits <- FALSE
    control@verbose <- FALSE
    sampler$setControl(control)

    samples <- sampler$run(0L, control@n.burn, FALSE)
    if (!is.null(samples$sigma)) burnInSigma <- samples$sigma
    if (!is.null(samples[["k"]])) burnInK <- samples[["k"]]
    
    if (length(oldX.test) > 0L)
      sampler$setTestPredictorAndOffset(oldX.test, oldOffset.test)
    control@keepTrainingFits <- oldKeepTrainingFits
    control@verbose <- oldVerbose
    if (keepTrees == TRUE) control@keepTrees <- TRUE
    sampler$setControl(control)

    samples <- sampler$run(0L, control@n.samples, updateState = FALSE)
  } else {
    samples <- sampler$run(updateState = FALSE)
  }

  result <- packageBartResults(sampler, samples, burnInSigma, burnInK, combineChains, keepSampler)
  # needed to extract ppd
  if (!is.null(sampler$data@weights) && length(sampler$data@weights) > 0L) {
    result$weights <- sampler$data@weights
    if (!is.null(sampler$data@weights.test) && length(sampler$data@weights.test) > 0L)
      result$weights.test <- sampler$data@weights.test
  }
  
  result
}

bart <- function(
  x.train, y.train, x.test = matrix(0.0, 0L, 0L),
  sigest = NA_real_, sigdf = 3.0, sigquant = 0.90, 
  k = 2.0,
  power = 2.0, base = 0.95, splitprobs = 1 / numvars,
  binaryOffset = 0.0, weights = NULL,
  ntree = 200L,
  ndpost = 1000L, nskip = 100L,
  printevery = 100L, keepevery = 1L, keeptrainfits = TRUE,
  usequants = FALSE, numcut = 100L, printcutoffs = 0L,
  verbose = TRUE, nchain = 1L, nthread = 1L, combinechains = TRUE,
  keeptrees = FALSE, keepcall = TRUE, sampleronly = FALSE,
  seed = NA_integer_,
  proposalprobs = NULL,
  keepsampler = keeptrees
)
{
  control <- dbartsControl(keepTrainingFits = as.logical(keeptrainfits), useQuantiles = as.logical(usequants),
                           keepTrees = FALSE,
                           n.burn = as.integer(nskip), n.trees = as.integer(ntree), n.chains = as.integer(nchain),
                           n.threads = as.integer(nthread), n.thin = as.integer(keepevery),
                           printEvery = as.integer(printevery), printCutoffs = as.integer(printcutoffs),
                           n.cuts = numcut, rngSeed = as.integer(seed))
  matchedCall <- if (keepcall) match.call() else call("NULL")
  control@call <- matchedCall
  control@n.burn <- control@n.burn %/% control@n.thin
  control@printEvery <- control@printEvery %/% control@n.thin
  keepsampler <- keepsampler || control@keepTrees
  if (control@n.burn == 0L && keeptrees == TRUE) control@keepTrees <- TRUE
  if (control@n.burn > 0L) control@keepTrees <- FALSE
  ndpost <- as.integer(ndpost) %/% control@n.thin

  tree.prior <- quote(cgm(power, base, split.probs))
  tree.prior[[2L]] <- power; tree.prior[[3L]] <- base
  if ("splitprobs" %in% names(matchedCall))
    tree.prior[[4L]] <- matchedCall$splitprobs
  else
    tree.prior[[4L]] <- formals(dbarts::bart)[["splitprobs"]]

  node.prior <- quote(normal(k))
  node.prior[[2L]] <- if (!is.null(matchedCall[["k"]])) matchedCall[["k"]] else k

  resid.prior <- quote(chisq(sigdf, sigquant))
  resid.prior[[2L]] <- sigdf; resid.prior[[3L]] <- sigquant
  
  args <- list(formula = x.train, data = y.train, test = x.test, subset = NULL, weights = weights,
               offset = binaryOffset, verbose = as.logical(verbose), n.samples = as.integer(ndpost),
               tree.prior = tree.prior, node.prior = node.prior, resid.prior = resid.prior,
               proposal.probs = proposalprobs, control = control, sigma = as.numeric(sigest))
  sampler <- do.call(dbarts::dbarts, args, envir = parent.frame(1L))
  
  if (sampleronly) return(sampler)
  
  control <- sampler$control
  
  burnInSigma <- NULL
  burnInK     <- NULL
  if (nskip > 0L) {
    oldX.test <- sampler$data@x.test
    oldOffset.test <- sampler$data@offset.test
    
    oldKeepTrainingFits <- control@keepTrainingFits
    oldVerbose <- control@verbose

    if (length(x.test) > 0) sampler$setTestPredictorAndOffset(NULL, NULL)
    control@keepTrainingFits <- FALSE
    control@verbose <- FALSE
    sampler$setControl(control)
    
    samples <- sampler$run(0L, control@n.burn, FALSE)
    if (!is.null(samples$sigma)) burnInSigma <- samples$sigma
    if (!is.null(samples[["k"]])) burnInK <- samples[["k"]]
    
    if (length(x.test) > 0) sampler$setTestPredictorAndOffset(oldX.test, oldOffset.test)
    control@keepTrainingFits <- oldKeepTrainingFits
    control@verbose <- oldVerbose
    if (keeptrees == TRUE) control@keepTrees <- TRUE
    sampler$setControl(control)

    samples <- sampler$run(0L, control@n.samples, updateState = FALSE)
  } else {
    samples <- sampler$run(updateState = FALSE)
  }

  result <- packageBartResults(sampler, samples, burnInSigma, burnInK, combinechains, keepsampler)
  
  result
}

makeind <- function(x, all = TRUE)
{
  ignored <- all ## for R check
  makeModelMatrixFromDataFrame(x, TRUE)
}

