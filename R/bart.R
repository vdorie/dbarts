# expects n.pars x n.samples x n.chains if n.chains > 1
# returns n.chains x n.samples x n.pars or (n.chains x n.samples) x n.pars
#
# for quantities such as yhat.train, n.pars is n.obs; for sigma, it is 1
# and the dimension is dropped
#
# preserves the per-parameter dimnames if they exist (for ranef)
convertSamplesFromDbartsToBart <-
  function(
    samples,
    n.chains = dim(samples)[length(dim(samples))],
    combineChains = FALSE
  ) {
    if (!combineChains) {
      ifelse_3(
        is.null(dim(samples)),
        length(dim(samples)) == 2L,
        samples,
        t(samples),
        aperm(samples, c(3L, 2L, 1L))
      )
    } else {
      ifelse_3(
        is.null(dim(samples)),
        length(dim(samples)) == 2L,
        samples,
        if (n.chains <= 1L) t(samples) else as.vector(t(samples)),
        {
          x <- NULL
          res <- evalx(dim(samples), t(matrix(samples, x[1L], prod(x[-1L]))))
          if (!is.null(dimnames(samples))) {
            colnames(res) <- dimnames(samples)[[1L]]
          }
          res
        }
      )
    }
  }

# expects n.chains x n.samples x n.pars or (n.chains x n.samples) x n.pars
# returns n.pars x n.samples x n.chains or n.pars x (n.samples x n.chains)
convertSamplesFromBartsToDbarts <- function(
  samples,
  n.chains,
  uncombineChains = FALSE
) {
  if (!uncombineChains) {
    ifelse_3(
      is.null(dim(samples)),
      is.matrix(samples),
      samples,
      t(samples),
      aperm(samples, c(3L, 2L, 1L))
    )
  } else {
    x <- NULL
    if (is.null(dim(samples))) {
      res <- if (n.chains == 1L) {
        samples
      } else {
        matrix(samples, length(samples) %/% n.chains, n.chains)
      }
      evalx(
        dimnames(samples),
        if (!is.null(x)) dimnames(res) <- list(x[[length(x)]], NULL)
      )
      res
    } else {
      res <- if (n.chains == 1L) {
        samples
      } else {
        array(
          t(samples),
          c(ncol(samples), nrow(samples) %/% n.chains, n.chains)
        )
      }
      evalx(
        dimnames(samples),
        if (!is.null(x)) dimnames(res) <- list(x[[length(x)]], NULL, NULL)
      )
      res
    }
  }
}

# input n.samples x n.chains x n.pars, or n.samples x n.pars when n.chains = 1
# output (n.samples * n.chains) x n.pars
combineChains <- function(samples) {
  ifelse_3(
    is.null(dim(samples)),
    length(dim(samples)) == 2L,
    samples,
    as.vector(samples),
    {
      x <- NULL
      res <- evalx(
        dim(samples),
        matrix(aperm(samples, c(2L, 1L, 3L)), prod(x[1L:2L]), x[3L])
      )
      if (!is.null(dimnames(samples))) {
        dimnames(res) <- evalx(dimnames(samples), list(NULL, x[[length(x)]]))
      }
      res
    }
  )
}

uncombineChains <- function(samples, n.chains) {
  if (is.null(dim(samples))) {
    if (n.chains == 1L) {
      samples
    } else {
      matrix(samples, n.chains, length(samples) %/% n.chains)
    }
  } else {
    res <- if (n.chains == 1L) {
      samples
    } else {
      aperm(
        array(
          samples,
          c(dim(samples)[1L] %/% n.chains, n.chains, ncol(samples))
        ),
        c(2L, 1L, 3L)
      )
    }
    if (!is.null(dimnames(samples))) {
      dimnames(res) <- list(NULL, NULL, dimnames(samples)[[2L]])
    }
    res
  }
}

packageBartResults <- function(
  fit,
  samples,
  burnInSigma,
  burnInK,
  combineChains,
  keepSampler
) {
  responseIsBinary <- fit$control@binary
  n.chains <- fit$control@n.chains

  yhat.train <- NULL
  yhat.train.mean <- NULL
  if (fit$control@keepTrainingFits) {
    yhat.train <- convertSamplesFromDbartsToBart(
      samples$train,
      n.chains,
      combineChains
    )
    if (!responseIsBinary) {
      yhat.train.mean <- apply(yhat.train, length(dim(yhat.train)), mean)
    }
  }

  yhat.test <- NULL
  yhat.test.mean <- NULL
  if (NROW(fit$data@x.test) > 0) {
    yhat.test <- convertSamplesFromDbartsToBart(
      samples$test,
      n.chains,
      combineChains
    )
    if (!responseIsBinary) {
      yhat.test.mean <- apply(yhat.test, length(dim(yhat.test)), mean)
    }
  }

  if (!responseIsBinary) {
    sigma <- convertSamplesFromDbartsToBart(
      samples$sigma,
      n.chains,
      combineChains
    )
  }

  varcount <- convertSamplesFromDbartsToBart(
    samples$varcount,
    n.chains,
    combineChains
  )
  if (!is.null(colnames(fit$data@x)) && !is.null(dim(varcount))) {
    dimnames(varcount) <- if (length(dim(varcount)) > 2L) {
      list(NULL, NULL, colnames(fit$data@x))
    } else {
      list(NULL, colnames(fit$data@x))
    }
  }

  varprobs <- NULL
  if (!is.null(samples[["varprobs"]])) {
    varprobs <- convertSamplesFromDbartsToBart(
      samples$varprobs,
      n.chains,
      combineChains
    )
    if (!is.null(colnames(fit$data@x)) && !is.null(dim(varprobs))) {
      dimnames(varprobs) <- if (length(dim(varprobs)) > 2L) {
        list(NULL, NULL, colnames(fit$data@x))
      } else {
        list(NULL, colnames(fit$data@x))
      }
    }
  }

  if (!is.null(burnInSigma)) {
    burnInSigma <- convertSamplesFromDbartsToBart(
      burnInSigma,
      n.chains,
      combineChains
    )
  }
  if (responseIsBinary) {
    result <- list(
      call = fit$control@call,
      family = fit$model@family,
      yhat.train = yhat.train,
      yhat.test = yhat.test,
      varcount = varcount,
      binaryOffset = fit$data@offset,
      # kept so residuals() works on binary fits (y - fitted); for gaussian
      # fits y rides the other branch
      y = fit$data@y
    )
  } else {
    result <- list(
      call = fit$control@call,
      family = fit$model@family,
      first.sigma = burnInSigma,
      sigma = sigma,
      sigest = fit$data@sigma,
      yhat.train = yhat.train,
      yhat.train.mean = yhat.train.mean,
      yhat.test = yhat.test,
      yhat.test.mean = yhat.test.mean,
      varcount = varcount,
      y = fit$data@y
    )
  }

  if (!is.null(varprobs)) {
    result$varprobs <- varprobs
  }

  if (keepSampler) {
    result$fit <- fit
  } else {
    result$n.chains <- n.chains
  }
  if (!is.null(samples[["k"]])) {
    result[["k"]] <- convertSamplesFromDbartsToBart(
      samples[["k"]],
      n.chains,
      combineChains
    )
    result[["first.k"]] <- convertSamplesFromDbartsToBart(
      burnInK,
      n.chains,
      combineChains
    )
  }

  class(result) <- "bart"
  invisible(result)
}

.kDefault <- quote(if (control@binary) quote(chi(1.5, 2.0)) else 2)

## Builds the quoted tree/node/resid prior calls that bart, bart2, and
## rbart_vi hand to dbarts. nodeK is the node prior's k argument exactly as
## it should enter the call - unevaluated for functions that redirect their
## matched call, evaluated for those that forward through do.call from
## internal frames - or NULL for no node prior. dart may be FALSE, TRUE, or
## a dbartsDartPrior spec; splitProbsName is the caller's argument spelling
## and splitProbsDefault its formal default.
buildSamplerPriors <- function(
  matchedCall,
  power,
  base,
  sigdf,
  sigquant,
  nodeK,
  dart = FALSE,
  splitProbsName = "split.probs",
  splitProbsDefault = NULL
) {
  if (inherits(dart, "dbartsDartPrior")) {
    # a full spec overrides the power/base arguments with its own
    tree.prior <- dart
  } else if (isTRUE(dart)) {
    if (splitProbsName %in% names(matchedCall)) {
      stop(
        "'",
        splitProbsName,
        "' cannot be combined with 'dart': a DART prior samples its split probabilities"
      )
    }
    tree.prior <- quote(dart(power, base))
    tree.prior[[2L]] <- power
    tree.prior[[3L]] <- base
  } else if (!isFALSE(dart)) {
    stop("'dart' must be TRUE, FALSE, or a prior created by dbartsPriors$dart")
  } else {
    tree.prior <- quote(cgm(power, base, split.probs))
    tree.prior[[2L]] <- power
    tree.prior[[3L]] <- base
    if (splitProbsName %in% names(matchedCall)) {
      tree.prior[[4L]] <- matchedCall[[splitProbsName]]
    } else {
      tree.prior[[4L]] <- splitProbsDefault
    }
  }

  if (!is.null(nodeK)) {
    node.prior <- quote(normal(k))
    node.prior[[2L]] <- nodeK
  } else {
    node.prior <- NULL
  }

  resid.prior <- quote(chisq(sigdf, sigquant))
  resid.prior[[2L]] <- sigdf
  resid.prior[[3L]] <- sigquant

  list(
    tree.prior = tree.prior,
    node.prior = node.prior,
    resid.prior = resid.prior
  )
}

bart2 <- function(
  formula,
  data,
  test,
  subset,
  weights,
  offset,
  offset.test = offset,
  sigest = NA_real_,
  sigdf = 3.0,
  sigquant = 0.90,
  k = NULL,
  power = 2.0,
  base = 0.95,
  split.probs = 1 / num.vars,
  dart = FALSE,
  n.trees = 75L,
  n.samples = 500L,
  n.burn = 500L,
  n.chains = 4L,
  n.threads = min(dbarts::guessNumCores(), n.chains),
  combineChains = TRUE,
  n.cuts = 100L,
  useQuantiles = FALSE,
  n.thin = 1L,
  keepTrainingFits = TRUE,
  printEvery = 100L,
  printCutoffs = 0L,
  verbose = TRUE,
  keepTrees = FALSE,
  keepCall = TRUE,
  samplerOnly = FALSE,
  seed = NA_integer_,
  proposal.probs = NULL,
  keepSampler = keepTrees,
  warm.start = NULL,
  n.grow.sweeps = 0L,
  factors = c("categorical", "indicators"),
  family = c("auto", "gaussian", "probit", "logistic", "aft"),
  missing = c("incorporate", "error"),
  ...
) {
  matchedCall <- match.call()
  callingEnv <- parent.frame()

  argNames <- names(matchedCall)[-1L]
  unknownArgs <- argNames %not_in%
    names(formals(dbarts::bart2)) &
    argNames %not_in% names(formals(dbarts::dbartsControl))
  if (any(unknownArgs)) {
    stop(
      "unknown arguments: '",
      paste0(argNames[unknownArgs], collapse = "', '"),
      "'"
    )
  }

  controlCall <- redirectCall(matchedCall, dbarts::dbartsControl)
  missingDefaultArgs <- names(formals(bart2))[
    names(formals(bart2)) %in%
      names(formals(dbarts::dbartsControl)) &
      names(formals(bart2)) %not_in% names(matchedCall)
  ]
  if (length(missingDefaultArgs) > 0L) {
    currentEnv <- sys.frame(sys.nframe())
    controlCall[missingDefaultArgs] <- lapply(
      formals(bart2)[missingDefaultArgs],
      eval,
      envir = currentEnv
    )
  }
  if (!is.na(seed)) {
    controlCall[["rngSeed"]] <- seed
  }
  control <- eval(controlCall, envir = callingEnv)

  control@call <- if (keepCall) matchedCall else call("NULL")
  control@n.burn <- control@n.burn %/% control@n.thin
  control@n.samples <- control@n.samples %/% control@n.thin
  control@printEvery <- control@printEvery %/% control@n.thin
  keepSampler <- keepSampler || control@keepTrees

  if (control@n.burn == 0L && keepTrees == TRUE) {
    control@keepTrees <- TRUE
  }
  if (control@n.burn > 0L) {
    control@keepTrees <- FALSE
  }

  # k enters unevaluated: bart2 redirects its matched call, so the stored
  # symbol resolves in the caller's frame
  priors <- buildSamplerPriors(
    matchedCall,
    power,
    base,
    sigdf,
    sigquant,
    nodeK = matchedCall[["k"]],
    dart = dart,
    splitProbsDefault = formals(dbarts::bart2)[["split.probs"]]
  )
  tree.prior <- priors$tree.prior
  node.prior <- priors$node.prior
  resid.prior <- priors$resid.prior

  samplerCall <- redirectCall(matchedCall, dbarts::dbarts)
  samplerCall$control <- control
  samplerCall$n.samples <- NULL
  samplerCall$tree.prior <- tree.prior
  samplerCall$node.prior <- node.prior
  samplerCall$resid.prior <- resid.prior
  samplerCall$sigma <- as.numeric(sigest)

  sampler <- eval(samplerCall, envir = callingEnv)
  if (samplerOnly == TRUE) {
    return(sampler)
  }

  control <- sampler$control

  # the initial forest: a warm-start donor, a grow-from-root warm start, or a
  # draw from the prior (the default, byte-identical to before this argument)
  n.grow.sweeps <- as.integer(n.grow.sweeps)[1L]
  if (is.na(n.grow.sweeps) || n.grow.sweeps < 0L) {
    stop("'n.grow.sweeps' must be a non-negative integer")
  }
  if (!is.null(warm.start) && n.grow.sweeps > 0L) {
    stop(
      "'warm.start' and 'n.grow.sweeps' both request an initialization; ",
      "supply at most one"
    )
  }

  if (!is.null(warm.start)) {
    sampler$installTrees(warm.start)
  } else if (n.grow.sweeps > 0L) {
    sampler$growFromRoot(n.grow.sweeps, updateState = FALSE)
  } else {
    sampler$sampleTreesFromPrior(updateState = FALSE)
  }

  burnInSigma <- NULL
  burnInK <- NULL
  if (n.burn > 0L) {
    oldX.test <- sampler$data@x.test
    oldOffset.test <- sampler$data@offset.test

    oldKeepTrainingFits <- control@keepTrainingFits
    oldVerbose <- control@verbose

    if (length(oldX.test) > 0L) {
      sampler$setTestPredictorAndOffset(NULL, NULL)
    }
    control@keepTrainingFits <- FALSE
    control@verbose <- FALSE
    sampler$setControl(control)

    samples <- sampler$run(0L, control@n.burn, updateState = FALSE)
    if (!is.null(samples$sigma)) {
      burnInSigma <- samples$sigma
    }
    if (!is.null(samples[["k"]])) {
      burnInK <- samples[["k"]]
    }

    if (length(oldX.test) > 0L) {
      sampler$setTestPredictorAndOffset(oldX.test, oldOffset.test)
    }
    control@keepTrainingFits <- oldKeepTrainingFits
    control@verbose <- oldVerbose
    if (keepTrees == TRUE) {
      control@keepTrees <- TRUE
    }
    sampler$setControl(control)

    samples <- sampler$run(0L, control@n.samples, updateState = FALSE)
  } else {
    samples <- sampler$run(updateState = FALSE)
  }

  result <- packageBartResults(
    sampler,
    samples,
    burnInSigma,
    burnInK,
    combineChains,
    keepSampler
  )
  # needed to extract ppd
  if (!is.null(sampler$data@weights) && length(sampler$data@weights) > 0L) {
    result$weights <- sampler$data@weights
    if (
      !is.null(sampler$data@weights.test) &&
        length(sampler$data@weights.test) > 0L
    ) {
      result$weights.test <- sampler$data@weights.test
    }
  }

  result
}

# Survival-probability draws from an AFT fit (docs/design/survival.md).
# Under the log-normal model log T = f(x) + sigma eps, S(t | x) =
# 1 - Phi((log t - f(x)) / sigma), evaluated at every posterior draw of f
# and sigma. Returns draws per the package's three-tier convention
# (extract = draws, fitted = mean, ci.level = interval): users take means
# and quantiles over the draw margin themselves. newdata predicts out of
# sample (requires a fit kept with keepTrees); NULL uses the training fits.
survivalProbabilities.bart <- function(
  object,
  times,
  newdata = NULL,
  combineChains = TRUE,
  ...
) {
  if (!identical(object[["family"]], "aft")) {
    stop("survivalProbabilities requires an aft (survival) fit")
  }
  times <- as.double(times)
  if (length(times) == 0L || any(!is.finite(times)) || any(times <= 0)) {
    stop("'times' must be finite and positive")
  }

  n.chains <- if (!is.null(object[["fit"]])) {
    object$fit$control@n.chains
  } else {
    object$n.chains
  }

  # work in the uncombined convention (chains x samples x observations),
  # where the sigma draws align with the fit draws unambiguously - the
  # loglik channel's approach; combine the chain margin at the end
  linearPredictor <- if (is.null(newdata)) {
    extract(object, type = "bart", sample = "train", combineChains = FALSE)
  } else {
    predict(object, newdata, type = "bart", combineChains = FALSE)
  }
  sigma <- object[["sigma"]]
  if (is.null(dim(sigma))) {
    sigma <- uncombineChains(as.vector(sigma), n.chains)
  }

  lpDims <- dim(linearPredictor)
  drawDims <- lpDims[-length(lpDims)]
  numObservations <- lpDims[length(lpDims)]
  numTimes <- length(times)
  if (prod(drawDims) != length(sigma)) {
    stop("the fit's draw count does not match its sigma draws")
  }

  result <- array(0, c(drawDims, numTimes, numObservations))
  # in column-major order the draw margins vary fastest, so the sigma draws
  # recycle across observations exactly as pointwiseLogLikelihood's do
  scale <- as.vector(sigma)
  if (length(drawDims) == 1L) {
    for (j in seq_len(numTimes)) {
      result[, j, ] <- pnorm(
        (log(times[j]) - linearPredictor) / scale,
        lower.tail = FALSE
      )
    }
  } else {
    for (j in seq_len(numTimes)) {
      result[,, j, ] <- pnorm(
        (log(times[j]) - as.vector(linearPredictor)) / scale,
        lower.tail = FALSE
      )
    }
  }

  if (combineChains && length(drawDims) > 1L) {
    # merge the chain and sample margins as combineChains() does, samples
    # varying fastest within each chain
    result <- aperm(result, c(2L, 1L, 3L, 4L))
    dim(result) <- c(prod(drawDims), numTimes, numObservations)
  }
  result
}

# Grouped (random-intercept) AFT fits are not yet reachable from R:
# rbart_vi has no 'family' argument, so an rbart fit is never an aft model,
# and a survival curve that dropped the fitted intercepts would be wrong.
survivalProbabilities.rbart <- function(object, ...) {
  stop(
    "survival probabilities are not available for rbart fits: rbart_vi ",
    "has no 'family' argument, so it cannot fit an aft (survival) model; ",
    "grouped survival support is planned"
  )
}

bart <- function(
  x.train,
  y.train,
  x.test = matrix(0.0, 0L, 0L),
  sigest = NA_real_,
  sigdf = 3.0,
  sigquant = 0.90,
  k = 2.0,
  power = 2.0,
  base = 0.95,
  splitprobs = 1 / numvars,
  binaryOffset = 0.0,
  weights = NULL,
  ntree = 200L,
  ndpost = 1000L,
  nskip = 100L,
  printevery = 100L,
  keepevery = 1L,
  keeptrainfits = TRUE,
  usequants = FALSE,
  numcut = 100L,
  printcutoffs = 0L,
  verbose = TRUE,
  nchain = 1L,
  nthread = 1L,
  combinechains = TRUE,
  keeptrees = FALSE,
  keepcall = TRUE,
  sampleronly = FALSE,
  seed = NA_integer_,
  proposalprobs = NULL,
  keepsampler = keeptrees
) {
  # coerce eagerly, naming the argument as the caller typed it - dbartsControl
  # re-coerces its own (already-integer) inputs and would otherwise blame its
  # internal slot names (n.burn/n.trees/...) for a bad ntree/nskip/... value
  ntree <- coerceOrError(ntree, "integer")
  nskip <- coerceOrError(nskip, "integer")
  nchain <- coerceOrError(nchain, "integer")
  nthread <- coerceOrError(nthread, "integer")
  keepevery <- coerceOrError(keepevery, "integer")
  printevery <- coerceOrError(printevery, "integer")
  printcutoffs <- coerceOrError(printcutoffs, "integer")
  numcut <- coerceOrError(numcut, "integer")
  ndpost <- coerceOrError(ndpost, "integer")

  control <- dbartsControl(
    keepTrainingFits = as.logical(keeptrainfits),
    useQuantiles = as.logical(usequants),
    keepTrees = FALSE,
    n.burn = nskip,
    n.trees = ntree,
    n.chains = nchain,
    n.threads = nthread,
    n.thin = keepevery,
    printEvery = printevery,
    printCutoffs = printcutoffs,
    n.cuts = numcut,
    rngSeed = as.integer(seed)
  )
  matchedCall <- if (keepcall) match.call() else call("NULL")
  control@call <- matchedCall
  control@n.burn <- control@n.burn %/% control@n.thin
  control@printEvery <- control@printEvery %/% control@n.thin
  keepsampler <- keepsampler || control@keepTrees
  if (control@n.burn == 0L && keeptrees == TRUE) {
    control@keepTrees <- TRUE
  }
  if (control@n.burn > 0L) {
    control@keepTrees <- FALSE
  }
  ndpost <- ndpost %/% control@n.thin

  priors <- buildSamplerPriors(
    matchedCall,
    power,
    base,
    sigdf,
    sigquant,
    nodeK = if (!is.null(matchedCall[["k"]])) matchedCall[["k"]] else k,
    splitProbsName = "splitprobs",
    splitProbsDefault = formals(dbarts::bart)[["splitprobs"]]
  )
  tree.prior <- priors$tree.prior
  node.prior <- priors$node.prior
  resid.prior <- priors$resid.prior

  # the frozen BayesTree-compatibility shim keeps dummy expansion
  args <- list(
    formula = x.train,
    data = y.train,
    test = x.test,
    subset = NULL,
    weights = weights,
    offset = binaryOffset,
    verbose = as.logical(verbose),
    n.samples = as.integer(ndpost),
    tree.prior = tree.prior,
    node.prior = node.prior,
    resid.prior = resid.prior,
    proposal.probs = proposalprobs,
    control = control,
    sigma = as.numeric(sigest),
    factors = "indicators",
    missing = "error"
  )
  sampler <- do.call(dbarts::dbarts, args, envir = parent.frame(1L))

  if (sampleronly) {
    return(sampler)
  }

  control <- sampler$control

  burnInSigma <- NULL
  burnInK <- NULL
  if (nskip > 0L) {
    oldX.test <- sampler$data@x.test
    oldOffset.test <- sampler$data@offset.test

    oldKeepTrainingFits <- control@keepTrainingFits
    oldVerbose <- control@verbose

    if (length(x.test) > 0) {
      sampler$setTestPredictorAndOffset(NULL, NULL)
    }
    control@keepTrainingFits <- FALSE
    control@verbose <- FALSE
    sampler$setControl(control)

    samples <- sampler$run(0L, control@n.burn, updateState = FALSE)
    if (!is.null(samples$sigma)) {
      burnInSigma <- samples$sigma
    }
    if (!is.null(samples[["k"]])) {
      burnInK <- samples[["k"]]
    }

    if (length(x.test) > 0) {
      sampler$setTestPredictorAndOffset(oldX.test, oldOffset.test)
    }
    control@keepTrainingFits <- oldKeepTrainingFits
    control@verbose <- oldVerbose
    if (keeptrees == TRUE) {
      control@keepTrees <- TRUE
    }
    sampler$setControl(control)

    samples <- sampler$run(0L, control@n.samples, updateState = FALSE)
  } else {
    samples <- sampler$run(updateState = FALSE)
  }

  result <- packageBartResults(
    sampler,
    samples,
    burnInSigma,
    burnInK,
    combinechains,
    keepsampler
  )

  result
}

makeind <- function(x, all = TRUE) {
  ignored <- all ## for R check # nolint: object_usage_linter.
  makeModelMatrixFromDataFrame(x, TRUE)
}
