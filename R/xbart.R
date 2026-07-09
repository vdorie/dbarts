xbart <- function(
  formula,
  data,
  subset,
  weights,
  offset,
  verbose = FALSE,
  n.samples = 200L,
  method = c("k-fold", "random subsample"),
  n.test = c(5, 0.2),
  n.reps = 40L,
  n.burn = c(200L, 150L, 50L),
  loss = c("rmse", "log", "mcr"),
  n.threads = dbarts::guessNumCores(),
  n.trees = 75L,
  k = NULL,
  power = 2,
  base = 0.95,
  split.probs = NULL,
  dart = FALSE,
  drop = TRUE,
  resid.prior = chisq,
  control = dbarts::dbartsControl(),
  sigma = NA_real_,
  seed = NA_integer_,
  factors = c("categorical", "indicators"),
  family = c("auto", "gaussian", "probit", "logistic"),
  missing = c("incorporate", "error"),
  node.prior = NULL
) {
  matchedCall <- match.call()

  currEnv <- sys.frame(sys.nframe())
  evalEnv <- parent.frame(1L)

  validateCall <- redirectCall(
    matchedCall,
    quoteInNamespace(validateArgumentsInEnvironment),
    control,
    verbose,
    n.samples,
    sigma
  )
  validateCall <- addCallArgument(validateCall, 1L, currEnv)
  validateCall <- addCallArgument(validateCall, 2L, xbart)
  validateCall <- addCallArgument(validateCall, 3L, "xbart")
  eval(validateCall, evalEnv, getNamespace("dbarts"))

  if (control@call != call("NA")[[1L]]) {
    control@call <- matchedCall
  }
  control@n.chains <- 1L
  control@n.threads <- 1L
  control@keepTrees <- FALSE
  control@keepTrainingFits <- FALSE
  control@updateState <- FALSE
  control@verbose <- FALSE

  dataCall <- redirectCall(
    matchedCall,
    quoteInNamespace(dbartsData),
    formula,
    data,
    subset,
    weights,
    offset,
    factors,
    missing
  )
  data <- eval(dataCall, evalEnv)
  data@n.cuts <- rep_len(control@n.cuts, ncol(data@x))
  data@sigma <- sigma

  family <- match.arg(family)
  uniqueResponses <- unique(data@y)
  responseIsBinary <- length(uniqueResponses) == 2L &&
    all(sort(uniqueResponses) == c(0, 1))
  if (family == "auto") {
    family <- if (responseIsBinary) "probit" else "gaussian"
  } else if (family != "gaussian" && !responseIsBinary) {
    # gaussian on a 0/1 response is a legitimate request; the binary
    # families need latent-variable coding
    stop("family \"", family, "\" requires a response coded 0/1")
  }
  control@binary <- family != "gaussian"

  if (is.na(data@sigma) && !control@binary) {
    data@sigma <- estimateSigmaFromLinearModel(data)
  }

  if (control@binary && is.null(matchedCall[["resid.prior"]])) {
    matchedCall[["resid.prior"]] <- quote(fixed(1))
  }

  if (
    !is.character(method) || method[1L] %not_in% eval(formals(xbart)$method)
  ) {
    stop(
      "method must be in '",
      paste0(eval(formals(xbart)$method), collapse = "', '"),
      "'"
    )
  }
  method <- method[1L]
  if (!is.null(matchedCall$method) && is.null(matchedCall$n.test)) {
    n.test <- eval(formals(xbart)$n.test)[match(
      method,
      eval(formals(xbart)$method)
    )]
  }
  n.test <- n.test[1L]

  if (is.null(matchedCall$loss)) {
    loss <- loss[if (!control@binary) 1L else 2L]
  } else if (is.function(loss)) {
    if (length(formals(loss)) != 3L) {
      stop("supplied loss function must take exactly three arguments")
    }
    loss <- list(loss, evalEnv)
  } else if (is.list(loss)) {
    if (!is.function(loss[[1L]])) {
      stop("first member of loss-list must be a function")
    }
    if (length(formals(loss[[1L]])) != 3L) {
      stop("supplied loss function must take exactly three arguments")
    }
    if (!is.environment(loss[[2L]])) {
      stop("second member of loss-list must be an environment")
    }
  }

  if (is.null(matchedCall$n.trees) && "n.trees" %not_in% names(matchedCall)) {
    n.trees <- control@n.trees
  } else {
    n.trees <- coerceOrError(n.trees, "integer")
  }
  if (anyNA(n.trees) || any(n.trees <= 0L)) {
    stop("'n.trees' must contain only positive integers")
  }

  # a supplied node.prior contributes the leaf model shape - normal(k), the
  # default, linear(columns, k), or gp(columns, k, ...), whose designated
  # covariate columns resolve against the model matrix; the k argument
  # drives the k grid as always, with a k inside the supplied prior standing
  # in for a missing k argument
  node.spec <- NULL
  if (!is.null(matchedCall[["node.prior"]])) {
    priorEnv <- new.env(parent = evalEnv)
    priorEnv[["normal"]] <- getNamespace("dbarts")[["normal"]]
    priorEnv[["linear"]] <- getNamespace("dbarts")[["linear"]]
    priorEnv[["gp"]] <- getNamespace("dbarts")[["gp"]]
    priorEnv[["chi"]] <- getNamespace("dbarts")[["chi"]]
    node.spec <- eval(matchedCall[["node.prior"]], priorEnv)
    if (is.function(node.spec)) {
      node.spec <- node.spec()
    }
    if (!is.null(node.spec) && !is(node.spec, "dbartsNodePrior")) {
      stop("'node.prior' must be a node prior specification; see ?dbartsPriors")
    }
  }

  if (is.null(matchedCall[["k"]])) {
    k <- if (!is.null(node.spec) && !is.null(node.spec@k)) {
      node.spec@k
    } else {
      eval(eval(quoteInNamespace(.kDefault)))
    }
  }
  kIsGrid <- is.numeric(k)
  if (kIsGrid && (anyNA(k) || any(k <= 0))) {
    stop("'k' must contain only positive values")
  }

  power <- coerceOrError(power, "numeric")
  base <- coerceOrError(base, "numeric")
  drop <- coerceOrError(drop, "logical")
  if (anyNA(power) || any(power <= 0)) {
    stop("'power' must contain only positive values")
  }
  if (anyNA(base) || any(base <= 0 | base >= 1)) {
    stop("'base' must contain only values in (0, 1)")
  }

  # the tree structure prior. DART (dart = TRUE or a dart() spec) samples its
  # own split probabilities; cgm takes a fixed split.probs. power and base are
  # xbart's grid axes, so the prior is built at the grid's first values and
  # every cell overrides them (see cellModel); a supplied dart() spec
  # contributes its Dirichlet hyperparameters, with power and base still swept.
  if (inherits(dart, "dbartsDartPrior")) {
    tree.prior <- dart
  } else if (isTRUE(dart)) {
    if ("split.probs" %in% names(matchedCall)) {
      stop(
        "'split.probs' cannot be combined with 'dart': a DART prior samples its split probabilities"
      )
    }
    tree.prior <- dbartsPriors$dart(power[1L], base[1L])
  } else if (!isFALSE(dart)) {
    stop("'dart' must be TRUE, FALSE, or a prior created by dbartsPriors$dart")
  } else {
    tree.prior <- cgm(power[1L], base[1L], split.probs)
  }
  tree.prior <- resolveSplitProbabilities(tree.prior, data)

  if (is.null(node.spec)) {
    node.prior <- quote(normal(k))
    node.prior[[1L]] <- quoteInNamespace(normal)
    node.prior[[2L]] <- ifelse_3(is.numeric(k), is.list(k), k[1L], k[[1L]], k)
    node.prior <- eval(node.prior)
  } else {
    kValue <- ifelse_3(is.numeric(k), is.list(k), k[1L], k[[1L]], k)
    if (is.call(kValue)) {
      kValue <- eval(kValue)
    }
    node.prior <- if (is(node.spec, "dbartsLinearPrior")) {
      resolveLeafCovariates(linear(node.spec@columns, kValue), data)
    } else if (is(node.spec, "dbartsGPPrior")) {
      resolveLeafCovariates(
        gp(
          node.spec@columns,
          kValue,
          node.spec@lengthscale,
          node.spec@max.leaf.size
        ),
        data
      )
    } else {
      normal(kValue)
    }
  }
  # xbart cells always run a fixed k (the default 2, not the binary
  # hyperprior), as before the priors became objects
  node.hyperprior <- resolveNodeHyperprior(node.prior@k, binary = FALSE)

  resid.prior <-
    if (
      !is.null(matchedCall$resid.prior) || "resid.prior" %in% names(matchedCall)
    ) {
      env <- new.env(parent = evalEnv)
      env[["chisq"]] <- getNamespace("dbarts")[["chisq"]]
      env[["fixed"]] <- getNamespace("dbarts")[["fixed"]]
      eval(matchedCall$resid.prior, env)
    } else {
      eval(formals(xbart)$resid.prior, getNamespace("dbarts"))()
    }
  if (is.call(resid.prior)) {
    resid.prior <- eval(resid.prior, getNamespace("dbarts"))
  }
  model <- newValidated(
    "dbartsModel",
    tree.prior,
    node.prior,
    node.hyperprior,
    resid.prior,
    family = family,
    # the logistic scale is probit's default widened by the logistic
    # latent's standard deviation, pi / sqrt(3)
    node.scale = switch(
      family,
      gaussian = 0.5,
      probit = 3.0,
      logistic = pi * sqrt(3.0)
    )
  )

  numObservations <- length(data@y)
  if (method == "k-fold") {
    n.test <- coerceOrError(n.test, "integer")
    if (n.test < 2L || n.test > numObservations) {
      stop(
        "for k-fold crossvalidation, n.test must be an integer in [2, ",
        numObservations,
        "]"
      )
    }
    foldSizes <- rep.int(numObservations %/% n.test, n.test) +
      rep.int(
        c(1L, 0L),
        c(numObservations %% n.test, n.test - numObservations %% n.test)
      )
    numTest <- 0L
  } else {
    n.test <- coerceOrError(n.test, "numeric")
    if (n.test > 1) {
      n.test <- n.test / numObservations
    }
    if (n.test <= 0 || n.test >= 1) {
      stop("for random subsample crossvalidation, n.test must be in (0, 1)")
    }
    numTest <- max(
      1L,
      min(numObservations - 1L, as.integer(round(n.test * numObservations)))
    )
    foldSizes <- integer()
  }

  n.reps <- coerceOrError(n.reps, "integer")
  if (is.na(n.reps) || n.reps <= 0L) {
    stop("'n.reps' must be a positive integer")
  }
  n.burn <- rep_len(coerceOrError(n.burn, "integer"), 3L)
  if (anyNA(n.burn) || any(n.burn < 0L)) {
    stop("'n.burn' must contain non-negative integers")
  }
  n.threads <- coerceOrError(n.threads, "integer")
  if (is.na(n.threads) || n.threads <= 0L) {
    stop("'n.threads' must be a positive integer")
  }

  # DART holds its Dirichlet updates until the forest is likelihood-informed;
  # as the fitting functions default it to half the burn-in, default here to
  # half the fresh-sampler burn-in each cell runs
  if (
    is(model@tree.prior, "dbartsDartPrior") &&
      is.na(model@tree.prior@update.delay)
  ) {
    model@tree.prior@update.delay <- as.numeric(n.burn[1L] %/% 2L)
  }

  lossFunction <- xbartLossFunction(loss, control, family)

  # a replication draws a data split and sweeps every parameter cell over
  # it. Chains warm-start only across cells - the training data is
  # identical there, so carrying trees is sound - and never across folds
  # or splits, whose held-out rows the previous training set contained.
  # (The retired implementation warm-started across splits; slow-mixing
  # cells never forgot the previous fold and scored optimistically.)
  # Tree counts are fixed at a sampler's creation, so they vary slowest
  # and each count gets a fresh fit per split.
  kLength <- if (kIsGrid) length(k) else 1L
  cells <- expand.grid(
    iBase = seq_along(base),
    iPower = seq_along(power),
    iK = seq_len(kLength),
    iTrees = seq_along(n.trees)
  )
  numCells <- nrow(cells)

  spec <- namedList(
    control,
    model,
    data,
    method,
    foldSizes,
    numTest,
    n.samples = control@n.samples,
    n.burn,
    n.trees,
    kValues = if (kIsGrid) k else NA_real_,
    kIsGrid,
    power,
    base,
    cells,
    lossFunction
  )

  numChunks <- max(1L, min(n.threads, n.reps))
  chunkIndices <- parallel::splitIndices(n.reps, numChunks)

  # a supplied seed drives every chunk deterministically and leaves the
  # caller's stream untouched; results are reproducible for a fixed
  # (seed, n.threads) pair, since the chunking changes with the latter
  if (!is.na(seed)) {
    oldSeed <- if (exists(".Random.seed", globalenv())) {
      get(".Random.seed", globalenv())
    } else {
      NULL
    }
    on.exit(
      if (!is.null(oldSeed)) {
        assign(".Random.seed", oldSeed, envir = globalenv())
      } else if (exists(".Random.seed", globalenv())) {
        rm(".Random.seed", envir = globalenv())
      },
      add = TRUE
    )
    set.seed(seed)
  }
  chunkSeeds <- sample.int(.Machine$integer.max, numChunks)

  if (verbose) {
    cat(
      "running ",
      numCells,
      " parameter combination",
      if (numCells > 1L) "s" else "",
      " x ",
      n.reps,
      " replication",
      if (n.reps > 1L) "s" else "",
      " on ",
      numChunks,
      " worker",
      if (numChunks > 1L) "s" else "",
      "\n",
      sep = ""
    )
  }

  if (numChunks == 1L) {
    chunkResults <- list(xbartRunChunk(
      spec,
      chunkIndices[[1L]],
      chunkSeeds[1L]
    ))
  } else {
    cluster <- parallel::makeCluster(numChunks)
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    # passing the namespace function itself serializes it by reference,
    # loading dbarts on the workers without shipping this frame
    chunkResults <- parallel::clusterMap(
      cluster,
      xbartRunChunk,
      repIndices = chunkIndices,
      chunkSeed = chunkSeeds,
      MoreArgs = list(spec = spec)
    )
  }
  # rep-major, cells within: chunks hold contiguous rep ranges
  lossValues <- do.call(rbind, chunkResults)
  numResults <- ncol(lossValues)

  # place by index so the array layout is independent of evaluation order
  dims <- c(n.reps, length(n.trees), kLength, length(power), length(base))
  repIndex <- rep(seq_len(n.reps), each = numCells)
  cellIndex <- rep.int(seq_len(numCells), n.reps)
  linearIndex <- repIndex +
    n.reps *
      ((cells$iTrees[cellIndex] - 1L) +
        length(n.trees) *
          ((cells$iK[cellIndex] - 1L) +
            kLength *
              ((cells$iPower[cellIndex] - 1L) +
                length(power) * (cells$iBase[cellIndex] - 1L))))
  result <- array(NA_real_, c(dims, numResults))
  cellCount <- prod(dims)
  for (i in seq_len(numResults)) {
    result[linearIndex + (i - 1L) * cellCount] <- lossValues[, i]
  }

  varNames <- c("n.trees", "k", "power", "base")
  dimIncluded <- c(
    TRUE,
    if (drop) length(n.trees) > 1L else TRUE,
    if (!kIsGrid) {
      FALSE
    } else if (drop) {
      length(k) > 1L
    } else {
      TRUE
    },
    if (drop) length(power) > 1L else TRUE,
    if (drop) length(base) > 1L else TRUE,
    numResults > 1L
  )
  newDims <- c(dims, numResults)[dimIncluded]
  if (length(newDims) == 1L) {
    result <- as.vector(result)
  } else {
    dim(result) <- newDims
    dimNames <- vector("list", length(newDims))
    names(dimNames) <- c(
      "rep",
      varNames[dimIncluded[2L:5L]],
      if (numResults > 1L) "loss"
    )
    for (varName in varNames[dimIncluded[2L:5L]]) {
      x <- get(varName)
      dimNames[[varName]] <- as.character(
        if (is.double(x)) signif(x, 2L) else x
      )
    }
    dimnames(result) <- dimNames
  }

  result
}

## Resolve the loss argument into function(y.test, testSamples, weights):
## testSamples is numTestObservations x numSamples, on the latent scale for
## binary responses. The built-in binary losses transform by the family's
## link.
xbartLossFunction <- function(loss, control, family) {
  if (is.list(loss)) {
    result <- loss[[1L]]
    environment(result) <- loss[[2L]]
    return(result)
  }

  if (!is.character(loss) || loss[1L] %not_in% c("rmse", "log", "mcr")) {
    stop("loss must be in 'rmse', 'log', 'mcr', or a function")
  }
  loss <- loss[1L]
  if (loss %in% c("log", "mcr") && !control@binary) {
    stop("loss '", loss, "' requires a binary response")
  }

  probFromLatent <- if (identical(family, "logistic")) plogis else pnorm

  switch(
    loss,
    rmse = function(y.test, testSamples, weights) {
      y.test.hat <- rowMeans(testSamples)
      if (is.null(weights)) {
        sqrt(mean((y.test - y.test.hat)^2))
      } else {
        sqrt(sum(weights * (y.test - y.test.hat)^2) / sum(weights))
      }
    },
    log = function(y.test, testSamples, weights) {
      p.test <- rowMeans(probFromLatent(testSamples))
      logLikelihood <- ifelse(y.test > 0, log(p.test), log1p(-p.test))
      if (is.null(weights)) {
        -mean(logLikelihood)
      } else {
        -sum(weights * logLikelihood) / sum(weights)
      }
    },
    mcr = function(y.test, testSamples, weights) {
      misclassified <- as.numeric(
        rowMeans(probFromLatent(testSamples)) > 0.5
      ) !=
        y.test
      if (is.null(weights)) {
        mean(misclassified)
      } else {
        sum(weights * misclassified) / sum(weights)
      }
    }
  )
}

## One worker's share of the replications. The predictor store (cuts +
## codes) is built once per chunk; each fold's sampler is a row-subset view
## over it, so every fold bins on the full data's cut grid and no fold
## re-quantizes the predictors. Per data split, every tree count gets a
## fresh sampler burned n.burn[1] iterations; the remaining parameter cells
## sweep warm off it with n.burn[2] iterations each, sound because the
## training data is unchanged. Chains never carry over between splits, whose
## held-out rows the previous training set contained; n.burn[3] is unused.
## Returns a (reps x cells) x numResults matrix, cells in spec$cells order.
xbartRunChunk <- function(spec, repIndices, chunkSeed) {
  set.seed(chunkSeed)

  data <- spec$data
  cells <- spec$cells
  numCells <- nrow(cells)
  numObservations <- length(data@y)
  hasWeights <- !is.null(data@weights)
  family <- if (spec$model@family == "auto") "" else spec$model@family

  handle <- bartcoreDataHandle(spec$control, data)

  cellModel <- function(cell) {
    result <- spec$model
    result@tree.prior@power <- spec$power[cells$iPower[cell]]
    result@tree.prior@base <- spec$base[cells$iBase[cell]]
    if (spec$kIsGrid) {
      result@node.hyperprior <- new(
        "dbartsFixedHyperprior",
        k = spec$kValues[cells$iK[cell]]
      )
    }
    result
  }

  numLossResults <- NULL

  # fit every cell against one split, fresh per tree count, warm across the
  # rest; cells are grouped by iTrees, so a single pass reuses each sampler
  # maximally. The view slices y/weights/offset by row and takes its test
  # offset from offset[testRows], as the retired data-slicing path did.
  sweepCells <- function(testRows) {
    trainRows <- seq_len(numObservations)[-testRows]
    y.test <- data@y[testRows]
    weights.test <- if (hasWeights) data@weights[testRows] else NULL

    lossValues <- NULL
    sampler <- NULL
    cellControl <- NULL
    currentTrees <- NA_integer_
    for (cell in seq_len(numCells)) {
      if (
        is.null(sampler) || spec$n.trees[cells$iTrees[cell]] != currentTrees
      ) {
        cellControl <- spec$control
        cellControl@n.trees <- spec$n.trees[cells$iTrees[cell]]
        sampler <- bartcoreSamplerFromHandle(
          handle,
          cellControl,
          cellModel(cell),
          data,
          trainRows,
          testRows,
          family
        )
        currentTrees <- spec$n.trees[cells$iTrees[cell]]
        numBurnIn <- spec$n.burn[1L]
      } else {
        bartcoreSetModel(sampler, cellModel(cell), cellControl, data)
        numBurnIn <- spec$n.burn[2L]
      }

      samples <- bartcoreRun(sampler, numBurnIn, spec$n.samples)
      lossValue <- spec$lossFunction(y.test, samples$test, weights.test)

      if (
        !is.numeric(lossValue) || length(lossValue) == 0L || anyNA(lossValue)
      ) {
        stop("loss function must return non-missing numeric values")
      }
      if (is.null(numLossResults)) {
        numLossResults <<- length(lossValue)
      } else if (length(lossValue) != numLossResults) {
        stop("loss function must always return the same number of values")
      }

      if (is.null(lossValues)) {
        lossValues <- matrix(0.0, numCells, numLossResults)
      }
      lossValues[cell, ] <- as.numeric(lossValue)
    }
    lossValues
  }

  results <- vector("list", length(repIndices))
  for (i in seq_along(repIndices)) {
    if (spec$method == "k-fold") {
      permutation <- sample.int(numObservations)
      repLoss <- NULL
      foldOffset <- 0L
      for (fold in seq_along(spec$foldSizes)) {
        testRows <- sort(permutation[
          foldOffset + seq_len(spec$foldSizes[fold])
        ])
        foldOffset <- foldOffset + spec$foldSizes[fold]
        foldLoss <- sweepCells(testRows)
        repLoss <- if (is.null(repLoss)) foldLoss else repLoss + foldLoss
      }
      results[[i]] <- repLoss / length(spec$foldSizes)
    } else {
      testRows <- sort(sample.int(numObservations, spec$numTest))
      results[[i]] <- sweepCells(testRows)
    }
  }

  do.call(rbind, results)
}
