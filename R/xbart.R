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
  n.burn = c(200L, 150L),
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
  sigest = NA_real_,
  seed = NA_integer_,
  factors = c("categorical", "indicators"),
  family = c("auto", "gaussian", "probit", "logistic"),
  missing = c("incorporate", "error"),
  node.prior = NULL,
  n.cuts = 100L,
  useQuantiles = FALSE,
  n.thin = 1L,
  storage = c("double", "single"),
  tree.prior = NULL
) {
  matchedCall <- match.call()

  currEnv <- sys.frame(sys.nframe())
  evalEnv <- parent.frame(1L)

  # the four flat knobs: shape-checked at the surface, mirroring
  # dbartsControl's own validity messages, then folded into the control
  # xbart builds for itself below
  n.cuts <- coerceOrError(n.cuts, "integer")
  if (is.na(n.cuts) || n.cuts <= 0L) {
    stop("'n.cuts' must be a positive integer")
  }
  useQuantiles <- coerceOrError(useQuantiles, "logical")
  if (is.na(useQuantiles)) {
    stop("'useQuantiles' must be TRUE/FALSE")
  }
  n.thin <- coerceOrError(n.thin, "integer")
  if (is.na(n.thin) || n.thin <= 0L) {
    stop("'n.thin' must be a positive integer")
  }
  storage <- match.arg(storage)

  # control = is no longer a formal: xbart builds its own control, the
  # knobs above landing where the flat formals' values landed, alongside
  # six fields xbart has always forced regardless of what a caller-supplied
  # control carried. n.trees/n.burn are xbart's own
  # grid axes (cellControl/cellModel overwrite them per cell below), seed is
  # handled by xbart's own RNG block, and printEvery/printCutoffs are inert
  # under verbose = FALSE, so control's copies of those five fields are left
  # at dbartsControl()'s own defaults.
  control <- dbartsControl(
    useQuantiles = useQuantiles,
    storage = storage,
    n.cuts = n.cuts,
    n.thin = n.thin,
    n.chains = 1L,
    n.threads = 1L,
    keepTrees = FALSE,
    keepTrainingFits = FALSE,
    updateState = FALSE,
    verbose = FALSE
  )

  validateCall <- redirectCall(
    matchedCall,
    quoteInNamespace(validateArgumentsInEnvironment),
    verbose,
    n.samples,
    sigest
  )
  validateCall <- addCallArgument(validateCall, 1L, currEnv)
  validateCall <- addCallArgument(validateCall, 2L, xbart)
  validateCall <- addCallArgument(validateCall, 3L, "xbart")
  validateCall <- addCallArgument(validateCall, "control", control)
  eval(validateCall, evalEnv, getNamespace("dbarts"))

  # the shared validator admits n.samples = 0 - a sampler is free to run
  # without keeping draws - but a cell scored on an empty sample matrix
  # would fault deeper, inside the loss
  if (control@n.samples <= 0L) {
    refuseZeroSamples("xbart")
  }

  if (control@call != call("NA")[[1L]]) {
    control@call <- matchedCall
  }

  # named ahead of the data build, matching bart2()/dbarts()/rbart_vi(), so
  # a bad family is refused before the response is ingested rather than after
  family <- match.arg(family)

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
  # a count-matrix data object declares the multinomial model, whose fitted
  # quantity is K probabilities per observation; every loss this function
  # evaluates is written against one location, so the fit would be scored as
  # though the slab were n rows. Refused ahead of the family resolution below,
  # which resolves counts to multinomial from "auto"
  refuseCountsCarryingData(formula, "xbart()")
  refuseResponseFreeFormula(formula, "xbart()")
  data <- eval(dataCall, evalEnv)
  data@n.cuts <- rep_len(control@n.cuts, ncol(data@x))
  data@sigma <- sigest

  # a factor/logical/character response is a classification; xbart cross-
  # validates the 2-level (probit) case only, never multinomial. A numeric
  # response takes the historic 0/1-vs-continuous path unchanged.
  family <- resolveClassificationFamily(data, family, "xbart", "gaussian")
  if (data@response.type == "numeric") {
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
  }
  control@binary <- isBinaryFamily(family)

  # the shared weight policy (R/spec.R's enforceWeightPolicy): a probit has
  # no tractable weighted latent-variable form and is refused (an all-ones
  # courtesy excepted), a logistic model requires positive integer count
  # weights, and a gaussian fit is unrestricted. xbart's own family is always
  # gaussian/probit/logistic, so the function's ordinal/nbinom branches never
  # fire here - the same function every other entry point reaches them with.
  data <- enforceWeightPolicy(data, family)

  if (is.na(data@sigma) && !control@binary) {
    data@sigma <- estimateStartingSigma(data)
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

  # a custom control could once supply an n.trees a caller left unnamed;
  # control = is gone, so this grid always comes from the formal itself
  n.trees <- coerceOrError(n.trees, "integer")
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

  # the grid default is a fixed k for every family, binary included: a
  # hyperprior k is held rather than swept and is drawn every sweep, so
  # defaulting binary fits to bart2's chi hyperprior would collapse the k
  # axis onto a single cell whose shrinkage moves within the fit
  if (is.null(matchedCall[["k"]])) {
    k <- if (!is.null(node.spec) && !is.null(node.spec@k)) {
      node.spec@k
    } else {
      2
    }
  }
  kIsGrid <- is.numeric(k)
  if (kIsGrid && (anyNA(k) || any(k <= 0))) {
    stop("'k' must contain only positive values")
  }
  # swept largest (most-shrunk) k first, so every warm start comes from a
  # simpler forest than the cell before it; kOrder un-permutes the reported
  # k axis back to the caller's order once the result array is final
  kOrder <- if (kIsGrid) order(k, decreasing = TRUE) else NULL
  if (kIsGrid) {
    k <- k[kOrder]
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

  # tree.prior (3.f, f4): follows the same grid-axis-overrides-the-object
  # rule as node.prior/k above - power and base are xbart's grid axes, so
  # cellModel overwrites them on the object every cell regardless of what is
  # supplied here, while the object's non-grid content (a cgm's split.probs,
  # a dart's Dirichlet hyperparameters) rides every cell unchanged. dart/
  # split.probs would only duplicate what a supplied tree.prior already
  # specifies, so they collide with it; power/base/k do not, since they are
  # grid axes rather than duplicates - this deliberately differs from
  # bart2's tree.prior, which does collide with power/base (R/bart.R's
  # buildSamplerPriors), because there they are ordinary scalars, not a grid.
  if (!is.null(matchedCall[["tree.prior"]])) {
    refuseColliding(matchedCall, "tree.prior", c("dart", "split.probs"))
    priorEnv <- new.env(parent = evalEnv)
    priorEnv[["cgm"]] <- getNamespace("dbarts")[["cgm"]]
    priorEnv[["dart"]] <- getNamespace("dbarts")[["dart"]]
    tree.prior <- eval(matchedCall[["tree.prior"]], priorEnv)
    if (is.function(tree.prior)) {
      tree.prior <- tree.prior()
    }
    if (!is(tree.prior, "dbartsTreePrior")) {
      stop(
        "'tree.prior' must be a tree prior specification; see ?dbartsPriors"
      )
    }
  } else {
    tree.prior <- resolveDartShorthand(
      dart,
      "split.probs" %in% names(matchedCall),
      "split.probs",
      function() dbartsPriors$dart(power[1L], base[1L]),
      function() cgm(power[1L], base[1L], split.probs)
    )
  }
  tree.prior <- resolveSplitProbabilities(tree.prior, data)

  if (is.null(node.spec)) {
    node.prior <- quote(normal(k))
    node.prior[[1L]] <- quoteInNamespace(normal)
    # take the first grid element: k[1L] for a numeric grid, k[[1L]] for a
    # list grid, otherwise k itself (a scalar or hyperprior spec)
    node.prior[[2L]] <- ifelse_3(is.numeric(k), is.list(k), k[1L], k[[1L]], k)
    node.prior <- eval(node.prior)
  } else {
    # first grid element, as above: k[1L] numeric, k[[1L]] list, else k
    kValue <- ifelse_3(is.numeric(k), is.list(k), k[1L], k[[1L]], k)
    if (is.call(kValue)) {
      kValue <- eval(kValue)
    }
    # the k argument replaces the supplied prior's own k, but its named
    # calibration is not a grid axis and rides every cell unchanged
    namedSd <- node.spec@prior.sd
    namedScale <- node.spec@prior.scale
    node.prior <- if (is(node.spec, "dbartsLinearPrior")) {
      resolveLeafCovariates(
        linear(node.spec@columns, kValue, namedSd, namedScale),
        data
      )
    } else if (is(node.spec, "dbartsGPPrior")) {
      resolveLeafCovariates(
        gp(
          node.spec@columns,
          kValue,
          node.spec@lengthscale,
          node.spec@max.leaf.size,
          namedSd,
          namedScale
        ),
        data
      )
    } else {
      normal(kValue, namedSd, namedScale)
    }
  }
  # xbart cells run a fixed k unless a hyperprior is named explicitly: the
  # grid default is 2 for every family, so the binary family default is never
  # taken here
  node.hyperprior <- resolveNodeHyperprior(node.prior@k, binary = FALSE)

  # a binary family runs on a fixed unit latent scale (R/spec.R's
  # fixedUnitScale rule): any supplied resid.prior is overridden, not just a
  # missing one, matching the shared resolver, so a caller cannot silently
  # fit an unfixed residual scale under a family that has none. The DEFAULT
  # value stays chisq rather than bart2's NULL-triggers-shorthand sentinel -
  # xbart has no sigdf/sigquant shorthands for a NULL to build from.
  resid.prior <- if (control@binary) {
    fixed(1)
  } else if (
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
    # a named calibration is held across every cell, created or re-modelled:
    # cellModel carries this model, and the setModel branch re-derives it
    prior.scale = resolvePriorScale(node.prior, node.hyperprior),
    node.scale = defaultNodeScale(family)
  )

  numObservations <- length(data@y)
  if (method == "k-fold") {
    n.test <- coerceOrError(n.test, "integer")
    if (n.test < 2L || n.test > numObservations) {
      stop(
        "for k-fold crossvalidation, 'n.test' must be an integer in [2, ",
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
      stop(
        "for random subsample crossvalidation, 'n.test' must be in (0, 1)"
      )
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
  n.burn <- rep_len(coerceOrError(n.burn, "integer"), 2L)
  if (anyNA(n.burn) || any(n.burn < 0L)) {
    stop("'n.burn' must contain non-negative integers")
  }
  n.threads <- coerceOrError(n.threads, "integer")
  if (length(n.threads) != 1L) {
    stop("'n.threads' must be of length 1")
  }
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
  # Restarting each split from a fresh forest keeps a slow-mixing cell from
  # remembering the previous fold and scoring optimistically on its own
  # held-out rows. Tree counts are fixed at a sampler's creation, so they
  # vary slowest and each count gets a fresh fit per split.
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
  chunkSeeds <- if (!is.na(seed)) {
    withFixedSeed(seed, sample.int(.Machine$integer.max, numChunks))
  } else {
    sample.int(.Machine$integer.max, numChunks)
  }

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

  # axis 3 is k, still in the decreasing sweep order; restore the caller's
  # order on both the array and the k vector before anything is reported
  if (kIsGrid && length(k) > 1L) {
    kOrderInv <- kOrder
    kOrderInv[kOrder] <- seq_along(kOrder)
    result <- result[,, kOrderInv, , , , drop = FALSE]
    k <- k[kOrderInv]
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
## held-out rows the previous training set contained.
## Returns a (reps x cells) x numResults matrix, cells in spec$cells order.
xbartRunChunk <- function(spec, repIndices, chunkSeed) {
  set.seed(chunkSeed)

  data <- spec$data
  cells <- spec$cells
  numCells <- nrow(cells)
  numObservations <- length(data@y)
  hasWeights <- !is.null(data@weights)
  family <- spec$model@family

  # linear and gp node priors read raw covariate values, fixed across cells;
  # the handle must own raw for them so each fold view can gather them
  nodePrior <- spec$model@node.prior
  leafCovariateColumns <-
    if (is(nodePrior, "dbartsLinearPrior") || is(nodePrior, "dbartsGPPrior")) {
      nodePrior@columns
    } else {
      NULL
    }
  handle <- bartcoreDataHandle(spec$control, data, leafCovariateColumns)

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
  # offset from offset[testRows], so each fold trains and scores on exactly
  # its own rows.
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
        bartcoreSetModel(sampler, cellModel(cell), data)
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
