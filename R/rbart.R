rbart.priors <- list(
  cauchy = function(x, rel.scale) dcauchy(x, 0, rel.scale * 2.5, TRUE),
  gamma = function(x, rel.scale) {
    dgamma(x, shape = 2.5, scale = rel.scale * 2.5, log = TRUE)
  }
)
cauchy <- NULL ## for R CMD check

## A seeded fit leaves the caller's stream where it found it, absent if it
## was absent.
readGlobalSeed <- function() {
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
}

writeGlobalSeed <- function(seed) {
  if (is.null(seed)) {
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  } else {
    assign(".Random.seed", seed, envir = .GlobalEnv)
  }
  invisible(NULL)
}

## The column of 'frame' named by 'name', or NULL when there is none; a name
## that is not a column must fall through to the caller's scope.
rbartColumn <- function(frame, name) {
  name <- as.character(name)
  if (!is.null(names(frame)) && name %in% names(frame)) {
    frame[[name]]
  } else {
    NULL
  }
}

rbart_vi <- function(
  formula,
  data,
  test,
  subset,
  weights,
  offset,
  offset.test = offset,
  group.by,
  group.by.test,
  prior = cauchy, ## can be a symbol in rbart.priors or a function; on log scale
  sigest = NA_real_,
  sigdf = 3.0,
  sigquant = 0.90,
  k = 2.0,
  power = 2.0,
  base = 0.95,
  n.trees = 75L,
  n.samples = 1500L,
  n.burn = 1500L,
  n.chains = 4L,
  n.threads = min(dbarts::guessNumCores(), n.chains),
  combineChains = FALSE,
  n.cuts = 100L,
  useQuantiles = FALSE,
  n.thin = 5L,
  keepTrainingFits = TRUE,
  printEvery = 100L,
  printCutoffs = 0L,
  verbose = TRUE,
  keepTrees = TRUE,
  keepCall = TRUE,
  seed = NA_integer_,
  keepSampler = keepTrees,
  keepTestFits = TRUE,
  callback = NULL,
  ...
) {
  matchedCall <- match.call()
  callingEnv <- parent.frame()
  warnOnce(
    "tombstone.rbart_vi",
    "'rbart_vi' is deprecated and is removed in dbarts ",
    tombstoneExpiry,
    "; grouped random effects live in stan4bart (stan4bart::stan4bart), whose ",
    "group-spread prior differs, so results move."
  )

  # the argument list is 0.9-34's, NA defaults included; NULL is the absent
  # spelling everywhere else, so it reads as NA here
  if (is.null(seed)) {
    seed <- NA_integer_
  }
  if (is.null(sigest)) {
    sigest <- NA_real_
  }

  # because we use a lot of trickery to redirect calls in the calling environment
  # (for example, to get the data), we replicate some base mechanisms like complaining
  # about unknown arguments
  argNames <- names(matchedCall)[-1L]
  unknownArgs <- argNames %not_in%
    names(formals(rbart_vi)) &
    argNames %not_in% names(formals(dbartsControl))
  if (any(unknownArgs)) {
    stop(
      "unknown arguments: '",
      paste0(argNames[unknownArgs], collapse = "', '"),
      "'"
    )
  }

  n.chains <- coerceOrError(n.chains, "integer")[1L]
  if (is.na(n.chains) || n.chains < 1L) {
    stop("n.chains must be a non-negative integer")
  }

  n.threads <- coerceOrError(n.threads, "integer")[1L]
  if (is.na(n.threads) || n.threads < 1L) {
    stop("n.threads must be a non-negative integer")
  }

  controlCall <- redirectCall(matchedCall, dbarts::dbartsControl)
  controlCall$seed <- NULL
  missingDefaults <- names(formals(rbart_vi))[
    names(formals(rbart_vi)) %in% names(formals(dbartsControl))
  ]
  missingDefaults <- missingDefaults[
    missingDefaults %not_in% c(names(controlCall), "...", "seed")
  ]
  controlCall[missingDefaults] <- formals(rbart_vi)[missingDefaults]
  if ("n.threads" %in% missingDefaults) {
    controlCall[["n.threads"]] <- eval(controlCall[["n.threads"]])
  }
  control <- eval(controlCall, envir = callingEnv)
  control@keepFits <- TRUE

  control@call <- if (keepCall) matchedCall else call("NULL")
  control@n.burn <- control@n.burn %/% control@n.thin
  control@n.samples <- control@n.samples %/% control@n.thin
  control@printEvery <- control@printEvery %/% control@n.thin
  if (control@n.samples == 0L) {
    stop("no posterior draws will be taken after thinning")
  }
  control@n.chains <- 1L
  control@n.threads <- max(control@n.threads %/% n.chains, 1L)
  if (n.chains > 1L && n.threads > 1L) {
    if (control@verbose) {
      warning("verbose output disabled for multiple threads")
    }
    control@verbose <- FALSE
  }

  keepSampler <- keepSampler || control@keepTrees

  tree.prior <- quote(cgm(power, base))
  tree.prior[[2L]] <- power
  tree.prior[[3L]] <- base

  if (!is.null(matchedCall[["k"]])) {
    leaf.prior <- quote(normal(k))
    leaf.prior[[2L]] <- k
  } else {
    leaf.prior <- NULL
  }

  family <- withResidPrior(
    newValidated("dbartsFamily", token = "auto"),
    chisq(sigdf, sigquant)
  )

  if (is.null(matchedCall[["group.by"]])) {
    stop("'group.by' must be specified to use rbart_vi")
  }

  group.by.literal <- NULL
  # look for group.by in data, if supplied, first
  if (is.symbol(matchedCall[["group.by"]]) && !missing(data)) {
    group.by.literal <- rbartColumn(data, matchedCall[["group.by"]])
  }

  if (is.null(group.by.literal)) {
    try(
      group.by.literal <- eval(matchedCall[["group.by"]], environment(formula)),
      silent = TRUE
    )
  }

  if (is.null(group.by.literal)) {
    try(group.by.literal <- group.by, silent = TRUE)
  }

  if (is.null(group.by.literal)) {
    stop("'group.by' not found")
  }
  group.by <- group.by.literal
  if (
    !is.numeric(group.by) && !is.factor(group.by) && !is.character(group.by)
  ) {
    stop("'group.by' must be coercible to factor type")
  }

  if (!is.null(matchedCall[["group.by.test"]])) {
    group.by.literal <- NULL
    if (is.symbol(matchedCall[["group.by.test"]]) && !missing(test)) {
      group.by.literal <- rbartColumn(test, matchedCall[["group.by.test"]])
    }

    if (is.null(group.by.literal)) {
      try(
        group.by.literal <- eval(
          matchedCall[["group.by.test"]],
          environment(formula)
        ),
        silent = TRUE
      )
    }

    if (is.null(group.by.literal)) {
      try(group.by.literal <- group.by.test, silent = TRUE)
    }

    if (
      is.null(group.by.literal) &&
        is.symbol(matchedCall[["group.by.test"]]) &&
        !missing(data)
    ) {
      group.by.literal <- rbartColumn(data, matchedCall[["group.by.test"]])
    }

    if (is.null(group.by.literal)) {
      stop("'group.by.test' not found")
    }

    group.by.test <- group.by.literal
    if (
      !is.numeric(group.by.test) &&
        !is.factor(group.by.test) &&
        !is.character(group.by.test)
    ) {
      stop("'group.by.test' must be coercible to factor type")
    }

    if (is.null(group.by.test)) {
      stop("'group.by.test' specified but not found")
    }

    if (
      !is.numeric(group.by.test) &&
        !is.factor(group.by.test) &&
        !is.character(group.by.test)
    ) {
      stop("'group.by.test' must be coercible to factor type")
    }
  }

  if (is.null(matchedCall$prior)) {
    matchedCall$prior <- formals(rbart_vi)$prior
  }

  if (
    is.symbol(matchedCall$prior) ||
      is.character(matchedCall$prior) &&
        any(names(rbart.priors) == matchedCall$prior)
  ) {
    prior <- rbart.priors[[which(names(rbart.priors) == matchedCall$prior)]]
  }

  dataCall <- redirectCall(matchedCall, dbarts::dbartsData)
  dataCall$factors <- "indicators"
  dataCall$na.action <- quote(stats::na.omit)
  data <- eval(dataCall, envir = callingEnv)

  if (!is.null(attr(data, "survivalStatus")) || !is.null(data@counts)) {
    stop("'rbart_vi' fits a continuous or binary response only")
  }
  if (
    all(data@y %in% c(0, 1)) &&
      !is.null(data@weights) &&
      !all(data@weights %in% c(0, 1))
  ) {
    stop(
      "'rbart_vi' takes only 0 and 1 weights with a binary response; ",
      "use stan4bart::stan4bart for weighted grouped fits"
    )
  }

  if (length(group.by) != length(data@y)) {
    stop(
      "'group.by' not of length equal to that of data; check for NAs in original data, and for name collisions with `data` argument and calling environment"
    )
  }
  group.by <- droplevels(as.factor(group.by))
  if (!is.null(matchedCall[["group.by.test"]])) {
    if (length(group.by.test) != nrow(data@x.test)) {
      stop("'group.by.test' not of length equal to that of data")
    }
    group.by.test <- droplevels(as.factor(group.by.test))
  } else if (!is.null(data@x.test)) {
    warning("'test' supplied by 'group.by.test' missing; recycling 'group.by'")
    group.by.test <- rep_len(group.by, nrow(data@x.test))
  } else {
    group.by.test <- NULL
  }

  if (!is.null(callback)) {
    if (!is.function(callback)) {
      stop("callback must be a function")
    }
    if (length(formals(callback)) != 5L) {
      stop("callback function must take exactly 5 arguments")
    }
  }

  rbartArgs <- namedList(
    group.by,
    prior,
    keepTrainingFits,
    keepTestFits,
    callback
  )

  samplerArgs <- namedList(
    formula = data,
    control,
    tree.prior,
    leaf.prior,
    family,
    sigest = if (is.na(sigest)) NULL else as.numeric(sigest)
  )
  if (is.null(leaf.prior)) {
    samplerArgs[["leaf.prior"]] <- NULL
  }
  if (!is.na(seed)) {
    oldSeed <- readGlobalSeed()
    on.exit(writeGlobalSeed(oldSeed), add = TRUE)
  }

  # any refusal of the chains' arguments surfaces here, once, rather than
  # inside a worker that then retries serially
  streamBefore <- readGlobalSeed()
  do.call(dbarts::dbarts, samplerArgs)
  writeGlobalSeed(streamBefore)

  chainResults <- vector("list", n.chains)
  runSingleThreaded <- n.threads <= 1L || n.chains <= 1L
  if (!runSingleThreaded) {
    tryResult <- tryCatch(
      cluster <- makeCluster(min(n.threads, n.chains), "PSOCK"),
      error = function(e) e
    )
    if (inherits(tryResult, "error")) {
      tryResult <- tryCatch(
        cluster <- makeCluster(min(n.threads, n.chains), "FORK"),
        error = function(e) e
      )
    }

    if (inherits(tryResult, "error")) {
      warning(
        "unable to multithread, defaulting to single: ",
        tryResult$message
      )
      runSingleThreaded <- TRUE
    } else {
      if (!is.na(seed)) {
        # one seed per chain, drawn sequentially from the given seed
        set.seed(seed)
        randomSeeds <- sample.int(.Machine$integer.max, n.chains)
      } else {
        randomSeeds <- rep.int(NA_integer_, n.chains)
      }

      clusterExport(
        cluster,
        c("rbart_vi_fit", "rbart_vi_run"),
        asNamespace("dbarts")
      )
      clusterEvalQ(cluster, require(dbarts))

      tryResult <- tryCatch(
        chainResults <- clusterMap(
          cluster,
          "rbart_vi_fit",
          seq_len(n.chains),
          randomSeeds,
          MoreArgs = namedList(samplerArgs, rbartArgs)
        ),
        error = function(e) e
      )

      stopCluster(cluster)

      if (inherits(tryResult, "error")) {
        warning(
          "error running multithreaded, defaulting to single: ",
          tryResult$message
        )
        runSingleThreaded <- TRUE
      }
    }
  }

  if (runSingleThreaded) {
    if (!is.na(seed)) {
      # run serially, every chain draws from the global generator
      set.seed(seed)
    }

    for (chainNum in seq_len(n.chains)) {
      chainResults[[chainNum]] <- rbart_vi_fit(
        1L,
        NA_integer_,
        samplerArgs,
        rbartArgs
      )
    }
  }
  packageRbartResults(
    control,
    data,
    group.by,
    group.by.test,
    chainResults,
    combineChains,
    seed,
    keepSampler
  )
}

rbart_vi_run <- function(
  sampler,
  data,
  state,
  prior,
  verbose,
  n.samples,
  isWarmup,
  rbartArgs
) {
  control <- sampler$control

  numRanef <- data$numRanef
  g.sel <- data$g.sel
  g <- data$g
  offset.orig <- data$offset.orig

  kIsModeled <- data$kIsModeled
  posteriorClosure <- prior$posteriorClosure
  evalEnv <- prior$evalEnv

  numObservations <- length(sampler$data@y)
  numTestObservations <- NROW(sampler$data@x.test)

  samples <- list(tau = rep(NA_real_, n.samples))
  if (!control@binary) {
    samples$sigma <- rep(NA_real_, n.samples)
  }
  if (kIsModeled) {
    samples$k <- rep(NA_real_, n.samples)
  }
  if (!isWarmup) {
    samples$ranef <- matrix(NA_real_, numRanef, n.samples)
    samples$yhat.train <- matrix(
      NA_real_,
      if (rbartArgs$keepTrainingFits) numObservations else 0L,
      n.samples
    )
    samples$yhat.test <- matrix(
      NA_real_,
      if (rbartArgs$keepTestFits) numTestObservations else 0L,
      n.samples
    )
    samples$varcount <- matrix(NA_integer_, ncol(sampler$data@x), n.samples)
  }

  # order of update matters - need to store a ranef that goes with a prediction
  # or else when they're added together they won't be consistent with `predict`
  for (i in seq_len(n.samples)) {
    # update ranef
    # row variance is sigma^2 / w_i, so a group's precision is its weight
    # total over sigma^2 and its data term is the weighted residual total
    resid <- with(state, y.st - treeFit.train)
    if (!is.null(data$weights)) {
      resid <- resid * data$weights
    }
    post.var <- 1.0 / (data$w.g / state$sigma^2.0 + 1.0 / state$tau^2.0)
    post.mean <- post.var *
      vapply(g.sel, function(sel) sum(resid[sel]), 0) /
      state$sigma^2.0
    ranef <- rnorm(numRanef, post.mean, sqrt(post.var))
    ranef.vec <- ranef[g]

    # update BART params
    sampler$setOffset(
      ranef.vec + if (!is.null(offset.orig)) offset.orig else 0,
      isWarmup
    )
    dbarts_samples <- sampler$run(0L, 1L)
    state$treeFit.train <- as.vector(dbarts_samples$train) - ranef.vec
    if (control@binary) {
      sampler$getLatents(state$y.st)
    }
    state$sigma <- dbarts_samples$sigma[1L]

    # update sd of ranef
    evalEnv$b.sq <- sum(ranef^2.0)
    state$tau <- sliceSample(
      posteriorClosure,
      state$tau,
      control@n.thin,
      boundary = c(0.0, Inf)
    )[control@n.thin]

    .Call(C_dbarts_assignInPlace, samples$tau, i, state$tau)
    if (!is.null(samples$sigma)) {
      .Call(C_dbarts_assignInPlace, samples$sigma, i, state$sigma)
    }
    if (!is.null(samples$ranef)) {
      .Call(C_dbarts_assignInPlace, samples$ranef, i, ranef)
    }
    if (!is.null(samples$yhat.train) && rbartArgs$keepTrainingFits) {
      .Call(C_dbarts_assignInPlace, samples$yhat.train, i, state$treeFit.train)
    }
    if (!is.null(samples$varcount)) {
      .Call(
        C_dbarts_assignInPlace,
        samples$varcount,
        i,
        dbarts_samples$varcount
      )
    }
    if (
      !is.null(samples$yhat.test) &&
        numTestObservations > 0L &&
        rbartArgs$keepTestFits
    ) {
      .Call(C_dbarts_assignInPlace, samples$yhat.test, i, dbarts_samples$test)
    }
    if (!is.null(samples$k)) {
      .Call(C_dbarts_assignInPlace, samples$k, i, dbarts_samples$k)
    }
    if (!isWarmup && !is.null(rbartArgs$callback)) {
      names(ranef) <- data$g.levels
      if (is.null(samples$callback)) {
        callback_i <- rbartArgs$callback(
          state$treeFit.train,
          dbarts_samples$test,
          ranef,
          state$sigma,
          state$tau
        )
        samples$callback <- matrix(
          NA_real_,
          length(callback_i),
          control@n.samples,
          dimnames = list(names(callback_i), NULL)
        )
        .Call(C_dbarts_assignInPlace, samples$callback, i, callback_i)
        rm(callback_i)
      } else {
        .Call(
          C_dbarts_assignInPlace,
          samples$callback,
          i,
          rbartArgs$callback(
            state$treeFit.train,
            dbarts_samples$test,
            ranef,
            state$sigma,
            state$tau
          )
        )
      }
    }

    if (verbose && i %% control@printEvery == 0L) {
      cat("iter: ", i, "\n", sep = "")
    }
  }

  list(state = state, samples = samples)
}

rbart_vi_fit <- function(chain.num, seed, samplerArgs, rbartArgs) {
  chain.num <- "ignored"

  if (!is.na(seed)) {
    set.seed(seed)
  }

  sampler <- do.call(dbarts::dbarts, samplerArgs)
  sampler$control@call <- samplerArgs$control@call

  oldUpdateState <- sampler$control@updateState
  verbose <- sampler$control@verbose
  control <- sampler$control
  control@updateState <- FALSE
  control@verbose <- FALSE
  control@keepTrainingFits <- TRUE
  control@keepFits <- TRUE
  sampler$setControl(control)

  # the loop sets each sweep's offset, intercepts included, and the test rows
  # must not inherit it
  oldTestUsesRegularOffset <- sampler$data@testUsesRegularOffset
  sampler$data@testUsesRegularOffset <- FALSE

  y <- sampler$data@y
  rel.scale <- if (!control@binary) sd(y) else 0.5

  g <- as.integer(rbartArgs$group.by)
  g.levels <- levels(rbartArgs$group.by)
  numRanef <- nlevels(rbartArgs$group.by)
  g.sel <- lapply(seq_len(numRanef), function(j) g == j)
  n.g <- sapply(g.sel, sum)
  offset.orig <- sampler$data@offset
  # a binary fit keeps its 0/1 weights as the active rows
  weights <- sampler$data@weights
  if (is.null(weights)) {
    weights <- sampler$activeRows
  }
  w.g <- if (is.null(weights)) {
    n.g
  } else {
    vapply(g.sel, function(sel) sum(weights[sel]), 0)
  }
  kIsModeled <- as.logical(sampler$getLeafPrior()[1L, "k.has.hyperprior"])
  data <- namedList(
    w.g,
    weights,
    numRanef,
    g.sel,
    g,
    g.levels,
    offset.orig,
    kIsModeled
  )

  evalEnv <- list2env(list(
    rel.scale = rel.scale,
    q = numRanef,
    prior = rbartArgs$prior
  ))
  b.sq <- NULL ## for R CMD check
  posteriorClosure <- function(x) {
    ifelse(
      x <= 0.0 | is.infinite(x),
      -.Machine$double.xmax * .Machine$double.eps,
      -q * base::log(x) - 0.5 * b.sq / x^2.0 + prior(x, rel.scale)
    )
  }
  environment(posteriorClosure) <- evalEnv
  prior <- namedList(posteriorClosure, evalEnv)

  sampler$sampleTreesFromPrior()
  state <- list(
    tau = rel.scale / 5.0,
    sigma = if (!control@binary) sampler$data@sigma else 1.0,
    y.st = if (!control@binary) y else sampler$getLatents()
  )
  # Sample from prior to get started
  ranef <- rnorm(numRanef, 0.0, state$tau)
  ranef.vec <- ranef[g]

  prior <- list(
    posteriorClosure = posteriorClosure,
    evalEnv = evalEnv
  )

  if (control@n.burn > 0L) {
    oldKeepTrees <- control@keepTrees
    control@keepTrees <- FALSE
    sampler$setControl(control)

    sampler$setOffset(
      ranef.vec + if (!is.null(offset.orig)) offset.orig else 0,
      TRUE
    )

    state$treeFit.train <- as.vector(sampler$run(0L, 1L)$train) - ranef.vec

    run_result <- rbart_vi_run(
      sampler,
      data,
      state,
      prior,
      FALSE,
      control@n.burn,
      TRUE,
      rbartArgs
    )
    state <- run_result$state

    firstTau <- run_result$samples$tau
    firstSigma <- run_result$samples$sigma
    firstK <- run_result$samples$k

    if (control@keepTrees != oldKeepTrees) {
      control@keepTrees <- TRUE
      sampler$setControl(control)
    }
  } else {
    sampler$setOffset(
      ranef.vec + if (!is.null(offset.orig)) offset.orig else 0,
      TRUE
    )

    if (control@keepTrees) {
      control@keepTrees <- FALSE
      sampler$setControl(control)
      state$treeFit.train <- as.vector(sampler$run(0L, 1L)$train) - ranef.vec
      control@keepTrees <- TRUE
      sampler$setControl(control)
    } else {
      state$treeFit.train <- as.vector(sampler$run(0L, 1L)$train) - ranef.vec
    }

    firstTau <- NULL
    firstSigma <- NULL
    firstK <- NULL
  }

  run_result <- rbart_vi_run(
    sampler,
    data,
    state,
    prior,
    verbose,
    control@n.samples,
    FALSE,
    rbartArgs
  )

  tau <- run_result$samples$tau
  sigma <- run_result$samples$sigma
  ranef <- run_result$samples$ranef
  yhat.train <- run_result$samples$yhat.train
  yhat.test <- run_result$samples$yhat.test
  k <- run_result$samples$k
  callback <- run_result$samples$callback
  varcount <- run_result$samples$varcount

  sampler$data@testUsesRegularOffset <- oldTestUsesRegularOffset
  sampler$setOffset(if (!is.null(offset.orig)) offset.orig else NULL, FALSE)

  control@updateState <- oldUpdateState
  sampler$setControl(control)
  # a sampler that leaves a worker, or is saved, is re-created from this
  sampler$storeState()

  rownames(ranef) <- g.levels

  result <- namedList(
    sampler,
    ranef,
    firstTau,
    firstSigma,
    tau,
    sigma,
    yhat.train,
    yhat.test,
    callback,
    varcount
  )
  if (kIsModeled) {
    result$firstK <- firstK
    result$k <- k
  }
  result
}

packageRbartResults <- function(
  control,
  data,
  group.by,
  group.by.test,
  chainResults,
  combineChains,
  seed,
  keepSampler
) {
  n.chains <- length(chainResults)

  responseIsBinary <- chainResults[[1L]]$sampler$control@binary

  result <- list(call = control@call, y = data@y, group.by = group.by)
  if (!responseIsBinary) {
    result$sigest <- chainResults[[1L]]$sampler$data@sigma
  }
  if (!is.null(group.by.test)) {
    result[["group.by.test"]] <- group.by.test
  }

  if (n.chains > 1L) {
    if (
      !is.null(group.by.test) &&
        any(unmeasuredLevels <- levels(group.by.test) %not_in% levels(group.by))
    ) {
      warning(
        "test includes random effect levels not present in training - ranef estimates default to draws from the ranef distribution parameterized by the posterior of its variance"
      )
      n.samples <- dim(chainResults[[1L]]$ranef)[2L]
      n.unmeasured <- sum(unmeasuredLevels)
      totalRanef <- sapply(seq_along(chainResults), function(k) {
        unmeasuredRanef <- matrix(
          rnorm(
            n.unmeasured * n.samples,
            0,
            rep(chainResults[[k]]$tau, each = n.unmeasured)
          ),
          n.unmeasured,
          n.samples,
          dimnames = list(levels(group.by.test)[unmeasuredLevels], NULL)
        )
        rbind(chainResults[[k]]$ranef, unmeasuredRanef)
      })
      ranefDim <- c(
        dim(chainResults[[1L]]$ranef)[1L] + n.unmeasured,
        n.samples,
        n.chains
      )
      ranefDimnames <- list(
        c(
          rownames(chainResults[[1L]]$ranef),
          levels(group.by.test)[unmeasuredLevels]
        ),
        NULL,
        NULL
      )
      ranef <- array(totalRanef, ranefDim, ranefDimnames)
      result$ranef <- convertSamplesFromDbartsToBart(
        ranef,
        n.chains,
        combineChains
      )
    } else {
      ranef <- array(
        sapply(chainResults, function(x) x$ranef),
        c(dim(chainResults[[1L]]$ranef), n.chains),
        list(rownames(chainResults[[1L]]$ranef), NULL, NULL)
      )
      result$ranef <- convertSamplesFromDbartsToBart(
        ranef,
        n.chains,
        combineChains
      )
    }
    result$first.tau <- convertSamplesFromDbartsToBart(
      sapply(chainResults, function(x) x$firstTau),
      n.chains,
      combineChains
    )
    if (!responseIsBinary) {
      result$first.sigma <- convertSamplesFromDbartsToBart(
        sapply(chainResults, function(x) x$firstSigma),
        n.chains,
        combineChains
      )
      result$sigma <- convertSamplesFromDbartsToBart(
        sapply(chainResults, function(x) x$sigma),
        n.chains,
        combineChains
      )
    }
    result$tau <- convertSamplesFromDbartsToBart(
      sapply(chainResults, function(x) x$tau),
      n.chains,
      combineChains
    )
    if (NROW(chainResults[[1L]]$yhat.train) <= 0L) {
      result$yhat.train <- NULL
    } else {
      result$yhat.train <- convertSamplesFromDbartsToBart(
        array(
          sapply(chainResults, function(x) x$yhat.train),
          c(dim(chainResults[[1L]]$yhat.train), n.chains)
        ),
        n.chains,
        combineChains
      )
    }
    if (NROW(chainResults[[1L]]$yhat.test) <= 0L) {
      result$yhat.test <- NULL
    } else {
      result$yhat.test <- convertSamplesFromDbartsToBart(
        array(
          sapply(chainResults, function(x) x$yhat.test),
          c(dim(chainResults[[1L]]$yhat.test), n.chains)
        ),
        n.chains,
        combineChains
      )
    }
    if (!is.null(chainResults[[1L]]$callback)) {
      result$callback <- convertSamplesFromDbartsToBart(array(
        sapply(chainResults, function(x) x$callback),
        c(dim(chainResults[[1L]]$callback), n.chains)
      ))
      dimnames(result$callback) <- list(
        NULL,
        NULL,
        dimnames(chainResults[[1L]]$callback)[[1L]]
      )
    }
    result$varcount <- convertSamplesFromDbartsToBart(
      array(
        sapply(chainResults, function(x) x$varcount),
        c(dim(chainResults[[1L]]$varcount), n.chains)
      ),
      n.chains,
      combineChains
    )
    if (!is.null(chainResults[[1L]]$firstK)) {
      result$first.k <- convertSamplesFromDbartsToBart(
        sapply(chainResults, function(x) x$firstK),
        n.chains,
        combineChains
      )
    }
    if (!is.null(chainResults[[1L]]$k)) {
      result$k <- convertSamplesFromDbartsToBart(
        sapply(chainResults, function(x) x$k),
        n.chains,
        combineChains
      )
    }
  } else {
    result$ranef <- t(chainResults[[1L]]$ranef)
    if (
      !is.null(group.by.test) &&
        any(unmeasuredLevels <- levels(group.by.test) %not_in% levels(group.by))
    ) {
      warning(
        "test includes random effect levels not present in training - ranef estimates default to draws from the ranef distribution parameterized by the posterior of its variance"
      )
      n.unmeasured <- sum(unmeasuredLevels)
      n.samples <- ncol(chainResults[[1L]]$ranef)
      unmeasuredRanef <-
        matrix(
          rnorm(
            n.samples * n.unmeasured,
            0,
            rep(chainResults[[1L]]$tau, n.samples)
          ),
          n.samples,
          n.unmeasured,
          dimnames = list(NULL, levels(group.by.test)[unmeasuredLevels])
        )
      result$ranef <- cbind(result$ranef, unmeasuredRanef)
    }
    result$first.tau <- chainResults[[1L]]$firstTau
    if (!responseIsBinary) {
      result$first.sigma <- chainResults[[1L]]$firstSigma
      result$sigma <- chainResults[[1L]]$sigma
    }
    result$tau <- chainResults[[1L]]$tau
    result$yhat.train <- if (NROW(chainResults[[1L]]$yhat.train) <= 0L) {
      NULL
    } else {
      t(chainResults[[1L]]$yhat.train)
    }
    result$yhat.test <- if (NROW(chainResults[[1L]]$yhat.test) <= 0L) {
      NULL
    } else {
      t(chainResults[[1L]]$yhat.test)
    }
    if (!is.null(chainResults[[1L]]$callback)) {
      result$callback <- t(chainResults[[1L]]$callback)
    }
    result$varcount <- chainResults[[1L]]$varcount
    if (!is.null(chainResults[[1L]]$firstK)) {
      result$first.k <- chainResults[[1L]]$firstK
    }
    if (!is.null(chainResults[[1L]]$k)) {
      result$k <- chainResults[[1L]]$k
    }
  }

  result$ranef.mean <- apply(result$ranef, length(dim(result$ranef)), mean)
  if (control@keepTrainingFits) {
    result$yhat.train.mean <- apply(
      result$yhat.train,
      length(dim(result$yhat.train)),
      mean
    )
  }
  if (!is.null(result$yhat.test)) {
    result$yhat.test.mean <- apply(
      result$yhat.test,
      length(dim(result$yhat.test)),
      mean
    )
  }

  if (keepSampler) {
    result$fit <- lapply(chainResults, function(x) x$sampler)
  } else {
    result$n.chains <- n.chains
  }

  if (!is.na(seed)) {
    oldSeed <- readGlobalSeed()
    set.seed(seed)
    result$seed <- .GlobalEnv$.Random.seed
    writeGlobalSeed(oldSeed)
  } else {
    if (!exists(".Random.seed", .GlobalEnv)) {
      runif(1L)
    }
    result$seed <- .GlobalEnv$.Random.seed
  }

  class(result) <- "rbart"
  result
}


## An rbart fit saved by dbarts 0.9-x holds samplers that predate the fields
## this version re-creates them from, so they cannot be used.
refuseLegacyRbart <- function(object) {
  if (
    !is.null(object$fit) &&
      !exists(
        "activeRows",
        envir = as.environment(object$fit[[1L]]),
        inherits = FALSE
      )
  ) {
    stop(
      "this fit was saved by dbarts 0.9-x; dbarts 1.0-0 cannot read its ",
      "trees; refit with this version",
      call. = FALSE
    )
  }
  invisible(NULL)
}

predict.rbart <- function(
  object,
  newdata,
  group.by,
  offset,
  type = c("ev", "ppd", "bart", "ranef"),
  combineChains = TRUE,
  ...
) {
  if (is.null(object$fit)) {
    stop("predict requires rbart to be called with 'keepTrees' == TRUE")
  }
  refuseLegacyRbart(object)

  dotsList <- list(...)
  if (!is.null(dotsList[["value"]])) {
    warning("argument 'value' has been deprecated; use 'type' instead")
    type <- dotsList[["value"]]
    dotsList[["value"]] <- NULL
  }

  if (is.character(type)) {
    if (type[1L] == "response") {
      type[1L] <- "ev"
    } else if (type[1L] == "link") {
      type[1L] <- "bart"
    }
  }
  if (is.character(type) && length(type) > 0L && type[1L] == "post-mean") {
    warning("type of 'post-mean' for predict deprecated; use 'ev' instead")
    type[1L] <- "ev"
  }
  if (
    !is.character(type) ||
      length(type) == 0L ||
      type[1L] %not_in% eval(formals(predict.rbart)$type)
  ) {
    stop(
      "type must be in '",
      paste0(eval(formals(predict.rbart)$type), collapse = "', '"),
      "'"
    )
  }
  type <- type[1L]

  if (missing(offset)) {
    offset <- NULL
  }

  n.chains <- if (is.null(object$n.chains)) {
    length(object$fit)
  } else {
    object$n.chains
  }
  n.samples <- object$fit[[1L]]$control@n.samples

  nonParametricPart <- 0
  # collects results in an array of n.obs x n.samples x n.chains, default for
  # internal sampler
  #
  # utilize bart stuff to get n.obs, since we would otherwise have to build
  # the test matrix
  if (type != "ranef") {
    if (n.chains > 1L) {
      n.obs <- NULL
      nonParametricPart <- array(
        sapply(seq_len(n.chains), function(i) {
          res <- object$fit[[i]]$predict(newdata, offset)
          if (is.null(n.obs)) {
            n.obs <<- dim(res)[1L]
          }
          res
        }),
        c(n.obs, n.samples, n.chains)
      )
    } else {
      nonParametricPart <- object$fit[[1L]]$predict(newdata, offset)
      n.obs <- nrow(nonParametricPart)
    }
    if (n.obs != length(group.by)) {
      stop("length of group.by not equal to number of rows in test")
    }

    nonParametricPart <- convertSamplesFromDbartsToBart(
      nonParametricPart,
      n.chains,
      combineChains
    )
  }

  if (type == "bart") {
    return(nonParametricPart)
  }

  ranef <- 0
  if (type != "bart") {
    ranefNames.test <- levels(group.by)
    ranefNames.train <- if (length(dim(object$ranef)) > 2L) {
      dimnames(object$ranef)[[3L]]
    } else {
      dimnames(object$ranef)[[2L]]
    }

    ranef <- object$ranef
    if (n.chains > 1L) {
      if (length(dim(ranef)) > 2L && combineChains) {
        ranef <- combineChains(ranef)
      } else if (length(dim(ranef)) == 2L && !combineChains && n.chains > 1L) {
        ranef <- uncombineChains(ranef, n.chains)
      }
    }

    if (!all(measuredLevels <- ranefNames.test %in% ranefNames.train)) {
      warning(
        "test includes random effect levels not present in training - ranef estimates default to draws from their latent distribution parameterized by the posterior of its variance; draws may not be the same across future calls to 'predict'"
      )
      n.unmeasured <- sum(!measuredLevels)
      if (n.chains > 1L) {
        # the draws' scales in the layout of the intercepts they fill: chain
        # fastest when split, chain-major when combined
        tauSplit <- if (is.null(dim(object$tau))) {
          uncombineChains(object$tau, n.chains)
        } else {
          object$tau
        }
        tauScale <- if (combineChains) {
          as.vector(t(tauSplit))
        } else {
          as.vector(tauSplit)
        }
        if (!combineChains) {
          unmeasuredRanef <- array(
            rnorm(
              n.chains * n.samples * n.unmeasured,
              0,
              rep.int(tauScale, n.unmeasured)
            ),
            c(n.chains, n.samples, n.unmeasured),
            dimnames = list(NULL, NULL, ranefNames.test[!measuredLevels])
          )
        } else {
          unmeasuredRanef <- matrix(
            rnorm(
              n.chains * n.samples * n.unmeasured,
              0,
              rep.int(tauScale, n.unmeasured)
            ),
            n.chains * n.samples,
            n.unmeasured,
            dimnames = list(NULL, ranefNames.test[!measuredLevels])
          )
        }
        if (length(dim(ranef)) == 2L) {
          ranef <- cbind(ranef, unmeasuredRanef)
        } else {
          # ranef are n.chains x n.samples x n.group
          ranef <- array(
            c(ranef, unmeasuredRanef),
            c(n.chains, n.samples, dim(ranef)[3L] + n.unmeasured),

            dimnames = list(
              NULL,
              NULL,
              c(dimnames(ranef)[[3L]], dimnames(unmeasuredRanef)[[3L]])
            )
          )
        }
      } else {
        unmeasuredRanef <- matrix(
          rnorm(n.samples * n.unmeasured, 0, rep.int(object$tau, n.unmeasured)),
          n.samples,
          n.unmeasured,
          dimnames = list(NULL, ranefNames.test[!measuredLevels])
        )
        ranef <- cbind(ranef, unmeasuredRanef)
      }
    }
  }

  if (type == "ranef") {
    ranef <- if (length(dim(ranef)) > 2L) {
      ranef[,, ranefNames.test, drop = FALSE]
    } else {
      ranef[, ranefNames.test, drop = FALSE]
    }
    return(rbartCombineOrUncombineChains(ranef, n.chains, combineChains))
  }

  ranef <- unname(
    if (length(dim(ranef)) > 2L) {
      ranef[,, as.character(group.by), drop = FALSE]
    } else {
      ranef[, as.character(group.by), drop = FALSE]
    }
  )
  ranef <- rbartCombineOrUncombineChains(ranef, n.chains, combineChains)

  if (
    length(dim(nonParametricPart)) != length(dim(ranef)) ||
      any(dim(nonParametricPart) != dim(ranef))
  ) {
    stop("internal error: fixed and random parts do not conform")
  }
  result <- nonParametricPart + ranef

  responseIsBinary <- is.null(object[["sigma"]])
  if (responseIsBinary) {
    result <- pnorm(result)
  }

  if (type == "ppd") {
    result <- sampleFromPPD(result, object, NULL, n.chains)
  }

  if (exists("unmeasuredRanef", inherits = FALSE)) {
    attr(result, "ranef") <- unmeasuredRanef
  }

  result
}

extract.rbart <- function(
  object,
  type = c("ev", "ppd", "bart", "ranef", "trees"),
  sample = c("train", "test"),
  combineChains = TRUE,
  ...
) {
  if (is.character(type)) {
    if (type[1L] == "response") {
      type[1L] <- "ev"
    } else if (type[1L] == "link") {
      type[1L] <- "bart"
    }
  }
  if (
    !is.character(type) || type[1L] %not_in% eval(formals(extract.rbart)$type)
  ) {
    stop(
      "type must be in '",
      paste0(eval(formals(extract.rbart)$type), collapse = "', '"),
      "'"
    )
  }
  type <- type[1L]

  n.chains <- if (is.null(object$n.chains)) {
    length(object$fit)
  } else {
    object$n.chains
  }

  if (type == "trees") {
    if (is.null(object$fit)) {
      stop(
        "extracting trees requires rbart to be called with 'keepTrees' == TRUE"
      )
    }
    refuseLegacyRbart(object)
    treesCall <- match.call()
    target <- quote(object$fit[[i]]$getTrees)
    target[[2L]][[2L]][[2L]] <- treesCall$object
    treesCall[[1L]] <- target
    treesCall$object <- NULL
    treesCall$type <- NULL
    treesCall$chainNums <- NULL
    evalEnv <- parent.frame()
    dotsList <- list(...)
    chainNums <- if ("chainNums" %in% names(dotsList)) {
      as.integer(dotsList[["chainNums"]])
    } else {
      seq_len(n.chains)
    }
    varOrder <- c("sample", "chain", "tree", "n", "var", "value")
    allTrees <- lapply(chainNums, function(i) {
      result_i <- eval(subTermInLanguage(treesCall, quote(i), i), evalEnv)
      if (n.chains > 1L) {
        result_i$chain <- i
      }
      result_i[, varOrder[varOrder %in% colnames(result_i)]]
    })
    if (length(allTrees) > 1L) {
      allTrees <- Reduce(rbind, allTrees)
    } else {
      allTrees <- allTrees[[1L]]
    }
    row.names(allTrees) <- as.character(seq_len(nrow(allTrees)))
    return(allTrees)
  }

  if (
    !is.character(sample) ||
      sample[1L] %not_in% eval(formals(extract.rbart)$sample)
  ) {
    stop(
      "sample must be in '",
      paste0(eval(formals(extract.rbart)$sample), collapse = "', '"),
      "'"
    )
  }
  sample <- sample[1L]

  if (sample == "test" && is.null(object[["yhat.test"]])) {
    stop(
      "cannot extract test sample predictions if no test data exists; use `predict` instead"
    )
  }

  if (type == "ranef") {
    ranefNames <- if (sample == "train") {
      levels(object$group.by)
    } else {
      levels(object$group.by.test)
    }
    ranef <- if (length(dim(object$ranef)) > 2L) {
      object$ranef[,, ranefNames, drop = FALSE]
    } else {
      object$ranef[, ranefNames, drop = FALSE]
    }
    if (n.chains > 1L) {
      if (length(dim(ranef)) > 2L && combineChains) {
        ranef <- combineChains(ranef)
      } else if (length(dim(ranef)) == 2L && !combineChains && n.chains > 1L) {
        ranef <- uncombineChains(ranef, n.chains)
      }
    }

    return(ranef)
  }

  result <- if (sample == "train") object$yhat.train else object$yhat.test
  # if necessary, recover chain information or throw it away
  if (n.chains > 1L) {
    if (length(dim(result)) > 2L && combineChains) {
      result <- combineChains(result)
    } else if (length(dim(result)) == 2L && !combineChains && n.chains > 1L) {
      result <- uncombineChains(result, n.chains)
    }
  }

  if (type == "bart") {
    return(result)
  }

  ranefNames <- if (sample == "train") {
    as.character(object$group.by)
  } else {
    as.character(object$group.by.test)
  }
  ranef <- unname(
    if (length(dim(object$ranef)) > 2L) {
      object$ranef[,, ranefNames, drop = FALSE]
    } else {
      object$ranef[, ranefNames, drop = FALSE]
    }
  )

  if (n.chains > 1L) {
    if (length(dim(ranef)) > 2L && combineChains) {
      ranef <- combineChains(ranef)
    } else if (length(dim(ranef)) == 2L && !combineChains && n.chains > 1L) {
      ranef <- uncombineChains(ranef, n.chains)
    }
  }

  result <- result + ranef

  responseIsBinary <- is.null(object[["sigma"]])
  if (responseIsBinary) {
    result <- pnorm(result)
  }

  if (type == "ppd") {
    result <- sampleFromPPD(result, object, NULL, n.chains)
  }

  result
}

fitted.rbart <- function(
  object,
  type = c("ev", "ppd", "bart", "ranef"),
  sample = c("train", "test"),
  ...
) {
  if (is.character(type)) {
    if (type[1L] == "response") {
      type[1L] <- "ev"
    } else if (type[1L] == "link") {
      type[1L] <- "bart"
    }
  }
  if (
    !is.character(type) || type[1L] %not_in% eval(formals(fitted.rbart)$type)
  ) {
    stop(
      "type must be in '",
      paste0(eval(formals(fitted.rbart)$type), collapse = "', '"),
      "'"
    )
  }
  type <- type[1L]

  if (
    !is.character(sample) ||
      sample[1L] %not_in% eval(formals(fitted.rbart)$sample)
  ) {
    stop(
      "sample must be in '",
      paste0(eval(formals(fitted.rbart)$sample), collapse = "', '"),
      "'"
    )
  }
  sample <- sample[1L]

  if (type == "ev") {
    ranefNames <- dimnames(object$ranef)
    ranefNames <- ranefNames[[length(ranefNames)]]
    if (sample == "train") {
      groupByMatch <- match(object$group.by, ranefNames)
      result <- rbartFittedMean(
        object$yhat.train,
        object$ranef,
        groupByMatch,
        is.null(object[["sigma"]])
      )
    } else {
      groupByMatch <- match(object$group.by.test, ranefNames)
      result <- rbartFittedMean(
        object$yhat.test,
        object$ranef,
        groupByMatch,
        is.null(object[["sigma"]])
      )
    }
  } else {
    result <- extract(object, type, sample, ...)

    result <- if (!is.null(dim(result))) {
      apply(result, length(dim(result)), mean)
    } else {
      mean(result)
    }
  }

  result
}

residuals.rbart <- function(object, ...) {
  object$y - fitted.rbart(object)
}
print.rbart <- function(x, ...) {
  cat(
    "\nCall:\n",
    paste(deparse(x$call), sep = "\n", collapse = "\n"),
    "\n\n",
    sep = ""
  )
  invisible(x)
}


rbartCombineOrUncombineChains <- function(x, n.chains, combineChains) {
  if (n.chains > 1L) {
    if (length(dim(x)) > 2L && combineChains) {
      x <- combineChains(x)
    } else if (length(dim(x)) == 2L && !combineChains) {
      x <- uncombineChains(x, n.chains)
    }
  }
  x
}

rbartFittedMean <- function(yhat, ranef, groupByMatch, responseIsBinary) {
  nd <- length(dim(yhat))
  n <- dim(yhat)[nd]
  yhat <- matrix(yhat, ncol = n)
  ranef <- matrix(ranef, ncol = dim(ranef)[length(dim(ranef))])
  eta <- yhat + ranef[, groupByMatch, drop = FALSE]
  if (responseIsBinary) {
    eta <- pnorm(eta)
  }
  colMeans(eta)
}

plot.rbart <- function(
  x,
  plquants = c(0.05, 0.95),
  cols = c("blue", "black"),
  ...
) {
  if ("sigma" %in% names(x)) {
    par(mfrow = c(1L, 2L))
    if (!is.null(dim(x$sigma))) {
      plot(
        NULL,
        type = "n",
        ylab = "sigma",
        xlim = c(1, ncol(x$first.sigma) + ncol(x$sigma)),
        ylim = range(x$first.sigma, x$sigma)
      )
      for (i in seq_len(nrow(x$sigma))) {
        lines(
          c(seq_len(ncol(x$first.sigma)), ncol(x$first.sigma) + 0.5),
          c(
            x$first.sigma[i, ],
            0.5 * (x$first.sigma[i, ncol(x$first.sigma)] + x$sigma[i, 1L])
          ),
          col = "red",
          lty = i
        )
        lines(
          c(
            ncol(x$first.sigma) + 0.5,
            seq.int(ncol(x$first.sigma) + 1, length.out = ncol(x$sigma))
          ),
          c(
            0.5 * (x$first.sigma[i, ncol(x$first.sigma)] + x$sigma[i, 1L]),
            x$sigma[i, ]
          ),
          lty = i
        )
      }
    } else {
      plot(
        c(x$first.sigma, x$sigma),
        col = rep(c("red", "black"), c(length(x$first.sigma), length(x$sigma))),
        ylab = "sigma",
        ...
      )
    }
  }

  if (length(dim(x$ranef)) > 2L) {
    ranef <- x$ranef[,, as.integer(x$group.by)]
  } else {
    ranef <- x$ranef[, as.integer(x$group.by)]
  }
  yhat.train <- x$yhat.train + ranef

  if ("sigma" %in% names(x)) {
    ql <- apply(
      yhat.train,
      length(dim(yhat.train)),
      quantile,
      probs = plquants[1L]
    )
    qm <- apply(yhat.train, length(dim(yhat.train)), quantile, probs = .5)
    qu <- apply(
      yhat.train,
      length(dim(yhat.train)),
      quantile,
      probs = plquants[2L]
    )
    plot(
      x$y,
      qm,
      ylim = range(ql, qu),
      xlab = "y",
      ylab = "posterior interval for E(Y | x)",
      ...
    )

    for (i in seq_along(qm)) {
      lines(rep(x$y[i], 2L), c(ql[i], qu[i]), col = cols[1L])
    }
    abline(0, 1, lty = 2L, col = cols[2L])
  } else {
    ## shouldn't happen for now
    pdrs <- pnorm(yhat.train) #draws of p(Y=1 | x)
    ql <- apply(pdrs, length(dim(pdrs)), quantile, probs = plquants[1L])
    qm <- apply(pdrs, length(dim(pdrs)), quantile, probs = .5)
    qu <- apply(pdrs, length(dim(pdrs)), quantile, probs = plquants[2L])
    plot(
      qm,
      qm,
      ylim = range(ql, qu),
      xlab = "median of p",
      ylab = "posterior interval for P(Y = 1|  x)",
      ...
    )
    for (i in seq_along(qm)) {
      lines(rep(qm[i], 2L), c(ql[i], qu[i]), col = cols[1L])
    }
    abline(0, 1, lty = 2L, col = cols[2L])
  }
}

rejectionSample <- function(
  target,
  dgenerator,
  rgenerator,
  constant,
  boundary,
  log = TRUE,
  maxIter = 100L
) {
  useLog <- log
  rm(log)
  numIters <- 0
  if (useLog) {
    while (TRUE) {
      u <- -rexp(1)
      x <- rgenerator()
      numIters <- numIters + 1
      if (x <= boundary[1] || x >= boundary[2]) {
        next
      }
      if (u < target(x) - constant - dgenerator(x)) {
        return(x)
      }
      if (numIters == maxIter) {
        stop("unable to obtain rejection sample after ", maxIter)
      }
    }
  } else {
    while (TRUE) {
      u <- runif(1)
      x <- rgenerator()
      numIters <- numIters + 1
      if (x <= boundary[1] || x >= boundary[2]) {
        next
      }
      if (u < (target(x) / (constant * dgenerator(x)))) {
        return(x)
      }
      if (numIters == maxIter) {
        stop("unable to obtain rejection sample after ", maxIter)
      }
    }
  }
}

sliceSample <- function(
  target,
  start,
  numSamples = 100L,
  width = NA,
  maxIter = 100L,
  boundary = c(-Inf, Inf),
  log = TRUE
) {
  useLog <- log
  rm(log)

  findMode <- function(target, start, boundary) {
    optimResult <- tryCatch(
      optim(
        start,
        target,
        method = "L-BFGS-B",
        lower = boundary[1L],
        upper = boundary[2L],
        hessian = TRUE,
        control = list(fnscale = -1)
      ),
      error = function(e) e
    )
    if (inherits(optimResult, "error")) {
      ## ||
      ##(is.finite(boundary[1]) && abs(optimResult$par - boundary[1]) < 1e-6) ||
      ##(is.finite(boundary[2]) && abs(optimResult$par - boundary[2]) < 1e-6))
      ## if optim fails, do own gradient ascent
      delta <- 1e-6
      while (
        start - 2 * delta <= boundary[1L] && start + 2 * delta >= boundary[2L]
      ) {
        delta <- delta / 2
      }
      if (start - delta <= boundary[1L]) {
        lh <- target(start)
        mh <- target(start + delta)
        rh <- target(start + 2 * delta)
        deriv <- -(3 * lh - 4 * mh + rh) / (2 * delta)
        hess <- (lh - 2 * mh + rh) / (delta^2)
      } else if (start + delta >= boundary[2L]) {
        lh <- target(start - 2 * delta)
        mh <- target(start - delta)
        rh <- target(start)
        deriv <- (3 * rh - 4 * mh + lh) / (2 * delta)
        hess <- (rh - 2 * mh + lh) / (delta^2)
      } else {
        rh <- target(start + delta)
        mh <- target(start)
        lh <- target(start - delta)
        deriv <- (rh - lh) / (2 * delta)
        hess <- (rh - 2 * mh + lh) / (delta^2)
      }
      step <- abs(deriv / hess) / 5

      lh <- start - step
      mh <- start
      rh <- start + step
      if (lh <= boundary[1L]) {
        lh <- if (is.finite(boundary[1L])) {
          boundary[1L] + delta
        } else {
          start - delta
        }
      }
      if (rh >= boundary[2L]) {
        rh <- if (is.finite(boundary[2L])) {
          boundary[2L] - delta
        } else {
          start + delta
        }
      }
      lf <- target(lh)
      mf <- target(mh)
      rf <- target(rh)

      ## keep going until we've got a middle that is higher than one side or the other
      while ((lf < mf && mf < rf) || (lf > mf && mf > rf)) {
        if (lf >= mf) {
          mh <- lh
          mf <- lf
          lh <- lh - step
          if (lh <= boundary[1L]) {
            lh <- if (is.finite(boundary[1L])) {
              boundary[1L] + delta
            } else {
              mh - delta
            }
          }
          lf <- target(lh)
        } else {
          mh <- rh
          mf <- rf
          rh <- rh + step
          if (rh >= boundary[2L]) {
            rh <- if (is.finite(boundary[2L])) {
              boundary[2L] - delta
            } else {
              mh + delta
            }
          }
          rf <- target(rh)
        }
      }

      optimResult <- tryCatch(
        optim(
          mh,
          target,
          method = "L-BFGS-B",
          lower = boundary[1L],
          upper = boundary[2L],
          hessian = TRUE,
          control = list(fnscale = -1)
        ),
        error = function(e) e
      )
    }
    optimResult
  }

  getInterval <- function(f, x, width, height, boundary) {
    r <- runif(1L)
    x.l <- x - r * width
    x.r <- x + (1 - r) * width

    if (is.finite(boundary[1L])) {
      while (x.l > boundary[1L] && f(x.l) > height) {
        x.l <- x.l - width
      }
      if (x.l < boundary[1L]) x.l <- boundary[1L]
    } else {
      while (f(x.l) > height) {
        x.l <- x.l - width
      }
    }
    if (is.finite(boundary[2L])) {
      while (x.r < boundary[2L] && f(x.r) > height) {
        x.r <- x.r + width
      }
      if (x.r > boundary[2L]) x.r <- boundary[2L]
    } else {
      while (f(x.r) > height) {
        x.r <- x.r + width
      }
    }

    c(x.l, x.r)
  }
  shrinkInterval <- function(x, x.p, int) {
    if (x.p > x) {
      int[2L] <- x.p
    } else {
      int[1L] <- x.p
    }
    int
  }

  f <- target
  if (is.na(width)) {
    optimResult <- findMode(target, start, boundary)

    if (!inherits(optimResult, "error")) {
      if (useLog == TRUE) {
        normalizingConstant <- NULL ## for R CMD check
        evalEnv <- list2env(list(
          target = target,
          normalizingConstant = optimResult$value
        ))
        f <- function(x) exp(target(x) - normalizingConstant)
        environment(f) <- evalEnv
        optimResult$value <- 1
        ## optimResult$hessian is this is theoretically unchanged by the transformation, since f'(x_0) = 0 && (h(x_0) - normConst) = 0, however it can be inaccurate numerically
        optimResult$hessian <- optimHess(
          optimResult$par,
          f,
          control = list(fnscale = -1)
        )
      }
      # width is derived from a normal approximation at the mode based on the equality of second derivatives
      # going out two standard deviations
      width <- 2 * abs(optimResult$hessian[1L] * sqrt(2 * pi))^(-1 / 3)
      if (is.nan(width) || is.infinite(width)) width <- 1000
    } else {
      width <- 1000
    }
  }

  result <- rep(NA_real_, numSamples)
  x <- start
  f.x <- f(x)
  if (f.x <= 1e-2) {
    ## if the starting point has a really low density, we grab a different one using
    ## rejection sampling
    mu <- NULL
    sigma <- NULL ## for R CMD check
    # us normal approximation with a slightly inflated standard deviation
    evalEnv <- list2env(list(mu = start, sigma = 1.15 * width / 2))
    r <- function() rnorm(1, mu, sigma)
    d <- function(x) dnorm(x, mu, sigma, log = TRUE)
    environment(r) <- evalEnv
    environment(d) <- evalEnv
    if (exists("optimResult") && useLog == TRUE) {
      if (inherits(optimResult, "error")) {
        stop("slice sampler failed: unable to determine initial curvature")
      }
      evalEnv$mu <- optimResult$par
      # special case for variance components
      if (
        is.finite(boundary[1L]) &&
          optimResult$par - boundary[1L] < evalEnv$sigma
      ) {
        evalEnv$sigma <- optimResult$par - boundary[1L]
      }
      if (
        is.finite(boundary[2L]) &&
          boundary[2L] - optimResult$par < evalEnv$sigma
      ) {
        evalEnv$sigma <- boundary[2L] - optimResult$par
      }
      c <- target(optimResult$par) - d(optimResult$par)
    } else {
      stop(
        "rejection start for case without optimization and/or not on log scale not yet implemented"
      )
    }

    tryResult <- tryCatch(
      x <- rejectionSample(target, d, r, c, boundary, maxIter = maxIter),
      error = function(e) e
    )
    if (inherits(tryResult, "error")) {
      warning(
        "rejection sample failed after ",
        maxIter,
        " iterations; dominating function may require hand-tuning"
      )
      x <- start
    } else {
      f.x <- f(x)
    }
  }
  for (i in seq_len(numSamples)) {
    u.p <- runif(1L, 0, f.x)
    int <- getInterval(f, x, width, u.p, boundary)
    for (j in seq_len(maxIter)) {
      x.p <- runif(1, int[1L], int[2L])
      f.x <- f(x.p)
      if (is.nan(f.x) || is.infinite(f.x)) {
        stop("slice sampler failed: likely due to underflow")
      }
      if (f.x > u.p) {
        break
      }
      int <- shrinkInterval(x, x.p, int)
    }
    if (j == maxIter) {
      stop("slice sampler failed: maxIter reached")
    }
    x <- x.p
    result[i] <- x
  }
  result
}
