rbart.priors <- list(
  cauchy = function(x, rel.scale) dcauchy(x, 0, rel.scale * 2.5, TRUE),
  gamma = function(x, rel.scale) {
    dgamma(x, shape = 2.5, scale = rel.scale * 2.5, log = TRUE)
  }
)
cauchy <- NULL ## for R CMD check

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
  k = NULL,
  prior.scale = NA_real_,
  power = 2.0,
  base = 0.95,
  split.probs = NULL,
  dart = FALSE,
  n.trees = 75L,
  n.samples = 1500L,
  n.burn = 1500L,
  n.chains = 4L,
  n.threads = min(dbarts::guessNumCores(), n.chains),
  combineChains = TRUE,
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
  factors = c("categorical", "indicators"),
  family = c("auto", "gaussian", "aft"),
  missing = c("incorporate", "error"),
  storage = c("double", "single"),
  updateState = TRUE,
  ...
) {
  matchedCall <- match.call()
  callingEnv <- parent.frame()
  family <- match.arg(family)

  # '...' is rejection-only - every dots name is diagnosed by name, never
  # forwarded.
  argNames <- names(matchedCall)[-1L]
  rejectUnknownDotsArgs(argNames, dbarts::rbart_vi)

  # factors/missing are forwarded formal defaults - redirectCall only
  # carries a name into the host dbartsData() call when the caller supplied
  # it, so an unsupplied one used to silently take dbartsData()'s own
  # default rather than the token this signature advertises. Resolve here,
  # in rbart_vi's own frame, and stamp the resolved value onto matchedCall
  # unconditionally so the call built from it below forwards it explicitly.
  # rbart_vi's formal defaults are kept textually identical to
  # dbartsData()'s, so this is draw-neutral.
  factors <- match.arg(factors)
  missing <- match.arg(missing)
  matchedCall$factors <- factors
  matchedCall$missing <- missing

  n.chains <- coerceOrError(n.chains, "integer")[1L]
  if (is.na(n.chains) || n.chains < 1L) {
    stop("'n.chains' must be a positive integer")
  }

  n.threads <- coerceOrError(n.threads, "integer")[1L]
  if (is.na(n.threads) || n.threads < 1L) {
    stop("'n.threads' must be a positive integer")
  }

  controlCall <- redirectCall(matchedCall, dbarts::dbartsControl)
  missingDefaults <- names(formals(rbart_vi))[
    names(formals(rbart_vi)) %in% names(formals(dbartsControl))
  ]
  missingDefaults <- missingDefaults[
    missingDefaults %not_in% names(controlCall)
  ]
  controlCall[missingDefaults] <- formals(rbart_vi)[missingDefaults]
  if ("n.threads" %in% missingDefaults) {
    controlCall[["n.threads"]] <- eval(controlCall[["n.threads"]])
  }
  control <- eval(controlCall, envir = callingEnv)

  control@call <- if (keepCall) matchedCall else call("NULL")
  control@n.burn <- control@n.burn %/% control@n.thin
  control@n.samples <- control@n.samples %/% control@n.thin
  control@printEvery <- control@printEvery %/% control@n.thin
  if (control@n.samples == 0L) {
    stop("no posterior draws will be taken after thinning")
  }

  keepSampler <- keepSampler || control@keepTrees

  # k enters EVALUATED, unlike bart2: rbart forwards its arguments to
  # dbarts through do.call from internal frames, where a stored symbol
  # cannot resolve. Without dbarts attached, pass hyperpriors as strings
  # (k = "chi(1.5)")
  priors <- buildSamplerPriors(
    matchedCall,
    power,
    base,
    sigdf,
    sigquant,
    nodeK = if (!is.null(matchedCall[["k"]])) k else NULL,
    priorScale = prior.scale,
    dart = dart,
    splitProbsDefault = formals(dbarts::rbart_vi)[["split.probs"]]
  )
  tree.prior <- priors$tree.prior
  node.prior <- priors$node.prior
  resid.prior <- priors$resid.prior

  if (is.null(matchedCall[["group.by"]])) {
    stop("'group.by' must be specified to use rbart_vi")
  }

  group.by.literal <- NULL
  # look for group.by in data, if supplied, first; by name, not position, so an
  # absent symbol falls through to the calling environment rather than silently
  # binding to the first column
  if (is.symbol(matchedCall[["group.by"]])) {
    try(
      group.by.literal <- data[[as.character(matchedCall[["group.by"]])]],
      silent = TRUE
    )
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
    if (is.symbol(matchedCall[["group.by.test"]])) {
      try(
        group.by.literal <- test[[as.character(matchedCall[[
          "group.by.test"
        ]])]],
        silent = TRUE
      )
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

    if (is.null(group.by.literal)) {
      try(
        group.by.literal <- data[[as.character(matchedCall[[
          "group.by.test"
        ]])]],
        silent = TRUE
      )
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
  }

  if (is.null(matchedCall[["prior"]])) {
    matchedCall[["prior"]] <- formals(rbart_vi)$prior
  }

  builtinTauPrior <- (is.symbol(matchedCall[["prior"]]) ||
    is.character(matchedCall[["prior"]])) &&
    any(names(rbart.priors) == as.character(matchedCall[["prior"]]))
  if (builtinTauPrior) {
    prior <- rbart.priors[[which(
      names(rbart.priors) == as.character(matchedCall[["prior"]])
    )]]
  }

  # survival (AFT) response ingestion (docs/design/survival.md, R surface):
  # a survival::Surv or two-column (time, status) response enters through the
  # formula's left-hand side, evaluated in the caller's data/environment as
  # group.by is. dbartsData's formula-path Surv refusal stays intact for every
  # other family; here the response is rewritten to the log event/censoring
  # time before dbartsData sees it, and the per-observation status rides the
  # control's bartcore.survival attribute, which the C bridge reads at ingestion
  # (the same attr-on-control convention applyGroupAttribute uses for
  # bartcore.groups). A Surv left-hand side auto-selects aft from "auto"; a
  # bare two-column response needs family = "aft" explicitly. rbart_vi's
  # matrix interface has no survival form; aft enters only through the formula.
  survivalStatus <- NULL
  survFormula <- NULL
  if (
    family %in% c("auto", "aft") && is.formula(formula) && length(formula) == 3L
  ) {
    responseExpr <- formula[[2L]]
    responseValue <- NULL
    if (!missing(data) && (is.list(data) || is.environment(data))) {
      try(
        responseValue <- eval(responseExpr, data, environment(formula)),
        silent = TRUE
      )
    }
    if (is.null(responseValue)) {
      try(
        responseValue <- eval(responseExpr, environment(formula)),
        silent = TRUE
      )
    }
    if (is.null(responseValue)) {
      try(responseValue <- eval(responseExpr, callingEnv), silent = TRUE)
    }

    responseIsSurv <- inherits(responseValue, "Surv")
    if (family == "aft" || (family == "auto" && responseIsSurv)) {
      survival <- extractSurvivalResponse(responseValue)
      if (is.null(survival)) {
        stop(
          "family \"aft\" needs a survival::Surv or two-column (time, ",
          "status) response through the formula's left-hand side"
        )
      }
      if (!is.null(matchedCall[["weights"]])) {
        stop("survival (aft) fits do not support 'weights' in this version")
      }
      if (!is.null(matchedCall[["subset"]])) {
        stop("survival (aft) fits do not support 'subset' in this version")
      }
      family <- "aft"
      survivalStatus <- survival$status
      # bind the log-time in a child of the formula's environment and point
      # the left-hand side at it, so model.frame builds the response from a
      # plain numeric while the right-hand side terms resolve as before
      survEnv <- new.env(parent = environment(formula))
      assign(".dbartsSurvivalResponse", survival$log.time, envir = survEnv)
      survFormula <- formula
      survFormula[[2L]] <- as.symbol(".dbartsSurvivalResponse")
      environment(survFormula) <- survEnv
    }
  }

  dataCall <- redirectCall(matchedCall, dbarts::dbartsData)
  if (!is.null(survFormula)) {
    # inject the rewritten formula via a symbol so eval does not rebind its
    # environment (a formula literal is a call to `~`, re-evaluated otherwise)
    survEvalEnv <- new.env(parent = callingEnv)
    survEvalEnv[[".dbartsSurvFormula"]] <- survFormula
    dataCall[["formula"]] <- quote(.dbartsSurvFormula)
    data <- eval(dataCall, envir = survEvalEnv)
  } else {
    data <- eval(dataCall, envir = callingEnv)
  }

  # both the R loop (predicts over data@x) and the in-core path would need
  # sparse-aware plumbing; reserved until a consumer appears
  if (predictorSourceIsSparse(data@x)) {
    stop(
      "rbart_vi does not support sparse predictor matrices; ",
      "dbarts() and bart2() do - use one of those, or pass a dense matrix"
    )
  }

  if (length(group.by) != length(data@y)) {
    stop(
      "'group.by' not of length equal to that of data; check for NAs in original data, and for name collisions with 'data' argument and calling environment"
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

  # a factor/logical/character response is a classification. rbart_vi's
  # random-effects model fits the 2-level (probit) case only; 3+ levels are
  # multinomial (bart2). Resolve here and pass an explicit family so the
  # per-chain dbarts() calls below do not each re-announce the verdict. A
  # numeric response is left alone - there is no message to deduplicate, so
  # the per-chain dbarts() calls resolve it themselves.
  family <- resolveClassificationFamily(
    data,
    family,
    "rbart_vi",
    c("gaussian", "aft")
  )
  # "auto" only survives resolveClassificationFamily on a numeric response
  # (aft's own is already forced explicit above, a Surv/two-column response);
  # a 0/1-coded numeric response is STILL a probit down in the per-chain
  # dbarts() call (resolveSamplerSpec's own binary test, R/spec.R), so redo
  # that test here rather than assume every unresolved "auto" is gaussian
  gatedFamily <- if (identical(family, "auto")) {
    uniqueY <- unique(data@y)
    if (length(uniqueY) == 2L && all(sort(uniqueY) == c(0, 1))) {
      "probit"
    } else {
      "gaussian"
    }
  } else {
    family
  }
  warnFamilyGatedArgs(argNames, gatedFamily)

  rbartArgs <- namedList(
    group.by,
    prior,
    keepTrainingFits,
    keepTestFits,
    callback
  )

  # the C bridge reads the per-observation status off the control's
  # bartcore.survival attribute at ingestion (the same attr-on-control
  # convention as bartcore.groups); both sampler paths carry the control, and
  # dbarts() permits the indirect aft path once the attribute is present
  if (!is.null(survivalStatus)) {
    attr(control, "bartcore.survival") <- survivalStatus
  }

  # in-core fast path: the built-in tau priors run rbart_vi's Gibbs blocks
  # inside the engine (one sampler, chains on worker threads, no
  # per-iteration R); custom prior functions and callbacks keep the R loop
  if (builtinTauPrior && is.null(callback)) {
    control@n.chains <- n.chains
    control@n.threads <- n.threads
    rel.scale <- if (!control@binary) sd(data@y) else 0.5
    attr(control, "bartcore.groups") <- list(
      indices = as.integer(group.by),
      n.groups = nlevels(group.by),
      prior = as.character(matchedCall[["prior"]]),
      rel.scale = rel.scale,
      n.steps = control@n.thin
    )

    samplerArgs <- namedList(
      formula = data,
      control,
      tree.prior,
      node.prior,
      resid.prior,
      sigma = as.numeric(sigest)
    )
    if (is.null(node.prior)) {
      samplerArgs[["node.prior"]] <- NULL
    }
    # the family stays off the call for numeric gaussian/binary (dbarts()
    # resolves it from the response); aft and a categorical-response probit
    # (resolved above) declare themselves so no per-chain re-announcement fires
    if (family %in% c("aft", "probit")) {
      samplerArgs[["family"]] <- family
    }

    fitResult <- rbart_vi_fit_bartcore(
      samplerArgs,
      rbartArgs,
      levels(group.by),
      seed
    )

    result <- packageRbartResults(
      control,
      data,
      group.by,
      group.by.test,
      fitResult$chainResults,
      combineChains,
      seed,
      keepSampler = FALSE
    )
    # one multi-chain sampler stands in for the per-chain fits; n.chains
    # stays on the object so the generics need not infer it from the list
    if (keepSampler) {
      result$fit <- list(fitResult$sampler)
    }
    return(result)
  }

  control@n.chains <- 1L
  control@n.threads <- max(control@n.threads %/% n.chains, 1L)
  if (n.chains > 1L && n.threads > 1L) {
    if (control@verbose) {
      warning("verbose output disabled for multiple threads")
    }
    control@verbose <- FALSE
  }

  samplerArgs <- namedList(
    formula = data,
    control,
    tree.prior,
    node.prior,
    resid.prior,
    sigma = as.numeric(sigest)
  )
  if (is.null(node.prior)) {
    samplerArgs[["node.prior"]] <- NULL
  }
  if (family %in% c("aft", "probit")) {
    samplerArgs[["family"]] <- family
  }
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
      # We draw sequentially from the given seed, one for each thread. To be polite
      # (more to match bart), we set the seed back when we're done.
      randomSeeds <- if (!is.na(seed)) {
        withFixedSeed(seed, sample.int(.Machine$integer.max, n.chains))
      } else {
        rep.int(NA_integer_, n.chains)
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
    # If the seed was passed in, a set.seed drives the chain seeds drawn
    # at sampler creation and any R-level draws; set the stream back when
    # done.
    fitChains <- function() {
      for (chainNum in seq_len(n.chains)) {
        chainResults[[chainNum]] <- rbart_vi_fit(
          1L,
          NA_integer_,
          samplerArgs,
          rbartArgs
        )
      }
      chainResults
    }
    chainResults <- if (!is.na(seed)) {
      withFixedSeed(seed, fitChains())
    } else {
      fitChains()
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

# One sampler run drives every kept sample through a per-sweep R closure
# instead of a run(0, 1) round trip apiece. The closure conditions each sweep
# on a fresh random-intercept draw (setOffset before the sweep) and, once the
# sweep lands, reads its fit back to update the residual and slice-sample tau -
# the identical operations in the identical order the round-trip loop ran, so
# the draws are unchanged. Under n.thin > 1 only the first sweep of each
# thinning block sets the offset and only the last is recorded, matching a
# run(0, 1) that internally thins.
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

  n.g <- data$n.g
  numRanef <- data$numRanef
  g.sel <- data$g.sel
  g <- data$g
  offset.orig <- data$offset.orig

  kIsModeled <- inherits(sampler$model@node.hyperprior, "dbartsChiHyperprior")
  # aft redraws its censored latents inside each sweep (against f + ranef the
  # offset carries), so the working response is pulled back like binary's
  familyIsAft <- sampler$model@family == "aft"
  posteriorClosure <- prior$posteriorClosure
  # used via evalEnv$b.sq in the nested postStep, which the linter misses
  evalEnv <- prior$evalEnv # nolint: object_usage_linter.

  numObservations <- length(sampler$data@y)
  numTestObservations <- NROW(sampler$data@x.test)
  numPredictors <- ncol(sampler$data@x)
  usesDart <- inherits(sampler$model@tree.prior, "dbartsDartPrior")

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
    samples$varcount <- matrix(NA_integer_, numPredictors, n.samples)
    if (usesDart) {
      samples$varprobs <- matrix(NA_real_, numPredictors, n.samples)
    }
  }

  # engine scratch the sweep callback reads: the same channels a run(0, 1)
  # fills, so the per-sample reads differ from the loop only by column index.
  # varcount aliases integer storage the engine writes as uint32 (see the
  # bridge entry point). Channels no read needs are left null so the engine
  # skips them.
  #
  # varcount is SINGLE-SLAB, numPredictors x n.samples, and stays so on a
  # multi-forest sampler: bartcore_runWithCallback pins the engine's declared
  # forest count to 1, so the layout here is a property of the call rather than
  # of who can reach it. Two guards keep a multi-forest sampler off rbart_vi
  # anyway, and neither is the same one - this R-loop path dies at the pre-run
  # rescale's setOffset(updateScale = TRUE), which refuseBCFMutation refuses,
  # while the in-core path (which owns no buffer) dies at the grouped x
  # multi-forest refusal in R/spec.R.
  raw <- list(
    sigma = numeric(n.samples),
    train = matrix(0.0, numObservations, n.samples)
  )
  if (kIsModeled) {
    raw$k <- numeric(n.samples)
  }
  if (!isWarmup) {
    raw$varcount <- matrix(0L, numPredictors, n.samples)
    if (numTestObservations > 0L) {
      raw$test <- matrix(0.0, numTestObservations, n.samples)
    }
    if (usesDart) {
      raw$varprobs <- matrix(0.0, numPredictors, n.samples)
    }
  }

  n.thin <- control@n.thin
  ranef <- NULL
  ranef.vec <- NULL
  callbackError <- NULL

  # drawn before the sweep it conditions
  preStep <- function() {
    resid <- state$y.st - state$treeFit.train
    post.var <- 1.0 / (n.g / state$sigma^2.0 + 1.0 / state$tau^2.0)
    post.mean <- (n.g / state$sigma^2.0) *
      sapply(seq_len(numRanef), function(j) mean(resid[g.sel[[j]]])) *
      post.var
    ranef <<- rnorm(numRanef, post.mean, sqrt(post.var))
    ranef.vec <<- ranef[g]
    sampler$setOffset(
      ranef.vec + if (!is.null(offset.orig)) offset.orig else 0,
      isWarmup
    )
  }

  # reads back the sweep recorded in column i, then updates tau and books the
  # sample
  postStep <- function(i) {
    state$treeFit.train <<- raw$train[, i] - ranef.vec
    if (control@binary || familyIsAft) {
      state$y.st <<- sampler$getLatents(state$y.st)
    }
    state$sigma <<- raw$sigma[i]

    evalEnv$b.sq <- sum(ranef^2.0)
    state$tau <<- sliceSample(
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
      .Call(C_dbarts_assignInPlace, samples$varcount, i, raw$varcount[, i])
    }
    if (!is.null(samples$varprobs)) {
      .Call(C_dbarts_assignInPlace, samples$varprobs, i, raw$varprobs[, i])
    }
    if (
      !is.null(samples$yhat.test) &&
        numTestObservations > 0L &&
        rbartArgs$keepTestFits
    ) {
      .Call(C_dbarts_assignInPlace, samples$yhat.test, i, raw$test[, i])
    }
    if (!is.null(samples$k)) {
      .Call(C_dbarts_assignInPlace, samples$k, i, raw$k[i])
    }
    if (!isWarmup && !is.null(rbartArgs$callback)) {
      names(ranef) <- data$g.levels
      callbackValue <- rbartArgs$callback(
        state$treeFit.train,
        raw$test[, i],
        ranef,
        state$sigma,
        state$tau
      )
      if (is.null(samples$callback)) {
        samples$callback <<- matrix(
          NA_real_,
          length(callbackValue),
          control@n.samples,
          dimnames = list(names(callbackValue), NULL)
        )
      }
      .Call(C_dbarts_assignInPlace, samples$callback, i, callbackValue)
    }

    if (verbose && i %% control@printEvery == 0L) {
      cat("iter: ", i, "\n", sep = "")
    }
  }

  # fires before every sweep; the offset is held across a thinning block, so
  # only its first sweep books the previous sample and draws the next intercept
  runCallback <- function(sweepIndex) {
    tryCatch(
      {
        if (sweepIndex %% n.thin == 0L) {
          group <- sweepIndex %/% n.thin
          if (group > 0L) {
            postStep(group)
          }
          preStep()
        }
        FALSE
      },
      error = function(e) {
        callbackError <<- e
        TRUE
      }
    )
  }

  .Call(
    C_dbarts_bartcore_runWithCallback,
    sampler$getPointer(),
    0L,
    as.integer(n.samples),
    raw,
    runCallback,
    environment()
  )
  if (!is.null(callbackError)) {
    stop(callbackError)
  }

  # the last block's sweep landed after the final callback; book it now
  postStep(n.samples)

  list(state = state, samples = samples)
}

# The in-core counterpart of rbart_vi_fit: the grouping rides an internal
# attribute on the control object, so the sampler runs every Gibbs block -
# including the random intercepts and tau - inside the engine. One
# multi-chain sampler runs warmup and sampling in two calls (warmup
# supplies first.tau/first.sigma with trees and training fits not kept),
# and the results split into the per-chain shapes rbart_vi_fit produces so
# the packaging is shared.
rbart_vi_fit_bartcore <- function(samplerArgs, rbartArgs, groupLevels, seed) {
  # a single chain draws through R's generator; several chains seed their
  # own generators from R's stream at creation - reproducible either way
  withSeedIfGiven <- function(expr) {
    if (is.na(seed)) expr else withFixedSeed(seed, expr)
  }

  withSeedIfGiven({
    sampler <- do.call(dbarts::dbarts, samplerArgs)
    sampler$control@call <- samplerArgs$control@call

    control <- sampler$control
    kIsModeled <- inherits(sampler$model@node.hyperprior, "dbartsChiHyperprior")

    sampler$sampleTreesFromPrior(updateState = FALSE)

    firstTau <- firstSigma <- firstK <- NULL
    if (control@n.burn > 0L) {
      warmupControl <- control
      warmupControl@keepTrees <- FALSE
      warmupControl@keepTrainingFits <- FALSE
      sampler$setControl(warmupControl)
      warmupSamples <- sampler$run(0L, control@n.burn, updateState = FALSE)
      firstTau <- warmupSamples$tau
      if (!control@binary) {
        firstSigma <- warmupSamples$sigma
      }
      if (kIsModeled) {
        firstK <- warmupSamples$k
      }
      sampler$setControl(control)
    }

    samples <- sampler$run(0L, control@n.samples)
  })

  # per-chain slices in the shapes rbart_vi_fit produces
  chainMatrix <- function(x, i) {
    if (is.null(x)) {
      NULL
    } else if (length(dim(x)) > 2L) {
      array(x[,, i], dim(x)[1L:2L])
    } else {
      x
    }
  }
  chainVector <- function(x, i) {
    if (is.null(x)) {
      NULL
    } else if (is.matrix(x)) {
      x[, i]
    } else {
      x
    }
  }

  n.chains <- control@n.chains
  chainResults <- vector("list", n.chains)
  for (i in seq_len(n.chains)) {
    ranef <- chainMatrix(samples$ranef, i)
    rownames(ranef) <- groupLevels
    chainResult <- list(
      sampler = sampler,
      ranef = ranef,
      firstTau = chainVector(firstTau, i),
      firstSigma = if (!control@binary) chainVector(firstSigma, i),
      tau = chainVector(samples$tau, i),
      sigma = if (!control@binary) chainVector(samples$sigma, i),
      yhat.train = if (rbartArgs$keepTrainingFits) {
        chainMatrix(samples$train, i)
      },
      yhat.test = if (rbartArgs$keepTestFits) chainMatrix(samples$test, i),
      callback = NULL,
      varcount = chainMatrix(samples$varcount, i)
    )
    if (!is.null(samples$varprobs)) {
      chainResult$varprobs <- chainMatrix(samples$varprobs, i)
    }
    if (kIsModeled) {
      chainResult$firstK <- chainVector(firstK, i)
      chainResult$k <- chainVector(samples$k, i)
    }
    chainResults[[i]] <- chainResult
  }

  namedList(sampler, chainResults)
}

rbart_vi_fit <- function(.chain.num.ignored, seed, samplerArgs, rbartArgs) {
  # clusterMap passes a chain index positionally as the first argument, but this
  # fit uses it only to occupy the slot; the value is never read.

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
  sampler$setControl(control)

  y <- sampler$data@y
  rel.scale <- if (!control@binary) sd(y) else 0.5

  g <- as.integer(rbartArgs$group.by)
  g.levels <- levels(rbartArgs$group.by)
  numRanef <- nlevels(rbartArgs$group.by)
  g.sel <- lapply(seq_len(numRanef), function(j) g == j)
  n.g <- sapply(g.sel, sum)
  offset.orig <- sampler$data@offset
  data <- namedList(n.g, numRanef, g.sel, g, g.levels, offset.orig)

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

  numObservations <- length(sampler$data@y) # nolint: object_usage_linter.
  numTestObservations <- NROW(sampler$data@x.test) # nolint: object_usage_linter.

  kIsModeled <- inherits(sampler$model@node.hyperprior, "dbartsChiHyperprior")

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

    # sweep the forest once before the first ranef draw so it absorbs the
    # response mean: a prior fit predicts the response midpoint, and for a
    # skewed response that offset would otherwise leak into the intercepts
    state$treeFit.train <- as.vector(sampler$run(0L, 1L)$train) - ranef.vec

    runResult <- rbart_vi_run(
      sampler,
      data,
      state,
      prior,
      FALSE,
      control@n.burn,
      TRUE,
      rbartArgs
    )
    state <- runResult$state

    firstTau <- runResult$samples$tau
    firstSigma <- runResult$samples$sigma
    firstK <- runResult$samples$k

    if (control@keepTrees != oldKeepTrees) {
      control@keepTrees <- TRUE
      sampler$setControl(control)
    }
  } else {
    oldKeepTrees <- control@keepTrees
    if (oldKeepTrees) {
      control@keepTrees <- FALSE
      sampler$setControl(control)
    }

    sampler$setOffset(
      ranef.vec + if (!is.null(offset.orig)) offset.orig else 0,
      TRUE
    )

    state$treeFit.train <- as.vector(sampler$run(0L, 1L)$train) - ranef.vec

    if (oldKeepTrees) {
      control@keepTrees <- TRUE
      sampler$setControl(control)
    }

    firstTau <- NULL
    firstSigma <- NULL
    firstK <- NULL
  }

  runResult <- rbart_vi_run(
    sampler,
    data,
    state,
    prior,
    verbose,
    control@n.samples,
    FALSE,
    rbartArgs
  )
  # the returned state is intentionally discarded; only the samples are kept

  tau <- runResult$samples$tau
  sigma <- runResult$samples$sigma
  ranef <- runResult$samples$ranef
  yhat.train <- runResult$samples$yhat.train
  yhat.test <- runResult$samples$yhat.test
  k <- runResult$samples$k
  callback <- runResult$samples$callback
  varcount <- runResult$samples$varcount
  varprobs <- runResult$samples$varprobs

  sampler$setOffset(if (!is.null(offset.orig)) offset.orig else NULL, FALSE)

  control@updateState <- oldUpdateState
  sampler$setControl(control)

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
  if (!is.null(varprobs)) {
    result$varprobs <- varprobs
  }
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

  result <- list(
    call = control@call,
    y = data@y,
    group.by = group.by,
    family = chainResults[[1L]]$sampler$model@family
  )
  if (!responseIsBinary) {
    result$sigest <- chainResults[[1L]]$sampler$data@sigma
  }
  if (!is.null(group.by.test)) {
    result[["group.by.test"]] <- group.by.test
  }
  # needed to extract ppd; mirrors bart2's packageBartResults
  if (!is.null(data@weights) && length(data@weights) > 0L) {
    result$weights <- data@weights
    if (!is.null(data@weights.test) && length(data@weights.test) > 0L) {
      result$weights.test <- data@weights.test
    }
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
    if (!is.null(chainResults[[1L]]$varprobs)) {
      result$varprobs <- convertSamplesFromDbartsToBart(
        array(
          sapply(chainResults, function(x) x$varprobs),
          c(dim(chainResults[[1L]]$varprobs), n.chains)
        ),
        n.chains,
        combineChains
      )
    }
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
    if (!is.null(chainResults[[1L]]$varprobs)) {
      result$varprobs <- chainResults[[1L]]$varprobs
    }
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

  result$seed <- if (!is.na(seed)) {
    withFixedSeed(seed, .GlobalEnv$.Random.seed)
  } else {
    if (!exists(".Random.seed", .GlobalEnv)) {
      runif(1L)
    }
    .GlobalEnv$.Random.seed
  }

  class(result) <- "rbart"
  result
}
