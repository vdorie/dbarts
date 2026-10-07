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
  sd = NULL,
  power = 2,
  base = 0.95,
  split.probs = NULL,
  drop = TRUE,
  sigest = NULL,
  seed = NULL,
  factors = c("categorical", "indicators"),
  family = c("auto", "gaussian", "probit", "logistic"),
  leaf.prior = NULL,
  n.cuts = 100L,
  useQuantiles = FALSE,
  n.thin = 1L,
  storage = c("double", "single"),
  tree.prior = NULL,
  parallel = getOption("dbarts.parallel", "auto"),
  cl = NULL,
  control = dbarts::dbartsControl(),
  sigma = NULL,
  ...
) {
  matchedCall <- match.call()
  # the creation-time estimate is 'sigest'; 'sigma' is 0.9-x's spelling of it,
  # accepted for one release
  sigmaSupplied <- !missing(sigma) && !is.null(sigma)
  sigest <- resolveRenamedSigma(
    !sigmaSupplied,
    missing(sigest),
    sigma,
    sigest,
    "xbart"
  )
  sigest <- if (sigmaSupplied) {
    resolveSigestArg(sigest, "xbart", "silent", "sigma")
  } else {
    resolveSigestArg(sigest, "xbart", "refuse")
  }
  if (sigmaSupplied) {
    matchedCall["sigest"] <- list(if (is.na(sigest)) NULL else sigest)
    matchedCall$sigma <- NULL
  }
  # '...' exists only so a retired argument name reaches a message naming
  # its successor; R refuses an unknown name before any body runs
  supplied <- dotNames(...)
  refuseForeignFrontDoorArgs(supplied, "xbart", names(formals(dbarts::xbart)))
  consolidated <- resolveConsolidatedArgs(
    matchedCall,
    supplied,
    "xbart",
    parent.frame(1L)
  )
  if (length(consolidated) > 0L) {
    matchedCall[names(consolidated)] <- NULL
  }

  currEnv <- sys.frame(sys.nframe())
  evalEnv <- parent.frame(1L)

  # the four flat knobs: shape-checked at the surface, mirroring
  # dbartsControl's own validity messages, then folded into the control
  # xbart builds for itself below
  n.cuts <- coerceOrError(n.cuts, "integer")
  if (length(n.cuts) == 0L || anyNA(n.cuts) || any(n.cuts <= 0L)) {
    stop("'n.cuts' must contain positive integers")
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

  # One precedence rule: a flat formal the caller named wins over the control's
  # slot, and a slot they did not name flat stands. n.burn is xbart's own grid
  # axis and n.threads its sweep width - neither is the control field of the
  # same name - so they are excluded. The fields below are forced after both:
  # a sweep runs one chain per cell, keeps no trees or training fits, stores
  # no state, prints nothing, and seeds each cell itself.
  control <- refuseFitStateControl(control, "xbart")
  resolved <- mergeFrontDoorControl(
    control,
    matchedCall,
    list(
      n.cuts = n.cuts,
      useQuantiles = useQuantiles,
      n.thin = n.thin,
      storage = storage,
      n.trees = n.trees
    )
  )
  n.cuts <- resolved$n.cuts
  useQuantiles <- resolved$useQuantiles
  n.thin <- resolved$n.thin
  storage <- resolved$storage
  # the grid and the RNG block read these; only the four knobs above are the
  # control's own copies
  n.trees <- resolved$n.trees
  # 'seed' is resolved from its value, not from whether the call named it: a
  # wrapper forwarding its own seed = NULL must still defer to the control.
  seed <- resolveSeedArg(seed, "xbart")
  if (is.na(seed)) {
    seed <- control@seed
  }
  control@n.cuts <- n.cuts
  control@useQuantiles <- useQuantiles
  control@n.thin <- n.thin
  control@storage <- storage
  control@n.chains <- 1L
  control@n.threads <- 1L
  control@keepTrees <- FALSE
  control@keepTrainingFits <- FALSE
  control@updateState <- FALSE
  control@verbose <- FALSE
  # the seed is the sweep's own: it is read above as 'seed' and drives the
  # per-replication splits and the per-sampler seeds every cell sampler is
  # created under, each handed to a fresh copy of this control's @seed below.
  # Left set here it would additionally seed every unit's first sampler
  # identically, off this one value instead of its own draw.
  control@seed <- NA_integer_

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
    control@call <- expandForwardedCall(matchedCall, evalEnv)
  }

  # named ahead of the data build, matching bart()/dbarts(), so
  # a bad family is refused before the response is ingested rather than after
  familySpec <- resolveFamily(
    matchedCall$family,
    eval(formals(dbarts::xbart)$family),
    "xbart",
    evalEnv
  )
  family <- familySpec@token

  dataCall <- redirectCall(
    matchedCall,
    quoteInNamespace(dbartsData),
    formula,
    data,
    subset,
    weights,
    offset,
    factors
  )
  # a count-matrix data object declares the multinomial model, whose fitted
  # quantity is K probabilities per observation; every loss this function
  # evaluates is written against one location, so the fit would be scored as
  # though the slab were n rows. Refused ahead of the family resolution below,
  # which resolves counts to multinomial from "auto"
  refuseCountsCarryingData(formula, "xbart()")
  refuseBasesCarryingData(
    formula,
    "xbart()",
    "xbart cross-validates single-forest models"
  )
  refuseResponseFreeFormula(formula, "xbart()")
  refuseForestTerm(formula, "xbart")
  data <- withMatrixResponseRestated(
    "xbart()",
    family,
    withBinaryResponsePrecision(family, eval(dataCall, evalEnv))
  )
  # a Surv formula response silently becomes log(time) with the censoring
  # status parked as an attribute (dbartsData()'s own short-circuit, which
  # has no family vocabulary to refuse it by) - xbart() reads neither the
  # attribute nor 'family' the way dbarts()'s aft/hazard blocks do, so left
  # unrefused this would cross-validate log(time) as an ordinary gaussian
  # response and silently discard every censored observation's status
  if (!is.null(attr(data, "survivalStatus"))) {
    stop(
      "xbart() does not cross-validate a survival (Surv) response; fit ",
      "with bart()/dbarts() using family = \"aft\" or \"hazard\" instead"
    )
  }
  data@n.cuts <- recycleNumCuts(control@n.cuts, ncol(data@x))
  data@sigma <- sigest

  # a factor/logical/character response is a classification; xbart cross-
  # validates the 2-level (probit) case only, never multinomial. A numeric
  # response takes the 0/1-vs-continuous path.
  autoDescription <- if (family == "auto") describeAutoResponse(data)
  family <- resolveClassificationFamily(
    data,
    family,
    "xbart",
    "gaussian"
  )
  if (data@response.type == "numeric") {
    uniqueResponses <- unique(data@y)
    responseIsBinary <- length(uniqueResponses) == 2L &&
      all(sort(uniqueResponses) == c(0, 1))
    if (family == "auto") {
      family <- if (responseIsBinary) "probit" else "gaussian"
    } else if (family != "gaussian" && !responseIsBinary) {
      # gaussian on a 0/1 response is a legitimate request; the binary
      # families need latent-variable coding
      refuseNonBinaryResponse(uniqueResponses, family)
    }
  }
  if (isBinaryFamily(family) && responseHasSingleClass(data@y)) {
    refuseNonBinaryResponse(unique(data@y[!is.na(data@y)]), family)
  }
  control@binary <- isBinaryFamily(family)

  # the shared weight policy (R/spec.R's enforceWeightPolicy): a probit has
  # no tractable weighted latent-variable form and is refused (an all-ones
  # courtesy excepted), a logistic model requires positive integer count
  # weights, and a gaussian fit is unrestricted. xbart's own family is always
  # gaussian/probit/logistic, so the function's ordinal/nbinom branches never
  # fire here - the same function every other entry point reaches them with.
  weightPolicy <- enforceWeightPolicy(data, family)
  # a probit 0/1 weight vector resolves to a row mask, which a sampler takes
  # and this does not: the folds partition the rows themselves, and each fit
  # is built and discarded inside the C loop with no channel to install one
  if (!is.null(weightPolicy$active)) {
    stop(
      "xbart does not accept weights of 0 and 1 under family \"",
      family,
      "\": they mark rows out of the likelihood, which cross-validation ",
      "already partitions; drop those rows before calling xbart"
    )
  }
  data <- weightPolicy$data

  # An unsupplied sigest is estimated per fold from the fold's training rows,
  # so no fold's residual prior reads its held-out responses. The estimate on
  # all rows still runs once, here, to raise any refusal or fallback once and
  # to choose the per-fold route: where it fell back to the marginal sd (a
  # sparse design, or no residual degrees of freedom), each fold takes its
  # own marginal sd and no fold attempts the linear fit again.
  sigmaPerFold <- NULL
  if (is.na(data@sigma) && !control@binary) {
    fellBack <- FALSE
    data@sigma <- withCallingHandlers(
      estimateStartingSigma(data),
      dbartsSigmaFallbackWarning = function(w) fellBack <<- TRUE
    )
    sigmaPerFold <- if (fellBack) "marginal" else "linear"
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

  # the grid comes from the formal where the caller named it and from the
  # control's own slot where they did not, one value either way
  n.trees <- coerceOrError(n.trees, "integer")
  if (anyNA(n.trees) || any(n.trees <= 0L)) {
    stop("'n.trees' must contain only positive integers")
  }

  # a supplied leaf.prior contributes the leaf model shape - normal(k), the
  # default, linear(columns, k), or gp(columns, k, ...), whose designated
  # covariate columns resolve against the model matrix; the k argument
  # drives the k grid as always, with a k inside the supplied prior standing
  # in for a missing k argument, and the sd argument drives the same axis in
  # absolute spreads, with a named sd inside the prior standing in for it
  leafSpec <- NULL
  if (!is.null(matchedCall[["leaf.prior"]])) {
    leafSpec <- evalInVocabulary(
      matchedCall[["leaf.prior"]],
      dbartsPriors[c("normal", "linear", "gp", "chi", "invchi")],
      evalEnv,
      resolvedAs(
        "leaf.prior",
        c("NULL", "dbartsLeafPrior"),
        "leaf prior specification"
      )
    )
  }

  # the k axis is 0.9-x's numeric vector or a list whose entries are numbers
  # and hyperprior objects, so one sweep scores fixed k against modelled k.
  # Read off the unevaluated argument in the prior vocabulary, as every other
  # prior argument here is, so k = list(1, chi()) resolves without the
  # constructor standing on the caller's search path. An absent k is ONE cell
  # at the front door's own default for the response type - fixed 2
  # continuous, chi(1.5, 2) binary - so a default xbart call scores the model
  # a default bart call fits; a k carried by a supplied leaf.prior stands in
  # for a missing argument.
  kGiven <- !is.null(matchedCall[["k"]])
  sdGiven <- !is.null(matchedCall[["sd"]])
  if (kGiven && sdGiven) {
    stop(
      "give either 'k' (relative to each fold's scale) or 'sd' (absolute ",
      "spreads on the family's scale) as the grid, not both"
    )
  }
  leafSd <- if (is.null(leafSpec)) NULL else leafSpec@prior.sd
  if (!is.null(leafSd) && (kGiven || sdGiven)) {
    grid <- if (kGiven) "k" else "sd"
    stop(
      "the leaf prior's 'sd' and the '",
      grid,
      "' grid both state the spread; drop the leaf prior's 'sd' and name ",
      "the spreads in the 'sd' grid"
    )
  }
  if (sdGiven && !is.null(leafSpec) && !is.null(leafSpec@k)) {
    stop(
      "the leaf prior's 'k' and the 'sd' grid both state the spread; drop ",
      "the leaf prior's 'k'"
    )
  }
  # the axis is one of k or sd; a named sd in the leaf prior is a one-cell sd
  # axis. Each cell is a leaf hyperprior plus the anchor it is relative to,
  # NA on the k axis, where the fold's own data fixes the anchor.
  sdAxis <- sdGiven || !is.null(leafSd)
  if (sdAxis) {
    sdSpec <- if (sdGiven) {
      evalInVocabulary(matchedCall[["sd"]], dbartsPriors, evalEnv)
    } else {
      leafSd
    }
    sdGrid <- resolveSdGrid(sdSpec)
    kGrid <- lapply(sdGrid, function(cell) cell$leaf.hyperprior)
    kAnchors <- vapply(sdGrid, function(cell) cell$prior.scale, 0.0)
    kLabels <- vapply(sdGrid, function(cell) cell$label, "")
    sortKeys <- vapply(sdGrid, function(cell) cell$sort.key, 0.0)
  } else {
    kSpec <- if (kGiven) {
      evalInVocabulary(matchedCall[["k"]], dbartsPriors, evalEnv)
    } else if (!is.null(leafSpec)) {
      leafSpec@k
    }
    kGrid <- resolveKGrid(kSpec, control@binary)
    kAnchors <- rep(NA_real_, length(kGrid))
    kLabels <- vapply(kGrid, kGridLabel, "")
    sortKeys <- vapply(kGrid, kGridSortKey, 0.0)
  }
  # swept most shrunk first - largest k, smallest sd - so every warm start
  # comes from a simpler forest than the cell before it; a modelled cell has
  # no fixed value to order by and sweeps last, in the order it was written.
  # kOrder un-permutes the reported axis back to the caller's order once the
  # result array is final
  kOrder <- order(sortKeys, decreasing = TRUE)
  kGrid <- kGrid[kOrder]
  kAnchors <- kAnchors[kOrder]

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
  # rule as leaf.prior/k above - power and base are xbart's grid axes, so
  # cellModel overwrites them on the object every cell regardless of what is
  # supplied here, while the object's non-grid content (a cgm's split.probs,
  # a dart's Dirichlet hyperparameters) rides every cell unchanged.
  # split.probs would only duplicate what a supplied tree.prior already
  # specifies, so it collides with it; power/base/k do not, since they are
  # grid axes rather than duplicates - this deliberately differs from
  # bart2's tree.prior, which does collide with power/base (R/bart.R's
  # buildSamplerPriors), because there they are ordinary scalars, not a grid.
  if (!is.null(matchedCall[["tree.prior"]])) {
    refuseColliding(matchedCall, "tree.prior", "split.probs")
    tree.prior <- evalInVocabulary(
      matchedCall[["tree.prior"]],
      dbartsPriors[c("cgm", "dart")],
      evalEnv,
      resolvedAs("tree.prior", "dbartsTreePrior", "tree prior specification")
    )
  } else {
    tree.prior <- cgm(power[1L], base[1L], split.probs)
  }
  tree.prior <- resolveSplitProbabilities(tree.prior, data)

  # the leaf model is built at the first cell's k; cellModel swaps the
  # hyperprior itself as the sweep moves along the axis
  kValue <- kGridValue(kGrid[[1L]])
  if (is.null(leafSpec)) {
    leafPrior <- quote(normal(k))
    leafPrior[[1L]] <- quoteInNamespace(normal)
    leafPrior[[2L]] <- kValue
    leafPrior <- eval(leafPrior)
  } else {
    # the grid replaces the supplied prior's own k or sd; the leaf model's
    # shape is all it keeps
    leafPrior <- if (is(leafSpec, "dbartsLinearPrior")) {
      resolveLeafCovariates(linear(leafSpec@columns, kValue), data)
    } else if (is(leafSpec, "dbartsGPPrior")) {
      resolveLeafCovariates(
        gp(
          leafSpec@columns,
          kValue,
          leafSpec@lengthscale,
          leafSpec@max.leaf.size
        ),
        data
      )
    } else {
      normal(kValue)
    }
  }
  leaf.hyperprior <- kGrid[[1L]]

  # a binary family runs on a fixed unit latent scale (R/spec.R's
  # fixedUnitScale rule): the residual prior is overridden where one is
  # given, not just where one is missing, matching the shared resolver, so a
  # caller cannot silently cross-validate an unfixed residual scale under a
  # family that has none. Otherwise the prior comes off the family object it
  # rides, or off the retired flat spelling this door still reads, which is
  # refused beside a family that named 'sigma' too unless the two agree.
  residPrior <- reconcileResidPrior(
    consolidatedResidPrior(consolidated),
    "resid.prior",
    familySpec
  )
  refuseSigestUnderFixedPrior(
    residPrior,
    sigest,
    if (sigmaSupplied) "sigma" else "sigest"
  )
  resid.prior <- if (control@binary) {
    fixed(1)
  } else if (!is.null(residPrior)) {
    residPrior
  } else {
    chisq()
  }
  # a fixed residual scale reads no estimate, so no fold fits one
  if (is(resid.prior, "dbartsFixedPrior")) {
    sigmaPerFold <- NULL
  }
  model <- newValidated(
    "dbartsModel",
    tree.prior,
    leafPrior,
    leaf.hyperprior,
    resid.prior,
    family = family,
    # an sd cell's anchor is held across folds, created or re-modelled:
    # cellModel swaps it per cell, and the setModel branch re-derives it
    prior.scale = kAnchors[[1L]],
    leaf.scale = defaultLeafScale(family)
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
  n.burn <- coerceOrError(n.burn, "integer")
  refuseThreeElementBurn(n.burn)
  n.burn <- rep_len(n.burn, 2L)
  if (anyNA(n.burn) || any(n.burn < 0L)) {
    stop("'n.burn' must contain non-negative integers")
  }
  n.threads <- coerceOrError(n.threads, "integer")
  if (length(n.threads) != 1L) {
    stop("'n.threads' must be of length 1")
  }
  if (is.na(n.threads)) {
    stop(naThreadsMessage)
  }
  if (n.threads <= 0L) {
    stop("'n.threads' must be a positive integer")
  }
  if (!is.null(cl) && (!inherits(cl, "cluster") || length(cl) == 0L)) {
    stop("'cl' must be NULL or a non-empty cluster from the parallel package")
  }
  parallel <- match.arg(parallel, c("auto", "fork", "socket"))
  if (parallel == "fork" && .Platform$OS.type == "windows") {
    stop("parallel = \"fork\" is not available on Windows; use \"socket\"")
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

  # announced once here, after the refusals and in the calling process,
  # however many replications or workers the cross-validation fans out to
  if (!is.null(autoDescription)) {
    announceAutoFamily(verbose, family, autoDescription)
  }

  # a replication draws a data split and sweeps every parameter cell over
  # it. Chains warm-start only across cells - the training data is
  # identical there, so carrying trees is sound - and never across folds
  # or splits, whose held-out rows the previous training set contained.
  # Restarting each split from a fresh forest keeps a slow-mixing cell from
  # remembering the previous fold and scoring optimistically on its own
  # held-out rows. Tree counts are fixed at a sampler's creation, so they
  # vary slowest and each count gets a fresh fit per split.
  kLength <- length(kGrid)
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
    n.samples = control@n.samples,
    n.burn,
    n.trees,
    kHyperpriors = kGrid,
    kAnchors,
    power,
    base,
    cells,
    lossFunction,
    sigmaPerFold,
    # a worker starts a fresh session; it is handed this one's warned-once
    # keys so a key already warned here stays silent there. A new key fires
    # once per worker, and the caller's deduplication reports it once
    onceKeys = warnedOnceKeys(),
    warn = getOption("warn")
  )

  # work is distributed over (replication, fold) UNITS rather than over
  # replication ranges: each fold of a replication is an independent fit, so
  # a k-fold call with a single replication fills every worker instead of
  # running the whole sweep on one. Cells still run in a fixed order within a
  # unit, which is what the k warm start inside a fold rides on.
  numFolds <- if (method == "k-fold") length(foldSizes) else 1L
  numUnits <- n.reps * numFolds
  numChunks <- max(
    1L,
    min(if (is.null(cl)) n.threads else min(n.threads, length(cl)), numUnits)
  )
  workerKind <- if (numChunks == 1L) {
    "session"
  } else if (!is.null(cl)) {
    "cluster"
  } else if (
    parallel == "fork" ||
      (parallel == "auto" &&
        .Platform$OS.type != "windows" &&
        !identical(Sys.getenv("RSTUDIO"), "1") &&
        !identical(Sys.getenv("POSITRON"), "1") &&
        !identical(.Platform$GUI, "AQUA"))
  ) {
    "fork"
  } else {
    "socket"
  }
  chunkIndices <- parallel::splitIndices(numUnits, numChunks)

  # each replication draws its data split from its own seed, and each unit's
  # sampler(s) their own: one seed per sampler a unit creates - one per
  # distinct tree count, since a fresh sampler is minted only when the tree
  # count changes (see xbartRunUnits) - all drawn here from the call's seed
  # alone, so a seed reproduces at any 'n.threads': no draw depends on which
  # worker ran a unit, on how many there were, or on RNGkind(), since a
  # unit's seed reaches its sampler through the control's seed slot rather
  # than a worker-side set.seed(). A supplied seed leaves the caller's stream
  # untouched; without one the seeds come off that stream, advancing it
  # exactly that far and no further at any thread count.
  numTreeCounts <- length(n.trees)
  seeds <- if (!is.na(seed)) {
    withFixedSeed(
      seed,
      sample.int(.Machine$integer.max, n.reps + numUnits * numTreeCounts)
    )
  } else {
    sample.int(.Machine$integer.max, n.reps + numUnits * numTreeCounts)
  }
  splitSeeds <- seeds[seq_len(n.reps)]
  # row i holds unit i's per-sampler seeds, in n.trees order
  unitSeeds <- matrix(
    seeds[n.reps + seq_len(numUnits * numTreeCounts)],
    numUnits,
    numTreeCounts,
    byrow = TRUE
  )

  # every stream below is one of those seeds, and at a single worker the
  # units run in THIS process, so the caller's own stream is saved across the
  # whole dispatch rather than left wherever the last fold stopped
  runUnits <- function() {
    # the splits are drawn here rather than on the worker that runs them, so
    # this process's RNGkind() governs them at every thread count regardless
    # of what a worker is set to
    unitRows <- vector("list", numUnits)
    for (replication in seq_len(n.reps)) {
      set.seed(splitSeeds[replication])
      if (method == "k-fold") {
        permutation <- sample.int(numObservations)
        foldOffset <- 0L
        for (fold in seq_len(numFolds)) {
          unitRows[[(replication - 1L) * numFolds + fold]] <- sort(
            permutation[foldOffset + seq_len(foldSizes[fold])]
          )
          foldOffset <- foldOffset + foldSizes[fold]
        }
      } else {
        unitRows[[replication]] <- sort(sample.int(numObservations, numTest))
      }
    }

    if (numChunks == 1L) {
      return(list(xbartRunChunk(spec, unitRows, unitSeeds)))
    }
    chunkRows <- lapply(chunkIndices, function(indices) unitRows[indices])
    chunkSeeds <- lapply(
      chunkIndices,
      function(indices) unitSeeds[indices, , drop = FALSE]
    )
    if (workerKind == "fork") {
      # a worker that fails hands back its condition rather than a try-error,
      # so a child that dies without one shows as a NULL result and mclapply's
      # own warning about it is left to signal
      results <- parallel::mclapply(
        seq_len(numChunks),
        function(i) {
          tryCatch(
            xbartRunChunk(spec, chunkRows[[i]], chunkSeeds[[i]]),
            error = function(e) {
              e$call <- NULL
              e
            }
          )
        },
        mc.cores = numChunks,
        mc.preschedule = FALSE
      )
      for (result in results) {
        if (inherits(result, "condition")) {
          stop(conditionMessage(result), call. = FALSE)
        }
        if (!is.list(result) || is.null(result$loss)) {
          stop(
            "a forked worker exited without a result; ",
            "try parallel = \"socket\""
          )
        }
      }
      return(results)
    }
    cluster <- if (workerKind == "cluster") {
      cl
    } else {
      made <- parallel::makeCluster(numChunks)
      on.exit(parallel::stopCluster(made), add = TRUE)
      made
    }
    # passing the namespace function itself serializes it by reference,
    # loading dbarts on the workers without shipping this frame; a worker's
    # own RNGkind() never matters, since every sampler it creates seeds off
    # control's seed slot rather than that worker's stream
    parallel::clusterMap(
      cluster[seq_len(numChunks)],
      xbartRunChunk,
      unitRows = chunkRows,
      unitSeeds = chunkSeeds,
      MoreArgs = list(spec = spec),
      SIMPLIFY = FALSE
    )
  }

  if (verbose) {
    cat(
      "running ",
      numCells,
      " parameter combination",
      if (numCells > 1L) "s" else "",
      " x ",
      numUnits,
      " (replication, fold) unit",
      if (numUnits > 1L) "s" else "",
      " on ",
      numChunks,
      " ",
      workerKind,
      " worker",
      if (numChunks > 1L) "s" else "",
      "\n",
      sep = ""
    )
  }

  # unit-major, cells within; the folds of one replication are contiguous, so
  # the reported loss is their average, as it was when one worker ran every
  # fold of a replication in sequence
  chunkResults <- withPreservedSeed(runUnits())
  unitLoss <- do.call(rbind, lapply(chunkResults, `[[`, "loss"))
  signalChunkWarnings(chunkResults)
  numResults <- ncol(unitLoss)
  lossValues <- matrix(
    apply(
      array(unitLoss, c(numCells, numFolds, n.reps, numResults)),
      c(1L, 3L, 4L),
      mean
    ),
    n.reps * numCells,
    numResults
  )

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
  # order on both the array and the k axis before anything is reported
  if (length(kGrid) > 1L) {
    kOrderInv <- kOrder
    kOrderInv[kOrder] <- seq_along(kOrder)
    result <- result[,, kOrderInv, , , , drop = FALSE]
  }
  # the k or sd axis labels its cells: a fixed cell by its value, at the two
  # significant digits every other axis prints, a modelled cell by the
  # constructor call that rebuilds it
  k <- sd <- kLabels

  varNames <- c("n.trees", if (sdAxis) "sd" else "k", "power", "base")
  dimIncluded <- c(
    TRUE,
    if (drop) length(n.trees) > 1L else TRUE,
    if (drop) length(k) > 1L else TRUE,
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
  # a supplied function keeps its own environment, so a closure keeps what it
  # captured; the list form calls it from the given environment
  if (is.function(loss)) {
    return(loss)
  }
  if (is.list(loss)) {
    lossFunction <- loss[[1L]]
    lossEnv <- loss[[2L]]
    return(function(y.test, testSamples, weights) {
      eval(as.call(list(lossFunction, y.test, testSamples, weights)), lossEnv)
    })
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

## One worker's share of the (replication, fold) units, as the rows each
## unit holds out and the seeds its fits run under (one column per distinct
## tree count, in n.trees order). The predictor store (cuts + codes) is built
## once per chunk; each unit's sampler is a row-subset view over it, so every
## fold bins on the full data's cut grid and no fold re-quantizes the
## predictors. Within a unit every tree count gets a fresh sampler, seeded
## from that unit's row of seeds through the control's seed slot, burned
## n.burn[1] iterations; the remaining parameter cells sweep warm off it with
## n.burn[2] iterations each on the same sampler, sound because the training
## data is unchanged. Chains never carry over between units, whose held-out
## rows the previous training set contained; seeding per sampler creation
## rather than per chunk is what keeps a result independent of how the units
## were distributed, and no worker ever calls set.seed().
## Returns a (units x cells) x numResults matrix, cells in spec$cells order.
xbartRunUnits <- function(spec, unitRows, unitSeeds) {
  data <- spec$data
  cells <- spec$cells
  numCells <- nrow(cells)
  numObservations <- length(data@y)
  hasWeights <- !is.null(data@weights)
  family <- spec$model@family

  # linear and gp leaf priors read raw covariate values, fixed across cells;
  # the handle must own raw for them so each fold view can gather them
  leafPrior <- spec$model@leaf.prior
  leafCovariateColumns <-
    if (is(leafPrior, "dbartsLinearPrior") || is(leafPrior, "dbartsGPPrior")) {
      leafPrior@columns
    } else {
      NULL
    }
  handle <- bartcoreDataHandle(spec$control, data, leafCovariateColumns)

  # the per-fold sigma's dense design is built once per chunk, not per fold
  sigmaDesign <- if (identical(spec$sigmaPerFold, "linear")) {
    sigmaDesignMatrix(data@x)
  }
  foldData <- function(trainRows) {
    if (is.null(spec$sigmaPerFold)) {
      return(data)
    }
    offset <- if (!is.null(data@offset)) data@offset[trainRows]
    y <- data@y[trainRows]
    residual <- if (!is.null(offset)) y - offset else y
    data@sigma <- if (is.null(sigmaDesign)) {
      floorMarginalSigma(sd(residual), residual)
    } else {
      floorSigmaEstimate(
        residualStandardError(
          y,
          sigmaDesign[trainRows, , drop = FALSE],
          if (hasWeights) data@weights[trainRows],
          offset
        ),
        residual
      )
    }
    data
  }

  cellModel <- function(cell) {
    result <- spec$model
    result@tree.prior@power <- spec$power[cells$iPower[cell]]
    result@tree.prior@base <- spec$base[cells$iBase[cell]]
    result@leaf.hyperprior <- spec$kHyperpriors[[cells$iK[cell]]]
    result@prior.scale <- spec$kAnchors[[cells$iK[cell]]]
    result
  }

  numLossResults <- NULL

  # fit every cell against one split, fresh per tree count, warm across the
  # rest; cells are grouped by iTrees, so a single pass reuses each sampler
  # maximally. The view slices y/weights/offset by row and takes its test
  # offset from offset[testRows], so each fold trains and scores on exactly
  # its own rows. treeSeeds holds one caller-drawn seed per distinct tree
  # count, in n.trees order; a fresh sampler consumes the next one through
  # its control's seed slot, which derives its chain's generator from a
  # dbarts generator rather than R's stream.
  sweepCells <- function(testRows, treeSeeds) {
    trainRows <- seq_len(numObservations)[-testRows]
    trainData <- foldData(trainRows)
    y.test <- data@y[testRows]
    weights.test <- if (hasWeights) data@weights[testRows] else NULL

    lossValues <- NULL
    sampler <- NULL
    cellControl <- NULL
    currentTrees <- NA_integer_
    treeIndex <- 0L
    for (cell in seq_len(numCells)) {
      if (
        is.null(sampler) || spec$n.trees[cells$iTrees[cell]] != currentTrees
      ) {
        treeIndex <- treeIndex + 1L
        cellControl <- spec$control
        cellControl@n.trees <- spec$n.trees[cells$iTrees[cell]]
        cellControl@seed <- treeSeeds[treeIndex]
        sampler <- bartcoreSamplerFromHandle(
          handle,
          cellControl,
          cellModel(cell),
          trainData,
          trainRows,
          testRows,
          family
        )
        currentTrees <- spec$n.trees[cells$iTrees[cell]]
        numBurnIn <- spec$n.burn[1L]
      } else {
        bartcoreSetModel(sampler, cellModel(cell), trainData, cellControl)
        numBurnIn <- spec$n.burn[2L]
      }

      samples <- bartcoreRun(sampler, numBurnIn, spec$n.samples)
      # bartcoreRun does not warn on its own (unlike bartcoreSamplerRun); each
      # call here is one complete cell/fold fit, so this is that fit's one
      # warning, not a per-sweep one
      warnOnGPFallback(samples)
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

  results <- vector("list", length(unitRows))
  for (i in seq_along(unitRows)) {
    results[[i]] <- sweepCells(unitRows[[i]], unitSeeds[i, ])
  }

  do.call(rbind, results)
}

## xbartRunUnits with every warning its fits raise captured rather than
## signalled, so a chunk run in this process and one run on a worker, whose
## own warnings would never reach the caller, report the same way. Returns
## the loss matrix, the distinct warnings in the order first raised, and the
## warned-once keys the chunk set. Under options(warn = 2) a warning is left
## to abort the run where it is raised, and an error signals the warnings
## captured before it ahead of propagating.
xbartRunChunk <- function(spec, unitRows, unitSeeds) {
  # a worker takes the caller's warning level, so 'warn = 2' escalates there
  oldWarn <- options(warn = spec$warn)
  on.exit(options(oldWarn), add = TRUE)
  for (key in spec$onceKeys) {
    onceWarnState[[key]] <- TRUE
  }
  captured <- list()
  keys <- character()
  loss <- withCallingHandlers(
    xbartRunUnits(spec, unitRows, unitSeeds),
    warning = function(w) {
      if (getOption("warn") >= 2L) {
        return()
      }
      # the call names this chunk's internals, not the caller's code, and can
      # carry a frame too large to send back from a worker
      w$call <- NULL
      key <- warningKey(w)
      if (key %not_in% keys) {
        keys <<- c(keys, key)
        captured[[length(captured) + 1L]] <<- w
      }
      invokeRestart("muffleWarning")
    },
    error = function(e) {
      for (w in captured) {
        warning(w)
      }
    }
  )
  list(
    loss = loss,
    warnings = captured,
    onceKeys = setdiff(warnedOnceKeys(), spec$onceKeys)
  )
}

## Re-signals the chunks' captured warnings once all units have finished, in
## unit order, each distinct (class, message) pair once, so a warning with a
## fixed message that recurs in every fit is reported once. The chunks'
## warned-once keys are marked in this session.
signalChunkWarnings <- function(chunkResults) {
  for (key in unlist(lapply(chunkResults, `[[`, "onceKeys"))) {
    onceWarnState[[key]] <- TRUE
  }
  captured <- unlist(lapply(chunkResults, `[[`, "warnings"), recursive = FALSE)
  keys <- vapply(captured, warningKey, "")
  for (w in captured[!duplicated(keys)]) {
    warning(w)
  }
  invisible(NULL)
}

warningKey <- function(w) {
  paste(c(class(w), conditionMessage(w)), collapse = "\n")
}

## The k axis, normalized to one leaf hyperprior per grid cell: a numeric
## vector is 0.9-x's fixed grid, a list mixes fixed values with hyperprior
## objects, a bare number or hyperprior is a one-cell grid, and NULL takes
## the response type's own front-door default.
resolveKGrid <- function(k, binary) {
  if (is.null(k)) {
    return(list(resolveLeafHyperprior(NULL, binary = binary)))
  }
  entries <- if (is.list(k)) {
    k
  } else if (is.numeric(k) || is.character(k)) {
    as.list(k)
  } else {
    list(k)
  }
  if (length(entries) == 0L) {
    stop("'k' must name at least one value")
  }
  lapply(entries, resolveKEntry)
}

## The sd axis: a numeric vector of absolute spreads, or a list mixing them
## with invchi() laws, each translated to the leaf hyperprior and anchor a
## leaf prior naming that sd reaches the engine with (resolveLeafPrior).
resolveSdGrid <- function(sd) {
  entries <- if (is.list(sd)) {
    sd
  } else if (is.numeric(sd)) {
    as.list(sd)
  } else {
    list(sd)
  }
  if (length(entries) == 0L) {
    stop("'sd' must name at least one value")
  }
  lapply(entries, function(entry) {
    if (is.function(entry)) {
      entry <- entry()
    }
    entry <- validateLeafSd(entry)
    if (is.null(entry)) {
      stop("'sd' must contain positive numbers and invchi() specifications")
    }
    translated <- resolveLeafPrior(
      new("dbartsNormalPrior", prior.sd = entry),
      binary = FALSE
    )
    c(
      translated,
      label = if (is.numeric(entry)) {
        as.character(signif(entry, 2L))
      } else {
        paste0("invchi(", format(entry@df), ", ", format(entry@scale), ")")
      },
      # most shrunk first: the smallest fixed sd sorts highest
      sort.key = if (is.numeric(entry)) -entry else -Inf
    )
  })
}

## One k grid entry: a positive number fixes k for its cell, a hyperprior
## object models it there. A character entry keeps normal()'s own string
## forms, so "2" and "chi(1.5)" read as they always did.
resolveKEntry <- function(entry) {
  if (is.character(entry)) {
    entry <- normal(entry)@k
  }
  if (is.numeric(entry)) {
    if (length(entry) != 1L || is.na(entry) || entry <= 0.0) {
      stop("'k' must contain only positive values")
    }
    return(newValidated("dbartsFixedHyperprior", k = as.numeric(entry)))
  }
  if (is.function(entry)) {
    entry <- entry()
  }
  if (!is(entry, "dbartsLeafHyperprior")) {
    stop(
      "'k' must contain positive numbers and hyperprior specifications; ",
      "see ?dbartsPriors"
    )
  }
  entry
}

## The k one grid cell builds its leaf model at: the value a fixed cell
## holds, or the hyperprior object a modelled one is drawn under.
kGridValue <- function(entry) {
  if (is(entry, "dbartsFixedHyperprior")) entry@k else entry
}

## Sort key for the sweep order: fixed cells sweep from the most shrunk
## down, and a modelled cell, having no fixed k to place, sweeps last.
kGridSortKey <- function(entry) {
  if (is(entry, "dbartsFixedHyperprior")) entry@k else -Inf
}

## One k axis label. A fixed cell prints its value at the two significant
## digits every grid axis prints; a modelled cell prints the constructor
## call that rebuilds it, so the two are told apart in the dimnames.
kGridLabel <- function(entry) {
  if (is(entry, "dbartsFixedHyperprior")) {
    return(as.character(signif(entry@k, 2L)))
  }
  if (!is(entry, "dbartsChiHyperprior")) {
    # a hyperprior class added without a label here would be reported under
    # some other constructor's name, which no reader could tell from a real
    # one; refused by class instead
    stop(
      "no k axis label for a hyperprior of class \"",
      class(entry),
      "\"; add one to kGridLabel"
    )
  }
  paste0(
    "chi(",
    format(entry@degreesOfFreedom),
    ", ",
    format(entry@scale),
    ")"
  )
}
