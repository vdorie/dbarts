setMethod("initialize", "dbartsControl", function(.Object, ...) {
  .Object <- callNextMethod()

  validObject(.Object)
  .Object
})

# Parse a survival response (docs/design/survival.md): a survival::Surv object
# (recognized by inherits(), so survival need not be imported; right-censoring
# only in v1) or a plain two-column (time, status) matrix or data frame.
# Returns the raw event/censoring time and the 0/1 status vector, or NULL when
# the value is not a survival response. Errors on a non-right Surv (with a
# factor-status hint for "mright"), a non-two-column matrix, non-positive
# times, or a status outside {0, 1}; a Surv-like object with no type attribute
# is treated as right-censored. Shared by the aft ingestion (which logs the
# time) and the discrete-time hazard expander (which keeps the raw time).
parseSurvivalResponse <- function(value) {
  if (inherits(value, "Surv")) {
    type <- attr(value, "type")
    if (identical(type, "mright")) {
      # survival::Surv codes a factor status as multi-state
      stop(
        "multi-state survival responses are not supported; the Surv status ",
        "must be 0/1 or logical, not a factor"
      )
    }
    if (!is.null(type) && type != "right") {
      stop("survival responses support only right-censoring in this version")
    }
    # unclass before extraction so [.Surv (or any classed-matrix method)
    # cannot re-wrap the columns
    value <- unclass(value)
    time <- as.double(value[, 1L])
    status <- as.double(value[, 2L])
  } else if ((is.matrix(value) || is.data.frame(value)) && NCOL(value) == 2L) {
    if (is.data.frame(value)) {
      if (is.factor(value[[2L]])) {
        stop("survival status must be 0/1 or logical, not a factor")
      }
      time <- as.double(value[[1L]])
      status <- as.double(value[[2L]])
    } else {
      value <- unclass(value)
      time <- as.double(value[, 1L])
      status <- as.double(value[, 2L])
    }
  } else {
    return(NULL)
  }
  if (any(!is.finite(time)) || any(time <= 0.0)) {
    stop("survival times must be finite and positive")
  }
  if (any(status != 0.0 & status != 1.0)) {
    stop("survival status must be 0 (censored) or 1 (event)")
  }
  list(time = time, status = status)
}

# Accelerated failure time ingestion: the log event/censoring time as the
# working response and the 0/1 status vector, or NULL. Wraps
# parseSurvivalResponse with the log() transform the AFT engine expects.
extractSurvivalResponse <- function(value) {
  survival <- parseSurvivalResponse(value)
  if (is.null(survival)) {
    return(NULL)
  }
  list(log.time = log(survival$time), status = survival$status)
}

# Discrete-time hazard ingestion (docs/design/survival.md, "Discrete-time
# hazard"): the RAW time and status, the AFT sibling that skips the log()
# transform, for the person-period expander.
extractSurvivalTimes <- parseSurvivalResponse

# Resolve the discrete-time grid and each subject's terminal period from the
# observed times (docs/design/survival.md section 1). `breaks` NULL (the
# default) uses the sorted distinct observed times (surv.bart's convention);
# a length-1 integer bins at the (1:K)/K quantiles (surv.bart's K); a longer
# numeric vector gives explicit interval boundaries b_0 < ... < b_K with
# right-closed intervals (b_{k-1}, b_k], the discSurv convention. Returns the
# representative period times (the right edges, sorted ascending) and each
# subject's terminal period index (1..K). Ties within a period are automatic:
# equal times share a period. findInterval(..., left.open = TRUE) counts grid
# points strictly below t, so a time exactly on grid point g_k lands in period
# k (its own interval's right edge).
resolveHazardGrid <- function(time, breaks) {
  if (is.null(breaks)) {
    periods <- sort(unique(time))
  } else {
    breaks <- as.double(breaks)
    if (anyNA(breaks)) {
      stop("'breaks' must not contain missing values")
    }
    if (length(breaks) == 1L) {
      K <- as.integer(breaks)
      if (is.na(K) || K < 1L) {
        stop(
          "'breaks' as a single value must be a positive integer period count"
        )
      }
      periods <- unique(as.double(
        quantile(time, probs = seq_len(K) / K, names = FALSE)
      ))
    } else {
      if (is.unsorted(breaks, strictly = TRUE)) {
        stop("'breaks' boundaries must be strictly increasing")
      }
      if (any(time <= breaks[1L]) || any(time > breaks[length(breaks)])) {
        stop(
          "every survival time must lie within the 'breaks' boundaries ",
          "(b_1, b_K]; widen the outer boundaries to cover the data"
        )
      }
      periods <- breaks[-1L]
    }
  }
  terminal <- findInterval(time, periods, left.open = TRUE) + 1L
  list(periods = periods, terminalPeriod = terminal)
}

# The person-period expander (docs/design/survival.md section 1): a subject
# i observed to time_i (event or censoring) becomes its at-risk rows, one per
# period k = 1..t_i, each carrying x_i, the ordinal period column (appended
# LAST), and the binary indicator y_ik = status_i * 1{k = t_i}. A censored
# subject's rows are all zero (right-censoring is pure data shape). Offsets and
# weights replicate per subject and follow the chosen binary family's policy
# downstream. Returns the expanded design, the binary response, the period
# grid (for the $periods marker), and the replicated offset/weights. The N'
# row guard (max.rows) refuses an over-fine grid, naming the coarsening levers.
expandDiscreteTimeHazard <- function(
  x,
  time,
  status,
  breaks = NULL,
  max.rows = 1e7,
  offset = NULL,
  weights = NULL
) {
  n <- length(time)
  grid <- resolveHazardGrid(time, breaks)
  periods <- grid$periods
  terminal <- grid$terminalPeriod

  Nprime <- sum(terminal)
  if (Nprime > max.rows) {
    stop(
      "person-period expansion would create ",
      Nprime,
      " rows, over the cap of ",
      max.rows,
      "; coarsen the time grid with 'breaks' (a boundary vector or an ",
      "integer period count) or raise 'max.rows'"
    )
  }

  # subject-major, period ascending within subject: the period-1 rows are the
  # subjects in order (survivalProbabilities reconstructs training covariates
  # from them), and sequence() supplies the within-subject period counter
  subjectOf <- rep.int(seq_len(n), terminal)
  periodOf <- sequence(terminal)
  y <- as.double(periodOf == terminal[subjectOf] & status[subjectOf] == 1.0)

  xExpanded <- if (is.data.frame(x)) {
    x[subjectOf, , drop = FALSE]
  } else {
    x <- as.matrix(x)
    x[subjectOf, , drop = FALSE]
  }
  xExpanded <- appendHazardPeriodColumn(xExpanded, periodOf)

  result <- list(x = xExpanded, y = y, periods = periods)
  if (!is.null(offset)) {
    offset <- as.double(offset)
    if (length(offset) == 1L) {
      offset <- rep_len(offset, n)
    }
    result$offset <- offset[subjectOf]
  }
  if (!is.null(weights)) {
    weights <- as.double(weights)
    if (length(weights) == 1L) {
      weights <- rep_len(weights, n)
    }
    result$weights <- weights[subjectOf]
  }
  result
}

# Append the ordinal period column (named "period") as the LAST column, the
# fixed convention both the hazard fit and its by-hand binary reduction target
# rely on. A named matrix keeps its names; an unnamed one stays unnamed so
# dbartsData assigns its usual defaults (the reduction target sees the same).
appendHazardPeriodColumn <- function(x, period) {
  if (is.data.frame(x)) {
    x[["period"]] <- period
    return(x)
  }
  named <- !is.null(colnames(x))
  out <- cbind(x, period)
  if (named) {
    colnames(out)[ncol(out)] <- "period"
  } else {
    colnames(out) <- NULL
  }
  out
}

## every slot below is passed explicitly to newValidated, so this
## function's own defaults are what apply, not A_class.R's prototype;
## only `binary` and `call`, which this constructor never sets, fall
## through to it. n.threads is deliberately one of the explicit ones:
## this default probes guessNumCores(), while the prototype's is the
## conservative n.threads = 1L for a bare new("dbartsControl").
dbartsControl <- function(
  verbose = FALSE,
  keepTrainingFits = TRUE,
  useQuantiles = FALSE,
  keepTrees = FALSE,
  storage = c("double", "single"),
  n.samples = NA_integer_,
  n.cuts = 100L,
  n.burn = 200L,
  n.trees = 75L,
  n.chains = 4L,
  n.threads = dbarts::guessNumCores(),
  n.thin = 1L,
  printEvery = 100L,
  printCutoffs = 0L,
  rngSeed = NA_integer_,
  updateState = TRUE
) {
  storage <- match.arg(storage)
  newValidated(
    "dbartsControl",
    verbose = as.logical(verbose),
    keepTrainingFits = as.logical(keepTrainingFits),
    useQuantiles = as.logical(useQuantiles),
    keepTrees = as.logical(keepTrees),
    storage = storage,
    n.samples = coerceOrError(n.samples, "integer"),
    n.cuts = coerceOrError(n.cuts, "integer"),
    n.burn = coerceOrError(n.burn, "integer"),
    n.trees = coerceOrError(n.trees, "integer"),
    n.chains = coerceOrError(n.chains, "integer"),
    n.threads = coerceOrError(n.threads, "integer"),
    n.thin = coerceOrError(n.thin, "integer"),
    printEvery = coerceOrError(printEvery, "integer"),
    printCutoffs = coerceOrError(printCutoffs, "integer"),
    rngSeed = coerceOrError(rngSeed, "integer"),
    updateState = as.logical(updateState)
  )
}

validateArgumentsInEnvironment <- function(
  envir,
  func,
  funcName,
  control,
  verbose,
  n.samples,
  sigma
) {
  controlIsMissing <- missing(control)

  if (!controlIsMissing) {
    if (!inherits(control, "dbartsControl")) {
      stop(
        "'control' argument must be of class dbartsControl; use ",
        "dbartsControl() function to create"
      )
    }
    envir$control <- control
  }

  if (!missing(verbose)) {
    if (!is.logical(verbose) || is.na(verbose)) {
      stop("'verbose' argument to ", funcName, " must be TRUE/FALSE")
    }
  } else if (!controlIsMissing) {
    envir$verbose <- control@verbose
  }

  if (!missing(n.samples)) {
    tryCatch(n.samples <- as.integer(n.samples), warning = function(e) {
      stop(
        "'n.samples' argument to ",
        funcName,
        " must be coercible to integer type"
      )
    })
    if (length(n.samples) != 1L) {
      stop("'n.samples' must be of length 1")
    }
    if (is.null(n.samples)) {
      stop("'n.samples' argument to ", funcName, " cannot be NULL")
    }
    if (is.na(n.samples) || n.samples < 0L) {
      stop(
        "'n.samples' argument to ",
        funcName,
        " must be a non-negative integer"
      )
    }
    envir$control@n.samples <- n.samples
  } else if (controlIsMissing || is.na(control@n.samples)) {
    envir$control@n.samples <- formals(func)[["n.samples"]]
  }

  if (!missing(sigma) && !is.na(sigma)) {
    tryCatch(sigma <- as.double(sigma), warning = function(e) {
      stop(
        "'sigma' argument to ",
        funcName,
        " must be coercible to numeric type"
      )
    })
    if (length(sigma) != 1L) {
      stop("'sigma' must be of length 1")
    }
    if (is.null(sigma) || sigma <= 0.0) {
      stop("'sigma' argument to ", funcName, " must be positive")
    }

    envir$sigma <- sigma
  }
}

dbarts <- function(
  formula,
  data,
  test,
  subset,
  weights,
  offset,
  offset.test = offset,
  verbose = FALSE,
  n.samples = 800L,
  tree.prior = cgm,
  node.prior = normal,
  resid.prior = chisq,
  resid.dist = gaussian,
  proposal.probs = c(birth_death = 0.5, swap = 0.1, change = 0.4, birth = 0.5),
  monotone = NULL,
  variance = NULL,
  n.trees.variance = 40L,
  power.variance = NULL,
  base.variance = NULL,
  control = dbarts::dbartsControl(),
  sigma = NA_real_,
  seed = NA_integer_,
  factors = c("categorical", "indicators"),
  family = c(
    "auto",
    "gaussian",
    "probit",
    "logistic",
    "aft",
    "ordinal",
    "nbinom",
    "hazard",
    "hazard.probit",
    "hazard.logistic",
    "hurdle.lognormal",
    "twopart"
  ),
  missing = c("incorporate", "error"),
  dispersion = NA_real_,
  breaks = NULL,
  max.rows = 1e7
) {
  matchedCall <- match.call()

  evalEnv <- parent.frame(1L)

  family <- match.arg(family)
  # hurdle.lognormal / twopart (docs/design/hurdle.md): the alias resolves to
  # the canonical token immediately, so every downstream message and any
  # packaged $family reads "hurdle.lognormal" regardless of which spelling
  # was requested
  if (identical(family, "twopart")) {
    family <- "hurdle.lognormal"
  }
  if (identical(family, "hurdle.lognormal")) {
    # a hurdle fit composes TWO independent samplers (an occupancy probit and
    # a positive-part gaussian, docs/design/hurdle.md section 2); dbarts()
    # returns exactly one sampler and cannot express that composition - only
    # bart2() (bart2Hurdle) builds it
    stop(
      "family \"hurdle.lognormal\" fits two component samplers and is only ",
      "available through bart2()"
    )
  }

  # survival response ingestion (docs/design/survival.md): a survival::Surv
  # object or an explicit family = "aft" with a two-column (time, status)
  # response fits the AFT log-normal model. The matrix (x.train, y.train)
  # interface is supported here, the response being the second positional
  # argument; the log event/censoring time replaces it and the status rides
  # the control attribute the bartcore survival family reads. Formula-LHS Surv
  # and a subset argument are deferred to a later surface pass.
  survivalStatus <- NULL
  directResponse <- !is.formula(formula) &&
    !inherits(formula, "dbartsData") &&
    !inherits(formula, "dgCMatrix") &&
    !missing(data)
  responseIsSurv <- directResponse && inherits(data, "Surv")
  hazardTokens <- c("hazard", "hazard.probit", "hazard.logistic")
  # a Surv response declares the model, so it auto-dispatches to aft from
  # "auto"; an explicit hazard token selects the discrete-time model instead
  # (the guard whitelist admits it, docs/design/survival.md section 2). Any
  # other explicit family with a Surv response is a conflict, never a silent
  # override.
  if (responseIsSurv && family %not_in% c("auto", "aft", hazardTokens)) {
    stop(
      "a survival (Surv) response cannot be fit with family \"",
      family,
      "\"; use family \"aft\", \"hazard\", or \"auto\""
    )
  }

  # discrete-time hazard ingestion (docs/design/survival.md, "Discrete-time
  # hazard"): person-period-expand (x, time, status) into an ordinary binary
  # (X', y') design and REMAP the hazard token to its underlying binary link
  # BEFORE any family-keyed switch runs (node.scale, control@binary,
  # fixedUnitScale, the weight policy). The engine, bridge, and ResponseModels
  # then see an ordinary probit/logistic fit; the hazard provenance survives
  # only as the period grid, parked on the control attribute the packaging
  # reads into the $periods marker (the bartcore.survival -> $status
  # precedent). No status vector or attribute reaches C++ - the censoring is
  # baked into y'.
  hazardPeriods <- NULL
  if (family %in% hazardTokens) {
    if (!directResponse) {
      stop(
        "discrete-time hazard fits currently use the matrix interface - ",
        "dbarts(x.train, y.train) or bart2(x.train, y.train) with a ",
        "survival::Surv or two-column (time, status) response"
      )
    }
    survival <- extractSurvivalTimes(data)
    if (is.null(survival)) {
      stop(
        "family \"",
        family,
        "\" needs a survival::Surv or two-column (time, status) response"
      )
    }
    if (!missing(subset)) {
      stop("survival responses do not support 'subset' in this version")
    }
    if (!missing(test)) {
      stop(
        "discrete-time hazard fits do not take a 'test' set; expand test ",
        "subjects with survivalProbabilities(fit, times, newdata = )"
      )
    }
    expansion <- expandDiscreteTimeHazard(
      formula,
      survival$time,
      survival$status,
      breaks = breaks,
      max.rows = max.rows,
      offset = if (missing(offset)) NULL else offset,
      weights = if (missing(weights)) NULL else weights
    )
    matchedCall$formula <- expansion$x
    matchedCall$data <- expansion$y
    matchedCall$test <- NULL
    matchedCall$offset.test <- NULL
    if (!is.null(expansion$offset)) {
      matchedCall$offset <- expansion$offset
    }
    if (!is.null(expansion$weights)) {
      matchedCall$weights <- expansion$weights
    }
    hazardPeriods <- expansion$periods
    # the remap: the engine-facing family is now an ordinary binary link
    family <- if (identical(family, "hazard.logistic")) "logistic" else "probit"
    # the survival response is consumed; do not let the aft block fire on it
    responseIsSurv <- FALSE
  }
  # aft is reachable through the direct-response form, or through the internal
  # rbart channel, which pre-sets the status on control@bartcore.survival and
  # passes a ready dbartsData; every other indirect route (the public formula
  # interface) is refused up front, before the response is materialized,
  # rather than failing hostilely downstream
  if (
    family == "aft" &&
      !directResponse &&
      is.null(attr(control, "bartcore.survival"))
  ) {
    stop(
      "survival (aft) fits currently use the matrix interface - ",
      "dbarts(x.train, y.train) or bart2(x.train, y.train) with a ",
      "survival::Surv or two-column (time, status) response"
    )
  }
  if (directResponse && (family == "aft" || responseIsSurv)) {
    survival <- extractSurvivalResponse(data)
    if (is.null(survival)) {
      stop(
        "family \"aft\" needs a survival::Surv or two-column ",
        "(time, status) response"
      )
    }
    if (!missing(subset)) {
      stop("survival responses do not support 'subset' in this version")
    }
    family <- "aft"
    matchedCall$data <- survival$log.time
    survivalStatus <- survival$status
  }

  validateCall <- redirectCall(
    matchedCall,
    quoteInNamespace(validateArgumentsInEnvironment)
  )
  validateCall <- addCallArgument(validateCall, 1L, sys.frame(sys.nframe()))
  validateCall <- addCallArgument(validateCall, 2L, dbarts::dbarts)
  validateCall <- addCallArgument(validateCall, 3L, "dbarts")
  eval(validateCall, evalEnv, getNamespace("dbarts"))

  if (length(control@call) == 1L && control@call == call("NA")) {
    control@call <- matchedCall
  }
  control@verbose <- verbose
  # a convenience mirror of dbartsControl(rngSeed = ), as the wrappers expose;
  # an explicit seed overrides the control's, NA leaves it untouched
  seed <- coerceOrError(seed, "integer")
  if (!is.na(seed)) {
    control@rngSeed <- seed
  }

  dataCall <- redirectCall(matchedCall, quoteInNamespace(dbartsData))
  data <- eval(dataCall, evalEnv)
  # cat(
  #   "x address after dbartsData call: ",
  #   .Call("dbarts_getPointerAddress", data@x),
  #   "\n",
  #   sep = ""
  # )

  data@n.cuts <- rep_len(control@n.cuts, ncol(data@x))
  data@sigma <- sigma

  # a factor/logical/character response declares a classification model. The
  # single-forest engine here fits only the 2-level (probit) case; 3+ levels
  # are multinomial, which only bart2(family = "multinomial") implements. A
  # numeric response takes the historic 0/1-vs-continuous path unchanged.
  # dbarts() is also reached anonymously through bart(), which never sets an
  # explicit family; see resolveClassificationFamily's doc comment for why
  # its auto-branch message lists every single-forest entry point instead of
  # naming itself. probit/logistic on a 2-level categorical response proceed
  # as binary.
  family <- resolveClassificationFamily(
    data,
    family,
    "dbarts()/rbart_vi()/bart()/xbart",
    c("gaussian", "aft", "nbinom"),
    splitMultinomialMessage = TRUE,
    allowOrdinal = TRUE
  )
  # ordinal (cumulative probit, docs/design/ordinal.md): a single-forest fixed-
  # unit-scale model like probit, but K-level. Recode the response to the
  # 1-based category codes the engine reads, and attach K on the control
  # attribute the bridge (C2) reads to select OrdinalResponse (the
  # bartcore.survival precedent below). The resolved ordered levels ride the
  # data object for the round-trip.
  if (identical(family, "ordinal")) {
    ordinal <- resolveOrdinalResponse(data)
    data@y <- ordinal$y
    data@response.levels <- ordinal$levels
    attr(control, "bartcore.n.categories") <- ordinal$K
  } else if (identical(family, "nbinom")) {
    # negative-binomial counts (docs/design/negative-binomial.md section 4): the
    # count response has no unambiguous class, so "nbinom" is never auto - only
    # explicit. y must be a non-negative integer count (the NB pmf has zero mass
    # off the integers and the grid kernel's count histogram presumes integer
    # y), validated here beside the binary 0/1 test.
    y <- data@y
    if (anyNA(y) || any(y < 0) || any(y != round(y))) {
      stop("family \"nbinom\" requires a non-negative integer (count) response")
    }
    # the dispersion r: NA (the default) estimates it on the capped integer grid;
    # a supplied value FIXES it and must be a positive integer (v1 ships the
    # exact integer envelope, section 2). The C bridge reads the resolved spec
    # off the control attribute the sampler build attaches below: a positive
    # value fixes r, a non-positive value estimates it on the grid.
    dispersionSpec <- resolveDispersion(dispersion)
    attr(control, "bartcore.dispersion") <- dispersionSpec
  } else if (data@response.type == "numeric") {
    uniqueResponses <- unique(data@y)
    responseIsBinary <- length(uniqueResponses) == 2 &&
      all(sort(uniqueResponses) == c(0, 1))
    if (family == "auto") {
      family <- if (responseIsBinary) "probit" else "gaussian"
    } else if (family != "gaussian" && family != "aft" && !responseIsBinary) {
      # gaussian on a 0/1 response is a legitimate request; the binary
      # families need latent-variable coding. aft fits continuous log-times.
      stop("family \"", family, "\" requires a response coded 0/1")
    }
  }
  # aft draws sigma and rescales like gaussian; only the binary families are
  # latent-variable models on a fixed unit scale
  control@binary <- family %in% c("probit", "logistic")
  # ordinal (cumulative probit) shares probit's fixed unit latent scale - sigma
  # fixed at 1, resid.prior fixed(1), no sigma estimate, node.scale 3.0 - but is
  # NOT binary: the bridge selects it by the bartcore.n.categories attribute
  # (not control@binary), and it reports K category levels. nbinom (counts) is
  # likewise a fixed-unit-scale family (sigma fixed at 1, the log-odds latent
  # entering kappa directly, docs/design/negative-binomial.md section 1),
  # selected by the bartcore.dispersion attribute. fixedUnitScale covers all
  # three families wherever the unit-scale handling matters.
  fixedUnitScale <- control@binary ||
    identical(family, "ordinal") ||
    identical(family, "nbinom")

  # binary weight policy, enforced here in the R layer (the bridge keeps the
  # same checks as a backstop for direct-API consumers): a probit has no
  # tractable weighted latent-variable form and is refused, except that
  # weights identically 1 are the unweighted likelihood and are treated as
  # absent (SuperLearner-style callers pass obsWeights = rep(1, n)
  # unconditionally); a logistic model treats weights as observation counts
  # (its Polya-Gamma latent is a sum of per-copy draws), so they must be
  # positive integers. Gaussian weights, including a gaussian fit of a 0/1
  # response, are unrestricted.
  if (!is.null(data@weights)) {
    if (family == "probit") {
      if (all(data@weights == 1)) {
        data@weights <- NULL
      } else {
        stop(
          "probit models do not support weights: a weighted probit has no ",
          "tractable latent-variable form. Integer count weights can be ",
          "fit exactly with family = \"logistic\"; for continuous weights, ",
          "model the latents directly."
        )
      }
    }
    if (family == "logistic") {
      w <- data@weights
      if (anyNA(w) || any(w <= 0) || any(w != round(w))) {
        stop(
          "logistic weights are observation counts and must be positive ",
          "integers; drop zero-count rows, and use a gaussian model for ",
          "continuous weights."
        )
      }
    }
    # ordinal refuses weights for probit's reason (a weighted truncated-normal
    # latent likelihood is not a coherent model), with the same all-ones-are-
    # absent courtesy the probit path extends to SuperLearner-style callers
    if (family == "ordinal") {
      if (all(data@weights == 1)) {
        data@weights <- NULL
      } else {
        stop(
          "ordinal models do not support weights: a weighted truncated-normal ",
          "latent likelihood is not a coherent model."
        )
      }
    }
    # nbinom refuses weights in v1 (docs/design/negative-binomial.md section 4):
    # the usual count "weight" is exposure, which belongs in the offset as a
    # log-exposure term, not in observation replication. The R refusal mirrors
    # the bridge's backstop, with the probit all-ones-are-absent courtesy.
    if (family == "nbinom") {
      if (all(data@weights == 1)) {
        data@weights <- NULL
      } else {
        stop(
          "nbinom (count) models do not support weights: exposure belongs in ",
          "the offset as a log-exposure term."
        )
      }
    }
  }

  if (is.na(data@sigma) && !fixedUnitScale) {
    tryResult <- tryCatch(
      data@sigma <- estimateSigmaFromLinearModel(data),
      error = function(e) e
    )
    if (inherits(tryResult, "error")) {
      stop("unable to obtain a starting estimate of sigma; provide one instead")
    }
  }

  # bart will passthrough with offset == something no matter what, which we
  # can NULL out; a latent-scale fixed-unit family (binary, ordinal) keeps its
  # zero offset as the meaningful reference, as probit always has
  if (!fixedUnitScale && !is.null(data@offset) && all(data@offset == 0.0)) {
    data@offset <- NULL
  }
  if (
    !fixedUnitScale &&
      !is.null(data@offset.test) &&
      all(data@offset.test == 0.0)
  ) {
    data@offset.test <- NULL
  }
  # cat(
  #   "x address after updating data: ",
  #   .Call("dbarts_getPointerAddress", data@x),
  #   "\n",
  #   sep = ""
  # )

  # resolve the monotone spec here, where its argument is a forced value (a
  # wrapper forwarding it through ... would otherwise reach parsePriors as an
  # unevaluated ...-reference); the resolved direction vector is injected below
  monotoneDirections <- resolveMonotone(monotone, data)

  parsePriorsCall <- redirectCall(matchedCall, quoteInNamespace(parsePriors))
  parsePriorsCall <- setDefaultsFromFormals(
    parsePriorsCall,
    formals(dbarts),
    "tree.prior",
    "node.prior",
    "resid.prior",
    "resid.dist"
  )
  parsePriorsCall$control <- control
  parsePriorsCall$data <- data
  parsePriorsCall$monotone <- monotoneDirections
  parsePriorsCall$parentEnv <- evalEnv

  if (fixedUnitScale) {
    if (any(names(parsePriorsCall) == "resid.prior")) {
      parsePriorsCall[[which(
        names(parsePriorsCall) == "resid.prior"
      )]] <- quote(fixed(1))
    } else {
      parsePriorsCall[[length(parsePriorsCall) + 1L]] <- quote(fixed(1))
      names(parsePriorsCall)[length(parsePriorsCall)] <- "resid.prior"
    }
  }
  priors <- eval(parsePriorsCall)

  # A monotone constraint restricts the forest to birth/death proposals (change
  # and swap would need a > 2-D constrained integral, monotone.md section 5): a
  # defaulted proposal.probs is forced to birth/death-only, an explicit
  # non-default one conflicts and errors.
  if (!is.null(monotoneDirections)) {
    defaultProbs <- c(birth_death = 0.5, swap = 0.1, change = 0.4, birth = 0.5)
    if (!isTRUE(all.equal(proposal.probs[names(defaultProbs)], defaultProbs))) {
      stop(
        "'monotone' forces birth/death-only proposals; a non-default ",
        "'proposal.probs' cannot be honored under the constraint"
      )
    }
    proposal.probs <- c(birth_death = 1, swap = 0, change = 0, birth = 0.5)
  }

  model <- newValidated(
    "dbartsModel",
    priors$tree.prior,
    priors$node.prior,
    priors$node.hyperprior,
    priors$resid.prior,
    proposal.probs = proposal.probs,
    family = family,
    # the logistic scale is probit's default widened by the logistic
    # latent's standard deviation, pi / sqrt(3)
    node.scale = switch(
      family,
      gaussian = 0.5,
      aft = 0.5,
      probit = 3.0,
      # ordinal reuses probit's latent scale (docs/design/ordinal.md section 2,
      # scheme C: the K = 2 anchor is probit exactly)
      ordinal = 3.0,
      # nbinom's psi is a log-odds, so it reuses logistic's node.scale
      # (docs/design/negative-binomial.md section 1)
      nbinom = pi * sqrt(3.0),
      logistic = pi * sqrt(3.0)
    )
  )

  # Student-t residuals (docs/design/robust-errors.md): only a continuous
  # gaussian response carries them (the binary families and aft have their own
  # latent scale), refused here R-side to match the C bridge's backstop. The
  # resolved degrees of freedom ride the model's resid.df attribute the bridge
  # reads - the bartcore.survival precedent above - absent for the Gaussian law.
  if (is(priors$resid.dist, "dbartsStudentDist") && family != "gaussian") {
    stop(
      "student residuals require a continuous gaussian response; family \"",
      family,
      "\" has its own fixed error scale"
    )
  }
  residDf <- residDistDf(priors$resid.dist)
  if (!is.null(residDf)) {
    attr(model, "resid.df") <- residDf
  }

  # the resolved per-column monotone directions ride the model attribute the C
  # bridge reads into SamplerOptions.monotoneDirections (the resid.df precedent)
  if (!is.null(monotoneDirections)) {
    attr(model, "monotone") <- monotoneDirections
  }

  # the AFT survival family reads its per-observation status off this control
  # attribute (the bartcore.groups precedent); the C bridge validates it
  if (!is.null(survivalStatus)) {
    if (length(survivalStatus) != length(data@y)) {
      stop("survival status length does not match the response")
    }
    attr(control, "bartcore.survival") <- survivalStatus
  }
  # the discrete-time hazard marker (docs/design/survival.md section 4): the
  # period grid, parked here for packageBartResults to read into $periods. The
  # C bridge never reads this attribute (unlike bartcore.survival), so a hazard
  # fit's draw stream is byte-identical to the by-hand binary fit's - the
  # reduction gate (benchmarks/R/hazard-reduction.R).
  if (!is.null(hazardPeriods)) {
    attr(control, "bartcore.hazard.periods") <- hazardPeriods
  }

  # the heteroscedastic variance forest (docs/design/heteroscedastic.md): a
  # `variance` selector installs a second forest modeling s^2(x); its config
  # rides the control attribute the C bridge reads. Gaussian + constant leaf
  # only (the C factory refuses otherwise; a friendly R check for the family).
  varianceColumns <- resolveVarianceColumns(variance, data)
  if (!is.null(varianceColumns)) {
    if (family != "gaussian") {
      stop(
        "a variance forest requires family = \"gaussian\"; the latent families ",
        "own the precision channel it routes through"
      )
    }
    if (!is.null(monotoneDirections)) {
      stop("a variance forest is not supported with monotone constraints")
    }
    allColumns <- setequal(varianceColumns, seq_len(ncol(data@x)))
    attr(control, "bartcore.variance") <- list(
      n.trees = coerceOrError(n.trees.variance, "integer"),
      base = if (is.null(base.variance)) {
        model@tree.prior@base
      } else {
        as.double(base.variance)
      },
      power = if (is.null(power.variance)) {
        model@tree.prior@power
      } else {
        as.double(power.variance)
      },
      columns = if (allColumns) NULL else as.integer(varianceColumns)
    )
  }

  result <- new("dbartsSampler", control, model, data)
  # cat(
  #   "x address after creating sampler: ",
  #   .Call("dbarts_getPointerAddress", data@x),
  #   "\n",
  #   sep = ""
  # )
  # cat(
  #   "x address in sampler$data: ",
  #   .Call("dbarts_getPointerAddress", result$data@x),
  #   "\n",
  #   sep = ""
  # )
  result
}

# Coerces a warm-start donor (a sampler, a bart fit with a kept sampler, or a
# raw state) to the stored "bartcoreState" its forests are read from.
warmStartState <- function(donor) {
  if (inherits(donor, "bartcoreState")) {
    return(donor)
  }
  sampler <-
    if (inherits(donor, "dbartsSampler")) {
      donor
    } else if (inherits(donor, "bart") && !is.null(donor$fit)) {
      donor$fit
    } else {
      stop(
        "'warm.start' must be a dbarts sampler, a bart fit made with ",
        "keepSampler = TRUE, or a bartcore state"
      )
    }
  sampler$storeState()
  if (is.null(sampler$state)) {
    stop("warm-start donor has no stored state")
  }
  sampler$state
}

# Draws from the BART prior (issue #31): repeatedly redraws trees and node
# parameters on a private sampler and evaluates the resulting forest at
# x.test, for calibrating priors before fitting - e.g. the prior
# distribution of a treatment effect via f(x1) - f(x0). Never touches the
# caller's sampler: a fresh construction gets its own external pointer, and
# while it borrows the caller's data object, nothing here mutates it.
samplePriorPredictive <- function(
  sampler,
  x.test = NULL,
  n.samples = 200L,
  type = c("ev", "ppd"),
  offset.test = NULL,
  n.threads = sampler$control@n.threads
) {
  if (!inherits(sampler, "dbartsSampler")) {
    stop("'sampler' must inherit from dbartsSampler")
  }
  type <- match.arg(type)
  n.samples <- as.integer(n.samples)[1L]
  if (is.na(n.samples) || n.samples <= 0L) {
    stop("'n.samples' must be a positive integer")
  }

  # a fresh sampler, not sampler$copy(): copy() installs the caller's saved
  # state - including the engine RNG - so successive calls would replay one
  # frozen stream. Fresh creation seeds the chain RNGs from R's stream when
  # control@rngSeed is NA (or pins them when it is set), giving independent
  # draws across calls by default with set.seed governing reproducibility.
  # The prior draws overwrite all tree state, so no donor state is needed.
  # keepTrees is forced off: it makes predict() serve saved posterior
  # samples instead of the live trees this function just drew.
  newControl <- sampler$control
  newControl@keepTrees <- FALSE
  draw <- dbartsSampler$new(newControl, sampler$model, sampler$data)

  xt <- if (is.null(x.test)) extract(draw, "predictors") else x.test
  responseIsBinary <- draw$control@binary

  sigmaDraws <- NULL
  if (type == "ppd" && !responseIsBinary) {
    residPrior <- draw$model@resid.prior
    if (inherits(residPrior, "dbartsChiSqPrior")) {
      # reported-scale scaled-inverse-chi-squared, matching the engine's own
      # sigma-prior calibration (P(sigma < sigest) == quantile)
      sigest <- draw$data@sigma
      df <- residPrior@df
      rawScale <- qchisq(1 - residPrior@quantile, df) / df
      sigmaDraws <- sqrt(df * sigest^2 * rawScale / rchisq(n.samples, df))
    } else if (inherits(residPrior, "dbartsFixedPrior")) {
      # no distributional uncertainty in sigma; getSigmas() already reports
      # the fixed value on the original scale
      sigmaDraws <- rep_len(draw$getSigmas()[1L], n.samples)
    } else {
      stop(
        "samplePriorPredictive does not support residual variance prior ",
        "class '",
        class(residPrior)[1L],
        "'"
      )
    }
  }

  results <- vector("list", n.samples)
  for (i in seq_len(n.samples)) {
    draw$sampleTreesFromPrior(updateState = FALSE)
    draw$sampleNodeParametersFromPrior(updateState = FALSE)
    fit <- draw$predict(xt, offset.test, n.threads)
    # multi-chain samplers draw an independent prior stream per chain; prior
    # draws are chain-free, so only the first chain's stream is kept
    if (length(dim(fit)) > 1L) {
      fit <- fit[, 1L]
    }
    results[[i]] <- fit
  }
  result <- do.call(rbind, results)

  if (responseIsBinary) {
    result <- probabilityFromLatents(result, list(family = draw$model@family))
  }

  if (type == "ppd") {
    if (responseIsBinary) {
      result <- matrix(
        rbinom(length(result), 1L, result),
        nrow(result),
        ncol(result)
      )
    } else {
      result <- result +
        matrix(
          rnorm(length(result), 0, rep(sigmaDraws, ncol(result))),
          nrow(result),
          ncol(result)
        )
    }
  }

  result
}


dbartsSampler <- setRefClass(
  "dbartsSampler",
  fields = list(
    pointer = "externalptr",
    control = "dbartsControl",
    model = "dbartsModel",
    data = "dbartsData",
    state = "ANY" # is either a list of states, or a promise to evaluate
  ),
  methods = list(
    initialize = function(control, model, data, ...) {
      if (!inherits(control, "dbartsControl")) {
        stop("'control' must inherit from dbartsControl")
      }
      if (!inherits(model, "dbartsModel")) {
        stop("'model' must inherit from dbartsModel")
      }
      if (!inherits(data, "dbartsData")) {
        stop("'data' must inherit from dbartsData")
      }
      .self$control <- control
      .self$model <- model
      .self$data <- data

      # "auto" (a hand-built model) keeps the bridge's own dispatch
      .self$pointer <- .Call(
        C_dbarts_bartcore_create,
        .self$control,
        .self$model,
        .self$data,
        if (model@family == "auto") "" else model@family
      )
      # materialized lazily on first access (forcing it before saveRDS
      # captures the sampler), or eagerly by storeState / updateState runs.
      # A deserialized object can force the promise after its pointer has
      # died; that must yield NULL - no stored state - not a C error
      delayedAssign(
        "state",
        {
          if (
            control@updateState &&
              .Call(C_dbarts_bartcore_isValidPointer, pointer)
          ) {
            .Call(C_dbarts_bartcore_storeState, pointer)
          } else {
            NULL
          }
        },
        eval.env = as.environment(.self),
        assign.env = as.environment(.self)
      )

      callSuper(...)
    },
    run = function(
      numBurnIn,
      numSamples,
      updateState = NA,
      n.threads = control@n.threads
    ) {
      "Runs the posterior sampler and returns a list with the results."
      if (missing(numBurnIn)) {
        numBurnIn <- NA_integer_
      }
      if (missing(numSamples)) {
        numSamples <- NA_integer_
      }

      samples <- bartcoreSamplerRun(.self, numBurnIn, numSamples)
      if (
        (is.na(updateState) && control@updateState == TRUE) ||
          identical(updateState, TRUE)
      ) {
        storeState()
      }
      if (is.null(samples)) {
        return(invisible(NULL))
      }
      samples
    },
    sampleTreesFromPrior = function(updateState = NA) {
      "Draws tree structure from prior"
      ptr <- getPointer()
      .Call(C_dbarts_bartcore_sampleTreesFromPrior, ptr)

      if (
        (is.na(updateState) && control@updateState == TRUE) ||
          identical(updateState, TRUE)
      ) {
        storeState(ptr)
      }

      invisible(NULL)
    },
    sampleNodeParametersFromPrior = function(updateState = NA) {
      "Draws end node parameters from prior; does not update tree structure."
      ptr <- getPointer()
      .Call(C_dbarts_bartcore_sampleNodeParametersFromPrior, ptr)

      if (
        (is.na(updateState) && control@updateState == TRUE) ||
          identical(updateState, TRUE)
      ) {
        storeState(ptr)
      }

      invisible(NULL)
    },
    growFromRoot = function(n.sweeps = 2L, updateState = NA) {
      "Builds an initial forest by XBART-style grow-from-root (He, Yalov and Hahn 2019) as a warm start, running n.sweeps grow sweeps in place; the exact MCMC sampler owns the forest once run() begins. Constant-leaf models only. See ?dbartsSampler."
      if (
        is(model@node.prior, "dbartsLinearPrior") ||
          is(model@node.prior, "dbartsGPPrior")
      ) {
        stop(
          "grow-from-root warm start is only available for the constant-leaf ",
          "model; linear and gp node priors initialize with ",
          "sampleTreesFromPrior instead"
        )
      }
      n.sweeps <- as.integer(n.sweeps)
      if (length(n.sweeps) != 1L || is.na(n.sweeps) || n.sweeps <= 0L) {
        stop("n.sweeps must be a single positive integer")
      }
      ptr <- getPointer()
      .Call(C_dbarts_bartcore_growFromRoot, ptr, n.sweeps)

      if (
        (is.na(updateState) && control@updateState == TRUE) ||
          identical(updateState, TRUE)
      ) {
        storeState(ptr)
      }

      invisible(NULL)
    },
    copy = function(shallow = FALSE) {
      "Creates a deep or shallow copy of the sampler."
      dupe <-
        if (shallow) {
          dbartsSampler$new(control, model, data)
        } else {
          newData <- data
          # only need to dupe things that can be changed internally, as the
          # rest will be simply swapped out
          newData@x <- .Call(C_dbarts_deepCopy, data@x)
          if (!is.null(data@x.test)) {
            newData@x.test <- .Call(C_dbarts_deepCopy, data@x.test)
          }
          dbartsSampler$new(control, model, newData)
        }

      # the stored state is opaque and never mutated in place (storeState
      # replaces it whole), so the copy can install the same object
      if (!is.null(state)) {
        dupe$setState(state)
      }
      dupe
    },
    show = function() {
      "Pretty prints the object."

      cat("dbarts sampler\n")
      cat("  call: ")
      writeLines(deparse(control@call))
      cat("\n")

      invisible(NULL)
    },
    predict = function(x.test, offset.test, n.threads = control@n.threads) {
      "Using existing sampler to predict for new data without re-running."
      ptr <- getPointer()

      x.test <- validateXTest(x.test, data@x)
      if (is.null(x.test)) {
        stop("x.test cannot be NULL")
      }
      # the engine's predict entry takes a dense matrix only; a frame/sparse
      # test set is coded above and materialized here, no resident sparse
      # predict store
      if (!is.matrix(x.test)) {
        x.test <- as.matrix(x.test)
      }

      if (missing(offset.test) || is.null(offset.test)) {
        offset.test <- NA_real_
      } else {
        offset.test <- as.double(offset.test)
        if (length(offset.test) == 1L) {
          offset.test <- rep_len(offset.test, nrow(x.test))
        }

        if (!identical(length(offset.test), nrow(x.test))) {
          stop(
            "length of test offset must be equal to number of rows in test matrix"
          )
        }
      }

      # the engine runs prediction serially; a missing offset is passed
      # as NULL rather than an NA sentinel
      if (length(offset.test) == 1L && is.na(offset.test)) {
        offset.test <- NULL
      }
      .Call(C_dbarts_bartcore_predict, ptr, x.test, offset.test)
    },
    setControl = function(newControl) {
      "Sets the control object for the sampler to a new one. Preserves the call() slot."
      if (!inherits(newControl, "dbartsControl")) {
        stop("'control' must inherit from dbartsControl")
      }

      selfEnv <- parent.env(environment())

      newControl@binary <- control@binary
      newControl@call <- control@call

      # settings fixed at creation: the generators and anything shaping
      # the cut grid
      for (slotName in c(
        "n.trees",
        "n.chains",
        "useQuantiles",
        "rngSeed"
      )) {
        if (!identical(slot(newControl, slotName), slot(control, slotName))) {
          stop(
            "the bartcore engine cannot change '",
            slotName,
            "' on an existing sampler"
          )
        }
      }
      if (newControl@keepTrees && is.na(newControl@n.samples)) {
        stop("keepTrees requires 'n.samples' to be specified")
      }

      ptr <- getPointer()
      selfEnv$control <- newControl
      .Call(C_dbarts_bartcore_setControl, ptr, control)

      invisible(NULL)
    },
    setModel = function(newModel) {
      "Sets the model object for the sampler to a new one."
      if (!inherits(newModel, "dbartsModel")) {
        stop("'model' must inherit from dbartsModel")
      }
      # the Dirichlet machinery is fixed at creation: a sampler cannot gain
      # or reconfigure it
      if (
        is(newModel@tree.prior, "dbartsDartPrior") ||
          is(model@tree.prior, "dbartsDartPrior")
      ) {
        stop("setModel cannot change a DART tree prior; recreate the sampler")
      }
      ptr <- getPointer()
      selfEnv <- parent.env(environment())

      newModel@family <- model@family
      oldModel <- model
      selfEnv$model <- newModel
      tryResult <- tryCatch(
        .Call(C_dbarts_bartcore_setModel, ptr, selfEnv$model, control, data),
        error = function(e) {
          selfEnv$model <- oldModel
          e$call <- quote(.Call(
            C_dbarts_bartcore_setModel,
            ptr,
            selfEnv$model,
            control,
            data
          ))
          e
        }
      )
      if (inherits(tryResult, "error")) {
        stop(tryResult)
      }

      invisible(NULL)
    },
    setData = function(newData, updateState = NA) {
      "Sets the data object for the sampler to a new one. Preserves the n.cuts and sigma slots. updateState is opt-in: only explicit TRUE stores state afterwards (NA/FALSE store nothing) - mutators are called per-sweep in Gibbs loops, so the default must stay free of that cost; contrast run()'s NA -> control@updateState convention."
      if (
        data@missing == "error" &&
          (anyNA(as.matrix(newData@x)) ||
            (!is.null(newData@x.test) && anyNA(newData@x.test)))
      ) {
        stop(
          "new predictors contain missing values and the sampler was built with missing = \"error\""
        )
      }
      bartcoreSamplerSetData(.self, newData)
      if (identical(updateState, TRUE)) {
        storeState()
      }
      invisible(NULL)
    },
    setResponse = function(y, updateScale = FALSE, updateState = NA) {
      "Changes the response against which the sampler is fitted. updateState is opt-in; see setData."
      bartcoreSamplerSetResponse(.self, y, updateScale)
      if (identical(updateState, TRUE)) {
        storeState()
      }
      invisible(NULL)
    },
    setOffset = function(offset, updateScale = FALSE, updateState = NA) {
      "Changes the offset slot used to adjust the response. updateState is opt-in; see setData."
      bartcoreSamplerSetOffset(.self, offset, updateScale)
      if (identical(updateState, TRUE)) {
        storeState()
      }
      invisible(NULL)
    },
    setWeights = function(weights, updateState = NA) {
      "Changes the weights with which the sampler is fitted. updateState is opt-in; see setData."
      weights <- as.double(weights)
      if (length(weights) != length(data@y)) {
        stop("'weights' must have length equal to that of 'y'")
      }
      if (anyNA(weights)) {
        stop("'weights' cannot be NA")
      }
      if (any(weights < 0.0)) {
        stop("'weights' must all be non-negative")
      }

      ptr <- getPointer()
      selfEnv <- parent.env(environment())

      oldWeights <- data@weights
      selfEnv$data@weights <- weights
      tryResult <- tryCatch(
        .Call(C_dbarts_bartcore_setWeights, ptr, data@weights),
        error = function(e) {
          selfEnv$data@weights <- oldWeights
          e
        }
      )
      if (inherits(tryResult, "error")) {
        stop(tryResult)
      }

      if (identical(updateState, TRUE)) {
        storeState(ptr)
      }
      invisible(NULL)
    },
    setSigma = function(sigma, updateState = NA) {
      "Changes the residual standard deviation parameter for each chain. updateState is opt-in; see setData."
      sigma <- as.double(sigma)
      if (length(sigma) != 1L) {
        stop("'sigma' must be of length 1")
      }
      if (is.na(sigma) || sigma <= 0.0) {
        stop("'sigma' must be positive")
      }

      ptr <- getPointer()
      .Call(C_dbarts_bartcore_setSigma, ptr, sigma)
      if (identical(updateState, TRUE)) {
        storeState(ptr)
      }
      invisible(NULL)
    },
    setPredictor = function(
      x,
      column,
      forceUpdate,
      updateCutPoints = FALSE,
      updateState = NA
    ) {
      "Changes a single column of the predictor matrix, or the entire matrix if column is missing. updateState is opt-in; see setData."

      checkMissingPolicy(data, anyNA(x), "predictors")
      result <- bartcoreSamplerSetPredictor(
        .self,
        x,
        column = if (missing(column)) NULL else column,
        forceUpdate = if (missing(forceUpdate)) NULL else forceUpdate,
        updateCutPoints = updateCutPoints
      )
      if (identical(updateState, TRUE)) {
        storeState()
      }
      # bartcoreSamplerSetPredictor returns invisible(NULL) or a visible
      # logical; preserve that (a bare 'result' would always be visible)
      if (is.null(result)) invisible(NULL) else result
    },
    setCutPoints = function(cuts, column, updateState = NA) {
      "Changes the cut points for the predictors in column, or the entire set itself if the column argument is missing. Forces the change by pruning any leaves that end up empty. updateState is opt-in; see setData."

      bartcoreSamplerSetCutPoints(
        .self,
        cuts,
        column = if (missing(column)) NULL else column
      )
      if (identical(updateState, TRUE)) {
        storeState()
      }
      invisible(NULL)
    },
    setTestPredictor = function(x.test, column) {
      "Changes a single column of the test predictor matrix."

      checkMissingPolicy(data, anyNA(x.test), "test predictors")
      bartcoreSamplerSetTestPredictor(
        .self,
        x.test,
        column = if (missing(column)) NULL else column
      )
    },
    setTestPredictorAndOffset = function(x.test, offset.test) {
      "Changes the test predictor matrix, and optionally the test offset."
      checkMissingPolicy(
        data,
        !is.null(x.test) && anyNA(x.test),
        "test predictors"
      )
      if (missing(offset.test)) {
        # predictors only; the engine keeps the current offset and the
        # bridge refuses if the row count would orphan its length
        return(bartcoreSamplerSetTestPredictor(.self, x.test, column = NULL))
      }

      x.test <- validateXTest(x.test, data@x)
      if (is.null(x.test) && !is.null(offset.test)) {
        stop("when test matrix is NULL, test offset must be as well")
      }
      if (!is.null(offset.test)) {
        offset.test <- as.double(offset.test)
        if (length(offset.test) == 1L) {
          offset.test <- rep_len(offset.test, nrow(x.test))
        }
        if (!identical(length(offset.test), nrow(x.test))) {
          stop(
            "length of test offset must be equal to number of rows in test matrix"
          )
        }
      }

      selfEnv <- parent.env(environment())
      oldTestUsesRegularOffset <- data@testUsesRegularOffset
      oldX.test <- data@x.test
      oldOffset.test <- data@offset.test

      selfEnv$data@testUsesRegularOffset <- FALSE
      selfEnv$data@x.test <- x.test
      selfEnv$data@offset.test <- offset.test
      tryResult <- tryCatch(
        .Call(
          C_dbarts_bartcore_setTestPredictorAndOffset,
          getPointer(),
          data@x.test,
          data@offset.test
        ),
        error = function(e) {
          selfEnv$data@testUsesRegularOffset <- oldTestUsesRegularOffset
          selfEnv$data@x.test <- oldX.test
          selfEnv$data@offset.test <- oldOffset.test
          e
        }
      )
      if (inherits(tryResult, "error")) {
        stop(tryResult)
      }
      invisible(NULL)
    },
    setTestOffset = function(offset.test) {
      "Changes the test offset."
      ptr <- getPointer()
      selfEnv <- parent.env(environment())

      selfEnv$data@testUsesRegularOffset <- FALSE
      if (!is.null(offset.test)) {
        if (is.null(data@x.test)) {
          stop("when test matrix is NULL, test offset must be as well")
        }
        if (length(offset.test) == 1L) {
          offset.test <- rep_len(offset.test, nrow(data@x.test))
        }
        if (length(offset.test) != nrow(data@x.test)) {
          stop(
            "length of test offset must be equal to number of rows in test matrix"
          )
        }
      }
      oldOffset.test <- data@offset.test
      selfEnv$data@offset.test <- offset.test
      tryResult <- tryCatch(
        .Call(C_dbarts_bartcore_setTestOffset, ptr, data@offset.test),
        error = function(e) {
          selfEnv$data@offset.test <- oldOffset.test
          e
        }
      )
      if (inherits(tryResult, "error")) {
        stop(tryResult)
      }

      invisible(NULL)
    },
    getLatents = function(result) {
      "For binary models, returns the current value of the latent variable representation."
      resultIsMissing <- missing(result)

      ptr <- getPointer()

      .Call(
        C_dbarts_bartcore_getLatents,
        ptr,
        if (resultIsMissing) NULL else result
      )
    },
    getSigmas = function(result) {
      "Return current residual error term on original, standard deviation scale."

      ptr <- getPointer()
      .Call(C_dbarts_bartcore_getSigmas, ptr)
    },
    getSumsOfSquaredResiduals = function(result) {
      "Return sum( (y - y.hat)^2 ) on original scale."
      ptr <- getPointer()
      .Call(C_dbarts_bartcore_getSumsOfSquaredResiduals, ptr)
    },
    getPointer = function() {
      "Returns the underlying reference pointer, checking for consistency first."
      selfEnv <- parent.env(environment())

      if (.Call(C_dbarts_bartcore_isValidPointer, pointer) == FALSE) {
        if (is.null(state)) {
          stop(
            "samplers cannot be re-created without a stored state; call ",
            "storeState() before serializing (see the Saving section of ",
            "?`dbartsSampler-class`)"
          )
        }
        selfEnv$pointer <- .Call(
          C_dbarts_bartcore_create,
          control,
          model,
          data,
          if (model@family == "auto") "" else model@family
        )
        # a same-spec continuation skips re-quantization; data@x serves any
        # cross-grid column (the engine keeps no predictor matrix)
        .Call(
          C_dbarts_bartcore_setState,
          pointer,
          state,
          rawPredictorMatrix(data@x)
        )
      }
      pointer
    },
    setState = function(newState) {
      "Sets the internal state from a cache."
      if (!inherits(newState, "bartcoreState")) {
        stop("'state' must inherit from bartcoreState")
      }
      selfEnv <- parent.env(environment())
      if (.Call(C_dbarts_bartcore_isValidPointer, pointer) == FALSE) {
        selfEnv$pointer <- .Call(
          C_dbarts_bartcore_create,
          control,
          model,
          data,
          if (model@family == "auto") "" else model@family
        )
      }
      .Call(
        C_dbarts_bartcore_setState,
        pointer,
        newState,
        rawPredictorMatrix(data@x)
      )
      selfEnv$state <- newState
      invisible(NULL)
    },
    storeState = function(ptr = getPointer()) {
      "Updates the cached internal state used for saving/loading."
      selfEnv <- parent.env(environment())
      selfEnv$state <- .Call(C_dbarts_bartcore_storeState, ptr)
      invisible(NULL)
    },
    installTrees = function(donor, samples = NULL) {
      "Warm-starts the forests from a donor sampler or bart fit over the same
       predictors. 'samples' maps each chain to a 1-based donor-sample index;
       NULL spreads the chains across the donor's kept samples."
      donorState <- warmStartState(donor)
      if (!is.null(samples)) {
        samples <- as.integer(samples)
      }
      ptr <- getPointer()
      .Call(C_dbarts_bartcore_installForests, ptr, donorState, samples)
      storeState(ptr)
      invisible(NULL)
    },
    printTrees = function(treeNums, chainNums, sampleNums) {
      "Produces an info dump of the internal state of the trees."
      matchedCall <- match.call()
      if (is.null(matchedCall$chainNums)) {
        chainNums <- seq_len(control@n.chains)
      }
      if (is.null(matchedCall$sampleNums)) {
        sampleNums <- if (control@keepTrees) {
          seq_len(control@n.samples)
        } else {
          NULL
        }
      } else {
        if (!control@keepTrees) {
          warning("sampleNums ignored if keepTrees is FALSE")
          sampleNums <- NULL
        } else {
          sampleNums <- as.integer(sampleNums)
        }
      }
      if (is.null(matchedCall$treeNums)) {
        treeNums <- seq_len(control@n.trees)
      }

      ptr <- getPointer()
      invisible(.Call(
        C_dbarts_bartcore_printTrees,
        ptr,
        as.integer(chainNums),
        sampleNums,
        as.integer(treeNums)
      ))
    },
    getTrees = function(
      treeNums,
      chainNums,
      sampleNums,
      current = FALSE,
      newdata = NULL
    ) {
      "Returns a data.frame containing the internal state of the trees."
      matchedCall <- match.call()
      current <- isTRUE(current)
      # live working trees have no sample dimension, so treat a current request
      # like a non-keepTrees sampler for sample handling
      useSaved <- control@keepTrees && !current
      if (is.null(matchedCall$chainNums)) {
        chainNums <- seq_len(control@n.chains)
      }
      if (is.null(matchedCall$sampleNums)) {
        sampleNums <- if (useSaved) {
          seq_len(control@n.samples)
        } else {
          NULL
        }
      } else {
        if (!useSaved) {
          warning(
            if (current) {
              "sampleNums ignored if current is TRUE"
            } else {
              "sampleNums ignored if keepTrees is FALSE"
            }
          )
          sampleNums <- NULL
        } else {
          sampleNums <- as.integer(sampleNums)
        }
      }
      if (is.null(matchedCall$treeNums)) {
        treeNums <- seq_len(control@n.trees)
      }

      chainNums <- as.integer(chainNums)
      treeNums <- as.integer(treeNums)

      if (any(chainNums <= 0 | chainNums > control@n.chains)) {
        stop("chainNums must be in [1, ", control@n.chains, "]")
      }
      if (
        useSaved &&
          any(sampleNums <= 0 | sampleNums > control@n.samples)
      ) {
        stop("sampleNums must be in [1, ", control@n.samples, "]")
      }
      if (any(treeNums <= 0 | treeNums > control@n.trees)) {
        stop("treeNums must be in [1, ", control@n.trees, "]")
      }

      # route new data through the trees so 'n' counts that data instead of the
      # training predictors; validated, coded, and materialized as for predict
      # (the tree-replay entry takes a dense matrix only)
      if (!is.null(newdata)) {
        newdata <- validateXTest(newdata, data@x)
        if (!is.matrix(newdata)) {
          newdata <- as.matrix(newdata)
        }
      }

      ptr <- getPointer()
      # saved-tree replay reads the current training predictors (the engine
      # keeps no matrix); a sparse data@x is skipped for a NULL replay source
      trees <- .Call(
        C_dbarts_bartcore_getTrees,
        ptr,
        chainNums,
        sampleNums,
        treeNums,
        current,
        newdata,
        rawPredictorMatrix(data@x),
        0L
      )
      # categorical rules report their split in 'directions' (value is NA);
      # when any column can hold one, pad the decode to the declared levels
      if (any(data@varTypes == CATEGORICAL_VARIABLE)) {
        trees <- decodeCategoricalSplits(trees, data@x, data@varTypes)
      }
      # rules on columns with missing values report their NA route
      if (!is.null(trees$missing)) {
        trees$missing <- c("L", "R")[trees$missing + 1L]
      }
      # linear leaves report one generically named slope column per
      # covariate; name them after the designated columns
      if (is(model@node.prior, "dbartsLinearPrior")) {
        covariateNames <- colnames(data@x)[model@node.prior@columns]
        if (!is.null(covariateNames)) {
          slopeColumns <- match(
            paste0("beta.", seq_along(covariateNames)),
            names(trees)
          )
          names(trees)[slopeColumns] <- paste0("beta.", covariateNames)
        }
      }
      trees
    },
    plotTree = function(
      treeNum,
      chainNum,
      sampleNum,
      treePlotPars = c(nodeHeight = 12, nodeWidth = 40, nodeGap = 8),
      ...
    ) {
      "Minimialist visualization of tree branching and contents."

      matchedCall <- match.call()
      if (is.null(matchedCall$chainNum)) {
        if (control@n.chains == 1L) {
          chainNum <- 1L
        } else {
          stop("chainNum required if more than one chain in sampler")
        }
      }
      if (is.null(matchedCall$sampleNum)) {
        sampleNum <- if (control@keepTrees) control@n.samples else 1L
      }

      tree <-
        if (control@keepTrees) {
          .self$getTrees(treeNum, chainNum, sampleNum)
        } else {
          .self$getTrees(treeNum, chainNum)
        }

      maxDepth <- getTreeDepthAndSize(tree)[["depth"]]

      tree <- cbind(
        tree,
        y = numeric(nrow(tree)),
        x = numeric(nrow(tree)),
        index = integer(nrow(tree))
      )
      tree <- fillPlotCoordinatesForNode(tree, maxDepth, 1L, 1L)
      numEndNodes <- tree$index[1L] - 1L
      #tree$index <- NULL

      plotHeight <- treePlotPars[["nodeHeight"]] *
        maxDepth +
        treePlotPars[["nodeGap"]] * (maxDepth - 1)
      dotsList <- list(...)
      dotsList$mar <- c(0, 0, 0, 0)
      par(dotsList)
      plot(
        NULL,
        type = "n",
        bty = "n",
        xaxt = "n",
        yaxt = "n",
        xlab = "",
        ylab = "",
        xlim = c(0, treePlotPars[["nodeWidth"]] * numEndNodes),
        ylim = c(0, plotHeight)
      )
      plotNode(tree, .self, treePlotPars)

      invisible(NULL)
    }
  )
)
