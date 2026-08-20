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
  seed = NA_integer_,
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
    seed = coerceOrError(seed, "integer"),
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
  sigma,
  sigest
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

  # Fitting functions (xbart) spell this 'sigest'; sampler constructors
  # (dbarts) keep 'sigma'.
  # Both route through this one validator under their own name, since the
  # two callers are matched here by argument name (redirectCall) and cannot
  # share a single formal without one silently dropping the other's value.
  if (!missing(sigest) && !is.na(sigest)) {
    tryCatch(sigest <- as.double(sigest), warning = function(e) {
      stop(
        "'sigest' argument to ",
        funcName,
        " must be coercible to numeric type"
      )
    })
    if (length(sigest) != 1L) {
      stop("'sigest' must be of length 1")
    }
    if (is.null(sigest) || sigest <= 0.0) {
      stop("'sigest' argument to ", funcName, " must be positive")
    }

    envir$sigest <- sigest
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
  interactions = NULL,
  blocks = NULL,
  variance = NULL,
  forests = NULL,
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

  # a forest() formula term declares an additional amplitude-coupled forest
  # (R/formulaTerms.R); checked against the requested family HERE, before any
  # family-specific remap or dispatch below can make the token unrecoverable
  # (hazard) or divert away from this function entirely (bart2's own
  # multinomial/hurdle.lognormal arcs, which carry the identical check)
  termIngestion <- ingestFormulaTerms(
    formula,
    family,
    if (missing(data)) NULL else data,
    missing(subset),
    matchedCall$subset,
    evalEnv
  )
  if (!is.null(termIngestion)) {
    if (!is.null(forests)) {
      stop(
        "'formula' declares a forest() term and 'forests' is also given; a ",
        "multi-forest model can only be declared one way - drop one"
      )
    }
    if (!missing(test)) {
      stop(
        "a forest() formula term does not support 'test'; drop the term or ",
        "fit a single-forest model"
      )
    }
    formula <- termIngestion$formula
    matchedCall$formula <- formula
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
  # a convenience mirror of dbartsControl(seed = ), as the wrappers expose;
  # an explicit seed overrides the control's, NA leaves it untouched
  seed <- coerceOrError(seed, "integer")
  if (!is.na(seed)) {
    control@seed <- seed
  }

  dataCall <- redirectCall(matchedCall, quoteInNamespace(dbartsData))
  # a forests = declaration's bases are conditioning DATA (docs/design/bcf.md),
  # so the columns they expand to ride the data object beside the weights they
  # mirror. Evaluated here, against this fit's own (pre-subset) data, since
  # dbartsData() has no index of its own to restrict a raw 'bases' entry by;
  # on the formula path, with 'subset' present, each one is then aligned to
  # the same rows the model frame will keep (resolveFormulaBasisSubset/
  # alignForestBasisToSubset, R/model.R).
  basisDeclarations <- forestBasisDeclarations(forests)
  # a basis declared on 'forests' has nowhere to ride once 'formula' is
  # already a built dbartsData: dbartsData() drops an unmatched 'bases'
  # argument in that case (its own ignored-args warning, R/data.R), so this
  # combination used to fit an ordinary single-forest model with the
  # declaration silently discarded. Refuse it by name instead. dbartsSpec() is
  # not touched - its first argument must already be a dbartsData, so this
  # predicate would be unconditionally true there and would refuse the
  # supported route (R/spec.R installs the declaration rather than dropping
  # it). The supported composition - a data object already carrying '@bases'
  # plus a knob-only 'forests' - is unaffected: forestBasisDeclarations()
  # returns NULL entries for a forest with no 'basis' of its own.
  if (
    !is.null(basisDeclarations) &&
      any(!vapply(basisDeclarations, is.null, logical(1L))) &&
      inherits(formula, "dbartsData")
  ) {
    stop(
      "'forests' declares a 'basis' but 'formula' is already a dbartsData; ",
      "a basis declaration cannot reach a pre-built data object through ",
      "dbarts() - use dbartsSpec(), or put the bases on the object with ",
      "dbartsData(bases = )"
    )
  }
  if (!is.null(basisDeclarations)) {
    # missing() reads this frame, so it is resolved here rather than inside the
    # per-forest closure below
    basisData <- if (missing(data)) NULL else data
    # NULL on the x/y interface or with no 'subset' - both leave every basis
    # untouched, below, exactly as before this rule existed
    subsetRows <- resolveFormulaBasisSubset(
      formula,
      basisData,
      matchedCall$subset
    )
    expanded <- lapply(seq_along(basisDeclarations), function(i) {
      basis <- expandForestBasis(
        evaluateForestBasis(basisDeclarations[[i]], basisData)
      )
      alignForestBasisToSubset(basis, i, subsetRows)
    })
    # a list in which no forest declares a basis is not a multi-forest
    # declaration at all - it names K ensembles with nothing to tell them
    # apart - so it falls through to resolveForests' own refusal by name
    if (any(!vapply(expanded, is.null, logical(1L)))) {
      dataCall$bases <- expanded
    } else {
      basisDeclarations <- NULL
    }
  }
  # a forest() term's bases were already evaluated against the model frame
  # (R/formulaTerms.R), post-subset, unlike the pre-subset raw-data channel
  # basisDeclarations above; 'forests' being NULL here is enforced by the
  # collision refusal above, so this never collides with it
  if (!is.null(termIngestion)) {
    # data@bases is read positionally against the forests it distinguishes
    # (forest 1 first); a term never speaks for forest 1, which has none
    dataCall$bases <- c(list(NULL), termIngestion$bases)
  }
  data <- if (is.null(basisDeclarations)) {
    eval(dataCall, evalEnv)
  } else {
    # the bases ride the data object's own 'bases' argument, which this caller
    # never wrote, so a refusal from it is restated in the word the caller
    # used; only validateForestBases names that argument. An unrelated error -
    # anything not naming 'bases' - is not this call's to relabel, so it keeps
    # its own condition class and call
    tryCatch(
      eval(dataCall, evalEnv),
      error = function(e) {
        message <- conditionMessage(e)
        if (!grepl("'bases'", message, fixed = TRUE)) {
          stop(e)
        }
        stop(gsub("'bases'", "'basis'", message, fixed = TRUE), call. = FALSE)
      }
    )
  }

  data@n.cuts <- rep_len(control@n.cuts, ncol(data@x))
  data@sigma <- sigma

  # a term's symbolic vars slot names design columns, which exist only now
  # (R/formulaTerms.R); forest 1 is the fit's own and reads every column, as
  # it does with no term at all
  if (!is.null(termIngestion)) {
    forests <- finalizeTermForests(termIngestion$pending, data)
  }

  spec <- resolveSamplerSpec(
    matchedCall,
    formals(dbarts),
    control,
    data,
    family,
    dispersion = dispersion,
    proposal.probs = proposal.probs,
    monotone = monotone,
    interactions = interactions,
    blocks = blocks,
    variance = variance,
    survivalStatus = survivalStatus,
    hazardPeriods = hazardPeriods,
    # the bases already ride the data object: dbartsData() is the one place
    # that knows which rows 'subset' kept
    bases = NULL,
    forests = forests,
    evalEnv = evalEnv
  )

  new("dbartsSampler", spec$control, spec$model, spec$data)
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
  # control@seed is NA (or pins them when it is set), giving independent
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
    # the noise a heteroscedastic prior predictive adds is s(x) eps, drawn
    # from the variance forest's own prior; the scalar draws below would
    # report a homoscedastic prior predictive instead. "ev" needs no noise
    # term and is unaffected.
    if (!is.null(attr(draw$control, "bartcore.variance"))) {
      stop(
        "samplePriorPredictive(type = \"ppd\") is not defined for a ",
        "heteroscedastic sampler: its residual scale is a variance forest, ",
        "not a single sigma"
      )
    }
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
    state = "ANY", # is either a list of states, or a promise to evaluate
    # Non-empty when this sampler is only the HOST SHELL of a fit whose model
    # lives elsewhere - bart2(family = "multinomial"), whose K-forest engine
    # rides the fit's $bc and reads nothing back from here. The string names
    # that fit and appears in the refusals below; empty for every sampler a
    # user builds, which is every sampler that owns its own model.
    hostFor = "character",
    # The per-forest, per-observation precision weight installed by
    # setForestWeights, mirrored here because it does not ride the engine's
    # saved state: forestWeights[[forest]] (1-based) holds the last vector
    # installed on that forest, NULL where none is. getPointer and setState
    # both re-apply it on every re-creation.
    forestWeights = "list"
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
      .self$hostFor <- character(0L)
      .self$forestWeights <- list()

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
    refuseHostMutation = function(what) {
      "Errors when this sampler is only the host shell of a fit whose model lives elsewhere; see the hostFor field."
      if (length(hostFor) == 0L) {
        return(invisible(NULL))
      }
      stop(
        what,
        " is not available on the host sampler of a ",
        hostFor,
        " fit; refit with bart2() to change the data"
      )
    },
    refuseHostRead = function(what) {
      "Errors when a read would answer from the host shell of a fit whose model lives elsewhere; the read-side sibling of refuseHostMutation. See the hostFor field."
      if (length(hostFor) == 0L) {
        return(invisible(NULL))
      }
      stop(
        what,
        " is not available on the host sampler of a ",
        hostFor,
        " fit; read the fit's own surface instead"
      )
    },
    run = function(
      numBurnIn,
      numSamples,
      updateState = NA,
      n.threads = control@n.threads
    ) {
      "Runs the posterior sampler and returns a list with the results."
      refuseHostMutation("$run")
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
      refuseHostMutation("$sampleTreesFromPrior")
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
      refuseHostMutation("$sampleNodeParametersFromPrior")
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
      refuseHostMutation("$growFromRoot")
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
        stop("'n.sweeps' must be a single positive integer")
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
      # forestWeights is a plain list field: assigning it shares the
      # underlying object, but R's copy-on-modify means a later
      # setForestWeights on either sampler duplicates before mutating, so
      # this is as alias-safe as the state install above. setState's own
      # reapply already ran against dupe's still-empty field, so reapply
      # again now that it carries .self's weights; dupe$pointer is valid
      # straight out of $new (and stays so after setState), which is what
      # lets this skip getPointer's re-creation branch
      dupe$forestWeights <- forestWeights
      dupe$reapplyForestWeights(dupe$pointer)
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
      # a sparse-backed test set rides to the engine as the container
      # validateXTest coded it; the engine routes its rows off that storage

      if (missing(offset.test) || is.null(offset.test)) {
        offset.test <- NA_real_
      } else {
        offset.test <- as.double(offset.test)
        if (length(offset.test) == 1L) {
          offset.test <- rep_len(offset.test, nrow(x.test))
        }

        if (!identical(length(offset.test), nrow(x.test))) {
          stop(
            "'offset.test' must have the same number of rows as 'x.test'"
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
      "Sets the control object for the sampler to a new one. Preserves the call() slot and any bartcore.* control attributes."
      refuseHostMutation("$setControl")
      if (!inherits(newControl, "dbartsControl")) {
        stop("'control' must inherit from dbartsControl")
      }

      selfEnv <- parent.env(environment())

      newControl@binary <- control@binary
      newControl@call <- control@call
      # bartcore.* attributes (the BCF, variance, grouped, survival and
      # ordinal/nbinom configuration resolveSamplerSpec attaches at creation,
      # docs/design/bcf.md) live outside the S4 slots newControl replaces
      # wholesale; a freshly built dbartsControl() never carries them, so
      # carry them forward exactly as binary/call are. Without this, a
      # legitimate setControl silently orphans the model configuration the
      # stored state's forest counts were built against, and the next
      # getPointer() re-creation refuses loudly instead (a BCF sampler: "the
      # data carry forest bases but no basis forest was configured"; a
      # heteroscedastic one: "state is not consistent with this sampler").
      for (attrName in grep(
        "^bartcore\\.",
        names(attributes(control)),
        value = TRUE
      )) {
        attr(newControl, attrName) <- attr(control, attrName)
      }

      # settings fixed at creation: the generators and anything shaping
      # the cut grid
      for (slotName in c(
        "n.trees",
        "n.chains",
        "useQuantiles",
        "seed"
      )) {
        if (!identical(slot(newControl, slotName), slot(control, slotName))) {
          stop(
            "changing '",
            slotName,
            "' is not available on an existing sampler"
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
      refuseHostMutation("$setModel")
      if (!inherits(newModel, "dbartsModel")) {
        stop("'model' must inherit from dbartsModel")
      }
      refuseBCFMutation(
        .self,
        "setModel",
        "every forest's node and tree priors are calibrated at creation; ",
        "make a new sampler instead"
      )
      # the Dirichlet machinery is fixed at creation: a sampler cannot gain
      # or reconfigure it
      if (
        is(newModel@tree.prior, "dbartsDartPrior") ||
          is(model@tree.prior, "dbartsDartPrior")
      ) {
        stop(
          "changing a DART tree prior is not available on an existing ",
          "sampler: recreate it instead"
        )
      }
      ptr <- getPointer()
      selfEnv <- parent.env(environment())

      newModel@family <- model@family
      oldModel <- model
      selfEnv$model <- newModel
      tryResult <- tryCatch(
        .Call(C_dbarts_bartcore_setModel, ptr, selfEnv$model, data),
        error = function(e) {
          selfEnv$model <- oldModel
          e$call <- quote(.Call(
            C_dbarts_bartcore_setModel,
            ptr,
            selfEnv$model,
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
      refuseHostMutation("$setData")
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
      refuseHostMutation("$setResponse")
      bartcoreSamplerSetResponse(.self, y, updateScale)
      if (identical(updateState, TRUE)) {
        storeState()
      }
      invisible(NULL)
    },
    setOffset = function(offset, updateScale = FALSE, updateState = NA) {
      "Changes the offset slot used to adjust the response. updateState is opt-in; see setData."
      refuseHostMutation("$setOffset")
      bartcoreSamplerSetOffset(.self, offset, updateScale)
      if (identical(updateState, TRUE)) {
        storeState()
      }
      invisible(NULL)
    },
    setWeights = function(weights, updateState = NA) {
      "Changes the weights with which the sampler is fitted. updateState is opt-in; see setData."
      refuseHostMutation("$setWeights")
      weights <- as.double(weights)
      if (length(weights) != length(data@y)) {
        stop("'weights' must have the same length as 'y'")
      }
      # NA/NaN and negative share one refusal, matching the bridge's
      # !(w >= 0.0) (is.na() is checked explicitly since NA_real_ >= 0.0
      # is itself NA in R, not FALSE)
      if (any(is.na(weights) | weights < 0.0)) {
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
    setActiveRows = function(active, updateState = NA) {
      "Sets the per-observation 0/1 mask of rows in the data set for this sampler. An inactive row leaves every sufficient statistic, every family-level parameter update and its own latent draw, but keeps its leaf occupancy and its fitted value. NULL clears, and an all-ones mask installs nothing. updateState is opt-in; see setData."
      refuseHostMutation("$setActiveRows")
      if (!is.null(active)) {
        active <- as.double(active)
        if (length(active) != length(data@y)) {
          stop("'active' must have the same length as 'y'")
        }
        if (anyNA(active)) {
          stop("'active' cannot be NA")
        }
        if (any(active != 0 & active != 1)) {
          stop("'active' must be all 0 or 1")
        }
      }

      ptr <- getPointer()
      .Call(C_dbarts_bartcore_setActiveRows, ptr, active)
      if (identical(updateState, TRUE)) {
        storeState(ptr)
      }
      invisible(NULL)
    },
    setForestWeights = function(forest, weights, updateState = NA) {
      "Sets a per-forest, per-observation weight: a multiplicative precision factor on the named forest's own leaf conditionals, composing with weights and active as (w_i * a_i) * m_f^2 * s_i rather than widening either channel. Only applies to a Bayesian causal forest built with forests = (see dbarts); forest indexes from 1, as with getCalibration/setCalibration (the basis forest is 2). The weight does not ride the sampler's saved state; it is mirrored on an R5 field that getPointer and setState both reinstall on every re-creation. updateState is opt-in; see setData."
      refuseHostMutation("$setForestWeights")
      weights <- as.double(weights)
      if (length(weights) != length(data@y)) {
        stop("'weights' must have the same length as 'y'")
      }
      # matches the bridge's !R_FINITE(w) || w < 0.0: a forest weight, unlike
      # a case weight, must also be finite, not merely non-negative
      if (!all(is.finite(weights)) || any(weights < 0.0)) {
        stop("'weights' must be finite and non-negative")
      }

      index <- resolveForestIndex(forest)
      ptr <- getPointer()
      selfEnv <- parent.env(environment())

      oldWeights <- forestWeights[index + 1L]
      selfEnv$forestWeights[index + 1L] <- list(weights)
      tryResult <- tryCatch(
        .Call(C_dbarts_bartcore_setForestWeights, ptr, index, weights),
        error = function(e) {
          selfEnv$forestWeights[index + 1L] <- oldWeights
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
    setForestBasis = function(forest, basis, updateState = NA) {
      "Changes the basis the named forest's amplitudes multiply, at any forest and any width. forest indexes from 1, as with setForestWeights and getCalibration/setCalibration (a Bayesian causal forest's basis forest is 2). A factor (or a one-sided formula naming one) expands to its level indicators, one amplitude per level, with no reference level dropped; a numeric vector or matrix is already those columns. This is the SOLE route by which a basis changes after creation, and the amplitudes are preserved and remapped: a width-preserving install leaves every one of them bitwise, and a width change carries each forest's block to its new offset and enters the added coordinates at 1. The matrix is mirrored into data@bases as setWeights mirrors weights, so it survives the sampler's re-creation. updateState is opt-in; see setData."
      refuseHostMutation("$setForestBasis")
      index <- resolveForestIndex(forest)
      if (is.null(basis)) {
        stop("'basis' cannot be NULL")
      }
      if (is.null(data@bases) || index >= length(data@bases)) {
        stop("forest index out of range")
      }
      values <- validateForestBases(
        list(expandForestBasis(evaluateForestBasis(basis))),
        length(data@y),
        argument = "basis"
      )[[1L]]

      ptr <- getPointer()
      selfEnv <- parent.env(environment())

      oldBases <- data@bases
      newBases <- oldBases
      newBases[[index + 1L]] <- values
      selfEnv$data@bases <- newBases
      tryResult <- tryCatch(
        .Call(C_dbarts_bartcore_setForestBasis, ptr, index, values),
        error = function(e) {
          selfEnv$data@bases <- oldBases
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
      refuseHostMutation("$setSigma")
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
      refuseHostMutation("$setPredictor")

      checkMissingPolicy(data, sourceAnyNA(x), "predictors")
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
      refuseHostMutation("$setCutPoints")

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
      refuseHostMutation("$setTestPredictor")

      checkMissingPolicy(data, sourceAnyNA(x.test), "test predictors")
      bartcoreSamplerSetTestPredictor(
        .self,
        x.test,
        column = if (missing(column)) NULL else column
      )
    },
    setTestPredictorAndOffset = function(x.test, offset.test) {
      "Changes the test predictor matrix, and optionally the test offset."
      refuseHostMutation("$setTestPredictorAndOffset")
      checkMissingPolicy(
        data,
        !is.null(x.test) && sourceAnyNA(x.test),
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
            "'offset.test' must have the same number of rows as 'x.test'"
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
      refuseHostMutation("$setTestOffset")
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
            "'offset.test' must have the same number of rows as 'x.test'"
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
      "Returns the current draw of the augmentation variable, whose meaning is per family and not uniform. A LOCATION, on the sampler's own latent scale, for probit (the truncated normal z), ordinal (the same z under the cut points) and aft (the imputed log survival time): these are regressed on directly. A PRECISION, one per observation, for logistic and nbinom (the Polya-Gamma omega) and Student-t (the scale-mixing lambda): these WEIGHT a working response and are not on the response scale at all. NULL for a plain gaussian sampler and for a multinomial one, neither of which augments. Note that a gaussian-family sampler built with resid.dist = student() DOES report latents, and they are precisions."
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
    getDispersion = function() {
      "Returns the negative-binomial dispersion r currently in force, one per chain, or NULL on every other family - the count analog of getSigmas(). It is the same scalar run()$dispersion records once per kept draw, read mid-sweep and without serializing state, so a host driving the sampler one sweep at a time reads it here instead of through storeState()$dispersion. Under a fixed dispersion it repeats the value the sampler was created with; otherwise it is that sweep's grid draw."
      refuseHostRead("$getDispersion")
      ptr <- getPointer()
      .Call(C_dbarts_bartcore_getDispersion, ptr)
    },
    getSumsOfSquaredResiduals = function(result) {
      "Return sum( (y - y.hat)^2 ) on original scale."
      ptr <- getPointer()
      .Call(C_dbarts_bartcore_getSumsOfSquaredResiduals, ptr)
    },
    getForestFits = function(forest) {
      "Returns a sampler's per-forest internal-scale fitted values (a Bayesian causal forest's 1 = prognostic, 2 = treatment; an ordinary sampler's only forest is 1), n.observations x n.chains. forest indexes from 1, as with setForestWeights/setForestBasis/getCalibration/setCalibration."
      ptr <- getPointer()
      .Call(C_dbarts_bartcore_getForestFits, ptr, resolveForestIndex(forest))
    },
    getFitsWithoutOffset = function() {
      "Returns the sampler's combined per-observation location on the RESPONSE scale and WITHOUT the installed offset, an n.observations x n.chains matrix; run()$train reports the same quantity with the offset folded in, so getFitsWithoutOffset() plus the installed offset is that value. This is the incremental read: getLatents() minus run()$train is biased, because the two are not on the same footing. Refused on a multinomial sampler, whose reported channels are per-category softmax probabilities rather than one additive location. Contrast getForestFits, which reports ONE forest's INTERNAL-scale totals."
      refuseHostRead("$getFitsWithoutOffset")
      ptr <- getPointer()
      .Call(C_dbarts_bartcore_getFitsWithoutOffset, ptr)
    },
    getForestAmplitudes = function(forest = NULL) {
      "Returns the named forest's amplitudes - the scalars its basis columns are multiplied by, one per column - as a q x n.chains matrix, or, at the default forest = NULL, every forest's stacked forest-major into a sum(q) x n.chains matrix, which is the row order the run's own glue channel carries. The vector is RAGGED, forest by forest, which is why a forest can be named: a Bayesian causal forest's forest 1 carries the single a on its implicit intercept and its forest 2 the pair (b0, b1) on its two level indicators, so the stacked read is its shipped (a, b0, b1). forest indexes from 1, as with setForestBasis/setForestWeights/getCalibration."
      ptr <- getPointer()
      .Call(
        C_dbarts_bartcore_getForestAmplitudes,
        ptr,
        if (is.null(forest)) NULL else resolveForestIndex(forest)
      )
    },
    getForestVariableCounts = function(forest) {
      "Returns a sampler's per-forest predictor split counts (a Bayesian causal forest's 1 = prognostic, 2 = treatment; an ordinary sampler's only forest is 1), n.predictors x n.chains, with rows named by the predictor columns when data@x carries colnames. forest indexes from 1, as with setForestWeights/setForestBasis/getCalibration/setCalibration."
      ptr <- getPointer()
      counts <- .Call(
        C_dbarts_bartcore_getForestVariableCounts,
        ptr,
        resolveForestIndex(forest)
      )
      predictorNames <- colnames(data@x)
      if (!is.null(predictorNames)) {
        rownames(counts) <- predictorNames
      }
      counts
    },
    getCalibration = function(forest = 1L) {
      "Returns the leaf-prior calibration in force, one row per chain and one column of prior.scale (the forest total's prior standard deviation at k = 1, in response units), prior.sd (prior.scale / k), prior.mean, k, k.has.hyperprior, response.scale, and response.shift, then the five multi-forest calibration-map quantities: amplitude.prior.variance and amplitude.prior.scale (exclusive - a forest carries a fixed amplitude variance or a half-Cauchy scale mixture, and the other reads NaN), node.scale.factor, node.scale.divisor, and basis.row.norm, which decompose prior.scale as factor * anchor / (divisor * row norm). All five are NaN on a forest whose scale the map does not own, and the two node.scale columns go NaN after a state install brings a foreign calibration, until setForestBasis re-imposes the map. The leaf model rides on a 'leaf.model' attribute and qualifies prior.sd: an equality only for the constant leaf. This is the authoritative reader of the calibration - model@prior.scale records the named intent, which a channel that re-anchors the response transform leaves untouched while moving what is in force."
      ptr <- getPointer()
      .Call(C_dbarts_bartcore_getCalibration, ptr, resolveForestIndex(forest))
    },
    setCalibration = function(
      prior.scale,
      prior.sd,
      prior.mean,
      forest = 1L,
      updateState = NA
    ) {
      "Restates a forest's leaf prior on every chain so that the forest total's prior standard deviation at k = 1 is prior.scale, in response units; prior.sd is the same statement at the current k and is refused when k is drawn from a hyperprior. Exactly one of the two is given. Nothing else moves - not k, not the response transform, not sigma, not the tree prior - and the write takes effect on the next sweep, reinterpreting no leaf value already drawn. updateState is opt-in; see setData."
      refuseHostMutation("$setCalibration")
      if (!missing(prior.mean)) {
        stop(
          "'prior.mean' is not writable: the leaf values it would shift are ",
          "already drawn and stored. The lever is the offset channel - ",
          "setOffset(rep_len(-getCalibration()[1L, \"prior.mean\"], n)) ",
          "re-centers the modelled quantity, and the getter then reports the ",
          "mean that is in force"
        )
      }
      if (missing(prior.scale) == missing(prior.sd)) {
        stop("give exactly one of 'prior.scale' and 'prior.sd'")
      }
      # after the argument checks above, so a malformed call is answered on its
      # own terms rather than by the refusal that would follow a well-formed one
      refuseBCFMutation(
        .self,
        "setCalibration",
        "every forest's leaf scale comes from the multi-forest calibration ",
        "map; make a new sampler instead"
      )

      ptr <- getPointer()
      index <- resolveForestIndex(forest)
      if (missing(prior.scale)) {
        prior.sd <- validateLiveScale(prior.sd, "prior.sd")
        calibration <- .Call(C_dbarts_bartcore_getCalibration, ptr, index)
        if (any(calibration[, "k.has.hyperprior"] != 0)) {
          stop(
            "'prior.sd' names a prior sd at the current 'k', but 'k' is ",
            "drawn every sweep under a hyperprior and the named sd would ",
            "drift: write 'prior.scale' instead (the prior sd at k = 1), or ",
            "pin 'k' at a number through setModel"
          )
        }
        # the sugar divides by a per-chain quantity, so it has no single
        # answer once the chains disagree about it; prior.scale is k-free and
        # remains available
        k <- unique(calibration[, "k"])
        if (length(k) != 1L) {
          stop(
            "the chains' 'k' have diverged (",
            paste0(format(k), collapse = ", "),
            "), so one 'prior.sd' names a different scale on each: write ",
            "'prior.scale', which does not divide by 'k'"
          )
        }
        prior.scale <- prior.sd * k
      } else {
        prior.scale <- validateLiveScale(prior.scale, "prior.scale")
      }

      .Call(C_dbarts_bartcore_setCalibration, ptr, index, prior.scale)
      if (identical(updateState, TRUE)) {
        storeState(ptr)
      }
      invisible(NULL)
    },
    reapplyForestWeights = function(ptr) {
      "Re-installs every forest weight mirrored on this sampler onto ptr, a freshly (re-)created or restated engine pointer that carries none of them. Called from getPointer and setState, never recursing through getPointer."
      for (forest in seq_along(forestWeights)) {
        weights <- forestWeights[[forest]]
        if (!is.null(weights)) {
          .Call(C_dbarts_bartcore_setForestWeights, ptr, forest - 1L, weights)
        }
      }
      invisible(NULL)
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
        reapplyForestWeights(pointer)
      }
      pointer
    },
    setState = function(newState) {
      "Sets the internal state from a cache."
      refuseHostMutation("$setState")
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
      reapplyForestWeights(pointer)
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
      refuseHostMutation("$installTrees")
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
        stop("'chainNums' must be in [1, ", control@n.chains, "]")
      }
      if (
        useSaved &&
          any(sampleNums <= 0 | sampleNums > control@n.samples)
      ) {
        stop("'sampleNums' must be in [1, ", control@n.samples, "]")
      }
      if (any(treeNums <= 0 | treeNums > control@n.trees)) {
        stop("'treeNums' must be in [1, ", control@n.trees, "]")
      }

      # route new data through the trees so 'n' counts that data instead of the
      # training predictors; validated and coded as for predict, and routed off
      # whatever storage it arrives in
      if (!is.null(newdata)) {
        newdata <- validateXTest(newdata, data@x)
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
