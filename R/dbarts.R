setMethod("initialize", "dbartsControl", function(.Object, ...) {
  .Object <- callNextMethod()

  validObject(.Object)
  .Object
})

# Parse a survival response: a survival::Surv object
# (recognized by inherits(), so survival need not be imported; right-censoring
# only in v1) or a plain two-column (time, status) matrix or data frame.
# Returns the raw event/censoring time and the 0/1 status vector, or NULL when
# the value is not a survival response. Errors on a non-right Surv (with a
# factor-status hint for "mright"), a non-two-column matrix, non-positive
# times, or a status outside {0, 1}; a Surv-like object with no type attribute
# is treated as right-censored. A missing time or status is a missing
# response, left NA for the caller's na.action, as survreg and coxph take it. Shared by the aft ingestion (which logs the
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
  observed <- !is.na(time)
  if (any(!is.finite(time[observed])) || any(time[observed] <= 0.0)) {
    stop("survival times must be finite and positive")
  }
  if (any(!is.na(status) & status != 0.0 & status != 1.0)) {
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
  # a row missing either part is a missing response
  log.time <- log(survival$time)
  log.time[is.na(survival$status)] <- NA_real_
  list(log.time = log.time, status = survival$status)
}

# Discrete-time hazard ingestion: the RAW time and status, the AFT sibling
# that skips the log() transform, for the person-period expander.
extractSurvivalTimes <- parseSurvivalResponse

# Resolve the discrete-time grid and each subject's terminal period from the
# observed times. `breaks` NULL (the
# default) uses the sorted distinct observed times (surv.bart's convention);
# a length-1 integer bins at the (1:K)/K quantiles (surv.bart's K); a longer
# numeric vector gives explicit interval boundaries b_0 < ... < b_K with
# right-closed intervals (b_{k-1}, b_k], the discSurv convention. Returns the
# representative period times (the right edges, sorted ascending) and each
# subject's terminal period index (1..K). Ties within a period are automatic:
# equal times share a period. findInterval(..., left.open = TRUE) counts grid
# points strictly below t, so a time exactly on grid point g_k lands in period
# k (its own interval's right edge). The grid comes from 'gridTime', the
# subjects the fit keeps, and a later time is placed in the last period.
resolveHazardGrid <- function(time, breaks, gridTime = time) {
  # a subject with no time has no period, and no place in the grid
  gridTime <- gridTime[!is.na(gridTime)]
  if (is.null(breaks)) {
    periods <- sort(unique(gridTime))
  } else {
    breaks <- as.double(breaks)
    if (anyNA(breaks)) {
      stop("'breaks' must not contain missing values")
    }
    if (length(breaks) == 1L) {
      K <- coerceOrError(breaks, "integer")
      if (is.na(K) || K < 1L) {
        stop(
          "'breaks' as a single value must be a positive integer period count"
        )
      }
      periods <- unique(as.double(
        quantile(gridTime, probs = seq_len(K) / K, names = FALSE)
      ))
    } else {
      if (is.unsorted(breaks, strictly = TRUE)) {
        stop("'breaks' boundaries must be strictly increasing")
      }
      if (
        any(time <= breaks[1L], na.rm = TRUE) ||
          any(time > breaks[length(breaks)], na.rm = TRUE)
      ) {
        stop(
          "every survival time must lie within the 'breaks' boundaries ",
          "(b_1, b_K]; widen the outer boundaries to cover the data"
        )
      }
      periods <- breaks[-1L]
    }
  }
  terminal <- pmin(
    findInterval(time, periods, left.open = TRUE) + 1L,
    length(periods)
  )
  # a subject with no time is at risk in the first period at least, as every
  # subject is; it keeps that one row, with a missing response
  terminal[is.na(terminal)] <- 1L
  list(periods = periods, terminalPeriod = terminal)
}

# The person-period expander: a subject
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
  weights = NULL,
  gridTime = time
) {
  n <- length(time)
  grid <- resolveHazardGrid(time, breaks, gridTime)
  periods <- grid$periods
  terminal <- grid$terminalPeriod

  Nprime <- sum(terminal)
  if (Nprime > max.rows) {
    stop(
      "person-period expansion would create ",
      Nprime,
      " rows, over the cap of ",
      max.rows,
      "; coarsen the time grid with family = hazard(breaks = ) (a boundary ",
      "vector or an integer period count), or raise the cap with ",
      "family = hazard(max.rows = )"
    )
  }

  # subject-major, period ascending within subject: the period-1 rows are the
  # subjects in order (survivalProbabilities reconstructs training covariates
  # from them), and sequence() supplies the within-subject period counter
  subjectOf <- rep.int(seq_len(n), terminal)
  periodOf <- sequence(terminal)
  y <- as.double(periodOf == terminal[subjectOf] & status[subjectOf] == 1.0)
  # a subject missing its time or its status is at risk with a missing
  # response, which the na.action then drops or refuses as it would any
  # other: to its time when it has one, and otherwise in the first period
  y[is.na(status[subjectOf]) | is.na(time[subjectOf])] <- NA_real_

  if ("period" %in% hazardPredictorNames(x)) {
    stop(
      "a hazard fit appends its own 'period' column; rename the predictor ",
      "'period'"
    )
  }
  xExpanded <- hazardRowSubset(x, subjectOf)
  xExpanded <- appendHazardPeriodColumn(xExpanded, periodOf)

  result <- list(x = xExpanded, y = y, periods = periods, subject = subjectOf)
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

# A formula-path hazard fit's na.action record, restated over person-period
# rows. The na.action ran on the subjects, before expansion, so the dropped
# subjects' rows are rebuilt from their own times on the kept subjects' grid,
# as the matrix interface, which expands first, would have dropped them; a
# subject with no time has its first-period row. Returns the record and the
# kept rows' make.unique names, taken over every subject so that they match
# that path.
hazardOmittedRows <- function(omitted, omittedTime, expansion, keptNames) {
  K <- length(expansion$periods)
  omittedSubjects <- unclass(omitted)
  # every kept subject has at least one row
  numSubjects <- max(expansion$subject) + length(omittedSubjects)
  keptSubjects <- seq_len(numSubjects)[-omittedSubjects]
  periodCounts <- integer(numSubjects)
  periodCounts[keptSubjects] <- tabulate(
    expansion$subject,
    length(keptSubjects)
  )
  omittedCounts <- pmin(
    findInterval(omittedTime, expansion$periods, left.open = TRUE) + 1L,
    K
  )
  # a subject with no time keeps the first period, as the expansion gives it
  omittedCounts[is.na(omittedCounts)] <- 1L
  periodCounts[omittedSubjects] <- omittedCounts
  subjectOf <- rep.int(seq_len(numSubjects), periodCounts)
  droppedRows <- which(subjectOf %in% omittedSubjects)
  names <- NULL
  if (!is.null(keptNames)) {
    subjectNames <- character(numSubjects)
    subjectNames[keptSubjects] <- keptNames
    subjectNames[omittedSubjects] <- names(omitted)
    names <- make.unique(subjectNames[subjectOf])
  }
  if (length(droppedRows) == 0L) {
    return(list(record = NULL, names = names))
  }
  record <- structure(
    droppedRows,
    names = names[droppedRows],
    class = class(omitted)
  )
  list(
    record = record,
    names = if (!is.null(names)) names[-droppedRows]
  )
}

# Names the person-period rows of both channels by R's make.unique over the
# subject names (dec-B34): a subject s1 at risk for three periods gives s1,
# s1.1, s1.2. Training rows are subject-major, test rows period-major (every
# test subject at each of the K periods).
hazardRowNames <- function(trainNames, subject, testNames, K) {
  list(
    train = if (!is.null(trainNames)) {
      make.unique(as.character(trainNames)[subject])
    },
    test = if (!is.null(testNames)) {
      make.unique(rep(as.character(testNames), times = K))
    }
  )
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
  kept <- attributes(x)[intersect(hazardDesignAttrs, names(attributes(x)))]
  if (inherits(x, "dbartsMixedMatrix")) {
    x$dense <- c(x$dense, list(as.double(period)))
    x$map <- c(x$map, length(x$dense))
    if (!is.null(x$columnNames)) {
      x$columnNames <- c(x$columnNames, "period")
    }
    out <- x
  } else {
    named <- !is.null(colnames(x))
    out <- cbind(x, period)
    if (named) {
      colnames(out)[ncol(out)] <- "period"
    } else {
      colnames(out) <- NULL
    }
  }
  # the builder attributes describe the design's columns, so the appended
  # ordinal column extends each of them
  if (!is.null(kept$term.labels)) {
    attr(out, "term.labels") <- c(kept$term.labels, "period")
  }
  if (!is.null(kept$drop)) {
    attr(out, "drop") <- c(kept$drop, list(period = FALSE))
  }
  if (!is.null(kept$varTypes)) {
    attr(out, "varTypes") <- c(kept$varTypes, ORDINAL_VARIABLE)
  }
  if (!is.null(kept$factor.levels)) {
    attr(out, "factor.levels") <- c(kept$factor.levels, list(NULL))
  }
  if (!is.null(kept$indicator.levels)) {
    attr(out, "indicator.levels") <- c(
      kept$indicator.levels,
      list(period = NULL)
    )
  }
  # the formula's terms cover the original columns; the period column is
  # read back by name
  if (!is.null(kept$terms)) {
    attr(out, "terms") <- kept$terms
  }
  out
}

hazardDesignAttrs <- c(
  "term.labels",
  "drop",
  "varTypes",
  "factor.levels",
  "indicator.levels",
  "terms"
)

# Row subset of a predictor set that keeps a container's columnar form and a
# matrix's builder attributes, which a bare subset drops.
hazardRowSubset <- function(x, rows) {
  if (is.data.frame(x) || inherits(x, "dbartsMixedMatrix")) {
    return(x[rows, , drop = FALSE])
  }
  x <- as.matrix(x)
  kept <- attributes(x)[intersect(hazardDesignAttrs, names(attributes(x)))]
  out <- x[rows, , drop = FALSE]
  for (a in names(kept)) {
    attr(out, a) <- kept[[a]]
  }
  out
}

hazardPredictorNames <- function(x) {
  if (is.data.frame(x)) names(x) else colnames(x)
}

## every slot below is passed explicitly to newValidated, so this
## function's own defaults are what apply, not A_class.R's prototype;
## only `binary` and `call`, which this constructor never sets, fall
## through to it. n.threads is deliberately one of the explicit ones:
## this default probes guessNumCores() capped to n.chains (dec-B115 -
## within-chain threading does not ship, so a budget above n.chains buys
## tree sampling nothing and is worth warning about, not defaulting to),
## while the prototype's is the conservative n.threads = 1L for a bare
## new("dbartsControl"). n.threads still keeps its own meaning, a total
## thread budget, distinct from n.chains; only the default is capped.
## dbartsControl(treeShift = ) in words, as the slot the bridge reads keeps
## it: "auto" is NA, the step taken where a forest's structure proposals are
## all zero, "always" TRUE and "never" FALSE.
resolveTreeShift <- function(treeShift) {
  choices <- c("auto", "always", "never")
  if (identical(treeShift, choices)) {
    return(NA)
  }
  index <- if (is.character(treeShift) && length(treeShift) == 1L) {
    pmatch(treeShift, choices)
  } else {
    NA_integer_
  }
  if (is.na(index)) {
    stop(
      "'treeShift' must be one of \"auto\", \"always\" or \"never\"",
      call. = FALSE
    )
  }
  c(NA, TRUE, FALSE)[index]
}

## The argument dbartsControl() would be called with to reproduce a control's
## slot: treeShift is spelled in words while its slot keeps the bridge's
## tri-state logical, and the slots hold NA where n.samples and seed take NULL.
controlArgumentFromSlot <- function(name, control) {
  if (identical(name, "treeShift")) {
    level <- control@levelGibbs
    return(
      if (is.na(level)) {
        "auto"
      } else if (level) {
        "always"
      } else {
        "never"
      }
    )
  }
  value <- methods::slot(control, name)
  if (name %in% c("n.samples", "seed") && is.na(value)) NULL else value
}

dbartsControl <- function(
  verbose = FALSE,
  keepTrainingFits = TRUE,
  keepFits = TRUE,
  useQuantiles = FALSE,
  treeShift = c("auto", "always", "never"),
  keepTrees = FALSE,
  storage = c("double", "single"),
  n.samples = NULL,
  n.cuts = 100L,
  n.burn = 200L,
  n.trees = 75L,
  n.chains = 4L,
  n.threads = min(dbarts::guessNumCores(), n.chains),
  n.thin = 1L,
  printEvery = 100L,
  printCutoffs = 0L,
  categoricalExhaustiveCap = 10L,
  testFitParallelCutoff = 65536L,
  predictParallelCutoff = 50000L,
  sparseDensityThreshold = 0.2,
  proposal.probs = c(
    birth_death = 0.6,
    swap = 0,
    change = 0.4,
    perturb = 0,
    rule_gibbs = 0,
    birth = 0.5
  ),
  seed = NULL,
  updateState = TRUE,
  ...
) {
  # the names THIS call carried, which the front doors' precedence rule reads
  # off the object: a slot set to the same value the constructor would have
  # chosen is otherwise indistinguishable from an untouched one
  namedHere <- names(match.call())[-1L]
  # '...' exists only so a retired argument name reaches a message naming
  # its successor; R refuses an unknown name before any body runs
  supplied <- dotNames(...)
  refuseForeignFrontDoorArgs(
    supplied,
    "dbartsControl",
    names(formals(dbarts::dbartsControl))
  )
  # every name this door carries on '...' is an ordinary value, so forcing
  # them here is safe (bart's are not: one of them is written in a
  # vocabulary that only resolves inside this package)
  seed <- resolveRenamedSeed(
    if ("rngSeed" %in% supplied) list(...)[["rngSeed"]] else NULL,
    "dbartsControl",
    seed
  )

  storage <- match.arg(storage)
  # the slot keeps the bridge's tri-state logical, whose NA is the automatic
  # mode; only this argument speaks in words
  levelGibbs <- resolveTreeShift(treeShift)
  # the slot keeps NA for "not set"; NULL is the argument's spelling of it
  if (is.null(n.samples)) {
    n.samples <- NA_integer_
  } else {
    refuseNaN(n.samples, "n.samples")
    if (isSingleNA(n.samples)) {
      warnNAForNull("n.samples", "dbartsControl")
    }
  }
  result <- newValidated(
    "dbartsControl",
    verbose = as.logical(verbose),
    keepTrainingFits = as.logical(keepTrainingFits),
    keepFits = as.logical(keepFits),
    useQuantiles = as.logical(useQuantiles),
    levelGibbs = levelGibbs,
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
    categoricalExhaustiveCap = coerceOrError(
      categoricalExhaustiveCap,
      "integer"
    ),
    testFitParallelCutoff = coerceOrError(testFitParallelCutoff, "integer"),
    predictParallelCutoff = coerceOrError(predictParallelCutoff, "integer"),
    sparseDensityThreshold = coerceOrError(sparseDensityThreshold, "numeric"),
    # the partial spellings are filled here rather than at the slot, so the
    # stored mixture is always the resolved six the bridge reads
    proposal.probs = resolveProposalProbs(proposal.probs),
    seed = resolveSeedArg(seed, "dbartsControl", refuse = TRUE),
    updateState = as.logical(updateState)
  )
  # a plain attribute, deliberately not a bartcore.* one (that prefix means
  # fit state, which setControl carries forward and the front doors refuse)
  # and deliberately not a slot (a slot would enter every validity check and
  # every stored control's printed form)
  attr(result, controlSuppliedAttr) <- unique(c(
    intersect(namedHere, names(formals(dbarts::dbartsControl))),
    # the retired spelling names the same slot
    if ("rngSeed" %in% supplied) "seed"
  ))
  result
}

## The attribute dbartsControl() stamps with the names its call carried.
controlSuppliedAttr <- "dbarts.supplied"

## The slots a control was CONSTRUCTED with, empty for one that carries no
## record (a bare new("dbartsControl"), or an object saved before the record
## existed), which then falls back to the differs-from-the-default test.
controlSuppliedSlots <- function(control) {
  named <- attr(control, controlSuppliedAttr)
  if (is.null(named)) character(0L) else named
}

## The one precedence rule bart() and xbart() apply to a supplied control: a
## flat formal the caller named wins over the control's slot of the same name,
## and a slot they did not name flat stands. 'flat' names each shared setting
## and holds its flat value; the resolved list comes back under the same names,
## for the door to read and to write onto the control where it keeps one.
## A slot the control never spoke for names nothing: each door carries its own
## defaults for the settings it also spells flat (bart's 500 samples and
## verbose = TRUE against the control's unset one and FALSE), so a control
## built to reach one engine setting must not silently move every other one
## with it. A slot SPEAKS when the control's own constructor call named it -
## which is what the dbarts.supplied record is for, since a value equal to the
## constructor's default is otherwise indistinguishable from an untouched slot
## - or when it differs from a fresh control's, which is how a post-
## construction edit (ctl@n.trees <- 200L) still speaks.
mergeFrontDoorControl <- function(control, matchedCall, flat) {
  supplied <- names(matchedCall)
  # no control named, nothing to merge: the door's own defaults stand, and the
  # reference control below (which probes the core count) is never built
  if ("control" %not_in% supplied) {
    return(flat)
  }
  fresh <- dbarts::dbartsControl()
  namedOnControl <- controlSuppliedSlots(control)
  for (name in names(flat)) {
    if (name %in% supplied) {
      next
    }
    if (
      name %in%
        namedOnControl ||
        !identical(methods::slot(control, name), methods::slot(fresh, name))
    ) {
      flat[name] <- list(controlArgumentFromSlot(name, control))
    }
  }
  flat
}

## The retired flat 'proposal.probs' beside a control whose own call NAMED the
## same slot: one setting written twice, refused by name exactly as a prior
## object supplied beside its shorthand is. A control whose slot merely DIFFERS
## from a fresh one - a post-construction edit - is not a collision, and the
## retired flat still wins there.
refuseCollidingMixture <- function(control) {
  if ("proposal.probs" %in% controlSuppliedSlots(control)) {
    stop(
      "'control' cannot be combined with 'proposal.probs': set the tree-move ",
      "mixture in one place, dbartsControl(proposal.probs = )"
    )
  }
  invisible(NULL)
}

## A control taken from a fitted sampler carries that fit's model configuration
## on bartcore.* attributes - the variance forest, the shape, the survival
## status, the forest map - which a new fit over new data has no claim to.
## bart and xbart, which build their own control, refuse it by name; dbarts and
## dbartsSpec, where passing a sampler's control on is a 0.9-x pattern, strip
## the attributes instead, and neither carries them into the new fit.
refuseFitStateControl <- function(control, caller) {
  if (!inherits(control, "dbartsControl")) {
    stop(
      "'control' argument to ",
      caller,
      " must be of class dbartsControl; use dbartsControl() to create one",
      call. = FALSE
    )
  }
  if (any(startsWith(names(attributes(control)), "bartcore."))) {
    stop(
      "'control' carries the model configuration of the fit it was taken ",
      "from; '",
      caller,
      "' builds its own - pass a fresh dbartsControl() instead",
      call. = FALSE
    )
  }
  control
}

validateArgumentsInEnvironment <- function(
  envir,
  func,
  funcName,
  control,
  verbose,
  n.samples,
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
    n.samples <- coerceOrError(n.samples, "integer")
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

  # One name for the residual-standard-deviation estimate supplied at
  # creation, on every entry point; the sampler's setSigma, which sets the
  # parameter rather than an estimate of it, is a different thing and keeps
  # its own name.
  # NULL and NA are "not given" and were resolved, with any warning, by the
  # entry point; the estimate the entry point holds is then NA_real_
  if (!missing(sigest) && !is.null(sigest) && !isSingleNA(sigest)) {
    envir$sigest <- validateSigest(sigest, funcName)
  }
}

validateSigest <- function(sigest, funcName) {
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
  if (is.na(sigest) || sigest <= 0.0) {
    stop("'sigest' argument to ", funcName, " must be positive")
  }
  sigest
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
  leaf.prior = normal,
  monotone = NULL,
  interactions = NULL,
  blocks = NULL,
  variance = NULL,
  forests = NULL,
  control = dbarts::dbartsControl(),
  sigest = NULL,
  seed = NULL,
  factors = c("categorical", "indicators"),
  family = c(
    "auto",
    "gaussian",
    "student",
    "probit",
    "logistic",
    "aft",
    "multinomial",
    "ordinal",
    "nbinom",
    "hazard",
    "hazard.probit",
    "hazard.logistic"
  ),
  na.action = dbarts::na.keepPredictors,
  sigma = NULL,
  node.prior = NULL,
  callback = NULL,
  ...
) {
  matchedCall <- match.call()

  evalEnv <- parent.frame(1L)

  # '...' carries the names dec-B98's consolidation moved onto the family and
  # prior objects, for one release; anything else is a caller mistake and is
  # refused by name rather than dropped without a word
  supplied <- dotNames(...)
  refuseForeignFrontDoorArgs(supplied, "dbarts", names(formals(dbarts::dbarts)))
  consolidated <- resolveConsolidatedArgs(
    matchedCall,
    supplied,
    "dbarts",
    evalEnv
  )
  # cleared from the matched call before anything is forwarded: the prior
  # resolver still carries a 'resid.prior' formal for the object to reach,
  # and a name left standing here would reach it without the reconciliation
  # below
  if (length(consolidated) > 0L) {
    matchedCall[names(consolidated)] <- NULL
  }
  # the tree-move mixture is a control setting now; the retired spelling is
  # honored where the control's own slot would otherwise stand, and refused
  # where the control's own call named that slot too
  proposal.probs <- consolidated[["proposal.probs"]]
  if (
    "proposal.probs" %in%
      names(consolidated) &&
      "control" %in% names(matchedCall)
  ) {
    refuseCollidingMixture(control)
  }

  # dbarts() never runs the sampler itself, so 'callback' has nothing to
  # drive here; it is validated anyway, ahead of the (possibly expensive)
  # sampler construction below, so a malformed pair fails at THIS call
  # rather than silently doing nothing until a later $run()
  validateCallback(callback)

  # the creation-time estimate is 'sigest' here as everywhere; the 0.9-x
  # spelling is folded in before the shared validator, which knows one name
  sigmaSupplied <- !missing(sigma) && !is.null(sigma)
  sigest <- resolveRenamedSigma(
    !sigmaSupplied,
    missing(sigest),
    sigma,
    sigest,
    "dbarts"
  )
  # an NA under the retired name is covered by that name's own warning
  sigest <- if (sigmaSupplied) {
    resolveSigestArg(sigest, "dbarts", "silent", "sigma")
  } else {
    resolveSigestArg(sigest, "dbarts", "refuse")
  }
  if (sigmaSupplied) {
    # the value, not the promise: it was evaluated once above, and NULL or NA
    # must reach the validator as the resolved estimate
    matchedCall["sigest"] <- list(if (is.na(sigest)) NULL else sigest)
    matchedCall$sigma <- NULL
  }

  # the leaf-value prior is 'leaf.prior' here as everywhere; 'node.prior'
  # is the 0.9-x spelling, accepted for one release. Both flags are read
  # before either name is assigned: an assignment makes missing() false.
  nodePriorSupplied <- !missing(node.prior)
  leafPriorSupplied <- !missing(leaf.prior)
  matchedCall <- resolveRenamedLeafPrior(
    matchedCall,
    nodePriorSupplied,
    leafPriorSupplied,
    "dbarts"
  )

  # 'family' is resolved from the caller's own unevaluated argument, so a
  # bare family constructor (student(3)) resolves in the family vocabulary;
  # nothing may force the argument before this. hurdle.lognormal is admitted
  # only to be refused by name below, with the reason.
  familySpec <- resolveFamily(
    matchedCall$family,
    eval(formals(dbarts::dbarts)$family),
    "dbarts",
    evalEnv,
    refused = "hurdle.lognormal"
  )
  family <- familySpec@token

  # a hurdle response is a composition of two independent samplers, which
  # this function cannot return; refused BY NAME, whose generic message would
  # name neither the front door nor why.
  if (identical(family, "hurdle.lognormal")) {
    stop(
      "dbarts() does not fit family = \"hurdle.lognormal\": it composes ",
      "two independent samplers (a zero-part probit and a positive-part ",
      "gaussian) and dbarts() returns one - use ",
      "bart(x, y, family = \"hurdle.lognormal\")"
    )
  }
  # the caller's own token, ahead of the student and hazard remaps below -
  # every downstream refusal that names 'family' echoes this, not the
  # resolved spelling, which is an implementation detail
  requestedFamily <- family

  # the family-only settings, read off the object rather than off formals
  # this signature no longer carries
  shape <- familySetting(familySpec, "shape", NA_real_)
  breaks <- familySetting(familySpec, "breaks", NULL)
  max.rows <- familySetting(familySpec, "max.rows", 1e7)
  # The residual scale's prior has one home, the family object it rides; the
  # retired flat spelling is still read for one release, and a flat spelling
  # beside a family that named 'sigma' too is refused unless the two agree.
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
  # Student-t is its own family token and its own engine family; on this
  # side of the bridge it is a gaussian response carrying a degrees-of-
  # freedom attribute, so the remap happens here, once
  residDf <- NULL
  if (identical(family, "student")) {
    residDf <- familySetting(familySpec, "df", NA_real_)
    family <- "gaussian"
  }

  # the forest constructors resolve by bare name inside their arguments, from
  # the caller's own unevaluated arguments; nothing may force them before this
  forestArguments <- resolveForestArguments(matchedCall, evalEnv)
  forests <- forestArguments$forests
  interactions <- forestArguments$interactions
  blocks <- forestArguments$blocks
  monotone <- forestArguments$monotone
  variance <- forestArguments$variance

  # a forest() formula term declares an additional amplitude-coupled forest
  # (R/formulaTerms.R); checked against the requested family HERE, before any
  # family-specific remap or dispatch below can make the token unrecoverable
  # (hazard) or divert away from this function entirely (bart2's own
  # multinomial/hurdle.lognormal arcs, which carry the identical check)
  termIngestion <- ingestFormulaTerms(
    formula,
    family,
    if (missing(data)) NULL else data
  )
  if (!is.null(termIngestion)) {
    if (!is.null(forests)) {
      stop(
        "'formula' declares a forest() term and 'forests' is also given; a ",
        "multi-forest model can only be declared one way - drop one"
      )
    }
    # a formula whose one forest() has no basis is a single-forest fit, which
    # takes 'test' as the same terms written plainly do
    if (
      !missing(test) &&
        !all(vapply(termIngestion$basisReads, is.null, logical(1L)))
    ) {
      stop(
        "a forest() formula term does not support 'test'; drop the term or ",
        "fit a single-forest model"
      )
    }
    formula <- termIngestion$formula
    matchedCall$formula <- formula
  }

  # survival response ingestion: a survival::Surv
  # object or an explicit family = "aft" with a two-column (time, status)
  # response fits the AFT log-normal model. The matrix (x.train, y.train)
  # interface is handled directly here, the response being the second
  # positional argument; the log event/censoring time replaces it and the
  # status rides the control attribute the bartcore survival family reads.
  # 'subset' is honoured (applied to the status vector alongside the
  # matchedCall$subset dbartsData() still applies to x/y itself, below). A
  # Surv left-hand side on 'formula' is ingested by dbartsData()'s own
  # formula branch (R/data.R) - it cannot be detected here, before the model
  # frame exists - so the matching conflict guard and auto-dispatch for that
  # route run again, once dbartsData() returns, further down.
  survivalStatus <- NULL
  directResponse <- !is.formula(formula) &&
    !inherits(formula, "dbartsData") &&
    !inherits(formula, "dgCMatrix") &&
    !missing(data)
  responseIsSurv <- directResponse && inherits(data, "Surv")
  # a caller-built dbartsData object (dbartsData(Surv(...) ~ ., data), called
  # directly rather than through dbarts()) can ALREADY carry the same
  # attributes the formula route stashes below - dbartsData() returns such
  # an object unchanged (R/data.R), so an explicit family = "aft"/hazard
  # token on it is a legitimate request, not the unsupported-interface case
  # the guards further down otherwise refuse
  survivalDataObject <- inherits(formula, "dbartsData") &&
    !is.null(attr(formula, "survivalStatus"))
  hazardTokens <- hazardFamilyTokens
  # a Surv response declares the model, so it auto-dispatches to aft from
  # "auto"; an explicit hazard token selects the discrete-time model instead
  # (the guard whitelist admits it). Any
  # other explicit family with a Surv response is a conflict, never a silent
  # override.
  if (responseIsSurv && family %not_in% c("auto", "aft", hazardTokens)) {
    stop(
      "a survival (Surv) response cannot be fit with family \"",
      family,
      "\"; use family \"aft\", \"hazard\", or \"auto\""
    )
  }

  # discrete-time hazard ingestion: person-period-expand (x, time, status)
  # into an ordinary binary
  # (X', y') design and REMAP the hazard token to its underlying binary link
  # BEFORE any family-keyed switch runs (leaf.scale, control@binary,
  # fixedUnitScale, the weight policy). The engine, bridge, and ResponseModels
  # then see an ordinary probit/logistic fit; the hazard provenance survives
  # only as the period grid, parked on the control attribute the packaging
  # reads into the $periods marker (the bartcore.survival -> $status
  # precedent). No status vector or attribute reaches C++ - the censoring is
  # baked into y'.
  hazardPeriods <- NULL
  # the person-period row names, set on the data object once it exists, and
  # whether the expansion ran before the na.action did
  hazardNames <- NULL
  hazardExpandedFirst <- FALSE
  hazardOffsetArgument <- NULL
  if (family %in% hazardTokens && directResponse) {
    survival <- extractSurvivalTimes(data)
    if (is.null(survival)) {
      stop(
        "family \"",
        family,
        "\" needs a survival::Surv or two-column (time, status) response"
      )
    }
    xForExpansion <- formula
    timeForExpansion <- survival$time
    statusForExpansion <- survival$status
    offsetForExpansion <- if (missing(offset)) NULL else offset
    weightsForExpansion <- if (missing(weights)) NULL else weights
    # the original row indices no longer mean anything once the design is
    # expanded to N' person-period rows (dec-B97), so 'subset' is applied
    # HERE, before expansion, rather than forwarded to dbartsData()
    if (!missing(subset)) {
      xForExpansion <- xForExpansion[subset, , drop = FALSE]
      timeForExpansion <- timeForExpansion[subset]
      statusForExpansion <- statusForExpansion[subset]
      if (!is.null(offsetForExpansion) && length(offsetForExpansion) > 1L) {
        offsetForExpansion <- offsetForExpansion[subset]
      }
      if (!is.null(weightsForExpansion) && length(weightsForExpansion) > 1L) {
        weightsForExpansion <- weightsForExpansion[subset]
      }
    }
    # the na.action drops a subject's rows together, and the grid must not
    # depend on them, as on the formula path; a missing time or status is a
    # missing response
    keptSubjects <- applyNaActionToXY(
      na.action,
      ifelse(is.na(statusForExpansion), NA_real_, timeForExpansion),
      xForExpansion
    )
    expansion <- expandDiscreteTimeHazard(
      xForExpansion,
      timeForExpansion,
      statusForExpansion,
      breaks = breaks,
      max.rows = max.rows,
      offset = offsetForExpansion,
      weights = weightsForExpansion,
      gridTime = if (!is.null(keptSubjects)) {
        timeForExpansion[keptSubjects$keep]
      } else {
        timeForExpansion
      }
    )
    matchedCall$formula <- expansion$x
    matchedCall$data <- expansion$y
    if (!missing(subset)) {
      matchedCall$subset <- NULL
    }
    if (!is.null(expansion$offset)) {
      # the per-subject 'offset' as written, which the expanded data object's
      # record of it replaces once that object exists (below)
      hazardOffsetArgument <- offsetArgumentFormula(
        matchedCall$offset,
        evalEnv,
        as.character(colnames(formula)),
        NROW(formula)
      )
      matchedCall$offset <- expansion$offset
    }
    if (!is.null(expansion$weights)) {
      matchedCall$weights <- expansion$weights
    }
    K <- length(expansion$periods)
    hazardExpandedFirst <- TRUE
    hazardNames <- hazardRowNames(
      observationRowNames(xForExpansion),
      expansion$subject,
      if (!missing(test)) observationRowNames(test),
      K
    )
    # a held-out subject has no event time to place it by, so 'test' expands
    # to every one of the SAME K training periods (the shape
    # hazardSurvivalProbabilities's own newdata expansion builds, reused here
    # via appendHazardPeriodColumn); survivalProbabilities then reads the
    # stored per-period draws straight off the fit (R/bart.R)
    if (!missing(test)) {
      n.test <- NROW(test)
      matchedCall$test <- if (is.data.frame(test)) {
        bigTest <- test[rep(seq_len(n.test), times = K), , drop = FALSE]
        bigTest[["period"]] <- rep(seq_len(K), each = n.test)
        bigTest
      } else {
        appendHazardPeriodColumn(
          as.matrix(test)[rep(seq_len(n.test), times = K), , drop = FALSE],
          rep(seq_len(K), each = n.test)
        )
      }
      if (!missing(offset.test)) {
        offsetTestForExpansion <- offset.test
        if (length(offsetTestForExpansion) == 1L) {
          offsetTestForExpansion <- rep_len(offsetTestForExpansion, n.test)
        }
        matchedCall$offset.test <- rep(offsetTestForExpansion, times = K)
      }
    } else {
      matchedCall$test <- NULL
      matchedCall$offset.test <- NULL
    }
    hazardPeriods <- expansion$periods
    # the remap: the engine-facing family is now an ordinary binary link
    family <- if (identical(family, "hazard.logistic")) "logistic" else "probit"
    # the survival response is consumed; do not let the aft block fire on it
    responseIsSurv <- FALSE
  } else if (
    family %in%
      hazardTokens &&
      !directResponse &&
      !is.formula(formula) &&
      !survivalDataObject
  ) {
    stop(
      "discrete-time hazard fits currently use the matrix interface - ",
      "dbarts(x, y) or bart(x, y) with a survival::Surv or two-column ",
      "(time, status) response, or a formula with a Surv left-hand side"
    )
  }
  # else: a Surv-formula hazard request (family %in% hazardTokens,
  # is.formula(formula)) OR one on a pre-built dbartsData object carrying
  # the same attributes (survivalDataObject) is expanded further down, once
  # dbartsData() has ingested and subsetted the response against the model
  # frame's own rows - it cannot be detected here, before that frame exists

  # aft is reachable through the direct-response form, through a Surv-formula
  # response or a pre-built dbartsData object carrying the same attributes
  # (dbartsData()'s own short-circuit, R/data.R; the matching conflict guard
  # and auto-dispatch run again below, once 'data' is built); every other
  # indirect route is refused up
  # front, before the response is materialized, rather than failing
  # hostilely downstream
  if (
    family == "aft" &&
      !directResponse &&
      !is.formula(formula) &&
      !survivalDataObject
  ) {
    stop(
      "survival (aft) fits currently use the matrix interface - ",
      "dbarts(x, y) or bart(x, y) with a survival::Surv or two-column ",
      "(time, status) response, or a formula with a Surv left-hand side"
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
    family <- "aft"
    matchedCall$data <- survival$log.time
    survivalStatus <- survival$status
    # 'subset' still reaches dbartsData() unchanged (matchedCall$subset) and
    # subsets x/y itself, reading the same value - aft's row count is
    # unchanged, unlike hazard's, so only the status vector needs its own
    # subsetting here
    if (!missing(subset)) {
      survivalStatus <- survivalStatus[subset]
    }
  }

  # multinomial (K-forest softmax): the response is
  # an n x K count matrix, which the response vector every other family takes
  # has no shape for, so it rides the data object's own 'counts' argument -
  # where 'subset' reaches it, where the validity method constrains it, and
  # where a re-created or reloaded sampler finds it without further discipline.
  # A factor, character or integer-code response is the single-trial special
  # case, one-hot expanded here. The matrix interface only, as the survival
  # families are: a formula LHS carries no count matrix, and data@y is the
  # DERIVED trials vector rather than anything the caller wrote.
  multinomialCounts <- NULL
  # family = "auto" reads a count matrix (3+ columns of non-negative whole
  # numbers) as multinomial; announced with the other resolutions, below
  if (family == "auto" && directResponse && isAutoCountMatrix(data)) {
    family <- "multinomial"
  }
  if (identical(family, "multinomial") && !inherits(formula, "dbartsData")) {
    if (!directResponse) {
      stop(
        "multinomial fits currently use the matrix interface - ",
        "dbarts(x, y) or bart(x, y) with a factor ",
        "or an n x K count-matrix response"
      )
    }
    multinomialCounts <- resolveMultinomialCounts(data)
    matchedCall$data <- NULL
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
    control@call <- expandForwardedCall(matchedCall, evalEnv)
  }
  control@verbose <- verbose
  # a convenience mirror of dbartsControl(seed = ), as the wrappers expose;
  # an explicit seed overrides the control's, NULL or NA leaves it untouched
  seed <- resolveSeedArg(seed, "dbarts", refuse = TRUE)
  if (!is.na(seed)) {
    control@seed <- seed
  }
  # dec-B115: within-chain threading does not ship, so tree sampling never
  # sees more than one thread per chain; a caller-supplied budget above
  # n.chains is not silently wasted - it still feeds the test-fit pool and
  # predict's fan-out - but is worth a word since it buys nothing for the
  # sampler's own sweep. Single site: bart() forwards here with its own
  # control already built, so this fires for both front doors, once per
  # fit. The row count it names is the control's own testFitParallelCutoff,
  # so a caller who moved that cutoff is told the number actually in force.
  if (control@n.threads > control@n.chains) {
    warning(
      sprintf(
        paste0(
          "n.threads (%d) exceeds n.chains (%d); tree sampling uses at ",
          "most one thread per chain, so the extra threads serve only ",
          "test-set fitting above %d rows and predict"
        ),
        control@n.threads,
        control@n.chains,
        control@testFitParallelCutoff
      ),
      call. = FALSE
    )
  }

  dataCall <- redirectCall(matchedCall, quoteInNamespace(dbartsData))
  # a data object written in the call has been built once already, to be
  # looked at above: it is that object the fit uses, not a second build of
  # it, whose 'subset' could draw other rows
  if (inherits(formula, "dbartsData")) {
    dataCall$formula <- formula
  }
  # a basis declared on 'forests' has nowhere to ride once 'formula' is
  # already a built dbartsData: dbartsData() drops an unmatched 'bases'
  # argument in that case (its own ignored-args warning, R/data.R), which
  # would silently fit an ordinary single-forest model with the declaration
  # discarded - refused here by name instead. dbartsSpec() is not touched -
  # its first argument must already be a dbartsData, so this
  # predicate would be unconditionally true there and would refuse the
  # supported route (R/spec.R installs the declaration rather than dropping
  # it). The supported composition - a data object already carrying '@bases'
  # plus a knob-only 'forests' - is unaffected: forestBasisDeclarations()
  # returns NULL entries for a forest with no 'basis' of its own.
  basisDeclarations <- forestBasisDeclarations(forests)
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
  # The bases of the model's forests, whichever door declared them, are read
  # by one function against this fit's data (readForestBasis): a term's in
  # the formula's ingestion, where which forest has one decides which is
  # which, and a list's here. A basis covers every row of the data, and the
  # data object is the one place that knows which rows the fit keeps: it
  # evaluates 'subset' and applies the na.action, once. So every basis rides
  # its 'bases' argument and is cut there (validateForestBases, R/data.R). A
  # value rides as itself. A basis written as code is built only on the rows
  # kept, so the numbers of the data's rows ride in its place and come back
  # as the kept rows, on which it is built below.
  basisRows <- NULL
  basisReads <- NULL
  if (!is.null(termIngestion)) {
    basisReads <- termIngestion$basisReads
    basisRows <- termIngestion$basisRows
  } else if (!is.null(basisDeclarations)) {
    # missing() reads this frame, so it is resolved here. A hazard fit on the
    # matrix interface has already been expanded to its person-period rows,
    # 'subset' applied, and those are the rows a basis covers
    basisRows <- if (hazardExpandedFirst) {
      list(data = NULL, full = NROW(matchedCall$formula))
    } else {
      fitBasisRows(formula, if (missing(data)) NULL else data)
    }
    declared <- readDeclaredBases(forests, basisRows$data, basisRows$full)
    forests <- declared$forests
    basisReads <- declared$reads
  }
  if (!is.null(basisReads) && all(vapply(basisReads, is.null, logical(1L)))) {
    # a list in which no forest declares a basis is not a multi-forest
    # declaration at all - it names K ensembles with nothing to tell them
    # apart - so it falls through to resolveForests' own refusal by name
    basisReads <- NULL
  }
  basisIsCode <- vapply(
    basisReads,
    function(read) !is.null(read$frame),
    logical(1L)
  )
  if (!is.null(basisReads)) {
    dataCall$bases <- lapply(basisReads, function(read) {
      if (is.null(read)) {
        NULL
      } else if (is.null(read$frame)) {
        expandValueBasis(read$value)
      } else {
        basisRowNumbers(basisRows$full)
      }
    })
  }
  if (!is.null(multinomialCounts)) {
    dataCall$counts <- multinomialCounts
  }
  data <- withMatrixResponseRestated(
    "bart()/dbarts()",
    requestedFamily,
    if (is.null(dataCall$bases)) {
      withBinaryResponsePrecision(family, eval(dataCall, evalEnv))
    } else {
      # the bases ride the data object's own 'bases' argument, which this caller
      # never wrote, so a refusal from it is restated in the word the caller
      # used; only validateForestBases names that argument. An unrelated error -
      # anything not naming 'bases' - is not this call's to relabel, so it keeps
      # its own condition class and call
      tryCatch(
        withBinaryResponsePrecision(family, eval(dataCall, evalEnv)),
        error = function(e) {
          message <- conditionMessage(e)
          if (!grepl("'bases'", message, fixed = TRUE)) {
            stop(e)
          }
          stop(gsub("'bases'", "'basis'", message, fixed = TRUE), call. = FALSE)
        }
      )
    }
  )

  # a Surv formula response was ingested and subsetted by dbartsData()'s own
  # short-circuit (R/data.R), which has no family vocabulary to dispatch on -
  # decode its stashed status/time now, against the SAME conflict guard and
  # auto-dispatch-to-aft the direct-response form applies above
  formulaSurvivalStatus <- attr(data, "survivalStatus")
  hazardFormulaRecord <- NULL
  if (!is.null(formulaSurvivalStatus)) {
    formulaSurvivalTime <- attr(data, "survivalTime")
    formulaSurvivalTimeOmitted <- attr(data, "survivalTimeOmitted")
    attr(data, "survivalStatus") <- NULL
    attr(data, "survivalTime") <- NULL
    attr(data, "survivalTimeOmitted") <- NULL
    if (family %not_in% c("auto", "aft", hazardTokens)) {
      stop(
        "a survival (Surv) response cannot be fit with family \"",
        family,
        "\"; use family \"aft\", \"hazard\", or \"auto\""
      )
    }
    if (family %in% hazardTokens) {
      # a transformed or indicator-coded period predictor is read back as the
      # appended column at prediction, so any term that reads it is refused
      termVars <- unlist(lapply(
        attr(data@x, "term.labels"),
        function(label) all.vars(str2lang(label))
      ))
      if ("period" %in% termVars) {
        stop(
          "a hazard fit appends its own 'period' column; rename the ",
          "predictor 'period'"
        )
      }
      expansion <- expandDiscreteTimeHazard(
        data@x,
        formulaSurvivalTime,
        formulaSurvivalStatus,
        breaks = breaks,
        max.rows = max.rows,
        offset = data@offset,
        weights = data@weights
      )
      # makeModelMatrix already typed the original columns (categorical vs
      # ordinal); the appended period column is ordinal by construction
      # (dec-B97) and rides last, so one more entry keeps the two aligned
      data@x <- expansion$x
      data@varTypes <- c(data@varTypes, ORDINAL_VARIABLE)
      data@y <- expansion$y
      data@offset <- expansion$offset
      data@weights <- expansion$weights
      hazardPeriods <- expansion$periods
      K <- length(expansion$periods)
      # the subject names are the ones dbartsData recorded for the
      # unexpanded rows
      hazardNames <- hazardRowNames(
        dataRowNames(data, "train"),
        expansion$subject,
        dataRowNames(data, "test"),
        K
      )
      if (!is.null(data@na.action)) {
        omittedRows <- hazardOmittedRows(
          data@na.action,
          formulaSurvivalTimeOmitted,
          expansion,
          dataRowNames(data, "train")
        )
        hazardFormulaRecord <- omittedRows$record
        hazardNames$train <- omittedRows$names
      }
      # dbartsData()'s own 'test' handling already built data@x.test (and
      # its offset/weights twins) family-agnostically, coded against the
      # SAME pre-expansion training columns data@x just was - a held-out
      # subject has no event time to place it by, so it expands to every
      # one of the SAME K periods, exactly as the matrix interface's own
      # hazard 'test' acceptance does (R/dbarts.R's directResponse block)
      if (!is.null(data@x.test)) {
        n.test <- nrow(data@x.test)
        data@x.test <- appendHazardPeriodColumn(
          hazardRowSubset(data@x.test, rep(seq_len(n.test), times = K)),
          rep(seq_len(K), each = n.test)
        )
        if (!is.null(data@offset.test)) {
          offsetTestForExpansion <- data@offset.test
          if (length(offsetTestForExpansion) == 1L) {
            offsetTestForExpansion <- rep_len(offsetTestForExpansion, n.test)
          }
          data@offset.test <- rep(offsetTestForExpansion, times = K)
        }
        if (!is.null(data@weights.test)) {
          weightsTestForExpansion <- data@weights.test
          if (length(weightsTestForExpansion) == 1L) {
            weightsTestForExpansion <- rep_len(weightsTestForExpansion, n.test)
          }
          data@weights.test <- rep(weightsTestForExpansion, times = K)
        }
      }
      # the remap: the engine-facing family is now an ordinary binary link
      family <- if (identical(family, "hazard.logistic")) {
        "logistic"
      } else {
        "probit"
      }
    } else {
      family <- "aft"
      survivalStatus <- formulaSurvivalStatus
    }
  } else if (
    is.formula(formula) && (family == "aft" || family %in% hazardTokens)
  ) {
    stop(
      "family \"",
      family,
      "\" needs a survival::Surv or two-column (time, status) response"
    )
  }

  # a subject's offset is re-evaluated on a new subject, as predict and
  # survivalProbabilities form it, never the person-period vector the
  # expansion turned it into
  if (!is.null(hazardOffsetArgument)) {
    attr(data, "offset.argument") <- hazardOffsetArgument
  }
  if (!is.null(hazardNames)) {
    # the matrix interface expands before the na.action runs, so any rows it
    # dropped are person-period rows
    omitted <- data@na.action
    if (hazardExpandedFirst) {
      if (!is.null(omitted) && !is.null(hazardNames$train)) {
        names(omitted) <- hazardNames$train[unclass(omitted)]
        data@na.action <- omitted
        hazardNames$train <- hazardNames$train[-unclass(omitted)]
      }
    } else {
      # the formula path's na.action ran on the SUBJECT-level model frame,
      # before expansion; hazardOmittedRows has restated its record over
      # the person-period rows the dropped subjects would have had
      data@na.action <- hazardFormulaRecord
    }
    data <- setDataRowNames(data, "train", hazardNames$train)
    data <- setDataRowNames(data, "test", hazardNames$test)
  }

  data@n.cuts <- recycleNumCuts(control@n.cuts, ncol(data@x))
  data@sigma <- sigest

  # A basis written as code is built now, on the rows the data object kept,
  # whose numbers it handed back in the basis's place. data@bases is
  # positional against the forests, the forest with no basis first; the
  # values it cut are already in their places.
  basisRecords <- NULL
  if (any(basisIsCode)) {
    built <- buildFitBases(basisReads, data@bases, length(data@y))
    data@bases <- built$bases
    basisRecords <- built$records
  }

  # a term's predictors name design columns, which exist only now
  # (R/formulaTerms.R)
  if (!is.null(termIngestion)) {
    forests <- finalizeTermForests(termIngestion$forests, data)
  }

  # the matrix interface's own status vector rides outside (x, y) - it was
  # cut by 'subset' alone, above - while dbartsData() has just applied
  # na.action to x and the log-time response internally; mirror the same
  # drop here, exactly as a term's basis is restricted above, or a
  # predictor-NA row leaves the status vector one longer than data@y and
  # the C bridge's own length check refuses the fit
  if (directResponse && !is.null(survivalStatus) && !is.null(data@na.action)) {
    survivalStatus <- survivalStatus[-unclass(data@na.action)]
  }

  spec <- resolveSamplerSpec(
    matchedCall,
    formals(dbarts),
    control,
    data,
    family,
    requestedFamily = requestedFamily,
    shape = shape,
    residDf = residDf,
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
    evalEnv = evalEnv,
    residPrior = residPrior,
    familySpec = familySpec,
    basisRecords = basisRecords
  )

  sampler <- new("dbartsSampler", spec$control, spec$model, spec$data)
  # a latent family's 0/1 case weights are membership, which the sampler
  # carries as its active-row mask rather than as weights: the spec has
  # already cleared the weights slot, so this is the only place the vector
  # lands. All-ones never reaches here, having resolved to no mask at all.
  # No store: a fresh sampler's state stays the promise that captures the
  # state current at first read.
  if (!is.null(spec$active)) {
    sampler$setActiveRows(spec$active, updateState = FALSE)
  }
  sampler
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
  n.samples <- coerceOrError(n.samples, "integer")[1L]
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
  # under the caller's prior, anchored where its model records
  draw$model <- sampler$model
  applyAnchor(draw$pointer, sampler$model, FALSE)
  reissueNamedLeafSd(draw, draw$pointer)

  xt <- if (is.null(x.test)) extract(draw, "predictors") else x.test
  responseIsBinary <- draw$control@binary

  # the noise a heteroscedastic prior predictive adds is s(x) eps, with s^2(x)
  # drawn from the variance forest's own prior once per sample and read at the
  # rows being predicted; the scalar draws below would report a homoscedastic
  # prior predictive instead. "ev" needs no noise term and is unaffected.
  drawsVariance <- type == "ppd" &&
    !responseIsBinary &&
    !is.null(attr(draw$control, "bartcore.variance"))
  if (drawsVariance) {
    # the surface is read through the draw sampler's own test rows, which is
    # where the variance forest evaluates off the training data; predict()
    # cannot serve it, keepTrees being off and these trees never recorded
    draw$setTestPredictorAndOffset(xt, NULL)
  }

  sigmaDraws <- NULL
  if (type == "ppd" && !responseIsBinary && !drawsVariance) {
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
  varianceResults <- if (drawsVariance) vector("list", n.samples) else NULL
  for (i in seq_len(n.samples)) {
    draw$sampleTreesFromPrior(updateState = FALSE)
    draw$sampleLeafParametersFromPrior(updateState = FALSE)
    fit <- draw$predict(xt, offset.test, n.threads)
    # multi-chain samplers draw an independent prior stream per chain; prior
    # draws are chain-free, so only the first chain's stream is kept
    if (length(dim(fit)) > 1L) {
      fit <- fit[, 1L]
    }
    results[[i]] <- fit
    if (drawsVariance) {
      draw$sampleVarianceForestFromPrior(updateState = FALSE)
      varianceResults[[i]] <- draw$getVariance(test = TRUE)[, 1L]
    }
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
    } else if (drawsVariance) {
      # s^2(x) is a variance on the response scale, one row per prior draw and
      # one column per predicted row, so the noise scale is its square root
      result <- result +
        sqrt(do.call(rbind, varianceResults)) *
          matrix(
            rnorm(length(result)),
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

## Whether a method that just changed the sampler should refresh its cached
## state (the copy storeState/setState carry): an explicit TRUE or FALSE
## wins, and NULL - every one of these methods' own default - resolves
## against control@updateState, exactly as run() has always resolved it.
## Anything else is refused by name. 0.9-x's default was NA, which is a
## missing value now and reads as NULL, after one warning.
checkUpdateState <- function(updateState) {
  if (is.null(updateState)) {
    return(NULL)
  }
  if (isSingleNA(updateState)) {
    warnNAForNull("updateState", "dbartsSampler")
    return(NULL)
  }
  if (!is.logical(updateState) || length(updateState) != 1L) {
    stop("'updateState' must be TRUE, FALSE or NULL", call. = FALSE)
  }
  updateState
}

## A run's burn-in or sample count: NULL takes the control's, which the
## engine layer reads as NA_integer_. 0.9-x documented NA for the same thing,
## and it reads that way for the release after one warning.
resolveRunCount <- function(count, argument) {
  if (is.null(count)) {
    return(NA_integer_)
  }
  refuseNaN(count, argument)
  if (isSingleNA(count)) {
    warnNAForNull(argument, "dbartsSampler")
    return(NA_integer_)
  }
  count
}

resolveUpdateState <- function(updateState, control) {
  updateState <- checkUpdateState(updateState)
  if (is.null(updateState)) control@updateState else updateState
}

## NULL is the one way to say "no test offset", as it is for lm's offset; an
## NA is a missing value, and nothing on these paths routes missing rows.
refuseMissingTestOffset <- function(offset.test) {
  if (anyNA(offset.test)) {
    stop(
      "'offset.test' contains missing values; use NULL for no offset",
      call. = FALSE
    )
  }
  invisible(NULL)
}

## The sampler's predict after validation: 'x.test' is already coded by
## validateXTest, so the predict methods that validate newdata themselves
## (and resolve its na.action) call this directly, and no warning fires twice.
predictCodedTest <- function(sampler, x.test, offset.test, n.threads) {
  # a sparse-backed test set rides to the engine as the container
  # validateXTest coded it; the engine routes its rows off that storage

  # A multinomial predict surface reports K probabilities per row, so its
  # offset takes the shape of that surface: the per-category matrix
  # entering the raw fits BEFORE the softmax, one row per PREDICTED row.
  # A flat vector stays refused there, and truthfully - after the blend it
  # would move the values off the simplex, and before it a common
  # per-observation shift is the softmax's own null direction. The rows are
  # the caller's, so a sampler holding either resident category offset
  # refuses a no-offset call rather than reporting the offset-free surface
  # (an all-zero matrix asks for that surface on purpose).
  counts <- dataCounts(sampler$data)
  if (!is.null(offset.test) && !is.null(counts)) {
    if (length(offset.test) == 1L) {
      offset.test <- matrix(
        as.double(offset.test),
        nrow(x.test),
        ncol(counts)
      )
    } else {
      offset.test <- as.matrix(offset.test)
      storage.mode(offset.test) <- "double"
    }
    refuseMissingTestOffset(offset.test)
    if (!identical(dim(offset.test), c(nrow(x.test), ncol(counts)))) {
      stop(
        "'offset.test' must be a per-category matrix with one row per ",
        "row of 'x.test' and ",
        ncol(counts),
        " categories"
      )
    }
  } else if (!is.null(offset.test)) {
    offset.test <- as.double(offset.test)
    refuseMissingTestOffset(offset.test)
    if (length(offset.test) == 1L) {
      offset.test <- rep_len(offset.test, nrow(x.test))
    }

    if (!identical(length(offset.test), nrow(x.test))) {
      stop(
        "'offset.test' must have the same number of rows as 'x.test'"
      )
    }
  }

  .Call(
    C_dbarts_bartcore_predict,
    sampler$getPointer(),
    x.test,
    offset.test,
    n.threads
  )
}

## The same for predictForests.
predictForestsCodedTest <- function(sampler, x.test, offset.test, n.threads) {
  .Call(
    C_dbarts_bartcore_predictPerForest,
    sampler$getPointer(),
    x.test,
    offset.test,
    n.threads
  )
}

## A named leaf-prior sd is absolute on the family's scale, but the engine
## holds the leaf scale against the response transform in force; a channel
## that re-anchors that transform, or an install that moves a re-created
## sampler into its recorded one, would carry the sd with it. The named scale
## is written back after each one, so the sd means what it did. Only a
## single-forest model can name one, and a write equal to what is in force is
## skipped inside the engine. An install passes the pointer it used, which
## getPointer has not yet bound.
reissueNamedLeafSd <- function(sampler, ptr = sampler$getPointer()) {
  anchor <- sampler$model@prior.scale
  if (is.na(anchor)) {
    return(invisible(NULL))
  }
  .Call(C_dbarts_bartcore_setLeafPrior, ptr, 0L, anchor)
  invisible(NULL)
}

## The response transform a sampler's leaf prior is anchored to and its
## chains hold their numbers in - (min, max) as a state's fit.scale holds it
## - is model, recorded on the model as the "response.range" attribute, where
## a model saved before the record existed reads NULL. A first creation and
## the re-anchoring channels write it from the engine.
recordAnchor <- function(model, ptr) {
  attr(model, "response.range") <- .Call(
    C_dbarts_bartcore_anchor,
    ptr,
    NULL,
    FALSE
  )
  model
}

## A re-creation from a sampler's own model takes its record as the engine's
## transform. The chains are moved there at once unless an install follows,
## whose state is converted into the record and moves them itself, exactly as
## the install moved a re-created sampler before the record existed.
applyAnchor <- function(ptr, model, installFollows) {
  record <- attr(model, "response.range", exact = TRUE)
  if (!is.null(record)) {
    .Call(C_dbarts_bartcore_anchor, ptr, record, installFollows)
  }
  invisible(ptr)
}

recreatePointer <- function(control, model, data, installFollows) {
  ptr <- .Call(
    C_dbarts_bartcore_create,
    control,
    model,
    data,
    if (model@family == "auto") "" else model@family
  )
  applyAnchor(ptr, model, installFollows)
}

## What prior.sd is the sd of, per leaf model.
priorSdOf <- function(leafModel) {
  switch(
    leafModel,
    linear = "coefficient",
    gp = "amplitude",
    "leaf value"
  )
}

## One forest's leaf prior, from the bridge's per-chain calibration: the
## specification in the terms it was named in, then the quantities the chains
## share. A fixed value is read off the engine, so it is what is in force; a
## law comes from the model, gated by the engine's per-forest flag, since a
## map forest pins k whatever the model says. The map entries are present only
## on a forest whose scale the map sets, and there the specification is the
## forest(sd = ) creation takes: the half-Cauchy median on a scale-mixture
## forest, the leaf-scale factor otherwise. Every chain runs under the
## sampler's one prior and transform, which no install moves, so the first
## chain's reading is the sampler's.
reportLeafPrior <- function(sampler, raw) {
  shared <- function(column) raw[[1L, column]]
  model <- sampler$model
  hyperprior <- model@leaf.hyperprior
  mapped <- !is.nan(raw[1L, "basis.row.norm"])
  drawn <- raw[1L, "k.has.hyperprior"] != 0
  anchor <- shared("prior.scale")
  spec <- model@leaf.prior
  sdNamed <- mapped || !is.na(model@prior.scale) || !is.null(spec@prior.sd)
  spec@k <- NULL
  spec@prior.sd <- NULL
  if (!sdNamed) {
    spec@k <- if (drawn) hyperprior else shared("k")
  } else if (!drawn) {
    spec@prior.sd <- anchor / shared("k")
  } else {
    spec@prior.sd <- invchi(
      hyperprior@degreesOfFreedom,
      if (is.finite(hyperprior@scale)) anchor / hyperprior@scale else 0
    )
  }
  leafModel <- attr(raw, "leaf.model")
  prior <- list(
    leaf.prior = spec,
    leaf.model = leafModel,
    prior.sd.of = priorSdOf(leafModel),
    prior.mean = shared("prior.mean"),
    k.scale = anchor,
    response.scale = shared("response.scale"),
    response.shift = shared("response.shift")
  )
  if (!mapped) {
    return(prior)
  }
  mixture <- is.nan(raw[1L, "amplitude.prior.variance"])
  amplitude <- if (mixture) {
    "amplitude.prior.scale"
  } else {
    "amplitude.prior.variance"
  }
  sd <- shared(if (mixture) "amplitude.prior.scale" else "leaf.scale.factor")
  prior$leaf.prior <- forest(sd = sd)
  prior$prior.sd.of <- if (mixture) "amplitude scale" else "forest total"
  for (column in c(amplitude, "leaf.scale.factor", "leaf.scale.divisor")) {
    prior[[column]] <- shared(column)
  }
  prior$basis.row.norm <- shared("basis.row.norm")
  prior
}

## Installs a restated leaf prior. A new anchor under the hyperprior in force
## is the engine's own leaf-scale write, which skips a value equal to what is
## in force; a new hyperprior, or a return to the data's scale, goes through
## the model, whose install re-pins a fixed sigma, so that is put back.
writeLeafPrior <- function(sampler, ptr, newModel) {
  oldModel <- sampler$model
  sameLaw <- identical(newModel@leaf.hyperprior, oldModel@leaf.hyperprior)
  if (sameLaw && !is.na(newModel@prior.scale)) {
    .Call(C_dbarts_bartcore_setLeafPrior, ptr, 0L, newModel@prior.scale)
    return(invisible(NULL))
  }
  if (sameLaw && is.na(oldModel@prior.scale)) {
    return(invisible(NULL))
  }
  sigmas <- .Call(C_dbarts_bartcore_getSigmas, ptr)
  .Call(
    C_dbarts_bartcore_setModel,
    ptr,
    newModel,
    sampler$data,
    sampler$control
  )
  if (
    !identical(.Call(C_dbarts_bartcore_getSigmas, ptr), sigmas) &&
      length(unique(sigmas)) == 1L
  ) {
    .Call(C_dbarts_bartcore_setSigma, ptr, sigmas[[1L]])
  }
  invisible(NULL)
}

## The prior vocabulary $setLeafPrior evaluates its specification in: the
## constructors' own, except that linear() and gp() default their columns to
## the sampler's, which fixes them, so a write may omit them. On a sampler of
## another leaf model the placeholder only lets the leaf-model refusal speak.
setLeafPriorVocabulary <- function(sampler) {
  current <- sampler$model@leaf.prior
  ownColumns <- if (
    is(current, "dbartsLinearPrior") || is(current, "dbartsGPPrior")
  ) {
    current@columns
  } else {
    1L
  }
  vocabulary <- dbartsPriors
  vocabulary$linear <- function(columns = ownColumns, k = NULL, sd = NULL) {
    linear(columns, k, sd)
  }
  vocabulary$gp <- function(
    columns = ownColumns,
    k = NULL,
    lengthscale = NULL,
    max.leaf.size = 256L,
    sd = NULL
  ) {
    gp(columns, k, lengthscale, max.leaf.size, sd)
  }
  vocabulary
}

## The model a $setLeafPrior write installs: the sampler's own, with the leaf
## prior's k or sd, and so its anchor and hyperprior, taken from the
## specification. The specification must name the sampler's leaf model, and
## any leaf-model detail it states must match, since only $setModel changes
## the leaf model.
restateLeafPrior <- function(sampler, spec, expr) {
  model <- sampler$model
  current <- model@leaf.prior
  constructor <- function(prior) {
    if (is(prior, "dbartsLinearPrior")) {
      "linear()"
    } else if (is(prior, "dbartsGPPrior")) {
      "gp()"
    } else {
      "normal()"
    }
  }
  refuseLeafModel <- function(what) {
    stop(
      "$setLeafPrior changes only the leaf prior's spread or its ",
      "hyperprior, and ",
      what,
      "; change the leaf model through $setModel",
      call. = FALSE
    )
  }
  if (!identical(class(spec), class(current))) {
    refuseLeafModel(paste0(
      "this sampler's leaf model is written ",
      constructor(current),
      ", not ",
      constructor(spec)
    ))
  }
  # a detail is compared only where the call states it; a prior passed as a
  # value states everything that is not at its constructor's default
  supplied <- if (is.call(expr) && is.name(expr[[1L]])) {
    constructorName <- as.character(expr[[1L]])
    if (constructorName %in% c("linear", "gp")) {
      names(match.call(dbartsPriors[[constructorName]], expr))[-1L]
    }
  }
  stated <- function(slot, default) {
    if (!is.null(supplied)) {
      return(slot %in% supplied)
    }
    !identical(methods::slot(spec, slot), default)
  }
  if (is(spec, "dbartsLinearPrior") || is(spec, "dbartsGPPrior")) {
    resolved <- resolveLeafCovariates(spec, sampler$data)
    if (!identical(resolved@columns, current@columns)) {
      refuseLeafModel("the leaf covariate columns it names differ")
    }
    if (
      is(spec, "dbartsGPPrior") &&
        stated("lengthscale", NULL) &&
        !identical(resolved@lengthscale, current@lengthscale)
    ) {
      refuseLeafModel("the lengthscale it names differs")
    }
    if (
      is(spec, "dbartsGPPrior") &&
        stated("max.leaf.size", 256L) &&
        spec@max.leaf.size != current@max.leaf.size
    ) {
      refuseLeafModel("the max.leaf.size it names differs")
    }
  }
  translated <- resolveLeafPrior(
    spec,
    drawsLeafKByDefault(sampler$model@family),
    monotone = !is.null(attr(model, "monotone"))
  )
  current@k <- spec@k
  current@prior.sd <- spec@prior.sd
  model@leaf.prior <- current
  model@leaf.hyperprior <- translated$leaf.hyperprior
  model@prior.scale <- translated$prior.scale
  model
}

## The spreads a $setLeafPrior(forests = ) call restates on a sampler whose
## forests carry amplitudes: forest(sd = ) in creation's positions, validated
## whole before anything is written. A NULL sd leaves its forest; every other
## knob is fixed at creation. Returns the per-forest sd, NA where none is
## stated.
resolveForestSpreads <- function(sampler, forests) {
  forestInfo <- attr(sampler$control, "bartcore.forests", exact = TRUE)
  numForests <- length(forestInfo$params)
  if (
    !is.list(forests) ||
      !all(vapply(forests, inherits, logical(1L), "dbartsForest"))
  ) {
    stop(
      "$setLeafPrior's 'forests' must be a list of forest() specifications, ",
      "as at creation"
    )
  }
  if (length(forests) == 0L || length(forests) > numForests) {
    stop(
      "$setLeafPrior's 'forests' names ",
      length(forests),
      " forests; this sampler has ",
      numForests
    )
  }
  given <- names(forests)
  if (!is.null(given)) {
    labels <- forestInfo$labels
    if (is.null(labels)) {
      labels <- rep("", numForests)
    }
    labels <- labels[seq_along(given)]
    mismatched <- which(nzchar(given) & given != labels)
    if (length(mismatched) > 0L) {
      index <- mismatched[[1L]]
      stop(
        "$setLeafPrior's 'forests' names forest ",
        index,
        " '",
        given[[index]],
        "', but it was created ",
        if (nzchar(labels[[index]])) {
          paste0("as '", labels[[index]], "'")
        } else {
          "unnamed"
        }
      )
    }
  }
  vapply(
    seq_along(forests),
    function(index) {
      spec <- forests[[index]]
      for (knob in setdiff(names(spec), "sd")) {
        if (!is.null(spec[[knob]])) {
          stop(
            "$setLeafPrior restates only a forest's 'sd': '",
            knob,
            "' is fixed at creation",
            if (knob == "basis") "; change it with $setForestBasis"
          )
        }
      }
      if (is.null(spec$sd)) {
        return(NA_real_)
      }
      validateForestSd(spec$sd)
    },
    numeric(1L)
  )
}

## Writes resolved per-forest spreads, then mirrors each into the control
## attribute creation reads, in the channel the forest was created in: the
## half-Cauchy median when it carries one, the leaf-scale factor otherwise.
## Every re-creation then builds with the write.
writeForestSpreads <- function(sampler, ptr, sds) {
  forestInfo <- attr(sampler$control, "bartcore.forests", exact = TRUE)
  for (index in which(!is.na(sds))) {
    .Call(C_dbarts_bartcore_setForestSd, ptr, index - 1L, sds[[index]])
    params <- forestInfo$params[[index]]
    params[[if (params[[7L]] > 0) 7L else 4L]] <- sds[[index]]
    forestInfo$params[[index]] <- params
  }
  newControl <- sampler$control
  attr(newControl, "bartcore.forests") <- forestInfo
  sampler$control <- newControl
  invisible(NULL)
}

## What a multinomial sampler's creation refuses, refused again on a write:
## its $setLeafPrior takes normal(k = ) with a fixed k, and nothing else.
refuseMultinomialLeafPrior <- function(spec) {
  reason <- if (!is.null(spec@prior.sd)) {
    paste0(
      "the softmax calibration map sets every category forest's leaf scale, ",
      "so a named 'sd' has nowhere to land"
    )
  } else if (is(spec@k, "dbartsLeafHyperprior")) {
    "a 'k' hyperprior is not supported on a multinomial sampler, at creation or after"
  }
  if (!is.null(reason)) {
    stop(multinomialLeafPriorMessage, reason, call. = FALSE)
  }
  invisible(NULL)
}

multinomialLeafPriorMessage <- paste0(
  "$setLeafPrior on a multinomial sampler takes normal(k = ) with a fixed k, ",
  "as its creation does: "
)

dbartsSampler <- setRefClass(
  "dbartsSampler",
  fields = list(
    pointer = "externalptr",
    control = "dbartsControl",
    model = "dbartsModel",
    data = "dbartsData",
    state = "ANY", # is either a list of states, or a promise to evaluate
    # The per-forest, per-observation precision weight installed by
    # setForestWeights, mirrored here because it does not ride the engine's
    # saved state: forestWeights[[forest]] (1-based) holds the last vector
    # installed on that forest, NULL where none is. getPointer, setState and
    # copy all re-apply it on every re-creation.
    forestWeights = "list",
    # The active-row mask installed by setActiveRows, mirrored here for the
    # same reason and re-applied on the same paths: NULL where no mask is in
    # force, an all-ones vector installing none. Without the mirror a masked
    # sampler re-created from its stored state - what getPointer does after a
    # save and load - would silently return every masked row to the
    # likelihood.
    activeRows = "ANY"
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
      # a model handed in directly is held to what the fitting functions hold
      # the one they resolve to
      refuseNoSplittableColumn(
        model@tree.prior@splitProbabilities,
        attr(model, "forest.columns", exact = TRUE)
      )
      .self$control <- control
      .self$model <- model
      .self$data <- data
      .self$forestWeights <- list()
      .self$activeRows <- NULL

      # "auto" (a hand-built model) keeps the bridge's own dispatch
      .self$pointer <- .Call(
        C_dbarts_bartcore_create,
        .self$control,
        .self$model,
        .self$data,
        if (model@family == "auto") "" else model@family
      )
      # a first creation anchors to its own data, whatever record the model
      # handed in carries; copy and the re-creations restate the saver's
      .self$model <- recordAnchor(model, .self$pointer)
      # the calibration map's anchor s, recorded at first creation so every
      # re-creation builds on it rather than on the response then in force; a
      # copy arrives with it already recorded
      forestInfo <- attr(control, "bartcore.forests", exact = TRUE)
      if (!is.null(forestInfo) && is.null(forestInfo$anchor)) {
        forestInfo$anchor <- .Call(
          C_dbarts_bartcore_getLeafPrior,
          .self$pointer,
          0L
        )[[1L, "map.anchor"]]
        attr(control, "bartcore.forests") <- forestInfo
        .self$control <- control
      }
      # after the call, so a refused creation never spends the warning's key;
      # re-creation from state (getPointer, setState) does not come through
      # here, and copy holds it off
      warnZeroTrials(dataCounts(data))
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
      numBurnIn = NULL,
      numSamples = NULL,
      updateState = NULL,
      ...,
      callback = NULL
    ) {
      "Runs the posterior sampler and returns a list with the results. NULL burn-in or sample counts take the control's."
      updateState <- checkUpdateState(updateState)
      ignoreRunThreadCount(...)
      numBurnIn <- resolveRunCount(numBurnIn, "numBurnIn")
      numSamples <- resolveRunCount(numSamples, "numSamples")

      samples <- bartcoreSamplerRun(.self, numBurnIn, numSamples, callback)
      if (resolveUpdateState(updateState, control)) {
        storeState()
      }
      if (is.null(samples)) {
        return(invisible(NULL))
      }
      samples
    },
    sampleTreesFromPrior = function(updateState = NULL) {
      "Draws tree structure from prior"
      updateState <- checkUpdateState(updateState)
      ptr <- getPointer()
      .Call(C_dbarts_bartcore_sampleTreesFromPrior, ptr)

      if (resolveUpdateState(updateState, control)) {
        storeState(ptr)
      }

      invisible(NULL)
    },
    sampleLeafParametersFromPrior = function(updateState = NULL) {
      "Draws leaf values from their prior; does not change tree structure."
      updateState <- checkUpdateState(updateState)
      ptr <- getPointer()
      .Call(C_dbarts_bartcore_sampleLeafParametersFromPrior, ptr)

      if (resolveUpdateState(updateState, control)) {
        storeState(ptr)
      }

      invisible(NULL)
    },
    sampleNodeParametersFromPrior = function(updateState = NULL) {
      "Retired: use $sampleLeafParametersFromPrior. Forwards for one release."
      updateState <- checkUpdateState(updateState)
      warnOnce(
        "tombstone.sampleNodeParametersFromPrior",
        "'$sampleNodeParametersFromPrior' is now ",
        "'$sampleLeafParametersFromPrior'; this call was forwarded. The old ",
        "name is removed in dbarts ",
        tombstoneExpiry,
        ".",
        class = "dbartsDeprecatedWarning"
      )
      sampleLeafParametersFromPrior(updateState)
    },
    sampleVarianceForestFromPrior = function(updateState = NULL) {
      "Draws the variance forest's tree structures and leaf factors from their priors; a no-op on a homoscedastic sampler."
      updateState <- checkUpdateState(updateState)
      ptr <- getPointer()
      .Call(C_dbarts_bartcore_sampleVarianceForestFromPrior, ptr)

      if (resolveUpdateState(updateState, control)) {
        storeState(ptr)
      }

      invisible(NULL)
    },
    growFromRoot = function(n.sweeps = 2L, updateState = NULL) {
      "Builds an initial forest by XBART-style grow-from-root (He, Yalov and Hahn 2019) as a warm start, running n.sweeps grow sweeps in place; the exact MCMC sampler owns the forest once run() begins. Constant-leaf models only. See ?dbartsSampler."
      updateState <- checkUpdateState(updateState)
      if (
        is(model@leaf.prior, "dbartsLinearPrior") ||
          is(model@leaf.prior, "dbartsGPPrior")
      ) {
        stop(
          "grow-from-root warm start is only available for the constant-leaf ",
          "model; linear and gp leaf priors initialize with ",
          "sampleTreesFromPrior instead"
        )
      }
      n.sweeps <- coerceOrError(n.sweeps, "integer")
      if (length(n.sweeps) != 1L || is.na(n.sweeps) || n.sweeps <= 0L) {
        stop("'n.sweeps' must be a single positive integer")
      }
      ptr <- getPointer()
      .Call(C_dbarts_bartcore_growFromRoot, ptr, n.sweeps)

      if (resolveUpdateState(updateState, control)) {
        storeState(ptr)
      }

      invisible(NULL)
    },
    copy = function(shallow = FALSE) {
      "Creates a deep or shallow copy of the sampler, keeping its model and installing its stored state."
      # a copy introduces no rows, so it does not raise the creation warning
      dupe <- withoutZeroTrialsWarning(
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
      )

      # a copy is a re-creation: it keeps this sampler's model, the record of
      # its transform included, and the install, or with none the move here,
      # puts the chains in it. The stored state is opaque and never mutated
      # in place (storeState replaces it whole), so the copy can install the
      # same object.
      dupe$model <- model
      applyAnchor(dupe$pointer, model, !is.null(state))
      if (!is.null(state)) {
        dupe$setState(state)
      } else {
        reissueNamedLeafSd(dupe, dupe$pointer)
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
      # the mask rides neither the state nor the data object either, so the
      # copy takes it from the same mirror by the same route
      dupe$activeRows <- activeRows
      dupe$reapplyActiveRows(dupe$pointer)
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
      "Using existing sampler to predict for new data without re-running. n.threads is a per-call worker count that does not persist, defaulting to the sampler's own: the replay is partitioned by (chain, saved draw), each partition writing its own rows, so the answer is identical bit for bit at every value."
      x.test <- validateXTest(x.test, data@x)
      if (is.null(x.test)) {
        stop("x.test cannot be NULL")
      }
      predictCodedTest(
        .self,
        x.test,
        if (!missing(offset.test)) offset.test,
        n.threads
      )
    },
    predictForests = function(
      x.test,
      offset.test,
      n.threads = control@n.threads
    ) {
      "Replays each forest separately at new data, without re-running: an n.new x n.forests x n.samples (x n.chains) array of each forest's own INTERNAL-scale total, the off-sample twin of getForestFits. Only a sampler that composes its forests through scalar amplitude glue reports per-forest fits; every other one, a multinomial sampler included, is refused by name. No glue, no response transform and no offset are folded in: the location an amplitude coupling reports is response.shift + sum_f (basis_f %*% glue_f) * (response.scale * f_f), and off the training rows the bases are the caller's, so the whole recombination is too. offset.test is refused for the same reason - a shift belongs to that recombination. Reports the saved samples under keepTrees, and otherwise the current trees, exactly as predict does. n.threads is predict's per-call worker count, with the same (chain, saved draw) partition and the same bitwise-identical result at every value."
      x.test <- validateXTest(x.test, data@x)
      if (is.null(x.test)) {
        stop("x.test cannot be NULL")
      }
      predictForestsCodedTest(
        .self,
        x.test,
        if (!missing(offset.test)) offset.test,
        n.threads
      )
    },
    setControl = function(newControl) {
      "Sets the control object for the sampler to a new one. Preserves the call() slot and any bartcore.* control attributes."
      if (!inherits(newControl, "dbartsControl")) {
        stop("'control' must inherit from dbartsControl")
      }

      selfEnv <- parent.env(environment())

      newControl@binary <- control@binary
      newControl@call <- control@call
      # a control taken from another sampler brings that sampler's attributes;
      # only this sampler's own are carried
      for (attrName in grep(
        "^bartcore\\.",
        names(attributes(newControl)),
        value = TRUE
      )) {
        attr(newControl, attrName) <- NULL
      }
      # bartcore.* attributes (the BCF, variance, survival and
      # ordinal/nbinom configuration resolveSamplerSpec attaches at creation)
      # live outside the S4 slots newControl replaces
      # wholesale; a freshly built dbartsControl() never carries them, so
      # carry them forward exactly as binary/call are. Without this, a
      # legitimate setControl silently orphans the model configuration the
      # stored state's forest counts were built against, and the next
      # getPointer() re-creation refuses loudly instead (a multi-forest
      # sampler: "the data carry forest bases but no basis forest was
      # configured"; a heteroscedastic one: "state is not consistent with this
      # sampler").
      for (attrName in grep(
        "^bartcore\\.",
        names(attributes(control)),
        value = TRUE
      )) {
        attr(newControl, attrName) <- attr(control, attrName)
      }

      # settings fixed at creation: the generators, anything shaping the cut
      # grid, and the four engine limits, which the sampler reads once when it
      # is created and never again. Accepting one here would leave the stored
      # control disagreeing with the engine, and a re-creation from that stored
      # control would then move the draws under categoricalExhaustiveCap.
      # proposal.probs is NOT among them: the mixture is installed with the
      # priors and has always been changeable mid-run through $setModel, so it
      # is honored here through that same install (below).
      for (slotName in c(
        "n.trees",
        "n.chains",
        "useQuantiles",
        "levelGibbs",
        "categoricalExhaustiveCap",
        "testFitParallelCutoff",
        "predictParallelCutoff",
        "sparseDensityThreshold",
        "seed"
      )) {
        if (
          !identical(
            methods::slot(newControl, slotName),
            methods::slot(control, slotName)
          )
        ) {
          stop(
            "changing '",
            # the slot keeps the bridge's name; the argument is treeShift
            if (slotName == "levelGibbs") "treeShift" else slotName,
            "' is not available on an existing sampler"
          )
        }
      }
      if (newControl@keepTrees && is.na(newControl@n.samples)) {
        stop("keepTrees requires 'n.samples' to be specified")
      }

      # a monotone sampler is birth/death-only, as creation forces: a defaulted
      # mixture is rewritten and one proposing other moves refused, before
      # anything is installed
      if (!is.null(attr(model, "monotone"))) {
        newControl@proposal.probs <- monotoneProposalProbs(
          newControl@proposal.probs,
          allowBirthDeath = TRUE
        )
      }

      mixtureMoved <- !identical(
        newControl@proposal.probs,
        control@proposal.probs
      )

      ptr <- getPointer()
      oldControl <- control
      # the engine takes the control first: a refusal there leaves the stored
      # control the one the engine still has
      .Call(C_dbarts_bartcore_setControl, ptr, newControl)
      selfEnv$control <- newControl
      # the engine reads the mixture off the control when the priors are
      # installed, so a changed one is pushed through the prior install and
      # meets every refusal that install already carries. A refusal rolls the
      # whole control back: the stored one must never name a mixture the
      # engine does not have.
      if (mixtureMoved) {
        tryCatch(
          .self$setModel(model),
          error = function(e) {
            selfEnv$control <- oldControl
            .Call(C_dbarts_bartcore_setControl, ptr, oldControl)
            stop(e)
          }
        )
      }

      invisible(NULL)
    },
    setModel = function(newModel) {
      "Sets the model object for the sampler to a new one."
      refuseCountsMutation(
        .self,
        "$setModel",
        "every category forest's prior is calibrated at creation from the ",
        "softmax map; make a new sampler instead"
      )
      if (!inherits(newModel, "dbartsModel")) {
        stop("'model' must inherit from dbartsModel")
      }
      refuseInvalidLeafPrior(newModel@leaf.prior)
      refuseAmplitudeMutation(
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
      # the columns a forest may split on are structure: its trees were drawn
      # under them, so a model stating others, or none, is another model
      if (
        !identical(
          attr(newModel, "forest.columns", exact = TRUE),
          attr(model, "forest.columns", exact = TRUE)
        )
      ) {
        stop(
          "$setModel cannot change the columns a forest may split on, its ",
          "'vars': they are fixed when a sampler is created; make a new ",
          "sampler instead"
        )
      }
      # split probabilities are a parameter, held to what creation holds them
      # to on a restricted forest
      refuseNoSplittableColumn(
        newModel@tree.prior@splitProbabilities,
        attr(model, "forest.columns", exact = TRUE)
      )
      ptr <- getPointer()
      selfEnv <- parent.env(environment())

      newModel@family <- model@family
      # the transform the prior is anchored to stays the sampler's
      attr(newModel, "response.range") <- attr(model, "response.range")
      oldModel <- model
      selfEnv$model <- newModel
      tryResult <- tryCatch(
        .Call(C_dbarts_bartcore_setModel, ptr, selfEnv$model, data, control),
        error = function(e) {
          selfEnv$model <- oldModel
          e$call <- quote(.Call(
            C_dbarts_bartcore_setModel,
            ptr,
            selfEnv$model,
            data,
            control
          ))
          e
        }
      )
      if (inherits(tryResult, "error")) {
        stop(tryResult)
      }

      invisible(NULL)
    },
    setData = function(newData, updateState = NULL) {
      "Sets the data object for the sampler to a new one. Preserves the n.cuts and sigma slots. updateState follows control@updateState: NULL, its default, resolves to the control's setting, and an explicit TRUE or FALSE overrides it - the same rule run() applies."
      updateState <- checkUpdateState(updateState)
      refuseCountsMutation(
        .self,
        "$setData",
        "its K category forests fix their data at creation, and the count ",
        "response is not a column of it; make a new sampler instead"
      )
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
      selfEnv <- parent.env(environment())
      selfEnv$model <- recordAnchor(model, getPointer())
      reissueNamedLeafSd(.self)
      if (resolveUpdateState(updateState, control)) {
        storeState()
      }
      invisible(NULL)
    },
    setResponse = function(
      y,
      updateScale = FALSE,
      updateState = NULL,
      status = NULL
    ) {
      "Changes the response against which the sampler is fitted, and, for an aft (survival) sampler given a non-null status, its censoring structure in the same call. updateState follows control@updateState; see setData."
      updateState <- checkUpdateState(updateState)
      # a caller porting $setResponse(y, updateState) from before the
      # updateScale/updateState reorder gets the same TRUE/FALSE/NA in the
      # same position, now meaning updateScale; sys.call() carries the raw,
      # unmatched call, so two unnamed arguments are caught before R's own
      # positional matching resolves them silently. Naming 'updateScale'
      # anywhere in the call - s$setResponse(updateScale = FALSE, y) - is
      # never the ported shape, whatever position its own unnamed argument
      # then falls into, so it is excluded rather than counted.
      rawCall <- sys.call()
      rawCallArgNames <- names(rawCall)[-1L]
      if (is.null(rawCallArgNames)) {
        rawCallArgNames <- character(length(rawCall) - 1L)
      }
      unnamedArgs <- sum(!nzchar(rawCallArgNames))
      if (unnamedArgs >= 2L && "updateScale" %not_in% rawCallArgNames) {
        # warnOnce's session-scoped key dedupes repeated calls inside one loop
        warnOnce(
          "setResponsePositionalUpdateScale",
          paste0(
            "the second argument to $setResponse is 'updateScale' in dbarts ",
            ">= 1.0-0, and 'updateState' has moved to third; pass both by ",
            "name to avoid depending on this order"
          )
        )
      }
      refuseCountsMutation(
        .self,
        "$setResponse",
        "its response is the n x K count matrix, which a flat vector cannot ",
        "express and which a length-n integer vector would only be guessed ",
        "into; replace it with $setCounts"
      )
      bartcoreSamplerSetResponse(.self, y, updateScale, status)
      if (isTRUE(updateScale)) {
        selfEnv <- parent.env(environment())
        selfEnv$model <- recordAnchor(model, getPointer())
        reissueNamedLeafSd(.self)
      }
      if (resolveUpdateState(updateState, control)) {
        storeState()
      }
      invisible(NULL)
    },
    setOffset = function(offset, updateScale = FALSE, updateState = NULL) {
      "Changes the offset slot used to adjust the response. updateState follows control@updateState; see setData."
      updateState <- checkUpdateState(updateState)
      refuseCountsMutation(
        .self,
        "$setOffset",
        "a common per-observation shift is the softmax's own null direction, ",
        "so a flat offset is inert; the shift that is not is the n x K ",
        "matrix $setCategoryOffset takes"
      )
      bartcoreSamplerSetOffset(.self, offset, updateScale)
      if (isTRUE(updateScale)) {
        selfEnv <- parent.env(environment())
        selfEnv$model <- recordAnchor(model, getPointer())
        reissueNamedLeafSd(.self)
      }
      if (resolveUpdateState(updateState, control)) {
        storeState()
      }
      invisible(NULL)
    },
    setWeights = function(weights, updateState = NULL) {
      "Changes the weights with which the sampler is fitted. A row of weight 0 leaves the likelihood but stays in the fit, occupying a leaf. A probit, ordinal or nbinom sampler carries no weight channel, and takes only weights of 0 and 1: those name the rows in its data set, so they install as the active-row mask (see setActiveRows) and the data object's weights slot stays empty. A Student-t sampler redraws the scale of each row whose weight leaves zero, from the sampler's own generators. updateState follows control@updateState; see setData."
      updateState <- checkUpdateState(updateState)
      refuseCountsMutation(
        .self,
        "$setWeights",
        "an integer case weight is already row-wise replication in its count ",
        "response, and a non-integer one has no exact augmentation sampler"
      )
      weights <- as.double(weights)
      if (length(weights) != length(data@y)) {
        stop("'weights' must have the same length as 'y'")
      }
      if (anyNA(weights)) {
        stop("'weights' cannot be NA")
      }
      if (any(weights < 0.0)) {
        stop("'weights' must all be non-negative")
      }
      if (!all(is.finite(weights))) {
        stop("'weights' must all be finite")
      }
      # the latent families that carry no weight at all but do carry the mask:
      # a 0/1 vector there is membership, not precision, so it goes to the
      # channel that means it - all-ones included, which the mask normalizes
      # to no mask and so also clears one already installed. The weights slot
      # is left empty, as creation leaves it.
      if (isMaskedWeightFamily(model@family)) {
        if (any(weights != 0 & weights != 1)) {
          stop(
            model@family,
            " models do not support case weights other than 0 and 1, which ",
            "mark rows in and out of the likelihood: such a vector installs ",
            "as the active-row mask, and ",
            if (model@family == "nbinom") {
              "an exposure belongs in the offset as a log-exposure term"
            } else {
              "a weighted truncated-normal latent likelihood is not a coherent model"
            }
          )
        }
        setActiveRows(weights, updateState = updateState)
        return(invisible(NULL))
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

      if (resolveUpdateState(updateState, control)) {
        storeState(ptr)
      }
      invisible(NULL)
    },
    setCounts = function(counts, updateState = NULL) {
      "Replaces a multinomial sampler's response: the n x K matrix of non-negative integer counts whose column k holds category k's successes, with trials n_i = sum_k counts[i, k] at least 0: a row with no trial enters no likelihood and still receives fitted probabilities, and the first such row in a session warns. n and K are fixed at creation - every combiner buffer is sized by n, and K is the forest count - so only the values change. The trees carry over, fitted to the previous counts exactly as setResponse leaves a single-forest sampler's, and the next run forms every category's working response against the new matrix. The matrix is mirrored into data@counts, and its row sums into data@y, so getPointer's transparent re-creation after save/load carries the current response rather than the one the sampler was created with. The sweep draws n_i Polya-Gamma variates per observation per category, so replacing single-trial labels with grouped counts multiplies sweep cost by mean(n_i). updateState follows control@updateState; see setData."
      updateState <- checkUpdateState(updateState)
      requireCountsCapability(.self, "$setCounts")
      ptr <- bartcoreSamplerSetCounts(.self, counts)
      if (resolveUpdateState(updateState, control)) {
        storeState(ptr)
      }
      invisible(NULL)
    },
    setCategoryOffset = function(offset, updateState = NULL) {
      "Installs, or at NULL clears, a multinomial sampler's n x K category offset: the latent becomes f_ik + o_ik, so the offset enters the log-sum-exp margins, every category's working response and the reported softmax probabilities, and never a leaf value. This is the response-side counterpart of setCounts rather than of setOffset, whose flat shift is added after the categories are blended - the wrong side of the nonlinearity - and is the softmax's own null direction besides. Only the row-centred part is identified: adding a constant to a whole row leaves every reported probability unchanged, and the entrance leaves the matrix as given rather than re-centring it. It shifts the TRAIN latent only; the test rows are other rows and carry their own (setCategoryTestOffset), and predict takes its own matrix per call. Mirrored into data@offset.category, so a re-created sampler carries it. updateState follows control@updateState; see setData."
      updateState <- checkUpdateState(updateState)
      requireCountsCapability(.self, "$setCategoryOffset")
      ptr <- bartcoreSamplerSetCategoryOffset(.self, offset)
      if (resolveUpdateState(updateState, control)) {
        storeState(ptr)
      }
      invisible(NULL)
    },
    setCategoryTestOffset = function(offset.test, updateState = NULL) {
      "Installs, or at NULL clears, a multinomial sampler's nTest x K category test offset: the recorded test channel becomes softmax(f_test + o_test), formed where the train blend forms softmax(f + o). The test fits enter no likelihood, so this moves the reported test probabilities and nothing else - no draw, no working response, no train channel. Its rows are the CURRENT test rows, so replacing those rows while it is installed is refused rather than silently reinterpreted; clear it first. Out-of-sample predict does not read it at all, taking its own matrix for the rows it is given. Mirrored into data@offset.category.test, so a re-created sampler carries it. updateState follows control@updateState; see setData."
      updateState <- checkUpdateState(updateState)
      requireCountsCapability(.self, "$setCategoryTestOffset")
      ptr <- bartcoreSamplerSetCategoryTestOffset(.self, offset.test)
      if (resolveUpdateState(updateState, control)) {
        storeState(ptr)
      }
      invisible(NULL)
    },
    setActiveRows = function(active, updateState = NULL) {
      "Sets the per-observation 0/1 mask of rows in the data set for this sampler. An inactive row leaves every sufficient statistic, every family-level parameter update and its own latent draw, but keeps its leaf occupancy and its fitted value. A row switched from inactive to active has its latent redrawn against the current fit before the call returns, from the sampler's own generators and not R's; a call that switches no row back in draws nothing. NULL clears, and an all-ones mask installs nothing. The mask does not ride the saved state; it is mirrored on an R5 field that getPointer, setState and copy reinstall on every re-creation. updateState follows control@updateState; see setData."
      updateState <- checkUpdateState(updateState)
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
      # mirrored only once the engine has taken it, and normalized the way the
      # engine normalizes: an all-ones vector installs nothing, so it records
      # as no mask rather than as a vector to re-apply
      selfEnv <- parent.env(environment())
      selfEnv$activeRows <- if (is.null(active) || all(active == 1)) {
        NULL
      } else {
        active
      }
      if (resolveUpdateState(updateState, control)) {
        storeState(ptr)
      }
      invisible(NULL)
    },
    setForestWeights = function(forest, weights, updateState = NULL) {
      "Sets a per-forest, per-observation weight: a multiplicative precision factor on the named forest's own leaf conditionals, composing with weights and active as (w_i * a_i) * m_f^2 * s_i rather than widening either channel. Only applies to a Bayesian causal forest built with forests = (see dbarts); forest indexes from 1, as with getLeafPrior/getK (the basis forest is 2). The weight does not ride the sampler's saved state; it is mirrored on an R5 field that getPointer and setState both reinstall on every re-creation. updateState follows control@updateState; see setData."
      updateState <- checkUpdateState(updateState)
      refuseCountsMutation(
        .self,
        "$setForestWeights",
        "its forests are its categories, whose margin is a log-sum-exp over ",
        "the other K - 1, so no forest carries a precision of its own"
      )
      weights <- as.double(weights)
      if (length(weights) != length(data@y)) {
        stop("'weights' must have the same length as 'y'")
      }
      # matches the bridge's !R_FINITE(w) || w < 0.0: a forest weight, unlike
      # a case weight, must also be finite, not merely non-negative
      if (!all(is.finite(weights)) || any(weights < 0.0)) {
        stop("forest weights must be finite and non-negative")
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

      if (resolveUpdateState(updateState, control)) {
        storeState(ptr)
      }
      invisible(NULL)
    },
    setForestBasis = function(forest, basis, updateState = NULL) {
      "Changes the basis the named forest's amplitudes multiply, at any forest and any width. forest indexes from 1, as with setForestWeights and getLeafPrior/getK (a Bayesian causal forest's basis forest is 2). basis is a value, or a one-sided formula, which is read as the basis of a forest() is, its names found where the formula was written. A factor expands to its level indicators, one amplitude per level, with no reference level dropped, and may leave a level empty (a swap can leave one momentarily unobserved); a numeric vector or matrix is already those columns, and one of all zeros is refused. Columns are taken by position: the names the basis was created with stay, whatever names the replacement has, and a replacement that has those names in another order is refused. A replacement of another width brings its own names. The forest's label does not change. This is the SOLE route by which a basis changes after creation, and the amplitudes are preserved and remapped: a width-preserving install leaves every one of them bitwise, and a width change carries each forest's block to its new offset and enters the added coordinates at 1. The matrix is mirrored into data@bases as setWeights mirrors weights, so it survives the sampler's re-creation. updateState follows control@updateState; see setData."
      updateState <- checkUpdateState(updateState)
      refuseCountsMutation(
        .self,
        "$setForestBasis",
        "its forests are its categories, which carry no amplitudes and so no ",
        "basis for any to multiply"
      )
      index <- resolveForestIndex(forest)
      if (is.null(basis)) {
        stop("'basis' cannot be NULL")
      }
      if (is.null(data@bases)) {
        stop(
          "$setForestBasis",
          " is not available on a sampler that carries no forest amplitudes: ",
          "amplitudes are fixed at creation; make a new sampler instead"
        )
      }
      if (index >= length(data@bases)) {
        stop("forest index out of range")
      }
      values <- validateForestBases(
        list(replacementForestBasis(basis, length(data@y))),
        length(data@y),
        argument = "basis"
      )[[1L]]
      # columns are taken by position and the names recorded at creation
      # stay; a replacement of another width brings its own
      current <- data@bases[[index + 1L]]
      refuseReorderedBasisNames(values, current, index + 1L)
      if (!is.null(current) && ncol(values) == ncol(current)) {
        dimnames(values) <- if (!is.null(colnames(current))) {
          list(NULL, colnames(current))
        }
      }

      forestInfo <- attr(control, "bartcore.forests", exact = TRUE)
      if (
        NCOL(values) == 1L &&
          !is.null(data@bases[[index + 1L]]) &&
          length(forestInfo$params) > index &&
          identical(forestInfo$params[[index + 1L]][8L], 0)
      ) {
        refuseHeldOneColumn(index + 1L)
      }

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

      if (resolveUpdateState(updateState, control)) {
        storeState(ptr)
      }
      invisible(NULL)
    },
    setSigma = function(sigma, updateState = NULL) {
      "Changes the residual standard deviation parameter for each chain; on a sampler that holds sigma fixed it rewrites the model's fixed value. updateState follows control@updateState; see setData."
      updateState <- checkUpdateState(updateState)
      refuseCountsMutation(
        .self,
        "$setSigma",
        "the softmax carries no residual scale to set"
      )
      sigma <- as.double(sigma)
      if (length(sigma) != 1L) {
        stop("'sigma' must be of length 1")
      }
      if (!is.finite(sigma) || sigma <= 0.0) {
        stop("'sigma' must be finite and positive")
      }

      ptr <- getPointer()
      .Call(C_dbarts_bartcore_setSigma, ptr, sigma)
      # a sigma the sampler holds fixed is model: the write rewrites the
      # fixed value, so a copy or a reload keeps it
      if (is(model@resid.prior, "dbartsFixedPrior")) {
        newModel <- model
        newModel@resid.prior@value <- sigma * sigma
        selfEnv <- parent.env(environment())
        selfEnv$model <- newModel
      }
      if (resolveUpdateState(updateState, control)) {
        storeState(ptr)
      }
      invisible(NULL)
    },
    setPredictor = function(
      x,
      column,
      forceUpdate,
      updateCutPoints = FALSE,
      updateState = NULL
    ) {
      "Changes a single column of the predictor matrix, or the entire matrix if column is missing. updateState follows control@updateState; see setData."
      updateState <- checkUpdateState(updateState)

      checkMissingPolicy(data, sourceAnyNA(x), "predictors")
      result <- withVisible(bartcoreSamplerSetPredictor(
        .self,
        x,
        column = if (missing(column)) NULL else column,
        forceUpdate = if (missing(forceUpdate)) NULL else forceUpdate,
        updateCutPoints = updateCutPoints
      ))
      if (resolveUpdateState(updateState, control)) {
        storeState()
      }
      # a forced update's TRUE comes back invisible and an unforced update's
      # verdict visible; a bare value here would always be visible
      if (result$visible) result$value else invisible(result$value)
    },
    setCutPoints = function(cuts, column, updateState = NULL) {
      "Changes the cut points for the predictors in column, or the entire set itself if the column argument is missing, when the entries of factor columns are not read. A grid out of order is sorted, and a point may appear only once, so for more splits near a value give a denser grid there; the one exception is the grid the column already holds, which is taken as it is, repeats included. Forces the change by pruning any leaves that end up empty. A later setData derives at most n.cuts cut points again, whatever grid was set. updateState follows control@updateState; see setData."
      updateState <- checkUpdateState(updateState)

      bartcoreSamplerSetCutPoints(
        .self,
        cuts,
        column = if (missing(column)) NULL else column
      )
      if (resolveUpdateState(updateState, control)) {
        storeState()
      }
      invisible(NULL)
    },
    setTestPredictor = function(x.test, column) {
      "Changes a single column of the test predictor matrix."

      checkMissingPolicy(data, sourceAnyNA(x.test), "test predictors")
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
        !is.null(x.test) && sourceAnyNA(x.test),
        "test predictors"
      )
      if (missing(offset.test)) {
        # predictors only; the engine keeps the current offset and the
        # bridge refuses if the row count would orphan its length
        return(bartcoreSamplerSetTestPredictor(.self, x.test, column = NULL))
      }

      testRowNames <- observationRowNames(x.test)
      x.test <- validateXTest(x.test, data@x)
      if (is.null(x.test) && !is.null(offset.test)) {
        stop("when test matrix is NULL, test offset must be as well")
      }
      if (!is.null(offset.test)) {
        offset.test <- as.double(offset.test)
        refuseMissingTestOffset(offset.test)
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
      selfEnv$data <- setDataRowNames(data, "test", testRowNames)
      invisible(NULL)
    },
    setTestOffset = function(offset.test) {
      "Changes the test offset."
      ptr <- getPointer()
      selfEnv <- parent.env(environment())

      # refused before the link to the regular offset is broken, so a refused
      # call leaves the sampler as it was
      if (!is.null(offset.test)) {
        refuseMissingTestOffset(offset.test)
      }
      selfEnv$data@testUsesRegularOffset <- FALSE
      if (!is.null(offset.test)) {
        if (is.null(data@x.test)) {
          stop("when test matrix is NULL, test offset must be as well")
        }
        offset.test <- as.double(offset.test)
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
      "Returns the current draw of the augmentation variable, whose meaning is per family and not uniform. A LOCATION, on the sampler's own latent scale, for probit (the truncated normal z), ordinal (the same z under the ordinal thresholds) and aft (the imputed log survival time): these are regressed on directly. A PRECISION, one per observation, for logistic and nbinom (the Polya-Gamma omega) and Student-t (the scale-mixing lambda): these WEIGHT a working response and are not on the response scale at all. NULL for a plain gaussian sampler and for a multinomial one, neither of which augments. Note that a gaussian-family sampler built with family = student() DOES report latents, and they are precisions."
      resultIsMissing <- missing(result)

      ptr <- getPointer()

      .Call(
        C_dbarts_bartcore_getLatents,
        ptr,
        if (resultIsMissing) NULL else result
      )
    },
    getSigmas = function(result) {
      "Returns each chain's current residual standard deviation on the original response scale, or NULL on a heteroscedastic sampler, whose scale is the surface getVariance() reports."

      # the formal is held so it cannot be quietly repurposed: this reader
      # allocates its own vector, and filling a caller's buffer in place is
      # getLatents' contract alone
      if (!missing(result)) {
        stop(
          "'result' is not used by getSigmas: this reader allocates its own vector"
        )
      }

      if (!is.null(attr(control, "bartcore.variance"))) {
        return(NULL)
      }
      ptr <- getPointer()
      .Call(C_dbarts_bartcore_getSigmas, ptr)
    },
    getShape = function() {
      "Returns the shape parameter of the sampler's family currently in force, one per chain, or NULL on a family with none (only nbinom has one today) - the count analog of getSigmas(). It is the same scalar run()$shape records once per kept draw, read mid-sweep and without serializing state, so a host driving the sampler one sweep at a time reads it here; the stored state holds it only where it is drawn. Under a fixed shape it repeats the value the sampler was created with; otherwise it is that sweep's grid draw."
      ptr <- getPointer()
      .Call(C_dbarts_bartcore_getShape, ptr)
    },
    getSumsOfSquaredResiduals = function(result) {
      "Return sum( (y - y.hat)^2 ) on original scale."
      # the formal is held so it cannot be quietly repurposed: this reader
      # allocates its own vector, and filling a caller's buffer in place is
      # getLatents' contract alone
      if (!missing(result)) {
        stop(
          "'result' is not used by getSumsOfSquaredResiduals: this reader ",
          "allocates its own vector"
        )
      }
      ptr <- getPointer()
      .Call(C_dbarts_bartcore_getSumsOfSquaredResiduals, ptr)
    },
    getForestFits = function(forest = NULL) {
      "Returns a sampler's per-forest internal-scale fitted values (a Bayesian causal forest's 1 = prognostic, 2 = treatment; an ordinary sampler's only forest is 1), n.observations x n.chains at one forest, or, at the default forest = NULL, every forest stacked with the forest margin between the observations and the chains, n.observations x n.forests x n.chains (a single-forest sampler's NULL read is bitwise its forest 1 read). forest indexes from 1, as with setForestWeights/setForestBasis/getLeafPrior/getK."
      ptr <- getPointer()
      if (!is.null(forest)) {
        return(.Call(
          C_dbarts_bartcore_getForestFits,
          ptr,
          resolveForestIndex(forest)
        ))
      }
      # the bridge counts forests from 0, as resolveForestIndex converts to
      numForests <- bartcoreNumForests(ptr)
      blocks <- lapply(
        seq_len(numForests),
        function(f) .Call(C_dbarts_bartcore_getForestFits, ptr, f - 1L)
      )
      if (numForests == 1L) {
        return(blocks[[1L]])
      }
      result <- array(
        0.0,
        c(nrow(blocks[[1L]]), numForests, ncol(blocks[[1L]]))
      )
      for (f in seq_len(numForests)) {
        result[, f, ] <- blocks[[f]]
      }
      result
    },
    getFitsWithoutOffset = function() {
      "Returns the sampler's combined per-observation location on the RESPONSE scale and WITHOUT the installed offset, an n.observations x n.chains matrix; run()$train reports the same quantity with the offset folded in, so getFitsWithoutOffset() plus the installed offset is that value. This is the incremental read: getLatents() minus run()$train is biased, because the two are not on the same footing. Refused on a multinomial sampler, whose reported channels are per-category softmax probabilities rather than one additive location. Contrast getForestFits, which reports ONE forest's INTERNAL-scale totals."
      refuseCountsMutation(
        .self,
        "$getFitsWithoutOffset",
        "its reported channels are per-category softmax probabilities rather ",
        "than one additive location; $predict(data@x) serves that read"
      )
      ptr <- getPointer()
      .Call(C_dbarts_bartcore_getFitsWithoutOffset, ptr)
    },
    getVariance = function(test = FALSE) {
      "Returns a heteroscedastic sampler's current variance surface s^2(x) on the ORIGINAL response scale, an n.observations x n.chains matrix at the default test = FALSE and an n.test x n.chains matrix at test = TRUE. This is the mid-sweep read of the channels run() records as 'variance' and 'varianceTest': at the state a recorded sweep left, the two agree exactly. It reports a VARIANCE, so a residual scale is its square root, and it is the surface analog of getSigmas(), which reports the scalar sigma a homoscedastic sampler carries. NULL where those channels report nothing: on a homoscedastic sampler, and at test = TRUE with no test rows installed. Unlike predict(), it reads the trees currently in force, so it answers after a prior draw and needs no keepTrees."
      ptr <- getPointer()
      .Call(C_dbarts_bartcore_getVariance, ptr, isTRUE(test))
    },
    getForestAmplitudes = function(forest = NULL) {
      "Returns the named forest's amplitudes - the scalars its basis columns are multiplied by, one per column - as a q x n.chains matrix, or, at the default forest = NULL, every forest's stacked forest-major into a sum(q) x n.chains matrix, which is the row order the run's own glue channel carries. The vector is RAGGED, forest by forest, which is why a forest can be named: a Bayesian causal forest's forest 1 carries the single a on its implicit intercept and its forest 2 the pair (b0, b1) on its two level indicators, so the stacked read is its shipped (a, b0, b1). forest indexes from 1, as with setForestBasis/setForestWeights/getLeafPrior."
      ptr <- getPointer()
      .Call(
        C_dbarts_bartcore_getForestAmplitudes,
        ptr,
        if (is.null(forest)) NULL else resolveForestIndex(forest)
      )
    },
    getForestVariableCounts = function(forest = NULL) {
      "Returns a sampler's per-forest predictor split counts (a Bayesian causal forest's 1 = prognostic, 2 = treatment; an ordinary sampler's only forest is 1), n.predictors x n.chains at one forest, or, at the default forest = NULL, every forest stacked with the forest margin between the predictors and the chains, n.predictors x n.forests x n.chains (a single-forest sampler's NULL read is bitwise its forest 1 read). Rows are named by the predictor columns when data@x carries colnames, on margin 1 in both shapes. forest indexes from 1, as with setForestWeights/setForestBasis/getLeafPrior/getK."
      ptr <- getPointer()
      if (is.null(forest)) {
        numForests <- bartcoreNumForests(ptr)
        blocks <- lapply(
          seq_len(numForests),
          function(f) {
            .Call(C_dbarts_bartcore_getForestVariableCounts, ptr, f - 1L)
          }
        )
        counts <- if (numForests == 1L) {
          blocks[[1L]]
        } else {
          result <- array(
            0L,
            c(nrow(blocks[[1L]]), numForests, ncol(blocks[[1L]]))
          )
          for (f in seq_len(numForests)) {
            result[, f, ] <- blocks[[f]]
          }
          result
        }
      } else {
        counts <- .Call(
          C_dbarts_bartcore_getForestVariableCounts,
          ptr,
          resolveForestIndex(forest)
        )
      }
      predictorNames <- colnames(data@x)
      if (!is.null(predictorNames)) {
        rownames(counts) <- predictorNames
      }
      counts
    },
    getLeafPrior = function(forest = NULL) {
      "Returns the leaf prior a forest runs under, alone, as a named list: leaf.prior, the specification in the terms it was named in - normal(), linear() or gp() carrying one of k (a number or a chi() law) or sd (a number or an invchi() law), the family default when none was named - which goes back into setLeafPrior or a fitting function's leaf.prior as is; leaf.model; prior.sd.of, what the sd is the sd of ('leaf value', 'coefficient' or 'amplitude'); prior.mean; k.scale, the value k is relative to, so the spread in force on each chain is k.scale / getK() - the data's scale under a k-named prior and under sd = invchi(df, 0), and otherwise, under an sd-named prior, twice the sd or invchi() scale in force; response.scale and response.shift. On a forest whose scale a multi-forest calibration map sets, k is pinned at 1, leaf.prior is the forest(sd = ) creation takes, which goes back into setLeafPrior(forests = ) - the half-Cauchy median on a forest created without a basis (prior.sd.of 'amplitude scale'), the leaf-scale factor otherwise ('forest total') - and the list adds basis.row.norm, leaf.scale.factor and leaf.scale.divisor, and one of amplitude.prior.variance or amplitude.prior.scale; they are absent elsewhere. Every chain runs under the sampler's one prior and response transform, which no state install moves, so every value is shared by the chains; an NA spread is refused on write. A drawn k is chain state, read by getK. At the default forest = NULL a multi-forest sampler returns an unnamed list of one prior per forest; a single-forest sampler's NULL read is bitwise its forest 1 read."
      ptr <- getPointer()
      read <- function(index) {
        reportLeafPrior(
          .self,
          .Call(C_dbarts_bartcore_getLeafPrior, ptr, index)
        )
      }
      if (!is.null(forest)) {
        return(read(resolveForestIndex(forest)))
      }
      numForests <- bartcoreNumForests(ptr)
      if (numForests == 1L) {
        return(read(0L))
      }
      lapply(seq_len(numForests) - 1L, read)
    },
    getK = function(forest = NULL) {
      "Returns each chain's current leaf-prior k, the value run()$k records per draw, read without running, as getSigmas reports sigma; after a run it is bitwise the last draw. A fixed k repeats per chain, and a forest whose scale a multi-forest calibration map sets reports 1. It is k whatever terms the prior was named in, relative to getLeafPrior()$k.scale. A vector of length n.chains at one forest, or, at the default forest = NULL on a multi-forest sampler, an n.forests x n.chains matrix; a single-forest sampler's NULL read is bitwise its forest 1 read."
      ptr <- getPointer()
      read <- function(index) {
        .Call(C_dbarts_bartcore_getLeafPrior, ptr, index)[, "k"]
      }
      if (!is.null(forest)) {
        return(read(resolveForestIndex(forest)))
      }
      numForests <- bartcoreNumForests(ptr)
      if (numForests == 1L) {
        return(read(0L))
      }
      do.call(rbind, lapply(seq_len(numForests) - 1L, read))
    },
    setLeafPrior = function(leaf.prior, forests = NULL, updateState = NULL) {
      "Restates the leaf prior's spread, or the hyperprior it is drawn under, on every chain, in the vocabulary a fitting function's leaf.prior takes: normal(sd = ), normal(k = ), an invchi() law on the sd, linear(sd = ) or gp(sd = ). The specification must name the sampler's own leaf model; leaf-model details such as a linear leaf's columns may be omitted and, if given, must match. Nothing else moves - not the tree prior, the response transform or sigma. Under a drawn k the engine keeps its current k across the write, so a change of k.scale - between the k and sd forms, or of an invchi() scale - scales the next sweep's spread by new k.scale / old k.scale, and getK and the spread in force jump with it until the law pulls k back. A multinomial sampler takes normal(k = ) with a fixed k, Inf included, applied to every category forest. A sampler whose forests carry amplitudes takes forests = list(forest(sd = ), ...) instead of leaf.prior, as its creation does: the same positions and names, a short list reaching the first forests, and a forest whose sd is not stated left as it is; normal() and normal(k = 2), which its creation also accepts, change nothing. Give exactly one of leaf.prior and forests. The write takes effect on the next sweep, reinterpreting no value already drawn; a write equal to what is in force is bitwise inert. The write is recorded on the model field, or for forests on the control, so a re-creation or a later re-anchoring channel restates it rather than the creation value. setModel changes everything else. updateState follows control@updateState; see setData."
      # a forest = index would otherwise match forests = partially
      if ("forest" %in% names(sys.call())) {
        stop(
          "$setLeafPrior takes no 'forest' index; a sampler whose forests ",
          "carry amplitudes restates them as forests = list(forest(sd = ), ...)"
        )
      }
      updateState <- checkUpdateState(updateState)
      multinomial <- samplerCarriesCounts(.self)
      amplitudes <- samplerCarriesAmplitudes(.self)
      if (!missing(forests)) {
        forests <- evalInForestVocabulary(
          substitute(forests),
          forestConstructors[FOREST_ARGUMENT_VOCABULARIES$forests],
          parent.frame()
        )
      }
      if (!is.null(forests)) {
        if (multinomial) {
          stop(
            multinomialLeafPriorMessage,
            "its forests are its categories; normal(k = ) states every one",
            call. = FALSE
          )
        }
        if (!amplitudes) {
          stop(
            "$setLeafPrior's 'forests' restates a forest's 'sd', which only a ",
            "sampler whose forests carry amplitudes has; state this sampler's ",
            "leaf prior as leaf.prior"
          )
        }
        if (!missing(leaf.prior)) {
          stop("give $setLeafPrior either 'leaf.prior' or 'forests', not both")
        }
        sds <- resolveForestSpreads(.self, forests)
        ptr <- getPointer()
        writeForestSpreads(.self, ptr, sds)
        if (resolveUpdateState(updateState, control)) {
          storeState(ptr)
        }
        return(invisible(NULL))
      }
      if (missing(leaf.prior)) {
        stop(
          "'leaf.prior' must be given: a leaf prior specification such as ",
          "normal(sd = 1)",
          if (amplitudes) ", or forests = list(forest(sd = ), ...)"
        )
      }
      expr <- substitute(leaf.prior)
      spec <- evalInVocabulary(
        expr,
        setLeafPriorVocabulary(.self),
        parent.frame(),
        resolvedAs(
          "leaf.prior",
          "dbartsLeafPrior",
          "leaf prior specification"
        )
      )
      # after the argument checks above, so a malformed call is answered on its
      # own terms rather than by the refusal that would follow a well-formed one
      if (amplitudes) {
        inert <- is(spec, "dbartsNormalPrior") &&
          is.null(spec@prior.sd) &&
          (is.null(spec@k) || identical(spec@k, 2) || identical(spec@k, 2L))
        if (!inert) {
          stop(
            "$setLeafPrior on a sampler whose forests carry amplitudes takes ",
            "forests = : the multi-forest calibration map sets every ",
            "forest's leaf scale; state a forest's spread as at creation, ",
            "forests = list(forest(sd = ), ...)",
            call. = FALSE
          )
        }
        if (resolveUpdateState(updateState, control)) {
          storeState()
        }
        return(invisible(NULL))
      }
      if (multinomial) {
        refuseMultinomialLeafPrior(spec)
      }
      newModel <- restateLeafPrior(.self, spec, expr)
      ptr <- getPointer()
      if (multinomial) {
        .Call(
          C_dbarts_bartcore_setForestK,
          ptr,
          as.double(newModel@leaf.hyperprior@k)
        )
      } else {
        writeLeafPrior(.self, ptr, newModel)
      }
      selfEnv <- parent.env(environment())
      selfEnv$model <- newModel
      if (resolveUpdateState(updateState, control)) {
        storeState(ptr)
      }
      invisible(NULL)
    },
    reapplyForestWeights = function(ptr) {
      "Re-installs every forest weight mirrored on this sampler onto ptr, a freshly (re-)created or restated engine pointer that carries none of them. Called from getPointer, setState and copy, never recursing through getPointer."
      for (forest in seq_along(forestWeights)) {
        weights <- forestWeights[[forest]]
        if (!is.null(weights)) {
          .Call(C_dbarts_bartcore_setForestWeights, ptr, forest - 1L, weights)
        }
      }
      invisible(NULL)
    },
    reapplyActiveRows = function(ptr) {
      "Re-installs the active-row mask mirrored on this sampler onto ptr, a freshly (re-)created or restated engine pointer that carries none. Called from getPointer, setState and copy, never recursing through getPointer."
      if (!is.null(activeRows)) {
        .Call(C_dbarts_bartcore_setActiveRows, ptr, activeRows)
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
        refuseLegacyState(state)
        ptr <- recreatePointer(control, model, data, TRUE)
        # a same-spec continuation skips re-quantization; data@x serves any
        # cross-grid column (the engine keeps no predictor matrix)
        # a store sized through the flat API is in no control, so the
        # re-created sampler takes the stored state's capacity
        .Call(
          C_dbarts_bartcore_setState,
          ptr,
          state,
          rawPredictorMatrix(data@x),
          TRUE
        )
        reapplyForestWeights(ptr)
        reapplyActiveRows(ptr)
        reissueNamedLeafSd(.self, ptr)
        # the replacement is bound only once it carries the state: a refused
        # install must leave the object exactly as it was rather than holding
        # a live but unfitted engine that the next run would silently sample
        # from and then store over the fitted state. The abandoned pointer
        # carries its own holder and finalizes itself.
        selfEnv$pointer <- ptr
      }
      pointer
    },
    setState = function(newState) {
      "Installs a stored state: the chains, never the model. A state in other response units is converted into the sampler's; a Gaussian-process leaf or forests with amplitudes refuse one under another response shift, saved draws or not. Invisibly returns TRUE when nothing had to be changed to install it, FALSE otherwise. See Saving and Value in ?dbartsSampler."
      refuseLegacyState(newState)
      if (!inherits(newState, "bartcoreState")) {
        stop("'state' must inherit from bartcoreState")
      }
      selfEnv <- parent.env(environment())
      ptr <- pointer
      if (.Call(C_dbarts_bartcore_isValidPointer, pointer) == FALSE) {
        ptr <- recreatePointer(control, model, data, TRUE)
      }
      exact <- .Call(
        C_dbarts_bartcore_setState,
        ptr,
        newState,
        rawPredictorMatrix(data@x),
        FALSE
      )
      reapplyForestWeights(ptr)
      reapplyActiveRows(ptr)
      reissueNamedLeafSd(.self, ptr)
      # as in getPointer: a re-created engine is bound only after the install
      # succeeds, so a refusal leaves a dead pointer dead instead of live and
      # unfitted, and leaves 'state' the one that is still installed
      selfEnv$pointer <- ptr
      selfEnv$state <- newState
      invisible(exact)
    },
    startThreads = function(n.threads = control@n.threads) {
      "Retired: threads are owned by each run. Does nothing."
      noOpThreadMethod("startThreads")
    },
    stopThreads = function() {
      "Retired: threads are owned by each run. Does nothing."
      noOpThreadMethod("stopThreads")
    },
    storeState = function(ptr = getPointer()) {
      "Updates the cached internal state used for saving/loading: the chains and the units they are stored in, no prior or fixed value."
      selfEnv <- parent.env(environment())
      selfEnv$state <- .Call(C_dbarts_bartcore_storeState, ptr)
      invisible(NULL)
    },
    installTrees = function(donor, samples = NULL) {
      "Warm-starts the forests from a donor sampler or bart fit over the same
       predictors, keeping this sampler's model; a donor in other response
       units is converted into this sampler's. 'samples' maps each chain to a
       1-based donor-sample index; NULL spreads the chains across the donor's
       kept samples. Single-forest samplers only."
      ptr <- getPointer()
      refuseMultiForestWarmStart(ptr, "$installTrees")
      donorState <- warmStartState(donor)
      if (!is.null(samples)) {
        samples <- coerceOrError(samples, "integer")
      }
      .Call(C_dbarts_bartcore_installForests, ptr, donorState, samples)
      reissueNamedLeafSd(.self, ptr)
      storeState(ptr)
      invisible(NULL)
    },
    printTrees = function(treeNums, chainNums, sampleNums) {
      "Produces an info dump of the internal state of the trees."
      matchedCall <- match.call()
      if (is.null(matchedCall$chainNums)) {
        chainNums <- seq_len(control@n.chains)
      }
      # NULL asks the engine for every RECORDED draw: a store still filling
      # holds fewer than n.samples, and only the engine knows how many
      if (is.null(matchedCall$sampleNums)) {
        sampleNums <- NULL
      } else {
        if (!control@keepTrees) {
          warning(
            "sampleNums ignored if keepTrees is FALSE",
            call. = FALSE
          )
          sampleNums <- NULL
        } else {
          sampleNums <- coerceOrError(sampleNums, "integer")
        }
      }
      if (is.null(matchedCall$treeNums)) {
        treeNums <- seq_len(control@n.trees)
      }

      ptr <- getPointer()
      invisible(.Call(
        C_dbarts_bartcore_printTrees,
        ptr,
        coerceOrError(chainNums, "integer"),
        sampleNums,
        coerceOrError(treeNums, "integer")
      ))
    },
    getTrees = function(
      treeNums,
      chainNums,
      sampleNums,
      current = FALSE,
      newdata = NULL,
      forest = NULL
    ) {
      "Returns a data.frame containing the internal state of the trees, one row per node. A sampler with several forests (a multinomial one, or one declared with forests =, whatever the count) puts a leading 'forest' column (indexed from 1) on it, and stacks every forest forest-major at the default forest = NULL, as the sampler's other per-forest readers stack at their own default; any other sampler returns no forest column, and accepts forest = 1. forest also takes a single index or a vector of them, each validated as getLeafPrior/getForestFits/getForestAmplitudes/getForestVariableCounts validate one. treeNums defaults to, and is validated against, EACH selected forest's own tree count, which need not match forest 1's."
      matchedCall <- match.call()
      current <- isTRUE(current)
      # live working trees have no sample dimension, so treat a current request
      # like a non-keepTrees sampler for sample handling
      useSaved <- control@keepTrees && !current
      if (is.null(matchedCall$chainNums)) {
        chainNums <- seq_len(control@n.chains)
      }
      # as for printTrees: NULL is every recorded draw, which only the engine
      # counts
      if (is.null(matchedCall$sampleNums)) {
        sampleNums <- NULL
      } else {
        if (!useSaved) {
          warning(
            if (current) {
              "sampleNums ignored if current is TRUE"
            } else {
              "sampleNums ignored if keepTrees is FALSE"
            },
            call. = FALSE
          )
          sampleNums <- NULL
        } else {
          sampleNums <- coerceOrError(sampleNums, "integer")
        }
      }
      # a later forest can carry its own n.trees (a forest() term or a
      # forests = entry), so a supplied treeNums is checked per forest below
      # rather than against control@n.trees (forest 1's count alone); left
      # unsupplied, it defaults to EACH forest's own seq_len
      treeNumsSupplied <- !is.null(matchedCall$treeNums)
      if (treeNumsSupplied) {
        treeNums <- coerceOrError(treeNums, "integer")
      }

      chainNums <- coerceOrError(chainNums, "integer")

      if (anyNA(chainNums)) {
        stop("'chainNums' contains missing values")
      }
      if (useSaved && anyNA(sampleNums)) {
        stop("'sampleNums' contains missing values")
      }
      if (treeNumsSupplied && anyNA(treeNums)) {
        stop("'treeNums' contains missing values")
      }
      if (any(chainNums <= 0 | chainNums > control@n.chains)) {
        stop("'chainNums' must be in [1, ", control@n.chains, "]")
      }
      if (
        useSaved &&
          any(sampleNums <= 0 | sampleNums > control@n.samples)
      ) {
        stop("'sampleNums' must be in [1, ", control@n.samples, "]")
      }

      # route new data through the trees so 'n' counts that data instead of the
      # training predictors; validated and coded as for predict, and routed off
      # whatever storage it arrives in
      if (!is.null(newdata)) {
        newdata <- validateXTest(newdata, data@x)
      }

      ptr <- getPointer()
      # NULL is every forest, forest-major, as getForestAmplitudes stacks.
      # The forest column appears only where the sampler has several forests
      # or was declared with forests =, so a single-forest sampler's table is
      # the one 0.9-x returned. The bound is
      # checked here, ahead of the .Call, so an out-of-range forest reads the
      # same "forest index out of range" the sibling readers raise rather
      # than getTrees' own bridge-side wording; an empty forest vector is
      # refused the same way a scalar one's length check would refuse it.
      forestIndices <- if (is.null(forest)) {
        seq_len(bartcoreNumForests(ptr)) - 1L
      } else if (length(forest) == 0L) {
        resolveForestIndex(forest)
      } else {
        indices <- vapply(forest, resolveForestIndex, 0L)
        if (any(indices >= bartcoreNumForests(ptr))) {
          stop("forest index out of range")
        }
        indices
      }
      # saved-tree replay reads the current training predictors (the engine
      # keeps no matrix); a sparse data@x is skipped for a NULL replay source
      trainingMatrix <- rawPredictorMatrix(data@x)
      hasForestColumn <- bartcoreNumForests(ptr) > 1L ||
        isTRUE(attr(control, "bartcore.forestsDeclared"))
      blocks <- lapply(forestIndices, function(forestIndex) {
        forestTreeCount <- bartcoreForestTreeCount(ptr, forestIndex)
        forestTreeNums <- if (treeNumsSupplied) {
          treeNums
        } else {
          seq_len(forestTreeCount)
        }
        if (any(forestTreeNums <= 0 | forestTreeNums > forestTreeCount)) {
          stop(
            "'treeNums' must be in [1, ",
            forestTreeCount,
            "] for forest ",
            forestIndex + 1L
          )
        }
        block <- .Call(
          C_dbarts_bartcore_getTrees,
          ptr,
          chainNums,
          sampleNums,
          forestTreeNums,
          current,
          newdata,
          trainingMatrix,
          forestIndex
        )
        # cbind's recycling refuses a length-1 scalar against a zero-row
        # block (an empty treeNums/sampleNums/chainNums selection), so the
        # forest column is sized explicitly rather than recycled
        if (hasForestColumn) {
          cbind(forest = rep(forestIndex + 1L, nrow(block)), block)
        } else {
          block
        }
      })
      trees <- if (length(blocks) == 1L) {
        blocks[[1L]]
      } else {
        do.call(rbind, blocks)
      }
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
      if (is(model@leaf.prior, "dbartsLinearPrior")) {
        covariateNames <- colnames(data@x)[model@leaf.prior@columns]
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
      forest = NULL,
      treePlotPars = c(nodeHeight = 12, nodeWidth = 40, nodeGap = 8),
      ...
    ) {
      "Minimialist visualization of tree branching and contents. forest, as with getTrees, defaults to the sampler's only forest and is required on a sampler with more than one."

      refusePlotTreeArgs(sys.call())
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
      if (is.null(forest)) {
        forest <- if (bartcoreNumForests(getPointer()) == 1L) {
          1L
        } else {
          stop("forest required if more than one forest in sampler")
        }
      } else {
        # a single tree's rows are ambiguous once more than one forest
        # contributes them; resolveForestIndex enforces a single positive
        # integer, its 0-based return unused since getTrees resolves forest
        # again on its own terms
        resolveForestIndex(forest)
      }

      tree <-
        if (control@keepTrees) {
          .self$getTrees(treeNum, chainNum, sampleNum, forest = forest)
        } else {
          .self$getTrees(treeNum, chainNum, forest = forest)
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
      oldpar <- par(no.readonly = TRUE)
      on.exit(par(oldpar), add = TRUE)
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
