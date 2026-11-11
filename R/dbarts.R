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
  placeOnHazardGrid(time, periods)
}

# Each subject's terminal period on a grid already made, at creation or a
# sampler's own: the period whose interval holds its time, a time past the last
# period in the last (expandDiscreteTimeHazard then leaves it censored there, at
# risk through every period with no event), and a subject with no time in the
# first, as every subject is at risk there; it keeps that one row, with a
# missing response.
placeOnHazardGrid <- function(time, periods) {
  terminal <- pmin(
    findInterval(time, periods, left.open = TRUE) + 1L,
    length(periods)
  )
  terminal[is.na(terminal)] <- 1L
  list(periods = periods, terminalPeriod = terminal)
}

# The times a hazard grid is read from: those of the subjects in the
# likelihood (weight positive). With none in, every subject's, under one
# warning. 'time' and 'weights' are per subject, the subjects the fit keeps.
hazardGridTimes <- function(time, weights) {
  if (is.null(weights)) {
    return(time)
  }
  if (length(weights) == 1L) {
    weights <- rep_len(weights, length(time))
  }
  subjectsIn <- rowsInLikelihood(weights)
  if (!any(subjectsIn)) {
    warnClassed(
      "dbartsScaleFallbackWarning",
      "no subject is in the likelihood (",
      countInLikelihood(0L, length(time), "subject"),
      "), so the hazard periods are taken from ",
      takenFromAll(length(time), "subject")
    )
    return(time)
  }
  time[subjectsIn]
}

# The person-period expander: a subject
# i observed to time_i (event or censoring) becomes its at-risk rows, one per
# period k = 1..t_i, each carrying x_i, the ordinal period column (appended
# LAST), and the binary indicator y_ik = status_i * 1{k = t_i}. A censored
# subject's rows are all zero (right-censoring is pure data shape). Offsets and
# weights replicate per subject and follow the chosen binary family's policy
# downstream. Returns the expanded design, the binary response, the period
# grid (for the $periods marker), and the replicated offset/weights. Once the
# grid gives the row count, and before any row is built, the fit's memory is
# projected and warned of (warnHazardMemory) from 'run' (hazardRunSettings),
# counting the 'n.test' held-out subjects at every period.
expandDiscreteTimeHazard <- function(
  x,
  time,
  status,
  run,
  breaks = NULL,
  offset = NULL,
  weights = NULL,
  gridTime = time,
  periods = NULL,
  n.test = 0L
) {
  n <- length(time)
  grid <- if (is.null(periods)) {
    resolveHazardGrid(time, breaks, gridTime)
  } else {
    placeOnHazardGrid(time, periods)
  }
  periods <- grid$periods
  terminal <- grid$terminalPeriod

  # in double precision, so a row count past the integer range is projected
  # rather than overflowing
  warnHazardMemory(
    sum(as.double(terminal)),
    as.double(n.test) * length(periods),
    NCOL(x),
    run
  )

  # subject-major, period ascending within subject: the period-1 rows are the
  # subjects in order (survivalProbabilities reconstructs training covariates
  # from them), and sequence() supplies the within-subject period counter
  subjectOf <- rep.int(seq_len(n), terminal)
  periodOf <- sequence(terminal)
  y <- as.double(periodOf == terminal[subjectOf] & status[subjectOf] == 1.0)
  # a subject past the last period (one out of the likelihood, whose time the
  # grid was not read from) is at risk in every period and has no event in it
  if (length(periods) > 0L) {
    y[which(time[subjectOf] > periods[length(periods)])] <- 0
  }
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

# The settings of the run a hazard fit's memory is projected for: its chains,
# its kept draws a chain, whether it keeps the training rows' fits and any
# row's, and 'n.trees', its trees a chain. 'setData' says the rows are handed
# to a sampler's $setData, whose periods are its own.
hazardRunSettings <- function(control, n.trees, setData = FALSE) {
  list(
    n.chains = control@n.chains,
    n.samples = control@n.samples,
    n.trees = n.trees,
    keepTrainingFits = control@keepTrainingFits,
    keepFits = control@keepFits,
    setData = setData
  )
}

# The peak memory, in bytes, of a hazard fit over 'rows' person-period
# training rows and 'testRows' held-out ones, with 'numPredictors' predictor
# columns besides the period, run as 'run' (hazardRunSettings) says. Measured
# on one build: a training row takes about 200 bytes and 34 a predictor
# column, which the chains share, and each chain about 20 bytes and 8 a tree,
# the fits of the trees it holds; a held-out row about 30 bytes and 36 a
# predictor column. A kept draw of a row's fit is 8 bytes a chain in the fit
# returned and about 16 at the peak, while it is copied out: the training
# rows' are kept under keepTrainingFits, and no row's without keepFits. Every
# figure counts doubles, integers and pointers, so it holds on any 64-bit
# machine.
hazardFitBytes <- function(rows, testRows, numPredictors, run) {
  keptDraws <- if (run$keepFits) run$n.chains * run$n.samples else 0
  trainingDraws <- if (run$keepTrainingFits) keptDraws else 0
  rows *
    (200 +
      34 * numPredictors +
      run$n.chains * (20 + 8 * run$n.trees) +
      16 * trainingDraws) +
    testRows * (30 + 36 * numPredictors + 16 * keptDraws)
}

# The machine's physical memory, in bytes, NA where the operating system does
# not report it. An unexported, undocumented option, dbarts.physicalMemory,
# read fresh per call, stands in for it, so a test can put a projection over
# it or make it unreadable with NA.
physicalMemory <- function() {
  forced <- getOption("dbarts.physicalMemory")
  if (!is.null(forced)) {
    return(as.double(forced)[1L])
  }
  .Call(C_dbarts_physicalMemory)
}

# A byte count in words, in megabytes or gigabytes of 2^20 and 2^30 bytes.
formatBytes <- function(bytes) {
  if (bytes < 2^30) {
    paste(format(signif(bytes / 2^20, 2L)), "MB")
  } else {
    paste(format(signif(bytes / 2^30, 2L)), "GB")
  }
}

# Warns, once per session, when a hazard fit's projected peak memory
# (hazardFitBytes) exceeds the machine's physical memory, naming the rows, the
# two figures and what lowers the first; silent where that memory cannot be
# read. Called before the rows are built, and the fit runs whatever it says.
warnHazardMemory <- function(rows, testRows, numPredictors, run) {
  available <- physicalMemory()
  if (is.na(available)) {
    return(invisible(NULL))
  }
  needed <- hazardFitBytes(rows, testRows, numPredictors, run)
  if (needed <= available) {
    return(invisible(NULL))
  }
  keepsTraining <- run$keepFits && run$keepTrainingFits
  keptDraws <- if (run$keepFits && (keepsTraining || testRows > 0)) {
    paste0(
      " keeping ",
      run$n.samples,
      " draws a chain of the ",
      if (keepsTraining) "training" else "test",
      " fits"
    )
  }
  # a sampler's grid and chain count are fixed at creation, so at $setData
  # those two levers take a new sampler
  fewerChains <- run$n.chains > 1L
  ways <- if (run$setData) {
    c(
      "hand over fewer subjects",
      if (keepsTraining) "set keepTrainingFits = FALSE",
      paste0(
        "make a new sampler with hazard(breaks = )",
        if (fewerChains) " or fewer chains"
      )
    )
  } else {
    c(
      "coarsen the grid with hazard(breaks = )",
      if (keepsTraining) "set keepTrainingFits = FALSE",
      if (fewerChains) "run fewer chains"
    )
  }
  numWays <- length(ways)
  if (numWays > 1L) {
    ways <- paste0(
      paste(ways[-numWays], collapse = ", "),
      if (numWays > 2L) ",",
      " or ",
      ways[numWays]
    )
  }
  warnOnce(
    "hazardMemory",
    class = "dbartsMemoryWarning",
    if (run$setData) "$setData: " else "",
    "the hazard fit expands to ",
    format(rows, scientific = FALSE),
    " person-period rows",
    if (testRows > 0) {
      paste0(" and ", format(testRows, scientific = FALSE), " test rows")
    },
    if (run$setData) " on the sampler's periods",
    "; at ",
    run$n.chains,
    if (run$n.chains == 1L) " chain" else " chains",
    keptDraws,
    " it needs about ",
    formatBytes(needed),
    " at its peak, and this machine has ",
    formatBytes(available),
    "; ",
    ways,
    ". This warning does not stop the fit and is shown once per session"
  )
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

# A hazard fit's held-out subjects at every one of the K periods, period-major,
# with their offsets and weights repeated alike. A held-out subject has no
# event time to place it by.
expandHazardTest <- function(data, K) {
  if (is.null(data@x.test)) {
    return(data)
  }
  n.test <- nrow(data@x.test)
  data@x.test <- appendHazardPeriodColumn(
    hazardRowSubset(data@x.test, rep(seq_len(n.test), times = K)),
    rep(seq_len(K), each = n.test)
  )
  if (!is.null(data@offset.test)) {
    offset <- data@offset.test
    if (length(offset) == 1L) {
      offset <- rep_len(offset, n.test)
    }
    data@offset.test <- rep(offset, times = K)
  }
  if (!is.null(data@weights.test)) {
    weights <- data@weights.test
    if (length(weights) == 1L) {
      weights <- rep_len(weights, n.test)
    }
    data@weights.test <- rep(weights, times = K)
  }
  data
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
## this default is guessNumCores() capped to n.chains, or one where that is NA
## (dec-B115 -
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
  useQuantiles = NULL,
  treeShift = c("auto", "always", "never"),
  keepTrees = FALSE,
  storage = c("double", "single"),
  n.samples = NULL,
  n.cuts = NULL,
  n.burn = 200L,
  n.trees = NULL,
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
  n.perturb.cuts = 1L,
  seed = NULL,
  updateState = TRUE,
  ...
) {
  # the names THIS call carried, which the front doors' precedence rule reads
  # off the object: a slot set to the same value the constructor would have
  # chosen is otherwise indistinguishable from an untouched one
  namedHere <- names(match.call())[-1L]
  # these three default to NULL, so a written NULL is the name left out
  namedHere <- setdiff(
    namedHere,
    names(Filter(
      is.null,
      list(
        n.trees = n.trees,
        n.cuts = n.cuts,
        useQuantiles = useQuantiles
      )
    ))
  )
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
  unstatedThreads <- missing(n.threads)
  n.threads <- fallBackToOneThread(n.threads, unstatedThreads)
  refuseNullArguments(proposal.probs)
  copies <- checkFlatCopies(n.trees, n.cuts, useQuantiles)
  result <- newValidated(
    "dbartsControl",
    verbose = as.logical(verbose),
    keepTrainingFits = as.logical(keepTrainingFits),
    keepFits = as.logical(keepFits),
    useQuantiles = if (is.null(useQuantiles)) {
      freshControlValue("useQuantiles")
    } else {
      copies$useQuantiles
    },
    levelGibbs = levelGibbs,
    keepTrees = as.logical(keepTrees),
    storage = storage,
    n.samples = coerceOrError(n.samples, "integer"),
    n.burn = coerceOrError(n.burn, "integer"),
    n.trees = if (is.null(n.trees)) {
      freshControlValue("n.trees")
    } else {
      copies$n.trees
    },
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
    n.perturb.cuts = coerceOrError(n.perturb.cuts, "numeric"),
    seed = resolveSeedArg(seed, "dbartsControl", refuse = TRUE),
    updateState = as.logical(updateState)
  )
  # a plain attribute, deliberately not a slot (a slot would enter every
  # validity check and every stored control's printed form)
  attr(result, controlSuppliedAttr) <- unique(c(
    intersect(namedHere, names(formals(dbarts::dbartsControl))),
    # the retired spelling names the same slot
    if ("rngSeed" %in% supplied) "seed"
  ))
  # the tombstoned cut count rides the attribute 0.9-34 kept it on, so
  # ctl@n.cuts reads what the caller wrote
  attr(result, controlCutsAttr) <- copies$n.cuts
  warnWrittenCopies(intersect(namedHere, controlWrittenNames))
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
## construction edit (ctl@n.burn <- 200L) still speaks. The tree count, the cut
## count and the quantile rule are not among these: readWrittenCopies has them.
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

## The class check of a control given to bart() and xbart(), which take any.
requireDbartsControl <- function(control, caller) {
  if (!inherits(control, "dbartsControl")) {
    stop(
      "'control' argument to ",
      caller,
      " must be of class dbartsControl; use dbartsControl() to create one",
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
    if (
      !(is.logical(verbose) || is.numeric(verbose)) ||
        length(verbose) != 1L ||
        is.na(verbose)
    ) {
      stop("'verbose' argument to ", funcName, " must be TRUE or FALSE")
    }
    envir$verbose <- as.logical(verbose)
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
  # entry point; the estimate the entry point holds is then NA_real_. So is
  # it for the name of a rule and for xbart's function, which the entry
  # point has read
  if (
    !missing(sigest) &&
      !is.null(sigest) &&
      !isSingleNA(sigest) &&
      !isSigestRule(sigest) &&
      !isSigestFunction(sigest)
  ) {
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
  power = 2,
  base = 0.95,
  split.probs = NULL,
  leaf.prior = normal,
  monotone = NULL,
  interactions = NULL,
  blocks = NULL,
  variance = NULL,
  forests = NULL,
  n.trees = NULL,
  control = dbarts::dbartsControl(),
  sigest = NULL,
  seed = NULL,
  factors = c("categorical", "indicators"),
  n.cuts = NULL,
  useQuantiles = NULL,
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
  weights.type = NULL,
  na.action = dbarts::na.keepPredictors,
  sigma = NULL,
  node.prior = NULL,
  tree.prior = NULL,
  callback = NULL,
  ...
) {
  matchedCall <- match.call()
  # a data object as handed over, which a forest's basis is read against
  dataObject <- NULL
  if (!missing(formula) && inherits(formula, "dbartsData")) {
    dataObject <- formula
    formula <- resolveDataObject(formula)
  }
  # a tree prior given as the tenth argument, where dbarts 0.9-34 took
  # 'tree.prior' and 'power' now stands, is the retired 'tree.prior'
  powerWritten <- c(
    "power" %in% names(matchedCall),
    "power" %in% names(sys.call())
  )
  if (
    powerWritten[1L] &&
      !powerWritten[2L] &&
      isPositionalTreePrior(matchedCall[["power"]], parent.frame())
  ) {
    names(matchedCall)[names(matchedCall) == "power"] <- "tree.prior"
  }
  # the argument names as the caller wrote them, before a retired one is
  # restated under its successor
  writtenArgs <- names(matchedCall)

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
  # the retired 'tree.prior' is read as the plain arguments it states, which
  # are what this call states of them from here on
  retiredTreePrior <- readRetiredTreePrior(matchedCall, evalEnv, "dbarts")
  writtenArgs <- unique(c(
    setdiff(writtenArgs, "tree.prior"),
    names(retiredTreePrior$args)
  ))
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
  # a number, or the name of the rule that estimates one; the retired
  # spelling takes a number only
  sigestRule <- "auto"
  if (sigmaSupplied) {
    if (!isSingleNA(sigest)) {
      sigest <- validateSigest(sigest, "dbarts")
    }
  } else {
    resolvedSigest <- resolveSigestRule(sigest, "dbarts")
    sigest <- resolvedSigest$value
    sigestRule <- resolvedSigest$rule
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
  # its default is NULL, so a written NULL is the name left out
  if (length(nullWritten(matchedCall, "node.prior", evalEnv)) > 0L) {
    nodePriorSupplied <- FALSE
    matchedCall$node.prior <- NULL
    writtenArgs <- setdiff(writtenArgs, "node.prior")
  }
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

  # the formula as written, kept for the call a sampler records
  writtenFormula <- NULL
  # judged without forcing: power may hold a retired tree prior in its
  # vocabulary
  for (name in nullWritten(matchedCall, c("power", "base"), evalEnv)) {
    refuseNullArgument(NULL, name)
  }
  # the fitting function's power and base are the forest with no basis's, so a
  # formula that writes that forest with one of its own states it twice
  givenPlain <- intersect(c("power", "base"), writtenArgs)
  if (length(givenPlain) > 0L && is.formula(formula)) {
    refuseGivenTwiceInFormula(formula, givenPlain)
  }
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
        (length(termIngestion$basisReads) > 1L ||
          !all(vapply(termIngestion$basisReads, is.null, logical(1L))))
    ) {
      stop(
        "a forest() formula term does not support 'test'; drop the term or ",
        "fit a single-forest model"
      )
    }
    writtenFormula <- matchedCall$formula
    formula <- termIngestion$formula
    matchedCall$formula <- formula
  }

  validateCall <- redirectCall(
    matchedCall,
    quoteInNamespace(validateArgumentsInEnvironment)
  )
  validateCall <- addCallArgument(validateCall, 1L, sys.frame(sys.nframe()))
  validateCall <- addCallArgument(validateCall, 2L, dbarts::dbarts)
  validateCall <- addCallArgument(validateCall, 3L, "dbarts")
  eval(validateCall, evalEnv, getNamespace("dbarts"))
  # the count and the cut settings as this call states them, flat or written
  # on the control; a wrapper has read them and passes them flat. Resolved
  # before any response is read, so that a hazard fit's memory can be
  # projected from them before its rows are built
  copies <- readDoorCopies(
    releaseCarriedControl(control, "dbarts"),
    list(n.trees = n.trees, n.cuts = n.cuts, useQuantiles = useQuantiles),
    writtenArgs
  )
  control <- copies$control

  # survival response ingestion: a survival::Surv
  # object or an explicit family = "aft" with a two-column (time, status)
  # response fits the AFT log-normal model. The matrix (x.train, y.train)
  # interface is handled directly here, the response being the second
  # positional argument; the log event/censoring time replaces it and the
  # status rides the data object, whose status slot the bartcore survival
  # family reads.
  # 'subset' is honored (applied to the status vector alongside the
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
  # a hazard fit with no subject in takes its periods from every subject;
  # held for the creation, whose leaves then read every row too, to say both
  # in one text
  heldPeriods <- NULL
  holdPeriods <- function(expr) {
    withCallingHandlers(expr, dbartsScaleFallbackWarning = function(w) {
      heldPeriods <<- w
      invokeRestart("muffleWarning")
    })
  }
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
  # only as the period grid, recorded on the model, which the packaging reads
  # into the $periods marker (as the data object's status becomes $status).
  # No status vector reaches C++ - the censoring is baked into y'.
  hazardPeriods <- NULL
  # the person-period row names, set on the data object once it exists, and
  # whether the expansion ran before the na.action did
  hazardNames <- NULL
  hazardExpandedFirst <- FALSE
  hazardOffsetArgument <- NULL
  # each row's subject, where a hazard fit repeats a subject's basis over its
  # periods after the data object is built
  basisExpansion <- NULL
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
    expansion <- holdPeriods(expandDiscreteTimeHazard(
      xForExpansion,
      timeForExpansion,
      statusForExpansion,
      hazardRunSettings(control, control@n.trees),
      breaks = breaks,
      offset = offsetForExpansion,
      weights = weightsForExpansion,
      gridTime = if (!is.null(keptSubjects)) {
        hazardGridTimes(
          timeForExpansion[keptSubjects$keep],
          if (length(weightsForExpansion) > 1L) {
            weightsForExpansion[keptSubjects$keep]
          } else {
            weightsForExpansion
          }
        )
      } else {
        hazardGridTimes(timeForExpansion, weightsForExpansion)
      },
      n.test = if (!missing(test)) NROW(test) else 0L
    ))
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
    # aft's row count is unchanged, unlike hazard's, so the data object
    # subsets x and the log times itself and only the status vector is cut
    # here: by the value of 'subset' that the data object is handed below,
    # so that the two are cut to the same rows
    if (!missing(subset)) {
      survivalStatus <- survivalStatus[subset]
    }
  }

  # multinomial (K-forest softmax): the response is
  # an n x K count matrix, which the response vector every other family takes
  # has no shape for, so it rides the data object's own 'counts' argument -
  # where 'subset' reaches it, where the validity method constrains it, and
  # where a re-created or reloaded sampler finds it without further discipline.
  # A factor or character response is the single-trial special case, one-hot
  # expanded here; a numeric vector is refused. The matrix interface only, as
  # the survival families are: a formula LHS carries no count matrix, and
  # data@y is the DERIVED trials vector rather than anything the caller wrote.
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
    multinomialAsFactor <- !is.data.frame(data) &&
      (!is.matrix(data) || (is.character(data) && ncol(data) == 1L))
    matchedCall$data <- NULL
  }
  # the binomial on the matrix interface: a two-column response under a
  # logistic family is successes and failures, which ride the data object's
  # counts as a multinomial count matrix does
  binomialCounts <- NULL
  # a hazard family has taken its (time, status) pair already, its family
  # remapped to the binary link
  if (
    directResponse &&
      !(requestedFamily %in% hazardTokens) &&
      is.matrix(data) &&
      ncol(data) > 1L
  ) {
    if (identical(family, "logistic")) {
      if (ncol(data) != 2L) {
        stop(binomialCountColumnsMessage(ncol(data), "the response is"))
      }
      binomialCounts <- validateBinomialCounts(data)
      matchedCall$data <- NULL
    } else if (identical(family, "probit") && ncol(data) == 2L) {
      stop(probitBinomialMessage)
    }
  }

  if (length(control@call) == 1L && control@call == call("NA")) {
    recordedCall <- matchedCall
    # the formula as written, forest() terms in it, not the one they were
    # taken out of
    if (!is.null(writtenFormula)) {
      recordedCall$formula <- writtenFormula
    }
    control@call <- expandForwardedCall(recordedCall, evalEnv)
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
  # session. The row count it names is the control's own testFitParallelCutoff,
  # so a caller who moved that cutoff is told the number actually in force.
  if (is.na(control@n.threads)) {
    stop(naThreadsControlMessage, call. = FALSE)
  }
  if (control@n.threads > control@n.chains) {
    warnOnce(
      excessThreadsKey,
      sprintf(
        paste0(
          "n.threads (%d) exceeds n.chains (%d); tree sampling uses at ",
          "most one thread per chain, so the extra threads serve only ",
          "test-set fitting above %d rows and predict. Set n.threads to at ",
          "most n.chains to avoid this warning, which is shown once per ",
          "session."
        ),
        control@n.threads,
        control@n.chains,
        control@testFitParallelCutoff
      )
    )
  }

  dataCall <- redirectCall(matchedCall, quoteInNamespace(dbartsData))
  # read above, and set once the data object has the columns it fits
  dataCall[c("n.cuts", "useQuantiles")] <- NULL
  dataCall$keep.data <- FALSE
  # Every argument that holds rows is evaluated once for a fit, and what this
  # function has read is handed to the data object as its value (handOn):
  # evaluated again there, an argument that draws its rows would give the fit
  # other rows than the ones read here, the predictors another draw than the
  # response, a basis or the censoring status another than either. The first
  # argument and 'data' have been read above, in that order, by this
  # function's own arguments: the formula, the predictors or a data object
  # written in the call, which is built once, and the data or the response,
  # for which a hazard or an aft fit has put the response it made of it in
  # the call. A formula's 'subset' is an expression over the data's columns,
  # which the data object evaluates, once, in the 'data' it is handed; with
  # the matrix interface it is an index, evaluated here. The arguments not
  # read here go on as written and are evaluated once there.
  handed <- list(
    formula = if (hazardExpandedFirst) dataCall$formula else formula
  )
  if ("data" %in% names(dataCall)) {
    handed["data"] <- list(
      if (is.language(dataCall$data)) data else dataCall$data
    )
  }
  if (
    !is.formula(formula) &&
      !inherits(formula, "dbartsData") &&
      "subset" %in% names(dataCall)
  ) {
    handed["subset"] <- list(subset)
  }
  # The bases of the model's forests, whichever door declared them, are read
  # by one function (readFitBasis): a term's in the formula's ingestion, where
  # which forest has one decides which is which, and a list's here, against
  # the data this call was given, a data object's against what it holds. A
  # value has one entry per row handed over and code is read over every row;
  # both are cut to the rows the data object keeps once it is built, below.
  # A hazard fit's has one value per subject.
  basisReads <- NULL
  basisSource <- NULL
  basisDeclarations <- forestBasisDeclarations(forests)
  if (!is.null(termIngestion)) {
    basisReads <- termIngestion$basisReads
    basisSource <- termIngestion$basisSource
    basisSource$hazard <- requestedFamily %in% hazardTokens
    basisSource$subset <- !missing(subset)
  } else if (!is.null(basisDeclarations)) {
    # missing() reads this frame, so the data are resolved here
    basisSource <- if (inherits(formula, "dbartsData")) {
      objectBasisSource(dataObject, forests)
    } else if (is.formula(formula)) {
      formulaBasisSource(formula, if (missing(data)) NULL else data)
    } else {
      matrixBasisSource(formula, forests)
    }
    basisSource$hazard <- requestedFamily %in% hazardTokens
    if (!inherits(formula, "dbartsData")) {
      basisSource$subset <- !missing(subset)
    }
    declared <- readDeclaredBases(forests, basisSource)
    forests <- declared$forests
    basisReads <- declared$reads
  }
  if (
    all(vapply(basisReads, is.null, logical(1L))) && length(basisReads) < 2L
  ) {
    # one forest with no basis is the single-forest fit, and two or more are
    # a model of several sharing one prior total
    basisReads <- NULL
  }
  if (!is.null(multinomialCounts)) {
    handed$counts <- multinomialCounts
  }
  if (!is.null(binomialCounts)) {
    handed$counts <- binomialCounts
  }
  # a formula's two-column response under a logistic family is the binomial's
  # successes and failures, which the data object reads off the model frame,
  # after 'subset' and the na.action
  if (
    identical(family, "logistic") &&
      !(requestedFamily %in% hazardTokens) &&
      is.formula(handed$formula)
  ) {
    attr(handed$formula, binomialResponseAttr) <- TRUE
  }
  dataCall <- handOn(dataCall, handed)
  data <- withMatrixResponseRestated(
    "bart()/dbarts()",
    requestedFamily,
    withBinaryResponsePrecision(family, eval(dataCall, evalEnv))
  )
  # the rows a basis is cut to: those the data object kept, the positions of
  # its rows among the rows handed over, or, for a hazard fit on the matrix
  # interface, which expanded its subjects first, each period row's subject
  # among the subjects handed over
  basisRows <- data@rows
  if (hazardExpandedFirst && !is.null(basisReads)) {
    subjects <- seq_len(NROW(formula))
    if (!missing(subset)) {
      subjects <- subjects[subset]
    }
    basisRows <- subjects[expansion$subject]
    if (!is.null(data@rows)) {
      basisRows <- basisRows[data@rows]
    }
  }

  if (!is.null(multinomialCounts)) {
    data@counts.as.factor <- multinomialAsFactor
  }

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
      expansion <- holdPeriods(expandDiscreteTimeHazard(
        data@x,
        formulaSurvivalTime,
        formulaSurvivalStatus,
        hazardRunSettings(control, control@n.trees),
        breaks = breaks,
        offset = data@offset,
        weights = data@weights,
        gridTime = hazardGridTimes(formulaSurvivalTime, data@weights),
        n.test = NROW(data@x.test)
      ))
      # makeModelMatrix already typed the original columns (categorical vs
      # ordinal); the appended period column is ordinal by construction
      # (dec-B97) and rides last, so one more entry keeps the two aligned
      data@x <- expansion$x
      data@varTypes <- c(data@varTypes, ORDINAL_VARIABLE)
      # the rows are the subjects' periods now, not rows handed over
      basisExpansion <- expansion$subject
      data@rows <- NULL
      data@rows.total <- NULL
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
      data <- expandHazardTest(data, K)
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

  data <- applyCutSettings(data, copies, !is.null(hazardPeriods))
  # A new fit calibrates from the response it is given: the response's scale
  # and sigma's starting estimate are read from these data's response, not
  # from the records a sampler wrote into a data object it held, unless the
  # caller states sigest. The count model's center and a survival fit's span
  # ride the range record and the calibration map's anchor and unit the other,
  # so they follow. The predictor-derived records (cut points, a leaf
  # covariate's standardization) are data and stay.
  data <- clearResponseCalibration(data)
  data@sigma <- sigest

  # Every basis is built now, on the rows the data object kept, a hazard
  # fit's repeated over each subject's periods (buildFitBases). data@bases is
  # positional against the forests, and a forests list states them in full:
  # what a data object handed in held is replaced, a held basis a forest
  # names kept.
  basisRecords <- NULL
  if (!is.null(basisReads) || inherits(formula, "dbartsData")) {
    built <- if (!is.null(basisReads)) {
      buildFitBases(
        basisReads,
        basisSource,
        basisRows,
        length(data@y),
        if (is.null(termIngestion)) forests else termIngestion$forests,
        expand = basisExpansion
      )
    }
    data <- installFitBases(data, built)
    basisRecords <- fitBasisRecords(built)
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

  # where a numeric column's scale is read from every row with the
  # response's, one warning names both; the response's alone is held until
  # the engine has said
  heldWarning <- NULL
  spec <- withCallingHandlers(
    resolveSamplerSpec(
      withRetiredTreePrior(matchedCall, retiredTreePrior),
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
      forests = forests,
      evalEnv = evalEnv,
      residPrior = residPrior,
      familySpec = familySpec,
      basisRecords = basisRecords,
      sigestRule = sigestRule,
      treePrior = retiredTreePrior$record,
      weightsType = validateWeightsType(weights.type)
    ),
    dbartsScaleFallbackWarning = function(w) {
      heldWarning <<- w
      invokeRestart("muffleWarning")
    }
  )

  # a latent family's 0/1 case weights are membership, which ride the data
  # object and which the engine composes with any mask installed later
  sampler <- new(
    "dbartsSampler",
    spec$control,
    spec$model,
    spec$data,
    periods = heldPeriods
  )
  # the sampler raised one text naming the response with whatever else was
  # read from every row, a multiplier's scale or a leaf's standardization
  if (
    !is.null(heldWarning) &&
      length(basisScaleFallbacks(sampler)) == 0L &&
      length(leafFromAllRows(
        sampler,
        .Call(C_dbarts_bartcore_derivedRecords, sampler$pointer, FALSE)
      )) ==
        0L
  ) {
    warning(heldWarning)
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

# The warm start behind the sampler's installTrees method and bart's
# 'warm.start': the donor's trees go into the sampler, which keeps its own
# model, cut points and data. Without force a donor that would have to be
# changed to fit is not installed and the sampler is left as it was, its
# stored state included; with it such a donor goes in with the changes. The
# predictors are handed over for a forced install of a donor recorded on
# other cut points, which re-quantizes them on the way. Returns the engine's
# report, a list whose 'installed' is FALSE only for a donor declined; the
# sampler's state is stored again only when the donor went in.
installDonorTrees <- function(sampler, donor, samples, force) {
  ptr <- sampler$getPointer()
  refuseMultiForestWarmStart(ptr, "$installTrees")
  donorState <- warmStartState(donor)
  if (!is.null(samples)) {
    samples <- coerceOrError(samples, "integer")
  }
  report <- .Call(
    C_dbarts_bartcore_installForests,
    ptr,
    donorState,
    samples,
    rawPredictorMatrix(sampler$data@x),
    force
  )
  if (report$installed) {
    sampler$storeState(ptr)
  }
  report
}

# Draws from the BART prior (issue #31): repeatedly redraws the hyperparameters,
# trees and leaf values on a private sampler and evaluates the resulting forest
# at x.test, for calibrating priors before fitting - e.g. the prior
# distribution of a treatment effect via f(x1) - f(x0). Draws under the prior
# the sampler runs under: k, a Student-t df, a negative-binomial shape, the
# ordinal cutpoints and DART's split weights are drawn from their priors for
# each sample, and each family's observation layer is its own. Never touches
# the caller's sampler: a fresh construction gets its own external pointer,
# and while it borrows the caller's data object, nothing here mutates it.
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
  if (samplerCarriesAmplitudes(sampler)) {
    stop(
      "samplePriorPredictive does not support a sampler of several forests: ",
      "their coefficients are not drawn from the prior here",
      call. = FALSE
    )
  }
  if (leafPriorIsImproper(sampler$model)) {
    stop(
      "samplePriorPredictive does not support an improper leaf prior ",
      "(invchi() with scale 0, or chi() with an infinite scale): it has no ",
      "prior draws",
      call. = FALSE
    )
  }
  type <- match.arg(type)
  n.samples <- coerceOrError(n.samples, "integer")[1L]
  if (is.na(n.samples) || n.samples <= 0L) {
    stop("'n.samples' must be a positive integer")
  }
  isMultinomial <- samplerCarriesCounts(sampler)
  if (
    isMultinomial &&
      is.null(offset.test) &&
      (!is.null(sampler$data@offset.category) ||
        !is.null(sampler$data@offset.category.test))
  ) {
    stop(
      "'offset.test' must be given on a sampler that carries a category ",
      "offset: one row per row of 'x.test' and one column per category (an ",
      "all-zero matrix for the offset-free surface)",
      call. = FALSE
    )
  }

  # a fresh sampler, not sampler$copy(): copy() installs the caller's saved
  # state - including the engine RNG - so successive calls would replay one
  # frozen stream. Fresh creation seeds the chain RNGs from R's stream when
  # control@seed is NA (or pins them when it is set), giving independent
  # draws across calls by default with set.seed governing reproducibility.
  # The prior draws overwrite all tree state, so no donor state is needed.
  # keepTrees is forced off: it makes predict() serve saved posterior
  # samples instead of the live trees this function just drew. Prior draws
  # are chain-free, so the sampler runs one chain.
  newControl <- sampler$control
  newControl@keepTrees <- FALSE
  newControl@n.chains <- 1L
  # under the caller's prior, on the grid, range and leaf standardization its
  # data object records
  draw <- dbartsSampler$new(newControl, sampler$model, sampler$data)

  xt <- if (is.null(x.test)) extract(draw, "predictors") else x.test
  responseIsBinary <- draw$control@binary
  family <- draw$model@family

  # the noise a heteroscedastic prior predictive adds is s(x) eps, with s^2(x)
  # drawn from the variance forest's own prior once per sample and read at the
  # rows being predicted; the scalar draws below would report a homoscedastic
  # prior predictive instead. "ev" needs no noise term and is unaffected.
  drawsVariance <- type == "ppd" &&
    !responseIsBinary &&
    !is.null(draw$model@variance)
  if (drawsVariance) {
    # the surface is read through the draw sampler's own test rows, which is
    # where the variance forest evaluates off the training data; predict()
    # cannot serve it, keepTrees being off and these trees never recorded
    withCodedPredictors(draw$setTestPredictorAndOffset(xt, NULL))
  }

  # sigma is a draw of its own prior, drawn here, ahead of the forests, as it
  # always has been, so a gaussian sampler's "ppd" is what it was
  scaleFamily <- !responseIsBinary &&
    !isMultinomial &&
    !family %in% c("nbinom", "ordinal")
  sigmaDraws <- NULL
  if (type == "ppd" && scaleFamily && !drawsVariance) {
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

  forests <- drawPriorForests(
    draw,
    xt,
    offset.test,
    n.samples,
    n.threads,
    drawsVariance
  )
  hyper <- forests$hyperparameters
  residualDf <- hyper[[1L]]$residual.df
  if (isMultinomial) {
    # n.samples x rows x K: each sample's softmax probabilities
    result <- aperm(simplify2array(forests$fits), c(3L, 1L, 2L))
  } else {
    result <- do.call(rbind, forests$fits)
  }

  if (responseIsBinary) {
    result <- probabilityFromLatents(result, list(family = family))
  } else if (identical(family, "nbinom")) {
    # predict serves the log mean, the center and the offset in it
    result <- exp(result)
  } else if (identical(family, "ordinal")) {
    # one cutpoint set per sample, shared by every row
    result <- aperm(
      simplify2array(lapply(seq_len(n.samples), function(i) {
        ordinalCategoryProbabilities(result[i, ], hyper[[i]]$thresholds)
      })),
      c(3L, 1L, 2L)
    )
  }

  if (type == "ppd") {
    if (responseIsBinary) {
      result <- matrix(
        rbinom(length(result), 1L, result),
        nrow(result),
        ncol(result)
      )
    } else if (isMultinomial) {
      # a draw in the shape the sampler's response was written in, whatever
      # the rows' trials: a category for a factor, a count vector of one
      # trial for a count matrix
      result <- if (isTRUE(dataCountsAsFactor(draw$data))) {
        multinomialPpdFromProbs(result)
      } else {
        multinomialCountPpdFromProbs(result, rep_len(1L, dim(result)[2L]))
      }
    } else if (identical(family, "ordinal")) {
      result <- multinomialPpdFromProbs(result)
    } else if (identical(family, "nbinom")) {
      result <- negbinPpd(
        result,
        vapply(hyper, function(h) h$shape, numeric(1L))
      )
    } else if (drawsVariance) {
      # s^2(x) is a variance on the response scale, one row per prior draw and
      # one column per predicted row, so the noise scale is its square root
      result <- result +
        sqrt(do.call(rbind, forests$variances)) *
          matrix(
            rnorm(length(result)),
            nrow(result),
            ncol(result)
          )
    } else if (!is.null(residualDf)) {
      # Student-t noise at each sample's df, sigma its conditional scale
      dfs <- vapply(hyper, function(h) h$residual.df, numeric(1L))
      result <- result +
        matrix(
          rep(sigmaDraws, ncol(result)) *
            rt(length(result), rep(dfs, ncol(result))),
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

# Whether a model's leaf prior is improper, a chi() k prior with an infinite
# scale (which invchi() with scale 0 becomes): it has no prior draws.
leafPriorIsImproper <- function(model) {
  hyperprior <- model@leaf.hyperprior
  inherits(hyperprior, "dbartsChiHyperprior") &&
    !is.finite(hyperprior@scale)
}

# The projected length, in seconds, of a call whose DART draws earn one
# warning. An unexported, undocumented option, read fresh per call, sets it
# instead, so a test can force the warning at any seed.
priorPredictiveWarnSeconds <- function() {
  forced <- getOption("dbarts.priorPredictiveWarnSeconds")
  if (is.null(forced)) {
    return(60)
  }
  seconds <- as.double(forced)[1L]
  if (is.na(seconds) || seconds < 0) {
    stop(
      "option 'dbarts.priorPredictiveWarnSeconds' must be a non-negative ",
      "number of seconds"
    )
  }
  seconds
}

# A projected length in words: seconds, minutes or hours.
formatPriorPredictiveSeconds <- function(seconds) {
  if (seconds < 120) {
    sprintf("%.0f seconds", seconds)
  } else if (seconds < 7200) {
    sprintf("%.0f minutes", seconds / 60)
  } else {
    sprintf("%.1f hours", seconds / 3600)
  }
}

# samplePriorPredictive's forests: n.samples independent draws of the prior
# the private sampler 'draw' runs under, each evaluated at xt. Each sample
# draws the hyperparameters first, then the trees, then the leaves. Under
# dart() the split weights and the trees are one draw, redrawn whole, with no
# cap, while any tree is one the sampler refuses. The engine makes the
# attempts in calls of a bounded count, sized to about a tenth of a second,
# and polls for an interrupt after each, so every call leaves the sampler in
# a state it holds; a call that accepts nothing has drawn k and the family's
# hyperparameters too, which the next redraws, as they do not bear on
# acceptance. Once a second (scaled with the threshold) has gone into the
# attempts, each call projects the whole: the time so far times n.samples
# over the samples accepted, or with none accepted yet, the time so far times
# n.samples as a lower bound. The first projection over
# priorPredictiveWarnSeconds() warns, once, printed when raised so that it
# arrives while the draw runs. Returns list(fits, variances, hyperparameters),
# one entry a sample.
drawPriorForests <- function(
  draw,
  xt,
  offset.test,
  n.samples,
  n.threads,
  drawsVariance
) {
  ptr <- draw$getPointer()
  isMultinomial <- samplerCarriesCounts(draw)
  warnSeconds <- priorPredictiveWarnSeconds()
  elapsed <- function() proc.time()[["elapsed"]]
  attemptsPerCall <- 1L
  attempts <- 0
  accepted <- 0L
  started <- NA_real_
  warned <- FALSE

  projectDuration <- function() {
    spent <- elapsed() - started
    if (warned || spent < warnSeconds / 60) {
      return(invisible(NULL))
    }
    if (accepted > 0L) {
      projected <- spent * n.samples / accepted
      if (projected < warnSeconds) {
        return(invisible(NULL))
      }
      message <- paste0(
        accepted,
        " of ",
        format(attempts, big.mark = ",", scientific = FALSE),
        " joint draws of the split weights and trees so far were ones the ",
        "sampler holds (acceptance ",
        signif(accepted / attempts, 2L),
        "), so the call is projected to take about ",
        formatPriorPredictiveSeconds(projected)
      )
    } else {
      projected <- spent * n.samples
      if (projected < warnSeconds) {
        return(invisible(NULL))
      }
      message <- paste0(
        "none of ",
        format(attempts, big.mark = ",", scientific = FALSE),
        " joint draws of the split weights and trees so far was one the ",
        "sampler holds (acceptance none in ",
        format(attempts, big.mark = ",", scientific = FALSE),
        "), so the call will take at least about ",
        formatPriorPredictiveSeconds(projected)
      )
    }
    warned <<- TRUE
    warning(
      "samplePriorPredictive: under dart() ",
      message,
      "; it can be interrupted",
      call. = FALSE,
      immediate. = TRUE
    )
  }

  # the entry's draw, repeated under dart() until the whole forest is held
  drawHyperparameters <- function() {
    repeat {
      if (is.na(started)) {
        started <<- elapsed()
      }
      callStarted <- elapsed()
      drawn <- .Call(
        C_dbarts_bartcore_sampleHyperparametersFromPrior,
        ptr,
        attemptsPerCall
      )
      if (is.null(drawn$split.probabilities)) {
        return(drawn)
      }
      callSeconds <- elapsed() - callStarted
      attempts <<- attempts + drawn$attempts
      if (drawn$forests.drawn) {
        accepted <<- accepted + 1L
      }
      if (callSeconds < 0.05) {
        attemptsPerCall <<- min(2L * attemptsPerCall, 1048576L)
      } else if (callSeconds > 0.2) {
        attemptsPerCall <<- max(attemptsPerCall %/% 2L, 1L)
      }
      projectDuration()
      if (drawn$forests.drawn) {
        return(drawn)
      }
    }
  }

  sampleOnce <- function() {
    drawn <- drawHyperparameters()
    if (!drawn$forests.drawn) {
      draw$sampleTreesFromPrior(updateState = FALSE)
    }
    draw$sampleLeafParametersFromPrior(updateState = FALSE)
    fit <- withCodedPredictors(draw$predict(xt, offset.test, n.threads))
    # one chain's fit, as a vector of rows; a multinomial sampler's rows
    # carry K probabilities, kept as rows x K
    d <- dim(fit)
    if (length(d) == 3L) {
      fit <- array(fit[,, 1L], d[1:2])
    } else if (length(d) == 2L && !isMultinomial) {
      fit <- fit[, 1L]
    }
    variance <- NULL
    if (drawsVariance) {
      draw$sampleVarianceForestFromPrior(updateState = FALSE)
      variance <- draw$getVariance(test = TRUE)[, 1L]
    }
    list(fit = fit, variance = variance, hyperparameters = drawn)
  }

  samples <- lapply(seq_len(n.samples), function(i) sampleOnce())
  list(
    fits = lapply(samples, `[[`, "fit"),
    variances = lapply(samples, `[[`, "variance"),
    hyperparameters = lapply(samples, `[[`, "hyperparameters")
  )
}

## Whether a method that just changed the sampler should refresh its cached
## state (the copy storeState/setState carry): an explicit TRUE or FALSE
## wins, and NULL - every one of these methods' own default - resolves
## against control@updateState, exactly as run() has always resolved it.
## A number reads as as.logical reads it; anything else is refused by name. 0.9-x's default was NA, which is a
## missing value now and reads as NULL, after one warning.
checkUpdateState <- function(updateState) {
  if (is.null(updateState)) {
    return(NULL)
  }
  if (isSingleNA(updateState)) {
    warnNAForNull("updateState", "dbartsSampler")
    return(NULL)
  }
  if (
    !(is.logical(updateState) || is.numeric(updateState)) ||
      length(updateState) != 1L
  ) {
    shown <- paste(deparse(updateState), collapse = " ")
    if (nchar(shown) > 40L) {
      shown <- paste0(substr(shown, 1L, 37L), "...")
    }
    stop(
      "'updateState' must be TRUE, FALSE or NULL, not ",
      shown,
      call. = FALSE
    )
  }
  as.logical(updateState)
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
  counts <- if (samplerCarriesCounts(sampler)) dataCounts(sampler$data)
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

## The chain the engine at ptr holds, as a stored state, with the level names
## of each forest a factor multiplies as the data object records them: the
## state's "basisLevels" attribute, one entry per forest, NULL for a forest
## of no factor. Written only where some forest is of a factor, so every other
## state is the object the engine wrote. Every store of a sampler's chain goes
## through here.
storedState <- function(ptr, data) {
  state <- .Call(C_dbarts_bartcore_storeState, ptr)
  levels <- dataBasisLevels(data)
  if (!all(vapply(levels, is.null, NA))) {
    attr(state, "basisLevels") <- levels
  }
  state
}

## How a state's coefficients go onto the levels a data object records: NULL
## where they go as stored, the state recording no levels (it holds no forest
## of a factor, or was stored before the record) or every factor's levels in
## the data's order; otherwise a list with, for each forest, NULL or the
## stored position of each of the data's levels. Each forest must multiply no
## factor on both sides, or a factor of the same level names on both, in any
## order; any other state is refused, whatever forceUpdate says.
stateLevelOrder <- function(state, data) {
  recorded <- attr(state, "basisLevels", exact = TRUE)
  if (is.null(recorded)) {
    return(NULL)
  }
  if (
    !is.list(recorded) ||
      !all(vapply(
        recorded,
        function(levels) is.null(levels) || is.character(levels),
        NA
      ))
  ) {
    stop(
      "the state's 'basisLevels' record must be a list of NULL or character ",
      "entries",
      call. = FALSE
    )
  }
  own <- dataBasisLevels(data)
  describe <- function(levels) {
    if (is.null(levels)) {
      "no factor"
    } else {
      paste(
        "a factor with the levels",
        quotedLevelList(levels, comma = TRUE)
      )
    }
  }
  order <- vector("list", length(recorded))
  moved <- FALSE
  for (forest in seq_len(max(length(recorded), length(own)))) {
    stored <- if (forest <= length(recorded)) recorded[[forest]]
    held <- if (forest <= length(own)) own[[forest]]
    if (is.null(stored) && is.null(held)) {
      next
    }
    position <- if (is.character(stored)) match(held, stored)
    if (
      is.null(position) ||
        length(stored) != length(held) ||
        anyNA(position) ||
        anyDuplicated(position) > 0L
    ) {
      stop(
        "the state's forest ",
        forest,
        " multiplies ",
        describe(stored),
        ", and this sampler's forest ",
        forest,
        " multiplies ",
        describe(held),
        ": a state's coefficients go only onto levels of the same names, ",
        "in any order",
        call. = FALSE
      )
    }
    if (!identical(position, seq_along(held))) {
      order[[forest]] <- position
      moved <- TRUE
    }
  }
  if (moved) order else NULL
}

## The install every re-creation shares - getPointer, $setState and copy: the
## state goes into the engine at ptr, then what rides no state is put back from
## the sampler's own mirrors, in this order. adopt is TRUE only where the
## engine was just created from a control that can name another saved-tree
## capacity than the state holds. force installs a state a tree of which has
## to be changed to fit, changing it; without it such a state is not
## installed, the engine and the mirrors are left alone and 0L is returned.
## The engine's cut points and weights are never written: a re-created engine
## holds the ones its data object records. Returns 1L for a state that went
## in, and 2L for one that went in onto cut points or weights other than the
## ones it records.
installStateOnto <- function(sampler, ptr, state, predictors, adopt, force) {
  levelOrder <- stateLevelOrder(state, sampler$data)
  status <- .Call(
    C_dbarts_bartcore_setState,
    ptr,
    state,
    predictors,
    adopt,
    force,
    levelOrder
  )
  if (status == 0L) {
    return(0L)
  }
  sampler$reapplyForestWeights(ptr)
  sampler$reapplyActiveRows(ptr)
  status
}

## The forests whose numeric column's scale the engine read from every row,
## too few being in the likelihood to read it from: the engine's status, raised
## only by a first creation, which takes the scale.
basisScaleFallbacks <- function(sampler) {
  ptr <- sampler$pointer
  numForests <- length(sampler$data@bases)
  Filter(
    function(forest) {
      raw <- .Call(C_dbarts_bartcore_getLeafPrior, ptr, forest - 1L)
      identical(raw[[1L, "basis.scale.all.rows"]], 1)
    },
    seq_len(numForests)
  )
}

## Refuses $setForestBasis(updateBasisScale = TRUE) on a forest whose default is
## not taken from its column, from the engine's raw leaf-prior row 'raw' of the
## forest called 'name': one that states its sd and one multiplied by a factor.
## The swap itself refuses a forest multiplied by nothing, and the engine
## refuses all three, under a generic text.
refuseBasisScaleUpdate <- function(raw, name) {
  stated <- raw[[1L, "named.sd"]]
  if (!is.nan(stated)) {
    stop(
      name,
      " states its sd, ",
      format(stated),
      ", so it has no default to update; restate it with ",
      "$setForest or leave updateBasisScale FALSE",
      call. = FALSE
    )
  }
  multiplier <- raw[[1L, "multiplier"]]
  if (identical(multiplier, 1)) {
    stop(
      name,
      " is multiplied by a factor, so its default is one unit and reads no ",
      "column; leave updateBasisScale FALSE",
      call. = FALSE
    )
  }
  invisible(NULL)
}

## Records on a sampler's data the scale the engine holds for the numeric
## forest at 0-based 'index', overwriting the entry there: after
## $setForestBasis(updateBasisScale = TRUE), whose scale every later creation
## from the object is to be handed.
recordBasisScale <- function(sampler, ptr, index) {
  raw <- .Call(C_dbarts_bartcore_getLeafPrior, ptr, index)
  data <- sampler$data
  scales <- dataBasisScale(data)
  if (is.null(scales)) {
    scales <- vector("list", length(data@bases))
  }
  length(scales) <- max(length(scales), length(data@bases))
  scales[[index + 1L]] <- raw[[1L, "basis.scale"]]
  data@basis.scale <- scales
  sampler$dataHeld <- data
  invisible(NULL)
}

## The warning that $setForestBasis(updateBasisScale = TRUE) kept a numeric
## forest's scale in force, on the bridge's answer: 1 where fewer than two of
## the 'n' rows were in the likelihood, 'rowsIn' of them, 2 where the forest
## 'name''s column was constant over them, 3 where their frequency counts
## totaled at most 1, each only on the sampler's first keep of any scale, and
## 0 otherwise.
warnBasisScaleKept <- function(answer, name, n, rowsIn) {
  if (identical(answer, 1L)) {
    warnClassed(
      "dbartsScaleFallbackWarning",
      "fewer than two observations are in the likelihood (",
      rowsIn,
      " of ",
      n,
      "), so the scale of ",
      name,
      "'s basis is kept as it was (said once per sampler)"
    )
  } else if (identical(answer, 2L)) {
    warnClassed(
      "dbartsScaleFallbackWarning",
      name,
      "'s basis is constant over the rows in the likelihood, so its scale ",
      "is kept as it was (said once per sampler)"
    )
  } else if (identical(answer, 3L)) {
    warnClassed(
      "dbartsScaleFallbackWarning",
      "the rows in the likelihood carry a total count of at most 1, so the ",
      "scale of ",
      name,
      "'s basis is kept as it was (said once per sampler)"
    )
  }
  invisible(NULL)
}

## Records on a sampler's data the scale the engine holds for each numeric
## forest's column, where the object has no entry for that forest: a first
## creation's, which every later creation from the object is handed. Raises
## the one warning where the engine read a scale from every row: under
## gaussian with the response's, whose scale was read from every row too;
## under probit and logistic counting the rows of positive weight the
## creation mask 'active' keeps, fewer than two.
recordBasisScales <- function(sampler, ptr, active = NULL) {
  data <- sampler$data
  numForests <- length(data@bases)
  scales <- dataBasisScale(data)
  if (is.null(scales)) {
    scales <- vector("list", numForests)
  }
  length(scales) <- max(length(scales), numForests)
  written <- FALSE
  for (forest in seq_len(numForests)) {
    raw <- .Call(C_dbarts_bartcore_getLeafPrior, ptr, forest - 1L)
    scale <- raw[[1L, "basis.scale"]]
    if (is.nan(scale) || !is.null(scales[[forest]])) {
      next
    }
    scales[[forest]] <- scale
    written <- TRUE
  }
  if (written) {
    data@basis.scale <- scales
    sampler$dataHeld <- data
  }
  fellBack <- basisScaleFallbacks(sampler)
  if (length(fellBack) == 0L) {
    return(invisible(NULL))
  }
  n <- length(data@y)
  counted <- fallbackCountsRows(sampler$model, data@weights, dataCounts(data))
  labels <- vapply(sampler$model@forests, methods::slot, "", "label")
  named <- vapply(fellBack, forestName, "", labels = labels)
  what <- if (length(named) == 1L) {
    paste0("the scale of ", named, "'s basis")
  } else {
    paste0(
      "the scales of ",
      paste(paste0(named, "'s"), collapse = " and "),
      " bases"
    )
  }
  if (sampler$model@family %in% c("probit", "logistic")) {
    warnClassed(
      "dbartsScaleFallbackWarning",
      fewRowsInLead(
        n,
        sum(rowsInLikelihood(data@weights, active, n)),
        byMask = !is.null(active) &&
          (is.null(data@weights) || all(data@weights > 0))
      ),
      ", so ",
      what,
      if (length(named) == 1L) " is" else " are",
      " taken from ",
      if (!isFALSE(counted)) {
        atTheirCounts(takenFromAll(n), counted)
      } else {
        takenFromAll(n)
      }
    )
  } else {
    warnClassed(
      "dbartsScaleFallbackWarning",
      "fewer than two distinct values of the response are in the ",
      "likelihood (",
      countInLikelihood(
        sum(rowsInLikelihood(data@weights, NULL, n)),
        n,
        "observation"
      ),
      if (fallbackCountsRows(sampler$model, data@weights)) {
        paste0(
          "), so the response's scale and the prior's center are taken from ",
          takenFromAll(n),
          ", and ",
          what,
          " from ",
          atTheirCounts(takenFromAll(n))
        )
      } else {
        paste0(
          "), so the response's scale, the prior's center and ",
          what,
          " are taken from ",
          takenFromAll(n)
        )
      }
    )
  }
  invisible(NULL)
}

## The response transform a sampler's leaf prior is anchored to and its
## chains hold their numbers in - (min, max) as a state's fit.scale holds it
## - is data, recorded on the data object as its response.range. A creation
## from an object that records one takes it as the engine's transform, and
## the chains are moved there at once: no install moves them, so a state that
## follows is read against the record only because they are already in it.
## The engine restates a named leaf-prior sd with the move.
applyAnchor <- function(ptr, data) {
  record <- dataResponseRange(data)
  if (!is.null(record)) {
    .Call(C_dbarts_bartcore_anchor, ptr, as.double(record))
  }
  invisible(ptr)
}

## Creates the engine for (control, model, data), on the values the data
## object records and deriving the ones it does not: the cut grids, the
## calibration map's anchor and unit and a numeric column's scale go in at
## construction, the range by the move applyAnchor makes. active is the mask
## handed to the creation door, handed only to a first creation, which reads
## a numeric column's scale over its rows; a re-creation is handed the
## recorded scale and installs its mask afterward. Returns the pointer and the
## engine's records.
createEngine <- function(
  control,
  model,
  data,
  active = NULL,
  tooFewUnits = FALSE
) {
  ptr <- .Call(
    C_dbarts_bartcore_create,
    control,
    model,
    data,
    bridgeFamily(control, model),
    if (!is.null(active)) as.double(active),
    if (tooFewUnits) TRUE
  )
  applyAnchor(ptr, data)
  records <- .Call(C_dbarts_bartcore_derivedRecords, ptr, TRUE)
  refuseCarriedOrderedGrid(data, records)
  list(pointer = ptr, records = records)
}

recreatePointer <- function(control, model, data) {
  createEngine(control, model, checkDerivedRecord(data, model))$pointer
}

## The refusal of a sampler whose engine is gone and that holds no stored
## state to rebuild it from, which a reload and a copy give alike.
stopWithoutStoredState <- function() {
  stop(
    "samplers cannot be re-created without a stored state; call ",
    "storeState() before serializing (see the Saving section of ",
    "?`dbartsSampler-class`)",
    call. = FALSE
  )
}

## The sd a model names for its leaf prior, read off its encoding: prior.scale
## over the reference k it rides, a fixed k or a chi prior's scale; NA where
## the prior is written with k. The reference is 2, so the quotient is the
## named sd exactly.
namedLeafSd <- function(model) {
  hyperprior <- model@leaf.hyperprior
  reference <- if (is(hyperprior, "dbartsChiHyperprior")) {
    hyperprior@scale
  } else {
    hyperprior@k
  }
  model@prior.scale / reference
}

## What prior.sd is the sd of, per leaf model.
priorSdOf <- function(leafModel) {
  switch(
    leafModel,
    linear = "coefficient",
    gp = "amplitude",
    "forest"
  )
}

## One forest's leaf prior, from the bridge's per-chain calibration: the
## specification in the terms it was named in, then the quantities the chains
## share. A fixed k is read off the engine, so it is what is in force; a
## prior on k comes from the model, gated by the engine's per-forest flag,
## since a map forest pins k whatever the model says; and an sd, fixed or
## under invchi(), comes from the model as named, the engine holding only its
## statement as k, which need not divide back to the same number. Every chain
## runs under the sampler's one prior and transform, which no install moves,
## so the first chain's reading is the sampler's.
##
## On a forest whose scale the calibration map sets, the specification is the
## forest's statement: forest(sd = ) with the number as written where one was
## stated, forest() where none was. `sd` is the number in force, in the
## response's units, and `sd.stated` whether it was stated; both are the
## engine's, and nothing is converted here. `multiplier` is what the engine
## was told the forest multiplies, and `amplitude` fixed() where its
## coefficient is held. A numeric forest's `sd` is per unit of its column and
## carries the column's name, and `basis.scale` is the scale the data object
## records for it. `index` is the forest's position from 0.
reportLeafPrior <- function(sampler, raw, index) {
  shared <- function(column) raw[[1L, column]]
  model <- sampler$model
  hyperprior <- model@leaf.hyperprior
  # a forest taking the fitting function's leaf prior runs under it and reads
  # as one forest does, its sd its part of the total
  sharer <- "share" %in% colnames(raw) && !is.nan(raw[1L, "share"])
  mapped <- !is.nan(raw[1L, "map.anchor"]) && !sharer
  drawn <- raw[1L, "k.has.hyperprior"] != 0
  anchor <- shared("prior.scale")
  spec <- model@leaf.prior
  sdNamed <- mapped || !is.na(model@prior.scale) || !is.null(spec@prior.sd)
  named <- spec@prior.sd
  spec@k <- NULL
  spec@prior.sd <- NULL
  if (!sdNamed) {
    spec@k <- if (drawn) hyperprior else shared("k")
  } else if (!is.null(named)) {
    spec@prior.sd <- named
  } else if (!mapped) {
    # a model whose encoding names an sd its specification does not carry
    sd <- namedLeafSd(model)
    spec@prior.sd <- if (drawn) invchi(hyperprior@degreesOfFreedom, sd) else sd
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
  # whether the fitting function's leaf prior is named by an sd the user
  # wrote, which is what a forest taking that prior takes its size from
  fittingSdStated <- !is.na(model@prior.scale) || !is.null(named)
  if (sharer) {
    return(sharerLeafPrior(prior, raw, fittingSdStated))
  }
  if (!mapped) {
    return(lonePriorRecord(prior, raw, sdNamed, drawn))
  }
  multiplier <- c("none", "factor", "numeric")[shared("multiplier") + 1]
  held <- shared("amplitude.held") != 0
  own <- !is.nan(shared("named.sd"))
  # a forest with no basis and a held coefficient, stating no sd of its own,
  # takes its size from the fitting function's leaf prior
  stated <- own || (multiplier == "none" && held && fittingSdStated)
  sd <- shared("map.sd")
  # a numeric forest's number is per unit of its column, and carries the
  # column's name where the basis has one
  column <- NULL
  if (multiplier == "numeric") {
    column <- colnames(sampler$data@bases[[index + 1L]])
    if (length(column) == 1L && !is.na(column) && nzchar(column)) {
      names(sd) <- column
    }
  }
  statement <- list(
    leaf.prior = if (own) forest(sd = shared("named.sd")) else forest(),
    sd = sd,
    sd.stated = stated,
    multiplier = multiplier
  )
  statement$amplitude <- statedCoefficientPrior(shared("amplitude.df"), held)
  prior <- c(statement, prior[-1L])
  # the size the sd is of: exact where the coefficient is held, and at a
  # coefficient of 1 where it is drawn
  prior$prior.sd.of <- switch(
    multiplier,
    none = "forest",
    factor = "level difference",
    numeric = "forest, per unit of basis"
  )
  if (!held) {
    prior$prior.sd.of <- paste(prior$prior.sd.of, "(at a coefficient of 1)")
  }
  # the scale the data object records for the column, which a default was
  # taken from
  if (multiplier == "numeric") {
    recorded <- dataBasisScale(sampler$data)
    if (index + 1L <= length(recorded) && !is.null(recorded[[index + 1L]])) {
      scale <- recorded[[index + 1L]]
      names(scale) <- names(sd)
      prior$basis.scale <- scale
    }
  }
  prior
}

## The reader's entry for a forest the engine's calibration map does not set,
## from its `prior`: the statement, then `sd`, the forest's size in force as one
## number in the response's units, and `sd.stated`, whether the statement is an
## sd. A fixed k gives k.scale over k, and a stated sd the number as stated; a
## drawn k, or an sd under invchi(), has no one number, and sd is NA as a
## sharer's is. The constant leaf's size is the sum of its trees', the forest's.
lonePriorRecord <- function(prior, raw, sdNamed, drawn) {
  sd <- if (drawn) {
    NA_real_
  } else if (sdNamed) {
    raw[[1L, "named.sd"]]
  } else {
    raw[[1L, "prior.scale"]] / raw[[1L, "k"]]
  }
  c(prior[1L], list(sd = sd, sd.stated = sdNamed), prior[-1L])
}

## The reader's entry for a forest that takes the fitting function's leaf
## prior, from the one-forest entry `prior` that prior gives: the leaf prior as
## stated, then sd, its part of the total in force where that is one number (a
## fixed k) and NA where k is drawn, sd.stated whether that prior is named by
## an sd the user wrote (`sdStated`), multiplier "none", its
## held coefficient as fixed() and, in prior.sd.of, that the size is a share
## where other forests take the rest. k.scale over k is its sd, k being the
## group's.
sharerLeafPrior <- function(prior, raw, sdStated) {
  drawn <- raw[[1L, "k.has.hyperprior"]] != 0
  statement <- list(
    leaf.prior = prior$leaf.prior,
    sd = if (drawn) NA_real_ else raw[[1L, "map.sd"]],
    sd.stated = sdStated,
    multiplier = "none",
    amplitude = fixed()
  )
  prior <- c(statement, prior[-1L])
  prior$prior.sd.of <- if (raw[[1L, "share"]] < 1) {
    "forest, shared"
  } else {
    "forest"
  }
  prior
}

## The 1-based positions of a sampler's forests that take the fitting
## function's leaf prior as a share of its total, empty where none does.
samplerSharers <- function(sampler, ptr = sampler$getPointer()) {
  numForests <- bartcoreNumForests(ptr)
  if (numForests < 2L || samplerCarriesCounts(sampler)) {
    return(integer(0L))
  }
  which(vapply(
    seq_len(numForests) - 1L,
    function(index) {
      raw <- .Call(C_dbarts_bartcore_getLeafPrior, ptr, index)
      !is.nan(raw[1L, "share"])
    },
    NA
  ))
}

## Installs a restated leaf prior. A named sd under the hyperprior in force
## is the engine's own named-sd write, which skips the sd it already holds; a
## new hyperprior, or a return to a prior written with k, goes through the
## model, whose install re-pins a fixed sigma, so that is put back. Neither
## moves a drawn k: k.scale is the data's under every prior, so the spread in
## force is kept until the next draw.
writeLeafPrior <- function(sampler, ptr, newModel) {
  oldModel <- sampler$model
  sameHyperprior <- identical(
    newModel@leaf.hyperprior,
    oldModel@leaf.hyperprior
  )
  if (sameHyperprior && !is.na(newModel@prior.scale)) {
    .Call(C_dbarts_bartcore_setLeafPrior, ptr, 0L, namedLeafSd(newModel))
    return(invisible(NULL))
  }
  if (sameHyperprior && is.na(oldModel@prior.scale)) {
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

## A model of several forests none of which takes the fitting function's leaf
## prior: the prior a creation with no leaf prior given writes there, k fixed at
## 2 whatever the family.
modelWithLeafPriorUnclaimed <- function(model) {
  model@leaf.prior@k <- NULL
  model@leaf.prior@prior.sd <- NULL
  model@leaf.hyperprior <- newValidated("dbartsFixedHyperprior", k = 2.0)
  model@prior.scale <- NA_real_
  model
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
    # the data the sampler runs on, always the block-holding kind, read-only
    # from outside: a write would change the R-side copy and not the engine,
    # so it stops and names the setter; the package writes dataHeld
    data = function(value) {
      if (!missing(value)) {
        stop(
          "the sampler's data are read-only here; change them with ",
          "s$setData(o2)",
          call. = FALSE
        )
      }
      dataHeld
    },
    dataHeld = "dbartsData",
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
    initialize = function(
      control,
      model,
      data,
      ...,
      active = NULL,
      periods = NULL
    ) {
      if (!inherits(control, "dbartsControl")) {
        stop("'control' must inherit from dbartsControl")
      }
      if (!inherits(model, "dbartsModel")) {
        stop("'model' must inherit from dbartsModel")
      }
      if (!inherits(data, "dbartsData")) {
        stop("'data' must inherit from dbartsData")
      }
      data <- resolveDataObject(data)
      # what holds of a sampler's data for its whole life, a copy coming
      # through here; never the class's validity, which refuses an object
      # saved before a slot existed
      checkBasisRecord(data)
      # a model handed in directly is held to what the fitting functions hold
      # the one they resolve to
      refuseNoSplittableColumn(
        model@tree.prior@splitProbabilities,
        attr(model, "forest.columns", exact = TRUE)
      )
      # the values a sampler derives from its data are data: the record the
      # data object carries is used and an absent one derived, and none rides
      # the model, where an object saved by an earlier build has the range
      data <- checkDerivedRecord(data, model)
      attr(model, "response.range") <- NULL
      # the model and data hold what the sampler runs, the tree count and
      # the quantile rule among them, and the control mirrors them
      taken <- takeControlReadings(control, model, data)
      model <- taken$model
      data <- taken$data
      # a model made by hand states the family's own reading as NA, which
      # resolves here as a fitting function resolves it
      if (is.na(model@weights.type)) {
        family <- bridgeFamily(control, model)
        model@weights.type <- familyWeightsType(
          NA_character_,
          if (identical(family, "")) "gaussian" else family
        )
      }
      handedControl <- control
      control <- mirrorControl(control, model, data)
      .self$control <- control
      .self$model <- model
      .self$dataHeld <- data
      .self$forestWeights <- list()
      .self$activeRows <- NULL

      # "auto" (a hand-built model) keeps the bridge's own dispatch; the mask
      # handed in, installed below and composed with the data's weights,
      # gives the rows a numeric column's scale is read from, and a hazard
      # fit's subjects in, which the engine cannot count, whether too few are
      # in
      units <- likelihoodUnits(.self$model, .self$data, active)
      created <- createEngine(
        .self$control,
        .self$model,
        .self$data,
        active,
        tooFewUnitsIn(units)
      )
      .self$pointer <- created$pointer
      # what the engine holds is the record from here on, written on the
      # sampler's own data object and never on the caller's
      recordDerived(.self, .self$pointer, created$records)
      if (length(model@forests) >= 2L) {
        recordBasisScales(.self, .self$pointer, active)
      }
      # after the call, so a refused creation never spends the warning's key;
      # re-creation from state (getPointer, setState) does not come through
      # here, and copy holds it off
      warnZeroTrials(dataCounts(data), identical(model@family, "logistic"))
      warnDirectDoorNotUsed(handedControl, model, data)
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
            storedState(pointer, data)
          } else {
            NULL
          }
        },
        eval.env = as.environment(.self),
        assign.env = as.environment(.self)
      )
      # the mask handed in, which the creation read the rows of a numeric
      # column's scale from, composed with the data's weights. No store: a fresh
      # sampler's state stays the promise that captures the state current at
      # first read.
      if (!is.null(active)) {
        installCreationMask(.self, active)
      }
      # once, at a creation: a copy is not one; verbose asks for the detail
      if (.self$control@verbose && !isTRUE(forestPrintState[["held"]])) {
        writeLines(samplerForestLines(.self, full = TRUE))
      }
      # a leaf covariate's center and scale read from every row, too few rows
      # being in the likelihood: under a family that reads the response, whose
      # scale was then read from every row too, one text names both
      fromAllRows <- leafFromAllRows(.self, created$records)
      if (length(fromAllRows) > 0L || !is.null(periods)) {
        warnRowsFallBack(
          units$n,
          units$numIn,
          model@family %in% c("gaussian", "nbinom"),
          fromAllRows,
          periods,
          byMask = model@family %in% c("logistic", "aft"),
          subjects = units$subjects,
          counted = fallbackCountsRows(model, data@weights, dataCounts(data)),
          countedRange = fallbackCountsRange(model, data@weights)
        )
      }

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
      "Creates a deep or shallow copy of the sampler, keeping its model, the response mapping included, and its data object, whose record of the cut points, their weights and the response range it is built on. A sampler whose engine exists in this R session gives the copy its current chain, even one stored since its last storeState(), and the copy's stored state is that chain. Otherwise, as after a save and load, the copy is given the stored state, which goes in as setState does with forceUpdate = TRUE: a split off the recorded cut points moves to the nearest point by value, and a tree that no longer fits the sampler is repaired, silently; see setState. Without a stored state such a copy is refused, as a reload is."
      shallow <- checkFlag(shallow, "shallow")
      live <- .Call(C_dbarts_bartcore_isValidPointer, pointer)
      if (!live && is.null(state)) {
        stopWithoutStoredState()
      }
      # a copy introduces no rows, so it does not raise the creation warning,
      # and it is not a creation, so it does not print the forests' block
      dupe <- withoutZeroTrialsWarning(withoutForestBlock(
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
      ))

      # a copy is a re-creation: it keeps this sampler's model and builds on
      # its data object's record, the range included, so the chains are in
      # the transform in force before any state goes in. A live engine's
      # chain is read now, drawing nothing; otherwise the stored state, which
      # is opaque and never mutated in place (storeState replaces it whole),
      # so the copy can install the same object
      installing <- if (live) {
        storedState(pointer, data)
      } else {
        state
      }
      # forestWeights is a plain list field: assigning it shares the
      # underlying object, but R's copy-on-modify means a later
      # setForestWeights on either sampler duplicates before mutating. The
      # mask rides neither the state nor the data object, so the copy takes
      # both mirrors before the install, which puts them back after the state.
      dupe$forestWeights <- forestWeights
      dupe$activeRows <- activeRows
      if (!is.null(installing)) {
        # as a reload installs it, not as a live setState: the state can hold
        # a saved-tree store of another capacity than this control names
        # (n.samples edited after the state was stored), and the copy takes
        # the state's
        refuseLegacyState(installing)
        status <- installStateOnto(
          dupe,
          dupe$pointer,
          installing,
          rawPredictorMatrix(data@x),
          TRUE,
          TRUE
        )
        # the copy's field is the chain it was made from, or, where that went
        # onto other cut points than it records, what the engine then holds,
        # so that a reload of the copy rebuilds what the copy was
        dupe$state <- if (status == 2L) {
          storedState(dupe$pointer, dupe$data)
        } else {
          installing
        }
      } else {
        dupe$reapplyForestWeights(dupe$pointer)
        dupe$reapplyActiveRows(dupe$pointer)
      }
      dupe
    },
    show = function() {
      "Pretty prints the object: its call and, for a sampler of several forests, a short block of what each forest resolved to."

      cat("dbarts sampler\n")
      cat("  call: ")
      writeLines(deparse(control@call))
      writeLines(samplerForestLines(.self))
      cat("\n")

      invisible(NULL)
    },
    predict = function(x.test, offset.test, n.threads = control@n.threads) {
      "Using existing sampler to predict for new data without re-running. n.threads is a per-call worker count that does not persist, defaulting to the sampler's own: the replay is partitioned by (chain, saved draw), each partition writing its own rows, so the answer is identical bit for bit at every value."
      refuseMatrixOfCodes(x.test, data@x, "'x.test'")
      x.test <- validateXTest(
        x.test,
        data@x,
        missing.seen = dataMissingSeen(data)
      )
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
      refuseMatrixOfCodes(x.test, data@x, "'x.test'")
      x.test <- validateXTest(
        x.test,
        data@x,
        missing.seen = dataMissingSeen(data)
      )
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
      "Sets the control object for the sampler to a new one. Preserves the call() slot; the tree count, the quantile rule and the binary flag are the sampler's, which the control mirrors."
      if (!inherits(newControl, "dbartsControl")) {
        stop("'control' must inherit from dbartsControl")
      }
      # a slot written after the constructor is not checked by it
      refusal <- perturbCutsRefusal(newControl@n.perturb.cuts)
      if (!is.null(refusal)) {
        stop(refusal, call. = FALSE)
      }

      selfEnv <- parent.env(environment())

      newControl@binary <- control@binary
      newControl@call <- control@call

      # the tree count, the quantile rule and the cut count are the sampler's:
      # a copy that differs from them is taken, the sampler's kept, the caller
      # told once the control is in; a value the caller wrote that differs is
      # refused
      held <- controlCopiesHeld(
        newControl,
        data@x,
        !is.null(model@hazard.periods)
      )
      takes <- heldCopies(model, data)
      for (name in names(held$written)) {
        if (!sameRecycled(held$written[[name]], takes[[name]])) {
          stop(
            "changing '",
            name,
            "' is not available on an existing sampler",
            call. = FALSE
          )
        }
      }

      # settings fixed at creation: the generators, anything shaping the cut
      # grid, and the four engine limits, which the sampler reads once when it
      # is created and never again. Accepting one here would leave the stored
      # control disagreeing with the engine, and a re-creation from that stored
      # control would then move the draws under categoricalExhaustiveCap.
      # The move settings are NOT among them: a sweep reads each forest's
      # mixture, window and tree-shift mode afresh, so a change is installed
      # (below) on every forest that states none of its own.
      for (slotName in c(
        "n.chains",
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
            slotName,
            "' is not available on an existing sampler"
          )
        }
      }
      if (newControl@keepTrees && is.na(newControl@n.samples)) {
        stop("keepTrees requires 'n.samples' to be specified")
      }
      # a fresh control states no draw count, and the sampler keeps its own,
      # as it keeps the copies above
      if (is.na(newControl@n.samples)) {
        newControl@n.samples <- control@n.samples
      }
      newControl <- mirrorControl(newControl, model, data)

      # a monotone sampler is birth/death-only, as creation forces: a defaulted
      # mixture is rewritten and one proposing other moves refused, before
      # anything is installed
      if (!is.null(attr(model, "monotone"))) {
        newControl@proposal.probs <- monotoneProposalProbs(
          newControl@proposal.probs,
          allowBirthDeath = TRUE
        )
      }

      refusal <- setControlMovesRefusal(newControl, control, model)
      if (!is.null(refusal)) {
        stop(refusal, call. = FALSE)
      }

      # the perturb window rides the mixture's install
      mixtureMoved <- !identical(
        newControl@proposal.probs,
        control@proposal.probs
      ) ||
        !identical(newControl@n.perturb.cuts, control@n.perturb.cuts)
      shiftMoved <- !identical(newControl@levelGibbs, control@levelGibbs)

      ptr <- getPointer()
      oldControl <- control
      # the engine takes the control first: a refusal there leaves the stored
      # control the one the engine still has
      .Call(C_dbarts_bartcore_setControl, ptr, newControl)
      selfEnv$control <- newControl
      # A changed mixture or window goes through the prior install where the
      # sampler's $setModel takes one, which writes the mode too. Everything
      # else takes the move install, which writes no prior: any change where
      # $setModel refuses, and a changed treeShift alone, since the prior
      # install would also put a fixed sigma back. Either way a refusal rolls
      # the whole control back: the stored one must never name settings the
      # engine does not have.
      if (mixtureMoved || shiftMoved) {
        tryCatch(
          if (mixtureMoved && !setModelRefused(.self)) {
            .self$setModel(model)
          } else {
            bartcoreSetMoveSettings(ptr, newControl, model)
          },
          error = function(e) {
            selfEnv$control <- oldControl
            .Call(C_dbarts_bartcore_setControl, ptr, oldControl)
            stop(e)
          }
        )
      }
      # once the control is in, so a refusal never spends the message
      warnWrittenCopies(names(held$written))
      warnNotUsed(held$carried, takes, notUsedWords("sampler"))

      invisible(NULL)
    },
    setModel = function(newModel) {
      "Sets the model object for the sampler to a new one, changing the numbers inside the model and refusing, by name and before anything is stored, a model of another structure. It takes power and base, the split probabilities, the leaf scale as k or sd (and a k drawn or held), the residual prior and a held residual sd, under the ordinary tree prior or DART alike; a split probability of zero on a column some chain's trees split on is refused, naming the column. It refuses a model of another response family (a Student-t error against a Gaussian one included), another kind of leaf or leaf covariate columns, a gp leaf's max.leaf.size, other interaction limits or blocks, DART against the ordinary tree prior, other columns a forest may split on, other monotone directions, a residual prior other than its own on a family with no residual scale (probit, logistic, ordinal, negative binomial, hazard), a chisq() residual prior of df 2 or below on a sampler whose variance forest has two or more trees, and a raw leaf specification whose k or sd does not resolve to the leaf hyperprior and prior scale the model carries; $setLeafPrior restates a leaf prior. It takes the Student-t degrees of freedom, a negative-binomial shape, a gp leaf's lengthscale and DART's own settings; a df, shape or DART concentration handed from held to drawn starts at the value in force, and a new lengthscale is refused on a sampler holding kept draws. A model stating no family, df or lengthscale keeps the sampler's. On a sampler of several forests or a multinomial one it takes the fitting function's level: a changed power or base reaches the forests with no basis that state none of their own, and the one such forest, or a multinomial model's category forests, those forests whatever they state; the leaf prior goes where setLeafPrior sends it; and the residual prior. A number a forest's record in the model states, a mean forest's own sd, power or base or the variance forest's power or base, is installed where it differs from the sampler's, through the code setForest or setVarianceForest runs on that forest, the later write winning; one the record does not state, as a variance forest's depth copied from the mean forest's at creation, is left as it is. A coefficient's prior that differs is refused naming setForest, and split probabilities other than equal by name. The sampler's own model, as after setResponse or setOffset, is always accepted. A new negative-binomial shape keeps the latents in force. This draws nothing; the next sweep draws the augmentation variables against the change."
      if (!inherits(newModel, "dbartsModel")) {
        stop("'model' must inherit from dbartsModel")
      }
      modelStructureRefusal(.self, newModel)
      # the reading is fixed at creation; a model stating none states the
      # sampler's
      if (
        !is.na(newModel@weights.type) &&
          !identical(
            newModel@weights.type,
            familyWeightsType(model@weights.type, model@family)
          )
      ) {
        stop(
          "$setModel cannot change weights.type, which is fixed when the ",
          "sampler is created; create a new sampler",
          call. = FALSE
        )
      }
      refuseInvalidLeafPrior(newModel@leaf.prior)
      # the variance forest recalibrates to the new residual prior
      if (!is.null(model@variance)) {
        refuseVarianceCalibration(newModel@resid.prior, model@variance@n.trees)
      }
      given <- newModel
      # the rest of the structure the model's records name, before anything
      # is stored
      refusal <- movedStructureRefusal(newModel, model)
      if (!is.null(refusal)) {
        stop(refusal, call. = FALSE)
      }
      # the variance forest's depth its record states, through
      # $setVarianceForest's code once the rest is installed; an unstated one,
      # a copy of the mean depth, is left
      varianceDepth <- statedVarianceDepth(.self, given)
      # split probabilities are a parameter, held to what creation holds them
      # to on a restricted forest
      refuseNoSplittableColumn(
        newModel@tree.prior@splitProbabilities,
        attr(model, "forest.columns", exact = TRUE)
      )
      ptr <- getPointer()
      selfEnv <- parent.env(environment())

      # a model stating no family, Student-t df or gp lengthscale states the
      # sampler's own, which a copy re-creates from
      if (!modelStatesFamily(newModel)) {
        attr(newModel, "family.spec") <- attr(
          model,
          "family.spec",
          exact = TRUE
        )
      }
      if (is.null(attr(newModel, "resid.df", exact = TRUE))) {
        attr(newModel, "resid.df") <- attr(model, "resid.df", exact = TRUE)
      }
      if (is(newModel@leaf.prior, "dbartsGPPrior")) {
        lengthscale <- newModel@leaf.prior@lengthscale
        newModel@leaf.prior@lengthscale <- if (is.null(lengthscale)) {
          model@leaf.prior@lengthscale
        } else {
          rep_len(as.double(lengthscale), length(newModel@leaf.prior@columns))
        }
      }
      newModel@family <- model@family
      newModel@weights.type <- model@weights.type
      # what the sampler holds for its life stays the sampler's
      newModel <- keepModelRecords(newModel, model)
      # a negative-binomial shape the model's family states is a parameter,
      # recorded on the model, which a copy reads
      familySpec <- attr(newModel, "family.spec", exact = TRUE)
      if (identical(model@family, "nbinom") && !is.null(familySpec)) {
        newModel@shape <- resolveShape(
          familySetting(familySpec, "shape", NA_real_)
        )
      }
      # the transform the prior is anchored to is data and stays the
      # sampler's, recorded on its data object
      attr(newModel, "response.range") <- NULL
      # several forests take the fitting function's level by its own route,
      # which stores only once the engine holds it
      if (bartcoreNumForests(ptr) > 1L) {
        selfEnv$model <- setModelOnForests(.self, newModel, given)
        return(invisible(NULL))
      }
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
      if (length(varianceDepth) > 0L) {
        selfEnv$model <- modelWithVarianceDepth(selfEnv$model, varianceDepth)
        installForestDepth(ptr, -1L, selfEnv$model@variance)
      }

      invisible(NULL)
    },
    setData = function(
      newData,
      updateState = NULL,
      updateCutPoints = NULL,
      updateScale = NULL,
      updateStandardization = NULL
    ) {
      "Sets the data object for the sampler to a new one; the data of a sampler of several forests are taken as rows, the bases and their records dropped, as dbarts() reads them. A hazard sampler also takes data built from a Surv response, spreading its subjects onto the sampler's own periods as dbarts() does at creation; a time past the last period leaves its subject at risk through every period with no event. A hazard sampler refuses a hand-built person-period matrix whose last column is not named period, at setData and at creation; give it survival data, or fit the rows with family = \"binomial\", the period a predictor. A cut count or quantile rule the object states wins, and the sampler derives by it from then on; one carrying cutPoints brings the count and rule that made them, an unstated rule reading as even; one stating neither keeps the sampler's, as updateCutPoints = \"none\" does whatever it states. Preserves the sigma slot and the level tables: columns are matched by name, the names being the sampler's in any order and any other refused, as predict refuses them, while a matrix with no column names is read by position with no message, as is any data handed to a sampler with none, whose record then takes the data's names, which predict needs; a factor is read by label against the levels the sampler was created with, so fewer levels, another order and declared levels no row holds are taken, while a row at a label the sampler never had, or a factor column of bare numbers, is refused. A factor whose labels were changed on purpose is relabeled to the sampler's labels with levels<- first, or given to a new sampler. The cut points, the response range and a linear or gp leaf's covariate standardization each follow updateCutPoints, updateScale and updateStandardization where given: \"position\", \"value\" or TRUE derive the value from the new data, and \"none\" or FALSE keep the sampler's. Left NULL, a value the data object records (its cutPoints, response.range and standardization slots, which the sampler writes into its own data) is used, and one it does not record is derived, so setData(sampler$data) keeps every derived value and a fresh dbartsData() derives them. A derived standardization is read over the new data's rows in the likelihood, over every row where fewer than two are in (subjects, for a hazard fit, whose rows must hold each subject's together at periods 1, 2, 3 and on, as dbarts() expands them) or the response's scale is read from every row, and refused where two or more rows in hold one value of a leaf covariate, or values that differ by rounding only, unless a gp leaf states its lengthscale, which takes the placeholder; FALSE keeps the sampler's. A probit, ordinal or nbinom sampler's 0/1 weights ride the data object, so setData(sampler$data) keeps the rows of weight 0 out; the active-row mask set by setActiveRows is cleared, the rows possibly changing. Saved draws are rewritten so that they predict what they did. updateState follows control@updateState: NULL, its default, resolves to the control's setting, and an explicit TRUE or FALSE overrides it - the same rule run() applies. It keeps the trees, the ordinal cutpoints, the Student-t degrees of freedom and the negative-binomial shape, and starts the latents as creation does (an ordinal or negative-binomial sampler from the kept cutpoints or shape). This draws nothing; the next sweep draws the augmentation variables against the change."
      updateState <- checkUpdateState(updateState)
      newData <- resolveDataObject(newData)
      # a missing test weight is refused here as at creation
      if (inherits(newData, "dbartsData")) {
        refuseMissingTestWeights(newData@weights.test)
      }
      # a binomial sampler's new data carry its successes and failures, the
      # number of rows free
      if (!samplerHoldsTrials(.self)) {
        refuseCountsMutation(
          .self,
          "$setData",
          "its K category forests fix their data at creation, and the count ",
          "response is not a column of it; make a new sampler instead"
        )
      }
      fellBack <- bartcoreSamplerSetData(
        .self,
        newData,
        updateCutPoints,
        updateScale,
        updateStandardization
      )
      if (resolveUpdateState(updateState, control)) {
        storeState()
      }
      # raised last, as warnScaleKept is: a warning turned into an error
      # leaves the record and the stored state the new data's
      units <- likelihoodUnits(model, data, activeRows)
      warnRowsFallBack(
        units$n,
        units$numIn,
        fellBack$range,
        fellBack$leaf,
        subjects = units$subjects,
        counted = fallbackCountsRows(model, data@weights, dataCounts(data)),
        countedRange = fallbackCountsRange(model, data@weights)
      )
      invisible(NULL)
    },
    setResponse = function(
      y,
      updateScale = FALSE,
      updateState = NULL,
      status = NULL
    ) {
      "Changes the response against which the sampler is fitted, and, for an aft (survival) sampler given a non-null status, its censoring structure in the same call. updateState follows control@updateState; see setData. This draws nothing; the next sweep draws the augmentation variables against the change. An aft sampler leaves a censored row at its new bound, a Student-t sampler keeps its degrees of freedom and scales, and an ordinal sampler keeps its cutpoints. A count response is written as its formula writes it: a binomial sampler takes an n x 2 matrix of successes and failures, or a vector of successes against each row's trials as they stand, a count above its row's trials refused; a multinomial sampler takes the n x K count matrix, and refuses a vector of numbers, which cannot say how each row's trials split among the categories. A factor, or a character vector read as one, is taken by a multinomial sampler whose rows each hold one trial, its labels matched by name to the categories the sampler records in any order (a label it lacks is refused by name, a factor is refused where any row holds more than one trial, naming the count matrix, and the response is then recorded as written as a factor, so draws come back as categories), and by a probit or logistic sampler whose response was written as a two-level factor, its labels matched by name to the two recorded levels, the first coded 0 and the second 1 whatever the new factor's level order (a label the sampler lacks is refused by name; a sampler whose response was written as numbers takes only 0 and 1). A binomial or multinomial sampler has the number of rows fixed, keeps its trees, reads no response scale and mirrors the counts into data@counts and their row sums into data@y; a probit or logistic sampler holds no counts and mirrors none."
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
      updateScale <- checkUpdateScale(updateScale)
      # a factor (or a character vector read as one) is read against the
      # levels the sampler recorded: a multinomial sampler's one-trial rows
      # as categories, a binary sampler's two as 0 and 1
      asFactor <- FALSE
      if (is.factor(y) || is.character(y)) {
        if (samplerCarriesCounts(.self)) {
          y <- factorAsCounts(.self, y)
          asFactor <- TRUE
        } else if (
          !samplerHoldsTrials(.self) &&
            model@family %in% c("probit", "logistic")
        ) {
          y <- factorAsBinaryCodes(.self, y)
        }
      }
      # a count response as its formula writes it: a binomial sampler's
      # successes and failures, or its successes against the trials held, and
      # a multinomial sampler's n x K counts. Neither reads a response scale.
      if (samplerCarriesCounts(.self) || samplerHoldsTrials(.self)) {
        if (!is.null(status)) {
          stop(
            "$setResponse: 'status' is an aft (survival) sampler's censoring ",
            "structure; a count sampler takes none"
          )
        }
        if (samplerHoldsTrials(.self)) {
          ptr <- bartcoreSamplerSetBinomialResponse(.self, y)
        } else {
          ptr <- bartcoreSamplerSetCounts(.self, y, asFactor)
        }
        if (resolveUpdateState(updateState, control)) {
          storeState(ptr)
        }
        return(invisible(NULL))
      }
      scaleKept <- bartcoreSamplerSetResponse(.self, y, updateScale, status)
      if (isTRUE(updateScale)) {
        recordResponseRange(.self, getPointer())
      }
      if (resolveUpdateState(updateState, control)) {
        storeState()
      }
      warnScaleKept(scaleKept)
      invisible(NULL)
    },
    setOffset = function(offset, updateScale = FALSE, updateState = NULL) {
      "Changes the offset slot used to adjust the response. updateState follows control@updateState; see setData. This draws nothing; the next sweep draws the augmentation variables against the change. A logistic sampler included."
      updateState <- checkUpdateState(updateState)
      refuseCountsMutation(
        .self,
        "$setOffset",
        "a common per-observation shift is the softmax's own null direction, ",
        "so a flat offset is inert; the shift that is not is the n x K ",
        "matrix $setCategoryOffset takes"
      )
      updateScale <- checkUpdateScale(updateScale)
      scaleKept <- bartcoreSamplerSetOffset(.self, offset, updateScale)
      if (isTRUE(updateScale)) {
        recordResponseRange(.self, getPointer())
      }
      if (resolveUpdateState(updateState, control)) {
        storeState()
      }
      warnScaleKept(scaleKept)
      invisible(NULL)
    },
    setWeights = function(
      weights,
      updateState = NULL,
      updateScale = FALSE,
      updateStandardization = FALSE
    ) {
      "Changes the weights with which the sampler is fitted. A row of weight 0 leaves the likelihood but stays in the fit, occupying a leaf. A probit or ordinal sampler takes only weights of 0 and 1, and an nbinom sampler whole-number counts: 0s and 1s name the rows in its data set, and are held in the data object's weights slot and composed with the active-row mask (see setActiveRows), a row being in when its weight is 1 and the mask keeps it. Under weights.type = \"frequency\", fixed at creation, a gaussian weight is a count of identical observations. A Student-t sampler keeps its scales, which the next sweep draws against the new weights. updateScale = TRUE then takes the response's scale (a count sampler's center) again from the rows the new weights leave in, as setOffset(updateScale = TRUE) does; it is refused on a sampler of several forests, and kept, with one warning per sampler, where those rows hold fewer than two distinct values. updateStandardization = TRUE then takes a linear or gp leaf's covariate center and scale, and a gp leaf's heuristic lengthscale, again from those rows, as setPredictor(updateStandardization = TRUE) does for every leaf covariate, keeping a covariate those rows hold at one value, or at values that differ by rounding only, with the same warning, and every covariate where updateScale = TRUE keeps the response's scale in the same call or, on a hazard fit, fewer than two subjects are in, the period among them, under a warning counting subjects; a gp sampler holding saved draws refuses it. Both are for burn-in: pass them TRUE while burning in and FALSE after, so every kept draw comes from one prior. updateState follows control@updateState; see setData. This draws nothing; the next sweep draws the augmentation variables against the change. A weight change is undone by setWeights and setState in either order, the sampler then continuing bit for bit as one that installed the stored state under the old weights would."
      updateState <- checkUpdateState(updateState)
      updateScale <- checkUpdateScale(updateScale)
      updateStandardization <- checkUpdateStandardization(
        updateStandardization
      )
      refuseCountsMutation(
        .self,
        "$setWeights",
        "an integer case weight is already row-wise replication in its count ",
        "response, and a non-integer one has no exact augmentation sampler"
      )
      if (updateScale) {
        refuseAmplitudeMutation(
          .self,
          "setWeights(updateScale = TRUE)",
          "every forest keeps its leaf calibration stated against the scale ",
          "fixed at creation; use updateScale = FALSE instead"
        )
      }
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
      # a latent family's weights are membership, 0 or 1, which the engine
      # composes with the active-row mask; nbinom's count identical rows
      if (identical(model@family, "nbinom")) {
        refuseNonCountFamilyWeights(weights, "nbinom")
      } else if (isMaskedWeightFamily(model@family)) {
        if (any(weights != 0 & weights != 1)) {
          stop(
            model@family,
            " models do not support case weights other than 0 and 1, which ",
            "mark rows in and out of the likelihood: a weighted ",
            "truncated-normal latent likelihood is not a coherent model"
          )
        }
      }

      ptr <- getPointer()
      selfEnv <- parent.env(environment())
      units <- rowSetterUnits(
        .self,
        weights,
        activeRows,
        updateStandardization,
        "$setWeights: "
      )

      oldWeights <- data@weights
      # all ones on a membership family is every row in, which the data
      # object records as no weights, as creation does
      selfEnv$dataHeld@weights <- if (
        isMaskedWeightFamily(model@family) && all(weights == 1)
      ) {
        NULL
      } else {
        weights
      }
      tryResult <- tryCatch(
        .Call(
          C_dbarts_bartcore_setWeights,
          ptr,
          weights,
          updateScale,
          updateStandardization,
          if (!is.null(units)) TRUE
        ),
        error = function(e) {
          selfEnv$dataHeld@weights <- oldWeights
          e
        }
      )
      if (inherits(tryResult, "error")) {
        stop(tryResult)
      }

      if (updateScale || updateStandardization) {
        recordResponseRange(.self, ptr, updateScale, updateStandardization)
      }
      if (resolveUpdateState(updateState, control)) {
        storeState(ptr)
      }
      warnRowsKept(.self, tryResult, units)
      invisible(NULL)
    },
    setCategoryOffset = function(offset, updateState = NULL) {
      "Installs, or at NULL clears, a multinomial sampler's n x K category offset: the latent becomes f_ik + o_ik, so the offset enters the log-sum-exp margins, every category's working response and the reported softmax probabilities, and never a leaf value. This is the response-side counterpart of setResponse's count matrix rather than of setOffset, whose flat shift is added after the categories are blended - the wrong side of the nonlinearity - and is the softmax's own null direction besides. Only the row-centered part is identified: adding a constant to a whole row leaves every reported probability unchanged, and the entrance leaves the matrix as given rather than re-centering it. It shifts the TRAIN latent only; the test rows are other rows and carry their own (setCategoryTestOffset), and predict takes its own matrix per call. Mirrored into data@offset.category, so a re-created sampler carries it. updateState follows control@updateState; see setData."
      updateState <- checkUpdateState(updateState)
      requireCountsCapability(.self, "$setCategoryOffset")
      ptr <- bartcoreSamplerSetCategoryOffset(.self, offset)
      if (resolveUpdateState(updateState, control)) {
        storeState(ptr)
      }
      invisible(NULL)
    },
    setCategoryTestOffset = function(offset.test) {
      "Installs, or at NULL clears, a multinomial sampler's nTest x K category test offset: the recorded test channel becomes softmax(f_test + o_test), formed where the train blend forms softmax(f + o). The test fits enter no likelihood, so this moves the reported test probabilities and nothing else - no draw, no working response, no train channel. Its rows are the CURRENT test rows, so replacing those rows while it is installed is refused rather than silently reinterpreted; clear it first. Out-of-sample predict does not read it at all, taking its own matrix for the rows it is given. Mirrored into data@offset.category.test, so a re-created sampler carries it."
      requireCountsCapability(.self, "$setCategoryTestOffset")
      bartcoreSamplerSetCategoryTestOffset(.self, offset.test)
      invisible(NULL)
    },
    setActiveRows = function(
      active,
      updateState = NULL,
      updateScale = FALSE,
      updateStandardization = FALSE
    ) {
      "Sets the per-observation 0/1 mask of rows in the data set for this sampler. An inactive row leaves every sufficient statistic, every family-level parameter update and its own latent draw, but keeps its leaf occupancy and its fitted value. The mask composes with the weights in either order: a row is in when its weight is positive and the mask keeps it, so on a probit, ordinal or nbinom sampler, whose 0/1 weights ride the data object, activeRows is only what this installed. Each sweep draws the latent of every row the mask holds in, so a row switched back in is drawn against the current fit with the rest. updateScale = TRUE then takes the response's scale (a count sampler's center) again from the rows the mask leaves in, as setOffset(updateScale = TRUE) does; a multinomial sampler takes it and derives nothing, its leaf scale and center being the same whatever rows are in, as a probit or logistic sampler does; it is refused on a sampler of several forests, and kept, with one warning per sampler, where those rows hold fewer than two distinct values. updateStandardization = TRUE then takes a linear or gp leaf's covariate standardization again from those rows, as setWeights does. NULL clears, and an all-ones mask installs nothing. The mask does not ride the saved state; it is mirrored on an R5 field that getPointer, setState and copy reinstall on every re-creation. updateState follows control@updateState; see setData. This draws nothing; the next sweep draws the augmentation variables against the change."
      updateState <- checkUpdateState(updateState)
      updateScale <- checkUpdateScale(updateScale)
      updateStandardization <- checkUpdateStandardization(
        updateStandardization
      )
      if (updateScale) {
        refuseAmplitudeMutation(
          .self,
          "setActiveRows(updateScale = TRUE)",
          "every forest keeps its leaf calibration stated against the scale ",
          "fixed at creation; use updateScale = FALSE instead"
        )
      }
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

      units <- rowSetterUnits(
        .self,
        data@weights,
        active,
        updateStandardization,
        "$setActiveRows: "
      )
      kept <- bartcoreSamplerSetActiveRows(
        .self,
        active,
        updateScale,
        updateStandardization,
        units
      )
      if (resolveUpdateState(updateState, control)) {
        storeState()
      }
      warnRowsKept(.self, kept, units)
      invisible(NULL)
    },
    setForestWeights = function(forest, weights, updateState = NULL) {
      "Sets a per-forest, per-observation weight: a multiplicative precision factor on the named forest's own leaf conditionals, composing with weights and active as (w_i * a_i) * m_f^2 * s_i rather than widening either channel. Only applies to a Bayesian causal forest built with forests = (see dbarts); forest is a number, the forest's position from 1, or a string, its label (forest<i> names position i), as with getLeafPrior/getK (the basis forest is 2). The weight does not ride the sampler's saved state; it is mirrored on an R5 field that getPointer and setState both reinstall on every re-creation. updateState follows control@updateState; see setData."
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

      index <- samplerForestIndex(forest, model, data, getPointer())
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
    setForestBasis = function(
      forest,
      basis,
      updateState = NULL,
      updateBasisScale = FALSE
    ) {
      "Changes the basis the named forest's amplitudes multiply. forest is a number, the forest's position from 1, or a string, its label (forest<i> names position i), as with setForestWeights and getLeafPrior/getK (a Bayesian causal forest's basis forest is 2). basis is a value, or a one-sided formula, which is read as the basis of a forest() is, its names found where the formula was written. It must be of the class the forest was created with, which decides what the forest multiplies and never changes. A forest of levels (created with a factor, a character vector or a logical vector) takes one of those three, matched to the forest's levels by name and not by position: the replacement may list its levels in another order or leave one out, a level of the forest that no row has is a column of zeros (a swap can leave one momentarily unobserved), and a row at a level the forest does not have, or at NA, is refused. Coefficient j multiplies level j's rows after any swap. Numbers are refused for such a forest, the level indicators the sampler stores among them; hand back the factor. A forest of numbers takes one numeric column, a vector or a matrix of one column, a column of all zeros among them (creation refuses one), and refuses several columns, a factor, a character vector and a logical vector; a forest created with no basis takes none. The sampler keeps the forest's block, in data@bases, and its levels, in data@basis.levels, each read by the forest's position (neither list has names, so a label does not index them), and not the factor it was handed: a caller that will hand a basis back keeps the factor it handed over. A state records the levels of each forest a factor multiplies, and setState puts each of its coefficients on this sampler's level of the same name, refusing a state whose names differ. No forest changes width, so the names the basis was created with stay, whatever names the replacement has. The forest's label does not change. This is the SOLE route by which a basis changes after creation, and the amplitudes are preserved: an install leaves every one of them bitwise. A swap moves the basis and nothing else: no forest's prior moves, so a numeric forest keeps the sd per unit its creation gave it, a default one the scale of the column it was created on (data@basis.scale). updateBasisScale = TRUE, for burn-in, also takes a numeric forest's default scale again from the column installed, over the rows in the likelihood at the call (a positive weight, kept by the active-row mask; a forest's own weights are not read), which may differ from the rows the response's scale was read from at creation, and restates that forest's default sd from it from the next sweep, leaving the trees and the coefficients as they are; the new scale replaces the forest's entry of data@basis.scale, so a copy or a reload runs the new prior. Re-installing the column in force with TRUE rescales it after the weights or the mask change. Where fewer than two rows are in, the column is constant over them, or the weights are counts that total at most 1 over them, the column is installed and the scale in force kept, with a warning on the sampler's first keep. TRUE is refused on a forest that states its sd and on one multiplied by a factor, neither having a default taken from a column. The matrix is mirrored into data@bases as setWeights mirrors weights, so it survives the sampler's re-creation. updateState follows control@updateState; see setData. This draws nothing; the next sweep draws the augmentation variables against the change."
      updateState <- checkUpdateState(updateState)
      updateBasisScale <- checkFlag(updateBasisScale, "updateBasisScale")
      refuseCountsMutation(
        .self,
        "$setForestBasis",
        "its forests are its categories, which carry no amplitudes and so no ",
        "basis for any to multiply"
      )
      index <- samplerForestIndex(forest, model, data, getPointer())
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
      # the forest takes the kind it was created with, read from the data
      # object, and nothing is stored before the value has been held to it
      current <- data@bases[[index + 1L]]
      labels <- vapply(model@forests, methods::slot, "", "label")
      recorded <- dataBasisLevels(data)
      recorded <- if (index < length(recorded)) recorded[[index + 1L]]
      values <- conformBasis(
        replacementForestBasis(basis, length(data@y)),
        forestBasisKind(current, recorded),
        recorded,
        "basis",
        forestName(index + 1L, labels)
      )
      if (nrow(values) != length(data@y)) {
        stop("'basis' must have the same length as 'y'")
      }
      # a forest keeps its width for life, so the names recorded at creation
      # stay, whatever names the replacement has
      dimnames(values) <- if (!is.null(colnames(current))) {
        list(NULL, colnames(current))
      }

      ptr <- getPointer()
      if (updateBasisScale) {
        refuseBasisScaleUpdate(
          .Call(C_dbarts_bartcore_getLeafPrior, ptr, index),
          forestName(index + 1L, labels)
        )
      }
      selfEnv <- parent.env(environment())

      oldData <- data
      newBases <- data@bases
      newBases[[index + 1L]] <- values
      selfEnv$dataHeld@bases <- newBases
      tryResult <- tryCatch(
        {
          kept <- .Call(
            C_dbarts_bartcore_setForestBasis,
            ptr,
            index,
            values,
            updateBasisScale
          )
          if (updateBasisScale) {
            recordBasisScale(.self, ptr, index)
          }
          kept
        },
        error = function(e) {
          selfEnv$dataHeld <- oldData
          e
        }
      )
      if (inherits(tryResult, "error")) {
        stop(tryResult)
      }

      if (resolveUpdateState(updateState, control)) {
        storeState(ptr)
      }
      # raised last, as warnScaleKept is: a warning turned into an error
      # leaves the record and the stored state the update's
      if (updateBasisScale) {
        warnBasisScaleKept(
          tryResult,
          forestName(index + 1L, labels),
          length(data@y),
          sum(rowsInLikelihood(data@weights, activeRows, length(data@y)))
        )
      }
      invisible(NULL)
    },
    setSigma = function(sigma, updateState = NULL) {
      "Changes the residual standard deviation parameter for each chain; on a sampler that holds sigma fixed it rewrites the model's fixed value. updateState follows control@updateState; see setData. This draws nothing; the next sweep draws the augmentation variables against the change."
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
      updateCutPoints = "none",
      updateState = NULL,
      updateStandardization = FALSE
    ) {
      "Changes a single column of the predictor matrix, or the entire matrix if column is missing, where x is a data frame, coded by label as a column is, or a numeric matrix, taken only where every predictor column is numeric and refused by name on a sampler with a factor column; the columns are matched by name, a name that is not the sampler's being refused, and a matrix with no column names is read by position with no message. Unforced, returns TRUE when the update was made and FALSE when it was refused and the sampler left as it was; forced, returns NULL invisibly. updateCutPoints is \"none\", which keeps each cut grid, or \"position\" or \"value\", which derive the grid a new sampler would have for the values and say where the splits on it go: each on its position, rescaled when the number of points changes, or each to the point nearest its threshold. updateStandardization = TRUE derives a linear or gp leaf's covariate standardization again for the covariates among the columns the update replaces, over their new values at the rows in the likelihood, keeping one those rows hold at one value, or at values that differ by rounding only, and every one, the period among them, on a hazard fit where fewer than two subjects are in, with a warning on the sampler's first keep, and keeps the others', a linear leaf's coefficients, current and saved, rewritten so that each leaf is the function it was; FALSE, the default, keeps it. updateState follows control@updateState; see setData. The accepted predictors reach the next sweep's latents, so a Metropolis step that scores proposals with the latents integrated out is exact. This draws nothing; the next sweep draws the augmentation variables against the change."
      updateState <- checkUpdateState(updateState)

      result <- bartcoreSamplerSetPredictor(
        .self,
        x,
        column = if (missing(column)) NULL else column,
        forceUpdate = if (missing(forceUpdate)) NULL else forceUpdate,
        updateCutPoints = updateCutPoints,
        updateStandardization = updateStandardization
      )
      if (resolveUpdateState(updateState, control)) {
        storeState()
      }
      # raised last: a warning turned into an error leaves the record and the
      # stored state the update's
      warnScaleKept(!is.null(result$kept), result$kept, units = result$units)
      # a forced update's NULL comes back invisible and an unforced update's
      # verdict visible; a bare value here would always be visible
      if (result$visible) result$value else invisible(result$value)
    },
    setCutPoints = function(
      cuts,
      column,
      updateState = NULL,
      splits = "position"
    ) {
      "Changes the cut points for the predictors in column, or the entire set itself if the column argument is missing, when the entries of factor columns are not read. A grid out of order is sorted, and a point may appear only once, so for more splits near a value give a denser grid there. splits says where the splits on a column go: \"position\" keeps each on its position, rescaled when the number of points changes, and \"value\" moves each to the point nearest its threshold. Always forced, as in 0.9-x: it has no unforced form, prunes any leaves that end up empty and returns NULL invisibly. It sets a grid, typically once; to have a grid follow values that change inside a loop, use setPredictor(updateCutPoints = ). The grid is recorded in the sampler's data object, so a later setData(sampler$data) keeps it, while a later setData that derives the cut points derives at most n.cuts again, whatever grid was set. updateState follows control@updateState; see setData."
      updateState <- checkUpdateState(updateState)

      bartcoreSamplerSetCutPoints(
        .self,
        cuts,
        column = if (missing(column)) NULL else column,
        splits = splits
      )
      if (resolveUpdateState(updateState, control)) {
        storeState()
      }
      invisible(NULL)
    },
    setTestPredictor = function(x.test, column) {
      "Changes a single column of the test predictor matrix, or the whole test set if column is missing, where x.test is a data frame or a numeric matrix, taken only where every predictor column is numeric."

      bartcoreSamplerSetTestPredictor(
        .self,
        x.test,
        column = if (missing(column)) NULL else column
      )
    },
    setTestPredictorAndOffset = function(x.test, offset.test) {
      "Changes the test predictor matrix, and optionally the test offset."
      refuseMatrixOfCodes(x.test, data@x, "'x.test'")
      if (missing(offset.test)) {
        # predictors only; the engine keeps the current offset and the
        # bridge refuses if the row count would orphan its length
        return(bartcoreSamplerSetTestPredictor(.self, x.test, column = NULL))
      }

      testRowNames <- observationRowNames(x.test)
      x.test <- validateXTest(
        x.test,
        data@x,
        missing.seen = dataMissingSeen(data)
      )
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

      selfEnv$dataHeld@testUsesRegularOffset <- FALSE
      selfEnv$dataHeld@x.test <- x.test
      selfEnv$dataHeld@offset.test <- offset.test
      tryResult <- tryCatch(
        .Call(
          C_dbarts_bartcore_setTestPredictorAndOffset,
          getPointer(),
          data@x.test,
          data@offset.test
        ),
        error = function(e) {
          restored <- restoreTestSet(getPointer(), oldX.test, oldOffset.test)
          selfEnv$dataHeld@testUsesRegularOffset <- oldTestUsesRegularOffset
          selfEnv$dataHeld@x.test <- restored$x.test
          selfEnv$dataHeld@offset.test <- restored$offset.test
          if (is.null(restored$x.test)) {
            selfEnv$dataHeld <- setDataRowNames(selfEnv$dataHeld, "test", NULL)
          }
          e
        }
      )
      if (inherits(tryResult, "error")) {
        stop(tryResult)
      }
      selfEnv$dataHeld <- setDataRowNames(data, "test", testRowNames)
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
      selfEnv$dataHeld@testUsesRegularOffset <- FALSE
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
      selfEnv$dataHeld@offset.test <- offset.test
      tryResult <- tryCatch(
        .Call(C_dbarts_bartcore_setTestOffset, ptr, data@offset.test),
        error = function(e) {
          selfEnv$dataHeld@offset.test <- oldOffset.test
          e
        }
      )
      if (inherits(tryResult, "error")) {
        stop(tryResult)
      }

      invisible(NULL)
    },
    getLatents = function(result) {
      "Returns the values the last sweep's trees were drawn from. Each sweep draws them first, after any change made between sweeps and before its trees, so between a setter and the next run they are the values from before the setter (setData starts them as creation does, an ordinal or negative-binomial sampler from the kept cutpoints or shape; an aft setResponse leaves a censored row at its new bound) and the next sweep draws them afresh. Their meaning is per family and not uniform. A LOCATION, on the sampler's own latent scale, for probit (the truncated normal z), ordinal (the same z under the ordinal thresholds) and aft (the imputed log survival time): these are regressed on directly. A PRECISION, one per observation, for logistic and nbinom (the Polya-Gamma omega) and Student-t (the scale-mixing lambda): these WEIGHT a working response and are not on the response scale at all. NULL for a plain gaussian sampler and for a multinomial one, neither of which augments. Note that a gaussian-family sampler built with family = student() DOES report latents, and they are precisions."
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

      if (!is.null(model@variance)) {
        return(NULL)
      }
      ptr <- getPointer()
      .Call(C_dbarts_bartcore_getSigmas, ptr)
    },
    getShape = function() {
      "Returns the shape parameter of the sampler's family currently in force, one per chain, or NULL on a family with none (only nbinom has one today) - the count analog of getSigmas(). It is the same scalar run()$shape records once per kept draw, read mid-sweep and without serializing state, so an outer sampler driving this one a sweep at a time reads it here; the stored state holds it only where it is drawn. Under a fixed shape it repeats the value the sampler was created with; otherwise it is that sweep's grid draw."
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
      "Returns a sampler's per-forest internal-scale fitted values (a Bayesian causal forest's 1 = prognostic, 2 = treatment; an ordinary sampler's only forest is 1), n.observations x n.chains at one forest, or, at the default forest = NULL, every forest stacked with the forest margin between the observations and the chains, n.observations x n.forests x n.chains (a single-forest sampler's NULL read is bitwise its forest 1 read). forest is a number, the forest's position from 1, or a string, its label (forest<i> names position i), as with setForestWeights/setForestBasis/getLeafPrior/getK."
      ptr <- getPointer()
      if (!is.null(forest)) {
        return(.Call(
          C_dbarts_bartcore_getForestFits,
          ptr,
          samplerForestIndex(forest, model, data, getPointer())
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
      .Call(C_dbarts_bartcore_getVariance, ptr, checkFlag(test, "test"))
    },
    getForestAmplitudes = function(forest = NULL) {
      "Returns the named forest's amplitudes - the scalars its basis columns are multiplied by, one per column - as a q x n.chains matrix, or, at the default forest = NULL, every forest's stacked forest-major into a sum(q) x n.chains matrix, which is the row order the run's own glue channel carries. A drawn one has prior median size 1, on a factor the difference between two levels' amplitudes, a level's own having 0.71. The vector is RAGGED, forest by forest, which is why a forest can be named: a Bayesian causal forest's forest 1 carries the single a on its implicit intercept and its forest 2 the pair (b0, b1) on its two level indicators, so the stacked read is its shipped (a, b0, b1). forest is a number, the forest's position from 1, or a string, its label (forest<i> names position i), as with setForestBasis/setForestWeights/getLeafPrior."
      ptr <- getPointer()
      .Call(
        C_dbarts_bartcore_getForestAmplitudes,
        ptr,
        if (is.null(forest)) {
          NULL
        } else {
          samplerForestIndex(forest, model, data, getPointer())
        }
      )
    },
    getForestVariableCounts = function(forest = NULL) {
      "Returns a sampler's per-forest predictor split counts (a Bayesian causal forest's 1 = prognostic, 2 = treatment; an ordinary sampler's only forest is 1), n.predictors x n.chains at one forest, or, at the default forest = NULL, every forest stacked with the forest margin between the predictors and the chains, n.predictors x n.forests x n.chains (a single-forest sampler's NULL read is bitwise its forest 1 read). Rows are named by the predictor columns when data@x carries colnames, on margin 1 in both shapes. forest is a number, the forest's position from 1, or a string, its label (forest<i> names position i), as with setForestWeights/setForestBasis/getLeafPrior/getK."
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
          samplerForestIndex(forest, model, data, getPointer())
        )
      }
      predictorNames <- colnames(data@x)
      if (!is.null(predictorNames)) {
        rownames(counts) <- predictorNames
      }
      counts
    },
    getLeafPrior = function(forest = NULL) {
      "Returns the leaf prior a forest runs under, alone, as a named list: leaf.prior, the specification in the terms it was named in - normal(), linear() or gp() carrying one of k (a number or a chi() prior) or sd (a number or an invchi() prior), the family default when none was named - which goes back into setLeafPrior or a fitting function's leaf.prior as is; leaf.model; prior.sd.of, what the sd is the sd of ('forest', 'coefficient' or 'amplitude'); prior.mean; k.scale, the value k is relative to, so the spread in force on each chain is k.scale / getK() - the data's scale, fixed by the family and the response mapping and the same however the prior is written; response.scale and response.shift. A forest the calibration map does not set also carries, after leaf.prior, sd, its size in force in the response's units - k.scale over k at a fixed k, the number as stated where an sd was, and NA where k is drawn or the sd is under invchi(), which have no one number - and sd.stated, whether the leaf prior is named by an sd, at creation or through setLeafPrior. On a forest whose scale a multi-forest calibration map sets, k is pinned at 1 and leaf.prior is the statement of the forest's size: forest(sd = v), v the number as it was written, where an sd was stated, and forest() where none was, and its sd is the one setForest states. The list then begins with leaf.prior, sd, the forest's sd in force, in the response's units - v itself, or the default's multiple of the response's scale times the unit the sampler derived at creation (see forest) - sd.stated, TRUE where the size comes from an sd the user wrote at either level, the forest's own forest(sd = ) or, on a forest with no basis that states none, the fitting function's normal(sd = ), multiplier, what the forest multiplies ('none', 'factor' or 'numeric'), and amplitude, fixed() where the forest's coefficient is held, student(df) or, on a forest with no basis, normal() where it is drawn under the prior its creation stated, and absent where it is drawn under its default. prior.sd.of says what sd is the size of: 'forest' for a forest with no basis, 'level difference' for a factor and 'forest, per unit of basis' for a numeric column, the size exactly where the coefficient is held and at a coefficient of 1, ' (at a coefficient of 1)', where it is drawn. A numeric forest's sd is per unit of its column and carries the column's name where the basis has one, and the list adds basis.scale, where the data object records one, the standard deviation of the column its default was taken from at creation, named as sd is. k.scale is the forest's own prior sd before its coefficient: the forest's sd, on every forest but a sharer. A drawn coefficient is a multiplier with prior median size 1 - on a factor the difference between two levels' coefficients, a level's own having 0.71 - and on a numeric column through bcf's 0.674. A forest with no basis stating neither sd nor a coefficient's prior, alone or sharing one prior total, reads under the fitting function's leaf prior, as one forest does: leaf.prior that prior, sd its part of the total in force, NA where k is drawn, sd.stated TRUE when the fitting function's leaf prior is named by an sd the user wrote, normal(sd = ), as it is on a sampler of one forest, and FALSE for the default and for k, amplitude fixed() and prior.sd.of 'forest' alone and 'forest, shared' in a group, one left of a group included. Every chain runs under the sampler's one prior and response transform, which no state install moves, so every value is shared by the chains; an NA spread is refused on write. A drawn k is chain state, read by getK. At the default forest = NULL a multi-forest sampler returns a list of one prior per forest, named by the forests' labels, a multinomial sampler's by the response's levels; a single-forest sampler's NULL read is bitwise its forest 1 read. forest is a number, the forest's position from 1, or a string, its label (forest<i> names position i)."
      ptr <- getPointer()
      read <- function(index) {
        reportLeafPrior(
          .self,
          .Call(C_dbarts_bartcore_getLeafPrior, ptr, index),
          index
        )
      }
      if (!is.null(forest)) {
        return(read(samplerForestIndex(forest, model, data, getPointer())))
      }
      numForests <- bartcoreNumForests(ptr)
      if (numForests == 1L) {
        return(read(0L))
      }
      priors <- lapply(seq_len(numForests) - 1L, read)
      names(priors) <- samplerForestLabels(model, data)
      priors
    },
    getK = function(forest = NULL) {
      "Returns each chain's current leaf-prior k, the value run()$k records per draw when k is drawn in that run (a run records no k while it is fixed), read without running, as getSigmas reports sigma; after a run it is bitwise the last draw. A fixed k is reported once per chain, and a forest whose scale a multi-forest calibration map sets reports 1; a forest with no basis stating neither sd nor a coefficient's prior reports the leaf prior's k, and those that share one prior total each the group's one. It is k whatever terms the prior was named in, relative to getLeafPrior()$k.scale. A vector of length n.chains at one forest, or, at the default forest = NULL on a multi-forest sampler, an n.forests x n.chains matrix; a single-forest sampler's NULL read is bitwise its forest 1 read. forest is a number, the forest's position from 1, or a string, its label (forest<i> names position i)."
      ptr <- getPointer()
      read <- function(index) {
        .Call(C_dbarts_bartcore_getLeafPrior, ptr, index)[, "k"]
      }
      if (!is.null(forest)) {
        return(read(samplerForestIndex(forest, model, data, getPointer())))
      }
      numForests <- bartcoreNumForests(ptr)
      if (numForests == 1L) {
        return(read(0L))
      }
      do.call(rbind, lapply(seq_len(numForests) - 1L, read))
    },
    setLeafPrior = function(leaf.prior, updateState = NULL) {
      "Restates the leaf prior's spread, or the hyperprior it is drawn under, on every chain, in the vocabulary a fitting function's leaf.prior takes: normal(sd = ), normal(k = ), an invchi() prior on the sd, linear(sd = ) or gp(sd = ). The specification must name the sampler's own leaf model; leaf-model details such as a linear leaf's columns may be omitted and, if given, must match. Nothing else moves - not the tree prior, the response transform or sigma. Into a drawn prior each chain's spread in force is kept at the call, from a fixed k or sd as from a drawn one: k.scale is the data's under every prior, so k is not touched and the new prior acts from the next draw of k. A fixed k or sd stated is taken as written, a fixed sd as k = k.scale / sd. Since the write keeps the chains' spread, on a sampler not yet run it is not the start creation under the new prior would make. A multinomial sampler takes normal(k = ) with a fixed k, Inf included, applied to every category forest. On a sampler whose forests carry amplitudes it is the fitting function's leaf prior, written on the forests it reaches as creation under it gives them. Where one forest has no basis it is written there whatever that forest states, a stated sd going; on a forest stating a coefficient's prior, which has no k, it sets the size the prior gives at the scale creation set, and no k is drawn there. Where several have none, the forests sharing one prior total take it, each keeping the part of it its creation gave it, so after a forest has left the group the sizes no longer add up to it, and those stating their own sd keep it; where none shares, as where every forest has a basis or every one with no basis states its sd, the default restated changes nothing and any other leaf.prior is refused, naming setForest(forest, sd = ). Into a drawn k a forest starts at the spread in force. A write changes the model. The write takes effect on the next sweep, reinterpreting no value already drawn; a write equal to what is in force is bitwise inert, that rounding aside. The write is recorded on the model field, so a re-creation or a later re-anchoring channel restates it rather than the creation value. setModel changes everything else. updateState follows control@updateState; see setData."
      updateState <- checkUpdateState(updateState)
      multinomial <- samplerCarriesCounts(.self)
      amplitudes <- samplerCarriesAmplitudes(.self)
      if (missing(leaf.prior)) {
        stop(
          "'leaf.prior' must be given: a leaf prior specification such as ",
          "normal(sd = 1)"
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
      selfEnv <- parent.env(environment())
      # the fitting function's leaf prior on several forests, by the writer
      # $setModel shares
      if (amplitudes || multinomial) {
        written <- setLeafPriorOnForests(.self, spec, expr)
        if (!is.null(written)) {
          selfEnv$model <- written
        }
        if (resolveUpdateState(updateState, control)) {
          storeState()
        }
        return(invisible(NULL))
      }
      newModel <- restateLeafPrior(.self, spec, expr)
      ptr <- getPointer()
      writeLeafPrior(.self, ptr, newModel)
      selfEnv$model <- newModel
      if (resolveUpdateState(updateState, control)) {
        storeState(ptr)
      }
      invisible(NULL)
    },
    getForest = function(forest = NULL) {
      "Returns every value in force on one forest as a plain named list, in setForest's order and with the names its kind takes: sd, NULL where k is drawn or no one positive number is in force; power; base; amplitude, the coefficient's prior as fixed(), normal() or student(df); proposal.probs, n.perturb.cuts and treeShift, a forest stating none reading the control's. do.call(s$setForest, c(list(f), s$getForest(f))) changes nothing. forest is a number, the forest's position from 1, or a string, its label (forest<i> names position i), left out on a sampler of one forest and on a multinomial one, whose list covers every category forest, whichever is named."
      ptr <- getPointer()
      forestInForce(.self, doorForestIndex(.self, forest, ptr), ptr)
    },
    setForest = function(
      forest = NULL,
      sd = NULL,
      power = NULL,
      base = NULL,
      amplitude = NULL,
      proposal.probs = NULL,
      n.perturb.cuts = NULL,
      treeShift = NULL,
      updateState = NULL,
      n.trees = NULL,
      vars = NULL,
      interactions = NULL,
      blocks = NULL,
      basis = NULL,
      ...
    ) {
      "Changes one forest, named by its position or label, in forest()'s names: sd, power, base, amplitude, proposal.probs, n.perturb.cuts and treeShift, each meaning what it means in forest() on that forest at creation. It changes only the names given, from the next sweep; no tree, leaf value or coefficient already drawn moves. A write changes the model. A name left out or NULL changes nothing, and so does a value equal to the one in force, which leaves the forest following what it followed, so the list getForest returns goes back as is. Each name is compared with the value in force before the call, whatever else the call changes, so an sd handed back beside amplitude = fixed() is no change; to pin a forest at its sd, state a different one or, where k is drawn, the same in a second call. A forest that states its sd or its depth is changed by stating another: nothing returns it to following the fitting function. On a model of one forest, sd is the size setLeafPrior(normal(sd = )) writes, and the later write wins between this door and setLeafPrior and setModel. Elsewhere sd is the forest's trees' scale: its size where the coefficient is held and its size at a coefficient of 1 where it is drawn, on a factor the difference between two levels' and on a numeric column per unit. A forest with no basis alone under the leaf prior leaves it at an sd stated for it, its k then pinned at 1. A forest sharing one prior total with others leaves the group the same way, and those still sharing keep their parts of the total, set once at creation, and their sizes; nothing puts a forest back in the group but a new sampler. amplitude moves the coefficient between held, fixed(), and drawn, normal() or student(df): a hold sets it to its held value at the call, a release draws it from there, a Student-t's mixing variance starts at the normal's, and no leaf or trees' scale moves, except that a forest the fitting function's leaf prior sizes takes creation's size when released and the last left of a group states its part as its sd. normal() writes a basis forest back to its default, fixed() one with no basis. vars, n.trees, interactions and blocks are refused by name whatever their value, being part of the model, and a basis goes through setForestBasis. What creation refuses of a statement is refused here for the same reason, and every check comes before any write. A multinomial model's forests take one statement: forest is left out, sd is refused, amplitude takes fixed() and nothing else, and a depth reaches every category. Every change is recorded on the model, so a copy and a reload run what the sampler runs, and setModel with that model installs the forest's stated sd, power and base through this code, the later write winning. A variance forest stays as it is. updateState follows control@updateState; see setData."
      vocabulary <- forestConstructors[FOREST_ARGUMENT_VOCABULARIES$forests]
      if (!missing(sd)) {
        sd <- evalInForestVocabulary(
          substitute(sd),
          vocabulary,
          parent.frame(),
          numbers = TRUE
        )
      }
      if (!missing(amplitude)) {
        amplitude <- evalInForestVocabulary(
          substitute(amplitude),
          vocabulary,
          parent.frame()
        )
      }
      writeForestDoor(
        .self,
        forest,
        list(
          sd = sd,
          power = power,
          base = base,
          amplitude = amplitude,
          proposal.probs = proposal.probs,
          n.perturb.cuts = n.perturb.cuts,
          treeShift = treeShift
        ),
        c(
          c("n.trees", "vars", "interactions", "blocks", "basis")[c(
            !missing(n.trees),
            !missing(vars),
            !missing(interactions),
            !missing(blocks),
            !missing(basis)
          )],
          dotNames(...)
        ),
        updateState
      )
    },
    setTreePrior = function(power = NULL, base = NULL, split.probs = NULL) {
      "Changes the fitting function's tree prior on a live sampler, as power, base and split.probs mean at creation, through setModel's route, so the two cannot differ. A name left out or NULL changes nothing, and so does a value equal to the one in force. power and base reach the forests with no basis that state none of their own and, where one forest has no basis, or on a multinomial sampler's category forests, those forests whatever they state, the later write winning; they make no forest state them, and nothing returns a forest that states its own to following, short of a new sampler. A name that reaches no forest is refused, naming the forests that state it and setForest. split.probs is a vector, which changes the values, 1 for equal probabilities, or on a DART sampler dart(...), which changes DART's settings, a setting left out taking creation's default; dart() on an ordinary sampler and a vector on a DART one change the model's structure and are refused, and a sampler of several forests takes equal probabilities only. A variance forest stays as it is; setVarianceForest changes it."
      splitExpr <- if (!missing(split.probs)) substitute(split.probs)
      writeTreePriorLevel(
        .self,
        list(power = power, base = base),
        splitExpr,
        parent.frame()
      )
    },
    getVarianceForest = function() {
      "Returns every value in force on the variance forest as a plain named list: power, base, proposal.probs and n.perturb.cuts, which setVarianceForest takes back and which then change nothing; NULL on a sampler with none."
      varianceForestInForce(.self)
    },
    setVarianceForest = function(
      power = NULL,
      base = NULL,
      proposal.probs = NULL,
      n.perturb.cuts = NULL,
      updateState = NULL,
      n.trees = NULL,
      vars = NULL,
      ...
    ) {
      "Changes the variance forest of a heteroscedastic sampler: power, base, proposal.probs and n.perturb.cuts, as varianceForest() takes them, by setForest's rules. vars and n.trees are refused by name, and a sampler with no variance forest refuses the call. The variance forest takes the mean forest's depth at creation and is on its own after that; a write is recorded as stated, and setModel installs a power or base a model's variance record states through this code, the later write winning. updateState follows control@updateState; see setData."
      writeVarianceForestDoor(
        .self,
        list(
          power = power,
          base = base,
          proposal.probs = proposal.probs,
          n.perturb.cuts = n.perturb.cuts
        ),
        c(
          c("n.trees", "vars")[c(!missing(n.trees), !missing(vars))],
          dotNames(...)
        ),
        updateState
      )
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
        .Call(
          C_dbarts_bartcore_setActiveRows,
          ptr,
          activeRows,
          FALSE,
          FALSE,
          NULL
        )
      }
      invisible(NULL)
    },
    getPointer = function() {
      "Returns the underlying reference pointer, checking for consistency first."
      selfEnv <- parent.env(environment())

      if (.Call(C_dbarts_bartcore_isValidPointer, pointer) == FALSE) {
        if (is.null(state)) {
          stopWithoutStoredState()
        }
        refuseLegacyState(state)
        # the new engine holds the cut points the data object records, and
        # the stored state goes onto them, forced, a split off them moved by
        # value; data@x serves a column whose splits move (the engine keeps
        # no predictor matrix). A store sized through the flat API is in no
        # control, so the re-created sampler takes the stored state's
        # capacity
        ptr <- recreatePointer(control, model, data)
        status <- installStateOnto(
          .self,
          ptr,
          state,
          rawPredictorMatrix(data@x),
          TRUE,
          TRUE
        )
        # the replacement is bound only once it carries the state: a refused
        # install must leave the object exactly as it was rather than holding
        # a live but unfitted engine that the next run would silently sample
        # from and then store over the fitted state. The abandoned pointer
        # carries its own holder and finalizes itself.
        selfEnv$pointer <- ptr
        # a state moved onto other cut points than it records is no longer
        # what the engine holds, so the field is the engine's own
        if (status == 2L) {
          storeState(ptr)
        }
      }
      pointer
    },
    setState = function(newState, forceUpdate = FALSE) {
      "Installs a stored state: the chains, never the model, and never the cut points or their weights, which stay this sampler's. A state records the level names of each forest a factor multiplies, and each of its coefficients goes onto this sampler's level of the same name, in whatever order the two list them; a state whose names differ from this sampler's, forest by forest, is refused under either forceUpdate, and one that records none, as from a sampler without such a forest, goes in by position. The state goes in as stored: its trees, leaf values, linear slopes, k, sigma and kept draws are numbers on the sampler's internal scale, read against this sampler's response mapping and its leaf covariates' centers and scales, which stay its own, and never converted; kept linear draws then predict under the sampler's standardization until replaced, and saved Gaussian-process draws made under another standardization or lengthscale are refused. A state stored under the mapping in force puts the chain back where it was stored, value for value; one stored before a re-anchor, or by another sampler, is legal and has no special meaning. Of the state's kept draws the newest that fit this sampler's store are taken. A state one of whose trees would have to be changed to fit this sampler - a split whose threshold is not one of this sampler's cut points, a leaf no row reaches, a split outside its interval or on a column, or in an interaction, its forest does not allow, a missing-value side on a column that has never held a missing value, monotone leaf values out of order, DART split weights of 0 on a column the forest may split on - is not installed: the sampler is left as it was and FALSE is returned, where TRUE says the state went in. With forceUpdate = TRUE such a state is installed, the trees repaired - a split off the cut points moved to the nearest point its tree still allows, or merged with what is beneath it into one leaf - and NULL is returned invisibly (DART weights of 0 on a splittable column start fresh, equal over the columns the forest may split on, with a drawn concentration and the delay counter restarted). The sampler's stored state is then the chain as installed, repairs included, not the state as handed in, so it installs again without forceUpdate, returning TRUE, and a copy or a reload starts from it. A state does not bring its cut points: to restore them too, hold the state in a variable st and call setCutPoints(attr(st, \"cutPoints\")), the whole set, before setState(st). A column whose cut points were changed and that attr(st, \"cutMass\") lists then has its points and not its weights: the restore still returns TRUE and the chain runs with that column's points equally likely; setPredictor(x[, j], j, updateCutPoints = \"value\") restores both while the column holds the values its cut points were derived from. A sampler whose engine is gone, as after a save and load, is first rebuilt on the cut points, weights, response range and leaf standardization its data object records, and from its own stored state where it has one. See Saving and Value in ?dbartsSampler. The latents it carries are the values its trees were drawn from. This draws nothing; the next sweep draws the augmentation variables against the change."
      refuseLegacyState(newState)
      if (!inherits(newState, "bartcoreState")) {
        stop("'state' must inherit from bartcoreState")
      }
      forceUpdate <- checkFlag(forceUpdate, "forceUpdate")
      # a state whose levels this sampler does not have is refused before a
      # sampler whose engine is gone is rebuilt, so it is left as it was
      stateLevelOrder(newState, data)
      selfEnv <- parent.env(environment())
      ptr <- pointer
      if (.Call(C_dbarts_bartcore_isValidPointer, pointer) == FALSE) {
        # a sampler whose engine is gone is first rebuilt on the cut points
        # its data object records, from its own saved state where it has one
        # that installs, as getPointer rebuilds it, so that the call below
        # gives what it gives a live sampler
        ptr <- recreatePointer(control, model, data)
        if (!is.null(state)) {
          rebuilt <- tryCatch(
            {
              refuseLegacyState(state)
              installStateOnto(
                .self,
                ptr,
                state,
                rawPredictorMatrix(data@x),
                TRUE,
                TRUE
              )
            },
            error = function(e) 0L
          )
          if (rebuilt == 0L) {
            ptr <- recreatePointer(control, model, data)
          } else {
            # this engine carries the sampler's own state: bound now, a state
            # then refused or declined leaves the sampler it would leave
            # live. Where that state went onto other cut points than it
            # records the field is the engine's own, and stands for the field
            # where the field is what was handed in
            selfEnv$pointer <- ptr
            if (rebuilt == 2L) {
              held <- storedState(ptr, data)
              if (identical(newState, state)) {
                newState <- held
              }
              selfEnv$state <- held
            }
          }
        }
      }
      status <- installStateOnto(
        .self,
        ptr,
        newState,
        rawPredictorMatrix(data@x),
        FALSE,
        forceUpdate
      )
      # as in getPointer: an engine created bare is bound only after the
      # install succeeds, so a refusal, or a state declined, leaves a dead
      # pointer dead instead of live and unfitted, and leaves 'state' the one
      # that is still installed
      if (status == 0L) {
        return(FALSE)
      }
      selfEnv$pointer <- ptr
      # a state that went onto other cut points or weights than it records,
      # or that was forced in, is no longer what the engine holds where a
      # tree was repaired: a copy or a reload made from it would take its
      # trees unmoved, so the engine's own chain is the one kept
      if (status == 2L || forceUpdate) {
        storeState(ptr)
      } else {
        selfEnv$state <- newState
      }
      if (!forceUpdate) TRUE else invisible(NULL)
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
      "Updates the cached internal state used for saving/loading: the chains and the units they are stored in, no prior or fixed value, with the cut points and weights in force, against which a reload moves the stored splits onto the ones its data object records. Call it before saving after the cut points change, so that a reload has the chain as it stands."
      selfEnv <- parent.env(environment())
      selfEnv$state <- storedState(ptr, data)
      invisible(NULL)
    },
    installTrees = function(donor, samples = NULL, forceUpdate = FALSE) {
      "Warm-starts the forests from a donor sampler or bart fit over the same
       predictors, keeping this sampler's model and cut points. The donor's
       trees and leaf values, and its k and sigma where this sampler draws
       them, go in as stored on the internal scale, not converted to this
       sampler's response; so do linear-leaf coefficients, numbers per unit
       of this sampler's own centered and scaled covariates, which stay this
       sampler's, so a seeded fit starts away from the donor's where the
       two standardize a covariate differently and the draws then move it. Where
       the donor's data could hold no missing value in a column this
       sampler's can, each chain draws the side of every rule on that column
       in the trees it takes. 'samples' maps each chain to a 1-based
       donor-sample index; NULL spreads the chains across the donor's kept
       samples. A donor every tree of which fits this sampler as stored is
       installed and TRUE is returned. One that would have to be changed - a
       split off this sampler's cut points, a leaf no row reaches, a split
       its forest does not allow, a missing-value side on a column that has
       never held a missing value, monotone leaf values out of order, linear
       slopes this sampler does not take, more trees than this sampler holds
       or, under blocks, trees the blocks cannot take as they stand - is not
       installed: the sampler is left as it was and FALSE is returned. With
       forceUpdate = TRUE such a donor is installed with the changes, the
       trees whose fits vary most going in where not all can, and NULL is
       returned invisibly. A donor of fewer trees gives all of them, the rest
       starting as single leaves. DART weights and a variance forest one side
       alone has are not judged: the donor's is left behind and this
       sampler's starts as a new fit's. Gaussian-process leaves go only to
       Gaussian-process leaves on as many rows. Single-forest samplers only.
       See Warm starts in ?dbartsSampler."
      forceUpdate <- checkFlag(forceUpdate, "forceUpdate")
      report <- installDonorTrees(.self, donor, samples, forceUpdate)
      if (!forceUpdate) report$installed else invisible(NULL)
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
        treeNums <- seq_len(model@forests[[1L]]@n.trees)
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
      "Returns a data.frame containing the internal state of the trees, one row per node. A sampler with more than one forest (a multinomial one, or one declared with forests = naming several) puts a leading 'forest' column (indexed from 1) on it, and stacks every forest forest-major at the default forest = NULL, as the sampler's other per-forest readers stack at their own default; any other sampler returns no forest column, and accepts forest = 1. forest also takes a single position or label, or a vector of them (positions, or labels as a string is read elsewhere; forest<i> names position i), each validated as getLeafPrior/getForestFits/getForestAmplitudes/getForestVariableCounts validate one. treeNums defaults to, and is validated against, EACH selected forest's own tree count, which need not match forest 1's."
      matchedCall <- match.call()
      current <- checkFlag(current, "current")
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
        newdata <- validateXTest(
          newdata,
          data@x,
          missing.seen = dataMissingSeen(data)
        )
      }

      ptr <- getPointer()
      # NULL is every forest, forest-major, as getForestAmplitudes stacks.
      # The forest column appears only where the sampler has more than one
      # forest, so a single-forest sampler's table is the one 0.9-x returned.
      # The bound is checked here, ahead of the .Call, so an out-of-range
      # forest reads the same "forest index out of range" the sibling readers
      # raise rather than getTrees' own bridge-side wording; an empty forest
      # vector is refused the same way a scalar one's length check would
      # refuse it.
      forestIndices <- if (is.null(forest)) {
        seq_len(bartcoreNumForests(ptr)) - 1L
      } else if (length(forest) == 0L) {
        samplerForestIndex(forest, model, data, getPointer())
      } else {
        indices <- if (is.character(forest)) {
          samplerForestIndex(forest, model, data, ptr, several = TRUE)
        } else {
          selectForest(forest, NULL, NA_integer_)
          vapply(forest, resolveForestIndex, 0L)
        }
        if (any(indices >= bartcoreNumForests(ptr))) {
          stop("forest index out of range")
        }
        indices
      }
      # saved-tree replay reads the current training predictors (the engine
      # keeps no matrix); a sparse data@x is skipped for a NULL replay source
      trainingMatrix <- rawPredictorMatrix(data@x)
      hasForestColumn <- bartcoreNumForests(ptr) > 1L
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
            "] for ",
            doorForestLabel(model, forestIndex + 1L)
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
      "Minimialist visualization of tree branching and contents. forest, as with getTrees, defaults to the sampler's only forest and is required on a sampler with more than one; it is a number, the forest's position from 1, or a string, its label (forest<i> names position i), and one only."

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
        samplerForestIndex(forest, model, data, getPointer())
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
