pdbart.getAndInitializeSampler <- function(bartCall, evalEnv) {
  # the two doors spell it differently; a stored call names the door it was
  # made through, so the spelling follows the call rather than the caller
  isLegacyDoor <- bartCall[[1L]] == quote(bartBT) ||
    bartCall[[1L]] == quote(dbarts::bartBT)
  samplerOnlyName <- if (isLegacyDoor) "sampleronly" else "samplerOnly"
  if (!is.null(bartCall[[samplerOnlyName]])) {
    stop(
      "'",
      samplerOnlyName,
      "' is set internally by pdbart/pd2bart and cannot be overridden"
    )
  }
  bartCall[[samplerOnlyName]] <- TRUE

  sampler <- eval(bartCall, evalEnv)

  control <- sampler$control
  verbose <- control@verbose
  keepTrainingFits <- control@keepTrainingFits
  control@verbose <- control@keepTrainingFits <- FALSE
  sampler$setControl(control)

  # a run of no sweeps at all is refused, as in 0.9-x, so nskip = 0 skips the
  # burn-in phase rather than asking for one
  samples <- if (sampler$control@n.burn > 0L) {
    sampler$run(0L, sampler$control@n.burn, updateState = FALSE)
  }
  fit <- list(first.sigma = samples[["sigma"]])
  control@verbose <- verbose
  control@keepTrainingFits <- keepTrainingFits
  sampler$setControl(control)
  # under keepTrees the callers predict from the saved trees instead of running,
  # so the sampling phase happens here rather than in their own branch. Burn-in
  # records nothing, so without this the store they read holds no draws at all.
  if (sampler$control@keepTrees) {
    invisible(sampler$run(0L, sampler$control@n.samples))
  }
  namedList(sampler, fit)
}

# Shared preamble for pdbart/pd2bart: obtain a sampler (and any accompanying
# bart fit) from whatever form of 'x.train' was supplied, for use predicting
# from or running with a total prediction matrix. 'name' is the caller name
# ("pdbart"/"pd2bart") used only in the diagnostic messages.
pdbart.prologue <- function(x.train, matchedCall, callingEnv, name) {
  sampler <- fit <- NULL
  # the formals of pdbart or pd2bart, whose call this is
  callerFormals <- formals(sys.function(sys.parent()))
  if (is.matrix(x.train) || is.data.frame(x.train) || is.formula(x.train)) {
    # pdbart/pd2bart carry the BayesTree spelling themselves (x.train,
    # y.train, and BayesTree names through '...'), so the fit they build is
    # the legacy door's, at 0.9-34's defaults
    bartCall <- redirectCall(
      matchedCall,
      dbarts::bartBT,
      callFormals = callerFormals
    )
    massign[sampler, fit] <- pdbart.getAndInitializeSampler(
      bartCall,
      callingEnv
    )
  } else if (inherits(x.train, "dbartsSampler")) {
    sampler <- x.train
    fit <- list()
    if (!sampler$control@keepTrees) {
      # shared with the bart-fit-object branch below: the input the caller
      # supplied cannot serve
      # the call as given, so one is substituted or regenerated. The thread
      # and draw fallbacks narrow this class instead of sharing it
      warning(warningCondition(
        paste0(
          "calling ",
          name,
          " with a sampler that does not have keepTrees set to TRUE will cause new samples to be generated and the state to be changed"
        ),
        class = c("dbartsFallbackWarning", "dbartsWarning")
      ))
    }
  } else if (inherits(x.train, "bart")) {
    fit <- x.train
    sampler <- fit$fit
    if (is.null(sampler)) {
      bartCall <- fit$call
      if (
        !is.call(bartCall) ||
          identical(bartCall, call("NA")) ||
          identical(bartCall, call("NULL"))
      ) {
        stop(
          "calling ",
          name,
          " with a bart fit object requires model to be fit with keepSampler == TRUE"
        )
      }
      warning(warningCondition(
        paste0(
          "calling ",
          name,
          " with a bart fit object requires model to be fit with keepSampler == TRUE; refitting using saved call"
        ),
        class = c("dbartsFallbackWarning", "dbartsWarning")
      ))
      massign[sampler, fit] <- pdbart.getAndInitializeSampler(
        bartCall,
        callingEnv
      )
    }
  } else if (
    inherits(
      x.train,
      c("bartMultinomial", "bartOrdinal", "bartNegbin", "bartHurdle")
    )
  ) {
    stop(name, " does not support a ", class(x.train)[1L], " fit")
  } else {
    stop(
      "'x.train' must be a matrix, data.frame, formula, fitted bart model, ",
      "or dbartsSampler"
    )
  }
  namedList(sampler, fit)
}

# Resolve the 'xind' argument (a formula-style expression, character column
# names, numeric indices, or NULL) down to integer column indices into the
# sampler's predictor matrix. 'xind' is forwarded as a promise so that the
# non-standard formula evaluation in the error branch behaves as if inlined.
pdbart.resolveXind <- function(xind, matchedCall, sampler) {
  tryResult <- tryCatch(xind, error = I)

  if (inherits(tryResult, "error")) {
    formula <- ~a
    formula[[2L]] <- matchedCall[["xind"]]
    terms <- terms(formula)

    xind <- attr(terms, "term.labels")
  } else if (
    !inherits(tryResult, "error") &&
      is.character(xind) &&
      length(xind) == 1L &&
      xind %not_in% colnames(sampler$data@x)
  ) {
    formula <- ~a
    formula[[2L]] <- parse(text = xind)[[1L]]
    terms <- terms(formula)

    xind <- attr(terms, "term.labels")
  } else if (is.null(xind)) {
    xind <- seq_len(ncol(sampler$data@x))
  }

  if (is.character(xind)) {
    if (is.null(colnames(sampler$data@x))) {
      stop("passing 'xind' by name requires 'x.train' to have column names")
    }
    unknownColumns <- xind %not_in% colnames(sampler$data@x)
    if (any(unknownColumns)) {
      stop(
        "unrecognized columns '",
        paste0(xind[unknownColumns], collapse = "', '"),
        "'"
      )
    }
    xind <- match(xind, colnames(sampler$data@x))
  }

  xind
}

# The level table of each selected predictor that is a factor (categorical or
# ordered), NULL for any other column.
pdbart.factorLevels <- function(sampler, xind) {
  factorLevels <- attr(sampler$data@x, "factor.levels")
  lapply(xind, function(j) {
    if (is.null(factorLevels)) NULL else factorLevels[[j]]
  })
}

# Default the 'levs' list: for each of the first 'numVariables' selected
# predictors, every level of a factor, by name; otherwise either the sorted
# unique values (when there are too few to bin) or the unique quantiles at
# 'levquants'. 'cmp' is the comparison deciding "too few": pdbart uses `<`,
# pd2bart uses `<=` (a long-standing difference in the two entry points,
# preserved here rather than reconciled).
pdbart.defaultLevs <- function(x, xind, levquants, numVariables, cmp, levels) {
  levs <- vector("list", numVariables)
  for (j in seq_len(numVariables)) {
    if (!is.null(levels[[j]])) {
      levs[[j]] <- levels[[j]]
      next
    }
    uniqueValues <- unique(x[, xind[j]])
    levs[[j]] <-
      if (cmp(length(uniqueValues), length(levquants))) {
        sort(uniqueValues)
      } else {
        unique(quantile(x[, xind[j]], probs = levquants))
      }
  }
  levs
}

# Validates user 'levs' against the factor columns: a factor's values are given
# by level name (a character vector or a factor), as the results report them.
pdbart.checkLevs <- function(levs, levels, xLabels) {
  for (j in seq_along(levs)) {
    if (is.null(levels[[j]])) {
      next
    }
    values <- levs[[j]]
    if (is.factor(values)) {
      values <- as.character(values)
    }
    if (!is.character(values)) {
      stop(
        "'levs' for factor predictor '",
        xLabels[j],
        "' must name its levels"
      )
    }
    unknown <- values[values %not_in% levels[[j]]]
    if (length(unknown) > 0L) {
      stop(
        "'levs' for factor predictor '",
        xLabels[j],
        "' names levels not present in training: ",
        paste0("'", unique(unknown), "'", collapse = ", ")
      )
    }
    levs[[j]] <- values
  }
  levs
}

# The value a predictor column takes at one 'levs' entry: a factor level's
# 0-based code, any other column's value itself.
pdbart.levelValues <- function(levs, levels) {
  lapply(seq_along(levs), function(j) {
    if (is.null(levels[[j]])) {
      levs[[j]]
    } else {
      match(levs[[j]], levels[[j]]) - 1
    }
  })
}

pdbart.xLabels <- function(sampler, xind) {
  if (is.null(colnames(sampler$data@x))) {
    paste0("x", xind)
  } else {
    colnames(sampler$data@x)[xind]
  }
}

# Per-draw predictions at each row of a prediction channel, as draws x rows
# with the chains in turn, the layout pdbart.drawMeans gives its means.
pdbart.drawsByRow <- function(pred) {
  t(matrix(pred, nrow = dim(pred)[1L]))
}

# The per-draw mean over the observation margin of a prediction channel: the
# observations are the channel's FIRST margin and the draw (and chain) margins
# the trailing ones, so each draw's observations are one contiguous slab and
# the reduction takes them without permuting the channel, which is n times the
# result. Several chains flatten in the array's own order - each chain's whole
# run in turn - which is the layout fdr's columns already held.
pdbart.drawMeans <- function(pred, n.chains) {
  if (n.chains > 1L) {
    as.vector(channelMeans(pred, 2L))
  } else {
    channelMeans(pred)
  }
}

# Assemble the returned pdbart/pd2bart result list. Identical between the two
# entry points except for the S3 class stamped on it ('className').
pdbart.buildResult <- function(sampler, fit, fdr, levs, xind, className) {
  xLabels <- pdbart.xLabels(sampler, xind)

  if (sampler$control@binary == FALSE) {
    result <- list(
      fd = fdr,
      levs = levs,
      xlbs = xLabels,
      bartcall = sampler$control@call,
      yhat.train = fit$yhat.train,
      first.sigma = fit$first.sigma,
      sigma = fit$sigma,
      yhat.train.mean = fit$yhat.train.mean,
      sigest = sampler$data@sigma,
      y = sampler$data@y,
      fit = sampler
    )
  } else {
    result <- list(
      fd = fdr,
      levs = levs,
      xlbs = xLabels,
      bartcall = fit$call,
      yhat.train = fit$yhat.train,
      y = sampler$data@y,
      fit = sampler
    )
  }
  class(result) <- className
  result
}

## create the contents to be used in partial dependence plots
pdbart <- function(
  x.train,
  y.train,
  xind = NULL,
  levs = NULL,
  levquants = c(0.05, seq(0.1, 0.9, 0.1), 0.95),
  pl = TRUE,
  plquants = c(0.05, 0.95),
  ...
) {
  matchedCall <- match.call()

  callingEnv <- parent.frame()

  sampler <- fit <- NULL ## for R CMD check (massign assigns these below)
  massign[sampler, fit] <- pdbart.prologue(
    x.train,
    matchedCall,
    callingEnv,
    "pdbart"
  )

  xind <- pdbart.resolveXind(xind, matchedCall, sampler)

  numVariables <- length(xind)

  # materialize the predictor codes once: a dense-frame/mixed container serves
  # them through as.matrix, a plain matrix (or dgCMatrix) is itself
  x <- extract(sampler, "predictors")
  levels <- pdbart.factorLevels(sampler, xind)

  if (is.null(levs)) {
    levs <- pdbart.defaultLevs(x, xind, levquants, numVariables, `<`, levels)
  } else if (length(levs) != numVariables) {
    stop("'levs' must have the same length as 'xind'")
  } else {
    levs <- pdbart.checkLevs(levs, levels, pdbart.xLabels(sampler, xind))
  }
  values <- pdbart.levelValues(levs, levels)

  numLevels <- sapply(levs, length)
  numSamples <- sampler$control@n.samples * sampler$control@n.chains

  if (sampler$control@keepTrees == TRUE) {
    fdr <- vector("list", numVariables)
    for (j in seq_len(numVariables)) {
      fdr[[j]] <- matrix(NA_real_, numSamples, numLevels[j])
      for (i in seq_len(numLevels[j])) {
        x.test <- x
        x.test[, xind[j]] <- values[[j]][i]

        pred <- pdbart.drawMeans(
          sampler$predict(x.test),
          sampler$control@n.chains
        )

        .Call(C_dbarts_assignInPlace, fdr[[j]], i, pred)
      }
    }
  } else {
    x.test <- NULL
    for (j in seq_len(numVariables)) {
      for (i in seq_len(numLevels[j])) {
        temp <- x
        temp[, xind[j]] <- values[[j]][i]
        x.test <- rbind(x.test, temp)
      }
    }
    sampler$setTestPredictor(x.test)

    samples <- sampler$run(0L, sampler$control@n.samples)
    if (is.null(fit[["call"]])) {
      fit <- packageBartResults(
        sampler,
        samples,
        fit$sigma,
        fit[["k"]],
        TRUE,
        TRUE
      )
      fit[["yhat.test"]] <- NULL
    }

    numObservations <- length(sampler$data@y)
    fdr <- vector("list", numVariables)
    offset <- 0
    for (j in seq_len(numVariables)) {
      fdr[[j]] <- matrix(NA_real_, numSamples, numLevels[j])
      for (i in seq_len(numLevels[j])) {
        indices <- seq.int(
          offset + (i - 1) * numObservations + 1,
          offset + i * numObservations
        )

        pred <- pdbart.drawMeans(
          if (sampler$control@n.chains > 1L) {
            samples$test[indices, , ]
          } else {
            samples$test[indices, ]
          },
          sampler$control@n.chains
        )

        .Call(C_dbarts_assignInPlace, fdr[[j]], i, pred)
      }
      offset <- offset + numObservations * numLevels[j]
    }
  }

  result <- pdbart.buildResult(sampler, fit, fdr, levs, xind, "pdbart")

  if (pl) {
    plot(result, plquants = plquants)
  }

  result
}

pd2bart <- function(
  x.train,
  y.train,
  xind = NULL,
  levs = NULL,
  levquants = c(0.05, seq(0.1, 0.9, 0.1), 0.95),
  pl = TRUE,
  plquants = c(0.05, 0.95),
  ...
) {
  matchedCall <- match.call()

  callingEnv <- parent.frame()

  sampler <- fit <- NULL ## for R CMD check (massign assigns these below)
  massign[sampler, fit] <- pdbart.prologue(
    x.train,
    matchedCall,
    callingEnv,
    "pd2bart"
  )

  xind <- pdbart.resolveXind(xind, matchedCall, sampler)

  # materialize the predictor codes once: a dense-frame/mixed container serves
  # them through as.matrix, a plain matrix (or dgCMatrix) is itself
  x <- extract(sampler, "predictors")
  levels <- pdbart.factorLevels(sampler, xind)

  if (is.null(levs)) {
    levs <- pdbart.defaultLevs(x, xind, levquants, 2L, `<=`, levels)
  } else {
    levs <- pdbart.checkLevs(levs, levels, pdbart.xLabels(sampler, xind))
  }
  values <- pdbart.levelValues(levs, levels)
  numSamples <- sampler$control@n.samples * sampler$control@n.chains

  xValues <- as.matrix(expand.grid(values[[1L]], values[[2L]]))
  numXValues <- nrow(xValues)

  # with two predictors each grid point is itself a whole row, so its
  # prediction needs no average over the training rows
  gridAsRows <- function() {
    x.test <- if (xind[1L] < xind[2L]) xValues else xValues[, c(2L, 1L)]
    colnames(x.test) <- colnames(x)
    x.test
  }

  if (sampler$control@keepTrees == TRUE) {
    if (ncol(sampler$data@x) == 2L) {
      fdr <- pdbart.drawsByRow(sampler$predict(gridAsRows()))
    } else {
      fdr <- matrix(NA_real_, numSamples, numXValues)
      for (i in seq_len(numXValues)) {
        x.test <- x
        x.test[, xind[1L]] <- xValues[i, 1L]
        x.test[, xind[2L]] <- xValues[i, 2L]

        pred <- pdbart.drawMeans(
          sampler$predict(x.test),
          sampler$control@n.chains
        )

        .Call(C_dbarts_assignInPlace, fdr, i, pred)
      }
    }
  } else {
    if (ncol(sampler$data@x) == 2L) {
      sampler$setTestPredictor(gridAsRows())
      samples <- sampler$run(0L, sampler$control@n.samples)
      fdr <- pdbart.drawsByRow(samples$test)
    } else {
      x.test <- NULL
      for (i in seq_len(numXValues)) {
        temp <- x
        temp[, xind[1L]] <- xValues[i, 1L]
        temp[, xind[2L]] <- xValues[i, 2L]
        x.test <- rbind(x.test, temp)
      }
      sampler$setTestPredictor(x.test)
      samples <- sampler$run(0L, sampler$control@n.samples)

      numObservations <- length(sampler$data@y)

      fdr <- matrix(NA_real_, numSamples, numXValues)
      for (i in seq_len(numXValues)) {
        indices <- seq.int((i - 1) * numObservations + 1, i * numObservations)
        pred <- pdbart.drawMeans(
          if (sampler$control@n.chains > 1L) {
            samples$test[indices, , ]
          } else {
            samples$test[indices, ]
          },
          sampler$control@n.chains
        )
        .Call(C_dbarts_assignInPlace, fdr, i, pred)
      }
    }
    if (is.null(fit[["call"]])) {
      fit <- packageBartResults(
        sampler,
        samples,
        fit$sigma,
        fit[["k"]],
        TRUE,
        TRUE
      )
      fit[["yhat.test"]] <- NULL
    }
  }

  result <- pdbart.buildResult(sampler, fit, fdr, levs, xind, "pd2bart")

  if (pl) {
    plot(result, plquants = plquants)
  }

  result
}
