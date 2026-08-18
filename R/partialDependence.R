pdbart.getAndInitializeSampler <- function(bartCall, evalEnv) {
  isBart2 <- bartCall[[1L]] == quote(bart2) ||
    bartCall[[1L]] == quote(dbarts::bart2)
  samplerOnlyName <- if (isBart2) "samplerOnly" else "sampleronly"
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

  samples <- sampler$run(0L, sampler$control@n.burn, updateState = FALSE)
  fit <- list(first.sigma = samples$sigma)
  control@verbose <- verbose
  control@keepTrainingFits <- keepTrainingFits
  sampler$setControl(control)
  namedList(sampler, fit)
}

# Shared preamble for pdbart/pd2bart: obtain a sampler (and any accompanying
# bart fit) from whatever form of 'x.train' was supplied, for use predicting
# from or running with a total prediction matrix. 'name' is the caller name
# ("pdbart"/"pd2bart") used only in the diagnostic messages.
pdbart.prologue <- function(x.train, matchedCall, callingEnv, name) {
  sampler <- fit <- NULL
  if (is.matrix(x.train) || is.data.frame(x.train) || is.formula(x.train)) {
    bartCall <- redirectCall(matchedCall, dbarts::bart)
    massign[sampler, fit] <- pdbart.getAndInitializeSampler(
      bartCall,
      callingEnv
    )
  } else if (inherits(x.train, "dbartsSampler")) {
    sampler <- x.train
    fit <- list()
    if (!sampler$control@keepTrees) {
      warning(
        "calling ",
        name,
        " with a sampler that does not have keepTrees set to TRUE will cause new samples to be generated and the state to be changed"
      )
    }
  } else if (inherits(x.train, "bart")) {
    fit <- x.train
    sampler <- fit$fit
    if (is.null(sampler)) {
      bartCall <- fit$call
      if (bartCall == call("NA") || bartCall == call("NULL")) {
        stop(
          "calling ",
          name,
          " with a bart fit object requires model to be fit with keepSampler == TRUE"
        )
      }
      warning(
        "calling ",
        name,
        " with a bart fit object requires model to be fit with keepSampler == TRUE; refitting using saved call"
      )
      massign[sampler, fit] <- pdbart.getAndInitializeSampler(
        bartCall,
        callingEnv
      )
    }
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

# Default the 'levs' list: for each of the first 'numVariables' selected
# predictors, either the sorted unique values (when there are too few to bin)
# or the unique quantiles at 'levquants'. 'cmp' is the comparison deciding
# "too few": pdbart uses `<`, pd2bart uses `<=` (a long-standing difference in
# the two entry points, preserved here rather than reconciled).
pdbart.defaultLevs <- function(x, xind, levquants, numVariables, cmp) {
  levs <- vector("list", numVariables)
  for (j in seq_len(numVariables)) {
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

# Assemble the returned pdbart/pd2bart result list. Identical between the two
# entry points except for the S3 class stamped on it ('className').
pdbart.buildResult <- function(sampler, fit, fdr, levs, xind, className) {
  if (is.null(colnames(sampler$data@x))) {
    xLabels <- paste0("x", xind)
  } else {
    xLabels <- colnames(sampler$data@x)[xind]
  }

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

  if (is.null(levs)) {
    levs <- pdbart.defaultLevs(x, xind, levquants, numVariables, `<`)
  } else if (length(levs) != numVariables) {
    stop("length of 'levs' must equal that of 'xind'")
  }

  numLevels <- sapply(levs, length)
  numSamples <- sampler$control@n.samples * sampler$control@n.chains

  if (sampler$control@keepTrees == TRUE) {
    fdr <- vector("list", numVariables)
    for (j in seq_len(numVariables)) {
      fdr[[j]] <- matrix(NA_real_, numSamples, numLevels[j])
      for (i in seq_len(numLevels[j])) {
        x.test <- x
        x.test[, xind[j]] <- levs[[j]][i]

        pred <-
          if (sampler$control@n.chains > 1L) {
            as.vector(apply(sampler$predict(x.test), c(2L, 3L), mean))
          } else {
            apply(sampler$predict(x.test), 2L, mean)
          }

        .Call(C_dbarts_assignInPlace, fdr[[j]], i, pred)
      }
    }
  } else {
    x.test <- NULL
    for (j in seq_len(numVariables)) {
      for (i in seq_len(numLevels[j])) {
        temp <- x
        temp[, xind[j]] <- levs[[j]][i]
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

        pred <-
          if (sampler$control@n.chains > 1L) {
            as.vector(apply(samples$test[indices, , ], c(2L, 3L), mean))
          } else {
            apply(samples$test[indices, ], 2L, mean)
          }

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

  if (is.null(levs)) {
    levs <- pdbart.defaultLevs(x, xind, levquants, 2L, `<=`)
  }
  numSamples <- sampler$control@n.samples * sampler$control@n.chains

  xValues <- as.matrix(expand.grid(levs[[1L]], levs[[2L]]))
  numXValues <- nrow(xValues)

  if (sampler$control@keepTrees == TRUE) {
    if (ncol(sampler$data@x) == 2L) {
      x.test <- if (xind[1L] < xind[2L]) xValues else xValues[, c(2L, 1L)]
      pred <- suppressWarnings(sampler$predict(x.test))
      fdr <- as.matrix(
        if (sampler$control@n.chains > 1L) {
          as.vector(apply(pred, c(2L, 3L), mean))
        } else {
          apply(pred, 2L, mean)
        }
      )
    } else {
      fdr <- matrix(NA_real_, numSamples, numXValues)
      for (i in seq_len(numXValues)) {
        x.test <- x
        x.test[, xind[1L]] <- xValues[i, 1L]
        x.test[, xind[2L]] <- xValues[i, 2L]

        pred <-
          if (sampler$control@n.chains > 1L) {
            as.vector(apply(sampler$predict(x.test), c(2L, 3L), mean))
          } else {
            apply(sampler$predict(x.test), 2L, mean)
          }

        .Call(C_dbarts_assignInPlace, fdr, i, pred)
      }
    }
  } else {
    if (ncol(sampler$data@x) == 2L) {
      x.test <- if (xind[1L] < xind[2L]) xValues else xValues[, c(2L, 1L)]
      sampler$setTestPredictor(x.test)
      samples <- sampler$run(0L, sampler$control@n.samples)
      fdr <- as.matrix(
        if (sampler$control@n.chains > 1L) {
          as.vector(apply(samples$test, c(2L, 3L), mean))
        } else {
          apply(samples$test, 2L, mean)
        }
      )
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
        pred <-
          if (sampler$control@n.chains > 1L) {
            as.vector(apply(samples$test[indices, , ], c(2L, 3L), mean))
          } else {
            apply(samples$test[indices, ], 2L, mean)
          }
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
