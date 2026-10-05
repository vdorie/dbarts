# The arguments pdbart and pd2bart keep for themselves; everything else in the
# call is bart's.
pdbart.ownArgs <- c("xind", "levs", "levquants", "pl", "plquants")

# bart's arguments pdbart sets itself, under either spelling, each mapped to
# the name the refusal gives.
pdbart.setInternally <- c(
  samplerOnly = "samplerOnly",
  sampleronly = "samplerOnly",
  test = "test",
  x.test = "test",
  offset.test = "offset.test"
)

# Families whose prediction is not one value per row per draw, and families
# this version does not yet serve. Refused by name before anything is fit.
pdbart.refuseFamily <- function(token, caller) {
  if (is.null(token) || length(token) != 1L || is.na(token)) {
    return(invisible(NULL))
  }
  if (token %in% c("multinomial", "ordinal")) {
    stop(
      "'",
      caller,
      "' does not serve ",
      if (token == "ordinal") "an " else "a ",
      token,
      " fit, whose prediction is a probability per category; predict on ",
      "new rows with the variable set gives each category's",
      call. = FALSE
    )
  }
  if (
    token %in%
      c(
        "nbinom",
        "hurdle.lognormal",
        "aft",
        "hazard",
        "hazard.probit",
        "hazard.logistic"
      )
  ) {
    stop(
      "'",
      caller,
      "' does not yet serve family = \"",
      token,
      "\"",
      call. = FALSE
    )
  }
  invisible(NULL)
}

pdbart.fitFamily <- function(fit) {
  if (inherits(fit, "bartMultinomial")) {
    "multinomial"
  } else if (inherits(fit, "bartOrdinal")) {
    "ordinal"
  } else if (inherits(fit, "bartNegbin")) {
    "nbinom"
  } else if (inherits(fit, "bartHurdle")) {
    "hurdle.lognormal"
  } else if (is.character(fit$family)) {
    fit$family[1L]
  }
}

# A hazard sampler's model reads as its binary link; the period grid on its
# control is what marks it.
pdbart.samplerFamily <- function(sampler) {
  if (!is.null(attr(sampler$control, "bartcore.hazard.periods"))) {
    "hazard"
  } else {
    sampler$model@family
  }
}

# The family a data call would fit, read before fitting: the caller's own,
# or under "auto" the one bart resolves from the response.
pdbart.dataFamily <- function(call, object, getData, callingEnv, caller) {
  token <- resolveFamily(
    call$family,
    eval(formals(dbarts::bart)$family),
    caller,
    callingEnv
  )@token
  if (token != "auto") {
    return(token)
  }
  data <- getData()
  dataMissing <- is.null(data)
  data <- if (dataMissing) NULL else data[[1L]]
  response <- autoRawResponse(object, data, dataMissing, callingEnv)
  if (inherits(response, "Surv")) {
    "aft"
  } else if (
    !is.null(detectAutoCounts(object, data, dataMissing, callingEnv)) ||
      !is.null(detectAutoMultinomial(object, data, dataMissing, callingEnv))
  ) {
    "multinomial"
  } else if (
    !is.null(detectAutoOrdinal(object, data, dataMissing, callingEnv))
  ) {
    "ordinal"
  } else {
    token
  }
}

pdbart.flag <- function(value, name) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop("'", name, "' must be TRUE or FALSE", call. = FALSE)
  }
  value
}

# A data call: the caller's call rewritten into a bart call, with trees and
# sampler kept, and evaluated where the caller wrote it, so an unevaluated
# argument resolves as bart would resolve it there.
pdbart.fitData <- function(object, getData, matchedCall, callingEnv, caller) {
  call <- matchedCall[names(matchedCall) %not_in% pdbart.ownArgs]
  argNames <- names(call)[-1L]
  refused <- intersect(argNames, names(pdbart.setInternally))
  if (length(refused) > 0L) {
    stop(
      "'",
      refused[1L],
      "' (bart's '",
      pdbart.setInternally[[refused[1L]]],
      "') is set internally by '",
      caller,
      "' and cannot be given",
      call. = FALSE
    )
  }
  for (name in intersect(argNames, c("keepTrees", "keeptrees"))) {
    if (!isTRUE(eval(call[[name]], callingEnv))) {
      stop(
        "'",
        caller,
        "' predicts from the saved trees, so '",
        name,
        "' can only be TRUE",
        call. = FALSE
      )
    }
  }
  legacy <- intersect(argNames, names(pdbartBayesTreeNames))
  call <- translatePdbartCall(call, legacy, callingEnv, caller)
  keepSampler <- if ("keepSampler" %in% names(call)) {
    pdbart.flag(eval(call$keepSampler, callingEnv), "keepSampler")
  } else {
    TRUE
  }
  pdbart.refuseFamily(
    pdbart.dataFamily(call, object, getData, callingEnv, caller),
    caller
  )
  call$keepTrees <- TRUE
  call$keepSampler <- TRUE
  call[[1L]] <- quote(dbarts::bart)
  if (length(legacy) == 0L) {
    notePdbartDefaults(callingEnv)
  }
  fit <- holdingBartNotices(eval(call, callingEnv))
  pdbart.refuseFamily(pdbart.fitFamily(fit), caller)
  namedList(fit, keepSampler)
}

# A fit kept without its trees or its sampler, refit through the function
# that made it, from its stored call, with both kept.
pdbart.refit <- function(fit, callingEnv, caller) {
  refitCall <- fit$call
  if (
    !is.call(refitCall) ||
      identical(refitCall, call("NA")) ||
      identical(refitCall, call("NULL"))
  ) {
    stop(
      "'",
      caller,
      "' needs a fit kept with keepTrees = TRUE, or one whose call was kept ",
      "so that it can be refit",
      call. = FALSE
    )
  }
  warning(warningCondition(
    paste0(
      "calling ",
      caller,
      " with a fit kept without its trees or sampler refits it from its ",
      "stored call; fit with keepTrees = TRUE to avoid this"
    ),
    class = c("dbartsFallbackWarning", "dbartsWarning")
  ))
  # the stored call is matched, so its names say which door made it, whatever
  # name the function was called by
  if ("x.train" %in% names(refitCall)) {
    refitCall[[1L]] <- quote(dbarts::bartBT)
    refitCall$keeptrees <- TRUE
    refitCall$keepsampler <- TRUE
  } else {
    refitCall[[1L]] <- quote(dbarts::bart)
    refitCall$keepTrees <- TRUE
    refitCall$keepSampler <- TRUE
  }
  holdingBartNotices(eval(refitCall, callingEnv))
}

# The sampler pdbart predicts from and the fit it reports, from whatever was
# passed first: data, a fit or a sampler. 'object' is list(value) or NULL
# when nothing was passed; 'getData' returns the data argument the same way.
pdbart.prologue <- function(object, getData, matchedCall, callingEnv, caller) {
  if (is.null(object)) {
    stop(
      "'formula' is required: a matrix, data frame or formula, a fit, or a ",
      "sampler",
      call. = FALSE
    )
  }
  object <- object[[1L]]
  isFit <- inherits(
    object,
    c("bart", "bartMultinomial", "bartOrdinal", "bartNegbin", "bartHurdle")
  )
  if (!isFit && !inherits(object, "dbartsSampler")) {
    massign[fit, keepSampler] <- pdbart.fitData(
      object,
      getData,
      matchedCall,
      callingEnv,
      caller
    )
    return(list(sampler = fit$fit, fit = fit, keepSampler = keepSampler))
  }

  if ("x.train" %in% names(matchedCall)) {
    warnClassed(
      "dbartsDeprecatedWarning",
      "'x.train' holds a ",
      if (isFit) "fit" else "sampler",
      "; pass it to '",
      caller,
      "' first, unnamed"
    )
  }
  matchedCall <- translatePdbartCall(
    matchedCall,
    intersect(names(matchedCall), "keepsampler"),
    callingEnv,
    caller
  )
  extra <- setdiff(
    names(matchedCall)[-1L],
    c("formula", "x.train", "keepSampler", pdbart.ownArgs)
  )
  if (length(extra) > 0L) {
    stop(
      "'",
      extra[1L],
      "' has no effect on a ",
      if (isFit) "fit" else "sampler",
      " passed to '",
      caller,
      "', which is not refit",
      call. = FALSE
    )
  }
  keepSampler <- if ("keepSampler" %in% names(matchedCall)) {
    pdbart.flag(eval(matchedCall$keepSampler, callingEnv), "keepSampler")
  } else {
    TRUE
  }

  if (!isFit) {
    pdbart.refuseFamily(pdbart.samplerFamily(object), caller)
    if (!object$control@keepTrees) {
      warning(warningCondition(
        paste0(
          "calling ",
          caller,
          " with a sampler that does not have keepTrees set to TRUE will ",
          "cause new samples to be generated and the state to be changed"
        ),
        class = c("dbartsFallbackWarning", "dbartsWarning")
      ))
    }
    return(list(sampler = object, fit = list(), keepSampler = keepSampler))
  }

  fit <- object
  pdbart.refuseFamily(pdbart.fitFamily(fit), caller)
  if (is.null(fit$fit) || !fit$fit$control@keepTrees) {
    fit <- pdbart.refit(fit, callingEnv, caller)
  }
  list(sampler = fit$fit, fit = fit, keepSampler = keepSampler)
}

# The value of an argument given under its own name or its BayesTree
# spelling, as list(value), or NULL when given under neither. 'value' is
# forced only when 'given'.
pdbart.argument <- function(given, value, legacyName, matchedCall, callingEnv) {
  if (given) {
    return(list(value))
  }
  if (legacyName %in% names(matchedCall)) {
    return(list(eval(matchedCall[[legacyName]], callingEnv)))
  }
  NULL
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
      stop(
        "passing 'xind' by name requires the predictors to have column names"
      )
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
# 'levquants', missing values left out of both. 'cmp' is the comparison
# deciding "too few": pdbart uses `<`, pd2bart uses `<=` (a long-standing
# difference in the two entry points, preserved here rather than reconciled).
pdbart.defaultLevs <- function(x, xind, levquants, numVariables, cmp, levels) {
  levs <- vector("list", numVariables)
  for (j in seq_len(numVariables)) {
    if (!is.null(levels[[j]])) {
      levs[[j]] <- levels[[j]]
      next
    }
    column <- x[, xind[j]]
    column <- column[!is.na(column)]
    uniqueValues <- unique(column)
    levs[[j]] <-
      if (cmp(length(uniqueValues), length(levquants))) {
        sort(uniqueValues)
      } else {
        unique(quantile(column, probs = levquants))
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

# The rows averaged over: the fit's own, less those it gives a 0 weight or
# masks out, each with its stored offset (NULL when the fit has none).
pdbart.averagedRows <- function(sampler, x) {
  data <- sampler$data
  keep <- rep_len(TRUE, nrow(x))
  if (length(data@weights) > 0L) {
    keep <- keep & data@weights > 0
  }
  if (!is.null(sampler$activeRows)) {
    keep <- keep & sampler$activeRows != 0
  }
  offset <- if (length(data@offset) > 0L) data@offset
  if (!all(keep)) {
    x <- x[keep, , drop = FALSE]
    offset <- offset[keep]
  }
  list(x = x, offset = offset)
}

# 'x' with each of a setting's columns set to its value.
pdbart.setColumns <- function(x, setting) {
  for (k in seq_along(setting$columns)) {
    x[, setting$columns[k]] <- setting$values[k]
  }
  x
}

# Per-draw averages over 'rows' at each of 'settings', draws x settings, a
# setting being list(columns, values). Each row's offset enters its
# prediction. A sampler with saved trees predicts from them; one without runs
# once over every setting's rows stacked, changing its state, and its samples
# come back for the result.
pdbart.drawsAt <- function(sampler, rows, settings) {
  n.chains <- sampler$control@n.chains
  fd <- matrix(
    NA_real_,
    sampler$control@n.samples * n.chains,
    length(settings)
  )
  if (sampler$control@keepTrees) {
    for (i in seq_along(settings)) {
      x.test <- pdbart.setColumns(rows$x, settings[[i]])
      pred <- if (is.null(rows$offset)) {
        sampler$predict(x.test)
      } else {
        sampler$predict(x.test, rows$offset)
      }
      .Call(C_dbarts_assignInPlace, fd, i, pdbart.drawMeans(pred, n.chains))
    }
    return(list(fd = fd, samples = NULL))
  }
  numRows <- nrow(rows$x)
  sampler$setTestPredictor(do.call(
    rbind,
    lapply(settings, pdbart.setColumns, x = rows$x)
  ))
  samples <- sampler$run(0L, sampler$control@n.samples)
  # averaging is linear on this scale, so the rows' offsets enter as their mean
  offsetMean <- if (is.null(rows$offset)) 0 else mean(rows$offset)
  for (i in seq_along(settings)) {
    indices <- seq.int((i - 1L) * numRows + 1L, i * numRows)
    pred <- if (n.chains > 1L) {
      samples$test[indices, , , drop = FALSE]
    } else {
      samples$test[indices, , drop = FALSE]
    }
    .Call(
      C_dbarts_assignInPlace,
      fd,
      i,
      pdbart.drawMeans(pred, n.chains) + offsetMean
    )
  }
  list(fd = fd, samples = samples)
}

# Whether a fit reports its draws split by chain, as its varcount does.
pdbart.chainsSplit <- function(fit, n.chains) {
  n.chains > 1L && length(dim(fit$varcount)) == 3L
}

# Draws x settings, the chains in turn, as chains x draws x settings.
pdbart.splitChains <- function(fd, n.chains) {
  aperm(
    array(fd, c(nrow(fd) %/% n.chains, n.chains, ncol(fd))),
    c(2L, 1L, 3L)
  )
}

# The draws of a sampler passed without saved trees, packaged as a fit is.
pdbart.packageRun <- function(sampler, fit, samples) {
  if (is.null(samples) || !is.null(fit[["call"]])) {
    return(fit)
  }
  fit <- packageBartResults(
    sampler,
    samples,
    fit$first.sigma,
    fit[["k"]],
    TRUE,
    TRUE
  )
  fit[["yhat.test"]] <- NULL
  fit
}

# Assemble the returned pdbart/pd2bart result list. Identical between the two
# entry points except for the S3 class stamped on it ('className').
pdbart.buildResult <- function(
  sampler,
  fit,
  fdr,
  levs,
  xind,
  keepSampler,
  className
) {
  xLabels <- pdbart.xLabels(sampler, xind)
  bartcall <- if (is.null(fit$call)) sampler$control@call else fit$call
  y <- if (is.null(fit$y)) sampler$data@y else fit$y
  n.chains <- sampler$control@n.chains

  result <- if (sampler$control@binary == FALSE) {
    list(
      fd = fdr,
      levs = levs,
      xlbs = xLabels,
      bartcall = bartcall,
      yhat.train = fit$yhat.train,
      first.sigma = fit$first.sigma,
      sigma = fit$sigma,
      yhat.train.mean = fit$yhat.train.mean,
      sigest = if (is.null(fit$sigest)) sampler$data@sigma else fit$sigest,
      y = y,
      n.chains = n.chains,
      fit = sampler
    )
  } else {
    list(
      fd = fdr,
      levs = levs,
      xlbs = xLabels,
      bartcall = bartcall,
      yhat.train = fit$yhat.train,
      y = y,
      n.chains = n.chains,
      fit = sampler
    )
  }
  if (!keepSampler) {
    result$fit <- NULL
  }
  class(result) <- className
  result
}

## create the contents to be used in partial dependence plots
pdbart <- function(
  formula,
  data,
  xind = NULL,
  levs = NULL,
  levquants = c(0.05, seq(0.1, 0.9, 0.1), 0.95),
  pl = TRUE,
  plquants = c(0.05, 0.95),
  ...
) {
  matchedCall <- match.call()
  callingEnv <- parent.frame()

  sampler <- fit <- keepSampler <- NULL ## for R CMD check (massign assigns)
  dataGiven <- !missing(data)
  massign[sampler, fit, keepSampler] <- pdbart.prologue(
    pdbart.argument(
      !missing(formula),
      formula,
      "x.train",
      matchedCall,
      callingEnv
    ),
    function() {
      pdbart.argument(dataGiven, data, "y.train", matchedCall, callingEnv)
    },
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

  rows <- pdbart.averagedRows(sampler, x)
  n.chains <- sampler$control@n.chains
  # every variable's settings in one pass, so that a sampler run without
  # saved trees draws them all from one run
  variable <- rep(seq_len(numVariables), lengths(values))
  settings <- unlist(
    lapply(seq_len(numVariables), function(j) {
      lapply(values[[j]], function(value) {
        list(columns = xind[j], values = value)
      })
    }),
    recursive = FALSE
  )
  draws <- pdbart.drawsAt(sampler, rows, settings)
  fit <- pdbart.packageRun(sampler, fit, draws$samples)
  split <- pdbart.chainsSplit(fit, n.chains)
  fdr <- lapply(seq_len(numVariables), function(j) {
    fd <- draws$fd[, variable == j, drop = FALSE]
    if (split) pdbart.splitChains(fd, n.chains) else fd
  })

  result <- pdbart.buildResult(
    sampler,
    fit,
    fdr,
    levs,
    xind,
    keepSampler,
    "pdbart"
  )

  if (pl) {
    plot(result, plquants = plquants)
  }

  result
}

pd2bart <- function(
  formula,
  data,
  xind = NULL,
  levs = NULL,
  levquants = c(0.05, seq(0.1, 0.9, 0.1), 0.95),
  pl = TRUE,
  plquants = c(0.05, 0.95),
  ...
) {
  matchedCall <- match.call()
  callingEnv <- parent.frame()

  sampler <- fit <- keepSampler <- NULL ## for R CMD check (massign assigns)
  dataGiven <- !missing(data)
  massign[sampler, fit, keepSampler] <- pdbart.prologue(
    pdbart.argument(
      !missing(formula),
      formula,
      "x.train",
      matchedCall,
      callingEnv
    ),
    function() {
      pdbart.argument(dataGiven, data, "y.train", matchedCall, callingEnv)
    },
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

  xValues <- as.matrix(expand.grid(values[[1L]], values[[2L]]))
  rows <- pdbart.averagedRows(sampler, x)
  n.chains <- sampler$control@n.chains

  # with two predictors and one offset for every row, each grid point is a
  # whole row and every averaged row is that row, so its prediction is the
  # average
  offset <- rows$offset
  if (ncol(x) == 2L && (is.null(offset) || all(offset == offset[1L]))) {
    gridRows <- if (xind[1L] < xind[2L]) xValues else xValues[, c(2L, 1L)]
    colnames(gridRows) <- colnames(x)
    if (sampler$control@keepTrees) {
      samples <- NULL
      fdr <- pdbart.drawsByRow(
        if (is.null(offset)) {
          sampler$predict(gridRows)
        } else {
          sampler$predict(gridRows, rep_len(offset[1L], nrow(gridRows)))
        }
      )
    } else {
      sampler$setTestPredictor(gridRows)
      samples <- sampler$run(0L, sampler$control@n.samples)
      fdr <- pdbart.drawsByRow(samples$test) +
        if (is.null(offset)) 0 else offset[1L]
    }
  } else {
    settings <- lapply(seq_len(nrow(xValues)), function(i) {
      list(columns = xind[1:2], values = xValues[i, ])
    })
    draws <- pdbart.drawsAt(sampler, rows, settings)
    samples <- draws$samples
    fdr <- draws$fd
  }
  fit <- pdbart.packageRun(sampler, fit, samples)
  if (pdbart.chainsSplit(fit, n.chains)) {
    fdr <- pdbart.splitChains(fdr, n.chains)
  }

  result <- pdbart.buildResult(
    sampler,
    fit,
    fdr,
    levs,
    xind,
    keepSampler,
    "pd2bart"
  )

  if (pl) {
    plot(result, plquants = plquants)
  }

  result
}
