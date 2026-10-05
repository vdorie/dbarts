# The arguments pdbart and pd2bart keep for themselves; everything else in the
# call is bart's.
pdbart.ownArgs <- c(
  "xind",
  "levs",
  "levquants",
  "pl",
  "plquants",
  "type",
  "newdata",
  "n.average.rows",
  "average.weights"
)

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
      c("aft", "hazard", "hazard.probit", "hazard.logistic")
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
    return(list(
      sampler = pdbart.rowSampler(fit),
      fit = fit,
      keepSampler = keepSampler,
      isSampler = FALSE,
      getData = getData
    ))
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
    caller,
    fits = FALSE
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
    return(list(
      sampler = object,
      fit = list(),
      keepSampler = keepSampler,
      isSampler = TRUE
    ))
  }

  fit <- object
  pdbart.refuseFamily(pdbart.fitFamily(fit), caller)
  sampler <- pdbart.rowSampler(fit)
  if (is.null(sampler) || !sampler$control@keepTrees) {
    fit <- pdbart.refit(fit, callingEnv, caller)
  }
  list(
    sampler = pdbart.rowSampler(fit),
    fit = fit,
    keepSampler = keepSampler,
    isSampler = FALSE,
    getData = pdbart.storedData(fit, callingEnv)
  )
}

# The data a fit's stored call names, re-evaluated where pdbart was called,
# as list(value); NULL when the call names none, and FALSE when there is no
# call or it cannot be evaluated.
pdbart.storedData <- function(fit, callingEnv) {
  function() {
    call <- fit$call
    if (
      !is.call(call) ||
        identical(call, call("NA")) ||
        identical(call, call("NULL"))
    ) {
      return(FALSE)
    }
    name <- if ("x.train" %in% names(call)) "y.train" else "data"
    if (is.null(call[[name]])) {
      return(NULL)
    }
    tryCatch(
      list(eval(call[[name]], callingEnv)),
      error = function(e) FALSE
    )
  }
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

# The sampler whose rows a fit was made on: a hurdle fit's zero part, whose
# rows are every observation, or any other fit's own.
pdbart.rowSampler <- function(fit) {
  if (inherits(fit, "bartHurdle")) fit$zero$fit else fit$fit
}

pdbart.isFormulaFit <- function(sampler) {
  !is.null(attr(sampler$data@x, "terms"))
}

# The variables of a formula fit's right-hand side, as its data names them.
pdbart.formulaVariables <- function(sampler) {
  labels <- attr(attr(sampler$data@x, "terms"), "term.labels")
  unique(unlist(lapply(labels, function(label) all.vars(str2lang(label)))))
}

# The values 'type' can take before the fit is known; each is checked against
# the fit once there is one.
pdbart.types <- c("auto", "bart", "link", "ev", "response", "ppd", "prob")
pdbart.types <- c(pdbart.types, "sigma")

# The type a fit's averages are taken on: "auto" is the link scale, except
# the mean response on a hurdle fit. A sampler carries no fit to transform
# through, so it takes the link scale only.
pdbart.resolveType <- function(type, fit, isSampler, caller) {
  type <- foldTypeAliases(type)
  if (isSampler) {
    if (type %not_in% c("auto", "bart")) {
      stop(
        "a sampler passed to '",
        caller,
        "' takes only type = \"bart\": it carries no fit to transform ",
        "through",
        call. = FALSE
      )
    }
    return("bart")
  }
  hurdle <- inherits(fit, "bartHurdle")
  if (type == "auto") {
    return(if (hurdle) "ev" else "bart")
  }
  allowed <- if (hurdle) {
    c("ev", "ppd", "prob", "bart")
  } else if (inherits(fit, "bartNegbin")) {
    c("ev", "ppd", "bart")
  } else {
    c("ev", "ppd", "bart", "sigma")
  }
  if (type %not_in% allowed) {
    stop(
      "'",
      caller,
      "' on a ",
      pdbart.fitFamily(fit),
      " fit does not take type = \"",
      type,
      "\"; it takes ",
      quotedNameList(allowed),
      call. = FALSE
    )
  }
  if (type == "sigma" && !fitIsHeteroscedastic(fit)) {
    stop(
      "type = \"sigma\" needs a fit with a variance forest, whose residual ",
      "spread varies by row",
      call. = FALSE
    )
  }
  type
}

# The averaging arguments, checked before anything is fit.
pdbart.checkAveraging <- function(type, newdata, n.average.rows, caller) {
  if (!is.character(type) || length(type) != 1L || is.na(type)) {
    stop("'type' must be a single string", call. = FALSE)
  }
  if (type == "forest") {
    stop(
      "'",
      caller,
      "' does not take type = \"forest\", which reports each forest of a ",
      "fit apart; a several-forest fit's partial dependence is on its ",
      "combined prediction",
      call. = FALSE
    )
  }
  if (type %not_in% pdbart.types) {
    stop(
      "'type' must be one of ",
      quotedNameList(pdbart.types),
      call. = FALSE
    )
  }
  if (!is.null(n.average.rows)) {
    if (!is.null(newdata)) {
      stop(
        "'n.average.rows' samples the fit's own rows, so it cannot be given ",
        "with 'newdata'; subsample 'newdata' instead",
        call. = FALSE
      )
    }
    if (
      !is.numeric(n.average.rows) ||
        length(n.average.rows) != 1L ||
        is.na(n.average.rows) ||
        n.average.rows < 1 ||
        n.average.rows != round(n.average.rows)
    ) {
      stop("'n.average.rows' must be a positive whole number", call. = FALSE)
    }
  }
  invisible(NULL)
}

# Resolve 'xind' to the varied predictors: in a formula fit the names of
# variables of the data, otherwise column indices into the predictor matrix.
# 'xind' is forwarded as a promise so that a formula-style expression, which
# names them as terms, is read unevaluated.
pdbart.resolveXind <- function(xind, matchedCall, sampler, formulaFit) {
  available <- if (formulaFit) {
    pdbart.formulaVariables(sampler)
  } else {
    colnames(sampler$data@x)
  }
  tryResult <- tryCatch(xind, error = I)

  if (inherits(tryResult, "error")) {
    formula <- ~a
    formula[[2L]] <- matchedCall[["xind"]]
    xind <- attr(terms(formula), "term.labels")
  } else if (
    is.character(xind) &&
      length(xind) == 1L &&
      xind %not_in% available
  ) {
    formula <- ~a
    formula[[2L]] <- parse(text = xind)[[1L]]
    xind <- attr(terms(formula), "term.labels")
  } else if (is.null(xind)) {
    xind <- if (formulaFit) available else seq_len(ncol(sampler$data@x))
  }

  if (formulaFit && !is.character(xind)) {
    stop(
      "in a formula fit 'xind' names variables of the data; a column number ",
      "of the model matrix is not taken",
      call. = FALSE
    )
  }
  if (is.character(xind)) {
    if (is.null(available)) {
      stop(
        "passing 'xind' by name requires the predictors to have column names"
      )
    }
    unknown <- xind %not_in% available
    if (any(unknown)) {
      stop(
        "unrecognized ",
        if (formulaFit) "variables" else "columns",
        " '",
        paste0(xind[unknown], collapse = "', '"),
        "'"
      )
    }
    if (!formulaFit) {
      xind <- match(xind, available)
    }
  }

  xind
}

# How the rows name each varied predictor: a formula fit's variable names, or
# a matrix fit's column names when it has them and its column indices
# otherwise.
pdbart.keys <- function(sampler, xind, formulaFit) {
  if (formulaFit || is.null(colnames(sampler$data@x))) {
    xind
  } else {
    colnames(sampler$data@x)[xind]
  }
}

pdbart.xLabels <- function(sampler, xind, formulaFit = FALSE) {
  if (formulaFit) {
    xind
  } else if (is.null(colnames(sampler$data@x))) {
    paste0("x", xind)
  } else {
    colnames(sampler$data@x)[xind]
  }
}

# The level table of each varied predictor that is a factor, NULL for any
# other: the fit's own when the predictor is one of its columns, otherwise
# the levels the grid's source rows carry.
pdbart.factorLevels <- function(sampler, xind, formulaFit, source) {
  factorLevels <- attr(sampler$data@x, "factor.levels")
  lapply(xind, function(key) {
    j <- if (formulaFit) match(key, colnames(sampler$data@x)) else key
    if (!is.na(j) && !is.null(factorLevels) && !is.null(factorLevels[[j]])) {
      return(factorLevels[[j]])
    }
    if (!formulaFit) {
      return(NULL)
    }
    column <- source[[key]]
    if (is.factor(column)) {
      levels(column)
    } else if (is.character(column)) {
      sort(unique(column[!is.na(column)]))
    }
  })
}

pdbart.column <- function(rows, key) {
  if (is.data.frame(rows)) {
    if (key %not_in% names(rows)) {
      stop("'newdata' has no column '", key, "'", call. = FALSE)
    }
    rows[[key]]
  } else {
    rows[, key]
  }
}

# Default the 'levs' list from the grid's source columns: for a factor, its
# levels, by name; otherwise either the sorted unique values (when there are
# too few to bin) or the unique quantiles at 'levquants', missing values left
# out of both. 'cmp' is the comparison deciding "too few": pdbart uses `<`,
# pd2bart uses `<=` (a long-standing difference in the two entry points,
# preserved here rather than reconciled).
pdbart.defaultLevs <- function(columns, levquants, cmp, levels) {
  lapply(seq_along(columns), function(j) {
    if (!is.null(levels[[j]])) {
      return(levels[[j]])
    }
    column <- columns[[j]]
    column <- column[!is.na(column)]
    uniqueValues <- unique(column)
    if (cmp(length(uniqueValues), length(levquants))) {
      sort(uniqueValues)
    } else {
      unique(quantile(column, probs = levquants))
    }
  })
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

# The value a predictor column of coded rows takes at one 'levs' entry: a
# factor level's 0-based code, any other column's value itself. Rows in a
# data frame take the level's name.
pdbart.levelValues <- function(levs, levels, rows) {
  lapply(seq_along(levs), function(j) {
    if (is.null(levels[[j]]) || is.data.frame(rows)) {
      levs[[j]]
    } else {
      match(levs[[j]], levels[[j]]) - 1
    }
  })
}

# 'rows' with each of 'keys' set to the matching entry of 'values', one value
# for every row or one per row. A factor keeps its levels, gaining the value
# when a subgroup's rows lack it.
pdbart.setVariables <- function(rows, keys, values) {
  for (k in seq_along(keys)) {
    value <- values[[k]]
    if (is.data.frame(rows)) {
      column <- rows[[keys[k]]]
      rows[[keys[k]]] <- if (is.factor(column)) {
        factor(
          rep_len(as.character(value), nrow(rows)),
          levels = union(levels(column), as.character(value)),
          ordered = is.ordered(column)
        )
      } else {
        rep_len(value, nrow(rows))
      }
    } else {
      rows[, keys[k]] <- value
    }
  }
  rows
}

# A formula fit's rows as the variables of its data, collected as
# get_all_vars collects them and cut to the rows the fit was made on by their
# names. 'getData' gives the data as list(value), NULL when the call named
# none, or FALSE when it cannot be had.
pdbart.trainingRows <- function(sampler, getData, caller) {
  refuse <- function() {
    stop(
      "'",
      caller,
      "' reads a formula fit's rows from the data its call names, which ",
      "cannot be evaluated here; give the rows as 'newdata'",
      call. = FALSE
    )
  }
  data <- getData()
  if (isFALSE(data)) {
    refuse()
  }
  terms <- attr(sampler$data@x, "terms")
  rows <- tryCatch(
    if (is.null(data)) {
      stats::get_all_vars(terms)
    } else {
      stats::get_all_vars(terms, as.data.frame(data[[1L]]))
    },
    error = function(e) NULL
  )
  rowNames <- sampler$data@rowNames$train
  index <- if (!is.null(rows)) match(rowNames, rownames(rows))
  if (
    is.null(rows) ||
      length(rowNames) != NROW(sampler$data@y) ||
      anyNA(index)
  ) {
    refuse()
  }
  rows[index, , drop = FALSE]
}

# 'average.weights' checked and kept as given, or NULL when absent.
pdbart.averageWeights <- function(weights, n, rows) {
  if (is.null(weights)) {
    return(NULL)
  }
  if (!is.numeric(weights) || length(weights) != n) {
    stop(
      "'average.weights' must have one number per row of ",
      rows,
      ", ",
      n,
      call. = FALSE
    )
  }
  if (anyNA(weights) || any(!is.finite(weights)) || any(weights < 0)) {
    stop("'average.weights' must be finite and non-negative", call. = FALSE)
  }
  if (all(weights == 0)) {
    stop("'average.weights' must not all be zero", call. = FALSE)
  }
  as.double(weights)
}

# The offset predict is to add to a subset of the fit's rows, or NULL when
# predict evaluates the fit's own on them. An offset argument that cannot be
# evaluated on these rows was given for the training rows alone, as a plain
# vector, and each row takes its own stored share; the offset() terms of a
# formula are evaluated by predict on the rows in either case.
pdbart.storedOffset <- function(sampler, all, rows, index) {
  argument <- attr(sampler$data, "offset.argument")
  if (is.null(argument) || !isFALSE(evaluateOffsetArgument(argument, rows))) {
    return(NULL)
  }
  share <- evaluateOffsetArgument(argument, all)
  if (isFALSE(share)) {
    termOffset <- if (pdbart.isFormulaFit(sampler)) {
      formulaTermOffset(sampler$data@x, all, "rows")
    }
    share <- sampler$data@offset - if (is.null(termOffset)) 0 else termOffset
  }
  if (length(share) == 1L) share else share[index]
}

# The rows averaged over and the weight each gets (NULL for equal weights):
# 'newdata' as given, or the fit's own rows less those it weights 0 or masks
# out, subsampled when asked. 'source' is where the default grid is read:
# 'newdata', or every row the fit was made on.
pdbart.frame <- function(
  sampler,
  formulaFit,
  getData,
  newdata,
  n.average.rows,
  average.weights,
  caller
) {
  if (!is.null(newdata)) {
    if (formulaFit && !is.data.frame(newdata)) {
      stop("'newdata' for a formula fit must be a data frame", call. = FALSE)
    }
    argument <- attr(sampler$data, "offset.argument")
    if (isFALSE(evaluateOffsetArgument(argument, newdata))) {
      stop(
        "the fit's 'offset' was given as ",
        describeOffsetArgument(argument),
        ", which cannot be evaluated on the rows of 'newdata'; write the ",
        "offset as a column of the data",
        call. = FALSE
      )
    }
    weights <- pdbart.averageWeights(
      average.weights,
      NROW(newdata),
      "'newdata'"
    )
    return(list(
      rows = newdata,
      weights = if (!is.null(weights)) weights / sum(weights),
      offset = NULL,
      source = newdata
    ))
  }

  all <- if (formulaFit) {
    pdbart.trainingRows(sampler, getData, caller)
  } else {
    extract(sampler, "predictors")
  }
  data <- sampler$data
  keep <- rep_len(TRUE, NROW(all))
  if (length(data@weights) > 0L) {
    keep <- keep & data@weights > 0
  }
  if (!is.null(sampler$activeRows)) {
    keep <- keep & sampler$activeRows != 0
  }
  weights <- pdbart.averageWeights(
    average.weights,
    NROW(all),
    "the fit"
  )
  if (!is.null(weights)) {
    keep <- keep & weights > 0
    if (!any(keep)) {
      stop(
        "every row 'average.weights' gives a positive weight has a fit ",
        "weight of 0",
        call. = FALSE
      )
    }
  }
  index <- which(keep)
  if (!is.null(n.average.rows) && n.average.rows < length(index)) {
    # kept in the fit's order, so a plain-vector offset stays aligned
    index <- index[sort(sample.int(length(index), n.average.rows))]
  }
  rows <- all[index, , drop = FALSE]
  if (!is.null(weights)) {
    weights <- weights[index] / sum(weights[index])
  }
  list(
    rows = rows,
    weights = weights,
    offset = pdbart.storedOffset(sampler, all, rows, index),
    source = all,
    index = index
  )
}

# Per-draw predictions of a fit at rows, draws x rows, the chains merged in
# turn.
pdbart.predictDraws <- function(fit, rows, type, offset) {
  pred <- predict(fit, rows, type = type, offset = offset)
  if (is.null(dim(pred))) matrix(pred, ncol = NROW(rows)) else pred
}

pdbart.average <- function(pred, weights) {
  if (is.null(weights)) rowMeans(pred) else drop(pred %*% weights)
}

# Per-draw averages over the frame's rows at each of 'settings', draws x
# settings, a setting being list(keys, values).
pdbart.fitDrawsAt <- function(fit, frame, type, settings) {
  fd <- NULL
  for (i in seq_along(settings)) {
    rows <- pdbart.setVariables(
      frame$rows,
      settings[[i]]$keys,
      settings[[i]]$values
    )
    draws <- pdbart.average(
      pdbart.predictDraws(fit, rows, type, frame$offset),
      frame$weights
    )
    if (is.null(fd)) {
      fd <- matrix(NA_real_, length(draws), length(settings))
    }
    fd[, i] <- draws
  }
  fd
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

# A sampler's rows averaged over: its own, less those it gives a 0 weight or
# masks out, each with its stored offset (NULL when it has none).
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

# Per-draw averages of a sampler's predictions over 'rows' at each of
# 'settings', draws x settings. Each row's offset enters its prediction. A
# sampler with saved trees predicts from them; one without runs once over
# every setting's rows stacked, changing its state, and its samples come back
# for the result.
pdbart.drawsAt <- function(sampler, rows, settings) {
  n.chains <- sampler$control@n.chains
  fd <- matrix(
    NA_real_,
    sampler$control@n.samples * n.chains,
    length(settings)
  )
  setRow <- function(setting) {
    pdbart.setVariables(rows$x, setting$keys, setting$values)
  }
  if (sampler$control@keepTrees) {
    for (i in seq_along(settings)) {
      x.test <- setRow(settings[[i]])
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
  sampler$setTestPredictor(do.call(rbind, lapply(settings, setRow)))
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
  varcount <- if (inherits(fit, "bartHurdle")) {
    fit$zero$varcount
  } else {
    fit$varcount
  }
  n.chains > 1L && length(dim(varcount)) == 3L
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

# What pdbart and pd2bart share once the fit is known: the scale, the varied
# predictors and the rows averaged over.
pdbart.setup <- function(
  prologue,
  xind,
  matchedCall,
  type,
  newdata,
  n.average.rows,
  average.weights,
  caller
) {
  sampler <- prologue$sampler
  isSampler <- prologue$isSampler
  if (
    isSampler &&
      (!is.null(newdata) ||
        !is.null(n.average.rows) ||
        !is.null(average.weights))
  ) {
    stop(
      "a sampler passed to '",
      caller,
      "' averages over its own rows; 'newdata', 'n.average.rows' and ",
      "'average.weights' take a fit",
      call. = FALSE
    )
  }
  formulaFit <- !isSampler && pdbart.isFormulaFit(sampler)
  xind <- pdbart.resolveXind(xind, matchedCall, sampler, formulaFit)
  list(
    sampler = sampler,
    fit = prologue$fit,
    isSampler = isSampler,
    type = pdbart.resolveType(type, prologue$fit, isSampler, caller),
    formulaFit = formulaFit,
    xind = xind,
    keys = pdbart.keys(sampler, xind, formulaFit),
    xLabels = pdbart.xLabels(sampler, xind, formulaFit),
    frame = if (isSampler) {
      list(
        rows = NULL,
        source = extract(sampler, "predictors")
      )
    } else {
      pdbart.frame(
        sampler,
        formulaFit,
        prologue$getData,
        newdata,
        n.average.rows,
        average.weights,
        caller
      )
    }
  )
}

# The grid of each varied predictor, given or by default, and the levels of
# those that are factors.
pdbart.grid <- function(setup, levs, levquants, numVariables, cmp) {
  source <- setup$frame$source
  levels <- pdbart.factorLevels(
    setup$sampler,
    setup$xind[seq_len(numVariables)],
    setup$formulaFit,
    source
  )
  labels <- setup$xLabels[seq_len(numVariables)]
  if (is.null(levs)) {
    columns <- lapply(setup$keys[seq_len(numVariables)], function(key) {
      pdbart.column(source, key)
    })
    levs <- pdbart.defaultLevs(columns, levquants, cmp, levels)
  } else if (length(levs) != numVariables) {
    stop("'levs' must have the same length as 'xind'")
  } else {
    levs <- pdbart.checkLevs(levs, levels, labels)
  }
  namedList(levs, levels)
}

# Per-draw averages at each setting, draws x settings, and the samples of a
# sampler run without saved trees.
pdbart.draws <- function(setup, settings) {
  if (setup$isSampler) {
    rows <- pdbart.averagedRows(setup$sampler, setup$frame$source)
    return(pdbart.drawsAt(setup$sampler, rows, settings))
  }
  list(
    fd = pdbart.fitDrawsAt(setup$fit, setup$frame, setup$type, settings),
    samples = NULL
  )
}

# Assemble the returned pdbart/pd2bart result list. Identical between the two
# entry points except for the S3 class stamped on it ('className').
pdbart.buildResult <- function(setup, fit, fdr, levs, keepSampler, className) {
  sampler <- setup$sampler
  bartcall <- if (is.null(fit$call)) sampler$control@call else fit$call
  hurdle <- inherits(fit, "bartHurdle")
  copied <- if (hurdle) {
    character()
  } else if (sampler$control@binary) {
    "yhat.train"
  } else {
    c("yhat.train", "first.sigma", "sigma", "yhat.train.mean", "sigest")
  }
  result <- c(
    list(
      fd = fdr,
      levs = levs,
      xlbs = setup$xLabels[seq_along(levs)],
      bartcall = bartcall
    ),
    fit[intersect(copied, names(fit))],
    list(
      y = if (is.null(fit$y)) sampler$data@y else fit$y,
      n.chains = sampler$control@n.chains,
      type = setup$type,
      family = if (setup$isSampler) {
        pdbart.samplerFamily(sampler)
      } else {
        pdbart.fitFamily(fit)
      },
      fit = if (hurdle) {
        list(zero = fit$zero$fit, positive = fit$positive$fit)
      } else {
        sampler
      }
    )
  )
  if (setup$isSampler && "sigest" %in% copied && is.null(result$sigest)) {
    result$sigest <- sampler$data@sigma
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
  type = "auto",
  newdata = NULL,
  n.average.rows = NULL,
  average.weights = NULL,
  ...
) {
  matchedCall <- match.call()
  callingEnv <- parent.frame()
  pdbart.checkAveraging(type, newdata, n.average.rows, "pdbart")

  dataGiven <- !missing(data)
  prologue <- pdbart.prologue(
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
  setup <- pdbart.setup(
    prologue,
    xind,
    matchedCall,
    type,
    newdata,
    n.average.rows,
    average.weights,
    "pdbart"
  )
  numVariables <- length(setup$xind)
  massign[levs, levels] <- pdbart.grid(
    setup,
    levs,
    levquants,
    numVariables,
    `<`
  )
  rows <- if (setup$isSampler) setup$frame$source else setup$frame$rows
  values <- pdbart.levelValues(levs, levels, rows)

  # every variable's settings in one pass, so that a sampler run without
  # saved trees draws them all from one run
  variable <- rep(seq_len(numVariables), lengths(values))
  settings <- unlist(
    lapply(seq_len(numVariables), function(j) {
      lapply(values[[j]], function(value) {
        list(keys = setup$keys[j], values = list(value))
      })
    }),
    recursive = FALSE
  )
  draws <- pdbart.draws(setup, settings)
  fit <- pdbart.packageRun(setup$sampler, prologue$fit, draws$samples)
  n.chains <- setup$sampler$control@n.chains
  split <- pdbart.chainsSplit(fit, n.chains)
  fdr <- lapply(seq_len(numVariables), function(j) {
    fd <- draws$fd[, variable == j, drop = FALSE]
    if (split) pdbart.splitChains(fd, n.chains) else fd
  })

  result <- pdbart.buildResult(
    setup,
    fit,
    fdr,
    levs,
    prologue$keepSampler,
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
  type = "auto",
  newdata = NULL,
  n.average.rows = NULL,
  average.weights = NULL,
  ...
) {
  matchedCall <- match.call()
  callingEnv <- parent.frame()
  pdbart.checkAveraging(type, newdata, n.average.rows, "pd2bart")

  dataGiven <- !missing(data)
  prologue <- pdbart.prologue(
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
  sampler <- prologue$sampler

  # with two predictors and one offset for every row, each grid point is a
  # whole row and every averaged row is that row, so its prediction is the
  # average; a posterior predictive draw is the exception, each row drawing
  # its own noise
  offset <- sampler$data@offset
  numPredictors <- if (!prologue$isSampler && pdbart.isFormulaFit(sampler)) {
    length(pdbart.formulaVariables(sampler))
  } else {
    ncol(sampler$data@x)
  }
  shortcut <- numPredictors == 2L &&
    (length(offset) == 0L || all(offset == offset[1L])) &&
    foldTypeAliases(type) != "ppd"
  if (shortcut) {
    unused <- c(
      if (!is.null(newdata)) "newdata",
      if (!is.null(n.average.rows)) "n.average.rows",
      if (!is.null(average.weights)) "average.weights"
    )
    if (length(unused) > 0L) {
      warning(
        quotedNameList(unused),
        " ha",
        if (length(unused) > 1L) "ve" else "s",
        " no effect on pd2bart over a fit's two predictors with one offset ",
        "for every row, where each grid point is a whole row",
        call. = FALSE
      )
    }
    newdata <- n.average.rows <- average.weights <- NULL
  }
  setup <- pdbart.setup(
    prologue,
    xind,
    matchedCall,
    type,
    newdata,
    n.average.rows,
    average.weights,
    "pd2bart"
  )
  massign[levs, levels] <- pdbart.grid(setup, levs, levquants, 2L, `<=`)
  rows <- if (setup$isSampler) setup$frame$source else setup$frame$rows
  values <- pdbart.levelValues(levs, levels, rows)
  grid <- expand.grid(values[[1L]], values[[2L]], stringsAsFactors = FALSE)
  keys <- setup$keys[1:2]
  n.chains <- sampler$control@n.chains

  samples <- NULL
  if (shortcut) {
    if (setup$isSampler) {
      gridRows <- pdbart.setVariables(
        setup$frame$source[rep_len(1L, nrow(grid)), , drop = FALSE],
        keys,
        grid
      )
      if (sampler$control@keepTrees) {
        fdr <- pdbart.drawsByRow(
          if (length(offset) == 0L) {
            sampler$predict(gridRows)
          } else {
            sampler$predict(gridRows, rep_len(offset[1L], nrow(gridRows)))
          }
        )
      } else {
        sampler$setTestPredictor(gridRows)
        samples <- sampler$run(0L, sampler$control@n.samples)
        fdr <- pdbart.drawsByRow(samples$test) +
          if (length(offset) == 0L) 0 else offset[1L]
      }
    } else {
      frame <- setup$frame
      first <- rep_len(1L, nrow(grid))
      gridRows <- pdbart.setVariables(
        frame$rows[first, , drop = FALSE],
        keys,
        grid
      )
      fdr <- pdbart.predictDraws(
        setup$fit,
        gridRows,
        setup$type,
        pdbart.storedOffset(
          sampler,
          frame$source,
          gridRows,
          frame$index[first]
        )
      )
    }
  } else {
    settings <- lapply(seq_len(nrow(grid)), function(i) {
      list(keys = keys, values = list(grid[[1L]][i], grid[[2L]][i]))
    })
    draws <- pdbart.draws(setup, settings)
    samples <- draws$samples
    fdr <- draws$fd
  }
  fit <- pdbart.packageRun(sampler, prologue$fit, samples)
  if (pdbart.chainsSplit(fit, n.chains)) {
    fdr <- pdbart.splitChains(fdr, n.chains)
  }

  result <- pdbart.buildResult(
    setup,
    fit,
    fdr,
    levs,
    prologue$keepSampler,
    "pd2bart"
  )

  if (pl) {
    plot(result, plquants = plquants)
  }

  result
}
