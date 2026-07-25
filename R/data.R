ORDINAL_VARIABLE <- 0L
CATEGORICAL_VARIABLE <- 1L

methods::setMethod(
  "initialize",
  "dbartsData",
  function(.Object, modelMatrices, n.cuts = 100L, sigma = NA_real_) {
    if (!missing(modelMatrices)) {
      .Object@y <- modelMatrices$y
      .Object@x <- modelMatrices$x
      # makeCategoricalModelMatrix types its columns; everything else is ordinal
      .Object@varTypes <- if (!is.null(attr(.Object@x, "varTypes"))) {
        as.integer(attr(.Object@x, "varTypes"))
      } else {
        rep.int(ORDINAL_VARIABLE, ncol(.Object@x))
      }
      .Object@x.test <- modelMatrices$x.test
      .Object@weights <- modelMatrices$weights
      .Object@weights.test <- modelMatrices$weights.test
      .Object@offset <- modelMatrices$offset
      .Object@offset.test <- modelMatrices$offset.test

      .Object@testUsesRegularOffset <- modelMatrices$testUsesRegularOffset
    }

    .Object@n.cuts <- rep_len(as.integer(n.cuts), ncol(.Object@x))
    .Object@sigma <- sigma

    validObject(.Object)
    .Object
  }
)

makeTestModelMatrix <- function(data, newdata) {
  validateXTest(newdata, data@x)
}

validateXTest <- function(x.test, x.train) {
  termLabels <- attr(x.train, "term.labels")
  numPredictors <- ncol(x.train)
  predictorNames <- colnames(x.train)
  drop <- attr(x.train, "drop")
  factorLevels <- attr(x.train, "factor.levels")

  if (is.null(x.test)) {
    return(x.test)
  }
  if (is.numeric(x.test) && is.null(dim(x.test)) && length(x.test) > 0L) {
    x.test <- matrix(x.test, ncol = length(x.test))
  }
  if (is.numeric(x.test) && NCOL(x.test) == 0L) {
    return(NULL)
  }
  if (is.data.frame(x.test)) {
    if (any(vapply(x.test, isSparseDataFrameColumn, FALSE))) {
      # sparse columns ride to the engine unexpanded, coded over the training
      # level table; the resulting container is preserved below (the model
      # frame replay takes no S4 columns, so it is skipped here)
      if (is.null(factorLevels)) {
        stop(
          "sparse test predictor columns require a categorical training ",
          "design; supply 'x' through the x/y interface"
        )
      }
      x.test <- mapFactorColumnsToTrainingLevels(
        x.test,
        predictorNames,
        factorLevels
      )
      x.test <- makeCategoricalModelMatrix(x.test)
    } else {
      if (!is.null(termLabels)) {
        testFormula <-
          as.formula(paste("~", paste(termLabels, collapse = " + ")))
        # model.frame resolves an absent term in the enclosing scope, so a
        # predictor missing from newdata that shares a name with a base object
        # (e.g. 'c') silently binds to it and fails with an opaque
        # "invalid type (builtin)"; name the missing variables up front instead
        neededVars <- all.vars(testFormula)
        missingVars <- neededVars[neededVars %not_in% names(x.test)]
        if (length(missingVars) > 0L) {
          stop(
            "'test' data is missing ",
            if (length(missingVars) > 1L) "variables" else "variable",
            " required by the model: '",
            toString(missingVars),
            "'"
          )
        }
        x.test <- model.frame(
          formula = testFormula,
          data = x.test,
          na.action = stats::na.pass
        )
      }
      if (!is.null(factorLevels)) {
        # trained with factors unexpanded: code against the training levels
        x.test <- mapFactorColumnsToTrainingLevels(
          x.test,
          predictorNames,
          factorLevels
        )
        x.test <- makeCategoricalModelMatrix(x.test)
      } else {
        x.test <- makeModelMatrixFromDataFrame(
          x.test,
          if (!is.null(drop)) drop else TRUE
        )
      }
    }
  }
  # a bare dgCMatrix test set takes the same resident path a mixed-container
  # test set's sparse columns already do (below), symmetric with a bare
  # dgCMatrix train set: wrap it as an all-sparse mixed container rather than
  # densifying. A bare numeric sparse matrix carries no factor levels, so it
  # cannot supply a categorical training column's values - refuse informatively
  # rather than let the bridge reject a malformed container.
  if (inherits(x.test, "dgCMatrix")) {
    if (!is.null(factorLevels)) {
      stop(
        "sparse (dgCMatrix) test predictors cannot supply values for a ",
        "categorical training column; supply 'test' as a dense matrix or ",
        "data frame instead"
      )
    }
    x.test <- wrapSparseTestMatrix(x.test)
  }
  # a sparse-backed container stays resident (the engine codes it against the
  # training cuts); everything else densifies as before
  xTestIsSparseContainer <-
    inherits(x.test, "dbartsMixedMatrix") && predictorSourceIsSparse(x.test)
  if (!is.matrix(x.test) && !xTestIsSparseContainer) {
    x.test <- as.matrix(x.test)
  }

  if (!xTestIsSparseContainer) {
    if (!is.numeric(x.test)) {
      stop("test matrix must be numeric")
    }

    if (is.integer(x.test)) {
      x.test <- matrix(as.double(x.test), nrow(x.test))
    }
  }

  if (!identical(NCOL(x.test), numPredictors)) {
    stop("number of columns in 'test' must be equal to that of 'x'")
  }
  if (numPredictors > 1) {
    xIsNamed <- !is.null(predictorNames)
    testIsNamed <- !is.null(colnames(x.test))

    columnIndices <- seq.int(numPredictors)
    if (xIsNamed && !testIsNamed) {
      # named fit, positional (unnamed) test: swapped columns would match
      # silently and return badly wrong numbers, so spell out the mapping the
      # positional match assumes rather than only noting that it happened
      shown <- min(numPredictors, 3L)
      mapping <- paste0(
        "column ",
        seq_len(shown),
        " = '",
        predictorNames[seq_len(shown)],
        "'",
        collapse = ", "
      )
      warning(
        "'test' is unnamed but 'x' had named predictors; columns of 'test' ",
        "are matched to 'x' by position (",
        mapping,
        if (numPredictors > shown) ", ..." else "",
        "). Supply 'test' with column names to match by name instead."
      )
    } else if (
      (!xIsNamed && testIsNamed) ||
        length(unique(predictorNames)) != length(predictorNames)
    ) {
      warning(
        "'x' and 'test' are not both named; columns of 'test' will be matched by position"
      )
    } else if (xIsNamed && testIsNamed) {
      matchIndices <- match(predictorNames, colnames(x.test))
      if (any(is.na(matchIndices))) {
        stop(
          "column names of 'test' do not match those of 'x': '",
          toString(predictorNames[is.na(matchIndices)]),
          "' present in 'x' but not in 'test' (whose columns are '",
          toString(colnames(x.test)),
          "')"
        )
      } else {
        columnIndices <- matchIndices
      }
    }

    if (xTestIsSparseContainer) {
      # reorder columns without densifying: the map re-points each predictor
      # position at its source, the names follow
      x.test$map <- x.test$map[columnIndices]
      if (!is.null(x.test$columnNames)) {
        x.test$columnNames <- x.test$columnNames[columnIndices]
      }
      if (xIsNamed) x.test$columnNames <- predictorNames
    } else {
      x.test <- x.test[, columnIndices, drop = FALSE]
      if (xIsNamed) colnames(x.test) <- predictorNames
    }
  }

  x.test
}

findTermInFormulaData <- function(formula, data, term) {
  dataIsMissing <- missing(data)
  matchedCall <- match.call()

  if (is.numeric(matchedCall$term)) {
    return(term)
  }

  if (!dataIsMissing) {
    if (is.symbol(matchedCall$term)) {
      if (any(names(data) == as.character(matchedCall$term))) {
        return(data[[as.character(matchedCall$term)]])
      }
    } else if (is.language(matchedCall$term)) {
      tryResult <- with(
        data,
        tryCatch(eval(matchedCall$term), error = function(e) e)
      )
      if (!inherits(tryResult, "error")) return(tryResult)
    }
  }
  if (is.symbol(matchedCall$term)) {
    if (any(ls(environment(formula)) == as.character(matchedCall$term))) {
      return(get(as.character(matchedCall$term), envir = environment(formula)))
    }
    tryResult <- tryCatch(
      get(as.character(matchedCall$term)),
      error = function(e) e
    )
    if (!inherits(tryResult, "error") && !is.null(tryResult)) return(tryResult)
  } else if (is.language(matchedCall$term)) {
    tryResult <- tryCatch(
      eval(matchedCall$term, environment(formula)),
      error = function(e) e
    )
    if (!inherits(tryResult, "error")) {
      return(tryResult)
    }
    tryResult <- tryCatch(eval(matchedCall$term), error = function(e) e)
    if (!inherits(tryResult, "error")) return(tryResult)
  }

  NULL
}

## this used to be a function evaluated in the caller's frame, but
## that causes warnings in R check so now it is just a block of code
getTestOffset <- quote({
  if (is.numeric(matchedCall$offset.test)) {
    return(list(offset.test = offset.test, testUsesRegularOffset = FALSE))
  }
  if (is.null(matchedCall$offset.test)) {
    return(list(offset.test = NULL, testUsesRegularOffset = FALSE))
  }

  if (is.symbol(matchedCall$offset.test)) {
    testOffsetName <- as.character(matchedCall$offset.test)

    if (identical(testOffsetName, "offset") && !is.null(offset)) {
      return(list(
        offset.test = if (offsetGivenAsScalar == TRUE) offset[1] else offset,
        testUsesRegularOffset = TRUE
      ))
    }

    if (is.formula(formula)) {
      if (!dataIsMissing && any(names(data) == testOffsetName)) {
        return(list(
          offset.test = data[[testOffsetName]],
          testUsesRegularOffset = FALSE
        ))
      }
      if (any(ls(environment(formula)) == testOffsetName)) {
        return(list(
          offset.test = get(testOffsetName, environment(formula)),
          testUsesRegularOffset = FALSE
        ))
      }
    }
    tryResult <- tryCatch(get(testOffsetName), error = function(e) e)
    if (!inherits(tryResult, "error") && !is.null(tryResult)) {
      return(list(offset.test = tryResult, testUsesRegularOffset = FALSE))
    }

    stop("cannot find test offset '", testOffsetName, "'")
  } else if (is.language(matchedCall$offset.test)) {
    ## offset.test could have been something like (offset + 0.5), or (offset + variable)
    baseOffset <- if (is.null(offset)) {
      NA_real_
    } else {
      if (offsetGivenAsScalar == TRUE) offset[1] else offset
    }

    if (identical(matchedCall$offset.test, quote(offset))) {
      return(list(offset.test = baseOffset, testUsesRegularOffset = TRUE))
    }

    testOffset <- subTermInLanguage(
      matchedCall$offset.test,
      quote(offset),
      baseOffset
    )

    if (is.formula(formula)) {
      if (!dataIsMissing) {
        tryResult <- with(
          data,
          tryCatch(eval(testOffset), error = function(e) e)
        )
        if (!inherits(tryResult, "error")) {
          return(list(offset.test = tryResult, testUsesRegularOffset = FALSE))
        }
      }
      tryResult <- tryCatch(
        eval(testOffset, environment(formula)),
        error = function(e) e
      )
      if (!inherits(tryResult, "error")) {
        return(list(offset.test = tryResult, testUsesRegularOffset = FALSE))
      }
    }
    tryResult <- tryCatch(
      eval(testOffset, parent.frame(3L)),
      error = function(e) e
    )
    if (!inherits(tryResult, "error")) {
      return(list(offset.test = tryResult, testUsesRegularOffset = FALSE))
    }
  }

  stop("cannot construct test offset")
})

# Classify a raw response for family routing. A factor/ordered/logical/
# character response declares a classification model; numeric passes through.
# The level count drives the 2 -> probit vs 3+ -> multinomial (bart2) / refusal
# (single-forest) split downstream. Shared by dbartsData's response coding and
# bart2's family = "auto" peek.
classifyResponse <- function(y) {
  if (is.factor(y)) {
    list(
      type = if (is.ordered(y)) "ordered factor" else "factor",
      n.levels = nlevels(y)
    )
  } else if (is.logical(y) && !is.matrix(y)) {
    list(type = "logical", n.levels = 2L)
  } else if (is.character(y) && !is.matrix(y)) {
    list(type = "character", n.levels = length(unique(y[!is.na(y)])))
  } else {
    list(type = "numeric", n.levels = NA_integer_)
  }
}

# Code a raw response to the doubles the engine reads and report its original
# type. A factor (or a character coerced with factor(), or a logical) becomes
# 0-based codes exactly as the historic x/y path did; a numeric response is
# passed through as.double, byte-identical to the previous handling. The
# original levels are returned alongside so an ordinal fit can round-trip them
# (docs/design/ordinal.md section 5); they are NULL for a numeric response
# (a numeric ordinal derives sort(unique(y)) itself) and were previously
# discarded.
codeResponse <- function(y) {
  info <- classifyResponse(y)
  coded <- if (info$type == "numeric" || info$type == "logical") {
    as.double(y)
  } else {
    as.double(as.integer(if (is.character(y)) factor(y) else y) - 1L)
  }
  levels <- if (is.factor(y)) {
    levels(y)
  } else if (info$type == "character") {
    levels(factor(y))
  } else {
    NULL
  }
  list(y = coded, type = info$type, n.levels = info$n.levels, levels = levels)
}

# Resolve family for a categorical (factor/logical/character) response ahead
# of a single-forest fit: a 2-level response is a binary classification
# (family = "auto" fits probit); 3+ levels is multinomial, which none of
# dbarts()/xbart()/rbart_vi() implement (only bart2(family = "multinomial")
# does). A numeric response is returned unchanged - a caller that also
# resolves the 0/1-vs-continuous ambiguity for numeric responses does so
# itself afterward (dbarts()/xbart() do; rbart_vi() defers it to the
# per-chain dbarts() call, since there is no message to deduplicate there).
#
# `caller` names the entry point for the 2-level conflict message and the
# non-split multinomial message ("CALLER does not fit a K-level ..."), used
# as-is by xbart() and rbart_vi(): each is reached only directly, so the
# auto/explicit distinction adds nothing at K >= 3 (every family choice is
# equally invalid). dbarts() passes splitMultinomialMessage = TRUE instead,
# because it is also reached anonymously through bart() (which never sets an
# explicit family): its auto-branch message cannot name a single caller and
# instead lists every single-forest entry point (passed via `caller`), while
# its explicit-family branch echoes the conflicting family like the 2-level
# message does.
resolveClassificationFamily <- function(
  data,
  family,
  caller,
  incompatibleFamilies,
  splitMultinomialMessage = FALSE,
  allowOrdinal = FALSE
) {
  responseType <- data@response.type
  K <- data@response.n.levels
  # ordinal (cumulative probit, docs/design/ordinal.md section 5): only
  # dbarts()/bart2() pass allowOrdinal = TRUE. family = "ordinal" is the
  # explicit primitive and forces the model on any response (numeric levels are
  # sort(unique(y)), resolved later); family = "auto" auto-dispatches an ORDERED
  # factor to it, announced. The other single-forest entries leave allowOrdinal
  # FALSE and fall through to their K >= 3 refusals below. is.ordered() is the
  # disjoint key: an unordered K >= 3 factor stays multinomial.
  if (allowOrdinal) {
    if (identical(family, "ordinal")) {
      return("ordinal")
    }
    # a 2-level ordered factor is binary (probit); only a 3+-level ordered
    # factor is a genuine ordinal scale worth auto-dispatching
    if (family == "auto" && responseType == "ordered factor" && K >= 3L) {
      announceAutoFamily(responseType, K, "ordinal")
      return("ordinal")
    }
  }
  if (responseType == "numeric") {
    return(family)
  }
  if (K >= 3L) {
    # a 3+-level ORDERED factor is ordinal (reached here only from an entry that
    # cannot fit it - xbart, rbart_vi; dbarts/bart2 route ordinal above); every
    # other 3+-level factor/character is unordered multinomial
    isOrdered <- identical(responseType, "ordered factor")
    model <- if (isOrdered) "ordinal" else "multinomial"
    suggestion <- if (isOrdered) {
      "bart2(family = \"ordinal\")"
    } else {
      "bart2(family = \"multinomial\")"
    }
    if (!splitMultinomialMessage) {
      stop(
        caller,
        " does not fit a ",
        K,
        "-level ",
        responseType,
        " response; ",
        model,
        " classification requires ",
        suggestion
      )
    }
    if (family == "auto") {
      stop(
        "a ",
        K,
        "-level ",
        responseType,
        " response is ",
        model,
        "; fit it with ",
        suggestion,
        " - ",
        caller,
        " fit only binary and continuous responses"
      )
    }
    stop(
      "family \"",
      family,
      "\" cannot fit a ",
      K,
      "-level ",
      responseType,
      " response; a 3+-level ",
      responseType,
      " is ",
      model,
      " (",
      suggestion,
      ")"
    )
  }
  if (family == "auto") {
    family <- "probit"
    announceAutoFamily(responseType, K, family)
  } else if (family %in% incompatibleFamilies) {
    stop(
      "family \"",
      family,
      "\" cannot fit a ",
      responseType,
      " response; a 2-level factor is a binary classification ",
      "(family = \"auto\" fits probit)"
    )
  }
  family
}

# Resolve the ordered category structure for an ordinal (cumulative-probit) fit
# from a dbartsData whose family has resolved to "ordinal" (docs/design/
# ordinal.md section 5). Returns the ONE-based category codes the engine reads
# (its y_ holds 1..K), the count K, and the ordered level labels for the
# round-trip. A factor/character response takes its stored level order - an
# UNORDERED one with an informational note, since the factor default
# (alphabetical) is rarely the intended scale; a numeric/integer/logical
# response uses sort(unique(y)) as the ordered levels. Called once, where the
# recoding happens (dbarts()), so the note fires once.
resolveOrdinalResponse <- function(data) {
  labels <- data@response.levels
  if (is.null(labels)) {
    # numeric/integer/logical: the distinct sorted values are the ordered levels
    levels <- sort(unique(data@y))
    codes <- match(data@y, levels)
    labels <- as.character(levels)
  } else {
    # data@y holds 0-based factor codes (codeResponse); the engine wants 1..K
    codes <- as.integer(data@y) + 1L
    if (data@response.type != "ordered factor") {
      message(
        "family = \"ordinal\": the ",
        data@response.type,
        " response is unordered; its category order is taken from the level ",
        "order (",
        paste(labels, collapse = " < "),
        ")"
      )
    }
  }
  K <- length(labels)
  if (K < 2L) {
    stop("family = \"ordinal\" requires a response with at least 2 categories")
  }
  list(y = as.double(codes), K = as.integer(K), levels = labels)
}

# Resolve the negative-binomial dispersion argument to the length-1 real the C
# bridge reads off the control's bartcore.dispersion attribute (docs/design/
# negative-binomial.md section 2). NA (the default) estimates r on the capped
# integer grid, encoded as a non-positive spec; a supplied value FIXES r. v1
# ships the exact integer envelope, so a fixed dispersion must be a single
# positive integer - a real fixed value is refused informatively, so admitting
# it later is a validation relaxation, not a signature change.
resolveDispersion <- function(dispersion) {
  if (length(dispersion) != 1L) {
    stop("'dispersion' must be a single value")
  }
  if (is.na(dispersion)) {
    return(-1) # a non-positive spec estimates r on the grid
  }
  if (!is.numeric(dispersion) || dispersion <= 0) {
    stop("'dispersion' must be a positive number")
  }
  if (dispersion != round(dispersion)) {
    stop(
      "family \"nbinom\" fits an integer dispersion; real dispersion is not ",
      "yet supported"
    )
  }
  as.double(dispersion)
}

dbartsData <- function(
  formula,
  data,
  test,
  subset,
  weights,
  offset,
  offset.test = offset,
  factors = c("categorical", "indicators"),
  missing = c("incorporate", "error")
) {
  dataIsMissing <- missing(data)
  testIsMissing <- missing(test)
  offsetIsMissing <- missing(offset)
  testOffsetIsMissing <- missing(offset.test)
  matchedCall <- match.call()

  # "indicators" dummy-expands factor columns as always; "categorical" keeps
  # them as single columns split by category subset, which only the bartcore
  # engine runs
  factors <- match.arg(factors)
  makeModelMatrix <- if (factors == "categorical") {
    makeCategoricalModelMatrix
  } else {
    makeModelMatrixFromDataFrame
  }
  # "incorporate" keeps NAs in the predictors - the trees route them by
  # learned per-rule directions; "error" rejects them
  missing <- match.arg(missing)

  offsetGivenAsScalar <- NA
  testUsesRegularOffset <- NA
  # the response's original type, recorded on the result so the fitters can
  # route family = "auto" and reject a categorical response an unsupported
  # family cannot fit; each y-producing branch below refreshes it
  responseInfo <- list(type = "numeric", n.levels = NA_integer_, levels = NULL)

  if (missing(formula)) {
    stop("first argument to dbartsData - 'formula'/'x.train' - must be present")
  }

  if (inherits(formula, "dbartsData")) {
    if (
      !dataIsMissing ||
        !testIsMissing ||
        !offsetIsMissing ||
        !testOffsetIsMissing
    ) {
      warning("if data supplied as dbartsData, remaining arguments are ignored")
    }
    return(formula)
  }

  if (is.formula(formula)) {
    if (
      !dataIsMissing &&
        !is.data.frame(data) &&
        !is.list(data) &&
        !is.environment(data)
    ) {
      stop(
        "for formula/data specification, data must be a data frame, list, or environment"
      )
    }

    modelFrameArgs <- c("formula", "data", "subset", "weights", "offset")

    ## extract offset prematurely, if necessary
    if (offsetIsMissing) {
      offset <- NULL
      modelFrameArgs <- c("formula", "data", "subset", "weights")
    } else {
      offsetCall <- matchedCall
      offsetCall <- offsetCall[c(
        1L,
        match(c("formula", "data", "offset"), names(offsetCall), nomatch = 0L)
      )]
      names(offsetCall)[which(names(offsetCall) == "offset")] <- "term"
      offsetCall[[1L]] <- quoteInNamespace(findTermInFormulaData)
      offset <- eval(offsetCall, parent.frame())

      if (!is.null(offset)) {
        offsetGivenAsScalar <- length(offset) == 1
        if (offsetGivenAsScalar) {
          modelFrameArgs <- c("formula", "data", "subset", "weights")
        }
      }
    }

    # pre-validate lengths against a known y/data length so a mismatch reads
    # as our own message rather than model.frame's "variable lengths differ
    # (found for '(weights)')"; a data.frame is the only case where the
    # eventual length is known this early without duplicating model.frame's
    # own work
    dataLength <- if (!dataIsMissing && is.data.frame(data)) {
      nrow(data)
    } else {
      NA_integer_
    }
    if (!is.na(dataLength)) {
      if (
        !is.null(offset) &&
          !isTRUE(offsetGivenAsScalar) &&
          length(offset) != dataLength
      ) {
        stop("length of 'offset' must equal length of 'y'")
      }
      if (!missing(weights)) {
        weightsCall <- matchedCall[c(
          1L,
          match(
            c("formula", "data", "weights"),
            names(matchedCall),
            nomatch = 0L
          )
        )]
        names(weightsCall)[names(weightsCall) == "weights"] <- "term"
        weightsCall[[1L]] <- quoteInNamespace(findTermInFormulaData)
        weightsValue <- tryCatch(
          eval(weightsCall, parent.frame()),
          error = function(e) NULL
        )
        if (!is.null(weightsValue) && length(weightsValue) != dataLength) {
          stop("length of 'weights' must equal length of 'y'")
        }
      }
    }

    modelFrameCall <- matchedCall
    modelFrameCall <- modelFrameCall[c(
      1L,
      match(modelFrameArgs, names(modelFrameCall), nomatch = 0L)
    )]
    modelFrameCall$drop.unused.levels <- FALSE
    # incomplete predictor rows stay; completeness is validated below
    # (previous versions silently na.omit-dropped them)
    modelFrameCall$na.action <- stats::na.pass
    modelFrameCall[[1L]] <- quote(stats::model.frame)
    ## this allows subset to be applied to offset, even if offset was a language construct (e.g. off + 0.1)
    if (identical(offsetGivenAsScalar, FALSE)) {
      modelFrameCall$offset <- offset
    }

    # a sparseFactor column would die inside model.frame with a bare S4
    # type error; refuse it explicitly first (the formula path takes plain
    # data-frame columns only, not S4 predictors)
    if (!dataIsMissing && (is.list(data) || is.environment(data))) {
      for (variableName in intersect(all.vars(formula), names(data))) {
        if (methods::is(data[[variableName]], "sparseFactor")) {
          stop(
            "sparse categorical predictors must be supplied through the ",
            "x/y interface; '",
            variableName,
            "' is a sparseFactor"
          )
        }
      }
    }

    modelFrame <- eval(modelFrameCall, parent.frame())
    if (NROW(modelFrame) == 0) {
      if (!is.null(matchedCall$subset)) {
        stop("empty 'subset' specified")
      }
      stop("cannot construct model matrices from formula")
    }

    ## pull out y - NO type coercion, so a factor response keeps its levels
    ## (model.response(., "numeric") would leave it a factor and warn, then
    ## trip "range not meaningful for factors" downstream); codeResponse then
    ## routes it exactly as the x/y path does
    y <- model.response(modelFrame)
    # a Surv response reaches here only through the formula interface, which
    # survival (aft) fits do not support; refuse with the supported surface
    # before any arithmetic trips survival's Ops.Surv guard
    if (inherits(y, "Surv")) {
      stop(
        "survival (Surv) responses are not supported by the formula ",
        "interface; use the matrix interface - dbarts(x.train, y.train) or ",
        "bart2(x.train, y.train) with a Surv or two-column (time, status) ",
        "response and family = \"aft\""
      )
    }
    if (is.null(y)) {
      y <- rep(0, NROW(modelFrame))
    }
    coded <- codeResponse(y)
    y <- coded$y
    responseInfo <- coded[c("type", "n.levels", "levels")]
    numObservations <- NROW(y)

    ## weights
    weights <- as.vector(model.weights(modelFrame))
    if (!is.null(weights)) {
      if (!is.numeric(weights)) {
        stop("'weights' must be of type numeric")
      }
      weights <- as.double(weights)
    }

    ## offset, when in data frame
    if (identical(offsetGivenAsScalar, FALSE)) {
      offset <- as.vector(model.offset(modelFrame))
    } else if (identical(offsetGivenAsScalar, TRUE)) {
      offset <- rep_len(offset, numObservations)
    }

    ## predictors
    modelTerms <- terms(modelFrame)
    if (is.empty.model(modelTerms)) {
      stop("predictors must be specified for regression tree analysis")
    }

    termLabels <- attr(modelTerms, "term.labels")
    badLabels <- grepl("`.* .*`", termLabels)
    if (sum(badLabels) > 0) {
      termLabels[badLabels] <- gsub("^`(.*)`$", "\\1", termLabels[badLabels])
    }

    x <- makeModelMatrix(modelFrame[termLabels])

    if (!testIsMissing) {
      testCall <- matchedCall
      testCall <- testCall[c(
        1L,
        match(c("formula", "data", "test"), names(testCall), nomatch = 0L)
      )]
      names(testCall)[which(names(testCall) == "test")] <- "term"
      testCall[[1L]] <- quoteInNamespace(findTermInFormulaData)

      temp <- eval(testCall, parent.frame())
      if (!is.null(temp)) test <- temp
    }
  } else if (inherits(formula, "dgCMatrix")) {
    ## sparse designs enter through the x/y interface only; columns are all
    ## ordinal and missing values are stored NaN entries (the Matrix
    ## convention), so no complete-case filtering applies
    if (dataIsMissing || is.null(data)) {
      data <- rep(0, nrow(formula))
    }
    if (
      !is.numeric(data) &&
        !is.factor(data) &&
        !is.logical(data) &&
        !is.character(data)
    ) {
      stop(
        "when 'formula' is a sparse matrix, 'data' must be numeric, a ",
        "factor, logical, or character"
      )
    }

    coded <- codeResponse(data)
    y <- coded$y
    responseInfo <- coded[c("type", "n.levels", "levels")]
    if (nrow(formula) != NROW(y)) {
      stop("'x' must have the same number of observations as 'y'")
    }
    initialNumObservations <- NROW(y)
    # an empty training set has no cut grid to quantize and would fault deeper
    # (a subscript-out-of-bounds subsetting the zero rows below); name it here,
    # as the formula path already rejects an empty model frame
    if (initialNumObservations == 0L) {
      stop("data has zero rows")
    }

    if (missing(subset) || is.null(subset)) {
      subset <- seq.int(length(y))
    }
    y <- y[subset]
    x <- formula[subset, , drop = FALSE]

    if (missing(weights)) {
      weights <- NULL
    }
    if (!is.null(weights)) {
      if (!is.numeric(weights)) {
        stop("'weights' must be a numeric vector")
      }
      weights <- as.double(weights)
      if (length(weights) == 1L) {
        weights <- rep_len(weights, initialNumObservations)
      }
      if (length(weights) != initialNumObservations) {
        stop("length of 'weights' must equal length of 'y'")
      }
      weights <- weights[subset]
    }

    if (offsetIsMissing) {
      offset <- NULL
    }
    if (!is.null(offset)) {
      if (!is.numeric(offset)) {
        stop("'offset' must be numeric")
      }
      if (length(offset) == 1L) {
        offset <- rep_len(offset, initialNumObservations)
        offsetGivenAsScalar <- TRUE
      } else {
        offsetGivenAsScalar <- FALSE
      }
      if (length(offset) != initialNumObservations) {
        stop("length of 'offset' must equal length of 'y'")
      }
      offset <- offset[subset]
    }
  } else if (
    is.numeric(formula) || is.data.frame(formula) || is.factor(formula)
  ) {
    ## backwards compatibility of bart(x.train, y.train, x.test)
    if (dataIsMissing || is.null(data)) {
      data <- rep(0, NROW(formula))
    }
    if (
      !is.numeric(data) &&
        !is.data.frame(data) &&
        !is.factor(data) &&
        !is.logical(data) &&
        !is.character(data)
    ) {
      stop(
        "when 'formula' is numeric, 'data' must be numeric, a factor, ",
        "logical, or character"
      )
    }

    coded <- codeResponse(data)
    y <- coded$y
    responseInfo <- coded[c("type", "n.levels", "levels")]
    if (NROW(formula) != NROW(y)) {
      stop("'x' must have the same number of observations as 'y'")
    }
    initialNumObservations <- NROW(y)
    # an empty training set has no cut grid to quantize and would fault deeper
    # (a subscript-out-of-bounds subsetting the zero rows below); name it here,
    # as the formula path already rejects an empty model frame
    if (initialNumObservations == 0L) {
      stop("data has zero rows")
    }

    if (missing(subset) || is.null(subset)) {
      subset <- seq.int(length(y))
    }
    y <- y[subset]

    if (is.data.frame(formula)) {
      formula <- makeModelMatrix(formula)
    }
    xIsMixed <- inherits(formula, "dbartsMixedMatrix")
    x <- if (is.matrix(formula) || xIsMixed) {
      formula[subset, , drop = FALSE]
    } else {
      formula[subset]
    }

    if (missing(weights)) {
      weights <- NULL
    }
    if (!is.null(weights)) {
      if (!is.numeric(weights)) {
        stop("'weights' must be a numeric vector")
      }
      weights <- as.double(weights)
      if (length(weights) == 1L) {
        weights <- rep_len(weights, initialNumObservations)
      }
      if (length(weights) != initialNumObservations) {
        stop("length of 'weights' must equal length of 'y'")
      }
      weights <- weights[subset]
    }

    if (offsetIsMissing) {
      offset <- NULL
    }
    if (!is.null(offset)) {
      if (!is.numeric(offset)) {
        stop("'offset' must be numeric")
      }
      if (length(offset) == 1L) {
        offset <- rep_len(offset, initialNumObservations)
        offsetGivenAsScalar <- TRUE
      } else {
        offsetGivenAsScalar <- FALSE
      }
      if (length(offset) != initialNumObservations) {
        stop("length of 'offset' must equal length of 'y'")
      }
      offset <- offset[subset]
    }

    # a mixed container keeps its rows and attributes: missing predictor
    # values are validated below, like the sparse-matrix branch above
    if (!xIsMixed) {
      # missing = "incorporate" must reach the shared NA handling at the end
      # of this function intact, exactly like the sparse and mixed-container
      # branches above - no row is dropped here for missingness in 'x' or
      # 'y'. That shared code unconditionally rejects a missing response
      # (anyNA(y) below) and only rejects a missing predictor when
      # missing = "error", so deferring to it makes this branch agree with
      # the formula interface instead of silently discarding incomplete
      # rows regardless of 'missing' (the bug: incorporation was only
      # reachable through the formula path, since this branch always
      # complete-cases-filtered first).
      completeCases <- if (missing == "error") {
        stats::complete.cases(x, y)
      } else {
        rep_len(TRUE, length(y))
      }

      # the check below happens before xHasNA is evaluated further down, so
      # missing = "error" must be checked against the pre-filter state here
      # or its request is silently defeated
      if (missing == "error") {
        if (anyNA(y)) {
          stop("response contains missing values")
        }
        if (anyNA(x)) {
          stop(
            "predictors contain missing values; use missing = \"incorporate\" to model them"
          )
        }
      }

      y <- y[completeCases]
      x <- if (!is.matrix(x)) {
        x[completeCases]
      } else {
        x[completeCases, , drop = FALSE]
      }
      if (length(attributes(formula)) > 0L) {
        for (attributeName in names(attributes(formula))) {
          if (attributeName == "dim") {
            next
          }
          if (attributeName == "dimnames" && !identical(dim(formula), dim(x))) {
            next
          }
          attr(x, attributeName) <- attr(formula, attributeName)
        }
      }
      if (!is.null(weights)) {
        weights <- weights[completeCases]
      }
      if (!is.null(offset)) offset <- offset[completeCases]
    }
  } else {
    stop(
      "unrecognized 'formula' type; must be coercible to numeric or a valid formula object"
    )
  }

  if (is.vector(x)) {
    x <- as.matrix(x)
  }
  if (is.data.frame(x)) {
    x <- makeModelMatrix(x)
  }

  x.test <- NULL
  if (!testIsMissing && !is.null(test)) {
    x.test <- validateXTest(test, x)
  }

  if (!is.null(x.test)) {
    if (testOffsetIsMissing) {
      ## default is offset.test = offset
      if (identical(offsetGivenAsScalar, TRUE)) {
        offset.test <- rep_len(offset[1L], nrow(x.test))
        testUsesRegularOffset <- TRUE
      } else if (identical(offsetGivenAsScalar, FALSE)) {
        if (nrow(x.test) != length(y)) {
          stop(
            "vectored 'offset' cannot be directly applied to test data of unequal length"
          )
        }
        offset.test <- offset
        testUsesRegularOffset <- TRUE
      }
    } else {
      testOffsetInfo <- eval(getTestOffset)

      offset.test <- testOffsetInfo$offset.test
      testUsesRegularOffset <- testOffsetInfo$testUsesRegularOffset

      if (!is.null(offset.test)) {
        offset.test <- rep_len(offset.test, nrow(x.test))
      }
    }
  } else {
    if (testOffsetIsMissing) offset.test <- NULL
  }

  weights.test <- NULL
  if (!is.null(x.test) && !is.null(matchedCall$weights)) {
    if (!is.formula(formula)) {
      warning(
        "'weights' are ignored for test data when model is not specified as a formula; this only impacts extracting samples from the posterior predictive distribution of the test data"
      )
    } else {
      testFormula <- formula
      lhs <- testFormula[[2L]]
      remainder <- testFormula
      remainder[[2L]] <- NULL
      testFormula <- as.formula(paste0(deparse(remainder), " - ", deparse(lhs)))
      environment(testFormula) <- environment(formula)
      modelFrameCall$formula <- testFormula
      modelFrameCall$data <- test
      tryResult <- tryCatch(
        testFrame <- eval(modelFrameCall, parent.frame()),
        error = function(e) e
      )
      if (inherits(tryResult, "error")) {
        warning("weights specified but not found in test data - ignoring")
      } else {
        weights.test <- testFrame[["(weights)"]]
      }
    }
  }

  # missingness is a predictor-only feature: rules route NAs in x, but the
  # response side must be complete. In a sparse x, NAs live only among the
  # stored entries and implicit zeros are observed values, so the checks
  # work off the slots without densifying.
  if (anyNA(y)) {
    stop("response contains missing values")
  }
  # Inf/-Inf survive the anyNA check (NaN does not), then poison the
  # precision-degeneracy ratio below into an NA condition; reject them here
  # with a named error instead.
  if (any(is.infinite(y))) {
    stop("response contains non-finite values")
  }

  # Precision-degenerate response: a response with a
  # large magnitude but a tiny spread quantizes to (near-)identical double
  # values before the engine ever sees it - e.g. y in [1e15, 1e15 + 1e-3]
  # rounds to a single representable double - so the engine fits an
  # apparently-fine model that can no longer tell observations apart. The
  # signature is range(y) tiny relative to scale(y): doubles near
  # magnitude s are spaced ~2.22e-16 * s apart, so 1e-10 is ~1e6x the ulp
  # spacing - real precision loss trips it with orders of magnitude of
  # headroom before a merely low-variance (but not degenerate) response
  # would. max(abs(y)) == 0 is guarded to avoid a 0/0; an all-zero
  # response has no precision loss to report (it is exactly, not
  # approximately, constant). No distinct-value-count check backs this up:
  # a rounding collapse to a handful of distinct doubles means the range
  # spans at most a few ulps (~1e-15 relative), which the ratio already
  # catches five orders of magnitude earlier, while low-cardinality data
  # the ratio does NOT flag has values separated by many ulps -
  # legitimately discrete responses (binary, counts, ordinal scales) that
  # standardization maps cleanly and the engine fits without degradation.
  yRange <- diff(range(y))
  yScale <- max(abs(y))
  if (yScale > 0 && yRange / yScale < 1e-10) {
    warning(
      "response values are indistinguishable, or nearly so, at double ",
      "precision (",
      length(unique(y)),
      " distinct value(s) among ",
      length(y),
      " observations); center and/or rescale the response before fitting"
    )
  }

  sparseAllMissingCheck <- function(x.sparse) {
    columnNnz <- diff(x.sparse@p)
    columnNumNA <- vapply(
      seq_len(ncol(x.sparse)),
      function(j) {
        sum(is.na(x.sparse@x[seq.int(
          x.sparse@p[j] + 1L,
          length.out = columnNnz[j]
        )]))
      },
      0L
    )
    if (any(columnNnz == nrow(x.sparse) & columnNumNA == nrow(x.sparse))) {
      stop("predictor columns cannot be entirely missing")
    }
  }
  if (is.matrix(x)) {
    xHasNA <- anyNA(x)
    if (xHasNA && any(colSums(!is.na(x)) == 0L)) {
      stop("predictor columns cannot be entirely missing")
    }
  } else if (inherits(x, "dbartsMixedMatrix")) {
    # both flavors hold dense columns as a per-column list; the mixed flavor
    # adds a sparse part
    denseColumns <- if (is.null(x$dense)) list() else x$dense
    denseHasNA <- any(vapply(denseColumns, anyNA, FALSE))
    sparseHasNA <- !is.null(x$sparse) && anyNA(x$sparse@x)
    xHasNA <- denseHasNA || sparseHasNA
    if (xHasNA) {
      denseAllMissing <- vapply(
        denseColumns,
        function(column) all(is.na(column)),
        FALSE
      )
      if (any(denseAllMissing)) {
        stop("predictor columns cannot be entirely missing")
      }
      if (!is.null(x$sparse)) {
        sparseAllMissingCheck(x$sparse)
      }
    }
  } else {
    xHasNA <- anyNA(x@x)
    if (xHasNA) sparseAllMissingCheck(x)
  }
  if (!is.null(offset) && anyNA(offset)) {
    stop("'offset' contains missing values")
  }
  if (!is.null(offset.test) && anyNA(offset.test)) {
    stop("'offset.test' contains missing values")
  }
  if (missing == "error") {
    if (xHasNA) {
      stop(
        "predictors contain missing values; use missing = \"incorporate\" to model them"
      )
    }
    if (!is.null(x.test) && anyNA(x.test)) {
      stop(
        "test predictors contain missing values; use missing = \"incorporate\" to model them"
      )
    }
  }

  result <- newValidated(
    "dbartsData",
    modelMatrices = namedList(
      y,
      x,
      x.test,
      weights,
      weights.test,
      offset,
      offset.test,
      testUsesRegularOffset
    ),
    n.cuts = NA_integer_,
    sigma = NA_real_
  )
  result@missing <- missing
  result@response.type <- responseInfo$type
  result@response.n.levels <- as.integer(responseInfo$n.levels)
  result@response.levels <- responseInfo$levels
  result
}
