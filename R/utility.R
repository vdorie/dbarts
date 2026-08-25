## validObject()'s default failure message is prefixed with S4 boilerplate
## ("invalid class \"...\" object: ") that leaks the internal class name;
## strip it so the sampler's user-facing validation errors read as plain
## sentences. Any other error is re-signaled unchanged.
rethrowValidityError <- function(e) {
  msg <- conditionMessage(e)
  cleaned <- sub("^invalid class [^ ]+ object: ", "", msg)
  if (identical(cleaned, msg)) {
    stop(e)
  }
  stop(cleaned, call. = FALSE)
}

newValidated <- function(Class, ...) {
  tryCatch(new(Class, ...), error = rethrowValidityError)
}

# One-line verdict when family = "auto" resolves a categorical response to a
# non-default family: probit for a 2-level response, multinomial for a 3+-level
# UNORDERED factor/character, ordinal for a 3+-level ORDERED factor (the level
# order is the category order - respected, not discarded). Set family
# explicitly to override.
announceAutoFamily <- function(responseType, nLevels, family) {
  message(
    "family = \"auto\": ",
    nLevels,
    "-level ",
    responseType,
    " response detected, fitting family = \"",
    family,
    "\"; set 'family' to override"
  )
}

# Family gating: an argument supplied by name whose only effect is on a
# response family this fit did not resolve to is diagnosed rather than
# silently dropped. Each row names the arguments it covers, the families
# under which they DO act ("liveIn"; every other family gates them), and why
# not; resolved once as data so each gated argument - resid.prior included -
# is a row here rather than a branch in warnFamilyGatedArgs.
familyGatingInventory <- list(
  list(
    names = c("sigest", "sigdf", "sigquant"),
    liveIn = c("gaussian", "aft", "hurdle.lognormal"),
    reason = "the residual scale is fixed, not estimated"
  ),
  list(
    # A fixed-unit-scale family overwrites a supplied resid.prior with
    # fixed(1) regardless (R/spec.R's fixedUnitScale branch), which was
    # silent before bart2 could forward a resid.prior object at all
    names = "resid.prior",
    liveIn = c("gaussian", "aft", "hurdle.lognormal"),
    reason = "the residual scale is fixed, not estimated"
  ),
  list(
    names = "dispersion",
    liveIn = "nbinom",
    reason = "only family = \"nbinom\" estimates a count dispersion"
  ),
  list(
    names = c("breaks", "max.rows"),
    liveIn = c("hazard", "hazard.probit", "hazard.logistic"),
    reason = "only a discrete-time hazard fit expands a time grid"
  )
)

# Warns once, naming every argument this call's resolved 'family' cannot act
# on, why, and the family itself (one warning per call, a classed
# condition). suppliedNames is the caller's own argNames snapshot - the
# matchedCall names taken before any family branch, so an unsupplied
# (defaulted) argument never appears here. A no-op when nothing is gated.
warnFamilyGatedArgs <- function(suppliedNames, family) {
  hits <- Filter(
    function(row) {
      family %not_in% row$liveIn && any(row$names %in% suppliedNames)
    },
    familyGatingInventory
  )
  if (length(hits) == 0L) {
    return(invisible(NULL))
  }
  clauses <- vapply(
    hits,
    function(row) {
      gatedNames <- row$names[row$names %in% suppliedNames]
      paste0(
        paste0("'", gatedNames, "'", collapse = ", "),
        " (",
        row$reason,
        ")"
      )
    },
    character(1L)
  )
  warning(warningCondition(
    paste0(
      "family = \"",
      family,
      "\" has no use for ",
      paste0(clauses, collapse = "; "),
      "; ignored"
    ),
    class = c("dbartsFamilyGatedWarning", "dbartsWarning")
  ))
  invisible(NULL)
}

# dbarts()/dbartsControl() accept n.samples %/% n.thin == 0 - a sampler
# meant to be driven by a host loop's own run() calls, never this entry
# point's. bart2/xbart/rbart_vi all return posterior draws, so the same
# zero would fault deeper (an empty-array reshape); refused here instead,
# by one message naming the entry point that was called.
refuseZeroSamples <- function(caller) {
  stop(
    "'n.samples' must leave at least one draw after thinning (n.samples ",
    "%/% n.thin = 0); dbarts() and dbartsControl() accept a zero-draw run ",
    "- a sampler driven by a host loop - but ",
    caller,
    "() returns posterior draws"
  )
}

# Missingness predicate for a setPredictor/setTestPredictor/
# setTestPredictorAndOffset argument. Never plain anyNA() on a
# dbartsMixedMatrix: it is a plain list, and base anyNA on a list is TRUE
# iff some TOP-LEVEL element is a length-one atomic NA - true of a
# single-sparse-column container's sparseReference sentinel (a false
# missing value) and false whenever a real NA sits in sparse@x behind a
# longer metadata vector (a missed one). recursive = TRUE fixes neither: it
# also descends into sparseReference. A sparseMatrix reads its own @x
# rather than densifying through the default method; a pattern matrix (no
# @x) carries no values to be missing.
sourceAnyNA <- function(x) {
  if (inherits(x, "dbartsMixedMatrix")) {
    return(
      any(vapply(x$dense, anyNA, NA)) ||
        (!is.null(x$sparse) && anyNA(x$sparse@x))
    )
  }
  if (methods::is(x, "sparseMatrix")) {
    return(methods::.hasSlot(x, "x") && anyNA(x@x))
  }
  anyNA(x)
}

# The dbartsSampler mutators (setPredictor/setTestPredictor/
# setTestPredictorAndOffset) refuse a new predictor matrix containing NA when
# the sampler was built with missing = "error"; `what` names the refused
# argument ("predictors"/"test predictors") in the message.
checkMissingPolicy <- function(data, hasMissing, what) {
  if (data@missing == "error" && hasMissing) {
    stop(
      "new ",
      what,
      " contain missing values and the sampler was built with missing = \"error\""
    )
  }
}

coerceOrError <- function(x, type) {
  mc <- match.call()

  if (is.null(x)) {
    stop("'", mc[[2L]], "' cannot be NULL")
  }

  func <- switch(
    type,
    logical = as.logical,
    integer = as.integer,
    numeric = as.numeric
  )
  result <- tryCatch(func(x), warning = function(e) e)
  if (inherits(result, "warning")) {
    stop("'", mc[[2L]], "' must be coercible to type: ", type)
  }

  result
}

"%not_in%" <- function(x, table) match(x, table, nomatch = 0L) <= 0L

## Bare name of the function a stored call invoked. A namespace-qualified call
## (dbarts::bart2) carries a `::` call in slot 1, which as.character() splits
## into c("::", "dbarts", "bart2") - comparing that to a name is a
## length-3 condition, an error since R 4.2.
callName <- function(call) {
  sub("^.*::", "", deparse(call[[1L]])[1L])
}

evalx.recurse <- function(x, e) {
  if (length(e) == 0L || typeof(e) == "symbol") {
    return(e)
  }

  for (i in seq_along(e)) {
    if (!is.language(e[[i]])) {
      next
    }

    e[[i]] <- if (e[[i]] == "x") x else evalx.recurse(x, e[[i]])
  }

  e
}

ifelse_3 <- function(cond1, cond2, then1, then2, else_) {
  mc <- match.call()
  env <- parent.frame()
  if (eval(mc[["cond1"]], env)) {
    then1
  } else if (eval(mc[["cond2"]], env)) {
    then2
  } else {
    else_
  }
}

## evaluates the expression 'e' after first replacing all instances of 'x' with
## the expression passed as x
evalx <- function(x, e) {
  mc <- match.call()
  callingEnv <- parent.frame()

  e <- evalx.recurse(mc$x, mc$e)
  eval(e, callingEnv)
}

redirectCall <- function(call, fn, ...) {
  matchedCall <- match.call()
  extraArgs <- if (length(matchedCall) > 3L) {
    as.character(matchedCall[-c(1L, 2L, 3L)])
  } else {
    character()
  }

  originalFn <- eval(call[[1L]])
  call[[1L]] <- if (is.function(fn)) matchedCall[[3L]] else fn
  if (length(extraArgs) == 0L) {
    fn <- if (is.function(fn)) fn else eval(fn)

    argsToKeep <- names(call)[-1L] %in% names(formals(fn))
    if (
      any(names(formals(originalFn)) == "...") &&
        any(names(formals(fn)) == "...")
    ) {
      argsToKeep <- argsToKeep |
        names(call)[-1L] %not_in% names(formals(originalFn))
    }

    call <- call[c(TRUE, argsToKeep)]
  } else {
    matchIndices <- match(extraArgs, names(call), nomatch = 0L)

    call <- call[c(1L, matchIndices)]
  }

  call
}

addCallArgument <- function(call, position, argument) {
  if (is.character(position)) {
    name <- position
    position <- length(call) + 1L
  } else {
    position <- as.integer(position) + 1L
    if (position <= length(call)) {
      for (i in seq.int(length(call), position)) {
        call[i + 1L] <- call[i]
        names(call)[[i + 1L]] <- names(call)[[i]]
      }
    }
    name <- ""
  }
  call[[position]] <- argument
  names(call)[[position]] <- name
  call
}

subTermInLanguage <- function(lang, oldTerm, newTerm) {
  if (length(lang) == 1L && is.symbol(lang)) {
    return(if (lang == oldTerm) newTerm else lang)
  }

  for (i in seq_along(lang)) {
    if (is.symbol(lang[[i]])) {
      if (lang[[i]] == oldTerm) lang[[i]] <- newTerm
    } else if (is.language(lang[[i]])) {
      lang[[i]] <- subTermInLanguage(lang[[i]], oldTerm, newTerm)
    }
  }
  lang
}

setDefaultsFromFormals <- function(call, formals, ...) {
  argsToReplace <- list(...)
  matchIndices <- match(argsToReplace, names(call), nomatch = 0L)
  missingFormals <- match(argsToReplace[matchIndices == 0L], names(formals))

  if (length(missingFormals) == 0L) {
    return(call)
  }

  newFormalIndices <- seq.int(length(missingFormals)) + length(call)
  call[newFormalIndices] <- formals[missingFormals]
  names(call)[newFormalIndices] <- names(formals)[missingFormals]
  call
}

is.formula <- function(x) is.language(x) && x[[1L]] == "~"

## from lme4
namedList <- function(...) {
  result <- list(...)
  substituteNames <- sapply(substitute(list(...)), deparse)[-1L]
  if (is.null(resultNames <- names(result))) {
    resultNames <- substituteNames
  }
  if (any(noNames <- resultNames == "")) {
    resultNames[noNames] <- substituteNames[noNames]
  }
  setNames(result, resultNames)
}

## Turns data.frame w/factors into matrices of indicator variables. Differs from
## model.matrix as it doesn't drop columns for co-linearity even with multiple
## factors
makeModelMatrixFromDataFrame <- function(x, drop = TRUE) {
  if (!is.data.frame(x)) {
    stop("x is not a dataframe")
  }
  if (is.logical(drop) && (length(drop) != 1L || is.na(drop))) {
    stop("when logical, drop must be TRUE or FALSE")
  }
  if (is.list(drop) && length(drop) != length(x)) {
    stop("when list, drop must have length equal to x")
  }

  characterCols <- sapply(x, typeof) == "character"
  if (any(characterCols)) {
    x[characterCols] <- lapply(x[characterCols], as.factor)
  }

  columnIsSparse <- vapply(x, isSparseDataFrameColumn, FALSE)
  if (!any(columnIsSparse)) {
    result <- .Call(C_dbarts_makeModelMatrixFromDataFrame, x, drop)
    attr(result, "term.labels") <- names(x)
    return(result)
  }

  # expand dense input columns one at a time - the C builder treats columns
  # independently, so the blocks and drop entries match a whole-frame call -
  # and splice sparse ones in place; their drop entries are all-FALSE so the
  # pattern replays over a fully dense test frame
  columns <- vector("list", length(x))
  blockNames <- vector("list", length(x))
  dropPattern <- if (is.list(drop) || isTRUE(drop)) {
    vector("list", length(x))
  } else {
    NULL
  }
  for (j in seq_along(x)) {
    if (columnIsSparse[j]) {
      # a sparse factor cannot be dummy-expanded without densifying it; the
      # categorical builder rides it to the engine unexpanded instead
      if (methods::is(x[[j]], "sparseFactor")) {
        stop(
          "sparse categorical predictors require factors = \"categorical\"; ",
          "column '",
          names(x)[j],
          "' is a sparseFactor"
        )
      }
      slices <- sparseColumnSlices(x[[j]], names(x)[j], nrow(x))
      columns[[j]] <- slices
      blockNames[[j]] <- slices$names
      if (!is.null(dropPattern)) {
        dropPattern[[j]] <- rep.int(FALSE, length(slices$i))
      }
    } else {
      block <- .Call(
        C_dbarts_makeModelMatrixFromDataFrame,
        x[j],
        if (is.list(drop)) drop[j] else drop
      )
      columns[[j]] <- block
      blockNames[[j]] <- colnames(block)
      if (!is.null(dropPattern)) {
        dropPattern[j] <- if (is.list(drop)) drop[j] else attr(block, "drop")
      }
    }
  }
  result <- assembleMixedMatrix(columns, columnIsSparse, blockNames, nrow(x))
  attr(result, "term.labels") <- names(x)
  if (!is.null(dropPattern)) {
    names(dropPattern) <- names(x)
    attr(result, "drop") <- dropPattern
  }
  result
}

## Turns a data.frame into a numeric matrix without expanding factors:
## unordered factors become categorical columns holding codes 0..K-1,
## ordered factors become ordinal codes, and matrix columns (e.g. from
## poly()) splice in as ordinal. Attaches "varTypes" (integer per column)
## and "factor.levels" (list per column, NULL for non-factors) for the
## sampler and for test-data mapping; only the bartcore engine accepts
## categorical columns.
makeCategoricalModelMatrix <- function(x) {
  if (!is.data.frame(x)) {
    stop("x is not a dataframe")
  }

  characterCols <- sapply(x, typeof) == "character"
  if (any(characterCols)) {
    x[characterCols] <- lapply(x[characterCols], as.factor)
  }

  columns <- vector("list", length(x))
  columnIsSparse <- logical(length(x))
  columnTypes <- vector("list", length(x))
  columnLevels <- vector("list", length(x))
  columnNames <- vector("list", length(x))
  for (j in seq_along(x)) {
    column <- x[[j]]
    name <- names(x)[j]
    if (isSparseDataFrameColumn(column)) {
      # sparse columns ride to the engine unexpanded: a sparseFactor as one
      # categorical column over its level table, every other kind ordinal
      slices <- sparseColumnSlices(column, name, nrow(x))
      columns[[j]] <- slices
      columnIsSparse[j] <- TRUE
      if (methods::is(column, "sparseFactor")) {
        columnTypes[[j]] <- CATEGORICAL_VARIABLE
        columnLevels[[j]] <- list(column@levels)
      } else {
        columnTypes[[j]] <- rep.int(ORDINAL_VARIABLE, length(slices$i))
        columnLevels[[j]] <- rep.int(list(NULL), length(slices$i))
      }
      columnNames[[j]] <- slices$names
    } else if (is.factor(column)) {
      if (!is.ordered(column) && nlevels(column) > 65535L) {
        stop(
          "factor '",
          name,
          "' has more than 65535 levels, the most a ",
          "categorical predictor supports"
        )
      }
      # the factor itself; its codes materialize as doubles only at use
      columns[[j]] <- column
      columnTypes[[j]] <- if (is.ordered(column)) {
        ORDINAL_VARIABLE
      } else {
        CATEGORICAL_VARIABLE
      }
      columnLevels[[j]] <- list(levels(column))
      columnNames[[j]] <- name
    } else if (is.matrix(column)) {
      columns[[j]] <- matrix(as.double(column), nrow(column))
      columnTypes[[j]] <- rep.int(ORDINAL_VARIABLE, ncol(column))
      columnLevels[[j]] <- rep.int(list(NULL), ncol(column))
      columnNames[[j]] <- paste(
        name,
        if (!is.null(colnames(column))) {
          colnames(column)
        } else {
          seq_len(ncol(column))
        },
        sep = "."
      )
    } else if (is.numeric(column) || is.logical(column)) {
      columns[[j]] <- as.double(column)
      columnTypes[[j]] <- ORDINAL_VARIABLE
      columnLevels[[j]] <- list(NULL)
      columnNames[[j]] <- name
    } else {
      stop("column '", name, "' cannot be converted to a predictor")
    }
  }
  if (any(columnIsSparse)) {
    # the mixed flavor keeps its dense columns as a per-column list too, so
    # factors ride unexpanded and the bridge codes them at assembly time
    result <- assembleMixedMatrix(columns, columnIsSparse, columnNames, nrow(x))
  } else {
    result <- assembleDenseColumnMatrix(columns, columnNames, nrow(x))
  }
  attr(result, "term.labels") <- names(x)
  attr(result, "varTypes") <- unlist(columnTypes)
  # c() keeps NULL elements where unlist would drop them, so the list stays
  # aligned with the columns
  attr(result, "factor.levels") <- do.call(c, columnLevels)
  result
}

## Whether the predictor source carries CSC-backed columns (a dgCMatrix, or
## a mixed container with sparse-mapped columns) - the designs whose raw
## values the engine retains as slices and whose mutation surface is fixed
## at creation.
predictorSourceIsSparse <- function(x) {
  inherits(x, "dgCMatrix") ||
    (inherits(x, "dbartsMixedMatrix") && any(x$map < 0L))
}

## Whether predictor column 'column' (1-based) of a predictor source is
## CSC-backed - the columns whose values live in the sparse block rather than a
## dense vector, and so the ones whose mutation is whole-column only. Every
## column of a dgCMatrix is; a container reads the sign of its map; a plain
## matrix has none. An out-of-range column reports FALSE, leaving the range
## error to the caller that knows the bound.
predictorColumnIsSparseBacked <- function(x, column) {
  if (inherits(x, "dgCMatrix")) {
    return(TRUE)
  }
  if (!inherits(x, "dbartsMixedMatrix")) {
    return(FALSE)
  }
  isTRUE(x$map[column] < 0L)
}

## The dense predictor matrix the re-quantize/replay bridges consume from the
## R side (setCutPoints, setState, saved-tree getTrees, and the low-level
## handle's tracked source). A plain matrix passes through; a dense-backed
## container materializes to its numeric codes; a genuinely sparse source (a
## CSC-bearing container or a dgCMatrix) yields NULL, which the bridge reads
## as "re-quantize from the retained slices" - the same signal a non-real
## source gave before the container existed.
rawPredictorMatrix <- function(x) {
  if (is.matrix(x)) {
    x
  } else if (inherits(x, "dbartsMixedMatrix") && !predictorSourceIsSparse(x)) {
    as.matrix(x)
  } else {
    NULL
  }
}

## Starting sigma estimate from a linear fit. NAs in the predictors are
## mean-imputed for this estimate only: complete cases can be scarce when
## missingness is scattered, and the estimate just anchors the residual
## variance prior.
estimateSigmaFromLinearModel <- function(data) {
  x <- data@x
  # a sparse design would densify under lm and is typically wide anyway;
  # the marginal estimate still anchors the residual variance prior. A dense
  # container (a frame with factors) still fits lm; only CSC-backed columns
  # fall back to the marginal estimate.
  if (predictorSourceIsSparse(x)) {
    warning(warningCondition(
      paste0(
        "starting sigma estimate falls back to the marginal response sd: ",
        "'x' has sparse-backed predictor columns, so a linear-model ",
        "estimate is not attempted"
      ),
      class = c("dbartsSparseSigmaFallbackWarning", "dbartsWarning")
    ))
    residual <- if (!is.null(data@offset)) data@y - data@offset else data@y
    return(sd(residual))
  }
  x <- as.matrix(x)
  if (anyNA(x)) {
    for (j in seq_len(ncol(x))) {
      column <- x[, j]
      if (anyNA(column)) {
        column[is.na(column)] <- mean(column, na.rm = TRUE)
        x[, j] <- column
      }
    }
  }
  # summary.lm's "essentially perfect fit" warning is a BLAS/architecture
  # property of this internal fit (residual variance can land within a
  # factor of two of its 1e-30 threshold), leaks an implementation detail
  # meaningless at the dbarts() call site, and duplicates the user-facing
  # signal dbarts already gives for a degenerate response (the
  # "response values are indistinguishable..." warning in R/data.R). Muffle
  # only that exact condition; everything else passes through.
  withCallingHandlers(
    summary(lm(data@y ~ x, weights = data@weights, offset = data@offset))$sigma,
    warning = function(w) {
      if (grepl("essentially perfect fit", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

## Recode a sparseFactor test column over the training level table, keeping
## it sparse: the stored entries' labels and the implicit reference level
## are re-coded in training level order, so its codes agree with a dense
## factor of the same values. A level unseen in training has no code and
## errors, exactly as the dense factor path does.
remapSparseFactorToTrainingLevels <- function(column, trainingLevels, name) {
  storedCodes <- match(column@levels[column@values], trainingLevels)
  referenceCode <- match(column@reference, trainingLevels)
  if (anyNA(storedCodes) || is.na(referenceCode)) {
    stop(
      "test data factor '",
      name,
      "' has levels not present in the ",
      "training data"
    )
  }
  # the stored entries stay non-reference (their labels differ from the
  # reference label), so the canonical structure carries over untouched
  newValidated(
    "sparseFactor",
    i = column@i,
    values = as.integer(storedCodes),
    levels = trainingLevels,
    reference = column@reference,
    length = column@length
  )
}

## Recode a test data.frame's factor, character, and sparseFactor columns
## against the training data's level tables (aligned with the training
## columns by name), so codes agree across the two; a sparseFactor stays
## sparse. A level unseen in training has no code and errors.
mapFactorColumnsToTrainingLevels <- function(
  x.test,
  predictorNames,
  factorLevels
) {
  for (name in names(x.test)) {
    j <- match(name, predictorNames)
    if (is.na(j)) {
      next
    }
    column <- x.test[[name]]
    columnIsSparseFactor <- methods::is(column, "sparseFactor")
    if (is.null(factorLevels[[j]])) {
      # the training column was numeric; a factor/character/sparseFactor test
      # column here would otherwise fall through to makeCategoricalModelMatrix,
      # which recodes it against the test set's own level order instead
      if (is.factor(column) || is.character(column) || columnIsSparseFactor) {
        stop(
          "test column '",
          name,
          "' is ",
          class(column)[1L],
          " but the training column '",
          name,
          "' is numeric"
        )
      }
      next
    }
    if (columnIsSparseFactor) {
      x.test[[name]] <- remapSparseFactorToTrainingLevels(
        column,
        factorLevels[[j]],
        name
      )
      next
    }
    if (!is.factor(column) && !is.character(column)) {
      next
    }
    refactored <- factor(as.character(column), levels = factorLevels[[j]])
    if (anyNA(refactored) && !anyNA(column)) {
      stop(
        "test data factor '",
        name,
        "' has levels not present in the ",
        "training data"
      )
    }
    x.test[[name]] <- refactored
  }
  x.test
}

## The indicators route stores no level table; it replays the training drop
## pattern positionally, and that pattern is sized by the table each training
## column carried - a factor's per-level instance counts, a matrix's per-column
## flags. A test column whose table is longer indexes off the end of it, so the
## disagreement is named here rather than left to the column count the replay
## happens to produce. The shorter direction indexes in bounds and is left to
## the column-count and name checks downstream: a training level the data never
## observed contributes no column, so a test frame that simply lacks it still
## aligns.
refuseWiderTestColumns <- function(x.test, drop) {
  for (j in seq_along(x.test)) {
    column <- x.test[[j]]
    trained <- length(drop[[j]])
    if (is.character(column)) {
      column <- as.factor(column)
    }
    if (is.factor(column)) {
      if (nlevels(column) > trained) {
        stop(
          "'test' factor '",
          names(x.test)[j],
          "' declares ",
          nlevels(column),
          " levels but the training design declared ",
          trained,
          "; use bart2() or dbarts(), which track levels across predict ",
          "by default"
        )
      }
    } else if (!is.null(dim(column)) && ncol(column) > trained) {
      stop(
        "'test' column '",
        names(x.test)[j],
        "' has ",
        ncol(column),
        " columns but the training design expanded ",
        trained
      )
    }
  }
  invisible(NULL)
}

## Re-code an assembled predictor container's factor columns against the
## training level tables, so a container built elsewhere - carrying whatever
## level order its own data implied - means by each code what the training
## design means. remapSparseFactorToTrainingLevels' treatment generalized from
## a loose sparseFactor to a container: a CSC-backed column's stored codes, its
## reference code, and its declared level count all lift together, and the
## stored pattern is left alone (re-coding is a bijection on levels, so a
## non-reference entry stays non-reference); a dense-backed column is a real
## factor in $dense and re-levels directly. Columns are matched to the training
## design the way validateXTest's reorder matches them: by name when both sides
## carry unique names, by position otherwise. A level the training data never
## saw has no code and errors, exactly as the dense factor path does.
alignContainerFactorLevels <- function(x.test, predictorNames, factorLevels) {
  containerLevels <- attr(x.test, "factor.levels")
  if (is.null(containerLevels)) {
    return(x.test)
  }
  containerNames <- x.test$columnNames
  positions <- if (
    !is.null(predictorNames) &&
      !is.null(containerNames) &&
      anyDuplicated(predictorNames) == 0L
  ) {
    match(predictorNames, containerNames)
  } else {
    seq_along(factorLevels)
  }

  for (j in seq_along(factorLevels)) {
    trainingLevels <- factorLevels[[j]]
    position <- positions[j]
    if (
      is.null(trainingLevels) ||
        is.na(position) ||
        position > length(x.test$map)
    ) {
      next
    }
    oldLevels <- containerLevels[[position]]
    if (is.null(oldLevels) || identical(oldLevels, trainingLevels)) {
      next
    }
    name <- if (!is.null(predictorNames)) predictorNames[j] else as.character(j)
    source <- x.test$map[position]
    if (source > 0L) {
      column <- x.test$dense[[source]]
      if (!is.factor(column)) {
        next
      }
      refactored <- factor(as.character(column), levels = trainingLevels)
      if (anyNA(refactored) && !anyNA(column)) {
        stop(
          "test data factor '",
          name,
          "' has levels not present in the training data"
        )
      }
      x.test$dense[[source]] <- refactored
    } else {
      rank <- -source
      pointers <- x.test$sparse@p
      entries <- seq.int(
        pointers[rank] + 1L,
        length.out = pointers[rank + 1L] - pointers[rank]
      )
      codes <-
        match(oldLevels[x.test$sparse@x[entries] + 1], trainingLevels) - 1L
      reference <-
        match(oldLevels[x.test$sparseReference[rank] + 1L], trainingLevels) - 1L
      if (anyNA(codes) || is.na(reference)) {
        stop(
          "test data factor '",
          name,
          "' has levels not present in the training data"
        )
      }
      # slot surgery, never Matrix's `[<-`: its drop0 is matrix-wide and would
      # strip the explicit zeros the container's other columns hold
      x.test$sparse@x[entries] <- as.double(codes)
      x.test$sparseReference[rank] <- as.integer(reference)
      x.test$sparseCategoryCount[rank] <- length(trainingLevels)
    }
    containerLevels[[position]] <- trainingLevels
  }
  attr(x.test, "factor.levels") <- containerLevels
  x.test
}

## Refuse a sparse test column's declared reference level when the training
## column it maps to is not categorical (attr(x.train, "varTypes")): a
## container's as.matrix (mixedMatrix.R) takes that column's implicit rows to
## be the reference whenever it is non-NA, but the engine takes them to be the
## reference only for a store-CATEGORICAL column and 0 otherwise - two
## different densifications of the same container. Positions match the way
## validateXTest's later reorder does: by name when both sides carry unique names, by position
## otherwise. The reverse mismatch (no reference against a store-categorical
## column) is refused downstream, inside resolveCscCategoricalReferences.
refuseSparseTestReferenceAgainstTrainTypes <- function(
  x.test,
  predictorNames,
  varTypes
) {
  if (is.null(x.test$sparse) || ncol(x.test$sparse) == 0L) {
    return(invisible(NULL))
  }
  containerNames <- x.test$columnNames
  positions <- if (
    !is.null(predictorNames) &&
      !is.null(containerNames) &&
      anyDuplicated(predictorNames) == 0L
  ) {
    match(predictorNames, containerNames)
  } else {
    seq_along(varTypes)
  }

  for (j in seq_along(varTypes)) {
    if (isTRUE(varTypes[j] == CATEGORICAL_VARIABLE)) {
      next
    }
    position <- positions[j]
    if (is.na(position) || position > length(x.test$map)) {
      next
    }
    source <- x.test$map[position]
    if (source > 0L || is.na(x.test$sparseReference[-source])) {
      next
    }
    stop(
      "a sparse predictor column may declare a reference level only for a ",
      "categorical predictor"
    )
  }
  invisible(NULL)
}

## use this to produce calls of the form
##  dbarts:::functionName
## so that we can evaluate non-exported functions in
## the user's environment
quoteInNamespace <- function(name, character.only = FALSE) {
  result <- quote(a + b)
  result[[1L]] <- as.symbol(":::")
  result[[2L]] <- as.symbol("dbarts")

  result[[3L]] <- if (character.only) name else match.call()[[2]]
  result
}
