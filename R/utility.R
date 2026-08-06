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

validateObject <- function(object) {
  tryCatch(validObject(object), error = rethrowValidityError)
  invisible(object)
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
  summary(lm(data@y ~ x, weights = data@weights, offset = data@offset))$sigma
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
