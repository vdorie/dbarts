coerceOrError <- function(x, type)
{
  mc <- match.call()
  
  if (is.null(x)) stop("'", mc[[2L]], "' cannot be NULL")
  
  func <- switch(type, logical = as.logical, integer = as.integer, numeric = as.numeric)
  result <- tryCatch(func(x), warning = function(e) e)
  if (inherits(result, "warning")) stop("'", mc[[2L]], "' must be coercible to type: ", type)
  
  result
}

"%not_in%" <- function(x, table) match(x, table, nomatch = 0L) <= 0L

evalx.recurse <- function(x, e) {
  if (length(e) == 0L || typeof(e) == "symbol") return(e)
  
  for (i in seq_along(e)) {
    if (!is.language(e[[i]])) next
    
    e[[i]] <- if (e[[i]] == "x") x else evalx.recurse(x, e[[i]])
  }
  
  e
}

ifelse_3 <- function(a, b, c, d, e) {
  mc <- match.call(); env <- parent.frame()
  if (eval(mc[["a"]], env)) {
    c
  } else if (eval(mc[["b"]], env)) {
    d
  } else {
    e
  }
}

ifelse_n <- function(n, ...) {
  mc <- match.call(); env <- parent.frame()
  
  for (i in seq_len(n - 1L))
    if (eval(mc[[i + 2L]], env)) return(eval(mc[[n + 1L + i]], env))
  eval(mc[[2L * n - 1L]], env)
}

## evaluates the expression 'e' by after first replacing all instances of 'x' with the expression x
#x <- NULL
evalx <- function(x, e) {
  mc <- match.call()
  callingEnv <- parent.frame()
  
  e <- evalx.recurse(mc$x, mc$e)
  eval(e, callingEnv)
}

redirectCall <- function(call, fn, ...)
{
  matchedCall <- match.call()
  extraArgs <- if (length(matchedCall) > 3L) as.character(matchedCall[-c(1L, 2L, 3L)]) else character()
  
  originalFn <- eval(call[[1L]])
  call[[1L]] <- if (is.function(fn)) matchedCall[[3L]] else fn
  if (length(extraArgs) == 0L) {
    fn <- if (is.function(fn)) fn else eval(fn)
    
    argsToKeep <- names(call)[-1L] %in% names(formals(fn))
    if (any(names(formals(originalFn)) == "...") && any(names(formals(fn)) == "..."))
      argsToKeep <- argsToKeep | names(call)[-1L] %not_in% names(formals(originalFn))
    
    call <- call[c(TRUE, argsToKeep)]
  } else {
    matchIndices <- match(extraArgs, names(call), nomatch = 0L)
    
    call <- call[c(1L, matchIndices)]
  }
  
  call
}

addCallArgument <- function(call, position, argument)
{
  if (is.character(position)) {
    name <- position
    position <- length(call) + 1L
  } else {
    position <- as.integer(position) + 1L
    if (position <= length(call)) for (i in seq.int(length(call), position)) {
      call[i + 1L] <- call[i]
      names(call)[[i + 1L]] <- names(call)[[i]]
    }
    name <- ""
  }
  call[[position]] <- argument
  names(call)[[position]] <- name
  call
}

subTermInLanguage <- function(lang, oldTerm, newTerm)
{
  if (length(lang) == 1L && is.symbol(lang))
    return(if (lang == oldTerm) newTerm else lang)
  
  for (i in seq_along(lang)) {
    if (is.symbol(lang[[i]])) {
      if (lang[[i]] == oldTerm) lang[[i]] <- newTerm
    } else if (is.language(lang[[i]])) {
      lang[[i]] <- subTermInLanguage(lang[[i]], oldTerm, newTerm)
    }
  }
  lang
}

setDefaultsFromFormals <- function(call, formals, ...)
{
  argsToReplace <- list(...)
  matchIndices <- match(argsToReplace, names(call), nomatch = 0L)
  missingFormals <- match(argsToReplace[matchIndices == 0L], names(formals))

  if (length(missingFormals) == 0L) return(call)
  
  newFormalIndices <- seq.int(length(missingFormals)) + length(call)
  call[newFormalIndices] <- formals[missingFormals]
  names(call)[newFormalIndices] <- names(formals)[missingFormals]
  call
}

is.formula <- function(x) is.language(x) && x[[1L]] == '~'

## from lme4
namedList <- function(...) {
  result <- list(...)
  substituteNames <- sapply(substitute(list(...)), deparse)[-1L]
  if (is.null(resultNames <- names(result))) resultNames <- substituteNames
  if (any(noNames <- resultNames == "")) resultNames[noNames] <- substituteNames[noNames]
  setNames(result, resultNames)
}

## Turns data.frame w/factors into matrices of indicator variables. Differs from
## model.matrix as it doesn't drop columns for co-linearity even with multiple
## factors 
makeModelMatrixFromDataFrame <- function(x, drop = TRUE) {
  if (!is.data.frame(x)) stop('x is not a dataframe')
  if (is.logical(drop) && (length(drop) != 1L || is.na(drop))) stop('when logical, drop must be TRUE or FALSE')
  if (is.list(drop) && length(drop) != length(x)) stop('when list, drop must have length equal to x')
  
  characterCols <- sapply(x, typeof) == "character"
  if (any(characterCols)) x[characterCols] <- lapply(x[characterCols], as.factor)
  
  result <- .Call(C_dbarts_makeModelMatrixFromDataFrame, x, drop)
  attr(result, "term.labels") <- names(x)
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
  if (!is.data.frame(x)) stop("x is not a dataframe")

  characterCols <- sapply(x, typeof) == "character"
  if (any(characterCols)) x[characterCols] <- lapply(x[characterCols], as.factor)

  columns      <- vector("list", length(x))
  columnTypes  <- vector("list", length(x))
  columnLevels <- vector("list", length(x))
  columnNames  <- vector("list", length(x))
  for (j in seq_along(x)) {
    column <- x[[j]]
    name <- names(x)[j]
    if (is.factor(column)) {
      if (!is.ordered(column) && nlevels(column) > 53L)
        stop("factor '", name, "' has more than 53 levels, the most a ",
             "categorical predictor supports")
      columns[[j]] <- matrix(as.double(as.integer(column) - 1L), ncol = 1L)
      columnTypes[[j]] <- if (is.ordered(column)) ORDINAL_VARIABLE else CATEGORICAL_VARIABLE
      columnLevels[[j]] <- list(levels(column))
      columnNames[[j]] <- name
    } else if (is.matrix(column)) {
      columns[[j]] <- matrix(as.double(column), nrow(column))
      columnTypes[[j]] <- rep.int(ORDINAL_VARIABLE, ncol(column))
      columnLevels[[j]] <- rep.int(list(NULL), ncol(column))
      columnNames[[j]] <- paste(
        name,
        if (!is.null(colnames(column))) colnames(column) else seq_len(ncol(column)),
        sep = "."
      )
    } else if (is.numeric(column) || is.logical(column)) {
      columns[[j]] <- matrix(as.double(column), ncol = 1L)
      columnTypes[[j]] <- ORDINAL_VARIABLE
      columnLevels[[j]] <- list(NULL)
      columnNames[[j]] <- name
    } else {
      stop("column '", name, "' cannot be converted to a predictor")
    }
  }
  result <- do.call(cbind, columns)
  colnames(result) <- unlist(columnNames)
  attr(result, "term.labels") <- names(x)
  attr(result, "varTypes") <- unlist(columnTypes)
  # c() keeps NULL elements where unlist would drop them, so the list stays
  # aligned with the columns
  attr(result, "factor.levels") <- do.call(c, columnLevels)
  result
}

## Recode a test data.frame's factor and character columns against the
## training data's level tables (aligned with the training columns by
## name), so codes agree across the two; a level unseen in training has no
## code and errors.
mapFactorColumnsToTrainingLevels <- function(x.test, predictorNames, factorLevels) {
  for (name in names(x.test)) {
    j <- match(name, predictorNames)
    if (is.na(j) || is.null(factorLevels[[j]])) next
    column <- x.test[[name]]
    if (!is.factor(column) && !is.character(column)) next
    refactored <- factor(as.character(column), levels = factorLevels[[j]])
    if (anyNA(refactored) && !anyNA(column))
      stop("test data factor '", name, "' has levels not present in the ",
           "training data")
    x.test[[name]] <- refactored
  }
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
