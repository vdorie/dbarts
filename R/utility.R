coerceOrError <- function(x, type)
{
  mc <- match.call()
  
  if (is.null(x)) stop("'", mc[[2L]], "' cannot be NULL")
  
  func <- switch(type, logical = as.logical, integer = as.integer, numeric = as.numeric)
  result <- tryCatch(func(x), warning = function(e) e)
  if (is(result, "warning")) stop("'", mc[[2L]], "' must be coercible to type: ", type)
  
  result
}

"%not_in%" <- function(x, table) match(x, table, nomatch = 0L) <= 0L

evalx <- function(x, e) {
  mc <- match.call()
  callingEnv <- parent.frame()
  evalEnv <- new.env(parent = callingEnv)
  evalEnv$x <- x
  eval(mc$e, evalEnv)
}

prepareCallWithArguments <- function(call, fn, ...)
{
  matchedCall <- match.call()
  argsToKeep <- as.character(matchedCall[-c(1L, 2L, 3L)])
  matchIndices <- match(argsToKeep, names(call), nomatch = 0L)
  
  call <- call[c(1L, matchIndices)]
  call[[1L]] <- if (is.function(fn)) matchedCall[[3L]] else fn
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
  for (i in seq_along(lang)) {
    if (is.symbol(lang[[i]])) {
      if (lang[[i]] == oldTerm) lang[[i]] <- newTerm
    } else if (is.language(lang[[i]])) {
      lang[[i]] <- subTermInLanguage(lang[[i]], oldTerm, newTerm)
    }
  }
  return(lang)
}

setDefaultsFromFormals <- function(call, formals, ...)
{
  argsToReplace <- list(...)
  matchIndices <- match(argsToReplace, names(call), nomatch = 0L)
  missingFormals <- match(argsToReplace[matchIndices == 0L], names(formals))

  if (length(missingFormals) == 0L) return(call)
  
  call[seq.int(length(missingFormals)) + length(call)] <- formals[missingFormals]
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
  if (is.logical(drop) && is.na(drop)) stop('when logical, drop must be TRUE or FALSE')
  if (is.list(drop) && length(drop) != length(x)) stop('when list, drop must have length equal to x')
  
  result <- .Call(C_dbarts_makeModelMatrixFromDataFrame, x, drop)
  attr(result, "term.labels") <- names(x)
  result
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
