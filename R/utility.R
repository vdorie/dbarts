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

evalx.recurse <- function(x, e) {
  if (length(e) == 0L || typeof(e) == "symbol") return(e)
  
  for (i in seq_along(e)) {
    if (!is.language(e[[i]])) next
    
    e[[i]] <- if (e[[i]] == "x") x else evalx.recurse(x, e[[i]])
  }
  
  e
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
  if (is.logical(drop) && (length(drop) != 1L || is.na(drop))) stop('when logical, drop must be TRUE or FALSE')
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
