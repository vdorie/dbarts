coerceOrError <- function(x, type)
{
  func <- switch(type, logical = as.logical, integer = as.integer)
  result <- tryCatch(func(x), warning = function(e) e)
  if (is(result, "warning")) stop(paste0("'", match.call()[[2]], "' must be coercible to type: ", type))
  
  result
}
  
prepareCallWithArguments <- function(call, name, ...)
{
  argsToKeep <- unlist(list(...))
  matchIndices <- match(argsToKeep, names(call), nomatch = 0L)
  
  call <- call[c(1L, matchIndices)]
  call[[1]] <- name
  call
}

addCallArgument <- function(call, position, argument)
{
  if (is.character(position)) {
    name <- position
    position <- length(call) + 1L
  } else {
    position <- as.integer(position) + 1L
    if (position <= length(call)) for (i in length(call):position) {
      call[[i + 1]] <- call[[i]]
      names(call)[[i + 1]] <- names(call)[[i]]
    }
    name <- ""
  }
  call[[position]] <- argument
  names(call)[[position]] <- name
  call
}

setDefaultsFromFormals <- function(call, formals, ...)
{
  argsToReplace <- list(...)
  matchIndices <- match(argsToReplace, names(call), nomatch = 0L)
  missingFormals <- match(argsToReplace[matchIndices == 0L], names(formals))

  if (length(missingFormals) == 0) return(call)
  
  call[seq.int(length(missingFormals)) + length(call)] <- formals[missingFormals]
  call
}

is.formula <- function(x) is.language(x) && x[[1]] == '~'

## from lme4
namedList <- function(...) {
  result <- list(...)
  substituteNames <- sapply(substitute(list(...)), deparse)[-1]
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
  
  result <- .Call("dbarts_makeModelMatrixFromDataFrame", x, drop)
  attr(result, "term.labels") <- names(x)
  result
}

## use this to produce calls of the form
##  dbarts:::functionName
## so that we can evaluate non-exported functions in
## the user's environment
quoteInNamespace <- function(name) {
  result <- quote(a + b)
  result[[1]] <- as.symbol(":::")
  result[[2]] <- as.symbol("dbarts")
  result[[3]] <- match.call()[[2]]
  result
}
