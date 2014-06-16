stripExtraArguments <- function(call, ...)
{
  argsToKeep <- unlist(list(...))
  matchIndices <- match(argsToKeep, names(call), nomatch = 0L)
  
  call <- call[c(1L, matchIndices)]
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
## factors. Also drops any non-numeric, non-factor columns.
makeModelMatrixFromDataFrame <- function(x) {
  if (!is.data.frame(x))
    stop("x in makeModelMatrixFromDataFrame is not a dataframe")
  

  factorToMatrix <- function(factor) {
    sapply(levels(factor), function(level) ifelse(factor == level, 1, 0))
  }

  numTerms <- length(x)
  termNames <- names(x)

  termMatrices <- lapply(seq.int(numTerms), function(index) {
    x_i <- x[[index]]
    if (is.factor(x_i)) {
      result <- factorToMatrix(x_i)
      colnames(result) <- paste0(termNames[index], ".", levels(x_i))
      return(result)
    }
    if (is.matrix(x_i)) {
      result <- x_i
      colnames(result) <- paste0(termNames[index], ".", seq.int(ncol(x_i)))
      return(result)
    }

    return(matrix(x_i, length(x_i),
                  dimnames = list(NULL, termNames[index])))
  })

  columnNames <- unlist(lapply(termMatrices, colnames))

  return(matrix(unlist(termMatrices), nrow(x), dimnames = list(NULL, columnNames)))
}
