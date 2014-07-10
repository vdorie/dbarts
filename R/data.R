ORDINAL_VARIABLE <- 0L
CATEGORICAL_VARIABLE <- 1L

setMethod("initialize", "dbartsData",
          function(.Object, modelMatrices, n.cuts = 100L, sigma = NA_real_)
{
  if (!missing(modelMatrices)) {
    .Object@y <- modelMatrices$y
    .Object@x <- modelMatrices$x
    .Object@varTypes <- rep.int(ORDINAL_VARIABLE, ncol(.Object@x))
    .Object@x.test <- modelMatrices$x.test
    .Object@weights <- modelMatrices$weights
    .Object@offset <- modelMatrices$offset
  }
  
  .Object@n.cuts <- rep_len(as.integer(n.cuts), ncol(.Object@x))
  .Object@sigma <- sigma

  validObject(.Object)
  .Object
})

validateXTest <- function(x.test, numPredictors, predictorNames)
{
  if (is.null(x.test)) return(x.test)
  
  if (is.data.frame(x.test)) x.test <- makeModelMatrixFromDataFrame(x.test)
  if (!is.matrix(x.test)) x.test <- as.matrix(x.test)

  if (!is.numeric(x.test))
    stop('test matrix must be numeric')

  if (is.integer(x.test)) x.test <- matrix(as.double(x.test), nrow(x.test))
  
  if (!identical(ncol(x.test), numPredictors))
    stop("number of columns in 'test' must be equal to that of 'x'")
  if (numPredictors > 1) {
    xIsNamed    <- !is.null(predictorNames)
    testIsNamed <- !is.null(colnames(x.test))
    
    columnIndices <- seq.int(numPredictors)
    if ((xIsNamed && !testIsNamed) || (!xIsNamed && testIsNamed)) {
      ## warning("'x' and 'test' are not both named; columns of test matrix will be selected by position")
    } else if (xIsNamed && testIsNamed) {
      matchIndices <- match(predictorNames, colnames(x.test))
      if (any(is.na(matchIndices))) {
        warning("column names of 'test' does not equal that of 'x': '", toString(predictorNames),
                "'; match will be made by position")
      } else {
        columnIndices <- matchIndices
      }
    }
    
    x.test <- x.test[, columnIndices]
    if (xIsNamed) colnames(x.test) <- predictorNames
  }
  
  x.test
}

parseData <- function(formula, data, test, subset, weights, offset)
{
  dataIsMissing <- missing(data)
  testIsMissing <- missing(test)
  matchedCall <- match.call()
  
  ## default case of fn(y ~ x, data, ...)
  if (is.language(formula) && formula[[1]] == '~') {

    modelFrameCall <- matchedCall

    ## remove "test"
    matchPositions <- match(c("formula", "data", "subset", "weights", "offset"),
                            names(modelFrameCall), nomatch = 0L)
    
    modelFrameCall <- modelFrameCall[c(1L, matchPositions)]
    modelFrameCall$drop.unused.levels <- TRUE
    modelFrameCall[[1L]] <- quote(stats::model.frame)
    
    modelFrame <- eval(modelFrameCall, parent.frame())

    y <- model.response(modelFrame, "numeric")
    if (is.null(y)) y <- rep(0, NROW(modelFrame))
    numObservations <- NROW(y)
    
    weights <- as.vector(model.weights(modelFrame))
    if (!is.null(weights)) {
      if (!is.numeric(weights)) stop("'weights' must be numeric vector")
      weights <- rep_len(weights, numObservations)
    }
    
    offset <- as.vector(model.offset(modelFrame))
    if (!is.null(offset) && length(offset) != numObservations) stop("length of offset must be equal to that of y")

    modelTerms <- terms(modelFrame)
    if (is.empty.model(modelTerms)) stop("covariates must be specified for regression tree analysis")

    attr(modelTerms, "intercept") <- 0L

    termIsFactor <- sapply(modelFrame, is.factor)
    numFactorTerms <- sum(termIsFactor)
    contrasts <-
      if (numFactorTerms == 0) NULL else lapply(modelFrame[,termIsFactor], contrasts, contrasts = FALSE)
    
    x <- model.matrix(modelTerms, modelFrame, contrasts)

    ## unless the test data has the same number of rows as the 
    if (!testIsMissing) {
      foundTest <- FALSE
      if (!dataIsMissing && is.list(data)) {
        tryResult <-
          tryCatch(temp <- eval(call("$", as.symbol("data"), matchedCall$test)), error = function(e) e)
        
        if (!is(tryResult, "error") && !is.null(temp)) {
          foundTest <- TRUE
          test <- temp
        }
      }
      if (!foundTest)
        test <- eval(matchedCall$test, environment(formula))
    }
    
  } else if (is.numeric(formula) || is.data.frame(formula)) {
    ## backwards compatibility of bart(x.train, y.train, x.test)
    if (dataIsMissing || is.null(data)) data <- rep(0, NROW(formula))
    if (!is.numeric(data) && !is.data.frame(data) && !is.integer(data) && !is.factor(data)) stop("when 'formula' is numeric, 'data' must be numeric as well")

    if (is.factor(data)) {
      y <- as.numeric(as.integer(data) - 1L)
    } else {
      y <- as.numeric(data)
    }
    if (NROW(formula) != NROW(y))
      stop("'x' must have the same number of observations as 'y'")
    initialNumObservations <- NROW(y)
    
    if (missing(subset) || is.null(subset)) subset <- seq.int(length(y))
    y <- y[subset]

    if (is.data.frame(formula)) formula <- makeModelMatrixFromDataFrame(formula)
    x <- if (!is.matrix(formula)) formula[subset] else formula[subset,]
    
    
    if (missing(weights)) weights <- NULL
    if (!is.null(weights)) {
      if (!is.numeric(weights)) stop("'weights' must be a numeric vector")
      weights <- rep_len(weights, initialNumObservations)[subset]
    }

    if (missing(offset)) offset <- NULL
    if (!is.null(offset)) {
      if (!is.numeric(offset)) stop("'offset' must be a numeric vector")
      if (length(offset) == 1) offset <- rep_len(offset, initialNumObservations)
      if (length(offset) != initialNumObservations) stop("length of 'offset' must equal length of 'y'")
      offset <- offset[subset]
    }
  } else {
    stop("unrecognized formula for data; must be coercible to numeric types")
  }

  if (is.vector(x)) x <- as.matrix(x)
  if (is.data.frame(x)) x <- makeModelMatrixFromDataFrame(x)

  
  x.test <- NULL
  if (!testIsMissing && !is.null(test) && NCOL(test) > 0)
    x.test <- validateXTest(test, ncol(x), colnames(x))
  
  list(y = y, x = x, x.test = x.test, weights = weights, offset = offset)
}
