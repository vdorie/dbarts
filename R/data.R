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

parseData <- function(formula, data, test, subset, weights, offset)
{
  ## default case of fn(y ~ x, data, ...)
  if (is.language(formula) && formula[[1]] == '~') {

    modelFrameCall <- match.call()

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
      if (!is.numeric(weights)) stop("'weights' must be numeric vector.")
      weights <- rep_len(weights, numObservations)
    }
    
    offset <- as.vector(model.offset(modelFrame))
    if (!is.null(offset) && length(offset) != numObservations) stop("Length of offset must be equal to that of y.")

    modelTerms <- terms(modelFrame)
    if (is.empty.model(modelTerms)) stop("Covariates must be specified for regression tree analysis.")

    attr(modelTerms, "intercept") <- 0L

    termIsFactor <- sapply(modelFrame, is.factor)
    numFactorTerms <- sum(termIsFactor)
    contrasts <-
      if (numFactorTerms == 0) NULL else lapply(modelFrame[,termIsFactor], contrasts, contrasts = FALSE)
    
    x <- model.matrix(modelTerms, modelFrame, contrasts)
    
  } else if (is.numeric(formula) || is.data.frame(formula)) {
    ## backwards compatibility of bart(x.train, y.train, x.test)
    if (missing(data) || is.null(data)) data <- rep(0, NROW(formula))
    if (!is.numeric(data) && !is.data.frame(data) && !is.integer(data) && !is.factor(data)) stop("When 'formula' is numeric, 'data' must be numeric as well.")

    if (is.factor(data)) {
      y <- as.numeric(as.integer(data) - 1L)
    } else {
      y <- as.numeric(data)
    }
    if (NROW(formula) != NROW(y))
      stop("'x' must have the same number of observations as 'y'.")
    initialNumObservations <- NROW(y)
    
    if (missing(subset) || is.null(subset)) subset <- seq.int(length(y))
    y <- y[subset]

    if (is.data.frame(formula)) formula <- makeModelMatrixFromDataFrame(formula)
    x <- if (!is.matrix(formula)) formula[subset] else formula[subset,]
    
    
    if (missing(weights)) weights <- NULL
    if (!is.null(weights)) {
      if (!is.numeric(weights)) stop("'weights' must be a numeric vector.")
      weights <- rep_len(weights, initialNumObservations)[subset]
    }

    if (missing(offset)) offset <- NULL
    if (!is.null(offset)) {
      if (!is.numeric(offset)) stop("Offset must be a numeric vector.")
      if (length(offset) == 1) offset <- rep_len(offset, initialNumObservations)
      if (length(offset) != initialNumObservations) stop("Length of 'offset' must equal length of 'y'.")
      offset <- offset[subset]
    }
  } else {
    stop("Unrecognized formula for data. Must be coercible to numeric types.")
  }

  if (is.vector(x)) x <- as.matrix(x)
  if (is.data.frame(x)) x <- makeModelMatrixFromDataFrame(x)

  
  x.test <- NULL
  if (!missing(test) && !is.null(test) && NCOL(test) > 0) {
    if (is.data.frame(test)) test <- makeModelMatrixFromDataFrame(test)
    if (!is.matrix(test)) test <- as.matrix(test)

    numPredictors <- ncol(x)
    
    if (!identical(ncol(test), numPredictors))
      stop("Number of columns in 'test' must be equal to that of 'x'.")
    if (numPredictors == 1) {
      x.test <- test
    } else {
      xIsNamed    <- !is.null(colnames(x))
      testIsNamed <- !is.null(colnames(test))

      columnIndices <- seq.int(numPredictors)
      if ((xIsNamed && !testIsNamed) || (!xIsNamed && testIsNamed)) {
        warning("'x' and 'test' are not both named. Columns of test matrix will be selected by position.")
      } else if (xIsNamed && testIsNamed){
        matchIndices <- match(colnames(x), colnames(test))
        if (any(is.na(matchIndices))) {
          warning("Column names of 'test' does not equal that of 'x': '", toString(colnames(x)),
                  "'. Match will be made by position.")
          columnIndices <- matchIndices
        }
      }
      
      x.test <- test[, columnIndices]
      if (xIsNamed) colnames(x.test) <- colnames(x)
    }
  }
  
  list(y = y, x = x, x.test = x.test, weights = weights, offset = offset)
}
