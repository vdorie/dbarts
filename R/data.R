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
    .Object@offset.test <- modelMatrices$offset.test

    .Object@testUsesRegularOffset <- modelMatrices$testUsesRegularOffset
  }
  
  .Object@n.cuts <- rep_len(as.integer(n.cuts), ncol(.Object@x))
  .Object@sigma <- sigma

  validObject(.Object)
  .Object
})

validateXTest <- function(x.test, termLabels, numPredictors, predictorNames, drop)
{
  if (is.null(x.test)) return(x.test)
  if (is.data.frame(x.test)) {
    if (!is.null(termLabels))
      x.test <- model.frame(formula = as.formula(paste("~", paste(termLabels, collapse = " + "))), data = x.test)
    x.test <- makeModelMatrixFromDataFrame(x.test, if (!is.null(drop)) drop else TRUE)
  }
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

parseData <- function(formula, data, test, subset, weights, offset, offset.test = offset)
{
  dataIsMissing <- missing(data)
  testIsMissing <- missing(test)
  offsetIsMissing <- missing(offset)
  matchedCall <- match.call()

  offsetGivenAsScalar <- NA
  testUsesRegularOffset <- NA

  ## default case of fn(y ~ x, data, ...)
  if (is.language(formula) && formula[[1]] == '~') {

    ## remove "test" and offset.test before calling model.frame
    modelFrameArgs <- c("formula", "data", "subset", "weights", "offset")

    ## pull out offset if it is a scalar, else model.frame will bug out
    if (!dataIsMissing && !offsetIsMissing) {
      offsetFound <- FALSE
      tempOffset <- NULL
      if (is.data.frame(data)) {
        if (is.symbol(matchedCall$offset) && any(names(data) == matchedCall$offset)) {
          ## model.frame can find it later
          offsetFound <- TRUE
        }
      } else {
        ## data should be a list, which means "offset" might be a member of it
        tryResult <- tryCatch(tempOffset <- eval(call("$", as.symbol("data"), matchedCall$offset)), error = function(e) e)
        if (!is(tryResult, "error") && !is.null(tempOffset)) offsetFound <- TRUE
      }
      if (!offsetFound) {
        tempOffset <- eval(matchedCall$offset, environment(formula))
        if (!is.null(tempOffset)) offsetFound <- TRUE
      }
      if (offsetFound && !is.null(tempOffset)) {
        offset <- tempOffset
        if (!is.numeric(offset)) stop("'offset' must be numeric")
        offsetGivenAsScalar <- length(offset) == 1
        modelFrameArgs <- c("formula", "data", "subset", "weights")
      }
    }
    modelFrameCall <- matchedCall

    matchPositions <- match(modelFrameArgs, names(modelFrameCall), nomatch = 0L)
    modelFrameCall <- modelFrameCall[c(1L, matchPositions)]
    modelFrameCall$drop.unused.levels <- FALSE
    modelFrameCall[[1L]] <- quote(stats::model.frame)
    
    modelFrame <- eval(modelFrameCall, parent.frame())
    if (nrow(modelFrame) == 0) {
      if (!is.null(matchedCall$subset)) stop("invalid 'subset'")
      stop("empty data argument")
    }

    y <- model.response(modelFrame, "numeric")
    if (is.null(y)) y <- rep(0, NROW(modelFrame))
    numObservations <- NROW(y)
    
    weights <- as.vector(model.weights(modelFrame))
    if (!is.null(weights)) {
      if (!is.numeric(weights)) stop("'weights' must be numeric vector")
      weights <- rep_len(weights, numObservations)
    }

    if (is.na(offsetGivenAsScalar)) {
      offset <- as.vector(model.offset(modelFrame))
      if (!is.null(offset)) {
        if (length(offset) != numObservations) stop("length of offset must be equal to that of y")
        offsetGivenAsScalar <- FALSE
      }
    } else {
      offset <- rep_len(offset, numObservations)
    }

    modelTerms <- terms(modelFrame)
    if (is.empty.model(modelTerms)) stop("covariates must be specified for regression tree analysis")
    
    x <- makeModelMatrixFromDataFrame(modelFrame[attr(modelTerms, "term.labels")])
    
    if (!testIsMissing) {
      foundTest <- FALSE
      if (!dataIsMissing && !is.data.frame(data)) {
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

    if (offsetIsMissing) offset <- NULL
    if (!is.null(offset)) {
      if (!is.numeric(offset)) stop("'offset' must be numeric")
      if (length(offset) == 1) {
        offset <- rep_len(offset, initialNumObservations)
        offsetGivenAsScalar <- TRUE
      } else {
        offsetGivenAsScalar <- FALSE
      }
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
    x.test <- validateXTest(test, attr(x, "term.labels"), ncol(x), colnames(x), attr(x, "drop"))

  if (!is.null(x.test)) {
    if (missing(offset.test)) {
      ## default is offset.test = offset
      if (identical(offsetGivenAsScalar, TRUE)) {
        offset.test <- rep_len(offset[1], nrow(x.test))
      } else if (identical(offsetGivenAsScalar, FALSE)) {
        if (nrow(x.test) != length(y)) stop("vectored 'offset' cannot be directly applied to test data of unequal length")
        offset.test <- offset
      }
      
      if (!is.na(offsetGivenAsScalar)) testUsesRegularOffset <- TRUE
    } else if (!is.null(matchedCall$offset.test)) {
      ## test.offset could have been something like (offset + 0.5), and we would just have redefined offset
      if (is.language(matchedCall$offset.test)) {


        ## we can also have wierdness such as (offset + variable), where offset is now in our
        ## our environment but the variable is in the caller
        testReferencesOffset <- any(as.list(matchedCall$offset.test) == "offset")
        
        evalEnv <- if (testReferencesOffset) {
          result <- new.env(parent = parent.frame())
          result$offset <- if (offsetGivenAsScalar == TRUE) offset[1] else offset
          result
        } else {
          parent.frame()
        }
        
        offset.test <- eval(matchedCall$offset.test, evalEnv)
      }
        
      if (length(offset.test) == 1) {
        offset.test <- rep_len(offset.test, nrow(x.test))
      }
      
      testUsesRegularOffset <- FALSE

      ## for when, for whatever reason, it is manually passed in "offset.test = offset"
      if (identical(matchedCall$offset.test, quote(offset))) {
        testUsesRegularOffset <- TRUE
      }
    } else {
      testUsesRegularOffset <- FALSE
    }
  } else {
    if (missing(offset.test)) offset.test <- NULL
  }
  
  list(y = y, x = x, x.test = x.test, weights = weights, offset = offset, offset.test = offset.test,
       testUsesRegularOffset = testUsesRegularOffset)
}
