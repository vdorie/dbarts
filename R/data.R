ORDINAL_VARIABLE <- 0L
CATEGORICAL_VARIABLE <- 1L

methods::setMethod("initialize", "dbartsData",
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
  if (is.numeric(x.test) && NCOL(x.test) == 0L) return(NULL)
  if (is.data.frame(x.test)) {
    if (!is.null(termLabels))
      x.test <- model.frame(formula = as.formula(paste("~", paste(termLabels, collapse = " + "))), data = x.test)
    x.test <- makeModelMatrixFromDataFrame(x.test, if (!is.null(drop)) drop else TRUE)
  }
  if (!is.matrix(x.test)) x.test <- as.matrix(x.test)

  if (!is.numeric(x.test))
    stop('test matrix must be numeric')

  if (is.integer(x.test)) x.test <- matrix(as.double(x.test), nrow(x.test))
  
  if (!identical(NCOL(x.test), numPredictors))
    stop("number of columns in 'test' must be equal to that of 'x'")
  if (numPredictors > 1) {
    xIsNamed    <- !is.null(predictorNames)
    testIsNamed <- !is.null(colnames(x.test))
    
    columnIndices <- seq.int(numPredictors)
    if ((xIsNamed && !testIsNamed) || (!xIsNamed && testIsNamed) || length(unique(predictorNames)) != length(predictorNames)) {
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
    
    x.test <- x.test[,columnIndices]
    if (xIsNamed) colnames(x.test) <- predictorNames
  }
  
  x.test
}

findTermInFormulaData <- function(formula, data, term)
{
  formulaIsMissing <- missing(formula)
  dataIsMissing <- missing(data)
  matchedCall <- match.call()
  
  if (is.numeric(matchedCall$term)) return(term)
  
  if (!dataIsMissing) {
    if (is.symbol(matchedCall$term)) {
      if (any(names(data) == as.character(matchedCall$term))) return(data[[as.character(matchedCall$term)]])
    } else if (is.language(matchedCall$term)) {
      #attach(data, warn.conflicts = FALSE, name = ".dbartsData_data")
      tryResult <- with(data, tryCatch(eval(matchedCall$term), error = function(e) e))
      #detach(data)
      if (!is(tryResult, "error")) return(tryResult)
    }
  }
  if (is.symbol(matchedCall$term)) {
    if (any(ls(environment(formula)) == as.character(matchedCall$term))) return(get(as.character(matchedCall$term), envir = environment(formula)))
    tryResult <- tryCatch(get(as.character(matchedCall$term)), error = function(e) e)
    if (!is(tryResult, "error") && !is.null(tryResult)) return(tryResult)
  } else if (is.language(matchedCall$term)) {
    tryResult <- tryCatch(eval(matchedCall$term, environment(formula)), error = function(e) e)
    if (!is(tryResult, "error")) return(tryResult)
    tryResult <- tryCatch(eval(matchedCall$term), error = function(e) e)
    if (!is(tryResult, "error")) return(tryResult)
  }
  
  NULL
}

## this used to be a function evaluated in the caller's frame, but
## that causes warnings in R check so now it is just a block of code
getTestOffset <- quote({
  if (is.numeric(matchedCall$offset.test))
    return(namedList(offset.test, testUsesRegularOffset = FALSE))
  if (is.null(matchedCall$offset.test))
    return(list(offset.test = NULL, testUsesRegularOffset = FALSE))
  
  if (is.symbol(matchedCall$offset.test)) {
    testOffsetName <- as.character(matchedCall$offset.test)
    
    if (identical(testOffsetName, "offset") && !is.null(offset))
      return(list(offset.test = if (offsetGivenAsScalar == TRUE) offset[1] else offset, testUsesRegularOffset = TRUE))
    
    if (is.formula(formula)) {
      if (!dataIsMissing && any(names(data) == testOffsetName))
        return(list(offset.test = data[[testOffsetName]], testUsesRegularOffset = FALSE))
      if (any(ls(environment(formula)) == testOffsetName))
        return(list(offset.test = get(testOffsetName, environment(formula)), testUsesRegularOffset = FALSE))
    }
    tryResult <- tryCatch(get(testOffsetName), error = function(e) e)
    if (!is(tryResult, "error") && !is.null(tryResult))
      return(list(offset.test = tryResult, testUsesRegularOffset = FALSE))
    
    stop("cannot find test offset '", testOffsetName, "'")
  } else if (is.language(matchedCall$offset.test)) {
    ## test.offset could have been something like (offset + 0.5), or (offset + variable)
    baseOffset <- if (is.null(offset)) NA_real_ else { if (offsetGivenAsScalar == TRUE) offset[1] else offset }
    
    if (identical(matchedCall$offset.test, quote(offset)))
      return(list(offset.test = baseOffset, testUsesRegularOffset = TRUE))
    
    testOffset <- subTermInLanguage(matchedCall$offset.test, quote(offset), baseOffset)

    if (is.formula(formula)) {
      if (!dataIsMissing) {
        #attach(data)
        tryResult <- with(data, tryCatch(eval(testOffset), error = function(e) e))
        #detach(data)
        if (!is(tryResult, "error")) return(list(offset.test = tryResult, testUsesRegularOffset = FALSE))
      }
      tryResult <- tryCatch(eval(testOffset, environment(formula)), error = function(e) e)
      if (!is(tryResult, "error")) return(list(offset.test = tryResult, testUsesRegularOffset = FALSE))
    }
    tryResult <- tryCatch(eval(testOffset, parent.frame(3L)), error = function(e) e)
    if (!is(tryResult, "error")) return(list(offset.test = tryResult, testUsesRegularOffset = FALSE))
  }
  
  stop("cannot construct test offset")
})

dbartsData <- function(formula, data, test, subset, weights, offset, offset.test = offset)
{
  dataIsMissing <- missing(data)
  testIsMissing <- missing(test)
  offsetIsMissing <- missing(offset)
  testOffsetIsMissing <- missing(offset.test)
  matchedCall <- match.call()
  
  offsetGivenAsScalar <- NA
  testUsesRegularOffset <- NA
  
  if (missing(formula)) stop("first argument to dbartsData - 'formula'/'x.train' - must be present")
  
  if (is(formula, "dbartsData")) {
    if (!dataIsMissing || !testIsMissing || !offsetIsMissing || !testOffsetIsMissing)
      warning("if data supplied as dbartsData, remaining arguments are ignored")
    return(formula)
  }
  
  if (is.formula(formula)) {
    if (!dataIsMissing && !is.data.frame(data) && !is.list(data) && !is.environment(data))
      stop("for formula/data specification, data must be a data frame, list, or environment")
    
    modelFrameArgs <- c("formula", "data", "subset", "weights", "offset")
    
    ## extract offset prematurely, if necessary
    if (offsetIsMissing) {
      offset <- NULL
      modelFrameArgs <- c("formula", "data", "subset", "weights")
    } else {
      offsetCall <- matchedCall
      offsetCall <- offsetCall[c(1L, match(c("formula", "data", "offset"), names(offsetCall), nomatch = 0L))]
      names(offsetCall)[which(names(offsetCall) == "offset")] <- "term"
      offsetCall[[1L]] <- quoteInNamespace(findTermInFormulaData)
      offset <- eval(offsetCall, parent.frame())
     
      if (!is.null(offset)) {
        offsetGivenAsScalar <- length(offset) == 1
        if (offsetGivenAsScalar) modelFrameArgs <- c("formula", "data", "subset", "weights")
      }      
      originalOffset <- offset
    }    
    modelFrameCall <- matchedCall
    modelFrameCall <- modelFrameCall[c(1L, match(modelFrameArgs, names(modelFrameCall), nomatch = 0L))]
    modelFrameCall$drop.unused.levels <- FALSE
    modelFrameCall$na.action <- stats::na.omit
    modelFrameCall[[1L]] <- quote(stats::model.frame)
    ## this allows subset to be applied to offset, even if offset was a language construct (e.g. off + 0.1)
    if (identical(offsetGivenAsScalar, FALSE)) modelFrameCall$offset <- offset
    
    modelFrame <- eval(modelFrameCall, parent.frame())
    if (NROW(modelFrame) == 0) {
      if (!is.null(matchedCall$subset)) stop("empty 'subset' specified")
      stop("cannot construct model matrices from formula")
    }
    
    ## pull out y
    y <- model.response(modelFrame, "numeric")
    if (is.null(y)) y <- rep(0, NROW(modelFrame))
    numObservations <- NROW(y)
    
    ## weights
    weights <- as.vector(model.weights(modelFrame))
    if (!is.null(weights) && !is.numeric(weights)) stop("'weights' must be of type numeric")
    
    ## offset, when in data frame
    if (identical(offsetGivenAsScalar, FALSE)) {
      offset <- as.vector(model.offset(modelFrame))
    } else if (identical(offsetGivenAsScalar, TRUE)) {
      offset <- rep_len(offset, numObservations)
    }
    
    ## predictors
    modelTerms <- terms(modelFrame)
    if (is.empty.model(modelTerms)) stop("predictors must be specified for regression tree analysis")
    
    termLabels <- attr(modelTerms, "term.labels")
    badLabels <- grepl("`.* .*`", termLabels)
    if (sum(badLabels) > 0)
      termLabels[badLabels] <- gsub("^`(.*)`$", "\\1", termLabels[badLabels])
    
    
    x <- makeModelMatrixFromDataFrame(modelFrame[termLabels])
    
    if (!testIsMissing) {
      testCall <- matchedCall
      testCall <- testCall[c(1L, match(c("formula", "data", "test"), names(testCall), nomatch = 0L))]
      names(testCall)[which(names(testCall) == "test")] <- "term"
      testCall[[1L]] <- quoteInNamespace(findTermInFormulaData)
      
      temp <- eval(testCall, parent.frame())
      if (!is.null(temp)) test <- temp
    }
  } else if (is.numeric(formula) || is.data.frame(formula) || is.factor(formula)) {
    ## backwards compatibility of bart(x.train, y.train, x.test)
    if (dataIsMissing || is.null(data)) data <- rep(0, NROW(formula))
    if (!is.numeric(data) && !is.data.frame(data) && !is.factor(data)) stop("when 'formula' is numeric, 'data' must be numeric as well")
    
    if (is.factor(data)) {
      y <- as.double(as.integer(data) - 1L)
    } else {
      y <- as.double(data)
    }
    if (NROW(formula) != NROW(y))
      stop("'x' must have the same number of observations as 'y'")
    initialNumObservations <- NROW(y)
    
    if (missing(subset) || is.null(subset)) subset <- seq.int(length(y))
    y <- y[subset]

    if (is.data.frame(formula)) formula <- makeModelMatrixFromDataFrame(formula)
    x <- if (!is.matrix(formula)) formula[subset] else formula[subset,,drop=FALSE]
    
    if (missing(weights)) weights <- NULL
    if (!is.null(weights)) {
      if (!is.numeric(weights)) stop("'weights' must be a numeric vector")
      weights <- rep_len(weights, initialNumObservations)[subset]
    }

    if (offsetIsMissing) offset <- NULL
    if (!is.null(offset)) {
      if (!is.numeric(offset)) stop("'offset' must be numeric")
      if (length(offset) == 1L) {
        offset <- rep_len(offset, initialNumObservations)
        offsetGivenAsScalar <- TRUE
      } else {
        offsetGivenAsScalar <- FALSE
      }
      if (length(offset) != initialNumObservations) stop("length of 'offset' must equal length of 'y'")
      originalOffset <- offset
      offset <- offset[subset]
    }
    
    completeCases <- stats::complete.cases(x, y)
    
    y <- y[completeCases]
    x <- if (!is.matrix(x)) x[completeCases] else x[completeCases,,drop=FALSE]
    if (length(attributes(formula)) > 0L) for (attributeName in names(attributes(formula))) {
      if (attributeName == "dim") next
      if (attributeName == "dimnames" && !identical(dim(formula), dim(x))) next
      attr(x, attributeName) <- attr(formula, attributeName)
    }
    if (!is.null(weights)) weights <- weights[completeCases]
    if (!is.null(offset)) offset <- offset[completeCases]
  } else {
    stop("unrecognized 'formula' type; must be coercible to numeric or a valid formula object")
  }
  
  if (is.vector(x)) x <- as.matrix(x)
  if (is.data.frame(x)) x <- makeModelMatrixFromDataFrame(x)
  
  x.test <- NULL
  if (!testIsMissing && !is.null(test))
    x.test <- validateXTest(test, attr(x, "term.labels"), ncol(x), colnames(x), attr(x, "drop"))
  
  if (!is.null(x.test)) {
    if (testOffsetIsMissing) {
      ## default is offset.test = offset
      if (identical(offsetGivenAsScalar, TRUE)) {
        offset.test <- rep_len(offset[1L], nrow(x.test))
        testUsesRegularOffset <- TRUE
      } else if (identical(offsetGivenAsScalar, FALSE)) {
        if (nrow(x.test) != length(y)) stop("vectored 'offset' cannot be directly applied to test data of unequal length")
        offset.test <- offset
        testUsesRegularOffset <- TRUE
      }
    } else {
      #environment(getTestOffset) <- sys.frame(sys.nframe())
      #testOffsetInfo <- getTestOffset()
      testOffsetInfo <- eval(getTestOffset)
      
      offset.test <- testOffsetInfo$offset.test
      testUsesRegularOffset <- testOffsetInfo$testUsesRegularOffset
      
      if (!is.null(offset.test)) offset.test <- rep_len(offset.test, nrow(x.test))
    }
  } else {
    if (testOffsetIsMissing) offset.test <- NULL
  }
  
  methods::new("dbartsData", modelMatrices = namedList(y, x, x.test, weights, offset, offset.test, testUsesRegularOffset), n.cuts = NA_integer_, sigma = NA_real_)
}

