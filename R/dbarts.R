setMethod("initialize", "dbartsControl",
          function(.Object, ...)
{
  .Object <- callNextMethod()

  validObject(.Object)
  .Object
})

## we don't actually use these defaults; see class definition
## this is only provided for UI hints. Exception is n.cuts, which
## isn't part of class
dbartsControl <-
  function(verbose = FALSE, keepTrainingFits = TRUE, useQuantiles = FALSE,
           keepTrees = FALSE, n.samples = NA_integer_, n.cuts = 100L,
           n.burn = 200L, n.trees = 75L, n.chains = 4L, n.threads = guessNumCores(),
           n.thin = 1L, printEvery = 100L, printCutoffs = 0L,
           rngKind = "default", rngNormalKind = "default", updateState = TRUE)
{
  result <- new("dbartsControl",
                verbose = as.logical(verbose),
                keepTrainingFits = as.logical(keepTrainingFits),
                useQuantiles = as.logical(useQuantiles),
                keepTrees = as.logical(keepTrees),
                n.samples = coerceOrError(n.samples, "integer"),
                n.burn = coerceOrError(n.burn, "integer"),
                n.trees = coerceOrError(n.trees, "integer"),
                n.chains = coerceOrError(n.chains, "integer"),
                n.threads = coerceOrError(n.threads, "integer"),
                n.thin = coerceOrError(n.thin, "integer"),
                printEvery = coerceOrError(printEvery, "integer"),
                printCutoffs = coerceOrError(printCutoffs, "integer"),
                rngKind = rngKind,
                rngNormalKind = rngNormalKind,
                updateState = as.logical(updateState))
  
  n.cuts <- coerceOrError(n.cuts, "integer")
  if (n.cuts <= 0L) stop("'n.cuts' must be a positive integer")
  attr(result, "n.cuts") <- n.cuts 

  result
}

validateArgumentsInEnvironment <- function(envir, func, control, verbose, n.samples, sigma)
{
  controlIsMissing <- missing(control)
  
  if (!controlIsMissing) {
    if (!inherits(control, "dbartsControl"))
      stop("'control' argument must be of class dbartsControl; use dbartsControl() function to create")
    envir$control <- control
  }

  if (!missing(verbose)) {
    if (!is.logical(verbose) || is.na(verbose)) stop("'verbose' argument to dbarts must be TRUE/FALSE")
  } else if (!controlIsMissing) {
    envir$verbose <- control@verbose
  }
  
  if (!missing(n.samples)) {
    tryCatch(n.samples <- as.integer(n.samples), warning = function(e)
             stop("'n.samples' argument to dbarts must be coerceable to integer type"))
    if (length(n.samples) != 1L) stop("'n.samples' must be of length 1")
    if (is.null(n.samples))
      stop("'n.samples' argument to dbarts cannot be NULL")
    if (is.na(n.samples) || n.samples < 0L)
      stop("'n.samples' argument to dbarts must be a non-negative integer")
    envir$control@n.samples <- n.samples
  } else if (controlIsMissing || is.na(control@n.samples)) {
    envir$control@n.samples <- formals(func)[["n.samples"]]
  }

  if (!missing(sigma) && !is.na(sigma)) {
    tryCatch(sigma <- as.double(sigma), warning = function(e)
             stop("'sigma' argument to dbarts must be coerceable to numeric type"))
    if (length(sigma) != 1L) stop("'sigma' must be of length 1")
    if (is.null(sigma) || sigma <= 0.0) stop("'sigma' argument to dbarts must be positive")
    
    envir$sigma <- sigma
  }
}

dbarts <- function(formula, data, test, subset, weights, offset, offset.test = offset,
                   verbose = FALSE, n.samples = 800L,
                   tree.prior = cgm, node.prior = normal, resid.prior = chisq,
                   control = dbartsControl(), sigma = NA_real_)
{
  matchedCall <- match.call()
  
  evalEnv <- parent.frame(1L)

  validateCall <- redirectCall(matchedCall, quoteInNamespace(validateArgumentsInEnvironment))
  validateCall <- addCallArgument(validateCall, 1L, sys.frame(sys.nframe()))
  validateCall <- addCallArgument(validateCall, 2L, dbarts::dbarts)
  eval(validateCall, evalEnv, getNamespace("dbarts"))

  if (length(control@call) == 1L && control@call == call("NA")) control@call <- matchedCall
  control@verbose <- verbose

  dataCall <- redirectCall(matchedCall, quoteInNamespace(dbartsData))
  data <- eval(dataCall, evalEnv)
  #cat("x address after dbartsData call: ", .Call("dbarts_getPointerAddress", data@x), "\n", sep = "")
  
  data@n.cuts <- rep_len(attr(control, "n.cuts"), ncol(data@x))
  data@sigma  <- sigma
  attr(control, "n.cuts") <- NULL
  
  
  uniqueResponses <- unique(data@y)
  if (length(uniqueResponses) == 2 && all(sort(uniqueResponses) == c(0, 1))) control@binary <- TRUE
  
  if (is.na(data@sigma) && !control@binary)
    data@sigma <- summary(lm(data@y ~ data@x, weights = data@weights, offset = data@offset))$sigma
  
  ## bart will passthrough with offset == something no matter what, which we can NULL out
  if (!control@binary && !is.null(data@offset) && all(data@offset == 0.0)) {
    data@offset <- NULL
  }
  if (!control@binary && !is.null(data@offset.test) && all(data@offset.test == 0.0)) {
    data@offset.test <- NULL
  }
  #cat("x address after updating data: ", .Call("dbarts_getPointerAddress", data@x), "\n", sep = "")

  parsePriorsCall <- redirectCall(matchedCall, quoteInNamespace(parsePriors))
  parsePriorsCall <- setDefaultsFromFormals(parsePriorsCall, formals(dbarts), "tree.prior", "node.prior", "resid.prior")
  parsePriorsCall$control <- control
  parsePriorsCall$data <- data
  parsePriorsCall$parentEnv <- evalEnv
  priors <- eval(parsePriorsCall)

  model <- new("dbartsModel", priors$tree.prior, priors$node.prior, priors$resid.prior)

  result <- new("dbartsSampler", control, model, data)
  #cat("x address after creating sampler: ", .Call("dbarts_getPointerAddress", data@x), "\n", sep = "")
  #cat("x address in sampler$data: ", .Call("dbarts_getPointerAddress", result$data@x), "\n", sep = "")
  result
}


setClassUnion("listOrNULL", c("list", "NULL"))

dbartsSampler <-
  setRefClass("dbartsSampler",
              fields = list(
                pointer = "externalptr",
                control = "dbartsControl",
                model   = "dbartsModel",
                data    = "dbartsData",
                state   = "listOrNULL"
                ),
              methods = list(
                initialize =
                  function(control, model, data, ...)
                {
                  if (!inherits(control, "dbartsControl")) stop("'control' must inherit from dbartsControl")
                  if (!inherits(model, "dbartsModel")) stop("'model' must inherit from dbartsModel")
                  if (!inherits(data, "dbartsData")) stop("'data' must inherit from dbartsData")
                  
                  .self$control <- control
                  .self$model   <- model
                  .self$data    <- data
                  
                  .self$state   <- NULL
                  .self$pointer <- .Call(C_dbarts_create, .self$control, .self$model, .self$data)
                  
                  callSuper(...)
                },
                run = function(numBurnIn, numSamples, updateState = NA) {
                  'Runs the posterior sampler for numBurnIn + numSamples iterations and
                   returns a list with the results.'
                  if (missing(numBurnIn))  numBurnIn  <- NA_integer_
                  if (missing(numSamples)) numSamples <- NA_integer_
                  
                  ptr <- getPointer()
                  samples <- .Call(C_dbarts_run, ptr, as.integer(numBurnIn), as.integer(numSamples))

                  if ((is.na(updateState) && control@updateState == TRUE) || identical(updateState, TRUE))
                    storeState(ptr)

                  if (is.null(samples)) return(invisible(NULL))
                  
                  samples
                },
                sampleTreesFromPrior = function(updateState = NA) {
                  'Draws tree structure from prior; does not update tree predictions, so sampler
                   will be in invalid state'
                  
                  ptr <- getPointer()
                  samples <- .Call(C_dbarts_sampleTreesFromPrior, ptr)

                  if ((is.na(updateState) && control@updateState == TRUE) || identical(updateState, TRUE))
                    storeState(ptr)

                  invisible(NULL)
                },
                copy = function(shallow = FALSE) {
                  'Creates a deep or shallow copy of the sampler.'
                  dupe <-
                    if (shallow) {
                      dbartsSampler$new(control, model, data)
                    } else {
                      newData <- data
                      ## only need to dupe things that can be changed internally, as the rest will be
                      ## simply swapped out
                      newData@x <- .Call(C_dbarts_deepCopy, data@x)
                      if (!is.null(data@x.test)) {
                        newData@x.test <- .Call(C_dbarts_deepCopy, data@x.test)
                      }
                      dbartsSampler$new(control, model, newData)
                    }
                  
                  if (!is.null(state)) {
                    newState <- state
                    for (chainNum in seq_len(control@n.chains)) {
                      newState[[chainNum]]@fit.tree <- .Call(C_dbarts_deepCopy, state[[chainNum]]@fit.tree)
                      newState[[chainNum]]@fit.total <- .Call(C_dbarts_deepCopy, state[[chainNum]]@fit.total)
                      if (!is.null(data@x.test)) {
                        newState[[chainNum]]@fit.test <- .Call(C_dbarts_deepCopy, state[[chainNum]]@fit.test)
                      }
                      newState[[chainNum]]@sigma <- .Call(C_dbarts_deepCopy, state[[chainNum]]@sigma)
                      newState[[chainNum]]@trees <- .Call(C_dbarts_deepCopy, state[[chainNum]]@trees)
                    }
                    attr(newState, "runningTime") <- .Call(C_dbarts_deepCopy, attr(state, "runningTime"))
                    
                    dupe$setState(newState)
                  }
                  
                  dupe
                },
                show = function() {
                  'Pretty prints the object.'
                  
                  cat("dbarts sampler\n")
                  cat("  call: ")
                  writeLines(deparse(control@call))
                  cat("\n")
                  
                  invisible(NULL)
                },
                predict = function(x.test, offset.test) {
                  'Using existing sampler to predict for new data without re-running.'
                  
                  selfEnv <- parent.env(environment())
                  
                  if (!control@keepTrees) stop("predict requires that keepTrees is TRUE")
                   
                  ptr <- getPointer()
                  
                  x.test <- validateXTest(x.test, attr(data@x, "term.labels"), ncol(data@x), colnames(data@x), attr(data@x, "drop"))
                  if (is.null(x.test)) stop("x.test cannot be NULL")
                  
                  if (missing(offset.test) || is.null(offset.test)) {
                    offset.test <- NA_real_
                  } else {
                    offset.test <- as.double(offset.test)
                    if (length(offset.test) == 1)
                      offset.test <- rep_len(offset, nrow(x.test))
                    if (!identical(length(offset.test), nrow(x.test)))
                      stop("length of test offset must be equal to number of rows in test matrix")
                  }
                  
                  .Call(C_dbarts_predict, ptr, x.test, offset.test)
                },
                setControl = function(newControl) {
                  'Sets the control object for the sampler to a new one. Preserves the call() slot.'
                  
                  if (!inherits(newControl, "dbartsControl")) stop("'control' must inherit from dbartsControl")
                  
                  selfEnv <- parent.env(environment())
                  
                  newControl@binary <- control@binary
                  newControl@call   <- control@call
                  
                  ptr <- getPointer()
                  
                  selfEnv$control <- newControl
                  .Call(C_dbarts_setControl, ptr, control)
                  
                  invisible(NULL)
                },
                setModel = function(newModel) {
                  'Sets the model object for the sampler to a new one.'
                  
                  if (!inherits(newModel, "dbartsModel")) stop("'model' must inherit from dbartsModel")
                  ptr <- getPointer()
                  selfEnv <- parent.env(environment())
                  
                  selfEnv$model <- newModel
                  .Call(C_dbarts_setModel, ptr, model)
                  
                  invisible(NULL)
                },
                setData = function(newData, updateState = NA) {
                  'Sets the data object for the sampler to a new one. Preserves the n.cuts and sigma slots.'
                  
                  if (!inherits(newData, "dbartsData")) stop("'data' must inherit from dbartsData")
                  
                  newData@n.cuts <- data@n.cuts
                  newData@sigma  <- data@sigma
                  
                  ptr <- getPointer()
                  selfEnv <- parent.env(environment())
                  
                  selfEnv$data <- newData
                  .Call(C_dbarts_setData, ptr, data)
                  
                  if ((is.na(updateState) && control@updateState == TRUE) || identical(updateState, TRUE))
                    storeState(ptr)
                  
                  invisible(NULL)
                },
                setResponse = function(y, updateState = NA) {
                  'Changes the response against which the sampler is fitted.'
                  ptr <- getPointer()
                  selfEnv <- parent.env(environment())
                  
                  selfEnv$data@y <- as.double(y)
                  .Call(C_dbarts_setResponse, ptr, data@y)

                  if ((is.na(updateState) && control@updateState == TRUE) || identical(updateState, TRUE))
                    storeState(ptr)

                  invisible(NULL)
                },
                setOffset = function(offset, updateState = NA) {
                  'Changes the offset slot used to adjust the response.'
                  ptr <- getPointer()
                  selfEnv <- parent.env(environment())
                  
                  offset.test <- NA
                  if (is.null(offset)) {
                    if (identical(data@testUsesRegularOffset, TRUE)) offset.test <- NULL
                  } else {
                    offset <- as.double(offset)

                    if (length(offset) == 1) {
                      if (identical(data@testUsesRegularOffset, TRUE)) {
                        if (!is.null(data@x.test)) {
                          offset.test <- rep_len(offset, nrow(data@x.test))
                        } else {
                          offset.test <- NULL
                        }
                      }
                      offset <- rep_len(offset, length(data@y))
                    } else {
                      if (!identical(length(offset), length(data@y)))
                        stop('length of replacement offset is not equal to number of observations')
                      if (identical(data@testUsesRegularOffset, TRUE)) {
                        if (!is.null(data@x.test) && length(offset) == nrow(data@x.test)) {
                          offset.test <- offset
                        } else {
                          offset.test <- NULL
                        }
                      }
                    }
                  }

                  selfEnv$data@offset <- offset
                  .Call(C_dbarts_setOffset, ptr, data@offset)
                  if (!identical(offset.test, NA)) {
                    selfEnv$data@offset.test <- offset.test
                    .Call(C_dbarts_setTestOffset, ptr, data@offset.test)
                  }

                  if ((is.na(updateState) && control@updateState == TRUE) || identical(updateState, TRUE))
                    storeState(ptr)
                  
                  invisible(NULL)
                },
                setPredictor = function(x, column, updateState = NA) {
                  'Changes a single column of the predictor matrix, or the entire matrix itself if the column argument is missing. TRUE/FALSE returned as to whether or not the operation was successful.'
                  
                  selfEnv <- parent.env(environment())
                  
                  if (control@keepTrees && (is.null(selfEnv$keepTreesWarnOnce) || self$keepTreesWarnOnce == FALSE)) {
                    warning("changing predictor with keepTrees == TRUE can render old predictions invalid")
                    selfEnv$keepTreesWarnOnce <- TRUE
                  }

                  columnIsMissing <- missing(column)

                  if (!columnIsMissing && is.character(column)) {
                    if (is.null(colnames(data@x))) stop("column names not specified at initialization, so cannot be replaced by name")
                    
                    column <- match(column, colnames(data@x))
                    if (is.na(column)) stop("column name not found in names of current X")
                  }
                  
                  x <- if (is.matrix(x)) matrix(as.double(x), nrow(x)) else as.double(x)
                  
                  ptr <- getPointer()
                  
                  updateSuccessful <-
                    if (columnIsMissing) {
                      selfEnv$data@x <- x
                      .Call(C_dbarts_setPredictor, ptr, data@x)
                    } else {
                      .Call(C_dbarts_updatePredictor, ptr, x, as.integer(column))
                    }
                  
                  if ((is.na(updateState) && control@updateState == TRUE) || identical(updateState, TRUE))
                    storeState(ptr)

                  updateSuccessful
                },
                setTestPredictor = function(x.test, column, updateState = NA) {
                  'Changes a single column of the test predictor matrix.'

                  columnIsMissing <- missing(column)
                  
                  if (!columnIsMissing && is.character(column)) {
                    if (is.null(colnames(data@x.test))) stop("column names not specified at initialization, so cannot be replaced by name")
                    
                    column <- match(column, colnames(data@x.test))
                    if (is.na(column)) stop("column name not found in names of current test predictor matrix")
                  }
                  
                  ptr <- getPointer()
                  selfEnv <- parent.env(environment())
                  
                  if (columnIsMissing) {
                    selfEnv$data@x.test <- validateXTest(x.test, attr(data@x, "term.labels"), ncol(data@x), colnames(data@x), attr(data@x, "drop"))
                    .Call(C_dbarts_setTestPredictor, ptr, data@x.test)
                  } else {
                    x.test <- if (is.matrix(x.test)) matrix(as.double(x.test), nrow(x.test)) else as.double(x.test)
                    .Call(C_dbarts_updateTestPredictor, ptr, x.test, as.integer(column))
                  }
                  
                  if ((is.na(updateState) && control@updateState == TRUE) || identical(updateState, TRUE))
                    storeState(ptr)

                  invisible(NULL)
                },
                setTestPredictorAndOffset = function(x.test, offset.test, updateState = NA) {
                  'Changes the test predictor matrix, and optionally the test offset.'
                  ptr <- getPointer()
                  selfEnv <- parent.env(environment())
                  
                  x.test <- validateXTest(x.test, attr(data@x, "term.labels"), ncol(data@x), colnames(data@x), attr(data@x, "drop"))
                  
                  if (!missing(offset.test)) {
                    if (is.null(x.test)) {
                      if (!is.null(offset.test)) stop("when test matrix is NULL, test offset must be as well")
                    } else {
                      if (!is.null(offset.test)) {
                        offset.test <- as.double(offset.test)
                        if (length(offset.test) == 1) {
                          offset.test <- rep_len(offset, nrow(x.test))
                        }
                        if (!identical(length(offset.test), nrow(x.test))) {
                          stop("length of test offset must be equal to number of rows in test matrix")
                        }
                      }
                    }
                    selfEnv$data@testUsesRegularOffset <- FALSE
                    selfEnv$data@x.test <- x.test
                    selfEnv$data@offset.test <- offset.test
                    .Call(C_dbarts_setTestPredictorAndOffset, ptr, data@x.test, data@offset.test)
                  } else {
                    selfEnv$data@x.test <- x.test
                    .Call(C_dbarts_setTestPredictorAndOffset, ptr, data@x.test, NA_real_)
                  }
                  if ((is.na(updateState) && control@updateState == TRUE) || identical(updateState, TRUE))
                    storeState(ptr)
                  invisible(NULL)
                },
                setTestOffset = function(offset.test, updateState = NA) {
                  'Changes the test offset.'
                  ptr <- getPointer()
                  selfEnv <- parent.env(environment())

                  selfEnv$data@testUsesRegularOffset <- FALSE
                  if (!is.null(offset.test)) {
                    if (is.null(data@x.test)) stop("when test matrix is NULL, test offset must be as well")
                    if (length(offset.test) == 1) offset.test <- rep_len(offset.test, nrow(data@x.test))
                    if (length(offset.test) != nrow(data@x.test))
                      stop("length of test offset must be equal to number of rows in test matrix")
                  }
                  selfEnv$data@offset.test <- offset.test
                  .Call(C_dbarts_setTestOffset, ptr, data@offset.test)
                  
                  if ((is.na(updateState) && control@updateState == TRUE) || identical(updateState, TRUE))
                    storeState(ptr)

                  invisible(NULL)
                },
                getPointer = function() {
                  'Returns the underlying reference pointer, checking for consistency first.'
                  selfEnv <- parent.env(environment())
                  
                  if (.Call(C_dbarts_isValidPointer, pointer) == FALSE) {
                    oldVerbose <- control@verbose
                    selfEnv$control@verbose <- FALSE
                    selfEnv$pointer <- .Call(C_dbarts_create, control, model, data)
                    if (!is.null(state)) {
                      selfEnv$state <- .Call(C_dbarts_restoreState, pointer, state)
                    }
                    if (control@verbose != oldVerbose) {
                      selfEnv$control@verbose <- oldVerbose
                      .Call(C_dbarts_setControl, pointer, control)
                    }
                  }
                  
                  pointer
                },
                setState = function(newState) {
                  'Sets the internal state from a cache.'
                  
                  if (!is.list(newState)) stop("'state' must be a list of dbartsState objects")
                  if (length(newState) != control@n.chains) stop("'state' length must equal number of chains")
                  for (chainNum in seq_along(newState))
                    if (!is(newState[[chainNum]], "dbartsState")) stop("'state' must inherit from dbartsState")
                  
                  selfEnv <- parent.env(environment())
                  if (.Call(C_dbarts_isValidPointer, pointer) == FALSE) {
                    oldVerbose <- control@verbose
                    selfEnv$control@verbose <- FALSE
                    selfEnv$pointer <- .Call(C_dbarts_create, control, model, data)
                    if (control@verbose != oldVerbose) {
                      selfEnv$control@verbose <- oldVerbose
                      .Call(C_dbarts_setControl, pointer, control)
                    }
                  }
                  if (!is.null(newState)) {
                    selfEnv$state <- .Call(C_dbarts_restoreState, pointer, newState)
                  }

                  invisible(NULL)
                },
                storeState = function(ptr = getPointer()) {
                  'Updates the cached internal state used for saving/loading.'
                  selfEnv <- parent.env(environment())
                  if (is.null(state)) {
                    selfEnv$state <- .Call(C_dbarts_createState, ptr)
                  } else {
                    .Call(C_dbarts_storeState, ptr, state)
                  }

                  invisible(NULL)
                },
                printTrees = function(chainNums, sampleNums, treeNums) {
                  'Produces an info dump of the internal state of the trees.'
                  matchedCall <- match.call()
                  if (is.null(matchedCall$chainNums)) chainNums <- seq_len(control@n.chains)
                  if (is.null(matchedCall$sampleNums)) {
                    sampleNums <- if (control@keepTrees) seq_len(control@n.samples) else NULL
                  } else {
                    if (!control@keepTrees) {
                      warning("sampleNums ignored if keepTrees is FALSE")
                      sampleNums <- NULL
                    } else {
                      sampleNums <- as.integer(sampleNums)
                    }
                  }
                  if (is.null(matchedCall$treeNums)) treeNums <- seq_len(control@n.trees)
                  
                  ptr <- getPointer()
                  invisible(.Call(C_dbarts_printTrees, ptr, as.integer(chainNums), sampleNums, as.integer(treeNums)))
                },
                plotTree = function(treeNum, treePlotPars = list(nodeHeight = 12, nodeWidth = 40, nodeGap = 8), ...) {
                  'Minimialist visualization of tree branching and contents.'
                  
                  cutPoints <- createCutPoints(.self)
  
                  tree <- buildTree(strsplit(gsub("\\.", "\\. ", state@trees[treeNum]), " ", fixed = TRUE)[[1]])
                  tree$remainder <- NULL
                  
                  tree$indices <- seq_len(nrow(data@x))
                  tree <- fillObservationsForNode(tree, .self, cutPoints)
                  
                  tree <- fillPlotInfoForNode(tree, .self, state@fit.tree[,treeNum])
                  
                  maxDepth <- getMaxDepth(tree)
                  
                  tree <- fillPlotCoordinatesForNode(tree, maxDepth, 1L, 1L)
                  numEndNodes <- tree$index - 1L
                  tree$index <- NULL
                  
                  plotHeight <- treePlotPars$nodeHeight * maxDepth + treePlotPars$nodeGap * (maxDepth - 1)
                  dotsList <- list(...)
                  dotsList$mar = c(0, 0, 0, 0)
                  par(dotsList)
                  plot(NULL, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "",
                       xlim = c(0, treePlotPars$nodeWidth * numEndNodes), ylim = c(0, plotHeight))
                  plotNode(tree, .self, cutPoints, treePlotPars)
                  
                  invisible(NULL)
                }
                #crossvalidate = function(K = 5L, n.reps = 200L, n.burn = c(control@n.burn, control@n.burn %/% 5L),
                #                         loss = if (control@binary) "mcr" else "rmse", n.threads = guessNumCores(),
                #                         n.trees, k, power, base, drop = TRUE) {
                #  'K-fold crossvalidates across sets of parameters.'
                #  
                #  if (missing(n.trees)) n.trees <- control@n.trees
                #  if (missing(k))       k <- model@node.prior@k
                #  if (missing(power))   power <- model@tree.prior@power
                #  if (missing(base))    base <- model@tree.prior@base
                #  
                #  if (is.function(loss)) loss <- list(loss, parent.frame(1L))
                #  
                #  kOrder <- order(k, decreasing = TRUE)
                #  kOrder.inv <- kOrder; kOrder.inv[kOrder] <- seq_along(kOrder)
                #  k <- k[kOrder] 
                #  
                #  n.trees <- as.integer(n.trees)
                #  k       <- as.double(k)
                #  power   <- as.double(power)
                #  base    <- as.double(base)
                #  drop    <- as.logical(drop)
                #  
                #  ptr <- getPointer()
                #  result <- .Call(C_dbarts_xbart, ptr, as.integer(K), as.integer(n.reps), as.integer(n.burn), loss, as.integer(n.threads),
                #                  n.trees, k, power, base, drop)
                #  
                #  if (is.null(result) || is.null(dim(result))) return(result)
                #  
                #  varNames <- c("n.trees", "k", "power", "base")
                #  if (identical(drop, TRUE))
                #    varNames <- varNames[sapply(varNames, function(varName) if (length(get(varName)) > 1L) TRUE else FALSE)]
                #  
                #  if ("k" %in% varNames && any(kOrder != seq_along(k))) {
                #    indices <- rep(list(bquote()), length(dim(result)))
                #    indices[[1L + which(varNames == "k")]] <- kOrder.inv
                #    indexCall <- as.call(c(list(as.name("["), quote(result)), indices))
                #    result <- eval(indexCall)
                #    k <- k[kOrder.inv]
                #  }
                #  
                #  dimNames <- vector("list", length(dim(result)))
                #  for (i in seq_along(varNames)) {
                #    x <- get(varNames[i])
                #    dimNames[[i + 1L]] <- as.character(if (is.double(x)) signif(x, 2L) else x)
                #  }
                #  names(dimNames) <- c("rep", varNames)
                #  dimnames(result) <- dimNames
                #  
                #  result
                #}
              ))
