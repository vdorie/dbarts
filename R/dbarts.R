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
           n.samples = NA_integer_, n.cuts = 100L,
           n.burn = 100L, n.trees = 200L, n.threads = 1L,
           n.thin = 1L, printEvery = 100L, printCutoffs = 0L, updateState = TRUE)
{
  result <- new("dbartsControl",
                verbose = as.logical(verbose),
                keepTrainingFits = as.logical(keepTrainingFits),
                useQuantiles = as.logical(useQuantiles),
                n.samples = coerceOrError(n.samples, "integer"),
                n.burn = coerceOrError(n.burn, "integer"),
                n.trees = coerceOrError(n.trees, "integer"),
                n.threads = coerceOrError(n.threads, "integer"),
                n.thin = coerceOrError(n.thin, "integer"),
                printEvery = coerceOrError(printEvery, "integer"),
                printCutoffs = coerceOrError(printCutoffs, "integer"),
                updateState = as.logical(updateState))
  
  n.cuts <- coerceOrError(n.cuts, "integer")
  if (n.cuts <= 0L) stop("'n.cuts' must be a positive integer")
  attr(result, "n.cuts") <- n.cuts 

  result
}

validateArgumentsInEnvironment <- function(envir, control, verbose, n.samples, sigma)
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
    if (n.samples < 0L) return("'n.samples' argument to dbarts must be a non-negative integer")
    envir$control@n.samples <- n.samples
  } else if (controlIsMissing || is.na(control@n.samples)) {
    envir$control@n.samples <- formals(dbarts)[["n.samples"]]
  }

  if (!missing(sigma) && !is.na(sigma)) {
    tryCatch(sigma <- as.double(sigma), warning = function(e)
             stop("'sigma' argument to dbarts must be coerceable to numeric type"))
    if (sigma <= 0.0) stop("'sigma' argument to dbarts must be positive")
    
    envir$sigma <- sigma
  }
}

dbarts <- function(formula, data, test, subset, weights, offset, offset.test = offset,
                   verbose = FALSE, n.samples = 1000L,
                   tree.prior = cgm, node.prior = normal, resid.prior = chisq,
                   control = dbartsControl(), sigma = NA_real_)
{
  matchedCall <- match.call()

  validateCall <- prepareCallWithArguments(matchedCall, quoteInNamespace(validateArgumentsInEnvironment), "control", "verbose", "n.samples", "sigma")
  validateCall <- addCallArgument(validateCall, 1L, sys.frame(sys.nframe()))
  eval(validateCall, parent.frame(1L), asNamespace("dbarts"))

  if (control@call[[1]] != call("NA")[[1]]) control@call <- matchedCall
  control@verbose <- verbose

  parseDataCall <- prepareCallWithArguments(matchedCall, quoteInNamespace(parseData), "formula", "data", "test", "subset", "weights", "offset", "offset.test")
  modelMatrices <- eval(parseDataCall, parent.frame(1L))
  
  data <- new("dbartsData", modelMatrices, attr(control, "n.cuts"), sigma)
  attr(control, "n.cuts") <- NULL
  if (is.na(data@sigma))
    data@sigma <- summary(lm(data@y ~ data@x, weights = data@weights))$sigma

  uniqueResponses <- unique(data@y)
  if (length(uniqueResponses) == 2 && all(sort(uniqueResponses) == c(0, 1))) control@binary <- TRUE
  ## bart will passthrough with offset == something no matter what, which we can NULL out
  if (!control@binary && !is.null(data@offset) && all(data@offset == 0.0)) {
    data@offset <- NULL
  }
  if (!control@binary && !is.null(data@offset.test) && all(data@offset.test == 0.0)) {
    data@offset.test <- NULL
  }

  parsePriorsCall <- prepareCallWithArguments(matchedCall, quoteInNamespace(parsePriors), "tree.prior", "node.prior", "resid.prior")
  parsePriorsCall <- setDefaultsFromFormals(parsePriorsCall, formals(dbarts), "tree.prior", "node.prior", "resid.prior")
  parsePriorsCall$control <- control
  parsePriorsCall$data <- data
  parsePriorsCall$parentEnv <- parent.frame(1L)
  priors <- eval(parsePriorsCall)

  model <- new("dbartsModel", priors$tree.prior, priors$node.prior, priors$resid.prior)

  new("dbartsSampler", control, model, data)
}

setClassUnion("dbartsStateOrNULL", c("dbartsState", "NULL"))

dbartsSampler <-
  setRefClass("dbartsSampler",
              fields = list(
                pointer = "externalptr",
                control = "dbartsControl",
                model   = "dbartsModel",
                data    = "dbartsData",
                state   = "dbartsStateOrNULL"
                ),
              methods = list(
                initialize =
                  function(control, model, data, ...)
                {
                  if (!inherits(control, "dbartsControl")) stop("'control' must inherit from dbartsControl")
                  if (!inherits(model, "dbartsModel")) stop("'model' must inherit from dbartsModel")
                  if (!inherits(data, "dbartsData")) stop("'data' must inherit from dbartsData")

                  control <<- control
                  model   <<- model
                  data    <<- data
                  state   <<- NULL
                  pointer <<- .Call("dbarts_create", control, model, data)
                  
                  callSuper(...)
                },
                run = function(numBurnIn, numSamples, updateState = NA) {
                  'Runs the posterior sampler for numBurnIn + numSamples iterations and
                   returns a list with the results.'
                  if (missing(numBurnIn))  numBurnIn  <- NA_integer_
                  if (missing(numSamples)) numSamples <- NA_integer_
                  
                  ptr <- getPointer()
                  samples <- .Call("dbarts_run", ptr, as.integer(numBurnIn), as.integer(numSamples))

                  if ((is.na(updateState) && control@updateState == TRUE) || identical(updateState, TRUE))
                    storeState(ptr)

                  if (is.null(samples)) return(invisible(NULL))
                  
                  samples
                },
                copy = function(shallow = FALSE)
                {
                  dupe <-
                    if (shallow) {
                      dbartsSampler$new(control, model, data)
                    } else {
                      newData <- data
                      ## only need to dupe things that can be changed internally, as the rest will be
                      ## simply swapped out
                      newData@x <- .Call("dbarts_deepCopy", data@x)
                      if (!is.null(data@x.test)) {
                        newData@x.test <- .Call("dbarts_deepCopy", data@x.test)
                      }
                      dupe <- dbartsSampler$new(control, model, newData)
                    }

                  if (!is.null(state)) {
                    newState <- state
                    newState@fit.tree <- .Call("dbarts_deepCopy", state@fit.tree)
                    newState@fit.total <- .Call("dbarts_deepCopy", state@fit.total)
                    if (!is.null(data@x.test)) {
                      newState@fit.test <- .Call("dbarts_deepCopy", state@fit.test)
                    }
                    newState@sigma <- .Call("dbarts_deepCopy", state@sigma)
                    newState@runningTime <- .Call("dbarts_deepCopy", state@runningTime)
                    newState@trees <- .Call("dbarts_deepCopy", state@trees)

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
                  
                  ##xDisplay <- NULL
                  ##stringConnection <- textConnection("xDisplay", "w", local = TRUE)
                  ##sink(stringConnection)
                    
                  ##cat("  x:")
                  ##if (!is.null(colnames(data@x))) {
                  ##  for (i in 1:min(ncol(data@x), 10)) cat(sprintf(" %8.8s", colnames(data@x)[i]))
                  ##  if (ncol(data@x) > 10) cat(" ...")
                  ##  cat("\n")
                  ##  for (i in 1:min(nrow(data@x), 3)) {
                  ##    cat("    ")
                  ##    for (j in 1:min(ncol(data@x), 10)) cat(sprintf("     %4.4s", data@x[i,j]))
                  ##    cat("\n")
                  ##  }
                  ##}
                  ##sink()
                  ##close(stringConnection)
                    
                  ##writeLines(xDisplay)
                  
                },
                setControl = function(newControl) {
                  'Sets the control object for the sampler to a new one. Preserves the call() slot.'
                  
                  if (!inherits(newControl, "dbartsControl")) stop("'control' must inherit from dbartsControl.")

                  newControl@binary <- control@binary
                  newControl@call   <- control@call
                  
                  ptr <- getPointer()
                  control <<- newControl
                  .Call("dbarts_setControl", ptr, control)
                  
                  invisible(NULL)
                },
                setResponse = function(y, updateState = NA) {
                  'Changes the response against which the sampler is fitted.'
                  ptr <- getPointer()
                  
                  data@y <<- as.double(y)
                  .Call("dbarts_setResponse", ptr, data@y)

                  if ((is.na(updateState) && control@updateState == TRUE) || identical(updateState, TRUE))
                    storeState(ptr)

                  invisible(NULL)
                },
                setOffset = function(offset, updateState = NA) {
                  'Changes the offset slot used to adjust the response.'
                  ptr <- getPointer()

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

                  data@offset <<- offset
                  .Call("dbarts_setOffset", ptr, data@offset)
                  if (!identical(offset.test, NA)) {
                    data@offset.test <<- offset.test
                    .Call("dbarts_setTestOffset", ptr, data@offset.test)
                  }

                  if ((is.na(updateState) && control@updateState == TRUE) || identical(updateState, TRUE))
                    storeState(ptr)
                  
                  invisible(NULL)
                },
                setPredictor = function(x, column, updateState = NA) {
                  'Changes a single column of the predictor matrix, or the entire matrix itself if the column argument is missing. TRUE/FALSE returned as to whether or not the operation was successful.'

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
                      data@x <<- x
                      .Call("dbarts_setPredictor", ptr, x)
                    } else {
                      .Call("dbarts_updatePredictor", ptr, x, as.integer(column))
                    }
                  
                  if ((is.na(updateState) && control@updateState == TRUE) || identical(updateState, TRUE))
                    storeState(ptr)

                  return(updateSuccessful)
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
                  if (columnIsMissing) {
                    data@x.test <<- validateXTest(x.test, attr(data@x, "term.labels"), ncol(data@x), colnames(data@x), attr(data@x, "drop"))
                    .Call("dbarts_setTestPredictor", ptr, data@x.test)
                  } else {
                    x.test <- if (is.matrix(x.test)) matrix(as.double(x.test), nrow(x.test)) else as.double(x.test)
                    .Call("dbarts_updateTestPredictor", ptr, x.test, as.integer(column))
                  }
                  
                  if ((is.na(updateState) && control@updateState == TRUE) || identical(updateState, TRUE))
                    storeState(ptr)

                  invisible(NULL)
                },
                setTestPredictorAndOffset = function(x.test, offset.test, updateState = NA) {
                   'Changes the test predictor matrix, and optionally the test offset.'
                   ptr <- getPointer()

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
                     data@testUsesRegularOffset <<- FALSE

                     data@x.test <<- x.test
                     data@offset.test <<- offset.test
                     .Call("dbarts_setTestPredictorAndOffset", ptr, data@x.test, data@offset.test)
                   } else {
                     data@x.test <<- x.test
                     .Call("dbarts_setTestPredictorAndOffset", ptr, data@x.test, NA_real_)
                   }

                   if ((is.na(updateState) && control@updateState == TRUE) || identical(updateState, TRUE))
                     storeState(ptr)

                   invisible(NULL)
                },
                setTestOffset = function(offset.test, updateState = NA) {
                  'Changes the test offset.'
                  ptr <- getPointer()

                  data@testUsesRegularOffset <<- FALSE
                  if (!is.null(offset.test)) {
                    if (is.null(data@x.test)) stop("when test matrix is NULL, test offset must be as well")
                    if (length(offset.test) == 1) offset.test <- rep_len(offset.test, nrow(data@x.test))
                    if (length(offset.test) != nrow(data@x.test))
                      stop("length of test offset must be equal to number of rows in test matrix")
                  }
                  data@offset.test <<- offset.test
                  .Call("dbarts_setTestOffset", ptr, data@offset.test)
                  
                  if ((is.na(updateState) && control@updateState == TRUE) || identical(updateState, TRUE))
                    storeState(ptr)

                  invisible(NULL)
                },
                getPointer = function() {
                  'Returns the underlying reference pointer, checking for consistency first.'
                  if (.Call("dbarts_isValidPointer", pointer) == FALSE) {
                    oldVerbose <- control@verbose
                    control@verbose <<- FALSE
                    pointer <<- .Call("dbarts_create", control, model, data)
                    if (!is.null(state)) {
                      state <<- .Call("dbarts_restoreState", pointer, state)
                    }
                    if (control@verbose != oldVerbose) {
                      control@verbose <<- oldVerbose
                      .Call("dbarts_setControl", pointer, control)
                    }
                  }
                  
                  pointer
                },
                setState = function(newState) {
                  'Sets the internal state from a cache.'
                  if (!is(newState, "dbartsState")) stop("'state' must inherit from dbartsState")
                  
                  if (.Call("dbarts_isValidPointer", pointer) == FALSE) {
                    oldVerbose <- control@verbose
                    control@verbose <<- FALSE
                    pointer <<- .Call("dbarts_create", control, model, data)
                    if (control@verbose != oldVerbose) {
                      control@verbose <<- oldVerbose
                      .Call("dbarts_setControl", pointer, control)
                    }
                  }
                  if (!is.null(newState)) {
                    state <<- .Call("dbarts_restoreState", pointer, newState)
                  }

                  invisible(NULL)
                },
                storeState = function(ptr = getPointer()) {
                  'Updates the cached internal state used for saving/loading.'
                  if (is.null(state)) {
                    state <<- .Call("dbarts_createState", ptr)
                  } else {
                    .Call("dbarts_storeState", ptr, state)
                  }

                  invisible(NULL)
                })
              )
