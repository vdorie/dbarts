setMethod("initialize", "dbartsControl",
          function(.Object, ...)
{
  .Object <- callNextMethod()
  ## keep around but sorta hide
##  attr(.Object, "n.cuts") <- as.integer(n.cuts)

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
  matchedCall <- match.call()
  
  argsToKeep <- names(formals(dbartsControl))
  argsToKeep <- argsToKeep[argsToKeep != "n.cuts"]
  
  newCall <- call("new", "dbartsControl")

  newRange <- seq_len(length(matchedCall) - 1L) + 2L
  oldRange <- 1L + seq_len(length(matchedCall) - 1)
  
  newCall[newRange] <- matchedCall[oldRange]
  names(newCall)[newRange] <- names(matchedCall)[oldRange]
  
  result <- eval(newCall, parent.frame())
  attr(result, "n.cuts") <- as.integer(n.cuts)

  if (attr(result, "n.cuts") <= 0L) stop("'n.cuts' must be a positive integer")
  
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
    if (n.samples <= 0L) return("'n.samples' argument to dbarts must be a positive integer")
    envir$control@n.samples <- n.samples
  } else {
    envir$control@n.samples <- formals(dbarts)[["n.samples"]]
  }

  if (!missing(sigma) && !is.na(sigma)) {
    tryCatch(sigma <- as.double(sigma), warning = function(e)
             stop("'sigma' argument to dbarts must be coerceable to numeric type"))
    if (sigma <= 0.0) stop("'sigma' argument to dbarts must be positive")
    
    envir$sigma <- sigma
  }
}

dbarts <- function(formula, data, test, subset, weights, offset,
                   verbose = FALSE, n.samples = 1000L,
                   tree.prior = cgm, node.prior = normal, resid.prior = chisq,
                   control = dbartsControl(), sigma = NA_real_)
{
  matchedCall <- match.call()

  validateCall <- prepareCallWithArguments(matchedCall, quoteInNamespace(validateArgumentsInEnvironment), "control", "verbose", "n.samples", "sigma")
  validateCall <- addCallArgument(validateCall, 1L, sys.frame(sys.nframe()))
  eval(validateCall, parent.frame(1L), asNamespace("dbarts"))
  
  if (!is.symbol(control@call[[1]])) control@call <- matchedCall
  control@verbose <- verbose

  parseDataCall <- prepareCallWithArguments(matchedCall, quoteInNamespace(parseData), "formula", "data", "test", "subset", "weights", "offset")
  modelMatrices <- eval(parseDataCall, parent.frame(1L))
  
  data <- new("dbartsData", modelMatrices, attr(control, "n.cuts"), sigma)
  attr(control, "n.cuts") <- NULL
  if (is.na(data@sigma))
    data@sigma <- summary(lm(data@y ~ data@x, weights = data@weights))$sigma

  uniqueResponses <- unique(data@y)
  if (length(uniqueResponses) == 2 && all(sort(uniqueResponses) == c(0, 1))) control@binary <- TRUE
  if (!control@binary && !is.null(data@offset) && length(unique(data@offset)) == 1) {
    data@offset <- NULL
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

                  samples
                },
                show = function() {
                  'Pretty prints the object.'
                  
                  cat("dbarts sampler\n")
                  cat("  call: ")
                  writeLines(deparse(control@call))
                  cat("\n")
                  
                  if (FALSE) {
                    xDisplay <- NULL
                    stringConnection <- textConnection("xDisplay", "w", local = TRUE)
                    sink(stringConnection)
                    
                    cat("  x:")
                    if (!is.null(colnames(data@x))) {
                      for (i in 1:min(ncol(data@x), 10)) cat(sprintf(" %8.8s", colnames(data@x)[i]))
                      if (ncol(data@x) > 10) cat(" ...")
                      cat("\n")
                      for (i in 1:min(nrow(data@x), 3)) {
                        cat("    ")
                        for (j in 1:min(ncol(data@x), 10)) cat(sprintf("     %4.4s", data@x[i,j]))
                        cat("\n")
                      }
                    }
                    sink()
                    close(stringConnection)
                    
                    writeLines(xDisplay)
                  }
                  
                },
                setControl = function(control) {
                  'Sets the control object for the sampler to a new one. Preserves the call() slot.'
                  if (!inherits(control, "dbartsControl")) stop("'control' must inherit from dbartsControl.")

                  oldCall <- .self$control@call
                  
                  ptr <- getPointer()
                  control <<- control
                  uniqueResponses <- unique(data@y)
                  if (length(uniqueResponses) == 2 && all(sort(uniqueResponses) == c(0, 1))) control@binary <<- TRUE
                  control@call <<- oldCall
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
                  data@offset <<- as.double(offset)
                  .Call("dbarts_setOffset", ptr, data@offset)

                  if ((is.na(updateState) && control@updateState == TRUE) || identical(updateState, TRUE))
                    storeState(ptr)
                  
                  invisible(NULL)
                },
                setPredictor = function(x, column, updateState = NA) {
                  'Changes a single column of the predictor matrix.'
                  if (is.character(column)) {
                    if (is.null(colnames(data@x))) stop("column names not specified at initialization, so cannot be replaced by name")
                    
                    column <- match(column, colnames(data@x))
                    if (is.na(column)) stop("column name not found in names of current X")
                  }

                  ptr <- getPointer()
                  .Call("dbarts_setPredictor", ptr, as.double(x), as.integer(column))

                  if ((is.na(updateState) && control@updateState == TRUE) || identical(updateState, TRUE))
                    storeState(ptr)

                  invisible(NULL)
                },
                setTestPredictor = function(x.test, column, updateState = NA) {
                  'changes either a single column or the entire test predictor matrix'
                  ptr <- getPointer()
                  
                  if (missing(column)) {
                    x.test <- validateXTest(x.test, ncol(data@x), colnames(data@x))

                    data@x.test <<- x.test
                    .Call("dbarts_setTestPredictors", ptr, data@x.test)
                  } else {
                    if (is.character(column)) {
                      if (is.null(colnames(data@x.test))) stop("column names not specified at initialization, so cannot be replaced by name")
                      
                      column <- match(column, colnames(data@x.test))
                      if (is.na(column)) stop("column name not found in names of current test predictor matrix")
                    }
                    
                    .Call("dbarts_setTestPredictor", ptr, as.double(x.test), as.integer(column))
                  }

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
