setMethod("initialize", "dbartsControl",
          function(.Object, n.cuts = 100L, ...)
{
  .Object <- callNextMethod()
  ## keep around but sorta hide
  attr(.Object, "n.cuts") <- as.integer(n.cuts)

  validObject(.Object)
  .Object
})

validateArgumentsInEnvironment <- function(control, verbose, n.samples, sigma, envir) {
  if (!missing(control) && !inherits(control, "dbartsControl"))
    stop("'control' argument must be of class dbartsControl; use dbartsControl() function to create.")
  envir$control <- control

  
  if (!is.logical(verbose) || is.na(verbose)) stop("'verbose' argument to dbarts must be TRUE/FALSE.")
  tryCatch(n.samples <- as.integer(n.samples), warning = function(e)
           stop("'n.samples' argument to dbarts must be coerceable to integer type."))
  if (n.samples <= 0L) return("'n.samples' argument to dbarts must be a positive integer.")
  envir$control@n.samples <- n.samples
  envir$control@verbose <- verbose

  if (!is.na(sigma)) {
    tryCatch(sigma <- as.double(sigma), warning = function(e)
             stop("'sigma' argument to dbarts must be coerceable to numeric type."))
    if (sigma <= 0.0) stop("'sigma' argument to dbarts must be positive.")
  }

  envir$sigma <- sigma
}

dbarts <- function(formula, data, test, subset, weights, offset,
                   verbose = TRUE, n.samples = 1000L,
                   tree.prior = cgm, node.prior = normal, resid.prior = chisq,
                   control = dbartsControl(), sigma = NA_real_)
{
  matchedCall <- match.call()

  do.call(validateArgumentsInEnvironment,
          list(control = control, verbose = verbose, n.samples = n.samples, sigma = sigma, envir = sys.frame(sys.nframe())),
          envir = parent.frame(1L))
  if (!is.symbol(control@call[[1]])) control@call <- matchedCall

  parseDataCall <- stripExtraArguments(matchedCall, "formula", "data", "test", "subset", "weights", "offset")
  parseDataCall[[1L]] <- quote(dbarts:::parseData)
  modelMatrices <- eval(parseDataCall, parent.frame(1L))
  
  data <- new("dbartsData", modelMatrices, attr(control, "n.cuts"), sigma)
  attr(control, "n.cuts") <- NULL
  if (is.na(data@sigma)) data@sigma <- summary(lm(data@y ~ data@x))$sigma

  uniqueResponses <- unique(data@y)
  if (length(uniqueResponses) == 2 && all(sort(uniqueResponses) == c(0, 1))) control@binary <- TRUE
  if (!control@binary && !is.null(data@offset) && length(unique(data@offset)) == 1) {
    data@offset <- NULL
  }

  parsePriorsCall <- stripExtraArguments(matchedCall, "tree.prior", "node.prior", "resid.prior")
  parsePriorsCall <- setDefaultsFromFormals(parsePriorsCall, formals(dbarts), "tree.prior", "node.prior", "resid.prior")
  parsePriorsCall[[1]] <- quote(dbarts:::parsePriors)
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
                  function(control, model, data)
                {
                  if (!inherits(control, "dbartsControl")) stop("'control' must inherit from dbartsControl.")
                  if (!inherits(model, "dbartsModel")) stop("'model' must inherit from dbartsModel.")
                  if (!inherits(data, "dbartsData")) stop("'data' must inherit from dbartsData.")

                  control <<- control
                  model   <<- model
                  data    <<- data
                  state   <<- NULL
                  pointer <<- .Call("dbarts_create", control, model, data)
                },
                run = function(numBurnIn, numSamples) {
                  if (missing(numBurnIn))  numBurnIn  <- NA_integer_
                  if (missing(numSamples)) numSamples <- NA_integer_
                  
                  ptr <- getPointer()
                  samples <- .Call("dbarts_run", ptr, as.integer(numBurnIn), as.integer(numSamples))

                  updateState(ptr)

                  return(samples)
                },
                setY = function(y) {
                  .Call("dbarts_setY", getPointer(), as.double(y))

                  invisible(updateState(ptr))
                },
                getPointer = function() {
                  if (.Call("dbarts_isValidPointer", pointer) == FALSE) {
                    oldVerbose <- control@verbose
                    control@verbose <<- FALSE
                    pointer <<- .Call("dbarts_create", control, model, data)
                    if (!is.null(state)) {
                      state <<- .Call("dbarts_restoreState", pointer, state)
                    }
                    control@verbose <<- oldVerbose
                  }
                  
                  return(pointer)
                },
                updateState = function(ptr) {
                  if (is.null(state)) {
                    state <<- .Call("dbarts_createState", ptr)
                  } else {
                    .Call("dbarts_storeState", ptr, state)
                  }
                })
              )
