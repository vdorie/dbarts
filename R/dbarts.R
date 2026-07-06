setMethod("initialize", "dbartsControl", function(.Object, ...) {
  .Object <- callNextMethod()

  validObject(.Object)
  .Object
})

## we don't actually use these defaults; see class definition
## this is only provided for UI hints. Exception is n.cuts, which
## isn't part of class
dbartsControl <- function(
  verbose = FALSE,
  keepTrainingFits = TRUE,
  useQuantiles = FALSE,
  keepTrees = FALSE,
  n.samples = NA_integer_,
  n.cuts = 100L,
  n.burn = 200L,
  n.trees = 75L,
  n.chains = 4L,
  n.threads = dbarts::guessNumCores(),
  n.thin = 1L,
  printEvery = 100L,
  printCutoffs = 0L,
  rngSeed = NA_integer_,
  updateState = TRUE
) {
  result <- new(
    "dbartsControl",
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
    rngSeed = coerceOrError(rngSeed, "integer"),
    updateState = as.logical(updateState)
  )

  n.cuts <- coerceOrError(n.cuts, "integer")
  if (anyNA(n.cuts) || any(n.cuts <= 0L)) {
    stop("'n.cuts' must contain only positive integers")
  }
  attr(result, "n.cuts") <- n.cuts

  result
}

validateArgumentsInEnvironment <- function(
  envir,
  func,
  funcName,
  control,
  verbose,
  n.samples,
  sigma
) {
  controlIsMissing <- missing(control)

  if (!controlIsMissing) {
    if (!inherits(control, "dbartsControl")) {
      stop(
        "'control' argument must be of class dbartsControl; use ",
        "dbartsControl() function to create"
      )
    }
    envir$control <- control
  }

  if (!missing(verbose)) {
    if (!is.logical(verbose) || is.na(verbose)) {
      stop("'verbose' argument to ", funcName, " must be TRUE/FALSE")
    }
  } else if (!controlIsMissing) {
    envir$verbose <- control@verbose
  }

  if (!missing(n.samples)) {
    tryCatch(n.samples <- as.integer(n.samples), warning = function(e) {
      stop(
        "'n.samples' argument to ",
        funcName,
        " must be coercible to integer type"
      )
    })
    if (length(n.samples) != 1L) {
      stop("'n.samples' must be of length 1")
    }
    if (is.null(n.samples)) {
      stop("'n.samples' argument to ", funcName, " cannot be NULL")
    }
    if (is.na(n.samples) || n.samples < 0L) {
      stop(
        "'n.samples' argument to ",
        funcName,
        " must be a non-negative integer"
      )
    }
    envir$control@n.samples <- n.samples
  } else if (controlIsMissing || is.na(control@n.samples)) {
    envir$control@n.samples <- formals(func)[["n.samples"]]
  }

  if (!missing(sigma) && !is.na(sigma)) {
    tryCatch(sigma <- as.double(sigma), warning = function(e) {
      stop(
        "'sigma' argument to ",
        funcName,
        " must be coercible to numeric type"
      )
    })
    if (length(sigma) != 1L) {
      stop("'sigma' must be of length 1")
    }
    if (is.null(sigma) || sigma <= 0.0) {
      stop("'sigma' argument to ", funcName, " must be positive")
    }

    envir$sigma <- sigma
  }
}

dbarts <- function(
  formula,
  data,
  test,
  subset,
  weights,
  offset,
  offset.test = offset,
  verbose = FALSE,
  n.samples = 800L,
  tree.prior = cgm,
  node.prior = normal,
  resid.prior = chisq,
  proposal.probs = c(birth_death = 0.5, swap = 0.1, change = 0.4, birth = 0.5),
  control = dbarts::dbartsControl(),
  sigma = NA_real_,
  seed = NA_integer_,
  factors = c("categorical", "indicators"),
  family = c("auto", "gaussian", "probit", "logistic"),
  missing = c("incorporate", "error")
) {
  matchedCall <- match.call()

  evalEnv <- parent.frame(1L)

  validateCall <- redirectCall(
    matchedCall,
    quoteInNamespace(validateArgumentsInEnvironment)
  )
  validateCall <- addCallArgument(validateCall, 1L, sys.frame(sys.nframe()))
  validateCall <- addCallArgument(validateCall, 2L, dbarts::dbarts)
  validateCall <- addCallArgument(validateCall, 3L, "dbarts")
  eval(validateCall, evalEnv, getNamespace("dbarts"))

  if (length(control@call) == 1L && control@call == call("NA")) {
    control@call <- matchedCall
  }
  control@verbose <- verbose
  # a convenience mirror of dbartsControl(rngSeed = ), as the wrappers expose;
  # an explicit seed overrides the control's, NA leaves it untouched
  seed <- coerceOrError(seed, "integer")
  if (!is.na(seed)) {
    control@rngSeed <- seed
  }

  dataCall <- redirectCall(matchedCall, quoteInNamespace(dbartsData))
  data <- eval(dataCall, evalEnv)
  # cat(
  #   "x address after dbartsData call: ",
  #   .Call("dbarts_getPointerAddress", data@x),
  #   "\n",
  #   sep = ""
  # )

  data@n.cuts <- rep_len(attr(control, "n.cuts"), ncol(data@x))
  data@sigma <- sigma
  attr(control, "n.cuts") <- NULL

  family <- match.arg(family)
  uniqueResponses <- unique(data@y)
  responseIsBinary <- length(uniqueResponses) == 2 &&
    all(sort(uniqueResponses) == c(0, 1))
  if (family == "auto") {
    family <- if (responseIsBinary) "probit" else "gaussian"
  } else if (family != "gaussian" && !responseIsBinary) {
    # gaussian on a 0/1 response is a legitimate request; the binary
    # families need latent-variable coding
    stop("family \"", family, "\" requires a response coded 0/1")
  }
  control@family <- family
  control@binary <- family != "gaussian"

  # binary weight policy, enforced here in the R layer (the bridge keeps the
  # same checks as a backstop for direct-API consumers): a probit has no
  # tractable weighted latent-variable form and is refused; a logistic model
  # treats weights as observation counts (its Polya-Gamma latent is a sum of
  # per-copy draws), so they must be positive integers. Gaussian weights,
  # including a gaussian fit of a 0/1 response, are unrestricted.
  if (!is.null(data@weights)) {
    if (control@family == "probit") {
      stop(
        "probit models do not support weights: a weighted probit has no ",
        "tractable latent-variable form. Use family = \"logistic\" for ",
        "weighted binary regression, or model the latents directly."
      )
    }
    if (control@family == "logistic") {
      w <- data@weights
      if (anyNA(w) || any(w <= 0) || any(w != round(w))) {
        stop(
          "logistic weights are observation counts and must be positive ",
          "integers; drop zero-count rows, and use a gaussian model for ",
          "continuous weights."
        )
      }
    }
  }

  if (is.na(data@sigma) && !control@binary) {
    tryResult <- tryCatch(
      data@sigma <- estimateSigmaFromLinearModel(data),
      error = function(e) e
    )
    if (inherits(tryResult, "error")) {
      stop("unable to obtain a starting estimate of sigma; provide one instead")
    }
  }

  # bart will passthrough with offset == something no matter what, which we
  # can NULL out
  if (!control@binary && !is.null(data@offset) && all(data@offset == 0.0)) {
    data@offset <- NULL
  }
  if (
    !control@binary &&
      !is.null(data@offset.test) &&
      all(data@offset.test == 0.0)
  ) {
    data@offset.test <- NULL
  }
  # cat(
  #   "x address after updating data: ",
  #   .Call("dbarts_getPointerAddress", data@x),
  #   "\n",
  #   sep = ""
  # )

  parsePriorsCall <- redirectCall(matchedCall, quoteInNamespace(parsePriors))
  parsePriorsCall <- setDefaultsFromFormals(
    parsePriorsCall,
    formals(dbarts),
    "tree.prior",
    "node.prior",
    "resid.prior"
  )
  parsePriorsCall$control <- control
  parsePriorsCall$data <- data
  parsePriorsCall$parentEnv <- evalEnv

  if (control@binary) {
    if (any(names(parsePriorsCall) == "resid.prior")) {
      parsePriorsCall[[which(
        names(parsePriorsCall) == "resid.prior"
      )]] <- quote(fixed(1))
    } else {
      parsePriorsCall[[length(parsePriorsCall) + 1L]] <- quote(fixed(1))
      names(parsePriorsCall)[length(parsePriorsCall)] <- "resid.prior"
    }
  }
  priors <- eval(parsePriorsCall)

  model <- new(
    "dbartsModel",
    priors$tree.prior,
    priors$node.prior,
    priors$node.hyperprior,
    priors$resid.prior,
    proposal.probs = proposal.probs,
    # the logistic scale is probit's default widened by the logistic
    # latent's standard deviation, pi / sqrt(3)
    node.scale = switch(
      control@family,
      gaussian = 0.5,
      probit = 3.0,
      logistic = pi * sqrt(3.0)
    )
  )

  result <- new("dbartsSampler", control, model, data)
  # cat(
  #   "x address after creating sampler: ",
  #   .Call("dbarts_getPointerAddress", data@x),
  #   "\n",
  #   sep = ""
  # )
  # cat(
  #   "x address in sampler$data: ",
  #   .Call("dbarts_getPointerAddress", result$data@x),
  #   "\n",
  #   sep = ""
  # )
  result
}


dbartsSampler <- setRefClass(
  "dbartsSampler",
  fields = list(
    pointer = "externalptr",
    control = "dbartsControl",
    model = "dbartsModel",
    data = "dbartsData",
    state = "ANY" # is either a list of states, or a promise to evaluate
  ),
  methods = list(
    initialize = function(control, model, data, ...) {
      if (!inherits(control, "dbartsControl")) {
        stop("'control' must inherit from dbartsControl")
      }
      if (!inherits(model, "dbartsModel")) {
        stop("'model' must inherit from dbartsModel")
      }
      if (!inherits(data, "dbartsData")) {
        stop("'data' must inherit from dbartsData")
      }
      .self$control <- control
      .self$model <- model
      .self$data <- data

      # "auto" (a hand-built control) keeps the bridge's own dispatch
      .self$pointer <- .Call(
        C_dbarts_bartcore_create,
        .self$control,
        .self$model,
        .self$data,
        if (control@family == "auto") "" else control@family
      )
      # materialized lazily on first access (forcing it before saveRDS
      # captures the sampler), or eagerly by storeState / updateState runs.
      # A deserialized object can force the promise after its pointer has
      # died; that must yield NULL - no stored state - not a C error
      delayedAssign(
        "state",
        {
          if (
            control@updateState &&
              .Call(C_dbarts_bartcore_isValidPointer, pointer)
          ) {
            .Call(C_dbarts_bartcore_storeState, pointer)
          } else {
            NULL
          }
        },
        eval.env = as.environment(.self),
        assign.env = as.environment(.self)
      )

      callSuper(...)
    },
    run = function(
      numBurnIn,
      numSamples,
      updateState = NA,
      numThreads = control@n.threads
    ) {
      "Runs the posterior sampler and returns a list with the results."
      if (missing(numBurnIn)) {
        numBurnIn <- NA_integer_
      }
      if (missing(numSamples)) {
        numSamples <- NA_integer_
      }

      samples <- bartcoreSamplerRun(.self, numBurnIn, numSamples)
      if (
        (is.na(updateState) && control@updateState == TRUE) ||
          identical(updateState, TRUE)
      ) {
        storeState()
      }
      if (is.null(samples)) {
        return(invisible(NULL))
      }
      samples
    },
    sampleTreesFromPrior = function(updateState = NA) {
      "Draws tree structure from prior"
      ptr <- getPointer()
      .Call(C_dbarts_bartcore_sampleTreesFromPrior, ptr)

      if (
        (is.na(updateState) && control@updateState == TRUE) ||
          identical(updateState, TRUE)
      ) {
        storeState(ptr)
      }

      invisible(NULL)
    },
    sampleNodeParametersFromPrior = function(updateState = NA) {
      "Draws end node parameters from prior; does not update tree structure."
      ptr <- getPointer()
      .Call(C_dbarts_bartcore_sampleNodeParametersFromPrior, ptr)

      if (
        (is.na(updateState) && control@updateState == TRUE) ||
          identical(updateState, TRUE)
      ) {
        storeState(ptr)
      }

      invisible(NULL)
    },
    copy = function(shallow = FALSE) {
      "Creates a deep or shallow copy of the sampler."
      dupe <-
        if (shallow) {
          dbartsSampler$new(control, model, data)
        } else {
          newData <- data
          # only need to dupe things that can be changed internally, as the
          # rest will be simply swapped out
          newData@x <- .Call(C_dbarts_deepCopy, data@x)
          if (!is.null(data@x.test)) {
            newData@x.test <- .Call(C_dbarts_deepCopy, data@x.test)
          }
          dbartsSampler$new(control, model, newData)
        }

      # the stored state is opaque and never mutated in place (storeState
      # replaces it whole), so the copy can install the same object
      if (!is.null(state)) {
        dupe$setState(state)
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
    predict = function(x.test, offset.test, n.threads = control@n.threads) {
      'Using existing sampler to predict for new data without re-running.'
      ptr <- getPointer()

      x.test <- validateXTest(x.test, data@x)
      if (is.null(x.test)) {
        stop("x.test cannot be NULL")
      }

      if (missing(offset.test) || is.null(offset.test)) {
        offset.test <- NA_real_
      } else {
        offset.test <- as.double(offset.test)
        if (length(offset.test) == 1L) {
          offset.test <- rep_len(offset.test, nrow(x.test))
        }

        if (!identical(length(offset.test), nrow(x.test))) {
          stop(
            "length of test offset must be equal to number of rows in test matrix"
          )
        }
      }

      # the engine runs prediction serially; a missing offset is passed
      # as NULL rather than an NA sentinel
      if (length(offset.test) == 1L && is.na(offset.test)) {
        offset.test <- NULL
      }
      .Call(C_dbarts_bartcore_predict, ptr, x.test, offset.test)
    },
    setControl = function(newControl) {
      'Sets the control object for the sampler to a new one. Preserves the call() slot.'
      if (!inherits(newControl, "dbartsControl")) {
        stop("'control' must inherit from dbartsControl")
      }

      selfEnv <- parent.env(environment())

      newControl@binary <- control@binary
      newControl@family <- control@family
      newControl@call <- control@call

      # settings fixed at creation: the generators and anything shaping
      # the cut grid
      for (slotName in c(
        "n.trees",
        "n.chains",
        "useQuantiles",
        "rngSeed"
      )) {
        if (!identical(slot(newControl, slotName), slot(control, slotName))) {
          stop(
            "the bartcore engine cannot change '",
            slotName,
            "' on an existing sampler"
          )
        }
      }
      if (newControl@keepTrees && is.na(newControl@n.samples)) {
        stop("keepTrees requires 'n.samples' to be specified")
      }

      ptr <- getPointer()
      selfEnv$control <- newControl
      .Call(C_dbarts_bartcore_setControl, ptr, control)

      invisible(NULL)
    },
    setModel = function(newModel) {
      'Sets the model object for the sampler to a new one.'
      if (!inherits(newModel, "dbartsModel")) {
        stop("'model' must inherit from dbartsModel")
      }
      # the Dirichlet machinery is fixed at creation: a sampler cannot gain
      # or reconfigure it
      if (
        is(newModel@tree.prior, "dbartsDartPrior") ||
          is(model@tree.prior, "dbartsDartPrior")
      ) {
        stop("setModel cannot change a DART tree prior; recreate the sampler")
      }
      ptr <- getPointer()
      selfEnv <- parent.env(environment())

      oldModel <- model
      selfEnv$model <- newModel
      tryResult <- tryCatch(
        .Call(C_dbarts_bartcore_setModel, ptr, selfEnv$model, control, data),
        error = function(e) {
          selfEnv$model <- oldModel
          e$call <- quote(.Call(
            C_dbarts_bartcore_setModel,
            ptr,
            selfEnv$model,
            control,
            data
          ))
          e
        }
      )
      if (inherits(tryResult, "error")) {
        stop(tryResult)
      }

      invisible(NULL)
    },
    setData = function(newData, updateState = NA) {
      'Sets the data object for the sampler to a new one. Preserves the n.cuts and sigma slots.'
      if (
        data@missing == "error" &&
          (anyNA(newData@x) ||
            (!is.null(newData@x.test) && anyNA(newData@x.test)))
      ) {
        stop(
          "new predictors contain missing values and the sampler was built with missing = \"error\""
        )
      }
      bartcoreSamplerSetData(.self, newData)
    },
    setResponse = function(y, updateState = NA) {
      'Changes the response against which the sampler is fitted.'
      bartcoreSamplerSetResponse(.self, y)
    },
    setOffset = function(offset, updateScale = FALSE, updateState = NA) {
      'Changes the offset slot used to adjust the response.'
      bartcoreSamplerSetOffset(.self, offset, updateScale)
    },
    setWeights = function(weights, updateState = NA) {
      'Changes the weights with which the sampler is fitted.'
      ptr <- getPointer()
      selfEnv <- parent.env(environment())

      oldWeights <- data@weights
      selfEnv$data@weights <- as.double(weights)
      tryResult <- tryCatch(
        .Call(C_dbarts_bartcore_setWeights, ptr, data@weights),
        error = function(e) {
          selfEnv$data@weights <- oldWeights
          e
        }
      )
      if (inherits(tryResult, "error")) {
        stop(tryResult)
      }

      invisible(NULL)
    },
    setSigma = function(sigma, updateState = NA) {
      'Changes the residual standard deviation parameter for each chain.'
      ptr <- getPointer()
      .Call(C_dbarts_bartcore_setSigma, ptr, as.double(sigma))
      invisible(NULL)
    },
    setPredictor = function(
      x,
      column,
      forceUpdate,
      updateCutPoints = FALSE,
      updateState = NA
    ) {
      "Changes a single column of the predictor matrix, or the entire matrix if column is missing."

      if (data@missing == "error" && anyNA(x)) {
        stop(
          "new predictors contain missing values and the sampler was built with missing = \"error\""
        )
      }
      bartcoreSamplerSetPredictor(
        .self,
        x,
        column = if (missing(column)) NULL else column,
        forceUpdate = if (missing(forceUpdate)) NULL else forceUpdate,
        updateCutPoints = updateCutPoints
      )
    },
    setCutPoints = function(cuts, column, updateState = NA) {
      'Changes the cut points for the predictors in column, or the entire set itself if the column argument is missing. Forces the change by pruning any leaves that end up empty.'

      bartcoreSamplerSetCutPoints(
        .self,
        cuts,
        column = if (missing(column)) NULL else column
      )
    },
    setTestPredictor = function(x.test, column) {
      'Changes a single column of the test predictor matrix.'

      if (data@missing == "error" && anyNA(x.test)) {
        stop(
          "new test predictors contain missing values and the sampler was built with missing = \"error\""
        )
      }
      bartcoreSamplerSetTestPredictor(
        .self,
        x.test,
        column = if (missing(column)) NULL else column
      )
    },
    setTestPredictorAndOffset = function(x.test, offset.test) {
      'Changes the test predictor matrix, and optionally the test offset.'
      if (data@missing == "error" && !is.null(x.test) && anyNA(x.test)) {
        stop(
          "new test predictors contain missing values and the sampler was built with missing = \"error\""
        )
      }
      if (missing(offset.test)) {
        # predictors only; the engine keeps the current offset and the
        # bridge refuses if the row count would orphan its length
        return(bartcoreSamplerSetTestPredictor(.self, x.test, column = NULL))
      }

      x.test <- validateXTest(x.test, data@x)
      if (is.null(x.test) && !is.null(offset.test)) {
        stop("when test matrix is NULL, test offset must be as well")
      }
      if (!is.null(offset.test)) {
        offset.test <- as.double(offset.test)
        if (length(offset.test) == 1L) {
          offset.test <- rep_len(offset.test, nrow(x.test))
        }
        if (!identical(length(offset.test), nrow(x.test))) {
          stop(
            "length of test offset must be equal to number of rows in test matrix"
          )
        }
      }

      selfEnv <- parent.env(environment())
      oldTestUsesRegularOffset <- data@testUsesRegularOffset
      oldX.test <- data@x.test
      oldOffset.test <- data@offset.test

      selfEnv$data@testUsesRegularOffset <- FALSE
      selfEnv$data@x.test <- x.test
      selfEnv$data@offset.test <- offset.test
      tryResult <- tryCatch(
        .Call(
          C_dbarts_bartcore_setTestPredictorAndOffset,
          getPointer(),
          data@x.test,
          data@offset.test
        ),
        error = function(e) {
          selfEnv$data@testUsesRegularOffset <- oldTestUsesRegularOffset
          selfEnv$data@x.test <- oldX.test
          selfEnv$data@offset.test <- oldOffset.test
          e
        }
      )
      if (inherits(tryResult, "error")) {
        stop(tryResult)
      }
      invisible(NULL)
    },
    setTestOffset = function(offset.test) {
      'Changes the test offset.'
      ptr <- getPointer()
      selfEnv <- parent.env(environment())

      selfEnv$data@testUsesRegularOffset <- FALSE
      if (!is.null(offset.test)) {
        if (is.null(data@x.test)) {
          stop("when test matrix is NULL, test offset must be as well")
        }
        if (length(offset.test) == 1L) {
          offset.test <- rep_len(offset.test, nrow(data@x.test))
        }
        if (length(offset.test) != nrow(data@x.test)) {
          stop(
            "length of test offset must be equal to number of rows in test matrix"
          )
        }
      }
      oldOffset.test <- data@offset.test
      selfEnv$data@offset.test <- offset.test
      tryResult <- tryCatch(
        .Call(C_dbarts_bartcore_setTestOffset, ptr, data@offset.test),
        error = function(e) {
          selfEnv$data@offset.test <- oldOffset.test
          e
        }
      )
      if (inherits(tryResult, "error")) {
        stop(tryResult)
      }

      invisible(NULL)
    },
    getLatents = function(result) {
      'For binary models, returns the current value of the latent variable representation.'
      resultIsMissing <- missing(result)

      ptr <- getPointer()

      .Call(
        C_dbarts_bartcore_getLatents,
        ptr,
        if (resultIsMissing) NULL else result
      )
    },
    getSigmas = function(result) {
      'Return current residual error term on original, standard deviation scale.'

      ptr <- getPointer()
      .Call(C_dbarts_bartcore_getSigmas, ptr)
    },
    getSumsOfSquaredResiduals = function(result) {
      'Return sum( (y - y.hat)^2 ) on original scale.'
      ptr <- getPointer()
      .Call(C_dbarts_bartcore_getSumsOfSquaredResiduals, ptr)
    },
    getPointer = function() {
      'Returns the underlying reference pointer, checking for consistency first.'
      selfEnv <- parent.env(environment())

      if (.Call(C_dbarts_bartcore_isValidPointer, pointer) == FALSE) {
        if (is.null(state)) {
          stop(
            "samplers cannot be re-created without a stored state; call ",
            "storeState() - or force $state, see the Saving section of ",
            "?bart - before serializing"
          )
        }
        selfEnv$pointer <- .Call(
          C_dbarts_bartcore_create,
          control,
          model,
          data,
          if (control@family == "auto") "" else control@family
        )
        .Call(C_dbarts_bartcore_setState, pointer, state)
      }
      pointer
    },
    setState = function(newState) {
      'Sets the internal state from a cache.'
      if (!inherits(newState, "bartcoreState")) {
        stop("'state' must inherit from bartcoreState")
      }
      selfEnv <- parent.env(environment())
      if (.Call(C_dbarts_bartcore_isValidPointer, pointer) == FALSE) {
        selfEnv$pointer <- .Call(
          C_dbarts_bartcore_create,
          control,
          model,
          data,
          if (control@family == "auto") "" else control@family
        )
      }
      .Call(C_dbarts_bartcore_setState, pointer, newState)
      selfEnv$state <- newState
      invisible(NULL)
    },
    storeState = function(ptr = getPointer()) {
      'Updates the cached internal state used for saving/loading.'
      selfEnv <- parent.env(environment())
      selfEnv$state <- .Call(C_dbarts_bartcore_storeState, ptr)
      invisible(NULL)
    },
    printTrees = function(treeNums, chainNums, sampleNums) {
      'Produces an info dump of the internal state of the trees.'
      matchedCall <- match.call()
      if (is.null(matchedCall$chainNums)) {
        chainNums <- seq_len(control@n.chains)
      }
      if (is.null(matchedCall$sampleNums)) {
        sampleNums <- if (control@keepTrees) {
          seq_len(control@n.samples)
        } else {
          NULL
        }
      } else {
        if (!control@keepTrees) {
          warning("sampleNums ignored if keepTrees is FALSE")
          sampleNums <- NULL
        } else {
          sampleNums <- as.integer(sampleNums)
        }
      }
      if (is.null(matchedCall$treeNums)) {
        treeNums <- seq_len(control@n.trees)
      }

      ptr <- getPointer()
      invisible(.Call(
        C_dbarts_bartcore_printTrees,
        ptr,
        as.integer(chainNums),
        sampleNums,
        as.integer(treeNums)
      ))
    },
    getTrees = function(
      treeNums,
      chainNums,
      sampleNums,
      current = FALSE,
      newdata = NULL
    ) {
      "Returns a data.frame containing the internal state of the trees."
      matchedCall <- match.call()
      current <- isTRUE(current)
      # live working trees have no sample dimension, so treat a current request
      # like a non-keepTrees sampler for sample handling
      useSaved <- control@keepTrees && !current
      if (is.null(matchedCall$chainNums)) {
        chainNums <- seq_len(control@n.chains)
      }
      if (is.null(matchedCall$sampleNums)) {
        sampleNums <- if (useSaved) {
          seq_len(control@n.samples)
        } else {
          NULL
        }
      } else {
        if (!useSaved) {
          warning(
            if (current) {
              "sampleNums ignored if current is TRUE"
            } else {
              "sampleNums ignored if keepTrees is FALSE"
            }
          )
          sampleNums <- NULL
        } else {
          sampleNums <- as.integer(sampleNums)
        }
      }
      if (is.null(matchedCall$treeNums)) {
        treeNums <- seq_len(control@n.trees)
      }

      chainNums <- as.integer(chainNums)
      treeNums <- as.integer(treeNums)

      if (any(chainNums <= 0 | chainNums > control@n.chains)) {
        stop("chainNums must be in [1, ", control@n.chains, "]")
      }
      if (
        useSaved &&
          any(sampleNums <= 0 | sampleNums > control@n.samples)
      ) {
        stop("sampleNums must be in [1, ", control@n.samples, "]")
      }
      if (any(treeNums <= 0 | treeNums > control@n.trees)) {
        stop("treeNums must be in [1, ", control@n.trees, "]")
      }

      # route new data through the trees so 'n' counts that data instead of the
      # training predictors; validated and coded as for predict
      if (!is.null(newdata)) {
        newdata <- validateXTest(newdata, data@x)
      }

      ptr <- getPointer()
      trees <- .Call(
        C_dbarts_bartcore_getTrees,
        ptr,
        chainNums,
        sampleNums,
        treeNums,
        current,
        newdata
      )
      # categorical rules report their raw direction mask in 'value'; when any
      # column can hold one, decode the masks into per-level L/R strings
      if (any(data@varTypes == CATEGORICAL_VARIABLE)) {
        trees <- decodeCategoricalSplits(trees, data@x, data@varTypes)
      }
      # rules on columns with missing values report their NA route
      if (!is.null(trees$missing)) {
        trees$missing <- c("L", "R")[trees$missing + 1L]
      }
      # linear leaves report one generically named slope column per
      # covariate; name them after the designated columns
      if (is(model@node.prior, "dbartsLinearPrior")) {
        covariateNames <- colnames(data@x)[model@node.prior@columns]
        if (!is.null(covariateNames)) {
          slopeColumns <- match(
            paste0("beta.", seq_along(covariateNames)),
            names(trees)
          )
          names(trees)[slopeColumns] <- paste0("beta.", covariateNames)
        }
      }
      trees
    },
    plotTree = function(
      treeNum,
      chainNum,
      sampleNum,
      treePlotPars = c(nodeHeight = 12, nodeWidth = 40, nodeGap = 8),
      ...
    ) {
      'Minimialist visualization of tree branching and contents.'

      matchedCall <- match.call()
      if (is.null(matchedCall$chainNum)) {
        if (control@n.chains == 1L) {
          chainNum <- 1L
        } else {
          stop("chainNum required if more than one chain in sampler")
        }
      }
      if (is.null(matchedCall$sampleNum)) {
        sampleNum <- if (control@keepTrees) control@n.samples else 1L
      }

      tree <-
        if (control@keepTrees) {
          .self$getTrees(treeNum, chainNum, sampleNum)
        } else {
          .self$getTrees(treeNum, chainNum)
        }

      maxDepth <- getTreeDepthAndSize(tree)[["depth"]]

      tree <- cbind(
        tree,
        y = numeric(nrow(tree)),
        x = numeric(nrow(tree)),
        index = integer(nrow(tree))
      )
      tree <- fillPlotCoordinatesForNode(tree, maxDepth, 1L, 1L)
      numEndNodes <- tree$index[1L] - 1L
      #tree$index <- NULL

      plotHeight <- treePlotPars[["nodeHeight"]] *
        maxDepth +
        treePlotPars[["nodeGap"]] * (maxDepth - 1)
      dotsList <- list(...)
      dotsList$mar <- c(0, 0, 0, 0)
      par(dotsList)
      plot(
        NULL,
        type = "n",
        bty = "n",
        xaxt = "n",
        yaxt = "n",
        xlab = "",
        ylab = "",
        xlim = c(0, treePlotPars[["nodeWidth"]] * numEndNodes),
        ylim = c(0, plotHeight)
      )
      plotNode(tree, .self, treePlotPars)

      invisible(NULL)
    },
    startThreads = function(n.threads = control@n.threads) {
      # the engine manages worker threads within each run
      invisible(NULL)
    },
    stopThreads = function() {
      invisible(NULL)
    }
  )
)
