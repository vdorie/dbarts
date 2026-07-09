## stupid file name to load first

methods::setClass("dbartsTreePrior")
methods::setClass(
  "dbartsCGMPrior",
  contains = "dbartsTreePrior",
  slots = list(
    power = "numeric",
    base = "numeric",
    splitProbabilities = "numeric",
    # the raw user specification (possibly named, referencing columns);
    # resolved against the data into splitProbabilities when a sampler is
    # built, and NULL thereafter
    splitProbabilitiesSpec = "ANY"
  ),
  prototype = list(splitProbabilitiesSpec = NULL)
)
methods::setValidity("dbartsCGMPrior", function(object) {
  if (object@power <= 0.0) {
    return("'power' must be positive")
  }
  if (object@base <= 0.0 || object@base >= 1.0) {
    return("'base' must be in (0, 1)")
  }
  if (
    length(object@splitProbabilities) > 0L &&
      (any(object@splitProbabilities < 0.0) ||
        abs(sum(object@splitProbabilities) - 1.0) > 1.0e-10)
  ) {
    return("'splitProbabilities' must form a simplex")
  }
  TRUE
})

# DART (Linero 2018): a Dirichlet prior over the split-variable
# probabilities on top of the CGM structure prior; only the bartcore engine
# runs it. rho of NA means the number of predictors; update.delay of NA
# resolves to half the control's burn-in when a sampler is built.
methods::setClass(
  "dbartsDartPrior",
  contains = "dbartsCGMPrior",
  slots = list(
    a = "numeric",
    b = "numeric",
    rho = "numeric",
    alpha = "numeric",
    update.alpha = "logical",
    update.delay = "numeric"
  ),
  prototype = list(
    a = 0.5,
    b = 1.0,
    rho = NA_real_,
    alpha = 1.0,
    update.alpha = TRUE,
    update.delay = NA_real_
  )
)
methods::setValidity("dbartsDartPrior", function(object) {
  if (object@a <= 0.0) {
    return("'a' must be positive")
  }
  if (object@b <= 0.0) {
    return("'b' must be positive")
  }
  if (!is.na(object@rho) && object@rho <= 0.0) {
    return("'rho' must be positive")
  }
  if (object@alpha <= 0.0) {
    return("'alpha' must be positive")
  }
  if (!is.na(object@update.delay) && object@update.delay < 0.0) {
    return("'update.delay' must be non-negative")
  }
  if (
    length(object@splitProbabilities) > 0L ||
      !is.null(object@splitProbabilitiesSpec)
  ) {
    return("a DART prior samples its split probabilities and cannot fix them")
  }
  TRUE
})

# this is a prior over k
methods::setClass("dbartsNodeHyperprior")
methods::setClass(
  "dbartsChiHyperprior",
  contains = "dbartsNodeHyperprior",
  slots = list(degreesOfFreedom = "numeric", scale = "numeric")
)
methods::setValidity("dbartsChiHyperprior", function(object) {
  if (object@degreesOfFreedom <= 0.0) {
    return("'degreesOfFreedom' must be positive")
  }
  if (object@scale <= 0.0) {
    return("'scale' must be positive")
  }
  TRUE
})
methods::setClass(
  "dbartsFixedHyperprior",
  contains = "dbartsNodeHyperprior",
  slots = list(k = "numeric"),
  prototype = list(k = 2)
)
methods::setValidity("dbartsFixedHyperprior", function(object) {
  if (object@k <= 0.0) {
    return("'k' must be positive")
  }
  TRUE
})


methods::setClass("dbartsNodePrior")
# k holds the raw user specification (positive scalar, hyperprior object, or
# NULL for the family-dependent default); it becomes the model's separate
# node.hyperprior when a sampler is built
methods::setClass(
  "dbartsNormalPrior",
  contains = "dbartsNodePrior",
  slots = list(k = "ANY"),
  prototype = list(k = NULL)
)
# each leaf fits an intercept plus a linear term in the designated
# continuous columns; columns holds the raw user designation (character
# names or numeric indices) until a fitting function resolves it against
# the model matrix, after which it is 1-based integer column indices
methods::setClass(
  "dbartsLinearPrior",
  contains = "dbartsNodePrior",
  slots = list(k = "ANY", columns = "ANY"),
  prototype = list(k = NULL, columns = NULL)
)
# each leaf fits a smooth Gaussian-process function of the designated
# continuous columns; columns resolves as the linear prior's does.
# lengthscale is NULL for the median-distance heuristic or per-column
# kernel lengthscales on the standardized scale (a scalar recycles when
# resolved); leaves larger than max.leaf.size fall back to constant fits
methods::setClass(
  "dbartsGPPrior",
  contains = "dbartsNodePrior",
  slots = list(
    k = "ANY",
    columns = "ANY",
    lengthscale = "ANY",
    max.leaf.size = "integer"
  ),
  prototype = list(
    k = NULL,
    columns = NULL,
    lengthscale = NULL,
    max.leaf.size = 256L
  )
)


methods::setClass("dbartsResidPrior")
methods::setClass(
  "dbartsChiSqPrior",
  contains = "dbartsResidPrior",
  slots = list(df = "numeric", quantile = "numeric")
)
methods::setValidity("dbartsChiSqPrior", function(object) {
  if (object@df <= 0.0) {
    return("'df' must be positive")
  }
  if (object@quantile <= 0.0) {
    return("'quantile' must be positive")
  }
  TRUE
})
methods::setClass(
  "dbartsFixedPrior",
  contains = "dbartsResidPrior",
  slots = list(value = "numeric")
)
methods::setValidity("dbartsFixedPrior", function(object) {
  if (object@value <= 0.0) {
    return("'value' must be positive")
  }
  TRUE
})

methods::setClass(
  "dbartsControl",
  slots = list(
    binary = "logical",
    verbose = "logical",
    keepTrainingFits = "logical",
    useQuantiles = "logical",
    keepTrees = "logical",
    n.samples = "integer",
    n.cuts = "integer",
    n.burn = "integer",
    n.trees = "integer",
    n.chains = "integer",
    n.threads = "integer",
    n.thin = "integer",
    printEvery = "integer",
    printCutoffs = "integer",
    rngSeed = "integer",
    updateState = "logical",
    call = "language"
  ),
  prototype = list(
    binary = FALSE,
    verbose = FALSE,
    keepTrainingFits = TRUE,
    useQuantiles = FALSE,
    keepTrees = FALSE,
    n.samples = NA_integer_,
    n.cuts = 100L,
    n.burn = 200L,
    n.trees = 75L,
    n.chains = 4L,
    n.threads = 1L,
    n.thin = 1L,
    printEvery = 100L,
    printCutoffs = 0L,
    rngSeed = NA_integer_,
    updateState = TRUE,
    call = quote(call("NA"))
  )
)

methods::setValidity("dbartsControl", function(object) {
  if (length(object@verbose) != 1L) {
    return("'verbose' must be of length 1")
  }
  if (length(object@keepTrainingFits) != 1L) {
    return("'keepTrainingFits' must be of length 1")
  }
  if (length(object@useQuantiles) != 1L) {
    return("'useQuantiles' must be of length 1")
  }
  if (length(object@keepTrees) != 1L) {
    return("'keepTrees' must be of length 1")
  }

  if (length(object@n.burn) != 1L) {
    return("'n.burn' must be of length 1")
  }
  if (length(object@n.trees) != 1L) {
    return("'n.trees' must be of length 1")
  }
  if (length(object@n.threads) != 1L) {
    return("'n.threads' must be of length 1")
  }
  if (length(object@n.thin) != 1L) {
    return("'n.thin' must be of length 1")
  }

  if (length(object@printEvery) != 1L) {
    return("'printEvery' must be of length 1")
  }
  if (length(object@printCutoffs) != 1L) {
    return("'printCutoffs' must be of length 1")
  }
  if (length(object@updateState) != 1L) {
    return("'updateState' must be of length 1")
  }
  if (length(object@n.samples) != 1L) {
    return("'n.samples' must be of length 1")
  }

  if (length(object@rngSeed) != 1L) {
    return("'rngSeed' must be of length 1")
  }

  if (is.na(object@verbose)) {
    return("'verbose' must be TRUE/FALSE")
  }
  if (is.na(object@keepTrainingFits)) {
    return("'keepTrainingFits' must be TRUE/FALSE")
  }
  if (is.na(object@useQuantiles)) {
    return("'useQuantiles' must be TRUE/FALSE")
  }
  if (is.na(object@keepTrees)) {
    return("'keepTrees' must be TRUE/FALSE")
  }

  if (is.na(object@n.burn) || object@n.burn < 0L) {
    return("'n.burn' must be a non-negative integer")
  }
  if (is.na(object@n.trees) || object@n.trees <= 0L) {
    return("'n.trees' must be a positive integer")
  }
  if (is.na(object@n.chains) || object@n.chains <= 0L) {
    return("'n.chains' must be a positive integer")
  }
  if (is.na(object@n.threads) || object@n.threads <= 0L) {
    return("'n.threads' must be a positive integer")
  }
  if (is.na(object@n.thin) || object@n.thin <= 0L) {
    return("'n.thin' must be a positive integer")
  }

  if (is.na(object@printEvery) || object@printEvery < 0L) {
    return("'printEvery' must be a non-negative integer")
  }
  if (is.na(object@printCutoffs) || object@printCutoffs < 0L) {
    return("'printCutoffs' must be a non-negative integer")
  }

  if (is.na(object@updateState)) {
    return("'updateState' must be TRUE/FALSE")
  }

  ## n.cuts may be length 1 (recycled per-predictor by dbarts()/dbartsData())
  ## or already one value per predictor, so no length-1 constraint here
  if (anyNA(object@n.cuts) || any(object@n.cuts <= 0L)) {
    return("'n.cuts' must contain only positive integers")
  }

  ## handle this in particular b/c it is set through dbarts, not
  ## standard initializer
  if (!is.na(object@n.samples) && object@n.samples < 0L) {
    return("'n.samples' must be a non-negative integer")
  }

  TRUE
})


methods::setClass(
  "dbartsModel",
  slots = list(
    p.birth_death = "numeric",
    p.swap = "numeric",
    p.change = "numeric",

    p.birth = "numeric",

    node.scale = "numeric",
    # "auto" until a fitting function resolves it against the response
    family = "character",

    tree.prior = "dbartsTreePrior",
    node.prior = "dbartsNodePrior",
    node.hyperprior = "dbartsNodeHyperprior",
    resid.prior = "dbartsResidPrior"
  ),
  prototype = list(
    p.birth_death = 1.0,
    p.swap = 0.0,
    p.change = 0.0,
    p.birth = 0.5,
    node.scale = 0.5,
    family = "auto",
    tree.prior = new("dbartsCGMPrior"),
    node.prior = new("dbartsNormalPrior"),
    node.hyperprior = new("dbartsFixedHyperprior"),
    resid.prior = new("dbartsChiSqPrior")
  )
)
methods::setValidity("dbartsModel", function(object) {
  proposalProbs <- c(object@p.birth_death, object@p.swap, object@p.change)
  if (any(proposalProbs < 0.0) || any(proposalProbs > 1.0)) {
    return("rule proposal probabilities must be in [0, 1]")
  }
  if (abs(sum(proposalProbs) - 1.0) >= 1e-10) {
    return("rule proposal probabilities must sum to 1")
  }

  if (object@p.birth <= 0.0 || object@p.birth >= 1.0) {
    return("birth probability for birth/death step must be in (0, 1)")
  }

  if (object@node.scale <= 0.0) {
    return("node.scale must be > 0")
  }

  if (
    length(object@family) != 1L ||
      is.na(object@family) ||
      !(object@family %in% c("auto", "gaussian", "probit", "logistic"))
  ) {
    return(
      "'family' must be \"auto\", \"gaussian\", \"probit\", or \"logistic\""
    )
  }

  TRUE
})

methods::setClassUnion("matrixOrNULL", c("matrix", "NULL"))
methods::setClassUnion("numericOrNULL", c("numeric", "NULL"))

methods::setClass(
  "dbartsData",
  slots = list(
    y = "numeric",
    # a dense matrix or a Matrix::dgCMatrix (validated below; the class
    # cannot appear in a slot union without Matrix at load time)
    x = "ANY",
    varTypes = "integer",
    x.test = "matrixOrNULL",
    weights = "numericOrNULL",
    weights.test = "numericOrNULL",
    offset = "numericOrNULL",
    offset.test = "numericOrNULL",
    n.cuts = "integer",
    sigma = "numeric",
    missing = "character",

    testUsesRegularOffset = "logical"
  ),
  prototype = list(
    y = numeric(0),
    x = matrix(0, 0, 0),
    varTypes = integer(0),
    x.test = NULL,
    weights = NULL,
    weights.test = NULL,
    offset = NULL,
    offset.test = NULL,
    n.cuts = integer(0),
    sigma = NA_real_,
    missing = "incorporate",

    testUsesRegularOffset = NA
  )
)
methods::setValidity("dbartsData", function(object) {
  numObservations <- length(object@y)
  if (
    !is.matrix(object@x) &&
      !inherits(object@x, "dgCMatrix") &&
      !inherits(object@x, "dbartsMixedMatrix")
  ) {
    return("'x' must be a matrix, a Matrix::dgCMatrix, or a mixed container")
  }
  if (nrow(object@x) != numObservations) {
    return("number of rows of 'x' must equal length of 'y'")
  }

  if (
    length(object@varTypes) > 0 &&
      any(
        object@varTypes != ORDINAL_VARIABLE &
          object@varTypes != CATEGORICAL_VARIABLE
      )
  ) {
    return("variable types must all be ordinal or categorical")
  }

  if (!is.null(object@weights)) {
    if (length(object@weights) != numObservations) {
      return("'weights' must be null or have length equal to that of 'y'")
    }
    if (anyNA(object@weights)) {
      return("'weights' cannot be NA")
    }
    if (any(object@weights < 0.0)) {
      return("'weights' must all be non-negative")
    }
    if (any(object@weights == 0.0)) {
      warning("'weights' of 0 will be ignored but increase computation time")
    }
  }
  if (!is.null(object@offset) && length(object@offset) != numObservations) {
    return("'offset' must be null or have length equal to that of 'y'")
  }
  if (!is.null(object@x.test)) {
    if (ncol(object@x.test) != ncol(object@x)) {
      return("'x.test' must be null or have number of columns equal to 'x'")
    }
    if (
      !is.null(object@weights.test) &&
        length(object@weights.test) != nrow(object@x.test)
    ) {
      return(
        "'weights.test' must be null or have the same number of rows as 'x.test'"
      )
    }
    if (
      !is.null(object@offset.test) &&
        length(object@offset.test) != nrow(object@x.test)
    ) {
      return(
        "'offset.test' must be null or have the same number of rows as 'x.test'"
      )
    }
  }
  if (!anyNA(object@n.cuts) && length(object@n.cuts) != ncol(object@x)) {
    return("length of 'n.cuts' must equal number of columns in 'x'")
  }

  if (!is.na(object@sigma) && object@sigma <= 0.0) {
    return("'sigma' must be positive")
  }
  if (
    length(object@missing) != 1L ||
      !(object@missing %in% c("incorporate", "error"))
  ) {
    return("'missing' must be \"incorporate\" or \"error\"")
  }

  TRUE
})
