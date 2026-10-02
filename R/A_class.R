## the A_ prefix forces this file first in the package's alphabetical
## load order, so these S4 class definitions exist before any other
## R/ file references them

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
methods::setClass("dbartsLeafHyperprior")
methods::setClass(
  "dbartsChiHyperprior",
  contains = "dbartsLeafHyperprior",
  slots = list(degreesOfFreedom = "numeric", scale = "numeric")
)
methods::setValidity("dbartsChiHyperprior", function(object) {
  if (object@degreesOfFreedom <= 0.0) {
    return("'df' must be positive")
  }
  if (object@scale <= 0.0) {
    return("'scale' must be positive")
  }
  TRUE
})
methods::setClass(
  "dbartsFixedHyperprior",
  contains = "dbartsLeafHyperprior",
  slots = list(k = "numeric"),
  prototype = list(k = 2)
)
methods::setValidity("dbartsFixedHyperprior", function(object) {
  if (object@k <= 0.0) {
    return("'k' must be positive")
  }
  TRUE
})

# A law on the leaf prior's sd itself, sd = scale / chi_df: not a law on k, so
# it is deliberately not a dbartsLeafHyperprior and 'k' refuses it. scale = 0
# is the improper sd^-(df + 1) limit.
methods::setClass(
  "dbartsSdHyperprior",
  slots = list(df = "numeric", scale = "numeric")
)
methods::setValidity("dbartsSdHyperprior", function(object) {
  if (
    length(object@df) != 1L ||
      is.na(object@df) ||
      !is.finite(object@df) ||
      object@df <= 0.0
  ) {
    return("invchi() 'df' must be a single positive finite number")
  }
  if (
    length(object@scale) != 1L ||
      is.na(object@scale) ||
      !is.finite(object@scale) ||
      object@scale < 0.0
  ) {
    return("invchi() 'scale' must be a single non-negative finite number")
  }
  TRUE
})


methods::setClass("dbartsLeafPrior")
# an NA spread states no prior; $getLeafPrior reports one where the chains
# disagree, and it must not write back as the unnamed spelling
methods::setValidity("dbartsLeafPrior", function(object) {
  for (name in c("k", "prior.sd")) {
    value <- methods::slot(object, name)
    if (is.numeric(value) && anyNA(value)) {
      return(paste0(
        "the leaf prior's '",
        if (name == "k") "k" else "sd",
        "' is NA, a missing value: $getLeafPrior() reports NA where the ",
        "chains disagree on it; name a value"
      ))
    }
  }
  TRUE
})
# k and prior.sd hold the raw user specification: k is a positive scalar, a
# dbartsLeafHyperprior, or NULL for the family-dependent default; prior.sd is
# NULL (unnamed), a positive scalar, or a dbartsSdHyperprior. At most one is
# non-NULL. Both become the model's leaf.hyperprior and prior.scale when a
# sampler is built.
methods::setClass(
  "dbartsNormalPrior",
  contains = "dbartsLeafPrior",
  slots = list(k = "ANY", prior.sd = "ANY"),
  prototype = list(k = NULL, prior.sd = NULL)
)
# each leaf fits an intercept plus a linear term in the designated
# continuous columns; columns holds the raw user designation (character
# names or numeric indices) until a fitting function resolves it against
# the model matrix, after which it is 1-based integer column indices
methods::setClass(
  "dbartsLinearPrior",
  contains = "dbartsLeafPrior",
  slots = list(k = "ANY", columns = "ANY", prior.sd = "ANY"),
  prototype = list(k = NULL, columns = NULL, prior.sd = NULL)
)
# each leaf fits a smooth Gaussian-process function of the designated
# continuous columns; columns resolves as the linear prior's does.
# lengthscale is NULL for the median-distance heuristic or per-column
# kernel lengthscales on the standardized scale (a scalar recycles when
# resolved); leaves larger than max.leaf.size fall back to constant fits
methods::setClass(
  "dbartsGPPrior",
  contains = "dbartsLeafPrior",
  slots = list(
    k = "ANY",
    columns = "ANY",
    lengthscale = "ANY",
    max.leaf.size = "integer",
    prior.sd = "ANY"
  ),
  prototype = list(
    k = NULL,
    columns = NULL,
    lengthscale = NULL,
    max.leaf.size = 256L,
    prior.sd = NULL
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

## A response family and the settings that ride only on it. `token` is the
## front-door spelling ("gaussian", "student", "hazard.probit", ...) and
## `settings` holds that family's own arguments - a Student-t df, a count
## shape, a hazard time grid - so that no family-specific name has to
## live on a fitting function's signature. The constructors are in
## R/family.R and the resolution into the engine's own family list is in
## R/spec.R.
methods::setClass(
  "dbartsFamily",
  slots = list(token = "character", settings = "list"),
  prototype = list(token = NA_character_, settings = list())
)
methods::setValidity("dbartsFamily", function(object) {
  if (length(object@token) != 1L || is.na(object@token)) {
    return("'token' must be a single family name")
  }
  if (is.null(names(object@settings)) && length(object@settings) > 0L) {
    return("'settings' must be named")
  }
  refused <- familyRefusedSettings(object@token, names(object@settings))
  if (length(refused) > 0L) {
    return(refused)
  }
  TRUE
})

methods::setClass(
  "dbartsControl",
  slots = list(
    binary = "logical",
    verbose = "logical",
    keepTrainingFits = "logical",
    keepFits = "logical",
    useQuantiles = "logical",
    levelGibbs = "logical",
    keepTrees = "logical",
    storage = "character",
    n.samples = "integer",
    n.cuts = "integer",
    n.burn = "integer",
    n.trees = "integer",
    n.chains = "integer",
    n.threads = "integer",
    n.thin = "integer",
    printEvery = "integer",
    printCutoffs = "integer",
    categoricalExhaustiveCap = "integer",
    testFitParallelCutoff = "integer",
    predictParallelCutoff = "integer",
    sparseDensityThreshold = "numeric",
    ## the tree-move mixture, resolved to its six canonical names: it selects
    ## the structure move a sweep proposes and, within a birth/death move,
    ## birth against death
    proposal.probs = "numeric",
    seed = "integer",
    updateState = "logical",
    call = "language"
  ),
  prototype = list(
    binary = FALSE,
    verbose = FALSE,
    keepTrainingFits = TRUE,
    keepFits = TRUE,
    useQuantiles = FALSE,
    levelGibbs = NA,
    keepTrees = FALSE,
    storage = "double",
    n.samples = NA_integer_,
    n.cuts = 100L,
    n.burn = 200L,
    n.trees = 75L,
    n.chains = 4L,
    n.threads = 1L,
    n.thin = 1L,
    printEvery = 100L,
    printCutoffs = 0L,
    categoricalExhaustiveCap = 10L,
    testFitParallelCutoff = 65536L,
    predictParallelCutoff = 50000L,
    sparseDensityThreshold = 0.2,
    proposal.probs = c(
      birth_death = 0.6,
      swap = 0,
      change = 0.4,
      perturb = 0,
      rule_gibbs = 0,
      birth = 0.5
    ),
    seed = NA_integer_,
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
  if (length(object@keepFits) != 1L) {
    return("'keepFits' must be of length 1")
  }
  if (length(object@useQuantiles) != 1L) {
    return("'useQuantiles' must be of length 1")
  }
  if (length(object@levelGibbs) != 1L) {
    return("'levelGibbs' must be of length 1")
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
  if (length(object@n.chains) != 1L) {
    return("'n.chains' must be of length 1")
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

  if (length(object@seed) != 1L) {
    return("'seed' must be of length 1")
  }

  if (is.na(object@verbose)) {
    return("'verbose' must be TRUE/FALSE")
  }
  if (is.na(object@keepTrainingFits)) {
    return("'keepTrainingFits' must be TRUE/FALSE")
  }
  if (is.na(object@keepFits)) {
    return("'keepFits' must be TRUE/FALSE")
  }
  if (is.na(object@useQuantiles)) {
    return("'useQuantiles' must be TRUE/FALSE")
  }
  # levelGibbs, the slot behind treeShift, alone reads NA as a value rather
  # than as a missing one: it is the automatic mode, which takes the step for
  # a forest exactly where that forest's structural mixture is frozen
  if (is.na(object@keepTrees)) {
    return("'keepTrees' must be TRUE/FALSE")
  }

  if (
    length(object@storage) != 1L ||
      is.na(object@storage) ||
      !(object@storage %in% c("double", "single"))
  ) {
    return("'storage' must be \"double\" or \"single\"")
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

  if (is.na(object@printEvery) || object@printEvery < 1L) {
    return("'printEvery' must be a positive integer")
  }
  if (is.na(object@printCutoffs) || object@printCutoffs < 0L) {
    return("'printCutoffs' must be a non-negative integer")
  }

  ## the four engine limits: each is read once, when the sampler is created
  if (
    length(object@categoricalExhaustiveCap) != 1L ||
      is.na(object@categoricalExhaustiveCap) ||
      object@categoricalExhaustiveCap < 2L
  ) {
    return("'categoricalExhaustiveCap' must be a single integer >= 2")
  }
  ## 2^(P - 1) - 1 candidates are enumerated at P present levels, so a cap
  ## past 30 asks for more candidates than an int can index
  if (object@categoricalExhaustiveCap > 30L) {
    return(
      "'categoricalExhaustiveCap' must be at most 30 (2^29 candidates at the cap; 8 to 14 is recommended)"
    )
  }
  if (
    length(object@testFitParallelCutoff) != 1L ||
      is.na(object@testFitParallelCutoff) ||
      object@testFitParallelCutoff < 1L
  ) {
    return("'testFitParallelCutoff' must be a single positive integer")
  }
  if (
    length(object@predictParallelCutoff) != 1L ||
      is.na(object@predictParallelCutoff) ||
      object@predictParallelCutoff < 1L
  ) {
    return("'predictParallelCutoff' must be a single positive integer")
  }
  if (
    length(object@sparseDensityThreshold) != 1L ||
      is.na(object@sparseDensityThreshold) ||
      object@sparseDensityThreshold < 0.0 ||
      object@sparseDensityThreshold > 1.0
  ) {
    return("'sparseDensityThreshold' must be a single number in [0, 1]")
  }

  ## the tree-move mixture, resolved by dbartsControl() before it lands here
  ## and read from this slot by the bridge
  if (
    length(object@proposal.probs) != 6L ||
      !identical(
        names(object@proposal.probs),
        c("birth_death", "swap", "change", "perturb", "rule_gibbs", "birth")
      )
  ) {
    return(paste0(
      "'proposal.probs' must name 'birth_death', 'swap', 'change', ",
      "'perturb', 'rule_gibbs' and 'birth', in that order"
    ))
  }
  proposalProbs <- object@proposal.probs[
    c("birth_death", "swap", "change", "perturb", "rule_gibbs")
  ]
  if (anyNA(proposalProbs) || any(proposalProbs < 0.0 | proposalProbs > 1.0)) {
    return("rule proposal probabilities must be in [0, 1]")
  }
  ## all five exactly zero is the frozen mixture: no structural proposal is
  ## made, so there is no share to normalize
  if (
    sum(proposalProbs) != 0.0 &&
      abs(sum(proposalProbs) - 1.0) >= sqrt(.Machine$double.eps)
  ) {
    return("rule proposal probabilities must sum to 1")
  }
  birth <- object@proposal.probs[["birth"]]
  if (is.na(birth) || birth <= 0.0 || birth >= 1.0) {
    return("birth probability for birth/death step must be in (0, 1)")
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
    leaf.scale = "numeric",
    # The anchor a named leaf-prior sd translates to, in response units: the
    # forest total's prior sd at k = 1, or NA to inherit leaf.scale's
    # family-keyed internal-unit default. It records the named intent, which
    # the sampler re-issues after every channel that re-anchors the response
    # transform and which setLeafPrior rewrites.
    prior.scale = "numeric",
    # "auto" until a fitting function resolves it against the response
    family = "character",

    tree.prior = "dbartsTreePrior",
    leaf.prior = "dbartsLeafPrior",
    leaf.hyperprior = "dbartsLeafHyperprior",
    resid.prior = "dbartsResidPrior"
  ),
  prototype = list(
    leaf.scale = 0.5,
    prior.scale = NA_real_,
    family = "auto",
    tree.prior = new("dbartsCGMPrior"),
    leaf.prior = new("dbartsNormalPrior"),
    leaf.hyperprior = new("dbartsFixedHyperprior"),
    resid.prior = new("dbartsChiSqPrior")
  )
)
methods::setValidity("dbartsModel", function(object) {
  if (object@leaf.scale <= 0.0) {
    return("leaf.scale must be > 0")
  }

  # NaN is not the unnamed spelling, though is.na() accepts it as one: it
  # names no intent and cannot serve as a divisor, so it is refused with every
  # other malformed value rather than read as an absent one
  if (
    length(object@prior.scale) != 1L ||
      is.nan(object@prior.scale) ||
      (!is.na(object@prior.scale) &&
        (!is.finite(object@prior.scale) || object@prior.scale <= 0.0))
  ) {
    return("prior.scale must be NA or a single positive number")
  }

  if (
    length(object@family) != 1L ||
      is.na(object@family) ||
      !(object@family %in%
        c(
          "auto",
          "gaussian",
          "probit",
          "logistic",
          "aft",
          "multinomial",
          "ordinal",
          "nbinom"
        ))
  ) {
    return(
      paste0(
        "'family' must be \"auto\", \"gaussian\", \"probit\", ",
        "\"logistic\", \"aft\", \"multinomial\", \"ordinal\", or \"nbinom\""
      )
    )
  }

  TRUE
})

methods::setClassUnion("matrixOrNULL", c("matrix", "NULL"))
methods::setClassUnion("numericOrNULL", c("numeric", "NULL"))
methods::setClassUnion("listOrNULL", c("list", "NULL"))

methods::setClass(
  "dbartsData",
  slots = list(
    y = "numeric",
    # a dense matrix or a Matrix::dgCMatrix (validated below; the class
    # cannot appear in a slot union without Matrix at load time)
    x = "ANY",
    varTypes = "integer",
    # a dense matrix, a Matrix::dgCMatrix, a mixed container, or null
    # (validated below; validateXTest wraps a bare dgCMatrix as an all-sparse
    # mixed container before storage, so this slot never actually holds a
    # bare dgCMatrix - the class union stays symmetric with 'x' regardless)
    x.test = "ANY",
    weights = "numericOrNULL",
    weights.test = "numericOrNULL",
    offset = "numericOrNULL",
    offset.test = "numericOrNULL",
    n.cuts = "integer",
    sigma = "numeric",
    # always "incorporate" from every entry point: missing predictor values
    # are modelled for whichever rows 'na.action' kept. The slot survives for
    # a host that drives a sampler directly and wants new predictors refused
    # rather than routed.
    missing = "character",
    # what the fit's 'na.action' removed, in the shape base R records it: a
    # named integer vector of dropped row numbers whose class ("omit" or
    # "exclude") decides whether training fits pad back to the caller's own
    # row count through stats::naresid. NULL when nothing was dropped.
    na.action = "ANY",
    # the row names every observation-indexed output of a fit carries:
    # list(train, test), each a character vector or NULL, or NULL when
    # neither set of rows is named. Recorded at entry from the raw inputs
    # (the model frame's rows, or rownames(x) and rownames(test)) after
    # 'subset' and the na.action, outside 'x' so that no container, builder
    # or setter has to carry it. A data object saved before the slot existed
    # lacks it, so every read goes through dataRowNames.
    rowNames = "ANY",
    # the original response's type before it was coded to the doubles the
    # engine reads: "numeric", "factor", "ordered factor", "logical", or
    # "character". The fitters key family = "auto" and the categorical-response
    # refusals off this and response.n.levels rather than re-inspecting the
    # already-coded y (0/1/2 codes are indistinguishable from integer data).
    response.type = "character",
    response.n.levels = "integer",
    # the original response's levels in order, kept for the round-trip an
    # ordinal (cumulative-probit) fit needs to label its K category-probability
    # columns and map an argmax code back to a level. Character for a
    # factor/character response, NULL otherwise (a
    # numeric ordinal fit derives sort(unique(y)) itself). The other families
    # never read it.
    response.levels = "ANY",
    # the per-forest amplitude bases a multi-forest fit combines its forests
    # through, null for an ordinary single-forest fit. Element f is forest
    # f's n x q_f numeric matrix, or NULL for a forest
    # whose basis is the implicit intercept its single amplitude scales. They
    # ride the data object rather than the control, exactly as `weights` does:
    # they are conditioning data, so the setter that replaces one mirrors into
    # this slot and a re-created sampler carries it without further
    # discipline - and, because they arrive at CREATION, a widening applied
    # after a state restore preserves the restored amplitudes rather than the
    # constructed ones.
    bases = "listOrNULL",
    # The multinomial (K-forest softmax) response: an n x K integer matrix
    # whose column k holds category k's counts, null
    # for every other family. It rides the data object rather than an R5 field
    # for the reason the bases do - it is conditioning data, so a re-created
    # sampler carries it without further discipline and a saved fit reloads
    # with it - and @y carries the trials n_i = sum_k counts[i, k], which is
    # what keeps every length(data@y) reader meaningful on such a sampler.
    counts = "matrixOrNULL",
    # The n x K category offset, added to the raw per-category fits BEFORE the
    # softmax; null for none. Not @offset, which is added after the blend and
    # is the softmax's own null direction besides.
    offset.category = "matrixOrNULL",
    # Its test twin, one row per test row rather than per training row.
    offset.category.test = "matrixOrNULL",

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
    na.action = NULL,
    rowNames = NULL,
    response.type = "numeric",
    response.n.levels = NA_integer_,
    response.levels = NULL,
    bases = NULL,
    counts = NULL,
    offset.category = NULL,
    offset.category.test = NULL,

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
    return("'x' must have the same length as 'y'")
  }

  if (
    length(object@varTypes) > 0 &&
      any(
        !object@varTypes %in%
          c(ORDINAL_VARIABLE, CATEGORICAL_VARIABLE, ORDERED_FACTOR_VARIABLE)
      )
  ) {
    return(
      "variable types must all be ordinal, categorical, or ordered factor"
    )
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
    if (
      any(object@weights == 0.0) &&
        !all(object@weights == 0.0 | object@weights == 1.0)
    ) {
      # a supplied value with no effect in the current context is the same
      # condition wherever it recurs - shared with sampleNums/verbose/weights
      # sites elsewhere that are likewise ignored rather than refused. A
      # vector of nothing but 0s and 1s is exempt, on every family: there the
      # zeros ARE the statement being made - which rows are in the data set -
      # rather than an inert value among real weights, and a probit or
      # ordinal fit reads such a vector as its active-row mask outright, so
      # calling the zeros ignored would be false as well as unwanted
      warning(
        "'weights' of 0 will be ignored but increase computation time",
        call. = FALSE
      )
    }
  }
  if (!is.null(object@offset) && length(object@offset) != numObservations) {
    return("'offset' must be null or have length equal to that of 'y'")
  }
  if (!is.null(object@bases)) {
    # a data object is a container, and one forest's basis is a well-formed
    # thing to carry; what refuses a one-forest model is the designed refusal at
    # spec resolution, which is the site both creation routes reach and the only
    # one that can name where the count came from
    if (length(object@bases) < 1L) {
      return("'bases' must be null or name at least one forest")
    }
    for (basis in object@bases) {
      if (is.null(basis)) {
        next
      }
      if (!is.numeric(basis) || NROW(basis) != numObservations) {
        return(paste0(
          "each 'bases' element must be null or a numeric matrix with as ",
          "many rows as 'y' has elements"
        ))
      }
      if (anyNA(basis) || !all(is.finite(basis))) {
        return("'bases' values must all be finite")
      }
    }
  }
  # the multinomial response and its two category offsets. Read bare,
  # unlike every read outside this function: a
  # dbartsData deserialized from a fit saved before the slots existed never
  # reaches here, because validObject rejects it for the missing slots before
  # any validity function runs. Such an object stays READABLE - that is what
  # the guarded accessor buys it - but it is not revalidatable.
  counts <- object@counts
  categoryOffset <- object@offset.category
  categoryTestOffset <- object@offset.category.test
  if (!is.null(counts)) {
    if (!is.integer(counts)) {
      return("'counts' must be null or an integer matrix")
    }
    if (nrow(counts) != numObservations) {
      return("'counts' must have the same number of rows as 'y' has elements")
    }
    # K = 1 names one category with nothing to be distinguished from, and the
    # engine sizes its forest count off K
    if (ncol(counts) < 2L) {
      return("'counts' must have at least two categories")
    }
    if (anyNA(counts)) {
      return("'counts' cannot be NA")
    }
    if (any(counts < 0L)) {
      return("'counts' must all be non-negative")
    }
    # a row with no trial is accepted: it enters no likelihood
    trials <- rowSums(counts)
    # G1a: 'y' is the trials vector, which is what keeps every length(y)
    # reader meaningful on a multinomial data object
    if (!isTRUE(all.equal(as.double(object@y), as.double(trials)))) {
      return(
        "'y' must hold the row sums of 'counts', the per-observation trials"
      )
    }
  } else if (!is.null(categoryOffset) || !is.null(categoryTestOffset)) {
    return(paste0(
      "'offset.category'/'offset.category.test' must be null on a data ",
      "object that carries no 'counts'"
    ))
  }
  if (!is.null(categoryOffset)) {
    if (
      !is.numeric(categoryOffset) ||
        nrow(categoryOffset) != numObservations ||
        ncol(categoryOffset) != ncol(counts)
    ) {
      return(paste0(
        "'offset.category' must be null or a numeric matrix with the ",
        "dimensions of 'counts'"
      ))
    }
    if (anyNA(categoryOffset) || !all(is.finite(categoryOffset))) {
      return("'offset.category' values must all be finite")
    }
  }
  if (!is.null(categoryTestOffset)) {
    # its rows are the TEST rows, so it has nothing to describe without them
    if (is.null(object@x.test)) {
      return("'offset.category.test' must be null when 'x.test' is null")
    }
    if (
      !is.numeric(categoryTestOffset) ||
        nrow(categoryTestOffset) != nrow(object@x.test) ||
        ncol(categoryTestOffset) != ncol(counts)
    ) {
      return(paste0(
        "'offset.category.test' must be null or a numeric matrix with as ",
        "many rows as 'x.test' and as many columns as 'counts'"
      ))
    }
    if (anyNA(categoryTestOffset) || !all(is.finite(categoryTestOffset))) {
      return("'offset.category.test' values must all be finite")
    }
  }
  if (!is.null(object@x.test)) {
    if (
      !is.matrix(object@x.test) &&
        !inherits(object@x.test, "dgCMatrix") &&
        !inherits(object@x.test, "dbartsMixedMatrix")
    ) {
      return(
        "'x.test' must be a matrix, a Matrix::dgCMatrix, a mixed container, or null"
      )
    }
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
    return(paste0("'n.cuts' must have length ", ncol(object@x)))
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

# An unordered-factor predictor in sparse form (the Matrix package's naming
# convention): rows listed in i carry the level coded in values, every
# other row the implicit reference level. A missing value is a stored entry
# whose value is NA, never the reference. Constructor and methods live in
# R/sparseFactor.R.
methods::setClass(
  "sparseFactor",
  slots = list(
    i = "integer", # 0-based rows of the stored entries, ascending
    values = "integer", # 1-based level codes of the stored entries
    levels = "character",
    reference = "character", # the implicit level of unstored rows
    length = "integer"
  ),
  prototype = list(
    i = integer(0),
    values = integer(0),
    levels = "0",
    reference = "0",
    length = 0L
  )
)
methods::setValidity("sparseFactor", function(object) {
  numLevels <- length(object@levels)
  if (anyNA(object@levels)) {
    return("'levels' must be a character vector without NAs")
  }
  if (anyDuplicated(object@levels) > 0L) {
    return("'levels' cannot contain duplicates")
  }
  if (numLevels == 0L) {
    # the one shape without a level: every row a stored missing value
    if (
      length(object@reference) != 1L ||
        !is.na(object@reference) ||
        length(object@length) != 1L ||
        length(object@i) != object@length ||
        !all(is.na(object@values))
    ) {
      return(
        "'levels' can be empty only when every row is a stored missing value"
      )
    }
  } else if (
    length(object@reference) != 1L ||
      is.na(object@reference) ||
      object@reference %not_in% object@levels
  ) {
    return("'reference' must be a single element of 'levels'")
  }
  if (
    length(object@length) != 1L ||
      is.na(object@length) ||
      object@length < 0L
  ) {
    return("'length' must be a single non-negative integer")
  }
  if (length(object@values) != length(object@i)) {
    return("'i' and 'values' must have equal length")
  }
  if (any(object@values < 1L | object@values > numLevels, na.rm = TRUE)) {
    return("'values' must be level codes in [1, length(levels)], or NA")
  }
  if (anyNA(object@i) || any(object@i < 0L | object@i >= object@length)) {
    return("'i' must hold 0-based rows in [0, length)")
  }
  if (length(object@i) > 1L && any(diff(object@i) <= 0L)) {
    return("'i' must be strictly increasing")
  }
  TRUE
})
