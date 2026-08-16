# expects n.pars x n.samples x n.chains if n.chains > 1
# returns n.chains x n.samples x n.pars or (n.chains x n.samples) x n.pars
#
# for quantities such as yhat.train, n.pars is n.obs; for sigma, it is 1
# and the dimension is dropped
#
# preserves the per-parameter dimnames if they exist (for ranef)
convertSamplesFromDbartsToBart <-
  function(
    samples,
    n.chains = dim(samples)[length(dim(samples))],
    combineChains = FALSE
  ) {
    if (!combineChains) {
      ifelse_3(
        is.null(dim(samples)),
        length(dim(samples)) == 2L,
        samples,
        t(samples),
        aperm(samples, c(3L, 2L, 1L))
      )
    } else {
      ifelse_3(
        is.null(dim(samples)),
        length(dim(samples)) == 2L,
        samples,
        if (n.chains <= 1L) t(samples) else as.vector(t(samples)),
        {
          x <- NULL
          res <- evalx(dim(samples), t(matrix(samples, x[1L], prod(x[-1L]))))
          if (!is.null(dimnames(samples))) {
            colnames(res) <- dimnames(samples)[[1L]]
          }
          res
        }
      )
    }
  }

# expects n.chains x n.samples x n.pars or (n.chains x n.samples) x n.pars
# returns n.pars x n.samples x n.chains or n.pars x (n.samples x n.chains)
convertSamplesFromBartsToDbarts <- function(
  samples,
  n.chains,
  uncombineChains = FALSE
) {
  if (!uncombineChains) {
    ifelse_3(
      is.null(dim(samples)),
      is.matrix(samples),
      samples,
      t(samples),
      aperm(samples, c(3L, 2L, 1L))
    )
  } else {
    x <- NULL
    if (is.null(dim(samples))) {
      res <- if (n.chains == 1L) {
        samples
      } else {
        matrix(samples, length(samples) %/% n.chains, n.chains)
      }
      evalx(
        dimnames(samples),
        if (!is.null(x)) dimnames(res) <- list(x[[length(x)]], NULL)
      )
      res
    } else {
      res <- if (n.chains == 1L) {
        samples
      } else {
        array(
          t(samples),
          c(ncol(samples), nrow(samples) %/% n.chains, n.chains)
        )
      }
      evalx(
        dimnames(samples),
        if (!is.null(x)) dimnames(res) <- list(x[[length(x)]], NULL, NULL)
      )
      res
    }
  }
}

# input n.samples x n.chains x n.pars, or n.samples x n.pars when n.chains = 1
# output (n.samples * n.chains) x n.pars
combineChains <- function(samples) {
  ifelse_3(
    is.null(dim(samples)),
    length(dim(samples)) == 2L,
    samples,
    as.vector(samples),
    {
      x <- NULL
      res <- evalx(
        dim(samples),
        matrix(aperm(samples, c(2L, 1L, 3L)), prod(x[1L:2L]), x[3L])
      )
      if (!is.null(dimnames(samples))) {
        dimnames(res) <- evalx(dimnames(samples), list(NULL, x[[length(x)]]))
      }
      res
    }
  )
}

uncombineChains <- function(samples, n.chains) {
  if (is.null(dim(samples))) {
    if (n.chains == 1L) {
      samples
    } else {
      matrix(samples, n.chains, length(samples) %/% n.chains)
    }
  } else {
    res <- if (n.chains == 1L) {
      samples
    } else {
      aperm(
        array(
          samples,
          c(dim(samples)[1L] %/% n.chains, n.chains, ncol(samples))
        ),
        c(2L, 1L, 3L)
      )
    }
    if (!is.null(dimnames(samples))) {
      dimnames(res) <- list(NULL, NULL, dimnames(samples)[[2L]])
    }
    res
  }
}

# One chain's sample-s column of a run channel; the engine drops the trailing
# chain margin when n.chains == 1. Shared by the ordinal and nbinom fit loops
# that read a bartcoreRun() channel column by column.
channelColumn <- function(channel, s, chain, n.chains) {
  if (n.chains == 1L) channel[, s] else channel[, s, chain]
}

# Reshapes a raw varcount channel to bart2's draws-first convention and, when
# the coded design carries column names, threads them onto the trailing
# predictor margin (2-D combined or 3-D by-chain). Pure reshape of an
# already-computed count channel; shared by the gaussian/binary, ordinal, and
# nbinom packagers.
nameVarcount <- function(raw, predictorNames, n.chains, combineChains) {
  varcount <- convertSamplesFromDbartsToBart(raw, n.chains, combineChains)
  if (!is.null(predictorNames) && !is.null(dim(varcount))) {
    dimnames(varcount) <- if (length(dim(varcount)) > 2L) {
      list(NULL, NULL, predictorNames)
    } else {
      list(NULL, predictorNames)
    }
  }
  varcount
}

packageBartResults <- function(
  fit,
  samples,
  burnInSigma,
  burnInK,
  combineChains,
  keepSampler
) {
  responseIsBinary <- fit$control@binary
  n.chains <- fit$control@n.chains

  yhat.train <- NULL
  yhat.train.mean <- NULL
  if (fit$control@keepTrainingFits) {
    yhat.train <- convertSamplesFromDbartsToBart(
      samples$train,
      n.chains,
      combineChains
    )
    if (!responseIsBinary) {
      yhat.train.mean <- apply(yhat.train, length(dim(yhat.train)), mean)
    }
  }

  yhat.test <- NULL
  yhat.test.mean <- NULL
  if (NROW(fit$data@x.test) > 0) {
    yhat.test <- convertSamplesFromDbartsToBart(
      samples$test,
      n.chains,
      combineChains
    )
    if (!responseIsBinary) {
      yhat.test.mean <- apply(yhat.test, length(dim(yhat.test)), mean)
    }
  }

  if (!responseIsBinary) {
    sigma <- convertSamplesFromDbartsToBart(
      samples$sigma,
      n.chains,
      combineChains
    )
  }

  # a multi-forest sampler's varcount channel carries a forest axis (raw
  # p x n.forests x n.samples (x n.chains), prognostic first), which the
  # single-forest reshape cannot see: it would fold the forest margin into the
  # draws margin when chains are combined and error inside aperm when they are
  # not. shapeMultinomialChannel is the K-margin reshape that already handles
  # all three cases, so the forest axis lands on the trailing margin exactly as
  # multinomial's category axis does, named in ENGINE vocabulary. The forest
  # count comes from the sampler rather than the array's rank, which is
  # ambiguous (a single-forest multi-chain channel is 3-D too); data@bases is
  # the same probe isBCFSampler uses.
  numForests <- if (!is.null(fit$data@bases)) length(fit$data@bases) else 1L
  varcount <- if (numForests > 1L) {
    shapeMultinomialChannel(
      samples$varcount,
      paste0("forest", seq_len(numForests)),
      n.chains,
      combineChains,
      leadNames = colnames(fit$data@x)
    )
  } else {
    nameVarcount(
      samples$varcount,
      colnames(fit$data@x),
      n.chains,
      combineChains
    )
  }

  # heteroscedastic variance surface s(x) = sqrt(s^2(x)), train and test, on the
  # original scale (docs/design/heteroscedastic.md); NULL for a homoscedastic fit
  s.train <- NULL
  s.test <- NULL
  if (!is.null(samples[["variance"]])) {
    s.train <- sqrt(convertSamplesFromDbartsToBart(
      samples$variance,
      n.chains,
      combineChains
    ))
    if (!is.null(samples[["varianceTest"]])) {
      s.test <- sqrt(convertSamplesFromDbartsToBart(
        samples$varianceTest,
        n.chains,
        combineChains
      ))
    }
  }

  varprobs <- NULL
  if (!is.null(samples[["varprobs"]])) {
    varprobs <- convertSamplesFromDbartsToBart(
      samples$varprobs,
      n.chains,
      combineChains
    )
    if (!is.null(colnames(fit$data@x)) && !is.null(dim(varprobs))) {
      dimnames(varprobs) <- if (length(dim(varprobs)) > 2L) {
        list(NULL, NULL, colnames(fit$data@x))
      } else {
        list(NULL, colnames(fit$data@x))
      }
    }
  }

  if (!is.null(burnInSigma)) {
    burnInSigma <- convertSamplesFromDbartsToBart(
      burnInSigma,
      n.chains,
      combineChains
    )
  }
  if (responseIsBinary) {
    result <- list(
      call = fit$control@call,
      family = fit$model@family,
      yhat.train = yhat.train,
      yhat.test = yhat.test,
      varcount = varcount,
      binaryOffset = fit$data@offset,
      # kept so residuals() works on binary fits (y - fitted); for gaussian
      # fits y rides the other branch
      y = fit$data@y
    )
  } else {
    result <- list(
      call = fit$control@call,
      family = fit$model@family,
      first.sigma = burnInSigma,
      sigma = sigma,
      sigest = fit$data@sigma,
      yhat.train = yhat.train,
      yhat.train.mean = yhat.train.mean,
      yhat.test = yhat.test,
      yhat.test.mean = yhat.test.mean,
      s.train = s.train,
      s.test = s.test,
      varcount = varcount,
      y = fit$data@y
    )
  }

  if (!is.null(varprobs)) {
    result$varprobs <- varprobs
  }

  # the per-observation censoring status an aft (survival) fit needs so its
  # log-likelihood can take the survival tail on censored rows; the survival
  # ingestion parks it on the control attribute, and it is NULL (so the
  # element stays absent) for every non-survival family
  result$status <- attr(fit$control, "bartcore.survival")

  # the discrete-time hazard marker (docs/design/survival.md section 4): the
  # ordered period grid, present only on hazard fits (parked on the control
  # attribute by dbarts()). $family reads the binary token on both, so every
  # link-keyed generic stays correct with zero change; survivalProbabilities
  # dispatches its hazard branch on THIS marker instead. NULL - and so the
  # element stays absent - for every non-hazard fit.
  result$periods <- attr(fit$control, "bartcore.hazard.periods")

  # the forest count the widened varcount channel's trailing margin carries,
  # the multi-forest analog of the multinomial packager's K; absent at one
  # forest, so every single-forest fit keeps exactly the elements it had.
  # fitSynopsis reads it to tell a forest margin from a chain margin, which
  # the packaged rank alone cannot do
  if (numForests > 1L) {
    result$n.forests <- numForests
  }

  if (keepSampler) {
    result$fit <- fit
  } else {
    result$n.chains <- n.chains
  }
  if (!is.null(samples[["k"]])) {
    result[["k"]] <- convertSamplesFromDbartsToBart(
      samples[["k"]],
      n.chains,
      combineChains
    )
    result[["first.k"]] <- convertSamplesFromDbartsToBart(
      burnInK,
      n.chains,
      combineChains
    )
  }

  class(result) <- "bart"
  invisible(result)
}

.kDefault <- quote(if (control@binary) quote(chi(1.5, 2.0)) else 2)

## Builds the quoted tree/node/resid prior calls that bart, bart2, and
## rbart_vi hand to dbarts. nodeK is the node prior's k argument exactly as
## it should enter the call - unevaluated for functions that redirect their
## matched call, evaluated for those that forward through do.call from
## internal frames - or NULL for no node prior. dart may be FALSE, TRUE, or
## a dbartsDartPrior spec; splitProbsName is the caller's argument spelling
## and splitProbsDefault its formal default.
buildSamplerPriors <- function(
  matchedCall,
  power,
  base,
  sigdf,
  sigquant,
  nodeK,
  priorScale = NA_real_,
  dart = FALSE,
  splitProbsName = "split.probs",
  splitProbsDefault = NULL
) {
  priorScale <- validateNamedScale(priorScale, "prior.scale")
  if (inherits(dart, "dbartsDartPrior")) {
    # a full spec overrides the power/base arguments with its own
    tree.prior <- dart
  } else if (isTRUE(dart)) {
    if (splitProbsName %in% names(matchedCall)) {
      stop(
        "'",
        splitProbsName,
        "' cannot be combined with 'dart': a DART prior samples its split probabilities"
      )
    }
    tree.prior <- quote(dart(power, base))
    tree.prior[[2L]] <- power
    tree.prior[[3L]] <- base
  } else if (!isFALSE(dart)) {
    stop("'dart' must be TRUE, FALSE, or a prior created by dbartsPriors$dart")
  } else {
    tree.prior <- quote(cgm(power, base, split.probs))
    tree.prior[[2L]] <- power
    tree.prior[[3L]] <- base
    if (splitProbsName %in% names(matchedCall)) {
      tree.prior[[4L]] <- matchedCall[[splitProbsName]]
    } else {
      tree.prior[[4L]] <- splitProbsDefault
    }
  }

  # a named prior.scale needs a node prior to ride even when k is left to the
  # family default, so it builds one; the k slot then drops out of the call and
  # dbarts() resolves k exactly as it would have with no node prior at all
  if (!is.null(nodeK) || !is.na(priorScale)) {
    node.prior <- quote(normal(k))
    node.prior[[2L]] <- nodeK
    if (!is.na(priorScale)) {
      node.prior[["scale"]] <- priorScale
    }
  } else {
    node.prior <- NULL
  }

  resid.prior <- quote(chisq(sigdf, sigquant))
  resid.prior[[2L]] <- sigdf
  resid.prior[[3L]] <- sigquant

  list(
    tree.prior = tree.prior,
    node.prior = node.prior,
    resid.prior = resid.prior
  )
}

# Builds the dbarts() host-sampler call shared by bart2's standard path and its
# alternate-family arcs (multinomial/ordinal/nbinom): redirect the matched call
# to dbarts(), install the resolved control, drop n.samples (the sampler is
# driven by run() rather than sampled at construction), and thread the prebuilt
# tree/node/resid priors. 'family', when supplied, overrides the forwarded
# family token - NULL removes it (the multinomial host, which dbarts() does not
# know), a string sets it (ordinal/nbinom); left missing it keeps the user's
# forwarded family (the normal path). 'sigest', when supplied, sets the
# sampler's sigma; left missing it is not touched (ordinal/nbinom read no
# sigma). The multinomial data override (the placeholder response) is applied
# by the caller afterward.
buildHostSamplerCall <- function(
  matchedCall,
  control,
  priors,
  family,
  sigest
) {
  samplerCall <- redirectCall(matchedCall, dbarts::dbarts)
  samplerCall$control <- control
  samplerCall$n.samples <- NULL
  samplerCall$tree.prior <- priors$tree.prior
  samplerCall$node.prior <- priors$node.prior
  samplerCall$resid.prior <- priors$resid.prior
  if (!missing(family)) {
    samplerCall$family <- family
  }
  if (!missing(sigest)) {
    samplerCall$sigma <- as.numeric(sigest)
  }
  samplerCall
}

# The two-phase burn-in/sample split shared by bart and bart2's standard path.
# When burn-in is requested the test set is dropped and keepTrainingFits/verbose
# forced FALSE for the burn run (a draw-neutral speedup - test fits do not feed
# the MCMC), then restored for the kept-sample run, with keepTrees re-enabled
# after burn when requested. The sequence of sampler$run() calls - and thus the
# draw stream - is identical to the inline form it replaces. Returns the kept
# samples plus the burn-in sigma/k channels.
runWithBurnIn <- function(sampler, control, keepTrees) {
  burnInSigma <- NULL
  burnInK <- NULL
  if (control@n.burn > 0L) {
    oldX.test <- sampler$data@x.test
    oldOffset.test <- sampler$data@offset.test

    oldKeepTrainingFits <- control@keepTrainingFits
    oldVerbose <- control@verbose

    if (length(oldX.test) > 0L) {
      sampler$setTestPredictorAndOffset(NULL, NULL)
    }
    control@keepTrainingFits <- FALSE
    control@verbose <- FALSE
    sampler$setControl(control)

    samples <- sampler$run(0L, control@n.burn, updateState = FALSE)
    if (!is.null(samples$sigma)) {
      burnInSigma <- samples$sigma
    }
    if (!is.null(samples[["k"]])) {
      burnInK <- samples[["k"]]
    }

    if (length(oldX.test) > 0L) {
      sampler$setTestPredictorAndOffset(oldX.test, oldOffset.test)
    }
    control@keepTrainingFits <- oldKeepTrainingFits
    control@verbose <- oldVerbose
    if (keepTrees == TRUE) {
      control@keepTrees <- TRUE
    }
    sampler$setControl(control)

    samples <- sampler$run(0L, control@n.samples, updateState = FALSE)
  } else {
    samples <- sampler$run(updateState = FALSE)
  }
  list(samples = samples, burnInSigma = burnInSigma, burnInK = burnInK)
}

# The alternate-family bart2 arcs (multinomial/ordinal/nbinom/hurdle.lognormal)
# all refuse samplerOnly, warm.start/n.grow.sweeps, and keepTrainingFits = FALSE,
# differing only in the family name and the keepTrainingFits reason (which
# completes "requires keepTrainingFits = TRUE (the default): ").
checkFamilyUnsupportedArgs <- function(
  family,
  samplerOnly,
  warm.start,
  n.grow.sweeps,
  control,
  reason
) {
  if (isTRUE(samplerOnly)) {
    stop("family = \"", family, "\" does not support 'samplerOnly'")
  }
  grownSweeps <- as.integer(n.grow.sweeps)[1L]
  if (!is.null(warm.start) || (!is.na(grownSweeps) && grownSweeps > 0L)) {
    stop(
      "family = \"",
      family,
      "\" does not support 'warm.start' or ",
      "'n.grow.sweeps'"
    )
  }
  if (!control@keepTrainingFits) {
    stop(
      "family = \"",
      family,
      "\" requires keepTrainingFits = TRUE (the default): ",
      reason
    )
  }
}

bart2 <- function(
  formula,
  data,
  test,
  subset,
  weights,
  offset,
  offset.test = offset,
  sigest = NA_real_,
  sigdf = 3.0,
  sigquant = 0.90,
  k = NULL,
  prior.scale = NA_real_,
  power = 2.0,
  base = 0.95,
  split.probs = 1 / num.vars,
  dart = FALSE,
  n.trees = 75L,
  n.samples = 500L,
  n.burn = 500L,
  n.chains = 4L,
  n.threads = min(dbarts::guessNumCores(), n.chains),
  combineChains = TRUE,
  n.cuts = 100L,
  useQuantiles = FALSE,
  n.thin = 1L,
  keepTrainingFits = TRUE,
  printEvery = 100L,
  printCutoffs = 0L,
  verbose = TRUE,
  keepTrees = FALSE,
  keepCall = TRUE,
  samplerOnly = FALSE,
  seed = NA_integer_,
  proposal.probs = c(birth_death = 0.5, swap = 0.1, change = 0.4, birth = 0.5),
  monotone = NULL,
  interactions = NULL,
  blocks = NULL,
  variance = NULL,
  n.trees.variance = 40L,
  power.variance = NULL,
  base.variance = NULL,
  keepSampler = keepTrees,
  warm.start = NULL,
  n.grow.sweeps = 0L,
  factors = c("categorical", "indicators"),
  family = c(
    "auto",
    "gaussian",
    "probit",
    "logistic",
    "aft",
    "multinomial",
    "ordinal",
    "nbinom",
    "hazard",
    "hazard.probit",
    "hazard.logistic",
    "hurdle.lognormal",
    "twopart"
  ),
  missing = c("incorporate", "error"),
  resid.dist = gaussian,
  dispersion = NA_real_,
  breaks = NULL,
  max.rows = 1e7,
  ...
) {
  matchedCall <- match.call()
  callingEnv <- parent.frame()
  family <- match.arg(family)
  # hurdle.lognormal / twopart (docs/design/hurdle.md): resolve the alias
  # before anything reads 'family', so the token prints and dispatches as
  # "hurdle.lognormal" regardless of which spelling was requested
  if (identical(family, "twopart")) {
    family <- "hurdle.lognormal"
  }

  # family = "auto" with a 3+-level UNORDERED factor/character response is
  # multinomial (docs/design/multinomial.md); a 3+-level ORDERED factor is
  # ordinal (docs/design/ordinal.md, the disjoint is.ordered() key). Route
  # either into its branch below and announce. 2-level and numeric responses
  # stay on the standard path, where dbarts() resolves probit/gaussian (and
  # announces a factor probit itself).
  if (family == "auto") {
    dataMissing <- missing(data)
    autoData <- if (dataMissing) NULL else data
    autoMultinomial <- detectAutoMultinomial(
      formula,
      autoData,
      dataMissing,
      callingEnv
    )
    if (!is.null(autoMultinomial)) {
      family <- "multinomial"
      announceAutoFamily(
        autoMultinomial$type,
        autoMultinomial$n.levels,
        family
      )
    } else {
      autoOrdinal <- detectAutoOrdinal(
        formula,
        autoData,
        dataMissing,
        callingEnv
      )
      if (!is.null(autoOrdinal)) {
        family <- "ordinal"
        announceAutoFamily(autoOrdinal$type, autoOrdinal$n.levels, family)
      }
    }
  }

  argNames <- names(matchedCall)[-1L]
  unknownArgs <- argNames %not_in%
    names(formals(dbarts::bart2)) &
    argNames %not_in% names(formals(dbarts::dbartsControl))
  if (any(unknownArgs)) {
    stop(
      "unknown arguments: '",
      paste0(argNames[unknownArgs], collapse = "', '"),
      "'"
    )
  }

  # d1/D6 (docs/plans/bart2-argument-consolidation.md 3.d): factors/missing/
  # proposal.probs are forwarded formal defaults - redirectCall only carries
  # a name into the host dbarts() call when the caller supplied it, so an
  # unsupplied one used to silently take dbarts()'s own default rather than
  # the token/value this signature advertises. Resolve here, in bart2's own
  # frame, and stamp the resolved value onto matchedCall unconditionally so
  # every call built from it below forwards it explicitly. T-B keeps the
  # default texts identical to dbarts()'s, so this is draw-neutral: the
  # resolved value is exactly what dbarts() would already have chosen.
  factors <- match.arg(factors)
  missing <- match.arg(missing)
  matchedCall$factors <- factors
  matchedCall$missing <- missing
  matchedCall$proposal.probs <- proposal.probs

  controlCall <- redirectCall(matchedCall, dbarts::dbartsControl)
  missingDefaultArgs <- names(formals(dbarts::bart2))[
    names(formals(dbarts::bart2)) %in%
      names(formals(dbarts::dbartsControl)) &
      names(formals(dbarts::bart2)) %not_in% names(matchedCall)
  ]
  if (length(missingDefaultArgs) > 0L) {
    currentEnv <- sys.frame(sys.nframe())
    controlCall[missingDefaultArgs] <- lapply(
      formals(dbarts::bart2)[missingDefaultArgs],
      eval,
      envir = currentEnv
    )
  }
  if (!is.na(seed)) {
    controlCall[["rngSeed"]] <- seed
  }
  control <- eval(controlCall, envir = callingEnv)

  control@call <- if (keepCall) matchedCall else call("NULL")
  control@n.burn <- control@n.burn %/% control@n.thin
  control@n.samples <- control@n.samples %/% control@n.thin
  control@printEvery <- control@printEvery %/% control@n.thin

  # a zero (or thinned-to-zero) sample count would otherwise fault deeper, in
  # the empty-array reshape (dim(X) has no positive length); the multinomial
  # branch below keeps its own family-named check
  if (family != "multinomial" && isTRUE(control@n.samples <= 0L)) {
    stop("'n.samples' must be a positive integer")
  }

  # multinomial (docs/design/multinomial.md): a K-forest softmax model over
  # a factor response or an n x K count matrix, validated and dispatched
  # here rather than threaded through the rest of bart2 - it bypasses the
  # standard single-forest dbarts()/run() path entirely
  # (bartcoreMultinomialSampler/bartcoreMultinomialCountSampler each build a
  # separate K-forest engine; see benchmarks/R/multinomial-equivalence.R for
  # the exact call sequences bart2Multinomial/bart2MultinomialCounts
  # mirror). Every refusal below names the limitation rather than silently
  # reshaping around it.
  if (family == "multinomial") {
    if (!missing(weights)) {
      stop(
        "family = \"multinomial\" does not support 'weights': a non-integer ",
        "weight has no exact augmentation sampler, and an integer weight is ",
        "already expressible as row-wise count replication in the response"
      )
    }
    # the n x K category offset (docs/design/multinomial.md "The surface"):
    # a matrix enters the raw fits before the softmax and is threaded to
    # bartcoreMultinomialSampler/bartcoreMultinomialCountSampler's own offset
    # argument, never to the host dbarts() call below (matchedCall$offset is
    # cleared so redirectCall does not forward it there, where it would be
    # read as a flat per-row offset). A flat vector stays refused - it points
    # exactly along the softmax's null direction (a common per-observation
    # shift) and is identically inert, exactly as the host data object's own
    # flat offset is (parseMultinomialData, refused separately, below).
    multinomialOffset <- NULL
    if (!missing(offset)) {
      if (!is.matrix(offset)) {
        stop(
          "family = \"multinomial\" does not support a flat 'offset': a ",
          "common per-observation shift is the softmax's own null direction ",
          "and is identically inert; pass an n x K numeric matrix instead"
        )
      }
      if (!is.numeric(offset)) {
        stop("family = \"multinomial\" 'offset' must be a numeric n x K matrix")
      }
      multinomialOffset <- offset
      matchedCall$offset <- NULL
    }
    # offset is train-side only: an explicit offset.test is caught HERE,
    # before it would otherwise fall through buildHostSamplerCall's
    # redirectCall into the host dbarts() call and be refused several layers
    # down by parseMultinomialData, several steps from where the caller wrote
    # it. yhat.test is always computed WITHOUT any category offset (see
    # docs/design/multinomial.md "The surface"), so a caller wanting one
    # needs the dbarts:::-only channel this message names.
    if (!missing(offset.test)) {
      stop(
        "family = \"multinomial\" does not support 'offset.test'; yhat.test ",
        "is always computed without any category offset. A category test ",
        "offset is an internal-channel capability only ",
        "(dbarts:::bartcoreSetCategoryTestOffset, or the internal creators' ",
        "own offset.test argument), not reachable from bart2"
      )
    }
    if (!missing(subset)) {
      stop(
        "family = \"multinomial\" does not support 'subset'"
      )
    }
    checkFamilyUnsupportedArgs(
      "multinomial",
      samplerOnly,
      warm.start,
      n.grow.sweeps,
      control,
      "there is no test surface to fall back on"
    )
    warnFamilyGatedArgs(argNames, "multinomial")
    if (control@n.samples <= 0L) {
      stop("family = \"multinomial\" requires a positive 'n.samples'")
    }
    if (inherits(formula, "dbartsData")) {
      stop(
        "family = \"multinomial\" does not support a pre-built dbartsData ",
        "object; use the formula interface or the matrix ",
        "interface: bart2(x.train, y.train, family = \"multinomial\")"
      )
    }
    if (is.formula(formula)) {
      if (missing(data) || is.null(data)) {
        stop(
          "family = \"multinomial\" requires 'data' when 'formula' is a formula"
        )
      }
      extracted <- extractMultinomialFormulaData(
        matchedCall,
        callingEnv,
        formula,
        data
      )
      y <- extracted$y
      # from here on, y drives the same factor-vs-count-matrix dispatch as
      # the matrix interface; the host sampler (bart2Multinomial/
      # bart2MultinomialCounts, below) is built from the term-labeled
      # predictor DATA FRAME, not the formula - dbartsData's own
      # data-frame-as-x.train branch codes it (makeModelMatrix, chosen by
      # 'factors' exactly as it would for any other family's formula fit),
      # which is what threads term.labels/factor.levels onto sampler$data@x
      # and makes predict's data.frame newdata coding (validateXTest) work
      matchedCall$formula <- extracted$x
    } else {
      if (missing(data) || is.null(data)) {
        stop(
          "family = \"multinomial\" requires a factor response as the ",
          "second argument (the matrix-interface y.train)"
        )
      }
      y <- data
    }
    # the count-matrix response form (docs/design/multinomial.md): an n x K
    # matrix of nonnegative integer trial counts, beside the factor path
    # below. Any numeric matrix is treated as an attempted count response, so
    # a too-narrow one gets its own error rather than falling through to the
    # factor message.
    if (is.matrix(y) && is.numeric(y)) {
      if (ncol(y) < 2L) {
        stop(
          "family = \"multinomial\" count response matrix must have at ",
          "least 2 columns (K >= 2 categories)"
        )
      }
      if (anyNA(y)) {
        stop(
          "family = \"multinomial\" does not support missing response values"
        )
      }
      if (any(y < 0)) {
        stop("family = \"multinomial\" count response must be nonnegative")
      }
      if (any(y != round(y))) {
        stop("family = \"multinomial\" count response must be whole numbers")
      }
      if (any(rowSums(y) < 1)) {
        stop(
          "family = \"multinomial\" count response requires every row to ",
          "have at least one trial (row sum >= 1)"
        )
      }

      levels <- colnames(y)
      if (is.null(levels)) {
        levels <- as.character(seq_len(ncol(y)))
      }

      return(bart2MultinomialCounts(
        matchedCall,
        callingEnv,
        control,
        y,
        levels,
        power,
        base,
        sigdf,
        sigquant,
        sigest,
        dart,
        combineChains,
        offset = multinomialOffset,
        prior.scale = prior.scale,
        keepSampler = keepSampler
      ))
    }
    if (is.character(y)) {
      y <- factor(y)
    }
    if (!is.factor(y)) {
      stop(
        "family = \"multinomial\" requires a factor (or character) response, ",
        "or an n x K count matrix"
      )
    }
    if (anyNA(y)) {
      stop(
        "family = \"multinomial\" does not support missing response values"
      )
    }
    if (anyNA(levels(y))) {
      stop("family = \"multinomial\" response levels must not include NA")
    }
    if (nlevels(y) < 2L) {
      stop(
        "family = \"multinomial\" requires a response with at least 2 levels"
      )
    }

    return(bart2Multinomial(
      matchedCall,
      callingEnv,
      control,
      y,
      power,
      base,
      sigdf,
      sigquant,
      sigest,
      dart,
      combineChains,
      offset = multinomialOffset,
      prior.scale = prior.scale,
      keepSampler = keepSampler
    ))
  }

  # ordinal (cumulative probit, docs/design/ordinal.md): a single-forest fixed-
  # scale latent model whose n x K category probabilities and per-draw cutpoints
  # are synthesized in R (the engine reports only the latent eta), dispatched
  # here rather than threaded through the standard single-forest path so its
  # K-widened fit object (class "bartOrdinal", never "bart") stays distinct.
  if (family == "ordinal") {
    checkFamilyUnsupportedArgs(
      "ordinal",
      samplerOnly,
      warm.start,
      n.grow.sweeps,
      control,
      "the category probabilities are built from the training latent fits"
    )
    warnFamilyGatedArgs(argNames, "ordinal")
    return(bart2Ordinal(
      matchedCall,
      callingEnv,
      control,
      power,
      base,
      sigdf,
      sigquant,
      dart,
      combineChains,
      prior.scale = prior.scale,
      keepSampler = keepSampler
    ))
  }

  # negative-binomial counts (docs/design/negative-binomial.md): a single-forest
  # fixed-scale latent model whose mean counts mu = r exp(f + o) and per-draw
  # dispersion r are synthesized in R (the engine reports the log-odds latent
  # psi; r lives only in its state block), dispatched here so its bartNegbin fit
  # object (never "bart") stays distinct. family = "nbinom" is always explicit -
  # a count response has no unambiguous class to auto-detect.
  if (family == "nbinom") {
    checkFamilyUnsupportedArgs(
      "nbinom",
      samplerOnly,
      warm.start,
      n.grow.sweeps,
      control,
      "the mean counts are built from the training latent fits"
    )
    warnFamilyGatedArgs(argNames, "nbinom")
    return(bart2Negbin(
      matchedCall,
      callingEnv,
      control,
      power,
      base,
      sigdf,
      sigquant,
      dart,
      combineChains,
      prior.scale = prior.scale,
      keepSampler = keepSampler
    ))
  }

  # hurdle.lognormal / twopart (docs/design/hurdle.md): a semicontinuous two-
  # part fit built from an occupancy probit on 1{y > 0} (all n) and a gaussian
  # on log(y) restricted to the y > 0 subset, glued at report time. Dispatched
  # here - not inside dbarts(), which returns a single sampler and cannot
  # express the two-sampler composition - so its bartHurdle fit object stays
  # distinct. family is always explicit: a semicontinuous response has no
  # unambiguous auto class.
  if (family == "hurdle.lognormal") {
    checkFamilyUnsupportedArgs(
      "hurdle.lognormal",
      samplerOnly,
      warm.start,
      n.grow.sweeps,
      control,
      "the combined mean is built from the two training fits"
    )
    warnFamilyGatedArgs(argNames, "hurdle.lognormal")
    if (!missing(weights)) {
      stop("family = \"hurdle.lognormal\" does not support 'weights'")
    }
    if (!missing(subset)) {
      stop("family = \"hurdle.lognormal\" does not support 'subset'")
    }
    if (!missing(offset) || !missing(offset.test)) {
      stop(
        "family = \"hurdle.lognormal\" does not support 'offset'/'offset.test'"
      )
    }
    if (!missing(test)) {
      stop(
        "family = \"hurdle.lognormal\" does not take a 'test' set; the ",
        "positive-part fit is given the full training x as its own x.test"
      )
    }
    return(bart2Hurdle(matchedCall, callingEnv, control, formula, data, seed))
  }

  keepSampler <- keepSampler || control@keepTrees

  if (control@n.burn == 0L && keepTrees == TRUE) {
    control@keepTrees <- TRUE
  }
  if (control@n.burn > 0L) {
    control@keepTrees <- FALSE
  }

  # k enters unevaluated: bart2 redirects its matched call, so the stored
  # symbol resolves in the caller's frame
  priors <- buildSamplerPriors(
    matchedCall,
    power,
    base,
    sigdf,
    sigquant,
    nodeK = matchedCall[["k"]],
    priorScale = prior.scale,
    dart = dart,
    splitProbsDefault = formals(dbarts::bart2)[["split.probs"]]
  )

  samplerCall <- buildHostSamplerCall(
    matchedCall,
    control,
    priors,
    sigest = sigest
  )

  sampler <- eval(samplerCall, envir = callingEnv)
  # the hazard tokens remap to their underlying binary link before the model
  # object is built (R/dbarts.R), so sampler$model@family can no longer say
  # "hazard" - keep the caller's own explicit token for that case, and let
  # sampler$model@family resolve everything family = "auto" left open
  # otherwise (numeric/binary responses are not settled until dbarts() runs)
  gatedFamily <- if (
    family %in% c("hazard", "hazard.probit", "hazard.logistic")
  ) {
    family
  } else {
    sampler$model@family
  }
  warnFamilyGatedArgs(argNames, gatedFamily)
  if (isTRUE(samplerOnly)) {
    return(sampler)
  }

  control <- sampler$control

  # the initial forest: a warm-start donor, a grow-from-root warm start, or a
  # draw from the prior (the default, byte-identical to before this argument)
  n.grow.sweeps <- as.integer(n.grow.sweeps)[1L]
  if (is.na(n.grow.sweeps) || n.grow.sweeps < 0L) {
    stop("'n.grow.sweeps' must be a non-negative integer")
  }
  if (!is.null(warm.start) && n.grow.sweeps > 0L) {
    stop(
      "'warm.start' and 'n.grow.sweeps' both request an initialization; ",
      "supply at most one"
    )
  }

  if (!is.null(warm.start)) {
    sampler$installTrees(warm.start)
  } else if (n.grow.sweeps > 0L) {
    sampler$growFromRoot(n.grow.sweeps, updateState = FALSE)
  } else {
    sampler$sampleTreesFromPrior(updateState = FALSE)
  }

  burn <- runWithBurnIn(sampler, control, keepTrees)
  samples <- burn$samples
  burnInSigma <- burn$burnInSigma
  burnInK <- burn$burnInK

  result <- packageBartResults(
    sampler,
    samples,
    burnInSigma,
    burnInK,
    combineChains,
    keepSampler
  )
  # needed to extract ppd
  if (!is.null(sampler$data@weights) && length(sampler$data@weights) > 0L) {
    result$weights <- sampler$data@weights
    if (
      !is.null(sampler$data@weights.test) &&
        length(sampler$data@weights.test) > 0L
    ) {
      result$weights.test <- sampler$data@weights.test
    }
  }

  result
}

# Formula ingestion for bart2's family = "multinomial" branch, above. y is
# pulled with model.response(modelFrame) - NO type coercion - so a factor
# LHS keeps its levels and a cbind(c1, ..., cK) ~ x LHS keeps its column
# names; dbartsData's own formula path cannot be reused for this because it
# calls model.response(modelFrame, "numeric"), which discards both. The
# right-hand side is returned as the term-labeled predictor DATA FRAME
# (model.frame -> terms, the data.R formula recipe short of the final
# makeModelMatrix call), not yet coded: the caller installs it as the
# matched call's 'formula', and the host sampler build below
# (bart2Multinomial/bart2MultinomialCounts) reuses dbartsData's own
# data-frame-as-x.train branch to code it - the same makeModelMatrix
# call, choosing the categorical or indicators builder the same way via
# 'factors' - so it carries the same term.labels/factor.levels/drop
# attributes onto sampler$data@x that let predict.bartMultinomial's
# validateXTest code a data.frame newdata. (A matrix coded here directly
# would not work: 'factors = "categorical"', the default, produces a
# dbartsMixedMatrix container that only dbartsData's data-frame branch
# knows how to route to a fresh sampler build - the matrix/vector top-level
# dispatch does not recognize it.)
extractMultinomialFormulaData <- function(
  matchedCall,
  callingEnv,
  formula,
  data
) {
  if (!is.data.frame(data) && !is.list(data) && !is.environment(data)) {
    stop(
      "for formula/data specification, 'data' must be a data frame, list, ",
      "or environment"
    )
  }

  # a sparseFactor predictor dies inside model.frame with a bare S4 type
  # error; refuse it explicitly first, as the standard formula path does
  if (is.list(data) || is.environment(data)) {
    for (variableName in intersect(all.vars(formula), names(data))) {
      if (methods::is(data[[variableName]], "sparseFactor")) {
        stop(
          "sparse categorical predictors must be supplied through the x/y ",
          "interface; '",
          variableName,
          "' is a sparseFactor"
        )
      }
    }
  }

  modelFrameCall <- matchedCall
  modelFrameCall <- modelFrameCall[c(
    1L,
    match(c("formula", "data"), names(modelFrameCall), nomatch = 0L)
  )]
  modelFrameCall$drop.unused.levels <- FALSE
  modelFrameCall$na.action <- stats::na.pass
  modelFrameCall[[1L]] <- quote(stats::model.frame)

  modelFrame <- eval(modelFrameCall, callingEnv)
  if (NROW(modelFrame) == 0) {
    stop("cannot construct model matrices from formula")
  }

  y <- model.response(modelFrame)
  if (is.null(y)) {
    stop("family = \"multinomial\" requires a response in 'formula'")
  }

  modelTerms <- terms(modelFrame)
  if (is.empty.model(modelTerms)) {
    stop("predictors must be specified for regression tree analysis")
  }
  termLabels <- attr(modelTerms, "term.labels")
  badLabels <- grepl("`.* .*`", termLabels)
  if (sum(badLabels) > 0) {
    termLabels[badLabels] <- gsub("^`(.*)`$", "\\1", termLabels[badLabels])
  }

  list(y = y, x = modelFrame[termLabels])
}

# family = "auto" peek for bart2: a 3+-level UNORDERED factor or character
# response selects multinomial (an ORDERED factor is ordinal, split out here and
# caught by detectAutoOrdinal instead - the disjoint is.ordered() key,
# docs/design/ordinal.md). Returns classifyResponse's descriptor for that case,
# NULL for anything that stays on the standard single-forest path (numeric,
# 2-level, logical, ordered, a count matrix, or a pre-built dbartsData). Type
# detection only - the multinomial branch re-extracts and validates the
# response, so this evaluates the formula LHS directly rather than building a
# model frame.
detectAutoResponse <- function(formula, data, dataIsMissing, callingEnv) {
  response <- if (is.formula(formula)) {
    if (length(formula) != 3L || dataIsMissing || is.null(data)) {
      return(NULL)
    }
    tryCatch(
      eval(formula[[2L]], envir = data, enclos = callingEnv),
      error = function(e) NULL
    )
  } else if (
    inherits(formula, "dbartsData") || inherits(formula, "dgCMatrix")
  ) {
    return(NULL)
  } else if (dataIsMissing) {
    return(NULL)
  } else {
    data
  }
  classifyResponse(response)
}

detectAutoMultinomial <- function(formula, data, dataIsMissing, callingEnv) {
  info <- detectAutoResponse(formula, data, dataIsMissing, callingEnv)
  if (
    !is.null(info) &&
      info$type %in% c("factor", "character") &&
      !is.na(info$n.levels) &&
      info$n.levels >= 3L
  ) {
    info
  } else {
    NULL
  }
}

# The ordinal analog of detectAutoMultinomial: a 3+-level ORDERED factor
# response selects family = "ordinal" (docs/design/ordinal.md section 5). Only
# is.ordered() responses match, so ordinal and unordered multinomial never
# overlap; everything else stays on the standard path.
detectAutoOrdinal <- function(formula, data, dataIsMissing, callingEnv) {
  info <- detectAutoResponse(formula, data, dataIsMissing, callingEnv)
  if (
    !is.null(info) &&
      identical(info$type, "ordered factor") &&
      !is.na(info$n.levels) &&
      info$n.levels >= 3L
  ) {
    info
  } else {
    NULL
  }
}

# The multinomial (softmax) fit path (docs/design/multinomial.md), reached
# from bart2's family = "multinomial" branch after ingestion validation. y is
# the validated factor response; labels (0-based, as the engine wants them)
# and K = nlevels(y) follow from it. The host sampler is built through
# bart2's usual tree.prior/node.prior/resid.prior/control machinery so that
# n.trees, n.chains, the tree prior, and k (the only host-model quantity the
# multinomial engine reads; see ?bart2) match what a non-multinomial bart2
# call would build - but its response is the label vector, never y itself,
# and its resolved family (whatever "auto" picks for that placeholder) is
# irrelevant and never surfaced. bartcoreMultinomialSampler then builds a
# separate K-forest engine; benchmarks/R/multinomial-equivalence.R exercises
# the identical sequence (host creation, then bartcoreMultinomialSampler,
# then one bartcoreRun call) for the bitwise reproduction gate in
# test-multinomial-surface.R. No warm start and no two-phase burn-in/sample
# split: both are skipped so the RNG stream matches that internal pattern
# exactly for a given seed. offset, when non-NULL, is the validated n x K
# category offset (bart2's own matrix-only surface, above); it is threaded to
# bartcoreMultinomialSampler's own offset argument, never to the host sampler
# call, whose 'offset' matchedCall has already had cleared.
bart2Multinomial <- function(
  matchedCall,
  callingEnv,
  control,
  y,
  power,
  base,
  sigdf,
  sigquant,
  sigest,
  dart,
  combineChains,
  offset = NULL,
  prior.scale = NA_real_,
  keepSampler = FALSE
) {
  K <- nlevels(y)
  labels <- as.integer(y) - 1L

  priors <- buildSamplerPriors(
    matchedCall,
    power,
    base,
    sigdf,
    sigquant,
    nodeK = matchedCall[["k"]],
    priorScale = prior.scale,
    dart = dart,
    splitProbsDefault = formals(dbarts::bart2)[["split.probs"]]
  )

  # the host sampler: identical machinery to bart2's normal path, but its
  # response is the integer label vector (as doubles - a placeholder the
  # multinomial engine never reads, since the category labels ride
  # separately) and no 'family' is forwarded (dbarts() does not know
  # "multinomial"; whatever "auto" resolves the placeholder to is immaterial)
  samplerCall <- buildHostSamplerCall(
    matchedCall,
    control,
    priors,
    family = NULL,
    sigest = sigest
  )
  samplerCall$data <- as.double(labels)

  sampler <- eval(samplerCall, envir = callingEnv)

  bc <- bartcoreMultinomialSampler(sampler, labels, K = K, offset = offset)
  # from here the host owns no model: bc does. Mark it, so a mutation through
  # the fit's $fit errors instead of silently changing a sampler nothing reads
  sampler$hostFor <- "bart2(family = \"multinomial\")"
  samples <- bartcoreRun(bc, control@n.burn, control@n.samples)

  result <- packageMultinomialResults(
    control,
    y,
    levels(y),
    K,
    samples,
    combineChains,
    predictorNames = colnames(sampler$data@x)
  )
  # keepTrees retains the handles predict.bartMultinomial replays through: bc
  # holds every one of the K forests' saved trees (the sampling sweeps wrote
  # them regardless), and the host sampler's coded design (sampler@data@x) codes
  # newdata to the training columns. Without keepTrees neither survives the call.
  # keepSampler retains $fit on its own, independent of keepTrees, so a caller
  # can see the host sampler without paying for the K forests' saved trees.
  if (control@keepTrees) {
    result$bc <- bc
  }
  if (control@keepTrees || keepSampler) {
    result$fit <- sampler
  }
  result
}

# The grouped-count analog of bart2Multinomial: y is the validated n x K
# count matrix (bart2's count-matrix branch above) and levels are the
# resolved category names (colnames(y), or the index fallback). Mirrors
# bart2Multinomial's host-sampler construction exactly - same priors,
# same redirectCall/override sequence - substituting
# bartcoreMultinomialCountSampler for bartcoreMultinomialSampler and a
# placeholder response column of y (as doubles, never read by the engine)
# for the label vector. A one-hot y with every row sum 1 is therefore the
# same draw stream as bart2Multinomial on the equivalent factor (the
# single-trial reduction; see benchmarks/R/multinomial-equivalence.R's
# k3counts scenario). offset, as for bart2Multinomial, is the validated n x K
# category offset threaded to bartcoreMultinomialCountSampler's own argument.
bart2MultinomialCounts <- function(
  matchedCall,
  callingEnv,
  control,
  y,
  levels,
  power,
  base,
  sigdf,
  sigquant,
  sigest,
  dart,
  combineChains,
  offset = NULL,
  prior.scale = NA_real_,
  keepSampler = FALSE
) {
  K <- ncol(y)

  priors <- buildSamplerPriors(
    matchedCall,
    power,
    base,
    sigdf,
    sigquant,
    nodeK = matchedCall[["k"]],
    priorScale = prior.scale,
    dart = dart,
    splitProbsDefault = formals(dbarts::bart2)[["split.probs"]]
  )

  samplerCall <- buildHostSamplerCall(
    matchedCall,
    control,
    priors,
    family = NULL,
    sigest = sigest
  )
  samplerCall$data <- as.double(y[, 1L])

  sampler <- eval(samplerCall, envir = callingEnv)

  bc <- bartcoreMultinomialCountSampler(sampler, y, K = K, offset = offset)
  # the host owns no model past this point; see bart2Multinomial
  sampler$hostFor <- "bart2(family = \"multinomial\")"
  samples <- bartcoreRun(bc, control@n.burn, control@n.samples)

  result <- packageMultinomialResults(
    control,
    y,
    levels,
    K,
    samples,
    combineChains,
    predictorNames = colnames(sampler$data@x)
  )
  if (control@keepTrees) {
    result$bc <- bc
  }
  if (control@keepTrees || keepSampler) {
    result$fit <- sampler
  }
  result
}

# Reshapes one K-carrying channel (lead x K x n.samples (x n.chains): the run's
# train/test channel, predict's replay, and the per-category varcount channel
# all share this raw shape - lead is the observations for the fits, the
# predictors for varcount) into the package's draws-first convention:
# (n.chains x) n.samples x lead x K, combining the chain margin as bart2's other
# families do, with levels named on the trailing K margin so every K-shaped
# output threads them. leadNames, when supplied, names the lead margin (the
# predictor names varcount threads); the fits leave it NULL.
shapeMultinomialChannel <- function(
  raw,
  levels,
  n.chains,
  combineChains,
  leadNames = NULL
) {
  out <- if (n.chains == 1L) {
    aperm(raw, c(3L, 1L, 2L))
  } else if (combineChains) {
    a <- aperm(raw, c(3L, 4L, 1L, 2L))
    d <- dim(a)
    dim(a) <- c(d[1L] * d[2L], d[3L], d[4L])
    a
  } else {
    aperm(raw, c(4L, 3L, 1L, 2L))
  }
  numDims <- length(dim(out))
  dn <- rep(list(NULL), numDims)
  dn[[numDims]] <- levels
  # assigning NULL through [[<- would drop the slot, so guard the lead margin
  if (!is.null(leadNames)) {
    dn[[numDims - 1L]] <- leadNames
  }
  dimnames(out) <- dn
  out
}

# Reshapes one bartcoreRun() result into a bart2(family = "multinomial") fit.
# samples$train is n.obs x K x n.samples (x n.chains); the K-carrying softmax
# probabilities are already the identified quantity the engine reports (the
# run's train channel is the identified deliverable; see
# docs/design/multinomial.md), so no probabilityFromLatents-style transform
# applies here, unlike the binary families. Reshaped to the package's
# draws-first convention (n.chains x) n.samples x n.obs x K - matching
# combineChains as bart2's other families do - with the resolved category
# levels named on the trailing K margin so every K-shaped output threads them.
# levels is passed in rather than derived from y, since y is a factor for the
# label path but the count matrix itself for the count path (levels(y.train)
# vs. colnames(y.train)/index fallback; see bart2Multinomial/
# bart2MultinomialCounts).
#
# varcount is the per-sample per-category split-usage channel: each category
# forest's per-draw variable counts, reshaped like yhat.train to (n.chains x)
# n.samples x p x K with levels on the K margin and predictorNames on the p
# margin, one K-margin array like every other K-shaped output here.
packageMultinomialResults <- function(
  control,
  y,
  levels,
  K,
  samples,
  combineChains,
  predictorNames = NULL
) {
  n.chains <- control@n.chains

  # both the train (n.obs x K x n.samples (x n.chains)) and the test channel
  # reshape identically to the package's draws-first convention with levels
  # named on the trailing K margin; the test channel (yhat.test) is the same
  # softmax blend on the held-out rows, present only when 'test' was supplied.
  # predict.bartMultinomial reshapes its replayed channel through the same map.
  shapeChannel <- function(raw) {
    shapeMultinomialChannel(raw, levels, n.chains, combineChains)
  }

  result <- list(
    call = control@call,
    family = "multinomial",
    levels = levels,
    K = K,
    n.chains = n.chains,
    n.trees = control@n.trees,
    y = y,
    yhat.train = shapeChannel(samples$train),
    # the per-category varcount channel shares the fits' raw p x K x n.samples
    # shape, so the same reshape applies with predictor names on the p margin
    varcount = shapeMultinomialChannel(
      samples$varcount,
      levels,
      n.chains,
      combineChains,
      leadNames = predictorNames
    )
  )
  if (!is.null(samples$test)) {
    result$yhat.test <- shapeChannel(samples$test)
  }
  class(result) <- "bartMultinomial"
  result
}

# The cumulative-probit category probabilities for one posterior draw
# (docs/design/ordinal.md section 1): given the latent means eta (length n) and
# the K-1 finite cutpoints gamma (gamma_1 < ... < gamma_{K-1}), returns the
# n x K matrix P(y = k) = Phi(gamma_k - eta) - Phi(gamma_{k-1} - eta), with
# gamma_0 = -Inf and gamma_K = +Inf supplying the boundary columns of 0 and 1.
# Shared by the fit-time reshape and predict.bartOrdinal's replay.
ordinalCategoryProbabilities <- function(eta, cutpoints) {
  n <- length(eta)
  K <- length(cutpoints) + 1L
  bounds <- matrix(0, n, K + 1L)
  bounds[, K + 1L] <- 1
  for (j in seq_len(K - 1L)) {
    bounds[, j + 1L] <- pnorm(cutpoints[j] - eta)
  }
  bounds[, 2L:(K + 1L), drop = FALSE] - bounds[, 1L:K, drop = FALSE]
}

# The ordinal (cumulative probit) fit path (docs/design/ordinal.md section 5),
# reached from bart2's family = "ordinal" branch. Unlike multinomial's K-forest
# engine, ordinal is a SINGLE forest whose kept samples report the latent
# eta = f(x) (like probit) and, in an ordinal-only run channel, the per-draw
# K-1 cutpoints; the n x K category probabilities are synthesized here from the
# two. A single run(n.burn, n.samples) drives the whole chain (the cutpoint
# channel is aligned with each kept sweep's latent draw), so no per-sample
# state read is needed. dbarts(family = "ordinal") does the 1-based recoding,
# the fixed unit scale, and attaches K, so this reuses the standard bart2
# host-build machinery unchanged. The fit is class "bartOrdinal", never "bart",
# so the K-widened arrays never flow through a single-forest method.
bart2Ordinal <- function(
  matchedCall,
  callingEnv,
  control,
  power,
  base,
  sigdf,
  sigquant,
  dart,
  combineChains,
  prior.scale = NA_real_,
  keepSampler = FALSE
) {
  priors <- buildSamplerPriors(
    matchedCall,
    power,
    base,
    sigdf,
    sigquant,
    nodeK = matchedCall[["k"]],
    priorScale = prior.scale,
    dart = dart,
    splitProbsDefault = formals(dbarts::bart2)[["split.probs"]]
  )

  samplerCall <- buildHostSamplerCall(
    matchedCall,
    control,
    priors,
    family = "ordinal"
  )

  sampler <- eval(samplerCall, envir = callingEnv)

  K <- attr(sampler$control, "bartcore.n.categories")
  levels <- sampler$data@response.levels
  n.chains <- control@n.chains
  n.obs <- length(sampler$data@y)
  n.test <- NROW(sampler$data@x.test)
  n.samples <- control@n.samples

  bc <- bartcoreSampler(sampler, family = "ordinal")
  # bc holds the fit's model; the host keeps only the design and priors it was
  # built from, so a mutation through the returned $fit would change a sampler
  # nothing reads. Mark it, as bart2Multinomial does
  sampler$hostFor <- "bart2(family = \"ordinal\")"

  probsTrain <- array(0, c(n.obs, K, n.samples, n.chains))
  latentTrain <- array(0, c(n.obs, n.samples, n.chains))
  probsTest <- if (n.test > 0L) {
    array(0, c(n.test, K, n.samples, n.chains))
  } else {
    NULL
  }
  latentTest <- if (n.test > 0L) {
    array(0, c(n.test, n.samples, n.chains))
  } else {
    NULL
  }
  cutpointsRaw <- array(0, c(K - 1L, n.samples, n.chains))

  # one run drives the whole chain; the cutpoints ride their own run channel
  # (r$cutpoints, present because the ordinal family carries them), aligned with
  # each kept sweep's latent draw, so no per-sample state read is needed
  r <- bartcoreRun(bc, control@n.burn, n.samples)

  varWidth <- if (n.chains == 1L) {
    nrow(as.matrix(r$varcount))
  } else {
    dim(r$varcount)[1L]
  }
  varcountRaw <- array(0, c(varWidth, n.samples, n.chains))

  for (chain in seq_len(n.chains)) {
    for (s in seq_len(n.samples)) {
      gamma <- channelColumn(r$cutpoints, s, chain, n.chains)
      cutpointsRaw[, s, chain] <- gamma
      etaTrain <- channelColumn(r$train, s, chain, n.chains)
      latentTrain[, s, chain] <- etaTrain
      probsTrain[,, s, chain] <-
        ordinalCategoryProbabilities(etaTrain, gamma)
      if (n.test > 0L) {
        etaTest <- channelColumn(r$test, s, chain, n.chains)
        latentTest[, s, chain] <- etaTest
        probsTest[,, s, chain] <- ordinalCategoryProbabilities(etaTest, gamma)
      }
      varcountRaw[, s, chain] <- channelColumn(r$varcount, s, chain, n.chains)
    }
  }

  # keep the full 3-D (K-1) x n.samples x n.chains cutpoints for predict, which
  # pairs each replayed latent draw with its own thresholds
  cutpointsRawFull <- cutpointsRaw

  # drop the trailing singleton chain margin so the reshapers see the
  # n.chains == 1 layout their multinomial/gaussian siblings emit
  if (n.chains == 1L) {
    probsTrain <- array(probsTrain, dim(probsTrain)[1:3])
    latentTrain <- matrix(latentTrain, n.obs, n.samples)
    if (!is.null(probsTest)) {
      probsTest <- array(probsTest, dim(probsTest)[1:3])
      latentTest <- matrix(latentTest, n.test, n.samples)
    }
    cutpointsRaw <- matrix(cutpointsRaw, K - 1L, n.samples)
    varcountRaw <- matrix(varcountRaw, dim(varcountRaw)[1L], n.samples)
  }

  result <- packageOrdinalResults(
    control,
    sampler,
    levels,
    K,
    probsTrain,
    latentTrain,
    probsTest,
    latentTest,
    cutpointsRaw,
    varcountRaw,
    combineChains
  )
  # keepTrees retains the handles predict.bartOrdinal replays through: bc holds
  # the saved trees (the sweeps wrote them), and the sampler codes newdata to
  # the training columns. cutpoints.raw supplies predict's per-draw thresholds
  # in the raw (K-1) x n.samples x n.chains layout that pairs with the replayed
  # latent draws, so no re-run is needed.
  if (control@keepTrees) {
    result$bc <- bc
    result$cutpoints.raw <- cutpointsRawFull
  }
  if (control@keepTrees || keepSampler) {
    result$fit <- sampler
  }
  result
}

# Assemble a bart2(family = "ordinal") fit from the synthesized channels
# (docs/design/ordinal.md section 5). yhat.train/test are the n x K category
# probabilities in the multinomial draws-first convention (levels named on the
# K margin); cutpoints are the per-draw K-1 thresholds, the ordinal analog of
# gaussian's sigma; y is rebuilt as the ordered factor of observed categories
# for residuals and reporting.
packageOrdinalResults <- function(
  control,
  sampler,
  levels,
  K,
  probsTrain,
  latentTrain,
  probsTest,
  latentTest,
  cutpointsRaw,
  varcountRaw,
  combineChains
) {
  n.chains <- control@n.chains

  varcount <- nameVarcount(
    varcountRaw,
    colnames(sampler$data@x),
    n.chains,
    combineChains
  )

  result <- list(
    call = control@call,
    family = "ordinal",
    levels = levels,
    K = K,
    n.chains = n.chains,
    n.trees = control@n.trees,
    # the observed category as an ordered factor over the resolved levels, for
    # residuals() and printing (sampler$data@y holds 1-based codes)
    y = factor(
      levels[as.integer(sampler$data@y)],
      levels = levels,
      ordered = TRUE
    ),
    cutpoints = convertSamplesFromDbartsToBart(
      cutpointsRaw,
      n.chains,
      combineChains
    ),
    yhat.train = shapeMultinomialChannel(
      probsTrain,
      levels,
      n.chains,
      combineChains
    ),
    # the latent eta = f(x) draws (type = "bart"/"link"), the single-column
    # probit-scale channel the K probabilities are formed from
    latent.train = convertSamplesFromDbartsToBart(
      latentTrain,
      n.chains,
      combineChains
    ),
    varcount = varcount
  )
  if (!is.null(probsTest)) {
    result$yhat.test <- shapeMultinomialChannel(
      probsTest,
      levels,
      n.chains,
      combineChains
    )
    result$latent.test <- convertSamplesFromDbartsToBart(
      latentTest,
      n.chains,
      combineChains
    )
  }
  class(result) <- "bartOrdinal"
  result
}

# The negative-binomial mean counts for one posterior draw (docs/design/
# negative-binomial.md section 1): given the log-odds latent psi = f(x) + o and
# the dispersion r, the mean is mu_i = r exp(psi_i). Shared by the fit-time
# reshape and predict.bartNegbin's replay.
negbinMeanCounts <- function(psi, r) {
  r * exp(psi)
}

# One posterior-predictive count per (draw, observation) of a negative-binomial
# fit: y ~ NB(size = r, mean = mu) with mu the draw's mean count and r its
# dispersion (docs/design/negative-binomial.md). mu is a (chains x) draws x obs
# array; r is the draws-shaped dispersion, broadcast across the observation
# margin so each column shares its draw's r. Only this and predict's ppd touch
# the RNG, so type = "ev"/"bart" stay draw-neutral.
negbinPpd <- function(mu, r) {
  d <- dim(mu)
  array(
    rnbinom(length(mu), size = as.vector(array(r, d)), mu = as.vector(mu)),
    d
  )
}

# The negative-binomial count fit path (docs/design/negative-binomial.md
# sections 4-5), reached from bart2's family = "nbinom" branch. A SINGLE forest
# fits the log-odds latent psi = f(x) + o (like logistic under the Polya-Gamma
# augmentation); the mean counts mu = r exp(psi) and the per-draw dispersion r
# are synthesized here. The run is driven one kept sample at a time (through the
# low-level bc) because mu = r exp(psi) pairs each sweep's latent draw with that
# sweep's r; r itself comes from the run's own per-draw dispersion channel, so no
# state is serialized per sweep. dbarts(family = "nbinom") does the count
# validation, the fixed
# unit scale, and attaches the dispersion spec, so this reuses the standard
# bart2 host-build machinery. The fit is class "bartNegbin", never "bart".
bart2Negbin <- function(
  matchedCall,
  callingEnv,
  control,
  power,
  base,
  sigdf,
  sigquant,
  dart,
  combineChains,
  prior.scale = NA_real_,
  keepSampler = FALSE
) {
  priors <- buildSamplerPriors(
    matchedCall,
    power,
    base,
    sigdf,
    sigquant,
    nodeK = matchedCall[["k"]],
    priorScale = prior.scale,
    dart = dart,
    splitProbsDefault = formals(dbarts::bart2)[["split.probs"]]
  )

  samplerCall <- buildHostSamplerCall(
    matchedCall,
    control,
    priors,
    family = "nbinom"
  )

  sampler <- eval(samplerCall, envir = callingEnv)

  n.chains <- control@n.chains
  n.obs <- length(sampler$data@y)
  n.test <- NROW(sampler$data@x.test)
  n.samples <- control@n.samples

  bc <- bartcoreSampler(sampler, family = "nbinom")
  # the host owns no model past this point; see bart2Ordinal. It matters twice
  # here: $getDispersion() on the shell would answer with the placeholder's own
  # r rather than the fit's
  sampler$hostFor <- "bart2(family = \"nbinom\")"

  latentTrain <- array(0, c(n.obs, n.samples, n.chains))
  meanTrain <- array(0, c(n.obs, n.samples, n.chains))
  latentTest <- if (n.test > 0L) {
    array(0, c(n.test, n.samples, n.chains))
  } else {
    NULL
  }
  meanTest <- if (n.test > 0L) {
    array(0, c(n.test, n.samples, n.chains))
  } else {
    NULL
  }
  # r is a scalar per (sample, chain), so it rides a sigma-shaped matrix
  dispersionRaw <- matrix(0, n.samples, n.chains)
  varcountRaw <- NULL

  for (s in seq_len(n.samples)) {
    # the first kept sample absorbs the burn-in, so every run keeps one sample
    r <- bartcoreRun(bc, if (s == 1L) control@n.burn else 0L, 1L)
    # sigma-shaped, so a single-sample run's channel is exactly this row
    dispersionRaw[s, ] <- r$dispersion
    if (is.null(varcountRaw)) {
      varWidth <- if (n.chains == 1L) {
        nrow(as.matrix(r$varcount))
      } else {
        dim(r$varcount)[1L]
      }
      varcountRaw <- array(0, c(varWidth, n.samples, n.chains))
    }
    for (chain in seq_len(n.chains)) {
      rDraw <- dispersionRaw[s, chain]
      psiTrain <- channelColumn(r$train, 1L, chain, n.chains)
      latentTrain[, s, chain] <- psiTrain
      meanTrain[, s, chain] <- negbinMeanCounts(psiTrain, rDraw)
      if (n.test > 0L) {
        psiTest <- channelColumn(r$test, 1L, chain, n.chains)
        latentTest[, s, chain] <- psiTest
        meanTest[, s, chain] <- negbinMeanCounts(psiTest, rDraw)
      }
      varcountRaw[, s, chain] <- channelColumn(r$varcount, 1L, chain, n.chains)
    }
  }

  # drop the trailing singleton chain margin so the reshapers see the
  # n.chains == 1 layout their gaussian siblings emit (dispersionRaw keeps its
  # matrix form, the sigma channel's shape)
  if (n.chains == 1L) {
    latentTrain <- matrix(latentTrain, n.obs, n.samples)
    meanTrain <- matrix(meanTrain, n.obs, n.samples)
    if (n.test > 0L) {
      latentTest <- matrix(latentTest, n.test, n.samples)
      meanTest <- matrix(meanTest, n.test, n.samples)
    }
    varcountRaw <- matrix(varcountRaw, dim(varcountRaw)[1L], n.samples)
  }

  result <- packageNegbinResults(
    control,
    sampler,
    latentTrain,
    meanTrain,
    latentTest,
    meanTest,
    dispersionRaw,
    varcountRaw,
    combineChains
  )
  # keepTrees retains the handles predict.bartNegbin replays through: bc holds
  # the saved trees, and dispersion.raw supplies predict's per-draw r in the raw
  # n.samples x n.chains layout that pairs with the replayed latent draws.
  if (control@keepTrees) {
    result$bc <- bc
    result$dispersion.raw <- dispersionRaw
  }
  if (control@keepTrees || keepSampler) {
    result$fit <- sampler
  }
  result
}

# Assemble a bart2(family = "nbinom") fit from the synthesized channels
# (docs/design/negative-binomial.md section 4). yhat.train/test are the mean
# counts mu = r exp(f + o), the reported deliverable; latent.train/test are the
# log-odds latent psi draws (type = "bart"/"link"); dispersion is the per-draw r,
# the count analog of gaussian's sigma; y is the observed counts.
packageNegbinResults <- function(
  control,
  sampler,
  latentTrain,
  meanTrain,
  latentTest,
  meanTest,
  dispersionRaw,
  varcountRaw,
  combineChains
) {
  n.chains <- control@n.chains

  varcount <- nameVarcount(
    varcountRaw,
    colnames(sampler$data@x),
    n.chains,
    combineChains
  )

  result <- list(
    call = control@call,
    family = "nbinom",
    n.chains = n.chains,
    n.trees = control@n.trees,
    y = sampler$data@y,
    # the per-draw dispersion r, the sigma-shaped count analog
    dispersion = convertSamplesFromDbartsToBart(
      if (n.chains == 1L) dispersionRaw[, 1L] else dispersionRaw,
      n.chains,
      combineChains
    ),
    # the mean counts mu = r exp(f + o) (type = "ev"/"response")
    yhat.train = convertSamplesFromDbartsToBart(
      meanTrain,
      n.chains,
      combineChains
    ),
    # the log-odds latent psi = f(x) + o draws (type = "bart"/"link")
    latent.train = convertSamplesFromDbartsToBart(
      latentTrain,
      n.chains,
      combineChains
    ),
    varcount = varcount
  )
  if (!is.null(meanTest)) {
    result$yhat.test <- convertSamplesFromDbartsToBart(
      meanTest,
      n.chains,
      combineChains
    )
    result$latent.test <- convertSamplesFromDbartsToBart(
      latentTest,
      n.chains,
      combineChains
    )
  }
  class(result) <- "bartNegbin"
  result
}

# Splits a hurdle response into its two ingested parts (docs/design/hurdle.md
# sections 0-1): the occupancy indicator z = 1{y > 0} over all n, and the
# subset mask {i : y_i > 0} with the positive part's working response
# log(y[S]). y must be finite and non-negative (the nbinom non-negative-count
# precedent, docs/design/negative-binomial.md section 4, extended to a
# continuous response) and carry at least one exact zero and one positive
# value - the two parts a hurdle needs in order to fit either part at all.
splitHurdleResponse <- function(y) {
  y <- as.double(y)
  if (anyNA(y) || any(!is.finite(y)) || any(y < 0)) {
    stop(
      "family = \"hurdle.lognormal\" requires a non-negative, finite ",
      "numeric response"
    )
  }
  positive <- y > 0
  if (!any(positive)) {
    stop(
      "family = \"hurdle.lognormal\" requires at least one positive response ",
      "value for the positive-part fit"
    )
  }
  if (all(positive)) {
    stop(
      "family = \"hurdle.lognormal\" requires at least one exact zero; a ",
      "response with no zeros is an ordinary continuous fit"
    )
  }
  list(
    z = as.double(positive),
    positive = positive,
    logPositive = log(y[positive])
  )
}

# The hurdle.lognormal / twopart fit path (docs/design/hurdle.md), reached
# from bart2's family = "hurdle.lognormal" branch. Composed R-side from two
# ordinary single-forest fits - never a coupled engine model, section 0 - so
# this simply calls bart2() twice, once per component, at two independently
# derived seeds (section 13 hardening b: a shared seed correlates the two
# chains and biases the combined credible interval). The positive fit's
# 'test' is forced to the FULL training x (hardening c) so its in-sample
# fitted()/extract() carries E[y | y > 0, x] at the zero rows it never
# trained on. Reusing bart2() itself - rather than building the sampler
# directly, as bart2Negbin does - makes the reduction to two standalone
# bart2() calls exact by construction (benchmarks/R/hurdle-reduction.R): the
# wrapper's component fits ARE those calls. family is always explicit and the
# matrix interface only: the response is split before any dbartsData
# ingestion runs (the discrete-time hazard precedent, R/dbarts.R's
# extractSurvivalTimes), which a formula LHS cannot supply.
bart2Hurdle <- function(matchedCall, callingEnv, control, formula, data, seed) {
  if (
    is.formula(formula) ||
      inherits(formula, "dbartsData") ||
      inherits(formula, "dgCMatrix") ||
      missing(data)
  ) {
    stop(
      "family = \"hurdle.lognormal\" fits currently use the matrix interface ",
      "- bart2(x.train, y.train, family = \"hurdle.lognormal\")"
    )
  }

  split <- splitHurdleResponse(data)
  xPositive <- formula[split$positive, , drop = FALSE]

  seeds <- c(NA_integer_, NA_integer_)
  if (!is.na(seed)) {
    # independent per-component seeds derived deterministically from the
    # user's seed (the per-chain seed derivation rbart_vi uses); politely
    # restore the caller's RNG stream afterward
    oldSeed <- .GlobalEnv[[".Random.seed"]]
    set.seed(seed)
    seeds <- sample.int(.Machine$integer.max, 2L)
    if (!is.null(oldSeed)) {
      .GlobalEnv$.Random.seed <- oldSeed
    }
  }

  # redirectCall forwards every bart2 formal the caller supplied, including
  # names this outer site already diagnosed against "hurdle.lognormal"
  # (3.c.4's table) - left alone, each component would re-diagnose them
  # against ITS OWN forced family, both firing a warning where the outer
  # site did not (dispersion/breaks/max.rows, inert on both components:
  # a second and third copy of the same diagnosis) and, on the occupancy
  # call, one the outer site correctly did NOT raise (sigest/sigdf/sigquant:
  # genuinely live on the positive half, so a "probit" warning about them
  # is a false diagnostic, not merely a redundant one). Strip what the
  # outer site already covers before either component call runs.
  gatedOnBoth <- c("dispersion", "breaks", "max.rows")
  gatedOnOccupancyOnly <- c("sigest", "sigdf", "sigquant")

  occupancyCall <- redirectCall(matchedCall, dbarts::bart2)
  occupancyCall[c(gatedOnBoth, gatedOnOccupancyOnly)] <- NULL
  occupancyCall$formula <- formula
  occupancyCall$data <- split$z
  occupancyCall$family <- "probit"
  occupancyCall$seed <- seeds[1L]
  occupancyCall$keepTrees <- control@keepTrees
  occupancy <- eval(occupancyCall, callingEnv)

  positiveCall <- redirectCall(matchedCall, dbarts::bart2)
  positiveCall[gatedOnBoth] <- NULL
  positiveCall$formula <- xPositive
  positiveCall$data <- split$logPositive
  positiveCall$test <- formula
  positiveCall$family <- "gaussian"
  positiveCall$seed <- seeds[2L]
  positiveCall$keepTrees <- control@keepTrees
  positive <- eval(positiveCall, callingEnv)

  result <- list(
    call = control@call,
    family = "hurdle.lognormal",
    # the original non-negative response over all n, so residuals() can take
    # y - E[y | x] on the natural scale (neither component stores it: the
    # occupancy fit keeps the 1{y > 0} indicator, the positive fit log(y[S]))
    y = as.double(data),
    occupancy = occupancy,
    positive = positive
  )
  class(result) <- "bartHurdle"
  result
}

# S(t | x) draws from an AFT linear predictor and its per-draw sigma, in the
# uncombined convention (chains x samples x observations) where the sigma
# draws align with the fit draws unambiguously - the loglik channel's
# approach. Shared by the bart and rbart methods; the caller supplies the
# linear predictor (BART component for bart, BART + drawn intercepts for
# rbart) and combines the chain margin at the end.
survivalProbabilitiesFromDraws <- function(
  linearPredictor,
  sigma,
  times,
  n.chains,
  combineChains
) {
  if (is.null(dim(sigma))) {
    sigma <- uncombineChains(as.vector(sigma), n.chains)
  }

  lpDims <- dim(linearPredictor)
  drawDims <- lpDims[-length(lpDims)]
  numObservations <- lpDims[length(lpDims)]
  numTimes <- length(times)
  if (prod(drawDims) != length(sigma)) {
    stop("the fit's draw count does not match its sigma draws")
  }

  result <- array(0, c(drawDims, numTimes, numObservations))
  # in column-major order the draw margins vary fastest, so the sigma draws
  # recycle across observations exactly as pointwiseLogLikelihood's do
  scale <- as.vector(sigma)
  if (length(drawDims) == 1L) {
    for (j in seq_len(numTimes)) {
      result[, j, ] <- pnorm(
        (log(times[j]) - linearPredictor) / scale,
        lower.tail = FALSE
      )
    }
  } else {
    for (j in seq_len(numTimes)) {
      result[,, j, ] <- pnorm(
        (log(times[j]) - as.vector(linearPredictor)) / scale,
        lower.tail = FALSE
      )
    }
  }

  if (combineChains && length(drawDims) > 1L) {
    # merge the chain and sample margins as combineChains() does, samples
    # varying fastest within each chain
    result <- aperm(result, c(2L, 1L, 3L, 4L))
    dim(result) <- c(prod(drawDims), numTimes, numObservations)
  }
  result
}

# Survival-probability draws from a discrete-time hazard fit
# (docs/design/survival.md section 4). The fit is an ordinary binary fit on
# the person-period-expanded rows, so the per-(subject, period) hazards are
# h(k | x) = g(f(x, k) + o) through the fit's link (probit/logistic), and
# S(t | x) = prod_{k : periods[k] <= t} (1 - h(k | x)). This ALWAYS re-expands
# its subjects to the full grid and predicts - the training design is ragged
# (subject i has only its at-risk rows), so a stored per-subject fit cannot
# supply hazards past its own event/censoring period - and therefore requires
# keepTrees unconditionally, where aft's training path never does. Requested
# `times` are horizons (default: the grid periods); S cumulates (1 - h)
# through each. Returns draws per the package's three-tier convention, the
# shape aft's method uses (draws x times x observations, a chain margin under
# combineChains = FALSE).
hazardSurvivalProbabilities <- function(object, times, newdata, combineChains) {
  if (is.null(object[["fit"]])) {
    stop(
      "survivalProbabilities on a discrete-time hazard fit requires the ",
      "trees; refit with keepTrees = TRUE"
    )
  }
  periods <- object$periods
  K <- length(periods)
  if (is.null(times)) {
    times <- periods
  }
  times <- as.double(times)
  if (length(times) == 0L || any(!is.finite(times)) || any(times <= 0)) {
    stop("'times' must be finite and positive")
  }

  periodCol <- ncol(object$fit$data@x)
  if (is.null(newdata)) {
    # reconstruct the per-subject covariates from the coded expanded design:
    # every subject is at risk in period 1, so the period-1 rows are the
    # subjects in order, their covariate columns the (coded) per-subject x
    fitX <- extract(object$fit, "predictors")
    subjectCov <- fitX[fitX[, periodCol] == 1L, -periodCol, drop = FALSE]
    n <- nrow(subjectCov)
    bigX <- cbind(
      subjectCov[rep(seq_len(n), times = K), , drop = FALSE],
      rep(seq_len(K), each = n)
    )
  } else if (is.data.frame(newdata)) {
    n <- nrow(newdata)
    bigX <- newdata[rep(seq_len(n), times = K), , drop = FALSE]
    bigX[["period"]] <- rep(seq_len(K), each = n)
  } else {
    newdata <- as.matrix(newdata)
    n <- nrow(newdata)
    bigX <- cbind(
      newdata[rep(seq_len(n), times = K), , drop = FALSE],
      rep(seq_len(K), each = n)
    )
  }

  # hazards through the correct link (type = "ev" keys on $family, the binary
  # token); predict codes bigX to the training columns and replays the trees
  haz <- predict(object, bigX, type = "ev", combineChains = FALSE)
  drawDims <- dim(haz)[-length(dim(haz))]
  D <- prod(drawDims)
  numTimes <- length(times)

  # split the (n * K) obs margin into (n fast, K slow) - the period-major order
  # bigX was built in - and cumulate (1 - h) across the period margin
  surv <- 1 - haz
  dim(surv) <- c(D, n, K)
  if (K >= 2L) {
    for (k in 2:K) {
      surv[,, k] <- surv[,, k - 1L] * surv[,, k]
    }
  }
  out <- array(0, c(D, numTimes, n))
  for (j in seq_len(numTimes)) {
    m <- sum(periods <= times[j])
    out[, j, ] <- if (m == 0L) 1 else surv[,, m]
  }

  dim(out) <- c(drawDims, numTimes, n)
  if (combineChains && length(drawDims) > 1L) {
    # merge the chain and sample margins as combineChains() does, samples
    # fastest within each chain (aft's survivalProbabilitiesFromDraws move)
    out <- aperm(out, c(2L, 1L, 3L, 4L))
    dim(out) <- c(D, numTimes, n)
  }
  out
}

# Survival-probability draws from an AFT fit (docs/design/survival.md).
# Under the log-normal model log T = f(x) + sigma eps, S(t | x) =
# 1 - Phi((log t - f(x)) / sigma), evaluated at every posterior draw of f
# and sigma. Returns draws per the package's three-tier convention
# (extract = draws, fitted = mean, ci.level = interval): users take means
# and quantiles over the draw margin themselves. newdata predicts out of
# sample (requires a fit kept with keepTrees); NULL uses the training fits.
#
# A discrete-time hazard fit carries the $periods marker (never $family, which
# reads the binary link token); it takes the hazard branch above instead.
survivalProbabilities.bart <- function(
  object,
  times,
  newdata = NULL,
  combineChains = TRUE,
  ...
) {
  if (!is.null(object[["periods"]])) {
    return(hazardSurvivalProbabilities(
      object,
      if (missing(times)) NULL else times,
      newdata,
      combineChains
    ))
  }
  if (!identical(object[["family"]], "aft")) {
    stop("survivalProbabilities requires an aft (survival) fit")
  }
  times <- as.double(times)
  if (length(times) == 0L || any(!is.finite(times)) || any(times <= 0)) {
    stop("'times' must be finite and positive")
  }

  n.chains <- if (!is.null(object[["fit"]])) {
    object$fit$control@n.chains
  } else {
    object$n.chains
  }

  linearPredictor <- if (is.null(newdata)) {
    extract(object, type = "bart", sample = "train", combineChains = FALSE)
  } else {
    predict(object, newdata, type = "bart", combineChains = FALSE)
  }

  survivalProbabilitiesFromDraws(
    linearPredictor,
    object[["sigma"]],
    times,
    n.chains,
    combineChains
  )
}

# Grouped (random-intercept) AFT survival curves (docs/design/survival.md,
# riAFTBART's model). The linear predictor is E[log T | x, group] = f(x) plus
# the drawn group intercept - sourced from the "ev" channel (extract for the
# training data, predict for newdata), NOT the bare BART component, so the
# intercepts enter the curve. "ev" is on the log scale here (aft carries
# sigma, so it is not probability-transformed). newdata needs group.by (an
# unseen group draws its intercept from N(0, tau), inherited from predict).
survivalProbabilities.rbart <- function(
  object,
  times,
  newdata = NULL,
  group.by,
  combineChains = TRUE,
  ...
) {
  if (!identical(object[["family"]], "aft")) {
    stop("survivalProbabilities requires an aft (survival) rbart fit")
  }
  times <- as.double(times)
  if (length(times) == 0L || any(!is.finite(times)) || any(times <= 0)) {
    stop("'times' must be finite and positive")
  }

  n.chains <- if (is.null(object$n.chains)) {
    length(object$fit)
  } else {
    object$n.chains
  }

  linearPredictor <- if (is.null(newdata)) {
    extract(object, type = "ev", sample = "train", combineChains = FALSE)
  } else {
    if (missing(group.by)) {
      stop("'group.by' must be supplied when 'newdata' is given")
    }
    # predict.rbart maps intercepts through the factor levels (an unseen
    # level draws from N(0, tau)); coerce as rbart_vi does its group.by
    predict(
      object,
      newdata,
      group.by = as.factor(group.by),
      type = "ev",
      combineChains = FALSE
    )
  }

  survivalProbabilitiesFromDraws(
    linearPredictor,
    object[["sigma"]],
    times,
    n.chains,
    combineChains
  )
}

bart <- function(
  x.train,
  y.train,
  x.test = matrix(0.0, 0L, 0L),
  sigest = NA_real_,
  sigdf = 3.0,
  sigquant = 0.90,
  k = 2.0,
  prior.scale = NA_real_,
  power = 2.0,
  base = 0.95,
  splitprobs = 1 / numvars,
  binaryOffset = 0.0,
  weights = NULL,
  ntree = 200L,
  ndpost = 1000L,
  nskip = 100L,
  printevery = 100L,
  keepevery = 1L,
  keeptrainfits = TRUE,
  usequants = FALSE,
  numcut = 100L,
  printcutoffs = 0L,
  verbose = TRUE,
  nchain = 1L,
  nthread = 1L,
  combinechains = TRUE,
  keeptrees = FALSE,
  keepcall = TRUE,
  sampleronly = FALSE,
  seed = NA_integer_,
  proposalprobs = NULL,
  keepsampler = keeptrees,
  resid.dist = gaussian
) {
  # forwarded to dbarts() unevaluated (as the prior expressions are), so a bare
  # gaussian()/student() resolves in dbarts()'s residual-distribution vocabulary
  residDist <- substitute(resid.dist)
  # coerce eagerly, naming the argument as the caller typed it - dbartsControl
  # re-coerces its own (already-integer) inputs and would otherwise blame its
  # internal slot names (n.burn/n.trees/...) for a bad ntree/nskip/... value
  ntree <- coerceOrError(ntree, "integer")
  nskip <- coerceOrError(nskip, "integer")
  nchain <- coerceOrError(nchain, "integer")
  nthread <- coerceOrError(nthread, "integer")
  keepevery <- coerceOrError(keepevery, "integer")
  printevery <- coerceOrError(printevery, "integer")
  printcutoffs <- coerceOrError(printcutoffs, "integer")
  numcut <- coerceOrError(numcut, "integer")
  ndpost <- coerceOrError(ndpost, "integer")

  # bart() is the frozen BayesTree shim and does not package an ordinal fit; an
  # ordered-factor response (which dbarts() would auto-dispatch to ordinal) is
  # refused up front, pointing to bart2 (docs/design/ordinal.md section 5). The
  # matrix interface names the response directly; a formula response is caught
  # by the family backstop after the sampler is built.
  if (
    !is.formula(x.train) &&
      is.factor(y.train) &&
      is.ordered(y.train) &&
      nlevels(y.train) >= 3L
  ) {
    stop(
      "bart() does not fit ordinal (ordered-factor) responses; use ",
      "bart2(x.train, y.train, family = \"ordinal\")"
    )
  }

  control <- dbartsControl(
    keepTrainingFits = as.logical(keeptrainfits),
    useQuantiles = as.logical(usequants),
    keepTrees = FALSE,
    n.burn = nskip,
    n.trees = ntree,
    n.chains = nchain,
    n.threads = nthread,
    n.thin = keepevery,
    printEvery = printevery,
    printCutoffs = printcutoffs,
    n.cuts = numcut,
    rngSeed = as.integer(seed)
  )
  matchedCall <- if (keepcall) match.call() else call("NULL")
  control@call <- matchedCall
  control@n.burn <- control@n.burn %/% control@n.thin
  control@printEvery <- control@printEvery %/% control@n.thin
  keepsampler <- keepsampler || control@keepTrees
  if (control@n.burn == 0L && keeptrees == TRUE) {
    control@keepTrees <- TRUE
  }
  if (control@n.burn > 0L) {
    control@keepTrees <- FALSE
  }
  ndpost <- ndpost %/% control@n.thin

  priors <- buildSamplerPriors(
    matchedCall,
    power,
    base,
    sigdf,
    sigquant,
    nodeK = if (!is.null(matchedCall[["k"]])) matchedCall[["k"]] else k,
    priorScale = prior.scale,
    splitProbsName = "splitprobs",
    splitProbsDefault = formals(dbarts::bart)[["splitprobs"]]
  )
  tree.prior <- priors$tree.prior
  node.prior <- priors$node.prior
  resid.prior <- priors$resid.prior

  # the frozen BayesTree-compatibility shim keeps dummy expansion
  args <- list(
    formula = x.train,
    data = y.train,
    test = x.test,
    subset = NULL,
    weights = weights,
    offset = binaryOffset,
    verbose = as.logical(verbose),
    n.samples = as.integer(ndpost),
    tree.prior = tree.prior,
    node.prior = node.prior,
    resid.prior = resid.prior,
    resid.dist = residDist,
    proposal.probs = proposalprobs,
    control = control,
    sigma = as.numeric(sigest),
    factors = "indicators",
    missing = "error"
  )
  sampler <- do.call(dbarts::dbarts, args, envir = parent.frame(1L))

  # formula-response backstop for the pre-check above: dbarts() auto-dispatched
  # an ordered-factor response to ordinal, which bart() cannot package
  if (identical(sampler$model@family, "ordinal")) {
    stop(
      "bart() does not fit ordinal (ordered-factor) responses; use ",
      "bart2(..., family = \"ordinal\")"
    )
  }

  if (sampleronly) {
    return(sampler)
  }

  control <- sampler$control

  burn <- runWithBurnIn(sampler, control, keeptrees)
  samples <- burn$samples
  burnInSigma <- burn$burnInSigma
  burnInK <- burn$burnInK

  result <- packageBartResults(
    sampler,
    samples,
    burnInSigma,
    burnInK,
    combinechains,
    keepsampler
  )

  result
}

makeind <- function(x, all = TRUE) {
  ignored <- all ## for R check # nolint: object_usage_linter.
  makeModelMatrixFromDataFrame(x, TRUE)
}
