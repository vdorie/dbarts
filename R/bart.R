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

  varcount <- convertSamplesFromDbartsToBart(
    samples$varcount,
    n.chains,
    combineChains
  )
  if (!is.null(colnames(fit$data@x)) && !is.null(dim(varcount))) {
    dimnames(varcount) <- if (length(dim(varcount)) > 2L) {
      list(NULL, NULL, colnames(fit$data@x))
    } else {
      list(NULL, colnames(fit$data@x))
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
  dart = FALSE,
  splitProbsName = "split.probs",
  splitProbsDefault = NULL
) {
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

  if (!is.null(nodeK)) {
    node.prior <- quote(normal(k))
    node.prior[[2L]] <- nodeK
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
  proposal.probs = NULL,
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
    "ordinal"
  ),
  missing = c("incorporate", "error"),
  resid.dist = gaussian,
  ...
) {
  matchedCall <- match.call()
  callingEnv <- parent.frame()
  family <- match.arg(family)

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

  controlCall <- redirectCall(matchedCall, dbarts::dbartsControl)
  missingDefaultArgs <- names(formals(bart2))[
    names(formals(bart2)) %in%
      names(formals(dbarts::dbartsControl)) &
      names(formals(bart2)) %not_in% names(matchedCall)
  ]
  if (length(missingDefaultArgs) > 0L) {
    currentEnv <- sys.frame(sys.nframe())
    controlCall[missingDefaultArgs] <- lapply(
      formals(bart2)[missingDefaultArgs],
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
        "family = \"multinomial\" does not support 'weights' this arc"
      )
    }
    if (!missing(offset)) {
      stop(
        "family = \"multinomial\" does not support 'offset' this arc"
      )
    }
    if (!missing(subset)) {
      stop(
        "family = \"multinomial\" does not support 'subset' this arc"
      )
    }
    if (isTRUE(samplerOnly)) {
      stop(
        "family = \"multinomial\" does not support 'samplerOnly' this arc"
      )
    }
    grownSweeps <- as.integer(n.grow.sweeps)[1L]
    if (!is.null(warm.start) || (!is.na(grownSweeps) && grownSweeps > 0L)) {
      stop(
        "family = \"multinomial\" does not support 'warm.start' or ",
        "'n.grow.sweeps' this arc"
      )
    }
    if (!control@keepTrainingFits) {
      stop(
        "family = \"multinomial\" requires keepTrainingFits = TRUE (the ",
        "default): there is no test surface to fall back on"
      )
    }
    if (control@n.samples <= 0L) {
      stop("family = \"multinomial\" requires a positive 'n.samples'")
    }
    if (inherits(formula, "dbartsData")) {
      stop(
        "family = \"multinomial\" does not support a pre-built dbartsData ",
        "object this arc; use the formula interface or the matrix ",
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
        combineChains
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
      combineChains
    ))
  }

  # ordinal (cumulative probit, docs/design/ordinal.md): a single-forest fixed-
  # scale latent model whose n x K category probabilities and per-draw cutpoints
  # are synthesized in R (the engine reports only the latent eta), dispatched
  # here rather than threaded through the standard single-forest path so its
  # K-widened fit object (class "bartOrdinal", never "bart") stays distinct.
  if (family == "ordinal") {
    if (isTRUE(samplerOnly)) {
      stop("family = \"ordinal\" does not support 'samplerOnly' this arc")
    }
    grownSweeps <- as.integer(n.grow.sweeps)[1L]
    if (!is.null(warm.start) || (!is.na(grownSweeps) && grownSweeps > 0L)) {
      stop(
        "family = \"ordinal\" does not support 'warm.start' or ",
        "'n.grow.sweeps' this arc"
      )
    }
    if (!control@keepTrainingFits) {
      stop(
        "family = \"ordinal\" requires keepTrainingFits = TRUE (the default): ",
        "the category probabilities are built from the training latent fits"
      )
    }
    return(bart2Ordinal(
      matchedCall,
      callingEnv,
      control,
      power,
      base,
      sigdf,
      sigquant,
      dart,
      combineChains
    ))
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
    dart = dart,
    splitProbsDefault = formals(dbarts::bart2)[["split.probs"]]
  )
  tree.prior <- priors$tree.prior
  node.prior <- priors$node.prior
  resid.prior <- priors$resid.prior

  samplerCall <- redirectCall(matchedCall, dbarts::dbarts)
  samplerCall$control <- control
  samplerCall$n.samples <- NULL
  samplerCall$tree.prior <- tree.prior
  samplerCall$node.prior <- node.prior
  samplerCall$resid.prior <- resid.prior
  samplerCall$sigma <- as.numeric(sigest)

  sampler <- eval(samplerCall, envir = callingEnv)
  if (samplerOnly == TRUE) {
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

  burnInSigma <- NULL
  burnInK <- NULL
  if (n.burn > 0L) {
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
# exactly for a given seed.
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
  combineChains
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
    dart = dart,
    splitProbsDefault = formals(dbarts::bart2)[["split.probs"]]
  )

  # the host sampler: identical machinery to bart2's normal path, but its
  # response is the integer label vector (as doubles - a placeholder the
  # multinomial engine never reads, since the category labels ride
  # separately) and no 'family' is forwarded (dbarts() does not know
  # "multinomial"; whatever "auto" resolves the placeholder to is immaterial)
  samplerCall <- redirectCall(matchedCall, dbarts::dbarts)
  samplerCall$control <- control
  samplerCall$n.samples <- NULL
  samplerCall$data <- as.double(labels)
  samplerCall$family <- NULL
  samplerCall$tree.prior <- priors$tree.prior
  samplerCall$node.prior <- priors$node.prior
  samplerCall$resid.prior <- priors$resid.prior
  samplerCall$sigma <- as.numeric(sigest)

  sampler <- eval(samplerCall, envir = callingEnv)

  bc <- bartcoreMultinomialSampler(sampler, labels, K = K)
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
  if (control@keepTrees) {
    result$bc <- bc
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
# k3counts scenario).
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
  combineChains
) {
  K <- ncol(y)

  priors <- buildSamplerPriors(
    matchedCall,
    power,
    base,
    sigdf,
    sigquant,
    nodeK = matchedCall[["k"]],
    dart = dart,
    splitProbsDefault = formals(dbarts::bart2)[["split.probs"]]
  )

  samplerCall <- redirectCall(matchedCall, dbarts::dbarts)
  samplerCall$control <- control
  samplerCall$n.samples <- NULL
  samplerCall$data <- as.double(y[, 1L])
  samplerCall$family <- NULL
  samplerCall$tree.prior <- priors$tree.prior
  samplerCall$node.prior <- priors$node.prior
  samplerCall$resid.prior <- priors$resid.prior
  samplerCall$sigma <- as.numeric(sigest)

  sampler <- eval(samplerCall, envir = callingEnv)

  bc <- bartcoreMultinomialCountSampler(sampler, y, K = K)
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
# probabilities are already the identified quantity the engine reports (Q3 of
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
# margin. Symmetric to mbart2's per-category varcount, as an array where every
# other K-output is one.
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
  combineChains
) {
  priors <- buildSamplerPriors(
    matchedCall,
    power,
    base,
    sigdf,
    sigquant,
    nodeK = matchedCall[["k"]],
    dart = dart,
    splitProbsDefault = formals(dbarts::bart2)[["split.probs"]]
  )

  samplerCall <- redirectCall(matchedCall, dbarts::dbarts)
  samplerCall$control <- control
  samplerCall$n.samples <- NULL
  samplerCall$family <- "ordinal"
  samplerCall$tree.prior <- priors$tree.prior
  samplerCall$node.prior <- priors$node.prior
  samplerCall$resid.prior <- priors$resid.prior

  sampler <- eval(samplerCall, envir = callingEnv)

  K <- attr(sampler$control, "bartcore.n.categories")
  levels <- sampler$data@response.levels
  n.chains <- control@n.chains
  n.obs <- length(sampler$data@y)
  n.test <- NROW(sampler$data@x.test)
  n.samples <- control@n.samples

  bc <- bartcoreSampler(sampler, family = "ordinal")

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

  # one chain's sample-s column of a run channel; the engine drops the trailing
  # chain margin when n.chains == 1
  colOf <- function(channel, s, chain) {
    if (n.chains == 1L) channel[, s] else channel[, s, chain]
  }

  varWidth <- if (n.chains == 1L) {
    nrow(as.matrix(r$varcount))
  } else {
    dim(r$varcount)[1L]
  }
  varcountRaw <- array(0, c(varWidth, n.samples, n.chains))

  for (chain in seq_len(n.chains)) {
    for (s in seq_len(n.samples)) {
      gamma <- colOf(r$cutpoints, s, chain)
      cutpointsRaw[, s, chain] <- gamma
      etaTrain <- colOf(r$train, s, chain)
      latentTrain[, s, chain] <- etaTrain
      probsTrain[,, s, chain] <-
        ordinalCategoryProbabilities(etaTrain, gamma)
      if (n.test > 0L) {
        etaTest <- colOf(r$test, s, chain)
        latentTest[, s, chain] <- etaTest
        probsTest[,, s, chain] <- ordinalCategoryProbabilities(etaTest, gamma)
      }
      varcountRaw[, s, chain] <- colOf(r$varcount, s, chain)
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
    result$fit <- sampler
    result$cutpoints.raw <- cutpointsRawFull
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

  varcount <- convertSamplesFromDbartsToBart(
    varcountRaw,
    n.chains,
    combineChains
  )
  predictorNames <- colnames(sampler$data@x)
  if (!is.null(predictorNames) && !is.null(dim(varcount))) {
    dimnames(varcount) <- if (length(dim(varcount)) > 2L) {
      list(NULL, NULL, predictorNames)
    } else {
      list(NULL, predictorNames)
    }
  }

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

# Survival-probability draws from an AFT fit (docs/design/survival.md).
# Under the log-normal model log T = f(x) + sigma eps, S(t | x) =
# 1 - Phi((log t - f(x)) / sigma), evaluated at every posterior draw of f
# and sigma. Returns draws per the package's three-tier convention
# (extract = draws, fitted = mean, ci.level = interval): users take means
# and quantiles over the draw margin themselves. newdata predicts out of
# sample (requires a fit kept with keepTrees); NULL uses the training fits.
survivalProbabilities.bart <- function(
  object,
  times,
  newdata = NULL,
  combineChains = TRUE,
  ...
) {
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

  burnInSigma <- NULL
  burnInK <- NULL
  if (nskip > 0L) {
    oldX.test <- sampler$data@x.test
    oldOffset.test <- sampler$data@offset.test

    oldKeepTrainingFits <- control@keepTrainingFits
    oldVerbose <- control@verbose

    if (length(x.test) > 0) {
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

    if (length(x.test) > 0) {
      sampler$setTestPredictorAndOffset(oldX.test, oldOffset.test)
    }
    control@keepTrainingFits <- oldKeepTrainingFits
    control@verbose <- oldVerbose
    if (keeptrees == TRUE) {
      control@keepTrees <- TRUE
    }
    sampler$setControl(control)

    samples <- sampler$run(0L, control@n.samples, updateState = FALSE)
  } else {
    samples <- sampler$run(updateState = FALSE)
  }

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
