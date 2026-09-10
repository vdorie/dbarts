# expects n.pars x n.samples x n.chains if n.chains > 1
# returns n.chains x n.samples x n.pars or (n.chains x n.samples) x n.pars
#
# for quantities such as yhat.train, n.pars is n.obs; for sigma, it is 1
# and the dimension is dropped
#
# when combineChains is TRUE and n.chains > 1, the chain margin collapses
# chain-major regardless of n.pars: chain 1's whole run, then chain 2's,
# and so on - so a combined draw index means the same thing for every
# field a fit returns (sigma, k, yhat.train, varcount, ...)
#
# preserves the per-parameter dimnames if they exist
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
        # chain-major: chain 1's whole run, then chain 2's, ... - the same
        # order the n.pars > 1 branch below produces, so a combined sigma[i]
        # pairs with a combined yhat.train[i, ]
        if (n.chains <= 1L) t(samples) else as.vector(samples),
        {
          x <- NULL
          # ONE permutation, not matrix() then t(). Rotating the parameter
          # margin to the end and flattening the rest with dim<- lands the
          # same (n.samples * n.chains) x n.pars layout the pair produced,
          # but through a single full-size allocation: dim<- rewrites an
          # attribute of an array nothing else references. On a channel the
          # size of yhat.train that second copy is the whole peak.
          res <- evalx(dim(samples), {
            permuted <- aperm(samples, c(seq_along(x)[-1L], 1L))
            dim(permuted) <- c(prod(x[-1L]), x[1L])
            permuted
          })
          if (!is.null(dimnames(samples))) {
            colnames(res) <- dimnames(samples)[[1L]]
          }
          res
        }
      )
    }
  }

# The posterior mean of a channel over its LAST 'trailing' margins - the
# reported ones (an observation, or an observation and a category), which the
# returned layout puts outermost, so each reported cell's draws are one
# contiguous slab. apply() would do the same reduction, but it aperms its
# argument first - even when the permutation is the identity - so on a
# channel the size of yhat.train it allocates a second full-size array while
# the engine's own is still bound. Reducing a slab at a time allocates the
# slab. The slab ORDER is load-bearing: mean() sums in the order it is handed
# and corrects over the same order, so reducing here rather than over the
# engine's own observation-major array is what keeps every value bit-for-bit
# what apply() returned. The result matches apply()'s shape too: a named vector
# over the last margin at trailing = 1, an array over the last 'trailing'
# margins otherwise.
channelMeans <- function(samples, trailing = 1L) {
  d <- dim(samples)
  numDims <- length(d)
  kept <- seq.int(numDims - trailing + 1L, numDims)
  blockLength <- prod(d[seq_len(numDims - trailing)])
  means <- vapply(
    seq_len(prod(d[kept])),
    function(i) {
      mean(samples[seq.int(i * blockLength - blockLength + 1, i * blockLength)])
    },
    0
  )
  if (trailing == 1L) {
    names(means) <- dimnames(samples)[[numDims]]
  } else {
    dim(means) <- d[kept]
    dimnames(means) <- dimnames(samples)[kept]
  }
  means
}

# input n.samples x n.chains x n.pars, or n.samples x n.pars when n.chains = 1
# output (n.samples * n.chains) x n.pars
combineChains <- function(samples) {
  ifelse_3(
    is.null(dim(samples)),
    length(dim(samples)) == 2L,
    samples,
    # chain-major, inverting uncombineChains' vector branch: samples is
    # (n.chains, n.samples), so transpose before flattening to get chain
    # 1's whole run first, then chain 2's, ...
    as.vector(t(samples)),
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
      # samples is chain-major (chain 1's whole run, then chain 2's, ...,
      # convertSamplesFromDbartsToBart's combined scalar-field order);
      # filling n.samples x n.chains column-major, then transposing, lands
      # each chain's contiguous block in its own row
      t(matrix(samples, length(samples) %/% n.chains, n.chains))
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

  # the channel's own presence, not fit$control@keepTrainingFits alone: the
  # bridge nulls samples$train whenever EITHER keepTrainingFits or keepFits
  # is FALSE, so testing the channel directly covers both switches without
  # restating their conjunction here
  yhat.train <- NULL
  yhat.train.mean <- NULL
  if (!is.null(samples$train)) {
    yhat.train <- convertSamplesFromDbartsToBart(
      samples$train,
      n.chains,
      combineChains
    )
    if (!responseIsBinary) {
      yhat.train.mean <- padOmittedRows(
        fit$data@na.action,
        channelMeans(yhat.train)
      )
    }
  }

  # likewise: NROW(fit$data@x.test) > 0 says a test set exists, but
  # samples$test is additionally null when keepFits dropped it, which
  # keepTrainingFits alone never did
  yhat.test <- NULL
  yhat.test.mean <- NULL
  if (NROW(fit$data@x.test) > 0 && !is.null(samples$test)) {
    yhat.test <- convertSamplesFromDbartsToBart(
      samples$test,
      n.chains,
      combineChains
    )
    if (!responseIsBinary) {
      yhat.test.mean <- channelMeans(yhat.test)
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
  # p x n.forests x n.samples (x n.chains), prognostic first); shapeMultinomialChannel
  # is the K-margin reshape that handles it, landing the forest axis on the
  # trailing margin exactly as multinomial's category axis does. The forest
  # count comes from the sampler rather than the array's rank, which is
  # ambiguous (a single-forest multi-chain channel is 3-D too); data@bases is
  # the same probe samplerCarriesAmplitudes uses.
  numForests <- if (!is.null(fit$data@bases)) length(fit$data@bases) else 1L
  forestNames <- if (numForests > 1L) {
    paste0("forest", seq_len(numForests))
  } else {
    NULL
  }
  varcount <- if (numForests > 1L) {
    shapeMultinomialChannel(
      samples$varcount,
      forestNames,
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

  # the per-forest in-sample channels: forestFits is each forest's own
  # RESPONSE-scale raw total (response.scale * f_k, no glue
  # folded in - $getForestFits' meaning up to that one scalar) and glue the
  # multiplier channel unchanged, so yhat = response.shift +
  # sum_k (basis_k %*% glue_k) * forestFits_k reconstructs train exactly. The
  # response transform is fixed at data creation and shared by every chain, so
  # one chain's reading of it is the whole fit's. Gated on the COUPLING, not
  # the forest count: forestReportingIsDefined() is true only for the
  # amplitude-coupled combiner, so a K-forest multinomial run has numForests
  # > 1 but samples$forestFits NULL (it packages through its own reshaper),
  # and samples$forestFits's presence is this R-side shadow of that predicate.
  forestFits <- NULL
  glue <- NULL
  hasForestReporting <- numForests > 1L && !is.null(samples$forestFits)
  if (hasForestReporting) {
    responseScale <- fit$getCalibration(1L)[1L, "response.scale"]
    forestFits <- shapeMultinomialChannel(
      samples$forestFits * responseScale,
      forestNames,
      n.chains,
      combineChains
    )
    glue <- convertSamplesFromDbartsToBart(
      samples$glue,
      n.chains,
      combineChains
    )
    # the ragged margin's forest-major width per forest (a basis-less forest's
    # implicit intercept is one column), the same layout $getForestAmplitudes()
    # documents
    widths <- vapply(
      fit$data@bases,
      function(basis) if (is.null(basis)) 1L else ncol(basis),
      integer(1L)
    )
    attr(glue, "forest") <- rep(forestNames, widths)
  }

  # heteroscedastic variance surface s(x) = sqrt(s^2(x)), train and test, on the
  # original scale; NULL for a homoscedastic fit
  #
  # hasVariance is whether the MODEL is heteroscedastic, which survives
  # keepFits = FALSE where !is.null(samples$variance) would not: the bridge
  # names the "variance" slot in samples whenever the model carries the
  # channel, even when keepFits nulls its VALUE, so a caller reading
  # object$hasVariance below can still tell a homoscedastic fit from a
  # heteroscedastic one whose s.train/s.test keepFits dropped - predict()'s
  # posterior-predictive guard needs exactly that (see predict.bart)
  hasVariance <- "variance" %in% names(samples)
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

  # the residual-distribution token every packaged fit carries: student()
  # residuals are refused unless family == "gaussian" (R/spec.R), so the
  # model's resid.df attribute (residDistDf, R/model.R) being non-NULL
  # already records which law was fit.
  residDist <- if (is.null(attr(fit$model, "resid.df"))) {
    "gaussian"
  } else {
    "student"
  }

  # what the fit's na.action dropped, in base R's own shape and under base
  # R's own name, so that fitted() and residuals() pad through naresid
  naOmitted <- fit$data@na.action

  if (responseIsBinary) {
    result <- list(
      call = fit$control@call,
      family = fit$model@family,
      resid.dist = residDist,
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
      resid.dist = residDist,
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

  # absent, not NULL, off a complete fit, so an ordinary fit object carries
  # exactly the elements it always did
  if (!is.null(naOmitted)) {
    result$na.action <- naOmitted
  }

  if (!is.null(varprobs)) {
    result$varprobs <- varprobs
  }

  # the per-draw residual degrees of freedom a student() fit conditions each
  # draw on - a fixed nu repeats, an estimated one is that sweep's grid draw -
  # in sigma's own layout, since it is one scalar per draw as sigma is. The
  # channel is absent (not NULL) off a Student-t error law, so a gaussian fit
  # keeps exactly the elements it had; pointwiseLogLikelihood reads it as the
  # t marginal's df.
  if (!is.null(samples[["resid.df"]])) {
    result[["resid.df"]] <- convertSamplesFromDbartsToBart(
      samples[["resid.df"]],
      n.chains,
      combineChains
    )
  }

  # the per-observation censoring status an aft (survival) fit needs so its
  # log-likelihood can take the survival tail on censored rows; the survival
  # ingestion parks it on the control attribute, and it is NULL (so the
  # element stays absent) for every non-survival family
  result$status <- attr(fit$control, "bartcore.survival")

  # the discrete-time hazard marker: the ordered period grid, present only
  # on hazard fits (parked on the control attribute by dbarts()).
  # survivalProbabilities dispatches its hazard branch on THIS marker
  # instead of $family. NULL - and so the element stays absent - for every
  # non-hazard fit.
  result$periods <- attr(fit$control, "bartcore.hazard.periods")

  # the forest count the widened varcount channel's trailing margin carries,
  # the multi-forest analog of the multinomial packager's K; absent at one
  # forest, so every single-forest fit keeps exactly the elements it had.
  # fitSynopsis reads it to tell a forest margin from a chain margin, which
  # the packaged rank alone cannot do
  if (numForests > 1L) {
    result$n.forests <- numForests
  }
  # absent (not FALSE) off a homoscedastic fit, the same "absent, not NULL"
  # convention the rest of this list uses
  if (hasVariance) {
    result$hasVariance <- TRUE
  }
  if (hasForestReporting) {
    result$forestFits <- forestFits
    result$glue <- glue
    # the expanded per-forest bases, no draw axis: what makes the
    # reconstruction identity evaluable from the fit alone, since
    # neither run() channel carries them
    result$bases <- fit$data@bases
    forestInfo <- attr(fit$control, "bartcore.forests")
    # the declaring formula and the fit-time levels of every basis written as a
    # forest() term, which is what lets the same basis be rebuilt at NEW rows;
    # the element stays absent on a fit whose bases arrived as values, whose
    # blend then needs the caller's own at those rows
    result$basis.terms <- forestInfo$basisTerms
    forestLabels <- forestInfo$labels
    if (!is.null(forestLabels)) {
      attr(result, "forest.labels") <- forestLabels
    }
  }

  if (keepSampler) {
    result$fit <- fit
  }
  # the Gaussian-process fallback census of the run that produced this fit:
  # the leaf evaluations that reached a GP entry point and the ones that found
  # the leaf over max.leaf.size and scored it as a constant leaf instead. NULL
  # (and so dropped below) on every other leaf model. The warning that fires
  # above a quarter is raised at the run, not here; these are the counts behind
  # it, so any other threshold can be checked by hand.
  result$gp.fallback <- attr(samples, "gp.fallback")
  result$n.chains <- n.chains
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

  # a component absent from this fit (no test set, no heteroscedastic
  # variance, ...) was still written into the list above as an explicit
  # NULL, which names(fit) then lists as present; drop those, order-
  # preserving, rather than reporting a component the fit does not have.
  # '[' on a plain list keeps only 'names' of the attributes present, so
  # "forest.labels" (attached above, when present) is carried across by hand
  extraAttrs <- attributes(result)[
    setdiff(names(attributes(result)), "names")
  ]
  result <- result[!vapply(result, is.null, NA)]
  attributes(result)[names(extraAttrs)] <- extraAttrs

  class(result) <- "bart"
  invisible(result)
}

## Builds the quoted tree/node/resid prior calls the entry points hand to
## dbarts. nodeK is the node prior's k argument exactly as
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

  # A caller-supplied tree.prior/node.prior/resid.prior object fully replaces
  # the flat build below and is forwarded UNEVALUATED, exactly as k already
  # is (nodeK), so a bare vocabulary name inside it (linear(), gp(), fixed(),
  # ...) resolves in the caller's own frame, not here. A shorthand that would
  # otherwise help build the same prior is a collision, refused by name
  # before anything is built. An explicit `tree.prior = NULL` is itself NULL
  # here (match.call() stores a literal NULL as NULL), indistinguishable
  # from not supplying one at all.
  treePriorObj <- matchedCall[["tree.prior"]]
  nodePriorObj <- matchedCall[["node.prior"]]
  residPriorObj <- matchedCall[["resid.prior"]]

  if (!is.null(treePriorObj)) {
    refuseColliding(
      matchedCall,
      "tree.prior",
      c("power", "base", splitProbsName, "dart")
    )
    tree.prior <- treePriorObj
  } else {
    tree.prior <- resolveDartShorthand(
      dart,
      splitProbsName %in% names(matchedCall),
      splitProbsName,
      function() {
        priorCall <- quote(dart(power, base))
        priorCall[[2L]] <- power
        priorCall[[3L]] <- base
        priorCall
      },
      function() {
        priorCall <- quote(cgm(power, base, split.probs))
        priorCall[[2L]] <- power
        priorCall[[3L]] <- base
        priorCall[[4L]] <- if (splitProbsName %in% names(matchedCall)) {
          matchedCall[[splitProbsName]]
        } else {
          splitProbsDefault
        }
        priorCall
      }
    )
  }

  if (!is.null(nodePriorObj)) {
    refuseColliding(matchedCall, "node.prior", c("k", "prior.scale"))
    node.prior <- nodePriorObj
    # a named prior.scale needs a node prior to ride even when k is left to the
    # family default, so it builds one; the k slot then drops out of the call and
    # dbarts() resolves k exactly as it would have with no node prior at all
  } else if (!is.null(nodeK) || !is.na(priorScale)) {
    node.prior <- quote(normal(k))
    node.prior[[2L]] <- nodeK
    if (!is.na(priorScale)) {
      node.prior[["scale"]] <- priorScale
    }
  } else {
    node.prior <- NULL
  }

  if (!is.null(residPriorObj)) {
    refuseColliding(
      matchedCall,
      "resid.prior",
      c("sigdf", "sigquant", "sigest")
    )
    resid.prior <- residPriorObj
  } else {
    resid.prior <- quote(chisq(sigdf, sigquant))
    resid.prior[[2L]] <- sigdf
    resid.prior[[3L]] <- sigquant
  }

  list(
    tree.prior = tree.prior,
    node.prior = node.prior,
    resid.prior = resid.prior
  )
}

# Builds the dbarts() host-sampler call shared by bart2's standard path and its
# alternate-family paths (multinomial/ordinal/nbinom): redirect the matched
# call to dbarts(), install the resolved control, drop n.samples (the sampler
# is driven by run() rather than sampled at construction), and thread the
# prebuilt tree/node/resid priors. 'family', when supplied, overrides the
# forwarded family token - NULL removes it, a string sets it (dbarts() accepts
# "multinomial", "ordinal", and "nbinom" as family tokens directly); left
# missing it keeps the user's forwarded family (the normal path). 'sigest',
# when supplied, sets the sampler's sigma; left missing it is not touched
# (ordinal/nbinom read no sigma). The multinomial call sites set
# samplerCall$data to the real response (a factor or count matrix, which
# dbarts()'s own resolveMultinomialCounts one-hot expands) after this call
# returns; it is not a placeholder.
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
    # the caller's own family object is already stamped on the matched call
    # and carries that family's settings; an override that names the same
    # family keeps it rather than flattening it back to a bare token
    stamped <- matchedCall$family
    samplerCall$family <- if (
      is(stamped, "dbartsFamily") && identical(stamped@token, family)
    ) {
      stamped
    } else {
      family
    }
  }
  if (!missing(sigest)) {
    samplerCall$sigest <- as.numeric(sigest)
  }
  samplerCall
}

# The two-phase burn-in/sample split shared by both doors' standard path.
# When burn-in is requested the test set is dropped and keepTrainingFits/verbose
# forced FALSE for the burn run (a draw-neutral speedup - test fits do not feed
# the MCMC), then restored for the kept-sample run, with keepTrees re-enabled
# after burn when requested. Returns the kept samples plus the burn-in
# sigma/k channels.
#
# 'callback' is installed on the KEPT-sample run only. Every burn-in sweep is
# a recorded draw at this layer (the burn phase is its own sampler$run(...)
# call, so nothing here tells the engine those draws are burn-in), so a hook
# installed for both calls would see - and average - burn-in draws into
# whatever it accumulates; the engine's own rule ("fires on saved draws") is
# left exactly as stated, the split just never hands it the burn call.
runWithBurnIn <- function(sampler, control, keepTrees, callback = NULL) {
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

    samples <- sampler$run(
      0L,
      control@n.samples,
      updateState = FALSE,
      callback = callback
    )
  } else {
    samples <- sampler$run(updateState = FALSE, callback = callback)
  }
  list(samples = samples, burnInSigma = burnInSigma, burnInK = burnInK)
}

# The alternate-family bart2 arcs (multinomial/ordinal/nbinom/hurdle.lognormal)
# all refuse warm.start/n.grow.sweeps, keepTrainingFits = FALSE, and
# keepFits = FALSE, differing only in the family name and the reason (which
# completes "requires keepTrainingFits = TRUE (the default): " /
# "requires keepFits = TRUE (the default): "). samplerOnly is refused by
# default too, but allow.samplerOnly lets a caller whose returned sampler is
# load-bearing (ordinal, nbinom) opt back in; hurdle.lognormal's $fit is a
# PAIR of samplers, so it stays refused.
#
# keepFits = FALSE reaches here two ways: the caller wrote it, or 'callback'
# set it AUTOMATICALLY (bart()'s own default expression). Both are refused
# alike - these families' packaging reads the training channel keepFits
# drops - but the message names 'callback' too, since the automatic path is
# the one a caller is least likely to expect.
checkFamilyUnsupportedArgs <- function(
  family,
  samplerOnly,
  warm.start,
  n.grow.sweeps,
  control,
  reason,
  allow.samplerOnly = FALSE
) {
  if (isTRUE(samplerOnly) && !allow.samplerOnly) {
    stop("family = \"", family, "\" does not support 'samplerOnly'")
  }
  grownSweeps <- coerceOrError(n.grow.sweeps, "integer")[1L]
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
  if (!control@keepFits) {
    stop(
      "family = \"",
      family,
      "\" requires keepFits = TRUE (the default; a supplied 'callback' ",
      "sets it FALSE automatically unless 'keepFits' is given explicitly): ",
      reason
    )
  }
}

bart <- function(
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
  split.probs = NULL,
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
  proposal.probs = c(
    birth_death = 0.6,
    swap = 0,
    change = 0.4,
    perturb = 0,
    rule_gibbs = 0,
    birth = 0.5
  ),
  monotone = NULL,
  interactions = NULL,
  blocks = NULL,
  variance = NULL,
  keepSampler = keepTrees,
  warm.start = NULL,
  n.grow.sweeps = 0L,
  factors = c("categorical", "indicators"),
  family = c(
    "auto",
    "gaussian",
    "student",
    "probit",
    "logistic",
    "aft",
    "multinomial",
    "ordinal",
    "nbinom",
    "hazard",
    "hazard.probit",
    "hazard.logistic",
    "hurdle.lognormal"
  ),
  na.action = dbarts::na.keepPredictors,
  tree.prior = NULL,
  node.prior = NULL,
  resid.prior = NULL,
  storage = c("double", "single"),
  updateState = TRUE,
  keepFits = is.null(callback),
  callback = NULL,
  ...
) {
  matchedCall <- match.call()
  callingEnv <- parent.frame()

  # Two doors, one name for one release. A call carrying a BayesTree
  # spelling is a 0.9-x call and is fit by the legacy door with 0.9-x's
  # defaults; everything else is fit here. Both checks precede every other
  # read of the arguments, so a forwarded call is never partly interpreted
  # under the modern vocabulary first.
  suppliedCall <- sys.call()
  forwarded <- forwardToLegacyDoor(
    suppliedCall,
    names(matchedCall)[-1L],
    callingEnv
  )
  if (!is.null(forwarded)) {
    return(forwarded$value)
  }
  refuseLegacyPositionalCall(suppliedCall)
  # the dots are inspected by NAME, never forced as a list: a retired
  # argument may be written in a vocabulary only this package holds
  # (resid.dist = student()), which would not resolve in the caller's frame
  supplied <- dotNames(...)
  refuseForeignFrontDoorArgs(supplied, "bart", names(formals(dbarts::bart)))
  # ahead of sampler construction below, so a malformed pair fails here
  # rather than after the (possibly expensive) sampler is already built
  validateCallback(callback)
  # the names dec-B98's consolidation moved onto the family and prior
  # objects: read once, then cleared from the matched call so no forwarding
  # can carry an old spelling on to dbarts()
  consolidated <- resolveConsolidatedArgs(
    matchedCall,
    supplied,
    "bart",
    callingEnv
  )
  if (length(consolidated) > 0L) {
    matchedCall[names(consolidated)] <- NULL
  }
  dart <- if (is.null(consolidated[["dart"]])) FALSE else consolidated[["dart"]]
  levelGibbs <- consolidated[["levelGibbs"]]
  # the collision the tree-prior shorthand ladder would otherwise catch by
  # name; 'dart' has already left the matched call it reads
  if (!is.null(matchedCall[["tree.prior"]]) && !isFALSE(dart)) {
    stop(
      "'tree.prior' cannot be combined with 'dart': supply the prior either ",
      "as an object or through its shorthand arguments, not both"
    )
  }
  if ("rngSeed" %in% supplied) {
    seed <- resolveRenamedSeed(
      eval(matchedCall$rngSeed, callingEnv),
      "bart",
      seed
    )
    # stamped onto the matched call as well: the control is built by
    # redirecting that call by name, which the old spelling never reaches
    matchedCall$rngSeed <- NULL
    matchedCall$seed <- seed
  }

  # 'family' is resolved from the caller's own unevaluated argument, so a
  # bare family constructor (student(3), hazard(breaks)) resolves in the
  # family vocabulary; nothing may force the argument before this. The
  # resolved object is stamped back on, so every forwarding below - the host
  # dbarts() call, the hurdle components - carries the settings rather than
  # a spelling that has to be re-resolved against a different environment.
  familySpec <- resolveFamily(
    matchedCall$family,
    eval(formals(dbarts::bart)$family),
    "bart",
    callingEnv
  )
  familySpec <- applyConsolidatedFamilyArgs(familySpec, consolidated)
  family <- familySpec@token
  # the caller's own 'family' expression, kept for the stored call: the
  # resolved object below is a forwarding detail, and a fit whose call was
  # never written with a family should not print one
  suppliedFamily <- matchedCall$family
  matchedCall$family <- familySpec

  # A data object carrying an n x K count matrix declares the multinomial
  # (softmax) model, whose fitted quantity is K probabilities per observation
  # rather than one location. Everything this function does after the fit -
  # the chain/sample reshaping, the sigma channel, the loss and interval
  # readers - is written against one location, so the slab would be reshaped
  # as if it were n rows and reported without a word. Refuse the object here,
  # ahead of every family branch, since family = "auto" resolves it downstream
  # and no branch below would see the token.
  refuseCountsCarryingData(formula, "bart()")
  refuseResponseFreeFormula(formula, "bart()")

  # family = "auto" with a 3+-level UNORDERED factor/character response is
  # multinomial; a 3+-level ORDERED factor is ordinal (the disjoint
  # is.ordered() key). Route
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

  # factors/missing/proposal.probs are forwarded formal defaults:
  # redirectCall only carries a name into the host dbarts() call when the
  # caller supplied it, so an unsupplied one must be resolved here, in
  # this function's own frame, and stamped onto matchedCall unconditionally,
  # or it would silently take dbarts()'s own default rather than the
  # token/value this signature advertises.
  factors <- match.arg(factors)
  matchedCall$factors <- factors
  matchedCall$proposal.probs <- proposal.probs

  controlCall <- redirectCall(matchedCall, dbarts::dbartsControl)
  # '...' is a formal of both and names no value, so it is excluded before
  # the shared names are evaluated as defaults
  sharedFormals <- setdiff(names(formals(dbarts::bart)), "...")
  missingDefaultArgs <- sharedFormals[
    sharedFormals %in%
      names(formals(dbarts::dbartsControl)) &
      sharedFormals %not_in% names(matchedCall)
  ]
  if (length(missingDefaultArgs) > 0L) {
    currentEnv <- sys.frame(sys.nframe())
    controlCall[missingDefaultArgs] <- lapply(
      formals(dbarts::bart)[missingDefaultArgs],
      eval,
      envir = currentEnv
    )
  }
  control <- eval(controlCall, envir = callingEnv)
  # the level Gibbs step is declared on the tree prior now; the retired
  # spelling lands where a declared one lands, on the control the bridge
  # reads (R/spec.R copies the prior's own declaration there)
  if (!is.null(levelGibbs)) {
    control@levelGibbs <- validateLevelGibbs(levelGibbs)
  }

  storedCall <- matchedCall
  storedCall$family <- suppliedFamily
  control@call <- if (keepCall) storedCall else call("NULL")
  control@n.burn <- control@n.burn %/% control@n.thin
  control@n.samples <- control@n.samples %/% control@n.thin
  # printEvery counts post-thinning samples and must stay positive: thinning
  # heavier than it would otherwise drive it to 0, which the sampler refuses
  control@printEvery <- max(1L, control@printEvery %/% control@n.thin)

  # the multinomial branch below keeps its own family-named check
  if (family != "multinomial" && isTRUE(control@n.samples <= 0L)) {
    refuseZeroSamples("bart")
  }

  # multinomial: a K-forest softmax model over a factor response or an n x K
  # count matrix, validated and dispatched here rather than threaded through
  # the rest of this function - it bypasses the standard single-forest
  # packaging entirely, though bart2Multinomial/bart2MultinomialCounts
  # build their
  # sampler through dbarts()'s own family = "multinomial" dispatch, the same
  # one bartcoreMultinomialSampler/bartcoreMultinomialCountSampler route
  # through. Every refusal below names the limitation rather than silently
  # reshaping around it.
  if (family == "multinomial") {
    # a K-forest softmax has no amplitude-coupled slot for a forest() term to
    # declare; caught here since multinomial never reaches the shared
    # dbarts() ingestion (R/formulaTerms.R) that catches every other family
    if (formulaHasForestTerm(formula)) {
      stop(
        "family = \"multinomial\" does not support a forest() formula term"
      )
    }
    if (!missing(weights)) {
      stop(
        "family = \"multinomial\" does not support 'weights': a non-integer ",
        "weight has no exact augmentation sampler, and an integer weight is ",
        "already expressible as row-wise count replication in the response"
      )
    }
    # the n x K category offset: a matrix enters the raw fits before the
    # softmax and is threaded to
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
        refuseFlatOffsetOnMultinomial(offset)
      }
      if (!is.numeric(offset)) {
        stop("family = \"multinomial\" 'offset' must be a numeric n x K matrix")
      }
      multinomialOffset <- offset
      matchedCall$offset <- NULL
    }
    # offset is train-side only: an explicit offset.test is refused here
    # rather than left to fall through to parseMultinomialData's later
    # refusal. yhat.test is always computed WITHOUT any category offset, so
    # a caller wanting one needs the sampler-level channel this message names.
    if (!missing(offset.test)) {
      stop(
        "family = \"multinomial\" does not support 'offset.test'; yhat.test ",
        "is always computed without any category offset. A category test ",
        "offset is a sampler-level capability only - a dbartsSampler's own ",
        "$setCategoryTestOffset method, reached through keepSampler = TRUE, ",
        "or the internal creators' own offset.test argument - not reachable ",
        "from bart()"
      )
    }
    if (!missing(subset)) {
      stop(
        "family = \"multinomial\" does not support 'subset'"
      )
    }
    # buildMultinomialForest fixes its own CGM tree prior and copies only
    # power/base/proposal-probability fields from the host sampler it briefly
    # builds, so none of the following reach the K-forest engine: a DART
    # flag on either formal ('tree.prior' is resolved with the same prior
    # vocabulary parsePriors, R/model.R, uses, so 'tree.prior = dart()' is
    # caught the same way 'tree.prior = dbartsPriors$dart()' is), fixed
    # split probabilities, the monotone direction constraints (only their
    # proposal-probability rewrite is copied), and a variance forest (the
    # host sampler that would resolve one is discarded).
    treePrior <- if (missing(tree.prior)) {
      NULL
    } else {
      eval(
        matchedCall[["tree.prior"]],
        list2env(dbartsPriors, parent = callingEnv)
      )
    }
    unsupported <- c(
      "'dart' or a DART 'tree.prior'" = !isFALSE(dart) ||
        inherits(treePrior, "dbartsDartPrior"),
      "'split.probs'" = !is.null(split.probs),
      "'monotone'" = !is.null(monotone),
      "'variance'" = !is.null(variance)
    )
    if (any(unsupported)) {
      stop(
        "family = \"multinomial\" does not support ",
        paste0(names(unsupported)[unsupported], collapse = ", "),
        ": the K-forest engine copies only power/base/proposal-probability ",
        "fields from the host sampler it briefly builds"
      )
    }
    checkFamilyUnsupportedArgs(
      "multinomial",
      samplerOnly,
      warm.start,
      n.grow.sweeps,
      control,
      "there is no test surface to fall back on",
      allow.samplerOnly = TRUE
    )
    warnFamilyGatedArgs(argNames, "multinomial")
    if (control@n.samples <= 0L) {
      stop("family = \"multinomial\" requires a positive 'n.samples'")
    }
    if (inherits(formula, "dbartsData")) {
      stop(
        "family = \"multinomial\" does not support a pre-built dbartsData ",
        "object; use the formula interface or the matrix ",
        "interface: bart(x, y, family = \"multinomial\")"
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
    # the count-matrix response form: an n x K matrix of nonnegative
    # integer trial counts, beside the factor path
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
        stop("family = \"multinomial\" count response must be non-negative")
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
        keepSampler = keepSampler,
        samplerOnly = samplerOnly
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
      keepSampler = keepSampler,
      samplerOnly = samplerOnly
    ))
  }

  # ordinal (cumulative probit): a single-forest fixed-scale latent model
  # whose n x K category probabilities and per-draw thresholds
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
      "the category probabilities are built from the training latent fits",
      allow.samplerOnly = TRUE
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
      keepSampler = keepSampler,
      samplerOnly = samplerOnly
    ))
  }

  # negative-binomial counts: a single-forest fixed-scale latent model whose
  # mean counts mu = r exp(f + o) and per-draw
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
      "the mean counts are built from the training latent fits",
      allow.samplerOnly = TRUE
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
      keepSampler = keepSampler,
      samplerOnly = samplerOnly
    ))
  }

  # hurdle.lognormal: a semicontinuous two-part fit built from an
  # occupancy probit on 1{y > 0} (all n) and a gaussian
  # on log(y) restricted to the y > 0 subset, glued at report time. Dispatched
  # here - not inside dbarts(), which returns a single sampler and cannot
  # express the two-sampler composition - so its bartHurdle fit object stays
  # distinct. family is always explicit: a semicontinuous response has no
  # unambiguous auto class.
  if (family == "hurdle.lognormal") {
    # the two-sampler composition has no amplitude-coupled slot for a forest()
    # term either; caught here for the same reason as the multinomial guard
    # above, and before the matrix-interface-only check below would otherwise
    # name neither the family's real reason nor the term
    if (formulaHasForestTerm(formula)) {
      stop(
        "family = \"hurdle.lognormal\" does not support a forest() formula ",
        "term"
      )
    }
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

  # k enters unevaluated: bart redirects its matched call, so the stored
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
    splitProbsDefault = formals(dbarts::bart)[["split.probs"]]
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
  n.grow.sweeps <- coerceOrError(n.grow.sweeps, "integer")[1L]
  if (is.na(n.grow.sweeps) || n.grow.sweeps < 0L) {
    stop("'n.grow.sweeps' must be a non-negative integer")
  }
  if (!is.null(warm.start) && n.grow.sweeps > 0L) {
    stop(
      "'warm.start' and 'n.grow.sweeps' both request an initialization; ",
      "supply at most one"
    )
  }
  # a data object carrying forest bases reaches here as an ordinary fit, so
  # this is where the multi-forest donor refusal is owed on the modelling
  # surface: named for the argument the caller wrote, ahead of the R5 method's
  # own wording and the bridge backstop under it. 'n.grow.sweeps' is not held
  # to it - grow-from-root composes through the combiner and is covered at two
  # forests
  if (!is.null(warm.start)) {
    refuseMultiForestWarmStart(sampler$getPointer(), "'warm.start'")
  }

  if (!is.null(warm.start)) {
    sampler$installTrees(warm.start)
  } else if (n.grow.sweeps > 0L) {
    sampler$growFromRoot(n.grow.sweeps, updateState = FALSE)
  } else {
    sampler$sampleTreesFromPrior(updateState = FALSE)
  }

  burn <- runWithBurnIn(sampler, control, keepTrees, callback)
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
# validateXTest code a data.frame newdata.
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

  # a sparseVector/dgCMatrix/sparseFactor predictor dies inside model.frame
  # with a bare S4 type error; refuse it explicitly first, as the standard
  # formula path does
  if (is.list(data) || is.environment(data)) {
    refuseSparseFormulaColumns(formula, data)
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
# response selects multinomial (an ORDERED factor is ordinal, split out here
# and caught by detectAutoOrdinal instead - the disjoint is.ordered() key).
# Returns classifyResponse's descriptor for that case,
# NULL for anything that stays on the standard single-forest path (numeric,
# 2-level, logical, ordered, a count matrix, or a pre-built dbartsData). Type
# detection only - the multinomial branch re-extracts and validates the
# response, so this evaluates the formula LHS directly rather than building a
# model frame.
# A dbartsData carrying the n x K count response is a multinomial
# declaration, and the single-forest packaging of the entry
# points below cannot read a K-location fit. The refusal is at the ENTRY, not
# at a family branch: family = "auto" resolves counts to multinomial during
# spec resolution, so a branch keyed on the token never sees it.
refuseCountsCarryingData <- function(formula, caller) {
  if (inherits(formula, "dbartsData") && !is.null(dataCounts(formula))) {
    stop(
      caller,
      " does not support a data object carrying an n x K count matrix; fit ",
      "it with dbarts(family = \"multinomial\"), whose sampler reports the K ",
      "category probabilities"
    )
  }
  invisible(NULL)
}

# A one-sided formula (~ x1 + x2) names no response. dbartsData() defaults a
# missing response to zeros, which is correct for a composed sampler whose
# response is set later, but a fitting entry point runs the chain and returns
# that all-zero fit as though it were a real answer, with no error or
# warning. Refuse it here, by name, before the formula reaches dbartsData().
refuseResponseFreeFormula <- function(formula, caller) {
  if (is.formula(formula) && length(formula) != 3L) {
    stop(caller, " requires a two-sided formula; 'formula' has no response")
  }
  invisible(NULL)
}

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
# response selects family = "ordinal". Only is.ordered() responses match,
# so ordinal and unordered multinomial never
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

# The multinomial (softmax) fit path, reached from bart2's family =
# "multinomial" branch after ingestion validation. y is the validated factor
# response; K = nlevels(y) follows from it. The sampler is built directly
# through dbarts()'s own family = "multinomial" dispatch, the same public
# construction path dbarts(x, y, family = "multinomial") reaches, so bart2's
# usual tree.prior/node.prior/resid.prior/control machinery resolves
# n.trees, n.chains, the tree prior and k exactly as it would for any other
# family. No warm start and no two-phase burn-in/sample split: both are
# skipped, as bartcoreRun's single call needs neither. offset, when
# non-NULL, is the validated n x K category offset (bart2's own matrix-only
# surface, above); it is installed on the constructed sampler through
# $setCategoryOffset, before the run - never through the host dbarts()
# call's own flat offset argument.
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
  keepSampler = FALSE,
  samplerOnly = FALSE
) {
  K <- nlevels(y)

  priors <- buildSamplerPriors(
    matchedCall,
    power,
    base,
    sigdf,
    sigquant,
    nodeK = matchedCall[["k"]],
    priorScale = prior.scale,
    dart = dart,
    splitProbsDefault = formals(dbarts::bart)[["split.probs"]]
  )

  # one bartcore_create, through the public multinomial dispatch: the factor
  # response IS the count-matrix declaration (dbarts()'s own
  # resolveMultinomialCounts one-hot expands it against the response y is
  # already validated to be)
  samplerCall <- buildHostSamplerCall(
    matchedCall,
    control,
    priors,
    family = "multinomial",
    sigest = sigest
  )
  samplerCall$data <- y

  sampler <- eval(samplerCall, envir = callingEnv)
  if (!is.null(offset)) {
    sampler$setCategoryOffset(offset)
  }
  if (isTRUE(samplerOnly)) {
    return(sampler)
  }

  samples <- bartcoreRun(
    list(ptr = sampler$getPointer()),
    control@n.burn,
    control@n.samples
  )

  result <- packageMultinomialResults(
    control,
    y,
    levels(y),
    K,
    samples,
    combineChains,
    predictorNames = colnames(sampler$data@x)
  )
  # keepTrees retains the saved trees predict.bartMultinomial replays
  # through (the sampling sweeps wrote them regardless), and the sampler's
  # coded design (sampler@data@x) codes newdata to the training columns.
  # keepSampler retains $fit on its own, independent of keepTrees, so a
  # caller can see the sampler without paying for the K forests' saved trees.
  if (control@keepTrees || keepSampler) {
    result$fit <- sampler
  }
  result
}

# The grouped-count analog of bart2Multinomial: y is the validated n x K
# count matrix (bart2's count-matrix branch above) and levels are the
# resolved category names (colnames(y), or the index fallback). Mirrors
# bart2Multinomial's direct construction exactly - same priors, same one
# bartcore_create through dbarts()'s own family = "multinomial" dispatch -
# with y itself, already an n x K matrix, as the response:
# resolveMultinomialCounts passes a matrix straight through. A one-hot y
# with every row sum 1 is therefore the same draw stream as bart2Multinomial
# on the equivalent factor (the single-trial reduction). offset, as for
# bart2Multinomial, is installed through $setCategoryOffset.
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
  keepSampler = FALSE,
  samplerOnly = FALSE
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
    splitProbsDefault = formals(dbarts::bart)[["split.probs"]]
  )

  samplerCall <- buildHostSamplerCall(
    matchedCall,
    control,
    priors,
    family = "multinomial",
    sigest = sigest
  )
  samplerCall$data <- y

  sampler <- eval(samplerCall, envir = callingEnv)
  if (!is.null(offset)) {
    sampler$setCategoryOffset(offset)
  }
  if (isTRUE(samplerOnly)) {
    return(sampler)
  }

  samples <- bartcoreRun(
    list(ptr = sampler$getPointer()),
    control@n.burn,
    control@n.samples
  )

  result <- packageMultinomialResults(
    control,
    y,
    levels,
    K,
    samples,
    combineChains,
    predictorNames = colnames(sampler$data@x)
  )
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

# Re-derives a K-margined packaged channel's requested combineChains shape
# from however it is currently stored: forestFits' trailing n/forest pair
# (trailing = 2) or glue's trailing ragged margin (trailing = 1), generalizing
# combineOrUncombineChains (a single trailing margin) to what
# shapeMultinomialChannel widens. Used by extract(type = "forest") to serve
# either convention regardless of the fit's own packaged combineChains, the
# same on-demand reshape yhat.train already gets there. n.chains <= 1 is a
# no-op, as it is for every other channel: one chain carries no separate axis
# to fold or split. Round-trips bitwise with shapeMultinomialChannel's own
# combine/uncombine forms.
reshapeChainedChannel <- function(x, n.chains, combine, trailing) {
  d <- dim(x)
  storedCombined <- length(d) == trailing + 1L
  if (n.chains <= 1L || storedCombined == combine) {
    return(x)
  }
  dn <- dimnames(x)
  trailIdx <- seq_len(trailing) + 2L
  if (combine) {
    a <- aperm(x, c(2L, 1L, trailIdx))
    dim(a) <- c(d[1L] * d[2L], d[trailIdx])
    if (!is.null(dn)) {
      dimnames(a) <- c(list(NULL), dn[trailIdx])
    }
  } else {
    a <- array(x, c(d[1L] %/% n.chains, n.chains, d[trailIdx - 1L]))
    a <- aperm(a, c(2L, 1L, trailIdx))
    if (!is.null(dn)) {
      dimnames(a) <- c(list(NULL, NULL), dn[trailIdx - 1L])
    }
  }
  a
}

# Reshapes one bartcoreRun() result into a bart2(family = "multinomial") fit.
# samples$train is n.obs x K x n.samples (x n.chains); the K-carrying softmax
# probabilities are already the identified quantity the engine reports (the
# run's train channel is the identified deliverable), so no
# probabilityFromLatents-style transform applies here, unlike the binary
# families. Reshaped to the package's
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

# The cumulative-probit category probabilities for one posterior draw: given
# the latent means eta (length n) and the K-1 finite thresholds gamma
# (gamma_1 < ... < gamma_{K-1}), returns the
# n x K matrix P(y = k) = Phi(gamma_k - eta) - Phi(gamma_{k-1} - eta), with
# gamma_0 = -Inf and gamma_K = +Inf supplying the boundary columns of 0 and 1.
# Shared by the fit-time reshape and predict.bartOrdinal's replay.
ordinalCategoryProbabilities <- function(eta, thresholds) {
  n <- length(eta)
  K <- length(thresholds) + 1L
  bounds <- matrix(0, n, K + 1L)
  bounds[, K + 1L] <- 1
  for (j in seq_len(K - 1L)) {
    bounds[, j + 1L] <- pnorm(thresholds[j] - eta)
  }
  bounds[, 2L:(K + 1L), drop = FALSE] - bounds[, 1L:K, drop = FALSE]
}

# The ordinal (cumulative probit) fit path, reached from bart2's family =
# "ordinal" branch. Unlike multinomial's K-forest
# engine, ordinal is a SINGLE forest whose kept samples report the latent
# eta = f(x) (like probit) and, in an ordinal-only run channel, the per-draw
# K-1 thresholds; the n x K category probabilities are synthesized here from the
# two. A single run(n.burn, n.samples) drives the whole chain (the threshold
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
  keepSampler = FALSE,
  samplerOnly = FALSE
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
    splitProbsDefault = formals(dbarts::bart)[["split.probs"]]
  )

  samplerCall <- buildHostSamplerCall(
    matchedCall,
    control,
    priors,
    family = "ordinal"
  )

  sampler <- eval(samplerCall, envir = callingEnv)
  if (isTRUE(samplerOnly)) {
    return(sampler)
  }

  K <- attr(sampler$control, "bartcore.n.categories")
  levels <- sampler$data@response.levels
  n.chains <- control@n.chains
  n.obs <- length(sampler$data@y)
  n.test <- NROW(sampler$data@x.test)
  n.samples <- control@n.samples

  bc <- bartcoreSampler(sampler, family = "ordinal")
  # bc's engine is the one that runs below; adopt it into sampler so $fit
  # becomes an R5 wrapper around the engine that actually ran, not the
  # abandoned first-created host.
  sampler$adoptPointer(bc$ptr)

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
  thresholdsRaw <- array(0, c(K - 1L, n.samples, n.chains))

  # one run drives the whole chain; the thresholds ride their own run channel
  # (r$thresholds, present because the ordinal family carries them), aligned with
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
      gamma <- channelColumn(r$thresholds, s, chain, n.chains)
      thresholdsRaw[, s, chain] <- gamma
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

  # keep the full 3-D (K-1) x n.samples x n.chains thresholds for predict, which
  # pairs each replayed latent draw with its own thresholds
  thresholdsRawFull <- thresholdsRaw

  # drop the trailing singleton chain margin so the reshapers see the
  # n.chains == 1 layout their multinomial/gaussian siblings emit
  if (n.chains == 1L) {
    probsTrain <- array(probsTrain, dim(probsTrain)[1:3])
    latentTrain <- matrix(latentTrain, n.obs, n.samples)
    if (!is.null(probsTest)) {
      probsTest <- array(probsTest, dim(probsTest)[1:3])
      latentTest <- matrix(latentTest, n.test, n.samples)
    }
    thresholdsRaw <- matrix(thresholdsRaw, K - 1L, n.samples)
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
    thresholdsRaw,
    varcountRaw,
    combineChains
  )
  # keepTrees retains the saved trees predict.bartOrdinal replays through (the
  # sweeps wrote them regardless), and the sampler codes newdata to the
  # training columns. thresholds.raw supplies predict's per-draw thresholds in
  # the raw (K-1) x n.samples x n.chains layout that pairs with the replayed
  # latent draws, so no re-run is needed.
  if (control@keepTrees) {
    result$thresholds.raw <- thresholdsRawFull
  }
  if (control@keepTrees || keepSampler) {
    result$fit <- sampler
  }
  result
}

# Assemble a bart2(family = "ordinal") fit from the synthesized channels.
# yhat.train/test are the n x K category probabilities in the multinomial
# draws-first convention (levels named on the
# K margin); thresholds is the per-draw K-1 vector, the ordinal analog of
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
  thresholdsRaw,
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
    thresholds = convertSamplesFromDbartsToBart(
      thresholdsRaw,
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

# The negative-binomial mean counts for one posterior draw: given the
# log-odds latent psi = f(x) + o and the dispersion r, the mean is
# mu_i = r exp(psi_i). Shared by the fit-time
# reshape and predict.bartNegbin's replay.
negbinMeanCounts <- function(psi, r) {
  r * exp(psi)
}

# One posterior-predictive count per (draw, observation) of a negative-binomial
# fit: y ~ NB(size = r, mean = mu) with mu the draw's mean count and r its
# dispersion. mu is a (chains x) draws x obs array; r is the draws-shaped
# dispersion, broadcast across the observation
# margin so each column shares its draw's r. Only this and predict's ppd touch
# the RNG, so type = "ev"/"bart" stay draw-neutral.
negbinPpd <- function(mu, r) {
  d <- dim(mu)
  array(
    rnbinom(length(mu), size = as.vector(array(r, d)), mu = as.vector(mu)),
    d
  )
}

# The negative-binomial count fit path, reached from bart2's family =
# "nbinom" branch. A SINGLE forest
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
  keepSampler = FALSE,
  samplerOnly = FALSE
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
    splitProbsDefault = formals(dbarts::bart)[["split.probs"]]
  )

  samplerCall <- buildHostSamplerCall(
    matchedCall,
    control,
    priors,
    family = "nbinom"
  )

  sampler <- eval(samplerCall, envir = callingEnv)
  if (isTRUE(samplerOnly)) {
    return(sampler)
  }

  n.chains <- control@n.chains
  n.obs <- length(sampler$data@y)
  n.test <- NROW(sampler$data@x.test)
  n.samples <- control@n.samples

  bc <- bartcoreSampler(sampler, family = "nbinom")
  # bc's engine is the one run below; adopt it into sampler so $fit becomes
  # an R5 wrapper around the engine that actually ran, not the abandoned
  # first-created host. It matters twice here: $getDispersion() would
  # otherwise answer with the abandoned host's own r rather than the fit's.
  sampler$adoptPointer(bc$ptr)

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
  # keepTrees retains the saved trees predict.bartNegbin replays through (the
  # sweeps wrote them regardless), and dispersion.raw supplies predict's
  # per-draw r in the raw n.samples x n.chains layout that pairs with the
  # replayed latent draws.
  if (control@keepTrees) {
    result$dispersion.raw <- dispersionRaw
  }
  if (control@keepTrees || keepSampler) {
    result$fit <- sampler
  }
  result
}

# Assemble a bart2(family = "nbinom") fit from the synthesized channels.
# yhat.train/test are the mean counts mu = r exp(f + o), the reported
# deliverable; latent.train/test are the
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

# Splits a hurdle response into its two ingested parts: the occupancy
# indicator z = 1{y > 0} over all n, and the subset mask {i : y_i > 0} with
# the positive part's working response log(y[S]). y must be finite and
# non-negative (the nbinom non-negative-count precedent, extended to a
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

# A column complete on the positive rows but NA anywhere else is NA only on
# the zero rows (the two row sets are exhaustive): the positive part trains
# on the positive rows alone, so it never sees that column's missing values,
# yet its OWN 'test' call - forced to the full design so the combine covers
# the zero rows too - evaluates it at every row. A column NA on both row
# sets stays constructible: the positive part saw the missing value in
# training and has a route for it.
refuseHurdlePositiveMissingness <- function(x, positive) {
  numColumns <- NCOL(x)
  offending <- integer(0L)
  for (j in seq_len(numColumns)) {
    column <- x[, j]
    if (anyNA(column[positive])) {
      next
    }
    if (!anyNA(column[!positive])) {
      next
    }
    offending <- c(offending, j)
  }
  if (length(offending) == 0L) {
    return(invisible(NULL))
  }
  predictorNames <- colnames(x)
  labels <- if (is.null(predictorNames)) {
    paste0("column ", offending)
  } else {
    paste0("'", predictorNames[offending], "'")
  }
  shown <- labels[seq_len(min(5L, length(labels)))]
  stop(
    "family = \"hurdle.lognormal\": ",
    toString(shown),
    if (length(labels) > 5L) {
      paste0(" and ", length(labels) - 5L, " more column(s)")
    },
    " carry missing values only on the zero (y == 0) rows: the positive ",
    "part trains on the positive rows alone, so it learns no route for ",
    "them, and the hurdle composition evaluates the positive part on ",
    "every row"
  )
}

# The hurdle.lognormal fit path, reached from the front door's family =
# "hurdle.lognormal" branch. Composed R-side from two ordinary single-forest
# fits - never a coupled engine model - so this simply calls the front door
# twice, once per component, at two independently derived seeds (a shared seed
# would correlate the two chains and bias the combined credible interval).
# The positive fit's 'test' is forced to the FULL training x so its
# in-sample fitted()/extract() carries E[y | y > 0, x] at the zero rows it
# never trained on. family is always explicit and the matrix interface
# only: the response is split before any dbartsData ingestion runs (the
# discrete-time hazard precedent, R/dbarts.R's extractSurvivalTimes), which
# a formula LHS cannot supply.
bart2Hurdle <- function(matchedCall, callingEnv, control, formula, data, seed) {
  if (
    is.formula(formula) ||
      inherits(formula, "dbartsData") ||
      inherits(formula, "dgCMatrix") ||
      missing(data)
  ) {
    stop(
      "family = \"hurdle.lognormal\" fits currently use the matrix interface ",
      "- bart(x, y, family = \"hurdle.lognormal\")"
    )
  }

  split <- splitHurdleResponse(data)
  refuseHurdlePositiveMissingness(formula, split$positive)
  xPositive <- formula[split$positive, , drop = FALSE]

  # independent per-component seeds derived deterministically from the
  # user's seed (the per-chain seed derivation); politely
  # restore the caller's RNG stream afterward
  seeds <- if (!is.na(seed)) {
    withFixedSeed(seed, sample.int(.Machine$integer.max, 2L))
  } else {
    c(NA_integer_, NA_integer_)
  }

  # redirectCall forwards every bart2 formal the caller supplied, including
  # names this outer site already diagnosed against "hurdle.lognormal";
  # strip them before either component call runs, since each component
  # would otherwise re-diagnose them against its own forced family
  # (sigest/sigdf/sigquant/resid.prior a false "probit" diagnostic on the
  # occupancy call, since they are genuinely live on the positive half).
  # tree.prior/node.prior are NOT stripped from either list: they are live on
  # both components, so they flow to both exactly as power/base/k/prior.scale
  # already do. The family-only settings ride the family object, which each
  # component call replaces with its own, so nothing is left to strip.
  gatedOnOccupancyOnly <- c("sigest", "sigdf", "sigquant", "resid.prior")

  occupancyCall <- redirectCall(matchedCall, dbarts::bart)
  occupancyCall[gatedOnOccupancyOnly] <- NULL
  occupancyCall$formula <- formula
  occupancyCall$data <- split$z
  occupancyCall$family <- "probit"
  occupancyCall$seed <- seeds[1L]
  occupancyCall$keepTrees <- control@keepTrees
  occupancy <- eval(occupancyCall, callingEnv)

  positiveCall <- redirectCall(matchedCall, dbarts::bart)
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
    # both components come from the same matchedCall, so they share n.chains
    n.chains = occupancy$n.chains,
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

# S(t | x) draws from an AFT linear predictor and its residual scale, in the
# uncombined convention (chains x samples x observations) where the scale
# draws align with the fit draws unambiguously - the loglik channel's
# approach. The caller supplies the linear predictor and combines the chain
# margin at the end.
#
# 'sigma' is one scalar per draw for a homoscedastic fit and one value per draw
# AND observation for a heteroscedastic one, laid out as the linear predictor
# is; the first recycles across the observation margin, the second pairs with
# it element for element, and both are the same as.vector() enumeration.
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
  numDraws <- prod(drawDims)
  if (!length(sigma) %in% c(numDraws, numDraws * numObservations)) {
    stop("the fit's draw count does not match its residual scale draws")
  }

  result <- array(0, c(drawDims, numTimes, numObservations))
  # in column-major order the draw margins vary fastest, so a scalar-per-draw
  # scale recycles across observations exactly as pointwiseLogLikelihood's does
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

# Survival-probability draws from a discrete-time hazard fit. The fit is an
# ordinary binary fit on the person-period-expanded rows, so the
# per-(subject, period) hazards are
# h(k | x) = g(f(x, k) + o) through the fit's link (probit/logistic), and
# S(t | x) = prod_{k : periods[k] <= t} (1 - h(k | x)). A subject that rode
# 'test' at fit time (dbarts()'s own hazard test acceptance, person-period-
# expanded on the SAME training grid) has its hazards already stored, read
# straight off the packaged fit - no trees, no re-expansion. Every OTHER
# 'newdata' still ALWAYS re-expands and predicts - the training design is
# ragged (subject i has only its at-risk rows), so a stored per-subject fit
# cannot supply hazards past its own event/censoring period - and therefore
# requires keepTrees, where aft's training path never does. Requested `times`
# are horizons (default: the grid periods); S cumulates (1 - h) through each.
# Returns draws per the package's three-tier convention, the shape aft's
# method uses (draws x times x observations, a chain margin under
# combineChains = FALSE).
hazardSurvivalProbabilities <- function(object, times, newdata, combineChains) {
  periods <- object$periods
  K <- length(periods)
  if (is.null(times)) {
    times <- periods
  }
  times <- as.double(times)
  if (length(times) == 0L || any(!is.finite(times)) || any(times <= 0)) {
    stop("'times' must be finite and positive")
  }

  usesStoredTest <- is.null(newdata) && !is.null(object[["yhat.test"]])
  if (usesStoredTest) {
    # extract() on the packaged fit, not object$fit (the dbartsSampler):
    # extract.dbartsSampler accepts only type = "predictors", while
    # extract.bart's sample = "test" arm reads object$yhat.test directly and
    # applies the same probability transform predict(type = "ev") does
    haz <- extract(object, type = "ev", sample = "test", combineChains = FALSE)
    n <- dim(haz)[length(dim(haz))] %/% K
  } else {
    if (is.null(object[["fit"]])) {
      stop(
        "survivalProbabilities on a discrete-time hazard fit requires the ",
        "trees; refit with keepTrees = TRUE"
      )
    }
    periodCol <- ncol(object$fit$data@x)
    if (is.null(newdata)) {
      # reconstruct the per-subject covariates from the coded expanded
      # design: every subject is at risk in period 1, so the period-1 rows
      # are the subjects in order, their covariate columns the (coded)
      # per-subject x
      fitX <- extract(object$fit, "predictors")
      subjectCov <- fitX[fitX[, periodCol] == 1L, -periodCol, drop = FALSE]
      n <- nrow(subjectCov)
      # name the appended column "period" under the same rule the training
      # design used, so a named fit's re-expanded design matches by name
      bigX <- appendHazardPeriodColumn(
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
      bigX <- appendHazardPeriodColumn(
        newdata[rep(seq_len(n), times = K), , drop = FALSE],
        rep(seq_len(K), each = n)
      )
    }

    # hazards through the correct link (type = "ev" keys on $family, the
    # binary token); predict codes bigX to the training columns and replays
    # the trees
    haz <- predict(object, bigX, type = "ev", combineChains = FALSE)
  }
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

# Survival-probability draws from an AFT fit. Under the log-normal model
# log T = f(x) + s(x) eps, S(t | x) =
# 1 - Phi((log t - f(x)) / s(x)), evaluated at every posterior draw of f
# and of the residual scale - the scalar sigma homoscedastically, the variance
# surface s(x) under a variance forest, where the scale is per draw AND per
# observation. Returns draws per the package's three-tier convention
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
  refuseUnusedGenericArgs(
    list(...),
    "survivalProbabilities",
    "bart",
    foreignArgsFor(
      survivalProbabilitiesForeignReasons,
      names(formals(survivalProbabilities.bart))
    )
  )
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

  # The scale the normal tail divides by. A homoscedastic fit's is the stored
  # per-draw sigma at any rows; a heteroscedastic fit's is the surface, which
  # the training rows read off s.train and 'newdata' off the replay predict
  # parks on its result. A sampler that replays no variance surface has no
  # scale at those rows, so the curves are refused rather than drawn at the
  # pinned sigma - the wording the ppd branch uses for the same gap.
  scale <- if (is.null(object[["s.train"]])) {
    object[["sigma"]]
  } else if (is.null(newdata)) {
    heteroscedasticScale(object[["s.train"]], n.chains)
  } else {
    replayed <- attr(linearPredictor, "s")
    if (is.null(replayed)) {
      stop(
        "survival probabilities at 'newdata' are not available on a ",
        "heteroscedastic fit whose sampler replays no variance surface; ",
        "refit with keepTrees = TRUE"
      )
    }
    heteroscedasticScale(replayed, n.chains)
  }

  survivalProbabilitiesFromDraws(
    linearPredictor,
    scale,
    times,
    n.chains,
    combineChains
  )
}

# The BayesTree-style door is a strict compatibility mode: 0.9-34's argument
# list, 0.9-34's defaults, no families beyond the numeric/binary pair the
# response itself declares. A factor response of three or more levels was
# fit as its integer level codes there, which is a model no one asks for on
# purpose, so it is refused here naming both remedies - the modern door's
# multinomial fit, or the explicit coding that reproduces the old numbers.
refuseLegacyFactorResponse <- function() {
  stop(
    "'bartBT' does not fit a factor response with three or more levels; ",
    "dbarts 0.9-x fit its integer level codes as numbers. Use ",
    "bart(x, y, family = \"multinomial\") for a multinomial ",
    "fit, or pass as.integer(y) - 1L to keep the old numeric behaviour",
    call. = FALSE
  )
}

bartBT <- function(
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
  keepsampler = keeptrees
) {
  # a count-matrix data object declares the multinomial model, which this
  # frozen BayesTree shim reports none of: its packaging is one location per
  # observation, so a K-location fit would be reshaped and returned without a
  # word. The dbartsData passthrough is the only route one can arrive by
  refuseCountsCarryingData(x.train, "bartBT()")
  refuseResponseFreeFormula(x.train, "bartBT()")

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
  seed <- coerceOrError(seed, "integer")

  # named ahead of dbartsControl(), whose own validity would otherwise
  # blame n.thin/n.burn - its slot names, not the formals these came in as
  if (isTRUE(keepevery <= 0L)) {
    stop("'keepevery' must be a positive integer")
  }
  if (isTRUE(nskip < 0L)) {
    stop("'nskip' must be a non-negative integer")
  }

  # the matrix interface names the response directly; a formula response
  # reaches the same refusal through the backstop below, after dbarts() has
  # resolved it
  if (
    !is.formula(x.train) &&
      !missing(y.train) &&
      is.factor(y.train) &&
      nlevels(y.train) >= 3L
  ) {
    refuseLegacyFactorResponse()
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
    seed = seed
  )
  matchedCall <- if (keepcall) match.call() else call("NULL")
  control@call <- matchedCall
  control@n.burn <- control@n.burn %/% control@n.thin
  # printEvery counts post-thinning samples and must stay positive: thinning
  # heavier than it would otherwise drive it to 0, which the sampler refuses
  control@printEvery <- max(1L, control@printEvery %/% control@n.thin)
  # control@keepTrees is still FALSE here regardless of 'keeptrees' (set
  # below, once burn-in is known); read the user's own argument, as the
  # modern door does through its differently-sequenced control construction
  keepsampler <- keepsampler || keeptrees
  if (control@n.burn == 0L && keeptrees == TRUE) {
    control@keepTrees <- TRUE
  }
  if (control@n.burn > 0L) {
    control@keepTrees <- FALSE
  }
  ndpost <- ndpost %/% control@n.thin
  # a zero (or thinned-to-zero) sample count would otherwise fault deeper, in
  # the empty-array reshape (dim(X) has no positive length); mirrors the
  # modern door's same-shaped guard on control@n.samples
  if (isTRUE(ndpost <= 0L)) {
    stop("'ndpost' must be a positive integer")
  }

  priors <- buildSamplerPriors(
    matchedCall,
    power,
    base,
    sigdf,
    sigquant,
    nodeK = if (!is.null(matchedCall[["k"]])) matchedCall[["k"]] else k,
    priorScale = NA_real_,
    splitProbsName = "splitprobs",
    splitProbsDefault = formals(dbarts::bartBT)[["splitprobs"]]
  )
  tree.prior <- priors$tree.prior
  node.prior <- priors$node.prior
  resid.prior <- priors$resid.prior

  # the frozen BayesTree-compatibility shim keeps dummy expansion. Every
  # setting 0.9-34 had no name for is left at dbarts()'s own default rather
  # than forwarded: this door adds no capability, so a numeric or binary
  # response is all it can be given.
  args <- list(
    formula = x.train,
    data = y.train,
    test = x.test,
    weights = weights,
    offset = binaryOffset,
    verbose = as.logical(verbose),
    n.samples = as.integer(ndpost),
    tree.prior = tree.prior,
    node.prior = node.prior,
    resid.prior = resid.prior,
    proposal.probs = proposalprobs,
    control = control,
    sigest = as.numeric(sigest),
    factors = "indicators",
    na.action = stats::na.omit
  )
  # 0.9-34's row rule, unchanged: a row with a missing value anywhere is
  # dropped rather than modelled. The modern door's na.action default keeps
  # missing predictors; this door adds no capability.
  sampler <- tryCatch(
    do.call(dbarts::dbarts, args, envir = parent.frame(1L)),
    error = function(e) {
      msg <- conditionMessage(e)
      # the formula path's categorical response is refused inside dbarts(),
      # whose message names the modern door's family tokens; this door has no
      # 'family' formal, so its own two-remedy message stands in
      if (grepl("response is multinomial; fit it with", msg, fixed = TRUE)) {
        refuseLegacyFactorResponse()
      }
      stop(e)
    }
  )

  # formula-response backstop for the pre-check above: dbarts() resolved an
  # ordered-factor response to ordinal, which this door does not package
  if (identical(sampler$model@family, "ordinal")) {
    refuseLegacyFactorResponse()
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
  # needed to extract ppd; mirrors the modern door's packageBartResults
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

makeind <- function(x, all = TRUE) {
  ignored <- all ## for R check # nolint: object_usage_linter.
  makeModelMatrixFromDataFrame(x, TRUE)
}
