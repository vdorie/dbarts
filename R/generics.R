# predict, extract, fitted, and residuals methods for bart and the
# multinomial, ordinal, nbinom, and hurdle fit objects

extract <- function(object, ...) UseMethod("extract")

plotTree <- function(object, ...) UseMethod("plotTree")

survivalProbabilities <- function(object, ...) {
  UseMethod("survivalProbabilities")
}

# What a fit is, read off its descriptors rather than off which draw
# channels a run happened to keep. $family is the family as specified; the
# engine family its link and likelihood follow is looked up from it below, so
# a family this table does not name stops rather than taking some default
# link. A fit saved by dbarts 0.9-x carries no family element, and in that
# release only a gaussian fit drew sigma, so the fallback in fitFamily is the
# one place a channel's presence still decides anything.
familyEngineTokens <- c(
  gaussian = "gaussian",
  student = "gaussian",
  probit = "probit",
  logistic = "logistic",
  aft = "aft",
  hazard.probit = "probit",
  hazard.logistic = "logistic",
  multinomial = "multinomial",
  ordinal = "ordinal",
  nbinom = "nbinom",
  hurdle.lognormal = "hurdle.lognormal"
)

fitFamily <- function(object) {
  family <- object[["family"]]
  if (!is.null(family)) {
    return(family)
  }
  if (is.null(object[["sigma"]])) "probit" else "gaussian"
}

fitEngineFamily <- function(object) {
  family <- fitFamily(object)
  engine <- familyEngineTokens[family]
  if (is.na(engine)) {
    stop("fit has an unknown family '", family, "'")
  }
  unname(engine)
}

fitIsBinary <- function(object) {
  fitEngineFamily(object) %in% c("probit", "logistic")
}

# a family with a residual law, and so the resid.scale descriptor
fitHasResidual <- function(object) {
  fitEngineFamily(object) %in% c("gaussian", "aft")
}

# the residual law's shape, which the family fixes
fitIsStudent <- function(object) {
  identical(fitFamily(object), "student")
}

fitIsHazard <- function(object) {
  fitFamily(object) %in% c("hazard.probit", "hazard.logistic")
}

fitIsHeteroscedastic <- function(object) {
  identical(object[["resid.scale"]], "forest")
}

fitNumForests <- function(object) {
  if (is.null(object[["n.forests"]])) 1L else object[["n.forests"]]
}

# extract's model- and predictor-level types (sigma, k, shape,
# thresholds, varcount) refuse a caller-supplied 'sample' by name
refuseSampleOnModelType <- function(type, sampleSupplied) {
  if (sampleSupplied) {
    stop(
      "'sample' is not used when type = \"",
      type,
      "\": it is not per-observation"
    )
  }
}

# the family as specified: the resolved family object, with its settings
family.bart <- function(object, ...) {
  familySpec <- object[["family.spec"]]
  if (is.null(familySpec)) {
    familySpec <- newValidated("dbartsFamily", token = fitFamily(object))
  }
  familySpec
}
family.bartMultinomial <- family.bart
family.bartOrdinal <- family.bart
family.bartNegbin <- family.bart
family.bartHurdle <- family.bart

# latent-scale draws to probabilities for a binary fit, by its link
probabilityFromLatents <- function(latents, object) {
  switch(
    fitEngineFamily(object),
    probit = pnorm(latents),
    logistic = plogis(latents),
    stop("fit's family '", fitFamily(object), "' has no binary link")
  )
}

# A heteroscedastic fit's residual scale is the per-observation surface s(x),
# stored - like the draw channels it is laid
# out as - either combined or split. This normalizes whichever storage the fit
# used to the split, chain-fastest layout, so as.vector() on it enumerates
# draws in the order as.vector() on the split fits does and the two pair
# element for element. s(x) is already on the response scale (the working
# surface times the response range), so it is the fit's whole residual scale:
# a heteroscedastic fit carries no scalar sigma to scale. NULL passes through,
# marking a homoscedastic fit, whose scale is its scalar sigma.
heteroscedasticScale <- function(s, n.chains) {
  if (is.null(s)) NULL else combineOrUncombineChains(s, n.chains, FALSE)
}

# per-draw, per-observation log-likelihood of the stored training response.
# ev enters with chains split ((n.chains x) n.samples x n.obs), so that
# as.vector(ev) enumerates draws chain-fastest. A scalar-per-draw field
# (sigma, resid.df) may be STORED combined (a flat, chain-major vector -
# chain 1's whole run, then chain 2's, ...) or split ((n.chains x)
# n.samples matrix, chain-fastest); chainFastest below normalizes either
# storage to the split matrix, so as.vector() on it always yields the
# chain-fastest order ev's own as.vector() does, and the two pair by plain
# recycling regardless of how the fit itself was combined. Dispatch is on
# object$family, not on the presence of sigma, so a new family cannot
# silently reuse a formula that does not fit it (an aft fit has non-null
# sigma but is not gaussian): gaussian evaluates the normal density with
# weights as precision (y | x ~ N(f(x), sigma^2 / w)), and a
# student() residual law the t marginal at the same location and scale; probit and
# logistic the bernoulli mass on the y scale, weights being trial counts for
# logistic (probit never stores weights); aft the log density for events and
# the log survival tail for right-censored rows, mirroring the engine's
# AFTResponse::computeLogLikelihood. Any other family errors rather than
# reporting a wrong number. A row an installed active-row mask takes out of
# the data set reports NaN whatever the family, as the engine's own channel
# does. A heteroscedastic gaussian or aft fit scores at its
# own per-observation s(x) instead of the scalar (heteroscedasticScale below).
pointwiseLogLikelihood <- function(object, ev) {
  y <- object[["y"]]
  if (is.null(y)) {
    stop(
      "cannot compute the log-likelihood; fit does not store the training response"
    )
  }
  family <- fitEngineFamily(object)
  weights <- object[["weights"]]
  n.draws <- length(ev) %/% length(y)
  y <- rep(y, each = n.draws)
  n.chains <- fitNChains(object)
  chainFastest <- function(x) {
    if (is.null(dim(x))) uncombineChains(as.vector(x), n.chains) else x
  }

  if (identical(family, "gaussian")) {
    # the family fixes the residual law. A student() fit scores the
    # MARGINAL t density - the observation-level likelihood loo/waic are
    # defined on, and the density the engine itself reports - rather than the
    # gaussian working likelihood conditional on the latent precisions, which
    # is a different quantity.
    isStudent <- fitIsStudent(object)
    # s(x) is one value per draw AND observation, so it pairs with ev directly
    # rather than recycling across the observation margin as sigma does; a
    # length mismatch means the two channels were not written by the same run
    s <- heteroscedasticScale(object[["s.train"]], n.chains)
    sd <- if (is.null(s)) {
      rep_len(as.vector(chainFastest(object$sigma)), length(ev))
    } else if (length(s) != length(ev)) {
      stop("the fit's 's.train' draws do not match its fitted draws")
    } else {
      as.vector(s)
    }
    if (!is.null(weights)) {
      sd <- sd / rep(sqrt(weights), each = n.draws)
    }
    if (isStudent) {
      # sigma is the CONDITIONAL scale under the scale mixture, so the marginal
      # is a location-scale t_nu with that scale: sqrt(w) (y - f(x)) / sigma ~
      # t_nu. The df is one scalar per draw, as sigma is, so it pairs by the
      # same recycling.
      df <- object[["resid.df"]]
      if (is.null(df)) {
        stop(
          "cannot compute the log-likelihood; fit does not store the per-draw residual degrees of freedom"
        )
      }
      df <- rep_len(as.vector(chainFastest(df)), length(ev))
      result <- dt((y - as.vector(ev)) / sd, df, log = TRUE) - log(sd)
    } else {
      result <- dnorm(y, as.vector(ev), sd, log = TRUE)
    }
    # a zero-weight row is not in the model, so the channel flags it as
    # unavailable rather than reporting the -Inf an infinite sd would give
    if (!is.null(weights)) {
      result[rep(weights, each = n.draws) == 0] <- NaN
    }
  } else if (identical(family, "probit") || identical(family, "logistic")) {
    result <- dbinom(y, 1L, as.vector(ev), log = TRUE)
    if (!is.null(weights)) {
      result <- rep(weights, each = n.draws) * result
    }
  } else if (identical(family, "aft")) {
    status <- object[["status"]]
    if (is.null(status)) {
      stop(
        "cannot compute the aft log-likelihood; fit does not store the censoring status"
      )
    }
    # the residual scale and y are on the log-time scale (y is log event time
    # for an event, log censoring time for a censored row); events keep the
    # normal density, censored rows take the log upper survival tail
    # log P(log T > log C). A heteroscedastic aft fit's scale is per draw AND
    # per observation, so it pairs with ev directly where the scalar sigma
    # recycles - the gaussian branch's split, under the same length check
    s <- heteroscedasticScale(object[["s.train"]], n.chains)
    sd <- if (is.null(s)) {
      rep_len(as.vector(chainFastest(object$sigma)), length(ev))
    } else if (length(s) != length(ev)) {
      stop("the fit's 's.train' draws do not match its fitted draws")
    } else {
      as.vector(s)
    }
    location <- as.vector(ev)
    result <- dnorm(y, location, sd, log = TRUE)
    censored <- rep(status, each = n.draws) == 0
    result[censored] <- pnorm(
      y[censored],
      location[censored],
      sd[censored],
      lower.tail = FALSE,
      log.p = TRUE
    )
  } else {
    stop(
      "family '",
      if (is.null(family)) "NULL" else family,
      "' does not support the log-likelihood"
    )
  }
  # a row the active-row mask takes out of the data set - what a probit or
  # ordinal fit's 0/1 case weights install - is not in the model and has no
  # likelihood to report, so the channel gives NaN there rather than the
  # finite value the row's fit would still yield. That is the engine's own
  # convention on this channel, and the gaussian branch's zero-weight rule
  # above is the same statement through the other channel.
  active <- object[["active"]]
  if (!is.null(active)) {
    result[rep(active, each = n.draws) == 0] <- NaN
  }
  array(result, dim(ev), dimnames(ev))
}

# per-observation posterior summary for the interval-returning generics: est
# (the posterior mean) plus a symmetric ci.level credible band from the draw
# quantiles, pooled over every margin except the trailing 'trailing' ones
# (observations are the sole trailing margin for the bart-family generics, as
# in the mean path). The interval KIND follows the caller's type: "ev" gives a
# credible interval for E[Y|x] (a probability for binary), "ppd" a prediction
# interval that also carries the residual noise, and "bart" a credible
# interval on the latent scale. A K-widened channel (a category-probability
# draw) keeps K as a second trailing margin (trailing = 2) instead of pooling
# across categories, which would average incomparable probabilities; the
# result is then an array with est/ci.lower/ci.upper on a new trailing margin
# rather than a 3-column matrix, since a plain matrix cannot carry both an
# observation and a category index.
posteriorInterval <- function(draws, ci.level, trailing = 1L) {
  if (
    !is.numeric(ci.level) ||
      length(ci.level) != 1L ||
      is.na(ci.level) ||
      ci.level <= 0 ||
      ci.level >= 1
  ) {
    stop("'ci.level' must be a single number in (0, 1)")
  }
  probs <- c((1 - ci.level) / 2, 1 - (1 - ci.level) / 2)
  if (is.null(dim(draws))) {
    result <- matrix(
      c(mean(draws), quantile(draws, probs, names = FALSE)),
      nrow = 1L
    )
    colnames(result) <- c("est", "ci.lower", "ci.upper")
    return(result)
  }
  d <- dim(draws)
  keepAxes <- seq.int(length(d) - trailing + 1L, length(d))
  est <- channelMeans(draws, trailing)
  bounds <- apply(draws, keepAxes, quantile, probs = probs, names = FALSE)
  if (trailing == 1L) {
    result <- cbind(est, bounds[1L, ], bounds[2L, ])
    colnames(result) <- c("est", "ci.lower", "ci.upper")
    return(result)
  }
  result <- array(
    c(as.vector(est), as.vector(bounds[1L, , ]), as.vector(bounds[2L, , ])),
    dim = c(dim(est), 3L)
  )
  dn <- dimnames(est)
  dimnames(result) <- c(
    if (is.null(dn)) rep(list(NULL), length(dim(est))) else dn,
    list(c("est", "ci.lower", "ci.upper"))
  )
  result
}

## Adds a length-1 leading margin to 'x' (dec-A79: extract and predict keep a
## chain dimension of length 1 on a one-chain fit under combineChains =
## FALSE, where a one-chain fit's own storage carries none to begin with).
## A vector becomes a 1 x length matrix; an array of any other rank gains a
## leading dimension. Either way 'x's own dimnames ride the old margins, a
## vector's own names among them, since assigning 'dim' would otherwise drop
## them.
addChainDimension <- function(x) {
  d <- dim(x)
  if (is.null(d)) {
    nms <- names(x)
    x <- array(x, c(1L, length(x)))
    if (!is.null(nms)) {
      dimnames(x) <- list(NULL, nms)
    }
    return(x)
  }
  dn <- dimnames(x)
  x <- array(x, c(1L, d))
  if (!is.null(dn)) {
    dimnames(x) <- c(list(NULL), dn)
  }
  x
}

## The inverse of addChainDimension: drops the leading length-1 margin 'x' is
## known to carry.
dropChainDimension <- function(x) {
  d <- dim(x)[-1L]
  dn <- dimnames(x)
  if (length(d) <= 1L) {
    nms <- if (!is.null(dn)) dn[[2L]] else NULL
    x <- as.vector(x)
    if (!is.null(nms)) {
      names(x) <- nms
    }
    return(x)
  }
  dim(x) <- d
  if (!is.null(dn)) {
    dimnames(x) <- dn[-1L]
  }
  x
}

combineOrUncombineChains <- function(x, n.chains, combine) {
  if (length(dim(x)) > 2L && combine) {
    x <- combineChains(x)
  } else if (length(dim(x)) == 2L && !combine) {
    x <- if (n.chains > 1L) {
      uncombineChains(x, n.chains)
    } else {
      addChainDimension(x)
    }
  }
  x
}

# combineOrUncombineChains for a scalar-per-draw field (sigma, k, shape):
# combined is a chain-major vector, uncombined a chains x samples matrix (a
# 1 x samples one at one chain, dec-A79)
reshapeScalarChannel <- function(x, n.chains, combine) {
  if (is.null(dim(x))) {
    if (combine) {
      x
    } else if (n.chains > 1L) {
      uncombineChains(x, n.chains)
    } else {
      addChainDimension(x)
    }
  } else {
    if (combine) as.vector(t(x)) else x
  }
}

## convertSamplesFromDbartsToBart, keeping a one-chain fit's chain margin
## under combineChains = FALSE (dec-A79), for predict's own per-call reshape
## of fresh engine output. convertSamplesFromDbartsToBart itself is
## untouched: R/bart.R's packaging calls pass the fit's own combineChains,
## and a fit's stored fields never gain this margin on their own.
convertSamplesForCaller <- function(samples, n.chains, combineChains) {
  x <- convertSamplesFromDbartsToBart(samples, n.chains, combineChains)
  if (combineChains || n.chains > 1L) x else addChainDimension(x)
}

# The per-call worker count for a saved-tree replay. The engine partitions by
# (chain, draw) and reduces nothing across workers, so this moves no value; it
# is refused rather than floored because a zero, a negative, an NA or a
# non-numeric is a caller mistake, and the offending value is echoed so which
# call carried it is visible. Every predict method takes it as its LAST
# positional formal, so the argument a caller is most likely to supply by
# position - 'type' - stays third on every one of them.
validatePredictThreads <- function(n.threads) {
  if (
    !is.numeric(n.threads) ||
      length(n.threads) != 1L ||
      is.na(n.threads) ||
      n.threads < 1L ||
      n.threads != round(n.threads)
  ) {
    stop(
      "'n.threads' must be a single positive integer, not ",
      deparse(n.threads)[1L]
    )
  }
  as.integer(n.threads)
}

# One offset spelling is live across every predict method - 'offset'. The
# fit-time channels keep 'offset.test' (dbartsData, bart2, and the
# sampler's own $predict), so a caller carrying that name here would otherwise
# vanish into '...' with the offset silently dropped instead of applied.
predictOffsetUnusedArgs <- list(
  offset.test = "this fit's out-of-sample offset argument is named 'offset'"
)

# predict, extract(type = "trees") and plotTree all read the fit's SAVED
# trees, so a fit kept without them has nothing to read. The message names
# the one argument that keeps them rather than restating the condition.
refuseWithoutTrees <- function(what, keepTrees = "keepTrees") {
  stop(
    what,
    " requires the fit's saved trees; refit with ",
    keepTrees,
    " = TRUE"
  )
}

# bartBT spells it 'keeptrees', bart 'keepTrees'. A fit kept with
# keepCall = FALSE stores no call and names neither, so it takes bart's
# spelling, which is the surface such a fit most likely came from.
bartKeepTreesArgument <- function(object) {
  if (callName(object[["call"]]) == "bartBT") "keeptrees" else "keepTrees"
}

# An offset shifts the latent at rows the sampler never saw, and a hurdle fit
# replays its trees with no offset channel at all, so either spelling would be
# dropped rather than applied.
noPredictOffsetReason <- paste0(
  "this fit has no out-of-sample offset channel; predict replays the ",
  "offset-free surface"
)
predictNoOffsetUnusedArgs <- list(
  offset = noPredictOffsetReason,
  offset.test = noPredictOffsetReason
)

# 'offset' occupies the same slot on all six predict methods so the argument
# order is uniform in position and not merely in relative order; the family
# with no offset channel takes it as a formal and refuses a non-NULL value with
# the same wording it would carry out of '...'.
refusePredictOffsetChannel <- function(offset, class) {
  if (!is.null(offset)) {
    refuseUnusedGenericArgs(
      list(offset = offset),
      "predict",
      class,
      predictNoOffsetUnusedArgs
    )
  }
  invisible(NULL)
}

# An offset is a number per row, or a matrix of them; a logical or character
# value would be coerced silently into one (TRUE an offset of 1), so it is
# refused by name. A missing scalar stays for the refusals that name it.
refuseNonNumericOffset <- function(offset) {
  if (
    is.null(offset) ||
      is.numeric(offset) ||
      (is.logical(offset) && all(is.na(offset))) ||
      (is.data.frame(offset) && all(vapply(offset, is.numeric, logical(1L))))
  ) {
    return(invisible(NULL))
  }
  stop("'offset' must be numeric", call. = FALSE)
}

# The offset predict applies at newdata, as predict.lm forms it: the fit's
# 'offset' argument and offset() terms evaluated there, plus the caller's
# 'offset'. An argument that cannot be evaluated there (a plain vector given
# for the training rows) is refused unless the caller gives 'offset' for
# these rows, which then stands in for it.
predictTermOffset <- function(data, newdata, offset, caller = "predict") {
  if (missing(newdata) || is.null(newdata)) {
    return(offset)
  }
  argument <- evaluateOffsetArgument(attr(data, "offset.argument"), newdata)
  if (isFALSE(argument)) {
    if (is.null(offset)) {
      stop(
        "the fit's 'offset' was given as ",
        describeOffsetArgument(attr(data, "offset.argument")),
        ", which cannot be evaluated on the rows of 'newdata'; give ",
        caller,
        " an 'offset' for them"
      )
    }
    argument <- NULL
  }
  offset <- addOffsetShares(argument, offset, "offset", "newdata")
  addFormulaTermOffset(data@x, newdata, offset, "offset", "newdata")
}

predict.bart <- function(
  object,
  newdata,
  type = c("ev", "ppd", "bart", "forest", "sigma"),
  offset = NULL,
  weights = NULL,
  combineChains = TRUE,
  ci.level = NULL,
  forest = NULL,
  bases = NULL,
  na.action = dbarts::na.keepPredictors,
  n.threads = object$fit$control@n.threads,
  ...
) {
  if (is.null(object[["fit"]])) {
    refuseWithoutTrees("predict", bartKeepTreesArgument(object))
  }

  refuseUnusedGenericArgs(
    list(...),
    "predict",
    "bart",
    c(
      predictOffsetUnusedArgs,
      foreignArgsFor(predictForeignReasons, names(formals(predict.bart)))
    )
  )
  warnUnusedDots(list(...), "predict", "bart")
  refuseNonNumericOffset(offset)
  type <- validateType(type, eval(formals(predict.bart)$type))
  # above the type = "forest" and amplitude-blend returns below, so every arm's
  # value is checked rather than only the one that reaches the sampler here
  n.threads <- validatePredictThreads(n.threads)
  refuseForestSelectionOutsideForestArm(
    type,
    forest,
    fitIsHeteroscedastic(object)
  )
  refuseDroppedForestChannel(object)
  if (type == "sigma") {
    if (!fitIsHeteroscedastic(object)) {
      stop(
        "type = \"sigma\" predicts a heteroscedastic fit's per-observation ",
        "scale; this fit's sigma is one scalar per draw, which ",
        "extract(type = \"sigma\") returns"
      )
    }
    if (!is.null(weights)) {
      stop(
        "type = \"sigma\" does not support 'weights': it reports the ",
        "variance forest's scale s(x), which a case weight does not change"
      )
    }
  }

  # both amplitude arms read the SAVED trees draw by draw, pairing each draw's
  # forests with that draw's own amplitudes; without the tree store only the
  # current trees replay, one set standing for every draw, and the pairing the
  # arms are defined by does not exist. A plain single-forest predict keeps its
  # long-standing keepTrees-free reading of the current trees.
  if (
    (type == "forest" || !is.null(object[["forestFits"]])) &&
      !object$fit$control@keepTrees
  ) {
    stop(
      "predict requires the fit's saved trees; refit with ",
      bartKeepTreesArgument(object),
      " = TRUE: an amplitude-coupled fit pairs each saved draw's forests ",
      "with that draw's own amplitudes, and without the tree store only the ",
      "current trees replay, one set for every draw"
    )
  }

  # without the tree store only the current trees replay: one chain's are the
  # long-standing keepTrees-free reading, but several chains' current trees
  # are one evaluation each, not a sequence of draws to report
  if (!object$fit$control@keepTrees && object$fit$control@n.chains > 1L) {
    stop(
      "predict requires the fit's saved trees; refit with ",
      bartKeepTreesArgument(object),
      " = TRUE: without the tree store only each chain's current trees ",
      "replay, one evaluation per chain rather than a draw per sample"
    )
  }

  # the per-forest arm answers off the sampler's own replay and shares none of
  # the combined arms' machinery below: there is no ci.level band, no latent
  # transform and no s(x) attribute on a raw per-forest total
  if (type == "forest") {
    if (!is.null(bases)) {
      stop(
        "'bases' does not apply to type = \"forest\": that arm reports each ",
        "forest's own total BEFORE any basis, which is what leaves the ",
        "recombination to the caller"
      )
    }
    if (!is.null(ci.level)) {
      stop(
        "type = \"forest\" does not support 'ci.level': that arm reports ",
        "each forest's own total before any basis"
      )
    }
  }
  if (type != "forest" && is.null(object[["forestFits"]]) && !is.null(bases)) {
    numForests <- fitNumForests(object)
    stop(
      "'bases' is only meaningful on an amplitude-coupled multi-forest fit; ",
      "this fit has ",
      numForests,
      if (numForests == 1L) " forest" else " forests"
    )
  }

  # the fit's offset argument and offset() terms are evaluated on newdata, as
  # predict.lm does, and added to an 'offset' given here; the per-forest arm
  # reports each forest's own total, with no offset folded in
  if (type != "forest") {
    offset <- predictTermOffset(object$fit$data, newdata, offset)
  }
  # validated once, here; the rows na.action keeps are what every arm below
  # predicts, and padPredictedRows puts them back on newdata's rows. A
  # missing offset or weight marks its row incomplete the same way (dec-A89).
  rows <- preparePredictRows(
    newdata,
    object$fit$data@x,
    na.action,
    list(offset = offset, weights = weights)
  )
  if (isTRUE(rows$placeholder)) {
    restoreSeed <- protectRandomSeed()
    on.exit(restoreSeed(), add = TRUE)
  }
  offset <- subsetPredictInput(offset, rows, "offset")
  weights <- subsetPredictInput(weights, rows, "weights")
  # a binary draw's weights at new rows are its trial counts, Binomial(w, p),
  # under either link, as a logistic fit's own weights are
  if (type == "ppd" && !is.null(weights) && fitIsBinary(object)) {
    refuseNonCountWeights(
      weights,
      what = paste0(
        "the posterior predictive 'weights' of a ",
        fitEngineFamily(object),
        " fit are trial counts"
      )
    )
  }

  if (type == "forest") {
    return(padPredictedRows(
      predictForest(
        object,
        rows$x,
        offset,
        combineChains,
        forest,
        n.threads,
        rows$keptNames
      ),
      rows,
      trailing = 1L
    ))
  }

  # an amplitude-coupled fit has no combined test surface in the engine - the
  # sampler holds no basis at the caller's rows - so the combination is done
  # here, from the per-forest replay and the fit's own glue
  if (!is.null(object[["forestFits"]])) {
    bases <- subsetPredictBases(
      bases,
      rows,
      lapply(object$bases, function(basis) {
        if (is.null(basis)) NULL else basis[1L, , drop = FALSE]
      })
    )
    return(padPredictedRows(
      predictBlend(
        object,
        rows$x,
        offset,
        weights,
        type,
        combineChains,
        ci.level,
        bases,
        n.threads,
        rows$keptNames,
        rows$newdata
      ),
      rows,
      first = !is.null(ci.level)
    ))
  }

  n.chains <- object$fit$control@n.chains
  rowNames <- rows$keptNames
  result <- predictCodedTest(object$fit, rows$x, offset, n.threads)
  # a heteroscedastic fit returns list(mean, variance); s(x) rides back as an
  # attribute on the returned yhat so plain predict callers are unaffected
  s <- NULL
  if (is.list(result)) {
    s <- sqrt(convertSamplesForCaller(
      result$variance,
      n.chains,
      combineChains
    ))
    s <- nameObservationMargin(s, rowNames)
    result <- result$mean
  }
  if (type == "sigma") {
    if (is.null(s)) {
      stop(
        "type = \"sigma\" is not available on a heteroscedastic fit whose ",
        "sampler replays no variance surface"
      )
    }
    if (!is.null(ci.level)) {
      return(padPredictedRows(
        posteriorInterval(s, ci.level),
        rows,
        first = TRUE
      ))
    }
    return(padPredictedRows(s, rows))
  }
  # result is n.obs x n.samples x n.chains
  result <- convertSamplesForCaller(result, n.chains, combineChains)
  result <- nameObservationMargin(result, rowNames)

  if (type != "bart") {
    if (fitIsBinary(object)) {
      result <- probabilityFromLatents(result, object)
    }

    if (type == "ppd") {
      # keepFits = FALSE drops a heteroscedastic fit's s.train, which its
      # resid.scale descriptor survives, so this refuses rather than silently
      # sampling as if homoscedastic
      if (
        is.null(s) &&
          is.null(object[["s.train"]]) &&
          fitIsHeteroscedastic(object)
      ) {
        stop(
          "posterior predictive sampling needs this heteroscedastic fit's ",
          "'s.train' draws to tell it apart from a homoscedastic one here, ",
          "and 'keepFits = FALSE' dropped them (automatically, when a ",
          "'callback' was supplied, unless overridden); refit with ",
          "'keepFits = TRUE'"
        )
      }
      # the replayed s(x) above IS the noise scale at these rows; a
      # heteroscedastic fit whose sampler replays none cannot be drawn from
      if (is.null(s) && !is.null(object[["s.train"]])) {
        stop(
          "posterior predictive sampling is not available on a ",
          "heteroscedastic fit whose sampler replays no variance surface"
        )
      }
      # ppd sampling pairs one draw of noise with one draw of ev, one per
      # posterior sample; without the tree store, object$fit$predict above
      # replayed only the current trees (one evaluation, not one per draw),
      # the same shape refuseWithoutTrees's callers see, so refuse the same
      # way rather than let the length mismatch reach rnorm/rep_len below.
      # Checked after the heteroscedastic-specific stops above, which name
      # the more informative cause (keepFits) when either one applies.
      if (!object$fit$control@keepTrees) {
        stop(
          "predict requires the fit's saved trees; refit with ",
          bartKeepTreesArgument(object),
          " = TRUE: posterior predictive sampling draws one sample per ",
          "posterior draw, and without the tree store only the current ",
          "trees replay, one evaluation for every draw"
        )
      }
      result <- sampleFromPPD(
        result,
        object,
        weights,
        n.chains,
        heteroscedasticScale(s, n.chains)
      )
    }
  }

  # ci.level opts into a per-observation est + credible band (kind follows type)
  if (!is.null(ci.level)) {
    interval <- posteriorInterval(result, ci.level)
    if (!is.null(s)) {
      attr(interval, "s") <- s
    }
    return(padPredictedRows(interval, rows, first = TRUE))
  }

  if (!is.null(s)) {
    attr(result, "s") <- s
  }
  padPredictedRows(result, rows)
}

# extract(type = "trees") rewrites the matched call onto the sampler's
# getTrees(treeNums, chainNums, sampleNums, current, newdata, forest); none of
# the extract method's own vocabulary but 'forest' reaches getTrees (bart.Rd's
# 'Extracting Trees' section documents chainNums/sampleNums/treeNums/newdata/
# forest as accepted there), so a caller-supplied argument that collides by
# name - sample, combineChains, contribution - is refused by name instead of
# being left to partial-match one of getTrees's differently-named formals
# (sample -> sampleNums) or fall through to a raw 'unused argument'.
refuseTreesArguments <- function(treesCall, ownNames) {
  supplied <- intersect(ownNames, names(treesCall))
  if (length(supplied) > 0L) {
    stop(
      "'",
      supplied[1L],
      "' is not used when type = \"trees\"; the sampler's getTrees accepts ",
      "'chainNums', 'sampleNums', 'treeNums', 'current', 'newdata', and ",
      "'forest' instead (see 'Extracting Trees' in ?bart)"
    )
  }
  invisible(NULL)
}

# Tree draws cannot be combined across chains, so extract(type = "trees")
# always reads in the combineChains = FALSE regime and always carries the chain
# margin; the sampler's own getTrees leaves the column off at one chain.
# Placed after the forest column where there is one, else first.
addTreesChainColumn <- function(trees) {
  if (!is.null(trees[["chain"]])) {
    return(trees)
  }
  position <- if (identical(names(trees)[1L], "forest")) 1L else 0L
  chain <- data.frame(chain = rep(1L, nrow(trees)))
  cbind(
    trees[seq_len(position)],
    chain,
    trees[seq_len(ncol(trees) - position) + position]
  )
}

# extract's four scalar-parameter types, served to every fit class that has
# them: a parameter the fit sampled comes back as its draws, in the layout the
# chain margin asks for, and one it held fixed as one number (fixed holds what
# the sampler held; the draw channel of a fixed sigma or shape repeats it).
# k is the sampler's own, and leaf.prior.sd the k.scale over it - the forest
# total's prior sd in the units the forest fits - so it is a fixed number or
# draws exactly as k is. A fit with several forests has one k and one sd per
# forest, named, or the one 'forest' selects. A fit saved before fits stored
# these descriptors answers from its channels where they suffice and refuses by
# name where they do not.
extractParameter <- function(
  object,
  type,
  combineChains,
  forest = NULL
) {
  if (type == "sigma" && !fitHasResidual(object)) {
    return(1)
  }
  n.chains <- fitNChains(object)
  fixed <- object[["fixed"]]
  if (type == "leaf.prior.sd") {
    prior <- object[["leaf.prior"]]
    if (is.null(prior)) {
      stop(
        "cannot extract 'leaf.prior.sd': this fit was saved before fits ",
        "recorded their leaf prior"
      )
    }
    anchor <- if (is.null(prior[["leaf.prior"]])) {
      vapply(prior, function(forestPrior) forestPrior$k.scale, 0)
    } else {
      prior$k.scale
    }
    k <- if (!is.null(fixed[["k"]])) {
      fixed[["k"]]
    } else {
      reshapeScalarChannel(object[["k"]], n.chains, combineChains)
    }
    return(selectForests(anchor / k, forest))
  }
  value <- fixed[[type]]
  if (!is.null(value)) {
    return(selectForests(value, forest))
  }
  channel <- object[[type]]
  if (is.null(channel)) {
    stop(
      "cannot extract '",
      type,
      "': this fit was saved before fits recorded a ",
      type,
      " they held fixed"
    )
  }
  reshapeScalarChannel(channel, n.chains, combineChains)
}

# one number per forest, named, or the forest selected by index or name
selectForests <- function(value, forest) {
  if (is.null(forest)) {
    return(value)
  }
  perForest <- if (is.matrix(value)) rownames(value) else names(value)
  index <- resolveForestSelection(forest, perForest)
  if (is.matrix(value)) {
    return(drop(value[index, , drop = FALSE]))
  }
  chosen <- value[index]
  if (length(chosen) == 1L) unname(chosen) else chosen
}

extract.bart <- function(
  object,
  type = c(
    "ev",
    "ppd",
    "bart",
    "loglik",
    "trees",
    "forest",
    "sigma",
    "k",
    "leaf.prior.sd",
    "varcount"
  ),
  sample = c("train", "test"),
  combineChains = TRUE,
  forest = NULL,
  contribution = FALSE,
  ...
) {
  type <- validateType(type, eval(formals(extract.bart)$type))
  sampleSupplied <- !missing(sample)

  if (type == "trees") {
    if (is.null(object$fit)) {
      refuseWithoutTrees(
        "extract(type = \"trees\")",
        bartKeepTreesArgument(object)
      )
    }
    treesCall <- match.call()
    refuseTreesArguments(
      treesCall,
      c("sample", "combineChains", "contribution")
    )
    target <- quote(object$fit$getTrees)
    target[[2L]][[2L]] <- treesCall$object
    treesCall[[1L]] <- target
    treesCall$object <- NULL
    treesCall$type <- NULL
    return(addTreesChainColumn(eval(treesCall, parent.frame())))
  }

  # below the type == "trees" branch and its own refuseTreesArguments, so
  # extract(type = "trees", newdata = ) keeps forwarding to getTrees instead
  # of being refused here for a name that arm alone accepts
  refuseUnusedGenericArgs(
    list(...),
    "extract",
    "bart",
    foreignArgsFor(extractForeignReasons, names(formals(extract.bart)))
  )

  refuseForestSelectionOutsideForestArm(
    type,
    forest,
    fitIsHeteroscedastic(object),
    fitNumForests(object)
  )
  if (type != "forest" && isTRUE(contribution)) {
    stop(
      "type = \"",
      type,
      "\" does not support 'contribution': the ",
      "per-observation decomposition applies to the per-forest channel alone"
    )
  }

  # a heteroscedastic fit's scale is per observation, so it reads like a
  # fitted channel: train or test, chains split or combined
  if (type == "sigma" && fitIsHeteroscedastic(object)) {
    sample <- validateSample(sample, eval(formals(extract.bart)$sample))
    s <- object[[if (sample == "train") "s.train" else "s.test"]]
    # only bart() fits a variance forest, and only its keepFits drops s(x):
    # keepTrainingFits leaves s.train in place
    if (is.null(s)) {
      stop(
        "cannot extract 'sigma' at the ",
        sample,
        " rows: this heteroscedastic fit stores no per-observation scale ",
        "draws there (",
        if (sample == "test") "no test rows, or ",
        "'keepFits = FALSE' dropped them, as a supplied 'callback' does ",
        "unless 'keepFits' is given)"
      )
    }
    return(combineOrUncombineChains(s, fitNChains(object), combineChains))
  }

  # served before any sample/test-channel check, so a fit kept with
  # keepTrainingFits = FALSE still serves sigma
  if (type %in% c("sigma", "k", "leaf.prior.sd", "varcount")) {
    refuseSampleOnModelType(type, sampleSupplied)
    if (type == "varcount") {
      trailing <- if (fitNumForests(object) > 1L) 2L else 1L
      return(reshapeChainedChannel(
        object$varcount,
        fitNChains(object),
        combineChains,
        trailing
      ))
    }
    return(extractParameter(object, type, combineChains, forest))
  }

  sample <- validateSample(sample, eval(formals(extract.bart)$sample))

  if (type == "forest") {
    return(extractForest(object, sample, combineChains, forest, contribution))
  }

  # the log-likelihood is against the stored training response; there is no
  # test response to evaluate
  if (type == "loglik" && sample == "test") {
    stop("cannot extract a test sample log-likelihood; no test response exists")
  }

  if (sample == "test" && is.null(object[["yhat.test"]])) {
    stop(
      "cannot extract test sample predictions: either no test data exists ",
      "(use 'predict' instead), or the fit was run with 'keepFits' == ",
      "FALSE (set automatically when 'callback' is supplied, unless ",
      "overridden), which drops the channel even when test data exists"
    )
  }
  if (sample == "train" && is.null(object[["yhat.train"]])) {
    if (callName(object$call) == "bartBT") {
      stop(
        "cannot extract train sample predictions; bartBT must be called with 'keeptrainfits' == TRUE"
      )
    } else {
      stop(
        "cannot extract train sample predictions; bart must be called with ",
        "'keepTrainingFits' == TRUE and 'keepFits' == TRUE (the latter set ",
        "FALSE automatically when 'callback' is supplied, unless overridden)"
      )
    }
  }

  n.chains <- if (!is.null(object[["fit"]])) {
    object$fit$control@n.chains
  } else {
    object$n.chains
  }

  if (type == "loglik") {
    ev <- extract.bart(
      object,
      type = "ev",
      sample = "train",
      combineChains = FALSE
    )
    result <- pointwiseLogLikelihood(object, ev)
    return(combineOrUncombineChains(result, n.chains, combineChains))
  }

  result <- if (sample == "train") object$yhat.train else object$yhat.test
  weights <- if (sample == "train") object$weights else object$weights.test

  result <- combineOrUncombineChains(result, n.chains, combineChains)

  if (type == "bart") {
    return(result)
  }

  if (fitIsBinary(object)) {
    result <- probabilityFromLatents(result, object)
  }

  if (type == "ppd") {
    s <- if (sample == "train") object[["s.train"]] else object[["s.test"]]
    # resid.scale survives keepFits = FALSE where object$s.train would not,
    # so this catches a heteroscedastic fit whose s.train/s.test keepFits
    # dropped before the narrower "no s.test at all" check below
    if (
      is.null(s) &&
        is.null(object[["s.train"]]) &&
        fitIsHeteroscedastic(object)
    ) {
      stop(
        "posterior predictive sampling needs this heteroscedastic fit's ",
        "'s.train'/'s.test' draws, which 'keepFits = FALSE' dropped ",
        "(automatically, when a 'callback' was supplied, unless ",
        "overridden); refit with 'keepFits = TRUE'"
      )
    }
    if (is.null(s) && !is.null(object[["s.train"]])) {
      stop(
        "posterior predictive sampling is not available at the test rows of a ",
        "heteroscedastic fit that stores no 's.test' draws"
      )
    }
    result <- sampleFromPPD(
      result,
      object,
      weights,
      n.chains,
      heteroscedasticScale(s, n.chains)
    )
  }

  result
}

# Selects one or more forests from a per-forest channel's trailing margin, by
# 1-based index or by the shipped forest1..forestK vocabulary; NULL
# selects every forest, in margin order. A declaration's own forest.labels are
# not a selector - they are a display attribute, not a second vocabulary.
resolveForestSelection <- function(forest, forestNames) {
  if (is.null(forest)) {
    return(seq_along(forestNames))
  }
  if (is.character(forest)) {
    idx <- match(forest, forestNames)
    if (anyNA(idx)) {
      stop(
        "'forest' must name one of '",
        paste0(forestNames, collapse = "', '"),
        "'"
      )
    }
    return(idx)
  }
  idx <- coerceOrError(forest, "integer")
  if (anyNA(idx) || any(idx < 1L | idx > length(forestNames))) {
    stop("'forest' index must be between 1 and ", length(forestNames))
  }
  idx
}

# keepFits = FALSE drops forestFits but not the n.forests descriptor, which is
# set off the fit's own bases: a fit with more than one forest and no
# forestFits is therefore an amplitude-coupled one whose per-forest channel
# was opted out. Every arm that
# reads the channel names that here, rather than falling through to "this fit
# has none" (which is wrong - it had one) or, on the combined arm, to the
# engine's own off-sample refusal, which names neither keepFits nor callback.
refuseDroppedForestChannel <- function(object) {
  if (is.null(object[["forestFits"]]) && fitNumForests(object) > 1L) {
    stop(
      "this amplitude-coupled fit's per-forest channel was dropped by ",
      "'keepFits' == FALSE (set automatically when 'callback' is supplied, ",
      "unless overridden); the combined surface is rebuilt from that ",
      "channel, so refit with 'keepFits = TRUE'"
    )
  }
}

# extract(type = "forest"): the packaged per-forest response-scale raw total
# by default (forestFits already carries response.scale), or its
# per-observation contribution under contribution = TRUE, computed on
# demand as (basis %*% glue) * raw rather than stored. The selected forests
# always keep the trailing forest margin, even at length one. Refuses by name
# on a fit without forest reporting (the amplitude coupling, not the forest
# count) and on sample = "test" (an amplitude-coupled fit has no test fits).
extractForest <- function(object, sample, combineChains, forest, contribution) {
  refuseDroppedForestChannel(object)
  if (is.null(object[["forestFits"]])) {
    stop(
      "type = \"forest\" is only available on a fit with per-forest ",
      "reporting (an amplitude-coupled multi-forest fit); this fit has none"
    )
  }
  if (sample == "test") {
    stop(
      "type = \"forest\" does not support sample = \"test\": no test-sample ",
      "per-forest channel is stored, since an amplitude-coupled fit takes no ",
      "test predictors; predict(type = \"forest\") replays the forests at new ",
      "rows instead"
    )
  }
  n.chains <- if (!is.null(object[["fit"]])) {
    object$fit$control@n.chains
  } else {
    object$n.chains
  }

  fits <- reshapeChainedChannel(object$forestFits, n.chains, TRUE, 2L)
  forestNames <- dimnames(fits)[[3L]]
  idx <- resolveForestSelection(forest, forestNames)

  if (!contribution) {
    result <- fits[,, idx, drop = FALSE]
    return(reshapeChainedChannel(result, n.chains, combineChains, 2L))
  }

  glue <- reshapeChainedChannel(object$glue, n.chains, TRUE, 1L)
  glueForest <- attr(object$glue, "forest")
  n.obs <- dim(fits)[2L]
  result <- array(
    0,
    c(dim(fits)[1L], n.obs, length(idx)),
    dimnames = list(NULL, dimnames(fits)[[2L]], forestNames[idx])
  )
  for (j in seq_along(idx)) {
    k <- idx[j]
    basis <- object$bases[[k]]
    if (is.null(basis)) {
      basis <- matrix(1, n.obs, 1L)
    }
    g <- glue[, glueForest == forestNames[k], drop = FALSE]
    result[,, j] <- (g %*% t(basis)) * fits[,, k]
  }
  reshapeChainedChannel(result, n.chains, combineChains, 2L)
}

# predict(type = "forest"): the out-of-sample twin of extract(type = "forest")'s
# raw slice - each selected forest's own RESPONSE-scale total at newdata,
# replayed from the saved trees, which predict.bart requires keepTrees for on
# this arm (the sampler method under it still reads the current trees when
# there is no store). The engine reports the forests on their internal
# scale, so response.scale is applied here exactly as packageBartResults applies
# it to the in-sample channel; the amplitude glue, the response shift and any
# offset are deliberately NOT folded in, because the recombination needs the
# caller's own bases at the new rows (man/bart.Rd states the idiom). Refuses by
# name on a fit without per-forest reporting, off the same stored channel
# extract reads, and there is no contribution = arm for the same reason.
predictForest <- function(
  object,
  newdata,
  offset,
  combineChains,
  forest,
  n.threads,
  rowNames = NULL
) {
  if (is.null(object[["forestFits"]])) {
    stop(
      "type = \"forest\" is only available on a fit with per-forest ",
      "reporting (an amplitude-coupled multi-forest fit); this fit has none"
    )
  }
  n.chains <- object$fit$control@n.chains
  responseScale <- object$fit$getLeafPrior(1L)$response.scale
  raw <- predictForestsCodedTest(object$fit, newdata, offset, n.threads) *
    responseScale
  # forestFits carries the fit's own combineChains shape (3-d combined, 4-d
  # split across chains), so the forest margin is always the LAST axis rather
  # than a fixed index
  forestNames <- dimnames(object$forestFits)[[length(dim(object$forestFits))]]
  idx <- resolveForestSelection(forest, forestNames)
  result <- shapeMultinomialChannel(
    raw,
    forestNames,
    n.chains,
    TRUE,
    leadNames = rowNames
  )
  reshapeChainedChannel(
    result[,, idx, drop = FALSE],
    n.chains,
    combineChains,
    2L
  )
}

# The per-forest bases at the PREDICTED rows. A caller's own 'bases' wins
# everywhere; failing that, a forest() term's stored formula is re-evaluated
# against newdata, and a fit whose bases arrived as raw values has nothing to
# re-evaluate and must be given them by name. A bare (non-list) value positions
# itself when exactly one forest carries a basis - the Bayesian causal forest
# call, bases = <arm at the new rows> - and is refused as ambiguous otherwise.
# Widths are checked against the FIT's own bases rather than left to %*% to
# recycle: amplitude j multiplies column j, so a width that drifts is a wrong
# answer rather than an error.
resolveForestBases <- function(object, bases, newdata, n.new) {
  storedBases <- object$bases
  numForests <- length(storedBases)
  carriers <- which(!vapply(storedBases, is.null, logical(1L)))
  if (is.null(bases)) {
    terms <- object[["basis.terms"]]
    bases <- if (is.null(terms)) {
      vector("list", numForests)
    } else {
      lapply(seq_len(numForests), function(k) {
        if (is.null(terms[[k]])) {
          NULL
        } else {
          replayForestBasis(terms[[k]], newdata, k)
        }
      })
    }
  } else {
    if (!is.list(bases) || is.data.frame(bases)) {
      if (length(carriers) != 1L) {
        stop(
          "'bases' takes a bare value only when exactly one forest carries a ",
          "basis; ",
          length(carriers),
          " of this fit's ",
          numForests,
          " forests do - give a length-",
          numForests,
          " list instead"
        )
      }
      value <- bases
      bases <- vector("list", numForests)
      bases[[carriers]] <- value
    }
    if (length(bases) != numForests) {
      stop(
        "'bases' must be a length-",
        numForests,
        " list, one entry per forest (NULL for a forest that declares none); ",
        "got ",
        length(bases)
      )
    }
  }
  bases <- lapply(bases, expandForestBasis, atPrediction = TRUE)
  bases <- validateForestBases(
    bases,
    n.new,
    argument = "bases",
    rows = "'newdata'"
  )
  for (k in seq_len(numForests)) {
    width <- if (is.null(storedBases[[k]])) 0L else ncol(storedBases[[k]])
    if (width == 0L) {
      if (!is.null(bases[[k]])) {
        stop(
          "'bases' gives forest ",
          k,
          " a basis, which it declares none of: its single amplitude ",
          "multiplies an implicit all-ones column"
        )
      }
    } else if (is.null(bases[[k]])) {
      stop(
        "forest ",
        k,
        " carries a basis, so the blend needs its ",
        width,
        if (width == 1L) " column" else " columns",
        " at the new rows: give them through 'bases =' (a length-",
        numForests,
        " list, or the bare value when only one forest carries a basis)"
      )
    } else if (ncol(bases[[k]]) != width) {
      stop(
        "'bases' gives forest ",
        k,
        " ",
        ncol(bases[[k]]),
        if (ncol(bases[[k]]) == 1L) " column" else " columns",
        "; its amplitudes take ",
        width
      )
    }
  }
  bases
}

# predict(type = "ev"/"ppd"/"bart") on an amplitude-coupled fit: the
# recombination predict(type = "forest") deliberately leaves out, performed
# here because at THIS level the bases at the predicted rows are available -
# the caller's own, or a forest() term's formula re-evaluated against newdata -
# where the sampler holds none. eta = response.shift +
# sum_k (glue_k %*% t(basis_k)) * forest_k + offset, the identity the packaged
# yhat.train satisfies in sample, with the family's link applied after (so
# "bart" is eta itself, as it is for a single forest) and "ppd" feeding the
# unchanged sampleFromPPD. The glue and the replay pair draw for draw only in
# the combined, chain-major layout both are stated in, so the whole
# accumulation runs there and the result is split at the end.
predictBlend <- function(
  object,
  newdata,
  offset,
  weights,
  type,
  combineChains,
  ci.level,
  bases,
  n.threads,
  rowNames = NULL,
  rawNewdata = newdata
) {
  n.chains <- object$fit$control@n.chains
  perForest <- predictForest(object, newdata, NULL, TRUE, NULL, n.threads)
  n.new <- dim(perForest)[2L]
  forestNames <- dimnames(perForest)[[3L]]
  bases <- resolveForestBases(object, bases, rawNewdata, n.new)

  # the caller's own offset and weights at those rows, read as the sampler's
  # own predict reads them: numeric, length-1 recycled or one per row
  if (!is.null(offset)) {
    offset <- as.double(offset)
    if (length(offset) == 1L) {
      offset <- rep_len(offset, n.new)
    }
    if (length(offset) != n.new) {
      stop("'offset' must have the same number of rows as 'newdata'")
    }
  }
  if (!is.null(weights)) {
    weights <- as.double(weights)
    if (length(weights) == 1L) {
      weights <- rep_len(weights, n.new)
    }
    if (length(weights) != n.new) {
      stop("'weights' must have the same number of rows as 'newdata'")
    }
  }

  glue <- reshapeChainedChannel(object$glue, n.chains, TRUE, 1L)
  # the ragged margin's forest key rides the STORED channel, which the reshape
  # above does not carry forward
  glueForest <- attr(object$glue, "forest")
  result <- matrix(
    object$fit$getLeafPrior(1L)$response.shift,
    nrow(glue),
    n.new
  )
  for (k in seq_along(forestNames)) {
    basis <- bases[[k]]
    if (is.null(basis)) {
      basis <- matrix(1, n.new, 1L)
    }
    g <- glue[, glueForest == forestNames[k], drop = FALSE]
    result <- result + (g %*% t(basis)) * perForest[,, k]
  }
  if (!is.null(offset)) {
    result <- result + rep(offset, each = nrow(result))
  }
  result <- nameObservationMargin(result, rowNames)
  result <- combineOrUncombineChains(result, n.chains, combineChains)

  if (type != "bart") {
    if (fitIsBinary(object)) {
      result <- probabilityFromLatents(result, object)
    }
    if (type == "ppd") {
      result <- sampleFromPPD(result, object, weights, n.chains)
    }
  }

  if (!is.null(ci.level)) {
    return(posteriorInterval(result, ci.level))
  }
  result
}

fitted.bart <- function(
  object,
  type = c("ev", "ppd", "bart"),
  sample = c("train", "test"),
  ci.level = NULL,
  ...
) {
  type <- validateType(type, eval(formals(fitted.bart)$type))
  sample <- validateSample(sample, eval(formals(fitted.bart)$sample))
  refuseUnusedGenericArgs(
    list(...),
    "fitted",
    "bart",
    c(
      bartUnusedArgs,
      foreignArgsFor(fittedForeignReasons, names(formals(fitted.bart)))
    )
  )

  result <- extract(object, type, sample)

  # a training-side quantity pads back to the caller's own row count through
  # whatever the fit's na.action recorded; the test side never lost a row
  padded <- if (identical(sample, "train")) {
    function(value) padOmittedRows(object[["na.action"]], value)
  } else {
    identity
  }

  # ci.level opts into a per-observation est + credible band instead of the
  # posterior mean; the interval kind follows type (see posteriorInterval)
  if (!is.null(ci.level)) {
    return(padded(posteriorInterval(result, ci.level)))
  }

  if (!is.null(dim(result))) {
    padded(channelMeans(result))
  } else {
    padded(mean(result))
  }
}

# residuals are always against the training response, so a caller-supplied
# 'sample' collides with the fixed sample = "train" residuals.* passes to
# fitted - refuse it by name rather than let it reach fitted's own 'sample'
# formal twice and raise a raw 'formal argument "sample" matched by multiple
# actual arguments'.
refuseResidualsSample <- function(dots) {
  if ("sample" %in% names(dots)) {
    stop(
      "'sample' is not used by residuals: residuals are always against the ",
      "training response"
    )
  }
  invisible(NULL)
}

residuals.bart <- function(object, type = "ev", ...) {
  # type flows to fitted so link-scale (type = "bart") residuals are reachable;
  # residuals are always against the training response, so sample is pinned
  refuseResidualsSample(list(...))
  refuseUnusedGenericArgs(
    list(...),
    "residuals",
    "bart",
    c(
      bartUnusedArgs,
      foreignArgsFor(residualsForeignReasons, names(formals(residuals.bart)))
    )
  )
  # the response the fit kept, padded the same way the fitted values are, so
  # the two line up row for row at the caller's own shape
  padOmittedRows(object[["na.action"]], object$y) -
    fitted.bart(object, type = type, sample = "train")
}

# bart2(family = "multinomial") generics. The
# fit object is class "bartMultinomial" - deliberately NOT "bart" - so it
# never falls through to the "bart" methods above: those assume an n x
# samples (x chains) shape with no K margin and would silently misread the
# K-widened arrays here rather than error.
#
# type = "ev" is the engine's own train channel: already softmax
# PROBABILITIES, so unlike the binary
# families there is no latent-to-probability transform. type = "bart" (the
# latent scale) is refused: the run records only the identified
# probabilities, and the raw per-category fits are non-identified and
# unrecorded. type = "ppd" draws one category per posterior draw from its
# probability vector, returned as integer codes (1-based, indexing
# object$levels) in an array shaped like "ev" minus the K margin - the same
# "ppd keeps ev's shape" convention the binary families use. type flows
# through validateType, as it does on a "bart" fit, so "response" and "link"
# are the predict.glm synonyms here too.

# The two latent-scale requests a multinomial fit refuses, by name and with
# the reason, from extract, fitted and predict alike: type = "bart", the raw
# per-category fits the run does not record, and type = "forest", the same
# fits replayed at new rows. Both are the one non-identification: the softmax
# is invariant to a common per-observation shift, so each row's level is free,
# and it is not noise either - the backfit reproduces it as a function of x -
# so a latent surface would be read as signal it is not. What IS identified is
# the log-ratio, which the logs of the reported probabilities carry exactly.
refuseMultinomialLatentType <- function(type) {
  if (type == "bart") {
    stop(
      "multinomial fits do not support type = \"bart\": the run records ",
      "only the identified softmax probabilities; the raw per-category ",
      "latent fits are non-identified and unrecorded"
    )
  }
  if (type == "forest") {
    stop(
      "multinomial fits do not support type = \"forest\": a category's ",
      "forest is a latent whose level is reproducibly structured yet not ",
      "identified, so a raw replay reads as signal; the identified content ",
      "is the log-ratio of the probabilities predict() reports"
    )
  }
  invisible(NULL)
}

# The posterior-mean n x K probability matrix of a K-widened draws array
# (observation margin next-to-last, category margin last in every chain
# layout), and its argmax as a factor over the fit's own levels - the class
# prediction fitted() and predict() share, so the two cannot drift.
meanCategoryProbabilities <- function(probs, levels) {
  meanProbs <- channelMeans(probs, 2L)
  dimnames(meanProbs) <- list(rownames(meanProbs), levels)
  meanProbs
}
categoryFromMeanProbabilities <- function(meanProbs, levels, ordered = FALSE) {
  nameObservationMargin(
    factor(
      levels[max.col(meanProbs, ties.method = "first")],
      levels = levels,
      ordered = ordered
    ),
    rownames(meanProbs)
  )
}

# 'forest'/'contribution' select among an amplitude-coupled fit's co-fit
# forests (extract.bart's own vocabulary); every own-class fit but
# bartMultinomial has a single forest per component to begin with, so both
# names refuse for the same reason there.
singleForestReason <- paste0(
  "this selects among an amplitude-coupled fit's co-fit forests; this fit ",
  "has a single forest"
)

# bart-family arguments this fit's K-widened shape has no room for: a single
# category's forest is not identified individually (refuseMultinomialLatentType
# already refuses type = "forest"; 'forest'/'contribution' are refused here so
# passing them does not silently vanish into '...' regardless of type).
multinomialUnusedArgs <- list(
  forest = paste0(
    "a multinomial fit's K category forests are not identified individually ",
    "- the identified content is the reported probabilities"
  ),
  contribution = paste0(
    "a multinomial fit's K category forests are not identified individually ",
    "- the identified content is the reported probabilities"
  )
)

extract.bartMultinomial <- function(
  object,
  type = c(
    "ev",
    "ppd",
    "bart",
    "forest",
    "loglik",
    "sigma",
    "k",
    "leaf.prior.sd",
    "varcount",
    "trees"
  ),
  sample = c("train", "test"),
  combineChains = TRUE,
  ...
) {
  type <- validateType(type, eval(formals(extract.bartMultinomial)$type))
  sampleSupplied <- !missing(sample)
  refuseMultinomialLatentType(type)

  # unlike 'bart'/'forest' above, a category's TREES are recorded and
  # identified (only the per-category level, not the structure, is
  # non-identified), so this forwards to the K-forest sampler's own getTrees
  # exactly as extract.bart does, ahead of multinomialUnusedArgs' blanket
  # 'forest' refusal below
  if (type == "trees") {
    if (is.null(object$fit)) {
      refuseWithoutTrees(
        "extract(type = \"trees\")",
        bartKeepTreesArgument(object)
      )
    }
    treesCall <- match.call()
    refuseTreesArguments(
      treesCall,
      c("sample", "combineChains", "contribution")
    )
    target <- quote(object$fit$getTrees)
    target[[2L]][[2L]] <- treesCall$object
    treesCall[[1L]] <- target
    treesCall$object <- NULL
    treesCall$type <- NULL
    return(addTreesChainColumn(eval(treesCall, parent.frame())))
  }

  refuseUnusedGenericArgs(
    list(...),
    "extract",
    "bartMultinomial",
    c(
      multinomialUnusedArgs,
      foreignArgsFor(
        extractForeignReasons,
        names(formals(extract.bartMultinomial))
      )
    )
  )

  if (type %in% c("sigma", "k", "leaf.prior.sd")) {
    refuseSampleOnModelType(type, sampleSupplied)
    return(extractParameter(object, type, combineChains))
  }

  if (type == "varcount") {
    refuseSampleOnModelType(type, sampleSupplied)
    return(reshapeChainedChannel(
      object$varcount,
      fitNChains(object),
      combineChains,
      2L
    ))
  }

  sample <- validateSample(
    sample,
    eval(formals(extract.bartMultinomial)$sample)
  )

  if (type == "loglik" && sample == "test") {
    stop("cannot extract a test sample log-likelihood; no test response exists")
  }

  probs <- if (sample == "test") {
    if (is.null(object$yhat.test)) {
      stop(
        "this multinomial fit carries no test channel; refit with 'test' ",
        "to report out-of-sample softmax probabilities"
      )
    }
    object$yhat.test
  } else {
    object$yhat.train
  }
  n.chains <- fitNChains(object)

  if (type == "loglik") {
    result <- multinomialLogLik(
      object,
      reshapeChainedChannel(probs, n.chains, FALSE, 2L)
    )
    return(combineOrUncombineChains(result, n.chains, combineChains))
  }

  probs <- reshapeChainedChannel(probs, n.chains, combineChains, 2L)
  if (type == "ev") {
    return(probs)
  }
  # a count-row fit's own rows draw what it modelled, a count vector of each
  # row's trials; test rows have no trial count, so they draw one category
  # per draw, as predict does
  if (sample == "train" && !is.factor(object[["y"]])) {
    return(multinomialCountPpdFromProbs(probs, rowSums(object[["y"]])))
  }
  multinomialPpdFromProbs(probs)
}

# ONE formula covers both response ingestions: with p[s,i,k] the reported
# probability and n_i = sum_k y_ik (= 1 for a labeled response), the log
# density of the observed row is the multinomial log-pmf including its
# coefficient (as dmultinom reports it), which reduces to log(p[s,i,y_i]) when
# n_i = 1. The likelihood unit is the observation ROW (n_i trials), not a
# single trial or a (row, category) cell, so loo/WAIC on this channel is
# leave-one-row-out. probs enters in the split (chains x) samples x obs x K
# layout; the result drops the K margin (dim(ev) minus its trailing margin).
multinomialLogLik <- function(object, probs) {
  y <- object[["y"]]
  levels <- object[["levels"]]
  d <- dim(probs)
  K <- d[length(d)]
  nObs <- d[length(d) - 1L]
  counts <- if (is.factor(y)) {
    indicator <- matrix(0, length(y), K)
    indicator[cbind(seq_along(y), match(y, levels))] <- 1
    indicator
  } else {
    y
  }
  n <- rowSums(counts)
  logCoef <- lgamma(n + 1) - rowSums(lgamma(counts + 1))
  n.draws <- length(probs) %/% (nObs * K)
  flat <- probs
  dim(flat) <- c(n.draws * nObs, K)
  idx <- rep(seq_len(nObs), each = n.draws)
  # a zero-count cell contributes 0 even where its probability underflows to
  # 0, as in dmultinom, so a row with no trial is exactly 0
  expanded <- counts[idx, , drop = FALSE]
  cells <- expanded * log(flat)
  cells[expanded == 0] <- 0
  term <- rowSums(cells)
  array(
    rep(logCoef, each = n.draws) + term,
    d[-length(d)],
    dimnames(probs)[-length(d)]
  )
}

# shared by extract.bartMultinomial (stored channels) and
# predict.bartMultinomial (freshly replayed channels), so both draw
# categories from probabilities the identical way: K rides the trailing
# dimension already, so reinterpreting the same flat storage as a
# (draws * obs) x K matrix needs no permutation. codes are 1-based, indexing
# 'levels', in an array shaped like probs minus the K margin.
multinomialPpdFromProbs <- function(probs) {
  d <- dim(probs)
  K <- d[length(d)]
  flat <- probs
  dim(flat) <- c(prod(d[-length(d)]), K)
  codes <- apply(flat, 1L, function(p) sample.int(K, 1L, prob = p))
  array(codes, d[-length(d)], dimnames(probs)[-length(d)])
}

# The count-row counterpart: a Multinomial(n_i, p) count vector per draw and
# row, laid out as probs is (K trailing), a row of zero trials all zeros. Drawn
# as rmultinom draws it, by sequential binomials - category k takes
# Binomial(the trials left, p_k / the probability left) - but one category at
# a time across every (draw, row) at once rather than one call per row.
multinomialCountPpdFromProbs <- function(probs, trials) {
  d <- dim(probs)
  K <- d[length(d)]
  nObs <- d[length(d) - 1L]
  n.draws <- length(probs) %/% (nObs * K)
  flat <- probs
  dim(flat) <- c(n.draws * nObs, K)
  remaining <- rep(as.integer(trials), each = n.draws)
  probabilityLeft <- rowSums(flat)
  counts <- matrix(0L, n.draws * nObs, K)
  for (k in seq_len(K - 1L)) {
    conditional <- ifelse(
      probabilityLeft > 0,
      pmin(1, flat[, k] / probabilityLeft),
      0
    )
    counts[, k] <- stats::rbinom(length(remaining), remaining, conditional)
    remaining <- remaining - counts[, k]
    probabilityLeft <- pmax(0, probabilityLeft - flat[, k])
  }
  counts[, K] <- remaining
  array(counts, d, dimnames(probs))
}

# The posterior-mean n x K probability matrix (colnames = levels(y)), or
# (type = "class") the argmax category of that mean as a factor over the
# original levels - the class-prediction convenience. ci.level opts into a
# per-(observation, category) credible band instead of the posterior mean,
# taken on the full probability draws before the class reduction so it is
# meaningful regardless of 'type'.
fitted.bartMultinomial <- function(
  object,
  type = c("ev", "class", "bart"),
  sample = c("train", "test"),
  ci.level = NULL,
  ...
) {
  type <- validateType(type, eval(formals(fitted.bartMultinomial)$type))
  sample <- validateSample(
    sample,
    eval(formals(fitted.bartMultinomial)$sample)
  )
  refuseMultinomialLatentType(type)
  refuseUnusedGenericArgs(
    list(...),
    "fitted",
    "bartMultinomial",
    c(
      multinomialUnusedArgs,
      foreignArgsFor(
        fittedForeignReasons,
        names(formals(fitted.bartMultinomial))
      )
    )
  )
  refuseClassCiLevel(type, ci.level)
  # a training-side quantity pads back to the caller's own row count through
  # whatever the fit's na.action recorded; the test side never lost a row
  padded <- if (identical(sample, "train")) {
    function(value) padOmittedRows(object[["na.action"]], value)
  } else {
    identity
  }
  probs <- extract.bartMultinomial(object, type = "ev", sample = sample)
  if (!is.null(ci.level)) {
    return(padded(posteriorInterval(probs, ci.level, trailing = 2L)))
  }
  meanProbs <- meanCategoryProbabilities(probs, object$levels)
  if (type == "ev") {
    return(padded(meanProbs))
  }
  padded(categoryFromMeanProbabilities(meanProbs, object$levels))
}

# residuals.bart is y - fitted() on the response scale; a multinomial fit has
# no single scalar response to subtract from, so the per-category analog is
# the observed proportion minus the fitted probability, an n x K matrix
# (columns named by 'levels'). For the labeled-response ingestion
# (bart2Multinomial) the observed proportion is the 1[y = k] indicator; for
# the grouped-count ingestion (bart2MultinomialCounts) it is y / rowSums(y),
# which reduces to the same indicator when every row is a single trial. There
# is no other residual to choose, so 'type' is refused by name.
multinomialResidualsTypeReason <- list(
  type = paste0(
    "the residual is the per-category observed proportion minus the ",
    "fitted probability"
  )
)

residuals.bartMultinomial <- function(object, ...) {
  refuseResidualsSample(list(...))
  refuseUnusedGenericArgs(
    list(...),
    "residuals",
    "bartMultinomial",
    c(
      multinomialUnusedArgs,
      multinomialResidualsTypeReason,
      foreignArgsFor(
        residualsForeignReasons,
        names(formals(residuals.bartMultinomial))
      )
    )
  )
  # phat is already padded (fitted's own train-side na.action rule); build
  # 'observed' at the fit's own (unpadded) length and pad it the same way,
  # so the two align before differencing and before phat's own dimnames -
  # already at the padded length - are copied onto the result
  phat <- fitted.bartMultinomial(object, type = "ev")
  y <- object$y
  observed <- if (is.factor(y)) {
    indicator <- matrix(0, length(y), length(object$levels))
    indicator[cbind(seq_along(y), match(y, object$levels))] <- 1
    indicator
  } else {
    # a row with no trial has no observed proportion; as glm's response
    # residual at a zero-weight row, it is observed 0, so its residual is -p
    trials <- rowSums(y)
    y / ifelse(trials == 0, 1, trials)
  }
  observed <- padOmittedRows(object[["na.action"]], observed)
  result <- observed - phat
  dimnames(result) <- dimnames(phat)
  result
}

# Out-of-sample softmax probabilities by replaying the K forests' saved
# trees. Requires a fit kept with keepTrees: a kept
# $fit alone is not enough, since a sampler kept ONLY via keepSampler carries
# no saved trees to replay. Returns a levels-named (n.chains x) n.samples x
# n.new x K probability array, the yhat.test/train convention. type = "bart"
# (the raw per-category latent scale) stays unavailable, as it is for
# extract: only the identified probabilities are recoverable. type = "ppd"
# draws one category per posterior draw from that probability vector via the
# exact same construction extract.bartMultinomial's ppd uses
# (multinomialPpdFromProbs), so the two agree on semantics and encoding; it
# is the only branch that touches the RNG, so the default type = "ev" is
# unchanged and draw-neutral. The replay reads through $fit's own pointer:
# $fit is the K-forest sampler that actually ran, so getPointer() can
# re-create it from stored state after a save/reload.
#
# offset is the per-category shift at the PREDICTED rows, the same name
# predict.bart uses for its own new-row shift (R/generics.R's predict.bart):
# an nNew x K matrix entering the raw fits before the softmax. It is never
# taken from the fit, because these rows are not the fit's rows, so a fit
# trained under a category offset requires one here rather than being served
# the offset-free surface by default - an all-zero matrix asks for that
# surface on purpose. Passing the training offset back at the training rows
# reproduces yhat.train.
predict.bartMultinomial <- function(
  object,
  newdata,
  type = c("ev", "ppd", "bart", "forest", "class"),
  offset = NULL,
  combineChains = TRUE,
  ci.level = NULL,
  na.action = dbarts::na.keepPredictors,
  n.threads = object$fit$control@n.threads,
  ...
) {
  type <- validateType(type, eval(formals(predict.bartMultinomial)$type))
  refuseMultinomialLatentType(type)
  refuseUnusedGenericArgs(
    list(...),
    "predict",
    "bartMultinomial",
    c(
      multinomialUnusedArgs,
      predictOffsetUnusedArgs,
      foreignArgsFor(
        predictForeignReasons,
        names(formals(predict.bartMultinomial))
      )
    )
  )
  warnUnusedDots(list(...), "predict", "bartMultinomial")
  refuseNonNumericOffset(offset)
  refuseClassCiLevel(type, ci.level)
  if (is.null(object[["fit"]]) || !object$fit$control@keepTrees) {
    refuseWithoutTrees("predict")
  }
  # after the fit check, whose absence the default here would otherwise report
  # as a missing slot
  n.threads <- validatePredictThreads(n.threads)
  offset <- alignCategoryColumns(
    offset,
    object$levels,
    "offset",
    identical(object$levels.source, "index")
  )
  # a missing offset row is incomplete the same way an unroutable predictor
  # is (dec-A89)
  rows <- preparePredictRows(
    newdata,
    object$fit$data@x,
    na.action,
    list(offset = offset)
  )
  if (isTRUE(rows$placeholder)) {
    restoreSeed <- protectRandomSeed()
    on.exit(restoreSeed(), add = TRUE)
  }
  newdata <- rows$x
  offset <- subsetPredictInput(
    offset,
    rows,
    "offset",
    if (!is.null(object$fit$data@offset.category)) matrix(0, 1L, object$K)
  )
  if (is.null(offset)) {
    if (!is.null(object$fit$data@offset.category)) {
      stop(
        "'offset' is required on a multinomial fit trained with a category ",
        "offset: the predicted rows are not the training rows, so pass ",
        "their own ",
        nrow(newdata),
        " x ",
        object$K,
        " matrix, all-zero for the offset-free surface"
      )
    }
  } else {
    offset <- validateCategoryOffset(
      offset,
      nrow(newdata),
      object$K,
      "'offset'"
    )
  }
  # raw is n.new x K x n.samples (x n.chains), the run's test-channel shape
  raw <- predictCodedTest(object$fit, newdata, offset, n.threads)
  probs <- shapeMultinomialChannel(
    raw,
    object$levels,
    object$n.chains,
    combineChains,
    leadNames = rows$keptNames
  )
  if (type == "ppd") {
    probs <- multinomialPpdFromProbs(probs)
  }
  if (!is.null(ci.level)) {
    return(padPredictedRows(
      posteriorInterval(
        probs,
        ci.level,
        trailing = if (type %in% c("ev", "class")) 2L else 1L
      ),
      rows,
      first = TRUE
    ))
  }
  if (type == "class") {
    meanProbs <- meanCategoryProbabilities(probs, object$levels)
    return(padPredictedRows(
      categoryFromMeanProbabilities(meanProbs, object$levels),
      rows
    ))
  }
  padPredictedRows(probs, rows, trailing = if (type == "ppd") 0L else 1L)
}

# Shared "Call:" preamble for the print and summary methods. A fit kept with
# keepCall = FALSE stores no call, and one saved by an earlier version the
# placeholder call("NULL"); either is omitted.
printCall <- function(x) {
  if (is.call(x[["call"]]) && !identical(x[["call"]], call("NULL"))) {
    cat(
      "\nCall:\n",
      paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n",
      sep = ""
    )
  }
  invisible(NULL)
}

print.bartMultinomial <- function(x, ...) {
  printCall(x)
  cat("family: multinomial\n")
  cat("levels: ", paste(x$levels, collapse = ", "), "\n", sep = "")
  cat("n.chains: ", x$n.chains, "\n", sep = "")
  cat("n.trees: ", x$n.trees, "\n", sep = "")
  d <- dim(x$yhat.train)
  # a 4-dim yhat.train (combineChains = FALSE) already separates chains, so
  # d[2L] is per-chain; a 3-dim one (single chain, or combineChains = TRUE,
  # the default) folds the chain margin into d[1L] and must be divided back
  # out
  n.kept <- if (length(d) == 4L) d[2L] else d[1L] %/% x$n.chains
  cat("kept draws (per chain): ", n.kept, "\n", sep = "")
  if (!is.null(x$yhat.test)) {
    dt <- dim(x$yhat.test)
    n.test <- if (length(dt) == 4L) dt[3L] else dt[2L]
    cat("test rows: ", n.test, "\n", sep = "")
  }
  invisible(x)
}

# bart2(family = "ordinal") generics. Like
# bartMultinomial, the fit object is class "bartOrdinal" - never "bart" - so the
# K-widened category-probability arrays never fall through to the single-forest
# "bart" methods. Unlike multinomial, ordinal DOES carry a latent scale: the
# cumulative-probit fits are formed from a single latent eta = f(x), so
# type = "bart"/"link" returns that latent (as it does for probit), while
# type = "ev"/"response" returns the n x K category probabilities computed from
# the latent and the sampled thresholds. type = "ppd" draws one category per
# posterior draw. The K-1 threshold draws ride the fit's $thresholds field.
# ordinal has a single forest, so 'forest'/'contribution' refuse for the same
# reason a bart-family single-forest fit does.
ordinalUnusedArgs <- list(
  forest = singleForestReason,
  contribution = singleForestReason
)

extract.bartOrdinal <- function(
  object,
  type = c(
    "ev",
    "ppd",
    "bart",
    "loglik",
    "thresholds",
    "sigma",
    "k",
    "leaf.prior.sd",
    "varcount"
  ),
  sample = c("train", "test"),
  combineChains = TRUE,
  ...
) {
  type <- validateType(type, eval(formals(extract.bartOrdinal)$type))
  sampleSupplied <- !missing(sample)
  refuseUnusedGenericArgs(
    list(...),
    "extract",
    "bartOrdinal",
    c(
      ordinalUnusedArgs,
      foreignArgsFor(extractForeignReasons, names(formals(extract.bartOrdinal)))
    )
  )
  n.chains <- fitNChains(object)

  if (type %in% c("sigma", "k", "leaf.prior.sd")) {
    refuseSampleOnModelType(type, sampleSupplied)
    return(extractParameter(object, type, combineChains))
  }

  if (type %in% c("thresholds", "varcount")) {
    refuseSampleOnModelType(type, sampleSupplied)
    channel <- if (type == "thresholds") object$thresholds else object$varcount
    return(reshapeChainedChannel(channel, n.chains, combineChains, 1L))
  }

  sample <- validateSample(sample, eval(formals(extract.bartOrdinal)$sample))

  if (type == "loglik") {
    if (sample == "test") {
      stop(
        "cannot extract a test sample log-likelihood; no test response exists"
      )
    }
    result <- ordinalLogLik(
      object,
      reshapeChainedChannel(object$yhat.train, n.chains, FALSE, 2L)
    )
    return(combineOrUncombineChains(result, n.chains, combineChains))
  }

  if (type == "bart") {
    latent <- if (sample == "test") {
      object$latent.test
    } else {
      object$latent.train
    }
    if (is.null(latent)) {
      stop(
        "this ordinal fit carries no test channel; refit with 'test' to ",
        "report out-of-sample latent fits"
      )
    }
    return(combineOrUncombineChains(latent, n.chains, combineChains))
  }
  probs <- if (sample == "test") {
    if (is.null(object$yhat.test)) {
      stop(
        "this ordinal fit carries no test channel; refit with 'test' to ",
        "report out-of-sample category probabilities"
      )
    }
    object$yhat.test
  } else {
    object$yhat.train
  }
  probs <- reshapeChainedChannel(probs, n.chains, combineChains, 2L)
  if (type == "ev") {
    return(probs)
  }
  multinomialPpdFromProbs(probs)
}

# log P(y_i = k | eta, gamma) IS the reported category probability at the
# observed level: the run already stores the cumulative-probit difference,
# so no recomputation from eta/thresholds is needed. probs enters in the split
# (chains x) samples x obs x K layout; the result drops the K margin, the same
# shape type = "ppd" already returns for this family.
ordinalLogLik <- function(object, probs) {
  y <- object[["y"]]
  levels <- object[["levels"]]
  d <- dim(probs)
  K <- d[length(d)]
  nObs <- d[length(d) - 1L]
  n.draws <- length(probs) %/% (nObs * K)
  flat <- probs
  dim(flat) <- c(n.draws * nObs, K)
  k <- match(y, levels)
  idx <- rep(seq_len(nObs), each = n.draws)
  result <- log(flat[cbind(seq_len(n.draws * nObs), k[idx])])
  # a row the active-row mask takes out of the data set has no likelihood to
  # report, as pointwiseLogLikelihood reports it
  active <- object[["active"]]
  if (!is.null(active)) {
    result[active[idx] == 0] <- NaN
  }
  array(result, d[-length(d)], dimnames(probs)[-length(d)])
}

# The posterior-mean n x K probability matrix (colnames = levels), or
# (type = "class") the argmax category as an ordered factor over the original
# levels, or (type = "bart") the posterior-mean latent eta per observation.
# ci.level opts into a credible band instead of the posterior mean, taken on
# the full draws before any mean/class reduction.
fitted.bartOrdinal <- function(
  object,
  type = c("ev", "class", "bart"),
  sample = c("train", "test"),
  ci.level = NULL,
  ...
) {
  type <- validateType(type, eval(formals(fitted.bartOrdinal)$type))
  sample <- validateSample(sample, eval(formals(fitted.bartOrdinal)$sample))
  refuseUnusedGenericArgs(
    list(...),
    "fitted",
    "bartOrdinal",
    c(
      ordinalUnusedArgs,
      foreignArgsFor(fittedForeignReasons, names(formals(fitted.bartOrdinal)))
    )
  )
  refuseClassCiLevel(type, ci.level)
  # a training-side quantity pads back to the caller's own row count through
  # whatever the fit's na.action recorded; the test side never lost a row
  padded <- if (identical(sample, "train")) {
    function(value) padOmittedRows(object[["na.action"]], value)
  } else {
    identity
  }
  if (type == "bart") {
    latent <- if (sample == "test") object$latent.test else object$latent.train
    if (is.null(latent)) {
      stop(
        "this ordinal fit carries no test channel; refit with 'test' to ",
        "report out-of-sample latent fits"
      )
    }
    if (!is.null(ci.level)) {
      return(padded(posteriorInterval(latent, ci.level, trailing = 1L)))
    }
    return(padded(channelMeans(latent)))
  }
  probs <- extract.bartOrdinal(object, type = "ev", sample = sample)
  if (!is.null(ci.level)) {
    return(padded(posteriorInterval(probs, ci.level, trailing = 2L)))
  }
  meanProbs <- meanCategoryProbabilities(probs, object$levels)
  if (type == "ev") {
    return(padded(meanProbs))
  }
  padded(categoryFromMeanProbabilities(
    meanProbs,
    object$levels,
    ordered = TRUE
  ))
}

ordinalResidualsTypeReason <- list(
  type = paste0(
    "the residual is the per-category observed-indicator minus the fitted ",
    "probability"
  )
)

# y - fitted() on the response scale has no single scalar for a categorical
# response, so the per-category analog is the observed 1[y = k] indicator minus
# the fitted probability, an n x K matrix (columns named by 'levels').
residuals.bartOrdinal <- function(object, ...) {
  refuseResidualsSample(list(...))
  refuseUnusedGenericArgs(
    list(...),
    "residuals",
    "bartOrdinal",
    c(
      ordinalUnusedArgs,
      ordinalResidualsTypeReason,
      foreignArgsFor(
        residualsForeignReasons,
        names(formals(residuals.bartOrdinal))
      )
    )
  )
  phat <- fitted.bartOrdinal(object, type = "ev")
  y <- object$y
  indicator <- matrix(0, length(y), length(object$levels))
  indicator[cbind(seq_along(y), match(y, object$levels))] <- 1
  indicator <- padOmittedRows(object[["na.action"]], indicator)
  dimnames(indicator) <- dimnames(phat)
  indicator - phat
}

# Out-of-sample category probabilities by replaying the saved forest's trees to
# the newdata latent, then differencing the cumulative probit at the STORED
# per-draw thresholds. Requires a fit kept with
# keepTrees. The latent is f + o, as probit's: the fit's offset argument and
# offset() terms evaluated on newdata plus the 'offset' given here.
# type = "bart" returns the replayed latent eta; type = "ppd" draws
# one category per posterior draw. Only ppd touches the RNG, so type = "ev" is
# draw-neutral. The replay reads through $fit's own pointer: $fit is the
# sampler whose engine actually ran, so getPointer() can re-create it from
# stored state after a save/reload. The presence gate re-points to
# thresholds.raw, which this function already reads below and rides the same
# keepTrees gate.
predict.bartOrdinal <- function(
  object,
  newdata,
  type = c("ev", "ppd", "bart", "class"),
  offset = NULL,
  combineChains = TRUE,
  ci.level = NULL,
  na.action = dbarts::na.keepPredictors,
  n.threads = object$fit$control@n.threads,
  ...
) {
  type <- validateType(type, eval(formals(predict.bartOrdinal)$type))
  refuseUnusedGenericArgs(
    list(...),
    "predict",
    "bartOrdinal",
    c(
      ordinalUnusedArgs,
      predictOffsetUnusedArgs,
      foreignArgsFor(predictForeignReasons, names(formals(predict.bartOrdinal)))
    )
  )
  warnUnusedDots(list(...), "predict", "bartOrdinal")
  refuseNonNumericOffset(offset)
  refuseClassCiLevel(type, ci.level)
  if (is.null(object[["thresholds.raw"]])) {
    refuseWithoutTrees("predict")
  }
  # after the store check, whose absence the default here would otherwise
  # report as a missing slot
  n.threads <- validatePredictThreads(n.threads)
  offset <- predictTermOffset(object$fit$data, newdata, offset)
  rows <- preparePredictRows(
    newdata,
    object$fit$data@x,
    na.action,
    list(offset = offset)
  )
  if (isTRUE(rows$placeholder)) {
    restoreSeed <- protectRandomSeed()
    on.exit(restoreSeed(), add = TRUE)
  }
  offset <- subsetPredictInput(offset, rows, "offset")
  rowNames <- rows$keptNames
  n.chains <- object$n.chains
  # raw is n.new x n.samples (x n.chains): the replayed latent eta + o, the
  # test channel's shape
  raw <- predictCodedTest(object$fit, rows$x, offset, n.threads)
  if (type == "bart") {
    result <- nameObservationMargin(
      convertSamplesForCaller(raw, n.chains, combineChains),
      rowNames
    )
    if (!is.null(ci.level)) {
      return(padPredictedRows(
        posteriorInterval(result, ci.level, trailing = 1L),
        rows,
        first = TRUE
      ))
    }
    return(padPredictedRows(result, rows))
  }
  K <- object$K
  thresholds <- object$thresholds.raw # (K-1) x n.samples x n.chains
  if (length(dim(raw)) == 2L) {
    dim(raw) <- c(dim(raw), 1L)
  }
  n.new <- dim(raw)[1L]
  n.samples <- dim(raw)[2L]
  probs <- array(0, c(n.new, K, n.samples, n.chains))
  for (s in seq_len(n.samples)) {
    for (chain in seq_len(n.chains)) {
      probs[,, s, chain] <-
        ordinalCategoryProbabilities(raw[, s, chain], thresholds[, s, chain])
    }
  }
  if (n.chains == 1L) {
    probs <- array(probs, dim(probs)[1:3])
  }
  probs <- shapeMultinomialChannel(
    probs,
    object$levels,
    n.chains,
    combineChains,
    leadNames = rowNames
  )
  if (type == "ppd") {
    probs <- multinomialPpdFromProbs(probs)
  }
  if (!is.null(ci.level)) {
    trailing <- if (type %in% c("ev", "class")) 2L else 1L
    return(padPredictedRows(
      posteriorInterval(probs, ci.level, trailing = trailing),
      rows,
      first = TRUE
    ))
  }
  if (type == "class") {
    meanProbs <- meanCategoryProbabilities(probs, object$levels)
    return(padPredictedRows(
      categoryFromMeanProbabilities(meanProbs, object$levels, ordered = TRUE),
      rows
    ))
  }
  padPredictedRows(probs, rows, trailing = if (type == "ppd") 0L else 1L)
}

print.bartOrdinal <- function(x, ...) {
  printCall(x)
  cat("family: ordinal (cumulative probit)\n")
  cat("levels: ", paste(x$levels, collapse = " < "), "\n", sep = "")
  cat("n.chains: ", x$n.chains, "\n", sep = "")
  cat("n.trees: ", x$n.trees, "\n", sep = "")
  d <- dim(x$yhat.train)
  n.kept <- if (length(d) == 4L) d[2L] else d[1L] %/% x$n.chains
  cat("kept draws (per chain): ", n.kept, "\n", sep = "")
  if (!is.null(x$yhat.test)) {
    dt <- dim(x$yhat.test)
    n.test <- if (length(dt) == 4L) dt[3L] else dt[2L]
    cat("test rows: ", n.test, "\n", sep = "")
  }
  invisible(x)
}

# bart2(family = "nbinom") generics. The fit object is class "bartNegbin" -
# never "bart" - so the count arrays never fall through to the single-forest
# "bart" methods. A single forest fits the log mean eta = f(x) + c + o, so
# type = "bart" returns eta, while type = "ev" returns the mean counts
# mu = exp(eta) (the reported posterior mean count) and type = "ppd" draws one
# count per posterior draw from NB(size = r, mu). The per-draw shape r
# rides the fit's $shape field, the count analog of gaussian's sigma, and
# a drawn leaf scale rides $k, as on a bart fit.
# nbinom has a single forest, so 'forest'/'contribution' refuse for the same
# reason a bart-family single-forest fit does.
negbinUnusedArgs <- list(
  forest = singleForestReason,
  contribution = singleForestReason
)

extract.bartNegbin <- function(
  object,
  type = c(
    "ev",
    "ppd",
    "bart",
    "loglik",
    "shape",
    "sigma",
    "k",
    "leaf.prior.sd",
    "varcount"
  ),
  sample = c("train", "test"),
  combineChains = TRUE,
  ...
) {
  type <- validateType(type, eval(formals(extract.bartNegbin)$type))
  sampleSupplied <- !missing(sample)
  refuseUnusedGenericArgs(
    list(...),
    "extract",
    "bartNegbin",
    c(
      negbinUnusedArgs,
      foreignArgsFor(extractForeignReasons, names(formals(extract.bartNegbin)))
    )
  )
  n.chains <- fitNChains(object)

  if (type %in% c("shape", "sigma", "k", "leaf.prior.sd", "varcount")) {
    refuseSampleOnModelType(type, sampleSupplied)
    if (type == "varcount") {
      return(reshapeChainedChannel(
        object$varcount,
        n.chains,
        combineChains,
        1L
      ))
    }
    return(extractParameter(object, type, combineChains))
  }

  sample <- validateSample(sample, eval(formals(extract.bartNegbin)$sample))

  if (type == "loglik" && sample == "test") {
    stop("cannot extract a test sample log-likelihood; no test response exists")
  }

  latent <- if (sample == "test") object$latent.test else object$latent.train
  mu <- if (sample == "test") object$yhat.test else object$yhat.train
  if (sample == "test" && is.null(mu)) {
    stop(
      "this nbinom fit carries no test channel; refit with 'test' to report ",
      "out-of-sample counts"
    )
  }
  if (type == "bart") {
    return(combineOrUncombineChains(latent, n.chains, combineChains))
  }
  if (type == "loglik") {
    result <- negbinLogLik(
      object,
      combineOrUncombineChains(mu, n.chains, FALSE),
      n.chains
    )
    return(combineOrUncombineChains(result, n.chains, combineChains))
  }
  if (type == "ev") {
    return(combineOrUncombineChains(mu, n.chains, combineChains))
  }
  # type == "ppd": pair mu with shape in a common split layout so the two
  # align regardless of either's own storage, then reshape the result to the
  # caller's request
  muSplit <- combineOrUncombineChains(mu, n.chains, FALSE)
  disp <- scalarDrawVec(object$shape, n.chains, length(muSplit))
  result <- array(
    rnbinom(length(muSplit), size = disp, mu = as.vector(muSplit)),
    dim(muSplit),
    dimnames(muSplit)
  )
  combineOrUncombineChains(result, n.chains, combineChains)
}

# l[s,i] = dnbinom(y_i, size = shape[s], mu = yhat.train[s,i]); the
# per-draw shape pairs with the draws the same chain-fastest way the
# gaussian arm pairs sigma (shape is already sigma-shaped). mu enters
# forced to the split (chains x) samples x obs layout so it aligns with
# scalarDrawVec's own normalization regardless of either's own storage.
negbinLogLik <- function(object, mu, n.chains) {
  y <- object[["y"]]
  n.draws <- length(mu) %/% length(y)
  disp <- scalarDrawVec(object[["shape"]], n.chains, length(mu))
  result <- dnbinom(
    rep(y, each = n.draws),
    size = disp,
    mu = as.vector(mu),
    log = TRUE
  )
  # a row the active-row mask takes out of the data set has no likelihood to
  # report, as pointwiseLogLikelihood reports it
  active <- object[["active"]]
  if (!is.null(active)) {
    result[rep(active, each = n.draws) == 0] <- NaN
  }
  array(result, dim(mu), dimnames(mu))
}

# The posterior-mean count per observation (type = "ev"), the posterior-mean
# log mean per observation (type = "bart"), or a Monte Carlo mean over
# ppd draws (type = "ppd"). The observation margin is the array's last
# dimension in every chain layout, so we take the mean over that observation
# margin. ci.level opts into a credible band instead, taken on the full draws
# before the mean.
fitted.bartNegbin <- function(
  object,
  type = c("ev", "ppd", "bart"),
  sample = c("train", "test"),
  ci.level = NULL,
  ...
) {
  type <- validateType(type, eval(formals(fitted.bartNegbin)$type))
  sample <- validateSample(sample, eval(formals(fitted.bartNegbin)$sample))
  refuseUnusedGenericArgs(
    list(...),
    "fitted",
    "bartNegbin",
    c(
      negbinUnusedArgs,
      foreignArgsFor(fittedForeignReasons, names(formals(fitted.bartNegbin)))
    )
  )
  if (sample == "test" && is.null(object$yhat.test)) {
    stop(
      "this nbinom fit carries no test channel; refit with 'test' to report ",
      "out-of-sample counts"
    )
  }
  channel <- switch(
    type,
    bart = if (sample == "test") object$latent.test else object$latent.train,
    ev = if (sample == "test") object$yhat.test else object$yhat.train,
    # the ppd arm is a draw, not a stored channel; extract pairs each mu with
    # its own draw's shape, and the mean over the observation margin
    # below is invariant to the chain layout it returns
    ppd = extract.bartNegbin(object, type = "ppd", sample = sample)
  )
  # a training-side quantity pads back to the caller's own row count through
  # whatever the fit's na.action recorded; the test side never lost a row
  padded <- if (identical(sample, "train")) {
    function(value) padOmittedRows(object[["na.action"]], value)
  } else {
    identity
  }
  if (!is.null(ci.level)) {
    return(padded(posteriorInterval(channel, ci.level, trailing = 1L)))
  }
  padded(channelMeans(channel))
}

negbinResidualsTypeReason <- list(
  type = "the residual is the observed count minus the posterior-mean count"
)

# y - fitted() on the count scale: the observed count minus the posterior-mean
# count, an n-vector (the gaussian residual, on counts).
residuals.bartNegbin <- function(object, ...) {
  refuseResidualsSample(list(...))
  refuseUnusedGenericArgs(
    list(...),
    "residuals",
    "bartNegbin",
    c(
      negbinUnusedArgs,
      negbinResidualsTypeReason,
      foreignArgsFor(
        residualsForeignReasons,
        names(formals(residuals.bartNegbin))
      )
    )
  )
  padOmittedRows(object[["na.action"]], object$y) -
    fitted.bartNegbin(object, type = "ev")
}

# Out-of-sample mean counts by replaying the saved forest's trees to the newdata
# log mean eta, then mu = exp(eta). A log-exposure offset enters eta
# additively, the fit-time convention. Requires a fit kept with keepTrees.
# type = "bart" returns the replayed log mean; type = "ppd" draws one count per
# posterior draw at that draw's STORED shape r; type = "bart" and "ev"
# read no r. Only ppd touches the RNG, so type = "ev" is draw-neutral. The
# replay reads through $fit's own pointer: $fit is the sampler whose engine
# actually ran, so getPointer() can re-create it from stored state after a
# save/reload. The presence gate re-points to shape.raw, which is read
# below and rides the same keepTrees gate.
predict.bartNegbin <- function(
  object,
  newdata,
  type = c("ev", "ppd", "bart"),
  offset = NULL,
  combineChains = TRUE,
  ci.level = NULL,
  na.action = dbarts::na.keepPredictors,
  n.threads = object$fit$control@n.threads,
  ...
) {
  type <- validateType(type, eval(formals(predict.bartNegbin)$type))
  refuseUnusedGenericArgs(
    list(...),
    "predict",
    "bartNegbin",
    c(
      negbinUnusedArgs,
      predictOffsetUnusedArgs,
      foreignArgsFor(predictForeignReasons, names(formals(predict.bartNegbin)))
    )
  )
  warnUnusedDots(list(...), "predict", "bartNegbin")
  refuseNonNumericOffset(offset)
  if (is.null(object[["shape.raw"]])) {
    refuseWithoutTrees("predict")
  }
  # after the store check, whose absence the default here would otherwise
  # report as a missing slot
  n.threads <- validatePredictThreads(n.threads)
  offset <- predictTermOffset(object$fit$data, newdata, offset)
  # a missing offset row is incomplete the same way an unroutable predictor
  # is (dec-A89)
  rows <- preparePredictRows(
    newdata,
    object$fit$data@x,
    na.action,
    list(offset = offset)
  )
  if (isTRUE(rows$placeholder)) {
    restoreSeed <- protectRandomSeed()
    on.exit(restoreSeed(), add = TRUE)
  }
  rowNames <- rows$keptNames
  offset <- subsetPredictInput(offset, rows, "offset")
  n.chains <- object$n.chains
  # raw is n.new x n.samples (x n.chains): the replayed log mean eta
  raw <- predictCodedTest(object$fit, rows$x, offset, n.threads)
  if (type == "bart") {
    result <- nameObservationMargin(
      convertSamplesForCaller(raw, n.chains, combineChains),
      rowNames
    )
    if (!is.null(ci.level)) {
      return(padPredictedRows(
        posteriorInterval(result, ci.level, trailing = 1L),
        rows,
        first = TRUE
      ))
    }
    return(padPredictedRows(result, rows))
  }
  if (length(dim(raw)) == 2L) {
    dim(raw) <- c(dim(raw), 1L)
  }
  disp <- object$shape.raw # n.samples x n.chains
  n.new <- dim(raw)[1L]
  n.samples <- dim(raw)[2L]
  means <- array(0, c(n.new, n.samples, n.chains))
  for (s in seq_len(n.samples)) {
    for (chain in seq_len(n.chains)) {
      means[, s, chain] <- negbinMeanCounts(raw[, s, chain])
    }
  }
  if (n.chains == 1L) {
    means <- matrix(means, n.new, n.samples)
  }
  means <- convertSamplesForCaller(means, n.chains, combineChains)
  means <- nameObservationMargin(means, rowNames)
  if (type == "ppd") {
    # each count is drawn with its own draw's shape: the shapes are
    # laid out as the means are and take the caller's layout with them,
    # whichever layout the fit stored its own in
    shapes <- array(
      rep(disp, each = n.new),
      c(n.new, n.samples, n.chains)
    )
    if (n.chains == 1L) {
      shapes <- matrix(shapes, n.new, n.samples)
    }
    means <- negbinPpd(
      means,
      convertSamplesForCaller(shapes, n.chains, combineChains)
    )
  }
  if (!is.null(ci.level)) {
    return(padPredictedRows(
      posteriorInterval(means, ci.level, trailing = 1L),
      rows,
      first = TRUE
    ))
  }
  padPredictedRows(means, rows)
}

print.bartNegbin <- function(x, ...) {
  printCall(x)
  cat("family: negative binomial (log link)\n")
  if (is.null(x[["fixed"]][["shape"]])) {
    cat(
      "posterior mean shape (r): ",
      format(mean(x$shape), digits = 4L),
      "\n",
      sep = ""
    )
  } else {
    cat(
      "shape (r): fixed at ",
      format(x[["fixed"]][["shape"]], digits = 4L),
      "\n",
      sep = ""
    )
  }
  cat("n.chains: ", x$n.chains, "\n", sep = "")
  cat("n.trees: ", x$n.trees, "\n", sep = "")
  d <- dim(x$yhat.train)
  n.kept <- if (length(d) == 3L) d[2L] else d[1L] %/% x$n.chains
  cat("kept draws (per chain): ", n.kept, "\n", sep = "")
  if (!is.null(x$yhat.test)) {
    dt <- dim(x$yhat.test)
    n.test <- if (length(dt) == 3L) dt[3L] else dt[2L]
    cat("test rows: ", n.test, "\n", sep = "")
  }
  invisible(x)
}

# bart2(family = "hurdle.lognormal") generics. The fit object is class
# "bartHurdle" - never "bart" - holding the two
# conditionally-independent component fits ($zero, a probit fit of
# 1{y > 0} over all n; $positive, a gaussian fit of log(y) over the y > 0
# subset whose x.test is the full-n x). The report-time combine glues their
# posterior draws by sample index (any pairing is a valid joint draw, the parts
# share no parameters) and retransforms the positive part to the natural scale.
#
# type = "prob" is pi(x) = P(y > 0 | x), the zero part's own probability;
# type = "bart"/"link"/"log"
# the positive part's log-scale linear predictor f(x); type = "ev"/"response"
# the combined natural-scale mean via posterior-predictive Monte Carlo,
# E[y | x]_s = pi_s exp(f_s + sigma_s^2 / 2) PER DRAW s then aggregated across
# draws - NOT the biased plug-in of posterior means into one exponential.
# type = "ppd" is the proper bimodal
# predictive the plain gaussian ppd cannot make: per draw a Bernoulli(pi_s)
# spike at zero, else a lognormal exp(f_s + sigma_s z), z ~ N(0, 1).

# Fold the "response"/"link" type aliases onto the canonical "ev"/"bart" (the
# non-hurdle predict/extract/fitted idiom). Validation of the folded value
# against each method's allowed set stays at the call site, since those vary
# (some also reject length-0 input).
foldTypeAliases <- function(type) {
  if (is.character(type)) {
    if (type[1L] == "response") {
      type[1L] <- "ev"
    } else if (type[1L] == "link") {
      type[1L] <- "bart"
    }
  }
  type
}

# Fold the response/link aliases, then validate the requested type against the
# method's allowed set (its own 'type' formal, evaluated once at the call site
# and passed in) and return the canonical scalar. Centralizes the fold +
# %not_in% + stop shared by the bart predict/extract/fitted methods.
validateType <- function(type, allowed) {
  type <- foldTypeAliases(type)
  if (!is.character(type) || length(type) == 0L || type[1L] %not_in% allowed) {
    stop("type must be in '", paste0(allowed, collapse = "', '"), "'")
  }
  type[1L]
}

# Validate a 'sample' argument (train/test) against the method's own allowed
# set and return the canonical scalar - one wording for every class instead of
# a bare match.arg's "'arg' should be one of ...".
validateSample <- function(sample, allowed) {
  if (!is.character(sample) || sample[1L] %not_in% allowed) {
    stop("sample must be in '", paste0(allowed, collapse = "', '"), "'")
  }
  sample[1L]
}

# The own-class extract/fitted/predict/residuals/summary methods share the
# bart-family generics' NAMES but not their whole vocabulary: a K-widened or
# two-part shape has no single forest to select, no per-observation
# contribution to decompose, no separate test-sample fitted values, and (for
# a fixed-formula residual or vars channel) no caller choice to make. Reading
# a caller-supplied name out of '...' and stopping on the first hit refuses
# these by name instead of letting it fall through silently and discarding
# them.
refuseUnusedGenericArgs <- function(dots, generic, class, reasons) {
  supplied <- intersect(names(reasons), names(dots))
  if (length(supplied) > 0L) {
    stop(
      "'",
      supplied[1L],
      "' is not used by ",
      generic,
      " on a ",
      class,
      " fit: ",
      reasons[[supplied[1L]]]
    )
  }
  # A positional extra is the same caller mistake as a named one, and the more
  # likely one: the sibling method that does take the name takes it in a slot
  # this method's own formals do not reach, so the value lands in '...' and
  # would otherwise be discarded without a word.
  unnamed <- if (is.null(names(dots))) {
    seq_along(dots)
  } else {
    which(!nzchar(names(dots)))
  }
  if (length(unnamed) > 0L) {
    stop(
      generic,
      " on a ",
      class,
      " fit does not support unnamed arguments: ",
      length(unnamed),
      " supplied, the first at position ",
      unnamed[1L],
      " of '...'"
    )
  }
  invisible(NULL)
}

# A name that is a formal on one method of this surface and not on another
# is a caller mistake wherever it is foreign, not an argument to discard.
# Deriving each method's list from its own formals - rather than listing
# the foreign names by hand per class - is what keeps a name added to one
# signature refused on every sibling that does not take it.
foreignArgsFor <- function(reasons, own) {
  reasons[setdiff(names(reasons), own)]
}

# 'forest' selects among the per-forest channels only the "forest" arm
# reports, and among the forests' own k and leaf.prior.sd on a fit that has
# several; every other arm has already recombined them into the reported
# location, so a selection there would silently choose nothing. The model
# parameters and varcount are no recombined location, so they get their own
# wording, and so does a heteroscedastic fit's sigma, which is the variance
# forest's surface rather than a model parameter.
refuseForestSelectionOutsideForestArm <- function(
  type,
  forest,
  heteroscedastic = FALSE,
  numForests = 1L
) {
  if (is.null(forest)) {
    return(invisible(NULL))
  }
  if (type == "sigma" && heteroscedastic) {
    stop(
      "type = \"sigma\" on a heteroscedastic fit is the variance forest's ",
      "per-observation scale, not a per-forest quantity of the mean"
    )
  }
  if (type %in% c("k", "leaf.prior.sd") && numForests > 1L) {
    return(invisible(NULL))
  }
  if (type %in% c("sigma", "k", "leaf.prior.sd", "shape", "thresholds")) {
    stop(
      "type = \"",
      type,
      "\" is a model parameter, not a per-forest quantity"
    )
  }
  if (type == "varcount") {
    stop(
      "type = \"varcount\" keeps every forest on its trailing margin; ",
      "subset that margin"
    )
  }
  if (type != "forest") {
    stop(
      "type = \"",
      type,
      "\" does not support 'forest': every forest is ",
      "already recombined into the location it reports"
    )
  }
  invisible(NULL)
}

# type = "class" is a label, not a quantity with a credible band; the class
# reduction below is otherwise unreachable whenever ci.level is supplied (the
# ci.level branch returns first), so the combination is refused by name
# instead of the ev band being silently returned in the class request's
# place.
refuseClassCiLevel <- function(type, ci.level) {
  if (type == "class" && !is.null(ci.level)) {
    stop(
      "type = \"class\" does not support 'ci.level': a class prediction is ",
      "a label rather than a quantity with a credible band"
    )
  }
  invisible(NULL)
}

# Every own-class family has its own list of names its K-widened or two-part
# shape has no room for; bart did not, so the two names a
# fit-reduction never selects among had nowhere to be refused.
bartUnusedArgs <- list(
  forest = paste0(
    "the reduction is over the combined location, in which every forest ",
    "is already included"
  ),
  contribution = paste0(
    "the per-observation contribution decomposes one forest's fit, and the ",
    "reduction here is over the combined location"
  )
)

# The derived reason tables the surface's predict/extract/fitted/residuals
# methods compose via foreignArgsFor above, one entry per name that is a
# formal on some method of the generic and foreign on another. Composed
# AFTER a method's own class list (multinomialUnusedArgs and siblings,
# bartUnusedArgs), so a class-specific reason for the same
# name still wins - refuseUnusedGenericArgs reports the first hit in
# 'reasons', and composition order is priority order.
predictForeignReasons <- list(
  sample = "the fit's stored train and test channels are extract's 'sample'",
  weights = "this family's posterior-predictive draw takes no per-observation weight",
  bases = "only an amplitude-coupled multi-forest fit takes 'bases' at the predicted rows",
  contribution = "the per-observation contribution decomposition belongs to extract(type = \"forest\")",
  value = "predict's channel argument is named 'type'"
)

# Every predict method's own dots are otherwise inert (handed nowhere but to
# refuseUnusedGenericArgs above), so a name that survives the blocklist -
# neither a foreign name from another method of this surface nor a
# class-specific one - would otherwise be silently discarded rather than
# pointing the caller at a typo. warnUnusedDots warns on it (never errors,
# so a subclass method forwarding its own extra formals through NextMethod's
# '...' is not refused for arguments its caller legitimately supplied). It is
# package-local rather than base's chkDots, which quotes with the locale's
# fancy quotes before R 4.6; this diagnosis must be matchable by message on
# every R the package supports.

extractReplaysNothingReason <- "extract reads stored channels and replays nothing"
extractForeignReasons <- list(
  ci.level = "extract returns the draws that fitted() and predict() take a band over",
  newdata = "predict(object, newdata) is the read at new rows",
  offset = extractReplaysNothingReason,
  weights = extractReplaysNothingReason,
  n.threads = extractReplaysNothingReason,
  bases = extractReplaysNothingReason
)

fittedSummarizesNothingReason <- "fitted summarizes stored channels and replays nothing"
fittedForeignReasons <- list(
  combineChains = "the per-chain draws are extract(object, combineChains = FALSE)",
  sample = "this fit carries no test channel; call predict on newdata",
  newdata = fittedSummarizesNothingReason,
  offset = fittedSummarizesNothingReason,
  weights = fittedSummarizesNothingReason,
  n.threads = fittedSummarizesNothingReason,
  bases = fittedSummarizesNothingReason
)

residualsSummarizeNothingReason <- "residuals summarize stored channels and replay nothing"
residualsForeignReasons <- list(
  ci.level = "residuals are the observed response minus the posterior-mean fit",
  combineChains = "the per-chain draws are extract(object, combineChains = FALSE)",
  newdata = residualsSummarizeNothingReason,
  offset = residualsSummarizeNothingReason,
  weights = residualsSummarizeNothingReason,
  n.threads = residualsSummarizeNothingReason,
  bases = residualsSummarizeNothingReason
)

survivalProbabilitiesDrawsReason <- "survivalProbabilities returns the draws of S(t | x) at 'times'"
survivalProbabilitiesOwnArgsReason <- "survivalProbabilities takes 'times', 'newdata' and 'offset' alone"
survivalProbabilitiesForeignReasons <- list(
  type = survivalProbabilitiesDrawsReason,
  sample = survivalProbabilitiesDrawsReason,
  ci.level = survivalProbabilitiesDrawsReason,
  weights = survivalProbabilitiesOwnArgsReason,
  n.threads = survivalProbabilitiesOwnArgsReason,
  forest = survivalProbabilitiesOwnArgsReason,
  contribution = survivalProbabilitiesOwnArgsReason,
  bases = survivalProbabilitiesOwnArgsReason
)

# Resolve a hurdle type argument: fold the "response"/"link"/"log" aliases onto
# the canonical "ev"/"bart" and validate against 'allowed' (the predict.bart
# idiom, so a mis-typed request errors rather than silently mis-reporting).
resolveHurdleType <- function(type, allowed) {
  if (is.character(type) && length(type) > 0L) {
    type <- foldTypeAliases(type)
    if (type[1L] == "log") {
      type[1L] <- "bart"
    }
  }
  if (!is.character(type) || length(type) == 0L || type[1L] %not_in% allowed) {
    stop("type must be in '", paste0(allowed, collapse = "', '"), "'")
  }
  type[1L]
}

hurdleNChains <- function(object) {
  zero <- object$zero
  if (!is.null(zero[["fit"]])) {
    zero$fit$control@n.chains
  } else {
    zero$n.chains
  }
}

# A scalar-per-draw field (sigma, shape, ...) as a flat vector aligned,
# draw for draw, with the fit draws' as.vector order (chain-fastest, then
# sample, then observation - the layout pointwiseLogLikelihood and
# sampleFromPPD pair sigma with fits in); the field may be stored combined
# (flat, chain-major) or split ((n.chains x) n.samples matrix, chain-fastest),
# so it is normalized to the split matrix first, exactly as chainFastest does
# in pointwiseLogLikelihood, THEN recycled across the n.obs draw-blocks with
# rep_len so it aligns with a channel of any shape (combined or split) whose
# trailing margin is the observations - regardless of the field's own
# storage, or of a caller-requested combineChains that differs from it.
scalarDrawVec <- function(x, n.chains, n.total) {
  if (is.null(dim(x))) {
    x <- uncombineChains(as.vector(x), n.chains)
  }
  rep_len(as.vector(x), n.total)
}

# Glue the flat, draw-aligned zero-part-probability, positive-log-mean, and
# positive-sigma vectors into the requested channel and reshape to the fit's
# uncombined draw layout ('shape'). Only "ppd" touches the RNG (Bernoulli then
# lognormal), so the default "ev" is draw-neutral.
combineHurdleChannel <- function(
  type,
  piVec,
  fVec,
  sigmaVec,
  shape,
  shapeNames = NULL
) {
  channel <- switch(
    type,
    prob = piVec,
    bart = fVec,
    ev = piVec * exp(fVec + 0.5 * sigmaVec^2),
    ppd = rbinom(length(piVec), 1L, piVec) *
      exp(fVec + sigmaVec * rnorm(length(fVec)))
  )
  array(channel, shape, shapeNames)
}

# A single-forest fit's draws at coded rows, uncombined, as
# predict(type = "ev") and predict(type = "bart") report them: a hurdle
# component's, or a discrete-time hazard fit's per-period hazards.
codedRowDraws <- function(
  component,
  x,
  type,
  n.threads,
  rowNames,
  offset = NULL
) {
  raw <- predictCodedTest(component$fit, x, offset, n.threads)
  if (is.list(raw)) {
    raw <- raw$mean
  }
  result <- convertSamplesFromDbartsToBart(
    raw,
    component$fit$control@n.chains,
    FALSE
  )
  if (type == "ev") {
    result <- probabilityFromLatents(result, component)
  }
  nameObservationMargin(result, rowNames)
}

# The zero part's pi(x), positive log-mean f(x), and positive per-observation
# sigma draws for the combine, each a flat vector in the fit's uncombined
# as.vector order, plus the uncombined 'shape' to fold back to. In-sample reads
# the stored channels - the zero fit's ev over all n, and the positive
# fit's log-scale (bart) fits at the FULL-n rows through its x.test channel (the
# zero rows it never trained on included); out-of-sample replays both saved
# forests at the rows preparePredictRows kept, coded once for each component.
hurdleParts <- function(object, rows = NULL, n.threads = 1L) {
  if (is.null(rows)) {
    pi <- extract(
      object$zero,
      type = "ev",
      sample = "train",
      combineChains = FALSE
    )
    f <- extract(
      object$positive,
      type = "bart",
      sample = "test",
      combineChains = FALSE
    )
  } else {
    pi <- codedRowDraws(
      object$zero,
      rows$x,
      "ev",
      n.threads,
      rows$keptNames
    )
    # the positive part codes the same rows against its own training design;
    # the caller has already been warned about anything the coding says
    positiveTrain <- object$positive$fit$data@x
    positiveX <- suppressTestMatchWarnings(validateXTest(
      if (isTRUE(rows$placeholder)) {
        positiveTrain[1L, , drop = FALSE]
      } else {
        rows$newdata
      },
      positiveTrain,
      refuseMissing = FALSE
    ))
    f <- codedRowDraws(
      object$positive,
      positiveX,
      "bart",
      n.threads,
      rows$keptNames
    )
  }
  sigmaVec <- scalarDrawVec(
    object$positive$sigma,
    hurdleNChains(object),
    length(f)
  )
  list(
    pi = as.vector(pi),
    f = as.vector(f),
    sigma = sigmaVec,
    shape = dim(f),
    # the positive part's names: its test rows are the full design
    names = dimnames(f)
  )
}

finishHurdle <- function(parts, type, n.chains, combineChains, ci.level) {
  channel <- combineHurdleChannel(
    type,
    parts$pi,
    parts$f,
    parts$sigma,
    parts$shape,
    parts$names
  )
  result <- combineOrUncombineChains(channel, n.chains, combineChains)
  if (!is.null(ci.level)) {
    return(posteriorInterval(result, ci.level))
  }
  result
}

# hurdle's two components are each a single forest, so 'forest'/'contribution'
# refuse for the same reason a bart-family single-forest fit does.
hurdleUnusedArgs <- list(
  forest = singleForestReason,
  contribution = singleForestReason
)

extract.bartHurdle <- function(
  object,
  type = c(
    "ev",
    "ppd",
    "prob",
    "bart",
    "loglik",
    "sigma",
    "k",
    "leaf.prior.sd",
    "varcount"
  ),
  sample = c("train", "test"),
  combineChains = TRUE,
  ...
) {
  type <- resolveHurdleType(type, eval(formals(extract.bartHurdle)$type))
  sampleSupplied <- !missing(sample)
  refuseUnusedGenericArgs(
    list(...),
    "extract",
    "bartHurdle",
    c(
      hurdleUnusedArgs,
      foreignArgsFor(extractForeignReasons, names(formals(extract.bartHurdle)))
    )
  )

  # sigma is positive$sigma, the only one the composition carries; k,
  # leaf.prior.sd and varcount are lists keyed zero/positive, each part's own
  if (type %in% c("sigma", "k", "leaf.prior.sd", "varcount")) {
    refuseSampleOnModelType(type, sampleSupplied)
    n.chains <- hurdleNChains(object)
    if (type == "sigma") {
      return(extractParameter(object$positive, "sigma", combineChains))
    }
    parts <- object[c("zero", "positive")]
    if (type == "varcount") {
      return(lapply(parts, function(part) {
        reshapeChainedChannel(part$varcount, n.chains, combineChains, 1L)
      }))
    }
    return(lapply(parts, extractParameter, type, combineChains))
  }

  sample <- validateSample(sample, eval(formals(extract.bartHurdle)$sample))
  if (sample == "test") {
    stop(
      "this hurdle fit carries no separate test channel; call predict on ",
      "newdata for out-of-sample combined draws"
    )
  }
  n.chains <- hurdleNChains(object)
  if (type == "loglik") {
    return(combineOrUncombineChains(
      hurdleLogLik(object),
      n.chains,
      combineChains
    ))
  }
  finishHurdle(
    hurdleParts(object),
    type,
    n.chains,
    combineChains,
    NULL
  )
}

# Reuses hurdleParts() verbatim: the pi/f/sigma draws the ev/ppd channels
# already glue, flat and draw-aligned, at ALL n rows (the zero part's own
# channel; the positive part's x.test channel, zero rows included). y == 0
# rows take the zero part's own log(1 - pi); y > 0 rows take the zero part's
# log(pi) plus the lognormal density of y on its NATURAL scale (a -log(y)
# Jacobian against the stored log-scale channel) - comparable to any other
# model of y, not of log(y); NO truncation, since the positive part's
# lognormal support is already (0, Inf) - a future truncated (count) hurdle
# would need one and must not reuse this formula unchanged. This is NOT the
# sum of the two components' own loglik channels: the positive fit's channel
# covers only its y > 0 rows, sits on the log scale, and carries no Jacobian.
hurdleLogLik <- function(object) {
  parts <- hurdleParts(object)
  y <- object[["y"]]
  n.draws <- length(parts$f) %/% length(y)
  yRep <- rep(y, each = n.draws)
  positive <- yRep > 0
  result <- numeric(length(parts$f))
  result[!positive] <- log1p(-parts$pi[!positive])
  result[positive] <- log(parts$pi[positive]) +
    dnorm(
      log(yRep[positive]),
      parts$f[positive],
      parts$sigma[positive],
      log = TRUE
    ) -
    log(yRep[positive])
  array(result, parts$shape, parts$names)
}

fitted.bartHurdle <- function(
  object,
  type = c("ev", "ppd", "prob", "bart"),
  ci.level = NULL,
  ...
) {
  type <- resolveHurdleType(type, eval(formals(fitted.bartHurdle)$type))
  refuseUnusedGenericArgs(
    list(...),
    "fitted",
    "bartHurdle",
    c(
      hurdleUnusedArgs,
      foreignArgsFor(fittedForeignReasons, names(formals(fitted.bartHurdle)))
    )
  )
  # a hurdle fit has no separate test channel (extract.bartHurdle refuses
  # sample = "test" unconditionally), so the read is always the training
  # rows, which pad back to the caller's own row count through whatever
  # the fit's na.action recorded
  draws <- extract(object, type = type, sample = "train", combineChains = TRUE)
  if (!is.null(ci.level)) {
    return(padOmittedRows(
      object[["na.action"]],
      posteriorInterval(draws, ci.level)
    ))
  }
  padOmittedRows(object[["na.action"]], channelMeans(draws))
}

residuals.bartHurdle <- function(object, type = "ev", ...) {
  # natural-scale residual against the stored original response over all n;
  # call the method by name (the residuals.bart idiom) so the package namespace
  # need not import the stats fitted generic
  refuseResidualsSample(list(...))
  refuseUnusedGenericArgs(
    list(...),
    "residuals",
    "bartHurdle",
    c(
      hurdleUnusedArgs,
      foreignArgsFor(
        residualsForeignReasons,
        names(formals(residuals.bartHurdle))
      )
    )
  )
  padOmittedRows(object[["na.action"]], object$y) -
    fitted.bartHurdle(object, type = type)
}

# Out-of-sample combined draws by replaying BOTH saved forests at newdata and
# gluing them the same way the in-sample channels are. Requires a fit kept
# with keepTrees (both components keep trees
# when the hurdle does).
predict.bartHurdle <- function(
  object,
  newdata,
  type = c("ev", "ppd", "prob", "bart"),
  offset = NULL,
  combineChains = TRUE,
  ci.level = NULL,
  na.action = dbarts::na.keepPredictors,
  n.threads = object$zero$fit$control@n.threads,
  ...
) {
  type <- resolveHurdleType(type, eval(formals(predict.bartHurdle)$type))
  refuseUnusedGenericArgs(
    list(...),
    "predict",
    "bartHurdle",
    c(
      hurdleUnusedArgs,
      predictNoOffsetUnusedArgs,
      foreignArgsFor(predictForeignReasons, names(formals(predict.bartHurdle)))
    )
  )
  warnUnusedDots(list(...), "predict", "bartHurdle")
  refusePredictOffsetChannel(offset, "bartHurdle")
  if (is.null(object$zero[["fit"]])) {
    refuseWithoutTrees("predict")
  }
  # after the zero fit check, whose absence the default here would
  # otherwise report as a missing slot
  n.threads <- validatePredictThreads(n.threads)
  n.chains <- hurdleNChains(object)
  # both parts share one routable set (refuseHurdlePositiveMissingness), so
  # the rows are resolved once, against the zero part's design
  rows <- preparePredictRows(newdata, object$zero$fit$data@x, na.action)
  if (isTRUE(rows$placeholder)) {
    restoreSeed <- protectRandomSeed()
    on.exit(restoreSeed(), add = TRUE)
  }
  padPredictedRows(
    finishHurdle(
      hurdleParts(object, rows, n.threads),
      type,
      n.chains,
      combineChains,
      ci.level
    ),
    rows,
    first = !is.null(ci.level)
  )
}

print.bartHurdle <- function(x, ...) {
  printCall(x)
  cat("family: hurdle.lognormal (probit zero part + lognormal positive part)\n")
  cat("zero n (all rows): ", length(x$zero$y), "\n", sep = "")
  cat("positive-part n (y > 0): ", length(x$positive$y), "\n", sep = "")
  invisible(x)
}

# this method's '...' forwards nowhere and nothing delegates through it, so a
# caller-supplied argument of any kind - named or positional - is refused
# rather than the generic's own "on a <class> fit" wording, which reads wrong
# for a sampler
refuseSamplerExtractArgs <- function(dots) {
  if (length(dots) > 0L) {
    named <- names(dots)
    named <- if (is.null(named)) character() else named[nzchar(named)]
    stop(
      if (length(named) > 0L) {
        paste0("'", named[1L], "'")
      } else {
        "a positional argument"
      },
      " is not used by extract on a dbartsSampler: this method returns the ",
      "sampler's coded predictor matrix"
    )
  }
  invisible(NULL)
}

# Materialize the sampler's predictor code matrix: factor columns as their
# integer codes, the form data@x holds and the matrix getTrees replays. A
# dense-frame/mixed container materializes through as.matrix; a
# plain matrix (or a sparse dgCMatrix held as such) is returned unchanged.
extract.dbartsSampler <- function(object, type = "predictors", ...) {
  refuseSamplerExtractArgs(list(...))
  if (!is.character(type) || length(type) == 0L || type[1L] != "predictors") {
    stop("'type' must be one of 'predictors'")
  }
  x <- object$data@x
  if (inherits(x, "dbartsMixedMatrix")) as.matrix(x) else x
}

# fit-level dispatch for the sampler's plotTree method, so a kept bart fit
# can be plotted directly instead of reaching into $fit; chainNum
# and sampleNum forward only when supplied, since the method detects them by
# their absence
plotTree.dbartsSampler <- function(object, ...) {
  refusePlotTreeArgs(sys.call())
  invisible(object$plotTree(...))
}

# do.call(object$fit$plotTree, args) below forwards whatever the caller wrote
# by name straight through, so a caller typing the extract/fitted vocabulary's
# 'sample'/'chain' - instead of this method's own 'sampleNum'/'chainNum' -
# partial-matches the wrong formal via R's own argument matching and silently
# draws a different tree than intended. Reading the RAW (unmatched) call
# catches the exact name the caller wrote, before that matching resolves it.
refusePlotTreeArgs <- function(rawCall) {
  supplied <- intersect(c("sample", "chain"), names(rawCall))
  if (length(supplied) > 0L) {
    stop(
      "'",
      supplied[1L],
      "' is not used by plotTree; the saved ",
      supplied[1L],
      " is '",
      supplied[1L],
      "Num'"
    )
  }
  invisible(NULL)
}

plotTree.bart <- function(
  object,
  treeNum = 1L,
  chainNum,
  sampleNum,
  forest = NULL,
  ...
) {
  refusePlotTreeArgs(sys.call())
  if (is.null(object[["fit"]])) {
    refuseWithoutTrees("plotTree", bartKeepTreesArgument(object))
  }
  args <- list(treeNum = treeNum, forest = forest, ...)
  if (!missing(chainNum)) {
    args$chainNum <- chainNum
  }
  if (!missing(sampleNum)) {
    args$sampleNum <- sampleNum
  }
  invisible(do.call(object$fit$plotTree, args))
}

# plotTree.bart reads the trees off object$fit; a K-widened or two-part
# own-class fit has no single sampler that reads that way (bartHurdle has
# two), so each refuses by name instead of falling through to "no applicable
# method", pointing at the sampler(s) that do carry the trees.
refusePlotTreeMethod <- function(class, hint) {
  stop(
    "plotTree is defined for bart and dbartsSampler fits; a ",
    class,
    " fit's trees live on its sampler - call ",
    hint
  )
}
plotTree.bartMultinomial <- function(object, ...) {
  refusePlotTreeMethod("bartMultinomial", "plotTree(object$fit, ...)")
}
plotTree.bartOrdinal <- function(object, ...) {
  refusePlotTreeMethod("bartOrdinal", "plotTree(object$fit, ...)")
}
plotTree.bartNegbin <- function(object, ...) {
  refusePlotTreeMethod("bartNegbin", "plotTree(object$fit, ...)")
}
plotTree.bartHurdle <- function(object, ...) {
  refusePlotTreeMethod(
    "bartHurdle",
    "plotTree(object$zero$fit, ...) or plotTree(object$positive$fit, ...)"
  )
}

# survivalProbabilities.bart dispatches on an aft or discrete-time hazard
# fit; none of the four own-class families is either, so each refuses by name
# instead of falling through to "no applicable method".
refuseSurvivalProbabilitiesMethod <- function(class) {
  stop(
    "survivalProbabilities applies to a discrete-time hazard fit ",
    "(bart(family = \"hazard\")); a ",
    class,
    " fit has no hazard channel"
  )
}
survivalProbabilities.bartMultinomial <- function(object, ...) {
  refuseSurvivalProbabilitiesMethod("bartMultinomial")
}
survivalProbabilities.bartOrdinal <- function(object, ...) {
  refuseSurvivalProbabilitiesMethod("bartOrdinal")
}
survivalProbabilities.bartNegbin <- function(object, ...) {
  refuseSurvivalProbabilitiesMethod("bartNegbin")
}
survivalProbabilities.bartHurdle <- function(object, ...) {
  refuseSurvivalProbabilitiesMethod("bartHurdle")
}

# The gaussian posterior predictive's noise scale, in the split layout's
# chain-fastest order sampleFromPPD draws in: the per-draw scalar sigma
# recycled across the observation margin, or a heteroscedastic fit's own
# per-observation s(x), which is already on the response scale and so stands
# in for sigma rather than scaling it. A case weight is a precision multiplier
# on whichever of the two it is, giving sd_i = scale_i / sqrt(w_i).
ppdNoiseScale <- function(sigma, s, weights, n.obs, n.draws) {
  sd <- if (is.null(s)) {
    rep_len(as.vector(sigma), n.obs * n.draws)
  } else if (length(s) != n.obs * n.draws) {
    stop("the fit's 's(x)' draws do not match its predicted draws")
  } else {
    as.vector(s)
  }
  if (!is.null(weights)) {
    sd <- sd * rep(sqrt(1 / weights), each = n.draws)
  }
  sd
}

# The posterior predictive noise at scale 'sd', laid out as ppdNoiseScale lays
# it out: gaussian, or, given a student() fit's per-draw degrees of freedom
# (chain-fastest, as sigma), t with each draw's own, sd * t_nu.
ppdNoise <- function(n, sd, df = NULL) {
  if (is.null(df)) {
    return(rnorm(n, 0, sd))
  }
  sd * stats::rt(n, rep_len(as.vector(df), n))
}

# the number of draws the noise scale spans: one per sigma draw, or, on a
# heteroscedastic fit, which carries no sigma, one per row of s(x)'s draws
ppdNumDraws <- function(sigma, s, n.obs) {
  if (is.null(s)) length(sigma) else length(s) %/% n.obs
}

# ev (expected value) enters in the caller's requested layout: chains split
# ((n.chains x) n.samples x n.obs, obs last) or chains combined ((n.chains *
# n.samples) x n.obs, chain-blocked rows - all of chain 1's samples, then
# chain 2's). Every family draws in the split layout's chain-fastest order,
# then reshapes to the caller's shape with the same combineChains() helper
# the stored draws go through, so a combined and a split ppd draw from the
# same seed agree bit-for-bit after accounting for row order. Gaussian draws
# noise in sigma's chain-fastest order, normalized below from whichever of
# its two storage layouts the fit used (combined: flat, chain-major; split:
# (n.chains x) n.samples matrix, already chain-fastest) - and adds it
# (reshaped when ev is combined); binary draws rbinom against the
# split-order probabilities and reshapes the outcome, since the draw
# depends on ev and cannot be reshaped after the fact. Single chain and
# already-split ev take the flat path unchanged. n.chains is needed only to
# perform that reshape. s carries a heteroscedastic fit's per-observation
# residual scale in that same split layout (heteroscedasticScale); it is NULL
# for a homoscedastic fit, whose scale is the per-draw scalar sigma.
sampleFromPPD <- function(ev, object, weights, n.chains = 1L, s = NULL) {
  oldSeed <- NULL
  if (!is.null(object[["seed"]])) {
    oldSeed <- .GlobalEnv$.Random.seed
    .GlobalEnv$.Random.seed <- object$seed
  }

  responseIsBinary <- fitIsBinary(object)
  sigma <- object$sigma
  if (!responseIsBinary && !is.null(sigma) && is.null(dim(sigma))) {
    sigma <- uncombineChains(as.vector(sigma), n.chains)
  }

  # a student() fit's noise is t with each draw's own degrees of freedom,
  # one scalar per draw as sigma is and paired with it the same way, as the
  # pointwise log-likelihood pairs them
  df <- NULL
  if (fitIsStudent(object)) {
    df <- object[["resid.df"]]
    if (is.null(df)) {
      stop(
        "posterior predictive sampling needs the fit's per-draw residual ",
        "degrees of freedom, which it does not store"
      )
    }
    if (is.null(dim(df))) {
      df <- uncombineChains(as.vector(df), n.chains)
    }
  }

  if (is.null(weights)) {
    if (responseIsBinary) {
      if (n.chains > 1L && length(dim(ev)) < 3L) {
        # ev is combined (chain-blocked rows). Draw in the split layout's
        # chain-fastest order and reshape with combineChains, so a combined
        # and a split draw from the same seed agree bit-for-bit (the gaussian
        # branch's guarantee), instead of consuming the RNG stream in the
        # combined layout's differing order.
        ev.split <- uncombineChains(ev, n.chains)
        draws <- rbinom(length(ev), 1L, as.vector(ev.split))
        result <- combineChains(array(draws, dim(ev.split)))
        dimnames(result) <- dimnames(ev)
      } else if (length(dim(ev)) > 2L) {
        result <- array(
          rbinom(length(ev), 1L, ev),
          dim(ev),
          dimnames = dimnames(ev)
        )
      } else {
        result <- matrix(
          rbinom(length(ev), 1L, ev),
          nrow(ev),
          ncol(ev),
          dimnames = list(rownames(ev), colnames(ev))
        )
      }
    } else {
      n.obs <- dim(ev)[length(dim(ev))]
      n.draws <- ppdNumDraws(sigma, s, n.obs)
      noise <- ppdNoise(
        n.obs * n.draws,
        ppdNoiseScale(sigma, s, NULL, n.obs, n.draws),
        df
      )
      if (n.chains > 1L && length(dim(ev)) < 3L) {
        noise <- combineChains(array(
          noise,
          c(n.chains, n.draws %/% n.chains, n.obs)
        ))
      }
      result <- ev + noise
    }
  } else {
    if (responseIsBinary) {
      # a weight-w row is w iid bernoulli trials; the coherent posterior
      # predictive draw is the number of successes, rbinom(, w, ev), not a
      # bernoulli draw scaled by w. size is recycled to match ev's own
      # column-major fill so each obs's weight lines up with its draws.
      if (n.chains > 1L && length(dim(ev)) < 3L) {
        # combined ev: draw in the split layout's chain-fastest order and
        # reshape, matching the unweighted binary and gaussian branches so a
        # combined and a split draw from the same seed agree bit-for-bit.
        ev.split <- uncombineChains(ev, n.chains)
        size <- rep(weights, each = prod(dim(ev.split)[1L:2L]))
        draws <- rbinom(length(ev), size, as.vector(ev.split))
        result <- combineChains(array(draws, dim(ev.split)))
        dimnames(result) <- dimnames(ev)
      } else if (length(dim(ev)) > 2L) {
        size <- rep(weights, each = prod(dim(ev)[1L:2L]))
        result <- array(
          rbinom(length(ev), size, ev),
          dim(ev),
          dimnames = dimnames(ev)
        )
      } else {
        size <- rep(weights, each = nrow(ev))
        result <- matrix(
          rbinom(length(ev), size, ev),
          nrow(ev),
          ncol(ev),
          dimnames = list(rownames(ev), colnames(ev))
        )
      }
    } else {
      n.obs <- dim(ev)[length(dim(ev))]
      n.draws <- ppdNumDraws(sigma, s, n.obs)
      sd <- ppdNoiseScale(sigma, s, weights, n.obs, n.draws)
      noise <- ppdNoise(n.obs * n.draws, sd, df)
      if (n.chains > 1L && length(dim(ev)) < 3L) {
        noise <- combineChains(array(
          noise,
          c(n.chains, n.draws %/% n.chains, n.obs)
        ))
      }
      result <- ev + noise
    }
  }
  if (!is.null(oldSeed)) {
    .GlobalEnv$.Random.seed <- oldSeed
  }

  result
}

# family/chain-count/tree-count/burn-in/kept-draws synopsis for print.bart,
# built only from fields that exist regardless of keepCall and keepSampler -
# so a fit created with keepCall = FALSE still prints something useful. n.trees and
# n.burn are only recoverable when the sampler itself was kept (keepTrees/
# keepSampler = TRUE); they are omitted otherwise, since the fit object
# does not retain them on its own.
fitSynopsis <- function(x) {
  fit <- x[["fit"]]

  n.chains <- if (!is.null(fit)) fit$control@n.chains else x$n.chains
  control <- if (!is.null(fit)) fit$control else NULL

  varcountDims <- dim(x[["varcount"]])
  # a multi-forest fit's varcount carries a trailing forest margin (the
  # shapeMultinomialChannel shape: draws x p x n.forests), so its rank is one
  # higher throughout and the single-forest arms below would read the predictor
  # count as the draw count. n.forests, not the rank, is what separates the two
  # - a single-forest uncombined varcount is rank 3 as well.
  n.forests <- fitNumForests(x)
  n.kept <- if (!is.null(control)) {
    control@n.samples
  } else if (is.null(varcountDims)) {
    NA_integer_
  } else if (n.forests > 1L) {
    if (length(varcountDims) == 4L) {
      varcountDims[2L]
    } else {
      varcountDims[1L] %/% n.chains
    }
  } else if (length(varcountDims) == 3L) {
    varcountDims[2L]
  } else if (n.chains > 1L) {
    varcountDims[1L] %/% n.chains
  } else {
    varcountDims[1L]
  }

  cat("family: ", x$family, "\n", sep = "")
  cat("n.chains: ", n.chains, "\n", sep = "")
  if (!is.null(control)) {
    cat("n.trees: ", control@n.trees, "\n", sep = "")
    cat("n.burn: ", control@n.burn, "\n", sep = "")
  }
  if (!is.na(n.kept)) {
    cat("kept draws (per chain): ", n.kept, "\n", sep = "")
  }
  if (!is.null(x[["monotone.prior"]])) {
    cat("monotone prior: ", x[["monotone.prior"]], "\n", sep = "")
  }
  invisible(NULL)
}

print.bart <- function(x, ...) {
  printCall(x)
  fitSynopsis(x)
  invisible(x)
}
