# predict, extract, fitted, and residuals methods for bart, rbart, and the
# multinomial, ordinal, nbinom, and hurdle fit objects

extract <- function(object, ...) UseMethod("extract")

plotTree <- function(object, ...) UseMethod("plotTree")

survivalProbabilities <- function(object, ...) {
  UseMethod("survivalProbabilities")
}

# latent-scale draws to probabilities for a binary fit; fits saved before
# the family element existed are all probit
probabilityFromLatents <- function(latents, object) {
  if (identical(object[["family"]], "logistic")) {
    plogis(latents)
  } else {
    pnorm(latents)
  }
}

# A heteroscedastic fit's residual scale is the per-observation surface s(x),
# stored - like the draw channels it is laid
# out as - either combined or split. This normalizes whichever storage the fit
# used to the split, chain-fastest layout, so as.vector() on it enumerates
# draws in the order as.vector() on the split fits does and the two pair
# element for element. s(x) is already on the response scale (the working
# surface times the response range), so it REPLACES the scalar sigma rather
# than scaling it: under this parameterization sigma is a fixed unit residual
# times that same range, a constant carrying no posterior content. NULL passes
# through, marking a homoscedastic fit, whose scale is that scalar.
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
# reporting a wrong number. A heteroscedastic gaussian fit scores at its own
# per-observation s(x) instead of the scalar (heteroscedasticScale below).
pointwiseLogLikelihood <- function(object, ev) {
  y <- object[["y"]]
  if (is.null(y)) {
    stop(
      "cannot compute the log-likelihood; fit does not store the training response"
    )
  }
  family <- object[["family"]]
  weights <- object[["weights"]]
  n.draws <- length(ev) %/% length(y)
  y <- rep(y, each = n.draws)
  n.chains <- fitNChains(object)
  chainFastest <- function(x) {
    if (is.null(dim(x))) uncombineChains(as.vector(x), n.chains) else x
  }

  if (identical(family, "gaussian")) {
    # resid.dist records the fitted residual law; a fit predating the field
    # carries no element at all and is read as gaussian, the historical
    # behavior, so old serialized fits keep working. "student" scores the
    # MARGINAL t density - the observation-level likelihood loo/waic are
    # defined on, and the density the engine itself reports - rather than the
    # gaussian working likelihood conditional on the latent precisions, which
    # is a different quantity. Any other token is refused rather than silently
    # scored against a formula that does not fit it.
    residDist <- object[["resid.dist"]]
    isStudent <- identical(residDist, "student")
    if (
      !is.null(residDist) && !isStudent && !identical(residDist, "gaussian")
    ) {
      stop(
        "pointwise log-likelihood does not support ",
        residDist,
        " residuals"
      )
    }
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
    # sigma and y are on the log-time scale (y is log event time for an event,
    # log censoring time for a censored row); events keep the normal density,
    # censored rows take the log upper survival tail log P(log T > log C)
    sd <- rep_len(as.vector(chainFastest(object$sigma)), length(ev))
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
  array(result, dim(ev))
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
  # 'sample' sat in the slot 'ci.level' now holds, so a "train"/"test" value
  # arriving here is a call written against the old order: name the argument
  # that moved rather than one the caller never wrote
  if (
    is.character(ci.level) &&
      length(ci.level) == 1L &&
      !is.na(ci.level) &&
      ci.level %in% c("train", "test")
  ) {
    stop(
      "'sample' is fitted's fourth argument and is matched by name; write ",
      "sample = \"",
      ci.level,
      "\""
    )
  }
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
  est <- apply(draws, keepAxes, mean)
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

combineOrUncombineChains <- function(x, n.chains, combine) {
  if (n.chains > 1L) {
    if (length(dim(x)) > 2L && combine) {
      x <- combineChains(x)
    } else if (length(dim(x)) == 2L && !combine) {
      x <- uncombineChains(x, n.chains)
    }
  }
  x
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
# fit-time channels keep 'offset.test' (dbartsData, bart2, rbart_vi, and the
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

# bart spells it 'keeptrees', bart2 and rbart_vi 'keepTrees'. A fit kept with
# keepCall = FALSE stores call("NULL") and names neither, so it takes bart's
# spelling, which is the surface such a fit most likely came from.
bartKeepTreesArgument <- function(object) {
  if (callName(object[["call"]]) == "bart2") "keepTrees" else "keeptrees"
}

# An offset shifts the latent at rows the sampler never saw, and these two
# families replay their trees with no offset channel at all, so either spelling
# would be dropped rather than applied.
noPredictOffsetReason <- paste0(
  "this fit has no out-of-sample offset channel; predict replays the ",
  "offset-free surface"
)
predictNoOffsetUnusedArgs <- list(
  offset = noPredictOffsetReason,
  offset.test = noPredictOffsetReason
)

# 'offset' occupies the same slot on all six predict methods so the argument
# order is uniform in position and not merely in relative order; the two
# families with no offset channel take it as a formal and refuse a non-NULL
# value with the same wording it would carry out of '...'.
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

predict.bart <- function(
  object,
  newdata,
  type = c("ev", "ppd", "bart", "forest"),
  offset = NULL,
  weights = NULL,
  combineChains = TRUE,
  ci.level = NULL,
  forest = NULL,
  bases = NULL,
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
  type <- validateType(type, eval(formals(predict.bart)$type))
  # above the type = "forest" and amplitude-blend returns below, so every arm's
  # value is checked rather than only the one that reaches the sampler here
  n.threads <- validatePredictThreads(n.threads)
  refuseForestSelectionOutsideForestArm(type, forest)

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
    return(predictForest(
      object,
      newdata,
      offset,
      combineChains,
      forest,
      n.threads
    ))
  }

  # an amplitude-coupled fit has no combined test surface in the engine - the
  # sampler holds no basis at the caller's rows - so the combination is done
  # here, from the per-forest replay and the fit's own glue
  if (!is.null(object[["forestFits"]])) {
    return(predictBlend(
      object,
      newdata,
      offset,
      weights,
      type,
      combineChains,
      ci.level,
      bases,
      n.threads
    ))
  }

  if (!is.null(bases)) {
    numForests <- if (is.null(object[["n.forests"]])) 1L else object$n.forests
    stop(
      "'bases' is only meaningful on an amplitude-coupled multi-forest fit; ",
      "this fit has ",
      numForests,
      if (numForests == 1L) " forest" else " forests"
    )
  }

  n.chains <- object$fit$control@n.chains
  result <- object$fit$predict(newdata, offset, n.threads)
  # a heteroscedastic fit returns list(mean, variance); s(x) rides back as an
  # attribute on the returned yhat so plain predict callers are unaffected
  s <- NULL
  if (is.list(result)) {
    s <- sqrt(convertSamplesFromDbartsToBart(
      result$variance,
      n.chains,
      combineChains
    ))
    result <- result$mean
  }
  # result is n.obs x n.samples x n.chains
  result <- convertSamplesFromDbartsToBart(result, n.chains, combineChains)

  if (type != "bart") {
    # nolint next: object_usage_linter. named for readability; value drives if.
    if ((responseIsBinary <- is.null(object[["sigma"]]))) {
      result <- probabilityFromLatents(result, object)
    }

    if (type == "ppd") {
      # the replayed s(x) above IS the noise scale at these rows; a
      # heteroscedastic fit whose sampler replays none cannot be drawn from
      if (is.null(s) && !is.null(object[["s.train"]])) {
        stop(
          "posterior predictive sampling is not available on a ",
          "heteroscedastic fit whose sampler replays no variance surface"
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
    return(interval)
  }

  if (!is.null(s)) {
    attr(result, "s") <- s
  }
  result
}

# extract(type = "trees") rewrites the matched call onto the sampler's
# getTrees(treeNums, chainNums, sampleNums, current, newdata); none of the
# extract method's own vocabulary reaches getTrees (bart.Rd's 'Extracting
# Trees' section documents only chainNums/sampleNums/treeNums/newdata as
# accepted there), so a caller-supplied argument that collides by name -
# sample, combineChains, forest, contribution - is refused by name instead of
# being left to partial-match one of getTrees's differently-named formals
# (sample -> sampleNums) or fall through to a raw 'unused argument'.
refuseTreesArguments <- function(treesCall, ownNames) {
  supplied <- intersect(ownNames, names(treesCall))
  if (length(supplied) > 0L) {
    stop(
      "'",
      supplied[1L],
      "' is not used when type = \"trees\"; the sampler's getTrees accepts ",
      "'chainNums', 'sampleNums', 'treeNums', 'current', and 'newdata' ",
      "instead (see 'Extracting Trees' in ?bart)"
    )
  }
  invisible(NULL)
}

extract.bart <- function(
  object,
  type = c("ev", "ppd", "bart", "loglik", "trees", "forest"),
  sample = c("train", "test"),
  combineChains = TRUE,
  forest = NULL,
  contribution = FALSE,
  ...
) {
  type <- validateType(type, eval(formals(extract.bart)$type))

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
      c("sample", "combineChains", "forest", "contribution")
    )
    target <- quote(object$fit$getTrees)
    target[[2L]][[2L]] <- treesCall$object
    treesCall[[1L]] <- target
    treesCall$object <- NULL
    treesCall$type <- NULL
    return(eval(treesCall, parent.frame()))
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

  sample <- validateSample(sample, eval(formals(extract.bart)$sample))
  refuseForestSelectionOutsideForestArm(type, forest)
  if (type != "forest" && isTRUE(contribution)) {
    stop(
      "type = \"",
      type,
      "\" does not support 'contribution': the ",
      "per-observation decomposition applies to the per-forest channel alone"
    )
  }

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
      "cannot extract test sample predictions if no test data exists; use 'predict' instead"
    )
  }
  if (sample == "train" && is.null(object[["yhat.train"]])) {
    if (callName(object$call) == "bart2") {
      stop(
        "cannot extract train sample predictions; bart2 must be called with 'keepTrainingFits' == TRUE"
      )
    } else {
      stop(
        "cannot extract train sample predictions; bart must be called with 'keeptrainfits' == TRUE"
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

  # nolint next: object_usage_linter. named for readability; value drives if.
  if ((responseIsBinary <- is.null(object[["sigma"]]))) {
    result <- probabilityFromLatents(result, object)
  }

  if (type == "ppd") {
    s <- if (sample == "train") object[["s.train"]] else object[["s.test"]]
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
  idx <- suppressWarnings(as.integer(forest))
  if (anyNA(idx) || any(idx < 1L | idx > length(forestNames))) {
    stop("'forest' index must be between 1 and ", length(forestNames))
  }
  idx
}

# extract(type = "forest"): the packaged per-forest response-scale raw total
# by default (forestFits already carries response.scale), or its
# per-observation contribution under contribution = TRUE, computed on
# demand as (basis %*% glue) * raw rather than stored. The selected forests
# always keep the trailing forest margin, even at length one. Refuses by name
# on a fit without forest reporting (the amplitude coupling, not the forest
# count) and on sample = "test" (an amplitude-coupled fit has no test fits).
extractForest <- function(object, sample, combineChains, forest, contribution) {
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
    dimnames = list(NULL, NULL, forestNames[idx])
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
  n.threads
) {
  if (is.null(object[["forestFits"]])) {
    stop(
      "type = \"forest\" is only available on a fit with per-forest ",
      "reporting (an amplitude-coupled multi-forest fit); this fit has none"
    )
  }
  n.chains <- object$fit$control@n.chains
  responseScale <- object$fit$getCalibration(1L)[1L, "response.scale"]
  raw <- object$fit$predictForests(newdata, offset, n.threads) * responseScale
  # forestFits carries the fit's own combineChains shape (3-d combined, 4-d
  # split across chains), so the forest margin is always the LAST axis rather
  # than a fixed index
  forestNames <- dimnames(object$forestFits)[[length(dim(object$forestFits))]]
  idx <- resolveForestSelection(forest, forestNames)
  result <- shapeMultinomialChannel(raw, forestNames, n.chains, TRUE)
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
  n.threads
) {
  n.chains <- object$fit$control@n.chains
  perForest <- predictForest(object, newdata, NULL, TRUE, NULL, n.threads)
  n.new <- dim(perForest)[2L]
  forestNames <- dimnames(perForest)[[3L]]
  bases <- resolveForestBases(object, bases, newdata, n.new)

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
    object$fit$getCalibration(1L)[1L, "response.shift"],
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
  result <- combineOrUncombineChains(result, n.chains, combineChains)

  if (type != "bart") {
    if (is.null(object[["sigma"]])) {
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
  ci.level = NULL,
  sample = c("train", "test"),
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

  # ci.level opts into a per-observation est + credible band instead of the
  # posterior mean; the interval kind follows type (see posteriorInterval)
  if (!is.null(ci.level)) {
    return(posteriorInterval(result, ci.level))
  }

  if (!is.null(dim(result))) {
    apply(result, length(dim(result)), mean)
  } else {
    mean(result)
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
  object$y - fitted.bart(object, type = type, sample = "train")
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
  d <- length(dim(probs))
  meanProbs <- apply(probs, c(d - 1L, d), mean)
  dimnames(meanProbs) <- list(NULL, levels)
  meanProbs
}
categoryFromMeanProbabilities <- function(meanProbs, levels, ordered = FALSE) {
  factor(
    levels[max.col(meanProbs, ties.method = "first")],
    levels = levels,
    ordered = ordered
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
  type = c("ev", "ppd", "bart", "forest", "loglik"),
  sample = c("train", "test"),
  combineChains = TRUE,
  ...
) {
  type <- validateType(type, eval(formals(extract.bartMultinomial)$type))
  sample <- validateSample(
    sample,
    eval(formals(extract.bartMultinomial)$sample)
  )
  refuseMultinomialLatentType(type)
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
# This does not contradict the engine's own per-observation channel staying
# undefined for this family: that channel cannot see the K-blend, while this
# is computed from the reported probabilities R already holds.
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
  term <- rowSums(counts[idx, , drop = FALSE] * log(flat))
  array(rep(logCoef, each = n.draws) + term, d[-length(d)])
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
  array(codes, d[-length(d)])
}

# fitted values are always the training rows (extract's own 'sample' formal
# reaches the test channel); a caller-supplied 'sample' is refused by name
# instead of vanishing into '...' unused, as it does today.
multinomialFittedSampleReason <- list(
  sample = paste0(
    "fitted values are always the training rows; use extract(object, ",
    "sample = \"test\")"
  )
)

# The posterior-mean n x K probability matrix (colnames = levels(y)), or
# (type = "class") the argmax category of that mean as a factor over the
# original levels - the class-prediction convenience. ci.level opts into a
# per-(observation, category) credible band instead of the posterior mean,
# taken on the full probability draws before the class reduction so it is
# meaningful regardless of 'type'.
fitted.bartMultinomial <- function(
  object,
  type = c("ev", "class", "bart"),
  ci.level = NULL,
  ...
) {
  type <- validateType(type, eval(formals(fitted.bartMultinomial)$type))
  refuseMultinomialLatentType(type)
  refuseUnusedGenericArgs(
    list(...),
    "fitted",
    "bartMultinomial",
    c(
      multinomialUnusedArgs,
      multinomialFittedSampleReason,
      foreignArgsFor(
        fittedForeignReasons,
        names(formals(fitted.bartMultinomial))
      )
    )
  )
  refuseClassCiLevel(type, ci.level)
  probs <- extract.bartMultinomial(object, type = "ev", sample = "train")
  if (!is.null(ci.level)) {
    return(posteriorInterval(probs, ci.level, trailing = 2L))
  }
  meanProbs <- meanCategoryProbabilities(probs, object$levels)
  if (type == "ev") {
    return(meanProbs)
  }
  categoryFromMeanProbabilities(meanProbs, object$levels)
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
  phat <- fitted.bartMultinomial(object, type = "ev")
  y <- object$y
  observed <- if (is.factor(y)) {
    indicator <- matrix(0, length(y), ncol(phat), dimnames = dimnames(phat))
    indicator[cbind(seq_along(y), match(y, object$levels))] <- 1
    indicator
  } else {
    y / rowSums(y)
  }
  observed - phat
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
  refuseClassCiLevel(type, ci.level)
  if (is.null(object[["fit"]]) || !object$fit$control@keepTrees) {
    refuseWithoutTrees("predict")
  }
  # after the fit check, whose absence the default here would otherwise report
  # as a missing slot
  n.threads <- validatePredictThreads(n.threads)
  newdata <- validateXTest(newdata, object$fit$data@x)
  if (is.null(newdata)) {
    stop("newdata cannot be NULL")
  }
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
  # raw is n.new x K x n.samples (x n.chains), the run's test-channel shape;
  # $fit$predict re-validates the coded newdata (idempotently) and shapes the
  # offset against the same rows
  raw <- object$fit$predict(newdata, offset, n.threads)
  probs <- shapeMultinomialChannel(
    raw,
    object$levels,
    object$n.chains,
    combineChains
  )
  if (type == "ppd") {
    probs <- multinomialPpdFromProbs(probs)
  }
  if (!is.null(ci.level)) {
    return(posteriorInterval(
      probs,
      ci.level,
      trailing = if (type %in% c("ev", "class")) 2L else 1L
    ))
  }
  if (type == "class") {
    meanProbs <- meanCategoryProbabilities(probs, object$levels)
    return(categoryFromMeanProbabilities(meanProbs, object$levels))
  }
  probs
}

# Shared "Call:" preamble for the print methods. A fit kept with
# keepCall = FALSE stores call("NULL") as a placeholder, which is suppressed.
printCall <- function(x) {
  if (!identical(x[["call"]], call("NULL"))) {
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
  # out - dividing by n.chains == 1 is a no-op, so this is correct either way
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
# the latent and the sampled cutpoints. type = "ppd" draws one category per
# posterior draw. The K-1 cutpoint draws ride the fit's $cutpoints field.
# ordinal has a single forest, so 'forest'/'contribution' refuse for the same
# reason a bart-family single-forest fit does.
ordinalUnusedArgs <- list(
  forest = singleForestReason,
  contribution = singleForestReason
)

extract.bartOrdinal <- function(
  object,
  type = c("ev", "ppd", "bart", "loglik"),
  sample = c("train", "test"),
  combineChains = TRUE,
  ...
) {
  type <- validateType(type, eval(formals(extract.bartOrdinal)$type))
  sample <- validateSample(sample, eval(formals(extract.bartOrdinal)$sample))
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
# verified to machine precision against a hand-built Phi difference, so no
# recomputation from eta/cutpoints is needed. probs enters in the split
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
  array(result, d[-length(d)])
}

ordinalFittedSampleReason <- list(
  sample = paste0(
    "fitted values are always the training rows; use extract(object, ",
    "sample = \"test\")"
  )
)

# The posterior-mean n x K probability matrix (colnames = levels), or
# (type = "class") the argmax category as an ordered factor over the original
# levels, or (type = "bart") the posterior-mean latent eta per observation.
# ci.level opts into a credible band instead of the posterior mean, taken on
# the full draws before any mean/class reduction.
fitted.bartOrdinal <- function(
  object,
  type = c("ev", "class", "bart"),
  ci.level = NULL,
  ...
) {
  type <- validateType(type, eval(formals(fitted.bartOrdinal)$type))
  refuseUnusedGenericArgs(
    list(...),
    "fitted",
    "bartOrdinal",
    c(
      ordinalUnusedArgs,
      ordinalFittedSampleReason,
      foreignArgsFor(fittedForeignReasons, names(formals(fitted.bartOrdinal)))
    )
  )
  refuseClassCiLevel(type, ci.level)
  if (type == "bart") {
    latent <- object$latent.train
    if (!is.null(ci.level)) {
      return(posteriorInterval(latent, ci.level, trailing = 1L))
    }
    return(apply(latent, length(dim(latent)), mean))
  }
  probs <- extract.bartOrdinal(object, type = "ev", sample = "train")
  if (!is.null(ci.level)) {
    return(posteriorInterval(probs, ci.level, trailing = 2L))
  }
  meanProbs <- meanCategoryProbabilities(probs, object$levels)
  if (type == "ev") {
    return(meanProbs)
  }
  categoryFromMeanProbabilities(meanProbs, object$levels, ordered = TRUE)
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
  indicator <- matrix(0, length(y), ncol(phat), dimnames = dimnames(phat))
  indicator[cbind(seq_along(y), match(y, object$levels))] <- 1
  indicator - phat
}

# Out-of-sample category probabilities by replaying the saved forest's trees to
# the newdata latent, then differencing the cumulative probit at the STORED
# per-draw cutpoints. Requires a fit kept with
# keepTrees. type = "bart" returns the replayed latent eta; type = "ppd" draws
# one category per posterior draw. Only ppd touches the RNG, so type = "ev" is
# draw-neutral. The replay reads through $fit's own pointer: $fit is the
# sampler whose engine actually ran, so getPointer() can re-create it from
# stored state after a save/reload. The presence gate re-points to
# cutpoints.raw, which this function already reads below and which rides the
# same keepTrees gate a deleted $bc field used to.
predict.bartOrdinal <- function(
  object,
  newdata,
  type = c("ev", "ppd", "bart", "class"),
  offset = NULL,
  combineChains = TRUE,
  ci.level = NULL,
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
      predictNoOffsetUnusedArgs,
      foreignArgsFor(predictForeignReasons, names(formals(predict.bartOrdinal)))
    )
  )
  refusePredictOffsetChannel(offset, "bartOrdinal")
  refuseClassCiLevel(type, ci.level)
  if (is.null(object[["cutpoints.raw"]])) {
    refuseWithoutTrees("predict")
  }
  # after the store check, whose absence the default here would otherwise
  # report as a missing slot
  n.threads <- validatePredictThreads(n.threads)
  newdata <- validateXTest(newdata, object$fit$data@x)
  if (is.null(newdata)) {
    stop("newdata cannot be NULL")
  }
  if (!is.matrix(newdata)) {
    newdata <- as.matrix(newdata)
  }
  n.chains <- object$n.chains
  # raw is n.new x n.samples (x n.chains): the replayed latent eta, the test
  # channel's shape
  raw <- bartcorePredict(
    list(ptr = object$fit$getPointer()),
    newdata,
    n.threads = n.threads
  )
  if (type == "bart") {
    result <- convertSamplesFromDbartsToBart(raw, n.chains, combineChains)
    if (!is.null(ci.level)) {
      return(posteriorInterval(result, ci.level, trailing = 1L))
    }
    return(result)
  }
  K <- object$K
  cutpoints <- object$cutpoints.raw # (K-1) x n.samples x n.chains
  if (length(dim(raw)) == 2L) {
    dim(raw) <- c(dim(raw), 1L)
  }
  n.new <- dim(raw)[1L]
  n.samples <- dim(raw)[2L]
  probs <- array(0, c(n.new, K, n.samples, n.chains))
  for (s in seq_len(n.samples)) {
    for (chain in seq_len(n.chains)) {
      probs[,, s, chain] <-
        ordinalCategoryProbabilities(raw[, s, chain], cutpoints[, s, chain])
    }
  }
  if (n.chains == 1L) {
    probs <- array(probs, dim(probs)[1:3])
  }
  probs <- shapeMultinomialChannel(
    probs,
    object$levels,
    n.chains,
    combineChains
  )
  if (type == "ppd") {
    probs <- multinomialPpdFromProbs(probs)
  }
  if (!is.null(ci.level)) {
    trailing <- if (type %in% c("ev", "class")) 2L else 1L
    return(posteriorInterval(probs, ci.level, trailing = trailing))
  }
  if (type == "class") {
    meanProbs <- meanCategoryProbabilities(probs, object$levels)
    return(categoryFromMeanProbabilities(
      meanProbs,
      object$levels,
      ordered = TRUE
    ))
  }
  probs
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
# never "bart" - so the count arrays
# never fall through to the single-forest "bart" methods. A single forest fits
# the log-odds latent psi = f(x) + o, so type = "bart" returns that latent (like
# probit/logistic), while type = "ev" returns the mean counts mu = r exp(psi)
# (the reported posterior mean count) and type = "ppd" draws one count per posterior draw
# from NB(r, plogis(psi)). The per-draw dispersion r rides the fit's $dispersion
# field, the count analog of gaussian's sigma.
# nbinom has a single forest, so 'forest'/'contribution' refuse for the same
# reason a bart-family single-forest fit does.
negbinUnusedArgs <- list(
  forest = singleForestReason,
  contribution = singleForestReason
)

extract.bartNegbin <- function(
  object,
  type = c("ev", "ppd", "bart", "loglik"),
  sample = c("train", "test"),
  combineChains = TRUE,
  ...
) {
  type <- validateType(type, eval(formals(extract.bartNegbin)$type))
  sample <- validateSample(sample, eval(formals(extract.bartNegbin)$sample))
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
  # type == "ppd": pair mu with dispersion in a common split layout so the two
  # align regardless of either's own storage, then reshape the result to the
  # caller's request
  muSplit <- combineOrUncombineChains(mu, n.chains, FALSE)
  disp <- scalarDrawVec(object$dispersion, n.chains, length(muSplit))
  result <- array(
    rnbinom(length(muSplit), size = disp, mu = as.vector(muSplit)),
    dim(muSplit)
  )
  combineOrUncombineChains(result, n.chains, combineChains)
}

# l[s,i] = dnbinom(y_i, size = dispersion[s], mu = yhat.train[s,i]); the
# per-draw dispersion pairs with the draws the same chain-fastest way the
# gaussian arm pairs sigma (dispersion is already sigma-shaped). mu enters
# forced to the split (chains x) samples x obs layout so it aligns with
# scalarDrawVec's own normalization regardless of either's own storage.
negbinLogLik <- function(object, mu, n.chains) {
  y <- object[["y"]]
  n.draws <- length(mu) %/% length(y)
  disp <- scalarDrawVec(object[["dispersion"]], n.chains, length(mu))
  result <- dnbinom(
    rep(y, each = n.draws),
    size = disp,
    mu = as.vector(mu),
    log = TRUE
  )
  array(result, dim(mu))
}

negbinFittedSampleReason <- list(
  sample = paste0(
    "fitted values are always the training rows; use extract(object, ",
    "sample = \"test\")"
  )
)

# The posterior-mean count per observation (type = "ev"), the posterior-mean
# log-odds latent per observation (type = "bart"), or a Monte Carlo mean over
# ppd draws (type = "ppd"). The observation margin is the array's last
# dimension in every chain layout, so we take the mean over that observation
# margin. ci.level opts into a credible band instead, taken on the full draws
# before the mean.
fitted.bartNegbin <- function(
  object,
  type = c("ev", "ppd", "bart"),
  ci.level = NULL,
  ...
) {
  type <- validateType(type, eval(formals(fitted.bartNegbin)$type))
  refuseUnusedGenericArgs(
    list(...),
    "fitted",
    "bartNegbin",
    c(
      negbinUnusedArgs,
      negbinFittedSampleReason,
      foreignArgsFor(fittedForeignReasons, names(formals(fitted.bartNegbin)))
    )
  )
  channel <- switch(
    type,
    bart = object$latent.train,
    ev = object$yhat.train,
    # the ppd arm is a draw, not a stored channel; extract pairs each mu with
    # its own draw's dispersion, and the mean over the observation margin
    # below is invariant to the chain layout it returns
    ppd = extract.bartNegbin(object, type = "ppd", sample = "train")
  )
  if (!is.null(ci.level)) {
    return(posteriorInterval(channel, ci.level, trailing = 1L))
  }
  apply(channel, length(dim(channel)), mean)
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
  object$y - fitted.bartNegbin(object, type = "ev")
}

# Out-of-sample mean counts by replaying the saved forest's trees to the newdata
# log-odds latent psi, then mu = r exp(psi) at the STORED per-draw
# dispersion. A log-exposure offset.test enters psi
# additively, the fit-time convention. Requires a fit kept with keepTrees.
# type = "bart" returns the replayed latent; type = "ppd" draws one count per
# posterior draw. Only ppd touches the RNG, so type = "ev" is draw-neutral. The
# replay reads through $fit's own pointer: $fit is the sampler whose engine
# actually ran, so getPointer() can re-create it from stored state after a
# save/reload. The presence gate re-points to dispersion.raw, which is read
# below and rides the same keepTrees gate a deleted $bc field used to.
predict.bartNegbin <- function(
  object,
  newdata,
  type = c("ev", "ppd", "bart"),
  offset = NULL,
  combineChains = TRUE,
  ci.level = NULL,
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
  if (is.null(object[["dispersion.raw"]])) {
    refuseWithoutTrees("predict")
  }
  # after the store check, whose absence the default here would otherwise
  # report as a missing slot
  n.threads <- validatePredictThreads(n.threads)
  newdata <- validateXTest(newdata, object$fit$data@x)
  if (is.null(newdata)) {
    stop("newdata cannot be NULL")
  }
  if (!is.matrix(newdata)) {
    newdata <- as.matrix(newdata)
  }
  n.chains <- object$n.chains
  # raw is n.new x n.samples (x n.chains): the replayed log-odds latent psi
  raw <- bartcorePredict(
    list(ptr = object$fit$getPointer()),
    newdata,
    offset,
    n.threads
  )
  if (type == "bart") {
    result <- convertSamplesFromDbartsToBart(raw, n.chains, combineChains)
    if (!is.null(ci.level)) {
      return(posteriorInterval(result, ci.level, trailing = 1L))
    }
    return(result)
  }
  if (length(dim(raw)) == 2L) {
    dim(raw) <- c(dim(raw), 1L)
  }
  disp <- object$dispersion.raw # n.samples x n.chains
  n.new <- dim(raw)[1L]
  n.samples <- dim(raw)[2L]
  means <- array(0, c(n.new, n.samples, n.chains))
  for (s in seq_len(n.samples)) {
    for (chain in seq_len(n.chains)) {
      means[, s, chain] <- negbinMeanCounts(raw[, s, chain], disp[s, chain])
    }
  }
  if (n.chains == 1L) {
    means <- matrix(means, n.new, n.samples)
  }
  means <- convertSamplesFromDbartsToBart(means, n.chains, combineChains)
  if (type == "ppd") {
    means <- negbinPpd(means, object$dispersion)
  }
  if (!is.null(ci.level)) {
    return(posteriorInterval(means, ci.level, trailing = 1L))
  }
  means
}

print.bartNegbin <- function(x, ...) {
  printCall(x)
  cat("family: negative binomial (log-odds latent)\n")
  cat(
    "posterior mean dispersion (r): ",
    format(mean(x$dispersion), digits = 4L),
    "\n",
    sep = ""
  )
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
# conditionally-independent component fits ($occupancy, a probit fit of
# 1{y > 0} over all n; $positive, a gaussian fit of log(y) over the y > 0
# subset whose x.test is the full-n x). The report-time combine glues their
# posterior draws by sample index (any pairing is a valid joint draw, the parts
# share no parameters) and retransforms the positive part to the natural scale.
#
# type = "prob" is the occupancy probability pi(x); type = "bart"/"link"/"log"
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
# %not_in% + stop shared by the bart/rbart predict/extract/fitted methods.
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
# a caller-supplied name out of '...' and stopping on the first hit - rather
# than letting it fall through silently, as it does today - refuses these by
# name instead of discarding them.
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
# reports; every other arm has already recombined them into the reported
# location, so a selection there would silently choose nothing.
refuseForestSelectionOutsideForestArm <- function(type, forest) {
  if (type != "forest" && !is.null(forest)) {
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
# shape has no room for; bart and rbart did not, so the two names a
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
rbartUnusedArgs <- list(
  forest = singleForestReason,
  contribution = singleForestReason
)

# The derived reason tables the surface's predict/extract/fitted/residuals
# methods compose via foreignArgsFor above, one entry per name that is a
# formal on some method of the generic and foreign on another. Composed
# AFTER a method's own class list (multinomialUnusedArgs and siblings,
# bartUnusedArgs, rbartUnusedArgs), so a class-specific reason for the same
# name still wins - refuseUnusedGenericArgs reports the first hit in
# 'reasons', and composition order is priority order.
predictForeignReasons <- list(
  sample = "the fit's stored train and test channels are extract's 'sample'",
  weights = "this family's posterior-predictive draw takes no per-observation weight",
  bases = "only an amplitude-coupled multi-forest fit takes 'bases' at the predicted rows",
  group.by = "'group.by' is the grouped (rbart_vi) fit's own predict argument",
  contribution = "the per-observation contribution decomposition belongs to extract(type = \"forest\")"
)

extractReplaysNothingReason <- "extract reads stored channels and replays nothing"
extractForeignReasons <- list(
  ci.level = "extract returns the draws that fitted() and predict() take a band over",
  newdata = "predict(object, newdata) is the read at new rows",
  offset = extractReplaysNothingReason,
  weights = extractReplaysNothingReason,
  n.threads = extractReplaysNothingReason,
  bases = extractReplaysNothingReason,
  group.by = "the stored channels already carry the fit's own grouping"
)

fittedSummarizesNothingReason <- "fitted summarizes stored channels and replays nothing"
fittedForeignReasons <- list(
  combineChains = "the per-chain draws are extract(object, combineChains = FALSE)",
  sample = "fitted values are always the fit's training rows",
  newdata = fittedSummarizesNothingReason,
  offset = fittedSummarizesNothingReason,
  weights = fittedSummarizesNothingReason,
  n.threads = fittedSummarizesNothingReason,
  bases = fittedSummarizesNothingReason,
  group.by = "the stored channels already carry the fit's own grouping"
)

residualsSummarizeNothingReason <- "residuals summarize stored channels and replay nothing"
residualsForeignReasons <- list(
  ci.level = "residuals are the observed response minus the posterior-mean fit",
  combineChains = "the per-chain draws are extract(object, combineChains = FALSE)",
  newdata = residualsSummarizeNothingReason,
  offset = residualsSummarizeNothingReason,
  weights = residualsSummarizeNothingReason,
  n.threads = residualsSummarizeNothingReason,
  bases = residualsSummarizeNothingReason,
  group.by = "the stored channels already carry the fit's own grouping"
)

survivalProbabilitiesDrawsReason <- "survivalProbabilities returns the draws of S(t | x) at 'times'"
survivalProbabilitiesOwnArgsReason <- "survivalProbabilities takes 'times' and 'newdata' alone"
survivalProbabilitiesForeignReasons <- list(
  group.by = "'group.by' is the grouped (rbart_vi) fit's own argument",
  type = survivalProbabilitiesDrawsReason,
  sample = survivalProbabilitiesDrawsReason,
  ci.level = survivalProbabilitiesDrawsReason,
  offset = survivalProbabilitiesOwnArgsReason,
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
  occupancy <- object$occupancy
  if (!is.null(occupancy[["fit"]])) {
    occupancy$fit$control@n.chains
  } else {
    occupancy$n.chains
  }
}

# A scalar-per-draw field (sigma, dispersion, ...) as a flat vector aligned,
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

# Glue the flat, draw-aligned occupancy-probability, positive-log-mean, and
# positive-sigma vectors into the requested channel and reshape to the fit's
# uncombined draw layout ('shape'). Only "ppd" touches the RNG (Bernoulli then
# lognormal), so the default "ev" is draw-neutral.
combineHurdleChannel <- function(type, piVec, fVec, sigmaVec, shape) {
  channel <- switch(
    type,
    prob = piVec,
    bart = fVec,
    ev = piVec * exp(fVec + 0.5 * sigmaVec^2),
    ppd = rbinom(length(piVec), 1L, piVec) *
      exp(fVec + sigmaVec * rnorm(length(fVec)))
  )
  array(channel, shape)
}

# The occupancy pi(x), positive log-mean f(x), and positive per-observation
# sigma draws for the combine, each a flat vector in the fit's uncombined
# as.vector order, plus the uncombined 'shape' to fold back to. In-sample reads
# the stored channels - the occupancy fit's ev over all n, and the positive
# fit's log-scale (bart) fits at the FULL-n rows through its x.test channel (the
# zero rows it never trained on included); out-of-sample replays both saved
# forests at newdata.
hurdleParts <- function(object, newdata, n.threads = 1L) {
  if (missing(newdata)) {
    pi <- extract(
      object$occupancy,
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
    pi <- predict(
      object$occupancy,
      newdata,
      type = "ev",
      combineChains = FALSE,
      n.threads = n.threads
    )
    f <- predict(
      object$positive,
      newdata,
      type = "bart",
      combineChains = FALSE,
      n.threads = n.threads
    )
  }
  sigmaVec <- scalarDrawVec(
    object$positive$sigma,
    hurdleNChains(object),
    length(f)
  )
  list(pi = as.vector(pi), f = as.vector(f), sigma = sigmaVec, shape = dim(f))
}

finishHurdle <- function(parts, type, n.chains, combineChains, ci.level) {
  channel <- combineHurdleChannel(
    type,
    parts$pi,
    parts$f,
    parts$sigma,
    parts$shape
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
  type = c("ev", "ppd", "prob", "bart", "loglik"),
  sample = c("train", "test"),
  combineChains = TRUE,
  ...
) {
  type <- resolveHurdleType(type, eval(formals(extract.bartHurdle)$type))
  sample <- validateSample(sample, eval(formals(extract.bartHurdle)$sample))
  refuseUnusedGenericArgs(
    list(...),
    "extract",
    "bartHurdle",
    c(
      hurdleUnusedArgs,
      foreignArgsFor(extractForeignReasons, names(formals(extract.bartHurdle)))
    )
  )
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
# already glue, flat and draw-aligned, at ALL n rows (the occupancy's own
# channel; the positive part's x.test channel, zero rows included). y == 0
# rows take the occupancy's own log(1 - pi); y > 0 rows take the occupancy
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
  array(result, parts$shape)
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
  # sample = "test" unconditionally), so the read is always the training rows
  draws <- extract(object, type = type, sample = "train", combineChains = TRUE)
  if (!is.null(ci.level)) {
    return(posteriorInterval(draws, ci.level))
  }
  apply(draws, length(dim(draws)), mean)
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
  object$y - fitted.bartHurdle(object, type = type)
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
  n.threads = object$occupancy$fit$control@n.threads,
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
  refusePredictOffsetChannel(offset, "bartHurdle")
  if (is.null(object$occupancy[["fit"]])) {
    refuseWithoutTrees("predict")
  }
  # after the occupancy fit check, whose absence the default here would
  # otherwise report as a missing slot
  n.threads <- validatePredictThreads(n.threads)
  n.chains <- hurdleNChains(object)
  finishHurdle(
    hurdleParts(object, newdata, n.threads),
    type,
    n.chains,
    combineChains,
    ci.level
  )
}

print.bartHurdle <- function(x, ...) {
  printCall(x)
  cat("family: hurdle.lognormal (probit occupancy + lognormal positive part)\n")
  cat("occupancy n (all rows): ", length(x$occupancy$y), "\n", sep = "")
  cat("positive-part n (y > 0): ", length(x$positive$y), "\n", sep = "")
  invisible(x)
}

# 'value' was predict.rbart's pre-1.0 name for 'type'. It is not accepted,
# only refused by name, since a supplied one would otherwise choose the
# default channel silently.
rbartPredictValueUnusedArgs <- list(
  value = "predict's channel argument is named 'type'"
)

predict.rbart <- function(
  object,
  newdata,
  type = c("ev", "ppd", "bart", "ranef"),
  offset = NULL,
  weights = NULL,
  combineChains = TRUE,
  ci.level = NULL,
  n.threads = object$fit[[1L]]$control@n.threads,
  ...,
  group.by
) {
  if (is.null(object$fit)) {
    refuseWithoutTrees("predict")
  }
  if (missing(group.by)) {
    stop(
      "'group.by' must be given by name: predict on an rbart fit needs the ",
      "test rows' grouping factor, and it is no longer the third positional ",
      "argument"
    )
  }
  n.threads <- validatePredictThreads(n.threads)

  type <- validateType(type, eval(formals(predict.rbart)$type))
  refuseUnusedGenericArgs(
    list(...),
    "predict",
    "rbart",
    c(
      predictOffsetUnusedArgs,
      rbartPredictValueUnusedArgs,
      rbartUnusedArgs,
      foreignArgsFor(predictForeignReasons, names(formals(predict.rbart)))
    )
  )

  n.chains <- if (is.null(object$n.chains)) {
    length(object$fit)
  } else {
    object$n.chains
  }
  n.samples <- object$fit[[1L]]$control@n.samples

  # coerce as rbart_vi does at fit time so numeric or character group.by carry
  # levels; unlike fitting, unused factor levels are kept, so out-of-sample
  # groups present only in the supplied factor still receive prior draws
  group.by <- as.factor(group.by)

  nonParametricPart <- 0
  # collects results in an array of n.obs x n.samples x n.chains, default for
  # internal sampler
  #
  # read n.obs off the sampler's prediction output (its first dimension),
  # since we would otherwise have to build the test matrix ourselves
  if (type != "ranef") {
    if (n.chains > 1L) {
      if (length(object$fit) == 1L) {
        # the in-core Gibbs path keeps one multi-chain sampler, whose
        # predictions already carry the chain dimension
        nonParametricPart <- object$fit[[1L]]$predict(
          newdata,
          offset,
          n.threads
        )
        n.obs <- dim(nonParametricPart)[1L]
      } else {
        n.obs <- NULL
        nonParametricPart <- array(
          sapply(seq_len(n.chains), function(i) {
            res <- object$fit[[i]]$predict(newdata, offset, n.threads)
            if (is.null(n.obs)) {
              n.obs <<- dim(res)[1L]
            }
            res
          }),
          c(n.obs, n.samples, n.chains)
        )
      }
    } else {
      nonParametricPart <- object$fit[[1L]]$predict(newdata, offset, n.threads)
      n.obs <- nrow(nonParametricPart)
    }
    if (n.obs != length(group.by)) {
      stop("length of group.by not equal to number of rows in test")
    }

    nonParametricPart <- convertSamplesFromDbartsToBart(
      nonParametricPart,
      n.chains,
      combineChains
    )
  }

  if (type == "bart") {
    if (!is.null(ci.level)) {
      return(posteriorInterval(nonParametricPart, ci.level))
    }
    return(nonParametricPart)
  }

  ranef <- 0
  if (type != "bart") {
    ranefNames.test <- levels(group.by)
    ranefNames.train <- if (length(dim(object$ranef)) > 2L) {
      dimnames(object$ranef)[[3L]]
    } else {
      dimnames(object$ranef)[[2L]]
    }

    ranef <- object$ranef
    ranef <- combineOrUncombineChains(ranef, n.chains, combineChains)

    if (!all(measuredLevels <- ranefNames.test %in% ranefNames.train)) {
      warning(
        "test includes random effect levels not present in training (",
        paste0(ranefNames.test[!measuredLevels], collapse = ", "),
        "); ranef estimates default to draws from their latent distribution parameterized by the posterior of its variance, and may not be the same across future calls to 'predict'"
      )
      n.unmeasured <- sum(!measuredLevels)
      if (n.chains > 1L) {
        # object$tau may be stored combined (flat, chain-major) or split
        # ((n.chains x) n.samples matrix, chain-fastest); normalize to the
        # split matrix once, then flatten to whichever order this call's
        # own (un)combined 'ranef' needs below
        tauMat <- if (is.null(dim(object$tau))) {
          uncombineChains(as.vector(object$tau), n.chains)
        } else {
          object$tau
        }
        if (!combineChains) {
          unmeasuredRanef <- array(
            rnorm(
              n.chains * n.samples * n.unmeasured,
              0,
              rep.int(as.vector(tauMat), n.unmeasured)
            ),
            c(n.chains, n.samples, n.unmeasured),
            dimnames = list(NULL, NULL, ranefNames.test[!measuredLevels])
          )
        } else {
          unmeasuredRanef <- matrix(
            rnorm(
              n.chains * n.samples * n.unmeasured,
              0,
              rep.int(as.vector(t(tauMat)), n.unmeasured)
            ),
            n.chains * n.samples,
            n.unmeasured,
            dimnames = list(NULL, ranefNames.test[!measuredLevels])
          )
        }
        # branch on ranef's current shape (already reshaped above to match
        # this call's combineChains, which is what unmeasuredRanef was just
        # built against) rather than object$ranef's stored shape - those two
        # can disagree whenever the fit's own combineChains default differs
        # from the one requested here
        if (length(dim(ranef)) == 2L) {
          ranef <- cbind(ranef, unmeasuredRanef)
        } else {
          # ranef are n.chains x n.samples x n.group
          ranef <- array(
            c(ranef, unmeasuredRanef),
            c(n.chains, n.samples, dim(ranef)[3L] + n.unmeasured),

            dimnames = list(
              NULL,
              NULL,
              c(dimnames(ranef)[[3L]], dimnames(unmeasuredRanef)[[3L]])
            )
          )
        }
      } else {
        unmeasuredRanef <- matrix(
          rnorm(n.samples * n.unmeasured, 0, rep.int(object$tau, n.unmeasured)),
          n.samples,
          n.unmeasured,
          dimnames = list(NULL, ranefNames.test[!measuredLevels])
        )
        ranef <- cbind(ranef, unmeasuredRanef)
      }
    }
  }

  if (type == "ranef") {
    ranef <- if (length(dim(ranef)) > 2L) {
      ranef[,, ranefNames.test, drop = FALSE]
    } else {
      ranef[, ranefNames.test, drop = FALSE]
    }
    ranef <- combineOrUncombineChains(ranef, n.chains, combineChains)
    if (!is.null(ci.level)) {
      return(posteriorInterval(ranef, ci.level))
    }
    return(ranef)
  }

  ranef <- unname(
    if (length(dim(ranef)) > 2L) {
      ranef[,, as.character(group.by), drop = FALSE]
    } else {
      ranef[, as.character(group.by), drop = FALSE]
    }
  )
  ranef <- combineOrUncombineChains(ranef, n.chains, combineChains)

  if (
    length(dim(nonParametricPart)) != length(dim(ranef)) ||
      any(dim(nonParametricPart) != dim(ranef))
  ) {
    stop(
      "internal error: fixed and random effect predictions have ",
      "mismatched dimensions"
    )
  }
  result <- nonParametricPart + ranef

  responseIsBinary <- is.null(object[["sigma"]])
  if (responseIsBinary) {
    result <- probabilityFromLatents(result, object)
  }

  if (type == "ppd") {
    result <- sampleFromPPD(result, object, weights, n.chains)
  }

  if (!is.null(ci.level)) {
    return(posteriorInterval(result, ci.level))
  }

  if (exists("unmeasuredRanef", inherits = FALSE)) {
    attr(result, "ranef") <- unmeasuredRanef
  }

  result
}

extract.rbart <- function(
  object,
  type = c("ev", "ppd", "bart", "loglik", "ranef", "trees"),
  sample = c("train", "test"),
  combineChains = TRUE,
  ...
) {
  type <- validateType(type, eval(formals(extract.rbart)$type))

  n.chains <- if (is.null(object$n.chains)) {
    length(object$fit)
  } else {
    object$n.chains
  }

  if (type == "trees") {
    if (is.null(object$fit)) {
      refuseWithoutTrees("extract(type = \"trees\")")
    }
    treesCall <- match.call()
    refuseTreesArguments(treesCall, c("sample", "combineChains"))
    # the in-core Gibbs path keeps one multi-chain sampler: route chain
    # selection through its chainNums argument instead of the fit list
    singleFit <- length(object$fit) == 1L && n.chains > 1L
    target <- quote(object$fit[[i]]$getTrees)
    target[[2L]][[2L]][[2L]] <- treesCall$object
    if (singleFit) {
      target[[2L]][[3L]] <- 1L
    }
    treesCall[[1L]] <- target
    treesCall$object <- NULL
    treesCall$type <- NULL
    treesCall$chainNums <- if (singleFit) quote(i) else NULL
    evalEnv <- parent.frame()
    dotsList <- list(...)
    chainNums <- if ("chainNums" %in% names(dotsList)) {
      as.integer(dotsList[["chainNums"]])
    } else {
      seq_len(n.chains)
    }
    varOrder <- c("chain", "sample", "tree", "n", "var", "value")
    allTrees <- lapply(chainNums, function(i) {
      result_i <- eval(subTermInLanguage(treesCall, quote(i), i), evalEnv)
      if (n.chains > 1L) {
        result_i$chain <- i
      }
      # varOrder's "chain" is absent for single-chain fits (getTrees omits
      # it); reorder only the columns that exist and keep any others
      # (directions/missing/beta.*) trailing in their original order
      knownOrder <- varOrder[varOrder %in% colnames(result_i)]
      result_i[, c(knownOrder, setdiff(colnames(result_i), knownOrder))]
    })
    if (length(allTrees) > 1L) {
      allTrees <- Reduce(rbind, allTrees)
    } else {
      allTrees <- allTrees[[1L]]
    }
    row.names(allTrees) <- as.character(seq_len(nrow(allTrees)))
    return(allTrees)
  }

  # below the type == "trees" branch and its own refuseTreesArguments, so
  # extract(type = "trees", newdata = ) keeps forwarding to getTrees instead
  # of being refused here for a name that arm alone accepts
  refuseUnusedGenericArgs(
    list(...),
    "extract",
    "rbart",
    c(
      rbartUnusedArgs,
      foreignArgsFor(extractForeignReasons, names(formals(extract.rbart)))
    )
  )

  sample <- validateSample(sample, eval(formals(extract.rbart)$sample))

  # the log-likelihood is against the stored training response; there is no
  # test response to evaluate
  if (type == "loglik" && sample == "test") {
    stop("cannot extract a test sample log-likelihood; no test response exists")
  }

  if (sample == "test" && is.null(object[["yhat.test"]])) {
    stop(
      "cannot extract test sample predictions if no test data exists; use 'predict' instead"
    )
  }

  if (type == "ranef") {
    ranefNames <- if (sample == "train") {
      levels(object$group.by)
    } else {
      levels(object$group.by.test)
    }
    ranef <- if (length(dim(object$ranef)) > 2L) {
      object$ranef[,, ranefNames, drop = FALSE]
    } else {
      object$ranef[, ranefNames, drop = FALSE]
    }
    ranef <- combineOrUncombineChains(ranef, n.chains, combineChains)

    return(ranef)
  }

  if (sample == "train" && is.null(object[["yhat.train"]])) {
    stop(
      "cannot extract train sample predictions; rbart_vi must be called with 'keepTrainingFits' == TRUE"
    )
  }

  # the "ev" draws are the fit's per-draw location (BART component plus the
  # drawn group intercepts), so the log-likelihood conditions on both
  if (type == "loglik") {
    ev <- extract.rbart(
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
  # if necessary, recover chain information or throw it away
  result <- combineOrUncombineChains(result, n.chains, combineChains)

  if (type == "bart") {
    return(result)
  }

  ranefNames <- if (sample == "train") {
    as.character(object$group.by)
  } else {
    as.character(object$group.by.test)
  }
  ranef <- unname(
    if (length(dim(object$ranef)) > 2L) {
      object$ranef[,, ranefNames, drop = FALSE]
    } else {
      object$ranef[, ranefNames, drop = FALSE]
    }
  )

  ranef <- combineOrUncombineChains(ranef, n.chains, combineChains)

  result <- result + ranef

  responseIsBinary <- is.null(object[["sigma"]])
  if (responseIsBinary) {
    result <- probabilityFromLatents(result, object)
  }

  if (type == "ppd") {
    result <- sampleFromPPD(result, object, weights, n.chains)
  }

  result
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
# integer codes, the historical form of data@x and the matrix getTrees
# replays. A dense-frame/mixed container materializes through as.matrix; a
# plain matrix (or a sparse dgCMatrix held as such) is returned unchanged.
extract.dbartsSampler <- function(object, type = "predictors", ...) {
  refuseSamplerExtractArgs(list(...))
  if (!is.character(type) || length(type) == 0L || type[1L] != "predictors") {
    stop("'type' must be one of 'predictors'")
  }
  x <- object$data@x
  if (inherits(x, "dbartsMixedMatrix")) as.matrix(x) else x
}

fitted.rbart <- function(
  object,
  type = c("ev", "ppd", "bart", "ranef"),
  ci.level = NULL,
  sample = c("train", "test"),
  ...
) {
  type <- validateType(type, eval(formals(fitted.rbart)$type))
  sample <- validateSample(sample, eval(formals(fitted.rbart)$sample))
  refuseUnusedGenericArgs(
    list(...),
    "fitted",
    "rbart",
    c(
      rbartUnusedArgs,
      foreignArgsFor(fittedForeignReasons, names(formals(fitted.rbart)))
    )
  )

  if (sample == "train" && type != "ranef" && is.null(object[["yhat.train"]])) {
    stop(
      "cannot extract train sample predictions; rbart_vi must be called with 'keepTrainingFits' == TRUE"
    )
  }

  # ci.level routes through the draws (extract) rather than the mean-only C
  # fast path below, then summarizes to est + credible band (kind follows type)
  if (!is.null(ci.level)) {
    return(posteriorInterval(extract(object, type, sample), ci.level))
  }

  if (type == "ev") {
    ranefNames <- dimnames(object$ranef)
    ranefNames <- ranefNames[[length(ranefNames)]]
    groupByName <- if (sample == "train") "group.by" else "group.by.test"
    groupByMatch <- match(object[[groupByName]], ranefNames)
    # C_rbart_fitted indexes the group dimension of 'ranef' by this vector
    # directly, so a label with no column - a fit whose 'ranef' dimnames were
    # stripped or do not cover every group - has to be refused here
    if (anyNA(groupByMatch)) {
      stop(
        "'",
        groupByName,
        "' must name groups present in the 'ranef' dimnames"
      )
    }
    result <- .Call(
      C_rbart_fitted,
      if (sample == "train") object$yhat.train else object$yhat.test,
      object$ranef,
      groupByMatch,
      is.null(object[["sigma"]])
    )
  } else {
    result <- extract(object, type, sample)

    result <- if (!is.null(dim(result))) {
      apply(result, length(dim(result)), mean)
    } else {
      mean(result)
    }
  }

  result
}

residuals.rbart <- function(object, type = "ev", ...) {
  # as residuals.bart: type reaches fitted for link-scale residuals, sample
  # is pinned to the training response
  refuseResidualsSample(list(...))
  refuseUnusedGenericArgs(
    list(...),
    "residuals",
    "rbart",
    c(
      rbartUnusedArgs,
      foreignArgsFor(residualsForeignReasons, names(formals(residuals.rbart)))
    )
  )
  object$y - fitted.rbart(object, type = type, sample = "train")
}

# fit-level dispatch for the sampler's plotTree method, so a kept bart or
# rbart fit can be plotted directly instead of reaching into $fit; chainNum
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

plotTree.bart <- function(object, treeNum = 1L, chainNum, sampleNum, ...) {
  refusePlotTreeArgs(sys.call())
  if (is.null(object[["fit"]])) {
    refuseWithoutTrees("plotTree", bartKeepTreesArgument(object))
  }
  args <- list(treeNum = treeNum, ...)
  if (!missing(chainNum)) {
    args$chainNum <- chainNum
  }
  if (!missing(sampleNum)) {
    args$sampleNum <- sampleNum
  }
  invisible(do.call(object$fit$plotTree, args))
}

plotTree.rbart <- function(
  object,
  treeNum = 1L,
  chainNum = 1L,
  sampleNum,
  ...
) {
  refusePlotTreeArgs(sys.call())
  if (is.null(object[["fit"]])) {
    refuseWithoutTrees("plotTree")
  }
  n.chains <- if (is.null(object$n.chains)) {
    length(object$fit)
  } else {
    object$n.chains
  }
  chainNum <- as.integer(chainNum)
  if (
    length(chainNum) != 1L ||
      is.na(chainNum) ||
      chainNum < 1L ||
      chainNum > n.chains
  ) {
    stop("'chainNum' must be a single chain index in [1, ", n.chains, "]")
  }
  # the in-core Gibbs path keeps one multi-chain sampler (select the chain
  # through its own chainNum); the R-loop path keeps one sampler per chain
  singleFit <- length(object$fit) == 1L && n.chains > 1L
  sampler <- if (singleFit) object$fit[[1L]] else object$fit[[chainNum]]
  args <- list(
    treeNum = treeNum,
    chainNum = if (singleFit) chainNum else 1L,
    ...
  )
  if (!missing(sampleNum)) {
    args$sampleNum <- sampleNum
  }
  invisible(do.call(sampler$plotTree, args))
}

# plotTree.bart/.rbart read the trees off object$fit; a K-widened or two-part
# own-class fit has no single sampler that reads that way (bartHurdle has
# two), so each refuses by name instead of falling through to "no applicable
# method", pointing at the sampler(s) that do carry the trees.
refusePlotTreeMethod <- function(class, hint) {
  stop(
    "plotTree is defined for bart, rbart_vi and dbartsSampler fits; a ",
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
    "plotTree(object$occupancy$fit, ...) or plotTree(object$positive$fit, ...)"
  )
}

# survivalProbabilities.bart/.rbart dispatch on an aft or discrete-time hazard
# fit; none of the four own-class families is either, so each refuses by name
# instead of falling through to "no applicable method".
refuseSurvivalProbabilitiesMethod <- function(class) {
  stop(
    "survivalProbabilities applies to a discrete-time hazard fit ",
    "(bart2(family = \"hazard\")); a ",
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

  responseIsBinary <- is.null(object$sigma)
  sigma <- object$sigma
  if (!responseIsBinary && is.null(dim(sigma))) {
    sigma <- uncombineChains(as.vector(sigma), n.chains)
  }

  # the noise added below is always gaussian (rnorm); resid.dist is absent
  # for a fit predating the field or a binary fit and reads as gaussian, the
  # historical behavior; a present non-"gaussian" token (student residuals)
  # means that noise is wrong, so the draw is refused rather than taken
  if (!responseIsBinary) {
    residDist <- object[["resid.dist"]]
    if (!is.null(residDist) && !identical(residDist, "gaussian")) {
      stop(
        "posterior predictive sampling does not support ",
        residDist,
        " residuals"
      )
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
      n.draws <- length(sigma)
      noise <- rnorm(
        n.obs * n.draws,
        0,
        ppdNoiseScale(sigma, s, NULL, n.obs, n.draws)
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
      n.draws <- length(sigma)
      sd <- ppdNoiseScale(sigma, s, weights, n.obs, n.draws)
      noise <- rnorm(n.obs * n.draws, 0, sd)
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

# family/chain-count/tree-count/burn-in/kept-draws synopsis for print.bart
# and print.rbart, built only from fields that exist regardless of
# keepCall, keepSampler, and (for rbart) whether one sampler is kept per
# chain or a single in-core Gibbs sampler handles them all - so a fit
# created with keepCall = FALSE still prints something useful. n.trees and
# n.burn are only recoverable when the sampler itself was kept (keepTrees/
# keepSampler = TRUE); they are omitted otherwise, since the fit object
# does not retain them on its own.
fitSynopsis <- function(x) {
  fit <- x[["fit"]]
  fitIsList <- !is.null(fit) && is.list(fit)

  n.chains <- if (fitIsList) {
    # rbart: object$n.chains, when present, is authoritative even when the
    # kept sampler list has only one element - the in-core Gibbs path keeps
    # one multi-chain sampler and stores n.chains alongside it (matches the
    # n.chains idiom in predict.rbart/extract.rbart/plotTree.rbart)
    if (is.null(x[["n.chains"]])) length(fit) else x$n.chains
  } else if (!is.null(fit)) {
    fit$control@n.chains
  } else {
    x$n.chains
  }

  control <- if (!is.null(fit)) {
    if (fitIsList) fit[[1L]]$control else fit$control
  } else {
    NULL
  }

  varcountDims <- dim(x[["varcount"]])
  # a multi-forest fit's varcount carries a trailing forest margin (the
  # shapeMultinomialChannel shape: draws x p x n.forests), so its rank is one
  # higher throughout and the single-forest arms below would read the predictor
  # count as the draw count. n.forests, not the rank, is what separates the two
  # - a single-forest uncombined varcount is rank 3 as well. The arithmetic is
  # print.bartMultinomial's, which is why a multinomial fit never needed this
  n.forests <- x[["n.forests"]]
  n.kept <- if (!is.null(control)) {
    control@n.samples
  } else if (is.null(varcountDims)) {
    NA_integer_
  } else if (!is.null(n.forests) && n.forests > 1L) {
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
  invisible(NULL)
}

print.bart <- function(x, ...) {
  printCall(x)
  fitSynopsis(x)
  invisible(x)
}

print.rbart <- function(x, ...) {
  printCall(x)
  fitSynopsis(x)
  invisible(x)
}
