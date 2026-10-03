# convergence diagnostics: a summary() method for bart/bart2 fits reporting
# per-variable mean/median/sd/mad/quantiles plus split-Rhat and bulk/tail
# effective sample size, computed in-package (no 'posterior' dependency),
# built over a plain (iteration, chain, variable) array with dimnames.

# n.chains survives on the object whether or not the sampler was kept (see
# packageBartResults); fit is a single dbartsSampler
fitNChains <- function(object) {
  if (!is.null(object[["n.chains"]])) {
    return(object[["n.chains"]])
  }
  fit <- object[["fit"]]
  if (inherits(fit, "dbartsSampler")) fit$control@n.chains else length(fit)
}

# maps one field's native bart-convention samples - (n.chains, n.samples[,
# n.vars]), collapsed to drop the chain dimension when n.chains == 1, or
# flattened to a vector/matrix when combineChains was requested at fit time
# - to posterior's (iteration, chain, variable) array. uncombineChains
# already knows how to invert the flattening; only the trailing transpose
# is new here. isScalar disambiguates the two shapes a combined,
# multi-chain, 2-D field can have: a scalar field (sigma/k) stores
# uncombined as (n.chains, n.samples); a per-variable field (varcount,
# yhat.train, ...) stores COMBINED (the default) as
# (n.chains * n.samples, n.vars) - dim length 2 either way.
toDrawsArray <- function(x, n.chains, isScalar) {
  d <- dim(x)
  if (n.chains <= 1L) {
    if (is.null(d)) {
      arr <- array(x, c(length(x), 1L, 1L))
      varNames <- NULL
    } else {
      arr <- array(x, c(d[1L], 1L, d[2L]))
      varNames <- dimnames(x)[[2L]]
    }
  } else if (is.null(d)) {
    mat <- uncombineChains(x, n.chains) # n.chains x n.samples
    arr <- array(t(mat), c(ncol(mat), n.chains, 1L))
    varNames <- NULL
  } else if (length(d) == 2L && isScalar) {
    arr <- array(t(x), c(d[2L], d[1L], 1L))
    varNames <- NULL
  } else if (length(d) == 2L) {
    arr <- aperm(uncombineChains(x, n.chains), c(2L, 1L, 3L))
    varNames <- dimnames(x)[[2L]]
  } else {
    arr <- aperm(x, c(2L, 1L, 3L))
    varNames <- dimnames(x)[[3L]]
  }
  dimnames(arr) <- list(
    NULL,
    NULL,
    if (is.null(varNames)) as.character(seq_len(dim(arr)[3L])) else varNames
  )
  arr
}

# fields with no per-variable axis; every other requested field (varcount,
# varprobs, yhat.train, yhat.test, ..., and nbinom's 'shape') has
# the same (n.chains-combined-or-not) scalar shape as sigma - one draws
# variable per column/observation, named "field[inner]"
scalarFields <- c(
  "sigma",
  "k",
  "leaf.prior.sd",
  "first.sigma",
  "first.k",
  "resid.df",
  "mean.s",
  "shape"
)

# One draws field by name: a stored channel, or a synthetic one. "mean.s" is a
# heteroscedastic fit's mean of s(x) over the training observations, one
# value per draw, in sigma's own (n.chains, n.samples) scalar-field layout.
# The variance surface has no scalar to summarize (summarizing every
# observation's draws would swamp the table), so its convergence is read off
# that pooled mean, as bartMultinomial's is off its pooled per-category prob.
# "leaf.prior.sd" is the k.scale over a drawn k, in k's own layout; with k
# fixed there are no draws and the value is the fixed line's.
drawsField <- function(object, v) {
  if (identical(v, "leaf.prior.sd")) {
    prior <- object[["leaf.prior"]]
    if (is.null(object[["k"]]) || is.null(prior[["leaf.prior"]])) {
      return(NULL)
    }
    return(prior$k.scale / object[["k"]])
  }
  if (!identical(v, "mean.s")) {
    return(object[[v]])
  }
  s <- object[["s.train"]]
  if (is.null(s)) {
    return(NULL)
  }
  # the split (never combined) layout, collapsed at one chain exactly as
  # every other summary here reads its own scalar field: this table's rows
  # are draws, not chains, so dec-A79's kept chain margin is not for it
  n.chains <- fitNChains(object)
  if (n.chains > 1L && length(dim(s)) == 2L) {
    s <- uncombineChains(s, n.chains)
  }
  apply(s, seq_len(length(dim(s)) - 1L), mean)
}

# "sigma" names the residual scale, which a heteroscedastic fit does not have
# as a scalar, and so does not carry. The token resolves to the fit's own
# scale channel there instead.
resolveDrawsVars <- function(object, vars) {
  if (is.null(object[["s.train"]])) {
    return(vars)
  }
  replace(vars, vars == "sigma", "mean.s")
}

# the requested fields this fit actually carries, in the requested order
presentDrawsVars <- function(object, vars) {
  vars <- resolveDrawsVars(object, vars)
  vars <- vars[
    !vapply(vars, function(v) is.null(drawsField(object, v)), logical(1L))
  ]
  if (numSampledThresholds(object) == 0L) {
    vars <- setdiff(vars, "thresholds")
  }
  vars
}

# bart(family = "ordinal")'s per-draw thresholds, the K - 1 gamma_1 < ... <
# gamma_{K-1}, in the (iteration, chain, variable) convention, labelled
# threshold[j]. They are stored like any per-column field, (n.samples [*
# n.chains]) x (K - 1) combined or chains x n.samples x (K - 1) not. The
# first is pinned at zero; sampledOnly leaves it out, as the table of
# parameters does.
ordinalThresholdsArray <- function(object, sampledOnly = FALSE) {
  arr <- toDrawsArray(object$thresholds, object$n.chains, isScalar = FALSE)
  first <- if (sampledOnly) 2L else 1L
  arr <- arr[,, seq.int(first, dim(arr)[3L]), drop = FALSE]
  dimnames(arr) <- list(
    NULL,
    NULL,
    paste0("threshold[", seq.int(first, length.out = dim(arr)[3L]), "]")
  )
  arr
}

# the thresholds there are to tabulate: none when only the pinned one exists
numSampledThresholds <- function(object) {
  d <- dim(object[["thresholds"]])
  if (is.null(d)) 0L else d[length(d)] - 1L
}

# gathers one or more chain-dimensioned fields off a bart/bart2 fit
# into a single (iteration, chain, variable) base array. 'thresholds'
# (bartOrdinal only) is special-cased to ordinalThresholdsArray for its
# threshold[j] labels.
bartDrawsArray <- function(object, vars) {
  n.chains <- fitNChains(object)
  present <- presentDrawsVars(object, vars)
  if (length(present) == 0L) {
    stop(
      "none of 'vars' (",
      paste0(vars, collapse = ", "),
      ") are present on this fit",
      call. = FALSE
    )
  }
  pieces <- lapply(present, function(v) {
    if (identical(v, "thresholds")) {
      return(ordinalThresholdsArray(object, sampledOnly = TRUE))
    }
    piece <- toDrawsArray(drawsField(object, v), n.chains, v %in% scalarFields)
    dimnames(piece)[[3L]] <- if (v %in% scalarFields) {
      v
    } else {
      paste0(v, "[", dimnames(piece)[[3L]], "]")
    }
    piece
  })
  varNames <- unlist(lapply(pieces, function(p) dimnames(p)[[3L]]))
  array(
    unlist(pieces),
    dim = c(dim(pieces[[1L]])[1L], n.chains, length(varNames)),
    dimnames = list(NULL, NULL, varNames)
  )
}

# The union of both components' present scalar fields, each labelled with
# a "zero."/"positive." prefix (a dot, not a bracket) - the same two
# blocks print.summary.bartHurdle prints under. Both components are driven
# by one n.chains/n.samples schedule and indexed draw for draw, so their
# (iteration, chain) margins match and the variable margins concatenate
# directly.
hurdleDrawsArray <- function(object, vars) {
  zero <- bartDrawsArray(object$zero, vars)
  pos <- bartDrawsArray(object$positive, vars)
  dimnames(zero)[[3L]] <- paste0("zero.", dimnames(zero)[[3L]])
  dimnames(pos)[[3L]] <- paste0("positive.", dimnames(pos)[[3L]])
  arr <- array(
    c(zero, pos),
    dim = c(dim(zero)[1:2], dim(zero)[3L] + dim(pos)[3L])
  )
  dimnames(arr) <- list(
    NULL,
    NULL,
    c(dimnames(zero)[[3L]], dimnames(pos)[[3L]])
  )
  arr
}

# ---- Rank-normalized split-Rhat and bulk/tail effective sample size, our
# own implementation of Vehtari, Gelman, Simpson, Carpenter, Burkner (2021,
# "Rank-normalization, folding, and localization"). Matched against the
# 'posterior' package's own internals (its exact constants and split/fold
# order, not merely the paper's prose) since summary()'s numbers must agree
# with a caller who separately has 'posterior' installed and runs it on the
# same (iteration, chain, variable) array bartDrawsArray builds.

# Splits one (iteration, chain) matrix into 2 * ncol(x) half-chains, each
# floor(nrow(x) / 2) draws long - the middle draw of an odd-length chain is
# dropped rather than assigned to either half. Column order is [chain 1's
# first half, ..., chain M's first half, chain 1's second half, ..., chain
# M's second half].
splitChainsMatrix <- function(x) {
  x <- as.matrix(x)
  n <- nrow(x)
  if (n == 1L) {
    return(x)
  }
  half <- n / 2
  cbind(
    x[seq_len(floor(half)), , drop = FALSE],
    x[ceiling(half + 1):n, , drop = FALSE]
  )
}

# Rank-normalizes every element of x as ONE POOL (ties averaged): the van
# der Waerden transform qnorm((rank - 3/8) / (S - 3/4 + 1)), S = length(x).
# Blom's constant c = 3/8 gives the S - 3/4 + 1 denominator (equivalently
# S + 1/4) - not the S - 1/4 a literal reading of the paper's rounded prose
# might suggest.
rankNormalizeMatrix <- function(x) {
  r <- rank(x, ties.method = "average")
  s <- length(r)
  z <- stats::qnorm((r - 3 / 8) / (s - 3 / 4 + 1))
  # rank()'s default na.last = TRUE still assigns NA a (last) rank, so
  # without this z would read as finite there; put the NA back so a
  # non-finite input propagates to NA instead of a number.
  z[is.na(x)] <- NA
  dim(z) <- dim(x)
  z
}

# Folds the whole pooled variable around its median before any split or
# rank-normalization - the tail-Rhat/tail-ESS input.
foldDraws <- function(x) abs(x - stats::median(x))

# A constant (or NA/Inf-containing) input carries no information; every
# statistic below reports NA for it rather than dividing by a zero
# variance.
diagnosticsReturnNA <- function(x) {
  any(!is.finite(x)) || (abs(max(x) - min(x)) < .Machine$double.eps)
}

# Gelman-Rubin Rhat over an ALREADY split (and, for bulk, rank-normalized)
# matrix of half-chains: sqrt(((n - 1) / n * W + B / n) / W), W the mean
# within-half-chain variance, B = n * var(half-chain means).
gelmanRubinRhat <- function(x) {
  if (diagnosticsReturnNA(x)) {
    return(NA_real_)
  }
  n <- nrow(x)
  chainMeans <- colMeans(x)
  chainVars <- apply(x, 2L, stats::var)
  varBetween <- n * stats::var(chainMeans)
  varWithin <- mean(chainVars)
  sqrt((varBetween / varWithin + n - 1) / n)
}

# Per-half-chain autocovariance at every lag via FFT (Geyer 1992's trick):
# zero-pad past twice the next highly composite length, multiply the
# transform by its own conjugate (the power spectrum), inverse-transform
# back, and rescale so lag 0 reads the ordinary sample variance.
autocovariance <- function(x) {
  n <- length(x)
  varX <- stats::var(x)
  if (varX == 0) {
    return(rep(0, n))
  }
  m <- stats::nextn(n)
  yc <- c(x - mean(x), rep(0, 2L * m - n))
  ac <- Re(stats::fft(abs(stats::fft(yc))^2, inverse = TRUE)[seq_len(n)])
  ac / ac[1L] * varX * (n - 1) / n
}

# The shared ESS estimator (Stan's, via Geyer's initial monotone sequence)
# over an already split (and, for bulk, rank-normalized; for tail, a raw
# 0/1 indicator) matrix of half-chains.
essFromHalfChains <- function(x) {
  nChains <- ncol(x)
  n <- nrow(x)
  if (n < 3L || diagnosticsReturnNA(x)) {
    return(NA_real_)
  }
  acov <- apply(x, 2L, autocovariance)
  acovMeans <- rowMeans(acov)
  meanVar <- acovMeans[1L] * n / (n - 1)
  varPlus <- meanVar * (n - 1) / n
  if (nChains > 1L) {
    varPlus <- varPlus + stats::var(colMeans(x))
  }

  rhoHatT <- rep(0, n)
  rhoHatEven <- 1
  rhoHatT[1L] <- rhoHatEven
  rhoHatOdd <- 1 - (meanVar - acovMeans[2L]) / varPlus
  rhoHatT[2L] <- rhoHatOdd
  t <- 0L
  while (
    t < nrow(acov) - 5L &&
      !is.nan(rhoHatEven + rhoHatOdd) &&
      (rhoHatEven + rhoHatOdd > 0)
  ) {
    t <- t + 2L
    rhoHatEven <- 1 - (meanVar - acovMeans[t + 1L]) / varPlus
    rhoHatOdd <- 1 - (meanVar - acovMeans[t + 2L]) / varPlus
    if (rhoHatEven + rhoHatOdd >= 0) {
      rhoHatT[t + 1L] <- rhoHatEven
      rhoHatT[t + 2L] <- rhoHatOdd
    }
  }
  maxT <- t
  if (rhoHatEven > 0) {
    rhoHatT[maxT + 1L] <- rhoHatEven
  }

  # Geyer's initial monotone sequence: smooth consecutive pair sums so they
  # never increase.
  t <- 0L
  while (t <= maxT - 4L) {
    t <- t + 2L
    if (rhoHatT[t + 1L] + rhoHatT[t + 2L] > rhoHatT[t - 1L] + rhoHatT[t]) {
      rhoHatT[t + 1L] <- (rhoHatT[t - 1L] + rhoHatT[t]) / 2
      rhoHatT[t + 2L] <- rhoHatT[t + 1L]
    }
  }

  ess <- nChains * n
  # 1:maxT, not seq_len(maxT): posterior's own .ess indexes rho_hat_t[1:max_t]
  # literally, and at max_t == 0 that 1:0 still selects element 1 (R drops
  # the 0), giving tau_hat = 2 - not the empty sum seq_len(0) would give,
  # which uncaps tau_hat toward 0 and inflates ESS by roughly log10(ess).
  tauHat <- -1 + 2 * sum(rhoHatT[1:maxT]) + rhoHatT[maxT + 1L]
  tauBound <- 1 / log10(ess)
  if (tauHat < tauBound) {
    tauHat <- tauBound
  }
  ess / tauHat
}

# Bulk Rhat: rank-normalize the pooled split-chain draws, then
# Gelman-Rubin. Tail (folded) Rhat: fold the raw pooled draws FIRST, split
# THAT, THEN rank-normalize; report the max of the two.
splitRhat <- function(x) {
  bulk <- gelmanRubinRhat(rankNormalizeMatrix(splitChainsMatrix(x)))
  tail <- gelmanRubinRhat(rankNormalizeMatrix(splitChainsMatrix(foldDraws(x))))
  max(bulk, tail)
}

# Bulk ESS: the same split-then-rank-normalize array Bulk Rhat uses.
essBulk <- function(x) {
  essFromHalfChains(rankNormalizeMatrix(splitChainsMatrix(x)))
}

# One quantile's ESS: an indicator on the RAW, UNSPLIT, UNRANKED pooled
# draws, split (never rank-normalized) and passed straight to the shared
# estimator. The NA/Inf/constant guard runs on the RAW draws, before
# stats::quantile ever sees them - an NA/NaN draw would otherwise error
# out of quantile() instead of propagating NA, and an Inf draw would
# still produce a (meaningless) finite indicator.
essQuantile <- function(x, prob) {
  if (diagnosticsReturnNA(x)) {
    return(NA_real_)
  }
  indicator <- x <= stats::quantile(x, probs = prob, names = FALSE)
  essFromHalfChains(splitChainsMatrix(indicator))
}

# Tail ESS: the smaller of the 5% and 95% quantile ESS values.
essTail <- function(x) min(essQuantile(x, 0.05), essQuantile(x, 0.95))

# The nine-column per-variable summary: mean/median/sd/mad/q5/q95 (ordinary
# pooled statistics, matching R's own mean/median/sd/mad/quantile) plus
# rhat/ess_bulk/ess_tail above - the columns 'posterior::summarise_draws'
# reports for a plain array, computed without it. Always run, unconditional
# on any package's availability (Decision 1).
summariseDraws <- function(arr) {
  varNames <- dimnames(arr)[[3L]]
  d <- dim(arr)[1:2]
  rows <- lapply(seq_along(varNames), function(i) {
    x <- arr[,, i, drop = TRUE]
    dim(x) <- d
    pooled <- as.vector(x)
    # posterior's own quantile2.default checks anyNA(x) before calling
    # quantile() and reports NA instead when it holds - quantile() itself
    # errors on an NA/NaN input rather than propagating NA the way
    # mean/median/sd/mad already do.
    qs <- if (anyNA(pooled)) {
      c(NA_real_, NA_real_)
    } else {
      stats::quantile(pooled, c(0.05, 0.95), names = FALSE)
    }
    data.frame(
      variable = varNames[i],
      mean = mean(pooled),
      median = stats::median(pooled),
      sd = stats::sd(pooled),
      mad = stats::mad(pooled),
      q5 = qs[1L],
      q95 = qs[2L],
      rhat = splitRhat(x),
      ess_bulk = essBulk(x),
      ess_tail = essTail(x),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# What summary names as fixed under the table instead of tabulating: each
# requested parameter the fit held fixed, by its label, with its value - one
# number, or one per forest - and the ordinal's first threshold, pinned at
# zero. A fit saved before fits recorded what they held fixed names none, and
# its constant channels tabulate as they always did.
fixedSummaryValues <- function(object, vars) {
  fixed <- object[["fixed"]]
  values <- list()
  for (v in vars) {
    if (v == "thresholds") {
      if (!is.null(object[["thresholds"]])) {
        values[["threshold[1]"]] <- object[["thresholds"]][1L]
      }
      next
    }
    value <- if (v == "leaf.prior.sd") {
      if (!is.null(fixed[["k"]])) extractParameter(object, v, TRUE)
    } else {
      fixed[[v]]
    }
    if (!is.null(value)) {
      values[[v]] <- value
    }
  }
  values
}

# the leaf-scale quantity a fit is reported in follows how it named its leaf
# prior: leaf.prior.sd for an sd, and k otherwise. A fit that records no
# naming keeps both.
defaultLeafScaleVars <- function(object, vars) {
  prior <- object[["leaf.prior"]]
  if (is.null(prior)) {
    return(vars)
  }
  if (is.null(prior[["leaf.prior"]])) {
    prior <- prior[[1L]]
  }
  named <- if (
    !is.null(prior[["basis.row.norm"]]) || !is.null(prior$leaf.prior@prior.sd)
  ) {
    "leaf.prior.sd"
  } else {
    "k"
  }
  setdiff(vars, setdiff(c("k", "leaf.prior.sd"), named))
}

# rhat > 1.01 is noted in the printed summary, not enforced: dbarts does not
# refuse to summarize a non-converged fit. Parameters the fit held fixed are
# named under the table, not tabulated as constants.
summary.bart <- function(
  object,
  vars = c("sigma", "k", "leaf.prior.sd", "resid.df"),
  ...
) {
  if (missing(vars)) {
    vars <- defaultLeafScaleVars(object, vars)
  }
  fixed <- fixedSummaryValues(object, vars)
  present <- presentDrawsVars(object, setdiff(vars, names(fixed)))
  stats <- if (length(present) == 0L) {
    NULL
  } else {
    summariseDraws(bartDrawsArray(object, present))
  }
  result <- list(
    call = object[["call"]],
    stats = stats,
    vars = vars,
    fixed = fixed
  )
  # the monotone prior, absent on a fit without a constraint
  result$monotone.prior <- object[["monotone.prior"]]
  structure(result, class = "summary.bart")
}

# bart2(family = "ordinal")'s scalar summary is the sampled thresholds beside
# whatever leaf scale the fit draws; its first threshold is pinned and named
# under the table, and sigma is not a parameter of the family.
summary.bartOrdinal <- function(
  object,
  vars = c("thresholds", "sigma", "k", "leaf.prior.sd"),
  ...
) {
  if (missing(vars)) {
    vars <- defaultLeafScaleVars(object, vars)
  }
  summary.bart(object, vars = vars, ...)
}

# bart2(family = "nbinom")'s per-draw shape r rides its own 'shape'
# field, the count analog of gaussian's sigma; scalarFields already gives it
# sigma's shape, so this is summary.bart with a widened default 'vars'.
summary.bartNegbin <- function(
  object,
  vars = c("shape", "sigma", "k", "leaf.prior.sd"),
  ...
) {
  if (missing(vars)) {
    vars <- defaultLeafScaleVars(object, vars)
  }
  summary.bart(object, vars = vars, ...)
}

# A hurdle fit is two ordinary bart2 fits under the
# hood - a zero-part probit on 1{y > 0} and a lognormal fit on the positive
# part - so each summarizes through summary.bart unchanged; only the
# packaging (both components, one call) and the print layout are new.
summary.bartHurdle <- function(
  object,
  vars = c("sigma", "k", "leaf.prior.sd", "resid.df"),
  ...
) {
  defaulted <- missing(vars)
  partSummary <- function(part) {
    if (defaulted) {
      summary.bart(part, ...)
    } else {
      summary.bart(part, vars = vars, ...)
    }
  }
  structure(
    list(
      call = object[["call"]],
      zero = partSummary(object$zero),
      positive = partSummary(object$positive)
    ),
    class = "summary.bartHurdle"
  )
}

# Collapses bartMultinomial's yhat.train (a (n.chains x) n.samples x n.obs x
# K probability array) over the observation margin into a per-category
# scalar channel shaped (iteration, chain, category) - the same
# (iteration, chain, variable) convention toDrawsArray produces for
# sigma/k, built directly here since this family has no such scalar
# field to reuse and its K-widened varcount/yhat shapes do not match the
# non-multinomial dims toDrawsArray assumes.
multinomialMeanProbArray <- function(object) {
  probs <- object$yhat.train
  d <- dim(probs)
  numDims <- length(d)
  obsMargin <- numDims - 1L
  means <- apply(probs, seq_len(numDims)[-obsMargin], mean)
  n.chains <- object$n.chains
  arr <- if (length(dim(means)) == 3L) {
    # already (chains, samples, K); reorder to (iteration, chain, K)
    aperm(means, c(2L, 1L, 3L))
  } else if (n.chains <= 1L) {
    array(means, c(dim(means)[1L], 1L, dim(means)[2L]))
  } else {
    # combineChains folds samples fastest within each chain (see
    # shapeMultinomialChannel), so splitting the leading margin back into
    # (samples, chains) in that order recovers the original layout
    array(means, c(dim(means)[1L] %/% n.chains, n.chains, dim(means)[2L]))
  }
  dimnames(arr) <- list(NULL, NULL, object$levels)
  arr
}

# multinomialMeanProbArray with its categories named on the variable margin -
# this family's only convergence instrument (no sigma/k scale), shared by
# summary alone so the label lives in one place.
multinomialDrawsArray <- function(object) {
  arr <- multinomialMeanProbArray(object)
  dimnames(arr)[[3L]] <- paste0("prob[", object$levels, "]")
  arr
}

multinomialSummaryVarsReason <- list(
  vars = paste0(
    "it pools the per-category pooled probability channel, which selects ",
    "nothing"
  )
)

# Convergence summary for a bart2(family = "multinomial") fit, mirroring
# summary.bart's shape (mean/sd/quantiles plus R-hat/ESS, unconditionally).
# This family has no sigma/k scale to summarize, so the scalar channel
# is each category's posterior mean predicted probability, pooled over the
# training observations per draw - enough to eyeball per-category
# convergence without dumping every observation's draws. There is no other
# channel to pool instead, so 'vars' is refused by name rather than
# silently ignored.
summary.bartMultinomial <- function(object, ...) {
  refuseUnusedGenericArgs(
    list(...),
    "summary",
    "bartMultinomial",
    multinomialSummaryVarsReason
  )
  arr <- multinomialDrawsArray(object)
  structure(
    list(
      call = object[["call"]],
      stats = summariseDraws(arr),
      vars = "prob",
      fixed = fixedSummaryValues(
        object,
        defaultLeafScaleVars(object, c("k", "leaf.prior.sd"))
      )
    ),
    class = "summary.bart"
  )
}

# one line under the table naming each parameter held fixed and its value, in
# summary.glm's manner; a value per forest is labelled by forest
printFixedLine <- function(fixed) {
  if (length(fixed) == 0L) {
    return(invisible(NULL))
  }
  show <- function(value) toString(format(unname(value), digits = 4L))
  labels <- unlist(Map(
    function(name, value) {
      if (is.matrix(value)) {
        paste0(name, "[", rownames(value), "] = ", apply(value, 1L, show))
      } else if (!is.null(names(value))) {
        paste0(
          name,
          "[",
          names(value),
          "] = ",
          format(unname(value), digits = 4L)
        )
      } else {
        paste0(name, " = ", show(value))
      }
    },
    names(fixed),
    fixed
  ))
  cat("(Fixed, not sampled: ", paste0(labels, collapse = ", "), ")\n", sep = "")
  invisible(NULL)
}

# the row-table body of a summary.bart object, with no Call: header - shared
# by print.summary.bart and print.summary.bartHurdle, which prints one header
# for the fit and this body once per component
printSummaryBartBody <- function(x, ...) {
  if (is.null(x$stats)) {
    # the requested set, not a fixed one: each family summarizes its own
    cat(
      "No scalar parameters (",
      paste0(x$vars, collapse = ", "),
      ") to summarize.\n",
      sep = ""
    )
    printFixedLine(x$fixed)
    return(invisible(NULL))
  }
  print(x$stats, ...)
  printFixedLine(x$fixed)
  if (any(x$stats$rhat > 1.01, na.rm = TRUE)) {
    cat(
      "\nNote: some R-hat values exceed 1.01; chains may not have converged.\n"
    )
  }
  invisible(NULL)
}

print.summary.bart <- function(x, ...) {
  printCall(x)
  if (!is.null(x$monotone.prior)) {
    cat("Monotone prior: ", x$monotone.prior, "\n\n", sep = "")
  }
  printSummaryBartBody(x, ...)
  invisible(x)
}

print.summary.bartHurdle <- function(x, ...) {
  printCall(x)
  cat("Zero-part component (probit, 1(y > 0)):\n")
  printSummaryBartBody(x$zero, ...)
  cat("\nPositive-part component (lognormal, y | y > 0):\n")
  printSummaryBartBody(x$positive, ...)
  invisible(x)
}
