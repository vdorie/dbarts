# convergence diagnostics: posterior-package draws accessors and a
# summary() method for the scalar parameters of bart/bart2 fits.
# posterior is Suggests-only; as_draws_array/as_draws_df are registered
# for it conditionally in hooks.R. summary() itself is a base generic and
# degrades to a plain quantile table when posterior is unavailable.

# indirection point so tests can simulate posterior's absence without
# uninstalling it
posteriorAvailable <- function() requireNamespace("posterior", quietly = TRUE)

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
# varprobs, yhat.train, yhat.test, ..., and nbinom's 'dispersion') has
# the same (n.chains-combined-or-not) scalar shape as sigma - one draws
# variable per column/observation, named "field[inner]"
scalarFields <- c(
  "sigma",
  "k",
  "tau",
  "first.sigma",
  "first.k",
  "first.tau",
  "resid.df",
  "mean.s",
  "dispersion"
)

# One draws field by name: a stored channel, or the synthetic "mean.s" of a
# heteroscedastic fit - the mean of s(x) over the training observations, one
# value per draw, in sigma's own (n.chains, n.samples) scalar-field layout.
# The variance surface has no scalar to summarize (summarizing every
# observation's draws would swamp the table), so its convergence is read off
# that pooled mean, as bartMultinomial's is off meanProb.
drawsField <- function(object, v) {
  if (!identical(v, "mean.s")) {
    return(object[[v]])
  }
  s <- object[["s.train"]]
  if (is.null(s)) {
    return(NULL)
  }
  s <- combineOrUncombineChains(s, fitNChains(object), FALSE)
  apply(s, seq_len(length(dim(s)) - 1L), mean)
}

# "sigma" names the residual scale, which a heteroscedastic fit does not have
# as a scalar: the sigma it stores is the variance forest parameterization's
# fixed unit residual times the response range, constant across draws and
# carrying no posterior content. The token resolves to the fit's own scale
# channel there instead, so nothing reports that constant as a parameter.
resolveDrawsVars <- function(object, vars) {
  if (is.null(object[["s.train"]])) {
    return(vars)
  }
  replace(vars, vars == "sigma", "mean.s")
}

# the requested fields this fit actually carries, in the requested order
presentDrawsVars <- function(object, vars) {
  vars <- resolveDrawsVars(object, vars)
  vars[!vapply(vars, function(v) is.null(drawsField(object, v)), logical(1L))]
}

# Reshapes bart2(family = "ordinal")'s per-draw thresholds - the K - 1
# gamma_1 < ... < gamma_{K-1} - into the (iteration, chain,
# variable) convention bartDrawsArray's other fields share. Combined storage
# (the default) is (n.samples [* n.chains]) x (K - 1), the same layout
# toDrawsArray already gives any per-column field, so that path is reused
# directly; storage kept per-chain (combineChains = FALSE) is (K - 1) x
# n.samples x n.chains - a different axis order than toDrawsArray's own 3-D
# case assumes (which is built for a per-OBSERVATION field, chain-first), so
# it is permuted here instead of routed through it.
ordinalThresholdsArray <- function(object) {
  thresholds <- object$thresholds
  arr <- if (length(dim(thresholds)) == 3L) {
    aperm(thresholds, c(2L, 3L, 1L))
  } else {
    toDrawsArray(thresholds, object$n.chains, isScalar = FALSE)
  }
  dimnames(arr) <- list(
    NULL,
    NULL,
    paste0("threshold[", seq_len(dim(arr)[3L]), "]")
  )
  arr
}

# gathers one or more chain-dimensioned fields off a bart/bart2 fit
# into a single (iteration, chain, variable) base array. 'thresholds'
# (bartOrdinal only) is the one field whose shape toDrawsArray cannot read
# directly and so is special-cased to ordinalThresholdsArray, already named.
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
      return(ordinalThresholdsArray(object))
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

# ---- Rank-normalized split-Rhat and bulk/tail effective sample size, our
# own implementation of Vehtari, Gelman, Simpson, Carpenter, Burkner (2021,
# "Rank-normalization, folding, and localization"). Matched against the
# 'posterior' package's own internals (its exact constants and split/fold
# order, not merely the paper's prose) since summary()'s numbers must agree
# with a caller who separately has 'posterior' installed and runs it on the
# same array draws() returns.

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
  tauHat <- -1 + 2 * sum(rhoHatT[seq_len(maxT)]) + rhoHatT[maxT + 1L]
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
# estimator.
essQuantile <- function(x, prob) {
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
    data.frame(
      variable = varNames[i],
      mean = mean(pooled),
      median = stats::median(pooled),
      sd = stats::sd(pooled),
      mad = stats::mad(pooled),
      q5 = stats::quantile(pooled, 0.05, names = FALSE),
      q95 = stats::quantile(pooled, 0.95, names = FALSE),
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

as_draws_array.bart <- function(x, vars = c("sigma", "k", "tau"), ...) {
  posterior::as_draws_array(bartDrawsArray(x, vars))
}

as_draws_df.bart <- function(x, vars = c("sigma", "k", "tau"), ...) {
  posterior::as_draws_df(bartDrawsArray(x, vars))
}

# matches summary.bartOrdinal's own default 'vars'; bartDrawsArray already
# special-cases "thresholds" to ordinalThresholdsArray.
as_draws_array.bartOrdinal <- function(
  x,
  vars = c("thresholds", "sigma", "k", "tau"),
  ...
) {
  posterior::as_draws_array(bartDrawsArray(x, vars))
}
as_draws_df.bartOrdinal <- function(
  x,
  vars = c("thresholds", "sigma", "k", "tau"),
  ...
) {
  posterior::as_draws_df(bartDrawsArray(x, vars))
}

# matches summary.bartNegbin's own default 'vars'; scalarFields already lists
# "dispersion" as a sigma-shaped field.
as_draws_array.bartNegbin <- function(
  x,
  vars = c("dispersion", "sigma", "k", "tau"),
  ...
) {
  posterior::as_draws_array(bartDrawsArray(x, vars))
}
as_draws_df.bartNegbin <- function(
  x,
  vars = c("dispersion", "sigma", "k", "tau"),
  ...
) {
  posterior::as_draws_df(bartDrawsArray(x, vars))
}

# This family has no scalar parameter at all: the pooled per-category mean
# predicted probability is its only convergence instrument, and 'summary'
# already chose it (multinomialDrawsArray), so the two agree. 'vars' is
# accepted for signature symmetry with the other classes but has only ever
# one meaning here - never object$yhat.train's n x K per-observation draws.
as_draws_array.bartMultinomial <- function(x, vars = "meanProb", ...) {
  posterior::as_draws_array(multinomialDrawsArray(x))
}
as_draws_df.bartMultinomial <- function(x, vars = "meanProb", ...) {
  posterior::as_draws_df(multinomialDrawsArray(x))
}

# The union of both components' present scalar fields, each labelled with an
# "occupancy."/"positive." prefix (a dot, not a bracket: posterior parses a
# bracket as an index) - the same two blocks print.summary.bartHurdle already
# prints under, so as_draws and summary agree. Both components are driven by
# one n.chains/n.samples schedule and indexed draw for draw, so their
# (iteration, chain) margins match and the variable margins concatenate
# directly.
hurdleDrawsArray <- function(object, vars) {
  occ <- bartDrawsArray(object$occupancy, vars)
  pos <- bartDrawsArray(object$positive, vars)
  dimnames(occ)[[3L]] <- paste0("occupancy.", dimnames(occ)[[3L]])
  dimnames(pos)[[3L]] <- paste0("positive.", dimnames(pos)[[3L]])
  arr <- array(
    c(occ, pos),
    dim = c(dim(occ)[1:2], dim(occ)[3L] + dim(pos)[3L])
  )
  dimnames(arr) <- list(NULL, NULL, c(dimnames(occ)[[3L]], dimnames(pos)[[3L]]))
  arr
}

as_draws_array.bartHurdle <- function(x, vars = c("sigma", "k", "tau"), ...) {
  posterior::as_draws_array(hurdleDrawsArray(x, vars))
}
as_draws_df.bartHurdle <- function(x, vars = c("sigma", "k", "tau"), ...) {
  posterior::as_draws_df(hurdleDrawsArray(x, vars))
}

# plain per-variable mean/sd/quantiles, pooling samples and chains; the
# degrade path when posterior is not installed
quantileSummary <- function(arr) {
  pooled <- matrix(arr, prod(dim(arr)[1L:2L]), dim(arr)[3L])
  data.frame(
    variable = dimnames(arr)[[3L]],
    mean = colMeans(pooled),
    sd = apply(pooled, 2L, sd),
    q2.5 = apply(pooled, 2L, quantile, probs = 0.025, names = FALSE),
    median = apply(pooled, 2L, quantile, probs = 0.5, names = FALSE),
    q97.5 = apply(pooled, 2L, quantile, probs = 0.975, names = FALSE),
    row.names = NULL
  )
}

# rhat > 1.01 is noted in the printed summary, not enforced: dbarts does not
# refuse to summarize a non-converged fit
summary.bart <- function(object, vars = c("sigma", "k", "tau"), ...) {
  present <- presentDrawsVars(object, vars)
  havePosterior <- length(present) > 0L && posteriorAvailable()
  stats <- if (length(present) == 0L) {
    NULL
  } else if (havePosterior) {
    posterior::summarise_draws(posterior::as_draws_array(bartDrawsArray(
      object,
      present
    )))
  } else {
    quantileSummary(bartDrawsArray(object, present))
  }
  structure(
    list(
      call = object[["call"]],
      stats = stats,
      vars = vars,
      posterior = havePosterior
    ),
    class = "summary.bart"
  )
}

# bart2(family = "ordinal")'s scalar summary is the K - 1 thresholds, the only
# parameters this family's outer fit carries beyond whatever mean-function
# scale it shares with 'vars': neither sigma nor k/tau is tracked on the
# ordinal fit object, so the summary is thresholds alone; any that are later
# tracked would be picked up automatically through 'vars'.
summary.bartOrdinal <- function(
  object,
  vars = c("thresholds", "sigma", "k", "tau"),
  ...
) {
  summary.bart(object, vars = vars, ...)
}

# bart2(family = "nbinom")'s per-draw dispersion r rides its own 'dispersion'
# field, the count analog of gaussian's sigma; scalarFields already gives it
# sigma's shape, so this is summary.bart with a widened default 'vars'.
summary.bartNegbin <- function(
  object,
  vars = c("dispersion", "sigma", "k", "tau"),
  ...
) {
  summary.bart(object, vars = vars, ...)
}

# A hurdle fit is two ordinary bart2 fits under the
# hood - an occupancy probit on 1{y > 0} and a lognormal fit on the positive
# part - so each summarizes through summary.bart unchanged; only the
# packaging (both components, one call) and the print layout are new.
summary.bartHurdle <- function(object, vars = c("sigma", "k", "tau"), ...) {
  structure(
    list(
      call = object[["call"]],
      occupancy = summary.bart(object$occupancy, vars = vars, ...),
      positive = summary.bart(object$positive, vars = vars, ...)
    ),
    class = "summary.bartHurdle"
  )
}

# Collapses bartMultinomial's yhat.train (a (n.chains x) n.samples x n.obs x
# K probability array) over the observation margin into a per-category
# scalar channel shaped (iteration, chain, category) - the same
# (iteration, chain, variable) convention toDrawsArray produces for
# sigma/k/tau, built directly here since this family has no such scalar
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
# this family's only convergence instrument (no sigma/k/tau scale), shared by
# summary and as_draws so the two report the same channel under the same
# names.
multinomialDrawsArray <- function(object) {
  arr <- multinomialMeanProbArray(object)
  dimnames(arr)[[3L]] <- paste0("meanProb[", object$levels, "]")
  arr
}

multinomialSummaryVarsReason <- list(
  vars = paste0(
    "it pools the per-category mean-probability channel, which selects ",
    "nothing"
  )
)

# Convergence summary for a bart2(family = "multinomial") fit, mirroring
# summary.bart's shape (posterior mean/sd/quantiles, R-hat/ESS when
# 'posterior' is installed). This family has no sigma/k/tau scale to
# summarize, so the scalar channel is each category's posterior mean
# predicted probability, pooled over the training observations per draw -
# enough to eyeball per-category convergence without dumping every
# observation's draws. There is no other channel to pool instead, so 'vars'
# is refused by name rather than silently ignored.
summary.bartMultinomial <- function(object, ...) {
  refuseUnusedGenericArgs(
    list(...),
    "summary",
    "bartMultinomial",
    multinomialSummaryVarsReason
  )
  arr <- multinomialDrawsArray(object)
  havePosterior <- posteriorAvailable()
  stats <- if (havePosterior) {
    posterior::summarise_draws(posterior::as_draws_array(arr))
  } else {
    quantileSummary(arr)
  }
  structure(
    list(
      call = object[["call"]],
      stats = stats,
      vars = "meanProb",
      posterior = havePosterior
    ),
    class = "summary.bart"
  )
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
    return(invisible(NULL))
  }
  print(x$stats, ...)
  if (!x$posterior) {
    cat("\nNote: install 'posterior' for R-hat/ESS convergence diagnostics.\n")
  } else if (any(x$stats$rhat > 1.01, na.rm = TRUE)) {
    cat(
      "\nNote: some R-hat values exceed 1.01; chains may not have converged.\n"
    )
  }
  invisible(NULL)
}

print.summary.bart <- function(x, ...) {
  cat(
    "\nCall:\n",
    paste(deparse(x$call), sep = "\n", collapse = "\n"),
    "\n\n",
    sep = ""
  )
  printSummaryBartBody(x, ...)
  invisible(x)
}

print.summary.bartHurdle <- function(x, ...) {
  cat(
    "\nCall:\n",
    paste(deparse(x$call), sep = "\n", collapse = "\n"),
    "\n\n",
    sep = ""
  )
  cat("Occupancy component (probit, 1(y > 0)):\n")
  printSummaryBartBody(x$occupancy, ...)
  cat("\nPositive-part component (lognormal, y | y > 0):\n")
  printSummaryBartBody(x$positive, ...)
  invisible(x)
}
