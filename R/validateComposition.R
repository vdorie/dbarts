# dbartsValidateComposition: simulation-based calibration over a one-sweep
# step the HOST supplies, reporting per scalar functional whether the composed
# kernel targets the posterior it claims. The host owns the prior, so nothing
# here draws from a sampler's own prior the way the in-repo SBC harness does;
# what is shared with it is the rank rule and the uniformity verdict below.

# Rank of theta0 among L draws WHEN THE LAW HAS ATOMS. #{draws < theta0} is
# uniform only for an atomless law: an atom parks all its mass on one rank.
# Tagging every draw AND theta0 with an iid Uniform(0, 1) and ranking the pairs
# lexicographically makes theta0's tag exchangeable with the tags of the tied
# draws, so an atom contributes a Uniform{0, ..., #ties} increment and the
# total stays uniform on {0, ..., L}. Applied to EVERY functional: with no ties
# it is #{draws < theta0} exactly and consumes no random numbers, and a
# host-supplied functional - an indicator, a discrete parameter, a probability
# that underflows to zero - ties more readily than a continuous draw does.
sbcDiscreteRank <- function(draws, theta0) {
  below <- sum(draws < theta0)
  ties <- sum(draws == theta0)
  if (ties == 0L) {
    return(below)
  }
  tag0 <- runif(1L)
  below + sum(runif(ties) < tag0)
}

# The band and the KS jitter run off a FIXED stream, so a verdict is
# reproducible from the ranks alone. An exported function must not leave that
# fix behind: the caller's .Random.seed is restored on exit, and removed again
# when this call is what created it.
withFixedSeed <- function(seed, expr) {
  hasSeed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  oldSeed <- if (hasSeed) get(".Random.seed", envir = globalenv()) else NULL
  on.exit({
    if (hasSeed) {
      assign(".Random.seed", oldSeed, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  })
  set.seed(seed)
  expr
}

# Uniformity verdict for one functional's ranks. The headline is the
# ecdf-difference statistic against a simulation-based simultaneous band (Talts
# et al. 2018, fig. 1), which is already corrected for the multiple looks
# across the rank grid and so is the robust primary test. A chi-square goodness
# of fit on equal-width bins and a KS test against the discrete uniform (the
# jitter handles rank ties) are secondary signals: at 20 bins a lone small
# chi-square p across many functionals is within multiple-comparison noise, so
# the verdict does not hinge on it.
rankUniformity <- function(
  ranks,
  L,
  nBins = 20L,
  nSim = 2000L,
  alpha = 0.05,
  seed = 1L
) {
  R <- length(ranks)
  nBins <- min(nBins, L + 1L)
  edges <- seq(0, L + 1L, length.out = nBins + 1L)
  counts <- as.integer(table(cut(
    ranks,
    breaks = edges,
    include.lowest = TRUE,
    right = FALSE
  )))
  expected <- R / nBins
  chisqStat <- sum((counts - expected)^2 / expected)
  chisqP <- pchisq(chisqStat, df = nBins - 1L, lower.tail = FALSE)

  # the ecdf of integer ranks is a cumulated tabulation, which is what makes a
  # Bonferroni'd alpha affordable: the band is a 1 - alpha quantile, so a small
  # alpha needs more null draws to place it stably (>= 20 in the tail)
  target <- seq_len(L + 1L) / (L + 1)
  ecdfDiff <- function(rk) {
    cumsum(tabulate(rk + 1L, L + 1L)) / length(rk) - target
  }
  observed <- max(abs(ecdfDiff(ranks)))
  nSim <- max(nSim, ceiling(20 / alpha))

  withFixedSeed(seed, {
    # KS against the discrete uniform on {0, ..., L}: rank -> (rank + U)/(L + 1)
    u <- (ranks + runif(R)) / (L + 1)
    ksP <- suppressWarnings(ks.test(u, "punif")$p.value)
    nullMax <- numeric(nSim)
    for (s in seq_len(nSim)) {
      nullMax[s] <- max(abs(
        ecdfDiff(sample.int(L + 1L, R, replace = TRUE) - 1L)
      ))
    }
    band <- as.numeric(quantile(nullMax, 1 - alpha))
    list(
      counts = counts,
      nBins = nBins,
      chisqP = chisqP,
      ksP = ksP,
      ecdfDiff = observed,
      ecdfBand = band,
      pass = observed <= band,
      mean = mean(ranks),
      meanTarget = L / 2
    )
  })
}

compositionCount <- function(x, name, minimum) {
  whole <- is.numeric(x) && length(x) == 1L && is.finite(x) && x == round(x)
  if (!whole || x < minimum) {
    stop(sprintf("'%s' must be a single whole number >= %d", name, minimum))
  }
  as.integer(x)
}

# Every functional value is checked against the one the first replication's
# init produced. A length or name disagreement means the ranked quantity is not
# the drawn quantity, which SBC cannot detect for itself: the ranks would come
# out of whatever the recycling happened to align.
compositionFunctionals <- function(value, reference, where) {
  if (!is.numeric(value) || length(value) == 0L || !all(is.finite(value))) {
    stop(sprintf(
      "'functionals' must return a finite, non-empty numeric vector (%s)",
      where
    ))
  }
  agrees <- is.null(reference) ||
    (length(value) == length(reference) &&
      identical(names(value), names(reference)))
  if (!agrees) {
    stop(sprintf(
      "'functionals' returned %d value(s) %s and %d at the initial state: %s",
      length(value),
      where,
      length(reference),
      "the ranked functionals must be the same ones, in the same order"
    ))
  }
  value
}

dbartsValidateComposition <- function(
  drawPrior,
  simulate,
  init,
  step,
  functionals,
  n.replications = 200L,
  n.draws = 200L,
  n.thin = 30L,
  n.burn = 200L,
  alpha = 0.05,
  seed = NULL
) {
  matchedCall <- match.call()
  callbacks <- list(
    drawPrior = drawPrior,
    simulate = simulate,
    init = init,
    step = step,
    functionals = functionals
  )
  notFunctions <- !vapply(callbacks, is.function, logical(1L))
  if (any(notFunctions)) {
    stop(sprintf("'%s' must be a function", names(callbacks)[notFunctions][1L]))
  }
  n.replications <- compositionCount(n.replications, "n.replications", 2L)
  n.draws <- compositionCount(n.draws, "n.draws", 2L)
  n.thin <- compositionCount(n.thin, "n.thin", 1L)
  n.burn <- compositionCount(n.burn, "n.burn", 0L)
  scalar <- is.numeric(alpha) && length(alpha) == 1L && is.finite(alpha)
  if (!scalar || alpha <= 0 || alpha >= 1) {
    stop("'alpha' must be a single number in (0, 1)")
  }
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed)) {
      stop("'seed' must be a single number, or NULL to use the stream as it is")
    }
    set.seed(seed)
  }

  L <- n.draws
  ranks <- NULL
  draws <- NULL
  reference <- NULL
  for (r in seq_len(n.replications)) {
    theta0 <- drawPrior()
    # the ranked quantity is the functional AT theta0, evaluated on the state
    # init builds there - never a name match against theta0 itself, which
    # carries no name contract and cannot supply a derived functional at all
    state <- init(theta0, simulate(theta0))
    t0 <- compositionFunctionals(functionals(state), reference, "at init")
    if (is.null(reference)) {
      reference <- t0
      labels <- if (is.null(names(t0))) {
        paste0("f", seq_along(t0))
      } else {
        names(t0)
      }
      ranks <- matrix(NA_integer_, n.replications, length(t0))
      colnames(ranks) <- labels
      draws <- matrix(NA_real_, L, length(t0))
    }
    for (b in seq_len(n.burn)) {
      state <- step(state)
    }
    for (l in seq_len(L)) {
      for (t in seq_len(n.thin)) {
        state <- step(state)
      }
      draws[l, ] <- compositionFunctionals(
        functionals(state),
        reference,
        "after a step"
      )
    }
    for (j in seq_len(ncol(ranks))) {
      ranks[r, j] <- sbcDiscreteRank(draws[, j], t0[[j]])
    }
  }

  # Bonferroni over the functionals ranked in THIS call, so that a whole call
  # passes with probability ~1 - alpha rather than each functional alarming
  # independently at the nominal rate
  adjusted <- alpha / ncol(ranks)
  uniformity <- lapply(
    seq_len(ncol(ranks)),
    function(j) rankUniformity(ranks[, j], L, alpha = adjusted)
  )
  field <- function(name) {
    vapply(uniformity, function(u) u[[name]], numeric(1L))
  }
  passed <- vapply(uniformity, function(u) u$pass, logical(1L))
  verdicts <- data.frame(
    functional = colnames(ranks),
    mean.rank = field("mean"),
    ecdf.diff = field("ecdfDiff"),
    band = field("ecdfBand"),
    chisq.p = field("chisqP"),
    ks.p = field("ksP"),
    verdict = ifelse(passed, "PASS", "FLAG"),
    stringsAsFactors = FALSE
  )
  structure(
    list(
      call = matchedCall,
      ranks = ranks,
      L = L,
      n.replications = n.replications,
      n.thin = n.thin,
      n.burn = n.burn,
      alpha = alpha,
      alpha.adjusted = adjusted,
      mean.rank.target = L / 2,
      verdicts = verdicts,
      pass = all(passed)
    ),
    class = "dbartsCompositionValidation"
  )
}

print.dbartsCompositionValidation <- function(x, ...) {
  cat(
    "\nCall:\n",
    paste(deparse(x$call), sep = "\n", collapse = "\n"),
    "\n\n",
    sep = ""
  )
  cat(sprintf(
    "%d replications, %d draws each (thin %d, burn %d)\n",
    x$n.replications,
    x$L,
    x$n.thin,
    x$n.burn
  ))
  cat(sprintf(
    "band alpha %.3g over %d functional(s) is %.3g; uniform mean rank %.1f\n\n",
    x$alpha,
    nrow(x$verdicts),
    x$alpha.adjusted,
    x$mean.rank.target
  ))
  rows <- x$verdicts
  rows$mean.rank <- round(rows$mean.rank, 1)
  for (column in c("ecdf.diff", "band", "chisq.p", "ks.p")) {
    rows[[column]] <- signif(rows[[column]], 3L)
  }
  print(rows, row.names = FALSE)
  if (!x$pass) {
    cat(
      "\nFLAG: the flagged functionals' ranks are not uniform, so the composed",
      "\nkernel does not target the posterior 'drawPrior' and 'simulate' imply.",
      "\n"
    )
  }
  invisible(x)
}
