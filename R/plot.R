# Shared left-hand sigma-trace panel for plot.bart, plot.rbart and
# plot.bartHurdle: splits the device into a 1x2 layout and draws the
# residual-scale trace. A matrix-shaped 'sigma' (multiple chains) is drawn as
# one line per chain bridging the burn-in ('first.sigma', red) into the
# sampling run; a vector is a single scatter. The posterior-interval panel is
# drawn by each caller afterward. setLayout = FALSE skips the mfrow call for
# a caller (bartHurdle's 2x2) that has already set its own multi-panel
# layout, so this panel does not reset it.
plotSigmaTrace <- function(first.sigma, sigma, ..., setLayout = TRUE) {
  if (setLayout) {
    par(mfrow = c(1L, 2L))
  }
  if (!is.null(dim(sigma))) {
    plot(
      NULL,
      type = "n",
      ylab = "sigma",
      xlim = c(1, ncol(first.sigma) + ncol(sigma)),
      ylim = range(first.sigma, sigma)
    )
    for (i in seq_len(nrow(sigma))) {
      lines(
        c(seq_len(ncol(first.sigma)), ncol(first.sigma) + 0.5),
        c(
          first.sigma[i, ],
          0.5 * (first.sigma[i, ncol(first.sigma)] + sigma[i, 1L])
        ),
        col = "red",
        lty = i
      )
      lines(
        c(
          ncol(first.sigma) + 0.5,
          seq.int(ncol(first.sigma) + 1, length.out = ncol(sigma))
        ),
        c(
          0.5 * (first.sigma[i, ncol(first.sigma)] + sigma[i, 1L]),
          sigma[i, ]
        ),
        lty = i
      )
    }
  } else {
    plot(
      c(first.sigma, sigma),
      col = rep(c("red", "black"), c(length(first.sigma), length(sigma))),
      ylab = "sigma",
      ...
    )
  }
}

plot.bart <- function(
  x,
  plquants = c(0.05, 0.95),
  cols = c("blue", "black"),
  ...
) {
  if (is.null(x[["yhat.train"]])) {
    if (callName(x$call) == "bart2") {
      stop("plot requires bart2 to be called with 'keepTrainingFits' == TRUE")
    } else {
      stop("plot requires bart to be called with 'keeptrainfits' == TRUE")
    }
  }

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)

  if ("sigma" %in% names(x)) {
    par(mfrow = c(1L, 2L))
    plotSigmaTrace(x$first.sigma, x$sigma, ..., setLayout = FALSE)
  }

  if ("sigma" %in% names(x)) {
    ql <- apply(
      x$yhat.train,
      length(dim(x$yhat.train)),
      quantile,
      probs = plquants[1]
    )
    qm <- apply(x$yhat.train, length(dim(x$yhat.train)), quantile, probs = .5)
    qu <- apply(
      x$yhat.train,
      length(dim(x$yhat.train)),
      quantile,
      probs = plquants[2]
    )
    plot(
      x$y,
      qm,
      ylim = range(ql, qu),
      xlab = "y",
      ylab = "posterior interval for E(Y|x)",
      ...
    )
    # nolint next: seq_linter. 1:length preserves the (empty qm) edge behavior.
    for (i in 1:length(qm)) {
      lines(rep(x$y[i], 2), c(ql[i], qu[i]), col = cols[1])
    }
    abline(0, 1, lty = 2, col = cols[2])
  } else {
    pdrs <- probabilityFromLatents(x$yhat.train, x) #draws of p(Y=1 | x)
    ql <- apply(pdrs, length(dim(pdrs)), quantile, probs = plquants[1])
    qm <- apply(pdrs, length(dim(pdrs)), quantile, probs = .5)
    qu <- apply(pdrs, length(dim(pdrs)), quantile, probs = plquants[2])
    plot(
      qm,
      qm,
      ylim = range(ql, qu),
      xlab = "median of p",
      ylab = "posterior interval for P(Y=1|x)",
      ...
    )
    # nolint next: seq_linter. 1:length preserves the (empty qm) edge behavior.
    for (i in 1:length(qm)) {
      lines(rep(qm[i], 2), c(ql[i], qu[i]), col = cols[1])
    }
    abline(0, 1, lty = 2, col = cols[2])
  }
}

plot.rbart <- function(
  x,
  plquants = c(0.05, 0.95),
  cols = c("blue", "black"),
  ...
) {
  if (is.null(x[["yhat.train"]])) {
    stop("plot requires rbart_vi to be called with 'keepTrainingFits' == TRUE")
  }

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)

  if ("sigma" %in% names(x)) {
    par(mfrow = c(1L, 2L))
    plotSigmaTrace(x$first.sigma, x$sigma, ..., setLayout = FALSE)
  }

  if (length(dim(x$ranef)) > 2L) {
    ranef <- x$ranef[,, as.integer(x$group.by)]
  } else {
    ranef <- x$ranef[, as.integer(x$group.by)]
  }
  yhat.train <- x$yhat.train + ranef

  if ("sigma" %in% names(x)) {
    ql <- apply(
      yhat.train,
      length(dim(yhat.train)),
      quantile,
      probs = plquants[1L]
    )
    qm <- apply(yhat.train, length(dim(yhat.train)), quantile, probs = .5)
    qu <- apply(
      yhat.train,
      length(dim(yhat.train)),
      quantile,
      probs = plquants[2L]
    )
    plot(
      x$y,
      qm,
      ylim = range(ql, qu),
      xlab = "y",
      ylab = "posterior interval for E(Y | x)",
      ...
    )

    for (i in seq_along(qm)) {
      lines(rep(x$y[i], 2L), c(ql[i], qu[i]), col = cols[1L])
    }
    abline(0, 1, lty = 2L, col = cols[2L])
  } else {
    ## shouldn't happen for now
    pdrs <- probabilityFromLatents(yhat.train, x) #draws of p(Y=1 | x)
    ql <- apply(pdrs, length(dim(pdrs)), quantile, probs = plquants[1L])
    qm <- apply(pdrs, length(dim(pdrs)), quantile, probs = .5)
    qu <- apply(pdrs, length(dim(pdrs)), quantile, probs = plquants[2L])
    plot(
      qm,
      qm,
      ylim = range(ql, qu),
      xlab = "median of p",
      ylab = "posterior interval for P(Y = 1|  x)",
      ...
    )
    for (i in seq_along(qm)) {
      lines(rep(qm[i], 2L), c(ql[i], qu[i]), col = cols[1L])
    }
    abline(0, 1, lty = 2L, col = cols[2L])
  }
}

# A draws-array (any number of leading chain/sample margins, observations
# last) flattened to a plain (draws x observations) matrix: apply(x, MARGIN =
# last, ...) already pools every other margin per observation, but pooling
# AND subsetting observations generically (own-class plot below, over a
# y > 0 subset) needs the matrix form instead.
lastMarginMatrix <- function(x) {
  matrix(x, ncol = dim(x)[length(dim(x))])
}

# median + plquants interval per column of a (draws x n) matrix, the
# posterior-interval panel plot.bart's own gaussian/binary branches compute by
# hand; shared by the three own-class plot methods below.
drawInterval <- function(m, plquants) {
  list(
    med = apply(m, 2L, quantile, probs = 0.5),
    lo = apply(m, 2L, quantile, probs = plquants[1L]),
    hi = apply(m, 2L, quantile, probs = plquants[2L])
  )
}

# A per-category trace of the training-mean predicted probability, pooling
# chains back into one draw sequence: the closest cheap analog of plot.bart's
# sigma trace for this family, which has no residual scale. P2 is the binary
# panel plot.bart draws (median vs median, an interval bar per point) on the
# predicted probability of each observation's OWN observed category - for a
# multi-trial count response, that per-observation probability does not
# summarize the n x K cell structure, so P2 becomes plot.bart's gaussian
# panel instead: the observed proportion y_ik / n_i (a fixed number, not a
# draw) against the interval of the drawn p_ik, one point per (row,
# category) cell.
plot.bartMultinomial <- function(
  x,
  plquants = c(0.05, 0.95),
  cols = NULL,
  ...
) {
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)
  par(mfrow = c(1L, 2L))
  arr <- multinomialMeanProbArray(x)
  d <- dim(arr)
  trace <- matrix(arr, d[1L] * d[2L], d[3L])
  if (is.null(cols)) {
    cols <- seq_len(ncol(trace))
  }
  plot(
    NULL,
    type = "n",
    xlim = c(1L, nrow(trace)),
    ylim = range(trace),
    xlab = "iteration",
    ylab = "mean predicted probability"
  )
  for (k in seq_len(ncol(trace))) {
    lines(seq_len(nrow(trace)), trace[, k], col = cols[k])
  }
  legend("topright", legend = x$levels, col = cols, lty = 1L, bty = "n")

  y <- x$y
  probs <- x$yhat.train # (n.chains x) n.samples x n x K
  K <- x$K
  n <- length(y) %/% if (is.factor(y)) 1L else K
  if (is.matrix(y) && any(rowSums(y) > 1)) {
    flat <- probs
    dim(flat) <- c(length(probs) %/% (n * K), n * K)
    observed <- as.vector(y / rowSums(y))
    band <- drawInterval(flat, plquants)
    plot(
      observed,
      band$med,
      ylim = range(band$lo, band$hi),
      xlab = "observed proportion",
      ylab = "posterior interval for p",
      ...
    )
    for (j in seq_along(observed)) {
      lines(rep(observed[j], 2L), c(band$lo[j], band$hi[j]), col = "blue")
    }
  } else {
    category <- if (is.factor(y)) match(y, x$levels) else max.col(y, "first")
    flat <- probs
    dim(flat) <- c(length(probs) %/% (n * K), n, K)
    selected <- vapply(
      seq_len(n),
      function(i) flat[, i, category[i]],
      numeric(dim(flat)[1L])
    )
    band <- drawInterval(selected, plquants)
    plot(
      band$med,
      band$med,
      ylim = range(band$lo, band$hi),
      xlab = "median of p(observed category)",
      ylab = "posterior interval for p(observed category)",
      ...
    )
    for (i in seq_along(band$med)) {
      lines(rep(band$med[i], 2L), c(band$lo[i], band$hi[i]), col = "blue")
    }
  }
  abline(0, 1, lty = 2L)
  invisible(x)
}

# P1: one trace per FREE cutpoint (gamma_1 is pinned at 0 and carries no
# information); at K = 2 there is none, so the plot degrades to the single
# full-device latent panel plot.bart's binary branch draws. P2: the
# per-observation latent posterior interval, observations ordered by median
# eta, coloured by observed level, with dashed reference lines at the
# posterior-median cutpoints - this shows whether the latent separates the
# observed levels at the fitted thresholds, which a per-category probability
# panel (the multinomial shape) would not, since it drops the order.
plot.bartOrdinal <- function(x, plquants = c(0.05, 0.95), cols = NULL, ...) {
  K <- x$K
  cpArr <- ordinalCutpointsArray(x) # (iteration, chain, cutpoint); [, , 1] == 0
  cutpointMedians <- apply(cpArr, 3L, quantile, probs = 0.5)

  if (K > 2L) {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar), add = TRUE)
    par(mfrow = c(1L, 2L))
    cd <- dim(cpArr)
    trace <- matrix(cpArr, cd[1L] * cd[2L], cd[3L])[, -1L, drop = FALSE]
    if (is.null(cols)) {
      cols <- seq_len(ncol(trace))
    }
    plot(
      NULL,
      type = "n",
      xlim = c(1L, nrow(trace)),
      ylim = range(trace),
      xlab = "iteration",
      ylab = "cutpoint"
    )
    for (j in seq_len(ncol(trace))) {
      lines(seq_len(nrow(trace)), trace[, j], col = cols[j])
    }
    legend(
      "topright",
      legend = paste0("gamma[", seq.int(2L, K - 1L), "]"),
      col = cols,
      lty = 1L,
      bty = "n"
    )
  }

  m <- lastMarginMatrix(x$latent.train)
  band <- drawInterval(m, plquants)
  ord <- order(band$med)
  levelCols <- as.integer(x$y)[ord]
  plot(
    seq_along(band$med),
    band$med[ord],
    ylim = range(band$lo, band$hi),
    col = levelCols,
    xlab = "observation (ordered by median eta)",
    ylab = "posterior interval for the latent eta",
    ...
  )
  for (i in seq_along(band$med)) {
    lines(rep(i, 2L), c(band$lo[ord[i]], band$hi[ord[i]]), col = levelCols[i])
  }
  abline(h = cutpointMedians, lty = 2L)
  invisible(x)
}

# P1: the dispersion trace. There is no burn-in channel (bart2 negbin drives
# one run(n.burn, n.samples), so there is no first.dispersion to bridge
# from, unlike plot.bart's sigma panel), and r is drawn on an integer grid,
# so the trace is a step plot rather than a scatter. P2: plot.bart's
# gaussian panel verbatim on counts - observed y vs the posterior interval
# of the mean count.
plot.bartNegbin <- function(
  x,
  plquants = c(0.05, 0.95),
  cols = c("blue", "black"),
  ...
) {
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)
  par(mfrow = c(1L, 2L))
  disp <- x$dispersion
  if (is.null(dim(disp))) {
    plot(disp, type = "s", xlab = "iteration", ylab = "dispersion (r)")
  } else {
    plot(
      NULL,
      type = "n",
      xlim = c(1L, ncol(disp)),
      ylim = range(disp),
      xlab = "iteration",
      ylab = "dispersion (r)"
    )
    for (i in seq_len(nrow(disp))) {
      lines(disp[i, ], type = "s", lty = i)
    }
  }

  band <- drawInterval(lastMarginMatrix(x$yhat.train), plquants)
  plot(
    x$y,
    band$med,
    ylim = range(band$lo, band$hi),
    xlab = "y",
    ylab = "posterior interval for E(Y|x)",
    ...
  )
  for (i in seq_along(band$med)) {
    lines(rep(x$y[i], 2L), c(band$lo[i], band$hi[i]), col = cols[1L])
  }
  abline(0, 1, lty = 2L, col = cols[2L])
  invisible(x)
}

# Two component fits and a composed model need four panels. P1 reuses
# plotSigmaTrace verbatim (setLayout = FALSE: the 2x2 grid is already set):
# the positive part is an ordinary gaussian bart fit and carries both
# channels. P2 the occupancy probability, plot.bart's binary panel. P3 the
# positive part on the scale it actually fit (log y over the y > 0 rows). P4
# the composed natural-scale mean over ALL n rows (zeros included) - the only
# panel that shows the model this family exists for.
plot.bartHurdle <- function(
  x,
  plquants = c(0.05, 0.95),
  cols = c("blue", "black"),
  ...
) {
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)
  par(mfrow = c(2L, 2L))
  plotSigmaTrace(x$positive$first.sigma, x$positive$sigma, setLayout = FALSE)

  piBand <- drawInterval(
    lastMarginMatrix(extract(x$occupancy, type = "ev", sample = "train")),
    plquants
  )
  plot(
    piBand$med,
    piBand$med,
    ylim = range(piBand$lo, piBand$hi),
    xlab = "median of p",
    ylab = "posterior interval for P(Y > 0 | x)"
  )
  for (i in seq_along(piBand$med)) {
    lines(rep(piBand$med[i], 2L), c(piBand$lo[i], piBand$hi[i]), col = cols[1L])
  }
  abline(0, 1, lty = 2L, col = cols[2L])

  positiveRows <- x$y > 0
  fMat <- lastMarginMatrix(extract(x$positive, type = "bart", sample = "test"))
  fBand <- drawInterval(fMat[, positiveRows, drop = FALSE], plquants)
  logY <- log(x$y[positiveRows])
  plot(
    logY,
    fBand$med,
    ylim = range(fBand$lo, fBand$hi),
    xlab = "log(y), y > 0 rows",
    ylab = "posterior interval for f(x)"
  )
  for (i in seq_along(fBand$med)) {
    lines(rep(logY[i], 2L), c(fBand$lo[i], fBand$hi[i]), col = cols[1L])
  }
  abline(0, 1, lty = 2L, col = cols[2L])

  evBand <- drawInterval(lastMarginMatrix(extract(x, type = "ev")), plquants)
  plot(
    x$y,
    evBand$med,
    ylim = range(evBand$lo, evBand$hi),
    xlab = "y",
    ylab = "posterior interval for E(Y | x)",
    ...
  )
  for (i in seq_along(evBand$med)) {
    lines(rep(x$y[i], 2L), c(evBand$lo[i], evBand$hi[i]), col = cols[1L])
  }
  abline(0, 1, lty = 2L, col = cols[2L])
  invisible(x)
}

plot.pdbart <- function(
  x,
  xind = seq_along(x$fd),
  plquants = c(0.05, 0.95),
  cols = c("blue", "black"),
  ...
) {
  rgy <- range(x$fd)
  for (i in xind) {
    tsum <- apply(
      x$fd[[i]],
      2,
      quantile,
      probs = c(plquants[1], .5, plquants[2])
    )
    plot(
      range(x$levs[[i]]),
      rgy,
      type = "n",
      xlab = x$xlbs[i],
      ylab = "partial-dependence",
      ...
    )
    lines(x$levs[[i]], tsum[2, ], col = cols[1], type = "b")
    lines(x$levs[[i]], tsum[1, ], col = cols[2], type = "b")
    lines(x$levs[[i]], tsum[3, ], col = cols[2], type = "b")
  }
}

plot.pd2bart <- function(
  x,
  plquants = c(0.05, 0.95),
  contour.color = "white",
  justmedian = TRUE,
  ...
) {
  pdquants <- apply(x$fd, 2, quantile, probs = c(plquants[1], .5, plquants[2]))
  qq <- vector("list", 3)
  for (i in 1:3) {
    qq[[i]] <- matrix(pdquants[i, ], nrow = length(x$levs[[1]]))
  }
  if (justmedian) {
    zlim <- range(qq[[2]])
    vind <- c(2)
  } else {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar), add = TRUE)
    par(mfrow = c(1, 3))
    zlim <- range(qq)
    vind <- 1:3
  }
  for (i in vind) {
    image(
      x = x$levs[[1]],
      y = x$levs[[2]],
      qq[[i]],
      zlim = zlim,
      xlab = x$xlbs[1],
      ylab = x$xlbs[2],
      ...
    )
    contour(
      x = x$levs[[1]],
      y = x$levs[[2]],
      qq[[i]],
      zlim = zlim,
      ,
      add = TRUE,
      method = "edge",
      col = contour.color
    )
    title(main = c("Lower quantile", "Median", "Upper quantile")[i])
  }
}
