# Shared left-hand sigma-trace panel for plot.bart and plot.bartHurdle:
# splits the device into a 1x2 layout and draws the
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
    if (callName(x$call) == "bartBT") {
      stop("plot requires bartBT to be called with 'keeptrainfits' == TRUE")
    } else {
      stop(
        "plot requires bart to be called with 'keepTrainingFits' == TRUE ",
        "and 'keepFits' == TRUE (the latter set FALSE automatically when ",
        "'callback' is supplied, unless overridden)"
      )
    }
  }

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)

  hasResidual <- fitHasResidual(x)
  # a heteroscedastic fit has no scalar sigma to trace, nor has a fit that
  # held sigma fixed
  if (
    hasResidual && !fitIsHeteroscedastic(x) && is.null(x[["fixed"]][["sigma"]])
  ) {
    par(mfrow = c(1L, 2L))
    plotSigmaTrace(x$first.sigma, x$sigma, ..., setLayout = FALSE)
  }

  if (hasResidual) {
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
  # a count row with no trial has no observed category or proportion, so the
  # second panel covers only the rows with trials, in either branch
  withTrials <- if (is.factor(y)) seq_len(n) else which(rowSums(y) > 0)
  if (length(withTrials) == 0L) {
    plot.new()
    title(main = "no rows with trials to compare against")
    return(invisible(x))
  }
  if (is.matrix(y) && any(rowSums(y) > 1)) {
    flat <- probs
    dim(flat) <- c(length(probs) %/% (n * K), n, K)
    flat <- flat[, withTrials, , drop = FALSE]
    dim(flat) <- c(dim(flat)[1L], length(withTrials) * K)
    kept <- y[withTrials, , drop = FALSE]
    observed <- as.vector(kept / rowSums(kept))
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
      withTrials,
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

# P1: one trace per FREE threshold (gamma_1 is pinned at 0 and carries no
# information); at K = 2 there is none, so the plot degrades to the single
# full-device latent panel plot.bart's binary branch draws. P2: the
# per-observation latent posterior interval, observations ordered by median
# eta, coloured by observed level, with dashed reference lines at the
# posterior-median thresholds - this shows whether the latent separates the
# observed levels at the fitted thresholds, which a per-category probability
# panel (the multinomial shape) would not, since it drops the order.
plot.bartOrdinal <- function(x, plquants = c(0.05, 0.95), cols = NULL, ...) {
  K <- x$K
  cpArr <- ordinalThresholdsArray(x) # (iteration, chain, threshold); [, , 1] == 0
  thresholdMedians <- apply(cpArr, 3L, quantile, probs = 0.5)

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
      ylab = "threshold"
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
  abline(h = thresholdMedians, lty = 2L)
  invisible(x)
}

# P1: the shape trace, left out when the shape was held fixed. There is no
# burn-in channel (bart2 negbin drives one run(n.burn, n.samples), so there is
# no first.shape to bridge from, unlike plot.bart's sigma panel), and r is
# drawn on an integer grid, so the trace is a step plot rather than a scatter.
# P2: plot.bart's gaussian panel verbatim on counts - observed y vs the
# posterior interval of the mean count.
plot.bartNegbin <- function(
  x,
  plquants = c(0.05, 0.95),
  cols = c("blue", "black"),
  ...
) {
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)
  if (is.null(x[["fixed"]][["shape"]])) {
    par(mfrow = c(1L, 2L))
    disp <- x$shape
    if (is.null(dim(disp))) {
      plot(disp, type = "s", xlab = "iteration", ylab = "shape (r)")
    } else {
      plot(
        NULL,
        type = "n",
        xlim = c(1L, ncol(disp)),
        ylim = range(disp),
        xlab = "iteration",
        ylab = "shape (r)"
      )
      for (i in seq_len(nrow(disp))) {
        lines(disp[i, ], type = "s", lty = i)
      }
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
# plotSigmaTrace verbatim (setLayout = FALSE: the grid is already set):
# the positive part is an ordinary gaussian bart fit and carries both
# channels; it is left out, and the grid made 1x3, when that part held
# sigma fixed. P2 the zero part's probability, plot.bart's binary panel. P3 the
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
  if (is.null(x$positive[["fixed"]][["sigma"]])) {
    par(mfrow = c(2L, 2L))
    plotSigmaTrace(x$positive$first.sigma, x$positive$sigma, setLayout = FALSE)
  } else {
    par(mfrow = c(1L, 3L))
  }

  piBand <- drawInterval(
    lastMarginMatrix(extract(x$zero, type = "ev", sample = "train")),
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

## A factor predictor's levels are names; they plot at positions 1..K.
pdAxisPositions <- function(levs) {
  if (is.character(levs)) seq_along(levs) else levs
}

# A chains x draws x values array, as pdbart gives under combineChains =
# FALSE, as one draws x values matrix.
pdMergedDraws <- function(fd) {
  if (length(dim(fd)) == 3L) combineChains(fd) else fd
}

# The response's name as the fit's call wrote it: a formula's left-hand side,
# or the data argument when it is a name. NULL when neither reads as one.
pdResponseName <- function(call) {
  if (!is.call(call)) {
    return(NULL)
  }
  formula <- if ("x.train" %in% names(call)) call$x.train else call$formula
  if (is.call(formula) && identical(formula[[1L]], as.name("~"))) {
    return(if (length(formula) == 3L) deparse1(formula[[2L]]))
  }
  data <- if ("y.train" %in% names(call)) call$y.train else call$data
  if (is.name(data)) as.character(data)
}

# The scale a pdbart or pd2bart result is on, by its type and family; a
# result that records no type is labelled as it always was.
pdScaleLabel <- function(x) {
  type <- x$type
  family <- if (is.null(x$family)) "" else x$family
  if (is.null(type)) {
    return("partial-dependence")
  }
  response <- pdResponseName(x$bartcall)
  if (is.null(response)) {
    response <- "response"
  }
  binary <- family %in% c("probit", "logistic")
  hurdle <- family == "hurdle.lognormal"
  switch(
    type,
    bart = if (family == "probit") {
      "probit scale"
    } else if (family == "logistic") {
      "logit scale"
    } else if (family == "nbinom") {
      paste0("log mean ", response)
    } else if (hurdle) {
      paste0("log ", response, " where positive")
    } else {
      response
    },
    ev = if (binary) {
      "probability"
    } else if (hurdle) {
      "mean response"
    } else {
      response
    },
    prob = paste0("probability ", response, " is positive"),
    ppd = paste0("predicted ", response),
    sigma = paste0("residual sd of ", response),
    response
  )
}

plot.pdbart <- function(
  x,
  xind = seq_along(x$fd),
  plquants = c(0.05, 0.95),
  cols = c("blue", "black"),
  ...
) {
  # type, xlab and ylab are this method's own for the frame it draws, so a
  # caller's go to the lines and the axes rather than colliding with them
  dots <- list(...)
  lineType <- if (is.null(dots$type)) "b" else dots$type
  xlab <- dots$xlab
  ylab <- if (is.null(dots$ylab)) pdScaleLabel(x) else dots$ylab
  dots$type <- dots$xlab <- dots$ylab <- NULL

  fd <- lapply(x$fd, pdMergedDraws)
  rgy <- range(fd)
  for (i in xind) {
    tsum <- apply(
      fd[[i]],
      2,
      quantile,
      probs = c(plquants[1], .5, plquants[2])
    )
    levs <- x$levs[[i]]
    isFactor <- is.character(levs)
    at <- pdAxisPositions(levs)
    do.call(
      plot,
      c(
        list(
          if (isFactor) c(0.5, length(levs) + 0.5) else range(levs),
          rgy,
          type = "n",
          xlab = if (is.null(xlab)) x$xlbs[i] else xlab,
          ylab = ylab,
          xaxt = if (isFactor) "n" else "s"
        ),
        dots
      )
    )
    if (isFactor) {
      # one point per level, no line: the levels carry no order to join
      axis(1, at = at, labels = levs)
      segments(at, tsum[1, ], at, tsum[3, ], col = cols[2])
      points(at, tsum[2, ], col = cols[1], pch = 19)
    } else {
      lines(levs, tsum[2, ], col = cols[1], type = lineType)
      lines(levs, tsum[1, ], col = cols[2], type = lineType)
      lines(levs, tsum[3, ], col = cols[2], type = lineType)
    }
  }
}

plot.pd2bart <- function(
  x,
  plquants = c(0.05, 0.95),
  contour.color = "white",
  justmedian = TRUE,
  ...
) {
  pdquants <- apply(
    pdMergedDraws(x$fd),
    2,
    quantile,
    probs = c(plquants[1], .5, plquants[2])
  )
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
  isFactor <- vapply(x$levs[1:2], is.character, FALSE)
  at <- lapply(x$levs[1:2], pdAxisPositions)
  for (i in vind) {
    image(
      x = at[[1]],
      y = at[[2]],
      qq[[i]],
      zlim = zlim,
      xlab = x$xlbs[1],
      ylab = x$xlbs[2],
      xaxt = if (isFactor[1]) "n" else "s",
      yaxt = if (isFactor[2]) "n" else "s",
      ...
    )
    if (isFactor[1]) {
      axis(1, at = at[[1]], labels = x$levs[[1]])
    }
    if (isFactor[2]) {
      axis(2, at = at[[2]], labels = x$levs[[2]])
    }
    # contours interpolate between neighbours, which factor levels are not
    if (!any(isFactor)) {
      contour(
        x = at[[1]],
        y = at[[2]],
        qq[[i]],
        zlim = zlim,
        ,
        add = TRUE,
        method = "edge",
        col = contour.color
      )
    }
    # a caller's main went to image with the other arguments
    if (is.null(list(...)$main)) {
      title(
        main = paste0(
          c("Lower quantile", "Median", "Upper quantile")[i],
          ", ",
          pdScaleLabel(x)
        )
      )
    }
  }
}
