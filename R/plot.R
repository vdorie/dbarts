# Shared left-hand sigma-trace panel for plot.bart and plot.rbart: splits the
# device into a 1x2 layout and draws the residual-scale trace. A matrix-shaped
# 'sigma' (multiple chains) is drawn as one line per chain bridging the burn-in
# ('first.sigma', red) into the sampling run; a vector is a single scatter. The
# posterior-interval panel is drawn by each caller afterward.
plotSigmaTrace <- function(first.sigma, sigma, ...) {
  par(mfrow = c(1L, 2L))
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

  if ("sigma" %in% names(x)) {
    plotSigmaTrace(x$first.sigma, x$sigma, ...)
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

  if ("sigma" %in% names(x)) {
    plotSigmaTrace(x$first.sigma, x$sigma, ...)
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

# A per-category trace of the training-mean predicted probability, pooling
# chains back into one draw sequence: the closest cheap analog of plot.bart's
# sigma trace for this family, which has no residual scale and a categorical
# y rather than a single scalar per observation to plot against a posterior
# interval.
plot.bartMultinomial <- function(x, cols = NULL, ...) {
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
    ylab = "mean predicted probability",
    ...
  )
  for (k in seq_len(ncol(trace))) {
    lines(seq_len(nrow(trace)), trace[, k], col = cols[k])
  }
  legend("topright", legend = x$levels, col = cols, lty = 1L, bty = "n")
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
