#!/usr/bin/env Rscript

# The same seeded fits on two builds, the base and the one that takes a
# prior's scale from the rows in the likelihood, compared as the slice says
# they must compare. Run `record` once on each build with the same seeds,
# then `compare` on the two files:
#
#   R_LIBS=<base lib>  Rscript benchmarks/R/rows-pair.R record base.rds
#   R_LIBS=<slice lib> Rscript benchmarks/R/rows-pair.R record slice.rds
#   Rscript benchmarks/R/rows-pair.R compare base.rds slice.rds
#
# Each scenario is recorded as stated ("asis") and, where rows are out with
# an extreme response, with the rows out given responses inside the range of
# the rows in ("moved"). The groups:
#   1 no row out (every family and composition, a mask after creation,
#     updateScale = TRUE, xbart): slice asis identical to base asis;
#   2 rows out with no extreme among them: identical;
#   3 gaussian and Student-t with an extreme on a row out: the slice's fit
#     asis equals, at the rows in, the base's fit moved;
#   4 a count fit masked at structural zeros, as many out as in, no offset:
#     the slice's fit asis equals, at the rows in, the base's with the masked
#     rows given the counts of the rows in;
#   5 several forests with rows out, hazard with an out subject off the
#     grid, and each other mover of the posterior: must differ, reported by
#     how much.
# Exit status 1 if a group 1 to 4 pair is not identical or a group 5 pair is.

suppressPackageStartupMessages(library(dbarts))

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) >= 1L) args[1L] else ""
nSweeps <- 20L

fitBart <- function(x, y, ..., weights = NULL, rows = NULL) {
  call <- list(
    x,
    y,
    n.trees = 10L,
    n.samples = nSweeps,
    n.burn = 10L,
    n.chains = 1L,
    n.threads = 1L,
    seed = 11L,
    verbose = FALSE,
    ...
  )
  if (!is.null(weights)) {
    call$weights <- weights
  }
  fit <- suppressWarnings(do.call(bart, call))
  train <- fit$yhat.train
  list(
    train = if (is.null(rows)) train else train[, rows, drop = FALSE],
    sigma = fit$sigma,
    periods = fit$periods
  )
}
fitSampler <- function(x, y, ..., after = NULL, rows = NULL) {
  control <- dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 10L,
    updateState = FALSE,
    seed = 11L
  )
  sampler <- suppressWarnings(dbarts(x, y, control = control, ...))
  if (!is.null(after)) {
    suppressWarnings(after(sampler))
  }
  out <- sampler$run(0L, nSweeps)
  list(
    train = if (is.null(rows)) out$train else out$train[rows, , drop = FALSE],
    sigma = out$sigma
  )
}

n <- 120L
set.seed(20261013L)
x <- matrix(runif(n * 4L), n, 4L)
mu <- 3 * x[, 1L] - 2 * x[, 2L]
y <- mu + rnorm(n, 0, 0.5)
rowsOut <- seq_len(n) %% 2L == 0L
rowsIn <- !rowsOut
wMask <- as.double(rowsIn)
yExtreme <- y
yExtreme[rowsOut] <- y[rowsOut] + 50
# the same response with the rows out inside the range of the rows in
yInside <- y
yInside[rowsOut] <- median(y[rowsIn])
yInside[rowsOut][1L] <- min(y[rowsIn])
yInside[rowsOut][2L] <- max(y[rowsIn])
wPos <- 0.5 + runif(n)
z <- as.double(seq_len(n) %% 2L)
latent <- as.vector(scale(mu)) + rnorm(n)
yBinary <- as.double(latent > 0)
yOrdinal <- cut(
  latent,
  c(-Inf, -0.8, 0.2, 1, Inf),
  labels = c("a", "b", "c", "d"),
  ordered_result = TRUE
)
yMulti <- factor(sample(c("a", "b", "c"), n, replace = TRUE))
counts <- as.double(rpois(n, exp(0.5 + 0.5 * x[, 1L])))
countsMasked <- counts
countsMasked[rowsOut] <- 0
# the masked rows given the counts of the rows in
countsInside <- counts
countsInside[rowsOut] <- counts[rowsIn]
timeS <- exp(0.3 * mu + rnorm(n, 0, 0.5))
statusS <- as.double(runif(n) < 0.7)
nh <- 40L
xh <- x[seq_len(nh), ]
timeH <- rep(1:5, length.out = nh)
statusH <- as.double(runif(nh) < 0.7)
wH <- rep(1, nh)
wH[1:2] <- 0
timeHIn <- timeH
timeHIn[1:2] <- c(3, 5)
timeHOff <- timeH
timeHOff[1:2] <- c(3.5, 20)
fitHazard <- function(time, w = NULL, family = "hazard") {
  fitBart(xh, cbind(time, statusH), family = family, weights = w)
}
twoForests <- function(response, weights = NULL) {
  fitSampler(
    x,
    response,
    weights = weights,
    forests = list(forest(), forest(basis = ~ factor(z)))
  )
}
xbartRun <- function(response, weights = NULL) {
  list(
    train = xbart(
      x,
      response,
      weights = weights,
      n.samples = 10L,
      n.burn = 5L,
      n.trees = 5L,
      n.reps = 2L,
      n.test = 4,
      n.threads = 1L,
      seed = 7L,
      verbose = FALSE
    )
  )
}
student <- dbarts:::student(df = 4)

# name -> list(group, asis = function(), moved = function() or NULL, rows)
scenarios <- list(
  gaussian = list(1L, function() fitBart(x, y)),
  weighted = list(1L, function() fitBart(x, y, weights = wPos)),
  student = list(1L, function() fitSampler(x, y, family = student)),
  variance = list(
    1L,
    function() fitSampler(x, y, variance = varianceForest(n.trees = 10L))
  ),
  gp = list(
    1L,
    function() fitSampler(data.frame(x), y, leaf.prior = gp("X1", k = 2))
  ),
  linear = list(
    1L,
    function() {
      fitSampler(
        data.frame(x),
        y,
        leaf.prior = dbartsPriors$linear("X2")
      )
    }
  ),
  twoforests = list(1L, function() twoForests(y)),
  aft = list(1L, function() fitBart(x, cbind(timeS, statusS), family = "aft")),
  probit = list(1L, function() fitBart(x, yBinary)),
  probitmask = list(
    1L,
    function() {
      fitSampler(x, yBinary, family = "probit", after = function(s) {
        s$setActiveRows(wMask)
      })
    }
  ),
  logistic = list(1L, function() fitBart(x, yBinary, family = "logistic")),
  ordinal = list(1L, function() fitBart(x, yOrdinal, family = "ordinal")),
  multinomial = list(1L, function() fitBart(x, yMulti, family = "multinomial")),
  nbinom = list(1L, function() fitSampler(x, counts, family = "nbinom")),
  nbinommask = list(
    1L,
    function() {
      fitSampler(x, counts, family = "nbinom", after = function(s) {
        s$setActiveRows(wMask)
      })
    }
  ),
  hazard = list(1L, function() fitHazard(timeH)),
  maskaftercreation = list(
    1L,
    function() fitSampler(x, y, after = function(s) s$setActiveRows(wMask))
  ),
  updatescale = list(
    1L,
    function() {
      fitSampler(x, y, after = function(s) {
        s$setOffset(rep(0, n), updateScale = TRUE)
        s$setResponse(y, updateScale = TRUE)
      })
    }
  ),
  xbart = list(1L, function() xbartRun(y)),
  rowsoutinside = list(
    2L,
    function() fitSampler(x, yInside, weights = wMask)
  ),
  hazardoutonthegrid = list(2L, function() fitHazard(timeHIn, wH)),
  gaussianextreme = list(
    3L,
    function() fitSampler(x, yExtreme, weights = wMask, rows = which(rowsIn)),
    function() fitSampler(x, yInside, weights = wMask, rows = which(rowsIn))
  ),
  studentextreme = list(
    3L,
    function() {
      fitSampler(
        x,
        yExtreme,
        weights = wMask,
        family = student,
        rows = which(rowsIn)
      )
    },
    function() {
      fitSampler(
        x,
        yInside,
        weights = wMask,
        family = student,
        rows = which(rowsIn)
      )
    }
  ),
  countsmasked = list(
    4L,
    function() {
      fitSampler(
        x,
        countsMasked,
        weights = wMask,
        family = "nbinom",
        rows = which(rowsIn)
      )
    },
    function() {
      fitSampler(
        x,
        countsInside,
        weights = wMask,
        family = "nbinom",
        rows = which(rowsIn)
      )
    }
  ),
  twoforestsout = list(5L, function() twoForests(yInside, wMask)),
  hazardoffgrid = list(5L, function() fitHazard(timeHOff, wH)),
  hazardbreaks = list(
    5L,
    function() fitHazard(timeHOff, wH, family = dbarts:::hazard(breaks = 3L))
  ),
  gaussianextrememover = list(
    5L,
    function() fitSampler(x, yExtreme, weights = wMask, rows = which(rowsIn))
  ),
  countmover = list(
    5L,
    function() fitSampler(x, counts, weights = wMask, family = "nbinom")
  ),
  xbartzero = list(5L, function() xbartRun(yExtreme, wMask)),
  standin = list(
    5L,
    function() {
      fitSampler(x, yExtreme, weights = c(1, 1, 1, rep(0, n - 3L)))
    }
  ),
  rederive = list(
    5L,
    function() {
      fitSampler(x, yExtreme, after = function(s) {
        s$setWeights(wMask)
        s$setOffset(rep(0, n), updateScale = TRUE)
      })
    }
  )
)

if (mode == "record") {
  out <- lapply(scenarios, function(sc) {
    list(
      group = sc[[1L]],
      asis = sc[[2L]](),
      moved = if (length(sc) >= 3L) sc[[3L]]()
    )
  })
  attr(out, "build") <- dbarts:::buildInfo()
  saveRDS(out, args[2L])
  cat("recorded", length(out), "scenarios to", args[2L], "\n")
} else if (mode == "compare") {
  base <- readRDS(args[2L])
  slice <- readRDS(args[3L])
  stopifnot(identical(names(base), names(slice)))
  worst <- function(a, b) {
    max(
      0,
      unlist(Map(
        function(u, v) {
          # Inf where the fits differ in shape (another hazard grid has
          # another number of person-period rows)
          if (!is.numeric(u)) {
            0
          } else if (!identical(dim(u), dim(v)) || length(u) != length(v)) {
            Inf
          } else {
            max(abs(u - v))
          }
        },
        a,
        b
      ))
    )
  }
  failed <- FALSE
  for (name in names(slice)) {
    group <- slice[[name]]$group
    reference <- if (group %in% c(3L, 4L)) {
      base[[name]]$moved
    } else {
      base[[name]]$asis
    }
    same <- identical(slice[[name]]$asis, reference)
    gap <- worst(slice[[name]]$asis, reference)
    expected <- if (group == 5L) "differ" else "identical"
    ok <- if (group == 5L) !same else same
    failed <- failed || !ok
    cat(sprintf(
      "%-22s group %d  %-9s max abs diff %-10.3g %s\n",
      name,
      group,
      if (same) "identical" else "differs",
      gap,
      if (ok) "ok" else paste0("<- FAIL (expected ", expected, ")")
    ))
  }
  if (failed) {
    cat("\nFAIL\n")
    quit(status = 1L)
  }
  cat("\nOK\n")
} else {
  stop("usage: rows-pair.R record <file> | compare <base> <slice>")
}
