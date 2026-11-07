#!/usr/bin/env Rscript

# Exact gates for frequency weights (docs/plans/frequency-weights.md, S7). A
# row of frequency weight w is w identical observations, so where the sampler
# consumes its random numbers the same way for a row of count w and for w
# copies of it, the weighted fit IS the fit on the rows written out, draw for
# draw; elsewhere it is that fit in distribution, or an exact posterior a
# single leaf admits. Arms:
#
#   gaussian     frequency gaussian against the expanded rows, 300 sweeps:
#                sigma and every row's fit within 1e-12 (the default sigma
#                estimate included, which reads the expanded rows)
#   linear       the same with a linear leaf, its covariate moments counted
#   forests      the same with two forests, the second on a factor basis,
#                whose default sd is the response's sd read at the counts;
#                within 1e-10, its coefficient's sums taken in another order
#   quantiles    the same with useQuantiles and fewer cuts than unique values:
#                the quantile grid reads unique values, which copies leave
#   hetero       a variance forest, each variance leaf counting w
#   nbinom       nbinom at a fixed and at a drawn shape: a row of count w
#                draws its Polya-Gamma latent as w unit draws, in row order
#   multinomial  the K softmax forests the same way
#   binomial     cbind(k, f) at row weight v against two adjacent rows of
#                logistic counts v k and v f, within 1e-10, 0 < k < n
#   sigma        one leaf (a constant predictor) under real weights: sigma's
#                draws against its exact posterior, the leaf integrated out
#                in closed form and sigma^2 on a grid, by batch-means z
#   nu           Student-t's degrees of freedom on fixed residuals (a leaf
#                prior too narrow to move the fit, sigma fixed) at integer,
#                real and zero counts against the grid's exact probabilities,
#                prior times the product of t densities to the counts
#   student      Student-t at integer weights against the expanded rows, the
#                two in distribution (one latent per row against one per
#                copy): Welch z on the posterior means of the fits, sigma and
#                nu, from batch means over independent chains
#
# The draw-for-draw arms are deterministic. The three in-distribution arms use
# |z| <= 4.5 on each statistic, their standard errors from batch means; at a
# fixed seed a given build passes or does not.
#
# Usage: Rscript frequency-exact.R [quick] [arm ...]

suppressPackageStartupMessages(library(dbarts))

args <- commandArgs(trailingOnly = TRUE)
quick <- "quick" %in% args
allArms <- c(
  "gaussian",
  "linear",
  "forests",
  "quantiles",
  "hetero",
  "nbinom",
  "multinomial",
  "binomial",
  "sigma",
  "nu",
  "student"
)
arms <- intersect(allArms, args)
if (length(arms) == 0L) {
  arms <- allArms
}

zBound <- 4.5
anyFailure <- FALSE
verdict <- function(label, fail, text) {
  cat(sprintf("%-12s %s%s\n", label, text, if (fail) "  <- FAIL" else ""))
  if (fail) {
    anyFailure <<- TRUE
  }
}

# ---- draw for draw against the rows written out ----

set.seed(20261108L)
n <- 40L
x <- matrix(runif(2L * n), n, 2L, dimnames = list(NULL, c("a", "b")))
w <- sample(1:4, n, replace = TRUE)
y <- sin(6 * x[, 1L]) + rnorm(n, 0, 0.4)
rows <- rep(seq_len(n), w)

drawControl <- dbartsControl(
  n.samples = 300L,
  n.burn = 0L,
  n.chains = 1L,
  n.threads = 1L,
  updateState = FALSE
)

# The largest difference between two runs' draws, every row-indexed draw
# read at rows 'rowsA' of the first and 'rowsB' of the second; a draw neither
# run indexes by row is compared whole. A draw one run reports and the other
# does not, or of another shape, is a failure (Inf).
drawGap <- function(a, b, rowsA, rowsB) {
  gap <- 0
  for (name in union(names(a), names(b))) {
    if (name == "test" || (is.null(a[[name]]) && is.null(b[[name]]))) {
      next
    }
    da <- a[[name]]
    db <- b[[name]]
    if (is.null(da) || is.null(db)) {
      return(Inf)
    }
    byRow <- !is.null(dim(da)) && dim(da)[1L] == n && name != "varcount"
    if (byRow) {
      da <- asplit(da, 1L)[rowsA]
      db <- asplit(db, 1L)[rowsB]
    }
    da <- unlist(da, use.names = FALSE)
    db <- unlist(db, use.names = FALSE)
    if (length(da) != length(db)) {
      return(Inf)
    }
    gap <- max(gap, abs(as.double(da) - as.double(db)))
  }
  gap
}

drawForDraw <- function(label, weighted, expanded, rowsA, rowsB, tol) {
  gap <- drawGap(weighted$run(), expanded$run(), rowsA, rowsB)
  verdict(
    label,
    !(gap <= tol),
    sprintf("max |difference| over 300 sweeps %.3g (tol %.0e)", gap, tol)
  )
}

# a weighted fit and the same call on the rows written out
expansionPair <- function(response, ..., weights.type = "frequency") {
  expandedResponse <- if (is.matrix(response)) {
    response[rows, , drop = FALSE]
  } else {
    response[rows]
  }
  list(
    weighted = dbarts(
      x,
      response,
      weights = w,
      weights.type = weights.type,
      control = drawControl,
      seed = 9L,
      ...
    ),
    expanded = dbarts(
      x[rows, , drop = FALSE],
      expandedResponse,
      control = drawControl,
      seed = 9L,
      ...
    )
  )
}

runExpansionArm <- function(label, pair, tol = 1e-12) {
  drawForDraw(
    label,
    pair$weighted,
    pair$expanded,
    rows,
    seq_along(rows),
    tol
  )
}

if ("gaussian" %in% arms) {
  runExpansionArm("gaussian", expansionPair(y, n.trees = 20L))
}
if ("linear" %in% arms) {
  runExpansionArm(
    "linear",
    expansionPair(y, n.trees = 20L, leaf.prior = linear("a"))
  )
}
if ("forests" %in% arms) {
  # a factor basis's default sd is one anchor, the response's sd read at the
  # counts; a forest's sd, stated or by default with no basis, is divided by
  # the anchor and multiplied back, so the anchor cancels there
  z <- factor(rep_len(c(0L, 1L), n))
  zRows <- z[rows]
  drawForDraw(
    "forests",
    dbarts(
      x,
      y,
      weights = w,
      weights.type = "frequency",
      control = drawControl,
      seed = 9L,
      forests = list(forest(n.trees = 10L), forest(n.trees = 10L, basis = z))
    ),
    dbarts(
      x[rows, , drop = FALSE],
      y[rows],
      control = drawControl,
      seed = 9L,
      forests = list(
        forest(n.trees = 10L),
        forest(n.trees = 10L, basis = zRows)
      )
    ),
    rows,
    seq_along(rows),
    1e-10
  )
}
if ("quantiles" %in% arms) {
  runExpansionArm(
    "quantiles",
    expansionPair(y, n.trees = 20L, useQuantiles = TRUE, n.cuts = 10L)
  )
}
if ("hetero" %in% arms) {
  runExpansionArm("hetero", expansionPair(y, n.trees = 20L, variance = TRUE))
}
if ("nbinom" %in% arms) {
  counts <- rnbinom(n, size = 3, mu = exp(1 + x[, 1L]))
  runExpansionArm(
    "nbinom",
    expansionPair(
      counts,
      n.trees = 20L,
      family = dbartsFamilies$nbinom(shape = 3),
      weights.type = NULL
    )
  )
  runExpansionArm(
    "nbinom r",
    expansionPair(counts, n.trees = 20L, family = "nbinom", weights.type = NULL)
  )
}
if ("multinomial" %in% arms) {
  categories <- factor(sample(c("p", "q", "r"), n, replace = TRUE))
  runExpansionArm(
    "multinomial",
    expansionPair(
      categories,
      n.trees = 10L,
      family = "multinomial",
      weights.type = NULL
    )
  )
}
if ("binomial" %in% arms) {
  trials <- sample(2:6, n, replace = TRUE)
  successes <- pmin(
    pmax(rbinom(n, trials, plogis(2 * x[, 1L] - 1)), 1L),
    trials - 1L
  )
  failures <- trials - successes
  rowWeights <- sample(1:3, n, replace = TRUE)
  pairRows <- rep(seq_len(n), each = 2L)
  drawForDraw(
    "binomial",
    dbarts(
      x,
      cbind(successes, failures),
      family = binomial,
      weights = rowWeights,
      control = drawControl,
      seed = 9L,
      n.trees = 20L
    ),
    dbarts(
      x[pairRows, , drop = FALSE],
      rep(c(1, 0), n),
      family = "logistic",
      weights = as.vector(rbind(
        rowWeights * successes,
        rowWeights * failures
      )),
      control = drawControl,
      seed = 9L,
      n.trees = 20L
    ),
    seq_len(n),
    seq.int(1L, 2L * n, by = 2L),
    1e-10
  )
}

# ---- batch means ----

# the batch-means standard error of a chain's mean, 'numBatches' batches
batchSe <- function(values, numBatches = 50L) {
  size <- length(values) %/% numBatches
  means <- colMeans(matrix(values[seq_len(size * numBatches)], size))
  sd(means) / sqrt(numBatches)
}

# ---- sigma's exact posterior on one leaf ----

if ("sigma" %in% arms) {
  set.seed(20261109L)
  nSigma <- 30L
  ySigma <- rnorm(nSigma, 1, 0.7)
  wSigma <- runif(nSigma, 0.3, 3)
  sigest <- 0.5
  nu <- 3
  quant <- 0.9
  leafSd <- 0.8
  lambda <- sigest^2 * qchisq(1 - quant, nu) / nu
  # the leaf prior is centered at the response's midrange, over the rows of
  # positive weight (all of them), in response units
  center <- min(ySigma) + 0.5 * (max(ySigma) - min(ySigma))
  total <- sum(wSigma)
  weightedMean <- sum(wSigma * ySigma) / total
  squares <- sum(wSigma * (ySigma - weightedMean)^2)
  # p(sigma^2 | y) on a grid in log sigma: the scaled inverse chi-squared
  # prior, the counted likelihood with the leaf integrated out against its
  # normal prior
  logSigma <- seq(log(0.05), log(5), length.out = 20001L)
  s2 <- exp(2 * logSigma)
  logDensity <- (-nu / 2 - 1) *
    log(s2) -
    nu * lambda / (2 * s2) +
    (-total / 2) * log(s2) +
    0.5 * log(s2 / total) -
    squares / (2 * s2) +
    dnorm(weightedMean, center, sqrt(leafSd^2 + s2 / total), log = TRUE) +
    2 * logSigma # d sigma^2 = 2 sigma^2 d log sigma
  mass <- exp(logDensity - max(logDensity))
  mass <- mass / sum(mass)
  sigmaGrid <- exp(logSigma)
  exactMean <- sum(mass * sigmaGrid)
  cdf <- cumsum(mass)
  probs <- c(0.1, 0.5, 0.9)
  exactQuantiles <- sigmaGrid[vapply(
    probs,
    function(p) which(cdf >= p)[1L],
    1L
  )]

  numSweeps <- if (quick) 200000L else 1000000L
  draws <- dbarts(
    matrix(0, nSigma, 1L, dimnames = list(NULL, "c")),
    ySigma,
    weights = wSigma,
    weights.type = "frequency",
    family = dbartsFamilies$gaussian(sigma = dbartsPriors$chisq(nu, quant)),
    sigest = sigest,
    leaf.prior = normal(sd = leafSd),
    n.trees = 1L,
    control = dbartsControl(
      n.samples = numSweeps,
      n.burn = 1000L,
      n.chains = 1L,
      n.threads = 1L,
      updateState = FALSE
    ),
    seed = 21L
  )$run()$sigma
  statistics <- list(mean = draws)
  for (j in seq_along(probs)) {
    statistics[[sprintf("P(<q%.1f)", probs[j])]] <- as.double(
      draws <= exactQuantiles[j]
    )
  }
  targets <- c(exactMean, probs)
  z <- mapply(
    function(values, target) (mean(values) - target) / batchSe(values),
    statistics,
    targets
  )
  verdict(
    "sigma",
    !all(abs(z) <= zBound),
    sprintf(
      "E[sigma] %.4f vs exact %.4f; z %s",
      mean(draws),
      exactMean,
      paste(sprintf("%.2f", z), collapse = " ")
    )
  )
}

# ---- Student-t's nu on fixed residuals ----

if ("nu" %in% arms) {
  set.seed(20261110L)
  nNu <- 40L
  sigmaNu <- 0.5
  # normal residuals, so the posterior spreads over the whole grid
  yNu <- sigmaNu * rnorm(nNu)
  counts <- c(0, sample(c(1, 2, 3, 0.5, 1.7, 2.4), nNu - 1L, replace = TRUE))
  # a row of count 0 is out of the likelihood, which the fit warns of
  # the fit is the leaf prior's center, the midrange of the rows in
  positive <- counts > 0
  center <- min(yNu[positive]) +
    0.5 * (max(yNu[positive]) - min(yNu[positive]))
  residuals <- yNu - center
  # the engine's grid and its prior, the gamma(2, 0.1) kernel
  grid <- c(3, 4, 5, 6, 8, 10, 12, 15, 20)
  logPrior <- log(grid * exp(-0.1 * grid))
  logPost <- logPrior +
    vapply(
      grid,
      function(df) {
        sum(
          counts[positive] *
            (dt(residuals[positive] / sigmaNu, df, log = TRUE) - log(sigmaNu))
        )
      },
      0
    )
  exact <- exp(logPost - max(logPost))
  exact <- exact / sum(exact)

  numSweeps <- if (quick) 200000L else 1000000L
  draws <- suppressWarnings(dbarts(
    matrix(0, nNu, 1L, dimnames = list(NULL, "c")),
    yNu,
    weights = counts,
    weights.type = "frequency",
    family = dbartsFamilies$student(
      sigma = dbartsPriors$fixed(sigmaNu^2)
    ),
    leaf.prior = normal(sd = 1e-9),
    n.trees = 1L,
    control = dbartsControl(
      n.samples = numSweeps,
      n.burn = 1000L,
      n.chains = 1L,
      n.threads = 1L,
      updateState = FALSE
    ),
    seed = 22L
  ))$run()$resid.df
  z <- vapply(
    seq_along(grid),
    function(k) {
      # a point the chain rarely visits has a batch-means error near 0;
      # the error of independent draws at the exact probability bounds it
      # from below
      hits <- as.double(draws == grid[k])
      se <- max(
        batchSe(hits),
        sqrt(exact[k] * (1 - exact[k]) / length(hits))
      )
      (mean(hits) - exact[k]) / se
    },
    0
  )
  verdict(
    "nu",
    !all(abs(z) <= zBound),
    sprintf(
      "exact %s; max |P - exact| %.4f; max |z| %.2f over %d grid points",
      paste(sprintf("%.3f", exact), collapse = " "),
      max(abs(
        tabulate(match(draws, grid), length(grid)) / length(draws) - exact
      )),
      max(abs(z)),
      length(grid)
    )
  )
}

# ---- Student-t against the copies, in distribution ----

if ("student" %in% arms) {
  set.seed(20261111L)
  nT <- 24L
  xT <- matrix(runif(nT), nT, 1L, dimnames = list(NULL, "a"))
  wT <- sample(1:3, nT, replace = TRUE)
  yT <- 2 * xT[, 1L] + 0.3 * rt(nT, df = 4)
  rowsT <- rep(seq_len(nT), wT)
  numChains <- if (quick) 4L else 8L
  numSweeps <- if (quick) 25000L else 100000L
  summaries <- function(xx, yy, weights, chainRows) {
    fits <- lapply(seq_len(numChains), function(chain) {
      run <- dbarts(
        xx,
        yy,
        weights = weights,
        weights.type = if (!is.null(weights)) "frequency",
        family = "student",
        n.trees = 10L,
        control = dbartsControl(
          n.samples = numSweeps,
          n.burn = 2000L,
          n.chains = 1L,
          n.threads = 1L,
          updateState = FALSE
        ),
        seed = 100L + chain + if (is.null(weights)) 50L else 0L
      )$run()
      rbind(run$train[chainRows, , drop = FALSE], run$sigma, run$resid.df)
    })
    # each statistic's mean over chains and the standard error of that mean,
    # the chains' batch-means variances pooled
    means <- sapply(fits, rowMeans)
    ses <- sapply(fits, function(f) apply(f, 1L, batchSe))
    list(
      mean = rowMeans(means),
      se = sqrt(rowSums(ses^2)) / numChains
    )
  }
  weighted <- summaries(xT, yT, wT, seq_len(nT))
  expanded <- summaries(
    xT[rowsT, , drop = FALSE],
    yT[rowsT],
    NULL,
    match(seq_len(nT), rowsT)
  )
  z <- (weighted$mean - expanded$mean) /
    sqrt(weighted$se^2 + expanded$se^2)
  verdict(
    "student",
    !all(abs(z) <= zBound),
    sprintf(
      "max |z| %.2f over %d fits, sigma (z %.2f) and nu (z %.2f)",
      max(abs(z)),
      nT,
      z[nT + 1L],
      z[nT + 2L]
    )
  )
}

if (anyFailure) {
  quit(status = 1L)
}
cat("\nOK: frequency weights match their expansions and exact posteriors\n")
