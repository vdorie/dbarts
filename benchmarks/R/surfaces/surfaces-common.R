# Shared generating processes and readouts for the published-surface battery.
# Each cell script in this directory sources this file, runs one problem from
# the literature against the SHIPPED default sampler, and reports its
# pre-registered primary statistic beside the number the source paper
# published for it.
#
# Every generating process below is transcribed from the primary source and
# carries the citation in a comment above it. Where a source does not print a
# function it uses, the reconstruction is labelled as such and the constraints
# it was built to are stated; a reconstruction cannot carry a reproduction
# verdict on its own.
#
# Not sourced directly: source() it from a cell script, or from R to reuse the
# generators.

suppressPackageStartupMessages(library(dbarts))

# --- output location -------------------------------------------------------

# Results land under a directory given as the first positional argument to a
# cell script. This default is a scratch path, not a repository path: nothing
# in this battery writes into the working tree.
surfacesDefaultOutputDir <- file.path(
  "/private/tmp/claude-501",
  "-Users-vdorie-Repositories-dbarts",
  "9a3b1fb9-76e6-4ef5-8d66-788d46a4c9a5",
  "scratchpad",
  "battery-pilot"
)

# First non-flag argument, else the default. Flags are the bare words a cell
# script recognizes (`quick`, arm names); anything else is read as a path.
surfacesOutputDir <- function(args, flags = character(0)) {
  paths <- setdiff(args, flags)
  dir <- if (length(paths) > 0L) paths[1L] else surfacesDefaultOutputDir
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
  dir
}

surfacesSave <- function(result, dir, name) {
  path <- file.path(dir, paste0(name, ".rds"))
  saveRDS(result, path)
  cat(sprintf("\nsaved %s\n", path))
  invisible(path)
}

# --- the matched-seed idiom ------------------------------------------------

# Data seed indexed by cell, design and replicate; sampler seed the replicate
# alone. Both are shared across arms, so within a replicate every arm sees the
# identical data set and the identical MCMC stream and every contrast between
# arms is paired. Distinct cells and designs get disjoint data streams because
# the cell and design names enter the seed.
surfacesDataSeed <- function(cell, design, replicate) {
  chars <- utf8ToInt(paste0(cell, ":", design))
  offset <- sum(chars * seq_along(chars)) %% 20000L
  as.integer(100000L + 1000L * offset + replicate)
}

surfacesSamplerSeed <- function(replicate) {
  as.integer(replicate)
}

# --- P2: the confounded step function --------------------------------------

# Pratola (2016), "Efficient Metropolis-Hastings proposal mechanisms for
# Bayesian regression tree models", Bayesian Analysis 11(3), section 2.3,
# taking the problem from Wu, Tjelmeland and West (2007). Transcribed from
# arXiv 1312.1895 equation (1) and the sentence following it:
#
#   y = 1 + N(0, 0.25)  if x1 <= 0.5 and x2 <= 0.5
#       3 + N(0, 0.25)  if x1 <= 0.5 and x2 >  0.5
#       5 + N(0, 0.25)  if x1 >  0.5
#
#   x1 ~ unif(0.1, 0.4), i =   1..200;  unif(0.6, 0.9), i = 201..300
#   x2 ~ unif(0.1, 0.4), i =   1..100;  unif(0.6, 0.9), i = 101..200;
#        unif(0.1, 0.9), i = 201..300
#   x3 ~ unif(0.6, 0.9), i =   1..200;  unif(0.1, 0.4), i = 201..300
#
# x1 and x3 are confounded by that block structure: on the 200/100 split that
# separates the mean levels, a cut on x3 partitions the same rows a cut on x1
# does, so root-on-x1 and root-on-x3 are exchangeable representations of one
# fitted function. N(0, 0.25) is read as variance 0.25 (the paper writes
# sigma^2 elsewhere in the same notation), so sd = 0.5. n = 300, p = 3.
surfacesConfoundedStep <- function() {
  n <- 300L
  x1 <- c(runif(200L, 0.1, 0.4), runif(100L, 0.6, 0.9))
  x2 <- c(runif(100L, 0.1, 0.4), runif(100L, 0.6, 0.9), runif(100L, 0.1, 0.9))
  x3 <- c(runif(200L, 0.6, 0.9), runif(100L, 0.1, 0.4))
  f <- ifelse(x1 > 0.5, 5, ifelse(x2 > 0.5, 3, 1))
  x <- cbind(x1 = x1, x2 = x2, x3 = x3)
  list(x = x, y = f + rnorm(n, 0, 0.5), f = f)
}

# The null control the battery requires alongside P2: two EXACTLY duplicated
# predictor columns. Their likelihood and prior ratios are identical, so the
# fraction of draws rooted on the first of the pair is exactly 1/2 by
# construction and switching between them must occur. An arm that cannot
# return 1/2 here, with a non-zero switch count in every chain, is measuring
# its own harness rather than the posterior.
#
# The columns are drawn on a FOUR-VALUE grid rather than continuously. A
# continuous duplicate pair is not a fair null: a rule-changing proposal must
# hit the twin column at the same cut point out of the whole grid, so the
# switch it is supposed to make is rare for a reason that has nothing to do
# with the kernel being tested. Four values leave three cuts per column and
# put the matching proposal within reach on every sweep.
surfacesDuplicateColumnNull <- function() {
  n <- 300L
  grid <- c(0.125, 0.375, 0.625, 0.875)
  x1 <- rep(grid, each = n %/% 4L)
  x3 <- sample(grid, n, replace = TRUE)
  x <- cbind(x1 = x1, x2 = x1, x3 = x3)
  f <- 1 + 4 * (x1 > 0.5)
  list(x = x, y = f + rnorm(n, 0, 0.5), f = f)
}

# --- P6: the diagonal shelf with targeted selection ------------------------

# Hahn, Murray and Carvalho (2020), "Bayesian regression tree models for
# causal inference", Bayesian Analysis 15(3), Example 1 (arXiv 1706.09523v4
# section 4.2): d = 2, n = 250, homogeneous effect, x1 and x2 iid U(0, 1),
# eps ~ N(0, 1), 200 replications, treatment effect -1.
#
# The propensity is the paper's, with one correction. The printed equation is
#
#   pi = 0.8 Phi(mu / (0.1(2 - x1 - x2) + 0.25)) + 0.025(x1 + x2) + 0.05,
#
# but the paper's own LaTeX source carries the generating expression on the
# line above it as a comment, and that expression divides mu by 3 first:
# `0.8*pnorm(m/3, 0, 0.1*(2-xtilde)+0.25) + 0.025*xtilde + 0.05`. Only the
# commented form reproduces the paper's Figure 4, which plots mu against pi
# for one realization: it puts pi near 0.21 at mu = -1 and near 0.10 at
# mu = -1.9, where the printed form puts pi near 0.08 and 0.05. The commented
# form is used here.
#
# RECONSTRUCTION. The paper never prints mu. Everything known about it is:
# a "shelf" at the line x1 = x2 (Figure 3 caption), values ranging from -3 to
# 3 (same caption), a "step function (or near-step function) along the
# diagonal" (Figure 5 caption), and a tight monotone mu-pi relation with
# pi = 0.5 near mu = 0 (Figure 4). The family below satisfies all of them:
#
#   mu(x1, x2) = 3 (2 Phi((x1 - x2) / shelfWidth) - 1)
#
# with shelfWidth setting how near a step the shelf is. It is NOT the paper's
# function, and the published bias and coverage cannot be attributed to a
# single width without it; the cell therefore sweeps the width.
surfacesShelfMu <- function(x1, x2, shelfWidth) {
  3 * (2 * pnorm((x1 - x2) / shelfWidth) - 1)
}

surfacesShelfPropensity <- function(mu, x1, x2) {
  xTilde <- x1 + x2
  0.8 * pnorm((mu / 3) / (0.1 * (2 - xTilde) + 0.25)) + 0.025 * xTilde + 0.05
}

surfacesDiagonalShelf <- function(n = 250L, shelfWidth = 0.15, tau = -1) {
  x1 <- runif(n)
  x2 <- runif(n)
  mu <- surfacesShelfMu(x1, x2, shelfWidth)
  pi <- surfacesShelfPropensity(mu, x1, x2)
  z <- rbinom(n, 1L, pi)
  list(
    x = cbind(x1 = x1, x2 = x2),
    z = z,
    y = mu + tau * z + rnorm(n),
    mu = mu,
    propensity = pi,
    tau = tau
  )
}

# --- P5: the checkerboard on an autocorrelated design ----------------------

# Zhu, Zeng and Kosorok (2015) scenario 3, transcribed from Feng and
# Baumgartner (2020), arXiv 2012.10737 section 4.1, which restates it:
#
#   X ~ N(0, Sigma),  Sigma_jk = 0.9^|j-k|
#   eps ~ N(0, 1)
#   f(X) = 2 x5 x10 + 2 x15 x20
#
# run at n in {800, 1600} and p in {20, 40}. Inclusion has an exactly right
# answer, {x5, x10, x15, x20}, and near-decoys: the immediate neighbours of
# each true column correlate 0.9 with it and the next ones out 0.81.
surfacesCheckerboardSignal <- c(5L, 10L, 15L, 20L)

surfacesCheckerboardNeighbours <- function(p) {
  nb <- c(surfacesCheckerboardSignal - 1L, surfacesCheckerboardSignal + 1L)
  sort(nb[nb >= 1L & nb <= p])
}

surfacesCheckerboard <- function(n = 1600L, p = 40L, nTest = 1000L) {
  lags <- abs(outer(seq_len(p), seq_len(p), "-"))
  sigmaRoot <- chol(0.9^lags)
  x <- matrix(rnorm((n + nTest) * p), n + nTest, p) %*% sigmaRoot
  colnames(x) <- paste0("x", seq_len(p))
  f <- 2 * x[, 5L] * x[, 10L] + 2 * x[, 15L] * x[, 20L]
  train <- seq_len(n)
  list(
    x = x[train, , drop = FALSE],
    y = f[train] + rnorm(n),
    f = f[train],
    xTest = x[-train, , drop = FALSE],
    fTest = f[-train]
  )
}

# --- C1: the He and Hahn factorial -----------------------------------------

# He and Hahn (2023), "Stochastic tree ensembles for regularized nonlinear
# regression", JASA, section 4.1 (arXiv 2002.03375v4, Tables 1 and 2).
#
# Mean functions, Table 1 verbatim:
#   Linear        x' gamma;  gamma_j = -2 + 4(j-1)/(d-1)
#   Single index  10 sqrt(a) + sin(5a);  a = sum_{j=1..10} (x_j - gamma_j)^2,
#                 gamma_j = -1.5 + (j-1)/3
#   Trig+poly     5 sin(3 x1) + 2 x2^2 + 3 x3 x4
#   Max           max(x1, x2, x3)
#
# Note that Single index carries its OWN gamma, different from Linear's.
#
# Correlated predictors with factor structure, section 4.1 verbatim: k = p/5
# factors, F (k x n) ~ N(0, 1); the loading matrix B (p x k) is 0/1 with
# exactly five ones in each column and a single 1 in each row, so B B' is
# block diagonal; X = (BF)' + eps with eps entries iid N(0, 0.01k); finally
# each column of X is scaled to standard deviation 1.
#
# Errors: eps_i ~ N(0, sigma^2) with sigma^2 = kappa^2 Var(f).
surfacesHeHahnMean <- function(x, which) {
  d <- ncol(x)
  if (which == "trigpoly") {
    return(5 * sin(3 * x[, 1L]) + 2 * x[, 2L]^2 + 3 * x[, 3L] * x[, 4L])
  }
  if (which == "singleindex") {
    gamma <- -1.5 + (seq_len(10L) - 1L) / 3
    a <- rowSums(
      (x[, seq_len(10L), drop = FALSE] -
        rep(gamma, each = nrow(x)))^2
    )
    return(10 * sqrt(a) + sin(5 * a))
  }
  if (which == "linear") {
    gamma <- -2 + 4 * (seq_len(d) - 1L) / (d - 1L)
    return(as.numeric(x %*% gamma))
  }
  if (which == "max") {
    return(pmax(x[, 1L], x[, 2L], x[, 3L]))
  }
  stop("unknown mean function '", which, "'")
}

surfacesFactorDesign <- function(n, p) {
  k <- p %/% 5L
  factors <- matrix(rnorm(k * n), k, n)
  loadings <- matrix(0, p, k)
  for (j in seq_len(k)) {
    loadings[((j - 1L) * 5L + 1L):(j * 5L), j] <- 1
  }
  x <- t(loadings %*% factors) +
    matrix(rnorm(n * p, 0, sqrt(0.01 * k)), n, p)
  x <- sweep(x, 2L, apply(x, 2L, sd), "/")
  colnames(x) <- paste0("x", seq_len(p))
  x
}

surfacesIndependentDesign <- function(n, p) {
  x <- matrix(rnorm(n * p), n, p)
  colnames(x) <- paste0("x", seq_len(p))
  x
}

# One replicate of C1: a training set and a held-out set drawn from the same
# design, with sigma calibrated on the TRAINING f so kappa is the training
# signal-to-noise ratio the paper defines. Section 5 of the paper, which
# carries the coverage table, does not say which of the two predictor arms
# it used, so the arm is a parameter here.
surfacesHeHahn <- function(n, nTest, p, which, kappa, design = "correlated") {
  x <- if (design == "correlated") {
    surfacesFactorDesign(n + nTest, p)
  } else {
    surfacesIndependentDesign(n + nTest, p)
  }
  f <- surfacesHeHahnMean(x, which)
  train <- seq_len(n)
  sigma <- kappa * sd(f[train])
  list(
    x = x[train, , drop = FALSE],
    y = f[train] + rnorm(n, 0, sigma),
    f = f[train],
    xTest = x[-train, , drop = FALSE],
    fTest = f[-train],
    sigma = sigma
  )
}

# --- P1: the low-noise Friedman emulator -----------------------------------

# Friedman (1991) five-dimensional test function on p iid uniform(0, 1)
# columns, as Chipman, George and McCulloch (2010) carry it (arXiv 0806.3286
# equations 26-27):
#
#   f(x) = 10 sin(pi x1 x2) + 20 (x3 - 0.5)^2 + 10 x4 + 5 x5
#
# Pratola (2016), "Efficient Metropolis-Hastings proposal mechanisms for
# Bayesian regression tree models", Bayesian Analysis 11(3), takes the same
# function as a deterministic simulator eta and observes it with noise,
# y(x) = eta(x) + eps, eps ~ N(0, sigma^2), at n = 5000 settings with m = 200
# trees. The section carrying that example is 2.2 in arXiv 1312.1895, which is
# the numbering the rest of this battery cites; the published version numbers
# it 2.3 and the confounded step function 2.2.
#
# His printed eta is 10 sin(2 pi x1 x2) + 20 (x3 - 0.5)^2 + 10 x4 + 5 x5, a
# full period in x1 x2 rather than Friedman's half, and the paper does not say
# which of the two his acceptance and coverage figures were taken on.
# `frequency` selects: 1 is Friedman's function, 2 is Pratola's as printed.
surfacesFriedmanMean <- function(x, frequency = 1) {
  10 *
    sin(frequency * pi * x[, 1L] * x[, 2L]) +
    20 * (x[, 3L] - 0.5)^2 +
    10 * x[, 4L] +
    5 * x[, 5L]
}

surfacesFriedman <- function(
  n,
  nTest = 1000L,
  p = 10L,
  sigma = 1,
  frequency = 1
) {
  x <- matrix(runif((n + nTest) * p), n + nTest, p)
  colnames(x) <- paste0("x", seq_len(p))
  f <- surfacesFriedmanMean(x, frequency)
  train <- seq_len(n)
  list(
    x = x[train, , drop = FALSE],
    y = f[train] + rnorm(n, 0, sigma),
    f = f[train],
    xTest = x[-train, , drop = FALSE],
    fTest = f[-train],
    sigma = sigma
  )
}

# --- readouts --------------------------------------------------------------

# Rows of a bart2 draw matrix belonging to one chain. The combined layout is
# chain-major: row (c - 1) * nSamples + s is chain c, sample s.
surfacesChainRows <- function(chain, nSamples) {
  ((chain - 1L) * nSamples + 1L):(chain * nSamples)
}

# Pointwise frequentist coverage of a known truth by equal-tailed posterior
# intervals. `draws` is draws x points, `truth` is length points.
surfacesCoverage <- function(draws, truth, level = 0.95) {
  alpha <- (1 - level) / 2
  lower <- apply(draws, 2L, quantile, probs = alpha, names = FALSE)
  upper <- apply(draws, 2L, quantile, probs = 1 - alpha, names = FALSE)
  mean(truth >= lower & truth <= upper)
}

surfacesRmse <- function(draws, truth) {
  sqrt(mean((colMeans(draws) - truth)^2))
}

# posterior::ess_basic on a draws x chains matrix, or on a plain vector for a
# single chain. Section 6.2 of the battery carries minimum ESS over a fixed
# set of points as the secondary statistic, so this is applied per point.
surfacesEss <- function(draws) {
  posterior::ess_basic(draws)
}

surfacesPointEss <- function(draws, points) {
  vapply(points, function(j) surfacesEss(draws[, j]), numeric(1L))
}

# varcount converted to per-draw inclusion proportions: each row of `counts`
# is one draw's split counts over the predictors, so dividing by the row total
# gives the share of that draw's splits falling on each column.
surfacesInclusion <- function(counts) {
  totals <- rowSums(counts)
  totals[totals == 0] <- NA_real_
  counts / totals
}

# The root split variable of every kept draw, from extract(type = "trees").
# Trees are listed in depth-first preorder, so the first row of each
# (chain, sample, tree) block is the root; var is the 1-based predictor index
# and -1 marks a leaf.
surfacesRootVariable <- function(trees) {
  key <- paste(trees$chain, trees$sample, trees$tree, sep = "/")
  first <- !duplicated(key)
  data.frame(
    chain = trees$chain[first],
    sample = trees$sample[first],
    tree = trees$tree[first],
    var = trees$var[first]
  )
}

# An acceptance-rate proxy that needs no engine hook: with one tree and no
# thinning, one kept draw is one sweep, so the tree structure differs from the
# previous draw's exactly when a structural move was accepted. Leaf values are
# redrawn every sweep regardless, so only the internal nodes' variables and
# cut points enter the comparison.
surfacesStructureChangeRate <- function(trees) {
  key <- paste(trees$chain, trees$sample, sep = "/")
  internal <- trees$var > 0L
  signature <- vapply(
    split(
      paste(trees$var, ifelse(internal, trees$value, NA_real_)),
      factor(key, levels = unique(key))
    ),
    function(rows) paste(rows, collapse = "|"),
    character(1L)
  )
  chain <- vapply(strsplit(names(signature), "/", fixed = TRUE), `[`, "", 1L)
  changed <- signature[-1L] != signature[-length(signature)]
  sameChain <- chain[-1L] == chain[-length(chain)]
  mean(changed[sameChain])
}

# --- printing --------------------------------------------------------------

# The battery reports a mean over seeds with the min-max range in
# parentheses, which is the shape the move-set tables already use.
surfacesRange <- function(x, digits = 3L) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) {
    return("      -")
  }
  fmt <- paste0("%.", digits, "f(%.", digits, "f-%.", digits, "f)")
  sprintf(fmt, mean(x), min(x), max(x))
}

# A paired difference against a control arm, replicate by replicate, in the
# shape the move-set tables quote a raw signal in: mean, standard deviation,
# the count of replicates on which the difference is positive, and the t
# statistic of the mean against zero. The t is what the battery's margins are
# read against, so it is printed beside the mean rather than left to be
# recomputed from it.
surfacesPairedDifference <- function(x, digits = 3L) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) {
    return("      -")
  }
  se <- sd(x) / sqrt(length(x))
  fmt <- paste0("%+.", digits, "f +/- %.", digits, "f (%d/%d) t %s")
  sprintf(
    fmt,
    mean(x),
    sd(x),
    sum(x > 0),
    length(x),
    if (se > 0) sprintf("%.2f", mean(x) / se) else "-"
  )
}

# The battery flags a regression only when BOTH of its conditions hold: the
# paired mean difference is worse than the margin, AND its one-sided 95%
# bound excludes the null, so that a noisy cell cannot flag on its point
# estimate alone. `worse` is the sign of a difference that counts as a
# regression: -1 where a smaller value is worse (coverage, an inclusion
# share) and +1 where a larger one is (any error).
surfacesMarginVerdict <- function(x, margin, worse) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) {
    return("-")
  }
  loss <- worse * mean(x)
  se <- sd(x) / sqrt(length(x))
  if (loss <= abs(margin)) {
    return("within margin")
  }
  if (loss - 1.645 * se > 0) {
    return("FLAG")
  }
  "past margin, not separated"
}

# The mirror image, and deliberately harder: "better on at least one
# pathology" means a paired improvement on the pre-registered primary
# exceeding FOUR times the paired standard error. `better` is the sign of an
# improving difference.
surfacesImprovementVerdict <- function(x, better) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) {
    return("-")
  }
  se <- sd(x) / sqrt(length(x))
  if (!(se > 0)) {
    return("no spread")
  }
  t <- better * mean(x) / se
  sprintf("%s (t %.2f)", if (t > 4) "clears 4x SE" else "below 4x SE", t)
}

surfacesHeader <- function(title) {
  cat(sprintf("\n== %s ==\n", title))
}

surfacesUptime <- function(label) {
  cat(sprintf("%s: %s\n", label, trimws(system("uptime", intern = TRUE))))
}
