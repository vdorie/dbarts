#!/usr/bin/env Rscript

# Exact-posterior gate for a zero split probability, which excludes its column
# as a column list does on every path, including where the columns of positive
# probability run out along it.
#
# One tree over three predictors: x1 and x2 are 0/1, so each holds one cut and
# is used up below a split on it; x3 is continuous and carries most of the
# signal, so a sampler that splits on it wherever it can is drawn there hard.
# Split probabilities 0.7, 0.3 and 0 for x1, x2, x3. Below a split on x1 and
# one on x2 no column of positive probability is left, so the node cannot
# grow; the tree space is the root leaf and the eight trees of one or two
# levels over x1 and x2, the split-variable prior over the available columns
# of positive probability, the cut prior trivial. Leaf values integrate in
# closed form, the residual sd held fixed.
#
# Two arms, one enumeration:
#   zero  split.probs = c(0.7, 0.3, 0);
#   vars  forest(vars = 1:2) with split.probs = c(0.7, 0.3, 0.5): the third
#         column left out by the list instead, the ratio among the allowed
#         columns the same.
# The enumeration sets the growth probability to 0 at a node with no column
# left, so the two arms agree with it only if the zero is read as the list is.
# Compared, per arm: the posterior probability of each split-count class
# (x1, x2) - (0, 0), (1, 0), (0, 1), (1, 1), (2, 1), (1, 2) - and the posterior
# mean of the fit in each of the four (x1, x2) cells, by batch-means z; any
# split on x3 fails outright.
#
# Sizes and what they detect. Measured on this design, a sampler that proposes
# the first available column where none of positive probability is left (the
# defect this gate holds closed) grows every tree to two levels and splits on
# x3 about five times a sweep below them: it misses the class probabilities by
# up to 0.52 and the cell fits by up to 1.24. Batch-means standard errors run
# near 1.0 / sqrt(N) (class) and 0.35 / sqrt(N) (fit) over N sweeps. Quick mode
# keeps 4 chains x 25 batches x 5000 sweeps, 5e5 sweeps in about a second,
# which puts that defect at hundreds of standard errors; full mode keeps
# 4 x 100 x 10000.
# A run of 2e7 sweeps an arm on two other seeds kept every |z| under 2. The
# bound is |z| <= 4.5 on each of 10 statistics an arm, their standard errors
# from 100 or 400 batch means: by a union bound with Student-t tails a correct
# sampler fails a fresh seed once in 5400 quick runs of an arm and once in
# 11000 full ones. The seed is fixed, so a given build either passes or does
# not.
#
# Usage: Rscript split-probs-exact.R [quick] [zero] [vars]; the zero arm alone
# when neither is named

suppressPackageStartupMessages(library(dbarts))
cgm <- dbartsPriors$cgm
normal <- dbartsPriors$normal
fixed <- dbartsPriors$fixed
gaussian <- dbartsFamilies$gaussian
forest <- dbartsForests$forest

args <- commandArgs(trailingOnly = TRUE)
quick <- "quick" %in% args
armNames <- intersect(c("zero", "vars"), args)
if (length(armNames) == 0L) {
  armNames <- "zero"
}

numChains <- 4L
batchSize <- if (quick) 5000L else 10000L
numBatches <- if (quick) 25L else 100L
engineSeed <- 20261029L
zBound <- 4.5

# ---- the design ----

x1 <- rep(c(0, 1, 0, 1), each = 4L)
x2 <- rep(c(0, 0, 1, 1), each = 4L)
x3 <- c(11, 62, 35, 94, 27, 78, 5, 51, 69, 18, 86, 42, 57, 2, 73, 31) / 100
noise <- c(6, -11, 2, 9, -4, 12, -7, 1, 5, -2, -9, 8, 3, -12, 7, -5) / 100
x <- cbind(x1 = x1, x2 = x2, x3 = x3)
cell <- 1L + x1 + 2L * x2
y <- 0.4 * x1 - 0.3 * x2 + 3 * (x3 > 0.5) + noise
n <- length(y)
splitProbabilities <- c(0.7, 0.3)
sigma <- 0.5
kLeaf <- 2
nodeScale <- 0.5
base <- 0.95
power <- 2

# ---- the exact posterior ----
#
# Each tree as its leaves, each a set of (x1, x2) cells, its split counts and
# its log prior. A node grows with probability base / (1 + depth)^power while
# some column of positive probability is available at it, and with
# probability 0 otherwise; the split variable is drawn in proportion to the
# probabilities of the available columns; a 0/1 column has one cut.
cellsWhere <- function(variable, value) {
  if (variable == 1L) {
    which(c(0, 1, 0, 1) == value)
  } else {
    which(c(0, 0, 1, 1) == value)
  }
}
growth <- function(depth) base / (1 + depth)^power

enumerateTrees <- function() {
  allCells <- 1:4
  trees <- list(list(
    leaves = list(allCells),
    counts = c(0L, 0L),
    logPrior = log(1 - growth(0))
  ))
  for (v in 1:2) {
    u <- 3L - v
    rootLog <- log(growth(0)) +
      log(splitProbabilities[v] / sum(splitProbabilities))
    # each child of the root: a leaf, or a split on the other column, below
    # which nothing is left
    childOptions <- function(cells) {
      list(
        list(leaves = list(cells), splits = 0L, logPrior = log(1 - growth(1))),
        list(
          leaves = list(
            intersect(cells, cellsWhere(u, 0)),
            intersect(cells, cellsWhere(u, 1))
          ),
          splits = 1L,
          logPrior = log(growth(1))
        )
      )
    }
    for (left in childOptions(cellsWhere(v, 0))) {
      for (right in childOptions(cellsWhere(v, 1))) {
        counts <- integer(2L)
        counts[v] <- 1L
        counts[u] <- left$splits + right$splits
        trees[[length(trees) + 1L]] <- list(
          leaves = c(left$leaves, right$leaves),
          counts = counts,
          logPrior = rootLog + left$logPrior + right$logPrior
        )
      }
    }
  }
  trees
}

# the split-count classes (x1, x2)
classes <- c("0,0", "1,0", "0,1", "1,1", "2,1", "1,2")
classOf <- function(x1Splits, x2Splits) {
  match(paste(x1Splits, x2Splits, sep = ","), classes)
}

exactPosterior <- function() {
  yRange <- max(y) - min(y)
  shift <- min(y) + 0.5 * yRange
  tau <- yRange * nodeScale / kLeaf
  y0 <- y - shift
  trees <- enumerateTrees()
  logWeight <- numeric(length(trees))
  cellMean <- matrix(0, length(trees), 4L)
  for (t in seq_along(trees)) {
    logWeight[t] <- trees[[t]]$logPrior
    for (leaf in trees[[t]]$leaves) {
      rows <- cell %in% leaf
      count <- sum(rows)
      s <- sum(y0[rows])
      ss <- sum(y0[rows]^2)
      precision <- 1 / tau^2 + count / sigma^2
      logWeight[t] <- logWeight[t] -
        0.5 * count * log(2 * pi * sigma^2) -
        0.5 * ss / sigma^2 +
        0.5 * log((1 / tau^2) / precision) +
        0.5 * (s / sigma^2)^2 / precision
      cellMean[t, leaf] <- shift + (s / sigma^2) / precision
    }
  }
  w <- exp(logWeight - max(logWeight))
  w <- w / sum(w)
  classProbability <- numeric(length(classes))
  for (t in seq_along(trees)) {
    k <- classOf(trees[[t]]$counts[1L], trees[[t]]$counts[2L])
    classProbability[k] <- classProbability[k] + w[t]
  }
  list(
    numTrees = length(trees),
    class = classProbability,
    fit = as.vector(w %*% cellMean)
  )
}

# ---- the sampler ----

makeSampler <- function(arm, seed) {
  control <- dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 1L,
    keepTrees = FALSE,
    n.samples = batchSize,
    n.burn = 0L,
    n.thin = 1L,
    updateState = FALSE,
    seed = seed
  )
  stated <- if (arm == "zero") {
    list(tree.prior = cgm(power, base, split.probs = c(splitProbabilities, 0)))
  } else {
    list(
      tree.prior = cgm(power, base, split.probs = c(splitProbabilities, 0.5)),
      forests = list(forest(vars = 1:2))
    )
  }
  do.call(
    dbarts,
    c(
      list(
        x,
        y,
        control = control,
        leaf.prior = normal(kLeaf),
        family = gaussian(sigma = fixed(sigma^2))
      ),
      stated
    )
  )
}

# batch means of each class's indicator and of the fit in each cell, one row
# per batch; the first batch of each chain is warm-up and dropped. Any split
# on x3 is counted.
runChain <- function(arm, seed) {
  sampler <- makeSampler(arm, seed)
  firstRow <- match(1:4, cell)
  classMeans <- matrix(0, numBatches, length(classes))
  fitMeans <- matrix(0, numBatches, 4L)
  splitsOnThird <- 0
  for (b in 0:numBatches) {
    draws <- sampler$run(0L, batchSize)
    counts <- draws$varcount
    splitsOnThird <- splitsOnThird + sum(counts[3L, ])
    k <- classOf(counts[1L, ], counts[2L, ])
    if (b > 0L) {
      classMeans[b, ] <- tabulate(k, length(classes)) / batchSize
      fitMeans[b, ] <- rowMeans(draws$train[firstRow, , drop = FALSE])
    }
  }
  list(class = classMeans, fit = fitMeans, splitsOnThird = splitsOnThird)
}

# ---- verdict ----

compare <- function(label, names, batches, exact) {
  estimate <- colMeans(batches)
  se <- apply(batches, 2L, sd) / sqrt(nrow(batches))
  z <- (estimate - exact) / se
  cat(sprintf("  %s\n", label))
  cat(sprintf(
    "    %-8s %9s %9s %9s %7s\n",
    "",
    "exact",
    "sampler",
    "MCse",
    "z"
  ))
  for (j in seq_along(exact)) {
    cat(sprintf(
      "    %-8s %9.5f %9.5f %9.5f %+7.1f%s\n",
      names[j],
      exact[j],
      estimate[j],
      se[j],
      z[j],
      if (is.na(z[j]) || abs(z[j]) > zBound) " <- FAIL" else ""
    ))
  }
  anyNA(z) || any(abs(z) > zBound)
}

exact <- exactPosterior()
anyFailure <- FALSE
for (armIndex in seq_along(armNames)) {
  arm <- armNames[armIndex]
  started <- proc.time()[3L]
  chains <- lapply(
    engineSeed + 100L * armIndex + seq_len(numChains),
    function(seed) runChain(arm, seed)
  )
  splitsOnThird <- sum(vapply(chains, function(chain) chain$splitsOnThird, 0))
  cat(sprintf(
    "\n%s arm: %d trees enumerated; %d chains x %d batches x %d sweeps in %.0f s; %d splits on x3\n",
    arm,
    exact$numTrees,
    numChains,
    numBatches,
    batchSize,
    proc.time()[3L] - started,
    as.integer(splitsOnThird)
  ))
  failedClass <- compare(
    "posterior probability of each split-count class (x1, x2)",
    classes,
    do.call(rbind, lapply(chains, function(chain) chain$class)),
    exact$class
  )
  failedFit <- compare(
    "posterior mean of the fit, by (x1, x2) cell",
    c("0,0", "1,0", "0,1", "1,1"),
    do.call(rbind, lapply(chains, function(chain) chain$fit)),
    exact$fit
  )
  anyFailure <- anyFailure || failedClass || failedFit || splitsOnThird > 0
}

if (anyFailure) {
  cat(
    "\nFAIL: the sampler deviates from the exact posterior under a zero split probability\n"
  )
  quit(status = 1L)
}
cat(
  "\nOK: the sampler matches the exact posterior under a zero split probability\n"
)
