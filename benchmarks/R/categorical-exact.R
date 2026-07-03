#!/usr/bin/env Rscript

# Exact-posterior gate for bartcore's categorical splits. A single-tree,
# one-categorical-predictor (4 levels) probit problem admits brute-force
# enumeration of the tree space - rules are direction assignments over the
# categories reachable at each node, uniform over the 2^R - 2 assignments
# that leave neither side empty, both orientations counted - with 1-D
# quadrature for the leaf marginals. Failure means the categorical rule
# prior, draw, or partition logic is wrong.
#
# Usage: Rscript categorical-exact.R [quick]

suppressPackageStartupMessages(library(dbarts))

args <- commandArgs(trailingOnly = TRUE)
quick <- "quick" %in% args
ndpost <- if (quick) 10000L else 50000L
n.seeds <- if (quick) 1L else 3L
tolerance <- if (quick) 0.01 else 0.005

set.seed(1001L)
numCategories <- 4L
nPerCategory <- 30L
category <- rep(0:(numCategories - 1L), each = nPerCategory)
x <- matrix(as.double(category), ncol = 1L)
categoryProbability <- c(0.25, 0.7, 0.45, 0.6)
y <- rbinom(length(category), 1L, categoryProbability[category + 1L])
n <- length(y)
offset <- qnorm(mean(y))
x.test <- matrix(as.double(0:(numCategories - 1L)), ncol = 1L)

base <- 0.5
power <- 2
k <- 2
tau <- 3.0 / k  # dbarts binary node scale over a single tree

successesByCategory <- as.vector(tapply(y, category, sum))
countsByCategory <- as.vector(table(category))

popcount <- function(mask) sum(as.integer(intToBits(mask)))

# trees over a reachable-category mask: leaf w.p. 1 - base / (1 + d)^power
# (forced when fewer than two categories reach the node), else an assignment
# sending a proper nonempty subset right, uniform over the 2^R - 2 of them
enumerate <- function(mask, depth) {
  numReachable <- popcount(mask)
  growth <- if (numReachable >= 2L) base / (1 + depth)^power else 0
  result <- list(list(leaves = list(mask), logPrior = log(1 - growth)))
  if (numReachable < 2L) return(result)

  subsets <- Filter(
    function(d) d != 0L && d != mask && bitwAnd(d, bitwNot(mask)) == 0L,
    seq_len(15L)
  )
  for (directions in subsets) {
    lefts <- enumerate(bitwAnd(mask, bitwNot(directions)), depth + 1L)
    rights <- enumerate(directions, depth + 1L)
    for (left in lefts) for (right in rights) {
      result[[length(result) + 1L]] <- list(
        leaves = c(left$leaves, right$leaves),
        logPrior = log(growth) - log(2^numReachable - 2) +
          left$logPrior + right$logPrior
      )
    }
  }
  result
}

trees <- enumerate(bitwShiftL(1L, numCategories) - 1L, 0L)
cat(sprintf("enumerated %d trees\n", length(trees)))

muGrid <- seq(-10, 10, by = 0.005)
muDensity <- dnorm(muGrid, 0, tau)
leafQuantities <- function(s, m) {
  w <- muDensity * exp(s * pnorm(muGrid + offset, log.p = TRUE) +
                         (m - s) * pnorm(-(muGrid + offset), log.p = TRUE))
  list(logMarginal = log(sum(w) * 0.005),
       meanProbability = sum(pnorm(muGrid + offset) * w) / sum(w))
}

logWeights <- numeric(length(trees))
predictions <- matrix(0, length(trees), numCategories)
for (t in seq_along(trees)) {
  logWeight <- trees[[t]]$logPrior
  for (leafMask in trees[[t]]$leaves) {
    inLeaf <- which(bitwAnd(bitwShiftL(1L, 0:(numCategories - 1L)),
                            leafMask) != 0L)
    q <- leafQuantities(sum(successesByCategory[inLeaf]),
                        sum(countsByCategory[inLeaf]))
    logWeight <- logWeight + q$logMarginal
    predictions[t, inLeaf] <- q$meanProbability
  }
  logWeights[t] <- logWeight
}
weights <- exp(logWeights - max(logWeights))
exact <- colSums(weights * predictions) / sum(weights)

fitBartcore <- function(seed) {
  set.seed(seed)
  control <- dbartsControl(n.chains = 1L, n.threads = 1L, n.trees = 1L,
                           updateState = FALSE)
  host <- dbarts(x, y, test = x.test, offset = offset, control = control,
                 node.prior = normal(k), tree.prior = cgm(power, base))
  # no public surface marks columns categorical yet; flip the type by hand
  host$data@varTypes[1L] <- 1L
  host$data@offset.test <- NULL
  bc <- dbarts:::bartcoreSampler(host)
  r <- dbarts:::bartcoreRun(bc, 5000L, ndpost)
  colMeans(pnorm(t(r$test) + offset))
}

fit <- colMeans(do.call(rbind, lapply(seq_len(n.seeds), fitBartcore)))

cat(sprintf("%9s %s\n", "exact", paste(sprintf("%.4f", exact), collapse = " ")))
cat(sprintf("%9s %s\n", "sampler", paste(sprintf("%.4f", fit), collapse = " ")))
gap <- max(abs(fit - exact))
cat(sprintf("max gap %.4f%s\n", gap, if (gap > tolerance) " <- FAIL" else ""))

if (gap > tolerance) quit(status = 1L)
cat("OK: categorical sampler matches the exact posterior\n")
