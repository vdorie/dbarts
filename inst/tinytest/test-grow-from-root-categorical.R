# grow-from-root now places real categorical split rules (the v1
# "categorical predictors are ordinal-only" contract inverted): a
# categorical-heavy design should see the init forest
# split on its signal factors, warm-start bart2 to a lower early RMSE than a
# cold start, and feed real split counts into a concurrent DART update.

makeFactorSignalData <- function(n = 400L, seed) {
  set.seed(seed)
  g1 <- factor(sample(letters[1:5], n, replace = TRUE))
  g2 <- factor(sample(letters[1:5], n, replace = TRUE))
  gNoise <- factor(sample(letters[1:5], n, replace = TRUE))
  x1 <- runif(n)
  x2 <- runif(n)
  effect1 <- c(a = -3, b = -1.5, c = 0, d = 1.5, e = 3)
  effect2 <- c(a = 2, b = 1, c = 0, d = -1, e = -2)
  mu <- unname(effect1[g1] + effect2[g2])
  y <- mu + rnorm(n, 0, 0.3)
  list(
    x = data.frame(g1 = g1, g2 = g2, gNoise = gNoise, x1 = x1, x2 = x2),
    y = y
  )
}

d <- makeFactorSignalData(seed = 20260812L)

## the init forest places categorical rules on the signal factors (g1, g2 are
## columns 1 and 2; gNoise, x1, x2 carry no signal)
set.seed(1L)
control <- dbarts::dbartsControl(
  n.trees = 50L,
  n.chains = 1L,
  updateState = FALSE,
  keepTrees = FALSE
)
sampler <- dbarts::dbarts(d$x, d$y, control = control)
sampler$growFromRoot(2L)
trees <- sampler$getTrees(current = TRUE)
isSignalRule <- trees$var %in% c(1L, 2L)
expect_true(any(trees$var == 1L))
expect_true(any(trees$var == 2L))
# a categorical rule's split decodes into 'directions'; ordinal rules and
# leaves report NA there
expect_true(!anyNA(trees$directions[isSignalRule]))

## a bart2 fit warm-started off two grow sweeps reaches a lower early training
## RMSE than a cold start, shaped like test-bart2-grow-from-root.R's assertion
earlyRMSE <- function(fit) {
  yh <- fit$yhat.train
  perSample <- sqrt(rowMeans(sweep(yh, 2L, d$y)^2))
  mean(perSample[seq_len(min(10L, length(perSample)))])
}
coldFit <- dbarts::bart2(
  d$x,
  d$y,
  n.trees = 50L,
  power = 1.5,
  n.samples = 20L,
  n.burn = 0L,
  n.chains = 1L,
  verbose = FALSE,
  seed = 2L
)
growFit <- dbarts::bart2(
  d$x,
  d$y,
  n.trees = 50L,
  power = 1.5,
  n.samples = 20L,
  n.burn = 0L,
  n.chains = 1L,
  verbose = FALSE,
  seed = 2L,
  n.grow.sweeps = 2L
)
expect_true(earlyRMSE(growFit) < 0.9 * earlyRMSE(coldFit))

## the DART assertion: the grow sweeps' own close-of-sweep DART update
## (chain.hpp) now sees real categorical split counts instead of structural
## zeros, so the signal factors carry more posterior split probability right
## after the grow sweeps (n.burn = 0, so the first sample is the first DART
## draw following them) than a matched no-grow arm gets from the same n.burn
## = 0 first draw, which never saw a categorical split
fitDart <- function(growSweeps) {
  dbarts::bart2(
    d$x,
    d$y,
    dart = TRUE,
    n.trees = 50L,
    n.samples = 5L,
    n.burn = 0L,
    n.chains = 1L,
    verbose = FALSE,
    seed = 4L,
    n.grow.sweeps = growSweeps
  )$varprobs
}
growVarprobs <- fitDart(2L)
noGrowVarprobs <- fitDart(0L)
growSignalMass <- growVarprobs[1L, "g1"] + growVarprobs[1L, "g2"]
noGrowSignalMass <- noGrowVarprobs[1L, "g1"] + noGrowVarprobs[1L, "g2"]
expect_true(growSignalMass > 0.5)
expect_true(growSignalMass > noGrowSignalMass)

## edge case: every categorical column has fewer than 2 present categories at
## the root (P < 2), so grow can never split on any of them; the fit must
## still complete cleanly and produce a legal, finite forest
set.seed(6L)
nDegenerate <- 200L
h1 <- factor(rep("a", nDegenerate), levels = c("a", "b", "c"))
h2 <- factor(rep("a", nDegenerate), levels = letters[1:4])
x1Degenerate <- runif(nDegenerate)
yDegenerate <- 3 * x1Degenerate + rnorm(nDegenerate, 0, 0.2)
xDegenerate <- data.frame(h1 = h1, h2 = h2, x1 = x1Degenerate)
degenerateSampler <- dbarts::dbarts(
  xDegenerate,
  yDegenerate,
  control = control
)
degenerateSampler$growFromRoot(2L)
degenerateFitted <- as.vector(degenerateSampler$predict(xDegenerate))
expect_true(all(is.finite(degenerateFitted)))
degenerateTrees <- degenerateSampler$getTrees(current = TRUE)
expect_true(nrow(degenerateTrees) > 0L)
# neither degenerate categorical column (var 1 = h1, var 2 = h2) is ever
# chosen, but the ordinal noise column still can be
expect_true(!any(degenerateTrees$var %in% c(1L, 2L)))
