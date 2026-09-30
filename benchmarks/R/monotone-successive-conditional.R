#!/usr/bin/env Rscript

# Successive-conditional check of the monotone sampler (Geweke 2004): the
# kernel leaves the joint law of (theta, y) invariant. Each replication draws
# theta0 = (tree, leaves, sigma) from the sampler's own prior, simulates y given
# theta0, installs theta0 as the chain's state and runs K sweeps given y. If the
# kernel targets the posterior under that prior, theta_K is again a prior draw,
# so every functional's paired difference g(theta_K) - g(theta0) has mean zero.
# The pairing cancels most of the replication variance, and the check needs no
# mixing: the chain starts at a posterior draw, so a one-tree birth/death chain
# that does not mix its structure on an informative design (the reason the
# one-tree SBC arm is only a diagnostic, docs/plans/monotone-exact-birth-death.md
# step 16) is still tested exactly.
#
# Arms: the "leaf" and "joint" monotone priors (x1 increasing) and the
# unconstrained birth/death-only twin, one tree, one predictor, n 100, under a
# tree prior of power 0.5 so trees are deep enough that most moves touch a leaf
# with a frozen constrained neighbor. Functionals: the leaf count, avg f, f at
# three test points, sigma, and the contrasts f(0.9) - f(0.1) and
# f(0.55) - f(0.45). An arm fails when any functional's paired |z| exceeds
# zBound.
#
# Power: with the move's old ratio restored (the d divisions back, the Z term
# dropped), quick mode fails "leaf" at |z| 9.0 and "joint" at 31 on the leaf
# count, in about a minute on arm64 macOS; under the default tree prior
# (power 2) "leaf" read only |z| 2.3.
#
# Usage: Rscript monotone-successive-conditional.R [quick] [arm ...]
#   arm  leaf, joint, twin (default: all three)

suppressPackageStartupMessages(library(dbarts))

args <- commandArgs(trailingOnly = TRUE)
quick <- "quick" %in% args
arms <- intersect(c("leaf", "joint", "twin"), args)
if (!length(arms)) {
  arms <- c("leaf", "joint", "twin")
}

nReps <- if (quick) 6000L else 20000L
nSweeps <- 20L
nObs <- 100L
# tree prior power 0.5 in place of 2: deeper trees, so most moves touch a leaf
# with a frozen constrained neighbor, where the move ratio is hardest
treePower <- 0.5
zBound <- 4.5
sigest <- 1
sigDf <- 3
sigQuant <- 0.9

set.seed(1L)
x <- matrix(runif(nObs), nObs, 1L, dimnames = list(NULL, "x1"))
xTest <- matrix(
  c(0.25, 0.5, 0.75, 0.1, 0.9, 0.45, 0.55),
  ncol = 1L,
  dimnames = list(NULL, "x1")
)
# a deterministic build response fixes the internal scale once; setResponse
# keeps it (updateScale = FALSE), so prior draw and posterior share one scale
yBuild <- seq(-2.5, 2.5, length.out = nObs)

# the engine's sigma prior on the reported scale: sigma^2 is scaled inverse
# chi-squared with P(sigma < sigest) = sigQuant
drawSigma <- function() {
  rawScale <- qchisq(1 - sigQuant, sigDf) / sigDf
  sqrt(sigDf * sigest^2 * rawScale / rchisq(1L, sigDf))
}

makeSampler <- function(arm, seed) {
  control <- dbartsControl(
    n.trees = 1L,
    n.chains = 1L,
    n.threads = 1L,
    n.samples = 1L,
    n.thin = 1L,
    updateState = FALSE,
    verbose = FALSE,
    keepTrainingFits = TRUE,
    seed = seed,
    proposal.probs = if (arm == "twin") {
      c(
        birth_death = 1,
        swap = 0,
        change = 0,
        perturb = 0,
        rule_gibbs = 0,
        birth = 0.5
      )
    } else {
      dbartsControl()@proposal.probs
    }
  )
  args <- list(
    x,
    yBuild,
    test = xTest,
    leaf.prior = dbartsPriors$normal(2),
    tree.prior = dbartsPriors$cgm(treePower, 0.95),
    sigest = sigest,
    control = control,
    family = dbartsFamilies$gaussian(
      sigma = dbartsPriors$chisq(sigDf, sigQuant)
    )
  )
  if (arm != "twin") {
    args$monotone <- dbartsForests$monotone(
      c(x1 = "increasing"),
      prior = arm
    )
  }
  do.call(dbarts, args)
}

functionals <- function(fTrain, fTest, sigma) {
  c(
    leaves = length(unique(round(fTrain, 9))),
    avg.f = mean(fTrain),
    f.25 = fTest[1L],
    f.50 = fTest[2L],
    f.75 = fTest[3L],
    sigma = sigma,
    wide = fTest[5L] - fTest[4L],
    local = fTest[7L] - fTest[6L]
  )
}

runArm <- function(arm, seed) {
  set.seed(seed)
  sampler <- makeSampler(arm, seed)
  nf <- length(functionals(yBuild, rep(0, nrow(xTest)), 1))
  before <- after <- matrix(NA_real_, nReps, nf)
  for (r in seq_len(nReps)) {
    sampler$sampleTreesFromPrior()
    sampler$sampleLeafParametersFromPrior()
    f0 <- as.numeric(sampler$predict(x))
    f0Test <- as.numeric(sampler$predict(xTest))
    sigma0 <- drawSigma()
    y <- f0 + sigma0 * rnorm(nObs)
    sampler$setSigma(sigma0)
    sampler$setResponse(y)
    res <- sampler$run(nSweeps - 1L, 1L)
    before[r, ] <- functionals(f0, f0Test, sigma0)
    after[r, ] <- functionals(
      res$train[, 1L],
      res$test[, 1L],
      as.numeric(res$sigma)[1L]
    )
  }
  colnames(before) <- colnames(after) <- names(functionals(
    yBuild,
    rep(0, nrow(xTest)),
    1
  ))
  list(before = before, after = after)
}

cat(sprintf(
  paste0(
    "successive-conditional: %d replications x %d sweeps, n %d, tree power ",
    "%.1f, |z| bound %.1f\n"
  ),
  nReps,
  nSweeps,
  nObs,
  treePower,
  zBound
))
failed <- FALSE
for (arm in arms) {
  started <- proc.time()[["elapsed"]]
  res <- runArm(arm, match(arm, c("leaf", "joint", "twin")))
  d <- res$after - res$before
  z <- colMeans(d) / (apply(d, 2L, sd) / sqrt(nrow(d)))
  # a functional the kernel never moves (sd 0) carries no evidence
  z[!is.finite(z)] <- 0
  worst <- max(abs(z))
  fail <- worst > zBound
  failed <- failed || fail
  cat(sprintf(
    "  %-5s %5.0f s  worst |z| %5.2f%s\n",
    arm,
    proc.time()[["elapsed"]] - started,
    worst,
    if (fail) "  <- FAIL" else ""
  ))
  cat(sprintf(
    "        %s\n",
    paste(sprintf("%s %.2f", names(z), z), collapse = ", ")
  ))
}

if (failed) {
  cat(
    "\nFAIL: the kernel does not leave the prior-predictive joint invariant\n"
  )
  quit(status = 1L)
}
cat("\nOK: every arm's paired differences are centred\n")
