#!/usr/bin/env Rscript

# Geweke (2004) joint-distribution test for the gaussian sampler at small
# ensemble size. Two simulators of the SAME joint p(theta, y), theta =
# (forest, sigma):
#
#   marginal-conditional  theta ~ p(theta), then y ~ p(y | theta); iid per draw
#   successive-conditional  theta_t ~ K(. | y_t-1), one dbarts sweep, then
#     y_t ~ p(y | theta_t) installed through the sampler's own between-sweep
#     setResponse
#
# The sweep leaves p(theta | y) invariant and the refresh draws y from the exact
# conditional, so the successive-conditional kernel preserves the joint: started
# at a marginal-conditional draw the chain is stationary FROM STEP 0 and every
# iterate is exactly joint-distributed. That is what makes the arm cheap where
# SBC is not - SBC buys near-independent posterior draws with a 72000-sweep
# burn-in per replication, this buys exactness at every step for nothing - and
# it is what makes the test valid without any convergence claim. Every
# functional of (theta, y) must then agree in distribution between the arms.
# Disagreement means the conditionals the sweep implements and the prior the
# marginal arm draws from are not two halves of one model.
#
# The recorded pair is (theta_t, y_t-1): AFTER the sweep, BEFORE the refresh.
# It is joint-distributed by the same argument and strictly more informative
# than the post-refresh pair, for which y | theta holds by construction. It is
# what lets two PIVOTS into the battery - sqrt(n) mean(y - f) / sigma is exactly
# N(0, 1) and sum((y - f)^2) / sigma^2 exactly chi-square(n) under the joint,
# whatever theta is - which read the sweep's draw of theta against the y that
# generated it, and mix at lag ~1 where every other functional does not.
#
# Self-consistency is the whole game. The forest prior is the ENGINE's own
# (sampleTreesFromPrior + sampleNodeParametersFromPrior, then predict), so a
# mismatch between the tree prior the initialiser grows from and the one the
# Metropolis ratio prices shows up rather than cancelling. sigma has no exposed
# prior draw and is transcribed from the model spec (sbc.R's sbcSigmaDraw:
# scaled-inverse-chi-squared on the REPORTED scale, calibrated so
# P(sigma < sigest) = quant); a wrong transcription would already have shown as
# a non-uniform sigma rank in the gaussian SBC arm, which is the corroboration
# this design leans on. The data scale is fixed at build and
# setResponse(updateScale = FALSE) keeps it, so both arms speak one scale and
# p(theta) does not move with y; the calibration is re-read after the run and
# must be bit-identical.
#
# Two tests per functional. The ecdf-difference simultaneous band (sbc.R
# rankUniformity, Bonferroni'd over the whole battery) reads the rank of each
# chain's FINAL iterate among that chain's OWN fresh block of iid
# marginal-conditional draws: exactly uniform and exactly independent across
# chains under the joint, so no thinning or mixing assumption enters, and
# discrete functionals go through sbc.R's tie-break. The z statistic compares
# marginal-conditional means against successive-conditional chain means over
# each chain's second half, with the standard error taken ACROSS the independent
# chains - no spectral estimate, valid whatever the within-chain correlation
# (reported as tau, the measured variance inflation over an iid run of the same
# length).
#
# A joint test names a disagreement, not a culprit, so a third arm does the
# adjudication the tree-shaped functionals need. Pinning sigma enormous
# (resid.prior = fixed) drives every structure move's marginal-likelihood ratio
# to 1 within ~1e-6, so the sweep's stationary tree law IS the forest prior the
# Metropolis ratio prices, with no response, no sigma draw and no leaf posterior
# in the way. Held against sampleTreesFromPrior's law on the same design, that
# is a prior-versus-prior statement between two ENGINE-side representations of
# one prior, and it separates the two readings of a tree-functional fire: a gap
# there puts the disagreement in the prior representations (the marginal arm's
# theta is then not the theta the sweep integrates over, which is an
# initialiser-consistency finding and not a statement about the posterior), and
# no gap there leaves it in the sweep.
#
# Sensitivity, measured by poisoning the successive-conditional refresh
# (poison=<name>, graded on the sixteen response-side functionals of the two
# arms and inverting the exit status, so a poison that fails to move them is
# itself a failure) at quick size: scaling the refresh noise by 1.02 fires six,
# at z = +23.2 on log.sigma and +15.0 on the log.resid.ss pivot; reading the
# refresh's sigma one step stale fires four, at z = +8.3 on log.resid.ss;
# centring the refresh on 0.99 f fires seven, at z = -14.2 on f.sd, and leaves
# the pivots alone, which is the shape an error in the mean alone has.
# poison=skip is a NEGATIVE control - refreshing only every other step still
# leaves the joint invariant, since a posterior kernel applied twice to one y is
# still a posterior kernel - and it fires nothing at all, so the arms are not
# alarming at any perturbation whatever.
#
# Not seen here: anything the fixed build scale hides (setResponse never
# rescales, so a scale-update path is out of scope); the k hyperprior, held
# fixed at k = 2 because no R-level k prior draw is exposed; weights, offsets
# and test predictors, all absent; families other than gaussian; and a
# state-dependent disagreement between predict() and the recorded training fits,
# which the two arms read separately and a precondition only compares at one
# state. The null self-check shares the marginal-conditional code with the arm
# it validates, so it certifies the rank machinery and the functionals but NOT
# the transcribed sigma prior.
#
# Usage: Rscript geweke-mc.R [quick] [poison=noise|sigmalag|fshrink|skip]

suppressPackageStartupMessages(library(dbarts))

gewekeScriptDir <- function() {
  fileArg <- grep("^--file=", commandArgs(), value = TRUE)
  if (length(fileArg) == 0L) {
    "benchmarks/R"
  } else {
    dirname(sub("^--file=", "", fileArg[1L]))
  }
}
# rankUniformity, sbcDiscreteRank, sbcSigmaDraw, sbcCheckSigmaPrior
source(file.path(gewekeScriptDir(), "sbc.R"))

args <- commandArgs(trailingOnly = TRUE)
quick <- "quick" %in% args
poisonArg <- grep("^poison=", args, value = TRUE)
poison <- if (length(poisonArg) == 0L) {
  "none"
} else {
  sub("^poison=", "", poisonArg[1L])
}
stopifnot(poison %in% c("none", "noise", "sigmalag", "fshrink", "skip"))

kLeaf <- 2
sigDf <- 3
sigQuant <- 0.9
sigest <- 1
zBound <- 4
fitTolerance <- 1e-12
poisonNoise <- 1.02 # refresh noise sd inflation under poison=noise
poisonShrink <- 0.99 # refresh centre under poison=fshrink
# the prior probe's pinned residual sd: 2000 internal-scale units against a
# forest prior sd of 0.25, so the structure moves see a flat likelihood
flatSigma <- 1e4
gewekeSeed <- 20260818L

# chains and steps per chain; the second half of each chain feeds the moment
# test and the final iterate the band. nStep is set at a few times the slowest
# functional's measured autocorrelation time (~70 sweeps for the forest level in
# arm A) so a per-step discrepancy has saturated by the time it is read.
nChain <- if (quick) 300L else 2000L
nStep <- if (quick) 100L else 200L
nPool <- if (quick) 50L else 100L # marginal-conditional draws per rank block
nNull <- if (quick) 300L else 1000L # replications of the null self-check
nProbe <- if (quick) 20000L else 200000L # draws per law in the prior probe
nBatch <- 100L # batches the prior probe's sweep series is read through

## The two configurations. Small n keeps the successive-conditional chain
## mixing (a tight posterior pins y to f and the pair random-walks) and is what
## gives the tree-shaped functionals their resolution: two forest priors can
## only be told apart where a prior-grown split is able to land a data-sparse
## child, and by 1600 observations at five trees the leaf-count gap is under
## 0.05% and invisible in 30000 draws. Small m is the mandated regime; arm B
## adds unused predictors and a design sparse enough that the cutpoint grid runs
## out inside the tree prior.
gewekeArms <- list(
  list(label = "A: n=100 p=3 m=8", n = 100L, p = 3L, m = 8L, designSeed = 11L),
  list(label = "B: n=30 p=5 m=5", n = 30L, p = 5L, m = 5L, designSeed = 12L)
)

gewekeFunctionalNames <- c(
  "f.mean",
  "f.row1",
  "f.sd",
  "log.sigma",
  "y.mean",
  "log.y.sd",
  "resid.mean",
  "log.resid.ss",
  "n.leaves",
  "leaf.depth",
  "splits.x1"
)
nFunctional <- length(gewekeFunctionalNames)
# every functional of every arm is one look; the band's level is corrected for
# all of them so a whole-battery pass has probability ~0.95
bandAlpha <- 0.05 / (nFunctional * length(gewekeArms))

## Fixed design and build response for one arm. The build response is
## deterministic: it fixes the internal scale (range 5, centre 0) once, and
## every later setResponse leaves it alone.
gewekeDesign <- function(arm) {
  set.seed(arm$designSeed)
  x <- matrix(runif(arm$n * arm$p), arm$n, arm$p)
  colnames(x) <- paste0("x", seq_len(arm$p))
  arm$x <- x
  arm$yBuild <- seq(-2.5, 2.5, length.out = arm$n)
  arm
}

## The sampler both arms share. flat = TRUE instead pins sigma at flatSigma for
## the prior probe, which needs the likelihood out of the structure moves.
gewekeSampler <- function(arm, flat = FALSE) {
  control <- dbartsControl(
    n.trees = arm$m,
    n.chains = 1L,
    n.threads = 1L,
    n.samples = 1L,
    n.thin = 1L,
    updateState = FALSE,
    verbose = FALSE,
    keepTrainingFits = TRUE
  )
  dbarts(
    arm$x,
    arm$yBuild,
    resid.prior = if (flat) {
      dbartsPriors$fixed(flatSigma^2)
    } else {
      dbartsPriors$chisq(sigDf, sigQuant)
    },
    node.prior = dbartsPriors$normal(kLeaf),
    sigma = if (flat) flatSigma else sigest,
    control = control
  )
}

## Depth of every node of a pre-order forest block, resetting at each tree.
## Leaves carry var == -1 and every internal node has exactly two children, so
## one pass with an explicit stack of owed children settles it.
gewekeDepth <- function(isLeaf, tree) {
  depths <- integer(length(isLeaf))
  owed <- integer(length(isLeaf))
  top <- 0L
  depth <- 0L
  current <- -1L
  for (i in seq_along(isLeaf)) {
    if (tree[i] != current) {
      current <- tree[i]
      top <- 0L
      depth <- 0L
    }
    depths[i] <- depth
    if (isLeaf[i]) {
      while (top > 0L) {
        owed[top] <- owed[top] - 1L
        if (owed[top] > 0L) {
          break
        }
        top <- top - 1L
        depth <- depth - 1L
      }
    } else {
      top <- top + 1L
      owed[top] <- 2L
      depth <- depth + 1L
    }
  }
  depths
}

## The tree-shaped tail of the battery, split out so the prior probe measures
## exactly the quantities the joint arms do.
gewekeTreeStats <- function(trees) {
  isLeaf <- trees$var == -1L
  depths <- gewekeDepth(isLeaf, trees$tree)
  c(sum(isLeaf), mean(depths[isLeaf]), sum(trees$var == 1L))
}

## The battery, computed by ONE function from the same primitives in both arms:
## the fit vector, sigma, the response, and the live forest. Scale functionals
## are logged where the scaled-inverse-chi-squared tail would otherwise cost the
## moment test its finite variance; f.sd is left raw because an all-stump forest
## makes it exactly 0.
gewekeFunctionals <- function(f, sigma, y, trees, n) {
  resid <- y - f
  c(
    mean(f),
    f[1L],
    sd(f),
    log(sigma),
    mean(y),
    log(sd(y)),
    sqrt(n) * mean(resid) / sigma,
    log(sum(resid^2) / sigma^2),
    gewekeTreeStats(trees)
  )
}

## One marginal-conditional draw: forest from the engine's own prior, sigma from
## the transcribed prior, y from the gaussian likelihood at that theta.
gewekeMarginalDraw <- function(sampler, arm, drawSigma) {
  sampler$sampleTreesFromPrior()
  sampler$sampleNodeParametersFromPrior()
  f <- as.numeric(sampler$predict(arm$x))
  sigma <- drawSigma(1L)
  y <- f + sigma * rnorm(arm$n)
  gewekeFunctionals(f, sigma, y, sampler$getTrees(current = TRUE), arm$n)
}

gewekeMarginalPool <- function(sampler, arm, drawSigma, nDraw) {
  pool <- matrix(NA_real_, nDraw, nFunctional)
  for (i in seq_len(nDraw)) {
    pool[i, ] <- gewekeMarginalDraw(sampler, arm, drawSigma)
  }
  pool
}

## One successive-conditional chain, initialised at a marginal-conditional draw
## so that it is stationary immediately. Returns the recorded second half; its
## last row is the chain's final iterate. Poisons act on the refresh only, never
## on the sweep, so a fire under one is the oracle's two arms disagreeing about
## the same engine rather than a changed engine.
gewekeChain <- function(sampler, arm, drawSigma) {
  sampler$sampleTreesFromPrior()
  sampler$sampleNodeParametersFromPrior()
  f <- as.numeric(sampler$predict(arm$x))
  sigma <- drawSigma(1L)
  y <- f + sigma * rnorm(arm$n)
  sampler$setSigma(sigma)
  sampler$setResponse(y)

  firstKept <- nStep %/% 2L + 1L
  kept <- matrix(NA_real_, nStep - firstKept + 1L, nFunctional)
  lastSigma <- sigma
  for (t in seq_len(nStep)) {
    result <- sampler$run(0L, 1L)
    f <- result$train[, 1L]
    sigma <- result$sigma[1L]
    if (t >= firstKept) {
      kept[t - firstKept + 1L, ] <- gewekeFunctionals(
        f,
        sigma,
        y,
        sampler$getTrees(current = TRUE),
        arm$n
      )
    }
    if (!(poison == "skip" && t %% 2L == 0L)) {
      noiseSd <- switch(
        poison,
        noise = poisonNoise * sigma,
        sigmalag = lastSigma,
        sigma
      )
      centre <- if (poison == "fshrink") poisonShrink * f else f
      y <- centre + noiseSd * rnorm(arm$n)
      sampler$setResponse(y)
    }
    lastSigma <- sigma
  }
  kept
}

## The deterministic precondition, which fires before any distributional arm
## could. The arms read f through different accessors - predict() in the
## marginal arm, which has no run to read, and the recorded training fits in the
## chain - so those must be the same map, or the refresh is not conditioning on
## the f the likelihood uses.
gewekeFitGap <- function(sampler, arm, drawSigma) {
  sampler$sampleTreesFromPrior()
  sampler$sampleNodeParametersFromPrior()
  f <- as.numeric(sampler$predict(arm$x))
  sampler$setResponse(f + drawSigma(1L) * rnorm(arm$n))
  result <- sampler$run(0L, 1L)
  max(abs(result$train[, 1L] - as.numeric(sampler$predict(arm$x))))
}

## Ranks of each chain's final iterate in that chain's own block of the pool,
## plus the across-chain z. Every rank must be taken before any rankUniformity
## call, which re-seeds the stream for its own null simulation. nAverage is the
## number of iterates behind each row of chainMeans, and turns the chain-mean
## variance into the reported tau.
gewekeCompare <- function(pool, chainMeans, finals, nAverage) {
  ranks <- matrix(NA_integer_, nrow(finals), nFunctional)
  for (i in seq_len(nrow(finals))) {
    block <- pool[((i - 1L) * nPool + 1L):(i * nPool), , drop = FALSE]
    for (j in seq_len(nFunctional)) {
      ranks[i, j] <- sbcDiscreteRank(block[, j], finals[i, j])
    }
  }
  poolVar <- apply(pool, 2L, var)
  chainVar <- apply(chainMeans, 2L, var)
  list(
    ranks = ranks,
    poolMean = colMeans(pool),
    chainMean = colMeans(chainMeans),
    z = (colMeans(chainMeans) - colMeans(pool)) /
      sqrt(chainVar / nrow(chainMeans) + poolVar / nrow(pool)),
    tau = chainVar * nAverage / poolVar
  )
}

## The adjudication arm: sampleTreesFromPrior's law against the law the sweep
## holds in stationarity once the likelihood is flat. Both are the engine's, so
## a gap is an inconsistency between two representations of the forest prior and
## nothing to do with the response. The sweep series is read through batches
## long against its ~20-sweep autocorrelation time, so the batch means are
## near-independent and their spread is the standard error.
gewekePriorProbe <- function(arm) {
  sampler <- gewekeSampler(arm, flat = TRUE)
  stopifnot(sampler$run(0L, 1L)$sigma[1L] == flatSigma)
  priorDraws <- matrix(NA_real_, nProbe, 3L)
  for (i in seq_len(nProbe)) {
    sampler$sampleTreesFromPrior()
    sampler$sampleNodeParametersFromPrior()
    priorDraws[i, ] <- gewekeTreeStats(sampler$getTrees(current = TRUE))
  }
  sweepDraws <- matrix(NA_real_, nProbe, 3L)
  for (i in seq_len(nProbe)) {
    sampler$run(0L, 1L)
    sweepDraws[i, ] <- gewekeTreeStats(sampler$getTrees(current = TRUE))
  }
  batches <- vapply(
    split(seq_len(nProbe), rep(seq_len(nBatch), each = nProbe %/% nBatch)),
    function(rows) colMeans(sweepDraws[rows, , drop = FALSE]),
    numeric(3L)
  )
  list(
    priorMean = colMeans(priorDraws),
    sweepMean = colMeans(sweepDraws),
    z = (colMeans(sweepDraws) - colMeans(priorDraws)) /
      sqrt(
        apply(batches, 1L, var) / nBatch + apply(priorDraws, 2L, var) / nProbe
      )
  )
}

gewekeColumns <- function(first, second) {
  cat(sprintf(
    "%-13s %9s %10s %7s %6s %8s %8s %7s\n",
    "functional",
    first,
    second,
    "z",
    "tau",
    "ecdfDiff",
    "band",
    "verdict"
  ))
}

gewekeReportProbe <- function(label, probe, elapsed) {
  cat(sprintf(
    "\n%s prior probe (sigma pinned at %g) | %d draws per law | %.1fs\n",
    label,
    flatSigma,
    nProbe,
    elapsed
  ))
  gewekeColumns("initialiser", "stationary")
  fired <- logical(3L)
  for (j in seq_len(3L)) {
    bad <- abs(probe$z[j]) > zBound
    fired[j] <- bad
    cat(sprintf(
      "%-13s %9.4f %10.4f %+7.2f %6s %8s %8s %7s\n",
      gewekeFunctionalNames[nFunctional - 3L + j],
      probe$priorMean[j],
      probe$sweepMean[j],
      probe$z[j],
      "-",
      "-",
      "-",
      if (bad) "FIRE" else "ok"
    ))
  }
  fired
}

gewekeReport <- function(label, comparison, steps, elapsed) {
  cat(sprintf(
    "\n%s | %d chains x %d steps, %d marginal draws | %.1fs\n",
    label,
    nrow(comparison$ranks),
    steps,
    nrow(comparison$ranks) * nPool,
    elapsed
  ))
  gewekeColumns("marginal", "successive")
  fired <- logical(nFunctional)
  for (j in seq_len(nFunctional)) {
    uniformity <- rankUniformity(
      comparison$ranks[, j],
      nPool,
      alpha = bandAlpha
    )
    bad <- !isTRUE(uniformity$pass) || abs(comparison$z[j]) > zBound
    fired[j] <- bad
    cat(sprintf(
      "%-13s %9.4f %10.4f %+7.2f %6.1f %8.4f %8.4f %7s\n",
      gewekeFunctionalNames[j],
      comparison$poolMean[j],
      comparison$chainMean[j],
      comparison$z[j],
      comparison$tau[j],
      uniformity$ecdfDiff,
      uniformity$ecdfBand,
      if (bad) "FIRE" else "ok"
    ))
  }
  fired
}

cat(sprintf(
  "Geweke marginal-conditional check | %s | band alpha %.5f, |z| < %g\n",
  if (quick) "quick" else "full",
  bandAlpha,
  zBound
))
if (poison != "none") {
  cat(sprintf(
    "POISON MODE (%s): an instrument, not a gate; here a FIRE is the pass\n",
    poison
  ))
}

set.seed(gewekeSeed)
drawSigma <- sbcSigmaDraw(sigest, sigDf, sigQuant)
sigmaCheck <- sbcCheckSigmaPrior(sigest, sigDf, sigQuant)
cat(sprintf(
  "sigma prior: P(sigma < sigest) %.4f vs %.4f; median %.4f vs %.4f -> %s\n",
  sigmaCheck$coverage,
  sigmaCheck$coverageTarget,
  sigmaCheck$medianEmpirical,
  sigmaCheck$medianTheory,
  if (isTRUE(sigmaCheck$pass)) "PASS" else "FAIL"
))
if (!isTRUE(sigmaCheck$pass)) {
  stop("transcribed sigma prior is miscalibrated; the marginal arm is invalid")
}

# a poisoned refresh is graded on the response side of the battery; the three
# tree-shaped functionals are the prior probe's business and the refresh is not
# in their path
responseIndex <- seq_len(nFunctional - 3L)
armFire <- FALSE
responseFire <- FALSE
probeFire <- FALSE
for (armIndex in seq_along(gewekeArms)) {
  arm <- gewekeDesign(gewekeArms[[armIndex]])
  set.seed(gewekeSeed + armIndex)
  generator <- gewekeSampler(arm)
  chainSampler <- gewekeSampler(arm)
  calibration <- chainSampler$getCalibration()
  # one scale, one leaf prior, no k hyperprior: what the marginal arm assumes
  stopifnot(
    identical(calibration, generator$getCalibration()),
    unname(calibration[1L, "k"]) == kLeaf,
    unname(calibration[1L, "k.has.hyperprior"]) == 0
  )

  fitGap <- gewekeFitGap(generator, arm, drawSigma)
  cat(sprintf(
    "\n%s precondition: predict vs recorded fits %.3g -> %s\n",
    arm$label,
    fitGap,
    if (fitGap < fitTolerance) "PASS" else "FAIL"
  ))
  if (!(fitGap < fitTolerance)) {
    stop(
      "predict() and the recorded fits disagree, so the two arms read ",
      "different f and the comparison is meaningless"
    )
  }

  started <- proc.time()[["elapsed"]]
  pool <- gewekeMarginalPool(generator, arm, drawSigma, nChain * nPool)
  chainMeans <- matrix(NA_real_, nChain, nFunctional)
  finals <- matrix(NA_real_, nChain, nFunctional)
  for (i in seq_len(nChain)) {
    kept <- gewekeChain(chainSampler, arm, drawSigma)
    chainMeans[i, ] <- colMeans(kept)
    finals[i, ] <- kept[nrow(kept), ]
  }
  elapsed <- proc.time()[["elapsed"]] - started

  # the scale both arms speak must not have moved under the run's setResponse
  if (!identical(chainSampler$getCalibration(), calibration)) {
    stop("the build calibration moved during the run; p(theta) is not fixed")
  }

  fired <- gewekeReport(
    arm$label,
    gewekeCompare(pool, chainMeans, finals, nStep %/% 2L),
    nStep,
    elapsed
  )
  armFire <- any(fired) || armFire
  responseFire <- any(fired[responseIndex]) || responseFire

  # the refresh is not in the probe's path, so a poisoned run has nothing to
  # learn from it
  if (poison == "none") {
    started <- proc.time()[["elapsed"]]
    probe <- gewekePriorProbe(arm)
    probeFire <- any(gewekeReportProbe(
      arm$label,
      probe,
      proc.time()[["elapsed"]] - started
    )) ||
      probeFire
  }
}

## Null self-check: the successive-conditional arm replaced by marginal-
## conditional draws from the OTHER sampler, which is the nStep = 0 limit of the
## chain. Ranks are uniform there by construction, so a fire is the rank
## machinery, the tie-break or the functionals, never the sweep.
arm <- gewekeDesign(gewekeArms[[1L]])
set.seed(gewekeSeed + 101L)
generator <- gewekeSampler(arm)
chainSampler <- gewekeSampler(arm)
started <- proc.time()[["elapsed"]]
nullPool <- gewekeMarginalPool(generator, arm, drawSigma, nNull * nPool)
nullDraws <- gewekeMarginalPool(chainSampler, arm, drawSigma, nNull)
nullFire <- any(gewekeReport(
  "null self-check (both arms marginal)",
  gewekeCompare(nullPool, nullDraws, nullDraws, 1L),
  0L,
  proc.time()[["elapsed"]] - started
))
if (nullFire) {
  stop(
    "null self-check fired: the rank machinery or the functionals are at ",
    "fault, so the arms above cannot be adjudicated"
  )
}

if (poison != "none") {
  # skip is the negative control: it leaves the joint invariant and must NOT
  # move the pivots
  if (responseFire == (poison != "skip")) {
    cat(sprintf(
      "\nOK: the response-side functionals %s under poison=%s, as they must\n",
      if (responseFire) "fired" else "held",
      poison
    ))
    quit(status = 0L)
  }
  cat(sprintf(
    "\nFAIL: the response-side functionals %s under poison=%s\n",
    if (responseFire) "fired" else "held",
    poison
  ))
  quit(status = 1L)
}
if (probeFire) {
  cat(
    "\nFAIL: sampleTreesFromPrior's forest law and the law the sweep holds\n",
    "under a flat likelihood are not the same prior, so a tree-shaped fire\n",
    "above is that inconsistency and not the sweep\n",
    sep = ""
  )
  quit(status = 1L)
}
if (armFire) {
  cat("\nFAIL: the two joint simulators disagree\n")
  quit(status = 1L)
}
cat("\nOK: marginal-conditional and successive-conditional agree\n")
