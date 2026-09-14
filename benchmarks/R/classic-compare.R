#!/usr/bin/env Rscript

# Classic comparison: does 1.0-0 compute what 0.9-34 computed, wherever both
# can fit the same model? This is the wide sibling of R/equivalence.R's
# nine-scenario cross-engine record, and it runs the SAME script under two
# installed libraries rather than recording from one and comparing from the
# other: 0.9-34 is a released package, not a build of this tree, so neither
# side can host the other.
#
# Every scenario is written against the 0.9-x API - bartBT's BayesTree
# spelling (which 0.9-34 spells 'bart'), dbarts() + $run(), xbart() - and
# pins every prior and control setting whose DEFAULT moved between the two
# releases: n.trees, n.burn, n.samples, n.chains, k, power, base, sigdf,
# sigquant, n.cuts, and the tree-move mixture, which 1.0-0 ships at
# birth/death 0.6, swap 0, change 0.4 where 0.9-x shipped 0.5 / 0.1 / 0.4.
# 'mixture=classic' sets the 0.9-x mixture on both sides, so the comparison
# isolates engine agreement; 'mixture=default' leaves it unset, which under
# 1.0-0 is the new kernel and under 0.9-34 is the old one, and so measures
# the mixture change instead.
#
# Data are fixed per scenario; only the MCMC seed varies (20 seeds), so the
# spread across seeds is pure Monte Carlo variability and a per-summary Welch
# z between the two libraries should look standard normal when the posteriors
# agree. Disjoint seed ranges are the exact, parameter-free complement, for
# the case where one side diverges and degenerates the z.
#
# Usage:
#   R_LIBS=<lib> Rscript classic-compare.R record out.rds [mixture=classic]
#   Rscript classic-compare.R compare a.rds b.rds
#   Rscript classic-compare.R merge out.rds part1.rds part2.rds ...
# Append 'quick' for a fast smoke run (not comparable to a full one).
# CLASSIC_COMPARE_SCENARIOS (comma-separated) restricts a run to a subset, so
# a long record can be taken in chunks and joined with 'merge'.
# CLASSIC_COMPARE_SEED_OFFSET shifts the seed block, for re-checking a
# marginal flag on fresh draws; CLASSIC_COMPARE_CORES sets the worker count;
# CLASSIC_COMPARE_TABLE names a csv for 'compare' to write its table to.

suppressPackageStartupMessages(library(dbarts))

args <- commandArgs(trailingOnly = TRUE)
quick <- "quick" %in% args
args <- setdiff(args, "quick")
mixture <- "classic"
mixtureArg <- grep("^mixture=", args, value = TRUE)
if (length(mixtureArg) > 0L) {
  mixture <- sub("^mixture=", "", mixtureArg[[1L]])
  args <- setdiff(args, mixtureArg)
}
mode <- if (length(args) >= 1L) args[[1L]] else "record"

if (!mixture %in% c("classic", "default")) {
  stop("mixture must be 'classic' or 'default'")
}

n.seeds <- if (quick) 3L else 20L
# A marginal flag on one summary out of a thousand is what a nominal-level
# statistical gate looks like when it is working; this shifts the seed block
# so a flagged scenario can be re-run on a fresh set of draws. It rides the
# recorded settings, so two result sets from different blocks refuse to
# compare.
seedOffset <- local({
  env <- Sys.getenv("CLASSIC_COMPARE_SEED_OFFSET", "0")
  as.integer(env)
})
ndpost <- if (quick) 200L else 1000L
nskip <- if (quick) 100L else 500L
ntree <- if (quick) 50L else 200L
n.test <- 25L
n.train.summ <- 25L

# the 0.9-x mixture, written the way both releases' front doors take it;
# NULL is each release's own default, which is the point of 'mixture=default'
classicProposalProbs <- c(
  birth_death = 0.5,
  swap = 0.1,
  change = 0.4,
  birth = 0.5
)
proposalProbs <- if (mixture == "classic") classicProposalProbs else NULL

# 0.9-34 spells the BayesTree door 'bart'; 1.0-0 spells it 'bartBT' and
# forwards a BayesTree-shaped 'bart' call to it. Calling the successor by
# name where it exists keeps the deprecation path out of the measurement.
bartFn <- if (
  exists("bartBT", where = asNamespace("dbarts"), inherits = FALSE)
) {
  dbarts::bartBT
} else {
  dbarts::bart
}
# 0.9-34's dbartsControl has no mixture slot (the model object carried it),
# so xbart - which exposes no flat knob on either release - can only be set
# through the control under 1.0-0.
controlTakesProposalProbs <- "proposal.probs" %in%
  names(formals(dbartsControl))

friedman <- function(x) {
  10 *
    sin(pi * x[, 1L] * x[, 2L]) +
    20 * (x[, 3L] - 0.5)^2 +
    10 * x[, 4L] +
    5 * x[, 5L]
}

# Warnings both releases raise for reasons that are part of the scenario,
# plus 1.0-0's rename notices, which are not a difference in what is computed.
muffleBenignWarning <- function(w) {
  msg <- conditionMessage(w)
  if (
    grepl("'weights' are ignored for test data", msg) ||
      grepl("'weights' of 0 will be ignored", msg) ||
      grepl("columns of 'test' will be matched by position", msg) ||
      grepl("weights specified but not found in test data", msg) ||
      grepl("deprecated", msg) ||
      grepl("proposal.probs", msg) ||
      grepl("bartBT", msg)
  ) {
    invokeRestart("muffleWarning")
  }
}

quietly <- function(expr) {
  withCallingHandlers(expr, warning = muffleBenignWarning)
}

makeScenarios <- function() {
  result <- list()
  p <- 10L

  # --- the BayesTree door ------------------------------------------------
  set.seed(6101L)
  x <- matrix(runif(500L * p), 500L)
  result$friedman <- list(
    kind = "bart",
    x = x,
    y = friedman(x) + rnorm(500L),
    x.test = matrix(runif(n.test * p), n.test),
    binary = FALSE
  )

  set.seed(6102L)
  x <- matrix(runif(1000L * p), 1000L)
  result$probit <- list(
    kind = "bart",
    x = x,
    y = rbinom(1000L, 1L, pnorm(scale(friedman(x)))),
    x.test = matrix(runif(n.test * p), n.test),
    binary = TRUE
  )

  # weights on the 0.9-x scale: a half-weight row is worth half an
  # observation, a double-weight row two
  set.seed(6103L)
  x <- matrix(runif(500L * p), 500L)
  weights <- sample(c(0.5, 1, 2), 500L, replace = TRUE)
  result$weighted <- list(
    kind = "bart",
    x = x,
    y = friedman(x) + rnorm(500L) / sqrt(weights),
    weights = weights,
    x.test = matrix(runif(n.test * p), n.test),
    binary = FALSE
  )

  # a zero weight is a row that is carried but contributes nothing
  set.seed(6104L)
  x <- matrix(runif(500L * p), 500L)
  weights <- rep(1, 500L)
  weights[sample.int(500L, 75L)] <- 0
  result$zeroweights <- list(
    kind = "bart",
    x = x,
    y = friedman(x) + rnorm(500L),
    weights = weights,
    x.test = matrix(runif(n.test * p), n.test),
    binary = FALSE
  )

  set.seed(6105L)
  x <- matrix(runif(1000L * p), 1000L)
  result$binaryoffset <- list(
    kind = "bart",
    x = x,
    y = rbinom(1000L, 1L, pnorm(scale(friedman(x)) + 0.4)),
    x.test = matrix(runif(n.test * p), n.test),
    binaryOffset = 0.4,
    binary = TRUE
  )

  set.seed(6106L)
  x <- matrix(runif(500L * p), 500L)
  result$fixedk <- list(
    kind = "bart",
    x = x,
    y = friedman(x) + rnorm(500L),
    x.test = matrix(runif(n.test * p), n.test),
    k = 3.0,
    binary = FALSE
  )

  set.seed(6107L)
  x <- matrix(runif(500L * p), 500L)
  y <- friedman(x) + rnorm(500L)
  x.test <- matrix(runif(n.test * p), n.test)
  result$cutssmall <- list(
    kind = "bart",
    x = x,
    y = y,
    x.test = x.test,
    numcut = 5L,
    binary = FALSE
  )
  result$cutslarge <- list(
    kind = "bart",
    x = x,
    y = y,
    x.test = x.test,
    numcut = 1000L,
    binary = FALSE
  )
  # quantile cut points over the same data, so the cut rule is the only change
  result$usequants <- list(
    kind = "bart",
    x = x,
    y = y,
    x.test = x.test,
    numcut = 20L,
    usequants = TRUE,
    binary = FALSE
  )

  # a factor predictor, twice: once as a data frame column (the BayesTree
  # door expands it to indicators on both releases) and once as the
  # indicator matrix written out by hand, which must reach the same model
  set.seed(6108L)
  n <- 600L
  xnum <- matrix(runif(n * 5L), n)
  g <- factor(sample(letters[1L:4L], n, replace = TRUE))
  gEffect <- c(a = -2, b = 0, c = 1.5, d = 4)[as.character(g)]
  yf <- friedman(cbind(xnum, xnum[, 1L:5L])) + gEffect + rnorm(n)
  frame <- data.frame(xnum, g = g)
  frameTest <- frame[seq_len(n.test), , drop = FALSE]
  indicators <- cbind(
    xnum,
    stats::model.matrix(~ g - 1)[,, drop = FALSE]
  )
  colnames(indicators) <- c(paste0("X", 1L:5L), paste0("g.", levels(g)))
  result$factorframe <- list(
    kind = "bart",
    x = frame,
    y = yf,
    x.test = frameTest,
    binary = FALSE
  )
  result$factorindicators <- list(
    kind = "bart",
    x = indicators,
    y = yf,
    x.test = indicators[seq_len(n.test), , drop = FALSE],
    binary = FALSE
  )

  set.seed(6109L)
  n <- 600L
  xnum <- matrix(runif(n * 5L), n)
  o <- factor(sample(1L:5L, n, replace = TRUE), ordered = TRUE)
  yo <- friedman(cbind(xnum, xnum[, 1L:5L])) +
    2 * as.integer(o) +
    rnorm(n)
  frame <- data.frame(xnum, o = o)
  result$orderedfactor <- list(
    kind = "bart",
    x = frame,
    y = yo,
    x.test = frame[seq_len(n.test), , drop = FALSE],
    binary = FALSE
  )

  set.seed(6110L)
  x <- matrix(runif(500L * p), 500L)
  y <- friedman(x) + rnorm(500L)
  x.test <- matrix(runif(n.test * p), n.test)
  result$trees20 <- list(
    kind = "bart",
    x = x,
    y = y,
    x.test = x.test,
    ntree = 20L,
    binary = FALSE
  )
  result$trees200 <- list(
    kind = "bart",
    x = x,
    y = y,
    x.test = x.test,
    ntree = 200L,
    binary = FALSE
  )
  result$singletree <- list(
    kind = "bart",
    x = x,
    y = y,
    x.test = x.test,
    ntree = 1L,
    binary = FALSE
  )

  set.seed(6111L)
  x <- matrix(runif(5000L * p), 5000L)
  result$largen <- list(
    kind = "bart",
    x = x,
    y = friedman(x) + rnorm(5000L),
    x.test = matrix(runif(n.test * p), n.test),
    binary = FALSE
  )

  set.seed(6112L)
  x <- matrix(runif(500L * p), 500L)
  result$chains4 <- list(
    kind = "bart",
    x = x,
    y = friedman(x) + rnorm(500L),
    x.test = matrix(runif(n.test * p), n.test),
    nchain = 4L,
    binary = FALSE
  )

  # a test set drawn wider than the training range, so the test rows are
  # answered by the edge leaves rather than interpolated
  set.seed(6113L)
  x <- matrix(runif(500L * p), 500L)
  result$testset <- list(
    kind = "bart",
    x = x,
    y = friedman(x) + rnorm(500L),
    x.test = matrix(runif(n.test * p, -0.25, 1.25), n.test),
    binary = FALSE
  )

  # --- the sampler door --------------------------------------------------
  set.seed(6114L)
  x <- matrix(runif(500L * p), 500L)
  x.test <- matrix(runif(n.test * p), n.test)
  offset <- 2 * x[, 6L] - 1
  result$offset <- list(
    kind = "sampler",
    x = x,
    y = friedman(x) + offset + rnorm(500L),
    offset = offset,
    offset.test = 2 * x.test[, 6L] - 1,
    x.test = x.test,
    binary = FALSE
  )

  # k drawn under the chi hyperprior at 0.9-x's OWN default shape
  # (degreesOfFreedom 1.25, scale Inf); 1.0-0's chi() defaults differ, so the
  # spec is written out rather than left to either release's formals. The
  # string form resolves through each release's own normal().
  set.seed(6115L)
  x <- matrix(runif(500L * p), 500L)
  result$chik <- list(
    kind = "sampler",
    x = x,
    y = friedman(x) + rnorm(500L),
    x.test = matrix(runif(n.test * p), n.test),
    nodePriorSpec = "chi(1.25, Inf)",
    binary = FALSE
  )

  # a 0.9-x embedded-Gibbs loop: the response and the offset are both
  # replaced between sweeps, which is the pattern dbartsSampler exists for
  set.seed(6116L)
  x <- matrix(runif(500L * p), 500L)
  base <- friedman(x)
  result$gibbsloop <- list(
    kind = "gibbs",
    x = x,
    y = base + rnorm(500L),
    x.test = matrix(runif(n.test * p), n.test),
    sweeps = lapply(seq_len(5L), function(i) {
      list(
        offset = rep(0.2 * i, 500L),
        y = base + 0.2 * i + rnorm(500L)
      )
    }),
    binary = FALSE
  )

  # a mid-chain predictor swap through the whole matrix
  set.seed(6117L)
  x <- matrix(runif(500L * p), 500L)
  x2 <- x
  x2[, 1L:3L] <- matrix(runif(500L * 3L), 500L)
  result$setpredictor <- list(
    kind = "setpredictor",
    x = x,
    y = friedman(x) + rnorm(500L),
    x.swap = x2,
    x.test = matrix(runif(n.test * p), n.test),
    binary = FALSE
  )

  # --- crossvalidation ---------------------------------------------------
  # 0.9-x's xbart carried one chain across replications and burned it in
  # again per replication (its n.burn's third element); 1.0-0 starts each
  # replication fresh and refuses that third element. The two arms separate
  # that documented change from the engine: 'xbart' runs several
  # replications, where the carry is in force on the old side only, and
  # 'xbart1rep' runs one, where there is no between-replication chain to
  # carry and the two calls describe the same procedure.
  set.seed(6118L)
  x <- matrix(runif(400L * p), 400L)
  y <- friedman(x) + rnorm(400L)
  result$xbart <- list(
    kind = "xbart",
    x = x,
    y = y,
    n.reps = if (quick) 2L else 8L,
    binary = FALSE
  )
  result$xbart1rep <- list(
    kind = "xbart",
    x = x,
    y = y,
    n.reps = 1L,
    binary = FALSE
  )
  # the same crossvalidation written out by hand over the same grid, each
  # fold its own fresh fit through the BayesTree door. It uses nothing but
  # the 0.9-x API, so it is an independent anchor for the two xbart rows:
  # whichever release's xbart agrees with its OWN hand-rolled CV is the one
  # whose crossvalidation loop reports the loss the fit actually earns.
  result$cvbyhand <- list(
    kind = "cv",
    x = x,
    y = y,
    binary = FALSE
  )

  # --- the tree-structure posterior ---------------------------------------
  # Three predictors with a hundred cut points against three with one, and a
  # response that is pure noise, so nothing in the data prefers one predictor
  # over another and the split counts read the tree-structure posterior
  # directly. That posterior is what the change move's repaired acceptance
  # ratio governs, and the unequal cut counts are the condition under which
  # the old ratio's missing proposal-density term fails to cancel. Every
  # other scenario here gives every predictor the same number of cut points,
  # where the term cancels and the repair is invisible.
  set.seed(6119L)
  n <- 500L
  xmixed <- cbind(
    matrix(runif(n * 3L), n),
    matrix(rbinom(n * 3L, 1L, 0.5), n)
  )
  colnames(xmixed) <- c("c1", "c2", "c3", "b1", "b2", "b3")
  xmixedTest <- cbind(
    matrix(runif(n.test * 3L), n.test),
    matrix(rbinom(n.test * 3L, 1L, 0.5), n.test)
  )
  colnames(xmixedTest) <- colnames(xmixed)
  result$mixedcuts <- list(
    kind = "bart",
    x = xmixed,
    y = rnorm(n),
    x.test = xmixedTest,
    binary = FALSE
  )

  # a cutoff for the nonlinear functional, fixed with the data: the fraction
  # of fitted values above it is a posterior quantity no linear summary of
  # the fits pins down
  for (name in names(result)) {
    scn <- result[[name]]
    if (is.null(scn[["cutoff"]])) {
      result[[name]]$cutoff <- if (isTRUE(scn[["binary"]])) {
        0
      } else {
        as.numeric(stats::median(scn[["y"]]))
      }
    }
  }
  result
}

# [d, S] or [d, S, C] -> (S * C) x d, pooling chains
poolChains <- function(a) {
  if (length(dim(a)) == 3L) t(matrix(a, nrow = dim(a)[1L])) else t(a)
}

makeControl <- function(scn, nChains = 1L) {
  dbartsControl(
    n.chains = nChains,
    n.threads = 1L,
    n.trees = if (!is.null(scn[["ntree"]])) scn[["ntree"]] else ntree,
    n.cuts = if (!is.null(scn[["numcut"]])) scn[["numcut"]] else 100L,
    n.thin = 1L,
    n.burn = 0L,
    useQuantiles = isTRUE(scn[["usequants"]]),
    keepTrainingFits = TRUE,
    keepTrees = FALSE,
    printEvery = ndpost + nskip + 1L,
    updateState = FALSE
  )
}

# Every prior is named, none left to a release's formals.
makeSampler <- function(scn, nChains = 1L) {
  nodePrior <- if (!is.null(scn[["nodePriorSpec"]])) {
    call("normal", scn[["nodePriorSpec"]])
  } else {
    call("normal", if (!is.null(scn[["k"]])) scn[["k"]] else 2.0)
  }
  samplerCall <- as.call(c(
    list(
      quote(dbarts),
      quote(scn[["x"]]),
      quote(scn[["y"]]),
      test = quote(scn[["x.test"]]),
      control = quote(makeControl(scn, nChains)),
      tree.prior = quote(cgm(2.0, 0.95)),
      node.prior = nodePrior,
      family = quote(gaussian(sigma = chisq(3.0, 0.90))),
      sigma = NA_real_
    ),
    # omitted, not passed as NULL, when the mixture is left to the release's
    # own default: 1.0-0 reads this name off '...', where a NULL is a value
    if (!is.null(proposalProbs)) list(proposal.probs = quote(proposalProbs)),
    if (!is.null(scn[["weights"]])) list(weights = quote(scn[["weights"]])),
    if (!is.null(scn[["offset"]])) list(offset = quote(scn[["offset"]])),
    if (!is.null(scn[["offset.test"]])) {
      list(offset.test = quote(scn[["offset.test"]]))
    }
  ))
  quietly(eval(samplerCall))
}

# The BayesTree door reads 'k' off its own matched call (both releases do),
# so the call is built with the value already in place rather than with an
# expression that would reach normal() unevaluated. keepcall = FALSE keeps
# the predictor matrices out of the stored call, and takes both releases down
# the same "no matched call, use the argument" branch for k.
fitViaBart <- function(scn) {
  bartArgs <- list(
    quote(bartFn),
    x.train = quote(scn[["x"]]),
    y.train = quote(scn[["y"]]),
    x.test = quote(scn[["x.test"]]),
    sigest = NA_real_,
    sigdf = 3.0,
    sigquant = 0.90,
    k = if (!is.null(scn[["k"]])) scn[["k"]] else 2.0,
    power = 2.0,
    base = 0.95,
    binaryOffset = if (!is.null(scn[["binaryOffset"]])) {
      scn[["binaryOffset"]]
    } else {
      0.0
    },
    ntree = if (!is.null(scn[["ntree"]])) scn[["ntree"]] else ntree,
    ndpost = ndpost,
    nskip = nskip,
    keepevery = 1L,
    keeptrainfits = TRUE,
    usequants = isTRUE(scn[["usequants"]]),
    numcut = if (!is.null(scn[["numcut"]])) scn[["numcut"]] else 100L,
    verbose = FALSE,
    nchain = if (!is.null(scn[["nchain"]])) scn[["nchain"]] else 1L,
    nthread = 1L,
    combinechains = TRUE,
    keeptrees = FALSE,
    keepcall = FALSE
  )
  if (!is.null(scn[["weights"]])) {
    bartArgs$weights <- quote(scn[["weights"]])
  }
  if (!is.null(proposalProbs)) {
    bartArgs$proposalprobs <- quote(proposalProbs)
  }
  fit <- quietly(eval(as.call(bartArgs)))
  list(
    train = fit$yhat.train,
    test = fit$yhat.test,
    sigma = fit$sigma,
    varcount = fit$varcount
  )
}

fitViaSampler <- function(scn) {
  sampler <- makeSampler(scn)
  r <- sampler$run(nskip, ndpost)
  list(
    train = poolChains(r$train),
    test = poolChains(r$test),
    sigma = as.vector(r$sigma),
    varcount = poolChains(r$varcount),
    k = if (!is.null(r$k)) as.vector(r$k)
  )
}

fitViaGibbsLoop <- function(scn) {
  sampler <- makeSampler(scn)
  sampler$run(nskip, 0L)
  perSweep <- ndpost %/% length(scn[["sweeps"]])
  parts <- lapply(scn[["sweeps"]], function(sweep) {
    sampler$setOffset(sweep[["offset"]])
    sampler$setResponse(sweep[["y"]])
    r <- sampler$run(0L, perSweep)
    list(
      train = poolChains(r$train),
      test = poolChains(r$test),
      sigma = as.vector(r$sigma),
      varcount = poolChains(r$varcount)
    )
  })
  list(
    train = do.call(rbind, lapply(parts, `[[`, "train")),
    test = do.call(rbind, lapply(parts, `[[`, "test")),
    sigma = unlist(lapply(parts, `[[`, "sigma")),
    varcount = do.call(rbind, lapply(parts, `[[`, "varcount"))
  )
}

fitViaSetPredictor <- function(scn) {
  sampler <- makeSampler(scn)
  half <- ndpost %/% 2L
  a <- sampler$run(nskip, half)
  sampler$setPredictor(scn[["x.swap"]], forceUpdate = TRUE)
  b <- sampler$run(0L, half)
  list(
    train = rbind(poolChains(a$train), poolChains(b$train)),
    test = rbind(poolChains(a$test), poolChains(b$test)),
    sigma = c(as.vector(a$sigma), as.vector(b$sigma)),
    varcount = rbind(poolChains(a$varcount), poolChains(b$varcount))
  )
}

# xbart reports a loss array rather than a fit, so it carries none of the
# channels below and returns its own summary vector whole.
fitViaXbart <- function(scn) {
  controlArgs <- list(
    n.cuts = 100L,
    n.thin = 1L,
    updateState = FALSE
  )
  if (controlTakesProposalProbs) {
    controlArgs$proposal.probs <- proposalProbs
  }
  frame <- data.frame(y = scn[["y"]], scn[["x"]])
  nTreesGrid <- c(20L, 100L)
  kGrid <- c(1, 2, 4)
  loss <- quietly(xbart(
    y ~ .,
    frame,
    n.samples = if (quick) 50L else 150L,
    method = "k-fold",
    n.test = 5L,
    n.reps = scn[["n.reps"]],
    n.burn = if (quick) c(50L, 25L) else c(150L, 75L),
    loss = "rmse",
    n.threads = 1L,
    n.trees = nTreesGrid,
    k = kGrid,
    power = 2.0,
    base = 0.95,
    family = gaussian(sigma = chisq(3.0, 0.90)),
    drop = TRUE,
    verbose = FALSE,
    control = do.call(dbartsControl, controlArgs)
  ))
  # n.reps x n.trees x k
  cellMean <- apply(loss, c(2L, 3L), mean)
  cellSd <- apply(loss, c(2L, 3L), sd)
  cellNames <- outer(
    paste0("t", nTreesGrid),
    paste0("k", kGrid),
    function(a, b) paste0(a, ".", b)
  )
  c(
    setNames(as.vector(cellMean), paste0("loss.", as.vector(cellNames))),
    setNames(as.vector(cellSd), paste0("losssd.", as.vector(cellNames))),
    loss.mean = mean(loss),
    loss.min.mean = mean(apply(loss, 1L, min))
  )
}

cvGrid <- list(n.trees = c(20L, 100L), k = c(1, 2, 4))

fitViaHandCV <- function(scn) {
  n <- length(scn[["y"]])
  nFolds <- 5L
  fold <- sample(rep_len(seq_len(nFolds), n))
  cells <- expand.grid(
    n.trees = cvGrid$n.trees,
    k = cvGrid$k,
    KEEP.OUT.ATTRS = FALSE
  )
  rmse <- vapply(
    seq_len(nrow(cells)),
    function(i) {
      squares <- 0
      for (f in seq_len(nFolds)) {
        train <- fold != f
        fit <- quietly(bartFn(
          x.train = scn[["x"]][train, , drop = FALSE],
          y.train = scn[["y"]][train],
          x.test = scn[["x"]][!train, , drop = FALSE],
          sigest = NA_real_,
          sigdf = 3.0,
          sigquant = 0.90,
          k = cells$k[i],
          power = 2.0,
          base = 0.95,
          ntree = cells$n.trees[i],
          ndpost = if (quick) 50L else 150L,
          nskip = if (quick) 50L else 150L,
          keepevery = 1L,
          keeptrainfits = FALSE,
          usequants = FALSE,
          numcut = 100L,
          verbose = FALSE,
          nchain = 1L,
          nthread = 1L,
          combinechains = TRUE,
          keeptrees = FALSE,
          keepcall = FALSE,
          proposalprobs = proposalProbs
        ))
        squares <- squares +
          sum((scn[["y"]][!train] - colMeans(fit$yhat.test))^2)
      }
      sqrt(squares / n)
    },
    0.0
  )
  names(rmse) <- paste0(
    "cv.t",
    cells$n.trees,
    ".k",
    cells$k
  )
  c(rmse, cv.mean = mean(rmse), cv.min = min(rmse))
}

summarizeFit <- function(scn, fit) {
  n <- nrow(scn[["x"]])
  idx <- unique(round(seq(1, n, length.out = n.train.summ)))
  vc <- fit$varcount
  result <- c(
    setNames(
      colMeans(fit$train[, idx, drop = FALSE]),
      paste0("fhat.train.", seq_along(idx))
    ),
    setNames(
      colMeans(fit$test),
      paste0("fhat.test.", seq_len(ncol(fit$test)))
    ),
    setNames(
      colMeans(vc),
      paste0("varcount.", seq_len(ncol(vc)))
    ),
    above.train = mean(rowMeans(fit$train > scn[["cutoff"]])),
    above.test = mean(rowMeans(fit$test > scn[["cutoff"]]))
  )
  if (!isTRUE(scn[["binary"]])) {
    result <- c(
      result,
      sigma.mean = mean(fit$sigma),
      sigma.sd = sd(fit$sigma)
    )
  }
  if (!is.null(fit$k)) {
    result <- c(result, k.mean = mean(fit$k), k.sd = sd(fit$k))
  }
  result
}

fitSummaries <- function(scn, seed) {
  set.seed(seed)
  if (scn[["kind"]] == "xbart") {
    return(fitViaXbart(scn))
  }
  if (scn[["kind"]] == "cv") {
    return(fitViaHandCV(scn))
  }
  fit <- switch(
    scn[["kind"]],
    bart = fitViaBart(scn),
    sampler = fitViaSampler(scn),
    gibbs = fitViaGibbsLoop(scn),
    setpredictor = fitViaSetPredictor(scn),
    stop("unknown scenario kind '", scn[["kind"]], "'")
  )
  summarizeFit(scn, fit)
}

numCores <- local({
  env <- Sys.getenv("CLASSIC_COMPARE_CORES", "")
  if (nzchar(env)) {
    as.integer(env)
  } else if (.Platform$OS.type == "windows") {
    1L
  } else {
    max(1L, parallel::detectCores() - 1L)
  }
})

runAll <- function(scenarios) {
  grid <- expand.grid(
    seed = seedOffset + seq_len(n.seeds),
    scenario = names(scenarios),
    stringsAsFactors = FALSE
  )
  fitOne <- function(i) {
    fitSummaries(scenarios[[grid$scenario[i]]], grid$seed[i])
  }
  rows <- if (numCores > 1L) {
    parallel::mclapply(
      seq_len(nrow(grid)),
      fitOne,
      mc.cores = numCores,
      mc.preschedule = FALSE
    )
  } else {
    lapply(seq_len(nrow(grid)), fitOne)
  }
  failed <- vapply(rows, inherits, TRUE, "try-error")
  if (any(failed)) {
    stop(
      "fit failed for ",
      paste(grid$scenario[failed], grid$seed[failed], collapse = ", "),
      ": ",
      conditionMessage(attr(rows[[which(failed)[1L]]], "condition"))
    )
  }
  lapply(setNames(names(scenarios), names(scenarios)), function(name) {
    do.call(rbind, rows[grid$scenario == name])
  })
}

scenarios <- makeScenarios()

scenarioFilter <- Sys.getenv("CLASSIC_COMPARE_SCENARIOS", "")
if (nzchar(scenarioFilter)) {
  wanted <- trimws(strsplit(scenarioFilter, ",", fixed = TRUE)[[1L]])
  unknown <- setdiff(wanted, names(scenarios))
  if (length(unknown) > 0L) {
    stop("unknown scenario(s): ", paste(unknown, collapse = ", "))
  }
  scenarios <- scenarios[wanted]
}

if (mode == "record") {
  out.file <- if (length(args) >= 2L) args[[2L]] else "classic-compare.rds"
  results <- runAll(scenarios)
  meta <- list(
    package.version = as.character(packageVersion("dbarts")),
    date = format(Sys.Date()),
    quick = quick,
    mixture = mixture,
    n.seeds = n.seeds,
    seed.offset = seedOffset,
    ndpost = ndpost,
    nskip = nskip,
    ntree = ntree
  )
  saveRDS(list(meta = meta, results = results), out.file)
  cat(
    "wrote",
    length(results),
    "scenarios x",
    n.seeds,
    "seeds (dbarts",
    meta$package.version,
    "mixture",
    mixture,
    ") to",
    out.file,
    "\n"
  )
} else if (mode == "merge") {
  # merge chunked record files (same meta) into one
  if (length(args) < 3L) {
    stop("usage: classic-compare.R merge out.rds in1.rds in2.rds ...")
  }
  parts <- lapply(args[-(1L:2L)], readRDS)
  meta <- parts[[1L]]$meta
  for (part in parts) {
    if (
      !identical(
        part$meta[setdiff(names(meta), "date")],
        meta[setdiff(names(meta), "date")]
      )
    ) {
      stop("parts were recorded with different settings")
    }
  }
  results <- do.call(c, lapply(parts, `[[`, "results"))
  saveRDS(list(meta = meta, results = results), args[[2L]])
  cat("merged", length(results), "scenarios to", args[[2L]], "\n")
} else if (mode == "compare") {
  if (length(args) < 3L) {
    stop("usage: classic-compare.R compare a.rds b.rds")
  }
  a.all <- readRDS(args[[2L]])
  b.all <- readRDS(args[[3L]])
  keys <- c("quick", "n.seeds", "seed.offset", "ndpost", "nskip", "ntree")
  if (!identical(a.all$meta[keys], b.all$meta[keys])) {
    stop("the two result sets were recorded with different settings")
  }
  cat(sprintf(
    "A: dbarts %s, mixture %s   B: dbarts %s, mixture %s\n",
    a.all$meta$package.version,
    a.all$meta$mixture,
    b.all$meta$package.version,
    b.all$meta$mixture
  ))
  cat(sprintf(
    "%d seeds, %d draws after %d burn-in\n\n",
    a.all$meta$n.seeds,
    a.all$meta$ndpost,
    a.all$meta$nskip
  ))
  cat(sprintf(
    "%-17s %5s  %7s  %-22s %5s %5s %5s\n",
    "scenario",
    "summ",
    "max |z|",
    "worst summary",
    ">2",
    ">3",
    ">4"
  ))
  anyFailure <- FALSE
  table <- NULL
  for (name in names(a.all$results)) {
    a <- a.all$results[[name]]
    b <- b.all$results[[name]]
    if (is.null(b)) {
      cat(sprintf("%-17s skipped (absent on side B)\n", name))
      next
    }
    if (identical(a, b)) {
      cat(sprintf("%-17s identical draws (same RNG stream)\n", name))
      next
    }
    z <- (colMeans(a) - colMeans(b)) /
      sqrt(apply(a, 2L, var) / nrow(a) + apply(b, 2L, var) / nrow(b))
    # The Welch z degenerates when one side diverges: a single blown-up seed
    # inflates the variance until |z| caps near sqrt(n.seeds). Disjoint seed
    # ranges are the exact complement - under exchangeability the chance is
    # 2 / choose(2n, n) per summary, negligible at 20 seeds a side.
    disjoint <- rep(FALSE, ncol(a))
    if (2 / choose(nrow(a) + nrow(b), nrow(a)) < 1e-8) {
      disjoint <- vapply(
        seq_len(ncol(a)),
        function(j) {
          isTRUE(max(a[, j]) < min(b[, j])) || isTRUE(max(b[, j]) < min(a[, j]))
        },
        TRUE
      )
    }
    finite <- is.finite(z)
    worst <- if (any(finite)) which.max(ifelse(finite, abs(z), -Inf)) else 1L
    row <- data.frame(
      scenario = name,
      summaries = length(z),
      max.abs.z = max(abs(z[finite]), 0),
      worst = names(z)[worst],
      n.gt2 = sum(abs(z) > 2, na.rm = TRUE),
      n.gt3 = sum(abs(z) > 3, na.rm = TRUE),
      n.gt4 = sum(abs(z) > 4, na.rm = TRUE),
      n.disjoint = sum(disjoint),
      stringsAsFactors = FALSE
    )
    table <- rbind(table, row)
    cat(sprintf(
      "%-17s %5d  %7.2f  %-22s %5d %5d %5d%s\n",
      name,
      row$summaries,
      row$max.abs.z,
      row$worst,
      row$n.gt2,
      row$n.gt3,
      row$n.gt4,
      if (row$n.disjoint > 0L) {
        sprintf("  %d disjoint <- FAIL", row$n.disjoint)
      } else {
        ""
      }
    ))
    if (row$n.gt4 > 0L || row$n.disjoint > 0L) {
      anyFailure <- TRUE
      bad <- which(abs(z) > 4)
      if (length(bad) > 0L) {
        cat(
          "    |z| > 4:",
          paste0(
            names(z)[bad],
            " (",
            round(z[bad], 2L),
            ")",
            collapse = ", "
          ),
          "\n"
        )
      }
      if (row$n.disjoint > 0L) {
        cat(
          "    disjoint:",
          paste0(
            names(z)[disjoint],
            " (",
            signif(colMeans(a)[disjoint], 4L),
            " vs ",
            signif(colMeans(b)[disjoint], 4L),
            ")",
            collapse = ", "
          ),
          "\n"
        )
      }
    }
  }
  if (!is.null(table)) {
    cat(sprintf(
      "\nover %d scenarios / %d summaries: max |z| = %.2f (%s, %s), %d with |z| > 3, %d with |z| > 4, %d disjoint\n",
      nrow(table),
      sum(table$summaries),
      max(table$max.abs.z),
      table$scenario[which.max(table$max.abs.z)],
      table$worst[which.max(table$max.abs.z)],
      sum(table$n.gt3),
      sum(table$n.gt4),
      sum(table$n.disjoint)
    ))
    out <- Sys.getenv("CLASSIC_COMPARE_TABLE", "")
    if (nzchar(out)) {
      utils::write.csv(table, out, row.names = FALSE)
    }
  }
  if (anyFailure) {
    quit(status = 1L)
  }
  cat("\nOK: posteriors statistically indistinguishable at |z| > 4\n")
} else {
  stop("unknown mode '", mode, "'; use record, merge or compare")
}
