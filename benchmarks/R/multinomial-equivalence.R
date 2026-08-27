#!/usr/bin/env Rscript

# Multinomial (softmax) bitwise-equivalence fixture. The single-forest anchor
# (equivalence.R) fits only through dbarts() - never a multinomial sampler - and
# byte-guards NEITHER the train nor the test SHAPE: it reads only
# test/varcount/sigma/k as summaries, never reads train, and its poolChains
# tolerates reshapes. So nothing in the existing anchors guards the multinomial
# engine's K-CHANNEL output seams (the n x K probability array, the per-category
# forest fits, the per-category variable counts). This fixture is that guard:
# record a baseline from build A, compare from a tree with build B installed, and
# every recorded channel must match to the bit (identical(), not tolerance). The
# heteroscedastic and hurdle models will re-touch exactly these multi-location
# output seams, so this fixture inherit-guards them the way bcf-equivalence.R
# guards the BCF report branches - the same lesson, verbatim.
#
# Each scenario drives an internal multinomial creation surface (R/bartcore.R;
# docs/design/multinomial.md) - the single-trial label path
# (bartcoreMultinomialSampler) for k3/k2, the grouped-count path
# (bartcoreMultinomialCountSampler) for k3counts and its mid-chain response
# swap for k3countsswap - at a fixed seed, single chain, one thread, and
# records:
#   - result$train, the K softmax-probability channels (n x K x n.samples),
#   - result$test, the same K softmax channels on a held-out x.test slice
#     (n.test x K x n.samples) - the C1 test-at-creation addition,
#   - the raw per-category forest fits (bartcoreForestFits 0 .. K-1),
#   - the per-category CUMULATIVE variable counts (bartcoreForestVariableCounts
#     0 .. K-1) - a final-state query, distinct from runVarcount below,
#   - result$varcount (runVarcount), the per-sample per-category run channel
#     (p x K x n.samples) - the widened storeSample varcount write.
#
# NEUTRALITY: supplying x.test consumes no rng, and widening the varcount write
# reads existing tree state (no draw), so the pre-existing channels (train, test,
# forestFits, varcount) reproduce the prior baselines bitwise. The compare loop
# checks only the channels the BASELINE recorded, so a compare against an older
# baseline lacking runVarcount checks exactly those and must pass; the
# re-recorded baseline additionally guards runVarcount.
#
# A getState -> setState continuation is deliberately NOT recorded: restore is
# structural by contract (omega is redrawn on the first restored sweep), so a
# bitwise assertion there would be a false gate.
#
# A scenario whose bitwise compare fails falls through to equivalence.R's
# statistical mode (docs/design/core-generalization.md): per-summary Welch
# z-scores over each statChannels entry's draws, against the baseline's
# recorded summaries, |z| >= 4 fails. A channel with no draws axis (a
# category's raw fit/varcount snapshot, a transaction's accept/reject
# verdict) has no statistical fallback and gates its mismatch directly, same
# as today. A baseline recorded before this mode existed carries no
# summaries and degrades loudly rather than silently passing.
#
# Usage:
#   Rscript multinomial-equivalence.R record [out.rds]
#   Rscript multinomial-equivalence.R compare baseline.rds
# Append 'quick' for a fast smoke pass (fewer draws; the settings guard refuses a
# mixed comparison against a full baseline).
# Append '--cross-host' to compare a baseline recorded on ANOTHER machine.
# Bitwise is unavailable there by construction: the scenario data run through
# exp()/plogis(), which are the platform libm, so the inputs already differ in
# the last ulp before any engine code runs. The flag exempts the point-in-time
# snapshot channels (one query of the live sampler after the last sweep, no
# draws axis to reduce over) and gates the draws-axis channels under a
# two-tier verdict - tier 1 a tight deviation bound that is the real gate,
# tier 2 a decoupled statistical fallback that adjudicates and never
# certifies. Compare only; a recording is host-local by definition.

source(
  system.file("common", "bartcoreHandle.R", package = "dbarts"),
  local = TRUE
)

suppressPackageStartupMessages(library(dbarts))

# drawMatrix, drawSummary, lag1Acf, toleranceRatio, essCompare,
# nonFiniteParts, compareCrossHost - shared with bcf-equivalence.R.
source(file.path(
  dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))),
  "equivalence-common.R"
))

args <- commandArgs(trailingOnly = TRUE)
quick <- "quick" %in% args
args <- setdiff(args, "quick")
crossHost <- "--cross-host" %in% args
args <- setdiff(args, "--cross-host")
mode <- if (length(args) >= 1L) args[[1L]] else "record"
if (crossHost && mode != "compare") {
  stop("--cross-host applies to compare only; a recording is host-local")
}

n.threads <- 1L
n.burn <- if (quick) 20L else 40L
n.samples <- if (quick) 20L else 40L
n.trees <- 40L

seeds <- c(
  k3.data = 6001L,
  k3.engine = 7001L,
  k2.data = 6002L,
  k2.engine = 7002L
)

makeControl <- function() {
  dbartsControl(
    n.chains = 1L,
    n.threads = n.threads,
    n.trees = n.trees,
    updateState = FALSE
  )
}

# The K softmax-probability channels (train and, when the fit carries x.test,
# test) plus every category forest's raw fits and split counts, at the sampler's
# current state, and the per-sample per-category varcount channel the run
# records (runVarcount, p x K x n.samples: the widened storeSample write). The
# test channel and runVarcount are additive over bb8855e/bcefa63; the earlier
# channels stay bitwise because widening the varcount write reads existing tree
# state (no draw) and supplying x.test consumes no rng. varcount below is the
# CUMULATIVE final-state per-category query (a different quantity), untouched.
#
# Channels a statistical compare can Welch-z: each carries a trailing draws
# axis (result$train/test/varcount, per bartcore_run's channel shape).
# forestFits and varcount are POINT-IN-TIME state queries with no repeated-
# draw axis to summarize, as are any accept/reject or install verdict.
statChannels <- c("train", "test", "runVarcount")

# The complement: every channel with no draws axis, and the only channels
# --cross-host exempts. These two vectors must partition the recorded names
# (minus summaries), so a channel added to either list has to be absent from
# the other; the compare loop re-checks the partition against what the
# baseline actually holds, which is where a channel appended at a call site
# rather than in recordChannels gets caught.
snapshotChannels <- c("forestFits", "varcount", "accepted", "installed")
stopifnot(length(intersect(snapshotChannels, statChannels)) == 0L)

# Tier 1's relative tolerance. Under a locked RNG stream two builds agree to
# rounding, so this is a stream-lock detector rather than a posterior test:
# loose enough that cross-host libm differences pass, ~1e6 times tighter than
# the tier-2 bar and so non-vacuous against every defect class the fixture
# exists to catch. The worst deviation an arm64-recorded baseline shows
# against an x86_64 Linux build sits about five orders of magnitude inside it.
crossHostRtol <- 1e-8

recordChannels <- function(bc, result, K) {
  ch <- list(
    train = result$train,
    test = result$test,
    forestFits = lapply(
      seq_len(K) - 1L,
      function(k) bartcoreForestFits(bc, k)
    ),
    varcount = lapply(
      seq_len(K) - 1L,
      function(k) bartcoreForestVariableCounts(bc, k)
    ),
    runVarcount = result$varcount
  )
  ch$summaries <- lapply(ch[statChannels], drawSummary)
  ch
}

runScenarios <- function() {
  n <- 200L
  p <- 4L
  result <- list()

  # (a) K = 3, covariate-dependent: a genuine softmax signal across three
  # categories, several predictors splitting.
  {
    set.seed(seeds[["k3.data"]])
    K <- 3L
    x <- matrix(runif(n * p), n, p)
    eta <- cbind(
      2 * (x[, 1L] - 0.5),
      x[, 2L] - x[, 3L],
      1.5 * (x[, 4L] - 0.5)
    )
    probs <- exp(eta) / rowSums(exp(eta))
    labels <- vapply(
      seq_len(n),
      function(i) sample.int(K, 1L, prob = probs[i, ]) - 1L,
      integer(1L)
    )
    # a held-out slice of x drives the additive test channel; slicing consumes
    # no rng, so labels/x above are byte-identical to the pre-test-channel run
    x.test <- x[seq_len(25L), , drop = FALSE]
    sampler <- dbarts(
      x,
      as.double(labels),
      test = x.test,
      control = makeControl()
    )
    set.seed(seeds[["k3.engine"]])
    bc <- dbarts:::bartcoreMultinomialSampler(sampler, labels, K = K)
    res <- bartcoreRun(bc, n.burn, n.samples)
    result$k3 <- recordChannels(bc, res, K)
  }

  # (b) K = 2, the logistic-equivalent two-category case.
  {
    set.seed(seeds[["k2.data"]])
    K <- 2L
    x <- matrix(runif(n * p), n, p)
    labels <- rbinom(n, 1L, plogis(2 * (x[, 1L] - 0.5) + x[, 2L]))
    x.test <- x[seq_len(25L), , drop = FALSE]
    sampler <- dbarts(
      x,
      as.double(labels),
      test = x.test,
      control = makeControl()
    )
    set.seed(seeds[["k2.engine"]])
    bc <- dbarts:::bartcoreMultinomialSampler(sampler, labels, K = K)
    res <- bartcoreRun(bc, n.burn, n.samples)
    result$k2 <- recordChannels(bc, res, K)
  }

  # (c) K = 3, GROUPED COUNTS (n_i > 1): a count matrix rather than single-trial
  # labels, driving the count-native combiner (the PG(n_i) summing draw and the
  # (y - n_i/2) working response) through bartcoreMultinomialCountSampler. Its
  # seeds are LITERALS kept out of the guarded `seeds` vector above so
  # settingsList() stays identical to the single-trial 5afb09a baseline and the
  # neutrality compare against it still runs (that compare checks only the k3/k2
  # scenarios it recorded). This scenario runs after k3/k2 with its own set.seed,
  # so it perturbs neither.
  {
    set.seed(6003L)
    K <- 3L
    x <- matrix(runif(n * p), n, p)
    eta <- cbind(
      2 * (x[, 1L] - 0.5),
      x[, 2L] - x[, 3L],
      1.5 * (x[, 4L] - 0.5)
    )
    probs <- exp(eta) / rowSums(exp(eta))
    trials <- sample(2:6, n, replace = TRUE)
    counts <- t(vapply(
      seq_len(n),
      function(i) rmultinom(1L, trials[i], probs[i, ])[, 1L],
      integer(K)
    ))
    x.test <- x[seq_len(25L), , drop = FALSE]
    sampler <- dbarts(
      x,
      as.double(counts[, 1L]),
      test = x.test,
      control = makeControl()
    )
    set.seed(7003L)
    bc <- dbarts:::bartcoreMultinomialCountSampler(sampler, counts, K = K)
    res <- bartcoreRun(bc, n.burn, n.samples)
    result$k3counts <- recordChannels(bc, res, K)
  }

  # (d) K = 3, FORCED whole-matrix setPredictor: run, replace the entire design
  # through the forced path - the only predictor mutation a multinomial sampler
  # accepts today, and the one that re-routes every category forest and
  # collapses any leaf the new codes empty - then run again and record the
  # post-swap state per category (docs/plans/multiforest-predictor-mutation.md,
  # "The harness this arc needs"). Nothing guarded that loop before. Its seeds
  # are LITERALS kept out of the guarded `seeds` vector above, as k3counts's
  # are, so settingsList() stays identical to the ec2a3d0 baseline and the
  # neutrality compare against it still runs; it runs last with its own
  # set.seed, so it perturbs none of the scenarios above. local() rather than
  # the bare block above: it scopes the same reused names without an own-line
  # brace.
  result$k3swap <- local({
    set.seed(6004L)
    K <- 3L
    x <- matrix(runif(n * p), n, p)
    eta <- cbind(
      2 * (x[, 1L] - 0.5),
      x[, 2L] - x[, 3L],
      1.5 * (x[, 4L] - 0.5)
    )
    probs <- exp(eta) / rowSums(exp(eta))
    labels <- vapply(
      seq_len(n),
      function(i) sample.int(K, 1L, prob = probs[i, ]) - 1L,
      integer(1L)
    )
    x2 <- matrix(runif(n * p), n, p)
    x.test <- x[seq_len(25L), , drop = FALSE]
    sampler <- dbarts(
      x,
      as.double(labels),
      test = x.test,
      control = makeControl()
    )
    set.seed(7004L)
    bc <- dbarts:::bartcoreMultinomialSampler(sampler, labels, K = K)
    bartcoreRun(bc, n.burn, n.samples)
    bartcoreSetPredictor(bc, x2, forceUpdate = TRUE)
    res <- bartcoreRun(bc, n.burn, n.samples)
    recordChannels(bc, res, K)
  })

  # (e)-(g) the TRANSACTIONAL predictor paths, which a K-forest multinomial
  # sampler began accepting when the two-phase revalidation widened across
  # forests_ (docs/plans/multiforest-predictor-mutation.md, S1). New streams:
  # these paths were refused at the bridge before this tip, so they have no
  # earlier baseline and become the regression floor from here. K = 3 makes
  # this the largest j-splitting tree count in scope, and the engine's verdict
  # rides beside the draws so a build that flipped an accept into a rollback
  # fails on the flag rather than only on the post-mutation state. Seeds are
  # LITERALS kept out of the guarded `seeds` vector, as k3counts's and
  # k3swap's are, so settingsList() stays identical to the 33f6fdc baseline;
  # each runs after the scenarios above with its own set.seed.
  makeK3 <- function(seed) {
    set.seed(seed)
    K <- 3L
    x <- matrix(runif(n * p), n, p)
    eta <- cbind(
      2 * (x[, 1L] - 0.5),
      x[, 2L] - x[, 3L],
      1.5 * (x[, 4L] - 0.5)
    )
    probs <- exp(eta) / rowSums(exp(eta))
    labels <- vapply(
      seq_len(n),
      function(i) sample.int(K, 1L, prob = probs[i, ]) - 1L,
      integer(1L)
    )
    list(x = x, labels = labels, K = K)
  }

  # (e) transactional whole-matrix setPredictor, sized to be ACCEPTED: every
  # tree of all three category forests re-routes and rebuilds. The jitter is
  # finer than the BCF fixture's, because the veto quantifies over three
  # forests of 40 trees rather than two of 50 and 25 (0.002 accepts at this
  # seed, 0.005 rolls back) - the count law the arc's measurement reports,
  # visible in the fixture itself.
  result$k3txn <- local({
    d <- makeK3(6005L)
    set.seed(6105L)
    x2 <- pmin(pmax(d$x + matrix(rnorm(n * p, 0, 0.002), n, p), 0), 1)
    sampler <- dbarts(
      d$x,
      as.double(d$labels),
      test = d$x[seq_len(25L), , drop = FALSE],
      control = makeControl()
    )
    set.seed(7005L)
    bc <- dbarts:::bartcoreMultinomialSampler(sampler, d$labels, K = d$K)
    bartcoreRun(bc, n.burn, n.samples)
    accepted <- bartcoreSetPredictor(bc, x2)
    res <- bartcoreRun(bc, n.burn, n.samples)
    c(recordChannels(bc, res, d$K), list(accepted = accepted))
  })

  # (f) transactional single-column updatePredictor: the subset the pruning
  # argument is stated over, summed over three forests rather than two.
  result$k3txncol <- local({
    d <- makeK3(6006L)
    set.seed(6106L)
    v <- pmin(pmax(d$x[, 2L] + rnorm(n, 0, 0.02), 0), 1)
    sampler <- dbarts(
      d$x,
      as.double(d$labels),
      test = d$x[seq_len(25L), , drop = FALSE],
      control = makeControl()
    )
    set.seed(7006L)
    bc <- dbarts:::bartcoreMultinomialSampler(sampler, d$labels, K = d$K)
    bartcoreRun(bc, n.burn, n.samples)
    accepted <- bartcoreUpdatePredictor(bc, v, 2L)
    res <- bartcoreRun(bc, n.burn, n.samples)
    c(recordChannels(bc, res, d$K), list(accepted = accepted))
  })

  # (g) a transactional proposal built to be ROLLED BACK, as the BCF fixture's
  # is: the veto reaches across every category forest, and the rollback has to
  # put every one of them back.
  result$k3reject <- local({
    d <- makeK3(6007L)
    sampler <- dbarts(
      d$x,
      as.double(d$labels),
      test = d$x[seq_len(25L), , drop = FALSE],
      control = makeControl()
    )
    set.seed(7007L)
    bc <- dbarts:::bartcoreMultinomialSampler(sampler, d$labels, K = d$K)
    bartcoreRun(bc, n.burn, n.samples)
    accepted <- bartcoreUpdatePredictor(
      bc,
      ifelse(seq_len(n) %% 2L == 0L, 0.25, 0.75),
      1L
    )
    res <- bartcoreRun(bc, n.burn, n.samples)
    c(recordChannels(bc, res, d$K), list(accepted = accepted))
  })

  # (h)-(i) the PER-OBSERVATION session, which a K-forest multinomial sampler
  # began accepting when the session's cell guard widened across forests_
  # (docs/plans/multiforest-predictor-mutation.md, S2). New streams, mirroring
  # the BCF fixture's pair. The verdict channel is the install MASK - the
  # session answers per row, not per transaction - and K = 3 makes this the
  # widest guarded set in the arc's scope: the veto quantifies over every tree
  # of all three category forests that splits on the column. Seeds are
  # LITERALS kept out of the guarded `seeds` vector.
  #
  # (h) the ACCEPT shape: a jitter nearly every row can take.
  result$k3perobs <- local({
    d <- makeK3(6008L)
    set.seed(6108L)
    v <- pmin(pmax(d$x[, 2L] + rnorm(n, 0, 0.02), 0), 1)
    sampler <- dbarts(
      d$x,
      as.double(d$labels),
      test = d$x[seq_len(25L), , drop = FALSE],
      control = makeControl()
    )
    set.seed(7008L)
    bc <- dbarts:::bartcoreMultinomialSampler(sampler, d$labels, K = d$K)
    bartcoreRun(bc, n.burn, n.samples)
    installed <- bartcoreUpdatePredictorPerObservation(bc, v, 2L)
    res <- bartcoreRun(bc, n.burn, n.samples)
    c(recordChannels(bc, res, d$K), list(installed = installed))
  })

  # (i) the DECLINE shape: the two-level replacement that rolls the whole
  # transaction back in (g) instead declines row by row here, leaving the
  # sampler on a partly installed column.
  result$k3perobspartial <- local({
    d <- makeK3(6009L)
    sampler <- dbarts(
      d$x,
      as.double(d$labels),
      test = d$x[seq_len(25L), , drop = FALSE],
      control = makeControl()
    )
    set.seed(7009L)
    bc <- dbarts:::bartcoreMultinomialSampler(sampler, d$labels, K = d$K)
    bartcoreRun(bc, n.burn, n.samples)
    installed <- bartcoreUpdatePredictorPerObservation(
      bc,
      ifelse(seq_len(n) %% 2L == 0L, 0.25, 0.75),
      1L
    )
    res <- bartcoreRun(bc, n.burn, n.samples)
    c(recordChannels(bc, res, d$K), list(installed = installed))
  })

  # (j) the CATEGORY OFFSET channel: an n x K TRAIN offset and an n.test x K
  # TEST offset, the response-side shifts that enter the latent on both sides
  # of the softmax (docs/plans/multinomial-counts-mutation.md). Both entrances
  # are exercised in the one scenario - both offsets taken at CREATION, then
  # both REPLACED mid-chain through their setters - so the recorded state
  # moves if either the creation lift or the mutation-time rematerialization
  # of the raw fits does. The train offset enters every category's working
  # response and rides the draws from the first sweep; the test offset enters
  # no likelihood and moves only the recorded test channel, which is why one
  # scenario can guard both. New stream: the channel did not exist before this
  # tip, so it has no earlier baseline and becomes the regression floor from
  # here. Its seeds are LITERALS kept out of the guarded `seeds` vector, as
  # (c)-(i)'s are, so settingsList() stays identical to the a825263 baseline
  # and the neutrality compare against it runs at all; it runs last with its
  # own set.seed, so it perturbs none of the scenarios above.
  result$k3offset <- local({
    d <- makeK3(6010L)
    n.test <- 25L
    x.test <- d$x[seq_len(n.test), , drop = FALSE]
    set.seed(6110L)
    offset <- matrix(rnorm(n * d$K, 0, 0.5), n, d$K)
    offset.test <- matrix(rnorm(n.test * d$K, 0, 0.5), n.test, d$K)
    offset2 <- matrix(rnorm(n * d$K, 0, 0.5), n, d$K)
    offset.test2 <- matrix(rnorm(n.test * d$K, 0, 0.5), n.test, d$K)
    sampler <- dbarts(
      d$x,
      as.double(d$labels),
      test = x.test,
      control = makeControl()
    )
    set.seed(7010L)
    bc <- dbarts:::bartcoreMultinomialSampler(
      sampler,
      d$labels,
      K = d$K,
      offset = offset,
      offset.test = offset.test
    )
    bartcoreRun(bc, n.burn, n.samples)
    bartcoreSetCategoryOffset(bc, offset2)
    bartcoreSetCategoryTestOffset(bc, offset.test2)
    res <- bartcoreRun(bc, n.burn, n.samples)
    recordChannels(bc, res, d$K)
  })

  # (k) K = 3 GROUPED COUNTS with a mid-chain RESPONSE swap: create on one n x K
  # count matrix, run, replace the whole response through bartcoreSetCounts -
  # the trees carry over, fitted to the previous counts, and the next run forms
  # every category's working response against the new ones - then run again and
  # record the post-swap state. The swap stream was only TRANSITIVELY pinned
  # before this scenario (test-multinomial-counts-mutation.R pins swap ==
  # rebuild bitwise, and k3counts pins the rebuild against the frozen
  # baseline), and that identity dissolves silently under a slice which
  # deliberately revises swap semantics and regenerates it as expected; the
  # recorded draws pin it directly. Its seeds are LITERALS kept out of the
  # guarded `seeds` vector, as (c)-(j)'s are, so settingsList() stays identical
  # to the 1027be5 baseline and the compare against it runs at all; it runs last
  # with its own set.seed, so it perturbs none of the scenarios above.
  result$k3countsswap <- local({
    set.seed(6011L)
    K <- 3L
    x <- matrix(runif(n * p), n, p)
    eta <- cbind(
      2 * (x[, 1L] - 0.5),
      x[, 2L] - x[, 3L],
      1.5 * (x[, 4L] - 0.5)
    )
    probs <- exp(eta) / rowSums(exp(eta))
    drawCounts <- function() {
      trials <- sample(2:6, n, replace = TRUE)
      t(vapply(
        seq_len(n),
        function(i) rmultinom(1L, trials[i], probs[i, ])[, 1L],
        integer(K)
      ))
    }
    counts <- drawCounts()
    counts2 <- drawCounts()
    sampler <- dbarts(
      x,
      as.double(counts[, 1L]),
      test = x[seq_len(25L), , drop = FALSE],
      control = makeControl()
    )
    set.seed(7011L)
    bc <- dbarts:::bartcoreMultinomialCountSampler(sampler, counts, K = K)
    bartcoreRun(bc, n.burn, n.samples)
    bartcoreSetCounts(bc, counts2)
    res <- bartcoreRun(bc, n.burn, n.samples)
    recordChannels(bc, res, K)
  })

  result
}

settingsList <- function() {
  list(
    quick = quick,
    n.threads = n.threads,
    n.burn = n.burn,
    n.samples = n.samples,
    n.trees = n.trees,
    seeds = seeds
  )
}

if (mode == "record") {
  out.file <- if (length(args) >= 2L) {
    args[[2L]]
  } else {
    "multinomial-equivalence-baseline.rds"
  }
  results <- runScenarios()
  meta <- c(
    list(
      rev = system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE),
      date = format(Sys.Date())
    ),
    settingsList()
  )
  saveRDS(list(meta = meta, results = results), out.file)
  cat(
    "wrote multinomial baseline for",
    length(results),
    "scenarios to",
    out.file,
    "\n"
  )
} else if (mode == "compare") {
  if (length(args) < 2L) {
    stop("usage: multinomial-equivalence.R compare baseline.rds")
  }
  baseline <- readRDS(args[[2L]])
  guarded <- names(settingsList())
  if (!identical(baseline$meta[guarded], settingsList())) {
    stop(
      "baseline was recorded with different settings: ",
      paste(
        guarded,
        vapply(
          baseline$meta[guarded],
          function(v) paste(v, collapse = ","),
          ""
        ),
        sep = "=",
        collapse = "; "
      )
    )
  }

  results <- runScenarios()
  anyFailure <- FALSE
  usedStatistical <- FALSE
  nCompared <- 0L
  nSkipped <- 0L
  for (name in names(baseline$results)) {
    a <- baseline$results[[name]]
    b <- results[[name]]
    # a gate that compares nothing must not pass: a renamed or dropped scenario
    # leaves every channel the baseline recorded for it unchecked
    if (is.null(b)) {
      anyFailure <- TRUE
      nSkipped <- nSkipped + 1L
      cat(sprintf("%-6s skipped (not produced this run) <- FAIL\n", name))
      next
    }
    nCompared <- nCompared + 1L
    # statChannels and snapshotChannels have to partition what the baseline
    # holds. A channel appended at a call site rather than in recordChannels is
    # invisible to a check placed there, so the partition is re-asserted here
    # against the recorded names: a channel added later stops the compare and
    # forces the taxonomy decision instead of silently defaulting to
    # gated-and-red off-host.
    unclassified <- setdiff(
      names(a),
      c("summaries", statChannels, snapshotChannels)
    )
    if (length(unclassified) > 0L) {
      stop(
        "unclassified channel(s) in ",
        name,
        ": ",
        paste(unclassified, collapse = ", ")
      )
    }
    if (crossHost) {
      verdict <- compareCrossHost(
        name,
        a,
        b,
        snapshotChannels,
        crossHostRtol,
        6L
      )
      anyFailure <- anyFailure || verdict$fail
      usedStatistical <- usedStatistical || verdict$decoupled
      next
    }
    channels <- setdiff(names(a), "summaries")
    ok <- vapply(
      channels,
      function(ch) identical(a[[ch]], b[[ch]]),
      logical(1L)
    )
    if (all(ok)) {
      cat(sprintf(
        "%-6s identical (all %d channels: %s)\n",
        name,
        length(channels),
        paste(channels, collapse = ", ")
      ))
      next
    }
    usedStatistical <- TRUE
    # summaries is a pure, deterministic reduction of the raw draws-axis
    # channels, so a baseline that predates the field but stores those channels
    # at full shape can still be compared: derive what it did not record.
    aSummaries <- a[["summaries"]]
    if (is.null(aSummaries) && all(statChannels %in% names(a))) {
      aSummaries <- lapply(a[statChannels], drawSummary)
    }
    # Nothing to Welch-z against at all: degrade loudly rather than silently
    # passing a possibly-real divergence.
    if (is.null(aSummaries)) {
      anyFailure <- TRUE
      cat(sprintf(
        "%-6s statistical compare unsupported by this baseline (recorded before summaries)\n",
        name
      ))
      next
    }
    bSummaries <- lapply(b[statChannels], drawSummary)
    z <- unlist(Map(
      function(sa, sb) {
        # A channel whose shape changed would otherwise recycle into a garbage
        # z rather than erroring.
        stopifnot(length(sa$mean) == length(sb$mean))
        (sa$mean - sb$mean) / sqrt(sa$var / sa$n + sb$var / sb$n)
      },
      aSummaries,
      bSummaries
    ))
    n.warn <- sum(abs(z) > 3, na.rm = TRUE)
    n.fail <- sum(abs(z) > 4, na.rm = TRUE)
    cat(sprintf(
      "%-6s %d summaries, max |z| = %.2f%s%s\n",
      name,
      length(z),
      max(abs(z), na.rm = TRUE),
      if (n.warn > 0L) sprintf(", %d with |z| > 3", n.warn) else "",
      if (n.fail > 0L) sprintf(", %d with |z| > 4 <- FAIL", n.fail) else ""
    ))
    if (n.fail > 0L) {
      anyFailure <- TRUE
      failed <- which(abs(z) > 4)
      cat(
        "  worst offenders:",
        paste0(
          names(z)[failed],
          " (z=",
          round(z[failed], 2L),
          ")",
          collapse = ", "
        ),
        "\n"
      )
    }
    # A NaN or Inf draw makes every z NaN and a max over NaNs report -Inf, so
    # the z-verdict above cannot see one at all. Count them directly.
    gated <- intersect(statChannels, names(a))
    nf <- vapply(
      gated,
      function(ch) c(sum(!is.finite(a[[ch]])), sum(!is.finite(b[[ch]]))),
      numeric(2L)
    )
    if (sum(nf) > 0) {
      anyFailure <- TRUE
      cat(sprintf(
        "%-6s NON-FINITE: %s <- FAIL\n",
        name,
        paste(nonFiniteParts(gated, nf), collapse = ", ")
      ))
    }
    # A mismatch outside statChannels (a point-in-time snapshot, or a
    # transaction's accept/reject verdict) has no draws to Welch-z, so it
    # gates independently of the z-verdict above - never let a clean z-score
    # silently pass a real divergence the z gate cannot see at all.
    pointMismatch <- setdiff(channels[!ok], statChannels)
    if (length(pointMismatch) > 0L) {
      anyFailure <- TRUE
      cat(sprintf(
        "%-6s also MISMATCH in %s (no statistical fallback) <- FAIL\n",
        name,
        paste(pointMismatch, collapse = ", ")
      ))
    }
  }

  cat(sprintf("%d compared / %d skipped\n", nCompared, nSkipped))
  if (anyFailure) {
    quit(status = 1L)
  }
  cat(
    if (crossHost && usedStatistical) {
      "\nOK: every gated multinomial channel passed cross-host tier 1, or tier 2 could not distinguish it - a weak bar, so report any decoupled scenario\n"
    } else if (crossHost) {
      "\nOK: every gated multinomial channel within the cross-host tier 1 bound (snapshot channels exempt)\n"
    } else if (usedStatistical) {
      "\nOK: every multinomial channel identical, or statistically indistinguishable (|z| < 4)\n"
    } else {
      "\nOK: every multinomial channel bitwise identical across every scenario\n"
    }
  )
} else {
  stop("unknown mode '", mode, "'; use record or compare")
}
