#!/usr/bin/env Rscript

# BCF bitwise-equivalence fixture. The single-forest anchor (equivalence.R)
# fits only through dbarts(), which never builds a two-forest BCF sampler, so
# it says nothing about the BCF path. This fixture is the BCF analogue: record
# a baseline from build A, compare from a tree with build B installed, and
# every recorded channel must match to the bit (identical(), not tolerance).
# A neutral refactor of the BCF path is proven neutral only if this reports
# identical on every scenario.
#
# Each scenario drives the internal bartcoreBCFSampler surface (R/bartcore.R;
# docs/design/bcf.md) at a fixed seed, single chain, one thread, and records:
#   - the raw per-forest fits of BOTH forests (bartcoreForestFits 0 and 1),
#   - the glue (bartcoreBCFGlue: a, b0, b1),
#   - sigma,
#   - the reported result$train and result$varcount channels,
#   - the treatment forest's own split counts
#     (bartcoreForestVariableCounts 1).
# train and varcount are the only bitwise guard on storeSample's BCF report
# branches; the per-forest fits and glue guard the combining math and the
# coupling draw. result$varcount is the reported channel, which for a BCF
# sampler is the PROGNOSTIC forest alone, so varcount.tau is the only guard on
# the treatment forest's tree structure - without it a change that moved tau's
# splits but not mu's would read as a leaf-value change.
#
# A getState -> setState continuation is deliberately NOT recorded: restore is
# structural by contract (the dropped accumulation history is not reproduced,
# test-bcf.R), so a bitwise assertion there would be a false gate.
#
# Usage:
#   Rscript bcf-equivalence.R record [out.rds]
#   Rscript bcf-equivalence.R compare baseline.rds
# Append 'quick' for a fast smoke pass (fewer draws; not comparable to a full
# baseline - the settings guard refuses the mixed comparison).

suppressPackageStartupMessages(library(dbarts))

args <- commandArgs(trailingOnly = TRUE)
quick <- "quick" %in% args
args <- setdiff(args, "quick")
mode <- if (length(args) >= 1L) args[[1L]] else "record"

# Global run knobs. These are pinned in the baseline meta and re-checked at
# compare; changing any invalidates a recorded baseline, so compare stops on a
# mismatch rather than silently comparing apples to oranges.
n.threads <- 1L
n.burn <- if (quick) 20L else 40L
n.samples <- if (quick) 20L else 40L
n.trees.mu <- 50L
n.trees.tau <- 25L

# Fixed per-scenario seeds. Pinned in the meta so a seed edit is caught as a
# settings change (the values themselves live here as ordinary code constants,
# as equivalence.R's scenario seeds do).
seeds <- c(
  default.data = 8001L,
  default.engine = 9001L,
  restricted.data = 8002L,
  restricted.engine = 9002L,
  glue.data = 8003L,
  glue.engine = 9003L,
  weighted.data = 8004L,
  weighted.weights = 8104L,
  weighted.engine = 9004L,
  treatment.data = 8005L,
  treatment.z2 = 8105L,
  treatment.engine = 9005L
)

# A modest two-forest problem: prognostic mu(x) plus a treatment effect
# z * tau(x), Gaussian noise. Small but nontrivial (both forests carry real
# structure, several predictors split).
makeData <- function(n, p, seed) {
  set.seed(seed)
  x <- matrix(runif(n * p), n, p)
  z <- rbinom(n, 1L, 0.5)
  mu <- 2 * sin(pi * x[, 1L]) + x[, 2L]
  tau <- 1 + 2 * x[, 3L]
  y <- mu + z * tau + rnorm(n, sd = 0.2)
  list(x = x, y = y, z = z)
}

makeControl <- function() {
  dbartsControl(
    n.chains = 1L,
    n.threads = n.threads,
    n.trees = n.trees.mu,
    updateState = FALSE
  )
}

# The full recorded channel set for one BCF sampler at its current state.
recordChannels <- function(bcSampler, result) {
  list(
    mu = dbarts:::bartcoreForestFits(bcSampler, 0L),
    tau = dbarts:::bartcoreForestFits(bcSampler, 1L),
    glue = dbarts:::bartcoreBCFGlue(bcSampler),
    sigma = result$sigma,
    train = result$train,
    varcount = result$varcount,
    varcount.tau = dbarts:::bartcoreForestVariableCounts(bcSampler, 1L)
  )
}

runScenarios <- function() {
  n <- 200L
  p <- 4L
  result <- list()

  # (a) default two-forest BCF
  {
    d <- makeData(n, p, seeds[["default.data"]])
    sampler <- dbarts(d$x, d$y, control = makeControl())
    set.seed(seeds[["default.engine"]])
    bc <- dbarts:::bartcoreBCFSampler(
      sampler,
      d$z,
      n.trees.treatment = n.trees.tau
    )
    res <- dbarts:::bartcoreRun(bc, n.burn, n.samples)
    result$default <- recordChannels(bc, res)
  }

  # (b) restricted-moderator BCF: the treatment forest reads only {x1, x3}
  {
    d <- makeData(n, p, seeds[["restricted.data"]])
    colnames(d$x) <- paste0("x", seq_len(p))
    sampler <- dbarts(d$x, d$y, control = makeControl())
    set.seed(seeds[["restricted.engine"]])
    bc <- dbarts:::bartcoreBCFSampler(
      sampler,
      d$z,
      n.trees.treatment = n.trees.tau,
      moderators = c("x1", "x3")
    )
    res <- dbarts:::bartcoreRun(bc, n.burn, n.samples)
    result$restricted <- recordChannels(bc, res)
  }

  # (c) asymmetric glue toggle: update.a on (a and its ridge move draw), b held
  # at (b0, b1) = (0, 1). Isolates the b-block-fixed routing the both-on
  # scenarios never exercise.
  {
    d <- makeData(n, p, seeds[["glue.data"]])
    sampler <- dbarts(d$x, d$y, control = makeControl())
    set.seed(seeds[["glue.engine"]])
    bc <- dbarts:::bartcoreBCFSampler(
      sampler,
      d$z,
      n.trees.treatment = n.trees.tau,
      update.a = TRUE,
      update.b = FALSE
    )
    res <- dbarts:::bartcoreRun(bc, n.burn, n.samples)
    result$glue_toggle <- recordChannels(bc, res)
  }

  # (d) weighted BCF: weights ~ U(0.5, 2) ride the data object into the
  # per-forest residual's w * m^2 weight channel and the glue draws, otherwise
  # unobserved by any gate.
  {
    d <- makeData(n, p, seeds[["weighted.data"]])
    set.seed(seeds[["weighted.weights"]])
    weights <- runif(n, 0.5, 2)
    sampler <- dbarts(d$x, d$y, weights = weights, control = makeControl())
    set.seed(seeds[["weighted.engine"]])
    bc <- dbarts:::bartcoreBCFSampler(
      sampler,
      d$z,
      n.trees.treatment = n.trees.tau
    )
    res <- dbarts:::bartcoreRun(bc, n.burn, n.samples)
    result$weighted <- recordChannels(bc, res)
  }

  # (e) setTreatment: run, swap the treatment vector, run again, record the
  # post-mutation state - the only bitwise guard on the setTreatment routing.
  {
    d <- makeData(n, p, seeds[["treatment.data"]])
    set.seed(seeds[["treatment.z2"]])
    z2 <- rbinom(n, 1L, 0.5)
    sampler <- dbarts(d$x, d$y, control = makeControl())
    set.seed(seeds[["treatment.engine"]])
    bc <- dbarts:::bartcoreBCFSampler(
      sampler,
      d$z,
      n.trees.treatment = n.trees.tau
    )
    dbarts:::bartcoreRun(bc, n.burn, n.samples)
    dbarts:::bartcoreSetTreatment(bc, z2)
    res <- dbarts:::bartcoreRun(bc, n.burn, n.samples)
    result$set_treatment <- recordChannels(bc, res)
  }

  # (f) forced whole-matrix setPredictor: run, replace the entire design
  # through the FORCED path - the only predictor mutation a BCF sampler accepts
  # today, and the one that re-routes both forests and collapses any leaf the
  # new codes empty - then run again and record the post-swap state per forest
  # (docs/plans/multiforest-predictor-mutation.md, "The harness this arc
  # needs"). Nothing guarded that loop before. Its seeds are LITERALS kept out
  # of the guarded `seeds` vector above so settingsList() stays identical to
  # the c820227 baseline and the neutrality compare against it still runs (that
  # compare checks only the five scenarios it recorded); the k3counts precedent
  # in multinomial-equivalence.R. This scenario runs after the five above with
  # its own set.seed, so it perturbs none of them. local() rather than the bare
  # block above: it scopes the same reused names without an own-line brace.
  result$set_predictor <- local({
    d <- makeData(n, p, 8006L)
    set.seed(8106L)
    x2 <- matrix(runif(n * p), n, p)
    sampler <- dbarts(d$x, d$y, control = makeControl())
    set.seed(9006L)
    bc <- dbarts:::bartcoreBCFSampler(
      sampler,
      d$z,
      n.trees.treatment = n.trees.tau
    )
    dbarts:::bartcoreRun(bc, n.burn, n.samples)
    dbarts:::bartcoreSetPredictor(bc, x2, forceUpdate = TRUE)
    res <- dbarts:::bartcoreRun(bc, n.burn, n.samples)
    recordChannels(bc, res)
  })

  # (g)-(i) the TRANSACTIONAL predictor paths, which a BCF sampler began
  # accepting when the two-phase revalidation widened across forests_
  # (docs/plans/multiforest-predictor-mutation.md, S1). New streams: these
  # paths were refused at the bridge before this tip, so they have no earlier
  # baseline and become the regression floor from here. Each records the
  # engine's verdict beside the draws, so a build that silently flipped an
  # accept into a rollback (or the reverse) fails on the flag rather than only
  # on the post-mutation state. Seeds are LITERALS kept out of the guarded
  # `seeds` vector, as (f)'s are, so settingsList() stays identical to the
  # 33f6fdc baseline and the neutrality compare against it still runs; each
  # runs after the scenarios above with its own set.seed and perturbs none.
  #
  # (g) transactional whole-matrix setPredictor: a jitter of the design sized
  # to be ACCEPTED, so every tree of BOTH forests re-routes against the new
  # codes and the fits rebuild from the recovered leaf parameters.
  result$set_predictor_txn <- local({
    d <- makeData(n, p, 8007L)
    set.seed(8107L)
    x2 <- pmin(pmax(d$x + matrix(rnorm(n * p, 0, 0.005), n, p), 0), 1)
    sampler <- dbarts(d$x, d$y, control = makeControl())
    set.seed(9007L)
    bc <- dbarts:::bartcoreBCFSampler(
      sampler,
      d$z,
      n.trees.treatment = n.trees.tau
    )
    dbarts:::bartcoreRun(bc, n.burn, n.samples)
    accepted <- dbarts:::bartcoreSetPredictor(bc, x2)
    res <- dbarts:::bartcoreRun(bc, n.burn, n.samples)
    c(recordChannels(bc, res), list(accepted = accepted))
  })

  # (h) transactional single-column updatePredictor: the subset the pruning
  # argument is stated over - only the trees splitting on column 3 can veto,
  # and only they are re-routed and rebuilt in the treatment forest. Column 3
  # is the moderator tau is a function of, so tau really does split on it.
  result$update_column_txn <- local({
    d <- makeData(n, p, 8008L)
    set.seed(8108L)
    v <- pmin(pmax(d$x[, 3L] + rnorm(n, 0, 0.02), 0), 1)
    sampler <- dbarts(d$x, d$y, control = makeControl())
    set.seed(9008L)
    bc <- dbarts:::bartcoreBCFSampler(
      sampler,
      d$z,
      n.trees.treatment = n.trees.tau
    )
    dbarts:::bartcoreRun(bc, n.burn, n.samples)
    accepted <- dbarts:::bartcoreUpdatePredictor(bc, v, 3L)
    res <- dbarts:::bartcoreRun(bc, n.burn, n.samples)
    c(recordChannels(bc, res), list(accepted = accepted))
  })

  # (i) a transactional proposal built to be ROLLED BACK: a two-level
  # replacement column collapses every observation onto two values of the
  # existing grid, which empties leaves in any tree splitting on it, so the
  # whole change reverts and the sampler continues from the pre-proposal
  # state. Gates that the rollback leaves a BCF run bitwise - the widened
  # revalidation touches both forests before it fails, and repartitionTrees
  # has to put both back.
  result$update_column_reject <- local({
    d <- makeData(n, p, 8009L)
    sampler <- dbarts(d$x, d$y, control = makeControl())
    set.seed(9009L)
    bc <- dbarts:::bartcoreBCFSampler(
      sampler,
      d$z,
      n.trees.treatment = n.trees.tau
    )
    dbarts:::bartcoreRun(bc, n.burn, n.samples)
    accepted <- dbarts:::bartcoreUpdatePredictor(
      bc,
      ifelse(seq_len(n) %% 2L == 0L, 0.25, 0.75),
      1L
    )
    res <- dbarts:::bartcoreRun(bc, n.burn, n.samples)
    c(recordChannels(bc, res), list(accepted = accepted))
  })

  # (j)-(k) the PER-OBSERVATION session, which a BCF sampler began accepting
  # when the session's cell guard widened across forests_
  # (docs/plans/multiforest-predictor-mutation.md, S2). New streams, on the
  # same terms as (g)-(i): the path was refused at the bridge before this tip.
  # The verdict channel here is the install MASK itself - the session's answer
  # is per row, not per transaction - so a build that moved one row's decision
  # fails on `installed` rather than only on the post-mutation draws. This is
  # also the only BCF scenario that consumes the engine's scan permutation.
  # Seeds are LITERALS kept out of the guarded `seeds` vector, as (f)-(i)'s
  # are, so settingsList() stays identical to the 938eb81 baseline.
  #
  # (j) the ACCEPT shape: a jitter of the moderator tau is a function of,
  # sized so nearly every row installs, then a second leg so the carried-over
  # tree state matters.
  result$per_observation <- local({
    d <- makeData(n, p, 8010L)
    set.seed(8110L)
    v <- pmin(pmax(d$x[, 3L] + rnorm(n, 0, 0.02), 0), 1)
    sampler <- dbarts(d$x, d$y, control = makeControl())
    set.seed(9010L)
    bc <- dbarts:::bartcoreBCFSampler(
      sampler,
      d$z,
      n.trees.treatment = n.trees.tau
    )
    dbarts:::bartcoreRun(bc, n.burn, n.samples)
    installed <- dbarts:::bartcoreUpdatePredictorPerObservation(bc, v, 3L)
    res <- dbarts:::bartcoreRun(bc, n.burn, n.samples)
    c(recordChannels(bc, res), list(installed = installed))
  })

  # (k) the DECLINE shape: the two-level replacement (i) rolls a whole
  # transaction back on instead declines row by row, so the recorded mask is
  # the per-row rollback and the sampler continues from a partly installed
  # column - the state the whole-transaction scenarios cannot reach.
  result$per_observation_partial <- local({
    d <- makeData(n, p, 8011L)
    sampler <- dbarts(d$x, d$y, control = makeControl())
    set.seed(9011L)
    bc <- dbarts:::bartcoreBCFSampler(
      sampler,
      d$z,
      n.trees.treatment = n.trees.tau
    )
    dbarts:::bartcoreRun(bc, n.burn, n.samples)
    installed <- dbarts:::bartcoreUpdatePredictorPerObservation(
      bc,
      ifelse(seq_len(n) %% 2L == 0L, 0.25, 0.75),
      1L
    )
    res <- dbarts:::bartcoreRun(bc, n.burn, n.samples)
    c(recordChannels(bc, res), list(installed = installed))
  })

  result
}

settingsList <- function() {
  list(
    quick = quick,
    n.threads = n.threads,
    n.burn = n.burn,
    n.samples = n.samples,
    n.trees.mu = n.trees.mu,
    n.trees.tau = n.trees.tau,
    seeds = seeds
  )
}

if (mode == "record") {
  out.file <- if (length(args) >= 2L) {
    args[[2L]]
  } else {
    "bcf-equivalence-baseline.rds"
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
  cat("wrote BCF baseline for", length(results), "scenarios to", out.file, "\n")
} else if (mode == "compare") {
  if (length(args) < 2L) {
    stop("usage: bcf-equivalence.R compare baseline.rds")
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
  for (name in names(baseline$results)) {
    a <- baseline$results[[name]]
    b <- results[[name]]
    if (is.null(b)) {
      cat(sprintf("%-14s skipped (not produced this run)\n", name))
      next
    }
    channels <- names(a)
    ok <- vapply(
      channels,
      function(ch) identical(a[[ch]], b[[ch]]),
      logical(1L)
    )
    if (all(ok)) {
      cat(sprintf(
        "%-14s identical (all %d channels: %s)\n",
        name,
        length(channels),
        paste(channels, collapse = ", ")
      ))
    } else {
      anyFailure <- TRUE
      cat(sprintf(
        "%-14s MISMATCH in: %s\n",
        name,
        paste(channels[!ok], collapse = ", ")
      ))
    }
  }

  if (anyFailure) {
    quit(status = 1L)
  }
  cat("\nOK: every BCF channel bitwise identical across every scenario\n")
} else {
  stop("unknown mode '", mode, "'; use record or compare")
}
