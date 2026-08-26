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
#   - the amplitudes (bartcoreForestAmplitudes: a, b0, b1),
#   - sigma,
#   - the reported result$train and result$varcount channels,
#   - the treatment forest's own split counts
#     (bartcoreForestVariableCounts 1).
# train and varcount are the only bitwise guard on storeSample's BCF report
# branches; the per-forest fits and glue guard the combining math and the
# coupling draw. result$varcount carries a FOREST AXIS (p x n.forests x
# n.samples here, single chain), slab 1 the prognostic forest and slab 2 the
# treatment forest, so it guards both forests' tree structure per draw;
# varcount.tau stays as the live-read cross-check against slab 2's last draw.
# Before bcf-bartcause-relocation D3 this channel was the prognostic forest
# alone and varcount.tau was the only guard on tau's splits.
#
# A getState -> setState continuation is deliberately NOT recorded: restore is
# structural by contract (the dropped accumulation history is not reproduced,
# test-bcf.R), so a bitwise assertion there would be a false gate.
#
# A scenario whose bitwise compare fails falls through to equivalence.R's
# statistical mode (docs/design/core-generalization.md): per-summary Welch
# z-scores over each statChannels entry's draws, against the baseline's
# recorded summaries, |z| >= 4 fails. A channel with no draws axis (a forest's
# raw fit/amplitude/varcount snapshot, a transaction's accept/reject verdict)
# has no statistical fallback and gates its mismatch directly, same as today.
# A baseline recorded before this mode existed carries no summaries and
# degrades loudly rather than silently passing.
#
# Usage:
#   Rscript bcf-equivalence.R record [out.rds]
#   Rscript bcf-equivalence.R compare baseline.rds
# Append 'quick' for a fast smoke pass (fewer draws; not comparable to a full
# baseline - the settings guard refuses the mixed comparison).
# Append '--cross-host' to compare a baseline recorded on ANOTHER machine.
# Bitwise is unavailable there by construction: makeData builds y through
# sin(), which is the platform libm, so the input data already differ in the
# last ulp before any engine code runs. The flag exempts the point-in-time
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

args <- commandArgs(trailingOnly = TRUE)
quick <- "quick" %in% args
args <- setdiff(args, "quick")
crossHost <- "--cross-host" %in% args
args <- setdiff(args, "--cross-host")
mode <- if (length(args) >= 1L) args[[1L]] else "record"
if (crossHost && mode != "compare") {
  stop("--cross-host applies to compare only; a recording is host-local")
}

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

# Channels a statistical compare can Welch-z: each carries a trailing draws
# axis (result$sigma/train/varcount, per bartcore_run's channel shape).
# Everything else recordChannels stores is a POINT-IN-TIME state query (a
# forest's raw fit/amplitude/varcount snapshot, a transaction's accept/reject
# verdict) with no repeated-draw axis to summarize.
statChannels <- c("sigma", "train", "varcount")

# The complement: every channel with no draws axis, and the only channels
# --cross-host exempts. These two vectors must partition the recorded names
# (minus summaries), so a channel added to either list has to be absent from
# the other; the compare loop re-checks the partition against what the
# baseline actually holds, which is where a channel appended at a call site
# rather than in recordChannels gets caught.
snapshotChannels <- c(
  "mu",
  "tau",
  "glue",
  "varcount.tau",
  "accepted",
  "installed"
)
stopifnot(length(intersect(snapshotChannels, statChannels)) == 0L)

# Tier 1's relative tolerance. Under a locked RNG stream two builds agree to
# rounding, so this is a stream-lock detector rather than a posterior test:
# loose enough that cross-host libm differences pass, ~1e6 times tighter than
# the tier-2 bar and so non-vacuous against every defect class the fixture
# exists to catch. The worst deviation an arm64-recorded baseline shows
# against an x86_64 Linux build sits about five orders of magnitude inside it.
crossHostRtol <- 1e-8

# A channel's draws laid out one row per cell, columns the trailing (draws)
# axis.
drawMatrix <- function(a) {
  d <- dim(a)
  matrix(a, ncol = if (is.null(d)) length(a) else d[length(d)])
}

# Posterior mean/var over a channel's trailing (draws) axis, one cell per
# leading-index combination - the reduction equivalence.R's fitSummaries
# performs by hand per channel, generalized over an arbitrary channel shape.
drawSummary <- function(a) {
  m <- drawMatrix(a)
  list(mean = rowMeans(m), var = apply(m, 1L, var), n = ncol(m))
}

# Per-cell lag-1 autocorrelation over the draws axis, 0 where a cell is
# constant.
lag1Acf <- function(m) {
  n <- ncol(m)
  if (n < 3L) {
    return(rep(0, nrow(m)))
  }
  cx <- m[, -n, drop = FALSE]
  cy <- m[, -1L, drop = FALSE]
  cx <- cx - rowMeans(cx)
  cy <- cy - rowMeans(cy)
  den <- sqrt(rowSums(cx^2) * rowSums(cy^2))
  ifelse(den > 0, rowSums(cx * cy) / den, 0)
}

# Tier 1 on a continuous channel: per-cell |a - b| against atol + rtol * |a|,
# atol taken from the channel's own scale in this scenario. Returns the ratio
# per cell, so <= 1 everywhere is the pass. A cell that agrees exactly scores 0
# even where the tolerance is 0; a non-finite difference scores Inf so it can
# never be read as a pass.
toleranceRatio <- function(x, y, rtol) {
  d <- abs(x - y)
  ratio <- ifelse(d == 0, 0, d / (rtol * max(abs(x)) + rtol * abs(x)))
  ratio[is.na(ratio)] <- Inf
  ratio
}

# Welch z per cell with an ESS-aware denominator. The draws axis here is ONE
# autocorrelated chain, not the independent seeds equivalence.R reduces over,
# so the nominal n understates the Monte Carlo error; the effective size comes
# from the pooled lag-1 autocorrelation, floored at 2 and capped at n so the
# denominator can only grow relative to the nominal statistic.
essCompare <- function(x, y) {
  mx <- drawMatrix(x)
  my <- drawMatrix(y)
  r <- 0.5 * (lag1Acf(mx) + lag1Acf(my))
  effective <- function(m) pmin(pmax(ncol(m) * (1 - r) / (1 + r), 2), ncol(m))
  ex <- effective(mx)
  ey <- effective(my)
  list(
    z = (rowMeans(mx) - rowMeans(my)) /
      sqrt(apply(mx, 1L, var) / ex + apply(my, 1L, var) / ey),
    ess = pmin(ex, ey)
  )
}

# One scenario's cross-host verdict, three name-prefixed lines: the exempt
# roster with a NON-GATING deviation diagnostic, the non-finite count, and the
# tier verdict. Tier 1 is the gate; tier 2 runs only after tier 1 fails and
# only adjudicates - did the posterior move, or did the stream merely decouple.
# Its bar is weak by construction (the line reports how weak), so a tier-2 pass
# is evidence the failure is not gross, never evidence the builds agree.
# Combinatorial channels are integer split counts: tier 1 compares them
# exactly, and tier 2 reports them diagnostically rather than through a z,
# because a Welch z on counts is a bitwise test in disguise and a structurally
# unsplit cell makes it 0/0.
compareCrossHost <- function(name, a, b) {
  exempt <- intersect(snapshotChannels, names(a))
  gated <- setdiff(names(a), c("summaries", exempt))
  for (ch in gated) {
    stopifnot(
      identical(dim(a[[ch]]), dim(b[[ch]])),
      length(a[[ch]]) == length(b[[ch]])
    )
  }

  n.differ <- 0L
  rel.dev <- 0
  for (ch in exempt) {
    if (!identical(a[[ch]], b[[ch]])) {
      n.differ <- n.differ + 1L
    }
    av <- unlist(a[[ch]])
    bv <- unlist(b[[ch]])
    if (is.double(av) && length(av) == length(bv) && max(abs(av)) > 0) {
      rel.dev <- max(rel.dev, max(abs(av - bv)) / max(abs(av)))
    }
  }
  cat(sprintf(
    "%-14s exempt (cross-host): %s - %d snapshot channels skipped by design [%d of %d differ, max rel dev %.1e]\n",
    name,
    paste(exempt, collapse = ", "),
    length(exempt),
    n.differ,
    length(exempt),
    rel.dev
  ))

  # Exempting the point-in-time channels deletes the only guard that saw a
  # NaN/Inf draw: every z built from one is NaN, and a max over NaNs with
  # na.rm reports -Inf and passes. Count them directly instead.
  nf <- vapply(
    gated,
    function(ch) c(sum(!is.finite(a[[ch]])), sum(!is.finite(b[[ch]]))),
    numeric(2L)
  )
  failed <- FALSE
  if (sum(nf) == 0) {
    cat(sprintf("%-14s NON-FINITE: none\n", name))
  } else {
    failed <- TRUE
    cat(sprintf(
      "%-14s NON-FINITE: %s <- FAIL\n",
      name,
      paste(nonFiniteParts(gated, nf), collapse = ", ")
    ))
  }

  is.continuous <- vapply(gated, function(ch) is.double(a[[ch]]), logical(1L))
  cont <- gated[is.continuous]
  comb <- gated[!is.continuous]
  ratios <- lapply(cont, function(ch) {
    toleranceRatio(a[[ch]], b[[ch]], crossHostRtol)
  })
  names(ratios) <- cont
  worst <- vapply(ratios, max, numeric(1L))
  cell <- vapply(ratios, which.max, integer(1L))
  comb.differ <- vapply(comb, function(ch) sum(a[[ch]] != b[[ch]]), numeric(1L))
  comb.n <- vapply(comb, function(ch) length(a[[ch]]), numeric(1L))

  bad.cont <- worst > 1
  bad.comb <- comb.differ > 0
  if (!any(bad.cont) && !any(bad.comb)) {
    cat(sprintf(
      "%-14s tier 1 PASS: max dev ratio %.1e (%s)%s%s\n",
      name,
      max(c(0, worst)),
      paste(sprintf("%s %.1e", cont, worst), collapse = ", "),
      if (length(comb) > 0L) {
        paste0(", ", paste(sprintf("%s identical", comb), collapse = ", "))
      } else {
        ""
      },
      if (all(worst == 0)) {
        sprintf(" - all %d GATED channels bitwise identical", length(gated))
      } else {
        ""
      }
    ))
    return(list(fail = failed, decoupled = FALSE))
  }

  cat(sprintf(
    "%-14s tier 1 FAIL: %s\n",
    name,
    paste(
      c(
        sprintf(
          "%s ratio %.1e (cell %d)",
          cont[bad.cont],
          worst[bad.cont],
          cell[bad.cont]
        ),
        sprintf(
          "%s %.0f of %.0f cells differ",
          comb[bad.comb],
          comb.differ[bad.comb],
          comb.n[bad.comb]
        )
      ),
      collapse = ", "
    )
  ))

  tier2 <- lapply(cont, function(ch) essCompare(a[[ch]], b[[ch]]))
  names(tier2) <- cont
  z <- unlist(lapply(tier2, function(s) s$z))
  ess <- unlist(lapply(tier2, function(s) s$ess))
  n.warn <- sum(abs(z) > 3, na.rm = TRUE)
  n.fail <- sum(abs(z) > 4, na.rm = TRUE)
  cat(sprintf(
    "%-14s decoupled: statistical - %d summaries, ESS-adjusted max |z| = %.2f (weak bar: tolerates %.2f posterior sd)%s%s%s\n",
    name,
    length(z),
    max(abs(z), na.rm = TRUE),
    4 * sqrt(2 / median(ess, na.rm = TRUE)),
    paste0(
      c(
        "",
        sprintf(
          "%s diagnostic-only, %.0f of %.0f cells differ, worst %.0f",
          comb[bad.comb],
          comb.differ[bad.comb],
          comb.n[bad.comb],
          vapply(
            comb[bad.comb],
            function(ch) max(abs(as.numeric(a[[ch]]) - as.numeric(b[[ch]]))),
            numeric(1L)
          )
        )
      ),
      collapse = "; "
    ),
    if (n.warn > 0L) sprintf(", %d with |z| > 3", n.warn) else "",
    if (n.fail > 0L) sprintf(", %d with |z| > 4 <- FAIL", n.fail) else ""
  ))
  if (n.fail > 0L) {
    failed <- TRUE
    offenders <- which(abs(z) > 4)
    cat(
      "  worst offenders:",
      paste0(
        names(z)[offenders],
        " (z=",
        round(z[offenders], 2L),
        ")",
        collapse = ", "
      ),
      "\n"
    )
  }
  list(fail = failed, decoupled = TRUE)
}

# Per-channel non-finite counts, zeros included so the line names every gated
# channel and a reader can see which one carries them.
nonFiniteParts <- function(gated, nf) {
  vapply(
    seq_along(gated),
    function(i) {
      n <- sum(nf[, i])
      if (n == 0) {
        return(sprintf("%s 0", gated[i]))
      }
      sprintf(
        "%s %.0f cells (%s)",
        gated[i],
        n,
        if (nf[1L, i] > 0 && nf[2L, i] > 0) {
          "baseline and this run"
        } else if (nf[1L, i] > 0) {
          "baseline"
        } else {
          "this run"
        }
      )
    },
    ""
  )
}

# The full recorded channel set for one BCF sampler at its current state,
# plus the draws-axis summaries a statistical compare needs (statChannels
# entries only).
recordChannels <- function(bcSampler, result) {
  ch <- list(
    mu = bartcoreForestFits(bcSampler, 0L),
    tau = bartcoreForestFits(bcSampler, 1L),
    glue = bartcoreForestAmplitudes(bcSampler),
    sigma = result$sigma,
    train = result$train,
    varcount = result$varcount,
    varcount.tau = bartcoreForestVariableCounts(bcSampler, 1L)
  )
  ch$summaries <- lapply(ch[statChannels], drawSummary)
  ch
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
    res <- bartcoreRun(bc, n.burn, n.samples)
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
    res <- bartcoreRun(bc, n.burn, n.samples)
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
    res <- bartcoreRun(bc, n.burn, n.samples)
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
    res <- bartcoreRun(bc, n.burn, n.samples)
    result$weighted <- recordChannels(bc, res)
  }

  # (e) setForestBasis: run, swap the treatment forest's basis, run again,
  # record the post-mutation state - the only bitwise guard on the sole
  # basis-mutation route.
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
    bartcoreRun(bc, n.burn, n.samples)
    bartcoreSetForestBasis(bc, 1L, cbind(1 - z2, z2))
    res <- bartcoreRun(bc, n.burn, n.samples)
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
    bartcoreRun(bc, n.burn, n.samples)
    bartcoreSetPredictor(bc, x2, forceUpdate = TRUE)
    res <- bartcoreRun(bc, n.burn, n.samples)
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
    bartcoreRun(bc, n.burn, n.samples)
    accepted <- bartcoreSetPredictor(bc, x2)
    res <- bartcoreRun(bc, n.burn, n.samples)
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
    bartcoreRun(bc, n.burn, n.samples)
    accepted <- bartcoreUpdatePredictor(bc, v, 3L)
    res <- bartcoreRun(bc, n.burn, n.samples)
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
    bartcoreRun(bc, n.burn, n.samples)
    accepted <- bartcoreUpdatePredictor(
      bc,
      ifelse(seq_len(n) %% 2L == 0L, 0.25, 0.75),
      1L
    )
    res <- bartcoreRun(bc, n.burn, n.samples)
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
    bartcoreRun(bc, n.burn, n.samples)
    installed <- bartcoreUpdatePredictorPerObservation(bc, v, 3L)
    res <- bartcoreRun(bc, n.burn, n.samples)
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
    bartcoreRun(bc, n.burn, n.samples)
    installed <- bartcoreUpdatePredictorPerObservation(
      bc,
      ifelse(seq_len(n) %% 2L == 0L, 0.25, 0.75),
      1L
    )
    res <- bartcoreRun(bc, n.burn, n.samples)
    c(recordChannels(bc, res), list(installed = installed))
  })

  # (l) active-row mask (docs/plans/latent-subset-mask.md) on the shared
  # gaussian response: a BCF chain composes the mask into its
  # GaussianResponse (the sigma df recount) and into composeForestWeights'
  # fan-out across BOTH forests, with no combiner-level code of its own -
  # the base ForestCombiner's setActiveRows is a no-op by design. Exercised
  # mid-chain through bartcoreSetActiveRows, same shape as (f)'s forced
  # setPredictor. Seeds are LITERAL, kept out of the guarded `seeds` vector
  # as (f)-(k)'s are, so settingsList() stays identical to the a825263
  # baseline and its neutrality compare still runs.
  result$masked <- local({
    d <- makeData(n, p, 8012L)
    set.seed(8112L)
    mask <- as.double(rbinom(n, 1L, 0.75))
    sampler <- dbarts(d$x, d$y, control = makeControl())
    set.seed(9012L)
    bc <- dbarts:::bartcoreBCFSampler(
      sampler,
      d$z,
      n.trees.treatment = n.trees.tau
    )
    bartcoreRun(bc, n.burn, n.samples)
    bartcoreSetActiveRows(bc, mask)
    res <- bartcoreRun(bc, n.burn, n.samples)
    recordChannels(bc, res)
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
  usedStatistical <- FALSE
  for (name in names(baseline$results)) {
    a <- baseline$results[[name]]
    b <- results[[name]]
    if (is.null(b)) {
      cat(sprintf("%-14s skipped (not produced this run)\n", name))
      next
    }
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
      verdict <- compareCrossHost(name, a, b)
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
        "%-14s identical (all %d channels: %s)\n",
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
        "%-14s statistical compare unsupported by this baseline (recorded before summaries)\n",
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
      "%-14s %d summaries, max |z| = %.2f%s%s\n",
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
        "%-14s NON-FINITE: %s <- FAIL\n",
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
        "%-14s also MISMATCH in %s (no statistical fallback) <- FAIL\n",
        name,
        paste(pointMismatch, collapse = ", ")
      ))
    }
  }

  if (anyFailure) {
    quit(status = 1L)
  }
  cat(
    if (crossHost && usedStatistical) {
      "\nOK: every gated BCF channel passed cross-host tier 1, or tier 2 could not distinguish it - a weak bar, so report any decoupled scenario\n"
    } else if (crossHost) {
      "\nOK: every gated BCF channel within the cross-host tier 1 bound (snapshot channels exempt)\n"
    } else if (usedStatistical) {
      "\nOK: every BCF channel identical, or statistically indistinguishable (|z| < 4)\n"
    } else {
      "\nOK: every BCF channel bitwise identical across every scenario\n"
    }
  )
} else {
  stop("unknown mode '", mode, "'; use record or compare")
}
