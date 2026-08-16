#!/usr/bin/env Rscript

# Stage 1 of change-move-fix (docs/plans/change-move-fix.md): measure, on stock
# fits, the no-op rate that variant (a) (propose the change rule from the
# UNRESTRICTED prior, auto-reject when a descendant becomes unsatisfiable) would
# incur. The current engine instead proposes from the descendant-valid menu, so
# no draw is ever wasted; the question is how much waste variant (a) trades for
# its clean birth/death cancellation.
#
# Per change proposal at a target node splitting on variable v the variant-(a)
# no-op probability is 1 - |Valid(v)| / |SI(v)|, where SI(v) is the ancestor-only
# split interval and Valid(v) the subset that keeps every same-variable
# descendant split satisfiable. Nodes whose two children are both leaves have
# Valid = SI exactly (nothing below to invalidate) -> zero waste; the waste
# lives at deeper nodes with a same-variable descendant split.
#
# The engine is instrumented (moves.hpp changefix_instr, log-only, reverts before
# commit): each change proposal appends a row and each move bumps a type counter,
# all read from data the move already computed - the ext_rng stream is untouched.
#   ordinal   : validCount = |Valid| (findGoodOrdinalRules good count),
#               siCount = |SI| (splitInterval size); the missing-direction coin
#               doubles both and cancels in the ratio, so it is dropped here.
#   categorical: attempts = rejection-loop iterations (truncated-geometric with
#               mean 1/validFraction), found, reachable-category count. The
#               per-proposal no-op estimate is 1 - found/attempts.
#
# Logging is gated on env vars set by this driver AFTER a silent burn-in, so only
# post-burn sweeps are counted; a fresh detail path per config resets the move
# counters. Requires the instrumented package to be installed.
#
# Usage: Rscript change-fix-instrumentation.R [quick]

suppressPackageStartupMessages(library(dbarts))

args <- commandArgs(trailingOnly = TRUE)
quick <- "quick" %in% args

nBurn <- if (quick) 100L else 500L
nSweeps <- if (quick) 100L else 1000L
seed <- 20260708L

# ---------------------------------------------------------------------------
# run one fit with logging confined to the post-burn sweeps, return the parsed
# per-proposal detail frame and the move-type counts
# ---------------------------------------------------------------------------

runLogged <- function(sampler, ntree) {
  detailPath <- tempfile("changefix-detail-", fileext = ".csv")
  summaryPath <- tempfile("changefix-summary-", fileext = ".csv")
  on.exit(
    {
      Sys.unsetenv("DBARTS_CHANGE_LOG")
      Sys.unsetenv("DBARTS_MOVE_SUMMARY")
      unlink(c(detailPath, summaryPath))
    },
    add = TRUE
  )

  Sys.unsetenv("DBARTS_CHANGE_LOG")
  Sys.unsetenv("DBARTS_MOVE_SUMMARY")
  invisible(sampler$run(nBurn, 0L)) # burn-in, NOT logged

  Sys.setenv(DBARTS_CHANGE_LOG = detailPath, DBARTS_MOVE_SUMMARY = summaryPath)
  invisible(sampler$run(0L, nSweeps)) # sampling, logged
  Sys.unsetenv("DBARTS_CHANGE_LOG")
  Sys.unsetenv("DBARTS_MOVE_SUMMARY")

  detail <- utils::read.csv(detailPath)
  summ <- utils::read.csv(summaryPath)
  list(detail = detail, counts = unlist(summ[1, ]), ntree = ntree)
}

# ---------------------------------------------------------------------------
# aggregate one config's log into the reported quantities
# ---------------------------------------------------------------------------

summarizeConfig <- function(name, res) {
  d <- res$detail
  d <- d[stats::complete.cases(d), , drop = FALSE] # defensive
  counts <- res$counts
  totalMoves <- sum(counts)
  nChange <- nrow(d)

  # per-proposal variant-(a) no-op probability
  ord <- d$kind == 0L
  cat <- d$kind == 1L
  noop <- numeric(nrow(d))
  # ordinals: exact 1 - |Valid|/|SI|; siCount >= 1 for an available variable,
  # validCount clamped at 0 (an empty good set means every prior draw is a no-op)
  vc <- pmax(d$validCount[ord], 0L)
  noop[ord] <- 1 - vc / d$siCount[ord]
  # categoricals: 1 - found/attempts (found=1,attempts=1 -> 0; aborts -> 1)
  noop[cat] <- 1 - d$found[cat] / d$attempts[cat]

  leafShare <- mean(d$childrenLeaves == 1L)
  deep <- d$childrenLeaves == 0L
  deepNoop <- noop[deep]

  q <- function(v, p) {
    if (length(v)) as.numeric(quantile(v, p, na.rm = TRUE)) else NA_real_
  }
  overallPerProp <- mean(noop) # includes zero-waste leaves
  changePerSweep <- nChange / nSweeps
  noopPerSweep <- overallPerProp * changePerSweep # exp. no-op changes / sweep

  # pooled categorical draw-level cross-check (successes / draws)
  catPooled <- if (any(cat)) {
    1 - sum(d$found[cat]) / sum(d$attempts[cat])
  } else {
    NA_real_
  }

  cat(sprintf("\n=== %s ===\n", name))
  cat(sprintf(
    "  moves: birth %.3f death %.3f swap %.3f change %.3f (n=%d)\n",
    counts["birth"] / totalMoves,
    counts["death"] / totalMoves,
    counts["swap"] / totalMoves,
    counts["change"] / totalMoves,
    totalMoves
  ))
  cat(sprintf(
    "  change proposals: %d  (%.1f per sweep, %.3f of moves)\n",
    nChange,
    changePerSweep,
    nChange / totalMoves
  ))
  cat(sprintf(
    "  frac ordinal / categorical: %.3f / %.3f\n",
    mean(ord),
    mean(cat)
  ))
  cat(sprintf("  leaf-children traffic share (zero-waste): %.3f\n", leafShare))
  cat(sprintf(
    "  deep-node no-op fraction: mean %.4f median %.4f p90 %.4f (n=%d)\n",
    if (length(deepNoop)) mean(deepNoop) else NA_real_,
    q(deepNoop, 0.5),
    q(deepNoop, 0.9),
    length(deepNoop)
  ))
  cat(sprintf(
    "  OVERALL variant-(a) no-op rate: %.4f per change proposal\n",
    overallPerProp
  ))
  cat(sprintf(
    "                                  %.4f no-op changes per sweep\n",
    noopPerSweep
  ))
  if (!is.na(catPooled)) {
    cat(sprintf("  [categorical pooled draw-level no-op: %.4f]\n", catPooled))
  }

  data.frame(
    config = name,
    ntree = res$ntree,
    moves = totalMoves,
    changeFracOfMoves = nChange / totalMoves,
    changePerSweep = changePerSweep,
    fracOrdinal = mean(ord),
    fracCategorical = mean(cat),
    leafChildrenShare = leafShare,
    deepNoopMean = if (length(deepNoop)) mean(deepNoop) else NA_real_,
    deepNoopMedian = q(deepNoop, 0.5),
    deepNoopP90 = q(deepNoop, 0.9),
    overallNoopPerProposal = overallPerProp,
    noopChangesPerSweep = noopPerSweep,
    catPooledNoop = catPooled,
    stringsAsFactors = FALSE
  )
}

baseControl <- function(ntree) {
  dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = ntree,
    n.samples = nSweeps,
    n.burn = nBurn,
    n.thin = 1L,
    keepTrees = FALSE,
    updateState = FALSE,
    seed = seed
  )
}

# ---------------------------------------------------------------------------
# data generators
# ---------------------------------------------------------------------------

makeFriedman <- function(n, p = 10L, s = 1L) {
  set.seed(s)
  x <- matrix(runif(n * p), n, p)
  y <- 10 *
    sin(pi * x[, 1] * x[, 2]) +
    20 * (x[, 3] - 0.5)^2 +
    10 * x[, 4] +
    5 * x[, 5] +
    rnorm(n)
  list(x = x, y = y)
}

makeMixed <- function(n, s = 2L) {
  set.seed(s)
  cont <- as.data.frame(matrix(runif(n * 5L), n, 5L))
  names(cont) <- paste0("c", 1:5)
  c4 <- factor(sample.int(4L, n, replace = TRUE))
  c8 <- factor(sample.int(8L, n, replace = TRUE))
  c16 <- factor(sample.int(16L, n, replace = TRUE))
  # continuous signal plus categorical effects that reward same-variable stacking
  eff8 <- (as.integer(c8) - 4.5) * 0.6
  eff16 <- (as.integer(c16) %% 4L) * 1.0
  y <- 10 *
    sin(pi * cont$c1 * cont$c2) +
    10 * cont$c3 +
    (as.integer(c4) == 1L) * 4 +
    eff8 +
    eff16 +
    rnorm(n)
  x <- cbind(cont, c4 = c4, c8 = c8, c16 = c16)
  list(x = x, y = y)
}

# ---------------------------------------------------------------------------
# run the config grid
# ---------------------------------------------------------------------------

rows <- list()
# treePriorExpr is an unevaluated call (quote(cgm(...))/quote(dart(...))) spliced
# inline so dbarts resolves the constructor in its own prior vocabulary env; the
# bare names cgm/dart are not exported and only bind there.
runOne <- function(name, x, y, ntree, treePriorExpr = NULL) {
  ctl <- baseControl(ntree)
  s <- if (is.null(treePriorExpr)) {
    dbarts(x, y, control = ctl)
  } else {
    eval(bquote(dbarts(x, y, control = ctl, tree.prior = .(treePriorExpr))))
  }
  res <- runLogged(s, ntree)
  rows[[length(rows) + 1L]] <<- summarizeConfig(name, res)
  invisible(NULL)
}

# (1) Friedman-1, continuous-only, defaults
for (n in c(1000L, 10000L)) {
  fr <- makeFriedman(n)
  for (m in c(75L, 200L)) {
    runOne(sprintf("Friedman n=%d m=%d", n, m), fr$x, fr$y, m)
  }
}

# (2) mixed 5 continuous + 3 categorical (4,8,16 levels), n=1e4, m=75
mx <- makeMixed(10000L)
runOne("Mixed n=1e4 m=75", mx$x, mx$y, 75L)

# (3) deep single-tree stress: ntree=1, power=1, n=1e3
fr3 <- makeFriedman(1000L)
runOne(
  "Deep single-tree n=1e3 ntree=1 power=1",
  fr3$x,
  fr3$y,
  1L,
  treePriorExpr = quote(cgm(power = 1, base = 0.95))
)

# (4) config (2) with DART enabled
runOne(
  "Mixed+DART n=1e4 m=75",
  mx$x,
  mx$y,
  75L,
  treePriorExpr = quote(dart(power = 2, base = 0.95))
)

# ---------------------------------------------------------------------------
# combined table
# ---------------------------------------------------------------------------

tab <- do.call(rbind, rows)
cat("\n================ STAGE 1 SUMMARY ================\n")
print(tab, digits = 4, row.names = FALSE)

outCsv <- file.path("benchmarks", "results", "change-fix-instrumentation.csv")
dir.create(dirname(outCsv), showWarnings = FALSE, recursive = TRUE)
utils::write.csv(tab, outCsv, row.names = FALSE)
cat(sprintf("\nwrote %s\n", outCsv))
