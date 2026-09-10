#!/usr/bin/env Rscript

# Measures LinearGaussianLeaf::maxNumCovariates (src/bartcore/model.hpp): the
# most columns a linear leaf regression may designate, eight, because the
# per-node sufficient-statistic scratch is a fixed-size stack array. The
# facade refuses a ninth (facade.hpp, "at most 8 leaf covariates are
# supported"), so a designation above the cap cannot be timed at all: this
# script sweeps 4 and 8, then records the refusal at 12 and 16 with the
# message the engine gives.
#
# FIT TIME is msec per iteration; the leaf draw is O(q^3) in the designated
# count, so the 4-to-8 slope says what a wider leaf would cost. LEAF
# CONDITIONING is the reason the cap is not purely a scratch-size question: a
# leaf holding fewer than q + 1 members has a rank-deficient design and leans
# entirely on the ridge the normal leaf prior supplies. Reported as the share
# of terminal nodes with n <= q + 1 and the median leaf size in units of
# q + 1; a cap that bound would show that share still near zero at eight.
#
# Usage: Rscript benchmarks/R/constant-linear-leaf-covariates.R [quick]
# No baseline, no pass/fail exit status: an informational sweep, run by hand
# when the cap is revisited.

suppressPackageStartupMessages(library(dbarts))

args <- commandArgs(trailingOnly = TRUE)
quick <- "quick" %in% args

n <- if (quick) 1000L else 4000L
p <- 16L
nTrees <- if (quick) 10L else 25L
nSamples <- if (quick) 40L else 200L
nBurn <- if (quick) 20L else 100L
counts <- if (quick) c(4L, 8L, 12L) else c(4L, 8L, 12L, 16L)

set.seed(11)
x <- matrix(runif(n * p), n, p)
colnames(x) <- paste0("x", seq_len(p))
y <- x[, 1L] + 2 * x[, 2L] * x[, 3L] + rnorm(n, 0, 0.3)
df <- data.frame(x, y = y)

cat(sprintf(
  "n = %d, %d trees, %d kept draws after %d burn\n",
  n,
  nTrees,
  nSamples,
  nBurn
))
cat(sprintf(
  "%8s %12s %14s %14s %s\n",
  "columns",
  "msec/iter",
  "leaves<=q+1",
  "median_n/(q+1)",
  "status"
))
for (q in counts) {
  set.seed(3)
  result <- tryCatch(
    {
      start <- proc.time()[["elapsed"]]
      fit <- bart(
        y ~ .,
        df,
        node.prior = dbartsPriors$linear(paste0("x", seq_len(q))),
        n.trees = nTrees,
        n.samples = nSamples,
        n.burn = nBurn,
        n.chains = 1L,
        n.threads = 1L,
        keepTrees = TRUE,
        seed = 7L,
        verbose = FALSE
      )
      elapsed <- proc.time()[["elapsed"]] - start
      trees <- extract(fit, type = "trees")
      leafSizes <- trees$n[trees$var < 0L]
      list(
        msec = elapsed / (nSamples + nBurn) * 1000,
        deficient = mean(leafSizes <= q + 1L),
        relative = median(leafSizes) / (q + 1L),
        status = "fit"
      )
    },
    error = function(e) {
      list(
        msec = NA_real_,
        deficient = NA_real_,
        relative = NA_real_,
        status = paste0("REFUSED: ", conditionMessage(e))
      )
    }
  )
  cat(sprintf(
    "%8d %12.3f %14.4f %14.2f %s\n",
    q,
    result$msec,
    result$deficient,
    result$relative,
    result$status
  ))
}
