# The bartMultinomial generics that used to fall through to a "bart"/default
# method (print's per-chain arithmetic,
# residuals, plot, summary) or had no type argument (predict).

parseKeptDraws <- function(printed) {
  line <- grep("kept draws", printed, value = TRUE)
  as.integer(sub(".*: ", "", line))
}

set.seed(2201)
n <- 150L
p <- 3L
x <- matrix(runif(n * p), n, p)
eta <- cbind(2 * (x[, 1L] - 0.5), x[, 2L] - x[, 3L], 1.5 * (x[, 3L] - 0.5))
probs <- exp(eta) / rowSums(exp(eta))
labels <- vapply(
  seq_len(n),
  function(i) sample.int(3L, 1L, prob = probs[i, ]) - 1L,
  integer(1L)
)
y <- factor(c("lo", "mid", "hi")[labels + 1L], levels = c("lo", "mid", "hi"))

n.trees <- 10L
n.burn <- 5L
n.samples <- 7L

# ---- fix 1: print.bartMultinomial's "kept draws (per chain)" arithmetic.
# Under the buggy version, combineChains = TRUE (the default) reported the
# COMBINED total (n.samples * n.chains) instead of the per-chain count.

set.seed(3301)
fitCombined <- bart2(
  x,
  y,
  family = "multinomial",
  n.trees = n.trees,
  n.chains = 3L,
  n.threads = 1L,
  n.burn = n.burn,
  n.samples = n.samples,
  verbose = FALSE,
  combineChains = TRUE
)
set.seed(3301)
fitSplit <- bart2(
  x,
  y,
  family = "multinomial",
  n.trees = n.trees,
  n.chains = 3L,
  n.threads = 1L,
  n.burn = n.burn,
  n.samples = n.samples,
  verbose = FALSE,
  combineChains = FALSE
)
# same fit (same seed, same settings) reported two ways: both must read the
# per-chain count, never the combined total
expect_equal(parseKeptDraws(capture.output(print(fitCombined))), n.samples)
expect_equal(parseKeptDraws(capture.output(print(fitSplit))), n.samples)
expect_false(
  parseKeptDraws(capture.output(print(fitCombined))) == n.samples * 3L
)

# single chain: division by n.chains == 1 is a no-op, so this path was never
# buggy, but must stay correct
fitSingle <- bart2(
  x,
  y,
  family = "multinomial",
  n.trees = n.trees,
  n.chains = 1L,
  n.threads = 1L,
  n.burn = n.burn,
  n.samples = 11L,
  verbose = FALSE
)
expect_equal(parseKeptDraws(capture.output(print(fitSingle))), 11L)

# ---- fix 2: residuals.bartMultinomial. Previously fell to residuals.default
# and silently returned NULL.

res <- residuals(fitCombined)
expect_false(is.null(res))
expect_equal(dim(res), c(n, 3L))
expect_equal(colnames(res), c("lo", "mid", "hi"))
# 1[y = k] - phat_k sums to 0 across categories, both indicator and phat rows
# being simplexes
expect_true(all(abs(rowSums(res)) < 1e-8))
# matches a hand-rolled indicator-minus-fitted computation exactly
phatCombined <- fitted(fitCombined)
indicator <- matrix(0, n, 3L, dimnames = dimnames(phatCombined))
indicator[cbind(seq_len(n), match(y, levels(y)))] <- 1
expect_equal(res, indicator - phatCombined)

# grouped-count ingestion: the observed-proportion analog (y / rowSums(y))
set.seed(2202)
n.c <- 120L
x.c <- matrix(runif(n.c * p), n.c, p)
eta.c <- cbind(
  2 * (x.c[, 1L] - 0.5),
  x.c[, 2L] - x.c[, 3L],
  1.5 * (x.c[, 3L] - 0.5)
)
probs.c <- exp(eta.c) / rowSums(exp(eta.c))
trials.c <- sample(2:6, n.c, replace = TRUE)
counts.c <- t(vapply(
  seq_len(n.c),
  function(i) rmultinom(1L, trials.c[i], probs.c[i, ])[, 1L],
  integer(3L)
))
colnames(counts.c) <- c("lo", "mid", "hi")
fitCounts <- bart2(
  x.c,
  counts.c,
  family = "multinomial",
  n.trees = n.trees,
  n.chains = 1L,
  n.threads = 1L,
  n.burn = n.burn,
  n.samples = n.samples,
  verbose = FALSE
)
resCounts <- residuals(fitCounts)
expect_equal(dim(resCounts), c(n.c, 3L))
expect_equal(resCounts, counts.c / rowSums(counts.c) - fitted(fitCounts))

# ---- fix 3: plot.bartMultinomial. Previously fell to plot.default with a
# "'x' is a list, but does not have components 'x' and 'y'" error.

pdf(NULL)
expect_silent(plot(fitCombined))
expect_silent(plot(fitSingle))
dev.off()

# ---- fix 4: summary.bartMultinomial. Previously fell to summary.default
# (a raw structure dump).

sm <- summary(fitCombined)
expect_equal(class(sm), "summary.bart")
expect_equal(sort(sm$stats$variable), sort(paste0("meanProb[", levels(y), "]")))
expect_true(is.data.frame(sm$stats) || inherits(sm$stats, "tbl_df"))
printedSummary <- capture.output(print(sm))
expect_true(any(grepl("meanProb", printedSummary, fixed = TRUE)))

havePosterior <- requireNamespace("posterior", quietly = TRUE)
if (havePosterior) {
  expect_true(sm$posterior)
  expect_true("rhat" %in% names(sm$stats))
} else {
  expect_false(sm$posterior)
  expect_true(is.data.frame(sm$stats))
}

# degrade path, simulated the same way test-convergence-diagnostics.R does
if (havePosterior) {
  oldPosteriorAvailable <- dbarts:::posteriorAvailable
  tryCatch(
    {
      assignInNamespace("posteriorAvailable", function() FALSE, ns = "dbarts")
      smDegraded <- summary(fitCombined)
      expect_false(smDegraded$posterior)
      expect_true(is.data.frame(smDegraded$stats))
      expect_false("rhat" %in% names(smDegraded$stats))
    },
    finally = assignInNamespace(
      "posteriorAvailable",
      oldPosteriorAvailable,
      ns = "dbarts"
    )
  )
  rm(oldPosteriorAvailable, smDegraded)
}

# ---- fix 5: predict.bartMultinomial's new 'type' argument.

x.test <- x[seq_len(20L), , drop = FALSE]
set.seed(3302)
fitKeep <- bart2(
  x,
  y,
  family = "multinomial",
  test = x.test,
  keepTrees = TRUE,
  n.trees = n.trees,
  n.chains = 1L,
  n.threads = 1L,
  n.burn = n.burn,
  n.samples = n.samples,
  verbose = FALSE
)

# the default type = "ev" is bitwise identical to a pre-change-style call
# (no type argument at all) and consumes no RNG state
seedBefore <- .Random.seed
predDefault <- predict(fitKeep, x.test)
predEv <- predict(fitKeep, x.test, type = "ev")
expect_identical(predDefault, predEv)
expect_identical(predEv, fitKeep$yhat.test)
expect_identical(.Random.seed, seedBefore)

# type = "ppd" draws valid categories, shaped like "ev" minus the K margin,
# and DOES consume RNG state (only when explicitly requested)
seedBeforePpd <- .Random.seed
predPpd <- predict(fitKeep, x.test, type = "ppd")
expect_false(identical(.Random.seed, seedBeforePpd))
expect_equal(dim(predPpd), dim(predEv)[-length(dim(predEv))])
expect_true(all(predPpd %in% seq_len(3L)))

# predict's ppd matches extract's ppd construction exactly: replaying the
# fit-time test rows reproduces the stored test channel bit for bit (the
# existing reproduction gate), so with the RNG aligned the same way, the two
# ppd draws - one built from the replayed channel, one from the stored one -
# must agree exactly as well
set.seed(4401)
ppdFromExtract <- extract(fitKeep, type = "ppd", sample = "test")
set.seed(4401)
ppdFromPredict <- predict(fitKeep, x.test, type = "ppd")
expect_identical(ppdFromExtract, ppdFromPredict)

# type is validated like every other type= argument in the package
expect_error(predict(fitKeep, x.test, type = "bart"), "'arg' should be one of")

rm(
  parseKeptDraws,
  n,
  p,
  x,
  eta,
  probs,
  labels,
  y,
  n.trees,
  n.burn,
  n.samples,
  fitCombined,
  fitSplit,
  fitSingle,
  res,
  phatCombined,
  indicator,
  n.c,
  x.c,
  eta.c,
  probs.c,
  trials.c,
  counts.c,
  fitCounts,
  resCounts,
  sm,
  printedSummary,
  havePosterior,
  x.test,
  fitKeep,
  seedBefore,
  predDefault,
  predEv,
  seedBeforePpd,
  predPpd,
  ppdFromExtract,
  ppdFromPredict
)
