# Internal bartcore engine surface (R/bartcore.R, src/bartcore/); the
# statistical-equivalence gates live outside the package in benchmarks/.

set.seed(99)
n <- 200L
p <- 5L
x <- matrix(runif(n * p), n, p)
f <- 10 * sin(pi * x[, 1L] * x[, 2L]) + 5 * x[, 4L]
y <- f + rnorm(n)
x.test <- matrix(runif(10L * p), 10L, p)

control <- dbartsControl(n.chains = 1L, n.threads = 1L,
                         n.trees = 50L, updateState = FALSE)
sampler <- dbarts(x, y, test = x.test, control = control)
bcSampler <- dbarts:::bartcoreSampler(sampler)

result <- dbarts:::bartcoreRun(bcSampler, 100L, 200L)

expect_equal(dim(result$train), c(n, 200L))
expect_equal(dim(result$test), c(10L, 200L))
expect_equal(dim(result$varcount), c(p, 200L))
expect_true(all(result$sigma > 0))

# the fit should explain most of the signal
fitMean <- rowMeans(result$train)
expect_true(mean((fitMean - f)^2) < 0.25 * mean((mean(y) - f)^2))

# embedded-Gibbs pattern: mutate offset between single draws
dbarts:::bartcoreSetOffset(bcSampler, rep(0.5, n))
result.offset <- dbarts:::bartcoreRun(bcSampler, 0L, 1L)
expect_equal(dim(result.offset$train), c(n, 1L))

dbarts:::bartcoreSetResponse(bcSampler, y + 1)
result.response <- dbarts:::bartcoreRun(bcSampler, 0L, 1L)
expect_true(mean(result.response$train) > mean(result.offset$train))

# no latents for a continuous response
expect_null(dbarts:::bartcoreGetLatents(bcSampler))

# binary response: probit latents match the response's signs; k pinned so
# this exercises the fixed-k path (the default binary chi hyperprior is
# tested below)
y.binary <- rbinom(n, 1L, pnorm(scale(f)))
sampler.binary <- dbarts(x, y.binary, control = control,
                         node.prior = normal(2))
bcSampler.binary <- dbarts:::bartcoreSampler(sampler.binary)
invisible(dbarts:::bartcoreRun(bcSampler.binary, 50L, 1L))

latents <- dbarts:::bartcoreGetLatents(bcSampler.binary)
expect_equal(length(latents), n)
expect_true(all(latents[y.binary == 1L] > 0) && all(latents[y.binary == 0L] < 0))

# fixed-k runs return no k samples
expect_null(result$k)

# the default binary spec uses the chi hyperprior on k
sampler.chik <- dbarts(x, y.binary, control = control)
bcSampler.chik <- dbarts:::bartcoreSampler(sampler.chik)
result.chik <- dbarts:::bartcoreRun(bcSampler.chik, 50L, 30L)
expect_equal(length(result.chik$k), 30L)
expect_true(all(result.chik$k > 0) && sd(result.chik$k) > 0)

# predictor mutation; on dedicated copies of x, since column and
# per-observation updates write into the borrowed matrix in place
x.mut <- x + 0
sampler.mut <- dbarts(x.mut, y, control = control)
bcSampler.mut <- dbarts:::bartcoreSampler(sampler.mut)
invisible(dbarts:::bartcoreRun(bcSampler.mut, 100L, 1L))

# an identity swap is always accepted; degenerate predictors would empty a
# leaf in some tree and roll back
expect_true(dbarts:::bartcoreSetPredictor(bcSampler.mut, x.mut + 0))
x.degenerate <- matrix(0.5, n, p)
expect_false(dbarts:::bartcoreSetPredictor(bcSampler.mut, x.degenerate))
expect_false(dbarts:::bartcoreSetPredictor(bcSampler.mut, x.degenerate,
                                           updateCutPoints = TRUE))
result.mut <- dbarts:::bartcoreRun(bcSampler.mut, 0L, 5L)
expect_true(all(is.finite(result.mut$train)))

# column-subset update: tiny jitter is accepted, degenerate column rejected
x.jitter <- x.mut[, 2L] + rnorm(n, 0, 1e-4)
expect_true(dbarts:::bartcoreUpdatePredictor(bcSampler.mut, x.jitter, 2L))
expect_false(dbarts:::bartcoreUpdatePredictor(bcSampler.mut, rep(0.5, n), 2L))
expect_error(dbarts:::bartcoreUpdatePredictor(bcSampler.mut, rep(0.5, n),
                                              p + 1L))

# per-observation update: extreme values install except where an observation
# is the last occupant of a leaf
installed <- dbarts:::bartcoreUpdatePredictorPerObservation(bcSampler.mut,
                                                            rep(10, n), 1L)
expect_equal(length(installed), n)
expect_true(any(installed) && any(!installed))
result.perobs <- dbarts:::bartcoreRun(bcSampler.mut, 0L, 5L)
expect_true(all(is.finite(result.perobs$train)))

# forced degenerate update collapses emptied splits instead of rolling back
expect_true(dbarts:::bartcoreSetPredictor(bcSampler.mut, x.degenerate,
                                          forceUpdate = TRUE,
                                          updateCutPoints = TRUE))
result.forced <- dbarts:::bartcoreRun(bcSampler.mut, 0L, 2L)
expect_true(all(is.finite(result.forced$train)))

# joint per-observation update: one mask, all-or-none across samplers that
# share an index-aligned column
x.jointA <- x + 0
x.jointB <- x + 0
x.jointB[, 2L:p] <- matrix(runif(n * (p - 1L)), n, p - 1L)
y.jointB <- -2 * x.jointB[, 1L] + rnorm(n)
bcSampler.jointA <- dbarts:::bartcoreSampler(dbarts(x.jointA, y,
                                                    control = control))
bcSampler.jointB <- dbarts:::bartcoreSampler(dbarts(x.jointB, y.jointB,
                                                    control = control))
invisible(dbarts:::bartcoreRun(bcSampler.jointA, 100L, 1L))
invisible(dbarts:::bartcoreRun(bcSampler.jointB, 100L, 1L))

installed.joint <- dbarts:::bartcoreUpdatePredictorPerObservationJointly(
  list(bcSampler.jointA, bcSampler.jointB), rep(10, n), c(1L, 1L))
expect_equal(length(installed.joint), n)
expect_true(any(installed.joint) && any(!installed.joint))
expect_true(all(is.finite(
  dbarts:::bartcoreRun(bcSampler.jointA, 0L, 1L)$train)))
expect_true(all(is.finite(
  dbarts:::bartcoreRun(bcSampler.jointB, 0L, 1L)$train)))

# quantile cut points and heterogeneous n.cuts
x.quants <- x + 0
x.quants[, 3L] <- round(x.quants[, 3L], 1L)  # 11 levels -> 10 quantile cuts
control.quants <- dbartsControl(n.chains = 1L,
                                n.threads = 1L,
                                n.trees = 50L, updateState = FALSE,
                                useQuantiles = TRUE,
                                n.cuts = c(100L, 50L, 100L, 100L, 25L))
sampler.quants <- dbarts(x.quants, y, control = control.quants)
bcSampler.quants <- dbarts:::bartcoreSampler(sampler.quants)
result.quants <- dbarts:::bartcoreRun(bcSampler.quants, 100L, 100L)
expect_equal(dim(result.quants$train), c(n, 100L))
fitMean.quants <- rowMeans(result.quants$train)
expect_true(mean((fitMean.quants - f)^2) < 0.25 * mean((mean(y) - f)^2))

# a coarser column cannot refresh quantile cuts: refused before any change
expect_error(
  dbarts:::bartcoreUpdatePredictor(bcSampler.quants,
                                   round(x.quants[, 3L] * 2) / 2, 3L,
                                   updateCutPoints = TRUE),
  pattern = "induced cut points")

# explicit cut points: installing a coarse grid collapses orphaned splits
dbarts:::bartcoreSetCutPoints(bcSampler.mut, list(c(0.25, 0.5, 0.75)), 1L)
result.cuts <- dbarts:::bartcoreRun(bcSampler.mut, 0L, 2L)
expect_true(all(is.finite(result.cuts$train)))
expect_error(
  dbarts:::bartcoreSetCutPoints(bcSampler.mut, list(c(0.5, 0.5)), 1L),
  pattern = "strictly increasing")

# multiple chains: per-chain slabs with a trailing chain dimension, run on
# worker threads; each chain gets its own generator seeded from R's stream
control.chains <- dbartsControl(n.chains = 2L,
                                n.threads = 2L, n.trees = 50L,
                                updateState = FALSE)
sampler.chains <- dbarts(x + 0, y, test = x.test, control = control.chains)
bcSampler.chains <- dbarts:::bartcoreSampler(sampler.chains)
result.chains <- dbarts:::bartcoreRun(bcSampler.chains, 100L, 60L)
expect_equal(dim(result.chains$sigma), c(60L, 2L))
expect_equal(dim(result.chains$train), c(n, 60L, 2L))
expect_equal(dim(result.chains$test), c(10L, 60L, 2L))
expect_equal(dim(result.chains$varcount), c(p, 60L, 2L))
expect_true(all(result.chains$sigma > 0))
expect_false(identical(result.chains$sigma[, 1L], result.chains$sigma[, 2L]))
fitMean.chains <- rowMeans(result.chains$train, dims = 1L)
expect_true(mean((fitMean.chains - f)^2) < 0.25 * mean((mean(y) - f)^2))

# transactions span chains; the sampler stays runnable afterward
expect_false(dbarts:::bartcoreSetPredictor(bcSampler.chains,
                                           matrix(0.5, n, p)))
installed.chains <- dbarts:::bartcoreUpdatePredictorPerObservation(
  bcSampler.chains, rep(10, n), 1L)
expect_true(any(installed.chains) && any(!installed.chains))
expect_true(all(is.finite(
  dbarts:::bartcoreRun(bcSampler.chains, 0L, 2L)$train)))

# multi-chain binary: latents and k gain the chain dimension
sampler.chains.binary <- dbarts(x + 0, y.binary, control = control.chains)
bcSampler.chains.binary <- dbarts:::bartcoreSampler(sampler.chains.binary)
result.chains.binary <- dbarts:::bartcoreRun(bcSampler.chains.binary, 50L, 10L)
expect_equal(dim(result.chains.binary$k), c(10L, 2L))
latents.chains <- dbarts:::bartcoreGetLatents(bcSampler.chains.binary)
expect_equal(dim(latents.chains), c(n, 2L))

# the standard dbartsSampler surface
control.engine <- dbartsControl(n.chains = 1L,
                                n.threads = 1L, n.trees = 50L,
                                n.burn = 100L, n.samples = 100L)
x.engine <- x + 0
colnames(x.engine) <- paste0("x", seq_len(p))
sampler.engine <- dbarts(x.engine, y, test = x.test,
                         control = control.engine)
r.engine <- sampler.engine$run()
expect_equal(dim(r.engine$train), c(n, 100L))
expect_equal(dim(r.engine$test), c(10L, 100L))
expect_true(all(r.engine$sigma > 0))
expect_true(mean((rowMeans(r.engine$train) - f)^2) <
              0.25 * mean((mean(y) - f)^2))

# standard mutation methods route to the new engine
sampler.engine$setOffset(rep(0.5, n))
r.offset.engine <- sampler.engine$run(0L, 1L)
expect_equal(dim(r.offset.engine$train), c(n, 1L))
sampler.engine$setResponse(y + 1)
sampler.engine$setSigma(1.0)
expect_equal(length(sampler.engine$getSigmas()), 1L)

# column updates default to the transaction and resolve names
expect_true(sampler.engine$setPredictor(
  sampler.engine$data@x[, 2L] + rnorm(n, 0, 1e-4), "x2"))
# full-matrix setPredictor defaults to forceUpdate = TRUE and returns
# nothing; like the classic engine it replaces data@x wholesale
expect_null(sampler.engine$setPredictor(matrix(runif(n * p), n, p)))
expect_true(all(is.finite(sampler.engine$run(0L, 2L)$train)))
# a degenerate matrix rolls back when not forced, and data@x is untouched
x.before <- sampler.engine$data@x
expect_false(sampler.engine$setPredictor(matrix(0.5, n, p),
                                         forceUpdate = FALSE))
expect_identical(sampler.engine$data@x, x.before)
# partial updates return the per-observation mask
installed.engine <- sampler.engine$setPredictor(rep(10, n), 1L,
                                                forceUpdate = "partial")
expect_equal(length(installed.engine), n)

sampler.engine$setCutPoints(list(c(0.25, 0.5, 0.75)), 1L)
sampler.engine$setTestPredictor(matrix(runif(10L * p), 10L, p))
sampler.engine$setTestPredictor(runif(10L), 2L)
expect_true(all(is.finite(sampler.engine$run(0L, 2L)$test)))

# test data can be removed (bart2's burn-in does this) and later restored
sampler.engine$setTestPredictorAndOffset(NULL, NULL)
expect_null(sampler.engine$data@x.test)
expect_null(sampler.engine$run(0L, 2L)$test)
sampler.engine$setTestPredictorAndOffset(x.test, rep(0.25, 10L))
expect_equal(sampler.engine$data@offset.test, rep(0.25, 10L))
expect_true(all(is.finite(sampler.engine$run(0L, 2L)$test)))

# printTrees dumps the live trees in the classic engine's format, one line
# per node
print.live <- capture.output(sampler.engine$printTrees(1L))
expect_equal(length(print.live),
             nrow(sampler.engine$getTrees(treeNums = 1L)))
expect_true(all(grepl(
  "^ *n: [0-9]+ TBN: [01]{3} Avail: [01]+ (var: [0-9]+ |ave: )", print.live)))

# weights swap like the classic engine's: a pointer swap reflected in the
# data slot, nothing else touched
w.engine <- runif(n, 0.5, 1.5)
sampler.engine$setWeights(w.engine)
expect_equal(sampler.engine$data@weights, w.engine)
expect_true(all(is.finite(sampler.engine$run(0L, 2L)$train)))
expect_error(sampler.engine$setWeights(rep(-1, n)),
             pattern = "non-negative")

# test offsets: created from a scalar offset the test offset stays synced
# with the regular one; setTestOffset breaks the sync
sampler.offset <- dbarts(x.engine + 0, y, test = x.test, offset = 0.5,
                         control = control.engine)
expect_equal(sampler.offset$data@offset.test, rep(0.5, 10L))
invisible(sampler.offset$run(20L, 1L))
sampler.offset$setOffset(1.5)
expect_equal(sampler.offset$data@offset.test, rep(1.5, 10L))
offset.explicit <- seq_len(10L) / 10
sampler.offset$setTestOffset(offset.explicit)
expect_identical(sampler.offset$data@testUsesRegularOffset, FALSE)
expect_equal(sampler.offset$data@offset.test, offset.explicit)
expect_error(sampler.offset$setTestOffset(rep(0, 3L)), pattern = "length")

# recorded test fits carry exactly the test offset: the same creation seed
# with and without an offset differs by it alone
set.seed(11)
sampler.to1 <- dbarts(x.engine + 0, y, test = x.test,
                      control = control.engine)
set.seed(11)
sampler.to2 <- dbarts(x.engine + 0, y, test = x.test,
                      control = control.engine)
sampler.to2$setTestOffset(rep(2, 10L))
r.to1 <- sampler.to1$run(20L, 3L)
r.to2 <- sampler.to2$run(20L, 3L)
expect_identical(r.to1$train, r.to2$train)
expect_identical(r.to1$test + 2, r.to2$test)

# replacing test predictors alone keeps the offset only when the row count
# still matches; the combined setter changes both
expect_silent(sampler.to2$setTestPredictor(matrix(runif(10L * p), 10L, p)))
expect_error(sampler.to2$setTestPredictor(matrix(runif(5L * p), 5L, p)),
             pattern = "together")
sampler.to2$setTestPredictorAndOffset(matrix(runif(5L * p), 5L, p),
                                      rep(1, 5L))
expect_equal(sampler.to2$data@offset.test, rep(1, 5L))
expect_equal(dim(sampler.to2$run(0L, 2L)$test), c(5L, 2L))

# thinning counts burn-in and samples at the kept rate
control.thin <- dbartsControl(n.chains = 1L,
                              n.threads = 1L, n.trees = 25L, n.thin = 2L)
sampler.thin <- dbarts(x + 0, y, control = control.thin)
r.thin <- sampler.thin$run(50L, 20L)
expect_equal(dim(r.thin$train), c(n, 20L))

# multiple chains through the flag: chain-dimensioned results and latents
control.engine2 <- dbartsControl(n.chains = 2L,
                                 n.threads = 2L, n.trees = 50L)
sampler.engine2 <- dbarts(x + 0, y.binary, control = control.engine2)
r.engine2 <- sampler.engine2$run(50L, 10L)
expect_equal(dim(r.engine2$k), c(10L, 2L))
expect_equal(dim(sampler.engine2$getLatents()), c(n, 2L))
expect_equal(length(sampler.engine2$getSigmas()), 2L)

# the package-level joint updater over full samplers
x.jointC <- x + 0
colnames(x.jointC) <- paste0("x", seq_len(p))
sampler.jointC <- dbarts(x.jointC, y, control = control.engine)
sampler.jointD <- dbarts(x.jointC + 0, y, control = control.engine)
invisible(sampler.jointC$run(50L, 1L))
invisible(sampler.jointD$run(50L, 1L))
installed.engine.joint <- updatePredictorPerObservationJointly(
  list(sampler.jointC, sampler.jointD), rep(10, n), "x1")
expect_equal(length(installed.engine.joint), n)

# whole-data replacement: cut points rebuild, splits remap, observation
# counts may change, and sigma holds on the original scale
control.setdata <- dbartsControl(n.chains = 1L,
                                 n.threads = 1L, n.trees = 50L)
sampler.setdata <- dbarts(x + 0, y, control = control.setdata)
invisible(sampler.setdata$run(100L, 1L))

sigma.before <- sampler.setdata$getSigmas()
sampler.setdata$setData(dbartsData(x + 0, y))
expect_equal(sampler.setdata$getSigmas(), sigma.before)
sampler.setdata$setData(dbartsData(x + 0, 2 * y + 3))
expect_equal(sampler.setdata$getSigmas(), sigma.before, tolerance = 1e-10)

n2 <- 150L
x2 <- matrix(runif(n2 * p), n2, p)
y2 <- 10 * sin(pi * x2[, 1L] * x2[, 2L]) + 5 * x2[, 4L] + rnorm(n2)
w2 <- runif(n2, 0.5, 1.5)
sampler.setdata$setData(dbartsData(x2, y2, test = x.test, weights = w2))
r.setdata <- sampler.setdata$run(0L, 5L)
expect_equal(dim(r.setdata$train), c(n2, 5L))
expect_equal(dim(r.setdata$test), c(10L, 5L))
expect_true(all(is.finite(r.setdata$train)))

# invalid replacements error and roll the data slot back
y.saved <- sampler.setdata$data@y
expect_error(sampler.setdata$setData(dbartsData(x2[, 1:2], y2)),
             pattern = "same predictors")
expect_identical(sampler.setdata$data@y, y.saved)
expect_true(all(is.finite(sampler.setdata$run(0L, 1L)$train)))

# whole-data replacement carries a test offset into the recorded test fits
sampler.setdata$setData(dbartsData(x2, y2 + 0.5, test = x.test,
                                   offset = 0.5))
expect_equal(sampler.setdata$data@offset.test, rep(0.5, 10L))
expect_true(all(is.finite(sampler.setdata$run(0L, 1L)$test)))

# binary samplers resize their latents with the data
sampler.setdata.binary <- dbarts(x + 0, y.binary, control = control.setdata)
invisible(sampler.setdata.binary$run(50L, 1L))
y2.binary <- rbinom(n2, 1L, 0.5)
sampler.setdata.binary$setData(dbartsData(x2, y2.binary))
expect_equal(length(sampler.setdata.binary$getLatents()), n2)
expect_true(all(is.finite(sampler.setdata.binary$run(0L, 2L)$train)))

# logistic via Polya-Gamma, reached through the internal helper's family
# argument; latents are the omega draws
sampler.logit.host <- dbarts(x, y.binary, control = control)
bcSampler.logit <- dbarts:::bartcoreSampler(sampler.logit.host,
                                            family = "logistic")
result.logit <- dbarts:::bartcoreRun(bcSampler.logit, 100L, 100L)
expect_equal(dim(result.logit$train), c(n, 100L))
expect_true(all(is.finite(result.logit$train)))
# log-odds fits classify well above the base rate
phat.logit <- plogis(rowMeans(result.logit$train))
expect_true(mean((phat.logit > 0.5) == (y.binary == 1L)) > 0.7)
omega <- dbarts:::bartcoreGetLatents(bcSampler.logit)
expect_equal(length(omega), n)
expect_true(all(omega > 0))

# family validation
expect_error(dbarts:::bartcoreSampler(sampler.logit.host, family = "cauchit"),
             pattern = "unrecognized response family")
expect_error(dbarts:::bartcoreSampler(sampler, family = "logistic"),
             pattern = "binary response")

# probit weights are rejected: a weighted probit has no tractable form
expect_error(dbarts(x, y.binary, weights = runif(n, 0.5, 1.5),
                    control = control),
             pattern = "probit models do not support weights")

# categorical predictors; no public surface marks matrix columns
# categorical (factors = "categorical" applies to data.frames), so the type
# is flipped on the data by hand
x.cat <- cbind(as.double(rep(0:3, length.out = n)), runif(n))
mu.cat <- c(2, -1, 3, 0)[x.cat[, 1L] + 1L]
y.cat <- mu.cat + 2 * x.cat[, 2L] + rnorm(n, 0, 0.5)
sampler.cat.host <- dbarts(x.cat, y.cat, control = control)
sampler.cat.host$data@varTypes[1L] <- 1L
bcSampler.cat <- dbarts:::bartcoreSampler(sampler.cat.host)
result.cat <- dbarts:::bartcoreRun(bcSampler.cat, 100L, 100L)
expect_true(all(is.finite(result.cat$train)))
# category means recovered after removing the known continuous effect
residual.means <- tapply(rowMeans(result.cat$train) - 2 * x.cat[, 2L],
                         x.cat[, 1L], mean)
expect_true(max(abs(residual.means - c(2, -1, 3, 0))) < 0.5)

# category codes outside the existing set are refused everywhere
expect_error(dbarts:::bartcoreUpdatePredictor(bcSampler.cat, rep(9, n), 1L),
             pattern = "existing category codes")
expect_error(dbarts:::bartcoreSetTestPredictor(
  bcSampler.cat, cbind(rep(7, 5L), runif(5L))),
  pattern = "existing category codes")
expect_error(dbarts:::bartcoreSetCutPoints(bcSampler.cat, list(0.5), 1L),
             pattern = "categorical predictor")

# valid categorical mutation routes through the mask logic
installed.cat <- dbarts:::bartcoreUpdatePredictorPerObservation(
  bcSampler.cat, as.double((x.cat[, 1L] + 1) %% 4), 1L)
expect_equal(length(installed.cat), n)
expect_true(all(is.finite(dbarts:::bartcoreRun(bcSampler.cat, 0L, 2L)$train)))

# non-integer codes are rejected at creation
x.cat.bad <- x.cat
x.cat.bad[1L, 1L] <- 0.5
sampler.cat.bad <- dbarts(x.cat.bad, y.cat, control = control)
sampler.cat.bad$data@varTypes[1L] <- 1L
expect_error(dbarts:::bartcoreSampler(sampler.cat.bad),
             pattern = "integer category codes")

# setData accepts categorical predictors (regression: a leftover ordinal-only
# check refused them) and refuses invalid codes; the internal handle skips
# the R5 slot fixups, so the type flip and cut counts carry over by hand
x.cat2 <- cbind(as.double(rep(c(1, 3, 0, 2), length.out = n)), runif(n))
data.cat2 <- dbartsData(x.cat2, mu.cat + 2 * x.cat2[, 2L])
data.cat2@n.cuts <- sampler.cat.host$data@n.cuts
data.cat2@sigma <- sampler.cat.host$data@sigma
data.cat2@varTypes[1L] <- 1L
dbarts:::bartcoreSetData(bcSampler.cat, data.cat2)
expect_true(all(is.finite(dbarts:::bartcoreRun(bcSampler.cat, 0L, 2L)$train)))
data.cat2@x[1L, 1L] <- 11
expect_error(dbarts:::bartcoreSetData(bcSampler.cat, data.cat2),
             pattern = "existing category codes")

# categories are capped at the code type's limit (65535); codes past 31
# exercise the high bits of the inline direction masks, and codes past 53
# the pooled tier (see test-data-categorical-wide.R)
n.wide <- 424L
x.wide <- matrix(as.double(rep(0:52, length.out = n.wide)), n.wide)
y.wide <- ifelse(x.wide[, 1L] >= 27, 2, 0) + rnorm(n.wide, 0, 0.3)
sampler.wide.host <- dbarts(x.wide, y.wide, control = control)
sampler.wide.host$data@varTypes[1L] <- 1L
bcSampler.wide <- dbarts:::bartcoreSampler(sampler.wide.host)
result.wide <- dbarts:::bartcoreRun(bcSampler.wide, 100L, 100L)
group.means <- tapply(rowMeans(result.wide$train), x.wide[, 1L] >= 27, mean)
expect_true(abs(group.means[[1L]]) < 0.3 && abs(group.means[[2L]] - 2) < 0.3)

x.over <- x.wide
x.over[1L, 1L] <- 65535
sampler.over.host <- dbarts(x.over, y.wide, control = control)
sampler.over.host$data@varTypes[1L] <- 1L
expect_error(dbarts:::bartcoreSampler(sampler.over.host),
             pattern = "codes in \\[0, 65535\\)")

# keepTrees through the flag: predictions from the saved trees reproduce the
# run's recorded test fits exactly
control.keep <- dbartsControl(n.chains = 1L,
                              n.threads = 1L, n.trees = 20L,
                              n.samples = 25L, n.burn = 50L,
                              keepTrees = TRUE, updateState = FALSE)
sampler.keep <- dbarts(x, y, test = x.test, control = control.keep)
result.keep <- sampler.keep$run()
predictions.keep <- sampler.keep$predict(x.test)
expect_equal(dim(predictions.keep), c(10L, 25L))
expect_identical(predictions.keep, unname(result.keep$test))

# getTrees: saved trees carry a sample dimension in the classic format
trees.saved <- sampler.keep$getTrees()
expect_equal(names(trees.saved), c("sample", "tree", "n", "var", "value"))
tree.first <- trees.saved[trees.saved$sample == 1L & trees.saved$tree == 1L, ]
expect_equal(tree.first$n[1L], n)
expect_equal(sum(tree.first$n[tree.first$var == -1L]), n)
expect_true(all(trees.saved$var == -1L |
                  (trees.saved$var >= 1L & trees.saved$var <= p)))

# live working trees have no sample dimension
trees.live <- sampler.keep$getTrees(current = TRUE)
expect_equal(names(trees.live), c("tree", "n", "var", "value"))

# newdata replays its rows through the trees for the counts
trees.newdata <- sampler.keep$getTrees(treeNums = 1L, sampleNums = 1L,
                                       newdata = x.test)
expect_equal(trees.newdata$n[1L], 10L)

# plotTree runs off the getTrees output
pdf(NULL)
expect_silent(sampler.keep$plotTree(1L))
dev.off()

# printTrees on a keepTrees sampler dumps the requested saved slots in the
# classic engine's saved-tree format
print.saved <- capture.output(
  sampler.keep$printTrees(treeNums = 1L, sampleNums = 1L))
expect_equal(length(print.saved),
             nrow(sampler.keep$getTrees(treeNums = 1L, sampleNums = 1L)))
expect_true(all(grepl(
  "^ *TBN: [01]{3}  (var: [0-9]+ ORDRule: |pred: )", print.saved)))

# prediction without keepTrees comes from the live trees, one set per chain
control.bc <- dbartsControl(n.chains = 1L,
                            n.threads = 1L, n.trees = 20L,
                            updateState = FALSE)
sampler.livepred <- dbarts(x, y, control = control.bc)
invisible(sampler.livepred$run(50L, 1L))
predictions.live <- sampler.livepred$predict(x.test)
expect_equal(length(predictions.live), 10L)
expect_equal(sampler.livepred$predict(x.test, offset.test = 2),
             predictions.live + 2)

# multiple chains add their dimension to predictions and a chain column to
# the trees
control.keep2 <- dbartsControl(n.chains = 2L,
                               n.threads = 1L, n.trees = 10L,
                               n.samples = 5L, keepTrees = TRUE,
                               updateState = FALSE)
sampler.keep2 <- dbarts(x, y, control = control.keep2)
invisible(sampler.keep2$run(30L, 5L))
expect_equal(dim(sampler.keep2$predict(x.test)), c(10L, 5L, 2L))
trees.keep2 <- sampler.keep2$getTrees()
expect_equal(names(trees.keep2),
             c("chain", "sample", "tree", "n", "var", "value"))

# multiple chains and samples get their header lines in the print dump
print.keep2 <- capture.output(
  sampler.keep2$printTrees(treeNums = 1L, sampleNums = 1:2))
expect_true(any(grepl("^chain [12]$", print.keep2)))
expect_true(any(grepl("sample [12]$", print.keep2)))

# state serialization: a restored sampler continues bitwise identically;
# multiple chains run on their own generators, so no seed sync is needed
control.state <- dbartsControl(n.chains = 2L,
                               n.threads = 1L, n.trees = 10L,
                               n.samples = 5L, updateState = FALSE)
sampler.state <- dbarts(x, y, control = control.state)
invisible(sampler.state$run(30L, 2L))
sampler.state$storeState()
state <- sampler.state$state
expect_inherits(state, "bartcoreState")

sampler.restored <- dbarts(x, y, control = control.state)
sampler.restored$setState(state)
expect_identical(sampler.state$run(0L, 3L), sampler.restored$run(0L, 3L))

# a single chain's state carries its generator like any other, so a restored
# sampler continues bitwise without any R-stream synchronization
sampler.state1 <- dbarts(x, y, control = control.bc)
invisible(sampler.state1$run(30L, 1L))
sampler.state1$storeState()
state1 <- sampler.state1$state
expect_true(length(sampler.state1$state[[1L]]$rng.state) > 0L)
sampler.restored1 <- dbarts(x, y, control = control.bc)
sampler.restored1$setState(state1)
result.a <- sampler.state1$run(0L, 2L)
result.b <- sampler.restored1$run(0L, 2L)
expect_identical(result.a, result.b)

# save/load: runs store state by default (updateState), and getPointer
# transparently re-creates the sampler from it after deserialization
control.us <- dbartsControl(n.chains = 2L,
                            n.threads = 1L, n.trees = 10L, n.samples = 5L)
sampler.us <- dbarts(x, y, control = control.us)
invisible(sampler.us$run(20L, 2L))
expect_inherits(sampler.us$state, "bartcoreState")
serialized <- tempfile(fileext = ".rds")
saveRDS(sampler.us, serialized)
result.a <- sampler.us$run(0L, 2L, updateState = FALSE)
sampler.loaded <- readRDS(serialized)
result.b <- sampler.loaded$run(0L, 2L, updateState = FALSE)
expect_identical(result.a, result.b)
unlink(serialized)

# malformed and inconsistent states are refused
expect_error(sampler.us$setState(list()), pattern = "bartcoreState")
state.bad <- sampler.us$state
internal.nodes <- which(state.bad[[1L]]$tree.vars > 0L)
expect_true(length(internal.nodes) > 0L)
state.bad[[1L]]$tree.values[internal.nodes[1L]] <-
  state.bad[[1L]]$tree.values[internal.nodes[1L]] + 1e-3
expect_error(sampler.us$setState(state.bad), pattern = "not consistent")

# prior sampling: tree structures and leaf parameters come from the CGM and
# node priors
control.prior.bc <- dbartsControl(n.chains = 1L,
                                  n.threads = 1L, n.trees = 50L,
                                  updateState = FALSE)
samplePrior <- function(sampler, numReplications) {
  leafCounts <- integer(0)
  leafValues <- numeric(0)
  for (r in seq_len(numReplications)) {
    sampler$sampleTreesFromPrior()
    sampler$sampleNodeParametersFromPrior()
    trees <- sampler$getTrees()
    leaves <- trees$var == -1L
    leafCounts <- c(leafCounts, tabulate(trees$tree[leaves], 50L))
    leafValues <- c(leafValues, trees$value[leaves])
  }
  list(counts = leafCounts, values = leafValues)
}
sampler.prior.bc <- dbarts(x, y, control = control.prior.bc)
set.seed(20)
prior.bc <- samplePrior(sampler.prior.bc, 10L)

# structures: the CGM(0.95, 2) prior's mean leaf count is about 2.5
expect_true(mean(prior.bc$counts) > 2.2 && mean(prior.bc$counts) < 2.8)

# parameters: leaf values match the node prior's spread
prior.sd <- 0.5 / (2 * sqrt(50))
expect_true(abs(sd(prior.bc$values) - prior.sd) < 0.15 * prior.sd)

# a sampler keeps running from a prior-drawn state
expect_true(all(is.finite(sampler.prior.bc$run(0L, 2L)$train)))

# setControl reconfigures a live sampler: the bart() dance of burning
# without training fits or stored trees, then flipping keepTrees on for
# the real samples
control.sc <- dbartsControl(n.chains = 1L,
                            n.threads = 1L, n.trees = 20L, n.samples = 8L,
                            n.burn = 50L, updateState = FALSE)
sampler.sc <- dbarts(x, y, test = x.test, control = control.sc)
control.burn <- control.sc
control.burn@keepTrainingFits <- FALSE
sampler.sc$setControl(control.burn)
r.burn <- sampler.sc$run(0L, 50L)
expect_null(r.burn$train)
expect_equal(length(r.burn$sigma), 50L)

control.keep <- control.sc
control.keep@keepTrees <- TRUE
sampler.sc$setControl(control.keep)
r.keep <- sampler.sc$run(0L, 8L)
expect_equal(dim(r.keep$train), c(n, 8L))
predictions.sc <- sampler.sc$predict(x.test)
expect_equal(dim(predictions.sc), c(10L, 8L))
expect_identical(predictions.sc, unname(r.keep$test))
expect_true(all(c("sample", "tree") %in% names(sampler.sc$getTrees())))

# turning storage back off frees it; predictions come from the live trees
sampler.sc$setControl(control.sc)
expect_equal(length(sampler.sc$predict(x.test)), 10L)

# settings fixed at creation refuse
control.bad <- control.sc
control.bad@n.trees <- 30L
expect_error(sampler.sc$setControl(control.bad), pattern = "n.trees")
control.bad <- control.sc
control.bad@n.chains <- 2L
expect_error(sampler.sc$setControl(control.bad), pattern = "n.chains")
control.bad <- control.sc
control.bad@rngSeed <- 71L
expect_error(sampler.sc$setControl(control.bad), pattern = "rngSeed")

# thinning through setControl matches creating with the rate
control.thin <- dbartsControl(n.chains = 1L,
                              n.threads = 1L, n.trees = 20L, n.samples = 5L,
                              n.thin = 3L, updateState = FALSE)
set.seed(31)
sampler.thin1 <- dbarts(x, y, control = control.thin)
control.nothin <- control.thin
control.nothin@n.thin <- 1L
set.seed(31)
sampler.thin2 <- dbarts(x, y, control = control.nothin)
sampler.thin2$setControl(control.thin)
r.thin1 <- sampler.thin1$run(5L, 5L)
r.thin2 <- sampler.thin2$run(5L, 5L)
expect_identical(r.thin1, r.thin2)

# setModel before running is indistinguishable from creating with the model
control.sm <- dbartsControl(n.chains = 1L,
                            n.threads = 1L, n.trees = 20L, n.samples = 5L,
                            updateState = FALSE)
set.seed(43)
sampler.sm1 <- dbarts(x, y, control = control.sm)
model.new <- sampler.sm1$model
model.new@tree.prior@power <- 1.5
model.new@tree.prior@base <- 0.8
model.new@node.hyperprior@k <- 3
sampler.sm1$setModel(model.new)
set.seed(43)
sampler.sm2 <- dbarts(x, y, control = control.sm,
                      tree.prior = cgm(1.5, 0.8), node.prior = normal(3))
r.sm1 <- sampler.sm1$run(20L, 5L)
r.sm2 <- sampler.sm2$run(20L, 5L)
expect_identical(r.sm1, r.sm2)

# a fixed residual prior holds sigma at sqrt(value), the documented
# variance semantics, at creation and through setModel alike
sampler.fix <- dbarts(x, y, resid.prior = fixed(4), control = control.sm)
r.fix <- sampler.fix$run(10L, 3L)
expect_equal(unique(r.fix$sigma), 2)
model.fixed <- sampler.sm1$model
model.fixed@resid.prior <- methods::new("dbartsFixedPrior", value = 9)
sampler.sm1$setModel(model.fixed)
expect_equal(unique(sampler.sm1$run(5L, 3L)$sigma), 3)
expect_equal(sampler.sm1$getSigmas(), 3)

# sums of squared residuals report on the original scale against the current
# (last recorded) fits
sampler.ssr <- dbarts(x, y, control = control.sm)
r.ssr <- sampler.ssr$run(20L, 1L)
expect_equal(sampler.ssr$getSumsOfSquaredResiduals(),
             sum((y - r.ssr$train[, 1L])^2))

# verbose keeps the classic engine's format: the creation summary, the loop
# header, the iteration counter, and the terminal summary
mkVerboseControl <- function() {
  dbartsControl(n.chains = 1L, n.threads = 1L,
                n.trees = 10L, n.samples = 5L, verbose = TRUE,
                printEvery = 2L, updateState = FALSE)
}
out.create.bc <- capture.output(
  sampler.vb <- dbarts(x, y, test = x.test, offset = 0.5,
                       control = mkVerboseControl()))
expect_true(any(out.create.bc == "Running BART with numeric y"))
expect_true(any(out.create.bc == "number of trees: 10"))
expect_true(any(grepl("^\tnumber of training observations: 200$",
                      out.create.bc)))
expect_true(any(grepl("^\treg : 0\\.50", out.create.bc)))

run.bc <- capture.output(invisible(sampler.vb$run(2L, 4L)))
expect_equal(sum(grepl("^iteration: ", run.bc)), 3L) # printEvery = 2, 6 total
expect_identical(run.bc[1L], "Running mcmc loop:")
expect_true(any(grepl("^total seconds in loop: ", run.bc)))
expect_true(any(run.bc == "Tree sizes, last iteration:"))
expect_true(any(run.bc == "DONE BART"))

# verbose toggles off through setControl
control.quiet <- mkVerboseControl()
control.quiet@verbose <- FALSE
sampler.vb$setControl(control.quiet)
expect_identical(capture.output(invisible(sampler.vb$run(0L, 2L))),
                 character(0))

# worker threads relay chain-prefixed progress through the main thread
control.threaded <- dbartsControl(n.chains = 2L,
                                  n.threads = 2L, n.trees = 10L,
                                  n.samples = 5L, verbose = TRUE,
                                  printEvery = 25L, updateState = FALSE)
out.threaded <- capture.output({
  sampler.vt <- dbarts(x, y, control = control.threaded)
  invisible(sampler.vt$run(25L, 25L))
})
expect_true(any(grepl("^\\[1\\] iteration: ", out.threaded)))
expect_true(any(grepl("^\\[2\\] iteration: ", out.threaded)))

# the state stores the response transform: an offset installed with
# updateScale moves it after creation, and a sampler restored over the same
# data must continue bitwise anyway (single-chain samplers share R's
# generator, so runs are seed-bracketed)
control.sc <- dbartsControl(n.chains = 1L, n.threads = 1L, n.trees = 25L,
                            updateState = FALSE)
sampler.sc <- dbarts(x, y, control = control.sc, sigma = 1)
bc.sc <- dbarts:::bartcoreSampler(sampler.sc)
set.seed(31)
invisible(dbarts:::bartcoreRun(bc.sc, 20L, 0L))
offset.sc <- 2 * sin(0.1 * seq_len(n))
dbarts:::bartcoreSetOffset(bc.sc, offset.sc, updateScale = TRUE)
invisible(dbarts:::bartcoreRun(bc.sc, 20L, 0L))
state.sc <- dbarts:::bartcoreStoreState(bc.sc)

# the same specification: an explicit sigma keeps the creation-time lm from
# folding the offset into a different variance prior
sampler.sc2 <- dbarts(x, y, offset = offset.sc, control = control.sc,
                      sigma = 1)
bc.sc2 <- dbarts:::bartcoreSampler(sampler.sc2)
dbarts:::bartcoreSetState(bc.sc2, state.sc)

set.seed(37)
run.sc <- dbarts:::bartcoreRun(bc.sc, 0L, 10L)
set.seed(37)
run.sc2 <- dbarts:::bartcoreRun(bc.sc2, 0L, 10L)
expect_identical(run.sc2$sigma, run.sc$sigma)
expect_identical(run.sc2$train, run.sc$train)

pred.sc <- dbarts:::bartcorePredict(bc.sc, x[1:5, , drop = FALSE])
pred.sc2 <- dbarts:::bartcorePredict(bc.sc2, x[1:5, , drop = FALSE])
expect_identical(pred.sc2, pred.sc)

# single-chain states gather into a multi-chain restore (the stan4bart
# pattern: fit chains separately, splice their states into one sampler for
# prediction); the R-stream-backed single-chain generators cannot restore
# into the dedicated multi-chain generators and are skipped
states.g <- lapply(1:3, function(i) {
  control.g <- dbartsControl(n.chains = 1L, n.threads = 1L, n.trees = 11L,
                             keepTrees = TRUE, n.samples = 6L,
                             updateState = FALSE)
  bc.g <- dbarts:::bartcoreSampler(dbarts(x, y, control = control.g,
                                          sigma = 1))
  invisible(dbarts:::bartcoreRun(bc.g, 7L, 6L))
  dbarts:::bartcoreStoreState(bc.g)
})
control.g3 <- dbartsControl(n.chains = 3L, n.threads = 1L, n.trees = 11L,
                            keepTrees = TRUE, n.samples = 6L,
                            updateState = FALSE)
bc.g3 <- dbarts:::bartcoreSampler(dbarts(x, y, control = control.g3,
                                         sigma = 1))
state.g <- states.g[[1L]]
state.g[[2L]] <- states.g[[2L]][[1L]]
state.g[[3L]] <- states.g[[3L]][[1L]]
expect_silent(dbarts:::bartcoreSetState(bc.g3, state.g))
pred.g <- dbarts:::bartcorePredict(bc.g3, x[1:4, , drop = FALSE])
expect_equal(dim(pred.g), c(4L, 6L, 3L))
# each restored chain predicts what its source sampler predicts
for (i in 1:3) {
  control.g1 <- dbartsControl(n.chains = 1L, n.threads = 1L, n.trees = 11L,
                              keepTrees = TRUE, n.samples = 6L,
                              updateState = FALSE)
  bc.g1 <- dbarts:::bartcoreSampler(dbarts(x, y, control = control.g1,
                                           sigma = 1))
  dbarts:::bartcoreSetState(bc.g1, states.g[[i]])
  expect_equal(dbarts:::bartcorePredict(bc.g1, x[1:4, , drop = FALSE]),
               pred.g[, , i, drop = FALSE][, , 1L])
}
