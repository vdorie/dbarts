# The public bart2(family = "multinomial") surface. The REPRODUCTION GATE
# below is the condition that matters most: bart2's fit path must reproduce,
# bit for bit, the internal bartcoreMultinomialSampler/bartcoreRun pattern
# benchmarks/R/multinomial-equivalence.R exercises, on the same data and
# seed. Everything else is level-threading, shape, and refusal coverage.

# The internal-path comparator: mirrors benchmarks/R/multinomial-equivalence.R
# exactly (host dbarts() sampler, then bartcoreMultinomialSampler, then one
# bartcoreRun call - no reseed in between, unlike the fixture script's
# cross-scenario isolation) so a single set.seed() before each of the two
# constructions is expected to agree bit for bit.
internalMultinomialFit <- function(
  x,
  labels,
  K,
  n.trees,
  n.burn,
  n.samples,
  test = NULL,
  offset = NULL
) {
  control <- dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = n.trees,
    updateState = FALSE
  )
  sampler <- if (is.null(test)) {
    dbarts(x, as.double(labels), control = control)
  } else {
    dbarts(x, as.double(labels), test = test, control = control)
  }
  bc <- dbarts:::bartcoreMultinomialSampler(
    sampler,
    labels,
    K = K,
    offset = offset
  )
  dbarts:::bartcoreRun(bc, n.burn, n.samples)
}

n.trees <- 20L
n.burn <- 12L
n.samples <- 12L

# --- K = 3, covariate-dependent (the reproduction gate) ---
set.seed(6301)
n <- 150L
p <- 4L
x3 <- matrix(runif(n * p), n, p)
eta <- cbind(
  2 * (x3[, 1L] - 0.5),
  x3[, 2L] - x3[, 3L],
  1.5 * (x3[, 4L] - 0.5)
)
probs3 <- exp(eta) / rowSums(exp(eta))
labels3 <- vapply(
  seq_len(n),
  function(i) sample.int(3L, 1L, prob = probs3[i, ]) - 1L,
  integer(1L)
)
y3 <- factor(c("lo", "mid", "hi")[labels3 + 1L], levels = c("lo", "mid", "hi"))

seed3 <- 7301L
set.seed(seed3)
fit3 <- bart2(
  x3,
  y3,
  family = "multinomial",
  n.trees = n.trees,
  n.chains = 1L,
  n.threads = 1L,
  n.burn = n.burn,
  n.samples = n.samples,
  verbose = FALSE
)
set.seed(seed3)
internal3 <- internalMultinomialFit(x3, labels3, 3L, n.trees, n.burn, n.samples)

probs3.fit <- fit3$yhat.train
dimnames(probs3.fit) <- NULL
expect_identical(probs3.fit, aperm(internal3$train, c(3L, 1L, 2L)))

# --- K = 2, the logistic-equivalent case ---
set.seed(6302)
n2 <- 120L
p2 <- 2L
x2 <- matrix(runif(n2 * p2), n2, p2)
labels2 <- rbinom(n2, 1L, plogis(2 * (x2[, 1L] - 0.5) + x2[, 2L]))
y2 <- factor(c("no", "yes")[labels2 + 1L], levels = c("no", "yes"))

seed2 <- 7302L
set.seed(seed2)
fit2 <- bart2(
  x2,
  y2,
  family = "multinomial",
  n.trees = n.trees,
  n.chains = 1L,
  n.threads = 1L,
  n.burn = n.burn,
  n.samples = n.samples,
  verbose = FALSE
)
set.seed(seed2)
internal2 <- internalMultinomialFit(x2, labels2, 2L, n.trees, n.burn, n.samples)

probs2.fit <- fit2$yhat.train
dimnames(probs2.fit) <- NULL
expect_identical(probs2.fit, aperm(internal2$train, c(3L, 1L, 2L)))

# --- shapes and probability sanity, both K = 3 and K = 2 fits ---
expect_equal(dim(fit3$yhat.train), c(n.samples, n, 3L))
expect_equal(dim(fit2$yhat.train), c(n.samples, n2, 2L))
expect_true(all(fit3$yhat.train >= 0 & fit3$yhat.train <= 1))
expect_true(all(fit2$yhat.train >= 0 & fit2$yhat.train <= 1))
expect_true(all(abs(apply(fit3$yhat.train, c(1L, 2L), sum) - 1) < 1e-8))
expect_true(all(abs(apply(fit2$yhat.train, c(1L, 2L), sum) - 1) < 1e-8))

# --- per-category varcount channel: shape, levels, reproduction ---
# varcount reshapes the run's per-sample per-category split-usage channel to
# n.samples x p x K with levels on the trailing K margin, mirroring yhat.train.
expect_equal(dim(fit3$varcount), c(n.samples, p, 3L))
expect_equal(dim(fit2$varcount), c(n.samples, p2, 2L))
expect_equal(dimnames(fit3$varcount)[[3L]], c("lo", "mid", "hi"))
expect_equal(dimnames(fit2$varcount)[[3L]], c("no", "yes"))
expect_true(is.integer(fit3$varcount))
# public == internal: fit$varcount strips to the internal bartcoreRun varcount
# channel (p x K x n.samples) bit for bit, the reproduction gate for varcount
vc3.fit <- fit3$varcount
dimnames(vc3.fit) <- NULL
expect_identical(vc3.fit, aperm(internal3$varcount, c(3L, 1L, 2L)))
vc2.fit <- fit2$varcount
dimnames(vc2.fit) <- NULL
expect_identical(vc2.fit, aperm(internal2$varcount, c(3L, 1L, 2L)))

# predictor names thread onto the p margin when x carries column names, as
# standard bart varcount does (an unnamed x leaves the margin NULL)
expect_null(dimnames(fit3$varcount)[[2L]])
set.seed(6303)
xn <- matrix(runif(n2 * p2), n2, p2, dimnames = list(NULL, c("aa", "bb")))
fitNamed <- bart2(
  xn,
  y2,
  family = "multinomial",
  n.trees = n.trees,
  n.chains = 1L,
  n.threads = 1L,
  n.burn = n.burn,
  n.samples = n.samples,
  verbose = FALSE
)
expect_equal(dimnames(fitNamed$varcount)[[2L]], c("aa", "bb"))
expect_equal(dimnames(fitNamed$varcount)[[3L]], c("no", "yes"))

# --- level threading: fitted colnames, extract dimnames, class factor levels ---
expect_equal(fit3$levels, c("lo", "mid", "hi"))
expect_equal(dimnames(fit3$yhat.train)[[3L]], c("lo", "mid", "hi"))

fitted3.ev <- fitted(fit3)
expect_equal(dim(fitted3.ev), c(n, 3L))
expect_equal(colnames(fitted3.ev), c("lo", "mid", "hi"))
expect_true(all(abs(rowSums(fitted3.ev) - 1) < 1e-8))

ev3 <- extract(fit3, type = "ev")
expect_equal(dimnames(ev3)[[3L]], c("lo", "mid", "hi"))

fitted3.class <- fitted(fit3, type = "class")
expect_true(is.factor(fitted3.class))
expect_equal(levels(fitted3.class), c("lo", "mid", "hi"))
expect_equal(length(fitted3.class), n)
# the class call is the argmax of the SAME posterior-mean matrix fitted()
# returns, so the two must agree observation-by-observation
expect_equal(
  as.character(fitted3.class),
  colnames(fitted3.ev)[max.col(fitted3.ev, ties.method = "first")]
)

# --- posterior-predictive class draws are valid categories ---
ppd3 <- extract(fit3, type = "ppd")
expect_equal(dim(ppd3), c(n.samples, n))
expect_true(all(ppd3 %in% seq_len(3L)))
ppd2 <- extract(fit2, type = "ppd")
expect_true(all(ppd2 %in% seq_len(2L)))

# --- multi-chain smoke: combineChains toggles the reported shape, both K ---
set.seed(431)
fitMulti <- bart2(
  x3,
  y3,
  family = "multinomial",
  n.trees = 10L,
  n.chains = 2L,
  n.threads = 2L,
  n.burn = 5L,
  n.samples = 5L,
  verbose = FALSE,
  combineChains = TRUE
)
expect_equal(dim(fitMulti$yhat.train), c(10L, n, 3L))
expect_equal(dim(fitMulti$varcount), c(10L, p, 3L))
set.seed(431)
fitMultiSplit <- bart2(
  x3,
  y3,
  family = "multinomial",
  n.trees = 10L,
  n.chains = 2L,
  n.threads = 2L,
  n.burn = 5L,
  n.samples = 5L,
  verbose = FALSE,
  combineChains = FALSE
)
expect_equal(dim(fitMultiSplit$yhat.train), c(2L, 5L, n, 3L))
expect_equal(dim(fitMultiSplit$varcount), c(2L, 5L, p, 3L))

# --- test data at creation: K = 3 reproduction gate, shape, levels ---
# test rows reuse train rows 1:20, so the test channel is well-defined and,
# by softmax invariance to the common level shift, matches those train columns
x3.test <- x3[seq_len(20L), , drop = FALSE]
set.seed(seed3)
fit3t <- bart2(
  x3,
  y3,
  family = "multinomial",
  test = x3.test,
  n.trees = n.trees,
  n.chains = 1L,
  n.threads = 1L,
  n.burn = n.burn,
  n.samples = n.samples,
  verbose = FALSE
)
set.seed(seed3)
internal3t <- internalMultinomialFit(
  x3,
  labels3,
  3L,
  n.trees,
  n.burn,
  n.samples,
  test = x3.test
)
# the public test channel reproduces the internal bartcoreRun test channel bit
# for bit (the reproduction gate extended to the test surface)
probs3t.fit <- fit3t$yhat.test
dimnames(probs3t.fit) <- NULL
expect_identical(probs3t.fit, aperm(internal3t$test, c(3L, 1L, 2L)))
expect_equal(dim(fit3t$yhat.test), c(n.samples, 20L, 3L))
expect_equal(dimnames(fit3t$yhat.test)[[3L]], c("lo", "mid", "hi"))
expect_true(all(fit3t$yhat.test >= 0 & fit3t$yhat.test <= 1))
expect_true(all(abs(apply(fit3t$yhat.test, c(1L, 2L), sum) - 1) < 1e-8))
# softmax invariance: test rows duplicate train rows 1:20, so their probabilities
# agree with those train columns to floating point
expect_true(
  max(abs(
    fit3t$yhat.test - fit3t$yhat.train[, seq_len(20L), , drop = FALSE]
  )) <
    1e-6
)
# extract(sample = "test") returns the test channel
expect_identical(
  extract(fit3t, type = "ev", sample = "test"),
  fit3t$yhat.test
)
# ppd over the test channel draws valid categories
ppd3t <- extract(fit3t, type = "ppd", sample = "test")
expect_equal(dim(ppd3t), c(n.samples, 20L))
expect_true(all(ppd3t %in% seq_len(3L)))

# --- test data at creation: K = 2 ---
x2.test <- x2[seq_len(15L), , drop = FALSE]
set.seed(seed2)
fit2t <- bart2(
  x2,
  y2,
  family = "multinomial",
  test = x2.test,
  n.trees = n.trees,
  n.chains = 1L,
  n.threads = 1L,
  n.burn = n.burn,
  n.samples = n.samples,
  verbose = FALSE
)
expect_equal(dim(fit2t$yhat.test), c(n.samples, 15L, 2L))
expect_equal(dimnames(fit2t$yhat.test)[[3L]], c("no", "yes"))
expect_true(all(fit2t$yhat.test >= 0 & fit2t$yhat.test <= 1))

# --- predict-on-newdata: keepTrees replays all K saved forests ---
set.seed(seed3)
fit3p <- bart2(
  x3,
  y3,
  family = "multinomial",
  test = x3.test,
  keepTrees = TRUE,
  n.trees = n.trees,
  n.chains = 1L,
  n.threads = 1L,
  n.burn = n.burn,
  n.samples = n.samples,
  verbose = FALSE
)
# THE REPRODUCTION GATE: predict at the fit-time test rows equals the recorded
# test channel BIT FOR BIT. Both replay the same saved trees through the same
# softmax, and neither carries the level-centering grand shift (afterCombine
# leaves totalTestFits and the saved leaves alone), so they agree exactly.
pred3p <- predict(fit3p, x3.test)
expect_identical(pred3p, fit3p$yhat.test)
# shape, levels, and simplex on a fresh newdata matrix
pred3new <- predict(fit3p, x3[21:40, , drop = FALSE])
expect_equal(dim(pred3new), c(n.samples, 20L, 3L))
expect_equal(dimnames(pred3new)[[3L]], c("lo", "mid", "hi"))
expect_true(all(pred3new >= 0 & pred3new <= 1))
expect_true(all(abs(apply(pred3new, c(1L, 2L), sum) - 1) < 1e-8))
# predict on the training rows reproduces the train channel to floating point:
# the train channel softmaxes totalFits WITH the common level shift, the replay
# softmaxes the saved (pre-shift) leaves; softmax invariance to that shift makes
# them agree up to rounding (not bitwise, unlike the shift-free test channel)
pred3train <- predict(fit3p, x3)
expect_equal(dim(pred3train), c(n.samples, n, 3L))
expect_true(max(abs(pred3train - fit3p$yhat.train)) < 1e-6)

# --- the retained $fit is a host shell: mutating it errors, reading does not ---
# it carries the design and priors the K-forest engine ($bc) was built from but
# none of the model, so every mutation through it used to be a silent no-op
hostRefusal <- "host sampler of a bart2"
expect_error(fit3p$fit$setResponse(rnorm(n)), hostRefusal)
expect_error(fit3p$fit$setOffset(rep(0.5, n)), hostRefusal)
expect_error(fit3p$fit$setSigma(2), hostRefusal)
expect_error(fit3p$fit$run(0L, 1L), hostRefusal)
expect_error(fit3p$fit$setPredictor(x3[, 1L], 1L), hostRefusal)
# the read surface predict() threads through is untouched
expect_equal(ncol(fit3p$fit$data@x), ncol(x3))
expect_identical(predict(fit3p, x3.test), pred3p)

# --- keepSampler retains $fit independent of keepTrees ---
fit3ks <- bart2(
  x3,
  y3,
  family = "multinomial",
  n.trees = n.trees,
  n.chains = 1L,
  n.threads = 1L,
  n.burn = n.burn,
  n.samples = n.samples,
  verbose = FALSE,
  keepSampler = TRUE,
  keepTrees = FALSE
)
expect_false(is.null(fit3ks$fit))
expect_true(is.null(fit3ks$bc))
rm(fit3ks)

# multi-chain predict threads the chain margin like the run channels
set.seed(431)
fit3pMulti <- bart2(
  x3,
  y3,
  family = "multinomial",
  keepTrees = TRUE,
  n.trees = 10L,
  n.chains = 2L,
  n.threads = 2L,
  n.burn = 5L,
  n.samples = 5L,
  verbose = FALSE
)
predMulti <- predict(fit3pMulti, x3.test)
expect_equal(dim(predMulti), c(10L, 20L, 3L))
predMultiSplit <- predict(fit3pMulti, x3.test, combineChains = FALSE)
expect_equal(dim(predMultiSplit), c(2L, 5L, 20L, 3L))

# --- count-matrix response: an n x K count matrix beside the factor
# path, both routed through bart2's multinomial branch. The internal-path
# comparator mirrors internalMultinomialFit above, substituting
# bartcoreMultinomialCountSampler for bartcoreMultinomialSampler.
internalMultinomialCountFit <- function(
  x,
  counts,
  K,
  n.trees,
  n.burn,
  n.samples,
  test = NULL,
  offset = NULL
) {
  control <- dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = n.trees,
    updateState = FALSE
  )
  sampler <- if (is.null(test)) {
    dbarts(x, as.double(counts[, 1L]), control = control)
  } else {
    dbarts(x, as.double(counts[, 1L]), test = test, control = control)
  }
  bc <- dbarts:::bartcoreMultinomialCountSampler(
    sampler,
    counts,
    K = K,
    offset = offset
  )
  dbarts:::bartcoreRun(bc, n.burn, n.samples)
}

# K = 3, grouped counts (n_i > 1): the reproduction gate extended to the
# count-matrix response form - a public count fit must reproduce the
# internal bartcoreMultinomialCountSampler channel bit for bit.
set.seed(6304)
n3c <- 150L
x3c <- matrix(runif(n3c * p), n3c, p)
eta3c <- cbind(
  2 * (x3c[, 1L] - 0.5),
  x3c[, 2L] - x3c[, 3L],
  1.5 * (x3c[, 4L] - 0.5)
)
probs3c <- exp(eta3c) / rowSums(exp(eta3c))
trials3c <- sample(2:6, n3c, replace = TRUE)
counts3c <- t(vapply(
  seq_len(n3c),
  function(i) rmultinom(1L, trials3c[i], probs3c[i, ])[, 1L],
  integer(3L)
))
colnames(counts3c) <- c("lo", "mid", "hi")

seed3c <- 7304L
set.seed(seed3c)
fit3c <- bart2(
  x3c,
  counts3c,
  family = "multinomial",
  n.trees = n.trees,
  n.chains = 1L,
  n.threads = 1L,
  n.burn = n.burn,
  n.samples = n.samples,
  verbose = FALSE
)
set.seed(seed3c)
internal3c <- internalMultinomialCountFit(
  x3c,
  counts3c,
  3L,
  n.trees,
  n.burn,
  n.samples
)

probs3c.fit <- fit3c$yhat.train
dimnames(probs3c.fit) <- NULL
expect_identical(probs3c.fit, aperm(internal3c$train, c(3L, 1L, 2L)))
vc3c.fit <- fit3c$varcount
dimnames(vc3c.fit) <- NULL
expect_identical(vc3c.fit, aperm(internal3c$varcount, c(3L, 1L, 2L)))

# colnames(Y) threads onto the K margin and onto $levels; the count matrix
# itself is stored as $y, as the factor path stores the factor as $y
expect_equal(fit3c$levels, c("lo", "mid", "hi"))
expect_equal(dimnames(fit3c$yhat.train)[[3L]], c("lo", "mid", "hi"))
expect_equal(dimnames(fit3c$varcount)[[3L]], c("lo", "mid", "hi"))
expect_identical(fit3c$y, counts3c)

# an unnamed count matrix falls back to as.character(seq_len(K)) index levels
unnamedCounts3c <- counts3c
colnames(unnamedCounts3c) <- NULL
set.seed(seed3c)
fitUnnamed3c <- bart2(
  x3c,
  unnamedCounts3c,
  family = "multinomial",
  n.trees = n.trees,
  n.chains = 1L,
  n.threads = 1L,
  n.burn = n.burn,
  n.samples = n.samples,
  verbose = FALSE
)
expect_equal(fitUnnamed3c$levels, c("1", "2", "3"))
expect_equal(dimnames(fitUnnamed3c$yhat.train)[[3L]], c("1", "2", "3"))

# generics work unchanged on a count fit: the softmax output is
# count-independent, so extract/fitted/predict need no special-casing
fitted3c.ev <- fitted(fit3c)
expect_equal(colnames(fitted3c.ev), c("lo", "mid", "hi"))
fitted3c.class <- fitted(fit3c, type = "class")
expect_true(is.factor(fitted3c.class))
expect_equal(levels(fitted3c.class), c("lo", "mid", "hi"))
ev3c <- extract(fit3c, type = "ev")
expect_identical(ev3c, fit3c$yhat.train)
ppd3c <- extract(fit3c, type = "ppd")
expect_true(all(ppd3c %in% seq_len(3L)))

set.seed(seed3c)
fit3cKeep <- bart2(
  x3c,
  counts3c,
  family = "multinomial",
  test = x3c[seq_len(10L), , drop = FALSE],
  keepTrees = TRUE,
  n.trees = n.trees,
  n.chains = 1L,
  n.threads = 1L,
  n.burn = n.burn,
  n.samples = n.samples,
  verbose = FALSE
)
pred3c <- predict(fit3cKeep, x3c[seq_len(10L), , drop = FALSE])
expect_identical(pred3c, fit3cKeep$yhat.test)

# print reports the count-matrix input the same way it reports a factor one
printed3c <- capture.output(print(fit3cKeep))
expect_true(any(grepl("levels: lo, mid, hi", printed3c, fixed = TRUE)))

# a one-hot count matrix (every row sum 1) reproduces the corresponding
# factor fit's probabilities bit for bit - the single-trial reduction,
# checked at the bart2 surface rather than the internal one above
onehot2 <- matrix(0L, n2, 2L, dimnames = list(NULL, c("no", "yes")))
onehot2[cbind(seq_len(n2), labels2 + 1L)] <- 1L
set.seed(seed2)
fitOneHot2 <- bart2(
  x2,
  onehot2,
  family = "multinomial",
  n.trees = n.trees,
  n.chains = 1L,
  n.threads = 1L,
  n.burn = n.burn,
  n.samples = n.samples,
  verbose = FALSE
)
expect_identical(fitOneHot2$yhat.train, fit2$yhat.train)
expect_identical(fitOneHot2$varcount, fit2$varcount)

# --- count-matrix validation: informative, R-side errors ---
badCounts <- onehot2 + 0L
negCounts <- badCounts
negCounts[1L, 1L] <- -1L
expect_error(
  bart2(x2, negCounts, family = "multinomial"),
  "non-negative"
)
fracCounts <- badCounts + 0.0
fracCounts[1L, 1L] <- 1.5
expect_error(
  bart2(x2, fracCounts, family = "multinomial"),
  "whole numbers"
)
zeroRowCounts <- badCounts
zeroRowCounts[1L, ] <- 0L
expect_error(
  bart2(x2, zeroRowCounts, family = "multinomial"),
  "row sum"
)
expect_error(
  bart2(x2, badCounts[, 1L, drop = FALSE], family = "multinomial"),
  "at least 2 columns"
)

# --- formula interface: a factor or cbind(...) count-matrix response
# through model.frame/model.response (no type coercion, so levels/column
# names survive), the right-hand side coded exactly as any other family's
# formula fit. A 3+-level factor response is also inferred under the default
# family = "auto" (checked below); a cbind(...) count matrix stays explicit.
x3Named <- x3
colnames(x3Named) <- c("x1", "x2", "x3", "x4")
df3 <- data.frame(x3Named, y3 = y3)

# factor LHS: reproduces the matrix-interface factor fit bit for bit at the
# same seed (the reproduction gate extended to the formula surface)
set.seed(seed3)
fit3Formula <- bart2(
  y3 ~ x1 + x2 + x3 + x4,
  data = df3,
  family = "multinomial",
  keepTrees = TRUE,
  n.trees = n.trees,
  n.chains = 1L,
  n.threads = 1L,
  n.burn = n.burn,
  n.samples = n.samples,
  verbose = FALSE
)
probs3Formula <- fit3Formula$yhat.train
dimnames(probs3Formula) <- NULL
probs3Matrix <- fit3$yhat.train
dimnames(probs3Matrix) <- NULL
expect_identical(probs3Formula, probs3Matrix)
expect_equal(fit3Formula$levels, c("lo", "mid", "hi"))

# predict on a data.frame newdata codes it through the retained terms
# (makeModelMatrix via sampler$data@x's attributes), matching predict on
# the equivalent already-coded matrix exactly
newdata3 <- as.data.frame(x3Named[seq_len(20L), , drop = FALSE])
predFromFrame <- predict(fit3Formula, newdata3)
predFromMatrix <- predict(fit3Formula, x3Named[seq_len(20L), , drop = FALSE])
expect_identical(predFromFrame, predFromMatrix)

# family = "auto" now detects a 3+-level factor response and fits multinomial:
# the default-family formula fit reproduces the explicit-multinomial formula
# fit bit for bit (the peek is RNG-neutral) and announces the verdict
set.seed(seed3)
expect_message(
  fit3Auto <- bart2(
    y3 ~ x1 + x2 + x3 + x4,
    data = df3,
    keepTrees = TRUE,
    n.trees = n.trees,
    n.chains = 1L,
    n.threads = 1L,
    n.burn = n.burn,
    n.samples = n.samples,
    verbose = FALSE
  ),
  "multinomial"
)
expect_identical(fit3Auto$yhat.train, fit3Formula$yhat.train)

# a factor with an unused level keeps it, as the matrix interface does
# (K = nlevels(y), never dropped)
y3Unused <- factor(y3, levels = c(levels(y3), "unused"))
df3Unused <- data.frame(x3Named, y3Unused = y3Unused)
set.seed(seed3)
fit3Unused <- bart2(
  y3Unused ~ x1 + x2 + x3 + x4,
  data = df3Unused,
  family = "multinomial",
  n.trees = n.trees,
  n.chains = 1L,
  n.threads = 1L,
  n.burn = n.burn,
  n.samples = n.samples,
  verbose = FALSE
)
expect_equal(fit3Unused$levels, c("lo", "mid", "hi", "unused"))
expect_equal(fit3Unused$K, 4L)

# cbind(...) count-matrix LHS: routes to the count path, levels from the
# cbind column names, and reproduces the matrix-interface count fit bit
# for bit at the same seed
x3cNamed <- x3c
colnames(x3cNamed) <- c("x1", "x2", "x3", "x4")
df3c <- data.frame(
  x3cNamed,
  lo = counts3c[, 1L],
  mid = counts3c[, 2L],
  hi = counts3c[, 3L]
)
set.seed(seed3c)
fit3cFormula <- bart2(
  cbind(lo, mid, hi) ~ x1 + x2 + x3 + x4,
  data = df3c,
  family = "multinomial",
  n.trees = n.trees,
  n.chains = 1L,
  n.threads = 1L,
  n.burn = n.burn,
  n.samples = n.samples,
  verbose = FALSE
)
expect_equal(fit3cFormula$levels, c("lo", "mid", "hi"))
probs3cFormula <- fit3cFormula$yhat.train
dimnames(probs3cFormula) <- NULL
probs3cMatrix <- fit3c$yhat.train
dimnames(probs3cMatrix) <- NULL
expect_identical(probs3cFormula, probs3cMatrix)

# character LHS coerces to a factor exactly as the matrix interface does
# (plain factor(), so default alphabetical levels)
df3Char <- data.frame(x3Named, y3Char = as.character(y3))
set.seed(seed3)
fit3Char <- bart2(
  y3Char ~ x1 + x2 + x3 + x4,
  data = df3Char,
  family = "multinomial",
  n.trees = n.trees,
  n.chains = 1L,
  n.threads = 1L,
  n.burn = n.burn,
  n.samples = n.samples,
  verbose = FALSE
)
expect_equal(fit3Char$levels, sort(levels(y3)))

# --- category offset: bart2's own n x K matrix 'offset', threaded to
# the internal creator's own offset argument, never to the host dbarts()
# call (whose flat offset stays refused separately, below). The
# reproduction gate extends to it on both response forms: a public offset
# fit must reproduce internalMultinomialFit/internalMultinomialCountFit's
# own offset-carrying channel bit for bit, and the offset must move the
# answer (non-vacuity). ---
offset3 <- cbind(
  0.6 * (x3[, 1L] - 0.5),
  -0.4 + 0.3 * x3[, 2L],
  0.2 * sin(4 * x3[, 3L])
)
set.seed(seed3)
fit3o <- bart2(
  x3,
  y3,
  family = "multinomial",
  offset = offset3,
  n.trees = n.trees,
  n.chains = 1L,
  n.threads = 1L,
  n.burn = n.burn,
  n.samples = n.samples,
  verbose = FALSE
)
set.seed(seed3)
internal3o <- internalMultinomialFit(
  x3,
  labels3,
  3L,
  n.trees,
  n.burn,
  n.samples,
  offset = offset3
)
probs3o.fit <- fit3o$yhat.train
dimnames(probs3o.fit) <- NULL
expect_identical(probs3o.fit, aperm(internal3o$train, c(3L, 1L, 2L)))
expect_true(all(fit3o$yhat.train >= 0 & fit3o$yhat.train <= 1))
expect_true(all(abs(apply(fit3o$yhat.train, c(1L, 2L), sum) - 1) < 1e-8))
expect_false(isTRUE(all.equal(fit3o$yhat.train, fit3$yhat.train)))

offset3c <- cbind(
  0.5 * (x3c[, 1L] - 0.5),
  -0.3 + 0.2 * x3c[, 2L],
  0.1 * sin(3 * x3c[, 4L])
)
set.seed(seed3c)
fit3co <- bart2(
  x3c,
  counts3c,
  family = "multinomial",
  offset = offset3c,
  n.trees = n.trees,
  n.chains = 1L,
  n.threads = 1L,
  n.burn = n.burn,
  n.samples = n.samples,
  verbose = FALSE
)
set.seed(seed3c)
internal3co <- internalMultinomialCountFit(
  x3c,
  counts3c,
  3L,
  n.trees,
  n.burn,
  n.samples,
  offset = offset3c
)
probs3co.fit <- fit3co$yhat.train
dimnames(probs3co.fit) <- NULL
expect_identical(probs3co.fit, aperm(internal3co$train, c(3L, 1L, 2L)))
expect_false(isTRUE(all.equal(fit3co$yhat.train, fit3c$yhat.train)))

# a wrong-shape matrix refuses naming n and K (validateCategoryOffset,
# reached through bartcoreMultinomialSampler)
expect_error(
  bart2(
    x3,
    y3,
    family = "multinomial",
    offset = offset3[, seq_len(2L), drop = FALSE],
    n.trees = 5L,
    n.burn = 2L,
    n.samples = 2L
  ),
  "150 x 3 matrix"
)
# a non-numeric matrix refuses by name, before it ever reaches validation
expect_error(
  bart2(
    x3,
    y3,
    family = "multinomial",
    offset = matrix("a", n, 3L),
    n.trees = 5L,
    n.burn = 2L,
    n.samples = 2L
  ),
  "numeric"
)

# --- honest refusals ---
# the improved message names the model reason (binding decision 4): an
# integer weight is already expressible as the count matrix's row-wise
# replication, a non-integer one has no exact augmentation sampler
expect_error(
  bart2(
    x2,
    y2,
    family = "multinomial",
    weights = rep(1, n2),
    n.trees = 5L,
    n.burn = 2L,
    n.samples = 2L
  ),
  "already expressible as row-wise count replication"
)
# a flat vector stays refused - it is the softmax's own null direction and
# identically inert (only an n x K matrix carries a meaningful offset)
expect_error(
  bart2(
    x2,
    y2,
    family = "multinomial",
    offset = rep(0, n2),
    n.trees = 5L,
    n.burn = 2L,
    n.samples = 2L
  ),
  "does not support a flat 'offset'"
)
# the internal constructor refuses a host offset for the same reason bart2
# does, rather than silently dropping it: the softmax is invariant to a common
# per-observation shift, so a flat offset is inert and a meaningful one is n x K
samplerOffset <- dbarts(
  x2,
  as.double(labels2),
  offset = rep(0.5, n2),
  control = dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 5L,
    updateState = FALSE
  )
)
expect_error(
  dbarts:::bartcoreMultinomialSampler(samplerOffset, labels2, K = 2L),
  "do not support a flat offset"
)
# offset is train-side only. bart2's own offset.test is caught at the R
# boundary, before it would otherwise fall through to the host dbarts() call
# (does not support 'offset.test'); the underlying host-object flat test
# offset refusal, reached only by building the internal sampler directly,
# names the internal channel instead of an unqualified "do not support"
expect_error(
  bart2(
    x2,
    y2,
    family = "multinomial",
    test = x2.test,
    offset.test = rep(0, 15L),
    n.trees = 5L,
    n.burn = 2L,
    n.samples = 2L
  ),
  "does not support 'offset.test'"
)
samplerTestOffset <- dbarts(
  x2,
  as.double(labels2),
  test = x2.test,
  offset.test = rep(0.5, 15L),
  control = dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 5L,
    updateState = FALSE
  )
)
expect_error(
  dbarts:::bartcoreMultinomialSampler(samplerTestOffset, labels2, K = 2L),
  "category test offset is an nTest x K matrix"
)
# setSigma is refused for the same class of reason: the softmax chain marks
# itself binary-sigma (sigmaScale() is 1 and the redraw is gated off), so a
# value written here would persist and silently rescale every category
# forest's leaf posterior precision
samplerSigma <- dbarts(
  x2,
  as.double(labels2),
  control = dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 5L,
    updateState = FALSE
  )
)
bcSigma <- dbarts:::bartcoreMultinomialSampler(samplerSigma, labels2, K = 2L)
expect_error(
  .Call(dbarts:::C_dbarts_bartcore_setSigma, bcSigma$ptr, 5),
  "response family fixes the residual standard deviation"
)
rm(samplerSigma, bcSigma)

# --- argument leaks (docs/design/multinomial.md): buildMultinomialSampler
# copies only power/base/proposal-probability fields into the K-forest
# engine, so DART (on either formal), split.probs, monotone, and variance
# are refused by name rather than silently dropped. The refusal sits ahead
# of the factor/count-matrix dispatch and of any sampler creation, so no
# n.trees/n.burn/n.samples are needed to reach it; 'dart' is checked on
# both entry shapes to confirm that, the rest on the factor shape only.
multinomialRefuses <- function(x, y, pattern, ...) {
  # nolint next: object_usage_linter. tinytest attaches expect_* at run time.
  expect_error(bart2(x, y, family = "multinomial", ...), pattern)
}
multinomialRefuses(x2, y2, "'dart'", dart = TRUE)
multinomialRefuses(x3c, counts3c, "'dart'", dart = TRUE)
multinomialRefuses(x2, y2, "'dart'", dart = dbartsPriors$dart())
# tree.prior = dart() is the second, formerly-unrefused route to DART; a
# direct call (not the '...'-forwarding helper above, which would itself
# break dart()'s deferred resolution)
expect_error(
  bart2(x2, y2, family = "multinomial", tree.prior = dart()),
  "tree.prior"
)
multinomialRefuses(x2, y2, "split.probs", split.probs = c(0.5, 0.5))
# monotone is half-applied upstream (its proposal.probs rewrite reaches the
# engine, the directions do not), so it is refused rather than silently fit
# as an unintended birth/death-only model
multinomialRefuses(x2, y2, "monotone", monotone = c(1, 0))
# K = 3 (not x2/y2): the host sampler's own auto-resolved family for a
# 2-level response is probit, which already refuses a variance forest for
# an unrelated reason; K = 3 isolates the multinomial-specific refusal
multinomialRefuses(x3, y3, "variance", variance = TRUE)

expect_error(extract(fit3, type = "bart"), "non-identified")
expect_error(fitted(fit3, type = "bart"), "'arg' should be one of")
# fit3 was built WITHOUT keepTrees, so it has no saved trees to replay
expect_error(predict(fit3, x3), "keepTrees")
expect_error(
  bart2(x2, as.double(labels2), family = "multinomial"),
  "factor"
)
expect_error(
  bart2(x2, addNA(y2), family = "multinomial"),
  "NA"
)
expect_error(
  bart2(x2, factor(rep("only", n2)), family = "multinomial"),
  "at least 2 levels"
)
expect_error(
  bart2(dbartsData(x2, as.double(labels2)), family = "multinomial"),
  "dbartsData"
)
expect_error(
  bart2(y2 ~ x2, family = "multinomial"),
  "'data'"
)
# a fit built WITHOUT test data has no test channel to extract
expect_error(extract(fit3, type = "ev", sample = "test"), "no test channel")
