# The public bart2(family = "multinomial") surface (docs/plans/multinomial.md
# C7). The REPRODUCTION GATE below is the abort condition the brief calls
# for: bart2's fit path must reproduce, bit for bit, the internal
# bartcoreMultinomialSampler/bartcoreRun pattern benchmarks/R/
# multinomial-equivalence.R exercises, on the same data and seed. Everything
# else is level-threading, shape, and refusal coverage.

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
  test = NULL
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
  bc <- dbarts:::bartcoreMultinomialSampler(sampler, labels, K = K)
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

# --- test data at creation (C1): K = 3 reproduction gate, shape, levels ---
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

# --- test data at creation (C1): K = 2 ---
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

# --- honest refusals ---
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
  "weights"
)
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
  "offset"
)
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
  bart2(
    y2 ~ x2,
    data = data.frame(y2 = y2, x2 = x2),
    family = "multinomial"
  ),
  "matrix interface"
)
# a fit built WITHOUT test data has no test channel to extract
expect_error(extract(fit3, type = "ev", sample = "test"), "no test channel")
