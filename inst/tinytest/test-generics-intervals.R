# Credible / prediction intervals through fitted() and predict(): a scalar
# ci.level opts into a per-observation est + ci.lower + ci.upper matrix, with
# the interval kind following type - "ev" is a credible interval for the mean
# (a probability for binary), "ppd" a prediction interval that adds residual
# noise (so it is wider), "bart" the latent scale. Default (NULL) is unchanged.

source(
  system.file("common", "friedmanData.R", package = "dbarts"),
  local = TRUE
)

x <- testData$x
y <- testData$y

fit <- bart2(
  x,
  y,
  n.samples = 100L,
  n.burn = 30L,
  n.trees = 25L,
  n.chains = 2L,
  n.threads = 1L,
  keepTrees = TRUE,
  verbose = FALSE
)

# default fitted is unchanged (a vector), ci.level returns the 3-column matrix
expect_null(dim(fitted(fit)))
cred <- fitted(fit, ci.level = 0.95)
expect_true(is.matrix(cred))
expect_equal(colnames(cred), c("est", "ci.lower", "ci.upper"))
expect_equal(nrow(cred), length(y))

# the est column is exactly the posterior mean that plain fitted() returns
expect_equal(unname(cred[, "est"]), unname(fitted(fit)))

# est lies within the band, and the band is ordered
expect_true(all(
  cred[, "ci.lower"] <= cred[, "est"] &
    cred[, "est"] <= cred[, "ci.upper"]
))

# a prediction interval (ppd) carries residual noise, so it is wider than the
# credible interval (ev) for the same level
pred <- fitted(fit, type = "ppd", ci.level = 0.95)
expect_true(
  mean(pred[, "ci.upper"] - pred[, "ci.lower"]) >
    mean(cred[, "ci.upper"] - cred[, "ci.lower"])
)

# a tighter level gives a narrower band
narrow <- fitted(fit, ci.level = 0.5)
expect_true(
  mean(narrow[, "ci.upper"] - narrow[, "ci.lower"]) <
    mean(cred[, "ci.upper"] - cred[, "ci.lower"])
)

# predict() takes ci.level too, on new data
pci <- predict(fit, x[1:5, ], ci.level = 0.9)
expect_true(is.matrix(pci) && nrow(pci) == 5L)
expect_equal(colnames(pci), c("est", "ci.lower", "ci.upper"))

# ci.level is validated
expect_error(fitted(fit, ci.level = 1.2), pattern = "must be a single number")
expect_error(
  fitted(fit, ci.level = c(0.9, 0.95)),
  pattern = "must be a single number"
)

# fitted()'s positional slot 3 is now ci.level (it used to be sample): the
# positional and named forms agree, and the old slot-3 usage - a train/test
# string in that position - is refused by a message naming 'sample' and where
# it went, instead of silently binding to it
expect_identical(
  fitted(fit, "ev", 0.9),
  fitted(fit, type = "ev", ci.level = 0.9)
)
expect_error(
  fitted(fit, "ev", "train"),
  pattern = "'sample' is fitted's fourth argument"
)

rm(fit, cred, pred, narrow, pci, x, y)
rm(testData)


source(system.file("common", "probitData.R", package = "dbarts"), local = TRUE)

X <- testData$X
Z <- testData$Z

fit <- bart2(
  X,
  Z,
  n.samples = 100L,
  n.burn = 30L,
  n.trees = 25L,
  n.chains = 1L,
  n.threads = 1L,
  keepTrees = TRUE,
  verbose = FALSE
)

# a binary ev interval is on the probability scale, entirely within [0, 1]
pci <- fitted(fit, ci.level = 0.95)
expect_true(all(pci >= 0 & pci <= 1))

# the latent-scale (bart) interval is not confined to [0, 1]
lci <- fitted(fit, type = "bart", ci.level = 0.95)
expect_true(any(lci[, "ci.lower"] < 0) || any(lci[, "ci.upper"] > 1))

rm(fit, pci, lci, X, Z)
rm(testData)


source(
  system.file("common", "friedmanData.R", package = "dbarts"),
  local = TRUE
)

g <- rep_len(seq_len(5L), length(testData$y))
rfit <- rbart_vi(
  testData$y ~ testData$x,
  group.by = g,
  n.samples = 40L,
  n.burn = 20L,
  n.thin = 1L,
  n.chains = 1L,
  n.trees = 20L,
  n.threads = 1L,
  keepTrees = TRUE,
  verbose = FALSE
)

# rbart fitted/predict carry ci.level; est matches the mean-only fitted
rci <- fitted(rfit, ci.level = 0.9)
expect_equal(colnames(rci), c("est", "ci.lower", "ci.upper"))
expect_equal(unname(rci[, "est"]), unname(fitted(rfit)), tolerance = 1.0e-6)
expect_true(is.matrix(predict(
  rfit,
  testData$x[1:4, ],
  group.by = g[1:4],
  ci.level = 0.9
)))

# as above: fitted.rbart's positional slot 3 is now ci.level too
expect_identical(
  fitted(rfit, "ev", 0.9),
  fitted(rfit, type = "ev", ci.level = 0.9)
)

rm(rfit, rci, g)
rm(testData)

# the remaining three classes' fitted() already had ci.level third before
# this slice, so the identity below asserts nothing this slice could break -
# a plain regression guard, kept for symmetry with the three above
n <- 40L
xSmall <- matrix(rnorm(n * 2L), n, 2L)

fitM <- bart2(
  xSmall,
  factor(sample(letters[1:3], n, replace = TRUE)),
  family = "multinomial",
  n.samples = 10L,
  n.burn = 5L,
  n.trees = 5L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)
expect_identical(
  fitted(fitM, "ev", 0.9),
  fitted(fitM, type = "ev", ci.level = 0.9)
)
rm(fitM)

fitO <- bart2(
  xSmall,
  ordered(
    sample(c("lo", "mid", "hi"), n, replace = TRUE),
    levels = c("lo", "mid", "hi")
  ),
  family = "ordinal",
  n.samples = 10L,
  n.burn = 5L,
  n.trees = 5L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)
expect_identical(
  fitted(fitO, "ev", 0.9),
  fitted(fitO, type = "ev", ci.level = 0.9)
)
rm(fitO)

fitN <- bart2(
  xSmall,
  rpois(n, 3),
  family = "nbinom",
  n.samples = 10L,
  n.burn = 5L,
  n.trees = 5L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)
expect_identical(
  fitted(fitN, "ev", 0.9),
  fitted(fitN, type = "ev", ci.level = 0.9)
)
rm(fitN)

# hurdle DOES discriminate: fitted.bartHurdle lost its 'sample' formal
# outright (section 8), so slot 3 is ci.level with no other candidate
fitH <- bart2(
  xSmall,
  ifelse(runif(n) < 0.5, 0, rlnorm(n)),
  family = "hurdle.lognormal",
  n.samples = 10L,
  n.burn = 5L,
  n.trees = 5L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)
expect_identical(
  fitted(fitH, "ev", 0.9),
  fitted(fitH, type = "ev", ci.level = 0.9)
)
rm(fitH, n, xSmall)
