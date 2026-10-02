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

fit <- bart(
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

# the band is the equal-tailed one: the (1 - ci.level) / 2 and
# (1 + ci.level) / 2 quantiles of the draws, per observation
draws <- predict(fit, x[1:5, ])
expect_equal(
  unname(pci[, c("ci.lower", "ci.upper")]),
  unname(t(apply(draws, 2L, quantile, probs = c(0.05, 0.95))))
)
expect_equal(unname(pci[, "est"]), unname(colMeans(draws)))
rm(draws)

# ci.level is validated
expect_error(fitted(fit, ci.level = 1.2), pattern = "must be a single number")
expect_error(
  fitted(fit, ci.level = c(0.9, 0.95)),
  pattern = "must be a single number"
)

# fitted()'s positional slot 3 is 'sample' (0.9-x's order, restored); ci.level
# is the fourth slot, and the positional and named forms agree at both
expect_identical(
  fitted(fit, "ev", "train"),
  fitted(fit, type = "ev", sample = "train")
)
expect_identical(
  fitted(fit, "ev", "train", 0.9),
  fitted(fit, type = "ev", sample = "train", ci.level = 0.9)
)

rm(fit, cred, pred, narrow, pci, x, y)
rm(testData)


source(system.file("common", "probitData.R", package = "dbarts"), local = TRUE)

X <- testData$X
Z <- testData$Z

fit <- bart(
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

# multinomial, ordinal and negbin gain 'sample' as their third argument too,
# matching every other family with a test channel; ci.level moves to the
# fourth slot on each
n <- 40L
xSmall <- matrix(rnorm(n * 2L), n, 2L)

fitM <- bart(
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
  fitted(fitM, "ev", "train"),
  fitted(fitM, type = "ev", sample = "train")
)
expect_identical(
  fitted(fitM, "ev", "train", 0.9),
  fitted(fitM, type = "ev", ci.level = 0.9)
)
rm(fitM)

fitO <- bart(
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
  fitted(fitO, "ev", "train"),
  fitted(fitO, type = "ev", sample = "train")
)
expect_identical(
  fitted(fitO, "ev", "train", 0.9),
  fitted(fitO, type = "ev", ci.level = 0.9)
)
rm(fitO)

fitN <- bart(
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
  fitted(fitN, "ev", "train"),
  fitted(fitN, type = "ev", sample = "train")
)
expect_identical(
  fitted(fitN, "ev", "train", 0.9),
  fitted(fitN, type = "ev", ci.level = 0.9)
)
rm(fitN)

# hurdle DOES discriminate: fitted.bartHurdle lost its 'sample' formal
# outright (section 8), so slot 3 is ci.level with no other candidate
fitH <- bart(
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
