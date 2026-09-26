# fitted's third positional argument is 'sample' (0.9-x's order, restored) on
# every class with a test channel: bart, and - newly, since they gain a
# 'sample' formal of their own - bartMultinomial, bartOrdinal and bartNegbin.
# fitted(fit, "ev", "test") must reach the same test-channel draws
# extract(fit, sample = "test") does, not silently stay on the training rows.
# n.chains = 1 throughout so extract's own chain-combination is a no-op and
# the independent computation below matches fitted's internals exactly.

set.seed(33)
n <- 60L
nTest <- 15L
x <- matrix(runif(n), n, 1L)
xTest <- matrix(runif(nTest), nTest, 1L)

# --- bart (gaussian) ---
y <- x[, 1L] + rnorm(n, 0, 0.2)
fit <- bart(
  x,
  y,
  test = xTest,
  n.samples = 20L,
  n.burn = 10L,
  n.trees = 10L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)
expect_identical(
  fitted(fit, "ev", "test"),
  dbarts:::channelMeans(extract(fit, type = "ev", sample = "test"))
)

# --- multinomial ---
yCat <- factor(sample(letters[1:3], n, replace = TRUE))
fitM <- bart(
  x,
  yCat,
  test = xTest,
  family = "multinomial",
  n.samples = 20L,
  n.burn = 10L,
  n.trees = 10L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)
expect_identical(
  fitted(fitM, "ev", "test"),
  dbarts:::meanCategoryProbabilities(
    extract(fitM, type = "ev", sample = "test"),
    fitM$levels
  )
)

# --- ordinal ---
yOrd <- ordered(
  sample(c("lo", "mid", "hi"), n, replace = TRUE),
  levels = c("lo", "mid", "hi")
)
fitO <- bart(
  x,
  yOrd,
  test = xTest,
  family = "ordinal",
  n.samples = 20L,
  n.burn = 10L,
  n.trees = 10L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)
expect_identical(
  fitted(fitO, "ev", "test"),
  dbarts:::meanCategoryProbabilities(
    extract(fitO, type = "ev", sample = "test"),
    fitO$levels
  )
)

# --- nbinom ---
yCount <- rpois(n, 3)
fitN <- bart(
  x,
  yCount,
  test = xTest,
  family = "nbinom",
  n.samples = 20L,
  n.burn = 10L,
  n.trees = 10L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)
expect_identical(
  fitted(fitN, "ev", "test"),
  dbarts:::channelMeans(extract(fitN, type = "ev", sample = "test"))
)

rm(fit, fitM, fitO, fitN, x, xTest, y, yCat, yOrd, yCount, n, nTest)
