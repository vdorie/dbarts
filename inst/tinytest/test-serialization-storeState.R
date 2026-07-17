source(
  system.file("common", "friedmanData.R", package = "dbarts"),
  local = TRUE
)

# Serializing a keepTrees fit for later predict requires capturing the tree
# state first with storeState() - it is not captured automatically because it
# duplicates the trees on the R side. The round-trip must then reproduce
# pre-serialization predictions exactly, and a fit saved without the capture
# must error on predict naming storeState().

roundTrip <- function(object) {
  tempFile <- tempfile()
  saveRDS(object, file = tempFile)
  on.exit(unlink(tempFile))
  readRDS(tempFile)
}

x <- testData$x
y <- testData$y

# --- classic bart ---------------------------------------------------------
set.seed(99L)
bartFit <- dbarts::bart(
  x,
  y,
  ntree = 5L,
  ndpost = 7L,
  nskip = 0L,
  keeptrees = TRUE,
  verbose = FALSE
)
preds.old <- predict(bartFit, x)
bartFit$fit$storeState()
preds.new <- predict(roundTrip(bartFit), x)
expect_equal(preds.old, preds.new)

# --- bart2 ----------------------------------------------------------------
set.seed(99L)
bart2Fit <- dbarts::bart2(
  x,
  y,
  n.trees = 5L,
  n.samples = 7L,
  n.burn = 0L,
  keepTrees = TRUE,
  n.threads = 1L,
  verbose = FALSE
)
preds.old <- predict(bart2Fit, x)

# without the capture, the reloaded fit refuses predict, naming the mechanism
expect_error(predict(roundTrip(bart2Fit), x), pattern = "storeState")

bart2Fit$fit$storeState()
expect_true(!is.null(bart2Fit$fit$state))
preds.new <- predict(roundTrip(bart2Fit), x)
expect_equal(preds.old, preds.new)

# --- rbart_vi (built-in tau prior: in-core / bartcore path) ---------------
n.g <- 5L
g <- factor(rep_len(seq_len(n.g), length(y)))

set.seed(0L)
rbartFit <- dbarts::rbart_vi(
  y ~ x,
  group.by = g,
  n.samples = 7L,
  n.burn = 0L,
  n.thin = 1L,
  n.chains = 2L,
  n.trees = 25L,
  n.threads = 1L,
  verbose = FALSE
)
preds.old <- predict(rbartFit, x, group.by = g)
invisible(lapply(rbartFit$fit, function(f) f$storeState()))
preds.new <- predict(roundTrip(rbartFit), x, group.by = g)
expect_equal(preds.old, preds.new)
