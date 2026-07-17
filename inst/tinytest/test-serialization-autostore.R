source(
  system.file("common", "friedmanData.R", package = "dbarts"),
  local = TRUE
)

# keepTrees fits store the sampler state automatically at fit time, so a
# save/load round-trip followed by predict works with no manual ritual and
# reproduces the pre-serialization predictions exactly.

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
preds.new <- predict(roundTrip(bart2Fit), x)
expect_equal(preds.old, preds.new)

# the automatically stored state is small relative to the posterior blocks it
# rides beside; on a save/load round-trip the state is what makes predict work
expect_true(!is.null(bart2Fit$fit$state))

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
preds.new <- predict(roundTrip(rbartFit), x, group.by = g)
expect_equal(preds.old, preds.new)

# --- rbart_vi (custom tau prior: per-chain R-loop path) -------------------
customCauchy <- function(x, rel.scale) dcauchy(x, 0, rel.scale * 2.5, TRUE)
set.seed(0L)
rbartFit.custom <- dbarts::rbart_vi(
  y ~ x,
  group.by = g,
  prior = customCauchy,
  n.samples = 7L,
  n.burn = 0L,
  n.thin = 1L,
  n.chains = 2L,
  n.trees = 25L,
  n.threads = 1L,
  verbose = FALSE
)
preds.old <- predict(rbartFit.custom, x, group.by = g)
preds.new <- predict(roundTrip(rbartFit.custom), x, group.by = g)
expect_equal(preds.old, preds.new)
