source(
  system.file("common", "friedmanData.R", package = "dbarts"),
  local = TRUE
)

# A bart() fit that keeps its trees stores the sampler's state at the end of the
# fit, so a fit saved and read back predicts as the original did; a sampler
# run further is stored again with storeState() before saving. A fit that keeps
# no trees has no sampler to store. A bartBT-style fit keeps its tree state as
# 0.9-x did.

roundTripStoreState <- function(object) {
  tempFile <- tempfile()
  saveRDS(object, file = tempFile)
  on.exit(unlink(tempFile))
  readRDS(tempFile)
}

x <- testData$x
y <- testData$y

# --- classic bart ---------------------------------------------------------
set.seed(99L)
bartFit <- dbarts::bartBT(
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
preds.new <- predict(roundTripStoreState(bartFit), x)
expect_equal(preds.old, preds.new)

# --- bart ----------------------------------------------------------------
set.seed(99L)
bart2Fit <- dbarts::bart(
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

# no manual step: the fit stored its state at the end of the fit
# (reading the state field would itself store it, so it is read after)
expect_equal(predict(roundTripStoreState(bart2Fit), x), preds.old)
expect_true(!is.null(bart2Fit$fit$state))
bart2Fit$fit$storeState()
expect_equal(predict(roundTripStoreState(bart2Fit), x), preds.old)

# a default hazard fit keeps its trees, so it reloads whole too
set.seed(3L)
nh <- 40L
h <- data.frame(
  time = sample(1:3, nh, TRUE),
  status = rbinom(nh, 1L, 0.7),
  x1 = runif(nh)
)
hazardFit <- bart(
  survival::Surv(time, status) ~ x1,
  data = h,
  family = "hazard",
  n.trees = 5L,
  n.samples = 6L,
  n.burn = 2L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)
expect_equal(
  survivalProbabilities(
    roundTripStoreState(hazardFit),
    1:3,
    newdata = h[1:5, ]
  ),
  survivalProbabilities(hazardFit, 1:3, newdata = h[1:5, ])
)

# a fit that keeps no trees keeps no sampler and is unchanged
noTrees <- dbarts::bart(
  x,
  y,
  n.trees = 5L,
  n.samples = 7L,
  n.burn = 0L,
  keepTrees = FALSE,
  n.threads = 1L,
  verbose = FALSE
)
expect_null(noTrees$fit)
expect_equal(roundTripStoreState(noTrees)$yhat.train, noTrees$yhat.train)
