source(
  system.file("common", "rbartGroupData.R", package = "dbarts"),
  local = TRUE
)
testData$g <- g
testData$b <- b
rm(b, sigma.b, g, n.g)

# test that extract works at baseline
x <- testData$x
y <- testData$y
g <- factor(testData$g)

## combineChains = FALSE pinned deliberately: the assertions below index
## rbartFit$ranef/$yhat.train as raw uncombined (n.chains x n.samples x ...)
## arrays, and are contrasted against rbartFit.2's explicit TRUE below - now
## that rbart_vi's own stored-object default is combined
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
  verbose = FALSE,
  combineChains = FALSE
)
expect_equal(
  dbarts::extract(rbartFit, type = "bart", combineChains = FALSE),
  rbartFit$yhat.train
)
expect_equal(
  dbarts::extract(rbartFit, type = "ranef", combineChains = FALSE),
  rbartFit$ranef
)
expect_equal(
  dbarts::extract(rbartFit, type = "ev", combineChains = FALSE),
  rbartFit$yhat.train + unname(rbartFit$ranef[,, as.character(g)])
)

ppd <- dbarts::extract(rbartFit, type = "ppd", combineChains = FALSE)
sigma.hat <- apply(
  ppd - dbarts::extract(rbartFit, type = "ev", combineChains = FALSE),
  c(1L, 2L),
  sd
)
# silly test with 7 samples but we need something
expect_true(
  cor(as.vector(sigma.hat), as.vector(rbartFit$sigma)) >= 0.85
)

set.seed(0L)
rbartFit.2 <- dbarts::rbart_vi(
  y ~ x,
  group.by = g,
  n.samples = 7L,
  n.burn = 0L,
  n.thin = 1L,
  n.chains = 2L,
  n.trees = 25L,
  n.threads = 1L,
  combineChains = TRUE,
  verbose = FALSE
)
expect_equal(
  dbarts::extract(rbartFit, type = "ev", combineChains = TRUE),
  dbarts::extract(rbartFit.2, type = "ev", combineChains = TRUE)
)
expect_equal(
  dbarts::extract(rbartFit, type = "ev", combineChains = FALSE),
  dbarts::extract(rbartFit.2, type = "ev", combineChains = FALSE)
)

rm(sigma.hat, ppd, rbartFit, g, y, x)


# test that fitted works correctly
x <- testData$x
y <- testData$y
g <- factor(testData$g)


## combineChains = FALSE pinned deliberately: the assertions below index
## rbartFit$ranef/$yhat.train/$yhat.test as raw uncombined 3-d arrays
rbartFit <- dbarts::rbart_vi(
  y ~ x,
  group.by = g,
  test = x,
  group.by.test = g,
  offset.test = 5,
  n.samples = 7L,
  n.burn = 0L,
  n.thin = 1L,
  n.chains = 2L,
  n.trees = 25L,
  n.threads = 1L,
  verbose = FALSE,
  combineChains = FALSE
)
expect_equal(
  as.vector(rbartFit$yhat.train),
  as.vector(rbartFit$yhat.test) - 5
)
expect_equal(
  apply(
    rbartFit$yhat.train + unname(rbartFit$ranef[,, as.character(g)]),
    3L,
    mean
  ),
  fitted(rbartFit)
)
expect_equal(
  apply(
    rbartFit$yhat.test + unname(rbartFit$ranef[,, as.character(g)]),
    3L,
    mean
  ),
  fitted(rbartFit, sample = "test")
)
expect_equal(fitted(rbartFit, type = "ranef"), rbartFit$ranef.mean)

rm(rbartFit, g, y, x)


# test that predict matches fitted
x <- testData$x
y <- testData$y
g <- factor(testData$g)

set.seed(0L)
rbartFit.0 <- dbarts::rbart_vi(
  y ~ x,
  group.by = g,
  n.samples = 14L,
  n.burn = 0L,
  n.thin = 1L,
  n.chains = 1L,
  n.trees = 25L,
  n.threads = 1L,
  keepTrees = TRUE,
  verbose = FALSE
)
expect_equal(
  fitted(rbartFit.0),
  apply(predict(rbartFit.0, x, group.by = g), 2L, mean)
)

set.seed(0L)
rbartFit.0 <- dbarts::rbart_vi(
  y ~ x,
  group.by = g,
  n.samples = 7L,
  n.burn = 0L,
  n.thin = 1L,
  n.chains = 2L,
  n.trees = 25L,
  n.threads = 1L,
  keepTrees = TRUE,
  combineChains = FALSE,
  verbose = FALSE
)
expect_equal(
  fitted(rbartFit.0),
  apply(predict(rbartFit.0, x, group.by = g, combineChains = FALSE), 3L, mean)
)

set.seed(0L)
rbartFit.1 <- dbarts::rbart_vi(
  y ~ x,
  group.by = g,
  n.samples = 7L,
  n.burn = 0L,
  n.thin = 1L,
  n.chains = 2L,
  n.trees = 25L,
  n.threads = 1L,
  keepTrees = TRUE,
  combineChains = TRUE,
  verbose = FALSE
)
expect_equal(
  fitted(rbartFit.1),
  apply(predict(rbartFit.1, x, group.by = g), 2L, mean)
)
expect_equal(
  predict(rbartFit.0, x, group.by = g),
  predict(rbartFit.1, x, group.by = g)
)
expect_equal(
  predict(rbartFit.0, x, group.by = g, combineChains = FALSE),
  predict(rbartFit.1, x, group.by = g, combineChains = FALSE)
)

# an rbart fit whose ranef dimnames do not name every group once crashed the
# session: fitted matched each label to NA and handed the NAs to C, which
# indexed the ranef array on NA_INTEGER
strippedFit <- rbartFit.1
dimnames(strippedFit$ranef) <- NULL
expect_error(
  fitted(strippedFit),
  "'group.by' must name groups present in the 'ranef' dimnames"
)
expect_error(
  residuals(strippedFit),
  "'group.by' must name groups present in the 'ranef' dimnames"
)
rm(strippedFit)

# the compiled entry point backstops a direct .Call that skips the guard above
n.obs <- length(rbartFit.1$group.by)
expect_error(
  .Call(
    dbarts:::C_rbart_fitted,
    rbartFit.1$yhat.train,
    rbartFit.1$ranef,
    rep_len(NA_integer_, n.obs),
    FALSE
  ),
  "group index must be a non-NA index into the 'ranef' group dimension"
)
expect_error(
  .Call(
    dbarts:::C_rbart_fitted,
    rbartFit.1$yhat.train,
    rbartFit.1$ranef,
    rep_len(ncol(rbartFit.1$ranef) + 1L, n.obs),
    FALSE
  ),
  "group index must be a non-NA index into the 'ranef' group dimension"
)
expect_error(
  .Call(
    dbarts:::C_rbart_fitted,
    rbartFit.1$yhat.train,
    rbartFit.1$ranef,
    rep_len(0L, n.obs),
    FALSE
  ),
  "group index must be a non-NA index into the 'ranef' group dimension"
)
expect_error(
  .Call(
    dbarts:::C_rbart_fitted,
    rbartFit.1$yhat.train,
    rbartFit.1$ranef,
    rep_len(1L, n.obs - 1L),
    FALSE
  ),
  "group index length must match the number of observations"
)
rm(n.obs)

rm(rbartFit.1, rbartFit.0, g, y, x)


rm(testData)
