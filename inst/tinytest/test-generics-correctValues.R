source(
  system.file("common", "friedmanData.R", package = "dbarts"),
  local = TRUE
)

# test that predict gives same result as x_train with linear data
bartFit <- dbarts::bartBT(
  testData$x,
  testData$y,
  ndpost = 20,
  nskip = 5,
  ntree = 5L,
  verbose = FALSE,
  keeptrees = TRUE
)
predictions <- predict(bartFit, testData$x, n.threads = 1L)
expect_equal(predictions, bartFit$yhat.train)

bartFit <- dbarts::bartBT(
  testData$x,
  testData$y,
  ndpost = 20,
  nskip = 5,
  ntree = 5L,
  nchain = 4L,
  nthread = 1L,
  verbose = FALSE,
  keeptrees = TRUE
)
predictions <- predict(bartFit, testData$x, n.threads = 1L)
expect_equal(predictions, bartFit$yhat.train)

rm(predictions, bartFit)

# test that extract and fitted give correct results
n.chains <- 4L
n.samples <- 20L
bartFit <- dbarts::bartBT(
  testData$x,
  testData$y,
  testData$x[1L:10L, ],
  ndpost = n.samples,
  nskip = 5L,
  ntree = 5L,
  nchain = n.chains,
  verbose = FALSE
)

expect_equal(extract(bartFit), bartFit$yhat.train)
expect_equal(fitted(bartFit), bartFit$yhat.train.mean)

expect_equal(extract(bartFit, sample = "test"), bartFit$yhat.test)
expect_equal(fitted(bartFit, sample = "test"), bartFit$yhat.test.mean)

extracted <- extract(bartFit, combineChains = FALSE)
for (i in seq_len(n.chains)) {
  expect_equal(
    extracted[i, , ],
    bartFit$yhat.train[seq_len(n.samples) + (i - 1L) * n.samples, ]
  )
}

bartFit <- dbarts::bartBT(
  testData$x,
  testData$y,
  testData$x[1L:10L, ],
  ndpost = n.samples,
  nskip = 5,
  ntree = 5L,
  nchain = n.chains,
  verbose = FALSE,
  combinechains = FALSE
)
extracted <- extract(bartFit)
for (i in seq_len(n.chains)) {
  expect_equal(
    extracted[seq_len(n.samples) + (i - 1L) * n.samples, ],
    bartFit$yhat.train[i, , ]
  )
}

rm(i, extracted, bartFit, n.samples, n.chains)

rm(testData)


source(system.file("common", "probitData.R", package = "dbarts"), local = TRUE)

# test that predict gives same result as x_train with binary data
bartFit <- dbarts::bartBT(
  y.train = testData$Z,
  x.train = testData$X,
  ndpost = 20L,
  nskip = 5L,
  ntree = 5L,
  k = 4.5,
  verbose = FALSE,
  keeptrees = TRUE
)
predictions <- predict(bartFit, testData$X, type = "bart", n.threads = 1L)
expect_equal(predictions, bartFit$yhat.train)

bartFit <- dbarts::bartBT(
  y.train = testData$Z,
  x.train = testData$X,
  ndpost = 20L,
  nskip = 5L,
  ntree = 5L,
  k = 4.5,
  nchain = 4L,
  nthread = 1L,
  verbose = FALSE,
  keeptrees = TRUE
)
predictions <- predict(bartFit, testData$X, type = "bart", n.threads = 1L)
expect_equal(predictions, bartFit$yhat.train)

# residuals reach fitted's type: the default is the response scale (y - phat),
# and type = "bart" gives link-scale residuals, distinct for a binary fit
resid.ev <- residuals(bartFit)
resid.link <- residuals(bartFit, type = "bart")
expect_equal(resid.ev, testData$Z - fitted(bartFit, type = "ev"))
expect_equal(resid.link, testData$Z - fitted(bartFit, type = "bart"))
expect_false(isTRUE(all.equal(resid.ev, resid.link)))

rm(predictions, bartFit, resid.ev, resid.link)

# sigma and nobs answer as for an lm or glm fit: the posterior mean of the
# residual sd, 1 with no residual scale, and the rows in the likelihood
set.seed(31L)
accessorData <- data.frame(a = rnorm(40L), b = rnorm(40L))
accessorData$y <- accessorData$a + rnorm(40L)
accessorData$z <- as.integer(accessorData$y > 0)
accessorData$w <- c(0, 0, rep(1, 38L))
accessorData$a[3L] <- NA
# the call is rebuilt so that weights = w is read in the data
accessorFit <- function(...) {
  fitCall <- match.call()
  fitCall[[1L]] <- quote(dbarts::bart)
  fitCall[c(
    "n.trees",
    "n.samples",
    "n.burn",
    "n.chains",
    "n.threads",
    "verbose"
  )] <- list(5L, 5L, 5L, 2L, 1L, FALSE)
  eval(fitCall, parent.frame())
}
gaussianFit <- accessorFit(
  y ~ a + b,
  accessorData,
  weights = w,
  na.action = na.omit
)
expect_identical(sigma(gaussianFit), mean(gaussianFit$sigma))
expect_identical(nobs(gaussianFit), 37L)
probitFit <- accessorFit(z ~ a + b, accessorData, subset = 4:40, weights = w)
expect_identical(sigma(probitFit), 1)
expect_identical(nobs(probitFit), 37L)
matrixFit <- accessorFit(
  as.matrix(accessorData[4:40, 1:2]),
  accessorData$y[4:40]
)
expect_identical(nobs(matrixFit), 37L)
expect_identical(sigma(matrixFit), mean(matrixFit$sigma))
hurdleFit <- accessorFit(
  ifelse(y > 0, y, 0) ~ a + b,
  accessorData,
  family = "hurdle.lognormal",
  weights = w,
  na.action = na.omit
)
expect_identical(sigma(hurdleFit), mean(extract(hurdleFit, type = "sigma")))
expect_identical(nobs(hurdleFit), 37L)
ordinalFit <- accessorFit(
  ordered(cut(y, 3L)) ~ a + b,
  accessorData,
  weights = w,
  na.action = na.omit
)
expect_identical(sigma(ordinalFit), 1)
expect_identical(nobs(ordinalFit), 37L)
expect_error(
  sigma(accessorFit(y ~ b, accessorData, variance = TRUE)),
  "a fit with a variance forest has no single residual sd; see extract(fit, type = \"sigma\")",
  fixed = TRUE
)
rm(
  accessorData,
  accessorFit,
  gaussianFit,
  probitFit,
  matrixFit,
  hurdleFit,
  ordinalFit
)

rm(testData)
