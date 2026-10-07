# fitted(type = "forest") is the posterior mean of extract(type = "forest"):
# a rows x forests matrix on a fit with several forests

set.seed(4101)
n <- 60L
d <- data.frame(x1 = rnorm(n), x2 = rnorm(n), z = rbinom(n, 1L, 0.5))
d$y <- d$x1 + d$z * (1 + d$x2) + rnorm(n, 0, 0.3)

fitForest <- function(data, n.chains = 1L, ...) {
  bart(
    y ~ x1 + x2 + forest(x1 + x2, basis = ~z),
    data,
    n.chains = n.chains,
    n.threads = 1L,
    n.trees = 8L,
    n.samples = 6L,
    n.burn = 4L,
    verbose = FALSE,
    keepTrees = TRUE,
    ...
  )
}

fit <- fitForest(d)
draws <- extract(fit, type = "forest")
expected <- dbarts:::channelMeans(draws, 2L)

mean.forest <- fitted(fit, type = "forest")
expect_equal(mean.forest, expected)
expect_equal(dim(mean.forest), dim(draws)[2:3])
expect_equal(dimnames(mean.forest), dimnames(draws)[2:3])
expect_equal(mean.forest, apply(draws, 2:3, mean))

# chains and draws average together
fit2 <- fitForest(d, n.chains = 2L)
draws2 <- extract(fit2, type = "forest")
expect_equal(dim(draws2)[1L], 12L)
expect_equal(
  fitted(fit2, type = "forest"),
  dbarts:::channelMeans(draws2, 2L)
)

# ci.level: est is the mean, and brackets it
band <- fitted(fit, type = "forest", ci.level = 0.9)
expect_equal(dim(band), c(dim(mean.forest), 3L))
expect_equal(
  dimnames(band),
  c(dimnames(mean.forest), list(c("est", "ci.lower", "ci.upper")))
)
expect_equal(band[,, "est"], mean.forest)
expect_true(all(band[,, "ci.lower"] <= band[,, "est"]))
expect_true(all(band[,, "est"] <= band[,, "ci.upper"]))

# rows the fit's na.action dropped come back as NA
dNA <- d
dNA$x1[c(3L, 10L)] <- NA
fitNA <- fitForest(dNA, na.action = na.exclude)
meanNA <- fitted(fitNA, type = "forest")
expect_equal(dim(meanNA), c(n, 2L))
expect_true(all(is.na(meanNA[c(3L, 10L), ])))
expect_false(anyNA(meanNA[-c(3L, 10L), ]))
bandNA <- fitted(fitNA, type = "forest", ci.level = 0.9)
expect_equal(dim(bandNA), c(n, 2L, 3L))
expect_true(all(is.na(bandNA[c(3L, 10L), , ])))
expect_false(anyNA(bandNA[-c(3L, 10L), , ]))
expect_equal(
  bandNA[-c(3L, 10L), , "est"],
  unname(meanNA[-c(3L, 10L), ]),
  check.attributes = FALSE
)

# what extract refuses, fitted refuses with extract's own message
single <- bart(
  y ~ x1 + x2,
  d,
  n.threads = 1L,
  n.trees = 8L,
  n.samples = 6L,
  n.burn = 4L,
  verbose = FALSE
)
singleMessage <- tryCatch(
  extract(single, type = "forest"),
  error = conditionMessage
)
expect_true(is.character(singleMessage))
expect_error(
  fitted(single, type = "forest"),
  singleMessage,
  fixed = TRUE
)
testMessage <- tryCatch(
  extract(fit, type = "forest", sample = "test"),
  error = conditionMessage
)
expect_true(is.character(testMessage))
expect_error(
  fitted(fit, type = "forest", sample = "test"),
  testMessage,
  fixed = TRUE
)

# the default is still one vector, and residuals does not take the type
expect_true(is.null(dim(fitted(fit))))
expect_equal(length(fitted(fit)), n)
expect_equal(eval(formals(dbarts:::fitted.bart)$type)[1L], "ev")
expect_error(
  residuals(fit, type = "forest"),
  "type = \"forest\" is not used by residuals",
  fixed = TRUE
)

# a multinomial fit is refused as extract refuses it
set.seed(4102)
xm <- matrix(runif(60L * 2L), 60L, 2L)
ym <- factor(sample(c("a", "b", "c"), 60L, replace = TRUE))
fitM <- bart(
  xm,
  ym,
  family = "multinomial",
  n.trees = 6L,
  n.threads = 1L,
  n.burn = 4L,
  n.samples = 6L,
  verbose = FALSE
)
multinomialMessage <- tryCatch(
  extract(fitM, type = "forest"),
  error = conditionMessage
)
expect_true(is.character(multinomialMessage))
expect_error(
  fitted(fitM, type = "forest"),
  multinomialMessage,
  fixed = TRUE
)
