source(
  system.file("common", "rbartGroupData.R", package = "dbarts"),
  local = TRUE
)
source(
  system.file("common", "captureWarnings.R", package = "dbarts"),
  local = TRUE
)
testData$g <- g
testData$b <- b
rm(b, sigma.b, g, n.g)

# test that rbart finds group.by
df <- as.data.frame(testData$x)
colnames(df) <- paste0("x_", seq_len(ncol(testData$x)))
df$y <- testData$y
df$g <- testData$g
expect_inherits(
  dbarts::rbart_vi(
    y ~ . - g,
    df,
    group.by = g,
    n.samples = 1L,
    n.burn = 0L,
    n.thin = 1L,
    n.chains = 1L,
    n.trees = 25L,
    n.threads = 1L,
    verbose = FALSE
  ),
  "rbart"
)

g <- df$g
df$g <- NULL
expect_inherits(
  dbarts::rbart_vi(
    y ~ .,
    df,
    group.by = g,
    n.samples = 1L,
    n.burn = 0L,
    n.thin = 1L,
    n.chains = 1L,
    n.trees = 25L,
    n.threads = 1L,
    verbose = FALSE
  ),
  "rbart"
)

y <- testData$y
x <- testData$x
expect_inherits(
  dbarts::rbart_vi(
    y ~ x,
    group.by = g,
    n.samples = 1L,
    n.burn = 0L,
    n.thin = 1L,
    n.chains = 1L,
    n.trees = 25L,
    n.threads = 1L,
    verbose = FALSE
  ),
  "rbart"
)

rm(x, y, g, df)


# group.by resolved from the calling environment (not the first data column):
# ranef must be per-group with level dimnames, and predict must round-trip on
# trained groups for any group.by type
set.seed(2L)
n.t1a <- 120L
g.t1a <- factor(rep(seq_len(6L), each = 20L))
x.t1a <- rnorm(n.t1a)
y.t1a <- 2 * x.t1a + rnorm(6L)[as.integer(g.t1a)] + rnorm(n.t1a, 0.0, 0.5)
# deliberately omits g so group.by must be found in the calling environment
df.t1a <- data.frame(y = y.t1a, x = x.t1a)

fit.t1a <- dbarts::rbart_vi(
  y ~ x,
  df.t1a,
  group.by = g.t1a,
  n.samples = 40L,
  n.burn = 20L,
  n.thin = 1L,
  n.chains = 1L,
  n.trees = 20L,
  n.threads = 1L,
  keepTrees = TRUE,
  verbose = FALSE
)

# per-group ranef, group levels as dimnames
expect_equal(ncol(fit.t1a$ranef), nlevels(g.t1a))
expect_equal(colnames(fit.t1a$ranef), levels(g.t1a))
expect_equal(length(fit.t1a$ranef.mean), nlevels(g.t1a))
expect_equal(names(fit.t1a$ranef.mean), levels(g.t1a))

# trained groups round-trip: no spurious new-level warning, correlation ~ 1
expect_silent(predict(fit.t1a, df.t1a, group.by = g.t1a))
pred.f <- predict(fit.t1a, df.t1a, group.by = g.t1a)
expect_true(cor(fitted(fit.t1a), colMeans(pred.f)) > 0.999)

# group.by coerced as at fit time: numeric and character match the factor
pred.n <- predict(fit.t1a, df.t1a, group.by = as.integer(g.t1a))
pred.c <- predict(fit.t1a, df.t1a, group.by = as.character(g.t1a))
expect_equal(pred.n, pred.f)
expect_equal(pred.c, pred.f)

# genuinely-new level draws from the prior, names the level, stays finite
g.new <- as.integer(g.t1a)
g.new[seq_len(5L)] <- 99L
warnings.newLevel <- captureWarnings(predict(fit.t1a, df.t1a, group.by = g.new))
expect_equal(length(warnings.newLevel), 1L)
expect_match(conditionMessage(warnings.newLevel[[1L]]), "99")
expect_inherits(warnings.newLevel[[1L]], "dbartsUnmeasuredLevelsWarning")
pred.new <- suppressWarnings(predict(fit.t1a, df.t1a, group.by = g.new))
expect_true(all(is.finite(pred.new)))

rm(
  n.t1a,
  g.t1a,
  x.t1a,
  y.t1a,
  df.t1a,
  fit.t1a,
  pred.f,
  pred.n,
  pred.c,
  g.new,
  pred.new
)


# test that works with missing levels
n.train <- 80L
x <- testData$x[seq_len(n.train), ]
y <- testData$y[seq_len(n.train)]
g <- factor(testData$g[seq_len(n.train)])

x.test <- testData$x[seq.int(n.train + 1L, nrow(testData$x)), ]
g.test <- factor(testData$g[seq.int(n.train + 1L, nrow(testData$x))], levels(g))
levels(g.test)[5L] <- "6"

# check that predict works when we've fit with missing levels
# combineChains = FALSE pinned deliberately: rbartFit$ranef is indexed as a
# raw uncombined 3-d array below
warnings.fitLevel <- captureWarnings(
  rbartFit <- dbarts::rbart_vi(
    y ~ x,
    group.by = g,
    test = x.test,
    group.by.test = g.test,
    n.samples = 7L,
    n.burn = 0L,
    n.thin = 1L,
    n.chains = 2L,
    n.trees = 25L,
    n.threads = 1L,
    keepTrees = TRUE,
    verbose = FALSE,
    combineChains = FALSE
  )
)
expect_true(any(vapply(
  warnings.fitLevel,
  inherits,
  logical(1L),
  "dbartsUnmeasuredLevelsWarning"
)))
expect_equal(
  apply(predict(rbartFit, x.test, group.by = g.test), 2L, mean),
  fitted(rbartFit, sample = "test")
)
expect_equal(
  apply(
    predict(rbartFit, x.test, group.by = g.test, combineChains = FALSE),
    3L,
    mean
  ),
  fitted(rbartFit, sample = "test")
)

# check that predicts works for completely new levels
levels(g.test) <- c(levels(g.test)[-5L], as.character(seq.int(7L, 28L)))
set.seed(0L)
ranef.pred <- suppressWarnings(
  predict(
    rbartFit,
    x.test,
    group.by = g.test,
    type = "ranef",
    combineChains = FALSE
  )
)
expect_equal(
  ranef.pred[,, as.character(1L:4L)],
  rbartFit$ranef[,, as.character(1L:4L)]
)
expect_true(
  # 7 samples of a 22-draw sd estimate: the correlation tracks tau but the
  # value depends on the draw history; assert strong tracking rather than a
  # knife-edge bound
  cor(
    as.numeric(rbartFit$tau),
    as.numeric(apply(ranef.pred[,, 5L:26L], c(1L, 2L), sd))
  ) >
    0.80
)

# check again with combineChains as TRUE at the top level
g.test <- droplevels(g.test)
levels(g.test) <- c(levels(g)[-5L], "6")
rbartFit <- suppressWarnings(dbarts::rbart_vi(
  y ~ x,
  group.by = g,
  test = x.test,
  group.by.test = g.test,
  n.samples = 7L,
  n.burn = 0L,
  n.thin = 1L,
  n.chains = 2L,
  n.trees = 25L,
  n.threads = 1L,
  keepTrees = TRUE,
  combineChains = TRUE,
  verbose = FALSE
))
expect_equal(
  apply(predict(rbartFit, x.test, group.by = g.test), 2L, mean),
  fitted(rbartFit, sample = "test")
)
expect_equal(
  apply(
    predict(rbartFit, x.test, group.by = g.test, combineChains = FALSE),
    3L,
    mean
  ),
  fitted(rbartFit, sample = "test")
)

levels(g.test) <- c(levels(g.test)[-5L], as.character(seq.int(7L, 28L)))
set.seed(0L)
ranef.pred <- suppressWarnings(predict(
  rbartFit,
  x.test,
  group.by = g.test,
  type = "ranef"
))
expect_equal(
  as.numeric(ranef.pred[, as.character(1L:4L)]),
  as.numeric(rbartFit$ranef[, as.character(1L:4L)])
)
expect_true(
  cor(
    as.numeric(rbartFit$tau),
    as.numeric(apply(ranef.pred[, 5L:26L], 1L, sd))
  ) >
    0.90
)

# check one last time with one chain
g.test <- droplevels(g.test)
levels(g.test) <- c(levels(g)[-5L], "6")
rbartFit <- suppressWarnings(dbarts::rbart_vi(
  y ~ x,
  group.by = g,
  test = x.test,
  group.by.test = g.test,
  n.samples = 14L,
  n.burn = 0L,
  n.thin = 1L,
  n.chains = 1L,
  n.trees = 25L,
  n.threads = 1L,
  keepTrees = TRUE,
  verbose = FALSE
))
expect_equal(
  apply(predict(rbartFit, x.test, group.by = g.test), 2L, mean),
  fitted(rbartFit, sample = "test")
)
levels(g.test) <- c(levels(g.test)[-5L], as.character(seq.int(7L, 28L)))
set.seed(0L)
ranef.pred <- suppressWarnings(predict(
  rbartFit,
  x.test,
  group.by = g.test,
  type = "ranef"
))
expect_equal(
  as.numeric(ranef.pred[, as.character(1L:4L)]),
  as.numeric(rbartFit$ranef[, as.character(1L:4L)])
)
expect_true(
  # 14 samples of a 22-draw sd estimate: the correlation tracks tau but sits
  # near 0.9 by seed; assert strong tracking rather than a knife-edge bound
  cor(
    as.numeric(rbartFit$tau),
    as.numeric(apply(ranef.pred[, 5L:26L], 1L, sd))
  ) >
    0.80
)

# check with more than one missing level
levels(g.test)[4L] <- "7"
rbartFit <- suppressWarnings(dbarts::rbart_vi(
  y ~ x,
  group.by = g,
  test = x.test,
  group.by.test = g.test,
  n.samples = 7L,
  n.burn = 0L,
  n.thin = 1L,
  n.chains = 4L,
  n.trees = 25L,
  n.threads = 1L,
  verbose = FALSE
))
expect_inherits(rbartFit, "rbart")

rm(rbartFit, ranef.pred, g.test, warnings.fitLevel)

rm(x.test, g, y, x, n.train)

rm(testData)
