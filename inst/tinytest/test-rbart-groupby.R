source(
  system.file("common", "friedmanData.R", package = "dbarts"),
  local = TRUE
)

# the port's deprecation warning is pinned in test-rbart-port.R
onceState <- dbarts:::onceWarnState
onceState[["tombstone.rbart_vi"]] <- TRUE
rm(onceState)

n.g <- 5L
if (getRversion() >= "3.6.0") {
  oldSampleKind <- RNGkind()[3L]
  suppressWarnings(RNGkind(sample.kind = "Rounding"))
}
g <- sample(n.g, length(testData$y), replace = TRUE)
if (getRversion() >= "3.6.0") {
  suppressWarnings(RNGkind(sample.kind = oldSampleKind))
  rm(oldSampleKind)
}

sigma.b <- 1.5
b <- rnorm(n.g, 0, sigma.b)

testData$y <- testData$y + b[g]
testData$g <- g
testData$b <- b
rm(b, sigma.b, g, n.g)

# test that rbart finds group.by
df <- as.data.frame(testData$x)
colnames(df) <- paste0("x_", seq_len(ncol(testData$x)))
df$y <- testData$y
df$g <- testData$g
set.seed(11L)
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
set.seed(22L)
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
set.seed(33L)
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


# test that works with missing levels
n.train <- 80L
x <- testData$x[seq_len(n.train), ]
y <- testData$y[seq_len(n.train)]
g <- factor(testData$g[seq_len(n.train)])

x.test <- testData$x[seq.int(n.train + 1L, nrow(testData$x)), ]
g.test <- factor(testData$g[seq.int(n.train + 1L, nrow(testData$x))], levels(g))
levels(g.test)[5L] <- "6"

# check that predict works when we've fit with missing levels
set.seed(44L)
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
  verbose = FALSE
))
expect_equal(
  apply(predict(rbartFit, x.test, g.test), 2L, mean),
  fitted(rbartFit, sample = "test")
)
expect_equal(
  apply(predict(rbartFit, x.test, g.test, combineChains = FALSE), 3L, mean),
  fitted(rbartFit, sample = "test")
)

# check that predicts works for completely new levels
levels(g.test) <- c(levels(g.test)[-5L], as.character(seq.int(7L, 28L)))
set.seed(0L)
ranef.pred <- suppressWarnings(
  predict(rbartFit, x.test, g.test, type = "ranef", combineChains = FALSE)
)
expect_equal(
  ranef.pred[,, as.character(1L:4L)],
  rbartFit$ranef[,, as.character(1L:4L)]
)
expect_true(
  cor(
    as.numeric(rbartFit$tau),
    as.numeric(apply(ranef.pred[,, 5L:26L], c(1L, 2L), sd))
  ) >
    0.90
)

# check again with combineChains as TRUE at the top level
g.test <- droplevels(g.test)
levels(g.test) <- c(levels(g)[-5L], "6")
set.seed(55L)
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
  apply(predict(rbartFit, x.test, g.test), 2L, mean),
  fitted(rbartFit, sample = "test")
)
expect_equal(
  apply(predict(rbartFit, x.test, g.test, combineChains = FALSE), 3L, mean),
  fitted(rbartFit, sample = "test")
)

levels(g.test) <- c(levels(g.test)[-5L], as.character(seq.int(7L, 28L)))
set.seed(0L)
ranef.pred <- suppressWarnings(predict(
  rbartFit,
  x.test,
  g.test,
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
set.seed(66L)
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
  apply(predict(rbartFit, x.test, g.test), 2L, mean),
  fitted(rbartFit, sample = "test")
)
levels(g.test) <- c(levels(g.test)[-5L], as.character(seq.int(7L, 28L)))
set.seed(0L)
ranef.pred <- suppressWarnings(predict(
  rbartFit,
  x.test,
  g.test,
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

# check with more than one missing level
levels(g.test)[4L] <- "7"
set.seed(77L)
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

rm(rbartFit, ranef.pred, g.test)

rm(x.test, g, y, x, n.train)

rm(testData)

# group.by, group.by.test and prior resolve through a forwarding wrapper: a
# column name is read as a column, a caller variable by its value
dfFwd <- data.frame(
  y = rnorm(40L),
  x = rnorm(40L),
  g = rep_len(1:4, 40L)
)
fwdRbart <- function(...) {
  dbarts::rbart_vi(
    ...,
    n.samples = 5L,
    n.burn = 0L,
    n.thin = 1L,
    n.chains = 1L,
    n.trees = 5L,
    n.threads = 1L,
    seed = 3L,
    verbose = FALSE
  )
}
fwdMiddle <- function(...) fwdRbart(...)
# a wrapper that renames: the parameter is called group.by-ish, and the data
# has a column of that name which must NOT be picked up
fwdRenamed <- function(d, grp) fwdRbart(y ~ x, d, group.by = grp)
# the group count is read off the fit's random effects
nGroups <- function(fit) ncol(fit$ranef)
direct <- dbarts::rbart_vi(
  y ~ x,
  dfFwd,
  group.by = g,
  n.samples = 5L,
  n.burn = 0L,
  n.thin = 1L,
  n.chains = 1L,
  n.trees = 5L,
  n.threads = 1L,
  seed = 3L,
  verbose = FALSE
)
expect_equal(nGroups(direct), 4L)
expect_equal(nGroups(fwdRbart(y ~ x, dfFwd, group.by = g)), 4L)
expect_equal(nGroups(fwdMiddle(y ~ x, dfFwd, group.by = g)), 4L)
expect_equal(fwdRbart(y ~ x, dfFwd, group.by = g)$ranef, direct$ranef)

# a local variable of the user's own name beats a global of that name
g <- rep_len(1:2, 40L)
localGroups <- (function() {
  g <- rep_len(1:5, 40L)
  fwdMiddle(y ~ x, dfFwd[c("y", "x")], group.by = g)
})()
expect_equal(nGroups(localGroups), 5L)
# and a column is read as a column, whatever the global says
expect_equal(nGroups(fwdMiddle(y ~ x, dfFwd, group.by = g)), 4L)
# the data has a column named like the intermediate wrapper's parameter: the
# column IS taken, forwarded or not, as lm does for weights and subset
dfGrp <- dfFwd
dfGrp$grp <- rep_len(1:3, 40L)
expect_equal(nGroups(fwdRenamed(dfGrp, rep_len(1:5, 40L))), 3L)
# a closure calling rbart_vi
closureGroups <- (function() {
  inner <- function(...) {
    dbarts::rbart_vi(
      ...,
      n.samples = 5L,
      n.burn = 0L,
      n.thin = 1L,
      n.chains = 1L,
      n.trees = 5L,
      n.threads = 1L,
      seed = 3L,
      verbose = FALSE
    )
  }
  gLocal <- rep_len(1:5, 40L)
  inner(y ~ x, dfFwd, group.by = gLocal)
})()
expect_equal(nGroups(closureGroups), 5L)
rm(g)

# group.by.test through a wrapper is the caller's local variable, not a
# same-named global: the fit carries the caller's test groups
gt <- rep_len(1:2, 40L)
gtLocal <- rep_len(c(4L, 3L, 2L, 1L), 40L)
testFit <- (function() {
  gt <- gtLocal
  fwdMiddle(y ~ x, dfFwd, test = dfFwd, group.by = g, group.by.test = gt)
})()
expect_equal(as.character(testFit$group.by.test), as.character(gtLocal))
rm(gt, gtLocal, testFit)

# prior = gamma is the gamma prior, not the default
fwdGamma <- fwdMiddle(y ~ x, dfFwd, group.by = g, prior = gamma)
expect_false(isTRUE(all.equal(fwdGamma$ranef, direct$ranef)))
expect_equal(
  fwdGamma$ranef,
  dbarts::rbart_vi(
    y ~ x,
    dfFwd,
    group.by = g,
    prior = gamma,
    n.samples = 5L,
    n.burn = 0L,
    n.thin = 1L,
    n.chains = 1L,
    n.trees = 5L,
    n.threads = 1L,
    seed = 3L,
    verbose = FALSE
  )$ranef
)
rm(
  fwdRbart,
  fwdMiddle,
  fwdRenamed,
  nGroups,
  direct,
  localGroups,
  dfGrp,
  closureGroups,
  testFit,
  fwdGamma,
  dfFwd
)


# predict takes a character or numeric group.by as the fit does: random
# effects by level, and an unseen level drawn from its distribution.
groupFrame <- data.frame(y = rnorm(60L), x = runif(60L))
characterGroups <- rep(c("a", "b", "c"), 20L)
characterFit <- suppressWarnings(rbart_vi(
  y ~ x,
  groupFrame,
  group.by = characterGroups,
  n.samples = 10L,
  n.burn = 10L,
  n.chains = 1L,
  n.trees = 5L,
  n.threads = 1L,
  keepTrees = TRUE,
  verbose = FALSE
))
ranef <- predict(
  characterFit,
  groupFrame[1:3, ],
  group.by = c("a", "b", "c"),
  type = "ranef"
)
expect_equal(
  ranef,
  predict(
    characterFit,
    groupFrame[1:3, ],
    group.by = factor(c("a", "b", "c")),
    type = "ranef"
  )
)
expect_equal(ncol(ranef), 3L)
expect_warning(
  unseen <- predict(
    characterFit,
    groupFrame[1:3, ],
    group.by = c("a", "b", "zz")
  ),
  "not present in training"
)
expect_equal(ncol(unseen), 3L)
numericFit <- suppressWarnings(rbart_vi(
  y ~ x,
  groupFrame,
  group.by = rep(1:3, 20L),
  n.samples = 10L,
  n.burn = 10L,
  n.chains = 1L,
  n.trees = 5L,
  n.threads = 1L,
  keepTrees = TRUE,
  verbose = FALSE
))
expect_warning(
  unseen <- predict(numericFit, groupFrame[1:3, ], group.by = c(1, 2, 9)),
  "not present in training"
)
expect_equal(ncol(unseen), 3L)
rm(groupFrame, characterGroups, characterFit, ranef, unseen, numericFit)
