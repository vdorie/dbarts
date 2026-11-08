# A hurdle fit takes a subset and weights of 0 and 1. A subset removes rows from
# both parts; a weight of 0 holds a row out of both likelihoods and keeps it in
# both designs. Each part is the standalone probit or gaussian fit it is called.

set.seed(5301L)
n <- 90L
d <- data.frame(a = runif(n), b = runif(n))
d$f <- factor(sample(c("u", "v", "w"), n, TRUE))
d$y <- ifelse(runif(n) < 0.5, 0, exp(rnorm(n)))
d$y[1L] <- 0
d$y[2L] <- 1.2
d$w <- 1
heldOut <- c(2L, 3L, 10L, 11L, 40L)
d$w[heldOut] <- 0
expect_true(any(d$y[heldOut] > 0) && any(d$y[heldOut] == 0))

seed <- 5302L
fitArgs <- list(
  n.samples = 6L,
  n.burn = 6L,
  n.trees = 8L,
  n.chains = 2L,
  n.threads = 1L,
  verbose = FALSE,
  keepTrees = TRUE
)
fitHurdle <- function(formula, ...) {
  # arguments are evaluated where bart reads them, so a column named in
  # weights or subset resolves in the data
  fitCall <- match.call()
  fitCall[[1L]] <- quote(bart)
  fitCall$family <- "hurdle.lognormal"
  fitCall[names(fitArgs)] <- fitArgs
  fitCall$seed <- seed
  eval(fitCall, parent.frame())
}
partSeeds <- dbarts:::withFixedSeed(seed, sample.int(.Machine$integer.max, 2L))
standalone <- function(formula, ..., family, seed) {
  do.call(bart, c(list(formula, ...), fitArgs, family = family, seed = seed))
}

# --- each part is the standalone fit with the same subset and weights ---------------
inSubset <- d$a > 0.15
fit <- fitHurdle(y ~ a + b + f, data = d, weights = w, subset = a > 0.15)
kept <- d[inSubset, ]
expect_identical(length(fit$y), nrow(kept))
zeroAlone <- standalone(
  as.double(y > 0) ~ a + b + f,
  data = kept,
  weights = kept$w,
  family = "probit",
  seed = partSeeds[1L]
)
expect_identical(fit$zero$yhat.train, zeroAlone$yhat.train)
positiveRows <- kept$y > 0
positiveAlone <- standalone(
  log(y) ~ a + b + f,
  data = kept[positiveRows, ],
  weights = quote(w),
  family = "gaussian",
  seed = partSeeds[2L],
  test = kept
)
expect_identical(fit$positive$yhat.train, positiveAlone$yhat.train)
expect_identical(fit$positive$yhat.test, positiveAlone$yhat.test)
expect_identical(fit$weights, kept$w)

# the matrix interface, with the subset an index, is the same fit
viaMatrix <- fitHurdle(
  d[, c("a", "b", "f")],
  d$y,
  weights = d$w,
  subset = which(inSubset)
)
expect_identical(unname(extract(viaMatrix, "ev")), unname(extract(fit, "ev")))
expect_identical(viaMatrix$weights, fit$weights)

# --- subset before na.action, as model.frame orders them ------------------------------
dGap <- d
dGap$a[8L] <- NA
for (action in list(na.omit, na.exclude)) {
  viaLm <- lm(y ~ a + b, dGap, subset = 2:40, na.action = action)
  viaHurdle <- fitHurdle(
    y ~ a + b,
    data = dGap,
    subset = 2:40,
    na.action = action
  )
  expect_identical(viaHurdle$na.action, viaLm$na.action)
  expect_identical(
    length(fitted(viaHurdle)),
    length(fitted(viaLm))
  )
}

# --- the weights ---------------------------------------------------------------------
expect_error(
  fitHurdle(y ~ a + b, data = d, weights = 2 * w),
  "weights of 0 and 1 only"
)
expect_error(
  fitHurdle(y ~ a + b, data = d, weights = 2 * w),
  "repeat each row its weight's number of times"
)
expect_error(
  fitHurdle(d[, c("a", "b")], d$y, weights = d$w / 2),
  "weights of 0 and 1 only"
)
held <- fitHurdle(y ~ a + b + f, data = d, weights = w)
# a 0 row stays in both designs and is out of both likelihoods
expect_identical(nrow(held$zero$fit$data@x), n)
expect_identical(ncol(held$positive$yhat.test), n)
expect_identical(held$weights, d$w)
# fitted and residuals are reported at a 0 row, as lm reports them
expect_false(anyNA(fitted(held)))
expect_equal(unname(residuals(held)), d$y - unname(fitted(held)))
# the predictive draws are NA there, the log-likelihood NaN
set.seed(1L)
ppd <- extract(held, type = "ppd")
expect_identical(unname(which(is.na(ppd[1L, ]))), heldOut)
expect_true(all(colSums(is.na(ppd)) == ifelse(d$w == 0, nrow(ppd), 0)))
loglik <- extract(held, type = "loglik")
expect_identical(unname(which(is.nan(loglik[1L, ]))), heldOut)
expect_false(anyNA(loglik[, d$w == 1]))
# and the draw is made before the row is set aside, so the stream and the other
# rows' draws are those of the fit without weights
set.seed(1L)
ppdAll <- extract(fitHurdle(y ~ a + b + f, data = d), type = "ppd")
expect_equal(dim(ppdAll), dim(ppd))
# weights.type: NULL and "frequency" are one reading, 0/1 weights being rows in
# and out of both parts, and "precision" is refused in the hurdle's words
stated <- fitHurdle(
  y ~ a + b + f,
  data = d,
  weights = w,
  weights.type = "frequency"
)
expect_identical(extract(stated, "ev"), extract(held, "ev"))
set.seed(1L)
expect_identical(extract(stated, type = "ppd"), ppd)
expect_error(
  fitHurdle(y ~ a + b, data = d, weights = w, weights.type = "precision"),
  "is not taken under family = \"hurdle.lognormal\"",
  fixed = TRUE
)
# a test frame's weights column gives the test rows' weights, a 0 an NA draw
te <- d[1:10, ]
te$w <- c(0, 1, 1, 1, 1, 0, 1, 1, 1, 1)
withTest <- fitHurdle(y ~ a + b + f, data = d, weights = w, test = te)
expect_identical(withTest$weights.test, te$w)
set.seed(2L)
testPpd <- extract(withTest, type = "ppd", sample = "test")
expect_identical(unname(which(is.na(testPpd[1L, ]))), c(1L, 6L))
expect_false(anyNA(extract(withTest, type = "ev", sample = "test")))
# the training channel is unchanged by a test set
expect_identical(extract(withTest, "ev"), extract(held, "ev"))

# --- a missing weight follows na.action on both routes (dec-A89) -----------------------
dNA <- d
dNA$w[5L] <- NA
for (action in list(na.omit, na.exclude)) {
  viaFormula <- fitHurdle(
    y ~ a + b + f,
    data = dNA,
    weights = w,
    na.action = action
  )
  viaMatrix <- fitHurdle(
    dNA[, c("a", "b", "f")],
    dNA$y,
    weights = dNA$w,
    na.action = action
  )
  expect_identical(viaMatrix$na.action, viaFormula$na.action)
  expect_identical(
    unname(extract(viaMatrix, "ev")),
    unname(extract(viaFormula, "ev"))
  )
  expect_identical(length(viaMatrix$y), n - 1L)
}
expect_identical(
  length(fitted(fitHurdle(
    dNA[, c("a", "b", "f")],
    dNA$y,
    weights = dNA$w,
    na.action = na.exclude
  ))),
  n
)
expect_error(
  fitHurdle(y ~ a + b, data = dNA, weights = w),
  "cannot be NA|missing"
)
expect_error(
  fitHurdle(dNA[, c("a", "b")], dNA$y, weights = dNA$w),
  "'weights' cannot be NA"
)

# --- predict's weights -------------------------------------------------------------------
# with no newdata a weight of 1 draws the 0 rows, and at the weight-1 rows the
# draws are those extract gives
set.seed(3L)
fresh <- predict(held, type = "ppd", weights = 1)
expect_false(anyNA(fresh))
set.seed(3L)
stored <- extract(held, type = "ppd")
expect_identical(fresh[, d$w == 1], stored[, d$w == 1])
set.seed(3L)
expect_identical(predict(held, type = "ppd", weights = rep(1, n)), fresh)
# a given weight replaces the fit's own, a 0 holding its row out
oneOut <- replace(rep(1, n), 7L, 0)
set.seed(3L)
expect_identical(
  unname(which(is.na(predict(held, type = "ppd", weights = oneOut)[1L, ]))),
  7L
)
expect_error(
  predict(held, type = "ppd", weights = c(1, 1)),
  "one number or one per training row"
)
expect_error(predict(held, type = "ppd", weights = 2), "0 and 1")
# types that draw nothing read no weight
expect_identical(
  predict(held, type = "ev", weights = oneOut),
  predict(held, type = "ev")
)
# a missing weight is its row's na.action's to answer (dec-B482)
withGap <- replace(rep(1, n), 4L, NA)
expect_identical(
  dim(predict(held, type = "ppd", weights = withGap, na.action = na.omit)),
  c(12L, n - 1L)
)
expect_identical(
  dim(predict(held, type = "ppd", weights = withGap, na.action = na.exclude)),
  c(12L, n)
)
expect_error(predict(held, type = "ppd", weights = withGap), "missing values")
# at newdata a weight of 1 draws a row as a new observation, 0 holds it out
newdata <- d[1:6, ]
set.seed(4L)
drawn <- predict(held, newdata, type = "ppd", weights = c(1, 0, 1, 1, 1, 1))
expect_identical(unname(which(is.na(drawn[1L, ]))), 2L)
# no weights stated: the weights column of newdata is read, a 0 holding its
# row out (dec-B648), and newdata without the column draws every row
set.seed(4L)
expect_identical(
  unname(which(is.na(predict(held, newdata, type = "ppd")[1L, ]))),
  2:3
)
set.seed(4L)
expect_false(anyNA(predict(held, newdata[c("a", "b", "f")], type = "ppd")))
expect_error(
  predict(held, newdata, type = "ppd", weights = 1:2),
  "same number of rows as 'newdata'"
)

# --- partial dependence leaves the 0 rows out ----------------------------------------------
grid <- c(0.3, 0.6)
pd <- pdbart(held, xind = "a", levs = list(grid), pl = FALSE)
manual <- vapply(
  grid,
  function(level) {
    rows <- d[d$w == 1, ]
    rows$a <- level
    rowMeans(predict(held, rows, type = "ev"))
  },
  numeric(2L * 6L)
)
expect_equal(unname(pd$fd[[1L]]), unname(manual), tolerance = 1e-8)
