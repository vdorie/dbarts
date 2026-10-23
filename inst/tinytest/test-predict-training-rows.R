# predict with no newdata returns every type at the training rows from the
# stored draws, needing no saved trees; an offset given there replaces the
# fit's own (dec-B341)

set.seed(21)
n <- 50L
x <- matrix(runif(n * 3L), n, dimnames = list(NULL, c("a", "b", "c")))
df <- data.frame(x)
df$y <- x[, 1L] + rnorm(n, 0, 0.3)
df$z <- as.integer(df$y + rnorm(n, 0, 0.3) > 0.5)
off <- rnorm(n, 0, 0.5)
new.off <- rnorm(n, 0, 0.5)

quick <- function(..., keepTrees = TRUE) {
  suppressWarnings(bart(
    ...,
    n.trees = 10L,
    n.burn = 10L,
    n.samples = 10L,
    n.chains = 2L,
    n.threads = 1L,
    keepTrees = keepTrees,
    verbose = FALSE
  ))
}

# gaussian: equal to extract, with and without saved trees
fit <- quick(y ~ a + b + c, df)
fit.bare <- quick(y ~ a + b + c, df, keepTrees = FALSE)
expect_null(fit.bare$fit)
for (type in c("ev", "bart")) {
  expect_equal(
    predict(fit, type = type),
    extract(fit, type, sample = "train")
  )
  expect_equal(
    predict(fit.bare, type = type, combineChains = FALSE),
    extract(fit.bare, type, sample = "train", combineChains = FALSE)
  )
}
set.seed(3)
a <- predict(fit.bare, type = "ppd")
set.seed(3)
b <- extract(fit.bare, "ppd", sample = "train")
expect_equal(a, b)
ci <- predict(fit.bare, ci.level = 0.9)
expect_equal(dim(ci), c(n, 3L))
# newdata = NULL is the same request
expect_equal(predict(fit.bare, NULL), predict(fit.bare))

# an offset replaces the fit's: none to start from
expect_equal(
  predict(fit, offset = new.off),
  predict(fit, df, offset = new.off)
)
expect_equal(
  predict(fit, offset = 2),
  predict(fit, df, offset = rep(2, n))
)
expect_error(predict(fit.bare, offset = 1:3), "one per training row")
expect_error(predict(fit.bare, offset = c(NA, rep(1, n - 1L))), "none missing")
expect_error(predict(fit.bare, weights = 1), "rows of 'newdata'")

# an offset given at fit time is replaced, not added to
fit.off <- quick(y ~ a + b + c + offset(off), df)
fit.off.bare <- quick(y ~ a + b + c + offset(off), df, keepTrees = FALSE)
expect_equal(predict(fit.off, type = "bart"), extract(fit.off, "bart"))
expect_equal(
  predict(fit.off, type = "bart", offset = new.off),
  sweep(extract(fit.off, "bart"), 2L, off - new.off, "-"),
  check.attributes = FALSE
)
# the fit's own offset given back changes nothing
expect_equal(
  predict(fit.off.bare, type = "bart", offset = off),
  predict(fit.off.bare, type = "bart")
)

# a binary fit shifts the latent before the link
fit.p <- quick(z ~ a + b + c, df, family = "probit")
expect_equal(
  predict(fit.p, offset = new.off),
  predict(fit.p, df, offset = new.off)
)
expect_equal(
  predict(fit.p, type = "bart", offset = new.off, combineChains = FALSE),
  predict(fit.p, df, type = "bart", offset = new.off, combineChains = FALSE)
)

# a count fit
df$cnt <- rpois(n, exp(0.5 + df$a))
fit.nb <- quick(cnt ~ a + b + c, df, family = "nbinom")
expect_equal(predict(fit.nb, type = "bart"), extract(fit.nb, "bart"))
expect_equal(
  predict(fit.nb, offset = new.off),
  predict(fit.nb, df, offset = new.off)
)

# an ordered response
df$ord <- cut(df$y, 3L, labels = c("lo", "mid", "hi"), ordered_result = TRUE)
fit.ord <- quick(ord ~ a + b + c, df)
expect_equal(predict(fit.ord), extract(fit.ord, "ev"))
expect_equal(
  predict(fit.ord, offset = new.off),
  predict(fit.ord, df, offset = new.off)
)
expect_equal(
  predict(fit.ord, type = "class"),
  predict(fit.ord, df, type = "class")
)

# an unordered response
df$cat <- factor(cut(df$y, 3L, labels = c("p", "q", "r")))
fit.mn <- quick(cat ~ a + b + c, df)
expect_equal(predict(fit.mn), extract(fit.mn, "ev"))
cat.off <- matrix(rnorm(n * 3L, 0, 0.5), n, 3L)
expect_equal(
  predict(fit.mn, offset = cat.off),
  predict(fit.mn, df, offset = cat.off),
  tolerance = 1e-6
)
# a fit trained with a category offset has it replaced, not added to
train.cat.off <- matrix(rnorm(n * 3L, 0, 0.8), n, 3L)
fit.mn.off <- suppressWarnings(bart(
  x,
  df$cat,
  offset = train.cat.off,
  n.trees = 10L,
  n.burn = 10L,
  n.samples = 10L,
  n.chains = 2L,
  n.threads = 1L,
  keepTrees = TRUE,
  verbose = FALSE
))
new.cat.off <- matrix(rnorm(n * 3L, 0, 0.8), n, 3L)
expect_equal(
  predict(fit.mn.off, offset = new.cat.off),
  predict(fit.mn.off, x, offset = new.cat.off),
  tolerance = 1e-6
)
# the fit's own offset given back changes nothing
expect_equal(
  predict(fit.mn.off, offset = train.cat.off),
  extract(fit.mn.off, "ev"),
  tolerance = 1e-8
)

# a hurdle fit has no offset channel
df$h <- ifelse(df$z == 1L, exp(df$y), 0)
fit.h <- quick(x, df$h, family = "hurdle.lognormal")
expect_equal(predict(fit.h), extract(fit.h, "ev"))
expect_error(predict(fit.h, offset = 1), "no out-of-sample offset channel")

# survivalProbabilities with an offset and no newdata applies it at the
# training rows off the stored draws (dec-B340), trees or none
tm <- exp(df$y + rnorm(n, 0, 0.2))
st <- cbind(tm, as.integer(runif(n) < 0.8))
fit.aft <- quick(x, st, family = "aft")
fit.aft.bare <- quick(x, st, family = "aft", keepTrees = FALSE)
times <- stats::quantile(tm, c(0.25, 0.75))
expect_equal(
  survivalProbabilities(fit.aft, times, offset = new.off),
  survivalProbabilities(fit.aft, times, newdata = x, offset = new.off)
)
expect_equal(
  survivalProbabilities(fit.aft, times, offset = 0),
  survivalProbabilities(fit.aft, times)
)
expect_equal(
  dim(survivalProbabilities(fit.aft.bare, times, offset = new.off)),
  c(20L, 2L, n)
)
expect_error(
  survivalProbabilities(fit.aft.bare, times, offset = 1:3),
  "one per training row"
)

# a row at weight 0 has NA posterior predictive draws by every route, with no
# warning: a weight of 0 holds the row out of the fit and gives it no noise
# scale (sd sigma / sqrt(w)); every other row's draw is what it was
w <- rep(c(0, 1, 4), length.out = n)
fit.w <- suppressWarnings(bart(
  x,
  df$y,
  weights = w,
  n.trees = 10L,
  n.burn = 10L,
  n.samples = 10L,
  n.chains = 2L,
  n.threads = 1L,
  keepTrees = TRUE,
  verbose = FALSE
))
counted <- function(expr) {
  seen <- 0L
  value <- withCallingHandlers(
    expr,
    warning = function(cond) {
      seen <<- seen + 1L
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = seen)
}
held <- which(w == 0)
kept <- which(w != 0)
# extract: the held rows are NA in every draw, and the rows the weights keep
# are the draws of the build before, which drew a weight-0 row at weight 1
set.seed(8)
ppd.w <- counted(extract(fit.w, "ppd"))
expect_identical(ppd.w$warnings, 0L)
expect_true(all(is.na(ppd.w$value[, held])))
expect_false(anyNA(ppd.w$value[, kept]))
fit.w1 <- fit.w
fit.w1$weights <- ifelse(w == 0, 1, w)
set.seed(8)
ppd.w1 <- extract(fit.w1, "ppd")
expect_identical(ppd.w$value[, kept], ppd.w1[, kept])
fitted.w <- counted(fitted(fit.w, type = "ppd"))
expect_identical(fitted.w$warnings, 0L)
expect_identical(which(is.na(fitted.w$value)), held)
ci.w <- counted(fitted(fit.w, type = "ppd", ci.level = 0.9))
expect_identical(ci.w$warnings, 0L)
expect_identical(which(is.na(ci.w$value[, "est"])), held)
expect_identical(which(is.na(ci.w$value[, "ci.lower"])), held)
expect_identical(which(is.na(ci.w$value[, "ci.upper"])), held)
# the other types are untouched by a weight of 0
expect_false(anyNA(extract(fit.w, "ev")))
expect_false(anyNA(fitted(fit.w, type = "bart")))
# a split-chain request holds the same rows out
expect_identical(
  which(apply(is.na(extract(fit.w, "ppd", combineChains = FALSE)), 3L, any)),
  held
)

# predict at newdata with weights = 0 is NA at those rows, where it was
# refused by name; a row stated at weight 1 is the way to draw it as a new
# observation
nd.w <- c(0, 1, 4, 0, 2)
set.seed(9)
pw <- counted(predict(fit.w, x[1:5, ], type = "ppd", weights = nd.w))
expect_identical(pw$warnings, 0L)
expect_identical(which(colSums(is.na(pw$value)) > 0L), c(1L, 4L))
set.seed(9)
pw1 <- predict(
  fit.w,
  x[1:5, ],
  type = "ppd",
  weights = ifelse(nd.w == 0, 1, nd.w)
)
expect_identical(pw$value[, c(2L, 3L, 5L)], pw1[, c(2L, 3L, 5L)])
expect_false(anyNA(pw1))
expect_true(all(is.na(predict(fit.w, x[1:5, ], type = "ppd", weights = 0))))
expect_identical(
  which(is.na(predict(
    fit.w,
    x[1:5, ],
    type = "ppd",
    weights = nd.w,
    ci.level = 0.9
  )[, "est"])),
  c(1L, 4L)
)
expect_false(anyNA(predict(fit.w, x[held, ], type = "ppd", weights = 1)))
# only a weight of 0 is held out: the rest of the check is as it was
expect_error(
  predict(fit.w, x[1:5, ], type = "ppd", weights = c(-1, 1, 1, 1, 1)),
  "finite and non-negative"
)
expect_error(
  predict(fit.w, x[1:5, ], type = "ppd", weights = c(Inf, 1, 1, 1, 1)),
  "finite and non-negative"
)

# the test sample's rows whose test-frame weight is 0
train.frame <- data.frame(x, y = df$y, weights = w)
test.rows <- 1:9
test.frame <- data.frame(
  x[test.rows, ],
  y = df$y[test.rows],
  weights = w[test.rows]
)
fit.wt <- suppressWarnings(bart(
  y ~ a + b + c,
  train.frame,
  test = test.frame,
  weights = weights,
  n.trees = 10L,
  n.burn = 10L,
  n.samples = 10L,
  n.chains = 2L,
  n.threads = 1L,
  keepTrees = TRUE,
  verbose = FALSE
))
expect_identical(fit.wt$weights.test, w[test.rows])
set.seed(4)
ppd.t <- counted(extract(fit.wt, "ppd", sample = "test"))
expect_identical(ppd.t$warnings, 0L)
heldTest <- which(w[test.rows] == 0)
keptTest <- which(w[test.rows] != 0)
expect_true(all(is.na(ppd.t$value[, heldTest])))
expect_false(anyNA(ppd.t$value[, keptTest]))
fit.wt1 <- fit.wt
fit.wt1$weights.test <- ifelse(w[test.rows] == 0, 1, w[test.rows])
set.seed(4)
expect_identical(
  ppd.t$value[, keptTest],
  extract(fit.wt1, "ppd", sample = "test")[, keptTest]
)

# Student-t and aft fits hold a weight of 0 out the same way, a heteroscedastic
# gaussian fit at its own scale s(x)
fit.t <- suppressWarnings(bart(
  x,
  df$y,
  weights = w,
  family = "student",
  n.trees = 10L,
  n.burn = 10L,
  n.samples = 10L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
))
expect_identical(which(is.na(fitted(fit.t, type = "ppd"))), held)
fit.s <- suppressWarnings(bart(
  x,
  df$y,
  weights = w,
  variance = varianceForest(n.trees = 5L),
  n.trees = 10L,
  n.burn = 10L,
  n.samples = 10L,
  n.chains = 1L,
  n.threads = 1L,
  keepTrees = TRUE,
  verbose = FALSE
))
expect_identical(which(is.na(fitted(fit.s, type = "ppd"))), held)
expect_true(all(is.na(predict(fit.s, x[1:3, ], type = "ppd", weights = 0))))
expect_identical(
  which(
    colSums(is.na(predict(fit.aft, x[1:5, ], type = "ppd", weights = nd.w))) >
      0L
  ),
  c(1L, 4L)
)

# a binary fit's weights are trial counts and a count of 0 stays refused; the
# rows a probit fit's 0/1 weights mask out draw as they did
fit.mask <- suppressWarnings(bart(
  x,
  df$z,
  weights = as.numeric(w > 0),
  n.trees = 10L,
  n.burn = 10L,
  n.samples = 10L,
  n.chains = 1L,
  n.threads = 1L,
  keepTrees = TRUE,
  verbose = FALSE
))
expect_false(anyNA(extract(fit.mask, "ppd")))
expect_error(
  predict(fit.mask, x[1:3, ], type = "ppd", weights = c(0, 1, 1)),
  "trial counts and must be positive integers"
)
