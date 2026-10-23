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

# weights given with no newdata replace the fit's own for the training rows'
# predictive draws, one value or one per training row, as they act at newdata
# for each family and type, worked from the stored draws with no trees
# (predict.lm with no newdata)
w.new <- runif(n, 0.5, 4)
for (given in list(w.new, 2.5)) {
  set.seed(31)
  a <- predict(fit, weights = given, type = "ppd")
  set.seed(31)
  b <- predict(fit, df, weights = given, type = "ppd")
  expect_equal(a, b, tolerance = 1e-8)
  # a fit kept without trees draws from its stored draws alone
  set.seed(31)
  a <- predict(fit.bare, weights = given, type = "ppd", combineChains = FALSE)
  set.seed(31)
  b <- dbarts:::sampleFromPPD(
    extract(fit.bare, "ev", combineChains = FALSE),
    fit.bare,
    rep_len(given, n),
    2L
  )
  expect_identical(a, b)
}
# the band over the replaced draws
set.seed(32)
a <- predict(fit, weights = w.new, type = "ppd", ci.level = 0.9)
set.seed(32)
b <- predict(fit, df, weights = w.new, type = "ppd", ci.level = 0.9)
expect_equal(a, b, tolerance = 1e-8)
# the weights only reach a posterior predictive draw, as at newdata
for (type in c("ev", "bart")) {
  expect_identical(
    predict(fit, weights = w.new, type = type),
    predict(fit, type = type)
  )
}
# they set the noise: a larger weight is a smaller spread
set.seed(33)
spread <- function(weights) {
  sd(predict(fit, weights = weights, type = "ppd") - predict(fit, type = "ev"))
}
expect_true(spread(100) < 0.5 * spread(1))
# the fit's own weights given back change nothing, and the replacement does not
# touch the fit
w.own <- rep(c(1, 4), length.out = n)
fit.own <- quick(x, df$y, weights = w.own)
set.seed(34)
a <- predict(fit.own, type = "ppd")
set.seed(34)
expect_identical(predict(fit.own, weights = w.own, type = "ppd"), a)
expect_identical(fit.own$weights, w.own)
# a student-t fit pairs each draw's own degrees of freedom, an aft fit draws
# on the log-time scale, a heteroscedastic fit at its own scale
fit.stu <- quick(x, df$y, family = "student")
fit.het <- quick(x, df$y, variance = varianceForest(n.trees = 5L))
for (family.fit in list(fit.stu, fit.het, fit.aft)) {
  set.seed(35)
  a <- predict(family.fit, weights = w.new, type = "ppd")
  set.seed(35)
  b <- predict(family.fit, x, weights = w.new, type = "ppd")
  expect_equal(a, b, tolerance = 1e-8, check.attributes = FALSE)
}
# trial counts on a binary fit
n.trials <- rep(c(1, 3, 5), length.out = n)
fit.lg <- quick(z ~ a + b + c, df, family = "logistic")
for (binary.fit in list(fit.p, fit.lg)) {
  set.seed(36)
  a <- predict(binary.fit, weights = n.trials, type = "ppd")
  set.seed(36)
  b <- predict(binary.fit, df, weights = n.trials, type = "ppd")
  expect_equal(a, b)
  expect_true(all(a >= 0 & a <= rep(n.trials, each = nrow(a))))
  expect_error(
    predict(binary.fit, weights = 0, type = "ppd"),
    "trial counts and must be positive integers"
  )
  expect_error(
    predict(binary.fit, weights = 1.5, type = "ppd"),
    "trial counts and must be positive integers"
  )
}
# checked as at newdata: one value or one per training row, none missing (the
# offset's rule, na.action having no say among the fit's own rows), and a
# posterior predictive draw alone looking at the values
expect_error(predict(fit, weights = 1:3), "one per training row \\(50\\)")
expect_error(predict(fit.bare, weights = c(NA, rep(1, n - 1L))), "none missing")
expect_error(predict(fit, weights = NA_real_), "none missing")
expect_error(
  predict(fit, weights = matrix(1, n, 2L), type = "ppd"),
  "one per training row"
)
expect_error(
  predict(fit, weights = c(-1, rep(1, n - 1L)), type = "ppd"),
  "finite and non-negative"
)
expect_error(
  predict(fit, weights = c(Inf, rep(1, n - 1L)), type = "ppd"),
  "finite and non-negative"
)
# a column of weights is a vector of them, as at newdata
set.seed(37)
a <- predict(fit, weights = matrix(w.new, n, 1L), type = "ppd")
set.seed(37)
expect_identical(a, predict(fit, weights = w.new, type = "ppd"))
# the two paths agree on what each type does with a weight that is negative,
# not a number, missing or the wrong length: a type that draws nothing never
# looks at the values, a posterior predictive draw refuses a negative or
# non-numeric one in the family's words, and a missing or misplaced one is
# refused on both
refused <- function(expr) {
  tryCatch(
    {
      force(expr)
      NA_character_
    },
    error = conditionMessage
  )
}
bad.weights <- list(
  negative = -1,
  negative.row = c(-1, rep(1, n - 1L)),
  text = "a",
  text.rows = rep("a", n),
  missing = NA_real_,
  missing.row = c(NA, rep(1, n - 1L)),
  short = 1:3
)
for (type in c("ev", "bart", "ppd")) {
  for (label in names(bad.weights)) {
    given <- bad.weights[[label]]
    here <- refused(predict(fit, weights = given, type = type))
    there <- refused(predict(fit, df, weights = given, type = type))
    info <- paste(type, label)
    expect_identical(is.na(here), is.na(there), info = info)
    if (label %in% c("negative", "negative.row", "text", "text.rows")) {
      # the family's check at a posterior predictive draw, in the same words
      expect_identical(
        here,
        if (type == "ppd") there else NA_character_,
        info = info
      )
    }
  }
}
for (type in c("ev", "bart")) {
  # weights that no draw reads change nothing, whatever they hold
  for (given in list(-1, "a", rep("a", n), c(-1, rep(1, n - 1L)))) {
    expect_identical(
      tryCatch(
        predict(fit, weights = given, type = type),
        error = conditionMessage
      ),
      predict(fit, type = type)
    )
  }
}
for (given in bad.weights) {
  # a heteroscedastic fit's scale takes no weight, whatever it holds
  expect_error(
    predict(fit.het, weights = given, type = "sigma"),
    "does not support 'weights'",
    info = "sigma, no newdata"
  )
  expect_error(
    predict(fit.het, x, weights = given, type = "sigma"),
    "does not support 'weights'",
    info = "sigma, newdata"
  )
}
# the other families take no weight in a posterior predictive draw, as at newdata
expect_error(predict(fit.nb, weights = 1), "takes no per-observation weight")
expect_error(predict(fit.ord, weights = 1), "takes no per-observation weight")

# an offset, weights or bases replace the fit's own on the stored training
# draws, so a fit kept without them is refused by name, the same for all three
fit.nofits <- quick(x, df$y, keepFits = FALSE)
fit.notrain <- quick(x, df$y, keepTrainingFits = FALSE)
for (dropped in list(fit.nofits, fit.notrain)) {
  for (request in list(
    list(offset = 1),
    list(offset = new.off),
    list(offset = 1:3),
    list(weights = 1, type = "ppd"),
    list(weights = w.new, type = "ppd")
  )) {
    expect_error(
      do.call(predict, c(list(dropped), request)),
      "cannot extract train sample predictions; bart must be called with 'keepTrainingFits' == TRUE and 'keepFits' == TRUE",
      fixed = TRUE,
      info = paste(names(request), collapse = ", ")
    )
  }
}
fit.bt <- suppressWarnings(bartBT(
  x,
  df$y,
  ntree = 10L,
  nskip = 10L,
  ndpost = 10L,
  keeptrainfits = FALSE,
  verbose = FALSE
))
expect_error(
  predict(fit.bt, offset = 1),
  "'keeptrainfits' == TRUE",
  fixed = TRUE
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
# and with no newdata: the fit's own weights give NA at the held rows, weights
# stated for them draw them as new observations, and weights of 0 given give NA
expect_identical(
  which(colSums(is.na(predict(fit.w, type = "ppd"))) > 0L),
  held
)
expect_false(anyNA(predict(fit.w, weights = 1, type = "ppd")))
expect_identical(
  which(
    colSums(is.na(predict(
      fit.w,
      weights = ifelse(w == 0, 1, 0),
      type = "ppd"
    ))) >
      0L
  ),
  kept
)
expect_true(all(is.na(predict(fit.w, weights = 0, type = "ppd"))))
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

# a logistic fit's test-frame weight of 0 is a count of no trials, whose draw
# is 0 successes and finite, where a precision weight of 0 has NA draws
count.train <- data.frame(x, z = df$z, count = rep(c(1, 2, 3), length.out = n))
count.test <- count.train[1:6, ]
count.test$count <- c(0, 1, 1, 0, 2, 1)
fit.count <- suppressWarnings(bart(
  z ~ a + b + c,
  count.train,
  test = count.test,
  weights = count,
  family = "logistic",
  n.trees = 10L,
  n.burn = 10L,
  n.samples = 10L,
  n.chains = 2L,
  n.threads = 1L,
  keepTrees = TRUE,
  verbose = FALSE
))
expect_identical(fit.count$weights.test, count.test$count)
ppd.count <- extract(fit.count, "ppd", sample = "test")
expect_false(anyNA(ppd.count))
expect_true(all(ppd.count[, c(1L, 4L)] == 0))
