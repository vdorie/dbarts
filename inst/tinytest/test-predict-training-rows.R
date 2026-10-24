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
# a missing offset is its row's na.action to answer, as at newdata: the default
# refuses, with the hint (the section at the end pins the rest)
expect_error(
  predict(fit.bare, offset = c(NA, rep(1, n - 1L))),
  "'offset' has missing values; use na.action = na.pass to return NA for them, or na.omit to drop them",
  fixed = TRUE
)

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
# checked as at newdata: one value or one per training row, a missing one being
# its row's na.action to answer (the section at the end), and a posterior
# predictive draw alone looking at the values
expect_error(predict(fit, weights = 1:3), "one per training row \\(50\\)")
expect_error(
  predict(fit.bare, weights = c(NA, rep(1, n - 1L))),
  "'weights' has missing values; use na.action = na.pass to return NA for them, or na.omit to drop them",
  fixed = TRUE
)
expect_error(
  predict(fit, weights = NA_real_),
  "'weights' has missing values",
  fixed = TRUE
)
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
# a column of weights is a vector of them, with no newdata and at it
set.seed(37)
a <- predict(fit, weights = matrix(w.new, n, 1L), type = "ppd")
set.seed(37)
expect_identical(a, predict(fit, weights = w.new, type = "ppd"))
set.seed(37)
a <- predict(fit, df, weights = matrix(w.new, n, 1L), type = "ppd")
set.seed(37)
expect_identical(a, predict(fit, df, weights = w.new, type = "ppd"))
# and a matrix of several columns is refused on both paths, where newdata had
# read the first column and let the others go
for (type in c("ev", "bart", "ppd")) {
  for (wide in list(
    matrix(1, n, 2L),
    matrix(1, 1L, 2L),
    matrix(w.new, n, 2L),
    array(1, c(n, 1L, 2L))
  )) {
    expect_error(
      predict(fit, weights = wide, type = type),
      "'weights' must be one number or one per training row",
      fixed = TRUE
    )
    expect_error(
      predict(fit, df, weights = wide, type = type),
      "'weights' must be one number or one per row of 'newdata'",
      fixed = TRUE
    )
  }
}
# so is a data frame or a list, which newdata had ignored or failed to coerce
for (type in c("ev", "bart", "ppd")) {
  for (notAtomic in list(data.frame(w = w.new), list(w.new), as.list(w.new))) {
    expect_error(
      predict(fit, weights = notAtomic, type = type),
      "'weights' must be one number or one per training row",
      fixed = TRUE
    )
    expect_error(
      predict(fit, df, weights = notAtomic, type = type),
      "'weights' must be one number or one per row of 'newdata'",
      fixed = TRUE
    )
  }
}
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
  short = 1:3,
  wide = matrix(1, n, 2L)
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

# a missing offset or weight at the training rows is its row's na.action to
# answer, exactly as at newdata (dec-B482): the default and na.fail refuse in
# newdata's words, na.pass and na.exclude return the row with NA draws, na.omit
# drops it. The row is computed at a stand-in and set aside afterward, so the
# random stream and every other row's draw are what a call without the missing
# value draws
na.w <- runif(n, 0.5, 3)
na.w[c(3L, 7L)] <- NA
na.o <- rnorm(n, 0, 0.5)
na.o[c(2L, 7L)] <- NA
na.rows <- list(weights = c(3L, 7L), offset = c(2L, 7L), both = c(2L, 3L, 7L))
na.given <- list(
  weights = list(weights = na.w),
  offset = list(offset = na.o),
  both = list(weights = na.w, offset = na.o)
)
seeded <- function(expr, seed = 41L) {
  set.seed(seed)
  expr
}
messageOf <- function(expr) {
  tryCatch(
    {
      force(expr)
      NA_character_
    },
    error = conditionMessage
  )
}
callOn <- function(object, args, newdata = NULL, ...) {
  do.call(
    predict,
    c(list(object), if (!is.null(newdata)) list(newdata), args, list(...))
  )
}
actions <- list(pass = na.pass, omit = na.omit, exclude = na.exclude)
# the same fit without its saved trees, which the training rows do not need
fit.stripped <- fit
fit.stripped[["fit"]] <- NULL

# the default refuses, and so does na.fail, in newdata's words, every type
for (type in c("ev", "bart", "ppd")) {
  for (label in names(na.given)) {
    for (action in list(NULL, na.fail)) {
      info <- paste(type, label, if (is.null(action)) "default" else "fail")
      here <- messageOf(callOn(
        fit,
        na.given[[label]],
        type = type,
        na.action = action
      ))
      there <- messageOf(callOn(
        fit,
        na.given[[label]],
        df,
        type = type,
        na.action = action
      ))
      expect_false(is.na(here), info = info)
      expect_identical(here, there, info = info)
    }
  }
}
expect_identical(
  messageOf(predict(fit, weights = na.w)),
  "'weights' has missing values; use na.action = na.pass to return NA for them, or na.omit to drop them"
)
expect_identical(
  messageOf(predict(fit, offset = na.o)),
  "'offset' has missing values; use na.action = na.pass to return NA for them, or na.omit to drop them"
)
expect_identical(
  messageOf(predict(fit, weights = na.w, na.action = na.fail)),
  "'weights' has missing values, which na.action = na.fail refuses"
)
expect_identical(
  messageOf(predict(
    fit,
    weights = na.w,
    offset = na.o,
    na.action = na.fail
  )),
  "'offset' and 'weights' have missing values, which na.action = na.fail refuses"
)
# the form of the one given is still checked, in words that do not say none is
# missing, which na.action may answer
expect_identical(
  messageOf(predict(fit, offset = 1:3)),
  "'offset' must be one number or one per training row (50)"
)
expect_identical(
  messageOf(predict(fit, weights = 1:3)),
  "'weights' must be one number or one per training row (50)"
)
# a lone NA has no row for na.action to answer, whatever it is
for (action in c(actions, list(default = NULL))) {
  for (given in list(list(weights = NA_real_), list(offset = NA_real_))) {
    here <- messageOf(callOn(fit, given, na.action = action))
    expect_identical(
      here,
      messageOf(callOn(fit, given, df, na.action = action))
    )
    expect_identical(here, paste0("'", names(given), "' has missing values"))
  }
}

# na.pass, na.omit and na.exclude: the numbers of the call at newdata, on the
# stored draws or on saved trees, with the rows set aside as newdata sets them
# aside
for (name in names(actions)) {
  action <- actions[[name]]
  for (label in names(na.given)) {
    aside <- na.rows[[label]]
    for (type in c("ev", "bart")) {
      info <- paste(name, label, type)
      there <- callOn(
        fit,
        na.given[[label]],
        df,
        type = type,
        na.action = action
      )
      expect_equal(
        callOn(fit, na.given[[label]], type = type, na.action = action),
        there,
        info = info
      )
      # no trees are needed
      expect_equal(
        callOn(
          fit.stripped,
          na.given[[label]],
          type = type,
          na.action = action
        ),
        there,
        info = info
      )
      expect_equal(
        callOn(
          fit.stripped,
          na.given[[label]],
          type = type,
          na.action = action,
          combineChains = FALSE
        ),
        callOn(
          fit,
          na.given[[label]],
          df,
          type = type,
          na.action = action,
          combineChains = FALSE
        ),
        info = info
      )
      # and the band
      expect_equal(
        callOn(
          fit.stripped,
          na.given[[label]],
          type = type,
          na.action = action,
          ci.level = 0.9
        ),
        callOn(
          fit,
          na.given[[label]],
          df,
          type = type,
          na.action = action,
          ci.level = 0.9
        ),
        info = info
      )
    }
    # the rows are set aside by NA, or cut
    got <- callOn(fit, na.given[[label]], na.action = action)
    expect_identical(
      dim(got),
      c(20L, n - if (name == "omit") length(aside) else 0L),
      info = paste(name, label)
    )
    if (name != "omit") {
      expect_identical(
        unname(which(colSums(is.na(got)) > 0L)),
        aside,
        info = name
      )
    } else {
      expect_false(anyNA(got), info = name)
    }
    # a posterior predictive draw reads and draws only the rows na.action
    # keeps, as at newdata, so under one seed it is the draw of the call at
    # newdata with the same action, and with the rows left out it is the call at
    # the kept rows of newdata
    given <- na.given[[label]]
    ppd <- seeded(callOn(fit, given, type = "ppd", na.action = action))
    expect_equal(
      ppd,
      seeded(callOn(fit, given, df, type = "ppd", na.action = action)),
      info = paste(name, label)
    )
    kept <- setdiff(seq_len(n), aside)
    atKept <- seeded(callOn(
      fit,
      lapply(given, function(value) value[kept]),
      df[kept, ],
      type = "ppd"
    ))
    if (name == "omit") {
      expect_equal(ppd, atKept, info = paste(name, label))
    } else {
      expect_equal(
        unname(ppd[, kept]),
        unname(atKept),
        info = paste(name, label)
      )
      expect_true(all(is.na(ppd[, aside])), info = paste(name, label))
    }
  }
}

# every row set aside: na.pass returns NA for all of them and na.omit none, a
# custom function is followed, and the band is NA over the rows set aside
all.w <- rep(NA_real_, n)
expect_true(all(is.na(predict(fit, weights = all.w, na.action = na.pass))))
expect_identical(
  dim(predict(fit, weights = all.w, na.action = na.pass)),
  dim(predict(fit, df, weights = all.w, na.action = na.pass))
)
expect_identical(
  dim(predict(fit, weights = all.w, na.action = na.omit)),
  dim(predict(fit, df, weights = all.w, na.action = na.omit))
)
expect_identical(
  dim(predict(fit, weights = all.w, na.action = na.omit)),
  c(20L, 0L)
)
# and with none kept the one row that stands in leaves R's random stream as it
# was found
for (action in list(na.pass, na.omit, na.exclude)) {
  set.seed(5L)
  stream <- .Random.seed
  invisible(predict(fit, weights = all.w, type = "ppd", na.action = action))
  expect_identical(.Random.seed, stream)
}
dropper <- function(object, ...) stats::na.omit(object, ...)
expect_identical(
  predict(fit, offset = na.o, na.action = dropper),
  predict(fit, offset = na.o, na.action = na.omit)
)
expect_identical(
  predict(fit, offset = na.o, na.action = "na.omit"),
  predict(fit, offset = na.o, na.action = na.omit)
)
band <- predict(fit, offset = na.o, ci.level = 0.9, na.action = na.pass)
expect_identical(unname(which(is.na(band[, "est"]))), na.rows$offset)
expect_identical(unname(which(is.na(band[, "ci.lower"]))), na.rows$offset)
expect_identical(unname(which(is.na(band[, "ci.upper"]))), na.rows$offset)
expect_identical(
  dim(predict(fit, offset = na.o, ci.level = 0.9, na.action = na.omit)),
  c(n - 2L, 3L)
)
# na.action is not looked at where nothing is missing, as predict.lm with no
# newdata does not look at it (at newdata it is resolved on every call): a call
# with a complete offset and weights is what it was, and a name that is no
# function is not asked about
expect_identical(
  predict(fit, offset = new.off, na.action = na.omit),
  predict(fit, offset = new.off)
)
expect_identical(predict(fit, na.action = na.omit), predict(fit))
expect_identical(predict(fit.h, na.action = na.omit), predict(fit.h))
expect_identical(predict(fit, na.action = "nonexistent"), predict(fit))
expect_identical(
  predict(fit, offset = new.off, weights = 2, na.action = "nonexistent"),
  predict(fit, offset = new.off, weights = 2)
)
expect_error(
  predict(fit, df, na.action = "nonexistent"),
  "nonexistent"
)
# and it is asked about once a value is missing
expect_error(
  predict(fit, offset = na.o, na.action = "nonexistent"),
  "nonexistent"
)
# a malformed weight is refused as malformed whatever na.action says of the
# missing value in it, in the words of the training rows
for (malformed in list(
  data.frame(w = c(NA, rep(1, n - 1L))),
  list(NA),
  matrix(c(NA, rep(1, 2L * n - 1L)), n, 2L)
)) {
  for (action in c(actions, list(default = NULL, fail = na.fail))) {
    expect_identical(
      messageOf(callOn(fit, list(weights = malformed), na.action = action)),
      "'weights' must be one number or one per training row (50)"
    )
  }
}
# a weight of 0 and a missing weight are different rows: the first has NA
# draws for want of a noise scale, the second is na.action's to answer
w.zero.na <- rep(1, n)
w.zero.na[4L] <- 0
w.zero.na[9L] <- NA
mixed <- predict(fit, weights = w.zero.na, type = "ppd", na.action = na.omit)
expect_identical(dim(mixed), c(20L, n - 1L))
expect_identical(unname(which(colSums(is.na(mixed)) > 0L)), 4L)
mixed <- predict(fit, weights = w.zero.na, type = "ppd", na.action = na.pass)
expect_identical(unname(which(colSums(is.na(mixed)) > 0L)), c(4L, 9L))

# the other families, against the call at newdata: a binary fit's weights are
# read by no type but a posterior predictive draw, whose trial counts a missing
# value leaves to na.action
for (action in actions) {
  expect_equal(
    predict(fit.p, offset = na.o, na.action = action),
    predict(fit.p, df, offset = na.o, na.action = action)
  )
  expect_equal(
    predict(fit.p, weights = na.w, type = "bart", na.action = action),
    predict(fit.p, df, weights = na.w, type = "bart", na.action = action)
  )
}
counts.na <- n.trials
counts.na[c(4L, 9L)] <- NA
for (binary.fit in list(fit.p, fit.lg)) {
  for (action in list(na.pass, na.omit, na.exclude)) {
    got <- seeded(predict(
      binary.fit,
      weights = counts.na,
      type = "ppd",
      na.action = action
    ))
    # the draw of the call at newdata, which draws the kept rows alone
    expect_equal(
      got,
      seeded(predict(
        binary.fit,
        df,
        weights = counts.na,
        type = "ppd",
        na.action = action
      ))
    )
    if (identical(action, na.omit)) {
      expect_equal(
        got,
        seeded(predict(
          binary.fit,
          df[-c(4L, 9L), ],
          weights = n.trials[-c(4L, 9L)],
          type = "ppd"
        ))
      )
    } else {
      expect_identical(unname(which(colSums(is.na(got)) > 0L)), c(4L, 9L))
    }
  }
}

# a count fit, an ordered response and an unordered one
for (action in actions) {
  expect_equal(
    predict(fit.nb, offset = na.o, na.action = action),
    predict(fit.nb, df, offset = na.o, na.action = action)
  )
  expect_equal(
    predict(
      fit.nb,
      offset = na.o,
      type = "bart",
      ci.level = 0.9,
      na.action = action
    ),
    predict(
      fit.nb,
      df,
      offset = na.o,
      type = "bart",
      ci.level = 0.9,
      na.action = action
    )
  )
  expect_equal(
    predict(fit.ord, offset = na.o, na.action = action),
    predict(fit.ord, df, offset = na.o, na.action = action)
  )
  expect_equal(
    predict(fit.ord, offset = na.o, type = "class", na.action = action),
    predict(fit.ord, df, offset = na.o, type = "class", na.action = action)
  )
  expect_equal(
    predict(fit.ord, offset = na.o, ci.level = 0.9, na.action = action),
    predict(fit.ord, df, offset = na.o, ci.level = 0.9, na.action = action)
  )
  na.cat <- cat.off
  na.cat[5L, ] <- NA
  na.cat[9L, 2L] <- NA
  expect_equal(
    predict(fit.mn, offset = na.cat, na.action = action),
    predict(fit.mn, df, offset = na.cat, na.action = action),
    tolerance = 1e-6
  )
  expect_equal(
    predict(fit.mn, offset = na.cat, type = "class", na.action = action),
    predict(fit.mn, df, offset = na.cat, type = "class", na.action = action)
  )
  # a posterior predictive draw of each of them reads and draws the kept rows
  # alone, as at newdata
  for (case in list(
    list(fit.nb, na.o),
    list(fit.ord, na.o),
    list(fit.mn, na.cat)
  )) {
    expect_equal(
      seeded(predict(
        case[[1L]],
        offset = case[[2L]],
        type = "ppd",
        na.action = action
      )),
      seeded(predict(
        case[[1L]],
        df,
        offset = case[[2L]],
        type = "ppd",
        na.action = action
      ))
    )
  }
}
# a count fit's posterior predictive draw reads and draws the kept rows alone,
# the same draws under each action
nb.pass <- seeded(predict(
  fit.nb,
  offset = na.o,
  type = "ppd",
  na.action = na.pass
))
nb.omit <- seeded(predict(
  fit.nb,
  offset = na.o,
  type = "ppd",
  na.action = na.omit
))
expect_identical(dim(nb.omit), c(20L, n - 2L))
expect_identical(unname(which(colSums(is.na(nb.pass)) > 0L)), na.rows$offset)
expect_identical(unname(nb.pass[, -na.rows$offset]), unname(nb.omit))
expect_false(anyNA(nb.omit))
expect_true(all(nb.omit == round(nb.omit)))
expect_true(all(is.na(
  predict(fit.mn, offset = na.cat, type = "class", na.action = na.pass)[c(
    5L,
    9L
  )]
)))
expect_length(
  predict(fit.mn, offset = na.cat, type = "class", na.action = na.omit),
  n - 2L
)
for (refused.fit in list(fit.p, fit.nb, fit.ord, fit.mn)) {
  given <- if (identical(refused.fit, fit.mn)) na.cat else na.o
  expect_identical(
    messageOf(predict(refused.fit, offset = given)),
    messageOf(predict(refused.fit, df, offset = given))
  )
  expect_true(startsWith(
    messageOf(predict(refused.fit, offset = given)),
    "'offset' has missing values; use na.action = na.pass"
  ))
}
# a heteroscedastic, a Student-t and an aft fit draw a row's noise from its
# weight: the rows set aside have NA draws, as at newdata, and under one seed
# the rest are the draws of the call at newdata (a heteroscedastic call at
# newdata also carries the scale s(x) as an attribute, which the training rows
# never did)
for (family.fit in list(fit.stu, fit.het, fit.aft)) {
  for (action in list(na.pass, na.omit, na.exclude)) {
    got <- seeded(predict(
      family.fit,
      weights = na.w,
      type = "ppd",
      na.action = action
    ))
    there <- seeded(predict(
      family.fit,
      x,
      weights = na.w,
      type = "ppd",
      na.action = action
    ))
    expect_identical(dim(got), dim(there))
    expect_equal(got, there, check.attributes = FALSE)
    expect_identical(
      unname(which(colSums(is.na(got)) > 0L)),
      unname(which(colSums(is.na(there)) > 0L))
    )
  }
}

# a multinomial fit on a count matrix draws, at its training rows, a vector of
# counts with each row's own trials: the rows kept are drawn with their own
# trials and the rows set aside are not drawn
trials.n <- rep_len(2:6, n)
counts.n <- t(vapply(
  seq_len(n),
  function(i) rmultinom(1L, trials.n[i], c(0.5, 0.3, 0.2))[, 1L],
  integer(3L)
))
fit.counts <- quick(x, counts.n, family = "multinomial")
off.counts <- matrix(0, n, 3L)
off.counts[c(4L, 9L), 2L] <- NA
omitted <- seeded(predict(
  fit.counts,
  offset = off.counts,
  type = "ppd",
  na.action = na.omit
))
expect_identical(dim(omitted)[2L], n - 2L)
expect_identical(
  unname(apply(omitted, c(1L, 2L), sum)),
  matrix(trials.n[-c(4L, 9L)], dim(omitted)[1L], n - 2L, byrow = TRUE)
)
passed <- seeded(predict(
  fit.counts,
  offset = off.counts,
  type = "ppd",
  na.action = na.pass
))
expect_identical(unname(which(apply(is.na(passed), 2L, any))), c(4L, 9L))
expect_identical(unname(passed[, -c(4L, 9L), ]), unname(omitted))
# its credible band keeps the category margin as the draws do, one est,
# ci.lower and ci.upper per row and category, whether or not a row is set
# aside; the mean counts of a row add up to its trials
band.counts <- predict(fit.counts, type = "ppd", ci.level = 0.9)
expect_identical(dim(band.counts), c(n, 3L, 3L))
expect_identical(
  dimnames(band.counts),
  dimnames(predict(fit.counts, type = "ev", ci.level = 0.9))
)
expect_equal(unname(rowSums(band.counts[,, 1L])), as.numeric(trials.n))
# a labeled fit draws one category a row, so its band has no category margin
expect_identical(dim(predict(fit.mn, type = "ppd", ci.level = 0.9)), c(n, 3L))
for (action in list(na.pass, na.exclude, na.omit)) {
  band.set <- seeded(predict(
    fit.counts,
    offset = off.counts,
    type = "ppd",
    ci.level = 0.9,
    na.action = action
  ))
  if (identical(action, na.omit)) {
    expect_identical(dim(band.set), c(n - 2L, 3L, 3L))
    expect_false(anyNA(band.set))
  } else {
    expect_identical(dim(band.set), c(n, 3L, 3L))
    expect_identical(
      unname(which(apply(is.na(band.set), 1L, any))),
      c(4L, 9L)
    )
    expect_equal(
      unname(rowSums(band.set[-c(4L, 9L), , 1L])),
      as.numeric(trials.n[-c(4L, 9L)])
    )
  }
}

# the case the review found: 40 rows, a missing weight at rows 3 and 7, one seed
df.40 <- df[1:40, c("a", "b", "c", "y")]
fit.40 <- quick(y ~ a + b + c, df.40)
w.40 <- runif(40L, 0.5, 3)
w.40[c(3L, 7L)] <- NA
ok.40 <- !is.na(w.40)
at.kept <- seeded(
  predict(
    fit.40,
    df.40[ok.40, ],
    weights = w.40[ok.40],
    type = "ppd"
  ),
  10L
)
expect_equal(
  seeded(
    predict(fit.40, weights = w.40, type = "ppd", na.action = na.omit),
    10L
  ),
  at.kept
)
expect_equal(
  seeded(
    predict(fit.40, df.40, weights = w.40, type = "ppd", na.action = na.omit),
    10L
  ),
  at.kept
)
for (action in list(na.pass, na.exclude)) {
  expect_equal(
    seeded(
      predict(fit.40, weights = w.40, type = "ppd", na.action = action),
      10L
    ),
    seeded(
      predict(fit.40, df.40, weights = w.40, type = "ppd", na.action = action),
      10L
    )
  )
  expect_equal(
    unname(seeded(
      predict(fit.40, weights = w.40, type = "ppd", na.action = action),
      10L
    )[, ok.40]),
    unname(at.kept)
  )
}

# the draws, the band and the classes of every family with every row set aside
# are NA for all of them or for none, shaped as when some are kept and leaving
# the random stream as it was; a call that read no row at all could not be
# shaped, which is why one stands in
all.o <- rep(NA_real_, n)
all.cat <- matrix(NA_real_, n, 3L)
categoryFits <- list(fit.ord, fit.mn, fit.counts)
for (case in list(
  list(fit, all.o),
  list(fit.p, all.o),
  list(fit.nb, all.o),
  list(fit.ord, all.o),
  list(fit.mn, all.cat),
  list(fit.counts, all.cat)
)) {
  object <- case[[1L]]
  for (type in c("ev", "ppd")) {
    for (ci in list(NULL, 0.9)) {
      info <- paste(class(object)[1L], type, if (!is.null(ci)) "band")
      set.seed(5L)
      stream <- .Random.seed
      passed <- predict(
        object,
        offset = case[[2L]],
        type = type,
        ci.level = ci,
        na.action = na.pass
      )
      omitted <- predict(
        object,
        offset = case[[2L]],
        type = type,
        ci.level = ci,
        na.action = na.omit
      )
      expect_true(all(is.na(passed)), info = info)
      expect_identical(.Random.seed, stream, info = info)
      # the observation margin: the first of a band, else the last, or the one
      # before the categories of the probabilities and of a count draw
      d <- dim(passed)
      margin <- if (!is.null(ci)) {
        1L
      } else if (
        any(vapply(categoryFits, identical, NA, object)) &&
          (type == "ev" || identical(object, fit.counts))
      ) {
        length(d) - 1L
      } else {
        length(d)
      }
      expect_identical(d[margin], n, info = info)
      expect_identical(dim(omitted)[margin], 0L, info = info)
    }
  }
}
for (object in list(fit.ord, fit.mn)) {
  given <- if (identical(object, fit.mn)) all.cat else all.o
  expect_true(all(is.na(
    predict(object, offset = given, type = "class", na.action = na.pass)
  )))
  expect_length(
    predict(object, offset = given, type = "class", na.action = na.omit),
    0L
  )
}
