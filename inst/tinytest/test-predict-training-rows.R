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
