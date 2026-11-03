# xbart() uses a grid its data object carries, as dbarts() does; a stated n.cuts
# or useQuantiles wins and derives it again.

set.seed(81)
n <- 200L
frame <- data.frame(x1 = runif(n), x2 = rexp(n))
frame$y <- 3 * frame$x1 + frame$x2 + rnorm(n, sd = 0.3)
run <- function() {
  dbartsControl(n.samples = 2L, n.burn = 0L, n.chains = 1L, n.threads = 1L)
}
form <- y ~ x1 + x2
xb <- function(object, ..., data = NULL, loss = "rmse") {
  # a data argument beside a dbartsData is ignored with a warning, so it goes
  # only where there is no dbartsData to carry the data
  source <- if (is.null(data)) list(object) else list(object, data)
  do.call(
    xbart,
    c(
      source,
      list(
        n.trees = 25L,
        ...,
        loss = loss,
        n.reps = 1L,
        n.samples = 15L,
        n.burn = c(15L, 5L),
        method = "k-fold",
        n.test = 3L,
        n.threads = 1L,
        seed = 1L,
        verbose = FALSE
      )
    )
  )
}
fitOf <- function(...) {
  sampler <- dbarts(form, frame, ..., control = run())
  invisible(sampler$run(0L, 2L))
  sampler
}

# a fit's own object with its predictors as they were gives the number the
# formula gives, to the bit, whatever made its grid
fit25 <- fitOf(n.cuts = 25L)
expect_identical(xb(fit25$data), xb(form, data = frame, n.cuts = 25L))
fitQuantile <- fitOf(n.cuts = 25L, useQuantiles = TRUE)
expect_identical(
  xb(fitQuantile$data),
  xb(form, data = frame, n.cuts = 25L, useQuantiles = TRUE)
)
fitPerColumn <- fitOf(n.cuts = c(x2 = 5L, x1 = 30L))
expect_identical(
  xb(fitPerColumn$data),
  xb(form, data = frame, n.cuts = c(30L, 5L))
)
# a column of few values, whose points carry weights
fewFrame <- frame
fewFrame$x2 <- sample(c(0, 1, 5), n, TRUE)
fitFew <- dbarts(form, fewFrame, n.cuts = 25L, control = run())
invisible(fitFew$run(0L, 2L))
expect_false(is.null(attr(fitFew$data@cutPoints[[2L]], "mass")))
expect_identical(xb(fitFew$data), xb(form, data = fewFrame, n.cuts = 25L))
# factor columns, ordered and not
factorFrame <- frame
factorFrame$x2 <- factor(sample(c("a", "b", "c"), n, TRUE))
factorFrame$o <- factor(
  sample(c("low", "mid", "high"), n, TRUE),
  ordered = TRUE
)
factorForm <- y ~ x1 + x2 + o
fitFactors <- dbarts(factorForm, factorFrame, n.cuts = 25L, control = run())
invisible(fitFactors$run(0L, 2L))
expect_identical(
  xb(fitFactors$data),
  xb(factorForm, data = factorFrame, n.cuts = 25L)
)

# a fold splits on the grid carried: at one point per column a draw's fits take
# at most four values over the held-out rows, as they do not on the grid derived
distinctFits <- function(y.test, y.test.hat, weights) {
  max(apply(y.test.hat, 2L, function(column) length(unique(column))))
}
handGrid <- function(...) {
  object <- fit25$data
  object@cutPoints <- list(...)
  object
}
atOnePoint <- handGrid(0.5, 1)
expect_true(all(xb(atOnePoint, loss = distinctFits) <= 4))
expect_true(any(xb(fit25$data, loss = distinctFits) > 4))
expect_identical(xb(atOnePoint), xb(handGrid(0.5, 1)))
expect_false(identical(xb(atOnePoint), xb(fit25$data)))
# on three points and two, twelve
atSome <- handGrid(c(0.2, 0.5, 0.8), c(0.5, 1.5))
expect_true(all(xb(atSome, loss = distinctFits) <= 12))
expect_true(any(xb(atSome, loss = distinctFits) > 4))
# dbarts() takes the same grid
sampler <- dbarts(atSome, control = run())
invisible(sampler$run(0L, 2L))
expect_identical(lengths(sampler$data@cutPoints), c(3L, 2L))

# a stated count or rule wins and derives again
expect_identical(
  xb(atOnePoint, n.cuts = 25L),
  xb(form, data = frame, n.cuts = 25L)
)
expect_identical(
  xb(atOnePoint, useQuantiles = TRUE),
  xb(form, data = frame, n.cuts = 25L, useQuantiles = TRUE)
)

# edited predictors stay on the carried grid, as in dbarts()
X <- as.matrix(frame[c("x1", "x2")])
fitMatrix <- dbarts(X, frame$y, n.cuts = 25L, control = run())
invisible(fitMatrix$run(0L, 2L))
edited <- fitMatrix$data
edited@x[, 1L] <- edited@x[, 1L] * 10
X[, 1L] <- X[, 1L] * 10
expect_false(identical(xb(edited), xb(dbartsData(X, frame$y, n.cuts = 25L))))

# a record is held as dbarts() holds it: put in order, its form refused by name
expect_identical(
  xb(handGrid(c(0.8, 0.2, 0.5), c(1.5, 0.5))),
  xb(atSome)
)
expect_error(
  xb(handGrid(0.5)),
  "'cutPoints' record must be a list with one entry"
)
expect_error(xb(handGrid(0.5, c(1, 1))), "repeats a point")
ownOrdered <- fitFactors$data
expect_identical(
  xb(ownOrdered),
  xb(factorForm, data = factorFrame, n.cuts = 25L)
)
ownOrdered@cutPoints[[3L]] <- c(0.25, 1.25)
expect_error(
  xb(ownOrdered),
  "an ordered factor, other cut points than those of its levels"
)
