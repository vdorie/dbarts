# A lone single-level factor leaves a root-only sampler with no available split
# variable anywhere. Birth/death must treat each tree's move as a no-op rather
# than force a birth and draw a rule for an invalid variable. Regression for the
# degenerate-root guard (src/bartcore/moves.hpp): unfixed, sampler$run segfaults
# here.

set.seed(0)
n <- 60L
y <- rnorm(n)

control <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 5L,
  updateState = FALSE
)
sampler <- dbarts(y ~ f, data.frame(f = factor(rep("a", n))), control = control)

res <- sampler$run(0L, 5L)

# the run completes with finite output instead of crashing
expect_true(all(is.finite(res$train)))
expect_true(all(is.finite(res$sigma)))

# every tree stays a lone root: one node per tree, all leaves (var == -1)
trees <- sampler$getTrees(current = TRUE)
expect_equal(nrow(trees), 5L)
expect_true(all(trees$var == -1L))

rm(sampler, res, trees)

# An ordinal column carries at least one cut point: with none, its own state
# validator refuses the store it sits in, and the summary printer indexes
# relative to a last cut that is not there. The entrance that sets a grid
# refuses an empty one rather than build that store.
x <- matrix(runif(n * 2L), n, 2L)
sampler <- dbarts(x, y, control = control)
expect_error(
  sampler$setCutPoints(list(numeric(0)), 1L),
  pattern = "at least one cut point"
)
expect_error(
  sampler$setCutPoints(list(c(0.25, 0.5), numeric(0)), 1:2),
  pattern = "at least one cut point"
)
# the grid the refusal left behind is the one the sampler still runs on
expect_true(all(is.finite(sampler$run(0L, 3L)$sigma)))
rm(sampler)

# A constant column induces no interior quantile, and the induced grid is
# floored to a single cut at the column's value, so the store stays
# restorable and the cutoff summary stays in bounds.
xConst <- cbind(rep(1.0, n), rnorm(n))
quantileControl <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 5L,
  useQuantiles = TRUE,
  updateState = FALSE
)
sampler <- dbarts(xConst, y, control = quantileControl)
sampler$storeState()
cutPoints <- attr(sampler$state, "cutPoints")
expect_equal(length(cutPoints[[1L]]), 1L)
expect_equal(cutPoints[[1L]], 1.0)
expect_true(length(cutPoints[[2L]]) > 1L)

# the state a zero-cut store could not round trip restores here
copied <- sampler$copy()
copied$storeState()
expect_equal(attr(copied$state, "cutPoints"), cutPoints)
expect_true(all(is.finite(copied$run(0L, 3L)$sigma)))
rm(sampler, copied, cutPoints)

# printing the cutoffs of that constant column completes, at every entrance
# that takes both a cutoff count and a quantile grid
verboseOutput <- capture.output(
  invisible(dbarts(
    xConst,
    y,
    control = dbartsControl(
      n.chains = 1L,
      n.threads = 1L,
      n.trees = 5L,
      useQuantiles = TRUE,
      printCutoffs = 10L,
      verbose = TRUE
    )
  ))
)
expect_true(any(grepl("x(1) cutoffs: 1.000000", verboseOutput, fixed = TRUE)))

invisible(capture.output(
  fit2 <- bart(
    xConst,
    y,
    n.samples = 5L,
    n.burn = 2L,
    n.chains = 1L,
    n.trees = 3L,
    n.threads = 1L,
    verbose = TRUE,
    printCutoffs = 5L,
    useQuantiles = TRUE
  )
))
expect_true(all(is.finite(fit2$yhat.train)))

invisible(capture.output(
  fit1 <- bart(
    xConst,
    y,
    ndpost = 5L,
    nskip = 2L,
    ntree = 3L,
    nthread = 1L,
    verbose = TRUE,
    printcutoffs = 5L,
    usequants = TRUE
  )
))
expect_true(all(is.finite(fit1$yhat.train)))

rm(verboseOutput, fit1, fit2, xConst, quantileControl, control, x, y, n)

# The uniform grid of a column whose range is one value repeats that value
# n.cuts times; a state carrying it restores, so a constant column (an all-zero
# dummy, a column made constant by subset) neither blocks a saved fit's predict
# nor the sampler's own copy and setState.
set.seed(1)
frame <- data.frame(y = rnorm(60L), a = runif(60L), b = 0, c = runif(60L))
frame$c[1:30] <- 0.5
fit <- bart(
  y ~ .,
  frame,
  subset = 1:30,
  n.samples = 5L,
  n.burn = 5L,
  n.chains = 2L,
  n.trees = 5L,
  n.threads = 1L,
  keepTrees = TRUE,
  verbose = FALSE
)
expected <- predict(fit, frame[1:3, ])
fit$fit$storeState()
reloaded <- unserialize(serialize(fit, NULL))
expect_identical(predict(reloaded, frame[1:3, ]), expected)

sampler <- dbarts(
  cbind(a = runif(40L), b = 0),
  rnorm(40L),
  control = dbartsControl(n.chains = 1L, n.threads = 1L, n.trees = 5L)
)
invisible(sampler$run(5L, 5L))
expect_true(inherits(sampler$copy(), "dbartsSampler"))
expect_silent(status <- sampler$setState(sampler$state))
expect_true(status)

# a column whose spread is below its own spacing collapses its uniform cuts
# too, and restores the same way
sampler <- dbarts(
  cbind(a = 1e15 + 4 * runif(40L)),
  rnorm(40L),
  control = dbartsControl(n.chains = 1L, n.threads = 1L, n.trees = 5L)
)
invisible(sampler$run(5L, 5L))
expect_silent(status <- sampler$setState(sampler$state))
expect_true(status)

# setCutPoints refuses a grid past the representable count, whose top bin
# would share the missing value's code, and a NaN cut
sampler <- dbarts(
  cbind(a = runif(40L)),
  rnorm(40L),
  control = dbartsControl(n.chains = 1L, n.threads = 1L, n.trees = 5L)
)
expect_error(
  sampler$setCutPoints(seq(0, 0.9, length.out = 65534L), 1L),
  "at most 65533 cut points"
)
expect_error(sampler$setCutPoints(c(0.2, NaN, 0.6), 1L), "NaN")
expect_error(sampler$setCutPoints(NaN, 1L), "NaN")

# The uniform grid spans a column's finite values, so one infinite value falls
# past an end cut instead of leaving the column unsplittable.
set.seed(2)
z <- runif(200L)
z[1L] <- Inf
y <- ifelse(z > 0.5, 2, 0) + rnorm(200L, sd = 0.2)
fit <- bart(
  data.frame(z = z),
  y,
  sigest = 0.2,
  n.samples = 50L,
  n.burn = 50L,
  n.chains = 1L,
  n.trees = 10L,
  n.threads = 1L,
  verbose = FALSE
)
expect_true(cor(fit$yhat.train.mean[-1L], y[-1L]) > 0.9)
# without a sigest, the starting estimate names the column it cannot fit
expect_error(
  bart(
    data.frame(z = z),
    y,
    n.samples = 5L,
    n.burn = 5L,
    n.chains = 1L,
    verbose = FALSE
  ),
  "predictor 'z' has infinite values"
)

# The quantile grid is over finite values too: a column holding only -Inf and
# Inf, or a few infinite values among finite ones, gives a grid whose state
# restores, so the saved fit predicts after reloading.
set.seed(3)
infiniteFrame <- data.frame(
  y = rnorm(40L),
  a = runif(40L),
  b = rep(c(-Inf, Inf), 20L),
  c = c(Inf, -Inf, runif(38L))
)
fit <- bart(
  y ~ .,
  infiniteFrame,
  useQuantiles = TRUE,
  sigest = 1,
  n.samples = 5L,
  n.burn = 5L,
  n.chains = 1L,
  n.trees = 5L,
  n.threads = 1L,
  keepTrees = TRUE,
  verbose = FALSE
)
expect_true(all(is.finite(unlist(attr(fit$fit$state, "cutPoints")))))
expected <- predict(fit, infiniteFrame[1:3, ])
fit$fit$storeState()
reloaded <- unserialize(serialize(fit, NULL))
expect_identical(predict(reloaded, infiniteFrame[1:3, ]), expected)
rm(infiniteFrame)

rm(frame, fit, expected, reloaded, sampler, z, y)
