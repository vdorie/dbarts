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
# floored to the single degenerate cut a uniform grid would place there, so
# the store stays restorable and the cutoff summary stays in bounds.
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
  fit2 <- bart2(
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
