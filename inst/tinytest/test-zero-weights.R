# Zero weights drop an observation from the leaf likelihood (the equivalence
# gate exercises the no-likelihood path bitwise against itself; this pins it
# against an oracle). The test: perturbing a zero-weight row's response must
# leave the fit at every other row bitwise unchanged, since that row
# contributes nothing to any leaf's sufficient statistics. The perturbed row's
# own reported fit is only equal up to rounding: the run loop keeps a running
# residual y - totalFits and rebuilds totalFits from it once per sweep, so the
# row's y enters and leaves its fit with one rounding apiece. The gaussian
# response transform is range-based (fitMin/fitMax), so the perturbation is
# kept strictly inside the response range - with the extremes held by
# non-zero-weight rows - to avoid moving the transform, which would
# legitimately shift every fit.

set.seed(20)
n <- 200L
x <- matrix(runif(n * 3L), n, 3L)
colnames(x) <- c("a", "b", "c")
w <- rep(1.0, n)
zeros <- sample(n, 40L)
w[zeros] <- 0.0
nz <- which(w > 0.0)

signal <- 2 * x[, 1L] - x[, 2L]
# non-zero-weight rows span [0, 10] exactly (so they hold the extremes)
ybase <- 10 * (signal - min(signal[nz])) / (max(signal[nz]) - min(signal[nz]))
ybase[zeros] <- 5.0 # zero-weight rows sit inside the range
ypert <- ybase
ypert[zeros] <- 5.5 # a different interior value

fitZeroWeight <- function(resp, nodePrior) {
  ctrl <- dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 25L,
    updateState = FALSE,
    rngSeed = 5L
  )
  args <- list(x, resp, weights = w, control = ctrl, sigma = 1.0)
  if (!is.null(nodePrior)) {
    args$node.prior <- nodePrior
  }
  suppressWarnings(do.call(dbarts, args))$run(50L, 50L)$train
}

# constant leaves: the fit is bitwise inert to the zero-weight rows' response
# on every retained row, equal to rounding on the perturbed rows themselves,
# recovers the signal on the retained rows, and is finite everywhere
const.base <- fitZeroWeight(ybase, NULL)
const.pert <- fitZeroWeight(ypert, NULL)
expect_identical(const.base[nz, ], const.pert[nz, ])
expect_equal(const.base[zeros, ], const.pert[zeros, ], tolerance = 1e-12)
expect_true(all(is.finite(const.base)))
expect_true(cor(rowMeans(const.base)[nz], signal[nz]) > 0.8)

# linear leaves: same likelihood-exclusion oracle over designated columns
lin.base <- fitZeroWeight(ybase, dbarts:::linear(c("a", "b")))
lin.pert <- fitZeroWeight(ypert, dbarts:::linear(c("a", "b")))
expect_identical(lin.base[nz, ], lin.pert[nz, ])
expect_equal(lin.base[zeros, ], lin.pert[zeros, ], tolerance = 1e-12)
expect_true(all(is.finite(lin.base)))

rm(
  x,
  w,
  zeros,
  nz,
  signal,
  ybase,
  ypert,
  fitZeroWeight,
  const.base,
  const.pert,
  lin.base,
  lin.pert,
  n
)

# The sigma posterior must ignore zero-weight rows in its degrees of freedom,
# as the weights validator documents. Paired design: fit A on n rows all
# w = 1; fit B on the same rows plus exact duplicates at w = 0. The duplicates
# cannot change the cut points, the response scaling, the leaf prior, or any
# leaf conditional, so the two sigma posteriors coincide - unless the df
# over-counts the zero-weight rows, which deflates fit B's sigma (the unfixed
# code gives a ratio near 0.13).
set.seed(20260708)
n <- 60L
x <- matrix(runif(n * 3L), n, 3L)
y <- rowSums(x) + rnorm(n, 0, 0.5)

x2 <- rbind(x, x)
y2 <- c(y, y)
w2 <- c(rep(1, n), rep(0, n))

set.seed(1L)
fitA <- bart(x, y, ntree = 50L, ndpost = 1500L, nskip = 500L, verbose = FALSE)
set.seed(1L)
fitB <- suppressWarnings(bart(
  x2,
  y2,
  weights = w2,
  ntree = 50L,
  ndpost = 1500L,
  nskip = 500L,
  verbose = FALSE
))

sigmaRatio <- mean(fitB$sigma) / mean(fitA$sigma)
expect_true(
  sigmaRatio > 0.7 && sigmaRatio < 1.4,
  info = paste0("zero-weight sigma ratio B/A = ", round(sigmaRatio, 3))
)

rm(n, x, y, x2, y2, w2, fitA, fitB, sigmaRatio)
