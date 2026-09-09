# Seeded-drift tripwire, not a correctness test: pins draws from a fixed
# seed so unintended RNG-affecting changes are caught. After an intentional
# shift, run tools/regenerate-snapshots.R and eyeball that the new values
# move by a plausible magnitude.

# The pinned values only mean anything on the reference build: its draw path
# is scalar and fixed-order, where the shipped build's is free to reassociate
# and lands elsewhere from the same seed. This exit guards test runs only;
# tools/regenerate-snapshots.R evaluates this file outside tinytest, where
# exit_file() returns its message and stops nothing, so that tool carries its
# own refusal.
if (!identical(dbarts:::buildInfo()$mode, "reference")) {
  exit_file("seeded-drift snapshots are pinned to the reference build")
}

source(
  system.file("common", "friedmanData.R", package = "dbarts"),
  local = TRUE
)

# xbart's fold/subsample splits are drawn via R's sample() when no seed is
# given; pin the sampling kind so results don't depend on what other test
# files left behind
oldSampleKind <- RNGkind()[3L]
suppressWarnings(RNGkind(sample.kind = "Rejection"))

x <- testData$x
y <- testData$y
k <- c(4, 8)

set.seed(0L)
xval.kf <- dbarts::xbart(
  x,
  y,
  method = "k-fold",
  n.reps = 4L,
  n.samples = 20L,
  n.burn = c(10L, 5L),
  n.test = 5,
  k = k,
  n.threads = 1L
)

set.seed(0L)
xval.rs <- dbarts::xbart(
  x,
  y,
  method = "random subsample",
  n.reps = 20L,
  n.samples = 20L,
  n.burn = c(10L, 5L),
  k = k,
  n.threads = 1L
)

res.kf <- apply(xval.kf, 2L, mean)
res.rs <- apply(xval.rs, 2L, mean)

reference <- list(
  kf = c(2.4366063501151, 4.5532985484437),
  rs = c(2.44200830568094, 4.64040529040152)
)

expect_equal(unname(res.kf), reference$kf)
expect_equal(unname(res.rs), reference$rs)

# the methods estimate the same quantity; at these tiny run lengths the
# Monte Carlo error leaves them agreeing only loosely
expect_true(all(abs(res.rs - res.kf) / res.kf < 0.15))

suppressWarnings(RNGkind(sample.kind = oldSampleKind))
rm(oldSampleKind, res.rs, res.kf, xval.rs, xval.kf, k, y, x, testData)
