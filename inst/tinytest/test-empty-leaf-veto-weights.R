# The empty-leaf veto counts POSITIVE-WEIGHT members, not members
# (docs/design/empty-leaf-veto.md). A weight of zero is absence, not
# downweighting, so a leaf all of whose rows carry weight zero contributes to no
# likelihood term and must be vetoed exactly as an unoccupied one is. The oracle
# is structural and needs no tree walk: routing ONLY the positive-weight rows
# through the live trees (getTrees(newdata = )) counts, per node, the rows the
# likelihood can see, so the law says every leaf reports n > 0. The zero-weight
# region is a half-space on one predictor so that whole leaves can fall inside
# it, which is what makes the check bite.

set.seed(20260812L)
n <- 400L
x <- matrix(runif(n * 3L), n, 3L)
w <- as.double(x[, 1L] <= 0.5) # the x1 > 0.5 half-space carries no weight
kept <- w > 0
signal <- 2 * x[, 1L] - x[, 2L]
y <- signal + rnorm(n, 0, 0.5)

control <- dbarts::dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 50L,
  updateState = FALSE,
  rngSeed = 7L
)

fitAndReport <- function(weights, ...) {
  args <- list(x, y, control = control, ...)
  if (!is.null(weights)) {
    args$weights <- weights
  }
  sampler <- suppressWarnings(do.call(dbarts::dbarts, args))
  fits <- sampler$run(100L, 50L)$train
  trees <- sampler$getTrees(
    current = TRUE,
    newdata = x[kept, , drop = FALSE]
  )
  list(leaves = trees[trees$var == -1L, ], fits = fits)
}

# gaussian: no live leaf may hold only zero-weight rows
gaussian <- fitAndReport(w)
expect_true(nrow(gaussian$leaves) > 50L)
expect_true(all(gaussian$leaves$n > 0L))

# non-vacuity: the same fixture under the count law - a sampler with no weights
# installed, whose veto is the member count the weighted path used to run -
# settles on leaves that hold no positive-weight row at all
counted <- fitAndReport(NULL)
expect_true(any(counted$leaves$n == 0L))

# the other half of the zero-weight contract stands: a zero-weight row is still
# partitioned and still reported a fit
expect_true(all(is.finite(gaussian$fits)))
expect_true(cor(rowMeans(gaussian$fits)[!kept], signal[!kept]) > 0.5)

# Student-t: the composed weight w * lambda is zero wherever w is, so the same
# law governs the family that carries the other shipped weight channel
student <- fitAndReport(w, resid.dist = dbarts:::student(df = 4))
expect_true(all(student$leaves$n > 0L))
expect_true(all(is.finite(student$fits)))

rm(
  n,
  x,
  w,
  kept,
  signal,
  y,
  control,
  fitAndReport,
  gaussian,
  counted,
  student
)
