# mutating a dense-backed column of a mixed dense/sparse container: the store
# owns its dense block, so cut-point requantization, a rolled-back
# transactional update, the setState cut-point replay, and a linear leaf's
# covariate regather all read the mutated values rather than the ones the
# container was built from

if (!requireNamespace("Matrix", quietly = TRUE)) {
  exit_file("Matrix not available")
}

set.seed(2718)
n <- 300L
u <- rnorm(n)
w <- rnorm(n)
nonzero <- runif(n) < 0.15
sv <- Matrix::sparseVector(
  x = 0.5 + runif(sum(nonzero)),
  i = which(nonzero),
  length = n
)
y <- 1.5 * u - w + as.double(sv) + rnorm(n, 0, 0.3)

x.frame <- data.frame(u = u, w = w)
x.frame$sv <- sv
expect_inherits(dbartsData(x.frame, y)@x, "dbartsMixedMatrix")

control <- dbartsControl(
  n.samples = 10L,
  n.burn = 0L,
  n.trees = 25L,
  n.chains = 1L,
  n.threads = 1L,
  updateState = FALSE
)

# the uniform grid the engine builds for a column: n.cuts interior points
# evenly spaced over the observed range
uniformGrid <- function(values, n.cuts) {
  min(values) + seq_len(n.cuts) * ((max(values) - min(values)) / (n.cuts + 1))
}

# CUT POINTS: handing the store back the grid it already holds must change
# nothing, exactly as it does for the dense equivalent
u.new <- u * 1.3 + 0.2

set.seed(11)
sampler.cuts <- dbarts(x.frame, y, control = control)
invisible(sampler.cuts$run(100L, 10L))
sampler.cuts$setPredictor(
  u.new,
  column = 1L,
  forceUpdate = TRUE,
  updateCutPoints = TRUE
)
ssr.mutated <- sampler.cuts$getSumsOfSquaredResiduals()
own.grid <- uniformGrid(u.new, sampler.cuts$data@n.cuts[1L])
sampler.cuts$setCutPoints(list(own.grid), 1L)
expect_equal(sampler.cuts$getSumsOfSquaredResiduals(), ssr.mutated)

x.dense <- as.matrix(dbartsData(x.frame, y)@x)
set.seed(11)
sampler.dense <- dbarts(x.dense, y, control = control)
invisible(sampler.dense$run(100L, 10L))
sampler.dense$setPredictor(
  u.new,
  column = 1L,
  forceUpdate = TRUE,
  updateCutPoints = TRUE
)
ssr.dense <- sampler.dense$getSumsOfSquaredResiduals()
sampler.dense$setCutPoints(list(own.grid), 1L)
expect_equal(sampler.dense$getSumsOfSquaredResiduals(), ssr.dense)

# ROLLBACK: a rejected transactional update leaves neither the refused values
# nor the creation-time ones behind for the next requantization to install
set.seed(55)
sampler.roll <- dbarts(x.frame, y, control = control)
invisible(sampler.roll$run(100L, 10L))
u.forced <- u * 1.15 - 0.1
sampler.roll$setPredictor(u.forced, column = 1L, forceUpdate = TRUE)
ssr.forced <- sampler.roll$getSumsOfSquaredResiduals()

# a constant column empties one side of every split on it
expect_false(sampler.roll$setPredictor(
  rep(0.5, n),
  column = 1L,
  forceUpdate = FALSE
))
expect_equal(sampler.roll$getSumsOfSquaredResiduals(), ssr.forced)
# the grid is still the creation-time one (the forced update did not refresh it)
sampler.roll$setCutPoints(
  list(uniformGrid(u, sampler.roll$data@n.cuts[1L])),
  1L
)
expect_equal(sampler.roll$getSumsOfSquaredResiduals(), ssr.forced)

# SET STATE: replaying a stored state reinstalls its cut points, which
# requantize the mutated column back to the codes the state was stored under
set.seed(77)
sampler.state <- dbarts(x.frame, y, control = control)
invisible(sampler.state$run(100L, 10L))
sampler.state$setPredictor(u * 0.8 + 0.4, column = 1L, forceUpdate = TRUE)
ssr.stored <- sampler.state$getSumsOfSquaredResiduals()
sampler.state$storeState()
stored <- sampler.state$state

sampler.state$setCutPoints(list(c(-1, 0, 1)), 1L)
sampler.state$setState(stored)
expect_equal(sampler.state$getSumsOfSquaredResiduals(), ssr.stored)

# LINEAR LEAVES: the leaf covariate is served from the store's raw, so a
# mutation that moves the signal into the covariate must be visible to the
# regather. The covariate carries no split weight, so the trees cannot reach
# the signal through the codes.
set.seed(303)
n.leaf <- 250L
covariate.old <- rnorm(n.leaf)
covariate.new <- rnorm(n.leaf)
nonzero.leaf <- runif(n.leaf) < 0.2
sv.leaf <- Matrix::sparseVector(
  x = 1 + runif(sum(nonzero.leaf)),
  i = which(nonzero.leaf),
  length = n.leaf
)
y.leaf <- 3 * covariate.new + rnorm(n.leaf, 0, 0.3)

leaf.frame <- data.frame(cv = covariate.old, noise = rnorm(n.leaf))
leaf.frame$sp <- sv.leaf

set.seed(404)
sampler.leaf <- dbarts(
  leaf.frame,
  y.leaf,
  node.prior = linear("cv"),
  tree.prior = cgm(split.probs = c(0, 1, 1)),
  control = control
)
invisible(sampler.leaf$run(150L, 10L))
ssr.noise <- sampler.leaf$getSumsOfSquaredResiduals()
sampler.leaf$setPredictor(covariate.new, column = 1L, forceUpdate = TRUE)
invisible(sampler.leaf$run(150L, 10L))
ssr.signal <- sampler.leaf$getSumsOfSquaredResiduals()

total.leaf <- sum((y.leaf - mean(y.leaf))^2)
expect_true(ssr.noise > 0.5 * total.leaf)
expect_true(ssr.signal < 0.1 * total.leaf)

# TRANSPOSED ARGUMENT: a whole-matrix mutation argument whose dim() disagrees
# with the design - same total length, dimensions swapped - is refused rather
# than reinterpreted column-major (is.matrix() is FALSE for every Matrix
# class, so the shape went unchecked and only the total-length check ran)
sampler.transposed <- dbarts(x.frame, y, control = control)
x.wrong.shape <- Matrix::t(as(matrix(rnorm(n * 3L), n, 3L), "CsparseMatrix"))
expect_error(
  sampler.transposed$setPredictor(x.wrong.shape, forceUpdate = TRUE),
  pattern = "dimension of x"
)

# ANYNA ON A CONTAINER, both directions: a container with exactly one sparse
# ordinal column carries sparseReference = NA_integer_ (the "no reference
# level" sentinel) as a length-one list element, which base anyNA() on a
# list misreads as a missing value that is not there
sampler.strict <- dbarts(x.frame, y, control = control, missing = "error")
test.one.sparse <- data.frame(u = rnorm(5L), w = rnorm(5L))
test.one.sparse$sv <- Matrix::sparseVector(
  x = c(0.6, 0.9),
  i = c(1L, 4L),
  length = 5L
)
container.no.na <- dbarts:::makeCategoricalModelMatrix(test.one.sparse)
expect_equal(length(container.no.na$sparseReference), 1L)
expect_true(is.na(container.no.na$sparseReference))
expect_silent(sampler.strict$setTestPredictor(container.no.na))

# with two or more sparse columns, sparseReference is a length>1 element -
# never inspected for NA by a top-level list anyNA() at all - so a real NA
# in sparse@x goes undetected in the other direction
sampler.strict2 <- dbarts(x.frame, y, control = control, missing = "error")
test.two.sparse <- data.frame(u = rnorm(5L))
test.two.sparse$w <- Matrix::sparseVector(
  x = c(NA_real_, 0.4),
  i = c(1L, 3L),
  length = 5L
)
test.two.sparse$sv <- Matrix::sparseVector(x = 0.7, i = 2L, length = 5L)
container.has.na <- dbarts:::makeCategoricalModelMatrix(test.two.sparse)
expect_true(length(container.has.na$sparseReference) > 1L)
expect_error(
  sampler.strict2$setTestPredictor(container.has.na),
  pattern = "missing values"
)
