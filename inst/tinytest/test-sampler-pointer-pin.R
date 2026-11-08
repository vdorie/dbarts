# The engine pointer holds exactly the R vectors its store points into: a CSC
# build's @i and @x, which R must neither free nor change while a column
# quantizes from them, until every CSC column owns its values, and nothing
# for a dense build, so a saved dense sampler writes its data once. A
# finalizer on an environment a data object carries says when the object is
# freed.

if (!requireNamespace("Matrix", quietly = TRUE)) {
  exit_file("Matrix not available")
}

control <- dbartsControl(
  n.samples = 5L,
  n.burn = 5L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)

attachProbe <- function(data, counter) {
  probe <- new.env()
  reg.finalizer(probe, function(e) counter$freed <- counter$freed + 1L)
  attr(data, "probe") <- probe
  data
}

set.seed(7L)
n <- 120L
p <- 4L
x <- matrix(rnorm(n * p), n, p)
y <- x[, 1L] + rnorm(n, 0, 0.5)

# a saved dense sampler writes its data once: the pointer serializes without it
sampler <- dbarts(x, y, n.trees = 10L, control = control)
invisible(sampler$run())
expect_true(length(serialize(sampler$data, NULL)) > 1024L)
expect_true(length(serialize(sampler$pointer, NULL)) < 1024L)
sampler$setData(dbartsData(x + 0.01, y))
expect_true(length(serialize(sampler$pointer, NULL)) < 1024L)

# the data a dense sampler was created from is let go at its first $setData,
# while the sampler lives
counter <- new.env()
counter$freed <- 0L
sampler.created <- local(
  dbarts(
    attachProbe(dbartsData(x, y), counter),
    n.trees = 10L,
    control = control
  )
)
invisible(gc())
invisible(gc())
expect_equal(counter$freed, 0L)
sampler.created$setData(dbartsData(x + 0.01, y))
invisible(gc())
invisible(gc())
expect_equal(counter$freed, 1L)
rm(sampler.created)

# a CSC build points into its dgCMatrix's @i and @x for the store's life, so
# the pointer holds those two vectors, of a bare dgCMatrix and of a mixed
# container's sparse part alike
slicesBytes <- function(sparse) {
  length(serialize(sparse@i, NULL)) + length(serialize(sparse@x, NULL))
}
x.dense <- matrix(0, n, p)
for (j in seq_len(p)) {
  nonzero <- runif(n) < 0.2
  x.dense[nonzero, j] <- 0.5 + runif(sum(nonzero))
}
y.sparse <- x.dense[, 1L] + rnorm(n, 0, 0.3)
sampler.sparse <- dbarts(
  methods::as(x.dense, "CsparseMatrix"),
  y.sparse,
  n.trees = 10L,
  control = control
)
expect_true(
  length(serialize(sampler.sparse$pointer, NULL)) >=
    slicesBytes(sampler.sparse$data@x)
)
expect_true(
  length(serialize(sampler.sparse$pointer, NULL)) <
    slicesBytes(sampler.sparse$data@x) + 1024L
)

# a whole $setPredictor leaves every column owning its values, so the
# pointer lets the creation vectors go; later updates and a run read only
# what the store owns
x.next <- x.dense
x.next[x.next == 0] <- 7
sampler.sparse$setPredictor(x.next, forceUpdate = TRUE)
expect_true(length(serialize(sampler.sparse$pointer, NULL)) < 1024L)
invisible(gc(full = TRUE))
sampler.sparse$setPredictor(x.dense[, 2L], column = 2L, forceUpdate = TRUE)
expect_equal(as.matrix(sampler.sparse$data@x)[, 2L], x.dense[, 2L])
expect_true(all(is.finite(sampler.sparse$run(0L, 5L)$sigma)))
rm(sampler.sparse)

set.seed(42L)
n.mixed <- 400L
x1 <- rnorm(n.mixed)
sv <- Matrix::sparseVector(
  x = 0.5 + runif(60L),
  i = sort(sample.int(n.mixed, 60L)),
  length = n.mixed
)
sm.dense <- matrix(0, n.mixed, 2L)
for (j in 1:2) {
  nonzero <- runif(n.mixed) < 0.1
  sm.dense[nonzero, j] <- 0.5 + runif(sum(nonzero))
}
x.frame <- data.frame(
  x1 = x1,
  f = factor(sample(c("a", "b", "c"), n.mixed, replace = TRUE))
)
x.frame$sv <- sv
x.frame$sm <- methods::as(sm.dense, "CsparseMatrix")
y.mixed <- 2 * as.double(sv) - 1.5 * x1 + rnorm(n.mixed, 0, 0.3)
sampler.mixed <- dbarts(
  x.frame,
  y.mixed,
  n.trees = 10L,
  control = control,
  sigest = 1
)
creationBytes <- slicesBytes(sampler.mixed$data@x$sparse)
expect_true(length(serialize(sampler.mixed$pointer, NULL)) >= creationBytes)
expect_true(
  length(serialize(sampler.mixed$pointer, NULL)) < creationBytes + 1024L
)

# a CSC column's update splices a new sparse part into data@x, leaving the
# creation vectors to the pointer alone; a collection then must not free
# them, which the next column's first update reads (an R-loaded sanitizer
# build reports the read of a freed vector)
sv.new <- as.double(sv)
sv.new[sv.new != 0] <- sv.new[sv.new != 0] * 1.5
sv.new[c(1L, 2L)] <- c(0.9, 1.1)
sampler.mixed$setPredictor(sv.new, column = 3L, forceUpdate = TRUE)
invisible(gc(full = TRUE))
sm1.new <- sm.dense[, 1L] * 2
sampler.mixed$setPredictor(sm1.new, column = 4L, forceUpdate = TRUE)
invisible(gc(full = TRUE))
expect_true(length(serialize(sampler.mixed$pointer, NULL)) >= creationBytes)
# column 5 is the last CSC column still reading its creation slice; once it
# owns its values the pointer holds nothing
sampler.mixed$setPredictor(
  sm.dense[, 2L] + 0.5,
  column = 5L,
  forceUpdate = TRUE
)
expect_true(length(serialize(sampler.mixed$pointer, NULL)) < 1024L)
invisible(gc(full = TRUE))
expect_equal(as.matrix(sampler.mixed$data@x)[, 3L], sv.new)
expect_equal(as.matrix(sampler.mixed$data@x)[, 4L], sm1.new)
expect_true(all(is.finite(sampler.mixed$run(0L, 5L)$sigma)))
rm(sampler.mixed)

# of k data objects handed to $setData, every one replaced is freed while the
# sampler lives
counter <- new.env()
counter$freed <- 0L
k <- 4L
for (i in seq_len(k)) {
  local({
    data <- dbartsData(x + rnorm(n * p, 0, 1e-3), y)
    sampler$setData(attachProbe(data, counter))
  })
}
invisible(gc())
invisible(gc())
invisible(gc())
expect_true(counter$freed >= k - 1L)
rm(sampler)
invisible(gc())
invisible(gc())
expect_equal(counter$freed, k)
