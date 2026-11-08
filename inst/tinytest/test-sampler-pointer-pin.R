# The engine pointer pins the data object it was built from for as long as
# the pointer lives, without writing that object a second time when a sampler
# is saved, and lets go of each object $setData replaces. A finalizer on an
# environment the data object carries says when the object is freed.

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

# a saved sampler writes its data once: the pointer serializes without it
sampler <- dbarts(x, y, n.trees = 10L, control = control)
invisible(sampler$run())
expect_true(length(serialize(sampler$data, NULL)) > 1024L)
expect_true(length(serialize(sampler$pointer, NULL)) < 1024L)
sampler$setData(dbartsData(x + 0.01, y))
expect_true(length(serialize(sampler$pointer, NULL)) < 1024L)

# a sparse build borrows the dgCMatrix's slices for the store's life, so data
# held by nothing but the pointer survives until the pointer goes
x.dense <- matrix(0, n, p)
for (j in seq_len(p)) {
  nonzero <- runif(n) < 0.2
  x.dense[nonzero, j] <- 0.5 + runif(sum(nonzero))
}
y.sparse <- x.dense[, 1L] + rnorm(n, 0, 0.3)
counter <- new.env()
counter$freed <- 0L
ptr <- local({
  sampler <- dbarts(
    methods::as(x.dense, "CsparseMatrix"),
    y.sparse,
    n.trees = 10L,
    control = control
  )
  dbarts:::recreatePointer(
    sampler$control,
    sampler$model,
    attachProbe(sampler$data, counter)
  )
})
invisible(gc())
invisible(gc())
expect_equal(counter$freed, 0L)
expect_true(.Call(dbarts:::C_dbarts_bartcore_isValidPointer, ptr))
expect_true(length(serialize(ptr, NULL)) < 1024L)
rm(ptr)
invisible(gc())
invisible(gc())
expect_equal(counter$freed, 1L)

# of k data objects handed to $setData, the k - 1 replaced are freed while
# the sampler lives
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
expect_equal(counter$freed, k - 1L)
rm(sampler)
invisible(gc())
invisible(gc())
expect_equal(counter$freed, k)
