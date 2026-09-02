# A test container's dense factor columns reach the replay as the integer
# level codes they already are, and the replay must route rows off them
# exactly as it does off the same values widened into doubles. The container
# keeps a sparse block, which is what keeps a test set resident rather than
# densified on the R side, so a CSC-backed column sits beside the coded ones.

if (!requireNamespace("Matrix", quietly = TRUE)) {
  exit_file("Matrix not available")
}

set.seed(2024L)
n.rr <- 200L
levels.rr <- c("a", "b", "c", "d") # "d" is declared but never observed
f.rr <- factor(sample(levels.rr[1:3], n.rr, replace = TRUE), levels = levels.rr)
f.rr[c(7L, 88L)] <- NA
o.rr <- ordered(sample(seq_len(4L), n.rr, replace = TRUE), levels = seq_len(4L))
o.rr[c(12L, 150L)] <- NA
z.rr <- rnorm(n.rr)
levels.sparse.rr <- c("p", "q", "r")
s.rr <- sparseFactor(
  factor(
    sample(levels.sparse.rr, n.rr, replace = TRUE, prob = c(8, 1, 1)),
    levels = levels.sparse.rr
  ),
  reference = "p"
)
y.rr <- ifelse(is.na(f.rr), 0, as.integer(f.rr)) +
  0.5 * ifelse(is.na(o.rr), 0, as.integer(o.rr)) +
  rnorm(n.rr, 0, 0.5)
df.rr <- data.frame(f = f.rr, o = o.rr, z = z.rr)
df.rr$s <- s.rr

sampler.rr <- dbarts(
  df.rr,
  y.rr,
  control = dbartsControl(
    n.trees = 15L,
    n.chains = 1L,
    n.threads = 1L,
    n.samples = 5L,
    n.burn = 0L,
    keepTrees = TRUE,
    updateState = FALSE
  )
)
expect_inherits(sampler.rr$data@x, "dbartsMixedMatrix")
expect_equal(sampler.rr$data@varTypes, c(1L, 2L, 0L, 1L))
invisible(sampler.rr$run(30L, 5L))

# the test rows carry the declared-but-unobserved level and missing values in
# both factor kinds
n.test.rr <- 40L
f.test.rr <- factor(rep_len(levels.rr, n.test.rr), levels = levels.rr)
f.test.rr[3L] <- NA
o.test.rr <- ordered(rep_len(seq_len(4L), n.test.rr), levels = seq_len(4L))
o.test.rr[5L] <- NA
df.test.rr <- data.frame(
  f = f.test.rr,
  o = o.test.rr,
  z = rnorm(n.test.rr)
)
df.test.rr$s <- sparseFactor(
  factor(rep_len(levels.sparse.rr, n.test.rr), levels = levels.sparse.rr),
  reference = "p"
)
expect_true(any(f.test.rr == "d", na.rm = TRUE))

# the container the read-only funnel parses: its dense columns are the factors
# themselves, so they cross the bridge as codes
container.codes.rr <- dbarts:::makeCategoricalModelMatrix(df.test.rr)
# the same container with those columns spelled as the doubles the bridge used
# to widen them into, which takes the other channel
container.doubles.rr <- container.codes.rr
for (k in seq_along(container.doubles.rr$dense)) {
  column <- container.doubles.rr$dense[[k]]
  if (is.factor(column)) {
    container.doubles.rr$dense[[k]] <- as.double(as.integer(column)) - 1
  }
}
expect_true(any(vapply(container.codes.rr$dense, is.factor, FALSE)))
expect_false(any(vapply(container.doubles.rr$dense, is.factor, FALSE)))

ptr.rr <- sampler.rr$getPointer()
predictOf <- function(container) {
  .Call(dbarts:::C_dbarts_bartcore_predict, ptr.rr, container, NULL, 1L)
}

# THE PIN: the channel changes how the rows cross the bridge, nothing about
# where they route
expect_identical(predictOf(container.codes.rr), predictOf(container.doubles.rr))

# and the user-facing entrance reaches the same answer, which is what makes
# the coded arm the one predict() actually takes
expect_identical(
  sampler.rr$predict(df.test.rr),
  predictOf(container.doubles.rr)
)

# the saved-tree replay is the second read-only entrance off the same funnel:
# its n column counts the rows each node routes, so a channel that mis-read a
# code would move it
treesOf <- function(newdata) {
  sampler.rr$getTrees(
    treeNums = 1:15,
    chainNums = 1L,
    sampleNums = 1L,
    newdata = newdata
  )
}
expect_identical(treesOf(df.test.rr)$n, treesOf(container.doubles.rr)$n)

rm(
  n.rr,
  levels.rr,
  f.rr,
  o.rr,
  z.rr,
  levels.sparse.rr,
  s.rr,
  y.rr,
  df.rr,
  sampler.rr,
  n.test.rr,
  f.test.rr,
  o.test.rr,
  df.test.rr,
  container.codes.rr,
  container.doubles.rr,
  ptr.rr,
  predictOf,
  treesOf
)
