# predict() frame/sparse parity: a data.frame/sparse test set is coded against
# the training levels (validateXTest) and then RIDES TO THE ENGINE AS IT IS -
# predict and getTrees route its rows off the container's own storage, with no
# n x p materialization anywhere. Every assertion here is bitwise against the
# dense twin, which is the whole claim.

if (!requireNamespace("Matrix", quietly = TRUE)) {
  exit_file("Matrix not available")
}

set.seed(70L)
n <- 200L
x1.tr <- rnorm(n)
levels.g <- paste0("L", 1:4)
codes.tr <- sample.int(4L, n, replace = TRUE)
codes.tr[1L] <- 1L
g.tr <- factor(levels.g[codes.tr], levels = levels.g)
y.tr <- 0.3 * codes.tr + rnorm(n)
train <- data.frame(x1 = x1.tr, g = g.tr)

sampler <- dbarts(
  train,
  y.tr,
  sigma = 1.0,
  control = dbartsControl(
    n.trees = 20L,
    n.chains = 1L,
    n.threads = 1L,
    updateState = FALSE
  )
)
sampler$run(10L, 20L)

n.test <- 30L
x1.te <- rnorm(n.test)
codes.te <- sample.int(4L, n.test, replace = TRUE)
codes.te[1L] <- 1L # the first level appears explicitly in the sparse column
g.te <- factor(levels.g[codes.te], levels = levels.g)

# training columns are ordered x1, g; the sparse test frame lists them in the
# opposite order, so validateXTest's column-reorder path runs before the
# resulting container reaches the engine
test.dense <- data.frame(x1 = x1.te, g = g.te)
test.sparse <- data.frame(g = g.te)
test.sparse$g <- sparseFactor(g.te, reference = "L2")
test.sparse$x1 <- x1.te

pred.dense <- sampler$predict(test.dense)
pred.sparse <- sampler$predict(test.sparse)
expect_identical(pred.dense, pred.sparse)

# offsets carry through the resident path identically as well
offset.te <- rnorm(n.test)
pred.dense.off <- sampler$predict(test.dense, offset.te)
pred.sparse.off <- sampler$predict(test.sparse, offset.te)
expect_identical(pred.dense.off, pred.sparse.off)

# a plain dense matrix test set takes the frozen dense entry, untouched
x.mat <- extract(sampler, "predictors")
pred.mat <- sampler$predict(x.mat[1:10L, , drop = FALSE])
expect_true(is.numeric(pred.mat))
expect_false(anyNA(pred.mat))

# predict.bart routes newdata through the same path; a mixed frame (dense +
# sparseFactor) predicts identically whichever way the factor column arrives
fit <- bart2(
  train,
  y.tr,
  n.samples = 20L,
  n.burn = 10L,
  n.trees = 20L,
  n.chains = 1L,
  n.threads = 1L,
  keepTrees = TRUE,
  verbose = FALSE
)
ev.dense <- predict(fit, test.dense)
ev.sparse <- predict(fit, test.sparse)
expect_identical(ev.dense, ev.sparse)

# getTrees' newdata takes the same coding boundary predict does, and routes
# off the same storage, so tree replay counts sparse and dense identically
trees.dense <- fit$fit$getTrees(1L, 1L, 1L, newdata = test.dense)
trees.sparse <- fit$fit$getTrees(1L, 1L, 1L, newdata = test.sparse)
expect_identical(trees.dense, trees.sparse)

# a test level absent from training errors identically for a sparse test
# column as it does for a dense one
unseenLabels <- c("L5", as.character(g.te[-1L]))
expect_error(
  sampler$predict(data.frame(x1 = x1.te, g = factor(unseenLabels))),
  pattern = "has levels not present in the"
)
test.unseen <- data.frame(x1 = x1.te)
test.unseen$g <- sparseFactor(
  unseenLabels,
  levels = c(levels.g, "L5"),
  reference = "L2"
)
expect_error(
  sampler$predict(test.unseen),
  pattern = "has levels not present in the"
)

# ---- a bare dgCMatrix test set against a numeric design ----
# validateXTest wraps it as an all-sparse container, so a numeric-only design
# (no factor levels to code against) predicts resident as well
set.seed(71L)
n.num <- 150L
x.num <- matrix(
  rnorm(n.num * 3L),
  n.num,
  3L,
  dimnames = list(NULL, c("a", "b", "c"))
)
y.num <- x.num[, 1L] - 0.5 * x.num[, 3L] + rnorm(n.num)
sampler.num <- dbarts(
  x.num,
  y.num,
  sigma = 1.0,
  control = dbartsControl(
    n.trees = 15L,
    n.chains = 1L,
    n.threads = 1L,
    updateState = FALSE
  )
)
sampler.num$run(10L, 15L)

n.sparse.test <- 40L
dense.entries <- matrix(0, n.sparse.test, 3L)
nonzero <- sample.int(n.sparse.test * 3L, 25L)
dense.entries[nonzero] <- rnorm(25L)
colnames(dense.entries) <- c("a", "b", "c")
sparse.test <- as(as(dense.entries, "CsparseMatrix"), "generalMatrix")
expect_identical(
  sampler.num$predict(sparse.test),
  sampler.num$predict(dense.entries)
)

# the bridge takes a bare dgCMatrix as the all-CSC container it is; the R
# surface wraps one before it gets there, so reach the entry point directly
expect_identical(
  .Call(
    dbarts:::C_dbarts_bartcore_predict,
    sampler.num$getPointer(),
    sparse.test,
    NULL
  ),
  sampler.num$predict(dense.entries)
)

# an offset rides the resident path identically, and a zero-row test set
# reaches the same refusal a dense one does
offset.sparse <- rnorm(n.sparse.test)
expect_identical(
  sampler.num$predict(sparse.test, offset.sparse),
  sampler.num$predict(dense.entries, offset.sparse)
)
expect_error(
  sampler.num$predict(sparse.test[integer(0L), , drop = FALSE]),
  pattern = "requires rows"
)
expect_error(
  sampler.num$predict(dense.entries[integer(0L), , drop = FALSE]),
  pattern = "requires rows"
)

# ---- a leaf covariate cannot arrive sparse ----
# the linear leaf reads contiguous raw values for its designated column, which
# CSC storage does not serve; predict refuses it exactly as setTestPredictor
sampler.leaf <- dbarts(
  x.num,
  y.num,
  sigma = 1.0,
  node.prior = linear("a"),
  control = dbartsControl(
    n.trees = 15L,
    n.chains = 1L,
    n.threads = 1L,
    updateState = FALSE
  )
)
sampler.leaf$run(10L, 5L)
expect_error(
  sampler.leaf$predict(sparse.test),
  pattern = "leaf covariate column cannot be a sparse test column"
)

# ---- a pooled categorical column (K > 63) ----
# past 63 levels a rule's mask leaves the rule word for the per-tree pool, and
# the flat replay reads it through the side channel; a sparse test column of
# such a factor routes exactly as its dense twin
set.seed(72L)
n.pool <- 400L
levels.pool <- paste0("P", 1:70)
codes.pool <- c(1:70, sample.int(70L, n.pool - 70L, replace = TRUE))
g.pool <- factor(levels.pool[codes.pool], levels = levels.pool)
train.pool <- data.frame(x1 = rnorm(n.pool), g = g.pool)
y.pool <- 0.05 * codes.pool + rnorm(n.pool)
sampler.pool <- dbarts(
  train.pool,
  y.pool,
  sigma = 1.0,
  control = dbartsControl(
    n.trees = 20L,
    n.chains = 1L,
    n.threads = 1L,
    updateState = FALSE
  )
)
sampler.pool$run(20L, 20L)

n.pool.test <- 50L
codes.pool.te <- sample.int(70L, n.pool.test, replace = TRUE)
g.pool.te <- factor(levels.pool[codes.pool.te], levels = levels.pool)
test.pool.dense <- data.frame(x1 = rnorm(n.pool.test), g = g.pool.te)
test.pool.sparse <- data.frame(x1 = test.pool.dense$x1)
test.pool.sparse$g <- sparseFactor(g.pool.te, reference = "P7")
expect_identical(
  sampler.pool$predict(test.pool.sparse),
  sampler.pool$predict(test.pool.dense)
)
expect_identical(
  sampler.pool$getTrees(1L, 1L, newdata = test.pool.sparse),
  sampler.pool$getTrees(1L, 1L, newdata = test.pool.dense)
)

# ---- the heteroscedastic list(mean, variance) route ----
# s^2(x) replays the variance forest over the same rows; both channels of the
# returned pair must agree with the dense twin, attributes included
set.seed(73L)
fit.var <- bart2(
  train,
  y.tr,
  n.samples = 15L,
  n.burn = 10L,
  n.trees = 20L,
  n.chains = 1L,
  n.threads = 1L,
  keepTrees = TRUE,
  verbose = FALSE,
  variance = TRUE,
  n.trees.variance = 10L
)
pred.var.dense <- predict(fit.var, test.dense)
pred.var.sparse <- predict(fit.var, test.sparse)
expect_identical(pred.var.dense, pred.var.sparse)
expect_false(is.null(attr(pred.var.sparse, "s")))

# and through the sampler's own entry, where the pair is still a list
raw.var.dense <- fit.var$fit$predict(test.dense)
raw.var.sparse <- fit.var$fit$predict(test.sparse)
expect_identical(names(raw.var.sparse), c("mean", "variance"))
expect_identical(raw.var.dense, raw.var.sparse)

rm(
  sampler,
  sampler.num,
  sampler.leaf,
  sampler.pool,
  fit,
  fit.var,
  train,
  train.pool,
  test.pool.dense,
  test.pool.sparse,
  x.num,
  y.num,
  y.pool,
  n.num,
  n.pool,
  n.pool.test,
  n.sparse.test,
  codes.pool,
  codes.pool.te,
  levels.pool,
  g.pool,
  g.pool.te,
  dense.entries,
  nonzero,
  sparse.test,
  offset.sparse,
  pred.var.dense,
  pred.var.sparse,
  raw.var.dense,
  raw.var.sparse,
  test.dense,
  test.sparse,
  test.unseen,
  pred.dense,
  pred.sparse,
  pred.dense.off,
  pred.sparse.off,
  ev.dense,
  ev.sparse,
  x.mat,
  pred.mat,
  trees.dense,
  trees.sparse,
  x1.tr,
  x1.te,
  levels.g,
  codes.tr,
  codes.te,
  g.tr,
  g.te,
  y.tr,
  n,
  n.test,
  offset.te,
  unseenLabels
)
