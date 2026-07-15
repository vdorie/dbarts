# predict() frame/sparse parity: the frozen dense predict entry point takes
# a data.frame/sparse test set through the same training-level coding as
# creation (validateXTest), then materializes it to a numeric matrix before
# the call - no resident sparse predict store.

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
# resulting container is materialized
test.dense <- data.frame(x1 = x1.te, g = g.te)
test.sparse <- data.frame(g = g.te)
test.sparse$g <- sparseFactor(g.te, reference = "L2")
test.sparse$x1 <- x1.te

pred.dense <- sampler$predict(test.dense)
pred.sparse <- sampler$predict(test.sparse)
expect_identical(pred.dense, pred.sparse)

# offsets carry through the materialized path identically as well
offset.te <- rnorm(n.test)
pred.dense.off <- sampler$predict(test.dense, offset.te)
pred.sparse.off <- sampler$predict(test.sparse, offset.te)
expect_identical(pred.dense.off, pred.sparse.off)

# a plain dense matrix test set is untouched by the new materialization step
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

# getTrees' newdata takes the same coding-and-materialization boundary as
# predict, so tree replay counts sparse and dense test sets identically
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

rm(
  sampler,
  fit,
  train,
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
