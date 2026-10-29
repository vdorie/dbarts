# a fit's data@x is a dbartsMixedMatrix, and as.data.frame() on it gives the
# predictors back in the form predict() and setTestPredictor take

set.seed(8)
nMD <- 60L
dfMD <- data.frame(
  y = rnorm(nMD),
  x1 = rnorm(nMD),
  f = factor(
    sample(c("a", "b", "c"), nMD, TRUE),
    levels = c("a", "b", "c", "d")
  )
)
fitMD <- function(formula) {
  dbarts::bart(
    formula,
    dfMD,
    n.trees = 10L,
    n.samples = 10L,
    n.burn = 10L,
    n.chains = 1L,
    n.threads = 1L,
    verbose = FALSE,
    keepTrees = TRUE,
    seed = 1L
  )
}

# a numeric and a factor column: the factor keeps its labels and its
# declared level that no row takes, and predicts the same draws
fitMixed <- fitMD(y ~ x1 + f)
xMixed <- fitMixed$fit$data@x
expect_inherits(xMixed, "dbartsMixedMatrix")
asDf <- as.data.frame(xMixed)
expect_inherits(asDf, "data.frame")
expect_identical(names(asDf), c("x1", "f"))
expect_identical(attr(asDf, "row.names"), seq_len(nMD))
expect_true(is.numeric(asDf$x1))
expect_identical(asDf$x1, dfMD$x1)
expect_true(is.factor(asDf$f))
expect_identical(levels(asDf$f), c("a", "b", "c", "d"))
expect_identical(as.character(asDf$f), as.character(dfMD$f))
expect_false("d" %in% as.character(asDf$f))
expect_identical(predict(fitMixed, asDf), predict(fitMixed, dfMD))
expect_identical(
  rownames(as.data.frame(xMixed, row.names = paste0("r", seq_len(nMD))))[1:2],
  c("r1", "r2")
)

# all numeric
fitNumeric <- fitMD(y ~ x1)
asDfNumeric <- as.data.frame(fitNumeric$fit$data@x)
expect_identical(names(asDfNumeric), "x1")
expect_identical(predict(fitNumeric, asDfNumeric), predict(fitNumeric, dfMD))

# base R's data frame constructors reach the method
expect_identical(data.frame(xMixed), asDf)
expect_identical(as.matrix(xMixed)[, "f"], as.double(as.integer(dfMD$f) - 1L))

# a sparse factor column rides through as a factor, and sparse numeric
# columns as numbers
set.seed(9)
sf <- factor(sample(c("a", "b", "c"), nMD, TRUE, prob = c(0.7, 0.2, 0.1)))
dfSparse <- data.frame(y = rnorm(nMD) + (sf == "b"), z = runif(nMD) + 1)
dfSparse$sf <- dbarts::sparseFactor(sf, reference = "b")
fitSparse <- dbarts::bart(
  y ~ sf + z,
  dfSparse,
  n.trees = 10L,
  n.samples = 10L,
  n.burn = 10L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE,
  keepTrees = TRUE,
  seed = 1L
)
asDfSparse <- as.data.frame(fitSparse$fit$data@x)
expect_identical(names(asDfSparse), c("z", "sf"))
expect_true(is.factor(asDfSparse$sf))
expect_identical(levels(asDfSparse$sf), c("a", "b", "c"))
expect_identical(as.character(asDfSparse$sf), as.character(sf))
expect_identical(asDfSparse$z, dfSparse$z)
expect_identical(
  predict(fitSparse, asDfSparse),
  predict(fitSparse, dfSparse)
)
