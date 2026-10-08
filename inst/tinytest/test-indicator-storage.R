# an indicator expansion's storage is invisible (dec-B370): each indicator
# column is built sparse or dense by the engine's own density threshold, the
# starting sigma is the linear-model estimate whichever the storage, no
# sparse-columns warning reaches a caller who gave no sparse column, and
# makeModelMatrixFromDataFrame returns a plain matrix unless the caller
# supplied a sparse column

if (!requireNamespace("Matrix", quietly = TRUE)) {
  exit_file("Matrix not available")
}

set.seed(71L)
n <- 400L
f.wide <- factor(sample.int(150L, n, replace = TRUE))
f.few <- factor(sample(c("a", "b", "c"), n, replace = TRUE))
z <- rnorm(n)
y <- 2 * z + as.integer(f.few) / 2 + rnorm(n, 0, 0.3)
d <- data.frame(f = f.wide, g = f.few, z = z)

withMode <- function(mode, expr) {
  old <- options(dbarts.sparseIndicators = mode)
  on.exit(options(old))
  force(expr)
}
# stand in for Matrix being absent, in the one place that asks
withoutMatrix <- function(expr) {
  original <- dbarts:::matrixAvailable
  assignInNamespace("matrixAvailable", function() FALSE, "dbarts")
  on.exit(assignInNamespace("matrixAvailable", original, "dbarts"))
  force(expr)
}
countWarnings <- function(expr) {
  seen <- character()
  value <- withCallingHandlers(
    expr,
    warning = function(w) {
      seen[[length(seen) + 1L]] <<- conditionMessage(w)
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = seen)
}

# --- the exported builder returns a plain matrix, with and without Matrix ---

mm <- makeModelMatrixFromDataFrame(d)
expect_true(is.matrix(mm))
expect_identical(dim(mm), c(n, 141L + 2L + 1L))
expect_identical(
  withoutMatrix(makeModelMatrixFromDataFrame(d)),
  mm
)
expect_identical(
  withMode("sparse", makeModelMatrixFromDataFrame(d)),
  mm
)
# the same through bayesTree's spelling
expect_true(is.matrix(makeind(d)))
# a column the caller made sparse is the one case a container comes back
d.sparse <- d
d.sparse$s <- Matrix::sparseVector(
  x = 1 + runif(40L),
  i = sort(sample.int(n, 40L)),
  length = n
)
expect_inherits(makeModelMatrixFromDataFrame(d.sparse), "dbartsMixedMatrix")
rm(d.sparse)

# --- the training design builds each indicator column by its density ---

mi <- dbarts:::makeIndicatorModelMatrix(d)
expect_inherits(mi, "dbartsMixedMatrix")
# the 150 level factor's columns are sparse; the three level factor's, with
# a density of a third, and z are dense
sparseColumns <- mi$map < 0L
expect_true(all(sparseColumns[grepl("^f\\.", colnames(mi))]))
expect_false(any(sparseColumns[grepl("^g\\.", colnames(mi))]))
expect_false(sparseColumns[colnames(mi) == "z"])
# the same values either way
expect_equal(
  unname(as.matrix(mi)),
  unname(mm),
  check.attributes = FALSE
)
expect_identical(colnames(mi), colnames(mm))
expect_identical(attr(mi, "drop"), attr(mm, "drop"))
# a factor all of whose indicator columns are dense leaves a plain matrix
expect_true(is.matrix(dbarts:::makeIndicatorModelMatrix(d["g"])))
# without Matrix, or forced dense, every column is dense
expect_true(is.matrix(withoutMatrix(dbarts:::makeIndicatorModelMatrix(d))))
expect_true(is.matrix(withMode("dense", dbarts:::makeIndicatorModelMatrix(d))))

# --- the starting sigma is the same estimate on either storage, and no
# sparse-columns warning is raised for storage the caller did not choose ---

fitSigest <- function(mode) {
  countWarnings(withMode(mode, bartBT(
    d,
    y,
    ntree = 5L,
    nskip = 5L,
    ndpost = 5L,
    nchain = 1L,
    nthread = 1L,
    keeptrees = TRUE,
    verbose = FALSE
  )))
}
auto <- fitSigest("auto")
dense <- fitSigest("dense")
sparse <- fitSigest("sparse")
expect_true(dbarts:::predictorSourceIsSparse(auto$value$fit$data@x))
expect_true(dbarts:::predictorSourceIsSparse(sparse$value$fit$data@x))
expect_false(dbarts:::predictorSourceIsSparse(dense$value$fit$data@x))
expect_equal(auto$value$sigest, dense$value$sigest, tolerance = 1e-10)
expect_equal(sparse$value$sigest, dense$value$sigest, tolerance = 1e-10)
expect_identical(auto$warnings, character())
expect_identical(sparse$warnings, character())
expect_identical(dense$warnings, character())
# and it is the linear model's, not the marginal sd's
expect_equal(
  dense$value$sigest,
  summary(lm(y ~ mm))$sigma,
  tolerance = 1e-6
)
no.matrix <- countWarnings(withoutMatrix(bartBT(
  d,
  y,
  ntree = 5L,
  nskip = 5L,
  ndpost = 5L,
  nchain = 1L,
  nthread = 1L,
  verbose = FALSE
)))
expect_equal(no.matrix$value$sigest, dense$value$sigest, tolerance = 1e-10)
expect_identical(no.matrix$warnings, character())

# a sparse column the caller supplied still falls back, with its warning
d.caller <- d["z"]
d.caller$s <- Matrix::sparseVector(
  x = 1 + runif(40L),
  i = sort(sample.int(n, 40L)),
  length = n
)
caller <- countWarnings(suppressMessages(dbarts::dbarts(
  y ~ .,
  d.caller,
  control = dbartsControl(
    n.trees = 5L,
    n.burn = 0L,
    n.samples = 5L,
    n.chains = 1L,
    n.threads = 1L,
    verbose = FALSE
  )
)))
expect_true(any(grepl("sparse-backed predictor columns", caller$warnings)))
