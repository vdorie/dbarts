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
  countWarnings(withMode(
    mode,
    bartBT(
      d,
      y,
      ntree = 5L,
      nskip = 5L,
      ndpost = 5L,
      nchain = 1L,
      nthread = 1L,
      keeptrees = TRUE,
      verbose = FALSE
    )
  ))
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

# --- the starting sigma by sparse QR equals the dense fit's, rank-deficient
# designs (the full indicator sets of two factors and an intercept), weights,
# an offset and missing values included ---

set.seed(5L)
n2 <- 300L
d2 <- data.frame(
  f = factor(sample.int(40L, n2, replace = TRUE)),
  g = factor(sample.int(30L, n2, replace = TRUE)),
  z = rnorm(n2)
)
y2 <- rnorm(n2) + d2$z
d2$g[c(4L, 90L)] <- NA
y2[7L] <- NA
w2 <- runif(n2, 0.5, 2)
w2[11L] <- 0
o2 <- rnorm(n2, 0, 0.2)
mi2 <- dbarts:::makeIndicatorModelMatrix(d2)
md2 <- dbarts:::makeIndicatorModelMatrix(d2, storage = "dense")
expect_true(dbarts:::predictorSourceIsSparse(mi2))
expect_true(is.matrix(md2))
for (case in list(
  list(NULL, NULL),
  list(w2, NULL),
  list(NULL, o2),
  list(w2, o2)
)) {
  expect_equal(
    dbarts:::sparseResidualStandardError(y2, mi2, case[[1L]], case[[2L]]),
    dbarts:::residualStandardError(
      y2,
      dbarts:::sigmaDesignMatrix(md2),
      case[[1L]],
      case[[2L]]
    ),
    tolerance = 1e-8
  )
}

# --- the sparse estimate matches the dense one where the pivot test and the
# QR's shape are tested: a numeric column collinear with another to a given
# relative size, a column in other units, and more columns than rows ---

sparseVsDense <- function(y, x, ...) {
  expect_equal(
    dbarts:::sparseResidualStandardError(
      y,
      dbarts:::makeIndicatorModelMatrix(x),
      ...
    ),
    dbarts:::residualStandardError(
      y,
      dbarts:::sigmaDesignMatrix(
        dbarts:::makeIndicatorModelMatrix(x, storage = "dense")
      ),
      ...
    ),
    tolerance = 1e-8
  )
}
set.seed(21)
n3 <- 120L
d3 <- data.frame(
  f = factor(sample.int(6L, n3, replace = TRUE)),
  a = rnorm(n3),
  b = rnorm(n3)
)
y3 <- rnorm(n3) + d3$a
# collinear to 1e-9 and 5e-8 (dropped by lm's 1e-7 tolerance), and to 2e-7,
# 1e-5, 3e-5 and 1e-3 (kept)
for (eps in c(0, 1e-9, 5e-8, 2e-7, 1e-5, 3e-5, 1e-3)) {
  d3c <- d3
  d3c$c <- d3$a + eps * rnorm(n3)
  sparseVsDense(y3, d3c, NULL, NULL)
}
# a column in other units
for (scale in c(1e9, 1e-9)) {
  d3s <- d3
  d3s$a <- scale * d3$a
  sparseVsDense(y3, d3s, NULL, NULL)
}
# more columns than rows: 200- and 60-level factors on 150 rows
n4 <- 150L
d4 <- data.frame(
  f = factor(sample.int(200L, n4, replace = TRUE)),
  g = factor(sample.int(60L, n4, replace = TRUE)),
  z = rnorm(n4)
)
y4 <- rnorm(n4) + d4$z
mi4 <- dbarts:::makeIndicatorModelMatrix(d4)
expect_true(dbarts:::predictorSourceIsSparse(mi4))
expect_true(ncol(mi4) > n4)
sparseVsDense(y4, d4, NULL, NULL)
# and the fit that reaches it raises no warning for a design the user did not
# make dependent
expect_silent(dbarts(
  y4 ~ .,
  d4,
  control = dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 5L,
    updateState = FALSE
  )
))
sparseVsDense(y4, d4, runif(n4, 0.5, 2), rnorm(n4, 0, 0.2))
# full rank at n columns leaves no residual degrees of freedom: not finite,
# as on the dense path
d5 <- data.frame(f = factor(seq_len(30L)), z = rnorm(30L))
expect_false(is.finite(dbarts:::sparseResidualStandardError(
  rnorm(30L),
  dbarts:::makeIndicatorModelMatrix(d5),
  NULL,
  NULL
)))

# --- the density threshold is inclusive: an indicator column at exactly
# sparseDensityThreshold of the rows is built sparse, one row more is dense ---

at <- factor(rep(c("a", "b", "c"), c(20L, 41L, 39L)))
above <- factor(rep(c("a", "b", "c"), c(21L, 40L, 39L)))
expect_true(dbarts:::factorMaySparse(at, "auto"))
expect_false(dbarts:::factorMaySparse(above, "auto"))
sparseMap <- function(m) {
  if (inherits(m, "dbartsMixedMatrix")) m$map < 0L else rep(FALSE, ncol(m))
}
mi.at <- dbarts:::makeIndicatorModelMatrix(data.frame(f = at))
expect_inherits(mi.at, "dbartsMixedMatrix")
expect_identical(sparseMap(mi.at), c(TRUE, FALSE, FALSE))
expect_true(is.matrix(dbarts:::makeIndicatorModelMatrix(data.frame(f = above))))
# a row coded missing is stored in every column, so it counts toward it
at.na <- factor(rep(c("a", "b", "c", NA), c(19L, 41L, 39L, 1L)))
expect_identical(
  sparseMap(dbarts:::makeIndicatorModelMatrix(data.frame(f = at.na))),
  c(TRUE, FALSE, FALSE)
)
above.na <- factor(rep(c("a", "b", "c", NA), c(20L, 40L, 39L, 1L)))
expect_true(is.matrix(
  dbarts:::makeIndicatorModelMatrix(data.frame(f = above.na))
))
