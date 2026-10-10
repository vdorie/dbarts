# The starting sigma is the residual sd of the linear regression of the
# response on the predictors, factors as indicators, by the routine 'sigest'
# names: lm.fit ("dense"), the sparse routine ("sparse"), or "auto", which
# follows what the caller passed. Each is held to lm.fit on the dense
# indicator design; where that leaves no residual degrees of freedom the
# estimate is the sd of the response, said under verbose and in the summary.

if (!requireNamespace("Matrix", quietly = TRUE)) {
  exit_file("Matrix not available")
}

control <- dbartsControl(
  n.trees = 3L,
  n.chains = 1L,
  n.threads = 1L,
  updateState = FALSE
)
# the dense indicator design of a frame, missing values at their column means
indicatorForm <- function(x) {
  for (name in names(x)) {
    if (isS4(x[[name]])) {
      x[[name]] <- if (methods::is(x[[name]], "sparseVector")) {
        as.numeric(x[[name]])
      } else {
        as.matrix(x[[name]])
      }
    }
  }
  dbarts:::sigmaDesignMatrix(makeModelMatrixFromDataFrame(x, drop = FALSE))
}
onRows <- function(v, rows) if (is.null(v) || is.null(rows)) v else v[rows]
lmSigma <- function(m, y, w = NULL, o = NULL, rows = NULL) {
  if (!is.null(rows)) {
    m <- m[rows, , drop = FALSE]
  }
  args <- list(onRows(y, rows), m, onRows(w, rows), onRows(o, rows))
  do.call(dbarts:::residualStandardError, args)
}
sparseSigma <- function(x, ...) {
  if (is.data.frame(x)) {
    x <- dbarts:::makeCategoricalModelMatrix(x)
  }
  dbarts:::sparseResidualStandardError(x = x, ...)
}
# the sparse routine on a frame against lm.fit on its indicator form; where
# lm.fit has no residual degrees of freedom the routine has no estimate
agrees <- function(x, y, w = NULL, o = NULL, rows = NULL, tolerance = 1e-10) {
  expected <- lmSigma(indicatorForm(x), y, w, o, rows)
  observed <- sparseSigma(x, y = y, weights = w, offset = o, rows = rows)
  if (is.finite(expected)) {
    expect_equal(observed, expected, tolerance = tolerance)
  } else {
    expect_true(is.na(observed))
  }
}
# counts the calls of the package functions named while expr runs, into
# calls$counts, which an error in expr leaves readable
calls <- new.env()
routines <- c(lm = "residualStandardError", D = "sparseResidualStandardError")
traced <- function(expr, what = routines) {
  calls$counts <- setNames(integer(length(what)), names(what))
  originals <- lapply(what, getFromNamespace, "dbarts")
  on.exit(
    for (name in names(what)) {
      assignInNamespace(what[[name]], originals[[name]], "dbarts")
    }
  )
  for (name in names(what)) {
    local({
      name <- name
      counted <- function(...) {
        calls$counts[[name]] <- calls$counts[[name]] + 1L
        originals[[name]](...)
      }
      assignInNamespace(what[[name]], counted, "dbarts")
    })
  }
  force(expr)
  calls$counts
}
withoutMatrix <- function(expr) {
  original <- dbarts:::matrixAvailable
  assignInNamespace("matrixAvailable", function() FALSE, "dbarts")
  on.exit(assignInNamespace("matrixAvailable", original, "dbarts"))
  force(expr)
}
# the messages and warnings expr raises, its printed output discarded
signals <- function(expr) {
  seen <- list(message = character(), warning = character())
  invisible(capture.output(withCallingHandlers(
    expr,
    message = function(m) {
      seen$message <<- c(seen$message, conditionMessage(m))
      invokeRestart("muffleMessage")
    },
    warning = function(w) {
      seen$warning <<- c(seen$warning, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )))
  seen
}
silent <- list(message = character(), warning = character())
sigmaOf <- function(x, y, ...) {
  dbarts(x, y, control = control, verbose = FALSE, ...)$data@sigma
}
bartOf <- function(x, y, verbose = FALSE, ...) {
  bart(
    x,
    y,
    verbose = verbose,
    n.trees = 3L,
    n.samples = 4L,
    n.burn = 2L,
    n.chains = 1L,
    n.threads = 1L,
    ...
  )
}
bartBTOf <- function(x, y, verbose = FALSE, ...) {
  bartBT(x, y, ntree = 3L, ndpost = 2L, nskip = 1L, verbose = verbose, ...)
}
xbartOf <- function(x, y, seed = 3L, ...) {
  xbart(
    x,
    y,
    n.samples = 3L,
    n.burn = c(2L, 1L),
    n.reps = 1L,
    n.test = 3L,
    n.trees = 3L,
    n.threads = 1L,
    seed = seed,
    ...
  )
}
# rbart_vi's deprecation has been shown this session
onceKeys <- dbarts:::onceWarnState
onceKeys[["tombstone.rbart_vi"]] <- TRUE

# --- every factor enters the regression as indicators, whatever the trees
# are given: the default frame, the indicators route and a matrix of codes ---

set.seed(31)
n <- 300L
x1 <- rnorm(n, 50)
factors <- list(
  five = factor(sample(letters[1:5], n, TRUE)),
  ordered = factor(sample(letters[1:5], n, TRUE), ordered = TRUE),
  two = factor(sample(c("u", "v"), n, TRUE)),
  forty = factor(sample.int(40L, n, TRUE)),
  incomplete = factor(sample(letters[1:5], n, TRUE), levels = letters[1:7])
)
for (f in factors) {
  y <- x1 + as.integer(f) / 2 + rnorm(n, sd = 0.5)
  expected <- summary(lm(y ~ f + x1))$sigma
  frame <- data.frame(f = f, x1 = x1)
  codes <- dbartsData(frame, y)
  codes@x <- structure(
    as.matrix(codes@x),
    varTypes = attr(codes@x, "varTypes"),
    factor.levels = attr(codes@x, "factor.levels")
  )
  for (rule in c("auto", "dense", "sparse")) {
    expect_equal(sigmaOf(frame, y, sigest = rule), expected, tolerance = 1e-10)
    expect_equal(
      sigmaOf(frame, y, sigest = rule, factors = "indicators"),
      expected,
      tolerance = 1e-10
    )
    expect_equal(
      dbarts:::estimateSigmaFromLinearModel(codes, rule),
      expected,
      tolerance = 1e-10
    )
  }
}

# --- "auto" follows what the caller passed; "dense" and "sparse" force a
# routine; a fold runs the routine all rows did ---

f6 <- factor(sample(letters[1:6], n, TRUE))
y <- x1 + as.integer(f6) / 2 + rnorm(n, sd = 0.5)
frame <- data.frame(f = f6, x1 = x1)
sparse <- Matrix::rsparsematrix(n, 4L, 0.3)
sparse[c(3L, 90L), 2L] <- NA
withVector <- frame["x1"]
withVector$s <- methods::as(sparse[, 1L], "sparseVector")
withFactor <- frame["x1"]
withFactor$f <- sparseFactor(f6)
doors <- list(dbarts = sigmaOf, bart = bartOf, xbart = xbartOf)
runs <- function(door, routine, ...) {
  counts <- suppressMessages(traced(doors[[door]](...)))
  # xbart: the estimate on all rows and each of three folds'
  expected <- if (door == "xbart") 4L else 1L
  expect_equal(counts[[routine]], expected, info = door)
  expect_equal(sum(counts), expected, info = door)
}
for (door in names(doors)) {
  runs(door, "lm", as.matrix(sparse[, 3:4]), y)
  runs(door, "lm", frame, y)
  for (mode in c("sparse", "dense")) {
    old <- options(dbarts.sparseIndicators = mode)
    runs(door, "lm", frame, y, factors = "indicators")
    options(old)
  }
  runs(door, "D", as.matrix(sparse[, 3:4]), y, sigest = "sparse")
  for (x in list(sparse, withVector, withFactor)) {
    runs(door, "D", x, y)
    runs(door, "lm", x, y, sigest = "dense")
  }
}
old <- options(dbarts.sparseIndicators = "sparse")
expect_true(dbarts:::predictorSourceIsSparse(
  dbartsData(frame, y, factors = "indicators")@x
))
storedSparse <- sigmaOf(frame, y, factors = "indicators")
options(dbarts.sparseIndicators = "dense")
expect_identical(storedSparse, sigmaOf(frame, y, factors = "indicators"))
options(old)
# no warning, and "dense" gives a sparse fit its dense twin's bits, missing
# entries included
expect_identical(signals(sigmaOf(sparse, y)), silent)
expect_equal(
  sigmaOf(sparse, y),
  sigmaOf(sparse, y, sigest = "dense"),
  tolerance = 1e-10
)
expect_identical(
  sigmaOf(sparse, y, sigest = "dense"),
  sigmaOf(as.matrix(sparse), y)
)
expect_identical(
  sigmaOf(withFactor, y, sigest = "dense"),
  sigmaOf(frame[c("x1", "f")], y)
)

# --- what each door takes for 'sigest' ---

binary <- as.numeric(y > median(y))
spec <- function(data, ...) {
  dbartsSpec(data, control = control, ...)$data@sigma
}
taking <- list(
  dbarts = sigmaOf,
  bart = bartOf,
  xbart = xbartOf,
  dbartsSpec = function(x, y, ...) spec(dbartsData(x, y), ...),
  bartBT = bartBTOf,
  rbart_vi = function(x, y, ...) {
    rbart_vi(y ~ x1, x, group.by = f, n.trees = 3L, verbose = FALSE, ...)
  }
)
aFunction <- function(x, y, weights, offset) 1
for (door in names(taking)) {
  expect_error(
    taking[[door]](frame, y, sigest = "sd"),
    "unknown 'sigest' rule \"sd\"; use \"auto\", \"dense\" or \"sparse\"",
    fixed = TRUE,
    info = door
  )
  expect_error(
    withoutMatrix(taking[[door]](frame, y, sigest = "sparse")),
    "sigest = \"sparse\" requires the Matrix package",
    fixed = TRUE,
    info = door
  )
  if (door != "xbart") {
    expect_error(
      taking[[door]](frame, y, sigest = aFunction),
      paste0(
        door,
        " must be a number, \"auto\", \"dense\" or \"sparse\"; only"
      ),
      fixed = TRUE
    )
  }
}
# a rule's name where sigest has no use is treated as a number is
expect_error(bartOf(frame, binary, sigest = "sd"), "unknown 'sigest' rule")
expect_warning(
  bartOf(frame, binary, sigest = "dense"),
  "has no use for 'sigest'"
)
expect_identical(signals(sigmaOf(frame, binary, sigest = "dense")), silent)
expect_identical(sigmaOf(frame, y, sigest = "1.5"), 1.5)
expect_error(
  suppressWarnings(dbarts(frame, y, control = control, sigma = "sparse")),
  "must be coercible to numeric type"
)
for (rule in c("auto", "dense")) {
  expect_equal(
    withoutMatrix(sigmaOf(frame, y, sigest = rule)),
    summary(lm(y ~ f6 + x1))$sigma
  )
}
# a data object's own sigma: dbartsSpec keeps it unless told how to replace
# it, and the other doors estimate over it
carried <- dbartsData(frame, y)
carried@sigma <- 7
expect_identical(spec(carried), 7)
expect_identical(spec(carried, sigest = "auto"), 7)
expect_identical(spec(carried, sigest = 2), 2)
expect_identical(spec(carried, sigest = "dense"), sigmaOf(frame, y))
expect_equal(spec(carried, sigest = "sparse"), sigmaOf(frame, y))
expect_identical(
  dbarts(carried, control = control, sigest = "auto")$data@sigma,
  sigmaOf(frame, y)
)

# --- the sparse routine against lm.fit: rows dropped beside a missing factor
# value, columns sharing their stored rows, and both beside the eliminated
# factor ---

set.seed(12)
m <- 400L
f30 <- factor(sample.int(30L, m, TRUE))
g12 <- factor(sample.int(12L, m, TRUE))
z1 <- rnorm(m, 50)
ym <- z1 + as.integer(f30) / 10 + rnorm(m, sd = 0.5)
wm <- rexp(m)
om <- rnorm(m)
# 150 levels, some of one row, 8 values missing
thin <- sample(c(seq_len(150L), sample.int(150L, m - 150L, TRUE)))
thin <- factor(replace(thin, sample.int(m, 8L), NA))
single <- which(thin %in% names(which(table(thin) == 1L)))[1:3]
fn <- replace(f30, sample.int(m, 25L), NA)
gn <- replace(g12, sample.int(m, 10L), NA)
two <- factor(sample.int(2L, m, TRUE))
twoN <- replace(two, sample.int(m, 12L), NA)
fold <- sample(rep_len(1:5, m))
thinFrame <- data.frame(f = thin, z1 = z1)
agrees(thinFrame, ym, replace(rep(1, m), single, 0))
agrees(thinFrame, replace(ym, single, NA))
for (k in 1:5) {
  agrees(thinFrame, ym, rows = which(fold != k))
}
agrees(thinFrame, ym, wm, om, which(fold != 2L))
both <- data.frame(f = fn, g = gn, z1 = z1)
agrees(both, ym, wm)
agrees(both, ym, wm, rows = which(fold != 3L))
agrees(data.frame(f = fn, f2 = fn, z1 = z1), ym)
twoFrame <- data.frame(a = twoN, z1 = z1)
agrees(twoFrame, ym)
agrees(data.frame(a = twoN, f = fn, z1 = z1), ym, wm)
agrees(twoFrame, ym, rows = which(two == "1" | is.na(twoN)))
agrees(data.frame(f = fn, z1 = z1), ym, rows = which(is.na(fn)))
agrees(data.frame(f = factor(seq_len(m)), z1 = z1), ym)
# nested and duplicated factors, a column constant within a level
codes <- as.integer(f30)
coarse <- factor(ifelse(codes <= 3L, codes, 4L + codes %/% 9L))
agrees(data.frame(h = coarse, f = f30, z1 = z1), ym)
agrees(data.frame(f = fn, z1 = z1, h = coarse), ym, wm)
agrees(data.frame(a = two, f = f30, f2 = f30, z1 = z1), ym)
agrees(data.frame(f = f30, xg = rnorm(30L)[f30], z1 = z1), ym)
agrees(data.frame(f = f30, xi = as.numeric(f30 == "3"), z1 = z1), ym)
# proportional columns, and times beside the indicator of their rows, alone
# and beside a factor with missing values
stored <- runif(m) < 0.4
a <- ifelse(stored, rexp(m, 1 / 30), 0)
start <- (1.7e9 + runif(m, 0, 3600)) * stored
end <- (start + runif(m, 0, 3600)) * stored
yt <- ym + a / 30 + (end - start) / 1200
columns <- function(...) methods::as(cbind(...), "CsparseMatrix")
beside <- function(x, ...) {
  x$s <- columns(...)
  x
}
for (x in list(data.frame(z1 = z1), data.frame(f = fn, z1 = z1))) {
  agrees(beside(x, a = a, b = 2.54 * a), yt)
  shifted <- beside(x, m = stored, a = a, b = 2.54 * a + 32 * stored)
  agrees(shifted, yt, tolerance = 1e-8)
  agrees(beside(x, start = start, end = end), yt, tolerance = 1e-8)
  timed <- beside(x, start = start, end = end, m = stored)
  agrees(timed, yt, wm, tolerance = 1e-8)
}
# more columns than rows: the rows' crossproduct, at rank below n and at n
wide <- Matrix::rsparsematrix(60L, 20L, 0.3)
mix <- Matrix::sparseMatrix(
  i = as.vector(replicate(90L, sample.int(20L, 2L))),
  j = rep(seq_len(90L), each = 2L),
  x = runif(180L, 0.5, 2),
  dims = c(20L, 90L)
)
yw <- rnorm(60L)
spanning <- methods::as(wide %*% mix, "CsparseMatrix")
expect_equal(
  sparseSigma(spanning, y = yw, weights = NULL, offset = NULL),
  lmSigma(as.matrix(spanning), yw),
  tolerance = 1e-10
)
full <- Matrix::rsparsematrix(60L, 90L, 0.2)
expect_true(is.na(sparseSigma(full, y = yw, weights = NULL, offset = NULL)))

# --- a factor's missing rows are never stored, and the Cholesky sees only
# the columns left beside the widest factor ---

big <- factor(replace(
  sample.int(500L, 2000L, TRUE),
  sample.int(2000L, 40L),
  NA
))
bigFrame <- data.frame(f = big, v = rnorm(2000L))
bigDesign <- dbarts:::sparseSigmaDesign(dbartsData(bigFrame, rnorm(2000L))@x)
expect_equal(length(bigDesign$X@x), 2L * 2000L - 40L)
expect_equal(length(bigDesign$imputed[["1"]]$rows), 40L)
sizes <- integer()
original <- dbarts:::sigmaPivotedCholesky
recording <- function(S, ...) {
  sizes <<- c(sizes, ncol(S))
  original(S, ...)
}
assignInNamespace("sigmaPivotedCholesky", recording, "dbarts")
observed <- try(sparseSigma(
  bigFrame,
  y = bigFrame$v,
  weights = NULL,
  offset = NULL
))
assignInNamespace("sigmaPivotedCholesky", original, "dbarts")
expect_identical(sizes, 2L)
expect_true(is.numeric(observed))
agrees(bigFrame, bigFrame$v + rnorm(2000L))
# a sparseFactor's reference level, first or not, gets its indicator
for (reference in c("a", "d")) {
  x <- data.frame(x1 = x1)
  x$f <- sparseFactor(f6, reference = reference)
  expect_equal(
    sparseSigma(x, y = y, weights = NULL, offset = NULL),
    summary(lm(y ~ f6 + x1))$sigma,
    tolerance = 1e-10
  )
}
# LAPACK never tests its first pivot: a factor alone leaves an intercept the
# factor already spans, and a value near 1.7e9 with a spread of 3600 on two
# levels' rows is no column at this tolerance
set.seed(12)
f4 <- factor(sample.int(4L, 400L, TRUE))
onTwo <- f4 %in% c("2", "3")
v <- ifelse(onTwo, 1.7e9 + runif(400L, 0, 3600), 0)
y4 <- rnorm(4L)[f4] +
  ifelse(onTwo, (v - 1.7e9) / 600, 0) +
  rnorm(400L, sd = 0.5)
alone <- data.frame(f = f4)
agrees(alone, y4, rexp(400L))
expect_equal(
  sparseSigma(beside(alone, v = v), y = y4, weights = NULL, offset = NULL),
  lmSigma(indicatorForm(alone), y4),
  tolerance = 1e-10
)

# --- stored zeros, a first level's rows, the rank band's edges, units ---

# a time beside the indicator of the rows it is missing on, with ten zeros
# stored in each: the pair is found by the nonzero entries
set.seed(7)
absent <- runif(m) < 0.1
time <- ifelse(absent, 0, 1.7e9 + runif(m, 0, 3600))
others <- Matrix::rsparsematrix(m, 5L, 0.2)
ys <- ifelse(absent, 0, (time - 1.7e9) / 600) + rnorm(m, sd = 0.5)
canonical <- cbind(columns(time = time, absent = absent), others)
triplets <- methods::as(canonical, "TsparseMatrix")
zeroRows <- c(which(absent)[1:10], which(!absent)[1:10])
withZeros <- Matrix::sparseMatrix(
  i = c(triplets@i + 1L, zeroRows),
  j = c(triplets@j + 1L, rep(1:2, each = 10L)),
  x = c(triplets@x, numeric(20L)),
  dims = dim(canonical)
)
expect_equal(sum(withZeros@x == 0), 20L)
for (x in list(canonical, withZeros)) {
  expect_equal(
    sparseSigma(x, y = ys, weights = NULL, offset = NULL),
    lmSigma(as.matrix(canonical), ys),
    tolerance = 1e-8
  )
}
onFirst <- ifelse(f30 == "1", 5e6 + rnorm(m, sd = 10), 0)
yf <- ym + (onFirst - 5e6 * (f30 == "1")) / 10
agrees(beside(data.frame(f = f30, z1 = z1), v = onFirst), yf, tolerance = 1e-8)
# a near-copy at 1e-4 of its norm stays, as in lm.fit; at 1e-6 it is dropped,
# where lm.fit keeps it and "dense" follows lm.fit
set.seed(51)
base <- Matrix::rsparsematrix(m, 6L, 0.4)
for (eps in c(1e-4, 1e-6)) {
  copy <- base[, 1L] + eps * rnorm(m) * (runif(m) < 0.4)
  yc <- (base[, 1L] - copy) / eps + rnorm(m, sd = 0.5)
  x <- cbind(base, columns(copy = copy))
  kept <- lmSigma(as.matrix(if (eps < 1e-5) base else x), yc)
  expect_equal(sigmaOf(x, yc, sigest = "sparse"), kept, tolerance = 1e-6)
  expect_identical(sigmaOf(x, yc, sigest = "dense"), lmSigma(as.matrix(x), yc))
  expect_true(
    abs(lmSigma(as.matrix(base), yc) / lmSigma(as.matrix(x), yc) - 1) > 0.1
  )
}
scaled <- cbind(base, columns(huge = 1e160 * base[, 1L]^2))
expect_equal(
  sparseSigma(scaled, y = ys, weights = NULL, offset = NULL),
  lmSigma(cbind(as.matrix(base), as.numeric(base[, 1L]^2)), ys),
  tolerance = 1e-8
)
# a column matching another's stored rows, or the rows it leaves, in count
# and in the sum of the row numbers but not row for row is not its partner
value <- replace(numeric(40L), c(2L, 3L, 6L), c(5, 7, 9))
near <- replace(numeric(40L), c(1L, 4L, 6L), 1)
unpaired <- columns(value = value, near = near, far = 1 - near)
yu <- value + rnorm(40L)
expect_equal(
  sparseSigma(unpaired, y = yu, weights = NULL, offset = NULL),
  lmSigma(as.matrix(unpaired), yu),
  tolerance = 1e-10
)

# --- no regression defined: the sd of the response, no warning, a line
# under verbose and in the summary, the fit recording it ---

levelPerRow <- data.frame(f = factor(seq_len(30L)), z = rnorm(30L))
y30 <- rnorm(30L)
fellBack <- "starting sigma is the sd of the response \\(.*\\): the linear model on"
estimating <- "estimating the starting sigma by a (dense|sparse) linear regression"
group <- rep(1:3, 10L)
grouped <- function(verbose, n.threads = 1L, ...) {
  rbart_vi(
    y30 ~ f + z,
    levelPerRow,
    group.by = group,
    n.trees = 3L,
    n.samples = 2L,
    n.burn = 1L,
    n.thin = 1L,
    n.chains = 2L,
    n.threads = n.threads,
    verbose = verbose,
    ...
  )$sigest
}
fitters <- list(
  dbarts = function(verbose, ...) {
    control@verbose <- verbose
    dbarts(levelPerRow, y30, control = control, ...)$data@sigma
  },
  bart = function(verbose, ...) bartOf(levelPerRow, y30, verbose, ...)$sigest,
  bartBT = function(verbose, ...) {
    bartBTOf(levelPerRow, y30, verbose, ...)$sigest
  },
  rbart_vi = grouped,
  xbart = function(verbose, ...) {
    xbartOf(levelPerRow, y30, verbose = verbose, ...)
  }
)
for (door in names(fitters)) {
  for (rule in c("dense", "sparse")) {
    loud <- signals(value <- fitters[[door]](TRUE, sigest = rule))
    # bartBT prints its lines with the rest of its output
    lines <- if (door == "bartBT") {
      capture.output(fitters[[door]](TRUE, sigest = rule))
    } else {
      loud$message
    }
    expect_identical(loud$warning, character(), info = door)
    expect_equal(sum(grepl(fellBack, lines)), 1L, info = door)
    expect_equal(sum(grepl(estimating, lines)), 1L, info = door)
    expect_equal(sum(grepl(paste("a", rule, "linear"), lines)), 1L, info = door)
    if (door == "bartBT") {
      expect_false(any(grepl("starting sigma", loud$message)))
    }
    if (door != "xbart") {
      expect_identical(value, sd(y30), info = door)
    }
    quiet <- signals(fitters[[door]](FALSE, sigest = rule))
    expect_identical(quiet, silent, info = door)
  }
}
if (at_home()) {
  threaded <- unlist(signals(grouped(TRUE, n.threads = 2L)))
  expect_false(any(grepl("starting sigma", threaded)))
}
# a caller's handler that muffles messages finds nothing it cannot muffle
expect_identical(
  withCallingHandlers(
    fitters$bart(FALSE),
    message = function(m) invokeRestart("muffleMessage")
  ),
  sd(y30)
)
recorded <- bartOf(levelPerRow, y30)
expect_true(recorded$sigest.fallback)
expect_true(bartBTOf(levelPerRow, y30)$sigest.fallback)
summaryLine <- paste0(
  "(Starting sigma: the sd of the response; the linear model had no ",
  "residual degrees of freedom)"
)
expect_true(summaryLine %in% capture.output(print(summary(recorded))))
partial <- pdbart(
  levelPerRow,
  y30,
  xind = "z",
  pl = FALSE,
  n.trees = 3L,
  n.samples = 2L,
  n.burn = 1L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)
expect_identical(partial$sigest, sd(y30))
for (fit in list(
  bartOf(frame, y),
  bartOf(levelPerRow, y30, sigest = 1),
  replace(recorded, "sigest.fallback", list(NULL))
)) {
  expect_null(fit$sigest.fallback)
  expect_false(summaryLine %in% capture.output(print(summary(fit))))
}

# --- an infinite entry is refused before any routine runs ---

infinite <- sparse
infinite[5L, 3L] <- Inf
colnames(infinite) <- paste0("c", 1:4)
infiniteVector <- withVector
infiniteVector$s[7L] <- Inf
building <- c(design = "startingSigmaDesign")
for (x in list(infinite, infiniteVector)) {
  for (rule in c("auto", "dense", "sparse")) {
    expect_error(
      traced(sigmaOf(x, y, sigest = rule), building),
      "starting estimate of sigma: predictor '(c3|s)' has infinite values"
    )
    expect_identical(calls$counts[["design"]], 0L)
  }
}

# --- under a fixed residual scale no estimate is made: the slot holds the
# fixed value and an infinite predictor value is no obstacle ---

estimating <- c(design = "startingSigmaDesign", all = "estimateStartingSigma")
for (rule in list(NULL, "dense", "sparse")) {
  none <- traced(
    value <- sigmaOf(
      infinite,
      y,
      family = gaussian(sigma = fixed(0.49)),
      sigest = rule
    ),
    estimating
  )
  expect_identical(value, sqrt(0.49))
  expect_identical(unname(none), c(0L, 0L))
}
fixedFit <- function(x = frame, response = y, verbose = FALSE, ...) {
  bartOf(
    x,
    response,
    verbose,
    family = gaussian(sigma = fixed(0.49)),
    seed = 11L,
    ...
  )
}
expect_identical(fixedFit()$sigest, sqrt(0.49))
expect_null(fixedFit(levelPerRow, y30)$sigest.fallback)
expect_identical(signals(fixedFit(levelPerRow, y30, TRUE)), silent)
expect_identical(
  signals(xbartOf(infinite, y, family = gaussian(sigma = fixed(0.49)))),
  silent
)
# a number beside it is accepted where the two agree, the draws those of a
# fit given none, and refused where they differ
onceKeys[["tombstone.sigest.fixed"]] <- NULL
expect_message(
  agreeing <- fixedFit(sigest = sqrt(0.49)),
  "has no effect under a fixed residual scale \\(it agrees with"
)
expect_identical(agreeing$yhat.train, fixedFit()$yhat.train)
expect_error(
  fixedFit(sigest = 0.5),
  "'sigest' has no effect under a fixed residual scale: sigma ="
)

# --- xbart: each fold runs the all-rows routine on its own rows, or the
# caller's function ---

# each fold's estimate before it is floored, and the responses it holds out
foldSigmas <- function(x, y, ...) {
  seen <- numeric()
  original <- dbarts:::floorSigmaEstimate
  recording <- function(sigma, residual) {
    seen <<- c(seen, sigma)
    original(sigma, residual)
  }
  assignInNamespace("floorSigmaEstimate", recording, "dbarts")
  on.exit(assignInNamespace("floorSigmaEstimate", original, "dbarts"))
  held <- list()
  loss <- function(y.test, testSamples, weights) {
    held[[length(held) + 1L]] <<- y.test
    0
  }
  xbartOf(x, y, loss = loss, ...)
  list(sigma = seen, held = held)
}
foldFrame <- both[c("f", "z1")]
for (rule in c("dense", "sparse")) {
  folds <- foldSigmas(foldFrame, ym, sigest = rule)
  # the estimate on all rows, then one per fold
  expect_equal(length(folds$sigma), 4L)
  for (k in 1:3) {
    rows <- which(!(ym %in% folds$held[[k]]))
    expect_equal(
      folds$sigma[k + 1L],
      lmSigma(indicatorForm(foldFrame), ym, rows = rows),
      tolerance = 1e-10
    )
  }
}
# where all rows leave no regression, no fold attempts one
expect_equal(suppressMessages(traced(xbartOf(levelPerRow, y30)))[["lm"]], 1L)
given <- list()
recorder <- function(x, y, weights, offset) {
  given[[length(given) + 1L]] <<- list(x, y, weights, offset)
  2 * sd(y)
}
functionFolds <- foldSigmas(frame, y, sigest = recorder)
expect_equal(length(functionFolds$sigma), 0L)
expect_equal(length(given), 3L)
for (k in 1:3) {
  rows <- which(!(y %in% functionFolds$held[[k]]))
  expect_identical(given[[k]][[2L]], y[rows])
  expect_equal(
    given[[k]][[1L]],
    indicatorForm(frame)[rows, ],
    check.attributes = FALSE
  )
  expect_null(given[[k]][[3L]])
  expect_null(given[[k]][[4L]])
}
given <- list()
xbartOf(
  withVector,
  y,
  weights = wm[1:300],
  offset = om[1:300],
  sigest = recorder
)
expect_inherits(given[[1L]][[1L]], "dgCMatrix")
expect_identical(lengths(given[[1L]][2:4]), rep(200L, 3L))
# the list form calls the function from the environment given
listed <- list(
  function(x, y, weights, offset) get("scale", parent.frame()) * sd(y),
  list2env(list(scale = 3))
)
tripled <- function(x, y, weights, offset) 3 * sd(y)
expect_identical(
  xbartOf(frame, y, sigest = listed),
  xbartOf(frame, y, sigest = tripled)
)
expect_false(identical(xbartOf(frame, y, sigest = tripled), xbartOf(frame, y)))
# a function that never reads x builds no design
expect_identical(
  traced(xbartOf(frame, y, sigest = tripled), building)[["design"]],
  0L
)
# one that draws random numbers draws on the stream its unit runs on
drawing <- function(x, y, weights, offset) sd(y) * (1 + sample.int(5L, 1L) / 10)
set.seed(8)
first <- xbartOf(frame, y, seed = NULL, sigest = drawing)
set.seed(8)
expect_identical(xbartOf(frame, y, seed = NULL, sigest = drawing), first)
set.seed(9)
before <- .Random.seed
xbartOf(frame, y, sigest = drawing)
expect_identical(.Random.seed, before)
# and no worker seeds it: under a seed its draws continue the stream the
# replication's split was drawn from
drawn <- integer()
drawer <- function(x, y, weights, offset) {
  drawn <<- c(drawn, sample.int(1000000L, 1L))
  sd(y)
}
xbartOf(frame, y, sigest = drawer)
set.seed(3L)
set.seed(sample.int(.Machine$integer.max, 4L)[1L])
split <- sample.int(n)
expect_identical(drawn, replicate(3L, sample.int(1000000L, 1L)))
for (value in list(NA_real_, 0, Inf, "1", c(1, 2))) {
  expect_error(
    xbartOf(frame, y, sigest = function(x, y, weights, offset) value),
    "'sigest' function must return one positive finite number"
  )
}
expect_error(
  xbartOf(frame, y, sigest = function(y) 1),
  "supplied sigest function must take exactly four arguments"
)
expect_error(
  xbartOf(frame, y, sigest = recorder, family = gaussian(sigma = fixed(0.49))),
  "'sigest' has no effect under a fixed residual scale"
)
given <- list()
xbartOf(frame, binary, sigest = recorder)
expect_equal(length(given), 0L)
