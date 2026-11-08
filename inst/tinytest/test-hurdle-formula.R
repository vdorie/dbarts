# The hurdle's front door reads its inputs once, into one data object, and cuts
# the two parts from it by rows: a formula, a data object and the matrix
# interface give the same fit, both parts see the same rows, columns and
# coding, and a row the na.action drops is gone from both.

set.seed(5101L)
n <- 90L
d <- data.frame(a = runif(n), b = runif(n))
d$f <- factor(sample(c("u", "v", "w"), n, TRUE), levels = c("u", "v", "w", "z"))
occupied <- rbinom(n, 1L, pnorm(0.8 * d$a - 0.3)) == 1L
d$y <- 0
d$y[occupied] <- exp(rnorm(sum(occupied), 0.5 + 0.5 * d$b[occupied], 0.4))
d$y[1L] <- 0
d$y[2L] <- 1.5
# a level that occurs on zero rows only, and one that occurs on none
zeroRows <- which(d$y == 0)
d$f[zeroRows[2:4]] <- "z"
d$f <- factor(as.character(d$f), levels = c("u", "v", "w", "z", "none"))
expect_true(all(d$y[d$f == "z"] == 0))
expect_equal(sum(d$f == "none"), 0L)

fitArgs <- list(
  n.samples = 8L,
  n.burn = 6L,
  n.trees = 6L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE,
  keepTrees = TRUE,
  seed = 5102L
)
fitHurdle <- function(formula, ...) {
  do.call(
    bart,
    c(list(formula, family = "hurdle.lognormal"), list(...), fitArgs)
  )
}
channels <- function(fit) {
  list(
    ev = extract(fit, type = "ev"),
    prob = extract(fit, type = "prob"),
    bart = extract(fit, type = "bart"),
    loglik = extract(fit, type = "loglik"),
    sigma = extract(fit, type = "sigma")
  )
}

# --- a formula, a data object and the matrix interface are one fit -----------
xFrame <- d[, c("a", "b", "f")]
fitMatrix <- fitHurdle(xFrame, d$y)
fitFormula <- fitHurdle(y ~ a + b + f, data = d)
fitObject <- fitHurdle(dbartsData(y ~ a + b + f, d))

expect_inherits(fitFormula, "bartHurdle")
expect_identical(channels(fitFormula), channels(fitMatrix))
expect_identical(channels(fitObject), channels(fitMatrix))
expect_identical(fitFormula$y, fitMatrix$y)

# each part is the standalone fit it is called: the same design, rows and draws
for (part in c("zero", "positive")) {
  expect_identical(
    fitFormula[[part]]$yhat.train,
    fitMatrix[[part]]$yhat.train
  )
}
expect_identical(
  fitFormula$positive$yhat.test,
  fitMatrix$positive$yhat.test
)

# an all-numeric formula gives the matrix fit on the same design
xNumeric <- as.matrix(d[, c("a", "b")])
values <- function(fit) lapply(channels(fit), unname)
expect_identical(
  values(fitHurdle(y ~ a + b, data = d)),
  values(fitHurdle(xNumeric, d$y))
)

# a formula that names no data reads its variables where it was written
local({
  a <- d$a
  b <- d$b
  y <- d$y
  expect_identical(
    values(fitHurdle(y ~ a + b)),
    values(fitHurdle(xNumeric, d$y))
  )
})

# --- both parts code the factor against the same levels ----------------------
designLevels <- function(part) {
  attr(part$fit$data@x, "factor.levels")
}
expect_identical(
  designLevels(fitFormula$positive),
  designLevels(fitFormula$zero)
)
expect_identical(designLevels(fitFormula$zero)[[3L]], levels(d$f))

# the zero-rows level is coded on the positive part's test rows: the stored
# draws at those rows are the replay of the trees at that level
zeroLevelRows <- which(d$f == "z")
stored <- extract(fitFormula, type = "bart")[, zeroLevelRows, drop = FALSE]
expect_false(anyNA(stored))
replayed <- predict(fitFormula, d[zeroLevelRows, ], type = "bart")
expect_equal(unname(replayed), unname(stored), tolerance = 1e-8)
expect_equal(
  unname(predict(fitFormula, d[zeroLevelRows, ], type = "prob")),
  unname(extract(fitFormula, type = "prob")[, zeroLevelRows, drop = FALSE]),
  tolerance = 1e-8
)

# --- row names ---------------------------------------------------------------
named <- d
rownames(named) <- paste0("r", seq_len(n))
fitNamed <- fitHurdle(y ~ a + b + f, data = named)
fitNamedMatrix <- fitHurdle(named[, c("a", "b", "f")], named$y)
expect_identical(names(fitted(fitNamed)), rownames(named))
expect_identical(names(fitted(fitNamed)), names(fitted(fitNamedMatrix)))
expect_identical(colnames(extract(fitNamed, type = "ev")), rownames(named))
expect_identical(names(residuals(fitNamed)), rownames(named))
expect_identical(fitNamed$row.names.train, fitNamedMatrix$row.names.train)
expect_identical(
  fitNamed$positive$row.names.train,
  rownames(named)[named$y > 0]
)

# --- na.action: the rows it drops are gone from both parts --------------------
dNA <- d
positiveRows <- which(dNA$y > 0)
dNA$a[positiveRows[3L]] <- NA
dNA$b[positiveRows[5L]] <- NA
dNA$y[7L] <- NA
dNA$f[which(dNA$y > 0)[8L]] <- NA
for (action in list(na.omit, na.exclude)) {
  viaFormula <- fitHurdle(y ~ a + b + f, data = dNA, na.action = action)
  viaMatrix <- fitHurdle(dNA[, c("a", "b", "f")], dNA$y, na.action = action)
  expect_identical(viaFormula$na.action, viaMatrix$na.action)
  expect_identical(
    class(viaFormula$na.action),
    class(attr(action(dNA), "na.action"))
  )
  expect_identical(channels(viaFormula), channels(viaMatrix))
  expect_identical(length(viaFormula$y), sum(stats::complete.cases(dNA)))
  expect_identical(
    length(viaFormula$zero$y) - length(viaFormula$positive$y),
    sum(viaFormula$y == 0)
  )
  expect_identical(
    length(fitted(viaFormula)),
    if (identical(action, na.exclude)) n else sum(stats::complete.cases(dNA))
  )
}
# the default keeps a missing predictor where the positive rows hold it
expect_inherits(
  fitHurdle(y ~ a + b + f, data = dNA[!is.na(dNA$y), ]),
  "bartHurdle"
)
expect_error(
  fitHurdle(y ~ a + b + f, data = dNA, na.action = na.fail),
  "missing values in object"
)

# a predictor missing on zero rows only cannot be routed on the positive part,
# whichever route the rows took
dZero <- d
dZero$a[which(dZero$y == 0)[6L]] <- NA
expect_error(
  fitHurdle(y ~ a + b + f, data = dZero),
  "missing values only on rows where y is zero"
)
expect_error(
  fitHurdle(dZero[, c("a", "b", "f")], dZero$y),
  "missing values only on rows where y is zero"
)

# --- a data object is the whole data: what is written beside it is ignored ----
warned <- tryCatch(
  fitHurdle(
    dbartsData(y ~ a + b + f, d),
    na.action = na.omit
  ),
  warning = function(w) conditionMessage(w)
)
expect_identical(
  warned,
  tryCatch(
    bart(
      dbartsData(y ~ a + b + f, d),
      n.samples = 4L,
      n.burn = 4L,
      n.chains = 1L,
      n.threads = 1L,
      verbose = FALSE,
      na.action = na.omit
    ),
    warning = function(w) conditionMessage(w)
  )
)
expect_identical(
  channels(suppressWarnings(
    fitHurdle(dbartsData(y ~ a + b + f, d), subset = 1:10)
  )),
  channels(fitObject),
  info = "a subset beside a data object is ignored"
)

# --- refusals that stand -----------------------------------------------------
expect_error(
  fitHurdle(y ~ a + b + forest(f), data = d),
  "forest\\(\\)"
)
expect_error(
  fitHurdle(
    Matrix::Matrix(as.matrix(d[, c("a", "b")]), sparse = TRUE),
    d$y
  ),
  "sparse"
)
expect_error(
  bart(xNumeric, family = "hurdle.lognormal"),
  "needs a response"
)
expect_error(
  fitHurdle(y ~ a + b, data = replace(d, "y", d$y + 1)),
  "at least one exact zero"
)

# --- the zero part's terms leave out the formula's offset() term ---------------
# (a fit with an offset is refused until offsets land; a data object reaches
# the parts' data now)
dOffset <- d
dOffset$o <- runif(n)
offsetData <- dbartsData(y ~ a + b + offset(o), dOffset, keep.data = FALSE)
parts <- dbarts:::hurdlePartData(
  offsetData,
  dbarts:::splitHurdleResponse(offsetData@y)
)
expect_null(attr(attr(parts$zero@x, "terms"), "offset"))
expect_null(parts$zero@offset)
expect_null(attr(parts$zero, "offset.argument"))
expect_false(is.null(attr(attr(parts$positive@x, "terms"), "offset")))
expect_identical(parts$positive@offset, offsetData@offset[offsetData@y > 0])
