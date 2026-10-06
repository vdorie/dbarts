# A forest's basis formula is rebuilt at new rows the way lm rebuilds a term of
# its formula: scale(), poly(), ns() and bs() keep the centre, scale or knots
# of the training rows (stats::makepredictcall), and any other expression is
# evaluated on the rows predict is given. The reference is lm's own machinery,
# model.frame() on the training rows and then on the terms object with the new
# rows, which evaluates the term's predvars.

set.seed(23)
n <- 60L
d <- data.frame(
  a = runif(n),
  b = runif(n),
  w = 50 + 10 * rnorm(n),
  v = rexp(n),
  g = factor(sample(c("lo", "mid", "hi"), n, replace = TRUE))
)
d$y <- d$a + rnorm(n, sd = 0.3)
nd <- data.frame(
  a = c(0.2, 0.5, 0.8),
  b = c(0.9, 0.1, 0.4),
  w = unname(stats::quantile(d$w, c(0.5, 0.75, 0.95))),
  v = c(0.1, 1, 2),
  g = factor(c("hi", "lo", "hi"), levels = levels(d$g))
)
stacked <- rbind(d[1:7, names(nd)], nd)

fitBasisTerms <- function(expr, data = d, ...) {
  f <- stats::as.formula(
    paste0("y ~ a + b + forest(a + b, basis = ~ ", expr, ", n.trees = 5L)")
  )
  bart(
    f,
    data,
    ...,
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 5L,
    n.burn = 0L,
    n.samples = 3L,
    keepTrees = TRUE,
    verbose = FALSE
  )
}
plain <- function(m) matrix(as.vector(m), nrow(m))
lmBasis <- function(expr, rows) {
  mf <- stats::model.frame(stats::as.formula(paste("~", expr)), d)
  tt <- stats::terms(mf)
  plain(as.matrix(stats::model.frame(tt, rows)[[1L]]))
}
replayAt <- function(fit, rows) {
  dbarts:::replayForestBasis(fit$basis.terms[[2L]], rows, 2L)
}

hasSplines <- requireNamespace("splines", quietly = TRUE)
exprs <- c(
  "scale(w)",
  "scale(w, scale = FALSE)",
  "poly(w, 2)",
  "poly(w, 2, raw = TRUE)",
  "scale(log(w))",
  "scale(cbind(w, v))"
)
if (hasSplines) {
  exprs <- c(exprs, "splines::ns(w, 3)", "splines::bs(w, df = 4)")
}

for (expr in exprs) {
  fit <- fitBasisTerms(expr)
  expect_equal(
    plain(replayAt(fit, nd)),
    lmBasis(expr, nd),
    info = paste(expr, "- several rows")
  )
  expect_equal(
    plain(replayAt(fit, nd[2L, ])),
    lmBasis(expr, nd[2L, ]),
    info = paste(expr, "- one row")
  )
  expect_equal(
    plain(replayAt(fit, stacked)),
    lmBasis(expr, stacked),
    info = paste(expr, "- among training rows")
  )
  # the rebuilt basis at the training rows is the fitted one
  expect_equal(
    plain(replayAt(fit, d)),
    plain(fit$bases[[2L]]),
    info = paste(expr, "- training rows")
  )
}

## predictions at a row do not depend on the rows beside it
fit <- fitBasisTerms("scale(w)")
together <- predict(fit, nd)
for (i in seq_len(nrow(nd))) {
  expect_equal(
    unname(predict(fit, nd[i, ])),
    unname(together[, i, drop = FALSE])
  )
}
expect_equal(
  unname(predict(fit, stacked)[, 8:10, drop = FALSE]),
  unname(together)
)

## a basis variable that is also a predictor: the partial dependence grid is a
## set of rows of one value, which scale() of the grid could not centre
fitA <- fitBasisTerms("scale(a)")
pd <- pdbart(
  fitA,
  xind = "a",
  levs = list(c(0.2, 0.5)),
  newdata = d,
  pl = FALSE
)
byHand <- vapply(
  c(0.2, 0.5),
  function(level) {
    rows <- d
    rows$a <- level
    rowMeans(predict(fitA, rows))
  },
  numeric(3L)
)
expect_equal(plain(pd$fd[[1L]]), plain(byHand))

## a fit saved with its sampler state and read back predicts the same
fit$fit$storeState()
path <- tempfile(fileext = ".rds")
saveRDS(fit, path)
expect_equal(predict(readRDS(path), nd), together)
unlink(path)

## what is stored: the training centre and scale, written into the call
stored <- fit$basis.terms[[2L]]$predcall
expect_equal(stored$center, mean(d$w))
expect_equal(stored$scale, stats::sd(d$w))

## subset: the centre and scale are those of the rows the fit used
keep <- d$w < 55
fitS <- fitBasisTerms("scale(w)", subset = keep)
expect_equal(fitS$basis.terms[[2L]]$predcall$center, mean(d$w[keep]))
expect_equal(
  plain(replayAt(fitS, nd)),
  matrix((nd$w - mean(d$w[keep])) / stats::sd(d$w[keep]), ncol = 1L)
)

## a fit stored without the rebuilt call evaluates the expression on the rows
## it is given
old <- fit
old$basis.terms[[2L]]$predcall <- NULL
expect_equal(
  plain(replayAt(old, nd)),
  matrix(as.vector(scale(nd$w)), ncol = 1L)
)

## an expression with no makepredictcall method is evaluated on the new rows,
## as lm evaluates it
fitM <- fitBasisTerms("I(w - mean(w))")
expect_equal(plain(replayAt(fitM, nd)), lmBasis("I(w - mean(w))", nd))
expect_equal(
  plain(replayAt(fitM, nd)),
  matrix(nd$w - mean(nd$w), ncol = 1L)
)
expect_equal(plain(replayAt(fitM, nd[1L, ])), matrix(0, 1L, 1L))

## a factor basis: a level the new rows lack is not an error, and the width is
## the fit's; a level the fit never saw is refused
fitG <- fitBasisTerms("g")
expect_equal(ncol(replayAt(fitG, nd)), nlevels(d$g))
expect_equal(dim(predict(fitG, nd)), c(3L, 3L))
bad <- nd
bad$g <- factor(c("hi", "new", "lo"))
expect_error(predict(fitG, bad), "new")

## a call under arithmetic or cbind() at the top of the expression is rebuilt
## from the training rows too, each operand as model.frame rebuilds it
d$z <- rep(c(0, 1), length.out = n)
nd$z <- c(1, 0, 1)
stackedZ <- rbind(d[1:7, names(nd)], nd)
lmVars <- function(vars, rows) {
  mf <- stats::model.frame(stats::reformulate(vars), d)
  tt <- stats::terms(mf)
  stats::model.frame(tt, rows)
}
descended <- list(
  list(
    "scale(w) + scale(v)",
    function(m) m[[1L]] + m[[2L]],
    c("scale(w)", "scale(v)")
  ),
  list("scale(w) * z", function(m) m[[1L]] * m$z, c("scale(w)", "z")),
  list("2 * scale(w)", function(m) 2 * m[[1L]], "scale(w)"),
  list(
    "cbind(scale(w), scale(v))",
    function(m) cbind(m[[1L]], m[[2L]]),
    c("scale(w)", "scale(v)")
  )
)
for (item in descended) {
  fit <- fitBasisTerms(item[[1L]])
  for (rows in list(nd, nd[2L, ], stackedZ)) {
    expect_equal(
      plain(replayAt(fit, rows)),
      plain(as.matrix(item[[2L]](lmVars(item[[3L]], rows)))),
      info = item[[1L]]
    )
  }
  expect_equal(
    plain(replayAt(fit, d)),
    plain(fit$bases[[2L]]),
    info = paste(item[[1L]], "- training rows")
  )
}

## two basis forests: each term stores the call of its own basis
fit2 <- bart(
  y ~ a +
    b +
    forest(a + b, basis = ~ scale(w), n.trees = 5L) +
    forest(a, basis = ~ scale(v), n.trees = 5L),
  d,
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 5L,
  n.burn = 0L,
  n.samples = 3L,
  keepTrees = TRUE,
  verbose = FALSE
)
expect_equal(fit2$basis.terms[[2L]]$predcall$center, mean(d$w))
expect_equal(fit2$basis.terms[[3L]]$predcall$center, mean(d$v))
expect_equal(unname(predict(fit2, d[1:5, ])), unname(predict(fit2, d)[, 1:5]))

## the stored call is evaluated in the formula's environment: a function that
## exists only where the formula was written is found at predict
fitLocal <- function() {
  half <- function(x) x / 2
  bart(
    y ~ a + b + forest(a + b, basis = ~ scale(half(w)), n.trees = 5L),
    d,
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 5L,
    n.burn = 0L,
    n.samples = 3L,
    keepTrees = TRUE,
    verbose = FALSE
  )
}
fitL <- fitLocal()
expect_equal(
  plain(replayAt(fitL, nd)),
  matrix((nd$w - mean(d$w)) / stats::sd(d$w), ncol = 1L)
)

## what the descent does not enter is evaluated on the rows given
fitI <- fitBasisTerms("I(scale(w))")
expect_equal(
  plain(replayAt(fitI, nd)),
  matrix(as.vector(scale(nd$w)), ncol = 1L)
)

## the rows a missing response drops still enter the centre, as in lm
dNA <- d
dNA$y[1:6] <- NA
fitNA <- fitBasisTerms("scale(w)", data = dNA)
expect_equal(fitNA$basis.terms[[2L]]$predcall$center, mean(d$w))

## pd2bart rebuilds the basis at its grid rows
pd2 <- pd2bart(
  fitA,
  xind = c("a", "b"),
  levs = list(c(0.2, 0.5), c(0.3, 0.6)),
  newdata = d,
  pl = FALSE
)
byHand <- vapply(
  list(c(0.2, 0.3), c(0.5, 0.3), c(0.2, 0.6), c(0.5, 0.6)),
  function(level) {
    rows <- d
    rows$a <- level[1L]
    rows$b <- level[2L]
    mean(predict(fitA, rows))
  },
  numeric(1L)
)
expect_equal(as.vector(colMeans(pd2$fd)), byHand)
