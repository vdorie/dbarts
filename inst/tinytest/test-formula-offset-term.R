# A formula's offset() terms count as in lm: added to the 'offset' argument
# at fit time, evaluated on 'test' and on predict's newdata. Data-dependent
# terms (poly(), ns(), scale()) are rebuilt on new rows from the training
# values, as predict.lm rebuilds them. An explicit offset.test is length-checked
# and a bare name in it is read from 'test' first.

fitArgs <- list(
  n.trees = 15L,
  n.samples = 20L,
  n.burn = 20L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE,
  seed = 3L,
  keepTrees = TRUE
)
# the arguments stay unevaluated, so a bare column name reaches bart as written
fitWith <- function(...) {
  call <- match.call()
  call[[1L]] <- quote(bart)
  for (name in setdiff(names(fitArgs), names(call))) {
    call[[name]] <- fitArgs[[name]]
  }
  eval(call, parent.frame())
}

set.seed(1L)
n <- 60L
d <- data.frame(a = runif(n), c = runif(n), o = 10 * rnorm(n))
d$y <- sin(3 * d$a) + d$o + rnorm(n, 0, 0.1)
te <- data.frame(a = c(0.2, 0.5, 0.8), c = c(0.1, 0.2, 0.3), o = c(50, -50, 0))

# ---- fit time: the term is the offset, alone or summed with the argument

f.term <- fitWith(formula = y ~ a + offset(o), data = d)
f.arg <- fitWith(formula = y ~ a, data = d, offset = d$o)
expect_identical(f.term$yhat.train, f.arg$yhat.train)
expect_identical(f.term$sigma, f.arg$sigma)
expect_identical(f.term$fit$data@offset, d$o)

f.scalar <- fitWith(formula = y ~ a + offset(o), data = d, offset = 2)
expect_identical(f.scalar$fit$data@offset, d$o + 2)
f.vector <- fitWith(formula = y ~ a + offset(o), data = d, offset = rep(1, n))
expect_identical(f.vector$fit$data@offset, d$o + 1)
expect_identical(dbartsData(y ~ a + offset(o), d)@offset, d$o)

# a missing value in the term is refused as one in the argument is
d.na <- d
d.na$o[3L] <- NA
expect_error(
  fitWith(formula = y ~ a + offset(o), data = d.na),
  pattern = "'offset' contains missing values"
)

# ---- test and predict: the term is evaluated on the new rows

f.test <- fitWith(formula = y ~ a + offset(o), data = d, test = te)
expect_identical(f.test$fit$data@offset.test, te$o)
p.test <- predict(f.test, te)
expect_equal(unname(p.test), unname(f.test$yhat.test))
expect_true(cor(colMeans(p.test), te$o) > 0.99)
# an offset given to predict adds to the term
expect_equal(
  unname(predict(f.test, te, offset = 10)),
  unname(p.test) + 10
)
# the fit's offset argument is applied at predict as well, as in lm
f.both <- fitWith(formula = y ~ a + offset(o), data = d, offset = 5, test = te)
expect_identical(f.both$fit$data@offset.test, te$o + 5)
expect_equal(unname(predict(f.both, te)), unname(f.both$yhat.test))
expect_error(
  predict(f.test, te[c("a", "c")]),
  pattern = "missing variable required by the formula's offset\\(\\) term: 'o'"
)
expect_error(
  predict(f.test, as.matrix(te)),
  pattern = "offset\\(\\) term, which is evaluated on 'newdata'"
)
expect_error(
  fitWith(formula = y ~ a + offset(o), data = d, test = te[c("a", "c")]),
  pattern = "missing variable required by the formula's offset\\(\\) term"
)

# ---- families: one refusing an offset refuses the term; others carry it

d$k <- factor(sample(c("p", "q", "r"), n, TRUE))
expect_error(
  fitWith(formula = k ~ a + offset(o), data = d, family = "multinomial"),
  pattern = "requires an n x K matrix \"offset\""
)
d$pos <- ifelse(d$a < 0.3, 0, exp(d$a))
expect_error(
  fitWith(formula = pos ~ a + offset(o), data = d, family = "hurdle.lognormal"),
  pattern = "does not support 'offset'/'offset.test'"
)
d$ex <- runif(n, 1, 5)
d$cnt <- rpois(n, d$ex * exp(d$a))
expect_identical(
  fitWith(
    formula = cnt ~ a + offset(log(ex)),
    data = d,
    family = "nbinom"
  )$yhat.train,
  fitWith(
    formula = cnt ~ a,
    data = d,
    offset = log(d$ex),
    family = "nbinom"
  )$yhat.train
)
if (requireNamespace("survival", quietly = TRUE)) {
  d$time <- sample(1:4, n, TRUE)
  d$status <- rbinom(n, 1L, 0.7)
  d$lo <- rnorm(n, 0, 0.3)
  # the person-period expansion carries the term as it carries the argument
  f.haz.term <- fitWith(
    formula = survival::Surv(time, status) ~ a + offset(lo),
    data = d,
    family = "hazard"
  )
  f.haz.arg <- fitWith(
    formula = survival::Surv(time, status) ~ a,
    data = d,
    offset = d$lo,
    family = "hazard"
  )
  expect_identical(f.haz.term$yhat.train, f.haz.arg$yhat.train)
  expect_identical(f.haz.term$fit$data@offset, f.haz.arg$fit$data@offset)
}

# ---- data-dependent terms rebuild from the training values

d$y2 <- 3 * d$c^2 + rnorm(n, 0, 0.1)
nd1 <- data.frame(c = c(0.5, 0.9, 0.1, 0.3))
nd2 <- data.frame(c = c(0.5, 0.51, 0.52, 0.53))
# a ':' in a label reads as an interaction, so ns is bound here rather than
# written splines::ns
ns <- splines::ns
for (rhs in c("poly(c, 2)", "scale(c)", "ns(c, 3)")) {
  f.basis <- fitWith(
    formula = as.formula(paste("y2 ~", rhs)),
    data = d,
    test = nd2
  )
  p1 <- predict(f.basis, nd1)[, 1L]
  p2 <- predict(f.basis, nd2)[, 1L]
  expect_equal(unname(p1), unname(p2))
  expect_equal(
    unname(predict(f.basis, nd1[1L, , drop = FALSE])[, 1L]),
    unname(p1)
  )
  expect_equal(unname(f.basis$yhat.test[, 1L]), unname(p2))
}
# a few training rows replay their own training fits
f.poly <- fitWith(formula = y2 ~ poly(c, 2), data = d)
expect_equal(
  unname(predict(f.poly, d[1:3, ])),
  unname(f.poly$yhat.train[, 1:3])
)

# ---- an explicit offset.test

te5 <- data.frame(a = runif(5L), c = runif(5L), o = 100)
expect_error(
  fitWith(
    formula = y ~ a,
    data = d,
    offset = o,
    test = te5,
    offset.test = c(7, 8, 9)
  ),
  pattern = "'offset.test' must have the same number of rows as 'test'"
)
expect_error(
  fitWith(formula = y ~ a, data = d, offset = o, test = te5, offset.test = d$o),
  pattern = "'offset.test' must have the same number of rows as 'test'"
)
# a bare name is the test set's own column before the training data's
f.named <- fitWith(
  formula = y ~ a,
  data = d,
  offset = o,
  test = te5,
  offset.test = o
)
expect_identical(f.named$fit$data@offset.test, rep(100, 5L))
f.expr <- fitWith(
  formula = y ~ a,
  data = d,
  offset = o,
  test = te5,
  offset.test = o / 2
)
expect_identical(f.expr$fit$data@offset.test, rep(50, 5L))
# A name in the caller's offset.test expression resolves
# where model.frame resolves one - in the data, then the formula's
# environment, or the caller's frame without a formula - and never in this
# package's own frames, whose locals ('x', 'data', 'offset') and namespace
# lie on the evaluator's own enclosure chain. Each case runs inside a
# function whose local takes a name those frames also bind, so the test
# needs no change to the global environment, which a conflicting global
# would otherwise be the way to show.
shadowedFit <- function() {
  x <- rep(7, 5L)
  fitWith(
    formula = y ~ a,
    data = d,
    offset = o,
    test = te5[, c("a", "o")],
    offset.test = o * 0 + x
  )
}
expect_identical(shadowedFit()$fit$data@offset.test, rep(7, 5L))
shadowedMatrixFit <- function() {
  x <- rep(3, 5L)
  bart(
    as.matrix(d["a"]),
    d$y,
    test = as.matrix(te5["a"]),
    offset = 1,
    offset.test = x + 0,
    n.trees = 5L,
    n.samples = 5L,
    n.burn = 5L,
    n.chains = 1L,
    n.threads = 1L,
    verbose = FALSE,
    keepTrees = TRUE
  )
}
expect_identical(shadowedMatrixFit()$fit$data@offset.test, rep(3, 5L))
rm(shadowedFit, shadowedMatrixFit)
expect_identical(
  fitWith(
    formula = y ~ a,
    data = d,
    offset = o,
    test = te5,
    offset.test = 3
  )$fit$data@offset.test,
  rep(3, 5L)
)

# ---- offset.test naming 'offset' beside a term: the argument's own share,
# evaluated on test, plus the term evaluated there - never the training term

f.named.term <- fitWith(
  formula = y ~ a + offset(o),
  data = d,
  test = te,
  offset.test = offset
)
expect_identical(f.named.term$fit$data@offset.test, te$o)
f.half <- fitWith(
  formula = y ~ a + offset(o / 2),
  data = d,
  offset = o / 2,
  test = te,
  offset.test = offset
)
expect_identical(f.half$fit$data@offset.test, te$o)
f.plus <- fitWith(
  formula = y ~ a + offset(o),
  data = d,
  test = te,
  offset.test = offset + 1
)
expect_identical(f.plus$fit$data@offset.test, te$o + 1)

# ---- a name in the term found in the formula's environment, as in lm

scale.o <- 2
gone <- 1
d.gone <- d
f.env <- fitWith(formula = y ~ a + offset(o * scale.o), data = d, test = te)
expect_identical(f.env$fit$data@offset.test, te$o * scale.o)
# at one x, rows differing only in o differ by scale.o times that
same.a <- data.frame(a = 0.5, c = 0.5, o = c(0, 1))
expect_equal(
  unname(apply(predict(f.env, same.a), 1L, diff)),
  rep(scale.o, nrow(predict(f.env, same.a)))
)
# a name that is no data column is frozen at its fit-time value, so it is not
# looked up again
before.gone <- predict(
  fitWith(formula = y ~ a + offset(o * gone), data = d),
  te
)
f.gone <- fitWith(formula = y ~ a + offset(o * gone), data = d.gone)
rm(gone)
expect_equal(predict(f.gone, te), before.gone)

# ---- the offset argument is re-evaluated on new rows, as predict.lm does

d$e <- runif(n, 1, 3)
te$e <- c(1, 2, 3)
f.log <- fitWith(formula = y ~ a, data = d, offset = log(e))
f.log.term <- fitWith(formula = y ~ a + offset(log(e)), data = d)
expect_identical(f.log$yhat.train, f.log.term$yhat.train)
expect_equal(predict(f.log, te), predict(f.log.term, te))
# at one x, rows differing only in e differ by lm's offset difference
same.x <- data.frame(a = 0.5, c = 0.5, o = 0, e = c(1, 4))
lm.log <- lm(y ~ a, d, offset = log(e))
expect_equal(
  unname(apply(predict(f.log, same.x), 1L, diff)),
  rep(unname(diff(predict(lm.log, same.x))), nrow(predict(f.log, same.x)))
)
# an offset given to predict adds to it
expect_equal(predict(f.log, te, offset = 1), predict(f.log, te) + 1)
# the test-set default reads the argument on test the same way
f.log.test <- fitWith(formula = y ~ a, data = d, offset = log(e), test = te)
expect_identical(f.log.test$fit$data@offset.test, log(te$e))
expect_false(f.log.test$fit$data@testUsesRegularOffset)
expect_equal(unname(predict(f.log.test, te)), unname(f.log.test$yhat.test))
# a plain vector for the training rows cannot be evaluated on others: predict
# asks for an offset, which then stands in for it
f.vector.arg <- fitWith(formula = y ~ a, data = d, offset = d$o)
expect_error(
  predict(f.vector.arg, te),
  pattern = "the fit's 'offset' was given as 'd\\$o', which cannot be evaluated"
)
expect_equal(
  predict(f.vector.arg, te, offset = te$o),
  predict(f.arg, te, offset = te$o)
)
expect_error(
  fitWith(formula = y ~ a, data = d, offset = d$o, test = te),
  pattern = "'offset' was given as 'd\\$o', which cannot be evaluated on the rows of 'test'"
)

# ---- ordinal predict and survivalProbabilities evaluate the offset on new
# rows as predict does

d$yo <- factor(cut(d$a + d$o / 10, 3L), ordered = TRUE)
f.ordinal <- fitWith(formula = yo ~ a + offset(o / 10), data = d)
# the training rows replay the fit's own offset
expect_equal(
  unname(predict(f.ordinal, d)),
  unname(extract(f.ordinal, type = "ev")),
  tolerance = 1e-12
)
# an offset given to predict adds to the term evaluated on newdata
expect_equal(
  predict(f.ordinal, te, offset = 1, type = "bart"),
  predict(f.ordinal, transform(te, o = o + 10), type = "bart")
)
if (requireNamespace("survival", quietly = TRUE)) {
  # a subject's offset() term is evaluated on newdata and applies at every
  # period, so the training subjects come back as the fit's own curves
  expect_equal(
    unname(survivalProbabilities(f.haz.term, newdata = d)),
    unname(survivalProbabilities(f.haz.term)),
    tolerance = 1e-12
  )
  expect_equal(
    survivalProbabilities(f.haz.term, newdata = d[1:3, ], offset = 0.5),
    survivalProbabilities(
      f.haz.term,
      newdata = transform(d[1:3, ], lo = lo + 0.5)
    )
  )
  # the argument given as a training vector applies at as many rows only
  expect_equal(
    unname(survivalProbabilities(f.haz.arg, newdata = d)),
    unname(survivalProbabilities(f.haz.arg)),
    tolerance = 1e-12
  )
  expect_error(
    survivalProbabilities(f.haz.arg, newdata = d[1:3, ]),
    pattern = "give survivalProbabilities an 'offset' for them"
  )
}

rm(
  f.named.term,
  f.half,
  f.plus,
  scale.o,
  f.env,
  same.a,
  d.gone,
  f.gone,
  before.gone,
  f.log,
  f.log.term,
  same.x,
  lm.log,
  f.log.test,
  f.vector.arg,
  f.ordinal
)

rm(
  fitArgs,
  fitWith,
  n,
  d,
  te,
  f.term,
  f.arg,
  f.scalar,
  f.vector,
  d.na,
  f.test,
  p.test,
  f.both,
  nd1,
  nd2,
  ns,
  rhs,
  f.basis,
  p1,
  p2,
  f.poly,
  te5,
  f.named,
  f.expr
)
if (exists("f.haz.term")) {
  rm(f.haz.term, f.haz.arg)
}
