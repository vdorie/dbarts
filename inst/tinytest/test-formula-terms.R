# The forests of a formula. A forest() call at the top of a formula's '+'
# chain is one forest of the model and its first argument the predictors it
# splits on; a multiplier is its 'basis' argument and nothing else.
# Block A: a formula with no forest() call is untouched. Block B: a forest()
# is a top-level term, and one crossed with another term is refused with the
# form to write. Block C: a forest's predictors are terms of its own. Block D:
# the forest with no basis. Block E: the tree count given twice. Block F:
# every forest multiplied. Block G: what else is refused.

n <- 60L
set.seed(0)
x1 <- runif(n)
x2 <- runif(n)
x3 <- runif(n)
a <- runif(n)
b <- runif(n)
z <- rbinom(n, 1L, 0.5)
zf <- factor(sample(c("u", "v", "w"), n, replace = TRUE))
g <- sample(c("p", "q", "r"), n, replace = TRUE)
o <- rnorm(n, 0, 0.3)
w <- runif(n, 0.5, 1.5)
weirdName <- runif(n)
y <- x1 + z * (1 + x2) + rnorm(n, 0, 0.2)
d <- data.frame(
  y = y,
  x1 = x1,
  x2 = x2,
  x3 = x3,
  a = a,
  b = b,
  z = z,
  zf = zf,
  g = g,
  o = o,
  w = w,
  stringsAsFactors = FALSE
)
d[["weird name"]] <- weirdName
# only the response, three predictors and a multiplier, for '.'
small <- d[c("y", "x1", "x2", "x3", "z")]

tinyArgs <- list(
  seed = 1L,
  n.samples = 3L,
  n.burn = 3L,
  n.trees = 3L,
  n.chains = 1L,
  n.threads = 1L,
  keepTrees = FALSE,
  keepSampler = TRUE,
  verbose = FALSE
)

fit <- function(formula, ..., data = d) {
  # an argument given as NULL is left out, for the fitting function's default
  settings <- utils::modifyList(tinyArgs, list(...))
  do.call(dbarts::bart, c(list(formula = formula, data = data), settings))
}
refuses <- function(formula, pattern, ..., fixed = TRUE) {
  # nolint next: object_usage_linter. tinytest attaches expect_* at run time.
  expect_error(
    fit(formula, ...),
    pattern,
    fixed = fixed,
    info = deparse(formula)
  )
}
forestInfo <- function(result) {
  attr(result$fit$control, "bartcore.forests", exact = TRUE)
}
predictors <- function(result) colnames(result$fit$data@x)
# the terms a fit stores, which predict reads new rows by
storedTerms <- function(result) attr(result$fit$data@x, "terms")
storedFormula <- function(result) {
  paste(deparse(stats::formula(storedTerms(result))), collapse = " ")
}
# the design without the terms it was built from: two spellings of one model
# may write those differently
design <- function(result) {
  x <- result$fit$data@x
  attr(x, "terms") <- NULL
  x
}
# the columns each forest splits on, by name; every column where unrestricted
splitsOn <- function(result) {
  lapply(forestInfo(result)$vars, function(columns) {
    if (is.null(columns)) predictors(result) else predictors(result)[columns]
  })
}
# the tree count of each forest: the first forest's is the control's
treeCounts <- function(result) {
  counts <- vapply(forestInfo(result)$params, function(p) p[[1L]], 0)
  counts[1L] <- result$fit$control@n.trees
  counts
}
# two spellings of one model: the same design, bases, forests and draws
expectSameForest <- function(formulaA, formulaB, ...) {
  fitA <- fit(formulaA, ...)
  fitB <- fit(formulaB, ...)
  info <- paste(deparse(formulaA), "against", deparse(formulaB))
  # nolint start: object_usage_linter. tinytest attaches expect_* at run time.
  expect_identical(design(fitA), design(fitB), info = info)
  expect_identical(
    attr(storedTerms(fitA), "term.labels"),
    attr(storedTerms(fitB), "term.labels"),
    info = info
  )
  expect_identical(fitA$fit$data@bases, fitB$fit$data@bases, info = info)
  expect_identical(forestInfo(fitA), forestInfo(fitB), info = info)
  expect_identical(fitA$yhat.train, fitB$yhat.train, info = info)
  expect_identical(fitA$sigma, fitB$sigma, info = info)
  # nolint end
  invisible(fitA)
}

## --- Block A: a formula with no forest() call is untouched -----------------
expect_true(is.null(dbarts:::walkFormulaTerms(y ~ a + b)))
expect_true(is.null(dbarts:::walkFormulaTerms(y ~ .)))
expect_true(is.null(dbarts:::walkFormulaTerms(y ~ . - z)))
# a ':' whose neither operand is a forest() call is not a term at all
expect_true(is.null(dbarts:::walkFormulaTerms(y ~ a + z:b)))
expect_true(is.null(dbarts:::walkFormulaTerms(y ~ a - 1)))
expect_true(is.null(dbarts:::walkFormulaTerms(y ~ a + offset(o))))
expect_true(is.null(dbarts:::walkFormulaTerms(as.formula("y ~ `weird name`"))))
# a call with an argument left empty is read without evaluating it
expect_true(is.null(dbarts:::walkFormulaTerms(y ~ m[, 1L] + a)))

expect_silent(fit(y ~ a + b))
expect_silent(fit(y ~ .))
expect_silent(fit(y ~ . - z))
expect_silent(fit(y ~ a - 1))
expect_silent(fit(y ~ a + offset(o)))
expect_silent(fit(as.formula("y ~ `weird name`")))
expect_silent(fit(y ~ a, subset = 1:40))
expect_silent(fit(y ~ a, weights = w, offset = o))
expect_silent(dbarts::bart(
  cbind(a, b),
  y,
  seed = 1L,
  n.samples = 3L,
  n.burn = 3L,
  n.trees = 3L,
  n.chains = 1L,
  n.threads = 1L,
  keepTrees = FALSE,
  verbose = FALSE
))

## --- Block B: a forest() is a top-level term --------------------------------
# a forest() crossed with another term, either way round and by ':' or '*', is
# refused with the forest() to write in its place; what it says to write is
# the model the crossing stood for, a column for each member of a sum
indicators <- function(f) outer(as.integer(f), seq_len(nlevels(f)), "==") * 1
crossed <- list(
  list(y ~ x1 + x2 + z:forest(x1 + x2), "forest(x1 + x2, basis = ~ z)", z),
  list(y ~ x1 + x2 + forest(x1 + x2):z, "forest(x1 + x2, basis = ~ z)", z),
  list(
    y ~ x1 + x2 + zf:forest(x1 + x2),
    "forest(x1 + x2, basis = ~ zf)",
    indicators(zf)
  ),
  list(
    y ~ x1 + x2 + factor(z):forest(x1 + x2),
    "forest(x1 + x2, basis = ~ factor(z))",
    indicators(factor(z))
  ),
  list(
    y ~ x1 + x2 + (a + b):forest(x1, sd = 2),
    "forest(x1, sd = 2, basis = ~ cbind(a, b))",
    cbind(a, b)
  ),
  list(
    y ~ x1 + x2 + (log(a) + b + x3):forest(x1),
    "forest(x1, basis = ~ cbind(log(a), b, x3))",
    cbind(log(a), b, x3)
  ),
  list(y ~ x1 + x2 + z * forest(x1 + x2), "forest(x1 + x2, basis = ~ z)", z),
  list(y ~ x1 + x2 + forest(x1) * z, "forest(x1, basis = ~ z)", z),
  list(
    y ~ x1 + x2 + scale(a):forest(x1),
    "forest(x1, basis = ~ scale(a))",
    scale(a)
  ),
  list(y ~ x1 + x2 + (z):forest(x1), "forest(x1, basis = ~ z)", z),
  list(y ~ x1 + x2 + z:forest(), "forest(basis = ~ z)", z)
)
for (case in crossed) {
  refuses(
    case[[1L]],
    paste0(
      "a forest() is not crossed with another term; a forest's multiplier ",
      "is its 'basis' argument: write ",
      case[[2L]]
    )
  )
  rewritten <- fit(stats::as.formula(paste("y ~ x1 + x2 +", case[[2L]])))
  expect_equal(
    rewritten$fit$data@bases[[2L]],
    matrix(as.double(case[[3L]]), n),
    check.attributes = FALSE,
    info = case[[2L]]
  )
}
# a sum with a member that is no column of numbers has no one basis to write:
# a factor() call, a factor column, a character column
for (formula in list(
  y ~ x1 + x2 + (a + factor(z)):forest(x1),
  y ~ x1 + x2 + (a + zf):forest(x1),
  y ~ x1 + x2 + (g + a):forest(x1)
)) {
  refuses(
    formula,
    paste0(
      "forest(x1)': a forest() is not crossed with another term; a forest's ",
      "multiplier is its 'basis' argument, as forest(x1 + x2, basis = ~ z)"
    )
  )
}
# and beside bart()'s own n.trees the crossing is still what is refused
refuses(y ~ x1 + x2 + (a + zf):forest(x1, n.trees = 4L), "is not crossed with")
# the crossing is quoted as written
refuses(y ~ x1 + x2 + z:forest(x1 + x2), "'z:forest(x1 + x2)': a forest() is")
# at any depth, ahead of where the crossing itself stands
refuses(
  y ~ x1 + x2 + I(z:forest(x1 + x2)),
  "'z:forest(x1 + x2)': a forest() is not crossed with another term"
)
refuses(
  y ~ (x1 + z:forest(x2))^2,
  "'z:forest(x2)': a forest() is not crossed with another term"
)
refuses(
  y ~ x1 - z:forest(x2),
  "'z:forest(x2)': a forest() is not crossed with another term"
)
# crossed more than once, or a forest() on both sides: no one forest() to write
for (formula in list(
  y ~ x1 + x2 + z:a:forest(x1),
  y ~ x1 + forest(x1):z:a,
  y ~ x1 + a:(b:forest(x1)),
  y ~ x1 + x2 + forest(x1):forest(x2),
  y ~ x1 + x2 + z:forest(x1):forest(x2)
)) {
  refuses(
    formula,
    paste0(
      "a forest() is not crossed with another term; a forest's multiplier ",
      "is its 'basis' argument, as forest(x1 + x2, basis = ~ z)"
    )
  )
}
# crossed, and stating a basis already
refuses(
  y ~ x1 + x2 + z:forest(x1, basis = ~z),
  "its multiplier is its 'basis' argument, which this one already states"
)
# reached only by evaluating an expression, in a removal, on the left
refuses(
  y ~ x1 + x2 + I(forest(x1)),
  "must appear as a top-level additive term, not inside 'I(forest(x1))'"
)
refuses(
  y ~ x1 + x2 + (forest(x1, basis = ~z)),
  "must appear as a top-level additive term, not inside"
)
refuses(
  y ~ x1 + x2 - forest(x1, basis = ~z),
  "top-level additive term, not inside 'x1 + x2 - forest(x1, basis = ~z)'"
)
refuses(forest(x1) ~ x2, "left-hand side")
# the constructor as it is written outside the arguments that resolve it is
# the same term
for (formula in list(
  y ~ dbartsForests$forest(x1 + x2) + forest(x1, basis = ~z),
  y ~ forest(x1 + x2) + dbarts::dbartsForests$forest(x1, basis = ~z),
  y ~ dbarts:::forest(x1 + x2) + dbarts:::dbartsForests$forest(x1, basis = ~z)
)) {
  expectSameForest(formula, y ~ x1 + x2 + forest(x1, basis = ~z))
}
refuses(
  y ~ x1 + x2 + z:dbartsForests$forest(x1),
  "write forest(x1, basis = ~ z)"
)
# forest() is not exported, so a dbarts::-qualified head is no term
expect_error(
  fit(y ~ x1 + x2 + dbarts::forest(x1 + x2, basis = ~z)),
  "not an exported object"
)

# an offset() and an intercept term are the fit's wherever they are written,
# and count as no forest: one model, six ways
placements <- list(
  y ~ offset(o) + forest(x1 + x2) + forest(x1, basis = ~z),
  y ~ forest(x1 + x2) + offset(o) + forest(x1, basis = ~z),
  y ~ forest(x1 + x2) + forest(x1, basis = ~z) + offset(o),
  y ~ 0 + offset(o) + forest(x1 + x2) + forest(x1, basis = ~z),
  y ~ offset(o) + forest(x1 + x2) + forest(x1, basis = ~z) - 1,
  y ~ x1 + x2 + offset(o) + forest(x1, basis = ~z) - 1
)
reference <- fit(y ~ x1 + x2 + offset(o) + forest(x1, basis = ~z))
expect_identical(reference$fit$data@offset, o)
expect_identical(splitsOn(reference), list(c("x1", "x2"), "x1"))
# the fit's stored terms carry the intercept term, wherever it is written
expect_identical(attr(storedTerms(reference), "intercept"), 1L)
expect_identical(
  vapply(
    placements,
    function(formula) attr(storedTerms(fit(formula)), "intercept"),
    0L
  ),
  c(1L, 1L, 1L, 0L, 0L, 0L)
)
expect_identical(
  storedFormula(fit(y ~ forest(x1 + x2) + forest(x1, basis = ~z) - 1)),
  "~x1 + x2 - 1"
)
# after a forest that names no term too
expect_identical(
  attr(
    storedTerms(fit(y ~ forest(basis = ~z) - 1 + forest(x1 + x2))),
    "intercept"
  ),
  0L
)
for (formula in placements) {
  placed <- fit(formula)
  expect_identical(predictors(placed), c("x1", "x2"), info = deparse(formula))
  expect_identical(placed$fit$data@offset, o, info = deparse(formula))
  expect_identical(
    placed$yhat.train,
    reference$yhat.train,
    info = deparse(formula)
  )
  expect_identical(placed$sigma, reference$sigma, info = deparse(formula))
}
# and the offset is in the fit
expect_false(identical(
  reference$yhat.train,
  fit(y ~ x1 + x2 + forest(x1, basis = ~z))$yhat.train
))

## --- Block C: a forest's predictors are terms of its own --------------------
# the fit's predictors are every forest's, and each forest splits on its own
separate <- dbarts::dbarts(
  y ~ x1 + x2 + forest(x3, basis = ~z),
  d,
  control = dbarts::dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 5L,
    n.samples = 1L,
    updateState = FALSE,
    verbose = FALSE,
    seed = 5L
  )
)
expect_identical(colnames(separate$data@x), c("x1", "x2", "x3"))
expect_identical(
  attr(separate$control, "bartcore.forests")$vars,
  list(1:2, 3L)
)
splitCounts <- matrix(0, 3L, 2L, dimnames = list(c("x1", "x2", "x3"), NULL))
for (sweep in seq_len(300L)) {
  separate$run(0L, 1L)
  splitCounts <- splitCounts + separate$getForestVariableCounts()[,, 1L]
}
expect_true(all(splitCounts[c("x1", "x2"), 1L] > 0))
expect_identical(splitCounts[["x3", 1L]], 0)
expect_identical(unname(splitCounts[c("x1", "x2"), 2L]), c(0, 0))
expect_true(splitCounts[["x3", 2L]] > 0)

# terms as a model formula writes them: a call, a factor
called <- fit(
  y ~ forest(log(x1) + factor(g) + x2) + forest(log(x1) + factor(g), basis = ~z)
)
expect_identical(predictors(called), c("log(x1)", "factor(g)", "x2"))
expect_identical(
  splitsOn(called),
  list(c("log(x1)", "factor(g)", "x2"), c("log(x1)", "factor(g)"))
)
expect_identical(called$fit$data@x[, "log(x1)"], log(x1))
# a term of several columns is all of them
several <- fit(y ~ x1 + forest(poly(x2, 2), basis = ~z))
expect_identical(predictors(several), c("x1", "poly(x2, 2).1", "poly(x2, 2).2"))
expect_identical(
  splitsOn(several),
  list("x1", c("poly(x2, 2).1", "poly(x2, 2).2"))
)
indicatorFit <- fit(
  y ~ x1 + x2 + zf + forest(zf, basis = ~x1),
  factors = "indicators"
)
expect_identical(predictors(indicatorFit)[3:5], c("zf.u", "zf.v", "zf.w"))
expect_identical(forestInfo(indicatorFit)$vars[[2L]], 3:5)

# '.' is every column of the data but the response, and '-' removes a term:
# among the plain terms
dotOutside <- expectSameForest(
  y ~ . - z + forest(x1, basis = ~z),
  y ~ x1 + x2 + x3 + forest(x1, basis = ~z),
  data = small
)
expect_identical(predictors(dotOutside), c("x1", "x2", "x3"))
# and inside a forest
dotInside <- expectSameForest(
  y ~ forest(. - z) + forest(. - z - x2, basis = ~z),
  y ~ x1 + x2 + x3 + forest(x1 + x3, basis = ~z),
  data = small
)
expect_identical(splitsOn(dotInside), list(c("x1", "x2", "x3"), c("x1", "x3")))
# a forest's own '.' brings its columns to the fit
dotBrings <- fit(y ~ x1 + forest(. - z, basis = ~z), data = small)
expect_identical(predictors(dotBrings), c("x1", "x2", "x3"))
expect_identical(splitsOn(dotBrings), list("x1", c("x1", "x2", "x3")))

# The forest with no basis has the terms of one formula: the right-hand side
# with that forest() replaced by its contents and every forest with a basis
# deleted. So the forest written out and the same forest left as plain terms
# are one model whatever stands beside them: a removal before, between or
# after the forests takes its term from that forest alone, a forest with a
# basis keeping its own.
spellings <- list(
  # between the forests, of a term the forest with a basis names
  list(
    y ~ forest(.) - x3 + forest(x3, basis = ~z),
    y ~ . - x3 + forest(x3, basis = ~z),
    list(c("x1", "x2", "z"), "x3"),
    small
  ),
  list(
    y ~ forest(x1 + x2) - x2 + forest(x2, basis = ~z),
    y ~ x1 + x2 - x2 + forest(x2, basis = ~z),
    list("x1", "x2"),
    d
  ),
  # after the forests
  list(
    y ~ forest(.) + forest(x1, basis = ~z) - z,
    y ~ . - z + forest(x1, basis = ~z),
    list(c("x1", "x2", "x3"), "x1"),
    small
  ),
  list(
    y ~ . + forest(x1, basis = ~z) - z,
    y ~ . - z + forest(x1, basis = ~z),
    list(c("x1", "x2", "x3"), "x1"),
    small
  ),
  # before them, where a model formula has nothing yet to remove it from
  list(
    y ~ -x2 + forest(x1 + x2) + forest(x1, basis = ~z),
    y ~ -x2 + x1 + x2 + forest(x1, basis = ~z),
    list(c("x1", "x2"), "x1"),
    d
  ),
  # between two forests with a basis, each of which keeps its own terms
  list(
    y ~ forest(x1 + x2 + x3) +
      forest(x1 + x3, basis = ~z) -
      x3 +
      forest(x3, basis = ~a),
    y ~ x1 +
      x2 +
      x3 +
      forest(x1 + x3, basis = ~z) -
      x3 +
      forest(x3, basis = ~a),
    list(c("x1", "x2"), c("x1", "x3"), "x3"),
    d
  ),
  # of a term nobody names, which is no predictor at all
  list(
    y ~ forest(x1 + x2) - b + forest(x1, basis = ~z),
    y ~ x1 + x2 - b + forest(x1, basis = ~z),
    list(c("x1", "x2"), "x1"),
    d
  ),
  # of a column that is only a basis
  list(
    y ~ forest(x1 + x2) + forest(x1, basis = ~z) - z,
    y ~ x1 + x2 + forest(x1, basis = ~z) - z,
    list(c("x1", "x2"), "x1"),
    d
  ),
  # of the one term a forest with a basis has: it is that forest's own
  list(
    y ~ forest(x1 + x2) + forest(x3, basis = ~z) - x3,
    y ~ x1 + x2 + forest(x3, basis = ~z) - x3,
    list(c("x1", "x2"), "x3"),
    d
  ),
  # with the forest's predictors given by name
  list(
    y ~ forest(c("x1", "x2", "x3")) - x3 + forest(x3, basis = ~z),
    y ~ x1 + x2 + x3 - x3 + forest(x3, basis = ~z),
    list(c("x1", "x2"), "x3"),
    d
  ),
  # and beside an intercept term and an offset
  list(
    y ~ forest(x1 + x2) - x2 + offset(o) + forest(x2, basis = ~z) - 1,
    y ~ x1 + x2 - x2 + offset(o) + forest(x2, basis = ~z) - 1,
    list("x1", "x2"),
    d
  )
)
for (spelling in spellings) {
  spelled <- expectSameForest(
    spelling[[1L]],
    spelling[[2L]],
    data = spelling[[4L]]
  )
  expect_identical(
    splitsOn(spelled),
    spelling[[3L]],
    info = deparse(spelling[[1L]])
  )
  expect_identical(
    attr(storedTerms(spelled), "intercept"),
    attr(storedTerms(fit(spelling[[2L]], data = spelling[[4L]])), "intercept"),
    info = deparse(spelling[[1L]])
  )
}
# what a removal names that no forest without a basis has changes nothing
expectSameForest(
  y ~ forest(x1 + x2) + forest(x3, basis = ~z) - x3 - z - b,
  y ~ x1 + x2 + forest(x3, basis = ~z)
)
# alone too
removedAlone <- fit(y ~ forest(x1 + x2) - x2)
expect_identical(predictors(removedAlone), "x1")
expect_identical(removedAlone$yhat.train, fit(y ~ x1)$yhat.train)
# a removal that leaves the forest with no basis nothing is refused by name,
# in either spelling
refuses(
  y ~ forest(x1 + x2) - x1 - x2 + forest(x3, basis = ~z),
  paste0(
    "'forest(x1 + x2)': what the formula removes beside it leaves the ",
    "forest no predictor to split on"
  )
)
refuses(
  y ~ forest(x1) - x1,
  "'forest(x1)': what the formula removes beside it leaves the forest no"
)
refuses(
  y ~ x1 + x2 - x1 - x2 + forest(x3, basis = ~z),
  paste0(
    "the formula removes every plain term it writes (x1 + x2), which ",
    "leaves the forest with no multiplier no predictor to split on"
  )
)
# plain terms that a removal takes away are still plain terms beside a forest
refuses(
  y ~ x1 + forest(x2) - x1 + forest(x3, basis = ~z),
  "the formula has plain terms (x1) and a forest() with no basis"
)

# plain terms are stored as they are written, so a fit whose forests add no
# term stores the terms of the same formula with no forest, and predict asks
# of new rows what it asks then: a column '.' brought in and '-' took out
withForest <- fit(
  y ~ . - x3 - z + forest(x1, basis = ~z),
  data = small,
  keepTrees = TRUE
)
withoutForest <- fit(y ~ . - x3 - z, data = small, keepTrees = TRUE)
expect_identical(storedTerms(withForest), storedTerms(withoutForest))
lacking <- small[1:5, c("x1", "x2", "z")]
expect_error(
  predict(withoutForest, lacking),
  "missing variable required by the model: 'x3'",
  fixed = TRUE
)
expect_error(
  predict(withForest, lacking),
  "missing variable required by the model: 'x3'",
  fixed = TRUE
)
expect_identical(dim(predict(withForest, small[1:5, ])), c(3L, 5L))
expect_identical(
  storedFormula(fit(y ~ 0 + x1 + x2 + forest(x1, basis = ~z))),
  "~0 + x1 + x2"
)
# the terms the forests add follow them
expect_identical(
  storedFormula(fit(y ~ . - x3 - z + forest(x3, basis = ~z), data = small)),
  "~(x1 + x2 + x3 + z) - x3 - z + x3"
)

# no first argument: every predictor of the fit
everyPredictor <- fit(y ~ x1 + x2 + forest(basis = ~z))
expect_identical(splitsOn(everyPredictor), list(c("x1", "x2"), c("x1", "x2")))
expect_null(forestInfo(everyPredictor)$vars[[2L]])
# naming every predictor is no restriction
expectSameForest(
  y ~ x1 + x2 + forest(basis = ~z),
  y ~ x1 + x2 + forest(x1 + x2, basis = ~z)
)
# 'vars' by name is the unnamed argument
expectSameForest(
  y ~ x1 + x2 + forest(x1, basis = ~z),
  y ~ x1 + x2 + forest(vars = x1, basis = ~z)
)
# the one unnamed argument may be written second
expectSameForest(
  y ~ x1 + x2 + forest(x1, basis = ~z),
  y ~ x1 + x2 + forest(basis = ~z, x1)
)
# a value in a term is names of the fit's predictors
named <- expectSameForest(
  y ~ x1 + x2 + x3 + forest(x1 + x3, basis = ~z),
  y ~ x1 + x2 + x3 + forest(c("x1", "x3"), basis = ~z)
)
expect_identical(splitsOn(named)[[2L]], c("x1", "x3"))
# and each that is a column of the data comes to the fit as the name written
# out does, so a forest's predictors can be named with no term beside it
broughtIn <- expectSameForest(
  y ~ forest(c("x1", "x2")) + forest(x1, basis = ~z),
  y ~ forest(x1 + x2) + forest(x1, basis = ~z)
)
expect_identical(splitsOn(broughtIn), list(c("x1", "x2"), "x1"))
broughtBySecond <- expectSameForest(
  y ~ x1 + forest(c("x2", "x3"), basis = ~z),
  y ~ x1 + forest(x2 + x3, basis = ~z)
)
expect_identical(predictors(broughtBySecond), c("x1", "x2", "x3"))
expect_identical(splitsOn(broughtBySecond), list("x1", c("x2", "x3")))
namedAlone <- fit(y ~ forest(c("x1", "x2")))
expect_identical(predictors(namedAlone), c("x1", "x2"))
expect_identical(namedAlone$yhat.train, fit(y ~ x1 + x2)$yhat.train)
expect_identical(
  predictors(fit(y ~ forest(c("weird name", "x1")) + forest(x1, basis = ~z))),
  c("weird name", "x1")
)
# a name that is no column of the data is a column of the design, one of a
# term's several among them
expectSameForest(
  y ~ x1 + x2 + zf + forest(zf, basis = ~x1),
  y ~ x1 + x2 + zf + forest(c("zf.u", "zf.v", "zf.w"), basis = ~x1),
  factors = "indicators"
)
refuses(
  y ~ x1 + x2 + forest(c("x1", "nosuch"), basis = ~z),
  paste0(
    "'nosuch' is not a predictor of this fit (x1, x2); here a forest()'s ",
    "first argument selects among them"
  )
)
refuses(
  y ~ x1 + x2 + forest(c("x1", "x1"), basis = ~z),
  "names 'x1' more than once"
)
# a column of the data hides a variable of the caller's with its name, a
# formula held in one among them
local({
  b <- ~ x1 + x2
  hidden <- fit(y ~ x1 + forest(b, basis = ~z))
  expect_identical(predictors(hidden), c("x1", "b"))
  expect_identical(splitsOn(hidden), list("x1", "b"))
})

# and not positions
positionText <- paste0(
  "in a formula a forest's predictors are written as terms, as ",
  "forest(x1 + x2), or named, as forest(c(\"x1\", \"x2\")); a position is ",
  "for a 'forests' list"
)
refuses(y ~ x1 + x2 + forest(1, basis = ~z), positionText)
refuses(y ~ x1 + x2 + forest(1:2, basis = ~z), positionText)
refuses(y ~ x1 + x2 + forest(vars = 2L, basis = ~z), positionText)
refuses(y ~ x1 + x2 + forest(1, basis = ~z), "'forest(1, basis = ~z)': in a")
# no predictors anywhere
noPredictors <- paste0(
  "the formula names no predictors: write them as plain terms or inside a ",
  "forest(), as forest(x1 + x2)"
)
refuses(y ~ forest(basis = ~z), noPredictors)
refuses(y ~ 0 + forest(basis = ~z), noPredictors)
refuses(y ~ offset(o) + forest(basis = ~z), noPredictors)
refuses(y ~ forest(basis = ~a) + forest(basis = ~b), noPredictors)
# a tilde on the predictors: the multiplier of another order of arguments
refuses(
  y ~ x1 + x2 + forest(~z),
  paste0(
    "forest()'s first argument is the predictors the forest splits on, ",
    "written without '~', as forest(x1 + x2); a multiplier is 'basis ='"
  )
)
# a formula held in a variable
heldTerms <- ~ x1 + x2
refuses(
  y ~ x1 + x2 + forest(heldTerms, basis = ~z),
  paste0(
    "forest()'s first argument, 'heldTerms', holds a formula; write its ",
    "terms in place, as forest(x1 + x2), or give the predictors' names as a ",
    "character vector"
  )
)
# an offset() and an intercept term are the fit's
refuses(
  y ~ forest(x1 + offset(o)) + forest(x1, basis = ~z),
  paste0(
    "'forest(x1 + offset(o))': an offset() is the fit's and no forest's; ",
    "write it beside the forests, as y ~ offset(o) + forest(...)"
  )
)
for (formula in list(
  y ~ forest(0 + x1) + forest(x1, basis = ~z),
  y ~ forest(x1 - 1) + forest(x1, basis = ~z),
  y ~ forest(1 + x1) + forest(x1, basis = ~z),
  y ~ forest(x1) + forest(x1 + 0, basis = ~z)
)) {
  refuses(
    formula,
    paste0(
      "an intercept term (1, 0 or - 1) is the fit's and no forest's; write ",
      "it beside the forests"
    )
  )
}
# terms that leave a forest nothing
refuses(
  y ~ x1 + forest(. - x1 - x2 - x3 - z, basis = ~z),
  "its terms leave the forest no predictor to split on",
  data = small
)
# an interaction is refused as it is among plain terms
refuses(
  y ~ forest(x1:x2) + forest(x1, basis = ~z),
  "':' and '*' terms are not supported in 'formula'"
)
refuses(
  y ~ x1 + x2 + forest(x1 * x2, basis = ~z),
  "':' and '*' terms are not supported in 'formula'"
)
# one unnamed argument
oneUnnamed <- paste0(
  "forest() takes one unnamed argument, the predictors the forest splits ",
  "on, joined by '+' as forest(x1 + x2); every other argument is given by ",
  "name: a multiplier is 'basis =' and a size is 'sd ='"
)
refuses(y ~ x1 + x2 + forest(x1, x2), oneUnnamed)
refuses(y ~ x1 + x2 + forest(x1, a, 30), oneUnnamed)
refuses(y ~ x1 + x2 + forest(x1, ~z), oneUnnamed)
refuses(y ~ x1 + x2 + forest(x1, x2, basis = ~z), oneUnnamed)

## --- Block D: the forest with no basis --------------------------------------
# alone it is the single-forest fit its terms written plainly are
plain <- fit(y ~ x1 + x2)
writtenAlone <- fit(y ~ forest(x1 + x2))
expect_identical(predictors(writtenAlone), c("x1", "x2"))
expect_identical(writtenAlone$yhat.train, plain$yhat.train)
expect_identical(writtenAlone$sigma, plain$sigma)
expect_null(forestInfo(writtenAlone))
expect_null(attr(writtenAlone$fit$control, "bartcore.forestsDeclared"))
expect_identical(writtenAlone$fit$control@n.trees, plain$fit$control@n.trees)
expect_identical(writtenAlone$fit$data@x, plain$fit$data@x)
# with 'test', 'subset', weights and an offset, as any single-forest fit
expect_identical(
  fit(y ~ forest(x1 + x2), test = d[1:7, ])$yhat.test,
  fit(y ~ x1 + x2, test = d[1:7, ])$yhat.test
)
expect_identical(
  fit(y ~ offset(o) + forest(x1 + x2), subset = 1:40, weights = w)$yhat.train,
  fit(y ~ offset(o) + x1 + x2, subset = 1:40, weights = w)$yhat.train
)

# with others it is forest 1, as the plain terms are, wherever it is written
expectSameForest(
  y ~ forest(x1 + x2) + forest(x1, basis = ~z),
  y ~ x1 + x2 + forest(x1, basis = ~z)
)
swapped <- expectSameForest(
  y ~ forest(x1, basis = ~z) + forest(x1 + x2),
  y ~ x1 + x2 + forest(x1, basis = ~z)
)
expect_identical(splitsOn(swapped), list(c("x1", "x2"), "x1"))
expect_null(swapped$fit$data@bases[[1L]])
expect_identical(dim(swapped$fit$data@bases[[2L]]), c(n, 1L))
# its terms are the first of the fit's predictors
expectSameForest(
  y ~ forest(x2, basis = ~z) + forest(x1 + x2),
  y ~ x1 + x2 + forest(x2, basis = ~z)
)
expect_identical(
  predictors(fit(y ~ forest(x1, basis = ~z) + forest(x2 + x1))),
  c("x2", "x1")
)
# between two multiplied forests, which keep the order written
between <- expectSameForest(
  y ~ forest(x1, basis = ~a) +
    forest(x1 + x2) +
    forest(x2, basis = ~ factor(z)),
  y ~ x1 + x2 + forest(x1, basis = ~a) + forest(x2, basis = ~ factor(z))
)
expect_identical(splitsOn(between), list(c("x1", "x2"), "x1", "x2"))
expect_identical(
  lapply(between$fit$data@bases, dim),
  list(NULL, c(n, 1L), c(n, 2L))
)
# predict rebuilds each forest's basis by the forest's position
for (formula in list(
  y ~ forest(x1, basis = ~ scale(a)) + forest(x1 + x2),
  y ~ forest(x1, basis = ~ scale(a)) + forest(x2, basis = ~ factor(z))
)) {
  kept <- fit(formula, keepTrees = TRUE)
  expect_equal(
    predict(kept, d),
    kept$yhat.train,
    check.attributes = FALSE,
    info = deparse(formula)
  )
}
# under 'subset', weights and an offset
expectSameForest(
  y ~ offset(o) + forest(x1 + x2) + forest(x1, basis = ~ factor(z)),
  y ~ offset(o) + x1 + x2 + forest(x1, basis = ~ factor(z)),
  subset = 5:55,
  weights = w
)
# the families that take a term take it written this way
binary <- d
binary$yb <- as.integer(y > stats::median(y))
for (family in c("probit", "logistic")) {
  expectSameForest(
    yb ~ forest(x1 + x2) + forest(x1, basis = ~z),
    yb ~ x1 + x2 + forest(x1, basis = ~z),
    family = family,
    data = binary
  )
}
# and alone it is the single-forest fit in every family, those that take no
# multiplied forest among them
familyData <- d
familyData$yo <- cut(
  y,
  3L,
  labels = c("lo", "mid", "hi"),
  ordered_result = TRUE
)
familyData$yc <- stats::rpois(n, exp(x1))
familyData$ym <- factor(rep_len(c("p", "q", "r"), n))
familyData$time <- stats::rexp(n, exp(x1))
familyData$status <- rep_len(c(1, 1, 0), n)
families <- list(
  ordinal = "yo",
  nbinom = "yc",
  multinomial = "ym",
  student = "y"
)
if (requireNamespace("survival", quietly = TRUE)) {
  familyData$surv <- survival::Surv(familyData$time, familyData$status)
  families <- c(families, list(aft = "surv", hazard = "surv"))
}
# everything a fit reports but the call it was made by and its sampler
reported <- function(result) {
  result <- unclass(result)
  result[setdiff(names(result), c("call", "fit"))]
}
for (family in names(families)) {
  response <- families[[family]]
  plainly <- fit(
    stats::as.formula(paste(response, "~ x1 + x2")),
    family = family,
    data = familyData
  )
  expect_true(length(reported(plainly)) > 3L, info = family)
  # its predictors written as terms, and given by name
  for (written in c("forest(x1 + x2)", "forest(c(\"x1\", \"x2\"))")) {
    alone <- fit(
      stats::as.formula(paste(response, "~", written)),
      family = family,
      data = familyData
    )
    info <- paste(family, written)
    expect_identical(class(alone), class(plainly), info = info)
    expect_identical(reported(alone), reported(plainly), info = info)
  }
}
# and nothing but that is: a multinomial fit reads its formula before a
# forest could be declared, and must not read a multiplied forest, a second
# forest or a forest with an argument of its own as plain terms
for (formula in list(
  ym ~ forest(x1 + x2, basis = ~z),
  ym ~ forest(x1 + x2) + forest(x1, basis = ~z),
  ym ~ forest(x1, basis = ~z) + forest(x1 + x2),
  ym ~ x1 + x2 + forest(x1, basis = ~z),
  ym ~ forest(x1 + x2, n.trees = 4L),
  ym ~ forest(x1) + forest(x2)
)) {
  refuses(
    formula,
    "family = \"multinomial\" does not support a forest() formula term",
    family = "multinomial",
    data = familyData,
    n.trees = NULL
  )
}

# beside plain predictor terms it is refused, the plain terms named and never
# an offset
refuses(
  y ~ x1 + x2 + forest(x1 + x2),
  paste0(
    "the formula has plain terms (x1 + x2) and a forest() with no basis ",
    "('forest(x1 + x2)'): each is the forest with no multiplier, and a model ",
    "has one. Write the plain terms inside that forest(), or give it a 'basis'"
  )
)
refuses(
  y ~ x3 + offset(o) + forest(x1 + x2) + forest(x1, basis = ~z),
  "the formula has plain terms (x3) and a forest() with no basis"
)
# and twice
refuses(
  y ~ forest(x1) + forest(x2),
  paste0(
    "the formula has 2 forest() terms with no basis; a model has one forest ",
    "with no multiplier, and every other forest states a 'basis'"
  )
)
refuses(
  y ~ forest(x1) + forest(x2) + forest(x3) + forest(x1, basis = ~z),
  "the formula has 3 forest() terms with no basis"
)
# a forest with a basis stands beside another forest
refuses(
  y ~ forest(x1 + x2, basis = ~z),
  paste0(
    "a multi-forest model needs at least two forests, and this call's ",
    "'basis' declarations resolve to 1: a forest with a 'basis' stands ",
    "beside another forest. Write the forest with no multiplier too, as ",
    "y ~ forest(x1 + x2) + forest(x1 + x2, basis = ~ z1) or forests = ",
    "list(forest(), forest(basis = ~ z1)), or use a single forest with ",
    "linear() leaves; otherwise drop the basis"
  )
)
# what one forest cannot state
refuses(
  y ~ forest(x1 + x2, sd = 2),
  "this model has one forest, so its size is the fitting function's"
)
refuses(
  y ~ forest(x1 + x2, amplitude = fixed()),
  "this model has one forest, which has no coefficient to hold"
)

# the written first forest's constraints are the fit's: served there, and
# refused when the fitting function is given them too
constrained <- fit(
  y ~ forest(x1 + x2, interactions = interactions(max.order = 1L)) +
    forest(x1, basis = ~z)
)
constrainedAtTop <- fit(
  y ~ x1 + x2 + forest(x1, basis = ~z),
  interactions = dbarts::dbartsForests$interactions(max.order = 1L)
)
expect_false(is.null(forestInfo(constrained)$interactions[[1L]]))
expect_identical(forestInfo(constrained), forestInfo(constrainedAtTop))
expect_identical(constrained$yhat.train, constrainedAtTop$yhat.train)
expect_false(identical(
  constrained$yhat.train,
  fit(y ~ x1 + x2 + forest(x1, basis = ~z))$yhat.train
))
refuses(
  y ~ forest(x1 + x2, interactions = interactions(max.order = 1L)) +
    forest(x1, basis = ~z),
  "'interactions' is declared both at the top level and on the first forest",
  interactions = dbarts::dbartsForests$interactions(max.order = 2L)
)
constrainedAlone <- fit(
  y ~ forest(x1 + x2, interactions = interactions(max.order = 1L))
)
expect_identical(
  constrainedAlone$yhat.train,
  fit(
    y ~ x1 + x2,
    interactions = dbarts::dbartsForests$interactions(max.order = 1L)
  )$yhat.train
)
expect_false(identical(constrainedAlone$yhat.train, plain$yhat.train))

## --- Block E: the tree count given twice ------------------------------------
# bart()'s own 'n.trees' is the count of the forest with no basis
refuses(
  y ~ forest(x1 + x2, n.trees = 7L) + forest(x1, basis = ~z),
  paste0(
    "'n.trees' is given to the fitting function and to the forest with no ",
    "basis ('forest(x1 + x2, n.trees = 7L)'), which are the same count; ",
    "give one"
  )
)
refuses(
  y ~ forest(x1 + x2, n.trees = 7L),
  "'n.trees' is given to the fitting function and to the forest with no basis"
)
# stated on the forest alone, it is served
expect_identical(
  treeCounts(fit(
    y ~ forest(x1 + x2, n.trees = 7L) + forest(x1, basis = ~z),
    n.trees = NULL
  )),
  c(7, 50)
)
expect_identical(
  fit(y ~ forest(x1 + x2, n.trees = 7L), n.trees = NULL)$fit$control@n.trees,
  7L
)
expect_identical(
  fit(y ~ forest(x1 + x2, n.trees = 7L), n.trees = NULL)$yhat.train,
  fit(y ~ x1 + x2, n.trees = 7L)$yhat.train
)
# another forest's count beside bart()'s own is that forest's
expect_identical(
  treeCounts(fit(
    y ~ forest(x1 + x2) + forest(x1, basis = ~z, n.trees = 7L),
    n.trees = 15L
  )),
  c(15, 7)
)

## --- Block F: every forest multiplied ---------------------------------------
# the forests keep the order written, the first taking the fitting function's
# tree count and tree prior, as a 'forests' list of the two does
allMultiplied <- function(...) {
  dbarts::dbarts(
    ...,
    control = dbarts::dbartsControl(
      n.chains = 1L,
      n.threads = 1L,
      n.trees = 15L,
      n.samples = 20L,
      updateState = FALSE,
      verbose = FALSE,
      seed = 9L
    )
  )
}
written <- allMultiplied(y ~ forest(x1, basis = ~a) + forest(x2, basis = ~b), d)
listed <- allMultiplied(
  y ~ x1 + x2,
  d,
  forests = list(forest(x1, basis = ~a), forest(x2, basis = ~b))
)
reversed <- allMultiplied(
  y ~ forest(x2, basis = ~b) + forest(x1, basis = ~a),
  d
)
expect_identical(colnames(written$data@x), c("x1", "x2"))
expect_identical(written$control@n.trees, 15L)
writtenInfo <- attr(written$control, "bartcore.forests")
expect_identical(writtenInfo$params[[2L]][[1L]], 50)
expect_identical(writtenInfo$vars, list(1L, 2L))
expect_identical(written$data@bases, list(matrix(a, n, 1L), matrix(b, n, 1L)))
expect_identical(
  writtenInfo$vars,
  attr(listed$control, "bartcore.forests")$vars
)
expect_identical(written$data@bases, listed$data@bases)
writtenRun <- written$run(0L, 20L)
listedRun <- listed$run(0L, 20L)
expect_identical(writtenRun$train, listedRun$train)
expect_identical(writtenRun$sigma, listedRun$sigma)
expect_identical(written$getForestAmplitudes(), listed$getForestAmplitudes())
# the two terms swapped are another model
expect_identical(colnames(reversed$data@x), c("x2", "x1"))
expect_false(identical(reversed$run(0L, 20L)$train, writtenRun$train))

## --- Block G: what else is refused ------------------------------------------
# an all-zero basis column
refuses(y ~ x1 + x2 + forest(x1 + x2, basis = ~ rep(0, n)), "all zeros")
# 'test' with a model of several forests
refuses(
  y ~ forest(x1 + x2) + forest(x1 + x2, basis = ~z),
  "forest() formula term does not support 'test'",
  test = d
)
# an argument forest() does not have
refuses(
  y ~ x1 + x2 + forest(x1, basis = ~z, by = z),
  "unused argument (by = z)"
)
# a 'forests' list reaches no data object built beforehand
dd <- dbarts::dbartsData(y ~ x1 + x2, d)
expect_error(
  dbarts:::dbarts(dd, forests = list(forest(), forest(basis = ~z))),
  "already a dbartsData"
)
# the forests of a model are written one way in one call
expect_error(
  dbarts:::dbarts(
    y ~ x1 + x2 + forest(x1 + x2, basis = ~z),
    d,
    forests = list(forest(), forest(basis = ~z)),
    control = dbarts::dbartsControl(
      n.chains = 1L,
      n.trees = 3L,
      n.burn = 3L,
      n.samples = 3L,
      n.threads = 1L,
      verbose = FALSE
    )
  ),
  "only be declared one way"
)
# the family matrix - refused at family resolution, naming the family
termFormula <- y ~ x1 + x2 + forest(x1 + x2, basis = ~z)
notDefined <- function(family, term = "forest(x1 + x2, basis = ~z)") {
  paste0(
    "family \"",
    family,
    "\" does not support a forest() formula term ('",
    term,
    "'): an amplitude-coupled fit is not defined for it"
  )
}
for (family in c("hazard", "hazard.probit", "hazard.logistic")) {
  refuses(termFormula, notDefined(family), family = family)
}
refuses(termFormula, "family = \"multinomial\"", family = "multinomial")
refuses(
  termFormula,
  "family = \"hurdle.lognormal\"",
  family = "hurdle.lognormal"
)
# "twopart" is simply unrecognized, refused at family resolution (through
# match.arg, whose list names only what dbarts() takes) before the term is
# looked at
refuses(termFormula, "'family' should be one of", family = "twopart")
for (family in c("aft", "ordinal", "nbinom")) {
  refuses(termFormula, notDefined(family), family = family)
  # with the first forest written out it is the multiplied term that is named
  refuses(
    y ~ forest(x1 + x2) + forest(x1, basis = ~z),
    notDefined(family, "forest(x1, basis = ~z)"),
    family = family
  )
}
