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
forestInfo <- function(result) attr(result$fit$control, "bartcore.forests")
predictors <- function(result) colnames(result$fit$data@x)
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
  expect_identical(fitA$fit$data@x, fitB$fit$data@x, info = info)
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
# refused with the forest() to write in its place
crossed <- list(
  list(y ~ x1 + x2 + z:forest(x1 + x2), "forest(x1 + x2, basis = ~ z)"),
  list(y ~ x1 + x2 + forest(x1 + x2):z, "forest(x1 + x2, basis = ~ z)"),
  list(y ~ x1 + x2 + zf:forest(x1 + x2), "forest(x1 + x2, basis = ~ zf)"),
  list(
    y ~ x1 + x2 + factor(z):forest(x1 + x2),
    "forest(x1 + x2, basis = ~ factor(z))"
  ),
  list(
    y ~ x1 + x2 + (a + b):forest(x1, sd = 2),
    "forest(x1, sd = 2, basis = ~ cbind(a, b))"
  ),
  list(y ~ x1 + x2 + z * forest(x1 + x2), "forest(x1 + x2, basis = ~ z)"),
  list(y ~ x1 + x2 + forest(x1) * z, "forest(x1, basis = ~ z)"),
  list(y ~ x1 + x2 + scale(a):forest(x1), "forest(x1, basis = ~ scale(a))"),
  list(y ~ x1 + x2 + z:forest(), "forest(basis = ~ z)")
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
}
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
# the term grammar names forest only, so a dbarts::-qualified head is not a
# term
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
refuses(
  y ~ x1 + x2 + forest(c("x1", "nosuch"), basis = ~z),
  "'vars' name not found in the design's column names: 'nosuch'"
)

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
for (family in c("hazard", "hazard.probit", "hazard.logistic")) {
  refuses(termFormula, paste0("family \"", family, "\""), family = family)
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
  refuses(termFormula, paste0("family \"", family, "\""), family = family)
}
