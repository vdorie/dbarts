# An unnamed matrix given as predictors is read in the sampler's column order
# with no message, at every door: the test argument at fit time, setData,
# setPredictor, setTestPredictor, the sampler's predict and predict on a fit.
# A matrix of another width is refused, columns with names are matched by
# name, and a named test matrix given to an unnamed fit keeps its warning.
# Every sampler is seeded before creation: its chain seeds are drawn then.

n.um <- 50L
set.seed(1L)
x.um <- matrix(
  runif(n.um * 3L),
  n.um,
  3L,
  dimnames = list(NULL, c("a", "b", "c"))
)
y.um <- drop(x.um %*% c(1, 2, 3)) + rnorm(n.um, 0, 0.1)
xt.um <- x.um[1:7, , drop = FALSE] + 0.1
xtUnnamed.um <- unname(xt.um)
control.um <- dbarts::dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.samples = 5L,
  n.burn = 0L,
  keepTrees = TRUE,
  verbose = FALSE
)
samplerOf <- function(x = x.um, test = NULL, ...) {
  set.seed(4L)
  sampler <- dbarts::dbarts(
    x,
    y.um,
    test,
    control = control.um,
    n.trees = 10L,
    ...
  )
  set.seed(5L)
  invisible(sampler$run(5L, 5L))
  sampler
}

# test at fit time
set.seed(6L)
named <- dbarts::dbartsData(x.um, y.um, test = xt.um)
set.seed(6L)
bare <- NULL
expect_silent(bare <- dbarts::dbartsData(x.um, y.um, test = xtUnnamed.um))
expect_identical(bare@x.test, named@x.test)
expect_identical(colnames(bare@x.test), c("a", "b", "c"))

# setData and setPredictor
sampler <- samplerOf()
twin <- samplerOf()
expect_silent(sampler$setData(dbarts::dbartsData(unname(x.um) + 1, y.um)))
twin$setData(dbarts::dbartsData(x.um + 1, y.um))
expect_identical(unname(sampler$data@x), unname(twin$data@x))
expect_silent(sampler$setPredictor(unname(x.um)))
twin$setPredictor(x.um)
expect_identical(unname(sampler$data@x), unname(twin$data@x))
set.seed(8L)
drawsSampler <- sampler$run(2L, 3L)
set.seed(8L)
expect_identical(drawsSampler, twin$run(2L, 3L))

# setTestPredictor and the sampler's predict
sampler <- samplerOf()
twin <- samplerOf()
expect_silent(sampler$setTestPredictor(xtUnnamed.um))
twin$setTestPredictor(xt.um)
expect_identical(sampler$data@x.test, twin$data@x.test)
expect_silent(drawsBare <- sampler$predict(xtUnnamed.um))
expect_identical(drawsBare, twin$predict(xt.um))

# a matrix of another width is refused, a named one is matched by name
expect_error(
  sampler$predict(xtUnnamed.um[, 1:2]),
  "number of columns in 'test' must be equal to that of 'x'"
)
expect_error(
  sampler$setTestPredictor(xtUnnamed.um[, 1:2]),
  "number of columns in 'test' must be equal to that of 'x'"
)
expect_silent(drawsSwapped <- sampler$predict(xt.um[, c("c", "a", "b")]))
expect_identical(drawsSwapped, drawsBare)

# a formula fit on an unnamed matrix records X.1, X.2, X.3: the same matrix
# handed back is read in the fit's order, with no message
X.um <- unname(x.um)
set.seed(4L)
formulaSampler <- dbarts::dbarts(
  y.um ~ X.um,
  control = control.um,
  n.trees = 10L
)
set.seed(5L)
invisible(formulaSampler$run(5L, 5L))
XT.um <- unname(xt.um)
Xnamed.um <- XT.um
colnames(Xnamed.um) <- colnames(formulaSampler$data@x)
expect_false(is.null(colnames(formulaSampler$data@x)))
expect_silent(formulaSampler$setTestPredictor(XT.um))
expect_identical(formulaSampler$data@x.test, {
  formulaSampler$setTestPredictor(Xnamed.um)
  formulaSampler$data@x.test
})
expect_silent(formulaDraws <- formulaSampler$predict(XT.um))
expect_identical(formulaDraws, formulaSampler$predict(Xnamed.um))

# predict on a fit
set.seed(7L)
fit <- dbarts::bart(
  x.um,
  y.um,
  keepTrees = TRUE,
  n.trees = 10L,
  n.burn = 5L,
  n.samples = 5L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)
expect_silent(predBare <- predict(fit, xtUnnamed.um))
expect_identical(predBare, predict(fit, xt.um))
expect_error(
  predict(fit, xtUnnamed.um[, 1:2]),
  "number of columns in 'test' must be equal to that of 'x'"
)
set.seed(7L)
formulaFit <- dbarts::bart(
  y.um ~ X.um,
  keepTrees = TRUE,
  n.trees = 10L,
  n.burn = 5L,
  n.samples = 5L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)
expect_silent(formulaPred <- predict(formulaFit, XT.um))
expect_identical(formulaPred, predict(formulaFit, Xnamed.um))

# the reverse keeps its warning: a named test matrix on an unnamed fit
expect_warning(
  dbarts::dbartsData(X.um, y.um, test = xt.um),
  "'x' and 'test' are not both named",
  strict = TRUE
)
