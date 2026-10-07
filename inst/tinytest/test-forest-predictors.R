# forest()'s first argument, the predictors a forest splits on. It is the one
# argument given unnamed, so that a later formal is an addition anywhere. In a
# 'forests' list it selects among the fit's predictors: as the terms of a
# model formula's right-hand side when a name in it is a predictor, and
# otherwise as a value of names or positions, found where forest() was
# called. Block A: one unnamed argument. Block B: one selection, six ways.
# Block C: what a selection refuses. Block D: a single forest.

forest <- dbartsForests$forest

set.seed(41)
n <- 120L
frame <- data.frame(
  x1 = runif(n),
  x2 = runif(n),
  x3 = runif(n),
  dose = runif(n, 0.5, 2),
  z = rbinom(n, 1L, 0.5)
)
frame$y <- with(frame, sin(pi * x1) + z * (1 + x3) + rnorm(n, sd = 0.3))
x <- as.matrix(frame[c("x1", "x2", "x3")])
y <- frame$y
z <- frame$z
dose <- frame$dose
# a value for every row, of the positions 1 and 2: a multiplier under the
# retired order of forest()'s arguments
z3 <- rep_len(c(1, 2), n)

predictorControl <- function() {
  dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 8L,
    n.samples = 5L,
    updateState = FALSE,
    verbose = FALSE,
    seed = 41L
  )
}
# a sampler of the formula's three predictors with `forests` written in the
# call itself, as a caller writes it
listed <- function(forests, formula = y ~ x1 + x2 + x3) {
  eval(
    bquote(dbarts(
      .(formula),
      frame,
      forests = .(forests),
      control = predictorControl()
    )),
    parent.frame()
  )
}
forestColumns <- function(sampler) {
  attr(sampler$control, "bartcore.forests")$vars
}
draws <- function(sampler) {
  run <- sampler$run(0L, 20L)
  list(run$train, run$sigma, sampler$getForestAmplitudes())
}

## --- Block A: one unnamed argument ------------------------------------------
oneUnnamed <- paste0(
  "forest() takes one unnamed argument, the predictors the forest splits ",
  "on, joined by '+' as forest(x1 + x2); every other argument is given by ",
  "name: a multiplier is 'basis =' and a size is 'sd ='"
)
w <- dose
for (forests in list(
  quote(list(forest(), forest(x1, dose))),
  quote(list(forest(), forest(x1, ~z))),
  quote(list(forest(), forest(x1, dose, 30))),
  quote(list(forest(), forest(~ factor(z), c("x1")))),
  quote(list(forest(), do.call(forest, list("x1", w)))),
  quote(list(forest(), forest("x1", basis = ~z, 2)))
)) {
  expect_error(
    listed(forests),
    oneUnnamed,
    fixed = TRUE,
    info = deparse(forests)
  )
}
# outside any argument too, and through a wrapper's dots
expect_error(forest("x1", w), oneUnnamed, fixed = TRUE)
expect_error(do.call(forest, list("x1", w)), oneUnnamed, fixed = TRUE)
forwarding <- function(...) forest(...)
expect_error(forwarding("x1", w), oneUnnamed, fixed = TRUE)
expect_identical(forwarding("x1", sd = 2)$vars, "x1")
expect_identical(forwarding("x1", sd = 2)$sd, 2)
# the one unnamed argument is the predictors wherever it is written
expect_identical(
  forestColumns(listed(quote(list(forest(), forest(basis = ~z, x1 + x3))))),
  list(NULL, c(1L, 3L))
)
expect_identical(forest(sd = 2, "x2")$vars, "x2")

# the constructor outside any argument builds what it built before the
# predictors moved first
expect_identical(
  forest(n.trees = 5L),
  structure(
    list(
      basis = NULL,
      vars = NULL,
      n.trees = 5L,
      base = NULL,
      power = NULL,
      sd = NULL,
      interactions = NULL,
      blocks = NULL,
      amplitude.prior.variance = NULL,
      amplitude = NULL
    ),
    class = "dbartsForest"
  )
)
expect_identical(forest(), forest(vars = NULL))
expect_identical(
  names(formals(forest)),
  c(
    "vars",
    "basis",
    "sd",
    "n.trees",
    "base",
    "power",
    "amplitude",
    "interactions",
    "blocks",
    "amplitude.prior.variance"
  )
)
# a value handed over is kept; code is kept as code, unevaluated
expect_identical(forest("x1")$vars, "x1")
expect_identical(forest(vars = 2L)$vars, 2L)
expect_identical(
  do.call(forest, list(vars = c("x1", "x3")))$vars,
  c("x1", "x3")
)
captured <- forest(x1 + noSuchName)
expect_true(inherits(captured$vars, "dbartsForestTerms"))
expect_identical(captured$vars$expr, quote(x1 + noSuchName))
expect_identical(captured$vars$env, environment())
# inside 'forests' too: with the caller's frame, not the frame laid over it
# in which the constructors resolve by bare name
capturedInside <- dbarts:::evalInForestVocabulary(
  quote(forest(x1 + noSuchName)),
  dbarts:::forestConstructors["forest"],
  environment()
)
expect_identical(capturedInside$vars$env, environment())

## --- Block B: one selection, six ways ---------------------------------------
reference <- listed(quote(list(
  forest(),
  forest(vars = c("x1", "x3"), basis = ~z)
)))
expect_identical(forestColumns(reference), list(NULL, c(1L, 3L)))
referenceDraws <- draws(reference)
heldNames <- c("x1", "x3")
heldPositions <- c(1L, 3L)
selections <- list(
  terms = quote(list(forest(), forest(x1 + x3, basis = ~z))),
  removal = quote(list(forest(), forest(. - x2, basis = ~z))),
  held = quote(list(forest(), forest(heldNames, basis = ~z))),
  heldByName = quote(list(forest(), forest(vars = heldNames, basis = ~z))),
  positions = quote(list(forest(), forest(c(1, 3), basis = ~z))),
  heldPositions = quote(list(forest(), forest(heldPositions, basis = ~z))),
  # a call built with the names in it, as a program builds one
  built = as.call(list(
    as.name("list"),
    quote(forest()),
    as.call(list(as.name("forest"), c("x1", "x3"), basis = ~z))
  )),
  handedOver = quote(list(
    forest(),
    do.call(forest, list(heldNames, basis = ~z))
  ))
)
for (name in names(selections)) {
  sampler <- listed(selections[[name]])
  expect_identical(forestColumns(sampler), list(NULL, c(1L, 3L)), info = name)
  expect_identical(draws(sampler), referenceDraws, info = name)
}
# and it is a selection: the forest restricted otherwise draws otherwise
expect_false(identical(
  draws(listed(quote(list(forest(), forest(x2, basis = ~z))))),
  referenceDraws
))

# '.' is every predictor of the fit, not every column of the data: z and dose
# are columns of `frame` and no predictor
expect_identical(
  forestColumns(listed(quote(list(forest(), forest(., basis = ~z))))),
  list(NULL, 1:3)
)
expect_identical(
  forestColumns(listed(quote(list(forest(. - x1), forest(basis = ~z))))),
  list(2:3, NULL)
)
# a predictor of the fit hides a variable of the caller's with its name
local({
  x1 <- c("x2", "x3")
  expect_identical(
    forestColumns(listed(quote(list(forest(), forest(x1, basis = ~z))))),
    list(NULL, 1L)
  )
})
# a term of the fit's formula is selected as it is written there
expect_identical(
  forestColumns(listed(
    quote(list(forest(), forest(log(x1) + x3, basis = ~z))),
    y ~ log(x1) + x2 + x3
  )),
  list(NULL, c(1L, 3L))
)
# found where forest() was called, inside a function that builds the list
buildForests <- function(selected) {
  list(forest(), forest(selected, basis = ~z))
}
built <- dbarts(
  y ~ x1 + x2 + x3,
  frame,
  forests = buildForests(c("x2", "x3")),
  control = predictorControl()
)
expect_identical(forestColumns(built), list(NULL, 2:3))
# with no data frame: the matrix interface's column names
expect_identical(
  forestColumns(dbarts(
    x,
    y,
    forests = list(forest(), forest(x1 + x3, basis = ~ factor(z))),
    control = predictorControl()
  )),
  list(NULL, c(1L, 3L))
)
# and through dbartsSpec
expect_identical(
  attr(
    dbartsSpec(
      dbartsData(x, y),
      predictorControl(),
      forests = list(forest(), forest(. - x1, basis = ~ factor(z)))
    )$control,
    "bartcore.forests"
  )$vars,
  list(NULL, 2:3)
)

## --- Block C: what a selection refuses --------------------------------------
refusesSelection <- function(forests, pattern) {
  # nolint next: object_usage_linter. tinytest attaches expect_* at run time.
  expect_error(listed(forests), pattern, fixed = TRUE, info = deparse(forests))
}
notPredictor <- function(name) {
  paste0(
    "'",
    name,
    "' is not a predictor of this fit (x1, x2, x3); here a forest()'s first ",
    "argument selects among them"
  )
}
# an unknown name, alone and among predictors
refusesSelection(
  quote(list(forest(), forest(nosuch, basis = ~z))),
  notPredictor("nosuch")
)
refusesSelection(
  quote(list(forest(), forest(x1 + nosuch, basis = ~z))),
  notPredictor("nosuch")
)
# a column of the data that is no predictor of the fit
refusesSelection(
  quote(list(forest(), forest(x1 + dose, basis = ~z))),
  notPredictor("dose")
)
refusesSelection(
  quote(list(forest(), forest(log(x1), basis = ~z))),
  notPredictor("log(x1)")
)
refusesSelection(
  quote(list(forest(), forest(mean, basis = ~z))),
  notPredictor("mean")
)
refusesSelection(
  quote(list(forest(), forest(c("x1", "nosuch"), basis = ~z))),
  "'vars' name not found in the design's column names: 'nosuch'"
)
# a repeat: a value for every row would otherwise be read as a few columns,
# a forest with no multiplier fitted in silence
repeated <- paste0(
  "forest()'s first argument selects predictors of the fit and names 'x1' ",
  "more than once; name each once. A multiplier is given as 'basis ='"
)
refusesSelection(quote(list(forest(z3), forest(basis = ~dose))), repeated)
refusesSelection(
  quote(list(forest(), forest(c("x1", "x1"), basis = ~z))),
  repeated
)
refusesSelection(
  quote(list(forest(), forest(vars = c(1L, 3L, 1L), basis = ~z))),
  repeated
)
refusesSelection(
  quote(list(forest(), do.call(forest, list(z3, basis = ~z)))),
  repeated
)
# a tilde, and a formula held in a variable
refusesSelection(
  quote(list(forest(), forest(~ factor(z)))),
  paste0(
    "forest()'s first argument is the predictors the forest splits on, ",
    "written without '~', as forest(x1 + x2); a multiplier is 'basis ='"
  )
)
heldFormula <- ~ x1 + x3
refusesSelection(
  quote(list(forest(), forest(heldFormula, basis = ~z))),
  paste0(
    "forest()'s first argument, 'heldFormula', holds a formula; write its ",
    "terms in place, as forest(x1 + x2), or give the predictors' names as a ",
    "character vector"
  )
)
# terms that leave nothing, and what is the fit's and no forest's
refusesSelection(
  quote(list(forest(), forest(. - x1 - x2 - x3, basis = ~z))),
  "forest()'s first argument, '. - x1 - x2 - x3', leaves the forest no predictor"
)
refusesSelection(
  quote(list(forest(), forest(x1 + offset(dose), basis = ~z))),
  "an offset() and an intercept term are the fit's and no forest's"
)
refusesSelection(
  quote(list(forest(), forest(0 + x1, basis = ~z))),
  "an offset() and an intercept term are the fit's and no forest's"
)
# as before: out of range, empty, missing
refusesSelection(
  quote(list(forest(), forest(4L, basis = ~z))),
  "'vars' column index out of range"
)
refusesSelection(
  quote(list(forest(), forest(character(0L), basis = ~z))),
  "'vars' is empty"
)
refusesSelection(
  quote(list(forest(), forest(NA, basis = ~z))),
  "'vars' contains missing values"
)

## --- Block D: a single forest -----------------------------------------------
single <- listed(quote(list(forest(x1 + x3))))
expect_identical(attr(single$model, "forest.columns"), c(1L, 3L))
expect_identical(
  single$run(0L, 20L)$train,
  listed(quote(list(forest(vars = c("x1", "x3")))))$run(0L, 20L)$train
)
# naming every predictor restricts nothing and stores nothing
expect_null(attr(listed(quote(list(forest(.))))$model, "forest.columns"))
refusesSelection(quote(list(forest(z3))), repeated)
