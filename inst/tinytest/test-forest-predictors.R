# forest()'s first argument, the predictors a forest splits on. It is the one
# argument given unnamed, so that a later formal is an addition anywhere. In a
# 'forests' list it selects among the fit's predictors: as the terms of a
# model formula's right-hand side when a name in it is a predictor, and
# otherwise as the names or positions it gave where forest() was called, at
# that moment. Block A: one unnamed argument. Block B: one selection, six
# ways. Block C: what a selection refuses. Block D: a single forest. Block E:
# a forest keeps the selection it was given at the call.

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
  attr(sampler$control, "bartcore.forests", exact = TRUE)$vars
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
# with the value the code has at the call beside it, when it has one
expect_null(captured$vars$evaluated)
expect_identical(captured$vars$unbound, c("x1", "noSuchName"))
heldAtCall <- c("x1", "x3")
taken <- forest(heldAtCall)
expect_identical(taken$vars$expr, quote(heldAtCall))
expect_true(taken$vars$evaluated)
expect_identical(taken$vars$value, c("x1", "x3"))
expect_true(forest(NULL[1L])$vars$evaluated)
expect_null(forest(NULL[1L])$vars$value)
# taken quietly: a warning does not escape, and code that stops is code that
# could not be evaluated, not an error of forest()'s
expect_silent(forest(as.integer("one")))
expect_silent(forest(log(-1) + heldAtCall))
stopped <- forest(stop("not now"))
expect_null(stopped$vars$evaluated)
expect_identical(stopped$vars$error, "not now")
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
  "forest()'s first argument, 'x1 + offset(dose)': an offset() is the fit's"
)
# an intercept term, as in a formula and with its text
for (forests in list(
  quote(list(forest(), forest(0 + x1, basis = ~z))),
  quote(list(forest(), forest(x1 - 1, basis = ~z))),
  quote(list(forest(), forest(x1 + 1, basis = ~z))),
  quote(list(forest(), forest(1 + x1, basis = ~z)))
)) {
  refusesSelection(
    forests,
    paste0(
      "'",
      deparse(forests[[3L]][-3L]),
      "': an intercept term (1, 0 or - 1) is the fit's and no forest's; ",
      "write it beside the forests"
    )
  )
}
# a removal names a term of the fit as the fit's formula writes it, and one
# that names none is refused, not passed by
removals <- list(
  list(y ~ log(x1) + x2 + x3, quote(. - log(x1)), 2:3),
  list(y ~ x1 + factor(z) + x2, quote(. - factor(z)), c(1L, 3L)),
  list(y ~ x1 + poly(x2, 2) + x3, quote(. - poly(x2, 2)), c(1L, 4L)),
  list(y ~ I(x1^2) + x2 + x3, quote(. - I(x1^2)), 2:3),
  list(y ~ log(x1) + x2 + x3, quote(. - (log(x1) + x3)), 2L),
  list(y ~ x1 + poly(x2, 2) + x3, quote(.), 1:4)
)
for (case in removals) {
  removing <- listed(
    bquote(list(forest(), forest(.(case[[2L]]), basis = ~z))),
    case[[1L]]
  )
  expect_identical(
    forestColumns(removing),
    list(NULL, case[[3L]]),
    info = deparse(case[[2L]])
  )
}
refusesSelection(
  quote(list(forest(), forest(. - nosuch, basis = ~z))),
  notPredictor("nosuch")
)
refusesSelection(
  quote(list(forest(), forest(. - dose, basis = ~z))),
  notPredictor("dose")
)
refusesSelection(
  quote(list(forest(), forest(x1 + x2 - log(x3), basis = ~z))),
  notPredictor("log(x3)")
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

## --- Block E: a forest keeps the selection it was given at the call --------
# The value of forest()'s first argument is taken once, where and when
# forest() is called. A forest built in a loop, by lapply() or by Map(), a
# variable that changes before the fit, and a forest saved and read back each
# fit the model written: by the columns every forest splits on, and draw for
# draw against the selection spelled out.
built <- function(forests) {
  dbarts(x, y, forests = forests, control = predictorControl())
}
selected <- function(sampler) {
  lapply(forestColumns(sampler), function(columns) {
    if (is.null(columns)) colnames(x) else colnames(x)[columns]
  })
}
expectAsWritten <- function(forests, explicit, columns, info) {
  sampler <- built(unname(forests))
  # nolint start: object_usage_linter. tinytest attaches expect_* at run time.
  expect_identical(selected(sampler), columns, info = info)
  expect_identical(draws(sampler), draws(built(explicit)), info = info)
  # nolint end
}
everyColumn <- colnames(x)
oneEach <- list(
  forest(),
  forest(vars = "x1", basis = ~z),
  forest(vars = "x3", basis = ~z)
)
oneEachColumns <- list(everyColumn, "x1", "x3")
twoBases <- list(
  forest(),
  forest(vars = "x1", basis = ~z),
  forest(vars = "x3", basis = ~dose)
)
both <- list(forest(), forest(vars = c("x1", "x3"), basis = ~z))
bothColumns <- list(everyColumn, c("x1", "x3"))

# (a) a for loop: over names, over positions, and indexing a vector
loopNames <- list(forest())
for (name in c("x1", "x3")) {
  loopNames[[length(loopNames) + 1L]] <- forest(name, basis = ~z)
}
expectAsWritten(loopNames, oneEach, oneEachColumns, "a loop over names")
loopPositions <- list(forest())
for (i in c(1L, 3L)) {
  loopPositions[[length(loopPositions) + 1L]] <- forest(i, basis = ~z)
}
expectAsWritten(loopPositions, oneEach, oneEachColumns, "a loop over positions")
loopIndexed <- list(forest())
wanted <- c("x1", "x3")
for (i in 1:2) {
  loopIndexed[[i + 1L]] <- forest(basis = ~z, vars = wanted[i])
}
expectAsWritten(
  loopIndexed,
  oneEach,
  oneEachColumns,
  "a loop indexing a vector"
)

# (b) lapply() and sapply() over the names
expectAsWritten(
  c(list(forest()), lapply(c("x1", "x3"), forest, basis = ~z)),
  oneEach,
  oneEachColumns,
  "lapply(names, forest)"
)
expectAsWritten(
  c(list(forest()), sapply(wanted, forest, basis = ~z, simplify = FALSE)),
  oneEach,
  oneEachColumns,
  "sapply(names, forest)"
)
# written in the call itself
expect_identical(
  selected(dbarts(
    x,
    y,
    forests = c(list(forest()), lapply(c("x1", "x3"), forest, basis = ~z)),
    control = predictorControl()
  )),
  oneEachColumns
)

# (d) Map() and mapply(), a basis each
expectAsWritten(
  c(list(forest()), Map(forest, wanted, basis = list(~z, ~dose))),
  twoBases,
  oneEachColumns,
  "Map(forest, names, basis = )"
)
expectAsWritten(
  c(
    list(forest()),
    mapply(forest, wanted, basis = list(~z, ~dose), SIMPLIFY = FALSE)
  ),
  twoBases,
  oneEachColumns,
  "mapply(forest, names, basis = )"
)

# (g) the variable changed, emptied, removed or replaced before the fit
chosen <- c("x1", "x3")
early <- list(forest(), forest(chosen, basis = ~z))
chosen <- "x2"
expectAsWritten(early, both, bothColumns, "the variable changed")
chosen <- NULL
expectAsWritten(early, both, bothColumns, "the variable set to NULL")
chosen <- ~x2
expectAsWritten(early, both, bothColumns, "the variable replaced by a formula")
rm(chosen)
expectAsWritten(early, both, bothColumns, "the variable removed")
# a variable that does not exist yet at the call is not looked up later
late <- list(forest(), forest(notYet, basis = ~z))
notYet <- c("x1", "x3")
expect_error(built(late), notPredictor("notYet"), fixed = TRUE)

# a forest saved and read back carries its selection with it
chosen <- c("x1", "x3")
saved <- unserialize(serialize(
  list(forest(), forest(chosen, basis = ~z)),
  NULL
))
chosen <- "x2"
expectAsWritten(saved, both, bothColumns, "saved, the variable changed")
rm(chosen)
expectAsWritten(saved, both, bothColumns, "saved, the variable gone")

# a frame that rewrites its variable on the way out, an argument not yet
# forced whose own variable then changes, a wrapper that reassigns afterwards
rewriting <- function(cols) {
  on.exit(cols <- "x2")
  list(forest(), forest(cols, basis = ~z))
}
expectAsWritten(rewriting(c("x1", "x3")), both, bothColumns, "on.exit")
lazily <- function(cols) list(forest(), forest(cols, basis = ~z))
k <- 1L
unforced <- lazily(list(c("x1", "x3"), "x2")[[k]])
k <- 2L
expectAsWritten(unforced, both, bothColumns, "an unforced argument")
reassigning <- function(cols) {
  out <- list(forest(), forest(cols, basis = ~z))
  cols <- "x2"
  out
}
expectAsWritten(reassigning(c("x1", "x3")), both, bothColumns, "reassigned")

# what was right already stays so. (c) a function called once per forest
expectAsWritten(
  c(list(forest()), lapply(wanted, function(name) forest(name, basis = ~z))),
  oneEach,
  oneEachColumns,
  "a function per forest"
)
# (e) wrappers: a formal handed on, and dots, with names, terms and a variable
handsOn <- function(cols, ...) forest(cols, ...)
passesDots <- function(...) forest(...)
local({
  mine <- c("x1", "x3")
  for (forests in list(
    list(forest(), handsOn(c("x1", "x3"), basis = ~z)),
    list(forest(), passesDots(c("x1", "x3"), basis = ~z)),
    list(forest(), passesDots(x1 + x3, basis = ~z)),
    list(forest(), passesDots(mine, basis = ~z))
  )) {
    expectAsWritten(forests, both, bothColumns, "a wrapper")
  }
})
expect_error(
  built(list(forest(), handsOn(x1 + x3, basis = ~z))),
  "forest()'s first argument, 'cols': object 'x1' not found",
  fixed = TRUE
)
# (f) do.call with a value, a quoted call and a name; a call built with the
# value in it
local({
  mine <- c("x1", "x3")
  for (forests in list(
    list(forest(), do.call(forest, list(mine, basis = ~z))),
    list(forest(), do.call("forest", list(mine, basis = ~z))),
    list(forest(), do.call(forest, list(quote(x1 + x3), basis = ~z))),
    list(forest(), do.call(forest, list(vars = quote(. - x2), basis = ~z))),
    list(forest(), do.call(forest, list(as.name("mine"), basis = ~z))),
    list(forest(), eval(call("forest", mine, basis = ~z))),
    list(forest(), eval(bquote(forest(.(mine), basis = ~z))))
  )) {
    expectAsWritten(forests, both, bothColumns, "do.call and built calls")
  }
})
# (h) built in a function that has returned, and in local()
returned <- (function(cols) list(forest(), forest(cols, basis = ~z)))(wanted)
expectAsWritten(returned, both, bothColumns, "a frame that returned")
expectAsWritten(
  local({
    mine <- c("x1", "x3")
    list(forest(), forest(mine, basis = ~z))
  }),
  both,
  bothColumns,
  "local()"
)
# (i) one list of terms on two designs: by name on each, refused where the
# names are no predictors
byTerms <- list(forest(), forest(x1 + x3, basis = ~z))
reordered <- frame[c("y", "x3", "x1", "x2", "z")]
expect_identical(
  forestColumns(dbarts(
    y ~ x3 + x1 + x2,
    reordered,
    forests = byTerms,
    control = predictorControl()
  )),
  list(NULL, 1:2)
)
renamed <- data.frame(y = y, a = x[, 1L], b = x[, 2L], c = x[, 3L], z = z)
expect_error(
  dbarts(
    y ~ a + b + c,
    renamed,
    forests = byTerms,
    control = predictorControl()
  ),
  "'x1' is not a predictor of this fit (a, b, c)",
  fixed = TRUE
)
# (k) values: positions, names, a variable, an expression on one
local({
  mine <- c("x9", "x3")
  for (forests in list(
    list(forest(), forest(c(1L, 3L), basis = ~z)),
    list(forest(), forest(c("x1", "x3"), basis = ~z)),
    list(forest(), forest(c("x1", mine[2L]), basis = ~z)),
    list(forest(), forest(colnames(x)[-2L], basis = ~z)),
    list(forest(), forest(grep("[13]$", colnames(x)), basis = ~z))
  )) {
    expectAsWritten(forests, both, bothColumns, "values")
  }
})
for (refused in list(
  list(quote(forest(-3, basis = ~z)), "'vars' column index out of range"),
  list(quote(forest(0, basis = ~z)), "'vars' column index out of range"),
  list(quote(forest(pi, basis = ~z)), "'vars' must be a whole number"),
  list(quote(forest(c, basis = ~z)), notPredictor("c")),
  list(quote(forest(x1:x3, basis = ~z)), notPredictor("x1:x3")),
  list(quote(forest(I(x1), basis = ~z)), notPredictor("I(x1)"))
)) {
  refusesSelection(bquote(list(forest(), .(refused[[1L]]))), refused[[2L]])
}
