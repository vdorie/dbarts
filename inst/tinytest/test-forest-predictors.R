# forest()'s first argument, the predictors a forest splits on. It is the one
# argument given unnamed, so that a later formal is an addition anywhere. In a
# 'forests' list it selects among the fit's predictors: as the terms of a
# model formula's right-hand side when a name in it is a predictor, and
# otherwise as the names or positions it gave where forest() was called, at
# that moment. Block A: one unnamed argument. Block B: one selection, six
# ways. Block C: what a selection refuses. Block D: a single forest. Block E:
# a forest keeps the selection it was given at the call. Block F: a predictor
# is named by its label as the fit holds it, on a formula and on a matrix.

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
# taken quietly: a warning does not escape at the call, where it is kept for
# the fit that uses the value, and code that stops is code that could not be
# evaluated, not an error of forest()'s
expect_identical(
  vapply(forest(as.integer("one"))$vars$warnings, conditionMessage, ""),
  "NAs introduced by coercion"
)
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
# what the evaluation at the call warned of is raised by the fit that uses
# the value, once and as the caller's own; a fit that reads the code as terms
# does not use the value and raises nothing
warningsOf <- function(expr) {
  raised <- character(0L)
  withCallingHandlers(expr, warning = function(w) {
    raised <<- c(raised, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  raised
}
careful <- function(names) {
  warning("check these names")
  names
}
warned <- NULL
expect_identical(
  warningsOf(warned <- list(forest(), forest(careful(wanted), basis = ~z))),
  character(0L)
)
expect_identical(warningsOf(built(warned)), "check these names")
usesValue <- NULL
expect_identical(
  warningsOf(usesValue <- built(warned)),
  "check these names"
)
expect_identical(selected(usesValue), bothColumns)
local({
  # the caller's vectors of these names add with a warning; as terms they
  # are the predictors, and the sum is never used
  x1 <- 1:2
  x3 <- 1:3
  asTerms <- NULL
  expect_identical(
    warningsOf(asTerms <- list(forest(), forest(x1 + x3, basis = ~z))),
    character(0L)
  )
  expect_identical(length(asTerms[[2L]]$vars$warnings), 1L)
  expect_identical(warningsOf(built(asTerms)), character(0L))
})

# a selection is names or positions: a logical and a factor are refused, not
# read through their codes, as the value of code and as a value handed over
notSelection <- function(kind) {
  paste0(
    "forest()'s first argument is a ",
    kind,
    "; a selection is names or positions. A multiplier is given as 'basis ='"
  )
}
codes <- factor(c("x3", "x1"))
for (refused in list(
  list(quote(forest(TRUE, basis = ~z)), "logical"),
  list(quote(forest(c(TRUE, FALSE, TRUE), basis = ~z)), "logical"),
  list(quote(forest(z == 1L, basis = ~z)), "logical"),
  list(quote(do.call(forest, list(TRUE, basis = ~z))), "logical"),
  list(quote(forest(codes, basis = ~z)), "factor"),
  list(quote(forest(factor(c("x3", "x1")), basis = ~z)), "factor"),
  list(quote(do.call(forest, list(vars = codes, basis = ~z))), "factor")
)) {
  refusesSelection(
    bquote(list(forest(), .(refused[[1L]]))),
    notSelection(refused[[2L]])
  )
}
refusesSelection(quote(list(forest(TRUE))), notSelection("logical"))

## --- Block F: a predictor is named by its label as the fit holds it ---------
# The term label on a formula fit and the column name on a matrix fit,
# written as code or as a backticked name. The arithmetic of '+', '-' and '.'
# is on the columns, so a column whose name is itself a call or an operator
# is not mistaken for that call, and a term that names no predictor is
# refused wherever it stands.
onDesign <- function(design, forests) {
  eval(
    bquote(dbarts(
      .(design),
      y,
      forests = .(forests),
      control = predictorControl()
    )),
    parent.frame()
  )
}
calls <- cbind("log(x1)" = log(x[, 1L]), x2 = x[, 2L], "I(x3^2)" = x[, 3L]^2)
crossedName <- cbind(x1 = x[, 1L], "x1:x2" = x[, 1L] * x[, 2L], x3 = x[, 3L])
# in a formula a column's name is written backticked; a name with ':' in it
# is no predictor of a formula fit at all, so the operator there is '*'
namedFrame <- data.frame(y = y, x, crossed = crossedName[, 2L], z = z)
names(namedFrame)[names(namedFrame) == "crossed"] <- "x1 * x2"
labelled <- list(
  # a column named like a call, removed as code and as a backticked name
  list(calls, quote(. - log(x1)), 2:3),
  list(calls, quote(. - `log(x1)`), 2:3),
  list(calls, quote(. - I(x3^2)), 1:2),
  list(calls, quote(log(x1) + I(x3^2)), c(1L, 3L)),
  list(calls, quote(`log(x1)` + x2), 1:2),
  # a column named like an operator on two others
  list(crossedName, quote(. - x1:x2), c(1L, 3L)),
  list(crossedName, quote(. - `x1:x2`), c(1L, 3L)),
  list(crossedName, quote(x1:x2 + x3), 2:3),
  # the same three on a formula fit, where the label is the term's
  list(y ~ log(x1) + x2 + x3, quote(. - log(x1)), 2:3),
  list(y ~ log(x1) + x2 + x3, quote(. - `log(x1)`), 2:3),
  list(y ~ log(x1) + x2 + x3, quote(`log(x1)` + x3), c(1L, 3L)),
  list(y ~ x1 + `x1 * x2` + x3, quote(. - x1 * x2), c(1L, 3L)),
  list(y ~ x1 + `x1 * x2` + x3, quote(. - `x1 * x2`), c(1L, 3L)),
  list(y ~ x1 + `x1 * x2` + x3, quote(x1 * x2 + x3), 2:3),
  # terms in the order a model formula reads them, and a group
  list(y ~ x1 + x2 + x3, quote(-x1 + .), 1:3),
  list(y ~ x1 + x2 + x3, quote(. - x1 + x1), 1:3),
  list(y ~ x1 + x2 + x3, quote(. - (x1 + x2)), 3L),
  list(y ~ x1 + x2 + x3, quote(. - (. - x2)), 2L)
)
for (case in labelled) {
  sampler <- if (inherits(case[[1L]], "formula")) {
    eval(bquote(dbarts(
      .(case[[1L]]),
      namedFrame,
      forests = list(forest(), forest(.(case[[2L]]), basis = ~z)),
      control = predictorControl()
    )))
  } else {
    onDesign(
      case[[1L]],
      bquote(list(forest(), forest(.(case[[2L]]), basis = ~z)))
    )
  }
  expect_identical(
    forestColumns(sampler),
    list(NULL, case[[3L]]),
    info = paste(
      if (inherits(case[[1L]], "formula")) "formula" else "matrix",
      deparse(case[[2L]])
    )
  )
}
# what names no predictor is refused, on a matrix as on a formula
expect_error(
  onDesign(calls, quote(list(forest(), forest(. - log(x9), basis = ~z)))),
  "'log(x9)' is not a predictor of this fit (log(x1), x2, I(x3^2))",
  fixed = TRUE
)
expect_error(
  onDesign(calls, quote(list(forest(), forest(. - x1, basis = ~z)))),
  "'x1' is not a predictor of this fit (log(x1), x2, I(x3^2))",
  fixed = TRUE
)
expect_error(
  onDesign(crossedName, quote(list(forest(), forest(. - x1:x3, basis = ~z)))),
  "'x1:x3' is not a predictor of this fit (x1, x1:x2, x3)",
  fixed = TRUE
)

# a column with no name is selected by position and is one of '.'; it is
# never named
unnamedOne <- x
colnames(unnamedOne) <- c("x1", "", "x3")
unnamedCases <- list(
  list(quote(forest(basis = ~z, vars = c("x1", "x3"))), c(1L, 3L)),
  list(quote(forest(wanted, basis = ~z)), c(1L, 3L)),
  list(quote(forest(1:2, basis = ~z)), 1:2),
  list(quote(forest(x1 + x3, basis = ~z)), c(1L, 3L)),
  list(quote(forest(., basis = ~z)), 1:3),
  list(quote(forest(. - x1, basis = ~z)), 2:3)
)
for (case in unnamedCases) {
  expect_identical(
    forestColumns(onDesign(unnamedOne, bquote(list(forest(), .(case[[1L]]))))),
    list(NULL, case[[2L]]),
    info = deparse(case[[1L]])
  )
}
expect_identical(
  attr(
    onDesign(unnamedOne, quote(list(forest(vars = wanted))))$model,
    "forest.columns"
  ),
  c(1L, 3L)
)
expect_error(
  onDesign(unnamedOne, quote(list(forest(), forest(c("x1", ""), basis = ~z)))),
  "forest()'s first argument has an empty name; a predictor with no name is",
  fixed = TRUE
)
# a column whose name is NA is unnamed too, alone and beside an empty one
for (unnamed in list(c("x1", NA, "x3"), c(NA, "", "x3"), c("", NA, "x3"))) {
  design <- x
  colnames(design) <- unnamed
  shape <- paste(unnamed, collapse = ",")
  for (case in list(
    list(quote(forest(basis = ~z, vars = "x3")), 3L),
    list(quote(forest(x3, basis = ~z)), 3L),
    list(quote(forest(., basis = ~z)), 1:3),
    list(quote(forest(. - x3, basis = ~z)), 1:2),
    list(quote(forest(2:3, basis = ~z)), 2:3)
  )) {
    expect_identical(
      forestColumns(onDesign(design, bquote(list(forest(), .(case[[1L]]))))),
      list(NULL, case[[2L]]),
      info = paste(shape, deparse(case[[1L]]))
    )
  }
  expect_identical(
    attr(
      onDesign(design, quote(list(forest(vars = "x3"))))$model,
      "forest.columns"
    ),
    3L,
    info = shape
  )
  expect_error(
    onDesign(design, quote(list(forest(), forest(c("x3", NA), basis = ~z)))),
    "'vars' contains missing values",
    fixed = TRUE,
    info = shape
  )
  expect_error(
    onDesign(design, quote(list(forest(), forest(c("x3", ""), basis = ~z)))),
    "forest()'s first argument has an empty name",
    fixed = TRUE,
    info = shape
  )
}
# with no column names at all, positions and '.'
expect_identical(
  forestColumns(onDesign(
    unname(x),
    quote(list(forest(), forest(., basis = ~z)))
  )),
  list(NULL, 1:3)
)
expect_identical(
  forestColumns(onDesign(
    unname(x),
    quote(list(forest(), forest(c(1, 3), basis = ~z)))
  )),
  list(NULL, c(1L, 3L))
)
expect_error(
  onDesign(unname(x), quote(list(forest(), forest(. - nosuch, basis = ~z)))),
  "'nosuch' is not a predictor of this fit",
  fixed = TRUE
)

# two columns with one name: '.' is every column, positions select either,
# and the name, which cannot tell them apart, is refused
twice <- x
colnames(twice) <- c("x1", "x1", "x3")
twiceCases <- list(
  list(quote(forest(., basis = ~z)), 1:3),
  list(quote(forest(. - x3, basis = ~z)), 1:2),
  list(quote(forest(c(1, 2), basis = ~z)), 1:2),
  list(quote(forest("x3", basis = ~z)), 3L),
  list(quote(forest(x3, basis = ~z)), 3L)
)
for (case in twiceCases) {
  expect_identical(
    forestColumns(onDesign(twice, bquote(list(forest(), .(case[[1L]]))))),
    list(NULL, case[[2L]]),
    info = deparse(case[[1L]])
  )
}
for (forests in list(
  quote(list(forest(), forest(x1, basis = ~z))),
  quote(list(forest(), forest(. - x1, basis = ~z))),
  quote(list(forest(), forest(c("x1", "x3"), basis = ~z))),
  quote(list(forest(), forest(basis = ~z, vars = "x1")))
)) {
  expect_error(
    onDesign(twice, forests),
    "'x1' is the name of 2 predictors of this fit; a name selects one",
    fixed = TRUE,
    info = deparse(forests)
  )
}

# a label is also matched as R writes the same code out, so a column named
# with other spacing is reached by its code; the exact name wins where both
# are there
spaced <- cbind(
  "log( x1 )" = log(x[, 1L]),
  "x2*x3" = x[, 2L] * x[, 3L],
  x3 = x[, 3L]
)
for (case in list(
  list(quote(. - log(x1)), 2:3),
  list(quote(log(x1) + x3), c(1L, 3L)),
  list(quote(. - x2 * x3), c(1L, 3L)),
  list(quote(x2 * x3 + x3), 2:3),
  list(quote(. - `log( x1 )`), 2:3)
)) {
  expect_identical(
    forestColumns(onDesign(
      spaced,
      bquote(list(forest(), forest(.(case[[1L]]), basis = ~z)))
    )),
    list(NULL, case[[2L]]),
    info = deparse(case[[1L]])
  )
}
bothSpacings <- cbind(spaced, "log(x1)" = x[, 1L])
expect_identical(
  forestColumns(onDesign(
    bothSpacings,
    quote(list(forest(), forest(. - log(x1), basis = ~z)))
  )),
  list(NULL, 1:3)
)

# a term's columns are those the term produced, by the design's own layout,
# and never those whose names begin as the term's does: a factor g under
# indicators beside predictors named g.total and g2, and beside a numeric
# predictor g.u, named as an indicator of g is
clash <- data.frame(frame["y"], x1 = x[, 1L], g = rep_len(c("u", "v", "w"), n))
clash$g.total <- x[, 2L]
clash$g2 <- x[, 3L]
clash$g.u <- dose
clash$z <- z
onClash <- function(forests) {
  eval(
    bquote(dbarts(
      y ~ x1 + g + g.total + g2 + g.u,
      clash,
      factors = "indicators",
      forests = .(forests),
      control = predictorControl()
    )),
    parent.frame()
  )
}
expect_identical(
  colnames(onClash(quote(list(forest(), forest(g, basis = ~z))))$data@x),
  c("x1", "g.u", "g.v", "g.w", "g.total", "g2", "g.u")
)
for (case in list(
  list(quote(g), 2:4),
  list(quote(. - g), c(1L, 5:7)),
  list(quote(g + g2), c(2:4, 6L)),
  list(quote(g.total), 5L),
  list(quote(. - g.total), c(1:4, 6:7)),
  list(quote(g2 + g.total), 5:6),
  # the numeric predictor is the term of that name; the indicator beside it
  # is the factor's
  list(quote(g.u), 7L),
  list(quote(. - g.u), 1:6),
  list(quote(. - g - g.u), c(1L, 5:6)),
  list(quote(g.v + x1), c(1L, 3L)),
  list(quote(c(2, 7)), c(2L, 7L)),
  list(quote(c("g.total", "g2")), 5:6)
)) {
  expect_identical(
    forestColumns(onClash(bquote(list(
      forest(),
      forest(.(case[[1L]]), basis = ~z)
    )))),
    list(NULL, case[[2L]]),
    info = deparse(case[[1L]])
  )
}
# by value the shared name tells neither column from the other
expect_error(
  onClash(quote(list(forest(), forest(c("g.u", "x1"), basis = ~z)))),
  "'g.u' is the name of 2 predictors of this fit",
  fixed = TRUE
)
# the single forest's constraints read a term the same way
expect_identical(
  attr(onClash(quote(list(forest(g + x1))))$model, "forest.columns"),
  1:4
)
# on a matrix a column is its own term, whatever its neighbours are called
neighbours <- cbind(g = x[, 1L], g.total = x[, 2L], g2 = x[, 3L])
expect_identical(
  forestColumns(onDesign(
    neighbours,
    quote(list(forest(), forest(. - g, basis = ~z)))
  )),
  list(NULL, 2:3)
)
expect_identical(
  forestColumns(onDesign(
    neighbours,
    quote(list(forest(), forest(g, basis = ~z)))
  )),
  list(NULL, 1L)
)

# a number alone is a position; among terms it would be a column's name, and
# is refused unless that name is what is written, backticked
numbered <- x
colnames(numbered) <- c("2", "1", "x3")
numeral <- paste0(
  "'2' is not a predictor; give positions as the whole argument, as c(2, 3)"
)
for (forests in list(
  quote(list(forest(), forest(. - 2, basis = ~z))),
  quote(list(forest(), forest(x3 + 2, basis = ~z)))
)) {
  expect_error(onDesign(numbered, forests), numeral, fixed = TRUE)
}
refusesSelection(quote(list(forest(), forest(x1 + 2, basis = ~z))), numeral)
for (case in list(
  list(quote(2), 2L),
  list(quote(c(2, 3)), 2:3),
  list(quote("2"), 1L),
  list(quote(`2` + x3), c(1L, 3L)),
  list(quote(. - `2`), 2:3)
)) {
  expect_identical(
    forestColumns(onDesign(
      numbered,
      bquote(list(forest(), forest(.(case[[1L]]), basis = ~z)))
    )),
    list(NULL, case[[2L]]),
    info = deparse(case[[1L]])
  )
}
