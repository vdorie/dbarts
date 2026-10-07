# forest()'s 'basis' is read as code. A term of a formula is read when the
# model is built, against the data and then where the formula was written, as
# lm() reads one. In a call of forest() the basis is bound to the call: a name
# of a column of the data is that column, and anything else is what it was
# where forest() was called, at that moment. Block A: where a name is found,
# at the four doors. Block B: what is code and what is a value. Block C: what
# is no basis. Block D: a forest keeps the basis its call was given. Block E:
# a number written beside a column. Block F: taken once, quietly. Block G: a
# variable called fixed.

forest <- dbartsForests$forest

set.seed(53)
n <- 90L
frame <- data.frame(
  x1 = runif(n),
  x2 = runif(n),
  x3 = runif(n),
  dose = runif(n, 0.5, 2),
  age = rnorm(n, 40, 10),
  z = rbinom(n, 1L, 0.5)
)
frame$y <- with(
  frame,
  sin(pi * x1) + z * (1 + x3) + 0.3 * dose * x2 + rnorm(n, sd = 0.3)
)
newRows <- data.frame(
  x1 = runif(7L),
  x2 = runif(7L),
  x3 = runif(7L),
  dose = runif(7L, 0.5, 2),
  age = rnorm(7L, 45, 8),
  z = rbinom(7L, 1L, 0.5)
)
x <- as.matrix(frame[c("x1", "x2", "x3")])
y <- frame$y

captureControl <- function() {
  dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 6L,
    n.samples = 4L,
    updateState = FALSE,
    verbose = FALSE,
    seed = 53L
  )
}
# the four doors: a term of the formula of dbarts() and of bart(), and a
# 'forests' list beside a formula and beside a matrix
termFit <- function(basis, data = frame, ...) {
  formula <- eval(bquote(y ~ x1 + x2 + x3 + forest(x1 + x3, basis = .(basis))))
  environment(formula) <- parent.frame()
  dbarts(formula, data, control = captureControl(), ...)
}
bartFit <- function(basis, data = frame, ...) {
  formula <- eval(bquote(y ~ x1 + x2 + x3 + forest(x1 + x3, basis = .(basis))))
  environment(formula) <- parent.frame()
  bart(
    formula,
    data,
    ...,
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 6L,
    n.samples = 4L,
    n.burn = 0L,
    keepTrees = TRUE,
    verbose = FALSE,
    seed = 53L
  )
}
listFit <- function(basis, data = frame, ...) {
  call <- bquote(dbarts(
    y ~ x1 + x2 + x3,
    data,
    forests = list(forest(), forest(x1 + x3, basis = .(basis))),
    control = captureControl()
  ))
  eval(as.call(c(as.list(call), list(...))), list(data = data), parent.frame())
}
matrixFit <- function(basis) {
  eval(
    bquote(dbarts(
      x,
      y,
      forests = list(forest(), forest(x1 + x3, basis = .(basis))),
      control = captureControl()
    )),
    parent.frame()
  )
}
# a list of forests built ahead of the fit
builtFit <- function(forests, data = NULL) {
  if (is.null(data)) {
    dbarts(x, y, forests = forests, control = captureControl())
  } else {
    dbarts(
      y ~ x1 + x2 + x3,
      data,
      forests = forests,
      control = captureControl()
    )
  }
}
# the model a value handed over fits, which reads no code at all
valueFit <- function(value, ...) {
  dbarts(
    x,
    y,
    forests = list(
      forest(),
      do.call(forest, list(vars = c("x1", "x3"), basis = value))
    ),
    control = captureControl(),
    ...
  )
}
basisOf <- function(sampler, index = 2L) sampler$data@bases[[index]]
labelsOf <- function(sampler) {
  attr(sampler$control, "bartcore.forests", exact = TRUE)$labels
}
draws <- function(sampler) {
  run <- sampler$run(0L, 12L)
  list(run$train, run$sigma, sampler$getForestAmplitudes())
}
column <- function(values, name) {
  matrix(as.double(values), ncol = 1L, dimnames = list(NULL, name))
}

## --- Block A: where a name is found -----------------------------------------
# a column of the data, with no variable of that name where the call is made
dataDraws <- draws(valueFit(frame$dose))
for (fit in list(termFit(quote(dose)), listFit(quote(dose)))) {
  expect_identical(basisOf(fit), column(frame$dose, "dose"))
  expect_identical(draws(fit), dataDraws)
}
expect_identical(
  bartFit(quote(dose))$bases[[2L]],
  column(frame$dose, "dose")
)
# with a second `dose` where the call is made, the data's column is the one
# used, at every door that has data
local({
  dose <- frame$dose * 100
  for (fit in list(termFit(quote(dose)), listFit(quote(dose)))) {
    expect_identical(basisOf(fit), column(frame$dose, "dose"))
    expect_identical(draws(fit), dataDraws)
  }
  expect_identical(
    bartFit(quote(dose))$bases[[2L]],
    column(frame$dose, "dose")
  )
  # with no data frame there is no column to name, and the caller's is found
  # where the call is made
  expect_identical(
    basisOf(matrixFit(quote(dose))),
    column(frame$dose * 100, "dose")
  )
  expect_identical(
    draws(matrixFit(quote(dose))),
    draws(valueFit(frame$dose * 100))
  )
  # handed over, the caller's value is used whatever the data's columns are
  handed <- dbarts(
    y ~ x1 + x2 + x3,
    frame,
    forests = list(
      forest(),
      do.call(forest, list(vars = c("x1", "x3"), basis = dose))
    ),
    control = captureControl()
  )
  expect_identical(unname(basisOf(handed)), unname(column(dose, "dose")))
  expect_null(colnames(basisOf(handed)))
  # and a column by its name, from a program
  named <- dbarts(
    y ~ x1 + x2 + x3,
    frame,
    forests = c(
      list(forest()),
      lapply(c("dose", "age"), function(name) {
        do.call(forest, list(basis = as.name(name)))
      })
    ),
    control = captureControl()
  )
  expect_identical(basisOf(named, 2L), column(frame$dose, "dose"))
  expect_identical(basisOf(named, 3L), column(frame$age, "age"))
})
# a name found nowhere says where it was looked for
expect_error(
  listFit(quote(nosuch)),
  "'basis' (nosuch): object 'nosuch' not found in 'data' or where forest() was called",
  fixed = TRUE
)
expect_error(
  termFit(quote(nosuch)),
  "'basis' (nosuch): object 'nosuch' not found in 'data' or where forest() was called",
  fixed = TRUE
)
expect_error(
  matrixFit(quote(nosuch)),
  "'basis' (nosuch): object 'nosuch' not found where forest() was called",
  fixed = TRUE
)
expect_error(
  listFit(quote(scale(nosuch) + dose)),
  "'basis' (scale(nosuch) + dose): object 'nosuch' not found in 'data' or where forest() was called",
  fixed = TRUE
)

## --- Block B: what is code and what is a value -------------------------------
twoColumns <- cbind(dose = frame$dose, age = frame$age)
twoDraws <- draws(valueFit(twoColumns))
# a formula held in a variable stands for that formula, and a tilde written
# in place is the same code
held <- ~ dose + age
for (fit in list(
  termFit(quote(held)),
  listFit(quote(held)),
  termFit(quote(~ dose + age)),
  listFit(quote(~ dose + age)),
  listFit(quote(dose + age))
)) {
  expect_identical(basisOf(fit), twoColumns)
  expect_identical(draws(fit), twoDraws)
  expect_identical(labelsOf(fit), c("forest1", "dose + age"))
}
# a formula handed over is that formula too, as is a call or a name
for (handedOver in list(held, quote(dose + age))) {
  fit <- builtFit(
    list(forest(), do.call(forest, list(c("x1", "x3"), basis = handedOver))),
    frame
  )
  expect_identical(basisOf(fit), twoColumns)
  expect_identical(draws(fit), twoDraws)
}
# a column of the data hides a variable of the caller's with its name, a
# formula it holds included; handed over, the formula is used
withB <- frame
withB$b <- frame$age
b <- ~dose
for (fit in list(termFit(quote(b), withB), listFit(quote(b), withB))) {
  expect_identical(basisOf(fit), column(frame$age, "b"))
}
handedFormula <- builtFit(
  list(forest(), do.call(forest, list(c("x1", "x3"), basis = b))),
  withB
)
expect_identical(basisOf(handedFormula), column(frame$dose, "dose"))
# a caller's variable named like a constructor is the caller's
local({
  forest <- frame$age
  fit <- dbarts(
    y ~ x1 + x2 + x3,
    frame,
    forests = list(
      dbartsForests$forest(),
      dbartsForests$forest(basis = forest)
    ),
    control = captureControl()
  )
  expect_identical(basisOf(fit), column(frame$age, "forest"))
})
# a value handed over has no code: its columns keep the names it has when
# every column has one and they differ, and otherwise have none
expect_identical(basisOf(valueFit(twoColumns)), twoColumns)
partlyNamed <- cbind(1 - frame$z, z = frame$z)
expect_identical(colnames(partlyNamed), c("", "z"))
expect_null(colnames(basisOf(valueFit(partlyNamed))))
expect_identical(basisOf(valueFit(partlyNamed)), unname(partlyNamed))
expect_error(
  valueFit(cbind(a = frame$dose, a = frame$age)),
  "'basis' has two columns named \"a\"; the columns of a basis are told apart by name",
  fixed = TRUE
)
# forwarded through a wrapper's dots, the code is read where the wrapper's
# caller wrote it; through a formal, it is what the formal holds
through <- function(...) forest(...)
expect_identical(
  basisOf(builtFit(
    list(forest(), through(x1 + x3, basis = dose + age)),
    frame
  )),
  twoColumns
)
byFormal <- function(multiplier) forest(basis = multiplier)
expect_identical(
  basisOf(builtFit(list(forest(), byFormal(frame$dose)))),
  column(frame$dose, "multiplier")
)

## --- Block C: what is no basis ----------------------------------------------
# NULL, held or computed, states none: the list is then a list of two forests
# with nothing to tell them apart
noBasis <- "forest 2 needs a 'basis'"
nothing <- NULL
useZ <- FALSE
expect_error(listFit(quote(nothing)), noBasis, fixed = TRUE)
expect_error(
  listFit(quote(if (useZ) z else NULL)),
  noBasis,
  fixed = TRUE
)
expect_error(matrixFit(quote(nothing)), noBasis, fixed = TRUE)
# and in a formula the term is then the forest with no multiplier
expect_error(
  termFit(quote(nothing)),
  "each is the forest with no multiplier",
  fixed = TRUE
)
useZ <- TRUE
expect_identical(
  unname(basisOf(listFit(quote(if (useZ) z else NULL)))),
  unname(column(frame$z, "z"))
)
# the slips of a caller who builds forests by program
name <- "dose"
expect_error(
  listFit(quote(name)),
  paste0(
    "'basis' is 'name', which holds the string \"dose\": a basis is a ",
    "column, not its name. Write basis = dose, or from a program ",
    "do.call(forest, list(basis = as.name(\"dose\"))), or hand over the ",
    "column itself"
  ),
  fixed = TRUE
)
expect_error(
  listFit("dose"),
  "'basis' is the string \"dose\": a basis is a column, not its name",
  fixed = TRUE
)
expect_error(
  termFit(quote(name)),
  "'basis' is 'name', which holds the string \"dose\"",
  fixed = TRUE
)
expect_error(
  valueFit("dose"),
  "'basis' is the string \"dose\"",
  fixed = TRUE
)
expect_error(
  listFit(2),
  "'basis' is the single value 2: a basis has a value for every observation",
  fixed = TRUE
)
two <- 2
expect_error(
  listFit(quote(two)),
  "'basis' is 'two', which holds the single value 2: a basis has a value for every observation",
  fixed = TRUE
)
expect_error(
  matrixFit(quote(two)),
  "'basis' is 'two', which holds the single value 2",
  fixed = TRUE
)
twoSided <- y ~ dose
for (basis in list(quote(twoSided), quote(y ~ dose))) {
  expect_error(
    listFit(basis),
    "a 'basis' formula must be one-sided, as ~ dose",
    fixed = TRUE
  )
  expect_error(
    termFit(basis),
    "a 'basis' formula must be one-sided, as ~ dose",
    fixed = TRUE
  )
}
# a tilde where the predictors go is still refused there
expect_error(
  forest(~dose),
  "forest()'s first argument is the predictors the forest splits on",
  fixed = TRUE
)

## --- Block D: a forest keeps the basis its call was given -------------------
w1 <- frame$dose
w2 <- frame$age
spelledOut <- function() {
  builtFit(list(
    forest(),
    do.call(forest, list(basis = w1)),
    do.call(forest, list(basis = w2))
  ))
}
spelledDraws <- draws(spelledOut())
expectEachOwn <- function(forests, info, name = NULL) {
  fit <- builtFit(forests)
  expect_identical(
    unname(basisOf(fit, 2L)),
    unname(column(w1, "")),
    info = info
  )
  expect_identical(
    unname(basisOf(fit, 3L)),
    unname(column(w2, "")),
    info = info
  )
  if (!is.null(name)) {
    expect_identical(colnames(basisOf(fit, 2L)), name, info = info)
    expect_identical(colnames(basisOf(fit, 3L)), name, info = info)
  }
  expect_identical(draws(fit), spelledDraws, info = info)
}
# a `for` loop over held columns: each forest has its own, not the last
looped <- list(forest())
for (w in list(w1, w2)) {
  looped[[length(looped) + 1L]] <- forest(basis = w)
}
expectEachOwn(looped, "for loop", "w")
# lapply over a closure, and with forest itself as the function
expectEachOwn(
  c(list(forest()), lapply(list(w1, w2), function(w) forest(basis = w))),
  "lapply over a closure",
  "w"
)
expectEachOwn(
  c(list(forest()), lapply(list(w1, w2), forest, vars = NULL)),
  "lapply(columns, forest, vars = )"
)
expectEachOwn(
  c(
    list(forest()),
    unname(Map(forest, list(NULL, NULL), basis = list(w1, w2)))
  ),
  "Map(forest, ..., basis = list(...))"
)
# a formula in a list of them, handed on by Map(), stands for that formula
mapped <- builtFit(
  c(
    list(forest()),
    unname(Map(forest, list(NULL, NULL), basis = list(~dose, ~age)))
  ),
  frame
)
expect_identical(basisOf(mapped, 2L), column(frame$dose, "dose"))
expect_identical(basisOf(mapped, 3L), column(frame$age, "age"))
expect_identical(
  draws(mapped),
  draws(builtFit(
    list(
      forest(),
      do.call(forest, list(basis = w1)),
      do.call(forest, list(basis = w2))
    ),
    frame
  ))
)
# the variable changed, set to NULL and removed between forest() and the fit
oneDraws <- draws(builtFit(list(forest(), do.call(forest, list(basis = w1)))))
moved <- w1
changed <- list(forest(), forest(basis = moved))
moved <- w2
expect_identical(basisOf(builtFit(changed)), column(w1, "moved"))
expect_identical(draws(builtFit(changed)), oneDraws)
moved <- w1
emptied <- list(forest(), forest(basis = scale(moved)))
moved <- NULL
expect_identical(
  basisOf(builtFit(emptied)),
  column(as.vector(scale(w1)), "scale(moved)")
)
moved <- w1
removed <- list(forest(), forest(basis = moved))
rm(moved)
expect_identical(basisOf(builtFit(removed)), column(w1, "moved"))
expect_identical(draws(builtFit(removed)), oneDraws)
# a forest saved and read back fits what it was given, and carries the names
# its basis uses and no frame of the function that built it
savedForests <- local({
  beside <- numeric(1e5)
  kept <- w1
  unserialize(serialize(
    list(forest(), forest(x1 + x3, basis = kept)),
    NULL
  ))
})
expect_true(length(serialize(savedForests, NULL)) < 20000L)
expect_identical(basisOf(builtFit(savedForests)), column(w1, "kept"))
expect_identical(draws(builtFit(savedForests)), draws(valueFit(w1)))
expect_null(savedForests[[2L]]$vars$env)
# what nothing binds at the call is not looked up later than the call
early <- list(forest(), forest(basis = I(dose / notYet)))
notYet <- 10
expect_error(
  builtFit(early, frame),
  "'basis' (I(dose/notYet)): object 'notYet' not found in 'data' or where forest() was called",
  fixed = TRUE
)
earlyValue <- list(forest(), forest(basis = alsoNotYet))
alsoNotYet <- w1
expect_error(
  builtFit(earlyValue),
  "'basis' (alsoNotYet): object 'alsoNotYet' not found where forest() was called",
  fixed = TRUE
)
# at the top level of a session the caller's frame is the workspace, whose
# variables are copied at the call like any frame's: one changed afterwards
# changes nothing, and one first made afterwards is not found
workspaceNames <- c("captureShrink", "captureLate", "captureColumn")
local({
  on.exit(rm(
    list = intersect(workspaceNames, ls(globalenv())),
    envir = globalenv()
  ))
  assign("captureShrink", 10, envir = globalenv())
  assign("captureColumn", w1, envir = globalenv())
  atTopLevel <- eval(
    quote(list(
      dbartsForests$forest(),
      dbartsForests$forest(basis = I(dose / captureShrink)),
      dbartsForests$forest(basis = captureColumn)
    )),
    globalenv()
  )
  lateAtTopLevel <- eval(
    quote(list(
      dbartsForests$forest(),
      dbartsForests$forest(basis = I(dose / captureLate))
    )),
    globalenv()
  )
  assign("captureShrink", 1000, envir = globalenv())
  assign("captureColumn", w2, envir = globalenv())
  assign("captureLate", 5, envir = globalenv())
  fit <- builtFit(atTopLevel, frame)
  expect_identical(
    basisOf(fit, 2L),
    column(frame$dose / 10, "I(dose/captureShrink)")
  )
  expect_identical(basisOf(fit, 3L), column(w1, "captureColumn"))
  expect_error(
    builtFit(lateAtTopLevel, frame),
    "'basis' (I(dose/captureLate)): object 'captureLate' not found in 'data' or where forest() was called",
    fixed = TRUE
  )
})
# one forest object fits the data it is given: where its name is a column it
# is that column, and where it is not it is what the call was given
dose <- frame$dose * 100
either <- list(forest(), forest(basis = dose))
expect_identical(basisOf(builtFit(either, frame)), column(frame$dose, "dose"))
expect_identical(basisOf(builtFit(either)), column(frame$dose * 100, "dose"))
rm(dose)

## --- Block E: a number written beside a column ------------------------------
# in a call of forest() it is the number at the call: a loop over it gives
# each forest its own, whatever the variable holds afterwards
rescaled <- list(forest())
for (k in c(10, 30)) {
  rescaled[[length(rescaled) + 1L]] <- forest(basis = I(dose / k))
}
expectRescaled <- function(info) {
  fit <- builtFit(rescaled, frame)
  expect_identical(
    basisOf(fit, 2L),
    column(frame$dose / 10, "I(dose/k)"),
    info = info
  )
  expect_identical(
    basisOf(fit, 3L),
    column(frame$dose / 30, "I(dose/k)"),
    info = info
  )
  expect_identical(labelsOf(fit), c("forest1", "I(dose/k)", "I(dose/k).1"))
}
expectRescaled("after the loop")
k <- 1000
expectRescaled("the variable changed")
rm(k)
expectRescaled("the variable removed")
expect_identical(
  draws(builtFit(rescaled, frame)),
  draws(builtFit(
    list(forest(), forest(basis = I(dose / 10)), forest(basis = I(dose / 30))),
    frame
  ))
)
# and predict uses the number of the call, not one found later
bartRescaled <- function() {
  bart(
    y ~ forest(x1 + x2 + x3) + forest(x1, basis = I(dose / k)),
    frame,
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 6L,
    n.samples = 4L,
    n.burn = 0L,
    keepTrees = TRUE,
    verbose = FALSE,
    seed = 53L
  )
}
# a term of a formula is read as lm() reads one: the number is looked up when
# the model is fitted, and again by predict
k <- 10
termRescaled <- bartRescaled()
expect_identical(
  termRescaled$bases[[2L]],
  column(frame$dose / 10, "I(dose/k)")
)
atTen <- dbarts:::replayForestBasis(termRescaled$basis.terms[[2L]], newRows, 2L)
expect_identical(atTen, column(newRows$dose / 10, "I(dose/k)"))
k <- 30
atThirty <- dbarts:::replayForestBasis(
  termRescaled$basis.terms[[2L]],
  newRows,
  2L
)
expect_identical(atThirty, column(newRows$dose / 30, "I(dose/k)"))
# as is a formula, in place or held: in a loop every forest reads the last
lazily <- list(forest())
for (k in c(10, 30)) {
  lazily[[length(lazily) + 1L]] <- forest(basis = ~ I(dose / k))
}
lazyFit <- builtFit(lazily, frame)
expect_identical(basisOf(lazyFit, 2L), column(frame$dose / 30, "I(dose/k)"))
expect_identical(basisOf(lazyFit, 3L), column(frame$dose / 30, "I(dose/k)"))
# the list's own record rebuilds the basis at new rows with the call's number
k <- 10
listRescaled <- builtFit(list(forest(), forest(basis = I(dose / k))), frame)
k <- 1000
listRecord <- attr(listRescaled$control, "bartcore.forests")$basisTerms[[2L]]
expect_identical(
  dbarts:::replayForestBasis(listRecord, newRows, 2L),
  column(newRows$dose / 10, "I(dose/k)")
)
rm(k)

## --- Block F: taken once, quietly -------------------------------------------
evaluations <- 0L
counted <- function(value) {
  evaluations <<- evaluations + 1L
  value
}
once <- list(forest(), forest(basis = counted(w1)))
expect_identical(evaluations, 1L)
expect_identical(basisOf(builtFit(once)), column(w1, "counted(w1)"))
invisible(builtFit(once))
expect_identical(evaluations, 1L)
# what the evaluation warns of is kept, and raised by each fit that uses the
# value; a fit that reads the code against its columns raises none
warns <- function(value) {
  warning("raised by the basis")
  value
}
expect_silent(warned <- list(forest(), forest(basis = warns(w1))))
countWarnings <- function(expr) {
  count <- 0L
  withCallingHandlers(expr, warning = function(w) {
    if (identical(conditionMessage(w), "raised by the basis")) {
      count <<- count + 1L
    }
    invokeRestart("muffleWarning")
  })
  count
}
expect_identical(countWarnings(builtFit(warned)), 1L)
expect_identical(countWarnings(builtFit(warned)), 1L)
# an argument a wrapper has not yet evaluated is evaluated at the call, once
# and as quietly
throughFormal <- function(multiplier) forest(basis = multiplier)
evaluations <- 0L
expect_silent(
  formal <- list(forest(), throughFormal(counted(warns(w1))))
)
expect_identical(evaluations, 1L)
formalFit <- NULL
expect_identical(countWarnings(formalFit <- builtFit(formal)), 1L)
expect_identical(evaluations, 1L)
expect_identical(basisOf(formalFit), column(w1, "multiplier"))
# and one that cannot be evaluated says why when the forest is fitted
expect_silent(unevaluable <- list(forest(), throughFormal(nosuchColumn)))
expect_error(
  builtFit(unevaluable, frame),
  "'basis' (multiplier): object 'nosuchColumn' not found",
  fixed = TRUE
)
# code that stops is code with no value, which a fit then says
stops <- function(value) stop("not this basis")
expect_silent(stopped <- list(forest(), forest(basis = stops(w1))))
expect_error(
  builtFit(stopped),
  "'basis' (stops(w1)): not this basis",
  fixed = TRUE
)
# what is decided from the code alone is said at the call
expect_error(
  forest(basis = dose * age),
  "'basis' does not take '*' between its terms ('dose * age')",
  fixed = TRUE
)
# a basis that draws random numbers moves the generator once, at the call
set.seed(7L)
random <- list(forest(), forest(basis = stats::rnorm(n)))
afterCall <- stats::runif(1L)
set.seed(7L)
expected <- stats::rnorm(n)
expect_identical(afterCall, stats::runif(1L))
expect_identical(
  unname(basisOf(builtFit(random))),
  unname(column(expected, ""))
)

## --- Block G: a variable called fixed ---------------------------------------
# inside 'forests' the name fixed is a constructor; a basis is the caller's
# code and reads it as the caller's variable or the data's column
local({
  fixed <- frame$z
  fit <- dbarts(
    x,
    y,
    forests = list(forest(), forest(x1 + x3, basis = ~ factor(fixed))),
    control = captureControl()
  )
  expect_identical(
    unname(basisOf(fit)),
    unname(cbind(1 - frame$z, frame$z) + 0)
  )
  expect_identical(
    colnames(basisOf(fit)),
    c("factor(fixed)0", "factor(fixed)1")
  )
  expect_identical(
    draws(fit),
    draws(valueFit(factor(frame$z)))
  )
})
withFixed <- frame
withFixed$fixed <- frame$age
expect_identical(
  basisOf(listFit(quote(fixed), withFixed)),
  column(frame$age, "fixed")
)
expect_identical(
  basisOf(termFit(quote(fixed), withFixed)),
  column(frame$age, "fixed")
)
# an attached column of that name is found where the call is made
attach(withFixed["fixed"])
attached <- tryCatch(
  dbarts(
    x,
    y,
    forests = list(forest(), forest(basis = fixed)),
    control = captureControl()
  ),
  finally = detach()
)
expect_identical(basisOf(attached), column(frame$age, "fixed"))
# and the constructor is still the constructor where it is called
held <- dbarts(
  x,
  y,
  forests = list(
    forest(amplitude = fixed()),
    forest(basis = factor(frame$z), amplitude = fixed())
  ),
  control = captureControl()
)
expect_identical(
  vapply(
    attr(held$control, "bartcore.forests")$params,
    function(params) params[8L],
    0
  ),
  c(0, 0)
)
