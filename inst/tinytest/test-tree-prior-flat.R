# The tree prior is three plain arguments at every door: power, base and
# split.probs, the last a vector or dart(). tree.prior is read on dbarts() and
# dbartsSpec() with a warning, and refused on bart() and xbart() with its
# translation. A setting given to the fitting function and to the forest with
# no basis is refused by name; split probabilities are judged by their value.

set.seed(11L)
n <- 60L
x <- matrix(runif(3L * n), n, dimnames = list(NULL, c("x1", "x2", "x3")))
y <- x[, 1L] + x[, 2L] * x[, 3L] + rnorm(n, sd = 0.3)
frame <- data.frame(x, y = y)
z <- rbinom(n, 1L, 0.5)

control <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.samples = 4L,
  n.burn = 6L,
  n.trees = 5L,
  updateState = FALSE
)
onBart <- function(...) {
  dbarts::bart(
    ...,
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 5L,
    n.burn = 6L,
    verbose = FALSE,
    samplerOnly = TRUE
  )
}
recordOf <- function(sampler) sampler$model@tree.prior
# the record a constructed cgm() or dart() resolves to against a sampler's data
resolved <- function(prior, sampler, delay = NA_real_) {
  record <- dbarts:::resolveSplitProbabilities(prior, sampler$data)
  if (is(record, "dbartsDartPrior") && is.na(record@update.delay)) {
    record@update.delay <- delay
  }
  record
}
resetKey <- function(key) {
  env <- dbarts:::onceWarnState
  env[[key]] <- NULL
  invisible(NULL)
}
warningsOf <- function(expr) {
  seen <- character(0L)
  withCallingHandlers(
    expr,
    warning = function(w) {
      seen <<- c(seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    message = function(m) invokeRestart("muffleMessage")
  )
  seen
}
errorOf <- function(expr) {
  tryCatch(
    {
      suppressWarnings(suppressMessages(expr))
      "created"
    },
    error = conditionMessage
  )
}

## --- every door writes the record --------------------------------------------
splitProbs <- c(x1 = 2, x2 = 1, x3 = 1)
cgmRecord <- dbartsPriors$cgm(power = 3, base = 0.5, split.probs = splitProbs)

atDbarts <- dbarts::dbarts(
  y ~ .,
  frame,
  control = control,
  power = 3,
  base = 0.5,
  split.probs = splitProbs
)
expect_identical(recordOf(atDbarts), resolved(cgmRecord, atDbarts))
atSpec <- dbarts::dbartsSpec(
  dbarts::dbartsData(y ~ ., frame),
  control = control,
  power = 3,
  base = 0.5,
  split.probs = splitProbs
)
expect_identical(atSpec$model@tree.prior, recordOf(atDbarts))
atBart <- onBart(y ~ ., frame, power = 3, base = 0.5, split.probs = splitProbs)
expect_identical(recordOf(atBart), recordOf(atDbarts))
atBT <- dbarts::bartBT(
  x,
  y,
  power = 3,
  base = 0.5,
  splitprobs = splitProbs,
  ntree = 5L,
  verbose = FALSE,
  sampleronly = TRUE
)
expect_identical(recordOf(atBT), resolved(cgmRecord, atBT))
atMonotone <- onBart(
  y ~ .,
  frame,
  power = 3,
  base = 0.5,
  monotone = c(x1 = 1)
)
expect_identical(
  recordOf(atMonotone),
  resolved(dbartsPriors$cgm(power = 3, base = 0.5), atMonotone)
)
atRbart <- suppressWarnings(dbarts::rbart_vi(
  y ~ .,
  frame,
  group.by = rep(1:3, 20L),
  power = 3,
  base = 0.5,
  n.samples = 10L,
  n.burn = 1L,
  n.thin = 1L,
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 5L,
  verbose = FALSE,
  keepTrees = TRUE
))
expect_identical(
  recordOf(atRbart$fit[[1L]]),
  resolved(dbartsPriors$cgm(power = 3, base = 0.5), atRbart$fit[[1L]])
)
atPdbart <- suppressMessages(dbarts::pdbart(
  x,
  y,
  xind = 1L,
  pl = FALSE,
  power = 3,
  base = 0.5,
  n.trees = 5L,
  n.samples = 2L,
  n.burn = 1L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
))
expect_identical(atPdbart$fit$model@tree.prior@power, 3)
expect_identical(atPdbart$fit$model@tree.prior@base, 0.5)
atHurdle <- dbarts::bart(
  x,
  ifelse(y > 0.5, exp(y), 0),
  family = "hurdle.lognormal",
  power = 3,
  base = 0.5,
  keepTrees = TRUE,
  n.trees = 5L,
  n.samples = 2L,
  n.burn = 2L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)
for (part in c("zero", "positive")) {
  record <- atHurdle[[part]]$fit$model@tree.prior
  expect_identical(c(record@power, record@base), c(3, 0.5), info = part)
}

# DART: the default update.delay is half the burn-in
dartBart <- onBart(
  y ~ .,
  frame,
  power = 3,
  base = 0.5,
  split.probs = dart(a = 1)
)
dartRecord <- dbartsPriors$dart(a = 1)
dartRecord@power <- 3
dartRecord@base <- 0.5
expect_identical(recordOf(dartBart), resolved(dartRecord, dartBart, 3))
dartDbarts <- dbarts::dbarts(
  y ~ .,
  frame,
  control = control,
  power = 3,
  base = 0.5,
  split.probs = dart(a = 1)
)
expect_identical(recordOf(dartDbarts), recordOf(dartBart))

# xbart: the grid's power and base on every cell's record, DART's settings
# kept; read where each cell's model reaches the engine
cellModels <- list()
keepModel <- function(model) {
  cellModels[[length(cellModels) + 1L]] <<- model
}
for (name in c("bartcoreSamplerFromHandle", "bartcoreSetModel")) {
  suppressMessages(trace(
    name,
    tracer = bquote(.(keepModel)(model)),
    where = asNamespace("dbarts"),
    print = FALSE
  ))
}
xbartCells <- function(...) {
  cellModels <<- list()
  dbarts::xbart(
    x,
    y,
    method = "k-fold",
    n.reps = 1L,
    n.samples = 4L,
    n.burn = c(6L, 2L),
    n.trees = 5L,
    n.threads = 1L,
    seed = 3L,
    ...
  )
  lapply(cellModels, function(model) model@tree.prior)
}
cells <- xbartCells(power = c(1.5, 3), base = 0.5, split.probs = dart(a = 1))
expect_true(length(cells) >= 2L)
for (cell in cells) {
  expect_inherits(cell, "dbartsDartPrior")
  expect_identical(c(cell@a, cell@update.delay, cell@base), c(1, 3, 0.5))
}
expect_identical(
  sort(unique(vapply(cells, function(cell) cell@power, 0))),
  c(1.5, 3)
)
cells <- xbartCells(power = 3, base = c(0.5, 0.8), split.probs = splitProbs)
for (cell in cells) {
  expect_inherits(cell, "dbartsCGMPrior")
  expect_identical(
    cell@splitProbabilities,
    recordOf(atDbarts)@splitProbabilities
  )
}
expect_identical(
  sort(unique(vapply(cells, function(cell) cell@base, 0))),
  c(0.5, 0.8)
)
uniform <- xbartCells(split.probs = c(2, 1, 1) / num.vars)
expect_equal(uniform[[1L]]@splitProbabilities, c(0.5, 0.25, 0.25))
for (name in c("bartcoreSamplerFromHandle", "bartcoreSetModel")) {
  suppressMessages(untrace(name, where = asNamespace("dbarts")))
}

## --- split.probs --------------------------------------------------------------
recordFor <- function(...) recordOf(onBart(y ~ ., frame, ...))
reference <- recordFor(split.probs = c(0.5, 0.25, 0.25))
expect_identical(recordFor(split.probs = c(2, 1, 1) / num.vars), reference)
expect_identical(
  unname(recordFor(split.probs = c(x1 = 2, .default = 1))@splitProbabilities),
  reference@splitProbabilities
)
expect_identical(recordFor(split.probs = dart), recordFor(split.probs = dart()))
expect_identical(
  recordFor(split.probs = dbarts::dbartsPriors$dart(a = 1)),
  recordFor(split.probs = dart(a = 1))
)
expect_identical(dbartsPriors$dart(0.3)@a, 0.3)
expect_error(dbartsPriors$dart(power = 3), "unused argument")
expect_error(
  onBart(y ~ ., frame, split.probs = cgm()),
  "'split.probs' takes a vector of split probabilities or dart()",
  fixed = TRUE
)
# through a wrapper's dots, in its caller's frame
wrapper <- function(...) onBart(y ~ ., frame, ...)
expect_identical(
  recordOf(wrapper(split.probs = c(2, 1, 1) / num.vars)),
  reference
)
expect_identical(recordOf(wrapper(split.probs = dart(a = 1)))@a, 1)

## --- tree.prior on dbarts() and dbartsSpec() -----------------------------------
retired <- function(...) {
  warnings <- character(0L)
  sampler <- withCallingHandlers(
    dbarts::dbarts(y ~ ., frame, control = control, ...),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(sampler = sampler, warnings = warnings)
}
plain <- function(...) {
  recordOf(dbarts::dbarts(y ~ ., frame, control = control, ...))
}
leftText <- paste0(
  "'tree.prior' has left 'dbarts': write the tree prior as power, base and ",
  "split.probs (split.probs = dart() for DART), as bart() takes them; ",
  "'tree.prior' is removed in dbarts 1.1-0. The value was used."
)
resetKey("tombstone.tree.prior.dbarts")
first <- retired(tree.prior = cgm(power = 3))
expect_identical(first$warnings, leftText)
expect_identical(recordOf(first$sampler), plain(power = 3))
# once a session
expect_identical(retired(tree.prior = cgm(power = 3))$warnings, character(0L))
resetKey("tombstone.tree.prior.dbartsSpec")
specWarnings <- warningsOf({
  dbarts::dbartsSpec(dbartsData(x, y), control = control, tree.prior = cgm())
  dbarts::dbartsSpec(dbartsData(x, y), control = control, tree.prior = cgm())
})
expect_identical(
  specWarnings,
  sub("'dbarts'", "'dbartsSpec'", leftText, fixed = TRUE)
)
# NULL is silent and is no tree prior
resetKey("tombstone.tree.prior.dbarts")
expect_identical(retired(tree.prior = NULL)$warnings, character(0L))
expect_null(dbarts:::onceWarnState[["tombstone.tree.prior.dbarts"]])

# written calls in any spelling, read against the constructors' old formals
expect_identical(
  recordOf(retired(tree.prior = cgm(3, 0.5))$sampler),
  plain(power = 3, base = 0.5)
)
expect_identical(
  recordOf(retired(tree.prior = dbartsPriors$cgm(power = 3))$sampler),
  plain(power = 3)
)
expect_identical(
  recordOf(
    retired(
      tree.prior = dbarts:::cgm(power = 3, split.probs = splitProbs)
    )$sampler
  ),
  plain(power = 3, split.probs = splitProbs)
)
expect_identical(
  recordOf(retired(tree.prior = dart(power = 3, a = 1))$sampler),
  plain(power = 3, split.probs = dart(a = 1))
)
expect_identical(recordOf(retired(tree.prior = cgm)$sampler), plain())
expect_identical(
  recordOf(retired(tree.prior = dart)$sampler),
  plain(split.probs = dart())
)
# a written cgm(power = 3) states power alone: base given beside it is taken
expect_identical(
  errorOf(dbarts::dbarts(
    y ~ forest(x1 + x2, base = 0.5) + forest(x3, basis = factor(z)),
    frame,
    control = control,
    tree.prior = cgm(power = 3)
  )),
  "created"
)
# a value states power and base and is handed through whole
constructed <- dbartsPriors$cgm(power = 3, base = 0.95)
expect_identical(
  recordOf(retired(tree.prior = constructed)$sampler),
  plain(power = 3, base = 0.95)
)
expect_match(
  errorOf(dbarts::dbarts(
    y ~ forest(x1 + x2, base = 0.5) + forest(x3, basis = factor(z)),
    frame,
    control = control,
    tree.prior = constructed
  )),
  "'base' is given to the fitting function and to the forest with no basis",
  fixed = TRUE
)
# a record read off a sampler: non-uniform probabilities, and DART with its
# resolved update.delay kept under another burn-in
expect_identical(
  recordOf(retired(tree.prior = recordOf(atDbarts))$sampler),
  recordOf(atDbarts)
)
longBurn <- control
longBurn@n.burn <- 40L
dartRead <- withCallingHandlers(
  dbarts::dbarts(
    y ~ .,
    frame,
    control = longBurn,
    tree.prior = recordOf(dartDbarts)
  ),
  warning = function(w) invokeRestart("muffleWarning")
)
expect_identical(recordOf(dartRead), recordOf(dartDbarts))
# an object saved by 0.9-34 has no splitProbabilitiesSpec slot
saved <- recordOf(atDbarts)
attr(saved, "splitProbabilitiesSpec") <- NULL
expect_false(methods::.hasSlot(saved, "splitProbabilitiesSpec"))
expect_identical(
  recordOf(retired(tree.prior = saved)$sampler),
  recordOf(atDbarts)
)
# positional: the tenth argument, where 0.9-34 took tree.prior
positional <- withCallingHandlers(
  dbarts::dbarts(
    y ~ .,
    frame,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    FALSE,
    4L,
    cgm(power = 3),
    control = control
  ),
  warning = function(w) invokeRestart("muffleWarning")
)
expect_identical(recordOf(positional), plain(power = 3))
# beside a plain argument it is refused
expect_identical(
  errorOf(dbarts::dbarts(
    y ~ .,
    frame,
    control = control,
    power = 3,
    tree.prior = cgm()
  )),
  "'tree.prior' and 'power' both set the tree prior of 'dbarts'; give 'power' alone"
)

## --- tree.prior on bart() and xbart() ------------------------------------------
expect_identical(
  errorOf(onBart(y ~ ., frame, tree.prior = cgm(power = 3, base = 0.5))),
  paste0(
    "'tree.prior' is not an argument of 'bart': write the tree prior as ",
    "plain arguments, power = 3, base = 0.5"
  )
)
expect_identical(
  errorOf(dbarts::xbart(x, y, tree.prior = dart(power = 3, a = 1))),
  paste0(
    "'tree.prior' is not an argument of 'xbart': write the tree prior as ",
    "plain arguments, power = 3, split.probs = dart(a = 1); xbart() takes ",
    "power and base as grid axes"
  )
)

# a value (not a written call) with uniform split probabilities shows no
# split.probs, and states power and base both
constructed <- dbartsPriors$cgm(power = 3)
expect_identical(
  errorOf(onBart(y ~ ., frame, tree.prior = constructed)),
  paste0(
    "'tree.prior' is not an argument of 'bart': write the tree prior as ",
    "plain arguments, power = 3, base = 0.95"
  )
)

## --- the retired tree.prior reaches the record at dbartsSpec() -----------------
expect_identical(
  suppressWarnings(dbarts::dbartsSpec(
    dbarts::dbartsData(x, y),
    control = control,
    tree.prior = dbartsPriors$cgm(power = 3, split.probs = splitProbs)
  ))$model@tree.prior,
  dbarts::dbartsSpec(
    dbarts::dbartsData(x, y),
    control = control,
    power = 3,
    split.probs = splitProbs
  )$model@tree.prior
)

## --- a whole-number power or base is taken at the plain doors -------------------
expect_identical(
  recordOf(dbarts::dbarts(y ~ ., frame, control = control, power = 3L)),
  plain(power = 3)
)
expect_identical(
  recordOf(onBart(y ~ ., frame, power = 3L)),
  recordOf(onBart(y ~ ., frame, power = 3))
)

## --- a setting given at both levels --------------------------------------------
twice <- function(name, where) {
  paste0(
    "'",
    name,
    "' is given to the fitting function and to the forest with no basis ",
    where,
    ", which are the same setting; give one"
  )
}
twoForests <- y ~ forest(x1 + x2, power = 1) + forest(x3, basis = factor(z))
expect_identical(
  errorOf(onBart(twoForests, frame, power = 3)),
  twice("power", "('forest(x1 + x2, power = 1)')")
)
# the default's own value, named, is stated
expect_identical(
  errorOf(onBart(twoForests, frame, power = 2)),
  twice("power", "('forest(x1 + x2, power = 1)')")
)
# a qualified written cgm() states only what it names
for (written in list(
  quote(dbarts:::cgm(power = 3)),
  quote(dbartsPriors$cgm(power = 3))
)) {
  expect_identical(
    errorOf(eval(bquote(dbarts::dbarts(
      y ~ forest(x1 + x2, base = 0.5) + forest(x3, basis = factor(z)),
      frame,
      control = control,
      tree.prior = .(written)
    )))),
    "created",
    info = deparse(written)
  )
}
# names one by one: base beside the forest's power is taken
expect_identical(errorOf(onBart(twoForests, frame, base = 0.5)), "created")
# a forest with a basis stating it is taken
expect_identical(
  errorOf(onBart(
    y ~ forest(x1 + x2) + forest(x3, basis = factor(z), power = 1),
    frame,
    power = 3
  )),
  "created"
)
# a forests list and a data object's bases name the forest by position
expect_identical(
  errorOf(dbarts::dbarts(
    x,
    y,
    control = control,
    forests = list(forest(basis = factor(z)), forest(base = 0.5)),
    base = 0.8
  )),
  twice("base", "(forest 2)")
)
# a lone forest term is the model's forest
expect_identical(
  errorOf(onBart(y ~ forest(x1 + x2, power = 1), frame)),
  "created"
)
expect_identical(
  errorOf(onBart(y ~ forest(x1 + x2, power = 1), frame, power = 3)),
  twice("power", "('forest(x1 + x2, power = 1)')")
)
# a variance forest never takes the fitting function's value instead
expect_identical(
  errorOf(onBart(
    y ~ forest(x1 + x2, power = 1),
    frame,
    power = 3,
    variance = ~x3
  )),
  twice("power", "('forest(x1 + x2, power = 1)')")
)
# with no forest with no basis, each name with its own remedy
noPlain <- function(name, remedy) {
  paste0(
    "'",
    name,
    "' given to the fitting function is the tree prior of the forest with no ",
    "basis, and every forest of this model has a basis; state it on a forest, ",
    "as forest(x1, basis = a, ",
    remedy,
    ")"
  )
}
onAllBasis <- function(...) {
  dbarts::bart(
    y ~ forest(x1, basis = factor(z)) + forest(x2, basis = x3),
    frame,
    ...,
    n.chains = 1L,
    n.threads = 1L,
    n.burn = 6L,
    verbose = FALSE,
    samplerOnly = TRUE
  )
}
expect_identical(
  errorOf(onAllBasis(power = 3)),
  noPlain("power", "power = 3")
)
expect_identical(
  errorOf(onAllBasis(base = 0.5)),
  noPlain("base", "base = 0.5")
)

## --- split probabilities judged by their value --------------------------------
severalText <- paste0(
  "a model of several forests does not support 'split.probs' (a vector that ",
  "is not uniform, or dart()); drop it or fit a single-forest model"
)
bcf <- y ~ forest(x1 + x2 + x3) + forest(x3, basis = factor(z))
expect_identical(
  errorOf(onBart(bcf, frame, split.probs = c(1, 2, 1))),
  severalText
)
expect_identical(errorOf(onBart(bcf, frame, split.probs = dart())), severalText)
expect_identical(
  errorOf(onBart(bcf, frame, split.probs = c(1, 1, 1))),
  "created"
)
expect_identical(
  errorOf(onBart(bcf, frame, split.probs = 1 / num.vars)),
  "created"
)
expect_identical(
  errorOf(onAllBasis(split.probs = c(x1 = 1, x2 = 1))),
  "created"
)
multinomialText <- paste0(
  "a multinomial (softmax) model does not support 'split.probs' (a vector ",
  "that is not uniform, or dart()); drop it or fit a single-forest model"
)
category <- factor(sample(c("a", "b", "c"), n, TRUE))
for (probs in list(quote(c(1, 2, 1)), quote(dart()))) {
  expect_identical(
    errorOf(eval(bquote(onBart(
      x,
      category,
      family = "multinomial",
      split.probs = .(probs)
    )))),
    multinomialText
  )
}

## --- the variance forest ------------------------------------------------------
varianceOf <- function(sampler) {
  unlist(attr(sampler$control, "bartcore.variance")[c("power", "base")])
}
expect_identical(
  varianceOf(onBart(y ~ ., frame, power = 3, base = 0.5, variance = ~x3)),
  c(power = 3, base = 0.5)
)
expect_identical(
  varianceOf(onBart(
    y ~ forest(x1 + x2 + x3, power = 1, base = 0.5),
    frame,
    variance = ~x3
  )),
  c(power = 1, base = 0.5)
)
expect_identical(
  varianceOf(onBart(
    y ~ .,
    frame,
    power = 3,
    variance = varianceForest(~x3, power = 1.5)
  )),
  c(power = 1.5, base = 0.95)
)
dartVariance <- onBart(
  y ~ .,
  frame,
  split.probs = dart(),
  variance = ~x3
)
expect_identical(
  names(attr(dartVariance$control, "bartcore.variance")),
  c("n.trees", "base", "power", "columns")
)

## --- pdbart -------------------------------------------------------------------
resetKey("tombstone.pdbart.splitprobs")
resetKey("tombstone.pdbart.power")
resetKey("tombstone.pdbart.base")
pdbartWarnings <- warningsOf(dbarts::pdbart(
  x,
  y,
  xind = 1L,
  pl = FALSE,
  power = 3,
  base = 0.5,
  splitprobs = c(2, 1, 1),
  n.trees = 5L,
  n.samples = 2L,
  n.burn = 1L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
))
expect_false(any(grepl("'power'|'base'", pdbartWarnings)))
expect_true(any(grepl("'splitprobs'", pdbartWarnings)))

## --- continuation ---------------------------------------------------------------
dartControl <- control
dartControl@n.chains <- 2L
dartControl@updateState <- TRUE
set.seed(5L)
continued <- dbarts::dbarts(
  y ~ .,
  frame,
  control = dartControl,
  split.probs = dart()
)
invisible(continued$run(10L, 2L))
continued$storeState()
copied <- continued$copy()
expect_identical(copied$state, continued$state)
expect_identical(recordOf(copied), recordOf(continued))
serialized <- tempfile(fileext = ".rds")
saveRDS(continued, serialized)
loaded <- readRDS(serialized)
unlink(serialized)
loaded$storeState()
expect_identical(loaded$state, continued$state)
expect_identical(recordOf(loaded), recordOf(continued))
# a copy and a reload continue the original's streams, to the last digits
nextRun <- continued$run(0L, 3L)
expect_equal(copied$run(0L, 3L), nextRun, tolerance = 1e-10)
expect_equal(loaded$run(0L, 3L), nextRun, tolerance = 1e-10)
# a one-forest $setModel changing power is taken
oneForest <- dbarts::dbarts(y ~ ., frame, control = control)
model <- oneForest$model
model@tree.prior@power <- 1
oneForest$setModel(model)
expect_identical(oneForest$model@tree.prior@power, 1)

rm(list = ls())
