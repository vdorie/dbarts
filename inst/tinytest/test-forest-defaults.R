# A forest's defaults go by its kind. The forest with no basis, the plain one,
# takes the fitting function's tree count, tree prior, 'interactions' and
# 'blocks' wherever it stands; a forest with a basis takes 50 trees and
# cgm(3, 0.25) wherever it stands. Block A: a plain forest anywhere. Block B:
# the defaults at every door. Block C: what counts as stated. Block D: with no
# plain forest, what the fitting function states for it is refused. Block E:
# 'interactions' and 'blocks' follow the plain forest. Block F: the count
# given twice at bart(). Block G: a control carried to another fit. Block H:
# print. Block I: which held coefficients are taken.

set.seed(71)
n <- 150L
frame <- data.frame(x1 = runif(n), x2 = runif(n), x3 = runif(n))
frame$dose <- runif(n, 5, 60)
frame$age <- rnorm(n, 50, 10)
frame$z <- rbinom(n, 1L, 0.5)
frame$y <- with(frame, x1 + x2 * x3 + 0.02 * dose * x3 + z * (1 + x1)) +
  rnorm(n, sd = 0.3)
x <- as.matrix(frame[c("x1", "x2", "x3")])
y <- frame$y
dose <- frame$dose
age <- frame$age
z <- frame$z
zBasis <- cbind(1 - z, z)
yBinary <- as.double(y > stats::median(y))

defaultsControl <- function(...) {
  settings <- list(n.chains = 1L, n.threads = 1L, n.samples = 5L, seed = 71L)
  do.call(dbartsControl, c(settings, list(updateState = FALSE, ...)))
}
# the fitting function's own and the multiplied forest's: count, base, power
fitDefault <- c(75, 0.95, 2)
fitStated <- c(15, 0.8, 3)
multiplied <- c(50, 0.25, 3)

forestInfo <- function(sampler) {
  attr(sampler$control, "bartcore.forests", exact = TRUE)
}
specSampler <- function(spec) {
  new("dbartsSampler", spec$control, spec$model, spec$data)
}
# the call of a fitting function written in 'text', given this file's
# settings, with 'stated' a tree count and a tree prior, and the arguments in
# '...', each code as text or a value; always a sampler
build <- function(text, stated = FALSE, ...) {
  call <- str2lang(text)
  door <- as.character(call[[1L]])
  if (door == "bart") {
    settings <- list(1L, 1L, 5L, 0L, FALSE, TRUE, 71L)
    names(settings) <- strsplit(
      "n.chains n.threads n.samples n.burn verbose samplerOnly seed",
      " "
    )[[1L]]
    call[names(settings)] <- settings
    if (stated) {
      call$n.trees <- 15L
    }
  } else {
    call$control <- call("defaultsControl")
    call$control$n.trees <- if (stated) 15L
  }
  if (stated) {
    call$tree.prior <- quote(cgm(3, 0.8))
  }
  extra <- lapply(list(...), function(argument) {
    if (is.character(argument)) str2lang(argument) else argument
  })
  call[names(extra)] <- extra
  built <- eval(call)
  if (door == "dbartsSpec") specSampler(built) else built
}
# each forest's tree count from the engine and its base and power from its
# record, one row a forest; the first forest's as the bridge reads them too
forestTable <- function(sampler) {
  info <- forestInfo(sampler)
  prior <- sampler$model@tree.prior
  onModel <- c(sampler$control@n.trees, prior@base, prior@power)
  if (is.null(info)) {
    return(rbind(onModel, deparse.level = 0L))
  }
  # nolint next: object_usage_linter. tinytest attaches expect_* at run time.
  expect_identical(info$params[[1L]][1:3], onModel)
  pointer <- sampler$getPointer()
  counts <- vapply(
    seq_along(info$params),
    function(index) dbarts:::bartcoreForestTreeCount(pointer, index - 1L),
    0L
  )
  cbind(
    counts,
    do.call(rbind, lapply(info$params, `[`, 2:3)),
    deparse.level = 0L
  )
}
tableOf <- function(...) rbind(..., deparse.level = 0L)
# what a call is refused with, or "created"
refusal <- function(text, ...) {
  tryCatch(
    {
      suppressWarnings(build(text, ...))
      "created"
    },
    error = conditionMessage
  )
}

## --- Block A: a plain forest anywhere ---------------------------------------
# a list on a formula and on a matrix, a specification and a data object
listOf <- function(forests, formula = FALSE) {
  paste0(
    if (formula) "dbarts(y ~ x1 + x2 + x3, frame" else "dbarts(x, y",
    ", forests = list(",
    forests,
    "))"
  )
}
specOf <- function(forests) {
  paste0("dbartsSpec(dbartsData(x, y), forests = list(", forests, "))")
}
dataOf <- function(bases, door = "dbarts") {
  paste0(door, "(dbartsData(x, y, bases = list(", bases, ")))")
}
doorsOf <- function(forests, bases) {
  c(listOf(forests, TRUE), listOf(forests), specOf(forests), dataOf(bases))
}
plainSecond <- doorsOf("forest(basis = dose), forest()", "dose, NULL")
plainFirst <- doorsOf("forest(), forest(basis = dose)", "NULL, dose")
# each is the reversed list's model, forest for forest
for (index in seq_along(plainSecond)) {
  second <- build(plainSecond[index], stated = TRUE)
  first <- build(plainFirst[index], stated = TRUE)
  info <- plainSecond[index]
  expect_identical(
    forestInfo(second)$params,
    rev(forestInfo(first)$params),
    info = info
  )
  expect_identical(forestTable(second), tableOf(multiplied, fitStated))
  expect_identical(forestInfo(second)$vars, rev(forestInfo(first)$vars))
  expect_identical(unname(second$data@bases), rev(unname(first$data@bases)))
  expect_null(second$data@bases[[2L]])
  expect_identical(dim(second$data@bases[[1L]]), c(n, 1L))
  # the half-Cauchy channel is on the plain forest wherever it stands
  expect_identical(second$getLeafPrior(2L)$prior.sd.of, "amplitude scale")
  expect_identical(second$getLeafPrior(1L)$prior.sd.of, "forest total")
  expect_true(all(is.finite(second$run(0L, 5L)$train)), info = info)
}
# it restores, copies and takes a new basis on either forest
second$storeState()
stored <- second$state
ahead <- second$run(0L, 3L)
second$setState(stored)
expect_equal(second$run(0L, 3L)$train, ahead$train)
copied <- second$copy()
expect_identical(forestInfo(copied), forestInfo(second))
expect_identical(dim(copied$run(0L, 3L)$forestFits), c(n, 2L, 3L))
second$setForestBasis(1L, age)
expect_identical(second$data@bases[[1L]][, 1L], age)
expect_null(second$data@bases[[2L]])
second$setForestBasis(2L, dose)
expect_identical(second$data@bases[[2L]][, 1L], dose)
expect_true(all(is.finite(second$run(0L, 3L)$train)))
# and a fit of it predicts: at the training rows, the training fit
plainSecondFit <- bart(
  dbartsData(x, y, bases = list(dose, NULL)),
  n.chains = 1L,
  n.samples = 5L,
  n.burn = 3L,
  verbose = FALSE,
  keepTrees = TRUE,
  seed = 71L
)
expect_equal(
  predict(plainSecondFit, x, bases = list(dose, NULL)),
  plainSecondFit$yhat.train,
  check.attributes = FALSE
)
# two forests with no basis are refused by their positions, wherever they are
twoPlain <- " have no 'basis': a model has one forest with no multiplier, and every other forest states a 'basis'"
for (entry in list(
  c(listOf("forest(), forest()"), "1 and 2"),
  c(listOf("forest(), forest(basis = z), forest()"), "1 and 3"),
  c(listOf("forest(basis = z), forest(), forest()", TRUE), "2 and 3"),
  c(specOf("forest(), forest(), forest()"), "1, 2 and 3"),
  c(dataOf("NULL, NULL"), "1 and 2"),
  c(dataOf("dose, NULL, NULL", "dbartsSpec"), "2 and 3")
)) {
  expect_identical(refusal(entry[1L]), paste0("forests ", entry[2L], twoPlain))
}
# one forest with a basis keeps its own refusal
expect_match(
  refusal(listOf("forest(basis = dose)")),
  "a multi-forest model needs at least two forests",
  fixed = TRUE
)

## --- Block B: the defaults at every door ------------------------------------
# each call with the position of its plain forest, 0 where every forest has a
# basis; the plain forest takes the fitting function's and every other forest
# the multiplied forest's, with nothing stated and with a count and a tree
# prior stated
twoTerms <- "(y ~ forest(x1, basis = dose) + forest(x2, basis = age), frame)"
allBasis <- c(
  paste0("dbarts", twoTerms),
  paste0("bart", twoTerms),
  doorsOf("forest(basis = dose), forest(basis = age)", "dose, age"),
  dataOf("dose, age", "dbartsSpec"),
  dataOf("dose, age", "bart")
)
shapes <- list(
  "dbarts(y ~ x1 + x2 + x3, frame)" = 1L,
  "dbarts(y ~ forest(x1 + x2 + x3), frame)" = 1L,
  "dbarts(x, y, forests = list(forest()))" = 1L,
  "bart(y ~ x1 + x2 + x3, frame)" = 1L,
  "dbartsSpec(dbartsData(x, y))" = 1L,
  "dbarts(y ~ x1 + x2 + forest(x3, basis = dose), frame)" = 1L,
  "bart(y ~ x1 + forest(x3, basis = dose) + forest(x1, basis = age), frame)" = 1L,
  "bart(dbartsData(x, y, bases = list(NULL, dose)))" = 1L,
  # written second or last in a formula it is the first forest still
  "dbarts(y ~ forest(x3, basis = dose) + forest(x1 + x2), frame)" = 1L,
  "bart(y ~ forest(x3, basis = dose) + x1 + x2, frame)" = 1L,
  "bart(dbartsData(x, y, bases = list(dose, NULL)))" = 2L
)
shapes[plainFirst] <- 1L
shapes[plainSecond] <- 2L
shapes[[listOf("forest(basis = dose), forest(basis = age), forest()")]] <- 3L
shapes[allBasis] <- 0L
for (shape in names(shapes)) {
  plain <- shapes[[shape]]
  unstated <- forestTable(build(shape))
  expected <- multiplied[col(unstated)]
  dim(expected) <- dim(unstated)
  expected[plain, ] <- fitDefault
  expect_identical(unstated, expected, info = shape)
  if (plain > 0L) {
    expected[plain, ] <- fitStated
    expect_identical(forestTable(build(shape, TRUE)), expected, info = shape)
  }
}
# a forest's own statement governs it and each number left out takes the
# default of its kind, one by one
expect_identical(
  forestTable(build(listOf(
    "forest(basis = dose, base = 0.5), forest(basis = age)"
  ))),
  tableOf(c(50, 0.5, 3), multiplied)
)
expect_identical(
  forestTable(build(
    listOf("forest(basis = dose, n.trees = 9L), forest()"),
    stated = TRUE
  )),
  tableOf(c(9, 0.25, 3), fitStated)
)
# with nothing stated a model of forests that all have a basis is, draw for
# draw, the model with the three numbers written on every forest
drawsOf <- function(text, ...) {
  sampler <- build(text, ...)
  run <- sampler$run(0L, 5L)
  list(run$train, run$sigma, sampler$getForestAmplitudes())
}
writtenOut <- function(text) {
  gsub("(basis = [a-z]+)", "\\1, n.trees = 50L, base = 0.25, power = 3", text)
}
for (text in c(
  allBasis[1:2],
  specOf("forest(basis = dose), forest(basis = age)"),
  "dbarts(x, yBinary, forests = list(forest(basis = dose), forest(basis = age)), family = \"probit\")"
)) {
  expect_identical(drawsOf(text), drawsOf(writtenOut(text)), info = text)
}
# and a plain forest second under a control naming 15 is the list that states
# every forest's three numbers under an untouched control
expect_identical(
  drawsOf(plainSecond[2L], control = "defaultsControl(n.trees = 15L)"),
  drawsOf(listOf(paste(
    "forest(basis = dose, n.trees = 50L, base = 0.25, power = 3),",
    "forest(n.trees = 15L, base = 0.95, power = 2)"
  )))
)

## --- Block C: what counts as stated -----------------------------------------
# stated means named, the default's own value included
noPlain <- " the forest with no basis, and every forest of this model has a basis; "
countText <- paste0(
  "'n.trees' given to the fitting function is the tree count of",
  noPlain,
  "state a count on a forest, as forest(x1, basis = a, n.trees = 100)"
)
controlText <- paste0(
  "the control's 'n.trees' is the tree count of",
  noPlain,
  "state a count on a forest, as forest(x1, basis = a, n.trees = 100), and ",
  "leave it out of dbartsControl()"
)
onBart <- allBasis[2L]
expect_identical(refusal(onBart), "created")
expect_identical(refusal(onBart, n.trees = 75L), countText)
expect_identical(refusal(onBart, n.trees = 100L), countText)
expect_identical(refusal(onBart, n.tr = 75L), countText)
# through a wrapper that forwards its own default, and through do.call
bartArguments <- list(
  y ~ forest(x1, basis = dose) + forest(x2, basis = age),
  frame,
  n.trees = 75L,
  n.chains = 1L,
  verbose = FALSE,
  samplerOnly = TRUE
)
forwarding <- function(n.trees = 75L) {
  bart(bartArguments[[1L]], frame, n.trees = n.trees, n.chains = 1L)
}
expect_error(forwarding(), countText, fixed = TRUE)
expect_error(do.call(bart, bartArguments), countText, fixed = TRUE)
edited <- function(count) {
  control <- defaultsControl()
  control@n.trees <- count
  control
}
blank <- new("dbartsControl")
blank@n.chains <- 1L
blank@n.threads <- 1L
blank@n.samples <- 5L
expect_identical(refusal(onBart, control = "dbartsControl()"), "created")
expect_identical(
  refusal(onBart, control = "dbartsControl(n.trees = 75L)"),
  controlText
)
expect_identical(refusal(onBart, control = "edited(30L)"), controlText)
for (text in allBasis[c(1L, 7L)]) {
  # with no control at all; nothing is run
  expect_identical(
    lapply(forestInfo(suppressWarnings(eval(str2lang(text))))$params, `[`, 1:3),
    list(multiplied, multiplied)
  )
  for (named in c(
    "defaultsControl(n.trees = 75L)",
    "defaultsControl(n.trees = 20L)",
    "do.call(defaultsControl, list(n.trees = 75L))",
    "edited(30L)"
  )) {
    expect_identical(refusal(text, control = named), controlText, info = named)
  }
  # the limit of what can be told: a slot edited back to the default and a
  # control made by new() read as not stated
  for (unnamed in c("defaultsControl()", "edited(75L)", "blank")) {
    expect_identical(
      refusal(text, control = unnamed),
      "created",
      info = unnamed
    )
  }
  # the defaults of the two priors, named
  expect_match(
    refusal(text, tree.prior = "cgm"),
    "'tree.prior' given to the fitting function is the tree prior of",
    fixed = TRUE
  )
  expect_match(
    refusal(text, leaf.prior = "normal"),
    "'leaf.prior' given to the fitting function is the leaf prior of",
    fixed = TRUE
  )
}
# bart()'s record of what its caller stated does not stay on the control
for (text in c(onBart, "bart(y ~ x1 + x2 + x3, frame, n.trees = 9L)")) {
  expect_null(attr(build(text)$control, "dbarts.stated", exact = TRUE))
}
# bartBT() builds its control and both priors too, so it records what its own
# caller stated, under its own names. The data object is its only door to
# several forests, and what it warns of there is not this file's
onBartBT <- function(given, ...) {
  withCallingHandlers(
    bartBT(
      dbartsData(x, y, bases = given),
      NULL,
      verbose = FALSE,
      sampleronly = TRUE,
      ...
    ),
    warning = function(condition) {
      if (startsWith(conditionMessage(condition), "if data supplied as")) {
        invokeRestart("muffleWarning")
      }
    }
  )
}
refusalBT <- function(...) {
  tryCatch(
    {
      onBartBT(list(dose, age), ...)
      "created"
    },
    error = conditionMessage
  )
}
unstatedBT <- onBartBT(list(dose, age))
expect_identical(forestTable(unstatedBT), tableOf(multiplied, multiplied))
expect_null(attr(unstatedBT$control, "dbarts.stated", exact = TRUE))
# the count a later fit inherits is bartBT's own default
expect_identical(forestInfo(unstatedBT)$control.n.trees, 200L)
expect_identical(refusalBT(keepcall = FALSE), "created")
countTextBT <- sub("'n.trees'", "'ntree'", countText, fixed = TRUE)
expect_identical(refusalBT(ntree = 200L), countTextBT)
expect_identical(refusalBT(ntree = 20L), countTextBT)
expect_identical(refusalBT(ntr = 20L), countTextBT)
expect_identical(refusalBT(ntree = 20L, keepcall = FALSE), countTextBT)
expect_error(
  do.call(onBartBT, list(list(dose, age), ntree = 200L)),
  countTextBT,
  fixed = TRUE
)
statedBT <- list(
  list("power", "tree", power = 2),
  list("base", "tree", base = 0.95),
  list("power", "tree", base = 0.95, power = 2),
  list("splitprobs", "tree", splitprobs = c(0.5, 0.25, 0.25)),
  list("k", "leaf", k = 2)
)
for (case in statedBT) {
  expect_match(
    do.call(refusalBT, case[-(1:2)]),
    paste0(
      "'",
      case[[1L]],
      "' given to the fitting function is the ",
      case[[2L]],
      " prior of the forest with no basis"
    ),
    fixed = TRUE
  )
}
# beside a plain forest, at either place, they are that forest's
fitBT <- c(200, 0.95, 2)
statedOnBT <- list(ntree = 20L, power = 3, base = 0.8)
plainBT <- c(20, 0.8, 3)
expect_identical(
  forestTable(onBartBT(list(dose, NULL))),
  tableOf(multiplied, fitBT)
)
expect_identical(
  forestTable(onBartBT(list(NULL, dose))),
  tableOf(fitBT, multiplied)
)
expect_identical(
  forestTable(do.call(onBartBT, c(list(list(dose, NULL)), statedOnBT))),
  tableOf(multiplied, plainBT)
)
expect_identical(
  forestTable(do.call(onBartBT, c(list(list(NULL, dose)), statedOnBT))),
  tableOf(plainBT, multiplied)
)

## --- Block D: no plain forest -----------------------------------------------
# each of the five, at every door, with its text; the retired names as written
given <- function(argument, what, remedy) {
  paste0(
    "'",
    argument,
    "' given to the fitting function is ",
    what,
    noPlain,
    remedy
  )
}
onForest <- "state it on a forest, as forest(x1, basis = a, "
treeRemedy <- paste0(onForest, "base = 0.25, power = 3)")
leafRemedy <- "state a forest's size on the forest, as forest(x1, basis = a, sd = 2)"
someBlocks <- "blocks(list(\"x1\", c(\"x2\", \"x3\")))"
for (text in allBasis[c(1L, 2L, 4L, 6L, 7L)]) {
  expect_identical(
    refusal(text, stated = TRUE),
    if (startsWith(text, "bart")) countText else controlText,
    info = text
  )
  expect_identical(
    refusal(text, tree.prior = "cgm(3, 0.8)"),
    given("tree.prior", "the tree prior of", treeRemedy),
    info = text
  )
  expect_identical(
    refusal(text, leaf.prior = "normal(k = 2)"),
    given("leaf.prior", "the leaf prior of", leafRemedy),
    info = text
  )
  expect_identical(
    refusal(text, interactions = "interactions(max.order = 1L)"),
    given(
      "interactions",
      "a constraint on",
      paste0(onForest, "interactions = interactions(max.order = 1))")
    ),
    info = text
  )
  expect_identical(
    refusal(text, blocks = someBlocks),
    given(
      "blocks",
      "a constraint on",
      "state it on a forest, as forest(x1 + x2, basis = a, blocks = blocks(list(\"x1\", \"x2\")))"
    ),
    info = text
  )
  # of two given, the first in order
  expect_identical(
    refusal(text, blocks = someBlocks, leaf.prior = "normal"),
    given("leaf.prior", "the leaf prior of", leafRemedy),
    info = text
  )
}
for (entry in list(
  list("power", 2, "the tree prior of", treeRemedy),
  list("base", 0.9, "the tree prior of", treeRemedy),
  list("split.probs", "c(x1 = 1, x2 = 1)", "the tree prior of", treeRemedy),
  list("k", 2, "the leaf prior of", leafRemedy),
  list("node.prior", "normal", "the leaf prior of", leafRemedy)
)) {
  door <- allBasis[if (entry[[1L]] == "node.prior") 1L else 2L]
  expect_identical(
    do.call(refusal, c(list(door), stats::setNames(entry[2L], entry[[1L]]))),
    given(entry[[1L]], entry[[3L]], entry[[4L]])
  )
}
# moved to a forest each is taken, and beside a plain forest none is refused
moved <- build(listOf(paste0(
  "forest(basis = dose, n.trees = 15L, base = 0.8, power = 3, sd = 2, ",
  "interactions = interactions(max.order = 1L), blocks = ",
  someBlocks,
  "), forest(basis = age)"
)))
expect_identical(forestTable(moved), tableOf(fitStated, multiplied))
expect_identical(forestInfo(moved)$params[[1L]][4L], 2)
besidePlain <- build(
  plainSecond[2L],
  stated = TRUE,
  leaf.prior = "normal",
  interactions = "interactions(max.order = 1L)",
  blocks = someBlocks
)
expect_identical(forestTable(besidePlain), tableOf(multiplied, fitStated))

## --- Block E: 'interactions' and 'blocks' follow the plain forest -----------
# given to the fitting function they are the plain forest's, at its position
# and with its tree count, and never the first forest's for being first
constraints <- forestInfo(besidePlain)
expect_null(constraints$interactions[[1L]])
expect_identical(constraints$interactions[[2L]]$max.order, 1L)
expect_null(attr(besidePlain$model, "interaction.max.order", exact = TRUE))
expect_null(constraints$blocks[[1L]])
expect_identical(sum(constraints$blocks[[2L]]$block.tree.counts), 15L)
expect_null(attr(besidePlain$model, "block.of.column", exact = TRUE))
# the largest number of predictors on one path of a tree, its nodes depth first
pathPredictors <- function(vars) {
  position <- 0L
  walk <- function(seen) {
    position <<- position + 1L
    var <- vars[position]
    if (var < 0L) {
      return(length(seen))
    }
    seen <- union(seen, var)
    max(walk(seen), walk(seen))
  }
  walk(integer())
}
additive <- build(
  listOf(
    "forest(basis = dose, n.trees = 10L, base = 0.95, power = 2), forest()"
  ),
  control = "defaultsControl(n.trees = 10L, keepTrees = TRUE)",
  interactions = "interactions(max.order = 1L)"
)
additive$run(0L, 300L)
orders <- lapply(1:2, function(index) {
  trees <- additive$getTrees(forest = index)
  vapply(
    split(trees$var, list(trees$sample, trees$tree), drop = TRUE),
    pathPredictors,
    0L
  )
})
expect_true(max(orders[[1L]]) > 1L)
expect_identical(max(orders[[2L]]), 1L)
# on the first, multiplied forest and at the top: two constraints, two forests
twoConstraints <- build(
  listOf(
    "forest(basis = dose, interactions = interactions(max.order = 2L)), forest()"
  ),
  interactions = "interactions(max.order = 1L)"
)
expect_identical(forestInfo(twoConstraints)$interactions[[1L]]$max.order, 2L)
expect_identical(forestInfo(twoConstraints)$interactions[[2L]]$max.order, 1L)
expect_identical(attr(twoConstraints$model, "interaction.max.order"), 2L)
# on the plain forest and at the top: one constraint twice, at either position
for (constraint in list(
  c("interactions", "interactions(max.order = 1L)"),
  c("blocks", someBlocks)
)) {
  own <- paste0("forest(", constraint[1L], " = ", constraint[2L], ")")
  atTop <- stats::setNames(constraint[2L], constraint[1L])
  for (forests in c(
    paste0(own, ", forest(basis = dose)"),
    paste0("forest(basis = dose), ", own)
  )) {
    expect_identical(
      do.call(refusal, c(list(listOf(forests)), atTop)),
      paste0(
        "'",
        constraint[1L],
        "' is given to the fitting function and to the forest with no basis, ",
        "which are the same constraint; give one"
      )
    )
  }
}

## --- Block F: the count given twice at bart() -------------------------------
# judged on the model, so a forest that states no basis in other words is the
# plain forest too
countTwice <- "'n.trees' is given to the fitting function and to the forest with no basis"
nothing <- NULL
expect_identical(
  refusal(
    "bart(y ~ forest(x1 + x2, n.trees = 7L) + forest(x3, basis = dose), frame)",
    n.trees = 15L
  ),
  paste0(
    countTwice,
    " ('forest(x1 + x2, n.trees = 7L)'), which are the same count; give one"
  )
)
for (text in c(
  "bart(y ~ forest(x1 + x2, n.trees = 7L, basis = NULL) + forest(x3, basis = dose), frame)",
  "bart(y ~ forest(x1 + x2, n.trees = 7L, basis = nothing) + forest(x3, basis = dose), frame)",
  "bart(y ~ forest(x1 + x2, n.trees = 7L, basis = NULL), frame)"
)) {
  expect_identical(
    refusal(text, n.trees = 15L),
    paste0(countTwice, ", which are the same count; give one")
  )
  expect_identical(forestTable(build(text))[1L, 1L], 7)
}
# a multiplied forest's count beside bart()'s is two counts
expect_identical(
  forestTable(build(
    "bart(y ~ x1 + x2 + forest(x3, basis = dose, n.trees = 7L), frame)",
    n.trees = 15L
  ))[, 1L],
  c(15, 7)
)

## --- Block G: a control carried to another fit ------------------------------
# the control's slot holds the first forest's count, and a multiplied forest's
# count is never the next fit's
carriedTo <- function(sampler, text = "dbarts(y ~ x1 + x2, frame)") {
  call <- str2lang(text)
  call$control <- quote(sampler$control)
  forestTable(eval(call))[, 1L]
}
allBasisSampler <- build(allBasis[1L])
expect_identical(allBasisSampler$control@n.trees, 50L)
expect_identical(carriedTo(allBasisSampler), 75)
expect_identical(carriedTo(allBasisSampler, allBasis[1L]), c(50, 50))
namedSecond <- build(plainSecond[2L], stated = TRUE)
expect_identical(namedSecond$control@n.trees, 50L)
expect_identical(carriedTo(namedSecond), 15)
ownSecond <- build(listOf("forest(basis = dose), forest(n.trees = 7L)"))
ownFirst <- build(listOf("forest(n.trees = 7L), forest(basis = dose)"))
expect_identical(carriedTo(ownSecond), 7)
expect_identical(carriedTo(ownFirst), 7)
# a control from a fit whose plain forest is first holds its own count already
expect_null(forestInfo(ownFirst)$control.n.trees)
# a count edited on a carried control is the caller's own: it is the next
# fit's plain forest's wherever that stands, and where every forest has a
# basis it is refused as any count on a control is
editedCarry <- function(sampler, count, text = "dbarts(y ~ x1 + x2, frame)") {
  control <- sampler$control
  control@n.trees <- count
  call <- str2lang(text)
  call$control <- control
  built <- eval(call)
  if (is.list(built)) {
    built <- specSampler(built)
  }
  forestTable(built)[, 1L]
}
for (carrier in list(allBasisSampler, namedSecond)) {
  expect_identical(editedCarry(carrier, 20L), 20)
  expect_identical(
    editedCarry(carrier, 20L, "dbartsSpec(dbartsData(x, y))"),
    20
  )
  expect_identical(editedCarry(carrier, 20L, plainFirst[2L]), c(20, 50))
  expect_identical(editedCarry(carrier, 20L, plainSecond[2L]), c(50, 20))
  expect_error(
    editedCarry(carrier, 20L, allBasis[1L]),
    controlText,
    fixed = TRUE
  )
}
# an edit to anything else leaves the carried count to go back
carriedOther <- allBasisSampler$control
carriedOther@n.samples <- 6L
expect_identical(
  forestTable(dbarts(y ~ x1 + x2, frame, control = carriedOther))[, 1L],
  75
)
# the limit of what can be told: a slot edited to the first forest's own
# count reads as untouched, and one edited to the default as not stated
expect_identical(editedCarry(allBasisSampler, 50L), 75)
expect_identical(editedCarry(allBasisSampler, 75L, allBasis[1L]), c(50, 50))
# and $setControl takes the control a sampler was created under, or its own,
# the slot keeping the first forest's count
allBasisSampler$setControl(defaultsControl(n.burn = 7L))
expect_identical(allBasisSampler$control@n.trees, 50L)
expect_identical(allBasisSampler$control@n.burn, 7L)
allBasisSampler$setControl(allBasisSampler$control)
namedSecond$setControl(defaultsControl(n.trees = 15L))
expect_identical(namedSecond$control@n.trees, 50L)
expect_identical(carriedTo(namedSecond), 15)
for (count in c(20L, 75L)) {
  expect_error(
    namedSecond$setControl(defaultsControl(n.trees = count)),
    "changing 'n.trees' is not available on an existing sampler",
    fixed = TRUE
  )
}

## --- Block H: print ---------------------------------------------------------
# a count for each forest, in the forests' order
treeLine <- function(fit) {
  printed <- utils::capture.output(print(fit))
  printed[startsWith(printed, "n.trees:")]
}
expect_identical(treeLine(plainSecondFit), "n.trees: 50, 75")
savedFit <- tempfile(fileext = ".rds")
saveRDS(plainSecondFit, savedFit)
expect_identical(treeLine(readRDS(savedFit)), "n.trees: 50, 75")
unlink(savedFit)
oneForestFit <- bart(
  y ~ x1 + x2 + x3,
  frame,
  n.trees = 9L,
  n.chains = 1L,
  n.samples = 3L,
  n.burn = 0L,
  verbose = FALSE,
  keepTrees = TRUE
)
expect_identical(treeLine(oneForestFit), "n.trees: 9")

## --- Block I: which held coefficients are taken -----------------------------
# Until a held coefficient has a value of its own for every shape, one is held
# only where the engine holds it at the value the help states, which goes by
# the width of the forest's basis and the forest's position.
three <- unname(stats::model.matrix(~ factor(rep_len(1:3, n)) - 1))
bases <- list(plain = NULL, one = dose, two = factor(z), three = three)
atDataDoor <- list(plain = NULL, one = dose, two = zBasis, three = three)
heldForests <- function(kinds, held, values = bases) {
  lapply(seq_along(kinds), function(index) {
    amplitude <- if (index == held) dbarts::dbartsPriors$fixed()
    basis <- values[[kinds[index]]]
    if (is.null(basis)) {
      dbarts::dbartsForests$forest(amplitude = amplitude)
    } else {
      dbarts::dbartsForests$forest(basis = basis, amplitude = amplitude)
    }
  })
}
heldDoors <- list(
  formula = function(kinds, held) {
    forests <- heldForests(kinds, held)
    dbarts(
      y ~ x1 + x2 + x3,
      frame,
      forests = forests,
      control = defaultsControl()
    )
  },
  matrix = function(kinds, held) {
    dbarts(
      x,
      y,
      forests = heldForests(kinds, held),
      control = defaultsControl()
    )
  },
  spec = function(kinds, held) {
    forests <- heldForests(kinds, held)
    specSampler(dbartsSpec(
      dbartsData(x, y),
      defaultsControl(),
      forests = forests
    ))
  },
  # two numeric columns reach the data door as a factor's two indicator
  # columns do, and held second they are taken: the assertion turns around
  # when a basis of several numeric columns is refused
  data = function(kinds, held) {
    dbarts(
      dbartsData(x, y, bases = unname(atDataDoor[kinds])),
      forests = heldForests(kinds, held, list()),
      control = defaultsControl()
    )
  }
)
fixedOn <- function(index, what) {
  paste0("forest ", index, ": amplitude = fixed() on ", what)
}
twoElsewhere <- "a basis of two columns is not supported yet unless the forest is the second; it would hold both coefficients at 1. Put the forest second, or let the coefficients be drawn"
plainHeldSecond <- "a forest with no basis is not supported yet where it is the second forest; it would hold the forest at zero. Put the forest with no basis first, or let the coefficient be drawn"
threeSecond <- "a basis of 3 columns is not supported; it would hold every column but the first at 1. Let the coefficients be drawn"
threeThird <- "a basis of 3 columns is not supported; it would hold every column at 1. Let the coefficients be drawn"
# the forests' kinds in order and which is held, with the values it is held
# at, or the refusal and whether the model is one to create with no hold
heldRows <- list(
  list(c("plain", "two"), 1L, values = 1),
  list(c("one", "two", "plain"), 3L, values = 1),
  list(c("plain", "two"), 2L, values = c(0, 1)),
  list(c("one", "two"), 2L, values = c(0, 1)),
  list(
    c("one", "plain"),
    2L,
    text = fixedOn(2L, plainHeldSecond),
    drawn = TRUE
  ),
  list(c("two", "one"), 1L, text = fixedOn(1L, twoElsewhere), drawn = TRUE),
  list(
    c("plain", "one", "two"),
    3L,
    text = fixedOn(3L, twoElsewhere),
    drawn = TRUE
  ),
  list(c("plain", "three"), 2L, text = fixedOn(2L, threeSecond), drawn = FALSE),
  list(
    c("plain", "one", "three"),
    3L,
    text = fixedOn(3L, threeThird),
    drawn = FALSE
  )
)
for (door in names(heldDoors)) {
  for (row in heldRows) {
    info <- paste(door, paste(row[[1L]], collapse = " "), "held", row[[2L]])
    if (is.null(row$text)) {
      sampler <- heldDoors[[door]](row[[1L]], row[[2L]])
      before <- sampler$getForestAmplitudes(row[[2L]])[, 1L]
      sampler$run(0L, 20L)
      after <- sampler$getForestAmplitudes(row[[2L]])[, 1L]
      expect_equal(before, row$values, info = info)
      expect_equal(after, row$values, info = info)
    } else {
      expect_error(
        heldDoors[[door]](row[[1L]], row[[2L]]),
        row$text,
        fixed = TRUE,
        info = info
      )
      # with the hold dropped the model is created
      if (row$drawn) {
        drawn <- heldDoors[[door]](row[[1L]], 0L)
        expect_identical(length(drawn$data@bases), length(row[[1L]]))
      }
    }
  }
}
# as formula terms, where a formula can write the shape
termHeld <- function(formula) {
  sampler <- dbarts(formula, frame, control = defaultsControl())
  sampler$getForestAmplitudes()[, 1L]
}
expect_equal(
  termHeld(
    y ~ forest(x1 + x2, amplitude = fixed()) +
      forest(x3, basis = factor(z), amplitude = fixed())
  ),
  c(1, 0, 1)
)
expect_error(
  termHeld(
    y ~ x1 +
      x2 +
      forest(x1, basis = dose) +
      forest(x3, basis = factor(z), amplitude = fixed())
  ),
  fixedOn(3L, twoElsewhere),
  fixed = TRUE
)
expect_error(
  termHeld(
    y ~ forest(x1, basis = factor(z), amplitude = fixed()) +
      forest(x3, basis = dose)
  ),
  fixedOn(1L, twoElsewhere),
  fixed = TRUE
)
