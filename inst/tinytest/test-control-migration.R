# The model facts a sampler runs ride its model and data objects, never its
# control: each forest's record, the variance forest's record, the category
# count, the shape setting and the hazard periods on the model; the
# censoring status and the quantile rule on the data. Each is read back at a
# value other than its default from the engine on the live sampler, after
# copy() and after a save and reload. The control's n.trees, useQuantiles and
# binary are copies the sampler writes from the model and data.

forest <- dbartsForests$forest
varianceForest <- dbartsForests$varianceForest

set.seed(23)
n <- 60L
frame <- data.frame(x1 = runif(n), x2 = runif(n), x3 = runif(n))
frame$z <- rbinom(n, 1L, 0.5)
frame$zn <- rnorm(n)
frame$y <- with(frame, x1 + z * (1 + x2) + rnorm(n, sd = 0.3))
x <- as.matrix(frame[c("x1", "x2", "x3")])
y <- frame$y
yb <- as.numeric(y > median(y))

migrationControl <- function(...) {
  dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.samples = 4L,
    n.burn = 2L,
    seed = 9L,
    ...
  )
}

# a copy and a reload of a sampler, each from its stored state
twins <- function(sampler, run = TRUE) {
  if (run) {
    sampler$run()
  }
  sampler$storeState()
  file <- tempfile(fileext = ".rds")
  on.exit(unlink(file))
  saveRDS(sampler, file)
  list(live = sampler, copy = sampler$copy(), reload = readRDS(file))
}
forestCounts <- function(sampler) {
  ptr <- sampler$getPointer()
  vapply(
    seq_len(dbarts:::bartcoreNumForests(ptr)),
    function(index) dbarts:::bartcoreForestTreeCount(ptr, index - 1L),
    0L
  )
}
varianceCount <- function(sampler) {
  sampler$storeState()
  length(sampler$state[[1L]]$variance.sizes)
}
countWarnings <- function(expr) {
  count <- 0L
  withCallingHandlers(expr, warning = function(w) {
    count <<- count + 1L
    invokeRestart("muffleWarning")
  })
  count
}

## ---- the variance forest's record -------------------------------------------

varianceSampler <- dbarts(
  y ~ x1 + x2 + x3,
  frame,
  variance = varianceForest(
    vars = ~ x1 + x2,
    n.trees = 7L,
    base = 0.8,
    power = 1.5
  ),
  control = migrationControl()
)
expect_null(attr(varianceSampler$control, "bartcore.variance"))
record <- varianceSampler$model@variance
expect_true(is(record, "dbartsVarianceRecord"))
expect_identical(record@n.trees, 7L)
expect_identical(record@columns, 1:2)
expect_equal(c(record@base, record@power), c(0.8, 1.5))
for (twin in twins(varianceSampler)) {
  expect_identical(varianceCount(twin), 7L)
  expect_identical(twin$model@variance, record)
}
# the record's prior is what the engine draws under: a default one draws
# otherwise
varianceOf <- function(...) {
  dbarts(
    y ~ x1 + x2 + x3,
    frame,
    variance = varianceForest(vars = ~ x1 + x2, n.trees = 7L, ...),
    control = migrationControl()
  )$run()$variance
}
expect_false(identical(
  varianceOf(),
  varianceOf(base = 0.8, power = 1.5)
))
rm(varianceSampler, varianceOf, record, twin)

## ---- the shape setting and the category count ---------------------------------

counts <- rpois(n, exp(x[, 1L]))
countSampler <- dbarts(
  x,
  counts,
  family = dbartsFamilies$nbinom(shape = 3),
  control = migrationControl()
)
expect_identical(countSampler$model@shape, 3)
for (twin in twins(countSampler)) {
  expect_equal(twin$getShape(), 3)
}
ordered <- cut(
  y,
  quantile(y, c(0, 0.25, 0.5, 0.75, 1)),
  include.lowest = TRUE,
  ordered_result = TRUE
)
ordinalSampler <- dbarts(
  x,
  ordered,
  family = "ordinal",
  control = migrationControl()
)
expect_identical(ordinalSampler$model@n.categories, 4L)
for (twin in twins(ordinalSampler)) {
  expect_identical(NROW(twin$run()$thresholds), 3L)
}
rm(counts, countSampler, ordered, ordinalSampler, twin)

## ---- the forests' records ------------------------------------------------------

moving <- c(birth_death = 0.5, swap = 0.1, change = 0.4, birth = 0.5)
threeForestsOf <- function(effectMoves) {
  dbarts(
    y ~ x1 + x2 + x3,
    frame,
    forests = list(
      plain = forest(n.trees = 9L, interactions = interactions(max.order = 1L)),
      effect = forest(
        x1 + x2,
        basis = ~ factor(z),
        n.trees = 6L,
        sd = 0.4,
        blocks = blocks(list("x1", "x2")),
        proposal.probs = effectMoves
      ),
      dose = forest(x3, basis = ~zn, n.trees = 4L, amplitude = fixed())
    ),
    control = migrationControl(keepTrees = TRUE)
  )
}
threeForests <- threeForestsOf(moving)
termless <- function(record) {
  record@basis.term <- NULL
  record
}
records <- threeForests$model@forests
expect_identical(length(records), 3L)
expect_true(all(vapply(records, is, NA, "dbartsForestRecord")))
expect_identical(vapply(records, function(r) r@n.trees, 0L), c(9L, 6L, 4L))
expect_identical(
  vapply(records, function(r) r@label, ""),
  c("plain", "effect", "dose")
)
expect_identical(records[[2L]]@vars, 1:2)
expect_identical(records[[3L]]@vars, 3L)
expect_equal(records[[2L]]@sd, 0.4)
expect_identical(
  vapply(records, function(r) r@held.amplitude, NA),
  c(FALSE, FALSE, TRUE)
)
expect_identical(records[[1L]]@interactions$max.order, 1L)
expect_false(is.null(records[[2L]]@blocks))
expect_identical(records[[2L]]@moves$proposal.probs[["swap"]], 0.1)
expect_identical(vapply(records, function(r) r@stated, ""), rep("n.trees", 3L))
expect_null(attr(threeForests$control, "bartcore.forests"))
for (twin in twins(threeForests)) {
  expect_identical(forestCounts(twin), c(9L, 6L, 4L))
  # a reload rebuilds a basis term's environment, which holds nothing
  expect_identical(
    lapply(twin$model@forests, termless),
    lapply(records, termless)
  )
  trees <- twin$getTrees(forest = "effect")
  expect_true(all(trees$var[trees$var > 0L] %in% 1:2))
  trees <- twin$getTrees(forest = 3L)
  expect_true(all(trees$var[trees$var > 0L] == 3L))
  expect_equal(c(twin$getForestAmplitudes(forest = "dose")), 1)
  expect_equal(unname(twin$getLeafPrior(forest = "effect")$sd), 0.4)
  expect_identical(
    twin$getLeafPrior(forest = "dose")$amplitude,
    dbartsPriors$fixed()
  )
}
# the moves install reads the record: the same forests stating none of their
# own draw otherwise
withoutMoves <- threeForestsOf(NULL)
expect_null(withoutMoves$model@forests[[2L]]@moves)
expect_false(identical(
  withoutMoves$run()$train,
  threeForestsOf(moving)$run()$train
))
rm(withoutMoves, trees, twin, threeForestsOf, termless)

# a student() coefficient's df rides its forest's record
studentFit <- function(df) {
  dbarts(
    y ~ x1 + x2 + x3,
    frame,
    forests = list(
      forest(),
      forest(basis = ~ factor(z), amplitude = student(df))
    ),
    control = migrationControl()
  )
}
heavy <- studentFit(3)
expect_equal(heavy$model@forests[[2L]]@amplitude.df, 3)
expect_true(is.na(heavy$model@forests[[1L]]@amplitude.df))
expect_false(identical(
  studentFit(3)$run()$train,
  studentFit(30)$run()$train
))
copies <- lapply(twins(heavy), function(twin) {
  set.seed(31L)
  twin$run()$train
})
expect_identical(copies$reload, copies$copy)
expect_equal(copies$copy, copies$live)
rm(studentFit, heavy, copies)

# a basis written as code is rebuilt at new rows from its record
codeFit <- bart(
  y ~ forest(x1 + x2 + x3) + forest(x1, basis = ~ factor(z)),
  frame,
  n.samples = 4L,
  n.burn = 2L,
  n.chains = 1L,
  n.threads = 1L,
  seed = 9L,
  verbose = FALSE,
  keepTrees = TRUE
)
expect_false(is.null(codeFit$fit$model@forests[[2L]]@basis.term))
expect_identical(
  codeFit$basis.terms[[2L]],
  codeFit$fit$model@forests[[2L]]@basis.term
)
expect_identical(dim(predict(codeFit, frame[1:5, ])), c(4L, 5L))
rm(codeFit)

# a multinomial model's one record sets every category forest's count
categories <- factor(sample(c("a", "b", "c"), n, TRUE))
multinomialSampler <- dbarts(
  x,
  categories,
  family = "multinomial",
  control = migrationControl(n.trees = 8L)
)
expect_identical(length(multinomialSampler$model@forests), 1L)
expect_identical(multinomialSampler$model@forests[[1L]]@n.trees, 8L)
for (twin in twins(multinomialSampler)) {
  expect_identical(forestCounts(twin), rep(8L, 3L))
}
rm(categories, multinomialSampler, twin)

## ---- the status and the hazard periods ---------------------------------------

time <- rexp(n, exp(-x[, 1L]))
status <- rbinom(n, 1L, 0.6)
aftOf <- function(status) {
  dbarts(x, cbind(time, status), family = "aft", control = migrationControl())
}
expected <- aftOf(rep(1, n))$run()$train
expect_false(identical(aftOf(status)$run()$train, expected))
setTo <- aftOf(status)
expect_identical(setTo$data@status, as.double(status))
expect_null(attr(setTo$control, "bartcore.survival"))
setTo$setResponse(setTo$data@y, status = rep(1, n))
expect_identical(setTo$data@status, rep(1, n))
copies <- twins(setTo, run = FALSE)
for (twin in copies[c("copy", "reload")]) {
  expect_identical(twin$data@status, rep(1, n))
  expect_identical(twin$run()$train, expected)
}
rm(aftOf, expected, setTo, copies, twin)

hazardFit <- bart(
  x,
  cbind(ceiling(time * 3), status),
  family = "hazard",
  n.samples = 4L,
  n.burn = 2L,
  n.chains = 1L,
  n.threads = 1L,
  seed = 9L,
  verbose = FALSE,
  keepSampler = TRUE
)
expect_false(is.null(hazardFit$fit$model@hazard.periods))
expect_identical(hazardFit$periods, hazardFit$fit$model@hazard.periods)
expect_null(attr(hazardFit$fit$control, "bartcore.hazard.periods"))
for (twin in twins(hazardFit$fit, run = FALSE)) {
  expect_identical(twin$model@hazard.periods, hazardFit$periods)
  expect_identical(dbarts:::pdbart.samplerFamily(twin), "hazard")
}
rm(hazardFit, time, status, twin)

## ---- the mirrors ---------------------------------------------------------------

mirrored <- function(sampler) {
  sampler$control@n.trees == sampler$model@forests[[1L]]@n.trees &&
    sampler$control@useQuantiles == isTRUE(sampler$data@useQuantiles) &&
    sampler$control@binary == dbarts:::isBinaryFamily(sampler$model@family)
}
for (family in c("gaussian", "probit", "logistic")) {
  response <- if (family == "gaussian") y else yb
  sampler <- dbarts(
    x,
    response,
    family = family,
    control = migrationControl(n.trees = 13L, useQuantiles = TRUE)
  )
  expect_identical(sampler$model@forests[[1L]]@n.trees, 13L)
  expect_true(sampler$data@useQuantiles)
  expect_identical(sampler$control@binary, family != "gaussian")
  expect_true(mirrored(sampler))
  for (twin in twins(sampler)) {
    expect_true(mirrored(twin))
  }
  # a sampler rebuilt from its stored triple, as stan4bart restores one
  expect_identical(
    countWarnings(
      restored <- new(
        "dbartsSampler",
        sampler$control,
        sampler$model,
        sampler$data
      )
    ),
    0L
  )
  expect_true(mirrored(restored))
}
rm(family, response, sampler, twin, restored)

# the family alone says binary: a model that records no forest, beside a
# control whose binary flag is FALSE, runs its binary family
for (family in c("probit", "logistic")) {
  fit <- dbarts(x, yb, family = family, control = migrationControl())
  bare <- fit$model
  bare@forests <- list()
  control <- migrationControl()
  expect_false(control@binary)
  handSampler <- new("dbartsSampler", control, bare, fit$data)
  expect_true(handSampler$control@binary)
  expect_identical(handSampler$model@forests[[1L]]@n.trees, control@n.trees)
  expect_identical(handSampler$run()$train, fit$run()$train)
}
rm(fit, bare, control, handSampler, family)

## ---- the quantile rule ---------------------------------------------------------

gridOf <- function(sampler) sampler$data@cutPoints[[1L]]
uniform <- dbarts(x, y, control = migrationControl(n.cuts = 7L))
quantile <- dbarts(
  x,
  y,
  control = migrationControl(n.cuts = 7L, useQuantiles = TRUE)
)
expect_false(isTRUE(uniform$data@useQuantiles))
expect_true(quantile$data@useQuantiles)
expect_false(identical(gridOf(uniform), gridOf(quantile)))
for (twin in twins(quantile)) {
  expect_identical(gridOf(twin), gridOf(quantile))
  expect_true(twin$data@useQuantiles)
}
spec <- dbartsSpec(
  dbartsData(x, y),
  migrationControl(n.cuts = 7L, useQuantiles = TRUE)
)
expect_true(spec$data@useQuantiles)
expect_true(spec$control@useQuantiles)
# xbart's folds bin on the handle the data object's rule builds
xbartUniform <- xbart(
  x,
  y,
  n.samples = 3L,
  n.burn = c(2L, 1L),
  n.reps = 1L,
  n.trees = c(3L, 5L),
  n.threads = 1L,
  seed = 4L,
  n.cuts = 3L
)
xbartQuantile <- xbart(
  x,
  y,
  n.samples = 3L,
  n.burn = c(2L, 1L),
  n.reps = 1L,
  n.trees = c(3L, 5L),
  n.threads = 1L,
  seed = 4L,
  n.cuts = 3L,
  useQuantiles = TRUE
)
expect_false(identical(xbartUniform, xbartQuantile))
rm(gridOf, uniform, quantile, spec, twin, xbartUniform, xbartQuantile)

## ---- what a sampler keeps across its setters -----------------------------------

keeper <- dbarts(
  y ~ x1 + x2 + x3,
  frame,
  variance = TRUE,
  control = migrationControl(n.trees = 11L, useQuantiles = TRUE)
)
otherModel <- dbarts(
  y ~ x1 + x2 + x3,
  frame,
  variance = TRUE,
  control = migrationControl(n.trees = 11L)
)$model
otherModel@tree.prior@power <- 1.5
keeper$setModel(otherModel)
expect_identical(keeper$model@forests[[1L]]@n.trees, 11L)
expect_false(is.null(keeper$model@variance))
expect_identical(forestCounts(keeper), 11L)
keeper$setData(dbartsData(y ~ x1 + x2 + x3, frame[sample(n), ]))
expect_true(keeper$data@useQuantiles)
keeper$setControl(
  migrationControl(n.trees = 11L, n.thin = 2L, useQuantiles = TRUE)
)
expect_true(mirrored(keeper))
expect_identical(varianceCount(keeper), 40L)
rm(keeper, otherModel)

## ---- a control copied from a fit states nothing --------------------------------

varianceSampler <- dbarts(
  y ~ x1 + x2 + x3,
  frame,
  variance = varianceForest(n.trees = 7L),
  control = migrationControl()
)

# the not-used warnings a call draws, and what else it draws, counted apart
warned <- function(expr) {
  notUsed <- character(0L)
  value <- withCallingHandlers(expr, dbartsControlNotUsedWarning = function(w) {
    notUsed <<- c(notUsed, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(value = value, notUsed = notUsed)
}
expect_notUsed <- function(call, name, value) {
  expect_identical(length(call$notUsed), length(name))
  for (index in seq_along(name)) {
    expect_true(
      startsWith(
        call$notUsed[index],
        paste0(
          "the control's ",
          name[index],
          " (",
          value[index],
          ") is not used"
        )
      )
    )
  }
}
fit50 <- dbarts(
  y ~ x1 + x2 + x3,
  frame,
  control = migrationControl(n.trees = 50L)
)
fit75 <- dbarts(y ~ x1 + x2 + x3, frame, control = migrationControl())
expect_identical(fit50$control@n.trees, 50L)
expect_false("n.trees" %in% dbarts:::controlSuppliedSlots(fit50$control))

# xbart(), with the count its grid takes
xbartOn <- function(control, ...) {
  xbart(
    x,
    y,
    n.samples = 3L,
    n.burn = c(2L, 1L),
    n.reps = 1L,
    n.threads = 1L,
    seed = 4L,
    ...,
    control = control
  )
}
copied <- warned(xbartOn(fit50$control, n.trees = 5L))
expect_notUsed(copied, "n.trees", 50L)
written <- warned(xbartOn(fit50$control, n.trees = 50L))
expect_identical(written$notUsed, character(0L))
expect_notUsed(
  warned(xbartOn(varianceSampler$control, n.trees = 5L)),
  "n.trees",
  75L
)
rm(copied, written, xbartOn)

## ---- the mirror reads the fitting function's level ------------------------------

plainSecond <- function(control = migrationControl()) {
  dbarts(
    y ~ x1 + x2 + x3,
    frame,
    forests = list(
      forest(x1, basis = ~ factor(z), n.trees = 7L),
      forest(n.trees = 9L)
    ),
    control = control
  )
}
secondSampler <- plainSecond()
expect_identical(forestCounts(secondSampler), c(7L, 9L))
expect_identical(secondSampler$control@n.trees, 9L)
copiedTwins <- twins(secondSampler)
expect_identical(copiedTwins$copy$control@n.trees, 9L)
expect_identical(copiedTwins$reload$control@n.trees, 9L)
expect_identical(copiedTwins$reload$getTrees()$tree[0L], integer(0L))
# given back to its own formula, none
handedBack <- warned(plainSecond(secondSampler$control))
expect_identical(handedBack$notUsed, character(0L))
expect_identical(forestCounts(handedBack$value), c(7L, 9L))
rm(plainSecond, secondSampler, copiedTwins, handedBack)
multinomialFit <- dbarts(
  x,
  factor(cut(y, 3L, labels = FALSE)),
  family = "multinomial",
  control = migrationControl(n.trees = 8L)
)
expect_identical(multinomialFit$control@n.trees, 8L)
expect_false(multinomialFit$control@binary)

## ---- $setControl ----------------------------------------------------------------

livePair <- function() {
  dbarts(y ~ x1 + x2 + x3, frame, control = migrationControl(n.trees = 12L))
}
# a carried count that differs: taken, the sampler's kept, the warning
sampler <- livePair()
carried <- fit50$control
carried@n.thin <- 2L
taken <- warned(sampler$setControl(carried))
expect_notUsed(taken, "n.trees", 50L)
expect_identical(sampler$control@n.trees, 12L)
expect_identical(sampler$control@n.thin, 2L)
expect_identical(forestCounts(sampler), 12L)
expect_identical(sampler$copy()$control@n.trees, 12L)
expect_true(grepl("keeps its 12", taken$notUsed, fixed = TRUE))
# a written one: refused, the sampler untouched, no warning
sampler <- livePair()
before <- sampler$control
expect_error(
  sampler$setControl(migrationControl(n.trees = 50L)),
  "changing 'n.trees' is not available on an existing sampler",
  fixed = TRUE
)
expect_identical(sampler$control, before)
edited <- fit50$control
edited@n.trees <- 20L
refusal <- warned(tryCatch(
  sampler$setControl(edited),
  error = conditionMessage
))
expect_identical(
  refusal$value,
  "changing 'n.trees' is not available on an existing sampler"
)
expect_identical(refusal$notUsed, character(0L))
expect_identical(sampler$control, before)
# the same value written: nothing to change
agrees <- warned(sampler$setControl(migrationControl(
  n.trees = 12L,
  n.thin = 3L
)))
expect_identical(agrees$notUsed, character(0L))
expect_identical(sampler$control@n.thin, 3L)
# a fresh control: taken, the count kept, silent
silent <- warned(sampler$setControl(migrationControl(n.thin = 4L)))
expect_identical(silent$notUsed, character(0L))
expect_identical(sampler$control@n.trees, 12L)
expect_identical(sampler$control@n.thin, 4L)
# the rule and the cut count read alike
sampler <- dbarts(
  y ~ x1 + x2 + x3,
  frame,
  control = migrationControl(useQuantiles = TRUE, n.cuts = 20L)
)
expect_error(
  sampler$setControl(migrationControl(useQuantiles = FALSE, n.cuts = 20L)),
  "changing 'useQuantiles' is not available on an existing sampler",
  fixed = TRUE
)
expect_error(
  sampler$setControl(migrationControl(useQuantiles = TRUE, n.cuts = 30L)),
  "changing 'n.cuts' is not available on an existing sampler",
  fixed = TRUE
)
expect_identical(
  warned(sampler$setControl(migrationControl()))$notUsed,
  character(0L)
)
expect_true(sampler$control@useQuantiles)
expect_identical(sampler$data@n.cuts, rep(20L, 3L))
taken <- warned(sampler$setControl(fit75$control))
expect_notUsed(taken, c("useQuantiles", "n.cuts"), c("FALSE", "100"))
expect_true(sampler$control@useQuantiles)
expect_identical(sampler$data@n.cuts, rep(20L, 3L))
rm(livePair, sampler, carried, taken, before, edited, refusal, agrees, silent)

## ---- $setModel: structure is refused, move settings are installed ---------------

structureRefused <- function(sampler, model, pattern) {
  before <- list(sampler$model, sampler$control)
  expect_error(sampler$setModel(model), pattern, fixed = TRUE)
  expect_identical(list(sampler$model, sampler$control), before)
  expect_true(is.numeric(sampler$run(0L, 1L)$train))
}
sampler <- dbarts(y ~ x1 + x2 + x3, frame, control = migrationControl())
structureRefused(sampler, fit50$model, "forest 1's tree count")
structureRefused(sampler, varianceSampler$model, "variance forest")
structureRefused(
  dbarts(
    y ~ x1 + x2 + x3,
    frame,
    variance = TRUE,
    control = migrationControl()
  ),
  fit75$model,
  "variance forest"
)
varianceOf <- function(...) {
  dbarts(
    y ~ x1 + x2 + x3,
    frame,
    variance = varianceForest(...),
    control = migrationControl()
  )
}
sampler <- varianceOf(n.trees = 7L)
structureRefused(
  sampler,
  varianceOf(n.trees = 8L)$model,
  "variance forest's tree count"
)
structureRefused(
  varianceOf(n.trees = 7L, vars = ~ x1 + x2),
  sampler$model,
  "columns of the variance forest"
)
# a variance forest's prior is a parameter: not refused, the sampler's kept
parameter <- varianceOf(n.trees = 7L, base = 0.8, power = 1.5)$model
sampler$setModel(parameter)
expect_identical(sampler$model@variance@base, 0.95)
expect_identical(sampler$model@variance@n.trees, 7L)
sampler$setModel(sampler$model)

ordinal <- function(k) {
  dbarts(
    x,
    cut(y, k, ordered_result = TRUE),
    family = "ordinal",
    control = migrationControl()
  )
}
structureRefused(ordinal(3L), ordinal(4L)$model, "number of categories")
hazardOf <- function(times) {
  dbarts(
    x,
    cbind(times, rep(1L, n)),
    family = "hazard",
    control = migrationControl()
  )
}
periods <- rep_len(1:3, n)
structureRefused(
  hazardOf(periods),
  hazardOf(rep_len(1:4, n))$model,
  "hazard periods"
)
sampler <- dbarts(y ~ x1 + x2 + x3, frame, control = migrationControl())
twoForest <- dbarts(
  y ~ x1 + x2 + x3,
  frame,
  forests = list(forest(), forest(basis = ~ factor(z))),
  control = migrationControl()
)
structureRefused(sampler, twoForest$model, "number of forests")

# each field of a record, by the check itself: no sampler of several forests
# reaches it yet. Only the fields named refuse: what else a record holds does not
refusedBy <- function(field, value) {
  changed <- twoForest$model
  methods::slot(changed@forests[[2L]], field) <- value
  dbarts:::movedStructureRefusal(changed, twoForest$model)
}
expect_null(dbarts:::movedStructureRefusal(twoForest$model, twoForest$model))
expect_true(grepl("forest 2's tree count", refusedBy("n.trees", 3L)))
expect_true(grepl("forest 2's 'vars'", refusedBy("vars", 1L)))
expect_true(grepl(
  "forest 2's 'interactions'",
  refusedBy("interactions", list(1:2))
))
expect_true(grepl("forest 2's 'blocks'", refusedBy("blocks", list(block = 1L))))
expect_true(grepl(
  "forest 2's 'basis'",
  refusedBy("basis.term", list(label = "b"))
))
expect_true(grepl("forest 2's label", refusedBy("label", "another")))
expect_null(refusedBy("stated", c("n.trees", "sd", "share")))
expect_null(refusedBy("base", 0.5))
expect_null(refusedBy("power", 0.5))
expect_null(refusedBy("sd", 0.5))
expect_null(refusedBy("held.amplitude", TRUE))
expect_null(refusedBy("amplitude.df", 3))
expect_null(refusedBy("moves", list(proposal.probs = 1)))
expect_null(refusedBy("n.trees", NULL))
# a model made by hand states none of it: the sampler's stand
handModel <- varianceSampler$model
handModel@forests <- list()
handModel@variance <- NULL
expect_null(dbarts:::movedStructureRefusal(handModel, varianceSampler$model))
live <- varianceOf(n.trees = 7L)
live$setModel(handModel)
expect_identical(live$model@variance@n.trees, 7L)
expect_identical(live$model@forests[[1L]]@n.trees, 75L)

# a model differing only in a forest's move settings is installed, the draws
# then those of a sampler created with them
moving <- c(birth_death = 0.3, swap = 0.2, change = 0.5, birth = 0.5)
movesOf <- function(moves) {
  dbarts(
    y ~ x1 + x2 + x3,
    frame,
    forests = list(forest(proposal.probs = moves)),
    variance = TRUE,
    control = migrationControl()
  )
}
stating <- movesOf(moving)
expect_false(is.null(stating$model@forests[[1L]]@moves))
installed <- dbarts(
  y ~ x1 + x2 + x3,
  frame,
  variance = TRUE,
  control = migrationControl()
)
withoutRun <- installed$copy()
installed$setModel(stating$model)
expect_identical(
  installed$model@forests[[1L]]@moves,
  stating$model@forests[[1L]]@moves
)
expect_identical(installed$run()$train, stating$run()$train)
expect_false(identical(
  withoutRun$copy()$run()$train,
  stating$copy()$run()$train
))
# what the model's records state is what is in force: none, the control's again
restated <- movesOf(moving)
restated$setModel(withoutRun$model)
expect_null(restated$model@forests[[1L]]@moves)
expect_identical(
  restated$run()$train,
  dbarts(
    y ~ x1 + x2 + x3,
    frame,
    variance = TRUE,
    control = migrationControl()
  )$run()$train
)
rm(structureRefused, sampler, varianceOf, parameter, ordinal, hazardOf, periods)
rm(twoForest, refusedBy, handModel, live, moving, movesOf, stating, installed)
rm(withoutRun, restated)

## ---- getTrees()' column follows the kind of model -------------------------------

columnOf <- function(...) {
  sampler <- dbarts(
    ...,
    control = migrationControl(keepTrees = TRUE)
  )
  sampler$run()
  "forest" %in% colnames(sampler$getTrees())
}
expect_false(columnOf(y ~ x1 + x2 + x3, frame))
expect_false(columnOf(y ~ x1 + x2 + x3, frame, forests = list(forest())))
expect_false(columnOf(y ~ forest(x1 + x2), frame))
expect_true(columnOf(
  y ~ x1 + x2 + x3,
  frame,
  forests = list(forest(), forest(basis = ~ factor(z)))
))
expect_true(columnOf(
  x,
  factor(cut(y, 3L, labels = FALSE)),
  family = "multinomial"
))
rm(columnOf, warned, expect_notUsed, fit50, fit75)
rm(multinomialFit, varianceSampler)
