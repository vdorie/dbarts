# A linear or gp leaf's covariate center and scale, and the lengthscale a gp
# leaf derives with them, are read from the rows in the likelihood wherever
# the sampler derives them: rows of positive weight that no mask switches
# off, each counting in full. Where those hold fewer than two distinct values
# of a covariate a derivation (creation, setData, an xbart fold) is refused,
# a gp leaf whose lengthscale is stated taking the placeholder, and a
# re-derivation (setPredictor with updateStandardization = TRUE) keeps that
# covariate's numbers, warning once per sampler. Where no row is in, every row
# is read under the creation warning. Every literal is taken in R over the
# rows named.

set.seed(20261023L, sample.kind = "Rejection")
n <- 200L
inRows <- seq_len(100L)
x1 <- c(runif(100L), runif(100L, 5, 15))
x2 <- runif(n)
y <- 2 * sin(3 * x1) + rnorm(n, 0, 0.3)
w <- rep(c(1, 0), each = 100L)
yBinary <- as.double(y > median(y))
yOrdered <- factor(cut(y, 3L, labels = FALSE), ordered = TRUE)
yCount <- rpois(n, 2)
times <- rep(1:5, length.out = n)
status <- as.double(seq_len(n) %% 4L != 0L)
frame <- data.frame(x1 = x1, x2 = x2)
x <- as.matrix(frame)

controlWith <- function(...) {
  settings <- list(
    n.chains = 2L,
    n.threads = 1L,
    n.trees = 10L,
    n.samples = 5L,
    n.burn = 0L,
    updateState = FALSE,
    seed = 7L
  )
  do.call(dbarts::dbartsControl, utils::modifyList(settings, list(...)))
}
# the response and family each family is created with
families <- list(
  gaussian = list(response = y, family = "gaussian"),
  student = list(response = y, family = "student"),
  probit = list(response = yBinary, family = "probit"),
  ordinal = list(response = yOrdered, family = "ordinal"),
  nbinom = list(response = yCount, family = "nbinom"),
  hazard.probit = list(
    response = cbind(times, status),
    family = "hazard.probit"
  ),
  logistic = list(response = yBinary, family = "logistic"),
  aft = list(response = cbind(exp(y), status), family = "aft")
)
make <- function(
  kind = "gaussian",
  weights = NULL,
  leaf.prior = dbartsPriors$linear("x1"),
  data = x,
  ...
) {
  entry <- families[[kind]]
  dbarts::dbarts(
    data,
    entry$response,
    weights = weights,
    family = entry$family,
    leaf.prior = leaf.prior,
    control = controlWith(...)
  )
}
# the record the sampler writes and the engine's own, unnamed
recordOf <- function(s) lapply(s$data@standardization[[1L]], unname)
readerOf <- function(s) {
  .Call(
    dbarts:::C_dbarts_bartcore_derivedRecords,
    s$getPointer(),
    FALSE
  )$standardization[[1L]]
}
# the median pairwise distance of standardized values, subsampled to at most
# 512 by a stride over the rows read, a missing value counted at 0
heuristic <- function(u) {
  u[is.na(u)] <- 0
  m <- length(u)
  values <- u[seq(1L, m, by = 1L + (m - 1L) %/% 512L)]
  if (length(values) < 2L) {
    return(1)
  }
  distances <- abs(outer(values, values, "-"))
  middle <- stats::median(distances[lower.tri(distances)])
  if (middle > 0) middle else 1
}
# a covariate's center, scale and heuristic lengthscale over the rows read
literal <- function(values, rows = seq_along(values), gp = FALSE) {
  read <- values[rows]
  center <- mean(read, na.rm = TRUE)
  scale <- stats::sd(read, na.rm = TRUE)
  result <- list(center = center, scale = scale)
  if (gp) {
    result$lengthscale <- heuristic((read - center) / scale)
  }
  result
}
expectRecord <- function(s, expected, info) {
  expect_equal(recordOf(s), expected, tolerance = 1e-12, info = info)
  expect_identical(recordOf(s), readerOf(s), info = info)
}
# a call's value and the creation and keep warnings it raised, by text
fallbacks <- function(expr) {
  caught <- character()
  value <- withCallingHandlers(
    expr,
    dbartsScaleFallbackWarning = function(w) {
      caught <<- c(caught, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = caught)
}
sweeps <- function(s, count = 5L) s$run(0L, count)$train

# ---- at creation, each family ----------------------------------------------

# rows out by weight 0 (gaussian, Student-t) or by 0/1 weights the mask takes
# (probit, ordinal, nbinom, hazard.probit): the rows in. logistic weights are
# counts and aft takes none, so every row is in
for (kind in names(families)) {
  rowsOut <- !kind %in% c("logistic", "aft")
  weights <- if (kind == "logistic") {
    rep(2, n)
  } else if (rowsOut) {
    w
  }
  for (gp in c(FALSE, TRUE)) {
    info <- paste(kind, if (gp) "gp" else "linear")
    s <- make(
      kind,
      weights,
      if (gp) dbartsPriors$gp("x1") else dbartsPriors$linear("x1")
    )
    if (kind == "hazard.probit") {
      # the person-period rows of the subjects in, as the expansion lays them
      rows <- s$activeRows != 0
      expected <- literal(s$data@x[, "x1"], rows, gp)
    } else {
      expected <- literal(x1, if (rowsOut) inRows else seq_len(n), gp)
    }
    expectRecord(s, expected, info)
  }
}
# under a family with rows out the every-row numbers are not the rows in's
expect_true(abs(literal(x1)$center - literal(x1, inRows)$center) > 1)

# fractional weights: membership, never size
fractional <- c(rep(c(0.001, 3), 50L), rep(0, 100L))
expect_warning(
  fractionalFit <- make(
    weights = fractional,
    leaf.prior = dbartsPriors$gp("x1")
  ),
  "rows with 'weights' of 0 are left out of the likelihood",
  strict = TRUE
)
expectRecord(
  fractionalFit,
  literal(x1, inRows, gp = TRUE),
  "fractional weights"
)
# the subset twin: rows at weight 0 take the numbers a fit on the rows in has
for (prior in list(dbartsPriors$linear("x1"), dbartsPriors$gp("x1"))) {
  expect_identical(
    recordOf(make(weights = w, leaf.prior = prior)),
    recordOf(dbarts::dbarts(
      x[inRows, ],
      y[inRows],
      leaf.prior = prior,
      control = controlWith()
    ))
  )
}
# every row in and every covariate spread: the every-row numbers, as before
for (kind in names(families)) {
  expectRecord(
    make(
      kind,
      if (kind == "logistic") rep(2, n),
      dbartsPriors$gp(c("x1", "x2"))
    ),
    if (kind == "hazard.probit") {
      s <- make(kind, NULL, dbartsPriors$linear("x1"))
      expanded <- s$data@x
      list(
        center = c(mean(expanded[, "x1"]), mean(expanded[, "x2"])),
        scale = c(stats::sd(expanded[, "x1"]), stats::sd(expanded[, "x2"])),
        lengthscale = c(
          heuristic(scale(expanded[, "x1"])[, 1L]),
          heuristic(scale(expanded[, "x2"])[, 1L])
        )
      )
    } else {
      list(
        center = c(mean(x1), mean(x2)),
        scale = c(stats::sd(x1), stats::sd(x2)),
        lengthscale = c(heuristic(scale(x1)[, 1L]), heuristic(scale(x2)[, 1L]))
      )
    },
    paste("every row in,", kind)
  )
}

# ---- no row in ---------------------------------------------------------------

# every row read, said once: one text naming the response's scale with the
# leaf's under the families that read the response, the leaf's alone under
# the rest
none <- rep(0, n)
joint <- paste0(
  "fewer than two distinct values of the response are in the likelihood ",
  "(0 of 200 observations are in it), so the response's scale, the ",
  "prior's center and the center and scale of \"x1\" are taken from all 200"
)
leafOnly <- paste0(
  "no observation is in the likelihood (0 of 200 have a positive weight), ",
  "so the center and scale of \"x1\" are taken from all 200"
)
for (kind in c("gaussian", "student", "probit", "ordinal", "nbinom")) {
  made <- fallbacks(make(kind, none))
  expect_identical(
    made$warnings,
    if (kind %in% c("probit", "ordinal")) leafOnly else joint,
    info = kind
  )
  expectRecord(made$value, literal(x1), paste("no row in,", kind))
}
# a hazard fit says first that its periods were taken from every subject, a
# count of subjects, and then the leaf's over its person-period rows
made <- fallbacks(make("hazard.probit", none))
periodRows <- length(made$value$data@y)
expect_identical(
  made$warnings,
  c(
    paste0(
      "no subject is in the likelihood (0 of 200 subjects are in it), so the ",
      "hazard periods are taken from all 200"
    ),
    gsub("200", periodRows, leafOnly)
  )
)
expectRecord(made$value, literal(made$value$data@x[, "x1"]), "hazard.probit")
# several covariates named together
made <- fallbacks(make("probit", none, dbartsPriors$linear(c("x1", "x2"))))
expect_identical(
  made$warnings,
  paste0(
    "no observation is in the likelihood (0 of 200 have a positive weight), ",
    "so the centers and scales of \"x1\" and \"x2\" are taken from all 200"
  )
)
# through bart, the one joint text; by hand, the spec's and the leaf's apart
viaBart <- fallbacks(bart(
  x,
  y,
  weights = none,
  leaf.prior = dbartsPriors$linear("x1"),
  n.trees = 5L,
  n.samples = 2L,
  n.burn = 0L,
  n.chains = 1L,
  n.threads = 1L,
  seed = 2L,
  verbose = FALSE
))
expect_identical(viaBart$warnings, joint)
byHand <- fallbacks({
  spec <- dbartsSpec(
    dbartsData(x, y, weights = none),
    control = controlWith(),
    leaf.prior = dbartsPriors$linear("x1")
  )
  new("dbartsSampler", spec$control, spec$model, spec$data)
})
expect_identical(length(byHand$warnings), 2L)
expect_identical(byHand$warnings[2L], joint)
# then a covariate constant over every row is refused
flat <- cbind(x1 = rep(0.5, n), x2 = x2)
expect_error(
  make(data = flat, weights = none),
  paste0(
    "the leaf covariate \"x1\" holds one value over all 200 observations, ",
    "which are read because none is in the likelihood"
  ),
  fixed = TRUE
)

# ---- refused at creation -------------------------------------------------------

oneValueIn <- x
oneValueIn[inRows, "x1"] <- 2
refusalOf <- function(holds, count, gp = FALSE) {
  paste0(
    "the leaf covariate \"x1\" ",
    holds,
    " over the rows in the likelihood (",
    count,
    " of 200 observations are in it), so its center and scale cannot be ",
    "taken from them; ",
    if (gp) {
      "state its lengthscale in gp(), drop it from gp(), "
    } else {
      "drop it from linear(), "
    },
    "or change the rows in (the weights, the 0/1 weights' mask, or subset)"
  )
}
missingIn <- x
missingIn[inRows, "x1"] <- NA
oneIn <- c(1, rep(0, n - 1L))
cases <- list(
  "rows in at one value" = list(
    data = oneValueIn,
    w = w,
    holds = "holds one value",
    count = 100L
  ),
  "constant over every row" = list(
    data = flat,
    w = NULL,
    holds = "holds one value",
    count = 200L
  ),
  "missing on every row in" = list(
    data = missingIn,
    w = w,
    holds = "holds no observed value",
    count = 100L
  ),
  "one row in" = list(
    data = x,
    w = oneIn,
    holds = "holds one value",
    count = 1L
  )
)
for (label in names(cases)) {
  case <- cases[[label]]
  for (gp in c(FALSE, TRUE)) {
    prior <- if (gp) dbartsPriors$gp("x1") else dbartsPriors$linear("x1")
    message <- refusalOf(case$holds, case$count, gp)
    # one row in also leaves the response one value, which is warned first
    expect_error(
      suppressWarnings(make(
        data = case$data,
        weights = case$w,
        leaf.prior = prior
      )),
      message,
      fixed = TRUE,
      info = label
    )
    expect_error(
      suppressWarnings(bart(
        case$data,
        y,
        weights = case$w,
        leaf.prior = prior,
        n.chains = 1L,
        n.threads = 1L,
        verbose = FALSE
      )),
      message,
      fixed = TRUE,
      info = label
    )
    spec <- suppressWarnings(dbartsSpec(
      dbartsData(case$data, yBinary, weights = case$w),
      control = controlWith(),
      leaf.prior = prior,
      family = "probit"
    ))
    expect_error(
      new(
        "dbartsSampler",
        spec$control,
        spec$model,
        spec$data,
        active = spec$active
      ),
      message,
      fixed = TRUE,
      info = label
    )
  }
  # a gp leaf whose lengthscale is stated takes the placeholder: the one
  # value, divided by 1
  stated <- suppressWarnings(make(
    data = case$data,
    weights = case$w,
    leaf.prior = dbartsPriors$gp("x1", lengthscale = 0.7)
  ))
  record <- recordOf(stated)
  expect_identical(names(record), c("center", "scale"), info = label)
  expect_true(all(is.finite(sweeps(stated))), info = label)
}
# its record: the value at its center, NA where no row has spread, and 1
# where rows out read off it
stated <- make(
  data = flat,
  leaf.prior = dbartsPriors$gp("x1", lengthscale = 0.7)
)
expect_identical(recordOf(stated), list(center = 0.5, scale = NA_real_))
stated <- make(
  data = oneValueIn,
  weights = w,
  leaf.prior = dbartsPriors$gp("x1", lengthscale = 0.7)
)
expect_identical(recordOf(stated), list(center = 2, scale = 1))
# rows out at one value, rows in spread: created, over the rows in
outAtOne <- x
outAtOne[-inRows, "x1"] <- 7
expectRecord(
  make(data = outAtOne, weights = w),
  literal(x1, inRows),
  "out at one"
)

# ---- missing values among the rows in -----------------------------------------

# the moments over the observed rows in, the heuristic counting a missing row
# in at 0
withNA <- x
withNA[c(3L, 40L, 150L), "x1"] <- NA
expectRecord(
  make(data = withNA, weights = w, leaf.prior = dbartsPriors$gp("x1")),
  literal(withNA[, "x1"], inRows, gp = TRUE),
  "missing among the rows in"
)

# ---- by hand -------------------------------------------------------------------

# new() handed the spec's mask is dbarts()'s; without it every row is read
for (kind in c("probit", "nbinom")) {
  entry <- families[[kind]]
  viaDbarts <- make(kind, w)
  spec <- dbartsSpec(
    dbartsData(x, entry$response, weights = w),
    control = controlWith(),
    leaf.prior = dbartsPriors$linear("x1"),
    family = entry$family
  )
  withActive <- new(
    "dbartsSampler",
    spec$control,
    spec$model,
    spec$data,
    active = spec$active
  )
  expect_identical(recordOf(withActive), recordOf(viaDbarts), info = kind)
  expect_identical(sweeps(withActive), sweeps(viaDbarts), info = kind)
  without <- new("dbartsSampler", spec$control, spec$model, spec$data)
  expectRecord(without, literal(x1), paste("by hand without active,", kind))
}

# ---- held ------------------------------------------------------------------------

# every call without updateStandardization holds the numbers in force
w2 <- rep(c(0, 1), each = 100L)
a <- as.double(seq_len(n) %% 3L != 0L)
calls <- list(
  "setWeights" = function(s) s$setWeights(w2),
  "setWeights, updateScale" = function(s) s$setWeights(w2, updateScale = TRUE),
  "setActiveRows" = function(s) s$setActiveRows(a),
  "setActiveRows, updateScale" = function(s) {
    s$setActiveRows(a, updateScale = TRUE)
  },
  "setResponse, updateScale" = function(s) {
    s$setResponse(y + 3, updateScale = TRUE)
  },
  "setOffset, updateScale" = function(s) {
    s$setOffset(rep(0.5, n), updateScale = TRUE)
  },
  "setPredictor, matrix" = function(s) {
    s$setPredictor(2 * x, forceUpdate = TRUE)
  },
  "setPredictor, column" = function(s) s$setPredictor(2 * x1, "x1", TRUE),
  "setPredictor, partial" = function(s) {
    s$setPredictor(x1 + 0.01, "x1", "partial")
  }
)
for (label in names(calls)) {
  s <- make(weights = w, leaf.prior = dbartsPriors$linear(c("x1", "x2")))
  before <- recordOf(s)
  invisible(suppressWarnings(calls[[label]](s)))
  expect_identical(recordOf(s), before, info = label)
  expect_identical(readerOf(s), before, info = label)
}

# ---- setPredictor(updateStandardization = TRUE) --------------------------------

liveFits <- function(s) as.vector(s$getFitsWithoutOffset())
x.probe <- cbind(x1 = runif(20L, 0, 15), x2 = runif(20L))
gap <- function(a, b) max(abs(a - b) / pmax(1, abs(a)))
stretched <- 3 * x1 + 1
# after setWeights(w2): over the rows w2 keeps, a second covariate held, the
# kept and live linear fits the functions they were
s <- make(
  weights = w,
  leaf.prior = dbartsPriors$linear(c("x1", "x2")),
  keepTrees = TRUE
)
invisible(s$run(10L, 5L))
s$setWeights(w2)
held <- recordOf(s)
kept <- s$predict(x.probe)
twin <- s$copy()
twin$setPredictor(stretched, "x1", TRUE)
s$setPredictor(stretched, "x1", TRUE, updateStandardization = TRUE)
expected <- literal(stretched, w2 == 1)
record <- recordOf(s)
expect_equal(record$center[[1L]], expected$center, tolerance = 1e-12)
expect_equal(record$scale[[1L]], expected$scale, tolerance = 1e-12)
expect_identical(record$center[[2L]], held$center[[2L]])
expect_identical(record$scale[[2L]], held$scale[[2L]])
expect_identical(recordOf(s), readerOf(s))
expect_true(gap(kept, s$predict(x.probe)) < 1e-12)
expect_true(gap(liveFits(twin), liveFits(s)) < 1e-12)
# after setActiveRows(a), gaussian and probit: the rows a keeps, under
# probit in place of the mask its creation installed
for (kind in c("gaussian", "probit")) {
  s <- make(kind, w)
  s$setActiveRows(a)
  s$setPredictor(stretched, "x1", TRUE, updateStandardization = TRUE)
  rows <- a == 1 & (kind == "probit" | w == 1)
  expectRecord(s, literal(stretched, rows), paste("after setActiveRows,", kind))
}
# a mask then setData to another row count: the new data's rows in, the mask
# in force read and none past it
s <- make("probit", w)
s$setActiveRows(a)
longer <- rbind(x, x[1:40, ])
s$setData(
  dbartsData(longer, c(yBinary, yBinary[1:40])),
  updateStandardization = FALSE
)
longerOther <- 2 * longer[, "x1"]
s$setPredictor(longerOther, "x1", TRUE, updateStandardization = TRUE)
expectRecord(s, literal(longerOther), "after setData to more rows")
s$setActiveRows(as.double(seq_len(240L) %% 2L == 0L))
s$setPredictor(longer[, "x1"], "x1", TRUE, updateStandardization = TRUE)
expectRecord(
  s,
  literal(longer[, "x1"], seq_len(240L) %% 2L == 0L),
  "a mask after setData"
)
# a copy of a sampler holding a mask holds it again
s <- make("probit", w)
s$setActiveRows(a)
copied <- s$copy()
copied$setPredictor(stretched, "x1", TRUE, updateStandardization = TRUE)
expectRecord(copied, literal(stretched, a == 1), "a copy's mask")
# too little keeps: a covariate the rows in hold at one value keeps its
# numbers, once warned per sampler in the leaf's words, the flag the response
# setters' too; a response keep's text is as it was; FALSE holds
keepText <- paste0(
  "fewer than two distinct values of the leaf covariate \"x1\" are in the ",
  "likelihood, so its center and scale are kept as they were (said once ",
  "per sampler)"
)
s <- make(weights = w, leaf.prior = dbartsPriors$linear(c("x1", "x2")))
before <- recordOf(s)
atOne <- ifelse(w == 1, 2, x1)
kept <- fallbacks(
  s$setPredictor(atOne, "x1", TRUE, updateStandardization = TRUE)
)
expect_identical(kept$warnings, keepText)
expect_identical(recordOf(s), before)
expect_identical(readerOf(s), before)
again <- fallbacks(
  s$setPredictor(atOne + 0, "x1", TRUE, updateStandardization = TRUE)
)
expect_identical(again$warnings, character())
# x2 replaced beside it re-derives
both <- cbind(x1 = atOne, x2 = 2 * x2)
s$setPredictor(both, forceUpdate = TRUE, updateStandardization = TRUE)
expect_identical(recordOf(s)$center[[1L]], before$center[[1L]])
expect_equal(
  recordOf(s)$center[[2L]],
  literal(2 * x2, inRows)$center,
  tolerance = 1e-12
)
# no row in: every covariate kept
s <- make(weights = w)
before <- recordOf(s)
s$setWeights(none)
kept <- fallbacks(
  s$setPredictor(stretched, "x1", TRUE, updateStandardization = TRUE)
)
expect_identical(kept$warnings, keepText)
expect_identical(recordOf(s), before)
# the flag is the response setters' too: a response keep first spends it,
# the rows in holding one value of the response less this offset
s <- make(weights = w, leaf.prior = dbartsPriors$linear("x2"))
responseKept <- fallbacks(
  s$setOffset(ifelse(w == 1, y - 1, 0), updateScale = TRUE)
)
expect_identical(
  responseKept$warnings,
  paste0(
    "fewer than two distinct values of the response are in the likelihood, ",
    "so the response's scale and the prior's center are kept as they were ",
    "(said once per sampler)"
  )
)
expect_identical(
  fallbacks(s$setPredictor(
    ifelse(w == 1, 0.5, x2),
    "x2",
    TRUE,
    updateStandardization = TRUE
  ))$warnings,
  character()
)
# several kept at once, named together
s <- make(weights = w, leaf.prior = dbartsPriors$linear(c("x1", "x2")))
expect_identical(
  fallbacks(s$setPredictor(
    cbind(x1 = atOne, x2 = ifelse(w == 1, 0.5, x2)),
    forceUpdate = TRUE,
    updateStandardization = TRUE
  ))$warnings,
  paste0(
    "fewer than two distinct values of each of the leaf covariates \"x1\" ",
    "and \"x2\" are in the likelihood, so their centers and scales are kept ",
    "as they were (said once per sampler)"
  )
)
# FALSE holds whatever the rows in hold
s <- make(weights = w)
before <- recordOf(s)
s$setPredictor(stretched, "x1", TRUE)
expect_identical(recordOf(s), before)

# ---- setData -------------------------------------------------------------------

# a fresh object: its rows in, by weight or under a mask family by the mask
# its 0/1 weights resolve to
s <- make()
s$setData(dbartsData(x, y, weights = w))
expectRecord(s, literal(x1, inRows), "setData, weights")
s <- make("nbinom")
s$setData(dbartsData(x, yCount, weights = w))
expectRecord(s, literal(x1, inRows), "setData, nbinom mask")
# a carried record is held, one without spread too; TRUE derives over the
# rows in; FALSE keeps the sampler's
s <- make()
carried <- dbartsData(x, y, weights = w)
carried@standardization <- list(list(center = c(x1 = 3), scale = c(x1 = NA)))
s$setData(carried)
# the placeholder 1 it is divided by, which these rows read off the center
expect_identical(recordOf(s), list(center = 3, scale = 1))
s$setData(carried, updateStandardization = TRUE)
expectRecord(s, literal(x1, inRows), "setData, TRUE")
s <- make()
before <- recordOf(s)
s$setData(dbartsData(x, y, weights = w), updateStandardization = FALSE)
expect_identical(recordOf(s), before)
# rows in at one value: refused at NULL and at TRUE, the sampler as it was,
# and served at FALSE
setDataRefusal <- paste0(
  "$setData: the leaf covariate \"x1\" holds one value over the rows in the ",
  "likelihood (100 of 200 observations are in it), so its center and scale ",
  "cannot be taken from them; drop it from linear(), change the rows in (the ",
  "weights, the 0/1 weights' mask, or subset), or keep the sampler's with ",
  "updateStandardization = FALSE"
)
s <- make(keepTrees = TRUE)
invisible(s$run(10L, 5L))
twin <- make(keepTrees = TRUE)
invisible(twin$run(10L, 5L))
s$storeState()
state <- s$state
data <- s$data
# another response, of the same length: a refused call must leave the
# response the sampler reads as it was
for (argument in list(NULL, TRUE)) {
  expect_error(
    s$setData(
      dbartsData(oneValueIn, 2 * y, weights = w),
      updateStandardization = argument
    ),
    setDataRefusal,
    fixed = TRUE
  )
  expect_identical(s$data, data)
  s$storeState()
  expect_identical(s$state, state)
}
# the response the engine reads is the old one: a re-derivation of the
# range reads it
s$setOffset(rep(0.1, n), updateScale = TRUE)
twin$setOffset(rep(0.1, n), updateScale = TRUE)
expect_identical(s$data@response.range, twin$data@response.range)
expect_identical(sweeps(s, 3L), sweeps(twin, 3L))
s$setData(
  dbartsData(oneValueIn, y, weights = w),
  updateStandardization = FALSE
)
expect_identical(recordOf(s), recordOf(twin))
# under probit the sampler's own data object carries no weights: setData of
# it drops the mask, and a derivation reads every row
s <- make("probit", w)
s$setData(s$data, updateStandardization = TRUE)
expect_null(s$activeRows)
expectRecord(s, literal(x1), "probit, its own data object")

# ---- the record through copies and reloads ---------------------------------

reloadOf <- function(s) {
  path <- tempfile(fileext = ".rds")
  saveRDS(s, path)
  reloaded <- readRDS(path)
  unlink(path)
  reloaded
}
for (prior in list(dbartsPriors$linear("x1"), dbartsPriors$gp("x1"))) {
  source <- make("probit", w, prior, keepTrees = TRUE, updateState = TRUE)
  invisible(source$run(5L, 5L))
  record <- recordOf(source)
  expect_equal(
    record,
    literal(x1, inRows, is(prior, "dbartsGPPrior")),
    tolerance = 1e-12
  )
  copied <- source$copy()
  routes <- list(
    copy = copied,
    "copy of a copy" = copied$copy(),
    reload = reloadOf(source),
    "reload of a copy" = reloadOf(copied)
  )
  # every route is made before the source draws on
  drawn <- sweeps(source)
  for (label in names(routes)) {
    route <- routes[[label]]
    expect_identical(recordOf(route), record, info = label)
    expect_identical(readerOf(route), record, info = label)
    # the live fits are rebuilt from the trees, to rounding
    expect_equal(sweeps(route), drawn, tolerance = 1e-10, info = label)
  }
  fromData <- new("dbartsSampler", source$control, source$model, source$data)
  expect_identical(recordOf(fromData), record)
  expect_identical(readerOf(fromData), record)
}

# ---- gp with kept draws ----------------------------------------------------

s <- make(weights = w, leaf.prior = dbartsPriors$gp("x1"), keepTrees = TRUE)
invisible(s$run(5L, 5L))
before <- recordOf(s)
expect_error(
  s$setPredictor(stretched, "x1", TRUE, updateStandardization = TRUE),
  "cannot be TRUE for a sampler with gp leaves that holds saved draws"
)
expect_error(
  s$setData(dbartsData(x, y, weights = w), updateStandardization = TRUE),
  "cannot change the covariate standardization"
)
expect_identical(recordOf(s), before)

# ---- xbart -------------------------------------------------------------------

# a fold derives over its training rows of positive weight
spec <- dbartsSpec(
  dbartsData(x, y, weights = w),
  control = controlWith(n.chains = 1L),
  leaf.prior = dbartsPriors$linear("x1")
)
handle <- dbarts:::bartcoreDataHandle(spec$control, spec$data, 1L)
trainRows <- seq(1L, n, by = 2L)
fold <- dbarts:::bartcoreSamplerFromHandle(
  handle,
  spec$control,
  spec$model,
  spec$data,
  trainRows,
  seq(2L, n, by = 2L)
)
foldRecord <- .Call(
  dbarts:::C_dbarts_bartcore_derivedRecords,
  fold$ptr,
  FALSE
)$standardization[[1L]]
expect_equal(
  foldRecord,
  literal(x1[trainRows], w[trainRows] == 1),
  tolerance = 1e-12
)
# a fold whose training rows in hold one value refuses the call
expect_error(
  xbart(
    oneValueIn,
    y,
    weights = w,
    leaf.prior = dbartsPriors$linear("x1"),
    n.reps = 1L,
    n.samples = 2L,
    n.burn = c(2L, 2L),
    n.trees = 5L,
    n.threads = 1L,
    seed = 1L
  ),
  paste0(
    "in a fold of xbart, the leaf covariate \"x1\" holds one value over the ",
    "fold's training rows in the likelihood, so its center and scale cannot ",
    "be taken from them; drop it from linear(), or change the rows (the ",
    "weights, subset or n.test)"
  ),
  fixed = TRUE
)
