# A linear or gp leaf's covariate center and scale, and the lengthscale a gp
# leaf derives with them, are read from the rows in the likelihood wherever
# the sampler derives them: rows of positive weight that no mask switches
# off, each counting in full. Where those hold fewer than two distinct values
# of a covariate a derivation (creation, setData, an xbart fold) is refused,
# a gp leaf whose lengthscale is stated taking the placeholder, and a
# re-derivation (setPredictor with updateStandardization = TRUE) keeps that
# covariate's numbers, warning once per sampler. Where fewer than two rows are
# in, a derivation reads every row instead, under the creation warning. Every
# literal is taken in R over the rows named.

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
# a hazard fit's periods are taken from every subject, counted in subjects,
# and its leaf's numbers from every subject's rows: one text says both
made <- fallbacks(make("hazard.probit", none))
expect_identical(
  made$warnings,
  paste0(
    "no subject is in the likelihood (0 of 200 subjects are in it), so the ",
    "hazard periods and the center and scale of \"x1\" are taken from all 200"
  )
)
expectRecord(made$value, literal(made$value$data@x[, "x1"]), "hazard.probit")
# a gp leaf says so as a linear one does
for (kind in c("gaussian", "probit")) {
  made <- fallbacks(make(kind, none, dbartsPriors$gp("x1")))
  expect_identical(
    made$warnings,
    if (kind == "probit") leafOnly else joint,
    info = paste("gp,", kind)
  )
}
# an unnamed column is the leaf covariate in its column
made <- fallbacks(make(
  data = unname(x),
  weights = none,
  leaf.prior = dbartsPriors$linear(1:2)
))
expect_identical(
  made$warnings,
  paste0(
    "fewer than two distinct values of the response are in the likelihood ",
    "(0 of 200 observations are in it), so the response's scale, the ",
    "prior's center and the centers and scales of the leaf covariate in ",
    "column 1 and the leaf covariate in column 2 are taken from all 200"
  )
)
# a setData that derives over no row in says so, last, as creation does
for (kind in c("gaussian", "probit")) {
  s <- make(kind)
  made <- fallbacks(s$setData(dbartsData(
    x,
    families[[kind]]$response,
    weights = none
  )))
  expect_identical(
    made$warnings,
    if (kind == "probit") leafOnly else joint,
    info = paste("setData,", kind)
  )
}
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
    "read because too few are in the likelihood"
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
  )
)
for (label in names(cases)) {
  case <- cases[[label]]
  for (gp in c(FALSE, TRUE)) {
    prior <- if (gp) dbartsPriors$gp("x1") else dbartsPriors$linear("x1")
    message <- refusalOf(case$holds, case$count, gp)
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
# one row in leaves the response one value, so its scale and the leaf's
# numbers are read from every row, said in one text (dec-B309), through
# dbarts() and bart(), under Student-t and under nbinom's mask, in the words
# it had before fewer than two rows in read every row under every family; a
# family that reads no response reads the leaf from every row in its own
# words (below)
oneRow <- paste0(
  "fewer than two distinct values of the response are in the likelihood ",
  "(1 of 200 observations are in it), so the response's scale, the prior's ",
  "center and the center and scale of \"x1\" are taken from all 200"
)
oneRowLeafOf <- function(count) {
  paste0(
    "one observation is in the likelihood (1 of ",
    count,
    " has a positive weight), so the center and scale of \"x1\" are taken ",
    "from all ",
    count
  )
}
oneRowLeaf <- oneRowLeafOf(200L)
for (gp in c(FALSE, TRUE)) {
  prior <- if (gp) dbartsPriors$gp("x1") else dbartsPriors$linear("x1")
  for (kind in c("gaussian", "student")) {
    made <- fallbacks(make(kind, weights = oneIn, leaf.prior = prior))
    expect_identical(made$warnings, oneRow, info = paste(kind, "gp", gp))
    expectRecord(
      made$value,
      literal(x1, gp = gp),
      paste("one row in,", kind, "gp", gp)
    )
  }
  viaBart <- fallbacks(bart(
    x,
    y,
    weights = oneIn,
    leaf.prior = prior,
    n.trees = 5L,
    n.samples = 2L,
    n.burn = 0L,
    n.chains = 1L,
    n.threads = 1L,
    seed = 2L,
    verbose = FALSE
  ))
  expect_identical(viaBart$warnings, oneRow, info = paste("bart, gp", gp))
  made <- fallbacks(make("probit", oneIn, prior))
  expect_identical(made$warnings, oneRowLeaf, info = paste("probit, gp", gp))
}
made <- fallbacks(make("nbinom", oneIn))
expect_identical(made$warnings, oneRow)
expectRecord(made$value, literal(x1), "one row in, nbinom")
# the value read is the response less its offset under gaussian and the count
# under nbinom: an offset that takes y to one value over the rows in falls
# back, and one beside a single count does not stop it
halfIn <- paste0(
  "fewer than two distinct values of the response are in the likelihood ",
  "(100 of 200 observations are in it), so the response's scale, the prior's ",
  "center and the center and scale of \"x1\" are taken from all 200"
)
made <- fallbacks(dbarts::dbarts(
  x,
  y,
  weights = w,
  offset = ifelse(w == 1, y, 0),
  leaf.prior = dbartsPriors$linear("x1"),
  control = controlWith()
))
expect_identical(made$warnings, halfIn)
expectRecord(made$value, literal(x1), "an offset taking y to one value")
oneCount <- ifelse(w == 1, 2, yCount)
made <- fallbacks(dbarts::dbarts(
  x,
  oneCount,
  weights = w,
  offset = log(1 + x2),
  family = "nbinom",
  leaf.prior = dbartsPriors$linear("x1"),
  control = controlWith()
))
expect_identical(made$warnings, halfIn)
expectRecord(made$value, literal(x1), "one count beside differing offsets")
# the same at a setData that derives the range; one that keeps it reads the
# leaf from every row too, one row being in, in the leaf's words
s <- make()
made <- fallbacks(s$setData(dbartsData(x, y, weights = oneIn)))
expect_identical(made$warnings, oneRow)
expectRecord(s, literal(x1), "setData, one row in")
s <- make()
range <- s$data@response.range
made <- fallbacks(s$setData(
  dbartsData(x, y, weights = oneIn),
  updateScale = FALSE
))
expect_identical(made$warnings, oneRowLeaf)
expectRecord(s, literal(x1), "setData keeping the range, one row in")
expect_identical(s$data@response.range, range)
# a multiplier forest takes no linear leaf, so one row in is refused before
# any warning
expect_identical(
  tryCatch(
    fallbacks(dbarts(
      y ~ forest(x1 + x2) + forest(x1 + x2, basis = x2),
      data.frame(y = y, x),
      weights = oneIn,
      leaf.prior = dbartsPriors$linear("x1"),
      control = controlWith()
    ))$warnings,
    error = function(e) conditionMessage(e)
  ),
  "a treatment forest does not support a linear leaf prior; drop it or fit a single-forest model"
)

# ---- one row in under every family ------------------------------------------

# fewer than two rows in leave no spread to read, so a family that reads no
# response reads the leaf from every row too, under one warning counting the
# rows in, its record flagged and the sampler running: at creation through
# dbarts(), bart() and new() by hand, and at a setData that derives. No fold
# of xbart reaches it under these families: xbart refuses probit's 0/1
# weights and fits neither ordinal nor hazard.probit
allRowsFlag <- function(s) {
  .Call(
    dbarts:::C_dbarts_bartcore_derivedRecords,
    s$getPointer(),
    FALSE
  )$standardization.all.rows
}
for (kind in c("probit", "ordinal", "hazard.probit")) {
  for (gp in c(FALSE, TRUE)) {
    info <- paste("one row in,", kind, if (gp) "gp" else "linear")
    prior <- if (gp) dbartsPriors$gp("x1") else dbartsPriors$linear("x1")
    made <- fallbacks(make(kind, oneIn, prior))
    expect_identical(made$warnings, oneRowLeaf, info = info)
    # a hazard fit's one subject in at time 1 takes one period, one row each
    expect_identical(length(made$value$data@y), n, info = info)
    expectRecord(made$value, literal(made$value$data@x[, "x1"], gp = gp), info)
    expect_identical(allRowsFlag(made$value), TRUE, info = info)
    expect_true(all(is.finite(sweeps(made$value))), info = info)
  }
  # a setData that derives over one row in says so, last, as creation does
  s <- make(kind)
  values <- s$data@x
  made <- fallbacks(s$setData(dbartsData(
    values,
    s$data@y,
    weights = as.double(seq_along(s$data@y) == 1L)
  )))
  expect_identical(
    made$warnings,
    oneRowLeafOf(length(s$data@y)),
    info = paste("setData,", kind)
  )
  expectRecord(s, literal(values[, "x1"]), paste("setData, one row in,", kind))
  expect_identical(allRowsFlag(s), TRUE, info = paste("setData,", kind))
}
# through bart(), and by hand: the spec's mask under probit, and a mask given
# to new() under logistic and aft, which take no 0/1 weights; the aft data
# object drops the record its fit wrote, so that the creation derives
viaBart <- fallbacks(bart(
  x,
  yBinary,
  weights = oneIn,
  leaf.prior = dbartsPriors$linear("x1"),
  n.trees = 5L,
  n.samples = 2L,
  n.burn = 0L,
  n.chains = 1L,
  n.threads = 1L,
  seed = 2L,
  verbose = FALSE
))
expect_identical(viaBart$warnings, oneRowLeaf)
byHand <- fallbacks({
  spec <- dbartsSpec(
    dbartsData(x, yBinary, weights = oneIn),
    control = controlWith(),
    leaf.prior = dbartsPriors$linear("x1"),
    family = "probit"
  )
  new(
    "dbartsSampler",
    spec$control,
    spec$model,
    spec$data,
    active = spec$active
  )
})
expect_identical(byHand$warnings, oneRowLeaf)
expectRecord(byHand$value, literal(x1), "one row in by hand, probit")
spec <- dbartsSpec(
  dbartsData(x, yBinary),
  control = controlWith(),
  leaf.prior = dbartsPriors$gp("x1"),
  family = "logistic"
)
made <- fallbacks(new(
  "dbartsSampler",
  spec$control,
  spec$model,
  spec$data,
  active = oneIn
))
expect_identical(made$warnings, oneRowLeaf)
expectRecord(made$value, literal(x1, gp = TRUE), "one row in, logistic")
aftFit <- make("aft")
aftData <- aftFit$data
aftData@standardization <- NULL
made <- fallbacks(new(
  "dbartsSampler",
  aftFit$control,
  aftFit$model,
  aftData,
  active = oneIn
))
expect_identical(made$warnings, oneRowLeaf)
expectRecord(made$value, literal(x1), "one row in, aft")
# a gp leaf whose lengthscale is stated reads every row as well, where the
# one value took the placeholder before
made <- fallbacks(make(
  "probit",
  oneIn,
  dbartsPriors$gp("x1", lengthscale = 0.7)
))
expect_identical(made$warnings, oneRowLeaf)
expectRecord(made$value, literal(x1), "one row in, stated lengthscale")
# the one row in missing the covariate, and several covariates named together
missingOne <- x
missingOne[1L, "x1"] <- NA
made <- fallbacks(make("probit", oneIn, data = missingOne))
expect_identical(made$warnings, oneRowLeaf)
expectRecord(made$value, literal(missingOne[, "x1"]), "one row in, missing")
made <- fallbacks(make("probit", oneIn, dbartsPriors$linear(c("x1", "x2"))))
expect_identical(
  made$warnings,
  paste0(
    "one observation is in the likelihood (1 of 200 has a positive weight), ",
    "so the centers and scales of \"x1\" and \"x2\" are taken from all 200"
  )
)
# a column constant over every row is refused, read over every row
expect_error(
  make("probit", oneIn, data = flat),
  paste0(
    "the leaf covariate \"x1\" holds one value over all 200 observations, ",
    "read because too few are in the likelihood"
  ),
  fixed = TRUE
)
# two rows in are read, and refused where they hold one value (dec-B465)
twoIn <- c(1, 1, rep(0, n - 2L))
made <- fallbacks(make("probit", twoIn))
expect_identical(made$warnings, character())
expectRecord(made$value, literal(x1, 1:2), "two rows in, probit")
twoAtOne <- x
twoAtOne[1:2, "x1"] <- 2
expect_error(
  make("probit", twoIn, data = twoAtOne),
  refusalOf("holds one value", 2L),
  fixed = TRUE
)
# a sampler of one row, that row in, has every row in: refused as one value
expect_error(
  suppressWarnings(dbarts::dbarts(
    x[1L, , drop = FALSE],
    y[1L],
    weights = 2,
    leaf.prior = dbartsPriors$linear("x1"),
    sigest = 1,
    control = controlWith()
  )),
  "holds one value over the rows in the likelihood (1 of 1 observations",
  fixed = TRUE
)

# rows out at one value, rows in spread: created, over the rows in
outAtOne <- x
outAtOne[-inRows, "x1"] <- 7
expectRecord(
  make(data = outAtOne, weights = w),
  literal(x1, inRows),
  "out at one"
)

# a long column name is named whole, the remedy after it
longName <- strrep("abcdefghij", 60L)
longFlat <- flat
colnames(longFlat)[1L] <- longName
expect_error(
  make(data = longFlat, leaf.prior = dbartsPriors$linear(longName)),
  paste0(
    '"',
    longName,
    '" holds one value over the rows in the likelihood (200 ',
    "of 200 observations are in it), so its center and scale cannot be ",
    "taken from them; drop it from linear(), or change the rows in (the ",
    "weights, the 0/1 weights' mask, or subset)"
  ),
  fixed = TRUE
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
# an unnamed column is the leaf covariate in its column
s <- make(data = unname(x), weights = w, leaf.prior = dbartsPriors$linear(1L))
expect_identical(
  fallbacks(s$setPredictor(
    atOne,
    1L,
    TRUE,
    updateStandardization = TRUE
  ))$warnings,
  paste0(
    "fewer than two distinct values of the leaf covariate in column 1 are in ",
    "the likelihood, so its center and scale are kept as they were (said ",
    "once per sampler)"
  )
)
# a covariate at one value over the rows in that the call does not replace is
# neither kept nor reported: replacing x2 alone warns of nothing
s <- make(weights = w, leaf.prior = dbartsPriors$linear(c("x1", "x2")))
s$setPredictor(atOne, "x1", TRUE)
expect_identical(
  fallbacks(s$setPredictor(
    2 * x2,
    "x2",
    TRUE,
    updateStandardization = TRUE
  ))$warnings,
  character()
)
expect_equal(
  recordOf(s)$center[[2L]],
  literal(2 * x2, inRows)$center,
  tolerance = 1e-12
)
# after the rows in move, a kept covariate leaves one the call does not
# replace held: x2 is not derived over the new rows
s <- make(weights = w, leaf.prior = dbartsPriors$linear(c("x1", "x2")))
before <- recordOf(s)
s$setWeights(w2)
invisible(fallbacks(s$setPredictor(
  ifelse(w2 == 1, 2, x1),
  "x1",
  TRUE,
  updateStandardization = TRUE
)))
expect_identical(recordOf(s), before)
expect_identical(readerOf(s), before)
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
# a carried record that gives x1 alone leaves x2 to derive, and refuses it
# before anything moves where the rows in hold one value of it
twoAt <- x
twoAt[inRows, "x2"] <- 0.25
s <- make(
  weights = w,
  leaf.prior = dbartsPriors$linear(c("x1", "x2")),
  n.burn = 5L
)
twin <- make(
  weights = w,
  leaf.prior = dbartsPriors$linear(c("x1", "x2")),
  n.burn = 5L
)
invisible(s$run())
invisible(twin$run())
partial <- dbartsData(twoAt, y, weights = w)
partial@standardization <- list(list(center = c(x1 = 0.5), scale = c(x1 = 0.3)))
expect_error(
  s$setData(partial),
  'leaf covariate "x2" holds one value',
  fixed = TRUE
)
expect_identical(readerOf(s), readerOf(twin))
expect_identical(sweeps(s, 3L), sweeps(twin, 3L))
# a covariate whose position among the leaf's is not its column is read from
# its column: x2, alone in its leaf, is refused
s <- make(weights = w, leaf.prior = dbartsPriors$linear("x2"))
before <- recordOf(s)
expect_error(
  s$setData(dbartsData(twoAt, y, weights = w)),
  'leaf covariate "x2" holds one value',
  fixed = TRUE
)
expect_identical(recordOf(s), before)
expect_identical(readerOf(s), before)
# a gp leaf whose lengthscale is stated takes the placeholder at a setData
# that derives, at NULL and at TRUE
for (argument in list(NULL, TRUE)) {
  s <- make(
    weights = w,
    leaf.prior = dbartsPriors$gp(c("x1", "x2"), lengthscale = c(1, 1))
  )
  s$setData(dbartsData(twoAt, y, weights = w), updateStandardization = argument)
  record <- recordOf(s)
  expect_identical(record$center[[2L]], 0.25)
  expect_identical(record$scale[[2L]], 1)
  expect_equal(record$center[[1L]], mean(x1[inRows]), tolerance = 1e-12)
}
# a refused setData leaves the mask in force: a later derivation reads the
# rows it keeps
s <- make("probit", w, dbartsPriors$linear(c("x1", "x2")))
expect_error(
  s$setData(dbartsData(twoAt, yBinary, weights = w)),
  "holds one value"
)
s$setPredictor(2 * x1, "x1", TRUE, updateStandardization = TRUE)
expectRecord(
  s,
  list(
    center = c(literal(2 * x1, inRows)$center, literal(x2, inRows)$center),
    scale = c(literal(2 * x1, inRows)$scale, literal(x2, inRows)$scale)
  ),
  "the mask after a refused setData"
)
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
# a fold whose training rows hold one row of positive weight reads its
# response and its leaf over every training row, as a creation does
# (the spec says the response's fallback, once)
specMade <- fallbacks(dbartsSpec(
  dbartsData(x, y, weights = oneIn),
  control = controlWith(n.chains = 1L),
  leaf.prior = dbartsPriors$linear("x1")
))
expect_identical(length(specMade$warnings), 1L)
spec <- specMade$value
fold <- dbarts:::bartcoreSamplerFromHandle(
  dbarts:::bartcoreDataHandle(spec$control, spec$data, 1L),
  spec$control,
  spec$model,
  spec$data,
  trainRows,
  seq(2L, n, by = 2L)
)
foldRecords <- .Call(dbarts:::C_dbarts_bartcore_derivedRecords, fold$ptr, FALSE)
expect_equal(
  foldRecords$standardization[[1L]],
  literal(x1[trainRows]),
  tolerance = 1e-12
)
expect_identical(foldRecords$standardization.all.rows, TRUE)
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

# ---- updateStandardization on setWeights and setActiveRows ---------------------

# rows that move: some out come in, some in go out
wCall <- as.double(seq_len(n) %% 3L != 0L)
leafPriorOf <- function(gp, columns = c("x1", "x2")) {
  if (gp) dbartsPriors$gp(columns) else dbartsPriors$linear(columns)
}
storedState <- function(s) {
  s$storeState()
  s$state
}
# the person-period rows of a hazard fit fill a gp leaf past its cap, which a
# run says; the draws are what is compared
quietSweeps <- function(s) {
  withCallingHandlers(
    sweeps(s),
    dbartsGPFallbackWarning = function(w) invokeRestart("muffleWarning")
  )
}
# a reference class's method by name, which `[[` does not install
setterOf <- function(s, setter) {
  if (setter == "setWeights") s$setWeights else s$setActiveRows
}
# the identity: the call with the argument leaves the sampler where the plain
# call and then setPredictor of every leaf column in force with the argument
# leave a twin, in the record, the stored state and the next sweeps. A gp
# leaf's fits are rebuilt by a predictor update in another order, so under gp
# the sampler is given the twin's plain update first
for (kind in c(
  "gaussian",
  "student",
  "probit",
  "ordinal",
  "nbinom",
  "hazard.probit"
)) {
  for (gp in c(FALSE, TRUE)) {
    info <- paste("identity,", kind, if (gp) "gp" else "linear")
    s <- make(kind, w, leafPriorOf(gp))
    twin <- make(kind, w, leafPriorOf(gp))
    rows <- if (kind == "hazard.probit") {
      as.double(seq_along(s$data@y) %% 3L != 0L)
    } else {
      wCall
    }
    leafValues <- s$data@x[, c("x1", "x2")]
    before <- recordOf(s)
    if (gp) {
      s$setPredictor(leafValues, c("x1", "x2"), TRUE)
    }
    s$setWeights(rows, updateStandardization = TRUE)
    twin$setWeights(rows)
    twin$setPredictor(
      leafValues,
      c("x1", "x2"),
      TRUE,
      updateStandardization = TRUE
    )
    expect_false(identical(recordOf(s), before), info = info)
    expect_identical(
      s$data@standardization,
      twin$data@standardization,
      info = info
    )
    expect_identical(recordOf(s), readerOf(s), info = info)
    expect_identical(storedState(s), storedState(twin), info = info)
    expect_identical(quietSweeps(s), quietSweeps(twin), info = info)
  }
}
# by setActiveRows on a gaussian sampler, over the weights in force and the
# new mask
s <- make(weights = w, leaf.prior = leafPriorOf(FALSE))
twin <- make(weights = w, leaf.prior = leafPriorOf(FALSE))
s$setActiveRows(wCall, updateStandardization = TRUE)
twin$setActiveRows(wCall)
twin$setPredictor(x, c("x1", "x2"), TRUE, updateStandardization = TRUE)
expectRecord(
  s,
  list(
    center = c(mean(x1[w == 1 & wCall == 1]), mean(x2[w == 1 & wCall == 1])),
    scale = c(sd(x1[w == 1 & wCall == 1]), sd(x2[w == 1 & wCall == 1]))
  ),
  "setActiveRows over the weights and the mask"
)
expect_identical(storedState(s), storedState(twin))
expect_identical(sweeps(s), sweeps(twin))
# clearing the mask with the argument re-derives over the rows the weights
# leave in
s <- make(weights = w, leaf.prior = leafPriorOf(FALSE))
twin <- make(weights = w, leaf.prior = leafPriorOf(FALSE))
s$setActiveRows(wCall, updateStandardization = TRUE)
twin$setActiveRows(wCall, updateStandardization = TRUE)
s$setActiveRows(NULL, updateStandardization = TRUE)
twin$setActiveRows(NULL)
twin$setPredictor(x, c("x1", "x2"), TRUE, updateStandardization = TRUE)
expectRecord(
  s,
  list(
    center = c(mean(x1[w == 1]), mean(x2[w == 1])),
    scale = c(sd(x1[w == 1]), sd(x2[w == 1]))
  ),
  "setActiveRows(NULL)"
)
expect_identical(sweeps(s), sweeps(twin))
# kept linear draws predict what they did
s <- make(weights = w, leaf.prior = leafPriorOf(FALSE), keepTrees = TRUE)
invisible(s$run(5L, 5L))
predicted <- s$predict(x.probe)
s$setWeights(wCall, updateStandardization = TRUE)
expectRecord(
  s,
  list(
    center = c(mean(x1[wCall == 1]), mean(x2[wCall == 1])),
    scale = c(sd(x1[wCall == 1]), sd(x2[wCall == 1]))
  ),
  "kept linear draws"
)
expect_true(gap(predicted, s$predict(x.probe)) < 1e-12)

# both arguments: the range and the standardization each over the rows in
s <- make(weights = w, leaf.prior = leafPriorOf(TRUE))
s$setWeights(wCall, updateScale = TRUE, updateStandardization = TRUE)
expect_equal(s$data@response.range, range(y[wCall == 1]), tolerance = 1e-12)
expectRecord(
  s,
  list(
    center = c(mean(x1[wCall == 1]), mean(x2[wCall == 1])),
    scale = c(sd(x1[wCall == 1]), sd(x2[wCall == 1])),
    lengthscale = c(
      literal(x1, wCall == 1, gp = TRUE)$lengthscale,
      literal(x2, wCall == 1, gp = TRUE)$lengthscale
    )
  ),
  "both arguments"
)
# a keep of both: the rows in hold one value of the response, so its scale is
# kept and every covariate with it, x2 though it has spread there, under one
# warning naming the response as the cause
lastRows <- seq_len(n) > 150L
wLast <- as.double(lastRows)
yLast <- replace(y, lastRows, 1)
xLast <- x
xLast[lastRows, "x1"] <- 9
s <- dbarts::dbarts(
  xLast,
  yLast,
  weights = w,
  leaf.prior = leafPriorOf(FALSE),
  control = controlWith()
)
range <- s$data@response.range
before <- recordOf(s)
both <- fallbacks(s$setWeights(
  wLast,
  updateScale = TRUE,
  updateStandardization = TRUE
))
jointText <- paste0(
  "fewer than two distinct values of the response are in the likelihood, so ",
  "the response's scale, the prior's center and the centers and scales of ",
  "\"x1\" and \"x2\" are kept as they were (said once per sampler)"
)
expect_identical(both$warnings, jointText)
expect_identical(s$data@response.range, range)
expect_identical(recordOf(s), before)
expect_identical(readerOf(s), before)
# a second joint keep, as each burn-in pass on the same rows makes, says nothing
expect_identical(
  fallbacks(s$setWeights(
    wLast,
    updateScale = TRUE,
    updateStandardization = TRUE
  ))$warnings,
  character()
)
expect_identical(recordOf(s), before)
# with updateScale FALSE nothing reads the response: each covariate keeps or
# re-derives over the rows in by its own values, x1 kept and x2 derived
s <- dbarts::dbarts(
  xLast,
  yLast,
  weights = w,
  leaf.prior = leafPriorOf(FALSE),
  control = controlWith()
)
before <- recordOf(s)
alone <- fallbacks(s$setWeights(wLast, updateStandardization = TRUE))
expect_identical(alone$warnings, keepText)
expect_identical(recordOf(s)$center[[1L]], before$center[[1L]])
expect_equal(
  recordOf(s)$center[[2L]],
  mean(x2[lastRows]),
  tolerance = 1e-12
)
expect_identical(recordOf(s), readerOf(s))
# the same rows give the same leaf numbers at creation and at a burn-in
# setter: a response at one value over the rows in reads every row at
# creation, under gaussian and nbinom, and a sampler created on every row
# keeps them at the setter, a covariate's spread there notwithstanding
in5 <- seq_len(n) <= 5L
zeros <- as.double(seq_len(n) %in% which(yCount == 0)[1:8])
tied <- list(
  gaussian = list(
    made = function(weights) {
      dbarts::dbarts(
        x,
        replace(y, in5, 1),
        weights = weights,
        leaf.prior = leafPriorOf(FALSE),
        control = controlWith()
      )
    },
    rows = as.double(in5)
  ),
  nbinom = list(
    made = function(weights) make("nbinom", weights, leafPriorOf(FALSE)),
    rows = zeros
  )
)
for (kind in names(tied)) {
  created <- fallbacks(tied[[kind]]$made(tied[[kind]]$rows))
  expect_identical(length(created$warnings), 1L, info = kind)
  expectRecord(
    created$value,
    list(
      center = c(mean(x1), mean(x2)),
      scale = c(sd(x1), sd(x2))
    ),
    paste(kind, "created on the tied rows")
  )
  s <- tied[[kind]]$made(NULL)
  kept <- fallbacks(s$setWeights(
    tied[[kind]]$rows,
    updateScale = TRUE,
    updateStandardization = TRUE
  ))
  expect_identical(kept$warnings, jointText, info = kind)
  expect_identical(recordOf(s), recordOf(created$value), info = kind)
  expect_identical(readerOf(s), recordOf(created$value), info = kind)
}

# the argument: after updateScale, absent and FALSE moving nothing; positional
# setWeights(w, TRUE) still names updateState
for (setter in c("setWeights", "setActiveRows")) {
  s <- make(weights = w, leaf.prior = leafPriorOf(TRUE))
  expect_identical(
    names(formals(setterOf(s, setter))),
    c(
      if (setter == "setWeights") "weights" else "active",
      "updateState",
      "updateScale",
      "updateStandardization"
    ),
    info = setter
  )
  twin <- make(weights = w, leaf.prior = leafPriorOf(TRUE))
  before <- recordOf(s)
  setterOf(s, setter)(wCall, updateStandardization = FALSE)
  setterOf(twin, setter)(wCall)
  expect_identical(recordOf(s), before, info = setter)
  expect_identical(recordOf(twin), before, info = setter)
  expect_identical(sweeps(s), sweeps(twin), info = setter)
  # NA and 1 are refused, nothing installed
  for (value in list(NA, 1)) {
    weights <- s$data@weights
    active <- s$activeRows
    expect_error(
      setterOf(s, setter)(a, updateStandardization = value),
      "'updateStandardization' must be TRUE or FALSE",
      fixed = TRUE,
      info = setter
    )
    expect_identical(s$data@weights, weights, info = setter)
    expect_identical(s$activeRows, active, info = setter)
  }
}
s <- make(weights = w, leaf.prior = leafPriorOf(FALSE), updateState = FALSE)
before <- recordOf(s)
state <- storedState(s)
invisible(s$run(0L, 1L))
s$setWeights(wCall, TRUE)
expect_identical(recordOf(s), before)
expect_false(identical(s$state, state))
# on a sampler of several forests and one with no leaf covariate the
# argument is taken and has nothing to do
noLeaf <- list(
  "several forests" = function() {
    dbarts::dbarts(
      y ~ forest(x1 + x2) + forest(x1 + x2, basis = b),
      data.frame(y = y, x, b = as.double(seq_len(n) %% 2L)),
      control = controlWith()
    )
  },
  "constant leaves" = function() {
    dbarts::dbarts(x, y, weights = w, control = controlWith())
  }
)
for (label in names(noLeaf)) {
  for (setter in c("setWeights", "setActiveRows")) {
    info <- paste(label, setter)
    s <- noLeaf[[label]]()
    twin <- noLeaf[[label]]()
    setterOf(s, setter)(wCall, updateStandardization = TRUE)
    setterOf(twin, setter)(wCall)
    expect_null(s$data@standardization, info = info)
    expect_identical(storedState(s), storedState(twin), info = info)
    expect_identical(sweeps(s), sweeps(twin), info = info)
  }
}

# the keep: the rows in hold one value of x1 and spread of x2, so x1 keeps
# its numbers and x2 re-derives, with one warning in the leaf's words; none on
# a second keep, by either setter
s <- dbarts::dbarts(
  xLast,
  y,
  weights = w,
  leaf.prior = leafPriorOf(FALSE),
  control = controlWith()
)
before <- recordOf(s)
kept <- fallbacks(s$setWeights(wLast, updateStandardization = TRUE))
expect_identical(kept$warnings, keepText)
expect_identical(recordOf(s)$center[[1L]], before$center[[1L]])
expect_identical(recordOf(s)$scale[[1L]], before$scale[[1L]])
expect_equal(
  recordOf(s)[c("center", "scale")],
  list(
    center = c(before$center[[1L]], mean(x2[lastRows])),
    scale = c(before$scale[[1L]], sd(x2[lastRows]))
  ),
  tolerance = 1e-12
)
expect_identical(recordOf(s), readerOf(s))
again <- fallbacks(s$setActiveRows(
  as.double(seq_len(n) > 160L),
  updateStandardization = TRUE
))
expect_identical(again$warnings, character())
# the flag is the response's keeps' too
s <- dbarts::dbarts(
  xLast,
  y,
  weights = w,
  leaf.prior = leafPriorOf(FALSE),
  control = controlWith()
)
invisible(fallbacks(s$setOffset(ifelse(w == 1, y - 1, 0), updateScale = TRUE)))
expect_identical(
  fallbacks(s$setWeights(wLast, updateStandardization = TRUE))$warnings,
  character()
)
# by setActiveRows on a mask family, in the same words
s <- dbarts::dbarts(
  xLast,
  yBinary,
  weights = w,
  family = "probit",
  leaf.prior = leafPriorOf(FALSE),
  control = controlWith()
)
expect_identical(
  fallbacks(s$setActiveRows(wLast, updateStandardization = TRUE))$warnings,
  keepText
)
# one row in: a re-derivation keeps, where creation reads every row under
# every family; with both arguments the response's scale and every
# covariate are kept under one warning, and a family that reads no response
# keeps every covariate in the leaf's words
s <- make(weights = w, leaf.prior = leafPriorOf(TRUE))
range <- s$data@response.range
before <- recordOf(s)
expect_identical(
  fallbacks(s$setWeights(
    oneIn,
    updateScale = TRUE,
    updateStandardization = TRUE
  ))$warnings,
  jointText
)
expect_identical(s$data@response.range, range)
expect_identical(recordOf(s), before)
expect_identical(readerOf(s), before)
s <- make("probit", w, leafPriorOf(FALSE))
before <- recordOf(s)
expect_identical(
  fallbacks(s$setActiveRows(oneIn, updateStandardization = TRUE))$warnings,
  paste0(
    "fewer than two distinct values of each of the leaf covariates \"x1\" ",
    "and \"x2\" are in the likelihood, so their centers and scales are kept ",
    "as they were (said once per sampler)"
  )
)
expect_identical(recordOf(s), before)

# gp with a kept draw: refused before the install, the weights, the mask and
# the engine as they were
for (kind in c("gaussian", "probit")) {
  for (setter in c("setWeights", "setActiveRows")) {
    info <- paste(kind, setter)
    s <- make(kind, w, dbartsPriors$gp("x1"), keepTrees = TRUE)
    twin <- make(kind, w, dbartsPriors$gp("x1"), keepTrees = TRUE)
    invisible(s$run(5L, 5L))
    invisible(twin$run(5L, 5L))
    weights <- s$data@weights
    active <- s$activeRows
    before <- recordOf(s)
    expect_error(
      setterOf(s, setter)(wCall, updateStandardization = TRUE),
      paste0(
        if (setter == "setWeights" && kind == "gaussian") {
          "$setWeights"
        } else {
          "$setActiveRows"
        },
        ": 'updateStandardization' cannot be TRUE for a sampler with gp ",
        "leaves that holds saved draws"
      ),
      fixed = TRUE,
      info = info
    )
    expect_identical(s$data@weights, weights, info = info)
    expect_identical(s$activeRows, active, info = info)
    expect_identical(recordOf(s), before, info = info)
    expect_identical(sweeps(s, 3L), sweeps(twin, 3L), info = info)
  }
}

# the burn-in loop of the help: during burn-in the record follows each pass's
# rows in, after it holds the last burn-in pass's, every kept draw replays
# under it, a twin written with the identity's calls is the same sampler, and
# a copy and a reload hold the record
set.seed(20261024L, sample.kind = "Rejection")
numBurn <- 3L
numKeep <- 3L
passes <- lapply(seq_len(numBurn + numKeep), function(i) {
  as.double(runif(n) < 0.6)
})
loopControl <- controlWith(keepTrees = TRUE, n.samples = numKeep)
s <- dbarts::dbarts(
  x,
  y,
  weights = w,
  leaf.prior = dbartsPriors$linear("x1"),
  control = loopControl
)
twin <- dbarts::dbarts(
  x,
  y,
  weights = w,
  leaf.prior = dbartsPriors$linear("x1"),
  control = loopControl
)
trained <- list()
for (i in seq_along(passes)) {
  burning <- i <= numBurn
  s$setWeights(
    passes[[i]],
    updateScale = burning,
    updateStandardization = burning
  )
  twin$setWeights(passes[[i]], updateScale = burning)
  if (burning) {
    twin$setPredictor(x1, "x1", TRUE, updateStandardization = TRUE)
  }
  read <- passes[[min(i, numBurn)]] == 1
  expectRecord(s, literal(x1, read), paste("setWeights loop, pass", i))
  expect_identical(s$data, twin$data, info = paste("pass", i))
  draw <- if (burning) s$run(1L, 0L) else s$run(0L, 1L)
  twinDraw <- if (burning) twin$run(1L, 0L) else twin$run(0L, 1L)
  expect_identical(draw, twinDraw, info = paste("pass", i))
  if (!burning) {
    trained[[length(trained) + 1L]] <- draw$train
  }
}
predicted <- s$predict(x)
for (k in seq_len(numKeep)) {
  expect_equal(
    predicted[, k, ],
    trained[[k]][, 1L, ],
    tolerance = 1e-12,
    info = paste("kept draw", k)
  )
}
record <- recordOf(s)
s$storeState()
for (route in list(s$copy(), reloadOf(s))) {
  expect_identical(recordOf(route), record)
  expect_identical(readerOf(route), record)
  expect_equal(route$predict(x), predicted, tolerance = 1e-12)
}
# by setActiveRows on a probit sampler with a gp leaf, the masks after
# burn-in moving rows the record does not follow: a copy and a reload are
# made on the record, and a saved gp draw replays only under the
# standardization it was drawn with
s <- make(
  "probit",
  w,
  dbartsPriors$gp("x1"),
  keepTrees = TRUE,
  n.samples = numKeep
)
for (i in seq_along(passes)) {
  burning <- i <= numBurn
  s$setActiveRows(
    passes[[i]],
    updateScale = burning,
    updateStandardization = burning
  )
  read <- passes[[min(i, numBurn)]] == 1
  expectRecord(
    s,
    literal(x1, read, gp = TRUE),
    paste("setActiveRows loop, pass", i)
  )
  invisible(if (burning) s$run(1L, 0L) else s$run(0L, 1L))
}
record <- recordOf(s)
predicted <- s$predict(x)
s$storeState()
for (route in list(s$copy(), reloadOf(s))) {
  expect_identical(route$activeRows, passes[[length(passes)]])
  expect_identical(recordOf(route), record)
  expect_identical(readerOf(route), record)
  expect_equal(route$predict(x), predicted, tolerance = 1e-12)
}
