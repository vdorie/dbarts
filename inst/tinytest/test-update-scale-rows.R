# updateScale = TRUE on setWeights and setActiveRows: the response's scale, or a
# count fit's center, is taken again from the rows the call leaves in the
# likelihood. The oracles are the identity (the one call is the plain setter
# followed by the offset entry called with the offset in force), the literals
# written from the rows in, and the keep (too little to read leaves the scale as
# the plain setter does). The fixtures put the rows out far from the rows in,
# have no missing predictor value and no test set but where one is named.

set.seed(20261013L)
n <- 120L
x <- matrix(runif(n * 2L), n, 2L, dimnames = list(NULL, c("x1", "x2")))
z <- as.double(seq_len(n) %% 2L)
out <- seq_len(n) %% 2L == 0L
rows.in <- !out
w01 <- as.double(rows.in)
# the rows in and twenty of the rows out, so that rows enter
w2 <- w01
w2[which(out)[1:20]] <- 1
y <- as.double(x[, 1L] + 2 * x[, 2L] + rnorm(n, 0, 0.3))
yA <- y
yA[out] <- y[out] + 50
off <- runif(n, 0, 0.5)
counts <- as.double(rpois(n, 3))
countsA <- counts
countsA[out] <- counts[out] + 30
ySurv <- y
ySurv[out] <- y[out] + 20
status <- as.double(seq_len(n) %% 3L != 0L)
offKeep <- rep(0.25, n)

ctl <- dbartsControl(
  n.chains = 2L,
  n.threads = 1L,
  n.trees = 10L,
  updateState = FALSE,
  n.samples = 10L,
  seed = 7L
)
ctlKeep <- dbartsControl(
  n.chains = 2L,
  n.threads = 1L,
  n.trees = 10L,
  updateState = FALSE,
  n.samples = 10L,
  seed = 7L,
  keepTrees = TRUE
)
ctlStore <- dbartsControl(
  n.chains = 2L,
  n.threads = 1L,
  n.trees = 10L,
  updateState = TRUE,
  n.samples = 10L,
  seed = 7L
)
mkG <- function(
  response = yA,
  weights = NULL,
  offset = off,
  control = ctl,
  ...
) {
  dbarts(
    x,
    response,
    weights = weights,
    offset = offset,
    control = control,
    ...
  )
}
mkC <- function(
  response = countsA,
  weights = NULL,
  offset = off,
  control = ctl,
  ...
) {
  dbarts(
    x,
    response,
    weights = weights,
    offset = offset,
    family = "nbinom",
    control = control,
    ...
  )
}
mkS <- function(
  logTime = ySurv,
  offset = off,
  control = ctl
) {
  dbarts(
    x,
    cbind(exp(logTime), status),
    offset = offset,
    family = "aft",
    control = control
  )
}
studentFamily <- dbarts:::student(df = 4)
readers <- function(s) s$getLeafPrior()
recordedRange <- function(s) s$data@response.range
noValuesPattern <- "fewer than two distinct values of the response are in the likelihood"
keepText <- paste0(
  noValuesPattern,
  ", so the response's scale and the prior's center are kept as they were"
)

caught <- function(expr) {
  warnings <- list()
  value <- withCallingHandlers(
    expr,
    warning = function(w) {
      warnings[[length(warnings) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }
  )
  scale <- vapply(
    warnings,
    function(w) inherits(w, "dbartsScaleFallbackWarning"),
    NA
  )
  list(
    value = value,
    warnings = warnings,
    scale = warnings[scale],
    others = warnings[!scale]
  )
}
# a sampler made where the plain creation of a keep case warns once
createQuietly <- function(builder) {
  made <- caught(builder())
  made$value
}
engineState <- function(s) {
  s$storeState()
  s$state
}
# the scale of the transform is kept: nothing the readers or the record carry
expectSame <- function(a, b, info) {
  expect_identical(readers(a), readers(b), info = info)
  expect_identical(a$model, b$model, info = info)
  expect_identical(engineState(a), engineState(b), info = info)
  expect_identical(a$run(0L, 5L), b$run(0L, 5L), info = info)
}

# --- the identity -----------------------------------------------------------
# One call against a twin given the plain setter and then the offset entry
# (the R method setOffset, which with no test set is the bridge's entry and the
# record) with the offset in force.
twinOf <- function(build, call) {
  b <- build()
  call(b, FALSE)
  b$setOffset(b$data@offset, updateScale = TRUE)
  b
}
expectIdentity <- function(build, call, label, moves = TRUE) {
  a <- build()
  before <- readers(a)
  beforeLatents <- a$getLatents()
  made <- caught(call(a, TRUE))
  b <- twinOf(build, call)
  expect_equal(length(made$warnings), 0L, info = label)
  if (moves) {
    expect_false(identical(readers(a), before), info = label)
  }
  expectSame(a, b, label)
  invisible(NULL)
}
setRows <- function(s, via, w, ...) {
  if (via == "setWeights") s$setWeights(w, ...) else s$setActiveRows(w, ...)
}
weightsCall <- function(w) function(s, u) s$setWeights(w, updateScale = u)
maskCall <- function(w) function(s, u) s$setActiveRows(w, updateScale = u)
# a range of rows in is not every row's: start from every row, leave
every <- NULL
gaussianCases <- list(
  gaussian = function() mkG(),
  "fixed sigma" = function() mkG(family = gaussian(sigma = fixed(1))),
  "student-t" = function() mkG(family = studentFamily),
  "variance forest" = function() mkG(variance = varianceForest(n.trees = 10L)),
  "linear leaves" = function() mkG(leaf.prior = dbartsPriors$linear("x2")),
  "gp leaves, no saved draw" = function() mkG(leaf.prior = gp("x1", k = 2)),
  "named sd" = function() mkG(leaf.prior = normal(sd = 0.65)),
  "drawn k" = function() mkG(leaf.prior = normal(k = chi(1.25, 3)))
)
for (label in names(gaussianCases)) {
  expectIdentity(
    gaussianCases[[label]],
    weightsCall(w01),
    paste(label, "weights")
  )
  expectIdentity(gaussianCases[[label]], maskCall(w01), paste(label, "mask"))
}
# rows entering at the call: the install draws their latents at the transform
# in force, and the re-derivation redraws nothing
expectEntering <- function(build, call, label) {
  a <- build()
  latents <- a$getLatents()
  expectIdentity(build, call, label)
  b <- build()
  call(b, FALSE)
  expect_false(identical(latents, b$getLatents()), info = label)
  c2 <- build()
  call(c2, TRUE)
  expect_identical(c2$getLatents(), b$getLatents(), info = label)
}
expectEntering(
  function() mkG(weights = w01, family = studentFamily),
  weightsCall(w2),
  "student-t entering"
)
expectIdentity(
  function() mkG(weights = w01, family = studentFamily),
  weightsCall(w2),
  "student-t entering, identity"
)
countBuild <- function() mkC(weights = w01)
expectEntering(countBuild, maskCall(w2), "count entering by a mask")
expectEntering(countBuild, weightsCall(w2), "count entering by weights")
expectIdentity(function() mkC(), maskCall(w01), "count leaving by a mask")
expectIdentity(function() mkC(), weightsCall(w01), "count leaving by weights")
survBuild <- function() {
  s <- mkS()
  s$setActiveRows(w01, updateScale = TRUE)
  s
}
expectEntering(survBuild, maskCall(w2), "survival entering")
expectIdentity(function() mkS(), maskCall(w01), "survival leaving")
# kept draws: the call rewrites them so they predict what they did
keptBuild <- function() {
  s <- mkG(control = ctlKeep)
  s$run(0L, 5L)
  s
}
expectIdentity(keptBuild, weightsCall(w01), "kept draws, weights")
expectIdentity(keptBuild, maskCall(w01), "kept draws, mask")
kept <- keptBuild()
predicted <- kept$predict(x)
kept$setWeights(w01, updateScale = TRUE)
expect_true(max(abs(kept$predict(x) - predicted)) < 1e-10)
# a gp sampler's saved draws refuse, and a call before one is saved does not
expectIdentity(
  function() mkG(control = ctlKeep, leaf.prior = gp("x1", k = 2)),
  weightsCall(w01),
  "gp leaves, keepTrees on, no draw yet"
)

# --- the test offset --------------------------------------------------------
# The new call leaves data@offset.test and the test fits where the bridge's
# offset entry leaves them; the R method setOffset moves both.
testRows <- c(1L, 3L, 5L, 7L, 9L, 11L, 13L)
mkT <- function() {
  dbarts(
    x,
    yA,
    test = x[testRows, ],
    offset = 0.7,
    control = ctl
  )
}
viaCall <- mkT()
viaBridge <- mkT()
viaMethod <- mkT()
viaCall$setWeights(w01, updateScale = TRUE)
viaBridge$setWeights(w01)
.Call(
  dbarts:::C_dbarts_bartcore_setOffset,
  viaBridge$getPointer(),
  viaBridge$data@offset,
  TRUE
)
dbarts:::recordDerived(viaBridge, viaBridge$getPointer())
viaMethod$setWeights(w01)
viaMethod$setOffset(viaMethod$data@offset, updateScale = TRUE)
expect_equal(viaCall$data@offset.test, rep(0.7, length(testRows)))
expect_identical(viaCall$data@offset.test, viaBridge$data@offset.test)
expect_identical(readers(viaCall), readers(viaBridge))
expect_identical(viaCall$run(0L, 5L)$test, viaBridge$run(0L, 5L)$test)
expect_null(viaMethod$data@offset.test)
expect_false(identical(viaMethod$run(0L, 5L)$test, viaCall$run(0L, 5L)$test))

# --- literals ---------------------------------------------------------------
inRange <- range((yA - off)[rows.in])
expectRange <- function(s, range, info) {
  expect_equal(recordedRange(s), range, tolerance = 1e-14, info = info)
  lp <- readers(s)
  expect_equal(lp$response.shift, mean(range), tolerance = 1e-14, info = info)
  expect_equal(lp$response.scale, diff(range), tolerance = 1e-14, info = info)
  expect_equal(lp$k.scale, diff(range) / 2, tolerance = 1e-14, info = info)
}
for (via in c("setWeights", "setActiveRows")) {
  s <- mkG()
  setRows(s, via, w01, updateScale = TRUE)
  expectRange(s, inRange, via)
  expect_equal(readers(s), readers(mkG(weights = w01)), info = via)
  expect_equal(
    recordedRange(s),
    recordedRange(mkG(weights = w01)),
    info = via
  )
  # back to every row, every row's
  setRows(s, via, rep(1, n), updateScale = TRUE)
  expectRange(s, range(yA - off), paste(via, "every row"))
  expect_equal(readers(s), readers(mkG()), info = via)
}
# a cleared mask re-derives too: every row's range again
cleared <- mkG()
cleared$setActiveRows(w01, updateScale = TRUE)
expectRange(cleared, inRange, "masked before the clear")
cleared$setActiveRows(NULL, updateScale = TRUE)
expectRange(cleared, range(yA - off), "setActiveRows(NULL, updateScale = TRUE)")
expect_null(cleared$activeRows)
# a named sd: k restated against the new scale
named <- mkG(leaf.prior = normal(sd = 0.65))
named$setWeights(w01, updateScale = TRUE)
expect_equal(
  as.numeric(named$getK()),
  rep(readers(named)$k.scale / 0.65, 2L),
  tolerance = 1e-12
)
expect_equal(
  readers(named)$k.scale,
  readers(mkG(weights = w01, leaf.prior = normal(sd = 0.65)))$k.scale,
  tolerance = 1e-14
)
# a student-t sampler, and weights and a mask together: the rows both keep
for (family in list(NULL, studentFamily)) {
  both <- if (is.null(family)) mkG() else mkG(family = family)
  keepBoth <- as.double(rows.in & seq_len(n) > 30L)
  both$setActiveRows(w01)
  both$setWeights(as.double(seq_len(n) > 30L), updateScale = TRUE)
  expectRange(both, range((yA - off)[as.logical(keepBoth)]), "weights and mask")
}
# a count fit's center: the log rate of the rows in (no offset, so a literal)
plain <- mkC(offset = NULL)
plain$setActiveRows(w01, updateScale = TRUE)
expect_equal(
  readers(plain)$response.shift,
  log(mean(countsA[rows.in])),
  tolerance = 1e-14
)
expect_equal(
  readers(plain),
  readers(mkC(offset = NULL, weights = w01)),
  tolerance = 1e-14
)
plain$setWeights(rep(1, n), updateScale = TRUE)
expect_equal(
  readers(plain)$response.shift,
  log(mean(countsA)),
  tolerance = 1e-14
)
# a survival sampler under a mask: the range of the observed log times less
# the offset
sv <- mkS()
sv$setActiveRows(w01, updateScale = TRUE)
subset <- dbarts(
  x[rows.in, ],
  cbind(exp(ySurv[rows.in]), status[rows.in]),
  offset = off[rows.in],
  family = "aft",
  control = ctl
)
expect_equal(
  recordedRange(sv),
  range((ySurv - off)[rows.in]),
  tolerance = 1e-14
)
expect_equal(recordedRange(sv), recordedRange(subset), tolerance = 1e-14)
expect_equal(readers(sv), readers(subset), tolerance = 1e-14)

# --- the argument -----------------------------------------------------------
held <- mkG(weights = w01)
before <- readers(held)
held$setWeights(rep(1, n))
held$setActiveRows(w01)
held$setWeights(w01, updateScale = FALSE)
held$setActiveRows(w2, updateScale = FALSE)
expect_identical(readers(held), before)
# positional: the second argument is updateState, as it was in 0.9-34
positional <- mkG(weights = w01)
before <- readers(positional)
positional$setWeights(rep(1, n), TRUE)
positional$setActiveRows(w2, TRUE)
expect_identical(readers(positional), before)
stored <- positional$state
positional$storeState()
expect_identical(positional$state, stored)
# a value that is not TRUE or FALSE stops with nothing installed
for (bad in list(NA, 1, "TRUE", c(TRUE, FALSE))) {
  s <- mkG()
  twin <- mkG()
  expect_error(
    s$setWeights(w01, updateScale = bad),
    "'updateScale' must be TRUE or FALSE"
  )
  expect_error(
    s$setActiveRows(w01, updateScale = bad),
    "'updateScale' must be TRUE or FALSE"
  )
  expect_null(s$data@weights)
  expect_null(s$activeRows)
  expectSame(s, twin, "bad updateScale")
}

# --- the keep ---------------------------------------------------------------
# From a range that is not every row's, a call whose rows in hold fewer than two
# distinct values installs the weights or the mask and keeps the transform,
# warning once per sampler.
keepRows <- list(
  none = rep(0, n),
  one = c(1, rep(0, n - 1L)),
  same = c(rep(1, 4L), rep(0, n - 4L))
)
keepVia <- function(builder, via) {
  function(s, u, w) {
    if (via == "setWeights") {
      s$setWeights(w, updateScale = u)
    } else {
      s$setActiveRows(w, updateScale = u)
    }
  }
}
expectKeep <- function(builder, via, case, label, count = FALSE) {
  info <- paste(label, via, case)
  wk <- keepRows[[case]]
  call <- keepVia(builder, via)
  s <- createQuietly(function() builder(case))
  twin <- createQuietly(function() builder(case))
  everyRow <- if (case == "same") NULL else builder(case, FALSE)
  before <- readers(s)
  first <- caught(call(s, TRUE, wk))
  second <- caught(call(s, TRUE, wk))
  call(twin, FALSE, wk)
  call(twin, FALSE, wk)
  expect_equal(length(first$scale), 1L, info = info)
  expect_equal(length(first$others), 0L, info = info)
  if (length(first$scale) == 1L) {
    expect_true(
      grepl(keepText, conditionMessage(first$scale[[1L]]), fixed = TRUE),
      info = info
    )
  }
  expect_equal(length(second$warnings), 0L, info = info)
  if (via == "setWeights" && !count) {
    expect_identical(s$data@weights, wk, info = info)
  } else {
    expect_identical(s$activeRows, wk, info = info)
  }
  expect_identical(readers(s), before, info = info)
  expect_identical(recordedRange(s), recordedRange(twin), info = info)
  expectSame(s, twin, info)
  if (!is.null(everyRow)) {
    expect_false(identical(readers(s), readers(everyRow)), info = info)
  }
  # a setter after a keep spends nothing more, in either order
  other <- caught(s$setOffset(s$data@offset, updateScale = TRUE))
  expect_equal(length(other$warnings), 0L, info = info)
}
yKeep <- function(case, response = yA, value = 1.5) {
  if (case == "same") {
    response[1:4] <- value
  }
  response
}
gaussianKeep <- function(family = NULL) {
  function(case, rowsOut = TRUE) {
    args <- list(
      response = yKeep(case),
      offset = offKeep,
      weights = if (rowsOut) w01 else NULL
    )
    if (!is.null(family)) {
      args$family <- family
    }
    do.call(mkG, args)
  }
}
countKeep <- function(case, rowsOut = TRUE) {
  mkC(
    response = yKeep(case, countsA, 2),
    offset = offKeep,
    weights = if (rowsOut) w01 else NULL
  )
}
survKeep <- function(case, rowsOut = TRUE) {
  s <- mkS(logTime = yKeep(case, ySurv, 1.5), offset = offKeep)
  if (rowsOut) {
    s$setActiveRows(w01, updateScale = TRUE)
  }
  s
}
for (case in names(keepRows)) {
  expectKeep(gaussianKeep(), "setWeights", case, "gaussian")
  expectKeep(gaussianKeep(), "setActiveRows", case, "gaussian")
  expectKeep(gaussianKeep(studentFamily), "setWeights", case, "student-t")
  expectKeep(countKeep, "setActiveRows", case, "count", count = TRUE)
  expectKeep(countKeep, "setWeights", case, "count", count = TRUE)
  expectKeep(survKeep, "setActiveRows", case, "survival")
}
# rows back in re-derive and warn nothing
back <- createQuietly(function() gaussianKeep()("none"))
caught(back$setWeights(keepRows$none, updateScale = TRUE))
returned <- caught(back$setWeights(w01, updateScale = TRUE))
expect_equal(length(returned$warnings), 0L)
expect_equal(
  readers(back),
  readers(gaussianKeep()("none")),
  tolerance = 1e-14
)
# the warning is once per sampler: not once per setter, across the four in
# either order, and a copy and a reload are new samplers
crossed <- createQuietly(function() gaussianKeep()("one"))
crossed$setWeights(keepRows$one)
expect_equal(
  length(caught(crossed$setOffset(offKeep, updateScale = TRUE))$scale),
  1L
)
expect_equal(
  length(caught(crossed$setWeights(keepRows$one, updateScale = TRUE))$warnings),
  0L
)
expect_equal(
  length(
    caught(crossed$setActiveRows(keepRows$one, updateScale = TRUE))$warnings
  ),
  0L
)
crossed2 <- createQuietly(function() gaussianKeep()("one"))
crossed2$setActiveRows(keepRows$one)
expect_equal(
  length(
    caught(crossed2$setActiveRows(keepRows$one, updateScale = TRUE))$scale
  ),
  1L
)
expect_equal(
  length(
    caught(crossed2$setResponse(yKeep("one"), updateScale = TRUE))$warnings
  ),
  0L
)
fresh <- createQuietly(function() gaussianKeep()("one"))
duplicate <- fresh$copy()
expect_equal(
  length(caught(duplicate$setWeights(keepRows$one, updateScale = TRUE))$scale),
  1L
)
expect_equal(
  length(
    caught(duplicate$setWeights(keepRows$one, updateScale = TRUE))$warnings
  ),
  0L
)
fresh$storeState()
file <- tempfile(fileext = ".rds")
saveRDS(fresh, file)
reloaded <- readRDS(file)
unlink(file)
expect_equal(
  length(
    caught(reloaded$setActiveRows(keepRows$one, updateScale = TRUE))$scale
  ),
  1L
)
# a warning turned into an error leaves the mirrors, the record and the stored
# state agreeing with the engine: sweeps not stored, then the keep with the
# store asked for, so that the call's own store is what brings it level
for (via in c("setWeights", "setActiveRows")) {
  behind <- createQuietly(function() {
    mkG(response = yKeep("same"), offset = offKeep, weights = w01)
  })
  behind$storeState()
  stored <- behind$state
  recorded <- recordedRange(behind)
  invisible(behind$run(0L, 3L, updateState = FALSE))
  old <- options(warn = 2L)
  result <- tryCatch(
    if (via == "setWeights") {
      behind$setWeights(keepRows$same, updateScale = TRUE, updateState = TRUE)
    } else {
      behind$setActiveRows(
        keepRows$same,
        updateScale = TRUE,
        updateState = TRUE
      )
    },
    error = function(e) e
  )
  options(old)
  expect_true(inherits(result, "error"), info = via)
  if (via == "setWeights") {
    expect_identical(behind$data@weights, keepRows$same, info = via)
  } else {
    expect_identical(behind$activeRows, keepRows$same, info = via)
  }
  expect_identical(recordedRange(behind), recorded, info = via)
  stateNow <- behind$state
  expect_false(identical(stateNow, stored), info = via)
  behind$storeState()
  expect_identical(behind$state, stateNow, info = via)
}

# --- refused, and nothing installed -----------------------------------------
binary <- as.double(y > median(y))
frame <- data.frame(y = yA, x, z = z)
mkTwoForests <- function() {
  dbarts(
    y ~ forest(x1 + x2) + forest(x1 + x2, basis = z),
    frame,
    control = ctl
  )
}
mkTwoProbit <- function() {
  dbarts(
    y ~ forest(x1 + x2) + forest(x1 + x2, basis = z),
    data.frame(y = binary, x, z = z),
    family = "probit",
    control = ctl
  )
}
mkMultinomial <- function() {
  dbarts(
    x,
    factor(rep(c("a", "b", "c"), length.out = n)),
    family = "multinomial",
    control = ctl
  )
}
mkGpSaved <- function() {
  s <- mkG(control = ctlKeep, leaf.prior = gp("x1", k = 2))
  s$run(0L, 3L)
  s
}
amplitudeText <- "does not support a sampler that carries forest amplitudes"
expectRefused <- function(build, call, pattern, label) {
  s <- build()
  twin <- build()
  expect_error(call(s), pattern, info = label)
  expect_identical(s$data@weights, twin$data@weights, info = label)
  expect_identical(s$activeRows, twin$activeRows, info = label)
  expect_identical(readers(s), readers(twin), info = label)
  expect_identical(s$run(0L, 3L), twin$run(0L, 3L), info = label)
}
for (build in list(mkTwoForests, mkTwoProbit)) {
  expectRefused(
    build,
    function(s) s$setWeights(w01, updateScale = TRUE),
    amplitudeText,
    "several forests, weights"
  )
  expectRefused(
    build,
    function(s) s$setActiveRows(w01, updateScale = TRUE),
    amplitudeText,
    "several forests, mask"
  )
}
expectRefused(
  mkMultinomial,
  function(s) s$setActiveRows(w01, updateScale = TRUE),
  "is not available on a multinomial sampler: .*; call it without updateScale$",
  "multinomial, mask"
)
expectRefused(
  mkMultinomial,
  function(s) s$setWeights(w01, updateScale = TRUE),
  "is not available on a multinomial sampler",
  "multinomial, weights"
)
expectRefused(
  mkGpSaved,
  function(s) s$setWeights(w01, updateScale = TRUE),
  "gp leaves that holds saved draws",
  "gp, weights"
)
expectRefused(
  mkGpSaved,
  function(s) s$setActiveRows(w01, updateScale = TRUE),
  "gp leaves that holds saved draws",
  "gp, mask"
)
expectRefused(
  function() mkG(),
  function(s) s$setActiveRows(replace(w01, 1L, 0.5), updateScale = TRUE),
  "'active' must be all 0 or 1",
  "mask holding 0.5"
)
logisticBuild <- function() {
  dbarts(x, binary, weights = rep(2, n), family = "logistic", control = ctl)
}
expectRefused(
  logisticBuild,
  function(s) s$setWeights(replace(rep(2, n), 1L, 0), updateScale = TRUE),
  "logistic weights are observation counts",
  "logistic count of 0"
)
expectRefused(
  function() mkS(),
  function(s) s$setWeights(w01, updateScale = TRUE),
  "do not support case weights",
  "survival weights"
)
# past R's refusal, the bridge's own: several forests refuse the argument
# before anything is installed
for (entry in c("setWeights", "setActiveRows")) {
  s <- mkTwoForests()
  twin <- mkTwoForests()
  symbol <- if (entry == "setWeights") {
    dbarts:::C_dbarts_bartcore_setWeights
  } else {
    dbarts:::C_dbarts_bartcore_setActiveRows
  }
  expect_error(
    .Call(symbol, s$getPointer(), w01, TRUE, FALSE, NULL),
    "updateScale",
    info = entry
  )
  expect_identical(s$run(0L, 3L), twin$run(0L, 3L), info = entry)
}

# --- nothing to read --------------------------------------------------------
# One-forest probit, ordinal, logistic and hazard samplers take the argument
# and ignore it, under a mask that brings rows back.
nothingBuilders <- list(
  probit = list(
    build = function() dbarts(x, binary, family = "probit", control = ctl),
    start = w01,
    call = w2
  ),
  ordinal = list(
    build = function() {
      dbarts(
        x,
        factor(rep(1:3, length.out = n), ordered = TRUE),
        family = "ordinal",
        control = ctl
      )
    },
    start = w01,
    call = w2
  ),
  logistic = list(
    build = function() dbarts(x, binary, family = "logistic", control = ctl),
    start = w01,
    call = w2
  ),
  hazard = list(
    build = function() {
      dbarts(
        x,
        cbind(rep(1:5, length.out = n), as.double(seq_len(n) %% 4L != 0L)),
        family = "hazard",
        control = ctl
      )
    },
    start = NULL,
    call = NULL
  )
)
for (label in names(nothingBuilders)) {
  entry <- nothingBuilders[[label]]
  build <- function() {
    s <- entry$build()
    rows <- length(s$data@y)
    startMask <- if (is.null(entry$start)) {
      as.double(seq_len(rows) %% 3L != 0L)
    } else {
      entry$start
    }
    s$setActiveRows(startMask)
    s
  }
  callMask <- if (is.null(entry$call)) {
    rep(1, length(entry$build()$data@y))
  } else {
    entry$call
  }
  for (via in c("setActiveRows", "setWeights")) {
    if (via == "setWeights" && label %in% c("logistic", "hazard")) {
      next
    }
    a <- build()
    b <- build()
    setRows(a, via, callMask, updateScale = TRUE)
    setRows(b, via, callMask, updateScale = FALSE)
    expectSame(a, b, paste(label, via))
  }
}

# --- by hand ----------------------------------------------------------------
# A count list becomes dbarts()'s sampler by the documented call.
countFrame <- data.frame(y = countsA, x)
countSampler <- function(weights = w01, ...) {
  dbarts(
    y ~ x1 + x2,
    countFrame,
    weights = weights,
    offset = off,
    family = "nbinom",
    control = ctl,
    ...
  )
}
countSpec <- function(weights = w01, response = countsA) {
  dbartsSpec(
    dbartsData(
      y ~ x1 + x2,
      data.frame(y = response, x),
      weights = weights,
      offset = off
    ),
    control = ctl,
    family = "nbinom"
  )
}
byHandOf <- function(spec, ...) {
  s <- new("dbartsSampler", spec$control, spec$model, spec$data)
  s$setActiveRows(spec$active, ...)
  s
}
expectByHand <- function(byHand, reference, info) {
  expect_identical(readers(byHand), readers(reference), info = info)
  expect_identical(byHand$model, reference$model, info = info)
  expect_identical(byHand$data, reference$data, info = info)
  expect_identical(byHand$activeRows, reference$activeRows, info = info)
  expect_identical(engineState(byHand), engineState(reference), info = info)
  expect_identical(byHand$run(0L, 20L), reference$run(0L, 20L), info = info)
}
spec <- countSpec()
expect_false(is.null(spec$active))
expect_equal(
  readers(new(
    "dbartsSampler",
    spec$control,
    spec$model,
    spec$data
  ))$response.shift,
  readers(countSampler(weights = NULL))$response.shift
)
expectByHand(
  byHandOf(spec, updateScale = TRUE),
  countSampler(),
  "count list"
)
# a probit and an ordinal list: the plain install, nothing to re-derive
probitSpec <- dbartsSpec(
  dbartsData(y ~ x1 + x2, data.frame(y = binary, x), weights = w01),
  control = ctl,
  family = "probit"
)
expectByHand(
  byHandOf(probitSpec),
  dbarts(
    y ~ x1 + x2,
    data.frame(y = binary, x),
    weights = w01,
    family = "probit",
    control = ctl
  ),
  "probit list"
)
ordinalFrame <- data.frame(
  y = factor(rep(1:3, length.out = n), ordered = TRUE),
  x
)
expectByHand(
  byHandOf(dbartsSpec(
    dbartsData(y ~ x1 + x2, ordinalFrame, weights = w01),
    control = ctl,
    family = "ordinal"
  )),
  dbarts(
    y ~ x1 + x2,
    ordinalFrame,
    weights = w01,
    family = "ordinal",
    control = ctl
  ),
  "ordinal list"
)
# a probit list with bases = carries amplitudes: the mask handed at creation
# makes it dbarts()'s, a numeric column's default scale read from the rows the
# mask keeps, and the argument is refused with the mask not installed
basesData <- function() {
  dbartsData(
    x,
    binary,
    weights = w01,
    bases = list(NULL, matrix(x[, 1L], n, 1L))
  )
}
basesSpec <- dbartsSpec(basesData(), control = ctl, family = "probit")
expect_false(is.null(basesSpec$active))
expectByHand(
  new(
    "dbartsSampler",
    basesSpec$control,
    basesSpec$model,
    basesSpec$data,
    active = basesSpec$active
  ),
  dbarts(basesData(), family = "probit", control = ctl),
  "probit bases list"
)
refusedHand <- new(
  "dbartsSampler",
  basesSpec$control,
  basesSpec$model,
  basesSpec$data
)
expect_error(
  refusedHand$setActiveRows(basesSpec$active, updateScale = TRUE),
  amplitudeText
)
expect_null(refusedHand$activeRows)
# a gaussian and a student-t list with zero weights are dbarts()'s with no call
for (family in list(NULL, studentFamily)) {
  gFrame <- data.frame(y = yA, x)
  args <- list(
    dbartsData(y ~ x1 + x2, gFrame, weights = w01, offset = off),
    control = ctl
  )
  if (!is.null(family)) {
    args$family <- family
  }
  gSpec <- do.call(dbartsSpec, args)
  expect_null(gSpec$active)
  hand <- new("dbartsSampler", gSpec$control, gSpec$model, gSpec$data)
  expectByHand(
    hand,
    if (is.null(family)) {
      dbarts(y ~ x1 + x2, gFrame, weights = w01, offset = off, control = ctl)
    } else {
      dbarts(
        y ~ x1 + x2,
        gFrame,
        weights = w01,
        offset = off,
        family = family,
        control = ctl
      )
    },
    "gaussian list"
  )
}
# the mask keeps three rows of one count: the creation warning, then the keep's,
# and a later keep silent by hand and audible on dbarts()'s sampler
three <- c(rep(1, 3L), rep(0, n - 3L))
threeCounts <- countsA
threeCounts[1:3] <- 2
specWarned <- caught(countSpec(weights = three, response = threeCounts))
expect_equal(length(specWarned$scale), 1L)
handWarned <- caught(byHandOf(specWarned$value, updateScale = TRUE))
expect_equal(length(handWarned$scale), 1L)
countFrame2 <- data.frame(y = threeCounts, x)
viaDbarts <- caught(
  dbarts(
    y ~ x1 + x2,
    countFrame2,
    weights = three,
    offset = off,
    family = "nbinom",
    control = ctl
  )
)
expect_equal(length(viaDbarts$scale), 1L)
expectByHand(handWarned$value, viaDbarts$value, "three rows of one count")
expect_equal(
  length(
    caught(handWarned$value$setActiveRows(three, updateScale = TRUE))$warnings
  ),
  0L
)
expect_equal(
  length(
    caught(viaDbarts$value$setActiveRows(three, updateScale = TRUE))$scale
  ),
  1L
)
# the loop: re-derived in burn-in, held after, and equal to the identity's calls
masks <- lapply(1:8, function(i) {
  set.seed(i)
  as.double(runif(n) < 0.6)
})
loopCount <- function(derive) {
  s <- countSampler(weights = NULL)
  centers <- numeric(8L)
  for (i in 1:8) {
    derive(s, masks[[i]], i <= 5L)
    centers[i] <- readers(s)$response.shift
    s$run(0L, 1L)
  }
  list(sampler = s, centers = centers)
}
byCall <- loopCount(function(s, m, u) s$setActiveRows(m, updateScale = u))
byIdentity <- loopCount(function(s, m, u) {
  s$setActiveRows(m)
  if (u) {
    s$setOffset(s$data@offset, updateScale = TRUE)
  }
})
expect_identical(byCall$centers, byIdentity$centers)
expect_identical(byCall$centers[6:8], rep(byCall$centers[5L], 3L))
for (i in 1:5) {
  m <- as.logical(masks[[i]])
  expect_equal(
    byCall$centers[i],
    log(sum(countsA[m]) / sum(exp(off[m]))),
    tolerance = 1e-14
  )
}
expect_identical(byCall$sampler$run(0L, 5L), byIdentity$sampler$run(0L, 5L))
loopGaussian <- function(derive) {
  s <- mkG()
  ranges <- numeric(8L)
  for (i in 1:8) {
    derive(s, masks[[i]], i <= 5L)
    ranges[i] <- readers(s)$response.scale
    s$run(0L, 1L)
  }
  list(sampler = s, ranges = ranges)
}
gByCall <- loopGaussian(function(s, m, u) s$setWeights(m, updateScale = u))
gByIdentity <- loopGaussian(function(s, m, u) {
  s$setWeights(m)
  if (u) {
    s$setOffset(s$data@offset, updateScale = TRUE)
  }
})
expect_identical(gByCall$ranges, gByIdentity$ranges)
expect_identical(gByCall$ranges[6:8], rep(gByCall$ranges[5L], 3L))
for (i in 1:5) {
  m <- as.logical(masks[[i]])
  expect_equal(
    gByCall$ranges[i],
    diff(range((yA - off)[m])),
    tolerance = 1e-14
  )
}
expect_identical(gByCall$sampler$run(0L, 5L), gByIdentity$sampler$run(0L, 5L))

# --- the record -------------------------------------------------------------
recorded <- mkG(control = ctlStore)
earlier <- recorded$state
recorded$setWeights(w01, updateScale = TRUE)
now <- readers(recorded)
copied <- recorded$copy()
expect_identical(readers(copied), now)
expect_equal(
  copied$run(0L, 5L)$train,
  recorded$run(0L, 5L)$train,
  tolerance = 1e-10
)
recorded$storeState()
file <- tempfile(fileext = ".rds")
saveRDS(recorded, file)
expect_identical(readers(readRDS(file)), now)
unlink(file)
# a state stored before the re-derivation does not undo it
restored <- mkG()
beforeState <- engineState(restored)
restored$setWeights(w01, updateScale = TRUE)
restoredReaders <- readers(restored)
restored$setState(beforeState)
expect_identical(readers(restored), restoredReaders)
# updateState: TRUE stores a state whose range is the new one, FALSE leaves
# the stored state the object it was
for (flag in c(TRUE, FALSE)) {
  s <- mkG()
  s$storeState()
  stateBefore <- s$state
  s$setActiveRows(w01, updateScale = TRUE, updateState = flag)
  stateAfter <- s$state
  if (flag) {
    expect_false(identical(stateAfter, stateBefore))
    expect_identical(engineState(s), stateAfter)
  } else {
    expect_identical(stateAfter, stateBefore)
    expect_false(identical(engineState(s), stateBefore))
  }
}
