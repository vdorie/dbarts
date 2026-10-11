# What the response gives a prior - its range, a count fit's center, the unit
# of a forest's default sd, sigest's stand-in and a hazard fit's periods - is
# taken from the rows in the likelihood when a sampler is created: weight
# positive and no mask switching the row off. The fixtures put the rows out
# far from the rows in, and hold no missing predictor value, so an install
# draws no direction and states compare bit for bit.

set.seed(20261011L)
n <- 120L
x <- matrix(runif(n * 2L), n, 2L, dimnames = list(NULL, c("x1", "x2")))
z <- as.double(seq_len(n) %% 2L)
out <- seq_len(n) %% 2L == 0L
rows.in <- !out
w01 <- as.double(rows.in)
y <- as.double(x[, 1L] + 2 * x[, 2L] + rnorm(n, 0, 0.3))
yA <- y
yA[out] <- y[out] + 50
yB <- y
yB[out] <- y[out] - 30
off <- runif(n, 0, 0.5)
df <- function(response) data.frame(y = response, x, z = z)

control <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 10L,
  updateState = FALSE,
  n.samples = 10L,
  seed = 7L
)
sampler <- function(response, ..., ctl = control) {
  dbarts(y ~ x1 + x2, df(response), control = ctl, ...)
}
verboseControl <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 10L,
  updateState = FALSE,
  verbose = TRUE,
  seed = 7L
)
readers <- function(s) s$getLeafPrior()
recordedRange <- function(s) attr(s$model, "response.range")

# Every warning and every message-free condition of a call, by class.
caught <- function(expr) {
  warnings <- list()
  sigmaFallbacks <- 0L
  value <- withCallingHandlers(
    expr,
    warning = function(w) {
      warnings[[length(warnings) + 1L]] <<- w
      invokeRestart("muffleWarning")
    },
    dbartsSigmaFallback = function(condition) {
      sigmaFallbacks <<- sigmaFallbacks + 1L
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
    others = warnings[!scale],
    sigmaFallbacks = sigmaFallbacks
  )
}
noValuesPattern <- "fewer than two distinct values of the response are in the likelihood"

w <- w01

counts <- as.double(rpois(n, 3))
countsA <- counts
countsA[out] <- counts[out] + 30
countsB <- counts
countsB[out] <- 0

# a constant offset keeps the rows in at one value
offKeep <- rep(0.25, n)
keepRows <- list(
  none = rep(0, n),
  one = c(1, rep(0, n - 1L)),
  same = c(rep(1, 4L), rep(0, n - 4L))
)
# the response a case reads: the rows in at one value for "same"
responseFor <- function(case, response, value) {
  if (case == "same") {
    response[1:4] <- value
  }
  response
}
# a sampler created in rule 4's case; its one creation warning is asserted
# where creation is under test
createQuietly <- function(...) {
  made <- caught(sampler(...))
  expect_equal(length(made$scale), 1L)
  expect_equal(length(made$others), 0L)
  made$value
}
rederive <- function(s, via, response, offset, updateScale = TRUE) {
  if (via == "setResponse") {
    s$setResponse(response, updateScale = updateScale)
  } else {
    s$setOffset(offset, updateScale = updateScale)
  }
}
keepText <- paste0(
  noValuesPattern,
  ", so the response's scale and the prior's "
)

creationText <- function(m, n) {
  paste0(
    noValuesPattern,
    " (",
    m,
    " of ",
    n,
    " observations are in it), so the ",
    "response's scale and the prior's center are taken from all ",
    n
  )
}
creationCases <- list(
  none = list(w = rep(0, n), m = 0L, n = n),
  one = list(w = c(1, rep(0, n - 1L)), m = 1L, n = n),
  same = list(w = c(rep(1, 60L), rep(0, n - 60L)), m = 60L, n = n)
)

startingLines <- function(expr) {
  lines <- character()
  utils::capture.output(withCallingHandlers(
    caught(expr),
    message = function(m) {
      lines <<- c(lines, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  ))
  lines
}

three <- c(rep(1, 3L), rep(0, n - 3L))

bcf <- function(...) {
  dbarts(
    y ~ forest(x1 + x2) + forest(x1 + x2, basis = z),
    control = control,
    ...
  )
}

# --- inert: the response of a row out reaches nothing -----------------------
# Two fits that differ only in the response of the rows out at creation are
# identical at the rows in, in sigma and in every reader, and agree at the rows
# out to rounding.
expectInert <- function(build, label) {
  a <- build(yA)
  b <- build(yB)
  expect_identical(readers(a), readers(b), info = label)
  expect_identical(recordedRange(a), recordedRange(b), info = label)
  da <- a$run(10L, 8L)
  db <- b$run(10L, 8L)
  expect_identical(da$train[rows.in, ], db$train[rows.in, ], info = label)
  expect_equal(da$train, db$train, tolerance = 1e-12, info = label)
  if (!is.null(da$sigma)) {
    expect_identical(da$sigma, db$sigma, info = label)
  }
}
expectInert(function(r) sampler(r, weights = w), "gaussian")
expectInert(
  function(r) sampler(r, weights = w, family = gaussian(sigma = fixed(1))),
  "fixed sigma"
)
expectInert(
  function(r) sampler(r, weights = w, family = dbarts:::student(df = 4)),
  "student-t"
)
expectInert(
  function(r) sampler(r, weights = w, variance = varianceForest(n.trees = 10L)),
  "variance forest"
)
expectInert(
  function(r) {
    sampler(r, weights = w, leaf.prior = dbartsPriors$linear("x2"))
  },
  "linear leaves"
)
expectInert(
  function(r) sampler(r, weights = w, leaf.prior = gp("x1", k = 2)),
  "gp leaves"
)
expectInert(
  function(r) sampler(r, weights = w, leaf.prior = normal(sd = 0.65)),
  "named sd"
)
expectInert(function(r) sampler(r, weights = w, offset = off), "offset")
expectInert(
  function(r) {
    dbarts(
      y ~ forest(x1 + x2) + forest(x1 + x2, basis = factor(z)),
      df(r),
      weights = w,
      control = control
    )
  },
  "two forests, a factor basis"
)
expectInert(
  function(r) {
    dbarts(
      y ~ forest(x1 + x2) + forest(x1 + x2, basis = z),
      df(r),
      weights = w,
      control = control
    )
  },
  "two forests, a numeric basis"
)
# the count arm takes its responses from the counts
expectCountInert <- function() {
  build <- function(r) sampler(r, weights = w, family = "nbinom")
  a <- build(countsA)
  b <- build(countsB)
  expect_identical(readers(a), readers(b))
  da <- a$run(10L, 8L)
  db <- b$run(10L, 8L)
  expect_identical(da$train[rows.in, ], db$train[rows.in, ])
}
expectCountInert()

fitBart <- function(r, ...) {
  bart(
    x,
    r,
    weights = w,
    n.trees = 10L,
    n.samples = 8L,
    n.burn = 10L,
    n.chains = 1L,
    n.threads = 1L,
    seed = 7L,
    verbose = FALSE,
    ...
  )
}
a <- fitBart(yA)
b <- fitBart(yB)
expect_identical(a$yhat.train[, rows.in], b$yhat.train[, rows.in])
expect_identical(a$sigma, b$sigma)
fitBartBT <- function(r) {
  bartBT(
    x,
    r,
    weights = w,
    ntree = 10L,
    ndpost = 8L,
    nskip = 10L,
    nchain = 1L,
    nthread = 1L,
    seed = 7L,
    verbose = FALSE
  )
}
a <- fitBartBT(yA)
b <- fitBartBT(yB)
expect_identical(a$yhat.train[, rows.in], b$yhat.train[, rows.in])
expect_identical(a$sigma, b$sigma)
fitXbart <- function(r) {
  xbart(
    x,
    r,
    weights = w,
    n.samples = 10L,
    n.burn = 5L,
    n.trees = 5L,
    n.reps = 2L,
    n.test = 4,
    n.threads = 1L,
    seed = 7L,
    verbose = FALSE
  )
}
expect_identical(fitXbart(yA), fitXbart(yB))

# --- literals, with 0/1 weights ---------------------------------------------
# The range, the center, the leaf scale and sigest are those of the rows in,
# and those of the fit on the subset.
inRange <- range(yA[rows.in])
s <- sampler(yA, weights = w)
sub <- dbarts(
  y ~ x1 + x2,
  df(yA)[rows.in, ],
  control = control
)
expect_equal(recordedRange(s), inRange, tolerance = 1e-14)
lp <- readers(s)
expect_equal(lp$response.shift, mean(inRange), tolerance = 1e-14)
expect_equal(lp$response.scale, diff(inRange), tolerance = 1e-14)
expect_equal(lp$prior.mean, mean(inRange), tolerance = 1e-14)
expect_equal(lp$k.scale, lp$response.scale / 2, tolerance = 1e-14)
expect_equal(lp, readers(sub), tolerance = 1e-14)
expect_equal(
  s$data@sigma,
  summary(lm(y ~ x1 + x2, df(yA), weights = w))$sigma,
  tolerance = 1e-12
)
expect_equal(s$data@sigma, sub$data@sigma, tolerance = 1e-12)
# an sd spelling: k is k.scale over the named sd
sd.s <- sampler(yA, weights = w, leaf.prior = normal(sd = 0.65))
expect_equal(
  readers(sd.s)$k.scale,
  readers(
    dbarts(
      y ~ x1 + x2,
      df(yA)[rows.in, ],
      control = control,
      leaf.prior = normal(sd = 0.65)
    )
  )$k.scale,
  tolerance = 1e-14
)
expect_equal(
  as.numeric(sd.s$getK()),
  readers(sd.s)$k.scale / 0.65,
  tolerance = 1e-12
)
# the unit of two forests' default sd
expect_equal(
  readers(bcf(data = df(yA), weights = w)),
  readers(bcf(data = df(yA)[rows.in, ])),
  tolerance = 1e-14
)
# the center of a count fit
cs <- sampler(countsA, weights = w, family = "nbinom")
expect_equal(
  readers(cs)$response.shift,
  log(sum(countsA[rows.in]) / sum(rows.in)),
  tolerance = 1e-14
)
# fractional weights in: the range, the center and the unit only; sigest is
# the weighted fit
wf <- ifelse(rows.in, 0.5 + runif(n), 0)
fr <- caught(sampler(yA, weights = wf))
expect_equal(recordedRange(fr$value), inRange, tolerance = 1e-14)
expect_equal(
  fr$value$data@sigma,
  summary(lm(y ~ x1 + x2, df(yA), weights = wf))$sigma,
  tolerance = 1e-12
)
expect_equal(length(fr$scale), 0L)
expect_equal(length(fr$others), 1L)
expect_true(grepl(
  "rows with 'weights' of 0 are left out of the likelihood",
  conditionMessage(fr$others[[1L]])
))

# --- held: no later call moves the numbers ----------------------------------
# Read from the engine (getLeafPrior), not from the record alone.
held <- sampler(yA, weights = w)
before <- readers(held)
held$setWeights(w * 0.5)
held$setActiveRows(as.double(seq_len(n) %% 3L != 0L))
expect_identical(readers(held), before)
copied <- held$copy()
expect_identical(readers(copied), before)
expect_equal(
  copied$run(2L, 5L)$train,
  held$run(2L, 5L)$train,
  tolerance = 1e-10
)
held$storeState()
file <- tempfile(fileext = ".rds")
saveRDS(held, file)
expect_identical(readers(readRDS(file)), before)
unlink(file)
rm(held, copied)

# --- updateScale = TRUE reads the rows in -----------------------------------
for (via in c("setResponse", "setOffset")) {
  for (rows in c("weights", "mask")) {
    s <- sampler(yA)
    if (rows == "weights") {
      s$setWeights(w)
    } else {
      s$setActiveRows(w)
    }
    reference <- sampler(yA, weights = w)
    result <- caught(
      if (via == "setResponse") {
        s$setResponse(yA, updateScale = TRUE)
      } else {
        s$setOffset(rep(0, n), updateScale = TRUE)
      }
    )
    expect_equal(length(result$warnings), 0L, info = paste(via, rows))
    expect_equal(
      recordedRange(s),
      range(yA[rows.in]),
      tolerance = 1e-14,
      info = paste(via, rows)
    )
    expect_equal(readers(s), readers(reference), info = paste(via, rows))
    expect_equal(
      recordedRange(s),
      recordedRange(reference),
      info = paste(via, rows)
    )
  }
}
two <- bcf(data = df(yA), weights = w)
expect_error(two$setOffset(rep(0, n), updateScale = TRUE), "forest amplitudes")
expect_error(two$setResponse(yA, updateScale = TRUE), "forest amplitudes")

# --- the keep: too little to read at a re-derivation ------------------------
# The transform and the center in force stay, the call completes as
# updateScale = FALSE does, and the R methods warn once per sampler.
for (family in c("gaussian", "nbinom")) {
  for (case in names(keepRows)) {
    for (via in c("setResponse", "setOffset")) {
      info <- paste(family, case, via)
      wk <- keepRows[[case]]
      yk <- if (family == "gaussian") {
        responseFor(case, yA, 1.5)
      } else {
        responseFor(case, countsA, 2)
      }
      kept <- createQuietly(yk, weights = wk, family = family)
      twin <- createQuietly(yk, weights = wk, family = family)
      before <- readers(kept)
      first <- caught(rederive(kept, via, yk, offKeep))
      second <- caught(rederive(kept, via, yk, offKeep))
      for (call in 1:2) {
        rederive(twin, via, yk, offKeep, updateScale = FALSE)
      }
      expect_equal(length(first$scale), 1L, info = info)
      expect_equal(length(first$others), 0L, info = info)
      expect_true(
        grepl(keepText, conditionMessage(first$scale[[1L]]), fixed = TRUE),
        info = info
      )
      expect_true(
        grepl(
          "kept as they were (said once per sampler)",
          conditionMessage(first$scale[[1L]]),
          fixed = TRUE
        ),
        info = info
      )
      expect_true(inherits(first$scale[[1L]], "dbartsWarning"), info = info)
      expect_equal(length(second$warnings), 0L, info = info)
      expect_identical(readers(kept), before, info = info)
      expect_identical(readers(kept), readers(twin), info = info)
      expect_identical(recordedRange(kept), recordedRange(twin), info = info)
      expect_identical(
        kept$run(3L, 5L)$train,
        twin$run(3L, 5L)$train,
        info = info
      )
      # the other setter does not warn again either
      other <- caught(rederive(
        kept,
        if (via == "setResponse") "setOffset" else "setResponse",
        yk,
        offKeep
      ))
      expect_equal(length(other$warnings), 0L, info = info)
    }
  }
}

# the flag survives setData, which derives as creation does and says so; a
# copy and a reload are new samplers and warn on their own first keep
yk <- responseFor("same", yA, 1.5)
kept <- createQuietly(yk, weights = keepRows$same)
expect_equal(
  length(caught(kept$setOffset(offKeep, updateScale = TRUE))$scale),
  1L
)
replacement <- kept$data
expect_equal(length(caught(kept$setData(replacement))$scale), 1L)
expect_equal(
  length(caught(kept$setOffset(offKeep, updateScale = TRUE))$warnings),
  0L
)
fresh <- createQuietly(yk, weights = keepRows$same)
duplicate <- fresh$copy()
expect_equal(
  length(caught(duplicate$setOffset(offKeep, updateScale = TRUE))$scale),
  1L
)
expect_equal(
  length(caught(duplicate$setOffset(offKeep, updateScale = TRUE))$warnings),
  0L
)
fresh$storeState()
file <- tempfile(fileext = ".rds")
saveRDS(fresh, file)
reloaded <- readRDS(file)
unlink(file)
expect_equal(
  length(caught(reloaded$setResponse(yk, updateScale = TRUE))$scale),
  1L
)
# a constant response with every row in is a window, not a keep
constant <- caught(sampler(rep(3, n)))
expect_equal(length(constant$scale), 0L)
expect_equal(readers(constant$value)$response.shift, 3)
expect_equal(readers(constant$value)$response.scale, 1)
expect_equal(
  length(caught(constant$value$setOffset(offKeep, updateScale = TRUE))$scale),
  0L
)
expect_equal(
  length(
    caught(constant$value$setResponse(rep(3, n), updateScale = TRUE))$scale
  ),
  0L
)

# a warning turned into an error leaves the object and its stored state
# agreeing with the engine
strict <- caught(sampler(
  yk,
  weights = keepRows$same,
  ctl = dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 10L,
    updateState = TRUE,
    n.samples = 10L,
    seed = 7L
  )
))$value
old <- options(warn = 2L)
result <- tryCatch(
  strict$setOffset(offKeep, updateScale = TRUE),
  error = function(e) e
)
options(old)
expect_true(inherits(result, "error"))
expect_identical(strict$data@offset, offKeep)
engineState <- strict$state
strict$storeState()
expect_identical(strict$state, engineState)

# --- rule 4 at creation -----------------------------------------------------
# The numbers are taken over every row, one warning, and the starting sigma is
# the regression on every row with no stand-in announced.
for (case in names(creationCases)) {
  cc <- creationCases[[case]]
  yc <- if (case == "same") responseFor("same", yA, 1.5) else yA
  if (case == "same") {
    yc[1:60] <- 1.5
  }
  made <- caught(sampler(yc, weights = cc$w))
  info <- paste("gaussian", case)
  expect_equal(length(made$scale), 1L, info = info)
  expect_equal(length(made$others), 0L, info = info)
  expect_identical(
    conditionMessage(made$scale[[1L]]),
    creationText(cc$m, cc$n),
    info = info
  )
  unweighted <- sampler(yc)
  expect_equal(readers(made$value), readers(unweighted), info = info)
  expect_equal(made$value$data@sigma, unweighted$data@sigma, info = info)
  expect_equal(recordedRange(made$value), range(yc), info = info)
  expect_equal(made$sigmaFallbacks, 0L, info = info)
  # a count fit under a mask of the same rows
  yk <- if (case == "same") responseFor("same", countsA, 2) else countsA
  if (case == "same") {
    yk[1:60] <- 2
  }
  madeCount <- caught(sampler(yk, weights = cc$w, family = "nbinom"))
  expect_equal(length(madeCount$scale), 1L, info = info)
  expect_identical(
    conditionMessage(madeCount$scale[[1L]]),
    creationText(cc$m, cc$n),
    info = info
  )
  expect_equal(
    readers(madeCount$value),
    readers(sampler(yk, family = "nbinom")),
    info = info
  )
}
# a bart fit records no stand-in, and a verbose fit counts every row
fit <- caught(bart(
  x,
  yA,
  weights = creationCases$none$w,
  n.trees = 5L,
  n.samples = 5L,
  n.burn = 5L,
  n.chains = 1L,
  n.threads = 1L,
  seed = 7L,
  verbose = FALSE
))
expect_equal(length(fit$scale), 1L)
expect_equal(fit$sigmaFallbacks, 0L)
expect_null(fit$value$sigest.fallback)
lines <- startingLines(sampler(
  yA,
  weights = creationCases$none$w,
  ctl = verboseControl
))
estimating <- grep("estimating the starting sigma", lines, value = TRUE)
expect_equal(length(estimating), 1L)
expect_true(grepl(paste0("on ", n, " rows"), estimating, fixed = TRUE))
expect_false(any(grepl("starting sigma is the sd of the response", lines)))
# nothing out, nothing warned: a constant response with every row in, and a
# probit under a mask of zeros
expect_equal(length(caught(sampler(rep(3, n)))$scale), 0L)
expect_equal(
  length(caught(sampler(rep(3, n), family = "nbinom"))$scale),
  0L
)
binary <- as.double(y > median(y))
probit <- caught(sampler(binary, weights = rep(0, n), family = "probit"))
expect_equal(length(probit$scale), 0L)
# fractional weights with a zero among them: the data object's own warning
# as well
fractional <- caught(sampler(yA, weights = c(0.5, rep(0, n - 1L))))
expect_equal(length(fractional$scale), 1L)
expect_equal(length(fractional$others), 1L)
expect_true(grepl(
  "rows with 'weights' of 0 are left out of the likelihood",
  conditionMessage(fractional$others[[1L]]),
  fixed = TRUE
))
# setData to such rows: every row, the same one warning
target <- sampler(y)
replacement <- target$data
replacement@weights <- creationCases$one$w
setDataCaught <- caught(target$setData(replacement))
expect_equal(length(setDataCaught$scale), 1L)
expect_identical(
  conditionMessage(setDataCaught$scale[[1L]]),
  creationText(1L, n)
)
expect_equal(recordedRange(target), range(y))
# dbartsSpec with one row in: the one warning, from the spec
specCaught <- caught(dbartsSpec(
  dbartsData(y ~ x1 + x2, df(yA), weights = creationCases$one$w)
))
expect_equal(length(specCaught$scale), 1L)

# --- sigest -----------------------------------------------------------------
sigestCases <- function(build, label) {
  made <- caught(build())
  expect_equal(length(made$scale), 0L, info = label)
  expect_equal(made$sigmaFallbacks, 1L, info = label)
  expect_equal(
    made$value$data@sigma,
    sd(yA[1:3]),
    tolerance = 1e-12,
    info = label
  )
}
sigestCases(function() sampler(yA, weights = three), "auto")
sigestCases(function() sampler(yA, weights = three, sigest = "dense"), "dense")
sigestCases(
  function() sampler(yA, weights = three, sigest = "sparse"),
  "sparse"
)
if (requireNamespace("Matrix", quietly = TRUE)) {
  sparseX <- Matrix::Matrix(x, sparse = TRUE)
  sigestCases(
    function() dbarts(sparseX, yA, weights = three, control = control),
    "dgCMatrix"
  )
  sigestCases(
    function() {
      dbarts(
        data.frame(f = factor(rep(letters[1:3], length.out = n)), x1 = x[, 1L]),
        yA,
        weights = three,
        control = control,
        sigest = "sparse"
      )
    },
    "indicators"
  )
}
fit <- caught(bart(
  x,
  yA,
  weights = three,
  n.trees = 5L,
  n.samples = 5L,
  n.burn = 5L,
  n.chains = 1L,
  n.threads = 1L,
  seed = 7L,
  verbose = FALSE
))
expect_true(isTRUE(fit$value$sigest.fallback))
lines <- startingLines(sampler(yA, weights = three, ctl = verboseControl))
expect_true(any(grepl("in 3 rows", lines, fixed = TRUE)))
# a response the predictors fit exactly among the rows in: the floor reads
# the rows in
exact <- 1 + 2 * x[, 1L] - x[, 2L]
exact[out] <- 100
floored <- sampler(exact, weights = w)
expect_equal(
  floored$data@sigma,
  sqrt(.Machine$double.eps) * max(1, max(abs(exact[rows.in]))),
  tolerance = 1e-12
)
# rows out with residual degrees of freedom left: the dense design's sigest
# over the rows in, by both routines
for (rule in c("dense", "sparse")) {
  expect_equal(
    sampler(yA, weights = w, sigest = rule)$data@sigma,
    summary(lm(y ~ x1 + x2, df(yA), weights = w))$sigma,
    tolerance = 1e-10,
    info = rule
  )
}
# a fixed residual sd and a supplied number are held, rows out or not
expect_equal(
  sampler(yA, weights = three, sigest = 2.5)$data@sigma,
  2.5
)
fixedSigma <- caught(
  sampler(yA, weights = three, family = gaussian(sigma = fixed(0.7)))
)
expect_equal(fixedSigma$value$data@sigma, sqrt(0.7))
expect_equal(fixedSigma$sigmaFallbacks, 0L)

# --- xbart in rule 4's case --------------------------------------------------
sameIn <- creationCases$same$w
ySame <- yA
ySame[1:60] <- 1.5
xb <- caught(xbart(
  x,
  ySame,
  weights = sameIn,
  n.samples = 10L,
  n.burn = 5L,
  n.trees = 5L,
  n.reps = 2L,
  n.test = 4,
  n.threads = 1L,
  seed = 7L,
  verbose = FALSE
))
expect_equal(length(xb$scale), 1L)
expect_equal(length(xb$others), 0L)
expect_equal(xb$sigmaFallbacks, 0L)
expect_true(all(is.finite(xb$value)))
# each fold's sigma is the unweighted fit over its training rows: the floor is
# handed a finite estimate and every training row's residual
seen <- list()
foldSigma <- function(sigma, residual) {
  seen[[length(seen) + 1L]] <<- c(sigma = sigma, rows = length(residual))
}
trace(
  "floorSigmaEstimate",
  tracer = bquote(.(foldSigma)(sigma, residual)),
  print = FALSE,
  where = asNamespace("dbarts")
)
invisible(caught(xbart(
  x,
  ySame,
  weights = sameIn,
  n.samples = 10L,
  n.burn = 5L,
  n.trees = 5L,
  n.reps = 1L,
  n.test = 4,
  n.threads = 1L,
  seed = 7L,
  verbose = FALSE
)))
untrace("floorSigmaEstimate", where = asNamespace("dbarts"))
seen <- do.call(rbind, seen)
folds <- seen[seen[, "rows"] < n, , drop = FALSE]
expect_true(nrow(folds) >= 4L)
expect_true(all(is.finite(folds[, "sigma"])) && all(folds[, "sigma"] > 0.1))
expect_true(all(folds[, "rows"] == n - n %/% 4L))

# --- creation routes for a count fit ----------------------------------------
rowsInCenter <- log(sum(countsA[rows.in]) / sum(rows.in))
everyRowCenter <- log(sum(countsA) / n)
viaDbarts <- sampler(countsA, weights = w, family = "nbinom")
expect_equal(readers(viaDbarts)$response.shift, rowsInCenter, tolerance = 1e-14)
viaSetData <- sampler(countsA, family = "nbinom")
withWeights <- viaSetData$data
withWeights@weights <- w
viaSetData$setData(withWeights)
expect_equal(
  readers(viaSetData)$response.shift,
  rowsInCenter,
  tolerance = 1e-14
)
spec <- dbartsSpec(
  dbartsData(y ~ x1 + x2, df(countsA), weights = w),
  control = control,
  family = "nbinom"
)
byHand <- new("dbartsSampler", spec$control, spec$model, spec$data)
byHand$setActiveRows(spec$active)
expect_equal(readers(byHand)$response.shift, everyRowCenter, tolerance = 1e-14)

# --- a mask over w against weights w * a ------------------------------------
wf <- 0.5 + runif(n)
a <- as.double(seq_len(n) %% 4L != 0L)
yMask <- y
yMask[a == 0] <- yMask[a == 0] + 40
masked <- sampler(yMask, weights = wf)
masked$setActiveRows(a)
reweighted <- sampler(yMask, weights = wf)
reweighted$setWeights(wf * a)
expect_identical(readers(masked), readers(reweighted))
expect_identical(masked$run(5L, 5L)$train, reweighted$run(5L, 5L)$train)
createdUnder <- caught(sampler(yMask, weights = wf * a))$value
expect_false(identical(readers(masked), readers(createdUnder)))

# --- a state across row sets ------------------------------------------------
# A state installs as stored whatever range its donor was created with.
wDonor <- as.double(seq_len(n) %% 3L != 0L)
yDonor <- y
yDonor[wDonor == 0] <- 90
stateArm <- function(make, compare) {
  donor <- make(y, w01)
  invisible(donor$run(2L, 3L))
  donor$storeState()
  recipient <- make(yDonor, wDonor)
  rangeBefore <- recordedRange(recipient)
  readersBefore <- readers(recipient)
  seen <- caught(recipient$setState(donor$state))
  expect_true(isTRUE(seen$value))
  expect_equal(length(seen$warnings), 0L)
  seenForced <- caught(recipient$setState(donor$state, forceUpdate = TRUE))
  expect_null(seenForced$value)
  expect_equal(length(seenForced$warnings), 0L)
  expect_identical(recordedRange(recipient), rangeBefore)
  expect_identical(readers(recipient), readersBefore)
  recipient$storeState()
  compare(donor$state[[1L]]$forests, recipient$state[[1L]]$forests)
}
sameForests <- function(a, b) expect_identical(a, b)
sameForestsButValues <- function(a, b) {
  for (f in seq_along(a)) {
    expect_identical(
      a[[f]][setdiff(names(a[[f]]), "tree.values")],
      b[[f]][setdiff(names(b[[f]]), "tree.values")]
    )
  }
}
stateArm(
  function(r, ww) sampler(r, weights = ww),
  sameForests
)
stateArm(
  function(r, ww) sampler(r, weights = ww, leaf.prior = gp("x1", k = 2)),
  sameForestsButValues
)
stateArm(
  function(r, ww) {
    dbarts(
      y ~ forest(x1 + x2) + forest(x1 + x2, basis = z),
      df(r),
      weights = ww,
      control = control
    )
  },
  sameForests
)

# --- hazard periods ---------------------------------------------------------
set.seed(20261012L)
nh <- 40L
xh <- matrix(runif(nh * 2L), nh, 2L, dimnames = list(NULL, c("x1", "x2")))
timeH <- rep(1:5, length.out = nh)
statusH <- as.double(runif(nh) < 0.7)
weightH <- rep(1, nh)
weightH[1:2] <- 0
timeH[1:2] <- c(3.5, 20)
statusH[1:2] <- 1
fitHazard <- function(xx, tt, ss, ww, family = "hazard") {
  caught(bart(
    xx,
    cbind(tt, ss),
    weights = ww,
    family = family,
    n.trees = 5L,
    n.samples = 5L,
    n.burn = 5L,
    n.chains = 1L,
    n.threads = 1L,
    seed = 3L,
    verbose = FALSE
  ))
}
without <- fitHazard(xh[-(1:2), ], timeH[-(1:2)], statusH[-(1:2)], NULL)
with <- fitHazard(xh, timeH, statusH, weightH)
expect_equal(length(with$warnings), 0L)
expect_identical(with$value$periods, without$value$periods)
expect_identical(with$value$periods, as.double(1:5))
# the subjects in expand to that fit's rows, covariates, period and response
timesIn <- timeH[weightH > 0]
expanded <- dbarts:::expandDiscreteTimeHazard(
  xh,
  timeH,
  statusH,
  gridTime = timesIn
)
reference <- dbarts:::expandDiscreteTimeHazard(
  xh[-(1:2), ],
  timeH[-(1:2)],
  statusH[-(1:2)]
)
inRows <- expanded$subject > 2L
expect_identical(expanded$x[inRows, ], reference$x)
expect_identical(expanded$y[inRows], reference$y)
# the subject past the last period is at risk in every period with no event
expect_identical(sum(expanded$subject == 2L), 5L)
expect_identical(expanded$y[expanded$subject == 2L], rep(0, 5L))
# a time within the same period of the grid changes nothing
timeH2 <- timeH
timeH2[1L] <- 3.7
moved <- fitHazard(xh, timeH2, statusH, weightH)
expect_identical(moved$value$yhat.train, with$value$yhat.train)
# breaks = K reads the quantiles of the subjects in
quantiled <- fitHazard(
  xh,
  timeH,
  statusH,
  weightH,
  family = dbarts:::hazard(breaks = 3L)
)
expect_identical(
  quantiled$value$periods,
  unique(as.double(quantile(timeH[-(1:2)], (1:3) / 3, names = FALSE)))
)
# stated boundaries are checked against every subject
expect_error(
  fitHazard(
    xh,
    timeH,
    statusH,
    weightH,
    family = dbarts:::hazard(breaks = c(0, 2, 4, 6))
  ),
  "every survival time must lie within"
)
# no subject in: every subject's periods, one warning
none <- fitHazard(xh, timeH, statusH, rep(0, nh))
expect_equal(length(none$scale), 1L)
expect_identical(
  conditionMessage(none$scale[[1L]]),
  paste0(
    "no subject is in the likelihood (0 of ",
    nh,
    " subjects are in it), so ",
    "the hazard periods are taken from all ",
    nh
  )
)
expect_identical(none$value$periods, sort(unique(timeH)))
# the logistic link still refuses a weight of zero
expect_error(
  fitHazard(xh, timeH, statusH, weightH, family = "hazard.logistic"),
  "positive integers"
)
