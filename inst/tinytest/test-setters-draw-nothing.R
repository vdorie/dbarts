# No setter draws. Each stores what it is given and rebuilds what follows from
# the latents in force; the next sweep draws the augmentation variables first,
# against the change. So between a setter and the next run, getLatents reads
# the values from before the call (an aft setResponse puts every row at its
# new time, a censored row at its new bound), no chain generator moves, and the
# drawn family parameters - the ordinal cutpoints, the Student-t degrees of
# freedom, the negative-binomial shape - are kept, by setData too.

nbinom <- dbartsFamilies$nbinom
student <- dbartsFamilies$student

set.seed(4401L, sample.kind = "Rejection")
n <- 90L
x <- matrix(runif(n * 2L), n, 2L)
xNew <- x
xNew[, 2L] <- runif(n)
f <- 2 * x[, 1L] - 1
yBinary <- as.double(rbinom(n, 1L, pnorm(f)))
yOrdinal <- as.double(cut(f + rnorm(n), c(-Inf, -0.3, 0.4, Inf)))
yCount <- as.double(rpois(n, exp(1 + x[, 1L])))
yHeavy <- f + 0.5 * rt(n, 4)
counts <- as.double(sample(1:3, n, replace = TRUE))
positive <- 0.5 + runif(n)
logTime <- f + rnorm(n, sd = 0.5)
bound <- as.double(quantile(logTime, 0.75))
status <- as.double(logTime <= bound)
logTime <- pmin(logTime, bound)
mask <- as.double(seq_len(n) %% 4L != 1L)
offset <- 0.2 * x[, 2L]

control <- dbartsControl(
  n.chains = 2L,
  n.threads = 1L,
  updateState = FALSE,
  seed = 31L
)
build <- function(y, ...) dbarts(x, y, ..., control = control, n.trees = 15L)
streamOf <- function(sampler) {
  sampler$storeState()
  lapply(sampler$state, function(chain) chain[["rng.state"]])
}
stateOf <- function(sampler, name) {
  sampler$storeState()
  lapply(sampler$state, function(chain) chain[[name]])
}
# a sampler some sweeps in, the call, then what moved
expectDrawsNothing <- function(make, call, info, latents = NULL) {
  sampler <- make()
  invisible(sampler$run(5L, 2L))
  held <- sampler$getLatents()
  stream <- streamOf(sampler)
  call(sampler)
  expect_identical(
    sampler$getLatents(),
    if (is.null(latents)) held else latents,
    info = info
  )
  expect_identical(streamOf(sampler), stream, info = info)
  expect_true(all(is.finite(sampler$run(0L, 2L)$train)), info = info)
}

families <- list(
  probit = list(
    make = function() build(yBinary, family = "probit"),
    y = 1 - yBinary
  ),
  ordinal = list(
    make = function() build(yOrdinal, family = "ordinal"),
    y = 4 - yOrdinal
  ),
  logistic = list(
    make = function() build(yBinary, family = "logistic", weights = counts),
    y = 1 - yBinary,
    weights = rev(counts)
  ),
  nbinom = list(
    make = function() build(yCount, family = nbinom()),
    y = rev(yCount)
  ),
  student = list(
    make = function() build(yHeavy, family = student(), weights = positive),
    y = rev(yHeavy),
    weights = rev(positive),
    sigma = TRUE
  ),
  aft = list(
    make = function() build(cbind(exp(logTime), status), family = "aft"),
    y = NULL
  )
)
for (name in names(families)) {
  arm <- families[[name]]
  if (!is.null(arm$y)) {
    expectDrawsNothing(
      arm$make,
      function(s) s$setResponse(arm$y),
      paste(name, "setResponse")
    )
  }
  expectDrawsNothing(
    arm$make,
    function(s) s$setOffset(offset),
    paste(name, "setOffset")
  )
  expectDrawsNothing(
    arm$make,
    function(s) s$setPredictor(xNew[, 2L], 2L),
    paste(name, "setPredictor")
  )
  expectDrawsNothing(
    arm$make,
    function(s) s$setActiveRows(mask),
    paste(name, "setActiveRows")
  )
  expectDrawsNothing(
    function() {
      s <- arm$make()
      s$setActiveRows(mask)
      s
    },
    function(s) s$setActiveRows(NULL),
    paste(name, "setActiveRows lifted")
  )
  if (!is.null(arm$weights)) {
    expectDrawsNothing(
      arm$make,
      function(s) s$setWeights(arm$weights),
      paste(name, "setWeights")
    )
  }
  if (isTRUE(arm$sigma)) {
    expectDrawsNothing(
      arm$make,
      function(s) s$setSigma(0.7),
      paste(name, "setSigma")
    )
  }
}

# aft setResponse: every row at its new log-time, a censored row at its new
# bound, the generators where they were
newLogTime <- logTime + 0.5
expectDrawsNothing(
  families$aft$make,
  function(s) s$setResponse(newLogTime),
  "aft setResponse",
  latents = matrix(newLogTime, n, 2L)
)

# a new negative-binomial shape keeps omega
expectDrawsNothing(
  function() build(yCount, family = nbinom(shape = 4)),
  function(s) s$setModel(build(yCount, family = nbinom(shape = 10))$model),
  "nbinom setModel shape"
)

# Student-t setResponse keeps nu, which the next run's draw starts from, and
# the scales
heavy <- families$student$make()
invisible(heavy$run(5L, 2L))
nu <- stateOf(heavy, "resid.df")
expect_true(all(
  unlist(nu) != stateOf(families$student$make(), "resid.df")[[1L]]
))
heavy$setResponse(rev(yHeavy))
expect_identical(stateOf(heavy, "resid.df"), nu)

# --- what the setters rebuild ------------------------------------------------
# setResponse recomputes what creation derives from the response - the ordinal
# proposal scale from the category counts, the negative-binomial count kernel
# from the histogram and, under updateScale, its shift - so a sampler created
# with one response, given a stored state and then a new response, runs as one
# created with the new response and given the same state. Each new response
# moves the category counts or the histogram; a permutation would not.
recomputed <- list(
  ordinal = list(
    make = function(y) build(y, family = "ordinal"),
    y = yOrdinal,
    yNew = pmin(yOrdinal + (seq_len(n) %% 3L == 0L), 3),
    updateScale = FALSE
  ),
  nbinom = list(
    make = function(y) build(y, family = nbinom()),
    y = yCount,
    yNew = pmax(yCount - 1, 0) + 2 * (seq_len(n) %% 3L == 0L),
    updateScale = TRUE
  )
)
for (name in names(recomputed)) {
  arm <- recomputed[[name]]
  donor <- arm$make(arm$y)
  invisible(donor$run(5L, 2L))
  donor$storeState()
  handed <- arm$make(arm$y)
  handed$setState(donor$state)
  handed$setResponse(arm$yNew, updateScale = arm$updateScale)
  created <- arm$make(arm$yNew)
  created$setState(donor$state)
  expect_identical(
    handed$run(0L, 4L)$train,
    created$run(0L, 4L)$train,
    info = paste(name, "setResponse")
  )
}

# the same under a mask: the count kernel is the rows in's, so a response
# handed to a masked sampler gives the draws of one handed it before the mask
# came, whichever order, the shift held by updateScale = FALSE; rows out are
# far from the rest, the shape would read them if the kernel did
maskedCount <- yCount
maskedCount[mask == 0] <- maskedCount[mask == 0] + 40
maskedNew <- pmax(maskedCount - 1, 0) + 2 * (seq_len(n) %% 3L == 0L)
maskedMake <- function(y, rows = NULL) {
  sampler <- build(y, family = nbinom())
  if (!is.null(rows)) {
    sampler$setActiveRows(rows)
  }
  sampler
}
donor <- maskedMake(maskedCount, mask)
invisible(donor$run(5L, 2L))
donor$storeState()
handed <- maskedMake(maskedCount, mask)
handed$setState(donor$state)
handed$setResponse(maskedNew, updateScale = FALSE)
late <- maskedMake(maskedCount)
late$setState(donor$state)
late$setResponse(maskedNew, updateScale = FALSE)
late$setActiveRows(mask)
expect_identical(
  handed$run(0L, 4L)$train,
  late$run(0L, 4L)$train,
  info = "nbinom setResponse under a mask"
)
# a mask cleared reads every row again
cleared <- maskedMake(maskedCount, mask)
cleared$setState(donor$state)
cleared$setActiveRows(NULL)
unmasked <- maskedMake(maskedCount)
unmasked$setState(donor$state)
expect_identical(
  cleared$run(0L, 4L)$train,
  unmasked$run(0L, 4L)$train,
  info = "nbinom mask cleared"
)
expect_false(
  identical(
    maskedMake(maskedCount, mask)$run(0L, 4L)$train,
    unmasked$run(0L, 4L)$train
  ),
  info = "the mask moves the draws"
)

# --- setData keeps the drawn family parameters ------------------------------
# and starts the latents as creation does, drawing nothing
keptBySetData <- list(
  ordinal = list(
    make = families$ordinal$make,
    data = dbartsData(xNew, 4 - yOrdinal),
    name = "thresholds"
  ),
  student = list(
    make = function() build(yHeavy, family = student()),
    data = dbartsData(xNew, rev(yHeavy)),
    name = "resid.df"
  ),
  nbinom = list(
    make = families$nbinom$make,
    data = dbartsData(xNew, rev(yCount)),
    name = "shape"
  )
)
for (name in names(keptBySetData)) {
  arm <- keptBySetData[[name]]
  sampler <- arm$make()
  invisible(sampler$run(5L, 2L))
  kept <- stateOf(sampler, arm$name)
  stream <- streamOf(sampler)
  sampler$setData(arm$data)
  expect_identical(stateOf(sampler, arm$name), kept, info = name)
  expect_identical(streamOf(sampler), stream, info = name)
  expect_true(all(is.finite(sampler$run(0L, 2L)$train)), info = name)
}

# an ordinal setData holding fewer categories than the sampler keeps all K - 1
# cutpoints
ordinal <- families$ordinal$make()
invisible(ordinal$run(5L, 2L))
cutpoints <- stateOf(ordinal, "thresholds")
expect_identical(length(cutpoints[[1L]]), 2L)
ordinal$setData(dbartsData(x, pmin(yOrdinal, 2)))
expect_identical(stateOf(ordinal, "thresholds"), cutpoints)
expect_true(all(is.finite(ordinal$run(0L, 2L)$train)))

rm(
  nbinom,
  student,
  n,
  x,
  xNew,
  f,
  yBinary,
  yOrdinal,
  yCount,
  yHeavy,
  counts,
  positive,
  logTime,
  bound,
  status,
  mask,
  offset,
  control,
  build,
  streamOf,
  stateOf,
  expectDrawsNothing,
  families,
  name,
  arm,
  newLogTime,
  heavy,
  nu,
  keptBySetData,
  recomputed,
  donor,
  handed,
  created,
  maskedCount,
  maskedNew,
  maskedMake,
  late,
  cleared,
  unmasked,
  sampler,
  kept,
  stream,
  ordinal,
  cutpoints
)
