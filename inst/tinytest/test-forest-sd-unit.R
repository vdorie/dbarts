# A forest's sd is a number of the response's units wherever it is written or
# read: forest(sd = ) at creation, $setLeafPrior(forests = ), what
# $getLeafPrior() returns and extract(type = "leaf.prior.sd"). The engine
# divides a stated number by a unit it derives once, at a first creation, and
# R records; a default never meets the unit. Each path below states an sd and
# reads it back three ways: the reader's `sd`, which is the number as written,
# `sd.stated`, and the forest's own scale, which is where a unit applied twice
# or not at all shows.

forest <- dbartsForests$forest
fixed <- dbartsPriors$fixed

set.seed(29)
n <- 150L
d <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
d$z <- rbinom(n, 1L, 0.5)
d$zf <- factor(d$z)
d$dose <- runif(n, 1, 30)
d$off <- 0.5 * d$x2
d$w <- runif(n, 0.5, 2)
d$w[seq(3L, n, length.out = 10L)] <- 0
d$y <- 30 + 6 * (d$x1 + d$z * (1 + d$x2)) + 0.2 * d$dose + rnorm(n)
d$yb <- rbinom(n, 1L, 0.4)
x <- as.matrix(d[c("x1", "x2")])
y <- d$y

unitControl <- function(...) {
  dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.samples = 20L,
    n.burn = 0L,
    updateState = FALSE,
    seed = 31L,
    ...
  )
}
forestInfo <- function(sampler) {
  attr(sampler$control, "bartcore.forests", exact = TRUE)
}
unitOf <- function(sampler) forestInfo(sampler)$unit
numWarnings <- 0L
counted <- function(expr, class = "warning") {
  withCallingHandlers(expr, warning = function(w) {
    if (inherits(w, class)) {
      numWarnings <<- numWarnings + 1L
    }
    invokeRestart("muffleWarning")
  })
}
# forest `f` of `sampler` states `value`: the reader returns the number
# itself, and the forest's own scale is the number's. A forest with a basis
# carries it in k.scale, over the half-normal median and per unit of row
# norm; one with none carries it as its half-Cauchy median, in units of the
# unit.
expectStated <- function(sampler, f, value, info) {
  prior <- sampler$getLeafPrior(f)
  expect_identical(prior$sd, value, info = info)
  expect_true(prior$sd.stated, info = info)
  inForce <- if (is.null(prior$amplitude.prior.scale)) {
    prior$k.scale * 0.674 * prior$basis.row.norm
  } else {
    prior$amplitude.prior.scale * unitOf(sampler)
  }
  expect_equal(inForce, value, tolerance = 1e-12, info = info)
}
twoStated <- list(forest(sd = 2), forest(basis = d$zf, sd = 0.5))
twoForests <- function(forests = twoStated, response = y, ...) {
  dbarts(x, response, forests = forests, control = unitControl(), ...)
}

# --- which model was fitted: one model at every door ---
doors <- list(
  term = function() {
    dbarts(
      y ~ forest(x1 + x2, sd = 2) + forest(x1 + x2, basis = zf, sd = 0.5),
      d,
      control = unitControl()
    )
  },
  listOnFormula = function() {
    dbarts(
      y ~ x1 + x2,
      d,
      forests = list(forest(sd = 2), forest(basis = ~zf, sd = 0.5)),
      control = unitControl()
    )
  },
  listOnMatrix = function() twoForests(),
  spec = function() {
    spec <- dbartsSpec(dbartsData(x, y), unitControl(), forests = twoStated)
    new("dbartsSampler", spec$control, spec$model, spec$data)
  },
  dataObject = function() {
    dbarts(
      dbartsData(x, y, bases = list(NULL, cbind(1 - d$z, d$z))),
      forests = list(forest(sd = 2), forest(sd = 0.5)),
      control = unitControl()
    )
  }
)
doorDraws <- list()
for (door in names(doors)) {
  sampler <- doors[[door]]()
  expectStated(sampler, 1L, 2, door)
  expectStated(sampler, 2L, 0.5, door)
  expect_equal(unitOf(sampler), sd(y), tolerance = 1e-12, info = door)
  doorDraws[[door]] <- sampler$run()$train
}
for (door in names(doors)[-1L]) {
  expect_identical(doorDraws[[door]], doorDraws[[1L]], info = door)
}
bartFit <- bart(
  y ~ forest(x1 + x2, sd = 2) + forest(x1 + x2, basis = zf, sd = 0.5),
  d,
  n.samples = 5L,
  n.burn = 0L,
  n.chains = 1L,
  n.threads = 1L,
  seed = 31L,
  verbose = FALSE,
  keepSampler = TRUE
)
expectStated(bartFit$fit, 1L, 2, "bart")
expectStated(bartFit$fit, 2L, 0.5, "bart")

# --- the unit written out: sd(y - offset) over the rows in the likelihood
# when the sampler is created, and held from then on ---
withMissing <- d
withMissing$y[5L] <- NA
keep <- 1:140
kept <- setdiff(keep, 5L)
kept <- kept[d$w[kept] > 0]
rowsIn <- sd((d$y - d$off)[kept])
threeForests <- function(sdOnFactor = NULL) {
  dbarts(
    y ~ forest(x1 + x2) +
      forest(x1, basis = zf, sd = sdOnFactor) +
      forest(x2, basis = dose),
    withMissing,
    subset = keep,
    weights = w,
    offset = off,
    control = unitControl()
  )
}
# the rows of weight zero raise creation's warning of them, and no fallback
numWarnings <- 0L
unstated <- counted(threeForests(), "dbartsScaleFallbackWarning")
expect_identical(numWarnings, 0L)
expect_equal(unitOf(unstated), rowsIn, tolerance = 1e-12)
onFactor <- unstated$getLeafPrior(2L)
expect_false(onFactor$sd.stated)
expect_equal(onFactor$sd, sqrt(2 / 3) * rowsIn, tolerance = 1e-12)
expect_identical(onFactor$sd, sqrt(2 / 3) * unitOf(unstated))
# the unit is creation's: neither the weights nor a mask moves it
unstated$setWeights(rep(1, length(unstated$data@y)))
expect_identical(unstated$getLeafPrior(2L), onFactor)
unstated$setActiveRows(rep(
  c(TRUE, FALSE),
  length.out = length(unstated$data@y)
))
expect_identical(unstated$getLeafPrior(2L)$sd, onFactor$sd)
# and a stated sd is divided by that unit, so one from the wrong rows shows
statedOnFactor <- suppressWarnings(threeForests(0.5))
expectStated(statedOnFactor, 2L, 0.5, "rows in the likelihood")
statedOnFactor$setWeights(rep(1, length(statedOnFactor$data@y)))
expectStated(statedOnFactor, 2L, 0.5, "after setWeights")
# where the rows in hold one value the unit is over every row, under
# creation's one warning
flat <- y
flat[d$w > 0] <- 7
numWarnings <- 0L
fallback <- counted(
  dbarts(
    x,
    flat,
    weights = d$w,
    forests = list(forest(), forest(basis = d$zf)),
    control = unitControl()
  ),
  "dbartsScaleFallbackWarning"
)
expect_identical(numWarnings, 1L)
expect_equal(unitOf(fallback), sd(flat), tolerance = 1e-12)
expect_equal(fallback$getLeafPrior(2L)$sd, sd(flat), tolerance = 1e-12)

# --- probit and logistic: the unit is the link's own error sd ---
for (family in c("probit", "logistic")) {
  latent <- twoForests(
    list(forest(sd = 0.5), forest(basis = d$zf, sd = 0.5)),
    d$yb,
    family = family
  )
  expectStated(latent, 1L, 0.5, family)
  expectStated(latent, 2L, 0.5, family)
  expect_equal(
    unitOf(latent),
    if (family == "probit") 1 else pi / sqrt(3),
    tolerance = 1e-15,
    info = family
  )
}

# --- not stated: the engine's own multiple, which never meets the unit ---
for (family in c("gaussian", "probit", "logistic")) {
  response <- if (family == "gaussian") y else d$yb
  for (numForests in 2:3) {
    forests <- list(forest(), forest(basis = d$zf), forest(basis = d$dose))
    sampler <- twoForests(
      forests[seq_len(numForests)],
      response,
      family = family
    )
    info <- paste(family, numForests)
    plain <- sampler$getLeafPrior(1L)
    median <- if (family == "gaussian") 2 else 1
    expect_identical(plain$amplitude.prior.scale, median, info = info)
    expect_false(plain$sd.stated, info = info)
    expect_identical(plain$sd, median * unitOf(sampler), info = info)
    for (f in seq_len(numForests)[-1L]) {
      multiplied <- sampler$getLeafPrior(f)
      expect_identical(
        multiplied$leaf.scale.factor,
        sqrt(2 / numForests),
        info = info
      )
      expect_false(multiplied$sd.stated, info = info)
      expect_identical(
        multiplied$sd,
        sqrt(2 / numForests) * unitOf(sampler),
        info = info
      )
    }
    expect_identical(forestInfo(sampler)$sd, rep(NA_real_, numForests))
  }
}
# a stated variance reaches the engine as written
withVariance <- twoForests(
  list(forest(), forest(basis = d$zf, amplitude.prior.variance = 2))
)
expect_identical(withVariance$getLeafPrior(2L)$amplitude.prior.variance, 2)
expect_identical(forestInfo(withVariance)$variance, c(NA, 2))

# --- one path at a time ---
created <- twoForests()
expectStated(created, 1L, 2, "creation")
expectStated(created, 2L, 0.5, "creation")

written <- twoForests(list(forest(), forest(basis = d$zf)))
written$setLeafPrior(forests = list(forest(sd = 3), forest(sd = 0.25)))
expectStated(written, 1L, 3, "the writer")
expectStated(written, 2L, 0.25, "the writer")

copied <- created$copy()
expectStated(copied, 1L, 2, "copy")
expectStated(copied, 2L, 0.5, "copy")

savedFile <- tempfile(fileext = ".rds")
created$storeState()
saveRDS(created, savedFile)
reloaded <- readRDS(savedFile)
unlink(savedFile)
invisible(reloaded$run(0L, 1L))
expectStated(reloaded, 1L, 2, "reload")
expectStated(reloaded, 2L, 0.5, "reload")

remade <- new("dbartsSampler", created$control, created$model, created$data)
expectStated(remade, 1L, 2, "new() from the sampler's own three")
expectStated(remade, 2L, 0.5, "new() from the sampler's own three")

restated <- twoForests()
restated$storeState()
restated$setState(restated$state)
expectStated(restated, 1L, 2, "setState")
expectStated(restated, 2L, 0.5, "setState")

# a response swapped without re-deriving the scale keeps the unit, on the
# sampler and on a copy made then
swapped <- twoForests()
priorBefore <- swapped$getLeafPrior()
swapped$setResponse(10 * y + 5, updateScale = FALSE)
expect_identical(swapped$getLeafPrior(), priorBefore)
expectStated(swapped, 2L, 0.5, "setResponse")
swappedCopy <- swapped$copy()
expect_identical(swappedCopy$getLeafPrior(), priorBefore)
expect_identical(unitOf(swappedCopy), unitOf(created))

# a write, a copy and a second write: the record's sd is what a copy states
writtenCopy <- written$copy()
expectStated(writtenCopy, 1L, 3, "a write, then a copy")
expectStated(writtenCopy, 2L, 0.25, "a write, then a copy")
writtenCopy$setLeafPrior(forests = list(forest(), forest(sd = 0.75)))
expectStated(writtenCopy, 1L, 3, "a second write")
expectStated(writtenCopy, 2L, 0.75, "a second write")
expect_identical(forestInfo(writtenCopy)$sd, c(3, 0.75))
expect_identical(forestInfo(written)$sd, c(3, 0.25))

# a swap to another factor of the same levels derives the forest's own scale
# again, from the stated number
swappedBasis <- twoForests()
swappedBasis$setForestBasis(2L, factor(1L - d$z))
expectStated(swappedBasis, 2L, 0.5, "setForestBasis")
# a column swapped for one in other units moves the row norm and leaves the
# stated number: the forest is the one created over the new column, on the
# sampler and on a copy
overColumn <- function(column) {
  twoForests(list(forest(sd = 2), forest(basis = column, sd = 0.5)))
}
swappedColumn <- overColumn(d$dose)
normBefore <- swappedColumn$getLeafPrior(2L)$basis.row.norm
swappedColumn$setForestBasis(2L, 3 * d$dose + 1)
afterSwap <- swappedColumn$getLeafPrior(2L)
expect_equal(normBefore, median(d$dose), tolerance = 1e-12)
expect_equal(
  afterSwap$basis.row.norm,
  3 * median(d$dose) + 1,
  tolerance = 1e-12
)
expectStated(swappedColumn, 2L, 0.5, "setForestBasis, a column in other units")
expect_identical(afterSwap$leaf.scale.factor, 0.5 / unitOf(swappedColumn))
expect_identical(afterSwap, overColumn(3 * d$dose + 1)$getLeafPrior(2L))
expect_identical(swappedColumn$copy()$getLeafPrior(2L), afterSwap)
# and a column given to the forest created with no basis leaves that forest
# the kind it was created as, on the sampler and on a copy
swappedBasis$setForestBasis(1L, d$dose)
expect_identical(swappedBasis$getLeafPrior(1L)$sd, 2)
expect_identical(swappedBasis$getLeafPrior(1L)$prior.sd.of, "amplitude scale")
expect_identical(
  swappedBasis$copy()$getLeafPrior(1L),
  swappedBasis$getLeafPrior(1L)
)

# --- read, then write ---
sdsOf <- function(sampler) {
  lapply(sampler$getLeafPrior(), function(prior) prior[c("sd", "sd.stated")])
}
for (forests in list(twoStated, list(forest(), forest(basis = d$zf)))) {
  twin <- twoForests(forests)
  readBack <- twoForests(forests)
  invisible(twin$run(0L, 5L))
  invisible(readBack$run(0L, 5L))
  readBack$setLeafPrior(
    forests = lapply(readBack$getLeafPrior(), function(prior) prior$leaf.prior)
  )
  expect_identical(sdsOf(readBack), sdsOf(twin))
  expect_identical(forestInfo(readBack)$sd, forestInfo(twin)$sd)
  expect_identical(readBack$run(0L, 20L), twin$run(0L, 20L))
}
# a default's number, stated: the forest is a stated one, and no bit moves,
# on the sampler and on a copy made from the record
for (numForests in 2:3) {
  forests <- list(forest(), forest(basis = d$zf), forest(basis = d$dose))
  forests <- forests[seq_len(numForests)]
  twin <- twoForests(forests)
  stated <- twoForests(forests)
  invisible(twin$run(0L, 5L))
  invisible(stated$run(0L, 5L))
  inForce <- lapply(stated$getLeafPrior(), function(prior) prior$sd)
  stated$setLeafPrior(forests = lapply(inForce, function(sd) forest(sd = sd)))
  expect_true(all(vapply(stated$getLeafPrior(), `[[`, NA, "sd.stated")))
  expect_identical(lapply(stated$getLeafPrior(), `[[`, "sd"), inForce)
  expect_identical(
    lapply(stated$getLeafPrior(), `[[`, "k.scale"),
    lapply(twin$getLeafPrior(), `[[`, "k.scale")
  )
  expect_identical(
    lapply(stated$copy()$getLeafPrior(), `[[`, "k.scale"),
    lapply(twin$getLeafPrior(), `[[`, "k.scale")
  )
  expect_identical(stated$run(0L, 20L), twin$run(0L, 20L))
}
# the same where the unit and back do not return sqrt(2 / 3): the response
# is stretched until its unit is such a one, and the default's number stated
# is the default's multiple at a write, on a copy made from the record and at
# a creation that states it
third <- sqrt(2 / 3)
threeDefaults <- list(forest(), forest(basis = d$zf), forest(basis = d$dose))
inexact <- NULL
for (stretch in seq(1, 4, by = 0.125)) {
  unit <- unitOf(twoForests(threeDefaults, stretch * y))
  if (third * unit / unit != third) {
    inexact <- stretch * y
    break
  }
}
expect_false(is.null(inexact))
twin <- twoForests(threeDefaults, inexact)
written <- twoForests(threeDefaults, inexact)
inForce <- lapply(written$getLeafPrior(), function(prior) prior$sd)
written$setLeafPrior(forests = lapply(inForce, function(sd) forest(sd = sd)))
channelsOf <- function(sampler) {
  lapply(sampler$getLeafPrior(), function(prior) {
    list(
      sd = prior$sd,
      k.scale = prior$k.scale,
      factor = prior$leaf.scale.factor,
      median = prior$amplitude.prior.scale
    )
  })
}
statingDefaults <- list(
  write = written,
  copy = written$copy(),
  creation = twoForests(
    list(
      forest(sd = inForce[[1L]]),
      forest(basis = d$zf, sd = inForce[[2L]]),
      forest(basis = d$dose, sd = inForce[[3L]])
    ),
    inexact
  )
)
twinChannels <- channelsOf(twin)
twinDraws <- twin$run(0L, 20L)
for (path in names(statingDefaults)) {
  sampler <- statingDefaults[[path]]
  expect_true(
    all(vapply(sampler$getLeafPrior(), `[[`, NA, "sd.stated")),
    info = path
  )
  expect_identical(channelsOf(sampler), twinChannels, info = path)
  expect_identical(
    sampler$getLeafPrior(2L)$leaf.scale.factor,
    third,
    info = path
  )
  expect_identical(sampler$run(0L, 20L), twinDraws, info = path)
}

# --- other data: the unit is not carried on a forests list or a control ---
wide <- dbarts(x, 10 * y, forests = twoStated, control = created$control)
expectStated(wide, 1L, 2, "other data")
expectStated(wide, 2L, 0.5, "other data")
expect_equal(unitOf(wide), 10 * sd(y), tolerance = 1e-12)

# --- names, and extract ---
termLabels <- c("forest1", "zf", "dose")
threeTerms <- y ~ forest(x1 + x2, sd = 2) +
  forest(x1, basis = zf, sd = 0.5) +
  forest(x2, basis = dose)
expect_identical(
  names(dbarts(threeTerms, d, control = unitControl())$getLeafPrior()),
  termLabels
)
namedList <- twoForests(
  list(prognostic = forest(sd = 2), effect = forest(basis = d$zf, sd = 0.5))
)
expect_identical(names(namedList$getLeafPrior()), c("prognostic", "effect"))
threeFit <- bart(
  threeTerms,
  d,
  n.samples = 5L,
  n.burn = 0L,
  n.chains = 1L,
  n.threads = 1L,
  seed = 31L,
  verbose = FALSE,
  keepSampler = TRUE
)
readerSds <- lapply(threeFit$fit$getLeafPrior(), function(prior) prior$sd)
expect_identical(names(readerSds), termLabels)
expect_identical(extract(threeFit, type = "leaf.prior.sd"), readerSds)
expect_identical(readerSds[1:2], list(forest1 = 2, zf = 0.5))
expect_identical(extract(threeFit, type = "leaf.prior.sd", forest = 2L), 0.5)
expect_identical(extract(threeFit, type = "leaf.prior.sd", forest = "zf"), 0.5)
expect_identical(
  extract(threeFit, type = "leaf.prior.sd", forest = "forest2"),
  0.5
)
expect_identical(
  extract(threeFit, type = "leaf.prior.sd", forest = c("dose", "forest1")),
  readerSds[c("dose", "forest1")]
)
# on one forest it is what it was
oneFit <- bart(
  y ~ x1 + x2,
  d,
  leaf.prior = normal(sd = 2),
  n.samples = 5L,
  n.burn = 0L,
  n.chains = 1L,
  n.threads = 1L,
  seed = 31L,
  verbose = FALSE
)
expect_identical(extract(oneFit, type = "leaf.prior.sd"), 2)

# --- one forest: its sd is the fitting function's normal(sd = ) ---
oneControl <- function() unitControl(n.trees = 12L)
for (family in c("gaussian", "probit")) {
  response <- if (family == "gaussian") y else d$yb
  frame <- d
  frame$response <- response
  reference <- dbarts(
    x,
    response,
    family = family,
    leaf.prior = normal(sd = 2),
    control = oneControl()
  )$run()
  expect_identical(
    dbarts(
      response ~ forest(x1 + x2, sd = 2),
      frame,
      family = family,
      control = oneControl()
    )$run(),
    reference,
    info = family
  )
  expect_identical(
    dbarts(
      response ~ x1 + x2,
      frame,
      family = family,
      forests = list(forest(sd = 2)),
      control = oneControl()
    )$run(),
    reference,
    info = family
  )
  expect_identical(
    dbarts(
      x,
      response,
      family = family,
      forests = list(forest(sd = 2)),
      control = oneControl()
    )$run(),
    reference,
    info = family
  )
}
bartArgs <- list(
  n.trees = 12L,
  n.samples = 5L,
  n.burn = 0L,
  n.chains = 1L,
  n.threads = 1L,
  seed = 31L,
  verbose = FALSE
)
expect_identical(
  do.call(bart, c(list(y ~ forest(x1 + x2, sd = 2), d), bartArgs))$yhat.train,
  do.call(
    bart,
    c(list(y ~ x1 + x2, d, leaf.prior = quote(normal(sd = 2))), bartArgs)
  )$yhat.train
)
statedTwice <- paste0(
  "this model has one forest and states its size twice: 'sd' on the forest ",
  "and 'leaf.prior' given to the fitting function are one statement; give one"
)
for (leafPrior in list(
  quote(normal(sd = 3)),
  quote(normal),
  quote(normal(k = 2)),
  quote(linear(columns = "x1"))
)) {
  expect_error(
    eval(bquote(dbarts(
      y ~ forest(x1 + x2, sd = 2),
      d,
      leaf.prior = .(leafPrior),
      control = oneControl()
    ))),
    statedTwice,
    fixed = TRUE,
    info = deparse(leafPrior)
  )
}
expect_error(
  do.call(bart, c(list(y ~ forest(x1 + x2, sd = 2), d, k = 3), bartArgs)),
  sub("'leaf.prior'", "'k'", statedTwice, fixed = TRUE),
  fixed = TRUE
)
expect_error(
  dbarts(
    x,
    y,
    forests = list(forest(amplitude = fixed())),
    control = oneControl()
  ),
  "this model has one forest, which has no coefficient to hold",
  fixed = TRUE
)

# --- a response with no spread has no unit to state an sd in ---
noSpread <- paste0(
  "a forest's 'sd' is in the response's units, and this response has no ",
  "spread to state one against (its standard deviation is 0); leave 'sd' ",
  "out, or give a response that varies"
)
constant <- rep(3, n)
onConstant <- suppressWarnings(
  twoForests(list(forest(), forest(basis = d$zf)), constant)
)
expect_identical(unitOf(onConstant), 0)
expect_error(
  suppressWarnings(twoForests(
    list(forest(), forest(basis = d$zf, sd = 0.7)),
    constant
  )),
  noSpread,
  fixed = TRUE
)
expect_error(
  onConstant$setLeafPrior(forests = list(forest(sd = 0.7))),
  noSpread,
  fixed = TRUE
)
# the fitting function's normal(sd = ) beside several forests says where the
# size is written
expect_error(
  twoForests(list(forest(), forest(basis = d$zf)), leaf.prior = normal(sd = 2)),
  paste0(
    "a model of several forests states each forest's size on the forest: ",
    "write the size of the forest with no basis as forest(x1 + x2, sd = 2), ",
    "which is in the response's units as leaf.prior = normal(sd = ) is, and ",
    "leave 'sd' out of 'leaf.prior'"
  ),
  fixed = TRUE
)

# --- a named single sd is the number, its name dropped ---
pars <- c(s = 0.5, t = 2)
unnamedDraws <- twoForests(
  list(forest(sd = 2), forest(basis = d$zf, sd = 0.5))
)$run()$train
namedOnList <- twoForests(
  list(forest(sd = pars["t"]), forest(basis = d$zf, sd = pars["s"]))
)
expect_identical(namedOnList$getLeafPrior(1L)$sd, 2)
expect_identical(namedOnList$getLeafPrior(2L)$sd, 0.5)
expect_identical(forestInfo(namedOnList)$sd, c(2, 0.5))
expect_identical(namedOnList$run()$train, unnamedDraws)
namedOnTerm <- dbarts(
  y ~ forest(x1 + x2, sd = pars["t"]) +
    forest(x1 + x2, basis = zf, sd = pars["s"]),
  d,
  control = unitControl()
)
expect_identical(namedOnTerm$getLeafPrior(2L)$sd, 0.5)
expect_identical(namedOnTerm$run()$train, unnamedDraws)
namedOnColumn <- twoForests(
  list(forest(), forest(basis = d$dose, sd = pars["s"]))
)
expect_identical(namedOnColumn$getLeafPrior(2L)$sd, 0.5)
namedWrite <- twoForests(list(forest(), forest(basis = d$zf)))
namedWrite$setLeafPrior(
  forests = list(forest(sd = pars["t"]), forest(sd = pars["s"]))
)
expect_identical(namedWrite$getLeafPrior(2L)$sd, 0.5)
expect_identical(forestInfo(namedWrite)$sd, c(2, 0.5))
expect_identical(namedWrite$run()$train, unnamedDraws)
expect_identical(
  dbarts(
    x,
    y,
    forests = list(forest(sd = pars["t"])),
    control = oneControl()
  )$run(),
  dbarts(x, y, leaf.prior = normal(sd = 2), control = oneControl())$run()
)
longSd <- paste0(
  "forest 'sd' must be a single number, not a vector of length 2: a forest ",
  "states one sd; to size a column differently, rescale it in 'basis', as ",
  "I(dose / 30)"
)
for (long in list(c(dose = 1, age = 2), c(1, 2))) {
  expect_error(
    twoForests(list(forest(), forest(basis = d$zf, sd = long))),
    longSd,
    fixed = TRUE
  )
  expect_error(
    dbarts(
      y ~ forest(x1 + x2) + forest(x1, basis = zf, sd = long),
      d,
      control = unitControl()
    ),
    longSd,
    fixed = TRUE
  )
  expect_error(
    namedWrite$setLeafPrior(forests = list(forest(), forest(sd = long))),
    longSd,
    fixed = TRUE
  )
}

# --- the size, from the engine's own prior draws under a flat likelihood: a
# forest that is drawn and states sd = v has size v, the prior median of the
# absolute multiplier times the prior sd of the forest's own fit ---
flatSampler <- function(forests, family) {
  sizeControl <- dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 20L,
    n.samples = 1L,
    n.burn = 0L,
    n.thin = 3L,
    updateState = FALSE,
    seed = 37L
  )
  if (family == "gaussian") {
    return(dbarts(
      x,
      y,
      family = gaussian(sigma = fixed(1e16)),
      forests = forests,
      control = sizeControl
    ))
  }
  sampler <- dbarts(
    x,
    d$yb,
    family = family,
    forests = forests,
    control = sizeControl
  )
  sampler$setActiveRows(rep(FALSE, n))
  sampler
}
priorSize <- function(sampler, f, multiplier, sweeps) {
  responseScale <- sampler$getLeafPrior(f)$response.scale
  fits <- matrix(NA_real_, sweeps, n)
  multipliers <- numeric(sweeps)
  for (sweep in seq_len(sweeps)) {
    sampler$run(0L, 1L)
    fits[sweep, ] <- sampler$getForestFits(f)[, 1L] * responseScale
    multipliers[sweep] <- multiplier(sampler$getForestAmplitudes(f)[, 1L])
  }
  median(abs(multipliers)) * mean(apply(fits, 2L, sd))
}
sizeArm <- function(values, sweeps, tolerance) {
  for (family in c("gaussian", "probit", "logistic")) {
    for (value in values) {
      # in the response's units: tens under gaussian, where sd(y) is near 10
      stated <- if (family == "gaussian") 10 * value else value
      sampler <- flatSampler(
        list(forest(sd = stated), forest(basis = d$zf, sd = stated)),
        family
      )
      info <- paste(family, stated)
      expect_equal(
        priorSize(sampler, 1L, function(a) a[1L], sweeps) / stated,
        1,
        tolerance = tolerance,
        info = info
      )
      expect_equal(
        priorSize(sampler, 2L, function(a) a[2L] - a[1L], sweeps) / stated,
        1,
        tolerance = tolerance,
        info = info
      )
    }
  }
}
sizeArm(0.7, 1500L, 0.1)
if (at_home()) {
  sizeArm(c(0.3, 0.7, 1.5), 5000L, 0.06)
}
