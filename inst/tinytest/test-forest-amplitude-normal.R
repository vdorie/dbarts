# The trees carry a forest's sd on every forest, held or drawn; a drawn
# coefficient is a multiplier with prior median size 1, whose prior does not
# read the sd; a live write sets the trees' scale and moves nothing drawn (dec-B550). And
# forest(amplitude = normal()) on a forest with no basis, taken where student(df)
# is and refused where it is (dec-B539).

forest <- dbartsForests$forest
fixed <- dbartsPriors$fixed
studentPrior <- dbartsPriors$student
normalPrior <- dbartsPriors$normal

set.seed(20261031)
n <- 100L
x1 <- runif(n)
x2 <- runif(n)
x3 <- runif(n)
z <- as.double(seq_len(n) %in% sample(n, 45L))
zf <- factor(z)
dose <- rnorm(n, 5, 2)
y <- 1 + sin(3 * x1) + z * (1 + x2) + 0.2 * dose * x2 + rnorm(n, sd = 0.5)
yb <- as.double(y > median(y))
x <- cbind(x1 = x1, x2 = x2)
d <- data.frame(y, yb, x1, x2, x3, zf, dose)

ctl <- function(chains = 1L, ...) {
  dbartsControl(
    n.chains = chains,
    n.threads = 1L,
    n.samples = 5L,
    n.burn = 0L,
    updateState = FALSE,
    ...
  )
}
make <- function(forests, ..., response = y, chains = 1L) {
  dbarts(
    x,
    response,
    forests = forests,
    ...,
    n.trees = 10L,
    n.cuts = 20L,
    control = ctl(chains, seed = 7L)
  )
}

# --- the creation and the reader, before anything else reads them: each is
# its own expectation, and an error is a failed one, so that a refusal of
# normal() or of the reader's answer does not stop the file ---
readAmplitude <- function(amplitude) {
  tryCatch(
    make(list(
      forest(sd = 1, amplitude = amplitude),
      forest(basis = zf)
    ))$getLeafPrior(1L)$amplitude,
    error = function(e) conditionMessage(e)
  )
}
expect_true(
  is(readAmplitude(normalPrior()), "dbartsNormalPrior"),
  info = "the reader gives back normal()"
)
expect_true(
  is(readAmplitude(normalPrior), "dbartsNormalPrior"),
  info = "a bare normal is normal()"
)

# the cases of a forest's kind and coefficient option: where it stands, and the
# model at a given sd
cases <- list(
  none.normal = list(1L, function(sd) {
    list(forest(sd = sd, amplitude = normalPrior()), forest(basis = zf))
  }),
  none.student3 = list(1L, function(sd) {
    list(forest(sd = sd, amplitude = studentPrior(3)), forest(basis = zf))
  }),
  none.student1 = list(1L, function(sd) {
    list(forest(sd = sd, amplitude = studentPrior(1)), forest(basis = zf))
  }),
  none.held = list(1L, function(sd) {
    list(forest(sd = sd, amplitude = fixed()), forest(basis = zf))
  }),
  factor.normal = list(2L, function(sd) {
    list(forest(), forest(basis = zf, sd = sd))
  }),
  factor.student3 = list(2L, function(sd) {
    list(forest(), forest(basis = zf, sd = sd, amplitude = studentPrior(3)))
  }),
  factor.held = list(2L, function(sd) {
    list(forest(), forest(basis = zf, sd = sd, amplitude = fixed()))
  }),
  numeric.normal = list(2L, function(sd) {
    list(forest(), forest(basis = dose, sd = sd))
  }),
  numeric.student3 = list(2L, function(sd) {
    list(forest(), forest(basis = dose, sd = sd, amplitude = studentPrior(3)))
  }),
  numeric.held = list(2L, function(sd) {
    list(forest(), forest(basis = dose, sd = sd, amplitude = fixed()))
  })
)
writeSd <- function(sampler, index, sd) {
  sampler$setForest(index, sd = sd)
  invisible(sampler)
}
coefficients <- function(sampler, index) sampler$getForestAmplitudes(index)
# the chain's glue as the state carries it: the coefficients and the mixture
# variances drawn beside them, which a write of sd leaves alone
auxiliaries <- function(sampler) {
  sampler$storeState()
  sampler$state[[1L]]$glue
}
engineSlot <- function(sampler, index, name) {
  .Call(
    dbarts:::C_dbarts_bartcore_getLeafPrior,
    sampler$getPointer(),
    index - 1L
  )[[1L, name]]
}

# --- the rule, on every kind: doubling sd doubles the trees' scale and only
# that, so no exception can come back (sd-carried) ---
for (name in names(cases)) {
  index <- cases[[name]][[1L]]
  build <- cases[[name]][[2L]]
  one <- make(build(0.6))$getLeafPrior(index)
  two <- make(build(1.2))$getLeafPrior(index)
  expect_equal(unname(one$sd), 0.6, info = name)
  expect_equal(unname(two$sd), 1.2, info = name)
  expect_equal(unname(one$k.scale), 0.6, info = name)
  expect_equal(unname(two$k.scale), 1.2, info = name)
  expect_identical(
    one[names(one) %in% c("amplitude", "multiplier", "prior.sd.of")],
    two[names(two) %in% c("amplitude", "multiplier", "prior.sd.of")],
    info = name
  )
}
# and the engine's own read of the coefficient's prior does not move with it
for (name in setdiff(
  names(cases),
  c("none.held", "factor.held", "numeric.held")
)) {
  index <- cases[[name]][[1L]]
  build <- cases[[name]][[2L]]
  one <- make(build(0.6))
  two <- make(build(1.2))
  for (column in c("amplitude.prior.scale", "amplitude.prior.variance")) {
    expect_identical(
      engineSlot(one, index, column),
      engineSlot(two, index, column),
      info = paste(name, column)
    )
  }
}

# --- a live write moves the trees' scale and nothing drawn ---
for (name in names(cases)) {
  index <- cases[[name]][[1L]]
  build <- cases[[name]][[2L]]
  # in mid-chain
  sampler <- make(build(0.7))
  invisible(sampler$run(0L, 8L))
  before <- coefficients(sampler, index)
  fits <- sampler$getForestFits(index)
  v <- auxiliaries(sampler)
  expect_true(is.numeric(v) && length(v) > 0L, info = name)
  writeSd(sampler, index, 2.1)
  expect_identical(coefficients(sampler, index), before, info = name)
  expect_identical(sampler$getForestFits(index), fits, info = name)
  expect_identical(auxiliaries(sampler), v, info = name)
  expect_equal(unname(sampler$getLeafPrior(index)$sd), 2.1, info = name)
  expect_equal(unname(sampler$getLeafPrior(index)$k.scale), 2.1, info = name)

  # before the first sweep and after a prior draw of the trees
  sampler <- make(build(0.7))
  before <- coefficients(sampler, index)
  writeSd(sampler, index, 2.1)
  expect_identical(coefficients(sampler, index), before, info = name)
  sampler <- make(build(0.7))
  sampler$sampleTreesFromPrior()
  before <- coefficients(sampler, index)
  fits <- sampler$getForestFits(index)
  writeSd(sampler, index, 2.1)
  expect_identical(coefficients(sampler, index), before, info = name)
  expect_identical(sampler$getForestFits(index), fits, info = name)

  # a copy and a reload hold the block and continue at the written sd
  sampler <- make(build(0.7))
  invisible(sampler$run(0L, 8L))
  writeSd(sampler, index, 2.1)
  sampler$storeState()
  copied <- sampler$copy()
  file <- tempfile(fileext = ".rds")
  saveRDS(sampler, file)
  reloaded <- readRDS(file)
  reloadedAgain <- readRDS(file)
  unlink(file)
  expect_identical(
    coefficients(copied, index),
    coefficients(sampler, index),
    info = name
  )
  expect_identical(
    coefficients(reloaded, index),
    coefficients(sampler, index),
    info = name
  )
  expect_equal(unname(copied$getLeafPrior(index)$sd), 2.1, info = name)
  expect_equal(unname(reloaded$getLeafPrior(index)$sd), 2.1, info = name)
  expect_equal(unname(reloaded$getLeafPrior(index)$k.scale), 2.1, info = name)
  expect_identical(
    copied$run(0L, 3L)$train,
    sampler$copy()$run(0L, 3L)$train,
    info = name
  )
  expect_identical(
    reloaded$run(0L, 3L)$train,
    reloadedAgain$run(0L, 3L)$train,
    info = name
  )

  # store, write, sweep, restore, write back
  sampler <- make(build(0.7))
  invisible(sampler$run(0L, 8L))
  before <- coefficients(sampler, index)
  v <- auxiliaries(sampler)
  stored <- sampler$state
  writeSd(sampler, index, 7)
  invisible(sampler$run(0L, 8L))
  expect_true(sampler$setState(stored), info = name)
  writeSd(sampler, index, 0.7)
  expect_identical(coefficients(sampler, index), before, info = name)
  expect_identical(auxiliaries(sampler), v, info = name)
}

# the block keeps creation's start on a written unswept sampler, and a chain on
# a forest whose start does not read sd is creation's, bit for bit
for (name in c(
  "factor.normal",
  "numeric.normal",
  "factor.student3",
  "none.held"
)) {
  index <- cases[[name]][[1L]]
  build <- cases[[name]][[2L]]
  written <- writeSd(make(build(0.7)), index, 2.1)
  expect_identical(
    written$run(0L, 3L)$train,
    make(build(2.1))$run(0L, 3L)$train,
    info = name
  )
}
# the sd in force, written, is no write
for (name in names(cases)) {
  index <- cases[[name]][[1L]]
  build <- cases[[name]][[2L]]
  written <- writeSd(make(build(0.7)), index, 0.7)
  expect_identical(
    written$run(0L, 3L)$train,
    make(build(0.7))$run(0L, 3L)$train,
    info = name
  )
}
# the restate of a numeric forest's default: its trees' scale and nothing drawn
restated <- make(list(forest(), forest(basis = dose)))
invisible(restated$run(0L, 8L))
before <- coefficients(restated, 2L)
fits <- restated$getForestFits(2L)
scale <- restated$getLeafPrior(2L)$k.scale
restated$setForestBasis(2L, 2 * dose, updateBasisScale = TRUE)
expect_identical(coefficients(restated, 2L), before)
expect_identical(restated$getForestFits(2L), fits)
expect_equal(unname(restated$getLeafPrior(2L)$k.scale), unname(scale) / 2)

# --- several chains: a write reaches every chain's trees and moves no chain's
# draws; at a tiny sd each chain's forest is tiny after a few sweeps ---
for (name in names(cases)) {
  index <- cases[[name]][[1L]]
  build <- cases[[name]][[2L]]
  sampler <- make(build(0.7), chains = 2L)
  invisible(sampler$run(0L, 8L))
  before <- coefficients(sampler, index)
  writeSd(sampler, index, 1e-4)
  expect_identical(ncol(before), 2L, info = name)
  expect_identical(coefficients(sampler, index), before, info = name)
  invisible(sampler$run(0L, 3L))
  expect_true(
    all(apply(abs(sampler$getForestFits(index)), 2L, max) < 0.01),
    info = name
  )
}

# --- an extreme sd: the start is kept finite, with no refusal ---
for (sd in c(1e-160, 1e300)) {
  for (prior in list(normalPrior(), studentPrior(3), studentPrior(1))) {
    result <- make(list(
      forest(sd = sd, amplitude = prior),
      forest(basis = zf)
    ))$run(
      0L,
      10L
    )
    expect_true(all(is.finite(result$train)), info = paste(sd, class(prior)))
    expect_true(all(is.finite(result$sigma)), info = paste(sd, class(prior)))
  }
}

# the start is 1 / m, m the multiple of the unit, held to 1e-100 through 1e100,
# so a multiple of 1e60 starts at 1e-60 and a multiple past the bound at its
# edge
unit <- unname(make(list(
  forest(),
  forest(basis = zf)
))$data@forest.calibration[[
  "unit"
]])
startAt <- function(multiple, prior) {
  as.vector(make(list(
    forest(sd = multiple * unit, amplitude = prior),
    forest(basis = zf)
  ))$getForestAmplitudes(1L))
}
for (prior in list(normalPrior(), studentPrior(3), studentPrior(1))) {
  info <- class(prior)
  expect_equal(startAt(1e60, prior) / 1e-60, 1, tolerance = 1e-12, info = info)
  expect_equal(startAt(1e-60, prior), 1e60, tolerance = 1e-12, info = info)
  expect_equal(
    startAt(1e100, prior) / 1e-100,
    1,
    tolerance = 1e-12,
    info = info
  )
  expect_equal(startAt(1e-100, prior), 1e100, tolerance = 1e-12, info = info)
  expect_equal(
    startAt(1e101, prior) / 1e-100,
    1,
    tolerance = 1e-12,
    info = info
  )
  expect_equal(startAt(1e-101, prior), 1e100, tolerance = 1e-12, info = info)
}

# --- forest.Rd's indexSd, as the help has it, against the product's spread ---
# The line is read from the installed help, so that a change to the Rd is what
# this reads, and evaluated at the example's fixture.
exampleLines <- local({
  db <- tools::Rd_db("dbarts")
  file <- tempfile(fileext = ".R")
  on.exit(unlink(file))
  tools::Rd2ex(db[["forest.Rd"]], file, commentDontrun = TRUE)
  readLines(file)
})
indexLine <- grep("^indexSd <- ", exampleLines, value = TRUE)
expect_identical(length(indexLine), 1L)
indexExpression <- parse(text = indexLine)[[1L]]
sdIndex <- 1.5
sampler <- dbarts(
  x,
  y,
  forests = list(forest(), forest(basis = zf, n.trees = 1L, sd = sdIndex)),
  n.trees = 10L,
  n.cuts = 20L,
  control = ctl(seed = 11L)
)
indexValue <- eval(
  indexExpression[[3L]],
  list(
    calibration = sampler$getLeafPrior(2L),
    basis = sampler$data@bases[[2L]]
  )
)
# the coefficients' draw from their prior times one tree at the prior's scale:
# each row's location is a level's coefficient times a normal with sd sd
drawnLevels <- rnorm(2e5, sd = sqrt(0.5) / 0.674)
drawnTrees <- rnorm(2e5, sd = sdIndex)
expect_equal(
  unique(round(indexValue, 8L)),
  round(sqrt(sdIndex^2 * 0.5 / 0.674^2), 8L)
)
expect_true(
  abs(sd(drawnLevels * drawnTrees) / unique(round(indexValue, 8L)) - 1) < 0.02
)

# --- the medians: normal() has student(3)'s no-sd size, under each family ---
families <- list(
  gaussian = list(y, "gaussian"),
  probit = list(yb, "probit"),
  logistic = list(yb, "logistic")
)
for (family in names(families)) {
  response <- families[[family]][[1L]]
  link <- families[[family]][[2L]]
  for (sd in list(NULL, 0.7)) {
    atPrior <- function(prior) {
      forests <- list(forest(amplitude = prior), forest(basis = zf))
      forests[[1L]]$sd <- sd
      dbarts(
        x,
        response,
        family = link,
        forests = forests,
        n.trees = 10L,
        n.cuts = 20L,
        control = ctl()
      )$getLeafPrior(1L)
    }
    normalRead <- atPrior(normalPrior())
    studentRead <- atPrior(studentPrior(3))
    expect_identical(normalRead$sd, studentRead$sd, info = family)
    expect_identical(
      normalRead$prior.sd.of,
      studentRead$prior.sd.of,
      info = family
    )
    expect_identical(
      normalRead$prior.sd.of,
      "forest (at a coefficient of 1)",
      info = family
    )
    expect_equal(
      unname(normalRead$k.scale),
      unname(normalRead$sd),
      tolerance = 1e-12,
      info = family
    )
  }
}
# a no-sd forest under probit and logistic, in the link's units
expect_equal(
  unname(
    dbarts(
      x,
      yb,
      family = "probit",
      forests = list(forest(amplitude = normalPrior()), forest(basis = zf)),
      n.trees = 10L,
      n.cuts = 20L,
      control = ctl()
    )$getLeafPrior(1L)$sd
  ),
  1.574,
  tolerance = 1e-3
)
expect_equal(
  unname(
    dbarts(
      x,
      yb,
      family = "logistic",
      forests = list(forest(amplitude = normalPrior()), forest(basis = zf)),
      n.trees = 10L,
      n.cuts = 20L,
      control = ctl()
    )$getLeafPrior(1L)$sd
  ),
  2.8547,
  tolerance = 1e-3
)

# --- the readers ---
sampler <- make(list(
  forest(sd = 0.7, amplitude = normalPrior()),
  forest(basis = zf)
))
prior <- sampler$getLeafPrior(1L)
expect_true(is(prior$amplitude, "dbartsNormalPrior"))
start <- coefficients(sampler, 1L)
invisible(sampler$run(0L, 5L))
expect_false(identical(coefficients(sampler, 1L), start))
printed <- capture.output(print(sampler))
expect_true(
  any(grepl("no basis; sd 0.7", printed, fixed = TRUE)),
  info = printed
)
expect_false(any(grepl("student", printed, fixed = TRUE)), info = printed)
fit <- bart(
  y ~ forest(x1, sd = 0.7, amplitude = normal()) + forest(x2, basis = zf),
  d,
  n.trees = 10L,
  n.samples = 5L,
  n.burn = 0L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)
summarized <- capture.output(summary(fit))
expect_true(
  any(grepl("coefficient normal", summarized, fixed = TRUE)),
  info = summarized
)
expect_null(fit$fixed$amplitude)
expect_true(any(grepl(
  "amplitude = normal()",
  capture.output(print(forest(sd = 2, amplitude = normalPrior()))),
  fixed = TRUE
)))

# --- the record: Inf on the forest stating it, NA beside it, absent otherwise,
# and on a basis forest normal() is its default, so the record is as left out ---
statementsOf <- function(forests) {
  sampler <- make(forests)
  vapply(
    seq_along(forests),
    function(index) engineSlot(sampler, index, "amplitude.df"),
    0
  )
}
expect_identical(
  statementsOf(list(
    forest(sd = 1, amplitude = normalPrior()),
    forest(basis = zf)
  )),
  c(Inf, NaN)
)
expect_identical(
  statementsOf(list(forest(), forest(basis = zf, amplitude = normalPrior()))),
  c(NaN, NaN)
)

# and on the record creation reads: Inf on the forest stating it, NA beside it
# and where none is stated, as on a basis forest whose normal() is its default
recordOf <- function(forests) {
  vapply(make(forests)$model@forests, slot, 0, "amplitude.df")
}
expect_identical(
  recordOf(list(forest(sd = 1, amplitude = normalPrior()), forest(basis = zf))),
  c(Inf, NA_real_)
)
expect_identical(
  recordOf(list(forest(), forest(basis = zf))),
  c(NA_real_, NA_real_)
)
expect_identical(
  recordOf(list(forest(), forest(basis = zf, amplitude = normalPrior()))),
  c(NA_real_, NA_real_)
)

# --- a basis forest's normal() is its default: same record, same draws ---
for (basis in list(zf, dose)) {
  written <- make(list(
    forest(),
    forest(basis = basis, amplitude = normalPrior())
  ))
  leftOut <- make(list(forest(), forest(basis = basis)))
  expect_identical(written$getLeafPrior(2L)$amplitude, NULL)
  expect_identical(written$run(0L, 5L)$train, leftOut$run(0L, 5L)$train)
}

# --- every door, a bare normal among them: one chain ---
reference <- dbarts(
  y ~ x1 + x2,
  d,
  forests = list(forest(sd = 0.7, amplitude = normal()), forest(basis = zf)),
  n.trees = 10L,
  n.cuts = 20L,
  control = ctl(seed = 13L)
)$run(0L, 5L)$train
doorDraws <- list(
  bare = dbarts(
    y ~ x1 + x2,
    d,
    forests = list(forest(sd = 0.7, amplitude = normal), forest(basis = zf)),
    n.trees = 10L,
    n.cuts = 20L,
    control = ctl(seed = 13L)
  ),
  term = dbarts(
    y ~ forest(x1 + x2, sd = 0.7, amplitude = normal()) +
      forest(x1 + x2, basis = zf),
    d,
    n.trees = 10L,
    n.cuts = 20L,
    control = ctl(seed = 13L)
  ),
  data = dbarts(
    dbartsData(x, y),
    forests = list(forest(sd = 0.7, amplitude = normal()), forest(basis = zf)),
    n.trees = 10L,
    n.cuts = 20L,
    control = ctl(seed = 13L)
  ),
  object = dbarts(
    y ~ x1 + x2,
    d,
    forests = list(
      forest(sd = 0.7, amplitude = normalPrior()),
      forest(basis = zf)
    ),
    n.trees = 10L,
    n.cuts = 20L,
    control = ctl(seed = 13L)
  )
)
for (door in names(doorDraws)) {
  expect_identical(doorDraws[[door]]$run(0L, 5L)$train, reference, info = door)
}
fromBart <- bart(
  y ~ forest(x1 + x2, sd = 0.7, amplitude = normal()) +
    forest(x1 + x2, basis = zf),
  d,
  n.trees = 10L,
  n.samples = 5L,
  n.burn = 0L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)
expect_true(is(fromBart$leaf.prior[[1L]]$amplitude, "dbartsNormalPrior"))

# --- beside what it stands beside ---
expect_silent(make(list(
  forest(sd = 0.7, amplitude = normalPrior()),
  forest(basis = zf),
  forest(basis = dose)
)))
weighted <- dbarts(
  y ~ x1 + x2,
  d,
  forests = list(forest(sd = 0.7, amplitude = normal()), forest(basis = zf)),
  weights = rep(c(1, 1, 0), length.out = n),
  n.trees = 10L,
  n.cuts = 20L,
  control = ctl()
)
expect_true(all(is.finite(weighted$run(0L, 5L)$train)))
# beside two sharers with sd: it stands outside the share, which keeps its parts
sharers <- function(first) {
  dbarts(
    y ~ forest(x1) + forest(x2) + forest(x3, sd = 0.7, amplitude = first),
    d,
    n.trees = 10L,
    n.cuts = 20L,
    control = ctl(seed = 17L)
  )
}
expect_identical(
  lapply(1:2, function(index) {
    sharers(normalPrior())$getLeafPrior(index)[c("sd", "prior.sd.of")]
  }),
  lapply(1:2, function(index) {
    sharers(studentPrior(3))$getLeafPrior(index)[c("sd", "prior.sd.of")]
  })
)

# --- what is refused, each counted, by exact words ---
refusal <- function(expr) {
  message <- NULL
  count <- 0L
  tryCatch(
    withCallingHandlers(
      expr,
      warning = function(w) {
        count <<- count + 1L
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) message <<- conditionMessage(e)
  )
  expect_identical(count, 0L)
  message
}
lone <- paste0(
  "this model has one forest, which has no coefficient to give a prior; a lone ",
  "forest draws its size through the fitting function's leaf prior, as ",
  "leaf.prior = normal(k = chi(1.5, 2)) or normal(sd = invchi(1.5, s)) with s a ",
  "size in the response's units"
)
expect_identical(refusal(make(list(forest(amplitude = normalPrior())))), lone)
expect_identical(
  refusal(make(list(forest(sd = 1, amplitude = normalPrior())))),
  lone
)
expect_identical(
  refusal(dbarts(
    y ~ forest(x1, amplitude = normal()) + forest(x2) + forest(x3),
    d,
    n.trees = 10L,
    n.cuts = 20L,
    control = ctl()
  )),
  paste0(
    "forest 1: amplitude = normal() draws one forest's coefficient; these ",
    "forests share one total: state each forest's 'sd' to give it its own"
  )
)
taken <- "'amplitude = normal()' takes no 'k' or 'sd': a drawn coefficient's size is the forest's own 'sd'"
for (written in list(
  normalPrior(k = 2),
  normalPrior(sd = 2),
  normalPrior(k = dbartsPriors$chi(1.5, 2))
)) {
  expect_identical(
    refusal(make(list(
      forest(sd = 1, amplitude = written),
      forest(basis = zf)
    ))),
    taken
  )
  expect_identical(
    refusal(make(list(forest(), forest(basis = zf, amplitude = written)))),
    taken
  )
  expect_identical(
    refusal(make(list(forest(), forest(basis = dose, amplitude = written)))),
    taken
  )
}
other <- paste0(
  "a forest's 'amplitude' must be fixed(), which holds its coefficient, normal() ",
  "or student(df), which draw it under that prior, or left out, which draws it ",
  "under its default prior"
)
expect_identical(
  refusal(make(list(forest(amplitude = "normal"), forest(basis = zf)))),
  other
)
expect_identical(
  refusal(studentPrior(Inf)),
  paste0(
    "a coefficient's student(df) takes one finite df of 1 or more; student() is ",
    "student(1), the half-Cauchy, and normal() is the normal"
  )
)
expect_identical(
  refusal(dbarts(
    y ~ x1 + x2,
    d,
    leaf.prior = normalPrior(k = dbartsPriors$chi(1.5, Inf)),
    forests = list(forest(amplitude = normal()), forest(basis = zf)),
    n.trees = 10L,
    n.cuts = 20L,
    control = ctl()
  )),
  paste0(
    "forest 1: amplitude = normal() with no 'sd' takes its size from the ",
    "fitting function's leaf prior, which here gives no size; state the ",
    "forest's 'sd'"
  )
)
expect_error(
  dbarts(
    y ~ forest(x1, amplitude = normal()),
    d,
    n.trees = 10L,
    n.cuts = 20L,
    control = ctl()
  ),
  "this model has one forest",
  fixed = TRUE
)
# live, normal() to student(3) is taken as creation takes it
live <- make(
  list(forest(sd = 1, amplitude = normalPrior()), forest(basis = zf)),
  response = y
)
live$setForest(1, amplitude = studentPrior(3))
expect_identical(
  live$getForest(1),
  make(
    list(forest(sd = 1, amplitude = studentPrior(3)), forest(basis = zf)),
    response = y
  )$getForest(1)
)

# the multinomial entry, whose category forests are lone ones
expect_identical(
  refusal(dbarts(
    x,
    factor(rep(c("a", "b", "c"), length.out = n)),
    family = "multinomial",
    forests = list(forest(amplitude = normalPrior())),
    n.trees = 10L,
    n.cuts = 20L,
    control = ctl()
  )),
  lone
)
# a variance forest refuses the whole model, in its own words
expect_identical(
  refusal(dbarts(
    y ~ forest(x1, sd = 1, amplitude = normal()) + forest(x2, basis = zf),
    d,
    variance = ~x1,
    n.trees = 10L,
    n.cuts = 20L,
    control = ctl()
  )),
  "a model of several forests does not support 'variance'; drop it or fit a single-forest model"
)
# --- a caller's own normal: in sd = and the other number arguments a function
# of that name is called; where an argument takes a prior it is the package's ---
local({
  normal <- function(v) 2 * v
  fit <- dbarts(
    y ~ forest(x1, sd = normal(0.5)) + forest(x2, basis = zf),
    d,
    n.trees = 10L,
    n.cuts = 20L,
    control = ctl()
  )
  expect_equal(fit$model@forests[[1L]]@sd, 1)
  fit <- dbarts(
    y ~ x1 + x2,
    d,
    forests = list(forest(sd = normal(0.5)), forest(basis = zf)),
    n.trees = 10L,
    n.cuts = 20L,
    control = ctl()
  )
  expect_equal(fit$model@forests[[1L]]@sd, 1)
})
local({
  expect_error(
    dbarts(
      y ~ x1 + x2,
      d,
      forests = list(forest(sd = normal(0.5)), forest(basis = zf)),
      n.trees = 10L,
      n.cuts = 20L,
      control = ctl()
    ),
    "forest 'sd' must be a number, not a dbartsNormalPrior",
    fixed = TRUE
  )
})
local({
  normal <- function() dbartsPriors$fixed()
  sampler <- dbarts(
    y ~ x1 + x2,
    d,
    forests = list(forest(sd = 0.7, amplitude = normal()), forest(basis = zf)),
    n.trees = 10L,
    n.cuts = 20L,
    control = ctl()
  )
  expect_true(is(sampler$getLeafPrior(1L)$amplitude, "dbartsNormalPrior"))
})
local({
  dn <- d
  dn$normal <- dose
  normal <- dose
  sampler <- dbarts(
    y ~ forest(x1 + normal) + forest(x2, basis = normal),
    dn,
    n.trees = 10L,
    n.cuts = 20L,
    control = ctl()
  )
  expect_identical(colnames(sampler$data@bases[[2L]]), "normal")
  expect_true("normal" %in% colnames(sampler$data@x))
})

# --- a constant column names the spelling that draws, and that spelling is taken
# where the call's advice meets the plain-terms refusal ---
one <- rep(1, n)
constantRefusal <- refusal(dbarts(
  y ~ forest(x1) + forest(x2, basis = one, sd = 2),
  d,
  n.trees = 10L,
  n.cuts = 20L,
  control = ctl()
))
expect_true(
  grepl(
    "for a drawn coefficient, forest(x2, sd = 2, amplitude = normal())",
    constantRefusal,
    fixed = TRUE
  ),
  info = constantRefusal
)
written <- dbarts(
  y ~ forest(x1) + forest(x2, sd = 2, amplitude = normal()),
  d,
  n.trees = 10L,
  n.cuts = 20L,
  control = ctl()
)
expect_true(is(written$getLeafPrior(2L)$amplitude, "dbartsNormalPrior"))
# with plain terms the clause's advice meets their own refusal, which names the
# spelling that is taken
plainTerms <- refusal(dbarts(
  y ~ x1 + forest(x2, sd = 2, amplitude = normal()),
  d,
  n.trees = 10L,
  n.cuts = 20L,
  control = ctl()
))
expect_true(
  grepl(
    "forest(x1) + forest(x2, sd = 2, amplitude = normal())",
    plainTerms,
    fixed = TRUE
  ),
  info = plainTerms
)
