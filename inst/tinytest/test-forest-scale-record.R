# A numeric multiplier's default is one unit of the response's scale for one
# standard deviation of its column (dec-B291), that standard deviation taken
# once, by the engine, over the rows in the likelihood when the sampler is
# first created (dec-B303), and recorded with the data: every later creation
# is handed the record and takes nothing again, so no swap, copy, reload or
# change of the weights or the mask moves the prior. A column with no standard
# deviation over those rows is refused where a default would need one, and
# created with an sd stated (dec-B307).

forest <- dbartsForests$forest
fixed <- dbartsPriors$fixed

set.seed(20261019)
n <- 150L
x1 <- runif(n)
x2 <- runif(n)
z <- as.double(seq_len(n) %in% sample(n, 60L))
zf <- factor(z)
dose <- rnorm(n, 0.05, 0.02)
age <- rnorm(n, 50, 17)
w <- rep(1, n)
out <- seq(5L, n, length.out = 10L)
w[out] <- 0
# the rows of weight zero hold a dose far from the rest
dose[out] <- 5
y <- 3 +
  2 * sin(3 * x1) +
  z * (1 + x2) +
  20 * dose * (dose < 1) +
  rnorm(n, sd = 0.6)
y[out] <- y[out] + 40
yb <- as.double(y > median(y))
d <- data.frame(y, yb, x1, x2, zf, dose, age, w)
x <- cbind(x1 = x1, x2 = x2)

recordControl <- function(...) {
  dbartsControl(
    n.chains = 2L,
    n.threads = 1L,
    n.samples = 10L,
    n.burn = 0L,
    updateState = FALSE,
    seed = 53L,
    ...
  )
}
unitOf <- function(sampler) {
  attr(sampler$control, "bartcore.forests", exact = TRUE)$unit
}
recordOf <- function(sampler, f = 2L) dbarts:::dataBasisScale(sampler$data)[[f]]
# what each path is held to: the record, the reader's sd and the forest's own
# scale, which is the sd over the half-normal median where it is drawn
readingOf <- function(sampler, f = 2L) {
  prior <- sampler$getLeafPrior(f)
  list(record = recordOf(sampler, f), sd = prior$sd, k.scale = prior$k.scale)
}
expectRule <- function(sampler, scale, f = 2L, info = NULL) {
  prior <- sampler$getLeafPrior(f)
  expect_equal(recordOf(sampler, f), scale, tolerance = 1e-12, info = info)
  expect_equal(unname(prior$basis.scale), scale, tolerance = 1e-12, info = info)
  expect_equal(
    unname(prior$sd),
    unitOf(sampler) / scale,
    tolerance = 1e-12,
    info = info
  )
  expect_equal(
    prior$k.scale * 0.674,
    unname(prior$sd),
    tolerance = 1e-12,
    info = info
  )
}
numWarnings <- 0L
counted <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    if (inherits(w, "dbartsScaleFallbackWarning")) {
      numWarnings <<- numWarnings + 1L
    }
    invokeRestart("muffleWarning")
  })
}

# --- which rows: subset, a dropped missing response, the weights; never a
# mask set after creation ---
keep <- seq_len(n) != 2L
missingY <- y
missingY[7L] <- NA
d$missingY <- missingY
kept <- keep & !is.na(missingY) & w > 0
held <- keep & !is.na(missingY)
for (numForests in 2:3) {
  forests <- if (numForests == 2L) {
    list(forest(), forest(basis = ~dose))
  } else {
    list(forest(), forest(basis = ~dose), forest(basis = ~zf))
  }
  numWarnings <- 0L
  sampler <- counted(dbarts(
    missingY ~ x1 + x2,
    d,
    subset = keep,
    weights = w,
    forests = forests,
    control = recordControl()
  ))
  info <- paste(numForests, "forests")
  expect_identical(numWarnings, 0L, info = info)
  expect_equal(
    recordOf(sampler),
    sd(dose[kept]),
    tolerance = 1e-12,
    info = info
  )
  expect_true(sd(dose[held]) > 10 * sd(dose[kept]), info = info)
  expect_equal(unitOf(sampler), sd(missingY[kept]), tolerance = 1e-12)
  expectRule(sampler, sd(dose[kept]), info = info)
  expect_identical(names(sampler$getLeafPrior(2L)$sd), "dose", info = info)
  before <- readingOf(sampler)
  sampler$setWeights(rep(1, sum(held)))
  expect_identical(readingOf(sampler), before, info = info)
  sampler$setActiveRows(rep(c(TRUE, FALSE), length.out = sum(held)))
  expect_identical(readingOf(sampler), before, info = info)
}
# created with every weight zero: one warning, and every row read
numWarnings <- 0L
none <- counted(dbarts(
  y ~ x1 + x2,
  d,
  weights = rep(0, n),
  forests = list(forest(), forest(basis = ~dose)),
  control = recordControl()
))
expect_identical(numWarnings, 1L)
expect_equal(recordOf(none), sd(dose), tolerance = 1e-12)
# a dose constant over the rows of positive weight, varying over the others:
# refused with no sd stated, created with one
d$flat <- ifelse(w > 0, 0.25, dose)
expect_error(
  dbarts(
    y ~ x1 + x2,
    d,
    weights = w,
    forests = list(forest(), forest(basis = ~flat)),
    control = recordControl()
  ),
  "forest 2 ('flat') states no sd and its basis is constant (0.25 on every row in the likelihood)",
  fixed = TRUE
)
expect_true(inherits(
  dbarts(
    y ~ x1 + x2,
    d,
    weights = w,
    forests = list(forest(), forest(basis = ~flat, sd = 2)),
    control = recordControl()
  ),
  "dbartsSampler"
))

# under probit the rows in are the ones the 0/1 weights keep, which install as
# the sampler's mask
latent <- dbarts(
  yb ~ x1 + x2,
  d,
  weights = w,
  family = "probit",
  forests = list(forest(), forest(basis = ~dose)),
  control = recordControl()
)
expectRule(latent, sd(dose[w > 0]), info = "probit")
expect_identical(latent$activeRows, w)
latentBefore <- readingOf(latent)
latent$setActiveRows(rep(TRUE, n))
expect_identical(readingOf(latent), latentBefore)
numWarnings <- 0L
latentNone <- counted(dbarts(
  yb ~ x1 + x2,
  d,
  weights = rep(0, n),
  family = "probit",
  forests = list(forest(), forest(basis = ~dose)),
  control = recordControl()
))
expect_identical(numWarnings, 1L)
expect_equal(recordOf(latentNone), sd(dose), tolerance = 1e-12)
expect_warning(
  dbarts(
    yb ~ x1 + x2,
    d,
    weights = rep(0, n),
    family = "probit",
    forests = list(forest(), forest(basis = ~dose)),
    control = recordControl()
  ),
  "no observation is in the likelihood (0 of 150 have a positive weight), so the scale of forest 2 ('dose')'s basis is taken from all 150",
  fixed = TRUE
)
expect_error(
  dbarts(
    yb ~ x1 + x2,
    d,
    weights = w,
    family = "probit",
    forests = list(forest(), forest(basis = ~flat)),
    control = recordControl()
  ),
  "states no sd and its basis is constant (0.25 on every row",
  fixed = TRUE
)

# --- once: a swap, a copy, a reload and a sampler made again from its own
# control, model and data take nothing again ---
twoDose <- function(...) {
  dbarts(
    y ~ x1 + x2,
    d,
    weights = w,
    forests = list(forest(), forest(basis = ~dose, ...)),
    control = recordControl()
  )
}
original <- twoDose()
original$run()
created <- readingOf(original)
amplitudes <- original$getForestAmplitudes()
original$setForestBasis(2L, dose / 10)
expect_identical(readingOf(original), created)
expect_identical(original$getForestAmplitudes(), amplitudes)
original$storeState()
file <- tempfile(fileext = ".rds")
saveRDS(original, file)
reloaded <- readRDS(file)
unlink(file)
invisible(reloaded$run(0L, 1L))
again <- new("dbartsSampler", original$control, original$model, original$data)
for (path in list(original$copy(), reloaded, again)) {
  expect_identical(readingOf(path), created)
}
# with the record cleared by hand, a creation takes the column in force
cleared <- original$data
cleared@basis.scale <- NULL
fresh <- new("dbartsSampler", original$control, original$model, cleared)
expect_equal(recordOf(fresh), sd(dose[w > 0]) / 10, tolerance = 1e-12)

# --- other data: a declaration that replaces a forest's basis takes its
# scale from the new column; a forests list and a control carry none ---
bigDose <- 100 * dose
replaced <- dbartsSpec(
  original$data,
  recordControl(),
  forests = list(forest(), forest(basis = bigDose))
)
replacedSampler <- new(
  "dbartsSampler",
  replaced$control,
  replaced$model,
  replaced$data
)
expect_equal(
  recordOf(replacedSampler),
  sd(100 * dose[w > 0]),
  tolerance = 1e-12
)
other <- d
other$dose <- 3 * age
otherFit <- dbarts(
  y ~ x1 + x2,
  other,
  weights = w,
  forests = list(forest(), forest(basis = ~dose)),
  control = original$control
)
expect_equal(recordOf(otherFit), sd(3 * age[w > 0]), tolerance = 1e-12)
handed <- dbartsData(x, y, weights = w, bases = list(NULL, dose))
first <- dbarts(handed, control = recordControl())
second <- dbarts(handed, control = recordControl())
expect_identical(recordOf(first), recordOf(second))
expect_null(dbarts:::dataBasisScale(handed))

# --- carried and wrong: the record is checked where it is handed on, to a
# forest that states no sd, and nowhere else ---
wrongData <- original$data
wrongData@basis.scale <- list(NULL, c(0.1, 0.2))
expect_error(
  new("dbartsSampler", original$control, original$model, wrongData),
  "'basis.scale' entry 2 must hold one positive finite number for each column of that forest's basis: it has 2 and the basis 1 column",
  fixed = TRUE
)
wrongData@basis.scale <- list(NULL, 0)
expect_error(
  new("dbartsSampler", original$control, original$model, wrongData),
  "'basis.scale' entry 2 must hold one positive finite number"
)
statedOriginal <- twoDose(sd = 30)
statedData <- statedOriginal$data
statedData@basis.scale <- list(NULL, c(0.1, 0.2))
expect_true(inherits(
  new(
    "dbartsSampler",
    statedOriginal$control,
    statedOriginal$model,
    statedData
  ),
  "dbartsSampler"
))

# --- stated, then put back: a write leaves the record ---
written <- twoDose()
written$setLeafPrior(forests = list(forest(), forest(sd = 0.5)))
written$storeState()
file <- tempfile(fileext = ".rds")
saveRDS(written, file)
writtenReloaded <- readRDS(file)
unlink(file)
for (path in list(written, written$copy(), writtenReloaded)) {
  prior <- path$getLeafPrior(2L)
  expect_identical(unname(prior$sd), 0.5)
  expect_true(prior$sd.stated)
  expect_identical(recordOf(path), created$record)
}

# --- not a path: new rows, a state and the response setters ---
fit <- bart(
  y ~ forest(x1 + x2) + forest(x1, basis = dose),
  d,
  weights = w,
  n.samples = 8L,
  n.burn = 4L,
  n.chains = 1L,
  n.threads = 1L,
  seed = 59L,
  keepTrees = TRUE,
  verbose = FALSE
)
nNew <- 20L
newRows <- data.frame(
  x1 = runif(nNew),
  x2 = runif(nNew),
  dose = 1000 * rnorm(nNew, 0.05, 0.02)
)
shift <- fit$fit$getLeafPrior(1L)$response.shift
glueForest <- attr(fit$glue, "forest")
source(
  system.file("common", "recombine.R", package = "dbarts"),
  local = TRUE
)
expect_equal(
  predict(fit, newRows),
  recombine(
    predict(fit, newRows, type = "forest"),
    list(NULL, matrix(newRows$dose)),
    nNew
  ),
  tolerance = 1e-12,
  check.attributes = FALSE
)
donor <- dbarts(
  y ~ x1 + x2,
  other,
  weights = w,
  forests = list(forest(), forest(basis = ~dose)),
  control = recordControl()
)
donor$run()
donor$storeState()
recipient <- twoDose()
recipientBefore <- readingOf(recipient)
recipient$setState(donor$state)
expect_identical(readingOf(recipient), recipientBefore)
setters <- twoDose()
settersBefore <- readingOf(setters)
setters$setResponse(y + 1)
setters$setOffset(rep(0.2, n))
setters$setWeights(rep(1, n))
setters$setPredictor(x2 / 2, 2L)
setters$setActiveRows(rep(c(TRUE, FALSE), length.out = n))
expect_identical(readingOf(setters), settersBefore)

# --- a column with no scale: refused at four doors where it states no sd,
# created where it states one ---
constants <- list(
  ones = rep(1, n),
  tenth = rep(0.1, n),
  odd = rep(61.7, n),
  large = rep(1e8 + 0.1, n)
)
d$kept <- seq_len(n) <= 100L
constantDoors <- list(
  term = function(column, ...) {
    d$col <- column
    formula <- if (length(list(...)) == 0L) {
      y ~ forest(x1 + x2) + forest(x1, basis = col)
    } else {
      y ~ forest(x1 + x2) + forest(x1, basis = col, sd = 2)
    }
    dbarts(formula, d, control = recordControl())
  },
  list.formula = function(column, ...) {
    d$col <- column
    dbarts(
      y ~ x1 + x2,
      d,
      forests = list(forest(), forest(basis = ~col, ...)),
      control = recordControl()
    )
  },
  list.matrix = function(column, ...) {
    dbarts(
      x,
      y,
      forests = list(forest(), forest(basis = column, ...)),
      control = recordControl()
    )
  },
  data = function(column, ...) {
    dbarts(
      dbartsData(x, y, bases = list(NULL, column)),
      forests = list(forest(), forest(...)),
      control = recordControl()
    )
  }
)
for (name in names(constants)) {
  column <- constants[[name]]
  value <- format(column[1L], digits = 15L)
  for (door in names(constantDoors)) {
    info <- paste(name, "at", door)
    expect_error(
      constantDoors[[door]](column),
      paste0(
        "states no sd and its basis is constant (",
        value,
        " on every row in the likelihood), so there is no standard ",
        "deviation to take a default from; state one with sd ="
      ),
      fixed = TRUE,
      info = info
    )
    expect_true(
      inherits(constantDoors[[door]](column, sd = 2), "dbartsSampler"),
      info = info
    )
  }
}
# constant on the rows subset keeps, varying over the rest
d$subsetFlat <- ifelse(d$kept, 2.5, dose)
expect_error(
  dbarts(
    y ~ x1 + x2,
    d,
    subset = kept,
    forests = list(forest(), forest(basis = ~subsetFlat)),
    control = recordControl()
  ),
  "forest 2 ('subsetFlat') states no sd and its basis is constant (2.5 on",
  fixed = TRUE
)
# and the same through a creation from a data object whose record is cleared
for (name in names(constants)) {
  edited <- original$data
  edited@bases[[2L]] <- matrix(constants[[name]], ncol = 1L)
  edited@basis.scale <- NULL
  expect_error(
    new("dbartsSampler", original$control, original$model, edited),
    "states no sd and its basis is constant",
    info = name
  )
}
# constant to rounding, and a column far from zero that varies
wiggle <- 61.7 + 1e-9 * sin(seq_len(n))
expect_error(
  constantDoors$list.matrix(wiggle),
  "states no sd and its basis is constant to rounding (its standard deviation is",
  fixed = TRUE
)
far <- 1e7 + seq_len(n)
expect_equal(
  recordOf(constantDoors$list.matrix(far)),
  sd(far),
  tolerance = 1e-12
)

# --- several forests, one prior-predictive door shut by name ---
expect_error(
  samplePriorPredictive(original),
  "samplePriorPredictive does not support a sampler of several forests: their coefficients are not drawn from the prior here",
  fixed = TRUE
)
