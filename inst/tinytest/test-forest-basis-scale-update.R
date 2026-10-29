# $setForestBasis(updateBasisScale = TRUE) takes a numeric forest's default
# scale again from the column it installs, over the rows in the likelihood at
# the call, and restates that forest's default sd from it, as a sampler created
# on that column over those rows would hold it; the trees and the coefficients
# are left as they are, and the record on the data carries the new scale to
# every re-creation. Too few rows in, or a column constant over them, installs
# the column and keeps the scale, with one warning per sampler (dec-B504). A
# forest that states its sd, or is multiplied by a factor, has no default taken
# from a column and is refused.

forest <- dbartsForests$forest

set.seed(20261029)
n <- 80L
x1 <- runif(n)
x2 <- runif(n)
z <- factor(rbinom(n, 1L, 0.5))
dose <- rnorm(n, 5, 2)
w <- rep(1, n)
out <- seq(3L, n, by = 5L)
w[out] <- 0
# the rows out hold a dose far from the rest
far <- ifelse(w > 0, dose, 500 + dose)
y <- 1 + 2 * sin(3 * x1) + 0.3 * dose * x2 + rnorm(n, sd = 0.5)
yb <- as.double(y > median(y))
x <- cbind(x1 = x1, x2 = x2)

updateControl <- function() {
  dbartsControl(
    n.chains = 2L,
    n.threads = 1L,
    n.samples = 5L,
    n.burn = 0L,
    updateState = FALSE,
    seed = 61L
  )
}
threeForests <- function(column, weights = NULL, family = "gaussian", ...) {
  dbarts(
    x,
    if (family == "gaussian") y else yb,
    weights = weights,
    family = family,
    forests = list(forest(), forest(basis = z), forest(basis = column, ...)),
    control = updateControl()
  )
}
recordOf <- function(sampler, f = 3L) {
  dbarts:::dataBasisScale(sampler$data)[[f]]
}
priorOf <- function(sampler, f = 3L) {
  prior <- sampler$getLeafPrior(f)
  list(sd = unname(prior$sd), basis.scale = unname(prior$basis.scale))
}
unitOf <- function(sampler) {
  sampler$data@forest.calibration[["unit"]]
}
keptWarnings <- 0L
counted <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    if (inherits(w, "dbartsScaleFallbackWarning")) {
      keptWarnings <<- keptWarnings + 1L
    }
    invokeRestart("muffleWarning")
  })
}

# --- the twin: rescaled before any sweep, the sampler created on the new
# column ---
swapped <- threeForests(dose)
twin <- threeForests(3 * dose)
swapped$setForestBasis(3L, 3 * dose, updateBasisScale = TRUE)
expect_identical(swapped$getLeafPrior(3L), twin$getLeafPrior(3L))
expect_identical(swapped$getLeafPrior(), twin$getLeafPrior())
expect_identical(recordOf(swapped), recordOf(twin))
expect_identical(swapped$run(0L, 5L)$train, twin$run(0L, 5L)$train)

# --- a swap without the argument keeps the scale, as before ---
plain <- threeForests(dose)
created <- priorOf(plain)
plain$setForestBasis(3L, 3 * dose)
expect_identical(priorOf(plain), created)
expect_identical(recordOf(plain), created$basis.scale)
plain$setForestBasis(3L, 3 * dose, updateBasisScale = FALSE)
expect_identical(priorOf(plain), created)

# --- after a run: the default over the new column, the trees and the
# coefficient as they were ---
running <- threeForests(dose)
invisible(running$run(0L, 5L))
fits <- running$getForestFits()
amplitudes <- running$getForestAmplitudes()
running$setForestBasis(3L, dose / 4, updateBasisScale = TRUE)
expect_identical(running$getForestFits(), fits)
expect_identical(running$getForestAmplitudes(), amplitudes)
expect_equal(recordOf(running), sd(dose / 4), tolerance = 1e-12)
expect_identical(priorOf(running)$basis.scale, recordOf(running))
expect_equal(
  priorOf(running)$sd,
  unitOf(running) / sd(dose / 4),
  tolerance = 1e-12
)
expect_false(running$getLeafPrior(3L)$sd.stated)
expect_identical(running$getLeafPrior(2L), plain$getLeafPrior(2L))

# --- the rows read: those in at the call ---
# created with rows of weight zero, then every row in: the scale over every
# row, where creation read the rows in
weighted <- threeForests(far, weights = w)
expect_equal(recordOf(weighted), sd(far[w > 0]), tolerance = 1e-12)
weighted$setWeights(rep(1, n))
weighted$setForestBasis(3L, far, updateBasisScale = TRUE)
expect_equal(recordOf(weighted), sd(far), tolerance = 1e-12, info = "weights")
# and the other way: created with every row in, rows of weight zero set
# later are not read
weightedBack <- threeForests(far)
weightedBack$setWeights(w)
weightedBack$setForestBasis(3L, far, updateBasisScale = TRUE)
expect_equal(
  recordOf(weightedBack),
  sd(far[w > 0]),
  tolerance = 1e-12,
  info = "weights back"
)
# under probit and logistic, a mask installed after creation
for (family in c("probit", "logistic")) {
  masked <- threeForests(far, family = family)
  expect_equal(recordOf(masked), sd(far), tolerance = 1e-12, info = family)
  masked$setActiveRows(w > 0)
  masked$setForestBasis(3L, far, updateBasisScale = TRUE)
  expect_equal(
    recordOf(masked),
    sd(far[w > 0]),
    tolerance = 1e-12,
    info = paste(family, "mask")
  )
}
# a forest's own weights are not read
ownWeights <- threeForests(dose)
ownWeights$setForestWeights(3L, rep(c(0, 1), length.out = n))
ownWeights$setForestBasis(3L, far, updateBasisScale = TRUE)
expect_equal(
  recordOf(ownWeights),
  sd(far),
  tolerance = 1e-12,
  info = "forest weights"
)

# --- a copy and a reload run the new prior ---
# after sweeps, the prior a copy and a reload hold
updated <- threeForests(dose)
invisible(updated$run(0L, 3L))
updated$setForestBasis(3L, 2 * dose + 1, updateBasisScale = TRUE)
updated$storeState()
file <- tempfile(fileext = ".rds")
saveRDS(updated, file)
reloaded <- readRDS(file)
unlink(file)
for (path in list(copy = updated$copy(), reload = reloaded)) {
  expect_identical(path$getLeafPrior(), updated$getLeafPrior())
  expect_identical(recordOf(path), recordOf(updated))
}
# and their draws, from an update before any sweep, where a re-creation holds
# the live sampler's fits to the bit
fresh <- threeForests(dose)
fresh$setForestBasis(3L, 2 * dose + 1, updateBasisScale = TRUE)
freshCopy <- fresh$copy()
fresh$storeState()
file <- tempfile(fileext = ".rds")
saveRDS(fresh, file)
freshReloaded <- readRDS(file)
unlink(file)
live <- fresh$run(0L, 5L)
expect_identical(freshCopy$run(0L, 5L), live, info = "copy")
expect_identical(freshReloaded$run(0L, 5L), live, info = "reload")

# --- refused by name, nothing moved ---
refusedOn <- function(sampler, f, basis, message, flag = TRUE) {
  before <- list(
    data = sampler$data,
    prior = sampler$getLeafPrior(),
    fits = sampler$getForestFits()
  )
  numErrors <- 0L
  keptWarnings <<- 0L
  counted(tryCatch(
    sampler$setForestBasis(f, basis, updateBasisScale = flag),
    error = function(e) {
      numErrors <<- numErrors + 1L
      expect_identical(conditionMessage(e), message)
    }
  ))
  expect_identical(numErrors, 1L)
  expect_identical(keptWarnings, 0L)
  expect_identical(sampler$data, before$data)
  expect_identical(sampler$getLeafPrior(), before$prior)
  expect_identical(sampler$getForestFits(), before$fits)
}
stated <- threeForests(dose, sd = 0.4)
invisible(stated$run(0L, 2L))
refusedOn(
  stated,
  3L,
  3 * dose,
  paste0(
    "forest 3 ('column') states its sd, 0.4, so it has no default to ",
    "update; restate ",
    "it with $setLeafPrior(forests = ) or leave updateBasisScale FALSE"
  )
)
statedLater <- threeForests(dose)
invisible(statedLater$run(0L, 2L))
statedLater$setLeafPrior(
  forests = list(forest(), forest(), forest(sd = 0.25))
)
refusedOn(
  statedLater,
  3L,
  3 * dose,
  paste0(
    "forest 3 ('column') states its sd, 0.25, so it has no default to ",
    "update; restate ",
    "it with $setLeafPrior(forests = ) or leave updateBasisScale FALSE"
  )
)
refusedOn(
  statedLater,
  2L,
  z,
  paste0(
    "forest 2 ('z') is multiplied by a factor, so its default is one unit and ",
    "reads no column; leave updateBasisScale FALSE"
  )
)
for (bad in list(NA, c(TRUE, TRUE), 1, "yes", NULL)) {
  refusedOn(
    running,
    3L,
    dose,
    "'updateBasisScale' must be TRUE or FALSE",
    flag = bad
  )
}
# a label names the forest in the refusal
labeled <- dbarts(
  y ~ forest(x1 + x2) + forest(x1, basis = dose, sd = 0.4),
  data.frame(y, x1, x2, dose),
  control = updateControl()
)
expect_error(
  labeled$setForestBasis(2L, dose, updateBasisScale = TRUE),
  "forest 2 ('dose') states its sd, 0.4, so it has no default",
  fixed = TRUE
)

# --- kept: too few rows in, a column constant over them; one warning per
# sampler, the column installed ---
keptText <- function(expr) {
  caught <- character()
  withCallingHandlers(
    expr,
    dbartsScaleFallbackWarning = function(w) {
      caught <<- c(caught, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  caught
}
oneIn <- threeForests(dose)
invisible(oneIn$run(0L, 2L))
before <- priorOf(oneIn)
oneIn$setWeights(c(1, rep(0, n - 1L)))
expect_identical(
  keptText(oneIn$setForestBasis(3L, far, updateBasisScale = TRUE)),
  paste0(
    "fewer than two observations are in the likelihood (1 of 80), so the ",
    "scale of forest 3 ('column')'s basis is kept as it was (said once per ",
    "sampler)"
  )
)
expect_identical(priorOf(oneIn), before)
expect_identical(recordOf(oneIn), before$basis.scale)
expect_identical(oneIn$data@bases[[3L]][, 1L], far)
# said once: a second keep, of either kind, is silent
expect_identical(
  keptText(oneIn$setForestBasis(3L, dose, updateBasisScale = TRUE)),
  character()
)
oneIn$setWeights(w)
expect_identical(
  keptText(oneIn$setForestBasis(
    3L,
    ifelse(w > 0, 2, dose),
    updateBasisScale = TRUE
  )),
  character()
)
expect_identical(priorOf(oneIn), before)
# no row in
noneIn <- threeForests(dose)
noneIn$setWeights(rep(0, n))
expect_match(
  keptText(noneIn$setForestBasis(3L, far, updateBasisScale = TRUE)),
  "(0 of 80), so the scale of forest 3 ('column')'s basis is kept",
  fixed = TRUE
)
# constant over the rows in, varying over the others
constant <- threeForests(dose, weights = w)
before <- priorOf(constant)
flat <- ifelse(w > 0, 2, dose)
expect_identical(
  keptText(constant$setForestBasis(3L, flat, updateBasisScale = TRUE)),
  paste0(
    "forest 3 ('column')'s basis is constant over the rows in the ",
    "likelihood, so its ",
    "scale is kept as it was (said once per sampler)"
  )
)
expect_identical(priorOf(constant), before)
expect_identical(constant$data@bases[[3L]][, 1L], flat)
# a column of zeros, as a plain swap takes it
zeros <- threeForests(dose)
keptWarnings <- 0L
counted(zeros$setForestBasis(3L, rep(0, n), updateBasisScale = TRUE))
counted(zeros$setForestBasis(3L, rep(0, n), updateBasisScale = TRUE))
expect_identical(keptWarnings, 1L)
# the same column with rows varying in again rescales, silently
keptWarnings <- 0L
counted(zeros$setForestBasis(3L, dose, updateBasisScale = TRUE))
expect_identical(keptWarnings, 0L)
expect_equal(recordOf(zeros), sd(dose), tolerance = 1e-12)
