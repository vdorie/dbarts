# What $getLeafPrior() and extract(type = "leaf.prior.sd") give for each kind
# of forest of a model of several: `sd` in the response's units, a numeric
# forest's per unit of its column and named by it, `sd.stated`, `multiplier`,
# `amplitude` and `prior.sd.of`, and the scale a numeric forest's default was
# taken from as `basis.scale`. Every number is computed here from the fixture
# and held against the number the engine draws under.

forest <- dbartsForests$forest
fixed <- dbartsPriors$fixed

set.seed(20261020)
n <- 150L
x1 <- runif(n)
x2 <- runif(n)
z <- as.double(seq_len(n) %in% sample(n, 60L))
dose <- rnorm(n, 0.05, 0.02)
age <- rnorm(n, 50, 17)
y <- 2 + sin(3 * x1) + z * x2 + 10 * dose + 0.01 * age + rnorm(n, sd = 0.5)
d <- data.frame(y, x1, x2, z, dose, age)
U <- sd(y)

readerControl <- function(...) {
  dbartsControl(
    n.chains = 2L,
    n.threads = 1L,
    n.samples = 10L,
    n.burn = 0L,
    updateState = FALSE,
    seed = 61L,
    ...
  )
}
fourForests <- function() {
  dbarts(
    y ~ forest(x1 + x2) +
      forest(x1, basis = dose) +
      forest(x1, basis = age) +
      forest(x1, basis = factor(z), amplitude = fixed()),
    d,
    control = readerControl()
  )
}
heldTwin <- function() {
  dbarts(
    y ~ forest(x1 + x2, amplitude = fixed()) +
      forest(x1, basis = dose, amplitude = fixed()) +
      forest(x1, basis = age, amplitude = fixed()) +
      forest(x1, basis = factor(z)),
    d,
    control = readerControl()
  )
}
engineMedian <- function(sampler, f) {
  .Call(
    dbarts:::C_dbarts_bartcore_getLeafPrior,
    sampler$getPointer(),
    f - 1L
  )[[1L, "amplitude.prior.scale"]]
}

# --- the literals of the table ---
expected <- list(
  forest1 = list(
    sd = 2 * U,
    multiplier = "none",
    held = FALSE,
    of = "forest (prior median)"
  ),
  dose = list(
    sd = c(dose = U / sd(dose)),
    multiplier = "numeric",
    held = FALSE,
    of = "forest, per unit of basis (prior median)"
  ),
  age = list(
    sd = c(age = U / sd(age)),
    multiplier = "numeric",
    held = FALSE,
    of = "forest, per unit of basis (prior median)"
  ),
  `factor(z)` = list(
    sd = U,
    multiplier = "factor",
    held = TRUE,
    of = "level difference"
  )
)
sampler <- fourForests()
priors <- sampler$getLeafPrior()
expect_identical(names(priors), names(expected))
for (label in names(expected)) {
  prior <- priors[[label]]
  want <- expected[[label]]
  expect_equal(prior$sd, want$sd, tolerance = 1e-12, info = label)
  expect_identical(names(prior$sd), names(want$sd), info = label)
  expect_false(prior$sd.stated, info = label)
  expect_identical(prior$multiplier, want$multiplier, info = label)
  expect_identical(is.null(prior[["amplitude"]]), !want$held, info = label)
  expect_identical(prior$prior.sd.of, want$of, info = label)
  expect_null(prior$amplitude.prior.scale, info = label)
  expect_null(prior$leaf.scale.factor, info = label)
  expect_null(prior$basis.row.norm, info = label)
}
expect_equal(priors$dose$basis.scale, c(dose = sd(dose)), tolerance = 1e-12)
expect_equal(priors$age$basis.scale, c(age = sd(age)), tolerance = 1e-12)
expect_null(priors$forest1$basis.scale)
expect_null(priors[["factor(z)"]]$basis.scale)
# a stated sd reads as stated, per unit of its column, before and after a write
sampler$setLeafPrior(
  forests = list(forest(), forest(sd = 30), forest(), forest())
)
expect_identical(sampler$getLeafPrior(2L)$sd, c(dose = 30))
expect_true(sampler$getLeafPrior(2L)$sd.stated)
expect_identical(
  sampler$getLeafPrior(2L)$leaf.prior,
  forest(sd = 30)
)
expect_false(sampler$getLeafPrior(3L)$sd.stated)
heldStated <- dbarts(
  y ~ forest(x1 + x2) + forest(x1, basis = dose, amplitude = fixed(), sd = 30),
  d,
  control = readerControl()
)
heldPrior <- heldStated$getLeafPrior(2L)
expect_identical(heldPrior$sd, c(dose = 30))
expect_identical(heldPrior$prior.sd.of, "forest, per unit of basis")
expect_true(is(heldPrior$amplitude, "dbartsFixedPrior"))

# --- through the engine: each forest's own scale is its sd, over the
# half-normal median where its coefficient is drawn, and a drawn forest with
# no basis carries its sd as its half-Cauchy median times the unit ---
expectEngine <- function(sampler, info) {
  unit <- sampler$data@forest.calibration[["unit"]]
  for (f in seq_along(sampler$data@bases)) {
    prior <- sampler$getLeafPrior(f)
    where <- paste(info, f)
    drawn <- is.null(prior[["amplitude"]])
    if (drawn && prior$multiplier == "none") {
      expect_equal(
        engineMedian(sampler, f) * unit,
        prior$sd,
        tolerance = 1e-12,
        info = where
      )
    } else {
      expect_equal(
        prior$k.scale * if (drawn) 0.674 else 1,
        unname(prior$sd),
        tolerance = 1e-12,
        info = where
      )
    }
  }
}
for (build in list(default = fourForests, held = heldTwin)) {
  sampler <- build()
  expectEngine(sampler, "creation")
  sampler$setLeafPrior(
    forests = list(
      forest(sd = 1.3),
      forest(sd = 25),
      forest(),
      forest(sd = 0.8)
    )
  )
  expectEngine(sampler, "write")
  sampler$setForestBasis(2L, dose / 10)
  expectEngine(sampler, "swap")
  expectEngine(sampler$copy(), "copy")
  sampler$storeState()
  file <- tempfile(fileext = ".rds")
  saveRDS(sampler, file)
  reloaded <- readRDS(file)
  unlink(file)
  invisible(reloaded$run(0L, 1L))
  expectEngine(reloaded, "reload")
}

# --- read then write: a default numeric forest written back its own number,
# the column's name left on it, is stated and draws as its untouched twin ---
twin <- fourForests()
written <- fourForests()
read <- written$getLeafPrior(2L)
written$setLeafPrior(
  forests = list(forest(), forest(sd = read$sd), forest(), forest())
)
expect_true(written$getLeafPrior(2L)$sd.stated)
expect_identical(written$getLeafPrior(2L)$sd, read$sd)
expect_identical(written$run(0L, 20L), twin$run(0L, 20L))

# --- extract: the list of the sampler's sd entries by label, and one forest's
# entry itself, a numeric forest's with its column's name ---
fit <- bart(
  y ~ forest(x1 + x2) +
    forest(x1, basis = dose) +
    forest(x1, basis = age) +
    forest(x1, basis = factor(z), amplitude = fixed()),
  d,
  n.samples = 10L,
  n.burn = 0L,
  n.chains = 1L,
  n.threads = 1L,
  seed = 67L,
  keepTrees = TRUE,
  verbose = FALSE
)
sds <- lapply(fit$fit$getLeafPrior(), `[[`, "sd")
expect_identical(extract(fit, "leaf.prior.sd"), sds)
expect_identical(extract(fit, "leaf.prior.sd", forest = "dose"), sds$dose)
expect_identical(extract(fit, "leaf.prior.sd", forest = 2L), sds$dose)
expect_identical(names(extract(fit, "leaf.prior.sd", forest = 3L)), "age")
expect_identical(extract(fit, "leaf.prior.sd", forest = 4L), sds[[4L]])
unkept <- bart(
  y ~ forest(x1 + x2) +
    forest(x1, basis = dose) +
    forest(x1, basis = age) +
    forest(x1, basis = factor(z), amplitude = fixed()),
  d,
  n.samples = 10L,
  n.burn = 0L,
  n.chains = 1L,
  n.threads = 1L,
  seed = 67L,
  verbose = FALSE
)
expect_null(unkept$fit)
expect_identical(extract(unkept, "leaf.prior.sd"), sds)
expect_identical(extract(unkept, "leaf.prior.sd", forest = "dose"), sds$dose)
