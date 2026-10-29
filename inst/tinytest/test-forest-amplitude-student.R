# A drawn coefficient under a Student-t prior, amplitude = student(df)
# (dec-B262, dec-B271): with a basis, a t whose median size is the normal
# default's, the forest's sd and leaf scale unmoved; with none, the forest's
# drawn coefficient at the median its sd states, student(1) there being the
# default restated (dec-B509). student() is student(1), and a df is finite and
# at least 1 (dec-B510); the prior is fixed at creation (dec-B511). The name
# student resolves as the prior only at a forest's 'amplitude': here it is
# left unbound, so every door below reads it through its own vocabulary.

forest <- dbartsForests$forest
studentPrior <- dbartsPriors$student

set.seed(20261029)
n <- 80L
x1 <- runif(n)
x2 <- runif(n)
z <- as.double(seq_len(n) %in% sample(n, 35L))
zf <- factor(z)
dose <- rnorm(n, 5, 2)
y <- 1 + sin(3 * x1) + z * (1 + x2) + 0.2 * dose * x2 + rnorm(n, sd = 0.5)
yb <- as.double(y > median(y))
x <- cbind(x1 = x1, x2 = x2)
d <- data.frame(y, yb, x1, x2, zf, dose)

ctl <- function(...) {
  dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.samples = 5L,
    n.burn = 0L,
    updateState = FALSE,
    ...
  )
}

# --- every door, family and kind: the reader says student(3) ---
doors <- list(
  term = function(response, family) {
    formula <- as.formula(paste(
      response,
      "~ forest(x1 + x2) + forest(x1, basis = zf, amplitude = student(3))",
      "+ forest(x2, basis = dose, amplitude = student(3))"
    ))
    dbarts(formula, d, family = family, control = ctl())
  },
  list = function(response, family) {
    dbarts(
      as.formula(paste(response, "~ x1 + x2")),
      d,
      family = family,
      forests = list(
        forest(),
        forest(basis = zf, amplitude = student(3)),
        forest(basis = dose, amplitude = student(3))
      ),
      control = ctl()
    )
  },
  data = function(response, family) {
    dbarts(
      dbartsData(x, d[[response]], bases = list(NULL, zf, dose)),
      family = family,
      forests = list(
        forest(),
        forest(amplitude = student(3)),
        forest(amplitude = student(3))
      ),
      control = ctl()
    )
  },
  spec = function(response, family) {
    spec <- dbartsSpec(
      dbartsData(x, d[[response]]),
      ctl(),
      family = family,
      forests = list(
        forest(),
        forest(basis = zf, amplitude = student(3)),
        forest(basis = dose, amplitude = student(3))
      )
    )
    new("dbartsSampler", spec$control, spec$model, spec$data)
  }
)
families <- list(
  gaussian = list(response = "y", family = "gaussian"),
  probit = list(response = "yb", family = "probit"),
  logistic = list(response = "yb", family = "logistic")
)
for (door in names(doors)) {
  for (family in names(families)) {
    info <- paste(door, family)
    arm <- families[[family]]
    sampler <- doors[[door]](arm$response, arm$family)
    expect_null(sampler$getLeafPrior(1L)[["amplitude"]], info = info)
    for (f in 2:3) {
      prior <- sampler$getLeafPrior(f)
      expect_identical(prior$amplitude, studentPrior(3), info = info)
      expect_identical(
        prior$multiplier,
        c("factor", "numeric")[f - 1L],
        info = info
      )
    }
    expect_identical(
      attr(sampler$control, "bartcore.forests", exact = TRUE)$amplitude.df,
      c(NA, 3, 3),
      info = info
    )
    expect_silent(sampler$run(0L, 2L))
  }
}

# --- the fit, its print and summary; a default fit's record is as it was ---
fit <- bart(
  y ~ forest(x1 + x2) +
    forest(x1, basis = zf, amplitude = student(3)) +
    forest(x2, basis = dose, amplitude = dbartsPriors$fixed()),
  d,
  n.samples = 5L,
  n.burn = 0L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)
expect_identical(fit$leaf.prior[[2L]]$amplitude, studentPrior(3))
expect_identical(fit$leaf.prior[[3L]]$amplitude, dbartsPriors$fixed())
expect_null(fit$leaf.prior[[1L]][["amplitude"]])
expect_identical(fit$fixed$amplitude, c(dose = 1))
full <- summary(fit)$forests[-1L]
expect_true(endsWith(full[[2L]], "(default); coefficient student(3)"))
expect_true(endsWith(full[[3L]], "; coefficient held"))
expect_false(grepl("coefficient", full[[1L]]))
short <- capture.output(print(fit))
expect_false(any(grepl("coefficient", short)))
plain <- dbarts(
  y ~ forest(x1 + x2) + forest(x1, basis = zf),
  d,
  control = ctl()
)
expect_null(
  attr(plain$control, "bartcore.forests", exact = TRUE)[["amplitude.df"]]
)
expect_null(plain$getLeafPrior(2L)[["amplitude"]])

# --- student() and a bare student are student(1) ---
for (written in c("student()", "student")) {
  sampler <- dbarts(
    as.formula(paste(
      "y ~ forest(x1 + x2) + forest(x1, basis = zf, amplitude =",
      written,
      ")"
    )),
    d,
    control = ctl()
  )
  expect_identical(
    sampler$getLeafPrior(2L)$amplitude,
    studentPrior(1),
    info = written
  )
}
expect_identical(studentPrior(), studentPrior(1))
expect_inherits(studentPrior(3L), "dbartsStudentPrior")
expect_identical(studentPrior(3L)@df, 3)

# --- refusals, each by its message and with no warning beside it ---
refused <- function(expr, pattern) {
  warnings <- 0L
  message <- withCallingHandlers(
    tryCatch(
      {
        expr
        NA_character_
      },
      error = conditionMessage
    ),
    warning = function(w) {
      warnings <<- warnings + 1L
      invokeRestart("muffleWarning")
    }
  )
  expect_true(grepl(pattern, message, fixed = TRUE), info = message)
  expect_identical(warnings, 0L, info = pattern)
}
dfRule <- "a coefficient's student(df) takes one finite df of 1 or more"
for (bad in c("0", "0.5", "NA", "NaN", "c(1, 3)", "Inf", "-1", "'3'")) {
  refused(
    dbarts(
      as.formula(paste0(
        "y ~ forest(x1) + forest(x2, basis = zf, amplitude = student(",
        bad,
        "))"
      )),
      d,
      control = ctl()
    ),
    dfRule
  )
}
refused(studentPrior(Inf), dfRule)
refused(
  dbarts(
    y ~ x1,
    d,
    forests = list(
      forest(),
      forest(basis = zf, amplitude = dbartsFamilies$student(3))
    ),
    control = ctl()
  ),
  "dbartsFamilies$student is the response's family"
)
refused(
  dbarts(
    y ~ x1,
    d,
    forests = list(forest(), forest(basis = zf, amplitude = "student")),
    control = ctl()
  ),
  "a forest's 'amplitude' is fixed(), which holds its coefficient, student(df)"
)
refused(
  dbarts(y ~ forest(x1, amplitude = student(3)), d, control = ctl()),
  paste0(
    "which has no coefficient to give a prior; a lone forest draws its size ",
    "through the fitting function's leaf prior, as k = chi(1.5, 2) or ",
    "leaf.prior = normal(sd = invchi(1.5, s)) with s a size in the ",
    "response's units"
  )
)
refused(
  dbarts(y ~ forest(x1, amplitude = dbartsPriors$fixed()), d, control = ctl()),
  "which has no coefficient to hold; 'amplitude' needs a model of several"
)
live <- dbarts(y ~ forest(x1 + x2) + forest(x1, basis = zf), d, control = ctl())
refused(
  live$setLeafPrior(forests = list(forest(), forest(amplitude = student(3)))),
  "'amplitude' is fixed at creation"
)
refused(
  bart(x, y, leaf.prior = student(3), verbose = FALSE),
  "'leaf.prior' must be a leaf prior specification"
)

# --- the family of the same name is still the family ---
expect_identical(
  bart(
    x,
    y,
    family = student(3),
    n.samples = 5L,
    n.burn = 0L,
    n.chains = 1L,
    verbose = FALSE
  )$family,
  "student"
)
expect_identical(
  dbarts(
    y ~ x1 + x2,
    d,
    family = student(3),
    control = ctl()
  )$model@family.spec@token,
  "student"
)
expect_identical(
  dbartsSpec(
    dbartsData(x, y),
    ctl(),
    family = student(3)
  )$model@family.spec@token,
  "student"
)

# --- a forest with no basis: student(1) is the drawn default restated ---
noBasis <- function(amplitude, family = "gaussian", response = "y", sd = NULL) {
  first <- if (is.null(amplitude)) {
    forest(sd = sd)
  } else {
    forest(sd = sd, amplitude = amplitude)
  }
  sampler <- dbarts(
    dbartsData(x, d[[response]], bases = list(NULL, zf)),
    family = family,
    forests = list(first, forest()),
    control = ctl(seed = 17L)
  )
  sampler
}
draws <- function(sampler) sampler$run(10L, 5L)
expect_identical(draws(noBasis(studentPrior(1))), draws(noBasis(NULL)))
expect_identical(
  draws(noBasis(studentPrior(1), "probit", "yb", sd = 1)),
  draws(noBasis(NULL, "probit", "yb"))
)
expect_identical(
  draws(noBasis(studentPrior(), "logistic", "yb")),
  draws(noBasis(NULL, "logistic", "yb"))
)
heavy <- noBasis(studentPrior(3))
expect_identical(heavy$getLeafPrior(1L)$amplitude, studentPrior(3))
expect_identical(heavy$getLeafPrior(1L)$sd, noBasis(NULL)$getLeafPrior(1L)$sd)
expect_false(heavy$getLeafPrior(1L)$sd.stated)
expect_true(all(is.finite(draws(heavy)$glue)))

# --- a copy, a reload and a state install carry the prior and its v ---
# outside the argument the name is unbound
expect_error(forest(amplitude = student(3)), "could not find function")
statefulForests <- list(
  forest(),
  forest(basis = zf, amplitude = studentPrior(3)),
  forest(basis = dose, amplitude = studentPrior(1))
)
stateful <- dbarts(x, y, forests = statefulForests, control = ctl())
set.seed(29L)
stateful$run(30L, 3L)
stateful$storeState()
duplicate <- stateful$copy()
file <- tempfile(fileext = ".rds")
saveRDS(stateful, file)
reloaded <- readRDS(file)
unlink(file)
continued <- lapply(list(stateful, duplicate, reloaded), function(s) {
  set.seed(31L)
  s$run(0L, 3L)
})
expect_identical(continued[[3L]], continued[[2L]])
expect_equal(continued[[2L]], continued[[1L]])
for (s in list(duplicate, reloaded)) {
  expect_identical(s$getLeafPrior(2L)$amplitude, studentPrior(3))
  expect_identical(s$getLeafPrior(3L)$amplitude, studentPrior(1))
}
# a state from a default sampler goes into a student one and runs
default <- dbarts(
  x,
  y,
  forests = list(forest(), forest(basis = zf), forest(basis = dose)),
  control = ctl()
)
default$run(10L, 2L)
default$storeState()
stateful$setState(default$state)
expect_identical(stateful$getForestAmplitudes(), default$getForestAmplitudes())
expect_true(all(is.finite(stateful$run(0L, 3L)$glue)))

# --- a coarse prior check: under a flat likelihood the numeric coefficient's
# size keeps the normal's median and only its tail moves ---
coefficientSizes <- function(amplitude) {
  numeric <- if (is.null(amplitude)) {
    forest(basis = dose, n.trees = 1L)
  } else {
    forest(basis = dose, n.trees = 1L, amplitude = amplitude)
  }
  sampler <- dbarts(
    x,
    y,
    family = gaussian(sigma = dbartsPriors$fixed(1e12)),
    forests = list(
      forest(n.trees = 1L, amplitude = dbartsPriors$fixed()),
      numeric
    ),
    control = ctl(seed = 23L)
  )
  abs(sampler$run(100L, 4000L)$glue[2L, ])
}
for (prior in list(NULL, studentPrior(3))) {
  sizes <- coefficientSizes(prior)
  quantiles <- quantile(sizes, c(0.5, 0.9), names = FALSE)
  expect_true(abs(quantiles[1L] - 0.674) < 0.1, info = quantiles[1L])
  ratio <- quantiles[2L] / quantiles[1L]
  if (is.null(prior)) {
    expect_true(ratio < 2.7, info = ratio)
  } else {
    expect_true(ratio > 2.8, info = ratio)
  }
}
