# A numeric multiplier in other units is the same forest. A row leaves a
# forest's update when its multiplier is at most 2^-26 of the largest that
# forest has in the sweep, so a column of values near 1e-9 is not lost to an
# absolute bound, and a default is one unit of the response's scale for one
# standard deviation of the column, so a default fit on dose, on 1024 dose and
# on 2^-40 dose draws the same, to the bit, drawn and held. And the size a
# default gives, read off the engine's own prior draws, is one unit per
# standard deviation of the column.

forest <- dbartsForests$forest
fixed <- dbartsPriors$fixed

set.seed(20261018)
n <- 120L
x1 <- runif(n)
x2 <- runif(n)
dose <- rnorm(n, 0.05, 0.02)
age <- rnorm(n, 50, 17)
treated <- as.double(seq_len(n) %in% sample(n, 24L))
y <- 2 + sin(3 * x1) + 25 * dose * (1 + x2) + rnorm(n, sd = 0.5)
yb <- as.double(y > median(y))
d <- data.frame(y, yb, x1, x2)

twinControl <- function(...) {
  dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.samples = 25L,
    n.burn = 0L,
    updateState = FALSE,
    seed = 47L,
    ...
  )
}
twinDraws <- function(basis, family = "gaussian", amplitude = NULL) {
  sampler <- dbarts(
    if (family == "gaussian") y ~ x1 + x2 else yb ~ x1 + x2,
    d,
    family = family,
    forests = list(forest(), forest(basis = basis, amplitude = amplitude)),
    control = twinControl()
  )
  samples <- sampler$run()
  list(train = samples$train, sigma = samples$sigma)
}

# --- the units twin, at the default, drawn and held, gaussian and probit ---
for (family in c("gaussian", "probit")) {
  for (amplitude in list(NULL, fixed())) {
    info <- paste(family, if (is.null(amplitude)) "drawn" else "held")
    reference <- twinDraws(dose, family, amplitude)
    expect_identical(
      twinDraws(1024 * dose, family, amplitude),
      reference,
      info = info
    )
    expect_identical(
      twinDraws(2^-40 * dose, family, amplitude),
      reference,
      info = info
    )
  }
}

# --- the size, from the engine's prior draws: a forest's own fit at a row,
# times the median absolute coefficient where it is drawn, is the sd per unit
# of its column, and one unit per standard deviation of it at the default ---
families <- list(
  gaussian = list(response = y, family = "gaussian"),
  probit = list(response = yb, family = "probit"),
  logistic = list(response = yb, family = "logistic")
)
priorSize <- function(sampler, f, row, draws) {
  scale <- sampler$getLeafPrior(f)$response.scale
  values <- vapply(
    seq_len(draws),
    function(draw) {
      sampler$sampleTreesFromPrior(updateState = FALSE)
      sampler$sampleLeafParametersFromPrior(updateState = FALSE)
      scale * sampler$getForestFits(f)[row, 1L]
    },
    0
  )
  sd(values)
}
sizeCheck <- function(
  column,
  family,
  draws,
  tolerance,
  sd = NULL,
  amplitude = NULL
) {
  sampler <- dbarts(
    cbind(x1 = x1, x2 = x2),
    family$response,
    family = family$family,
    forests = list(
      forest(),
      forest(basis = column, sd = sd, amplitude = amplitude),
      forest(basis = factor(treated))
    ),
    control = twinControl(n.trees = 20L)
  )
  unit <- attr(sampler$control, "bartcore.forests", exact = TRUE)$unit
  set.seed(5)
  measured <- priorSize(sampler, 2L, 1L, draws) *
    if (is.null(amplitude)) qnorm(0.75) else 1
  expected <- if (is.null(sd)) unit / stats::sd(column) else sd
  abs(measured / expected - 1) < tolerance
}
for (name in names(families)) {
  expect_true(sizeCheck(dose, families[[name]], 1500L, 0.1), info = name)
}
if (at_home()) {
  columns <- list(
    dose = dose,
    wide = 1000 * dose,
    age = age,
    treated = treated
  )
  for (name in names(families)) {
    for (column in names(columns)) {
      for (amplitude in list(NULL, fixed())) {
        info <- paste(name, column, if (is.null(amplitude)) "drawn" else "held")
        expect_true(
          sizeCheck(
            columns[[column]],
            families[[name]],
            5000L,
            0.06,
            amplitude = amplitude
          ),
          info = info
        )
        stated <- 0.3 / stats::sd(columns[[column]])
        expect_true(
          sizeCheck(
            columns[[column]],
            families[[name]],
            5000L,
            0.06,
            sd = stated,
            amplitude = amplitude
          ),
          info = paste(info, "stated")
        )
      }
    }
  }
}
