# A numeric multiplier in other units is the same forest. A row leaves a
# forest's update when its multiplier is at most 2^-26 of the largest that
# forest has in the sweep, so a column of values near 1e-9 is not lost to an
# absolute bound: a fit on dose and on 2^-40 dose draws the same, to the bit.

forest <- dbartsForests$forest

set.seed(20261018)
n <- 120L
x1 <- runif(n)
x2 <- runif(n)
dose <- rnorm(n, 0.05, 0.02)
y <- 2 + sin(3 * x1) + 25 * dose * (1 + x2) + rnorm(n, sd = 0.5)
d <- data.frame(y, x1, x2)

twinControl <- function() {
  dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.samples = 25L,
    n.burn = 0L,
    updateState = FALSE,
    seed = 47L
  )
}
twinDraws <- function(basis) {
  sampler <- dbarts(
    y ~ x1 + x2,
    d,
    forests = list(forest(), forest(basis = basis)),
    control = twinControl()
  )
  samples <- sampler$run()
  list(train = samples$train, sigma = samples$sigma)
}

# --- the units twin, a drawn coefficient ---
expect_identical(
  twinDraws(dose),
  twinDraws(2^-40 * dose),
  info = "dose and 2^-40 dose draw the same"
)
