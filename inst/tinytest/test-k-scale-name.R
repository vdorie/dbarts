# The value k is relative to is named k.scale on the leaf-prior reader, on a
# single-forest sampler, a multi-forest sampler and a fit; no entry is named
# anchor, and the sampler's recorded response transform is response.range.

set.seed(71L)
n <- 60L
x <- matrix(rnorm(2L * n), n, 2L)
y <- x[, 1L] + rnorm(n)
labels <- factor(sample(0:2, n, replace = TRUE))

control <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 5L,
  n.samples = 4L,
  seed = 71L
)

hasName <- function(prior) {
  expect_true("k.scale" %in% names(prior))
  expect_false("anchor" %in% names(prior))
}

single <- dbarts(x, y, control = control)
hasName(single$getLeafPrior())
expect_true(is.numeric(single$getLeafPrior()$k.scale))
expect_true(is.null(attr(single$model, "response.anchor", exact = TRUE)))
expect_equal(length(attr(single$model, "response.range", exact = TRUE)), 2L)

multi <- dbarts(x, labels, family = "multinomial", control = control)
for (prior in multi$getLeafPrior()) {
  hasName(prior)
}

fit <- suppressWarnings(suppressMessages(bart(
  x,
  y,
  n.trees = 5L,
  n.samples = 4L,
  n.burn = 2L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)))
hasName(fit$leaf.prior)
