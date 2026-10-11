# run()$k records k per draw whenever k is drawn during that run, however the
# sampler was created: setLeafPrior moves a forest between a fixed k and a
# drawn one, and the next run follows.

set.seed(3L)
n <- 80L
x <- matrix(runif(n * 2L), n, 2L)
y <- 3 * x[, 1L] + rnorm(n)
makeControl <- function() {
  dbartsControl(
    n.chains = 2L,
    n.threads = 1L,
    n.samples = 15L,
    n.burn = 10L,
    seed = 7L,
    n.trees = 20L
  )
}

# a fixed k turned drawn
sampler <- dbarts(
  x,
  y,
  leaf.prior = normal(k = 2),
  control = makeControl()
)
expect_null(sampler$run(0L, 5L)$k)
sampler$setLeafPrior(normal(k = chi(1.5, 2)))
run <- sampler$run(5L, 8L)
expect_identical(dim(run$k), c(8L, 2L))
expect_false(anyNA(run$k))
expect_true(all(run$k > 0))
expect_true(length(unique(run$k[, 1L])) > 1L)
expect_identical(sampler$getK(), run$k[8L, ])

# a drawn k turned fixed records no k: nothing is drawn in that run, and the
# value in force is getK's
sampler$setLeafPrior(normal(k = 3))
expect_null(sampler$run(0L, 6L)$k)
expect_identical(sampler$getK(), c(3, 3))

# and drawn again
sampler$setLeafPrior(normal(k = chi(1.5, 2)))
run <- sampler$run(2L, 4L)
expect_identical(dim(run$k), c(4L, 2L))
expect_identical(sampler$getK(), run$k[4L, ])

# a sampler created with a drawn k, turned fixed, likewise
created <- dbarts(
  x,
  y,
  leaf.prior = normal(k = chi(1.5, 2)),
  control = makeControl()
)
expect_identical(dim(created$run(0L, 4L)$k), c(4L, 2L))
created$setLeafPrior(normal(k = 3))
expect_null(created$run(0L, 4L)$k)
rm(sampler, created, run, x, y, n, makeControl)
