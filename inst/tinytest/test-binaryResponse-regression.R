source(system.file("common", "probitData.R", package = "dbarts"), local = TRUE)

# test that basic probit example passes regression test
n.burn <- 10L
n.sims <- 100L

set.seed(99)
bartFit <- dbarts::bart(
  y.train = testData$Z, x.train = testData$X,
  ntree = 50L, ndpost = n.sims, nskip = n.burn,
  k = 4.5, verbose = FALSE
)

expect_equal(
  bartFit$yhat.train[n.sims, 1L:5L],
  c(
    0.0773149992140284, 0.400502985797696, 0.151936745406016,
    0.362320505596542, -0.130398793675396
  )
)
expect_null(bartFit$yhat.test)
expect_equal(bartFit$varcount[n.sims,], c(25L, 24L, 23L))

expect_equal(extract(bartFit), pnorm(bartFit$yhat.train))
rm(bartFit, n.sims, n.burn)

# test that basic probit example with offset regression test
n.burn <- 10L
n.sims <- 100L

set.seed(99)
bartFit <- dbarts::bart(
  y.train = testData$Z, x.train = testData$X,
  ntree = 50L, ndpost = n.sims, nskip = n.burn,
  k = 4.5, binaryOffset = 0.1, verbose = FALSE
)

n.sims <- nrow(bartFit$yhat.train)

expect_equal(
  bartFit$yhat.train[n.sims, 1L:5L],
  c(
    0.170250482240226, 0.305633294712975, 0.636605516298103,
    0.391682523723046, -0.166274808432018
  )
)
expect_null(bartFit$yhat.test)
expect_equal(bartFit$varcount[n.sims,], c(27L, 28L, 18L))
rm(bartFit, n.sims, n.burn)

rm(testData)

