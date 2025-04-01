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
    0.360083720993859, 0.213898385154795, 0.514888279642085,
    0.402547682652614, 0.0641376173491096
  )
)
expect_null(bartFit$yhat.test)
expect_equal(bartFit$varcount[n.sims,], c(19L, 26L, 24L))

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
    0.157043723005439, 0.649674546901119, 0.392826725618914,
    0.510142912804732, -0.27263185358599
  )
)
expect_null(bartFit$yhat.test)
expect_equal(bartFit$varcount[n.sims,], c(32L, 21L, 21L))
rm(bartFit, n.sims, n.burn)

rm(testData)

