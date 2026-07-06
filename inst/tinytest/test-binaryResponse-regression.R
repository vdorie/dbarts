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
    -0.135254389899302, 0.305296648320596, 0.667004345507013,
    0.627483394905248, -0.152297510506963
  )
)
expect_null(bartFit$yhat.test)
expect_equal(bartFit$varcount[n.sims,], c(17L, 21L, 27L))

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
    0.28314561499688, 0.00876684849456905, 0.471386423757719,
    0.484111685661704, 0.0461271758169942
  )
)
expect_null(bartFit$yhat.test)
expect_equal(bartFit$varcount[n.sims,], c(21L, 27L, 26L))
rm(bartFit, n.sims, n.burn)

rm(testData)

