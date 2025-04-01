source(system.file("common", "friedmanData.R", package = "dbarts"), local = TRUE)

# test that posterior predictive distribution samples use correct sigma
n.samples <- 7L
n.chains  <- 2L
n.obs <- length(testData$y)
bartFit <- dbarts::bart(
  testData$x, testData$y, verbose = FALSE,
  ndpost = n.samples, nskip = 0L, nchain = n.chains,
  ntree = 25L, nthread = 1L
)
set.seed(0L)
samples.ppd <- extract(bartFit, type = "ppd")

set.seed(0L)
samples.pm  <- extract(bartFit)
for (i in seq_len(n.obs)) {
  expect_equal(
    samples.pm[,i] + rnorm(n.samples * n.chains, 0, bartFit$sigma),
    samples.ppd[,i]
  )
}

set.seed(0L)
samples.ppd <- extract(bartFit, type = "ppd", combineChains = FALSE)

set.seed(0L)
samples.pm  <- extract(bartFit, combineChains = FALSE)
for (i in seq_len(n.obs)) {
  expect_equal(
    samples.pm[,,i] + matrix(rnorm(n.samples * n.chains, 0, bartFit$sigma), nrow = n.chains),
    samples.ppd[,,i]
  )
}

rm(i, samples.pm, samples.ppd, bartFit, n.obs, n.chains, n.samples)

rm(testData)

