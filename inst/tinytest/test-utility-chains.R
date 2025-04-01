source(system.file("common", "friedmanData.R", package = "dbarts"), local = TRUE)

# test that combine/uncombine chains and convert from/to bart style works correctly
n.chains <- 3L
n.samples <- 5L
n.obs <- 7L

dbartsSamples <- array(
  seq_len(n.chains * n.samples * n.obs),
  c(n.obs, n.samples, n.chains),
  dimnames = list(as.character(seq_len(n.obs)), NULL, NULL)
)

bartSamples <- dbarts:::convertSamplesFromDbartsToBart(dbartsSamples, n.chains)

expect_equal(dim(bartSamples), c(n.chains, n.samples, n.obs))
expect_equal(
  dimnames(bartSamples)[c(3L, 2L, 1L)],
  dimnames(dbartsSamples)
)

for (k in seq_len(n.chains))
  expect_equal(t(bartSamples[k,,]), dbartsSamples[,,k])

bartSamples.cc <- dbarts:::convertSamplesFromDbartsToBart(
  dbartsSamples, n.chains, combineChains = TRUE
)
expect_equal(dim(bartSamples.cc), c(n.chains * n.samples, n.obs))
expect_equal(colnames(bartSamples.cc), dimnames(dbartsSamples)[[1L]])

for (k in seq_len(n.chains)) {
  expect_equal(
    bartSamples[k,,],
    bartSamples.cc[seq_len(n.samples) + (k - 1L) * n.samples,]
  )
}


expect_equal(
  dbarts:::uncombineChains(bartSamples.cc, n.chains),
  bartSamples
)
expect_equal(dbarts:::combineChains(bartSamples), bartSamples.cc)

expect_equal(
  dbarts:::convertSamplesFromBartsToDbarts(bartSamples, n.chains),
  dbartsSamples
)
expect_equal(
  dbarts:::convertSamplesFromBartsToDbarts(
    bartSamples.cc, n.chains, uncombineChains = TRUE
  ),
  dbartsSamples
)

rm(k, bartSamples.cc, bartSamples, dbartsSamples, n.obs, n.samples, n.chains)

rm(testData)

