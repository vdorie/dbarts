# Seeded-drift tripwire, not a correctness test: pins draws from a fixed
# seed so unintended RNG-affecting changes are caught. After an intentional
# shift, run tools/regenerate-snapshots.R and eyeball that the new values
# move by a plausible magnitude.

source(
  system.file("common", "friedmanData.R", package = "dbarts"),
  local = TRUE
)

n.g <- 5L
if (getRversion() >= "3.6.0") {
  oldSampleKind <- RNGkind()[3L]
  suppressWarnings(RNGkind(sample.kind = "Rounding"))
}
g <- sample(n.g, length(testData$y), replace = TRUE)
if (getRversion() >= "3.6.0") {
  suppressWarnings(RNGkind(sample.kind = oldSampleKind))
  rm(oldSampleKind)
}

sigma.b <- 1.5
b <- rnorm(n.g, 0, sigma.b)

testData$y <- testData$y + b[g]
testData$g <- g
rm(b, sigma.b, g, n.g)

df <- as.data.frame(testData$x)
colnames(df) <- paste0("x_", seq_len(ncol(testData$x)))
df$y <- testData$y
df$g <- testData$g

set.seed(99L)
rbartFit <- dbarts::rbart_vi(
  y ~ . - g,
  df,
  group.by = g,
  n.samples = 1L,
  n.burn = 5L,
  n.thin = 1L,
  n.chains = 1L,
  n.trees = 25L,
  n.threads = 1L,
  verbose = FALSE
)

reference <- list(
  ranef = c(
    -0.774543316213674,
    -0.535798600197457,
    -0.946375382431218,
    -0.383713620127667,
    0.255301267190684
  )
)

expect_equal(as.numeric(rbartFit$ranef), reference$ranef)

rm(rbartFit, df, testData)
