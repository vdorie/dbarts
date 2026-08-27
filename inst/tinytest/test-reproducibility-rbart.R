# Seeded-drift tripwire, not a correctness test: pins draws from a fixed
# seed so unintended RNG-affecting changes are caught. After an intentional
# shift, run tools/regenerate-snapshots.R and eyeball that the new values
# move by a plausible magnitude.

source(
  system.file("common", "rbartGroupData.R", package = "dbarts"),
  local = TRUE
)
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
    1.00660488039067,
    -0.198294267268686,
    -0.785368620164461,
    0.830880629071865,
    2.51966063263758
  )
)

expect_equal(as.numeric(rbartFit$ranef), reference$ranef)

rm(rbartFit, df, testData)
