source(
  system.file("common", "friedmanData.R", package = "dbarts"),
  local = TRUE
)

# the port's deprecation warning is pinned in test-rbart-port.R
onceState <- dbarts:::onceWarnState
onceState[["tombstone.rbart_vi"]] <- TRUE
rm(onceState)

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
testData$b <- b
rm(b, sigma.b, g, n.g)

expect_error(
  dbarts::rbart_vi(y ~ x, testData, group.by = g, n.threads = NA_integer_),
  "'n.threads' must be a positive integer, not NA; leave it out"
)

# test that rbart fails with invalid group.by
expect_error(
  dbarts::rbart_vi(y ~ x, testData, group.by = NA, n.threads = 1L),
  "'group.by' must be coercible to factor type"
)
expect_error(
  dbarts::rbart_vi(y ~ x, testData, group.by = not_a_symbol, n.threads = 1L),
  "'group.by' not found"
)
expect_error(
  dbarts::rbart_vi(y ~ x, testData, group.by = testData$g[-1L], n.threads = 1L),
  "group.by' not of length equal to that of data"
)
expect_error(
  dbarts::rbart_vi(y ~ x, testData, group.by = "not a factor", n.threads = 1L),
  "'group.by' not of length equal to that of data"
)

# a NULL power or base is refused: the defaults are 2 and 0.95, not NULL
rbartDefault <- function(...) {
  dbarts::rbart_vi(
    y ~ x,
    testData,
    group.by = g,
    n.trees = 3L,
    n.samples = 2L,
    n.burn = 1L,
    n.chains = 1L,
    n.thin = 1L,
    n.threads = 1L,
    verbose = FALSE,
    seed = 1L,
    ...
  )
}
expect_error(rbartDefault(power = NULL), "'power' cannot be NULL; leave it out")
expect_error(rbartDefault(base = NULL), "'base' cannot be NULL; leave it out")
# predict's and extract's combineChains is TRUE or FALSE; a missing type gets
# the list of types
rbartSmall <- rbartDefault()
expect_error(
  predict(
    rbartSmall,
    testData[1:2, ],
    group.by = testData$g[1:2],
    combineChains = NA
  ),
  "'combineChains' must be TRUE/FALSE",
  fixed = TRUE
)
expect_error(
  extract(rbartSmall, combineChains = "yes"),
  "'combineChains' must be TRUE/FALSE",
  fixed = TRUE
)
expect_error(
  predict(
    rbartSmall,
    testData[1:2, ],
    group.by = testData$g[1:2],
    combineChains = "yes"
  ),
  "'combineChains' must be TRUE/FALSE",
  fixed = TRUE
)
expect_error(
  extract(rbartSmall, combineChains = NA),
  "'combineChains' must be TRUE/FALSE",
  fixed = TRUE
)
# a number is read as as.logical reads it
expect_identical(
  extract(rbartSmall, combineChains = 0),
  extract(rbartSmall, combineChains = FALSE)
)
rbartKept <- rbartDefault(keepTrees = TRUE)
keptRows <- data.frame(x = I(testData$x[1:2, ]))
expect_identical(
  predict(
    rbartKept,
    keptRows,
    group.by = testData$g[1:2],
    combineChains = 1
  ),
  predict(
    rbartKept,
    keptRows,
    group.by = testData$g[1:2],
    combineChains = TRUE
  )
)
rm(rbartKept, keptRows)
for (type in list(NA_character_, character(0))) {
  expect_error(
    extract(rbartSmall, type = type),
    "type must be in '",
    fixed = TRUE
  )
  expect_error(
    fitted(rbartSmall, type = type),
    "type must be in '",
    fixed = TRUE
  )
  expect_error(
    predict(
      rbartSmall,
      testData[1:2, ],
      group.by = testData$g[1:2],
      type = type
    ),
    "type must be in '",
    fixed = TRUE
  )
}
rm(rbartDefault, rbartSmall, type)

rm(testData)
