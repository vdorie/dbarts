source(system.file("common", "hillData.R", package = "dbarts"), local = TRUE)

fitTinyModel <- function(testData, n.samples) {  
  train <- data.frame(y = testData$y, x = testData$x, z = testData$z)
  test <- data.frame(x = testData$x, z = 1 - testData$z)
  
  control <- dbarts::dbartsControl(
    n.burn = 0L, n.samples = 1L, n.thin = 5L,
    n.chains = 1L, n.threads = 1L, n.trees = 5L,
    updateState = FALSE, verbose = FALSE
  )

  sampler <- dbarts::dbarts(y ~ x + z, train, test, control = control)
  
  n <- testData$n
  y <- testData$y
  z <- testData$z
  p <- testData$p
  
  set.seed(0L)
  for (i in seq_len(n.samples)) {
    samples <- sampler$run()
    
    mu0 <- ifelse(z == 0L, samples$train[,1L], samples$test[,1L])
    mu1 <- ifelse(z == 1L, samples$train[,1L], samples$test[,1L])
    
    p0 <- dnorm(y, mu0, samples$sigma[1L]) * (1 - p)
    p1 <- dnorm(y, mu1, samples$sigma[1L]) * p
    p.z <- p1 / (p0 + p1)

    new.z <- rbinom(n, 1L, p.z)
    sampler$setPredictor(x = new.z, column = 2L, forceUpdate = TRUE)
    z <- new.z
    sampler$setTestPredictor(x = 1L - z, column = 2L)

    n1 <- sum(z); n0 <- n - n1
    p <- rbeta(1L, 1L + n0, 1L + n1)
  }
  
  dbarts:::namedList(sampler, z, p, samples)
}

# test that dbarts sampler runs and matches previous versions
tinyModel <- fitTinyModel(testData, 5L)
sampler <- tinyModel$sampler
z <- tinyModel$z
p <- tinyModel$p
samples <- tinyModel$samples

expect_equal(
  z[1L:20L],
  c(0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0
))
expect_equal(p, 0.939263201710153)
expect_equal(
  samples$train[1L:5L],
  c(
    92.2954347050941, 92.2954347050941, 93.2401525161323,
    92.2954347050941, 93.2708029591151
  )
)
expect_equal(
  samples$test[1L:5L],
  c(
    93.6109481580171, 93.6109481580171, 94.5556659690554,
    93.6109481580171, 94.5863164120381
  )
)

tree <- sampler$getTrees(5L)
expect_equal(tree$n, c(120L, 61L, 59L, 54L, 5L))
expect_equal(tree$var, c(1L, -1L, 1L, -1L, -1L))
expect_equal(
  tree$value,
  c(
    34.8679947872694, -0.0436339769593335, 57.3007502897521,
    -0.0220356180367407, 0.0966367066202154
  )
)

rm(tree, samples, p, z, sampler, tinyModel)


# test that sampler sets cut points correctly
sampler <- fitTinyModel(testData, 5L)$sampler

tree.orig <- sampler$getTrees(1L)

sampler$storeState()
sampler$setCutPoints(attr(sampler$state, "cutPoints")[[1L]][1L:87L], 1L)
tree.new <- sampler$getTrees(1L)

expect_equal(nrow(tree.orig), 5L)
expect_equal(nrow(tree.new), 3L)
expect_equal(tree.orig$value[4L], tree.new$value[3L])


sampler <- fitTinyModel(testData, 5L)$sampler

output.orig <- capture.output(sampler$printTrees(1L))

sampler$storeState()
sampler$setCutPoints(attr(sampler$state, "cutPoints")[[1L]][1L:88L], 1L)

output.new <- capture.output(sampler$printTrees(1L))

expect_equal(length(output.orig), length(output.new))
expect_true(grepl("Avail: 01", output.new[5L]))

rm(output.new, output.orig, sampler, tree.new, tree.orig)


# test that sampler set data collapses nodes correctly
sampler <- fitTinyModel(testData, 15L)$sampler

tree.orig <- sampler$getTrees(2L)

x.new <- sampler$data@x

valuesToSwap <- which(
  x.new[,1L] <= tree.orig$value[1L] & x.new[,2L] <= tree.orig$value[2L]
)
targetValues <- which(
  x.new[,1L] > tree.orig$value[1L]
  & x.new[,2L] > tree.orig$value[2L]
)[seq_along(valuesToSwap)]

temp <- x.new[valuesToSwap,2L]
x.new[valuesToSwap,2L] <- x.new[targetValues,2L]
x.new[targetValues,2L] <- temp

data.new <- sampler$data
data.new@x <- x.new
sampler$setData(data.new)

tree.new <- sampler$getTrees(2L)

expect_equal(tree.new$value[2L], tree.orig$value[4L])
rm(tree.new, temp, targetValues, valuesToSwap, x.new, tree.orig, sampler)

rm(fitTinyModel, testData)

