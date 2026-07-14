source(system.file("common", "hillData.R", package = "dbarts"), local = TRUE)

control <- dbarts::dbartsControl(
  n.burn = 0L,
  n.samples = 1L,
  n.thin = 5L,
  n.chains = 1L,
  n.threads = 1L,
  updateState = FALSE,
  verbose = FALSE,
)

# test that dbarts sampler updates predictors correctly
train <- data.frame(y = testData$y, x = testData$x, z = testData$z)
test <- data.frame(x = testData$x, z = 1 - testData$z)

sampler <- dbarts::dbarts(y ~ x + z, train, test, control = control)

n <- testData$n
z <- testData$z

sampler$setOffset(numeric(n))
expect_equal(sampler$data@offset, numeric(n))
sampler$setOffset(NULL)
expect_null(sampler$data@offset)

invisible(sampler$setPredictor(numeric(n), 2L))
expect_equal(as.numeric(sampler$data@x[, 2L]), numeric(n))

invisible(sampler$setPredictor(x = z, column = 2L))
expect_equal(as.numeric(sampler$data@x[, 2L]), z)

sampler$setTestPredictor(x = 1 - z, column = 2L)
expect_equal(as.numeric(sampler$data@x.test[, 2L]), 1 - z)

sampler$setTestPredictor(x = z, column = "z")
expect_equal(as.numeric(sampler$data@x.test[, "z"]), z)

sampler$setTestPredictor(NULL)
expect_null(sampler$data@x.test)

sampler$setTestPredictor(test)
expect_equal(sampler$data@x.test, as.matrix(test))

set.seed(0L)
new.x <- rnorm(n)
new.z <- as.double(rbinom(n, 1L, 0.5))
new.data <- cbind(new.x, new.z)
invisible(sampler$setPredictor(new.data))

expect_equal(as.numeric(sampler$data@x), as.numeric(new.data))

rm(new.data, new.z, new.x, z, n, sampler)


# test that dbarts sampler shallow/deep copies
## train, test defined above
sampler <- dbarts::dbarts(y ~ x + z, train, test, control = control)

shallowCopy <- sampler$copy(shallow = TRUE)

n <- testData$n

# a shallow copy shares the creation-time predictor object, but the engine
# keeps no matrix to write through: a mutation on one sampler is maintained
# R-side (copy-on-write) and no longer propagates to the copy (design plan 1)
x.shared <- shallowCopy$data@x
invisible(sampler$setPredictor(numeric(n), 2L))
expect_equal(shallowCopy$data@x, x.shared)
expect_equal(as.numeric(sampler$data@x[, 2L]), numeric(n))

# extract() is the on-demand materializer of the current data@x, on both
# the mutated original and the copy it no longer follows (data@x contract,
# design/data-ownership.md plan 3)
expect_equal(extract(sampler, "predictors"), as.matrix(sampler$data@x))
expect_equal(extract(shallowCopy, "predictors"), as.matrix(x.shared))
expect_false(isTRUE(all.equal(
  extract(sampler, "predictors"),
  extract(shallowCopy, "predictors")
)))

rm(shallowCopy)
gc(verbose = FALSE)

deepCopy <- sampler$copy(shallow = FALSE)

invisible(sampler$setPredictor(1 - train$z, 2L))
expect_false(all(sampler$data@x[, 2L] == deepCopy$data@x[, 2L]))

invisible(sampler$setPredictor(deepCopy$data@x[, 2L], 2L))
expect_equal(sampler$data@x, deepCopy$data@x)

rm(deepCopy, n, sampler)


# extract() after a per-observation partial mutation, and a shallow copy
# taken before it: extends the full-column case above to the partial path
# (does not duplicate it - the divergence there is already covered)
sampler <- dbarts::dbarts(y ~ x + z, train, test, control = control)
n <- testData$n

partialCopy <- sampler$copy(shallow = TRUE)
origZ <- as.numeric(sampler$data@x[, 2L])

installed <- sampler$setPredictor(1 - origZ, 2L, forceUpdate = "partial")

predictors <- extract(sampler, "predictors")
expect_equal(predictors, as.matrix(sampler$data@x))
expect_equal(predictors[installed, 2L], (1 - origZ)[installed])
expect_equal(predictors[!installed, 2L], origZ[!installed])

expect_equal(extract(partialCopy, "predictors")[, 2L], origZ)
expect_false(isTRUE(all.equal(
  extract(sampler, "predictors"),
  extract(partialCopy, "predictors")
)))

rm(sampler, n, partialCopy, origZ, installed, predictors)
gc(verbose = FALSE)


# a deep copy of a sampler with a stored state preserves the fitted trees and is
# fully independent of the source (regression: copy() referenced removed state slots)
stateControl <- dbarts::dbartsControl(
  n.burn = 0L,
  n.samples = 1L,
  n.thin = 5L,
  n.chains = 2L,
  n.threads = 1L,
  updateState = TRUE,
  verbose = FALSE
)
sampler <- dbarts::dbarts(y ~ x + z, train, control = stateControl)
invisible(sampler$run(10L, 1L))

stateCopy <- sampler$copy(shallow = FALSE)
expect_equal(sampler$getTrees(), stateCopy$getTrees()) # trees carried over

treesBefore <- sampler$getTrees()
invisible(stateCopy$run(10L, 1L)) # mutate only the copy
expect_equal(sampler$getTrees(), treesBefore) # source unaffected
expect_false(isTRUE(all.equal(stateCopy$getTrees(), treesBefore))) # copy diverged

rm(stateControl, sampler, stateCopy, treesBefore)


# the same, with keepTrees = TRUE, to exercise the savedTrees deep-copy branch
keepControl <- dbarts::dbartsControl(
  n.burn = 0L,
  n.samples = 3L,
  n.thin = 1L,
  n.chains = 1L,
  n.threads = 1L,
  keepTrees = TRUE,
  updateState = TRUE,
  verbose = FALSE
)
sampler <- dbarts::dbarts(y ~ x + z, train, control = keepControl)
invisible(sampler$run(5L, 3L))

keepCopy <- sampler$copy(shallow = FALSE)
expect_equal(sampler$getTrees(), keepCopy$getTrees()) # saved trees carried over

rm(keepControl, sampler, keepCopy)


rm(test, train)


# test that setPredictor with matrix specification doesn't change variables in parent frame
x.train <- dbarts::makeModelMatrixFromDataFrame(data.frame(
  x = testData$x,
  z = testData$z
))
y.train <- testData$y

sampler <- dbarts::dbarts(x.train, y.train, control = control)

n <- testData$n
z <- testData$z

invisible(sampler$setPredictor(numeric(n), 2L))
expect_equal(as.numeric(sampler$data@x[, 2L]), numeric(n))
expect_equal(as.numeric(x.train[, 2L]), z)

rm(z, n, sampler, y.train, x.train)


rm(testData)
