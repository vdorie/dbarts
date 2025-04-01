source(system.file("common", "hillData.R", package = "dbarts"), local = TRUE)

control <- dbarts::dbartsControl(
  n.burn = 0L, n.samples = 1L, n.thin = 5L,
  n.chains = 1L, n.threads = 1L,
  updateState = FALSE, verbose = FALSE,
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
expect_equal(as.numeric(sampler$data@x[,2L]), numeric(n))

invisible(sampler$setPredictor(x = z, column = 2L))
expect_equal(as.numeric(sampler$data@x[,2L]), z)

sampler$setTestPredictor(x = 1 - z, column = 2L)
expect_equal(as.numeric(sampler$data@x.test[,2L]), 1 - z)

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

invisible(sampler$setPredictor(numeric(n), 2L))
expect_equal(sampler$data@x, shallowCopy$data@x)

rm(shallowCopy)
gc(verbose = FALSE)

deepCopy <- sampler$copy(shallow = FALSE)

invisible(sampler$setPredictor(1 - train$z, 2L))
expect_false(all(sampler$data@x[,2L] == deepCopy$data@x[,2L]))

invisible(sampler$setPredictor(deepCopy$data@x[,2L], 2L))
expect_equal(sampler$data@x, deepCopy$data@x)

rm(deepCopy, n, sampler)


rm(test, train)



# test that setPredictor with matrix specification doesn't change variables in parent frame
x.train <- dbarts::makeModelMatrixFromDataFrame(data.frame(x = testData$x, z = testData$z))
y.train <- testData$y

sampler <- dbarts::dbarts(x.train, y.train, control = control)
  
n <- testData$n
z <- testData$z
  
invisible(sampler$setPredictor(numeric(n), 2L))
expect_equal(as.numeric(sampler$data@x[,2L]), numeric(n))
expect_equal(as.numeric(x.train[,2L]), z)

rm(z, n, sampler, y.train, x.train)



rm(testData)
