context("simd instructions")

source(system.file("common", "friedmanData.R", package = "dbarts"))

test_that("basic Friedman example gets same result with various SIMD sets", {
  maxSIMDLevel <- .Call(dbarts:::C_dbarts_getMaxSIMDInstructionSet)
  
  n.burn <- 0L
  n.sims <- 10L
  
  .Call(dbarts:::C_dbarts_setSIMDInstructionSet, 0L)
  set.seed(99)
  bartFit.0 <- bart(testData$x, testData$y, ndpost = n.sims, nskip = n.burn, ntree = 50L, verbose = FALSE)
  
  if (maxSIMDLevel >= 2L) {
    .Call(dbarts:::C_dbarts_setSIMDInstructionSet, 2L)
    set.seed(99)
    bartFit <- bart(testData$x, testData$y, ndpost = n.sims, nskip = n.burn, ntree = 50L, verbose = FALSE)
    
    expect_equal(bartFit.0$yhat.test, bartFit$yhat.test)
  }
  if (maxSIMDLevel >= 4L) {
    .Call(dbarts:::C_dbarts_setSIMDInstructionSet, 4L)
    set.seed(99)
    bartFit <- bart(testData$x, testData$y, ndpost = n.sims, nskip = n.burn, ntree = 50L, verbose = FALSE)
    
    expect_equal(bartFit.0$yhat.test, bartFit$yhat.test)
  }
  if (maxSIMDLevel >= 6L) {
    .Call(dbarts:::C_dbarts_setSIMDInstructionSet, 6L)
    set.seed(99)
    bartFit <- bart(testData$x, testData$y, ndpost = n.sims, nskip = n.burn, ntree = 50L, verbose = FALSE)
    
    expect_equal(bartFit.0$yhat.test, bartFit$yhat.test)
  }
  if (maxSIMDLevel >= 7L) {
    .Call(dbarts:::C_dbarts_setSIMDInstructionSet, 7L)
    set.seed(99)
    bartFit <- bart(testData$x, testData$y, ndpost = n.sims, nskip = n.burn, ntree = 50L, verbose = FALSE)
    
    expect_equal(bartFit.0$yhat.test, bartFit$yhat.test)
  }
  .Call(dbarts:::C_dbarts_setSIMDInstructionSet, maxSIMDLevel)
})

source(system.file("common", "multithreadData.R", package = "dbarts"))

test_that("long data gets same result with various SIMD sets", {
  maxSIMDLevel <- .Call(dbarts:::C_dbarts_getMaxSIMDInstructionSet)
  
  n.burn <- 0L
  n.sims <- 10L
  
  .Call(dbarts:::C_dbarts_setSIMDInstructionSet, 0L)
  set.seed(99)
  bartFit.0 <- bart(testData$x, testData$y, ndpost = n.sims, nskip = n.burn, ntree = 50L, verbose = FALSE)
  
  if (maxSIMDLevel >= 2L) {
    .Call(dbarts:::C_dbarts_setSIMDInstructionSet, 2L)
    set.seed(99)
    bartFit <- bart(testData$x, testData$y, ndpost = n.sims, nskip = n.burn, ntree = 50L, verbose = FALSE)
    
    expect_equal(bartFit.0$yhat.test, bartFit$yhat.test)
  }
  if (maxSIMDLevel >= 4L) {
    .Call(dbarts:::C_dbarts_setSIMDInstructionSet, 4L)
    set.seed(99)
    bartFit <- bart(testData$x, testData$y, ndpost = n.sims, nskip = n.burn, ntree = 50L, verbose = FALSE)
    
    expect_equal(bartFit.0$yhat.test, bartFit$yhat.test)
  }
  if (maxSIMDLevel >= 6L) {
    .Call(dbarts:::C_dbarts_setSIMDInstructionSet, 6L)
    set.seed(99)
    bartFit <- bart(testData$x, testData$y, ndpost = n.sims, nskip = n.burn, ntree = 50L, verbose = FALSE)
    
    expect_equal(bartFit.0$yhat.test, bartFit$yhat.test)
  }
  if (maxSIMDLevel >= 7L) {
    .Call(dbarts:::C_dbarts_setSIMDInstructionSet, 7L)
    set.seed(99)
    bartFit <- bart(testData$x, testData$y, ndpost = n.sims, nskip = n.burn, ntree = 50L, verbose = FALSE)
    
    expect_equal(bartFit.0$yhat.test, bartFit$yhat.test)
  }
  .Call(dbarts:::C_dbarts_setSIMDInstructionSet, maxSIMDLevel)
})

