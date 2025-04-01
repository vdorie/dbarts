source(system.file("common", "friedmanData.R", package = "dbarts"), local = TRUE)

architectureIsArm <- R.version$arch %in% c("aarch64", "arm", "arm64")

# test that basic Friedman example gets same result with various SIMD sets
maxSIMDLevel <- .Call(dbarts:::C_dbarts_getMaxSIMDInstructionSet)

if (maxSIMDLevel > 0L) {
  n.burn <- 0L
  n.sims <- 10L

  .Call(dbarts:::C_dbarts_setSIMDInstructionSet, 0L)
  set.seed(99L)
  bartFit.0 <- dbarts::bart(
    testData$x, testData$y,
    ndpost = n.sims, nskip = n.burn, ntree = 50L,
    verbose = FALSE
  )

  if (architectureIsArm) {
    if (maxSIMDLevel >= 1L) {
      .Call(dbarts:::C_dbarts_setSIMDInstructionSet, 1L) # NEON
      set.seed(99L)
      bartFit <- dbarts::bart(
        testData$x, testData$y,
        ndpost = n.sims, nskip = n.burn,
        ntree = 50L, verbose = FALSE
      )
    
      expect_equal(bartFit.0$yhat.test, bartFit$yhat.test)
    }
  } else {
    if (maxSIMDLevel >= 2L) {
      .Call(dbarts:::C_dbarts_setSIMDInstructionSet, 2L) # SSE2
      set.seed(99L)
      bartFit <- dbarts::bart(
        testData$x, testData$y,
        ndpost = n.sims, nskip = n.burn, ntree = 50L,
        verbose = FALSE
      )
      
      expect_equal(bartFit.0$yhat.test, bartFit$yhat.test)
    }
    if (maxSIMDLevel >= 5L) {
      .Call(dbarts:::C_dbarts_setSIMDInstructionSet, 5L) # SSE4_1
      set.seed(99L)
      bartFit <- dbarts::bart(
        testData$x, testData$y,
        ndpost = n.sims, nskip = n.burn, ntree = 50L,
        verbose = FALSE
      )
      
      expect_equal(bartFit.0$yhat.test, bartFit$yhat.test)
    }
    if (maxSIMDLevel >= 7L) {
      .Call(dbarts:::C_dbarts_setSIMDInstructionSet, 7L) # AVX
      set.seed(99L)
      bartFit <- dbarts::bart(
        testData$x, testData$y,
        ndpost = n.sims, nskip = n.burn, ntree = 50L,
        verbose = FALSE
      )
      
      expect_equal(bartFit.0$yhat.test, bartFit$yhat.test)
    }
    if (FALSE && maxSIMDLevel >= 8L) {
      .Call(dbarts:::C_dbarts_setSIMDInstructionSet, 8L) # AVX2
      set.seed(99L)
      bartFit <- dbarts::bart(
        testData$x, testData$y,
        ndpost = n.sims, nskip = n.burn, ntree = 50L,
        verbose = FALSE
      )
      
      expect_equal(bartFit.0$yhat.test, bartFit$yhat.test)
    }
  }
  .Call(dbarts:::C_dbarts_setSIMDInstructionSet, maxSIMDLevel)

  rm(bartFit, n.sims, n.burn)
}

rm(testData)

source(system.file("common", "multithreadData.R", package = "dbarts"), local = TRUE)

# test that long data gets same result with various SIMD sets

if (maxSIMDLevel > 0L) {
  n.burn <- 0L
  n.sims <- 10L

  .Call(dbarts:::C_dbarts_setSIMDInstructionSet, 0L)

  set.seed(99L)
  bartFit.0 <- dbarts::bart(
    testData$x, testData$y,
    ndpost = n.sims, nskip = n.burn, ntree = 10L,
    verbose = FALSE
  )

  if (architectureIsArm) {
    if (maxSIMDLevel >= 1L) {
      .Call(dbarts:::C_dbarts_setSIMDInstructionSet, 1L) # NEON
      set.seed(99L)
      bartFit <- dbarts::bart(
        testData$x, testData$y,
        ndpost = n.sims, nskip = n.burn, ntree = 10L,
        verbose = FALSE
      )
    
      expect_equal(bartFit.0$yhat.test, bartFit$yhat.test)
    }
  } else {
    if (maxSIMDLevel >= 2L) {
      .Call(dbarts:::C_dbarts_setSIMDInstructionSet, 2L) # SSE2
      set.seed(99L)
      bartFit <- dbarts::bart(
        testData$x, testData$y,
        ndpost = n.sims, nskip = n.burn, ntree = 10L,
        verbose = FALSE
      )
      
      expect_equal(bartFit.0$yhat.test, bartFit$yhat.test)
    }
    if (maxSIMDLevel >= 5L) {
      .Call(dbarts:::C_dbarts_setSIMDInstructionSet, 5L) # SSE4_1
      set.seed(99L)
      bartFit <- dbarts::bart(
        testData$x, testData$y,
        ndpost = n.sims, nskip = n.burn, ntree = 10L,
        verbose = FALSE
      )
      
      expect_equal(bartFit.0$yhat.test, bartFit$yhat.test)
    }
    if (maxSIMDLevel >= 7L) {
      .Call(dbarts:::C_dbarts_setSIMDInstructionSet, 7L) # AVX
      set.seed(99L)
      bartFit <- dbarts::bart(
        testData$x, testData$y,
        ndpost = n.sims, nskip = n.burn, ntree = 10L,
        verbose = FALSE
      )
      
      expect_equal(bartFit.0$yhat.test, bartFit$yhat.test)
    }
    if (maxSIMDLevel >= 8L) {
      .Call(dbarts:::C_dbarts_setSIMDInstructionSet, 8L) # AVX2
      set.seed(99L)
      bartFit <- dbarts::bart(
        testData$x, testData$y,
        ndpost = n.sims, nskip = n.burn, ntree = 10L,
        verbose = FALSE
      )
      
      expect_equal(bartFit.0$yhat.test, bartFit$yhat.test)
    }
  }
  .Call(dbarts:::C_dbarts_setSIMDInstructionSet, maxSIMDLevel)

  rm(bartFit, n.sims, n.burn)
}

rm(maxSIMDLevel, architectureIsArm, testData)


