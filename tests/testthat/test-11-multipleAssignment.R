context("multiple assignment")

test_that("multiple assignment works with missing arguments", {
  rh <- c(2, 5)
  massign[a,b] <- rh
  expect_equal(a, 2)
  expect_equal(b, 5)
  rm(a, b)
  
  massign[a,] <- rh
  expect_equal(a, 2)
  expect_error(b)
  rm(a)
  massign[,b] <- rh
  expect_equal(b, 5)
  expect_error(a)
  rm(b)
  
  expect_warning(massign[a = b,] <- rh)
  rm(a)
})

test_that("multiple assignment works with named arguments", {
  rh <- c(a = 2, b = 5)
  massign[a, c = b] <- rh
  expect_equal(a, 2)
  expect_equal(c, 5)
  expect_error(b)
  rm(a, c)
  
  massign[c = a,] <- rh
  expect_equal(c, 2)
  
  expect_error(massign[c = d,] <- rh)
  
  massign[c = a, d = a] <- rh
  expect_equal(c, 2)
  expect_equal(d, 2)
  
  expect_warning(massign[c = a, c = b] <- rh)
  rm(c)
  
  rh <- c(a = 2, a = 5)
  expect_warning(massign[b = a, c] <- rh)
  expect_equal(b, 2)
  expect_equal(c, 5)
  rm(b, c)
})

test_that("unpack works with different orderings, missing names", {
  rh <- list(a = 2, b = 5)
  unpack[a, b] <- rh
  expect_equal(a, 2)
  expect_equal(b, 5)
  rm(a, b)
  
  unpack[b, a] <- rh
  expect_equal(a, 2)
  expect_equal(b, 5)
  rm(a, b)
  
  rh <- c(a = 2, c = 5)
  unpack[a, b] <- rh
  expect_equal(a, 2)
  expect_error(b)
  rm(a)
})

