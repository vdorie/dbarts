## large but simple
generateMultithreadData <- function() {
  set.seed(99)
  
  n <- 500001
  x <- matrix(rnorm(n), n, 1)
  y <- sin(x) + rnorm(n)

  return (list(y = y, x = x))
}
testData <- generateMultithreadData()
rm(generateMultithreadData)
