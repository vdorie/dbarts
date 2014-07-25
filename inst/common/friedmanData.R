generateFriedmanData <- function() {
  f <- function(x) {
    10 * sin(pi * x[,1] * x[,2]) + 20 * (x[,3] - 0.5)^2 + 10 * x[,4] + 5 * x[,5]
  }
  
  set.seed(99)
  sigma <- 1.0
  n     <- 100
  
  x <- matrix(runif(n * 10), n, 10)
  y <- rnorm(n, f(x), sigma)
  
  list(x = x, y = y)
}
testData <- generateFriedmanData()
rm(generateFriedmanData)
