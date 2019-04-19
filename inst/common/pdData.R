generatePDData <- function() {
  f <- function(x) 0.5 * x[,1] + 2 * x[,2] * x[,3]
  
  sigma <- 0.2
  n     <- 100
  
  set.seed(27)
  x <- matrix(2 * runif(n * 3) -1, ncol = 3)
  colnames(x) <- c('rob', 'hugh', 'ed')
  
  Ey <- f(x)
  y  <- rnorm(n, Ey, sigma)
  
  list(x = x, y = y, mu = Ey)
}
testData <- generatePDData()
rm(generatePDData)
