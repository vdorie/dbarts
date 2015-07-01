generateHillData <- function() {
  f0 <- function(x) 90 + exp(0.06 * x)
  f1 <- function(x) 72 + 3 * sqrt(x)
  
  set.seed(2793)
  
  ## generate true values
  n <- 120
  p.0 <- 0.5
  
  z.0 <- rbinom(n, 1, p.0)
  n1 <- sum(z.0); n0 <- n - n1
  
  x <- numeric(n)
  x[z.0 == 0] <- rnorm(n0, 20, 10)
  x[z.0 == 1] <- rnorm(n1, 40, 10)
  y <- numeric(n)
  y[z.0 == 0] <- f0(x[z.0 == 0])
  y[z.0 == 1] <- f1(x[z.0 == 1])
  y <- y + rnorm(n)

  p <- rbeta(1, 1, 1)
  z <- rbinom(n, 1, p)
  
  return(list(x = x, y = y, n = n, z.0 = z.0, p.0 = p.0, p = p, z = z))
}
testData <- generateHillData()
rm(generateHillData)
