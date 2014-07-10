generateProbitData <- function() {
  n <- 800
  beta <- c(0.12, -0.05, 0.3)
  p <- 3
  
  set.seed(0)
  X <- matrix(rnorm(p * n), n, p)
  
  mu <- pnorm(X %*% beta)
  Z <- rbinom(n, 1, mu)

  return(list(X = X, Z = Z, p = mu))
}
testData <- generateProbitData()
rm(generateProbitData)
