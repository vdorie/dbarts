generateAlmostLinearBinaryData <- function() {
  f <- function(x) {
    res <- 0.5 * x[,1] - 0.3 * x[,2] + 0.3 * x[,3]^2 - 0.1 * x[,1] * x[,2]
  }

  set.seed(99)
  sigma <- 1.0
  n     <- 200

  offset <- 0.5
  x <- matrix(rnorm(n * 3), n, 3)
  mu <- 1.25 * f(x) + offset
  eta <- pnorm(mu)
  y <- rbinom(n, 1, eta)

  list(x = x, y = y, mu = mu, eta = eta, offset = offset)
}
testData <- generateAlmostLinearBinaryData()
rm(generateAlmostLinearBinaryData)
