source(system.file("common", "friedmanData.R", package = "dbarts"), local = TRUE)

# test that works with custom loss
x <- testData$x
y <- testData$y

n.reps  <- 3L
n.trees <- c(5L, 7L)
k       <- c(1, 2, 4)
power   <- c(1.5, 2)
base    <- c(0.75, 0.8, 0.95)

mad <- function(y.train, y.train.hat, weights) 
  mean(abs(y.train - apply(y.train.hat, 1L, mean)))

xval <- dbarts::xbart(
  x, y, n.samples = 6L, n.burn = c(5L, 3L, 1L), method = "k-fold", n.test = 5,
  n.reps = n.reps,
  n.trees = n.trees,
  k = k,
  power = power,
  base = base, loss = mad,
  n.threads = 1L
)

expect_inherits(xval, "array")
expect_equal(
  dim(xval),
  c(n.reps, length(n.trees), length(k), length(power), length(base))
)
expect_true(!anyNA(xval))
expect_equal(
  dimnames(xval),
  list(
    rep     = NULL,
    n.trees = as.character(n.trees),
    k       = as.character(k),
    power   = as.character(power),
    base    = as.character(base)
  )
)


xval <- dbarts::xbart(
  x, y, n.samples = 6L, n.burn = c(5L, 3L, 1L), method = "k-fold", n.test = 5,
  n.reps = n.reps,
  n.trees = n.trees,
  k = k,
  power = power,
  base = base, loss = mad,
  n.threads = 2L
)

expect_inherits(xval, "array")
expect_equal(
  dim(xval),
  c(n.reps, length(n.trees), length(k), length(power), length(base))
)
expect_true(!anyNA(xval))
expect_equal(
  dimnames(xval),
  list(
    rep     = NULL,
    n.trees = as.character(n.trees),
    k       = as.character(k),
    power   = as.character(power),
    base    = as.character(base)
  )
)

rm(xval, mad, base, power, k, n.trees, n.reps, y, x)

rm(testData)

