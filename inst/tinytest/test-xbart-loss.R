source(
  system.file("common", "friedmanData.R", package = "dbarts"),
  local = TRUE
)
source(
  system.file("common", "checkXvalShape.R", package = "dbarts"),
  local = TRUE
)

# test that works with custom loss
x <- testData$x
y <- testData$y

n.reps <- 3L
n.trees <- c(5L, 7L)
k <- c(1, 2, 4)
power <- c(1.5, 2)
base <- c(0.75, 0.8, 0.95)

mad <- function(y.train, y.train.hat, weights) {
  mean(abs(y.train - apply(y.train.hat, 1L, mean)))
}

xval <- dbarts::xbart(
  x,
  y,
  n.samples = 6L,
  n.burn = c(5L, 3L, 1L),
  method = "k-fold",
  n.test = 5,
  n.reps = n.reps,
  n.trees = n.trees,
  k = k,
  power = power,
  base = base,
  loss = mad,
  n.threads = 1L
)

checkXvalShape(xval, n.reps, n.trees, k, power, base)


xval <- dbarts::xbart(
  x,
  y,
  n.samples = 6L,
  n.burn = c(5L, 3L, 1L),
  method = "k-fold",
  n.test = 5,
  n.reps = n.reps,
  n.trees = n.trees,
  k = k,
  power = power,
  base = base,
  loss = mad,
  n.threads = 2L
)

checkXvalShape(xval, n.reps, n.trees, k, power, base)

# the reported cell IS the loss function's value, averaged over the folds
# rather than summed: a loss that always returns 3 makes every cell 3
constantLoss <- function(y.test, y.test.hat, weights) 3.0
xval <- dbarts::xbart(
  x,
  y,
  n.samples = 6L,
  n.burn = c(5L, 3L, 1L),
  method = "k-fold",
  n.test = 5,
  n.reps = 2L,
  n.trees = n.trees,
  loss = constantLoss,
  n.threads = 1L
)
expect_equal(as.vector(xval), rep_len(3.0, length(xval)))

rm(xval, constantLoss, mad, base, power, k, n.trees, n.reps, y, x)

rm(testData)
