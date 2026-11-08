# A data frame given as test data to a fit made from a matrix is read as a
# frame fit reads one: every training name among its columns is taken by name,
# extra columns are ignored, a missing or renamed column is refused, and a
# matrix keeps its exact width. The same at predict, at test = and at
# $setTestPredictor.

set.seed(19L)
n <- 60L
df <- data.frame(
  x1 = rnorm(n),
  x2 = rnorm(n),
  g = factor(sample(c("a", "b"), n, replace = TRUE)),
  k = 1
)
df$y <- df$x1 - df$x2 + rnorm(n, 0, 0.3)
train <- as.matrix(df[, c("x1", "x2")])

fitMatrix <- function(...) {
  dbarts::bart(
    train,
    df$y,
    n.trees = 5L,
    n.samples = 10L,
    n.burn = 10L,
    n.chains = 1L,
    n.threads = 1L,
    keepTrees = TRUE,
    seed = 3L,
    verbose = FALSE,
    ...
  )
}
fit <- fitMatrix()
rows <- df[1:7, ]

# extras and a reordering are taken by name
byName <- predict(fit, rows[, c("x1", "x2")])
expect_identical(predict(fit, rows), byName)
expect_identical(predict(fit, rows[, c("k", "x2", "g", "x1")]), byName)
expect_identical(predict(fit, as.matrix(rows[, c("x1", "x2")])), byName)

# at test =
withTest <- fitMatrix(test = rows)
expect_identical(
  withTest$yhat.test,
  fitMatrix(test = rows[, c("x1", "x2")])$yhat.test
)

# at $setTestPredictor
sampler <- dbarts::dbarts(
  train,
  df$y,
  n.trees = 5L,
  control = dbarts::dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.samples = 2L,
    n.burn = 2L,
    verbose = FALSE
  )
)
expect_silent(sampler$setTestPredictor(rows))
expect_equal(
  unname(sampler$data@x.test),
  unname(as.matrix(rows[, c("x1", "x2")]))
)

# a missing or renamed column is still refused, and a matrix keeps its width
expect_error(predict(fit, rows[, c("x1", "g")]), "x2")
expect_error(predict(fit, setNames(rows[, c("x1", "x2")], c("x1", "z"))), "x2")
expect_error(
  predict(fit, cbind(as.matrix(rows[, c("x1", "x2")]), k = 1)),
  "number of columns in 'test'"
)
