# bart() and dbarts() take the same multinomial response shapes: labels (a
# factor, a character vector or a one-column character matrix) and counts (a
# numeric matrix or a data frame of numeric columns). A wider character matrix
# and a data frame with a column that is not numeric are refused at both with
# one message.

set.seed(7601)
n <- 30L
x <- matrix(runif(n * 2L), n, 2L)
labels <- sample(c("u", "v", "w"), n, TRUE)
counts <- sapply(c("u", "v", "w"), function(k) as.integer(labels == k))
control <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.samples = 5L,
  updateState = FALSE
)

viaBart <- function(y) {
  bart(
    x,
    y,
    family = "multinomial",
    n.samples = 5L,
    n.burn = 2L,
    n.trees = 5L,
    n.chains = 1L,
    n.threads = 1L,
    seed = 3L,
    verbose = FALSE
  )
}
viaDbarts <- function(y) {
  dbarts(x, y, family = "multinomial", control = control, n.trees = 5L)
}
message.of <- function(expr) {
  tryCatch(
    {
      expr
      NA_character_
    },
    error = function(e) conditionMessage(e)
  )
}

# a data frame of counts is its numeric matrix
frame <- as.data.frame(counts)
expect_identical(
  viaBart(frame)$yhat.train,
  viaBart(counts)$yhat.train
)
expect_identical(viaDbarts(frame)$data@counts, viaDbarts(counts)$data@counts)

# a one-column character matrix is its labels
labelMatrix <- matrix(labels, ncol = 1L)
expect_identical(
  viaBart(labelMatrix)$yhat.train,
  viaBart(labels)$yhat.train
)
expect_identical(
  viaDbarts(labelMatrix)$data@counts,
  viaDbarts(labels)$data@counts
)
expect_true(viaDbarts(labelMatrix)$data@counts.as.factor)

# a wider character matrix and a data frame with a column that is not numeric
# are refused, alike at both
wide <- matrix(as.character(counts), n, 3L)
mixed <- data.frame(u = counts[, 1L], v = as.character(counts[, 2L]))
expected <- paste0(
  "a multinomial count-matrix response must be numeric, not character; ",
  "give labels as a vector, or as a one-column matrix of labels"
)
for (response in list(wide, mixed)) {
  expect_identical(message.of(viaBart(response)), expected)
  expect_identical(message.of(viaDbarts(response)), expected)
}

# a one-column data frame of labels, character or factor, is refused too
for (response in list(
  data.frame(l = labels),
  data.frame(l = factor(labels))
)) {
  expect_identical(message.of(viaBart(response)), expected)
  expect_identical(message.of(viaDbarts(response)), expected)
}
