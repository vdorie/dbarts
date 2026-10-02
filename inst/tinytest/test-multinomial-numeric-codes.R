# a numeric multinomial response is a vector of 0-based category codes: code c
# is category c + 1, so the count matrix one-hot expands each row to its own
# column, exactly as the same labels given as a factor do

set.seed(5901)
n <- 12L
x <- matrix(runif(n * 2L), n, 2L)
codes <- rep(c(0, 2, 1), 4L)
control <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 5L,
  n.samples = 5L,
  updateState = FALSE
)

expected <- outer(codes, 0:2, "==") * 1L
sampler <- dbarts(x, codes, family = "multinomial", control = control)
expect_equal(unname(sampler$data@counts), expected)
expect_equal(colnames(sampler$data@counts), c("0", "1", "2"))

twin <- dbarts(x, factor(codes), family = "multinomial", control = control)
expect_identical(
  unname(sampler$data@counts),
  unname(twin$data@counts)
)

rm(sampler, twin, expected, control, codes, x, n)
