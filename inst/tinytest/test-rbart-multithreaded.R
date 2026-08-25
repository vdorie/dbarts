source(
  system.file("common", "friedmanData.R", package = "dbarts"),
  local = TRUE
)

n.g <- 5L
if (getRversion() >= "3.6.0") {
  oldSampleKind <- RNGkind()[3L]
  suppressWarnings(RNGkind(sample.kind = "Rounding"))
}
g <- sample(n.g, length(testData$y), replace = TRUE)
if (getRversion() >= "3.6.0") {
  suppressWarnings(RNGkind(sample.kind = oldSampleKind))
  rm(oldSampleKind)
}

sigma.b <- 1.5
b <- rnorm(n.g, 0, sigma.b)

testData$y <- testData$y + b[g]
testData$g <- g
testData$b <- b
rm(b, sigma.b, g, n.g)

# test that works with multiple threads
x <- testData$x
y <- testData$y
g <- factor(testData$g)

set.seed(0)
fit <- dbarts::rbart_vi(
  y ~ x,
  group.by = g,
  n.samples = 7L,
  n.burn = 0L,
  n.thin = 1L,
  n.chains = 2L,
  n.trees = 25L,
  n.threads = 2L,
  verbose = FALSE
)
expect_inherits(fit, "rbart")

# both chains ran and both are reported: the shapes carry the requested
# chain and sample counts, the draws are finite, and the two chains hold
# different streams rather than one duplicated result
expect_equal(dim(fit$yhat.train), c(2L * 7L, length(y)))
expect_equal(dim(fit$ranef), c(2L * 7L, nlevels(g)))
expect_equal(length(fit$tau), 2L * 7L)
expect_true(all(is.finite(fit$yhat.train)))
expect_false(isTRUE(all.equal(
  fit$yhat.train[seq_len(7L), ],
  fit$yhat.train[7L + seq_len(7L), ]
)))
# the group labels ride through to the reported intercepts
expect_equal(names(fit$ranef.mean), levels(g))
expect_equal(colnames(fit$ranef), levels(g))

rm(fit)

rm(testData)
