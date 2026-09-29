# a data frame holding a sparseFactor column behaves like a dense-factor
# frame: predict and test take the same frames, and rows subset and print

set.seed(7)
n <- 150L
f <- factor(sample(c("a", "b", "c"), n, TRUE, prob = c(0.7, 0.2, 0.1)))
d <- data.frame(y = rnorm(n) + (f == "b"), z = runif(n) + 1)
d$sf <- sparseFactor(f, reference = "a")
exact <- data.frame(sf = d$sf, z = d$z)
dExtra <- cbind(d, w = runif(n))

fitArgs <- list(
  n.samples = 30L,
  n.burn = 10L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE,
  sigest = 1,
  seed = 1L,
  keepTrees = TRUE
)
fit <- do.call(bart, c(list(y ~ sf + z, d), fitArgs))
base <- predict(fit, newdata = exact)
expect_equal(dim(base), c(30L, n))
expect_equal(predict(fit, newdata = d), base)
expect_equal(predict(fit, newdata = dExtra), base)
expect_equal(predict(fit, newdata = dExtra[c(5L, 2L), ]), base[, c(5L, 2L)])

# a transformed dense term is replayed, the sparse term rides along
fitLog <- do.call(bart, c(list(y ~ sf + log(z), d), fitArgs))
expect_equal(
  predict(fitLog, newdata = d),
  predict(fitLog, newdata = data.frame(sf = d$sf, z = d$z))
)
expect_error(predict(fitLog, newdata = d["sf"]), "missing variable")

# test = at fit time takes the same frames
fitTest <- do.call(bart, c(list(y ~ sf + z, d, test = dExtra), fitArgs))
expect_equal(dim(fitTest$yhat.test), c(30L, n))
expect_equal(fitTest$yhat.test, predict(fitTest, newdata = exact))
fitTestLog <- do.call(bart, c(list(y ~ sf + log(z), d, test = d), fitArgs))
expect_equal(dim(fitTestLog$yhat.test), c(30L, n))

# only sparse terms
fitOnly <- do.call(bart, c(list(y ~ sf, d), fitArgs))
expect_equal(dim(predict(fitOnly, newdata = d)), c(30L, n))

# sparseFactor methods follow the factor's
sf <- d$sf
expect_equal(as.character(sf), as.character(f))
expect_equal(as.character(sf[c(3L, 1L, 3L)]), as.character(f[c(3L, 1L, 3L)]))
expect_equal(as.character(sf[-(1:140)]), as.character(f[-(1:140)]))
expect_equal(as.character(sf[f == "b"]), as.character(f[f == "b"]))
expect_equal(as.character(sf[c(TRUE, FALSE)]), as.character(f[c(TRUE, FALSE)]))
expect_equal(as.character(sf[0L]), character(0L))
expect_equal(as.character(sf[]), as.character(f))
expect_inherits(sf[1:10], "sparseFactor")
expect_equal(sf[1:10]@levels, sf@levels)
expect_equal(sf[1:10]@reference, sf@reference)
expect_equal(length(sf[1:10]), 10L)
expect_true(length(sf[1:10]@i) <= 10L)
expect_equal(format(sf), format(as.character(f)))
expect_error(sf["a"], "position")
expect_error(sf[n + 1L], "missing")
expect_error(sf[1L, 2L], "dimensions")

# data frame operations
expect_silent(print(head(d)))
expect_equal(as.character(d[d$y > 1, "sf"]), as.character(f[d$y > 1]))
sub <- d[d$y > 0.5, ]
expect_equal(as.character(sub$sf), as.character(f[d$y > 0.5]))
expect_equal(rownames(d[c(3L, 1L, 3L), ]), c("3", "1", "3.1"))
expect_equal(as.character(d[-1L, ]$sf), as.character(f[-1L]))
expect_silent(str(d))
expect_equal(as.character(as.data.frame(sf)$sf), as.character(f))

# a fit on a subset frame equals a fit on the same rows built directly
keep <- which(d$y > 0)
direct <- data.frame(y = d$y[keep], z = d$z[keep])
direct$sf <- sparseFactor(f[keep], reference = "a")
fitSub <- do.call(bart, c(list(y ~ sf + z, d[keep, ]), fitArgs))
fitDirect <- do.call(bart, c(list(y ~ sf + z, direct), fitArgs))
expect_equal(unname(fitSub$yhat.train), unname(fitDirect$yhat.train))
