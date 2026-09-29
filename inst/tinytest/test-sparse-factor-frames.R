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
expect_error(sf[n + 1L], "cannot hold NA")
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

# ---- aft survival curves, the sampler's getTrees, na.exclude, '.' ----
if (requireNamespace("survival", quietly = TRUE)) {
  dA <- data.frame(time = exp(d$y), status = rep(c(0, 1, 1, 1), length.out = n))
  dA$z <- d$z
  dA$sf <- d$sf
  fitA <- do.call(
    bart,
    c(
      list(survival::Surv(time, status) ~ sf + z, dA, family = "aft"),
      fitArgs
    )
  )
  times <- c(0.5, 1, 2)
  expect_equal(
    survivalProbabilities(fitA, times, newdata = dA),
    survivalProbabilities(fitA, times, newdata = dA[c("sf", "z")])
  )
}

sampler <- dbarts(
  y ~ sf + z,
  d,
  control = dbartsControl(
    keepTrees = TRUE,
    verbose = FALSE,
    n.chains = 1L,
    n.threads = 1L,
    n.samples = 10L,
    n.burn = 5L
  ),
  sigest = 1
)
sampler$run()
expect_equal(
  sampler$getTrees(newdata = dExtra),
  sampler$getTrees(newdata = exact)
)

fitDot <- do.call(bart, c(list(y ~ ., d), fitArgs))
expect_equal(
  predict(fitDot, newdata = dExtra),
  predict(fitDot, newdata = exact)
)
dNA <- dExtra
dNA$z[3L] <- NA
excluded <- predict(fitDot, newdata = dNA, na.action = na.exclude)
expect_true(all(is.na(excluded[, 3L])))
expect_equal(
  unname(excluded[, -3L]),
  unname(predict(fitDot, newdata = exact[-3L, ]))
)

# ---- factor-like behavior wherever a data frame operation reaches it ----
sf <- d$sf
ff <- f
# as.data.frame follows the vector method's row names
expect_equal(
  rownames(as.data.frame(sf, row.names = paste0("r", seq_len(n)))),
  paste0("r", seq_len(n))
)
expect_equal(
  rownames(as.data.frame(sf, row.names = seq_len(n) + 5L)),
  as.character(seq_len(n) + 5L)
)
expect_error(as.data.frame(sf, row.names = c("a", "b")), "row.names")

# drop = TRUE drops the levels no row takes; an unused reference is replaced
expect_equal(levels(sf[f == "b", drop = TRUE]), "b")
expect_equal(
  as.character(sf[f != "a", drop = TRUE]),
  as.character(ff[f != "a", drop = TRUE])
)
expect_equal(
  levels(sf[f != "a", drop = TRUE]),
  levels(ff[f != "a", drop = TRUE])
)
expect_equal(levels(sf[1:3]), levels(ff[1:3]))
expect_error(sf[NA_integer_], "cannot hold NA")
expect_error(sf[n + 1L], "cannot hold NA")

# replacement
sfr <- sf
sfr[2L] <- "c"
ffr <- ff
ffr[2L] <- "c"
expect_equal(as.character(sfr), as.character(ffr))
expect_inherits(sfr, "sparseFactor")
sfr[c(1L, 3L)] <- factor(c("b", "c"))
ffr[c(1L, 3L)] <- factor(c("b", "c"))
expect_equal(as.character(sfr), as.character(ffr))
sfr[f == "a"] <- "a"
expect_equal(length(sfr@i), sum(ffr != "a" & f != "a"))
expect_error(sfr[1L] <- "zzz", "cannot hold")
expect_error(sfr[n + 2L] <- "a", "gap")
sfx <- sf
sfx[n + 1L] <- "b"
expect_equal(as.character(sfx), c(as.character(ff), "b"))
expect_error(length(sfx) <- 3L, "cannot be set")
dr <- d
dr[2L, "sf"] <- "c"
expect_equal(as.character(dr$sf)[2L], "c")

# combining
expect_equal(as.character(c(sf, sf)), c(as.character(ff), as.character(ff)))
expect_inherits(c(sf, sf), "sparseFactor")
comb <- c(sf, factor(c("z", "a")))
expect_equal(levels(comb), c("a", "b", "c", "z"))
expect_equal(as.character(comb), c(as.character(ff), "z", "a"))
sfb <- sparseFactor(f, reference = "b")
expect_equal(as.character(c(sf, sfb)), c(as.character(ff), as.character(ff)))
expect_error(c(sf, 1), "combined")
expect_equal(nrow(rbind(d, d)), 2L * n)
expect_equal(
  as.character(rbind(d, d)$sf),
  c(as.character(ff), as.character(ff))
)
pieces <- split(d, d$y > 0)
expect_equal(
  sort(as.character(do.call(rbind, pieces)$sf)),
  sort(as.character(ff))
)

# conversion, comparison, ordering, uniqueness, counting
expect_equal(as.factor(sf), ff)
expect_equal(as.vector(sf), as.character(ff))
expect_equal(levels(sf), levels(ff))
expect_equal(nlevels(sf), nlevels(ff))
expect_equal(is.na(sf), is.na(ff))
expect_equal(sf == "b", ff == "b")
expect_equal(sf != "b", ff != "b")
expect_equal("b" == sf, "b" == ff)
expect_equal(sf == sf, ff == ff)
expect_equal(sf == ff, ff == ff)
expect_warning(sf > "a", "not meaningful")
expect_equal(order(sf), order(ff))
expect_equal(as.character(sort(sf)), as.character(sort(ff)))
expect_equal(xtfrm(sf), xtfrm(ff))
expect_equal(as.character(unique(sf)), as.character(unique(ff)))
expect_equal(duplicated(sf), duplicated(ff))
expect_equal(duplicated(d), duplicated(as.data.frame(lapply(d, as.vector))))
expect_equal(nrow(unique(rbind(d, d))), nrow(unique(d)))
expect_equal(c(table(sf)), c(table(ff)))
expect_equal(summary(sf), summary(ff))
expect_equal(as.character(rev(sf)), as.character(rev(ff)))
expect_equal(as.character(factor(sf)), as.character(ff))
expect_equal(nrow(na.omit(d)), n)

# str reads as a factor's line
capture <- capture.output(str(d))
expect_true(any(grepl("Factor w/ 3 levels \"a\",\"b\",\"c\"", capture)))
expect_false(any(grepl("@", capture, fixed = TRUE)))
