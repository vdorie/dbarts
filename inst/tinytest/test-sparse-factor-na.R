# a sparseFactor holds a missing value as an explicit stored entry, and every
# method answers it as a factor does

ff <- factor(
  c("a", NA, "b", "a", NA, "c", "a"),
  levels = c("a", "b", "c", "d")
)
sf <- sparseFactor(ff)
n <- length(ff)

# outcome and warnings of an expression
run <- function(expr) {
  warnings <- character()
  value <- withCallingHandlers(
    tryCatch(expr, error = function(e) {
      structure(conditionMessage(e), class = "err")
    }),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = warnings)
}
same <- function(sparse, dense) {
  expect_equal(as.character(sparse), as.character(dense))
  expect_equal(levels(sparse), levels(dense))
}

# ---- storage ----
expect_equal(sf@i, c(1L, 2L, 4L, 5L))
expect_equal(sf@values, c(NA, 2L, NA, 3L))
expect_equal(sf@reference, "a")
expect_equal(sf@i[is.na(sf@values)], c(1L, 4L))
sfi <- sparseFactor(
  c("b", NA, "c"),
  levels = c("a", "b", "c"),
  i = c(1L, 3L, 6L),
  length = 7L
)
expect_equal(as.character(sfi), c("b", "a", NA, "a", "a", "c", "a"))
sfn <- sparseFactor(c(1L, NA, 2L, NA), levels = c("a", "b"))
expect_equal(as.character(sfn), c("a", NA, "b", NA))
expect_equal(as.character(sparseFactor(c("a", NA, "b"))), c("a", NA, "b"))
expect_equal(levels(sparseFactor(c("b", NA, "a"))), c("a", "b"))
expect_error(sparseFactor(c("a", "z"), levels = "a"), "absent from 'levels'")
expect_error(
  sparseFactor(c("a", NA), levels = c("a", NA)),
  "'levels' cannot contain NA"
)
expect_error(sparseFactor(addNA(ff)), "NA level")
expect_error(
  new(
    "sparseFactor",
    i = 0L,
    values = 9L,
    levels = "a",
    reference = "a",
    length = 2L
  ),
  "level codes"
)
expect_error(
  new(
    "sparseFactor",
    i = 0L,
    values = 1L,
    levels = character(),
    reference = NA_character_,
    length = 2L
  ),
  "can be empty only"
)

# every row missing: no level is needed, as for a factor
allNA <- sparseFactor(c(NA_character_, NA))
expect_equal(levels(allNA), levels(factor(c(NA, NA))))
expect_equal(length(levels(allNA)), 0L)
expect_true(is.na(allNA@reference))
expect_equal(is.na(allNA), c(TRUE, TRUE))
expect_equal(as.character(sparseFactor(c(NA, NA))), c(NA_character_, NA))
expect_equal(
  as.character(sparseFactor(factor(c(NA, NA)))),
  c(NA_character_, NA)
)
expect_equal(length(levels(sparseFactor(character(0)))), 0L)
expect_equal(
  as.character(sparseFactor(c(NA, NA), levels = c("a", "b"))),
  c(NA_character_, NA)
)
same(droplevels(sf[c(2L, 5L)]), droplevels(ff[c(2L, 5L)]))
same(sf[c(2L, 5L), drop = TRUE], ff[c(2L, 5L), drop = TRUE])
expect_true(is.na(droplevels(sf[c(2L, 5L)])@reference))
allA <- sf
levels(allA) <- c(NA, NA, NA, NA)
fA <- ff
levels(fA) <- c(NA, NA, NA, NA)
same(allA, fA)
expect_true(is.na(allA@reference))
# a stored-NA vector takes new levels
fN <- factor(c(NA, NA))
sN <- allNA
levels(fN) <- c("x", "y")
levels(sN) <- c("x", "y")
same(sN, fN)
expect_equal(sN@reference, "x")
same(c(allNA, sf), c(factor(c(NA, NA)), ff))
expect_equal(as.character(c(allNA, allNA)), rep(NA_character_, 4L))
r <- run({
  f <- factor(c(NA, NA))
  f[1L] <- "a"
  f
})
r2 <- run({
  s <- allNA
  s[1L] <- "a"
  s
})
expect_equal(as.character(r2$value), as.character(r$value))
expect_equal(r2$warnings, r$warnings)
expect_equal(length(rep(allNA, 3L)), 6L)
expect_equal(capture.output(print(allNA))[2L], "  levels: ")
expect_equal(
  as.character(summary(allNA)),
  as.character(summary(factor(c(NA, NA))))
)
expect_equal(capture.output(str(allNA)), capture.output(str(factor(c(NA, NA)))))

# ---- reading ----
expect_equal(is.na(sf), is.na(ff))
expect_true(anyNA(sf))
expect_false(anyNA(sparseFactor(factor(c("a", "b")))))
expect_equal(as.character(sf), as.character(ff))
expect_equal(as.integer(sf), as.integer(ff))
expect_equal(as.vector(sf), as.vector(ff))
expect_equal(xtfrm(sf), xtfrm(ff))
expect_equal(order(sf), order(ff))
expect_equal(as.character(sort(sf)), as.character(sort(ff)))
expect_equal(order(sf, na.last = NA), order(ff, na.last = NA))
expect_equal(format(sf), format(ff))
expect_equal(summary(sf), summary(ff))
expect_equal(capture.output(str(sf)), capture.output(str(ff)))
expect_equal(as.character(unique(sf)), as.character(unique(ff)))
expect_equal(duplicated(sf), duplicated(ff))
expect_equal(as.character(factor(sf)), as.character(factor(ff)))
expect_equal(levels(factor(sf)), levels(factor(ff)))
expect_equal(
  as.vector(table(sf, useNA = "ifany")),
  as.vector(table(factor(ff), useNA = "ifany"))
)
if (getRversion() >= "4.6.0") {
  expect_equal(match(sf, "a"), match(ff, "a"))
  expect_equal(sf %in% "a", ff %in% "a")
}
shown <- capture.output(show(sf))
expect_true(any(grepl("(2 missing)", shown, fixed = TRUE)))
expect_false(any(grepl(
  "missing",
  capture.output(show(sparseFactor(ff[c(1L, 3L)])))
)))

# ---- indexing ----
same(sf[c(1L, NA, 9L)], ff[c(1L, NA, 9L)])
expect_equal(as.character(sf[c(1L, NA, 9L)]), c("a", NA, NA))
same(sf[c(TRUE, NA)], ff[c(TRUE, NA)])
same(sf[-10L], ff[-10L])
same(sf[c(2L, 5L, 2L)], ff[c(2L, 5L, 2L)])
same(sf[NA], ff[NA])
expect_error(sf["a"], "position")
expect_error(sf[[NA]], "subscript out of bounds")
expect_error(sf[[9L]], "subscript out of bounds")
expect_error(ff[[9L]], "subscript out of bounds")
expect_equal(as.character(sf[[2L]]), as.character(ff[[2L]]))
same(sf[, drop = TRUE], ff[, drop = TRUE])
same(droplevels(sf), droplevels(ff))
same(sf[c(1L, 2L), drop = TRUE], ff[c(1L, 2L), drop = TRUE])
same(sf[c(1L, 2L, 3L), drop = TRUE], ff[c(1L, 2L, 3L), drop = TRUE])

# ---- assignment ----
assign_same <- function(index, value, expected = NULL) {
  f <- ff
  s <- sf
  rf <- run(f[index] <- value)
  rs <- run(s[index] <- value)
  if (inherits(rf$value, "err")) {
    expect_true(inherits(rs$value, "err"))
    expect_equal(unclass(rs$value), unclass(rf$value))
  } else {
    same(s, f)
  }
  expect_equal(rs$warnings, rf$warnings)
}
assign_same(2L, "b")
assign_same(1L, NA)
assign_same(10L, "b")
assign_same(c(rep(FALSE, 8L), TRUE), "b")
assign_same(c(1L, NA), "b")
assign_same(c(TRUE, NA), "b")
assign_same(NA, "a")
assign_same(c(1L, NA), c("a", "b"))
assign_same(1L, "z")
assign_same(c(1L, 3L), c("z", "b"))
assign_same(c(1L, 2L), factor(c("b", NA)))
assign_same(-1L, "c")
assign_same(0L, "c")
assign_same(3L, character(0))
assign_same(integer(0), character(0))
assign_same(c(3L, 3L), c("a", "c"))
assign_same(11L, NA)
f <- ff
s <- sf
is.na(f) <- 3L
is.na(s) <- 3L
same(s, f)
f <- ff
s <- sf
f[[3L]] <- NA
s[[3L]] <- NA
same(s, f)
f <- ff
s <- sf
length(f) <- 9L
length(s) <- 9L
same(s, f)
expect_equal(length(s@i), 6L)
length(f) <- 2L
length(s) <- 2L
same(s, f)
expect_inherits(s, "sparseFactor")
expect_error(length(s) <- -1L, "invalid value")
f <- ff
s <- sf
levels(f) <- c("a", NA, "c", "d")
levels(s) <- c("a", NA, "c", "d")
same(s, f)
expect_equal(as.character(s)[3L], NA_character_)
f <- ff
s <- sf
levels(f) <- c(NA, "b", "c", "d")
levels(s) <- c(NA, "b", "c", "d")
same(s, f)
expect_false(is.na(s@reference))
expect_true(s@reference %in% s@levels)
f <- ff
s <- sf
levels(f) <- c("x", "x", "y", "z")
levels(s) <- c("x", "x", "y", "z")
same(s, f)
same(rep(sf, 2L), rep(ff, 2L))
same(c(sf, factor(c(NA, "z"))), c(ff, factor(c(NA, "z"))))
same(c(sf, sf), c(ff, ff))
expect_inherits(c(sf, sf), "sparseFactor")
same(
  c(sf, sparseFactor(factor(c(NA, "a")), reference = "a")),
  c(ff, factor(c(NA, "a")))
)

# ---- ops ----
expect_equal(sf == "a", ff == "a")
expect_equal(sf != "a", ff != "a")
missing <- NA_character_
expect_equal(sf == missing, ff == missing)
expect_equal(sf == sf, ff == ff)
expect_equal(
  run(sf < "a")$warnings,
  run(ff < "a")$warnings
)

# ---- frames ----
d <- data.frame(y = seq_len(n), f = sf)
dd <- data.frame(y = seq_len(n), f = ff)
r <- run(capture.output(print(d)))
expect_equal(r$warnings, character())
expect_equal(r$value, capture.output(print(dd)))
expect_equal(
  capture.output(print(head(d, 3L))),
  capture.output(print(head(dd, 3L)))
)
expect_equal(as.character(d[is.na(d$f), ]$f), as.character(dd[is.na(dd$f), ]$f))
expect_equal(d[is.na(d$f), "y"], dd[is.na(dd$f), "y"])
expect_equal(as.character(d[c(1L, NA), ]$f), as.character(dd[c(1L, NA), ]$f))
expect_equal(is.na(d), is.na(dd))
expect_true(anyNA(d))
expect_false(anyNA(d[!is.na(d$f), ]))
# base limits no method reaches, pinned so a change surfaces
expect_equal(nrow(na.omit(d)), n)
expect_equal(nrow(na.omit(dd)), sum(complete.cases(dd)))
expect_error(complete.cases(d), "invalid 'type'")
expect_error(na.fail(d), "invalid 'type'")
if (getRversion() >= "4.6.0") {
  expect_equal(
    as.character(rbind(d, d)$f),
    as.character(rbind(dd, dd)$f)
  )
}

# ---- fits ----
set.seed(11)
nn <- 120L
fd <- factor(sample(c("a", "b", "c"), nn, TRUE, prob = c(0.6, 0.3, 0.1)))
fd[sample(nn, 15L)] <- NA
z <- runif(nn)
y <- rnorm(nn) + (fd %in% "b") + z
xd <- data.frame(f = fd, z = z)
xs <- data.frame(f = sparseFactor(fd, reference = "a"), z = z)
dfd <- cbind(y = y, xd)
dfs <- cbind(y = y, xs)
dfs$f <- xs$f
fitArgs <- list(
  n.samples = 20L,
  n.burn = 10L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE,
  sigest = 1,
  seed = 3L,
  keepTrees = TRUE
)
fitDense <- do.call(bart, c(list(xd, y), fitArgs))
fitSparse <- do.call(bart, c(list(xs, y), fitArgs))
expect_equal(fitSparse$yhat.train, fitDense$yhat.train)
expect_equal(fitSparse$sigma, fitDense$sigma)
fitFD <- do.call(bart, c(list(y ~ z + f, dfd), fitArgs))
fitFS <- do.call(bart, c(list(y ~ z + f, dfs), fitArgs))
expect_equal(fitFS$yhat.train, fitFD$yhat.train)
expect_equal(fitFS$sigma, fitFD$sigma)

# predict on a test frame with NA, and one with another level order
testF <- factor(c("b", NA, "a", "c", NA, "a"), levels = c("a", "b", "c"))
tz <- runif(6L)
testDense <- data.frame(f = testF, z = tz)
testSparse <- data.frame(f = sparseFactor(testF, reference = "a"), z = tz)
expect_equal(predict(fitSparse, testSparse), predict(fitDense, testDense))
expect_equal(predict(fitFS, testSparse), predict(fitFD, testDense))
testF2 <- factor(c("b", NA, "a", "c", NA, "a"), levels = c("c", "b", "a"))
testSparse2 <- data.frame(f = sparseFactor(testF2, reference = "c"), z = tz)
expect_equal(
  predict(fitSparse, testSparse2),
  predict(fitDense, data.frame(f = testF2, z = tz))
)
# a reference unseen in training that no row takes is fine
unseenRef <- sparseFactor(
  c("b", "a", "c"),
  levels = c("zz", "a", "b", "c"),
  reference = "zz",
  i = 1:3,
  length = 3L
)
expect_equal(
  predict(fitSparse, data.frame(f = unseenRef, z = tz[1:3])),
  predict(
    fitDense,
    data.frame(
      f = factor(c("b", "a", "c"), levels = c("a", "b", "c")),
      z = tz[1:3]
    )
  )
)
# one that some row takes is refused
takenRef <- sparseFactor(
  "b",
  levels = c("zz", "a", "b"),
  reference = "zz",
  i = 1L,
  length = 2L
)
expect_error(
  predict(fitSparse, data.frame(f = takenRef, z = tz[1:2])),
  "not present in the training"
)

# a test NA in a column with none in training is refused, or answered NA
xd0 <- xd[!is.na(fd), ]
xs0 <- xs[!is.na(fd), ]
xs0$f <- sparseFactor(fd[!is.na(fd)], reference = "a")
y0 <- y[!is.na(fd)]
fit0d <- do.call(bart, c(list(xd0, y0), fitArgs))
fit0s <- do.call(bart, c(list(xs0, y0), fitArgs))
expect_equal(fit0s$yhat.train, fit0d$yhat.train)
r1 <- run(predict(fit0d, testDense))
r2 <- run(predict(fit0s, testSparse))
expect_true(inherits(r1$value, "err"))
expect_true(inherits(r2$value, "err"))
expect_true(grepl("f", unclass(r2$value), fixed = TRUE))
expect_equal(unclass(r2$value), unclass(r1$value))
p1 <- predict(fit0d, testDense, na.action = na.pass)
p2 <- predict(fit0s, testSparse, na.action = na.pass)
expect_equal(p2, p1)
expect_true(all(is.na(p2[, is.na(testF)])))

# an all-missing column is refused
allMissing <- data.frame(f = sparseFactor(rep(NA_character_, nn)), z = z)
expect_error(
  do.call(bart, c(list(allMissing, y), fitArgs)),
  "predictor columns cannot be entirely missing"
)
allMissingLevels <- data.frame(
  f = sparseFactor(rep(NA_character_, nn), levels = c("a", "b")),
  z = z
)
expect_error(
  do.call(bart, c(list(allMissingLevels, y), fitArgs)),
  "predictor columns cannot be entirely missing"
)

# ---- formula na.action ----
naArgs <- fitArgs
naArgs$n.samples <- 10L
naArgs$n.burn <- 5L
keptRows <- function(fit) fit$fit$data@y
actions <- list(na.omit = na.omit, na.exclude = na.exclude)
for (nm in names(actions)) {
  fs <- do.call(
    bart,
    c(list(y ~ z + f, dfs, na.action = actions[[nm]]), naArgs)
  )
  fd_ <- do.call(
    bart,
    c(list(y ~ z + f, dfd, na.action = actions[[nm]]), naArgs)
  )
  expect_equal(fs$yhat.train, fd_$yhat.train, info = nm)
  expect_equal(fs$na.action, fd_$na.action, info = nm)
  expect_equal(dim(fs$yhat.train), c(10L, nn - 15L), info = nm)
}
fx <- do.call(bart, c(list(y ~ z + f, dfs, na.action = na.exclude), naArgs))
expect_equal(length(fitted(fx)), nn)
expect_equal(unname(is.na(fitted(fx))), unname(is.na(fd)))
rs <- run(do.call(bart, c(list(y ~ z + f, dfs, na.action = na.fail), naArgs)))
rd <- run(do.call(bart, c(list(y ~ z + f, dfd, na.action = na.fail), naArgs)))
expect_true(inherits(rs$value, "err"))
expect_equal(unclass(rs$value), unclass(rd$value))
fdef <- do.call(bart, c(list(y ~ z + f, dfs), naArgs))
expect_equal(dim(fdef$yhat.train), c(10L, nn))
# subset with na.omit keeps the same rows as the dense fit
sub <- which(seq_len(nn) %% 3L != 0L)
fsub_s <- do.call(
  bart,
  c(list(y ~ z + f, dfs, subset = sub, na.action = na.omit), naArgs)
)
fsub_d <- do.call(
  bart,
  c(list(y ~ z + f, dfd, subset = sub, na.action = na.omit), naArgs)
)
expect_equal(fsub_s$yhat.train, fsub_d$yhat.train)
expect_equal(fsub_s$na.action, fsub_d$na.action)

# sparseVector and dgCMatrix columns, against their dense twins
if (requireNamespace("Matrix", quietly = TRUE)) {
  zv <- z
  zv[c(4L, 20L)] <- NA
  zv[zv < 0.7 & !is.na(zv)] <- 0
  dnum <- data.frame(y = y, f = fd, zn = zv)
  dnum <- dnum[!is.na(dnum$f), ]
  yy <- dnum$y
  dsv <- dnum
  dsv$zn <- as(dnum$zn, "sparseVector")
  for (act in c("na.omit", "na.fail")) {
    a <- run(do.call(
      bart,
      c(list(y ~ z + fn, dsv, na.action = get(act)), naArgs)
    ))
    b <- run(do.call(
      bart,
      c(list(y ~ z + fn, dnum, na.action = get(act)), naArgs)
    ))
    if (inherits(b$value, "err")) {
      expect_true(inherits(a$value, "err"), info = act)
      expect_equal(unclass(a$value), unclass(b$value), info = act)
    } else {
      expect_equal(a$value$yhat.train, b$value$yhat.train, info = act)
      expect_equal(a$value$na.action, b$value$na.action, info = act)
    }
  }
  mat <- cbind(zv, w = ifelse(is.na(zv), 0, zv^2))
  mat[mat < 0.5 & !is.na(mat)] <- 0
  dmat <- data.frame(y = y, f = fd)
  dmat$m <- mat
  dmat <- dmat[!is.na(dmat$f), ]
  dgc <- dmat
  dgc$m <- as(mat[!is.na(fd), ], "CsparseMatrix")
  for (act in c("na.omit", "na.fail")) {
    a <- run(do.call(
      bart,
      c(list(y ~ f + m, dgc, na.action = get(act)), naArgs)
    ))
    b <- run(do.call(
      bart,
      c(list(y ~ f + m, dmat, na.action = get(act)), naArgs)
    ))
    if (inherits(b$value, "err")) {
      expect_true(inherits(a$value, "err"), info = act)
      expect_equal(unclass(a$value), unclass(b$value), info = act)
    } else {
      expect_equal(a$value$yhat.train, b$value$yhat.train, info = act)
      expect_equal(a$value$na.action, b$value$na.action, info = act)
    }
  }
}

# test shorter than training, weights and weights.test, a sparse NA column
wts <- runif(nn) + 0.5
dW <- dfs
dW$w <- wts
dWd <- dfd
dWd$w <- wts
tst <- data.frame(y = 0, z = tz, f = testSparse$f)
tstD <- data.frame(y = 0, z = tz, f = testF)
tst$w <- wts[1:6]
tstD$w <- wts[1:6]
dataS <- dbartsData(y ~ z + f, dW, test = tst, weights = w)
dataD <- dbartsData(y ~ z + f, dWd, test = tstD, weights = w)
expect_equal(dataS@weights.test, wts[1:6])
expect_equal(dataS@weights.test, dataD@weights.test)
expect_equal(dataS@weights, dataD@weights)
