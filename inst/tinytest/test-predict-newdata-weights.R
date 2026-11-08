# predict at new data reads the weights column of a fit that took its weights
# from a column of its data, with the fit's meaning, as a test frame does, and
# from the new data alone (dec-B648, dec-B595, dec-B596, dec-B638, dec-B641)

source(system.file("common", "countWarnings.R", package = "dbarts"))

set.seed(71)
n <- 40L
d <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
d$w <- as.double(sample(1:3, n, replace = TRUE))
d$m <- rep(c(0, 1, 1, 1), length.out = n)
d$y <- d$x1 + rnorm(n) / sqrt(d$w)
d$yb <- as.double(d$y > 0)
d$yc <- rpois(n, exp(0.5 * d$x1))
d$yo <- factor(cut(d$y, 3L), ordered = TRUE)
d$k <- rbinom(n, 5L, 0.5)
d$f <- 5L - d$k

quick <- function(..., n.chains = 2L) {
  call <- match.call()
  call[[1L]] <- quote(bart)
  call$n.samples <- 20L
  call$n.chains <- n.chains
  call$keepTrees <- TRUE
  call$n.trees <- 10L
  call$n.burn <- 10L
  call$n.threads <- 1L
  call$seed <- 3L
  call$verbose <- FALSE
  eval(call, parent.frame())
}
drawn <- function(fit, newdata, ..., seed = 9L) {
  set.seed(seed)
  predict(fit, newdata, type = "ppd", ...)
}
naRows <- function(draws) {
  unname(which(apply(draws, length(dim(draws)), function(x) all(is.na(x)))))
}

# ---- a precision: the column is the stated weights ----------------------------

fit <- quick(y ~ x1 + x2, d, weights = w)
fromColumn <- drawn(fit, d)
expect_identical(fromColumn, drawn(fit, d, weights = d$w))
expect_false(identical(fromColumn, drawn(fit, d, weights = rep(1, n))))
# the noise follows the weights: every weight times ten, a tenth of the
# variance at the same rows
tenfold <- d
tenfold$w <- 10 * d$w
noise <- function(draws, fit, rows) {
  sd(draws[, rows] - as.vector(extract(fit, "ev")[, rows]))
}
set.seed(9)
ev <- predict(fit, d, type = "ev")
scale <- function(draws) {
  tapply(apply(draws - ev, 2L, sd), d$w, mean)
}
unit <- scale(fromColumn)
big <- scale(drawn(fit, tenfold))
expect_true(all(big < unit * 0.5))
expect_true(all(abs(unit[["1"]] / unit[["3"]] - sqrt(3)) < 0.6))

# a weight of 0 gives NA draws, no warning, and every other draw as it was
zeroed <- d[1:12, ]
zeroed$w[c(3L, 8L)] <- 0
reference <- drawn(fit, d[1:12, ])
warned <- countWarnings(withZero <- drawn(fit, zeroed))
expect_identical(warned, 0L)
expect_identical(naRows(withZero), c(3L, 8L))
expect_identical(withZero[, -c(3L, 8L)], reference[, -c(3L, 8L)])
rm(warned)

# ---- the column is read from newdata alone ------------------------------------

# no column in newdata: draws at weight 1, as before
noColumn <- d[, setdiff(names(d), "w")]
expect_identical(drawn(fit, noColumn), drawn(fit, noColumn, weights = 1))
# a same-named variable in the caller's frame is not used
w <- rep(10, n)
expect_identical(drawn(fit, noColumn), drawn(fit, noColumn, weights = 1))
rm(w)
# weights stated at predict keep their meaning, the column ignored
expect_identical(
  drawn(fit, d, weights = rep(2, n)),
  drawn(fit, noColumn, weights = rep(2, n))
)
# a fit whose weights were a vector draws at weight 1 whatever newdata holds
vectorFit <- quick(y ~ x1 + x2, d, weights = d$w)
expect_null(attr(vectorFit$fit$data, "weights.argument"))
expect_identical(drawn(vectorFit, d), drawn(vectorFit, noColumn))
callerWeights <- d$w
callerFit <- quick(y ~ x1 + x2, d, weights = callerWeights)
asColumn <- d
asColumn$callerWeights <- 7
expect_identical(drawn(callerFit, asColumn), drawn(callerFit, noColumn))
rm(callerWeights)
# an expression of columns is evaluated on newdata
exprFit <- quick(y ~ x1 + x2, d, weights = 2 * w)
expect_identical(drawn(exprFit, d), drawn(exprFit, d, weights = 2 * d$w))
# keepCall = FALSE keeps the record: it rides the data, not the call
noCall <- quick(y ~ x1 + x2, d, weights = w, keepCall = FALSE)
expect_identical(drawn(noCall, d), drawn(noCall, d, weights = d$w))
# the other types read no weight
expect_identical(
  predict(fit, d, type = "ev"),
  predict(fit, noColumn, type = "ev")
)

# ---- a missing weight is answered as a stated one is --------------------------

gap <- d[1:6, ]
gap$w[2L] <- NA
expect_error(
  predict(fit, gap, type = "ppd"),
  "the fit's weights column in 'newdata' has missing values; use na.action"
)
expect_identical(naRows(drawn(fit, gap, na.action = na.pass)), 2L)
expect_identical(dim(drawn(fit, gap, na.action = na.omit))[2L], 5L)
expect_identical(dim(drawn(fit, gap, na.action = na.exclude))[2L], 6L)
# and a stated weights wins, the column unread
expect_identical(dim(drawn(fit, gap, weights = rep(1, 6L)))[2L], 6L)

# ---- a frequency fit and the binary fits: one observation, 0 held out --------

freq <- quick(y ~ x1 + x2, d, weights = w, weights.type = "frequency")
atNew <- d[1:12, ]
atNew$w[5L] <- 0
frequency <- drawn(freq, atNew)
expect_identical(naRows(frequency), 5L)
expect_identical(
  frequency[, -5L],
  drawn(freq, d[1:12, ], weights = rep(1, 12L))[, -5L]
)
probit <- quick(yb ~ x1 + x2, d, weights = m)
atNew <- d[1:8, ]
expect_identical(naRows(drawn(probit, atNew)), c(1L, 5L))
expect_true(all(drawn(probit, atNew)[, c(2L:4L, 6L:8L)] %in% c(0, 1)))
logistic <- quick(yb ~ x1 + x2, d, weights = w, family = "logistic")
atNew <- d[1:8, ]
atNew$w[4L] <- 0
expect_identical(naRows(drawn(logistic, atNew)), 4L)
# a weight of 3 is three copies of a row, not three trials
expect_true(all(drawn(logistic, atNew)[, -4L] %in% c(0, 1)))

# ---- a Student-t fit: a precision --------------------------------------------

student <- quick(y ~ x1 + x2, d, weights = w, family = "student")
expect_identical(drawn(student, d), drawn(student, d, weights = d$w))

# ---- the count families --------------------------------------------------------

nbinom <- quick(yc ~ x1 + x2, d, weights = w, family = "nbinom", n.chains = 1L)
atNew <- d[1:8, ]
atNew$w[c(2L, 6L)] <- 0
counts <- drawn(nbinom, atNew)
expect_identical(naRows(counts), c(2L, 6L))
# one observation per row, not w summed
unit <- drawn(nbinom, d[1:8, ], weights = 1)
expect_identical(counts[, -c(2L, 6L)], unit[, -c(2L, 6L)])
expect_identical(drawn(nbinom, d[1:8, ]), unit)

ordinal <- quick(yo ~ x1 + x2, d, weights = m, family = "ordinal")
atNew <- d[1:8, ]
categories <- drawn(ordinal, atNew)
expect_identical(naRows(categories), c(1L, 5L))
expect_identical(
  categories[, -c(1L, 5L)],
  drawn(ordinal, atNew[, setdiff(names(atNew), "m")])[, -c(1L, 5L)]
)

# ---- a hurdle fit: 0 or 1 -------------------------------------------------------

h <- data.frame(
  a = runif(n),
  b = runif(n),
  y = ifelse(runif(n) < 0.5, 0, 1 + d$w)
)
h$w <- 1
h$w[c(3L, 4L)] <- 0
hurdle <- quick(y ~ a + b, h, weights = w, family = "hurdle.lognormal")
atNew <- h[1:8, ]
atNew$w[6L] <- 0
held <- drawn(hurdle, atNew)
expect_identical(naRows(held), c(3L, 4L, 6L))
expect_identical(held, drawn(hurdle, atNew, weights = atNew$w))

# ---- a binomial fit: the trials come from the response columns ----------------

binomial <- quick(cbind(k, f) ~ x1 + x2, d, weights = m, family = "binomial")
atNew <- d[1:12, ]
trialled <- drawn(binomial, atNew)
expect_identical(naRows(trialled), c(1L, 5L, 9L))
expect_true(max(trialled, na.rm = TRUE) > 1)
expect_true(all(trialled <= 5, na.rm = TRUE))
# the same draws as the trials stated with the weight 0 rows held out by hand
stated <- drawn(binomial, atNew, weights = ifelse(atNew$m == 0, 1, 5))
expect_identical(trialled[, -c(1L, 5L, 9L)], stated[, -c(1L, 5L, 9L)])
# without the response columns a row is one trial
noCounts <- atNew[, setdiff(names(atNew), c("k", "f"))]
expect_true(all(drawn(binomial, noCounts) %in% c(0, 1, NA)))
expect_identical(naRows(drawn(binomial, noCounts)), c(1L, 5L, 9L))
