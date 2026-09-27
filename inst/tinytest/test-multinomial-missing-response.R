# family = "multinomial" routes a missing response through na.action exactly
# as every other family does (dec-A70's standing ruling, dec-A70's neighbor
# dec-B108): the default and na.omit drop the row, na.exclude drops it and
# pads fitted() back with NA (across the category margin), na.fail errors,
# and na.pass keeps the row and trips the generic missing-response check.
# Both response forms - a factor and an n x K count matrix - go through the
# same rule; a count-matrix row is missing when ANY of its cells is NA, not
# only when every cell is.

set.seed(408L)
n <- 60L
x <- matrix(runif(n * 2L), n, 2L)
labels <- sample(c("a", "b", "c"), n, replace = TRUE)
missingRow <- 7L

quick <- list(
  n.trees = 5L,
  n.samples = 6L,
  n.burn = 5L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE,
  seed = 11L
)
fitMultinomial <- function(y, ...) {
  do.call(dbarts::bart, c(list(x, y, family = "multinomial"), quick, list(...)))
}

yFactor <- factor(labels)
yFactor[missingRow] <- NA

counts <- matrix(0L, n, 3L, dimnames = list(NULL, levels(yFactor)))
counts[cbind(seq_len(n), as.integer(yFactor))] <- 1L
# every cell NA, not just one: a fully unobserved row
countsFactorNA <- counts
countsFactorNA[missingRow, ] <- NA_integer_
# ONE cell NA is enough to make the row missing, even though the row still
# sums to a legal single trial if the NA were ignored
countsOneCellNA <- counts
countsOneCellNA[missingRow, 1L] <- NA_integer_
countsOneCellNA[missingRow, 2L] <- 1L

kept <- n - 1L

# --- factor response: default, na.omit, na.exclude, na.fail, na.pass -------

fitDefault <- fitMultinomial(yFactor)
expect_equal(length(fitDefault$y), kept)
expect_inherits(fitDefault[["na.action"]], "exclude")
expect_equal(unname(unclass(fitDefault[["na.action"]])), missingRow)

fitOmit <- fitMultinomial(yFactor, na.action = na.omit)
expect_equal(length(fitOmit$y), kept)

fitExclude <- fitMultinomial(yFactor, na.action = na.exclude)
expect_equal(length(fitExclude$y), kept)
fittedExclude <- fitted(fitExclude, type = "ev")
expect_equal(dim(fittedExclude), c(n, 3L))
expect_true(all(is.na(fittedExclude[missingRow, ])))
expect_true(!anyNA(fittedExclude[-missingRow, ]))
residExclude <- residuals(fitExclude)
expect_equal(dim(residExclude), c(n, 3L))
expect_true(all(is.na(residExclude[missingRow, ])))

expect_error(
  fitMultinomial(yFactor, na.action = na.fail),
  "missing values in object"
)
expect_error(
  fitMultinomial(yFactor, na.action = na.pass),
  "response contains missing values"
)

# --- count-matrix response: the same five, plus the any-cell-NA rule -------

fitCountsDefault <- fitMultinomial(countsFactorNA)
expect_equal(nrow(fitCountsDefault$y), kept)
expect_inherits(fitCountsDefault[["na.action"]], "exclude")

fitCountsOmit <- fitMultinomial(countsFactorNA, na.action = na.omit)
expect_equal(nrow(fitCountsOmit$y), kept)

fitCountsExclude <- fitMultinomial(countsFactorNA, na.action = na.exclude)
fittedCountsExclude <- fitted(fitCountsExclude, type = "ev")
expect_equal(dim(fittedCountsExclude), c(n, 3L))
expect_true(all(is.na(fittedCountsExclude[missingRow, ])))
expect_true(!anyNA(fittedCountsExclude[-missingRow, ]))

expect_error(
  fitMultinomial(countsFactorNA, na.action = na.fail),
  "missing values in object"
)
expect_error(
  fitMultinomial(countsFactorNA, na.action = na.pass),
  "response contains missing values"
)

# a single missing cell drops the row exactly as a fully-missing row does
fitOneCellNA <- fitMultinomial(countsOneCellNA)
expect_equal(nrow(fitOneCellNA$y), kept)
expect_inherits(fitOneCellNA[["na.action"]], "exclude")

# --- a response with no missing values is unaffected -----------------------

fitComplete <- fitMultinomial(factor(labels))
expect_null(fitComplete[["na.action"]])
expect_equal(length(fitComplete$y), n)
