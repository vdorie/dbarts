# a bare number is not a category label: a multinomial response is a factor, a
# character vector or an n x K count matrix, and a numeric vector is refused
# at dbarts() with the message bart() gives, whatever its values - 1, 2, 3
# would otherwise build an empty category 0

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

errorMessage <- function(expr) {
  tryCatch(
    {
      expr
      NA_character_
    },
    error = function(e) conditionMessage(e)
  )
}

bartMessage <- errorMessage(bart(
  x,
  codes,
  family = "multinomial",
  n.samples = 5L,
  n.burn = 2L,
  n.trees = 5L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
))
expect_identical(
  bartMessage,
  paste0(
    "family = \"multinomial\" requires a factor (or character) response, ",
    "or an n x K count matrix"
  )
)

numeric.responses <- list(
  zero.based = codes,
  one.based = codes + 1,
  integer = as.integer(codes),
  with.NA = replace(codes, 3L, NA),
  all.NA = rep(NA_real_, n),
  zero.one = rep(c(0, 1), n / 2L),
  one.two = rep(c(1, 2), n / 2L),
  fractional = rep(c(0.5, 2, 1), 4L),
  negative = rep(c(-1, 2, 1), 4L)
)
for (label in names(numeric.responses)) {
  expect_identical(
    errorMessage(dbarts(
      x,
      numeric.responses[[label]],
      family = "multinomial",
      control = control
    )),
    bartMessage,
    info = label
  )
}

# a count matrix, a factor and a character response stay accepted, and the
# last two give the counts the same labels do
expected <- outer(codes, 0:2, "==") * 1L
fromFactor <- dbarts(
  x,
  factor(codes),
  family = "multinomial",
  control = control
)
expect_equal(unname(fromFactor$data@counts), expected)
expect_equal(colnames(fromFactor$data@counts), c("0", "1", "2"))
fromCharacter <- dbarts(
  x,
  as.character(codes),
  family = "multinomial",
  control = control
)
expect_identical(fromCharacter$data@counts, fromFactor$data@counts)
fromMatrix <- dbarts(x, expected, family = "multinomial", control = control)
expect_equal(unname(fromMatrix$data@counts), expected)

# the other doors a numeric vector could reach: a data object's counts takes a
# matrix and refuses a vector by name, a live multinomial sampler's
# $setResponse is unavailable on it, and $setCounts takes a matrix and refuses
# a vector, a factor and a character vector as a category count that is not K
expect_error(
  dbartsData(x, rep(1, n), counts = codes),
  "'counts' must be a numeric matrix"
)
expect_error(
  fromFactor$setResponse(codes),
  "\\$setResponse is not available on a multinomial sampler"
)
expect_error(fromFactor$setCounts(codes), "'counts' must have 3 categories")
expect_error(
  fromFactor$setCounts(factor(codes)),
  "'counts' must have 3 categories"
)
expect_error(
  fromFactor$setCounts(as.character(codes)),
  "'counts' must have 3 categories"
)

# a formula never carries a multinomial response at dbarts()
d <- data.frame(y = codes, x1 = x[, 1L], x2 = x[, 2L])
expect_error(
  dbarts(y ~ x1 + x2, d, family = "multinomial", control = control),
  "multinomial fits currently use the matrix interface"
)

# the other families read a numeric response as they always did: 0/1 is binary
# under "auto", and ordinal keeps whole numbers as its ordered levels
expect_identical(
  dbarts(x, rep(c(0, 1), n / 2L), control = control)$model@family,
  "probit"
)
expect_identical(
  dbarts(x, codes + 1, family = "ordinal", control = control)$model@family,
  "ordinal"
)

rm(
  fromMatrix,
  fromCharacter,
  fromFactor,
  d,
  expected,
  numeric.responses,
  bartMessage,
  errorMessage,
  control,
  codes,
  x,
  n
)
