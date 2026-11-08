# $setResponse takes a factor (or a character vector read as one) where the
# sampler's response was written as one: a multinomial sampler whose rows each
# hold one trial, and a probit or logistic sampler built from a two-level
# factor. Labels are matched to the levels the sampler recorded by name and
# never by position, so a releveled factor does not change the model.

set.seed(5521L)
n <- 90L
x <- matrix(runif(n * 3L), n, 3L)
categories <- c("a", "b", "c")
labelsOne <- sample(categories, n, replace = TRUE)
labelsTwo <- sample(categories, n, replace = TRUE)

ctl <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  updateState = FALSE,
  seed = 7L
)
build <- function(y, ...) {
  dbarts(x, y, control = ctl, n.trees = 10L, verbose = FALSE, ...)
}
asCounts <- function(labels) {
  counts <- matrix(0L, n, 3L, dimnames = list(NULL, categories))
  counts[cbind(seq_len(n), match(labels, categories))] <- 1L
  counts
}

# ---- multinomial ----------------------------------------------------------

fresh <- build(factor(labelsTwo), family = "multinomial")
for (response in list(
  factor(labelsTwo),
  factor(labelsTwo, levels = c("c", "a", "b")),
  factor(labelsTwo, levels = c("c", "b", "a", "unused")),
  labelsTwo
)) {
  sampler <- build(factor(labelsOne), family = "multinomial")
  sampler$setResponse(response)
  expect_identical(sampler$data@counts, fresh$data@counts)
  expect_identical(sampler$data@counts, asCounts(labelsTwo))
  expect_identical(sampler$data@y, fresh$data@y)
  expect_identical(
    sampler$run(2L, 3L),
    build(
      factor(labelsTwo),
      family = "multinomial"
    )$run(2L, 3L)
  )
}

# the same response as the matching count matrix
viaMatrix <- build(factor(labelsOne), family = "multinomial")
viaMatrix$setResponse(asCounts(labelsTwo))
viaFactor <- build(factor(labelsOne), family = "multinomial")
viaFactor$setResponse(factor(labelsTwo, levels = c("b", "c", "a")))
expect_identical(viaFactor$run(2L, 3L), viaMatrix$run(2L, 3L))

# a label the sampler lacks is refused by name, the sampler untouched
sampler <- build(factor(labelsOne), family = "multinomial")
twin <- build(factor(labelsOne), family = "multinomial")
foreign <- replace(labelsTwo, 3L, "zzz")
expect_error(sampler$setResponse(factor(foreign)), "'zzz'", fixed = TRUE)
expect_identical(sampler$data@counts, twin$data@counts)
expect_identical(sampler$run(2L, 3L), twin$run(2L, 3L))
expect_error(sampler$setResponse(foreign), "not among them", fixed = TRUE)
expect_error(
  sampler$setResponse(factor(replace(labelsTwo, 1L, NA))),
  "missing values"
)
expect_error(
  sampler$setResponse(factor(labelsTwo[-1L])),
  "one label per row"
)

# a numeric vector of codes stays refused
expect_error(
  sampler$setResponse(match(labelsTwo, categories)),
  "a vector cannot say how each row's trials split among the categories",
  fixed = TRUE
)

# a row of more than one trial refuses a factor, naming the count matrix; a
# sampler that records no category names has no labels to match
grouped <- asCounts(labelsOne)
grouped[1L, ] <- c(2L, 0L, 1L)
fromGrouped <- build(grouped, family = "multinomial")
expect_error(
  fromGrouped$setResponse(factor(labelsTwo)),
  "give the n x K count matrix",
  fixed = TRUE
)
unnamed <- build(unname(asCounts(labelsOne)), family = "multinomial")
expect_error(
  unnamed$setResponse(factor(labelsTwo)),
  "records no category names",
  fixed = TRUE
)

# the record of the written form follows the last response: draws come back as
# categories after a factor and as count vectors after a matrix
fromMatrix <- build(asCounts(labelsOne), family = "multinomial")
expect_false(isTRUE(fromMatrix$data@counts.as.factor))
expect_identical(
  dim(samplePriorPredictive(fromMatrix, n.samples = 2L, type = "ppd")),
  c(2L, n, 3L)
)
fromMatrix$setResponse(factor(labelsTwo))
expect_true(fromMatrix$data@counts.as.factor)
expect_identical(
  dim(samplePriorPredictive(fromMatrix, n.samples = 2L, type = "ppd")),
  c(2L, n)
)
fromMatrix$setResponse(asCounts(labelsOne))
expect_false(isTRUE(fromMatrix$data@counts.as.factor))

# ---- probit and logistic --------------------------------------------------

binaryLabels <- sample(c("no", "yes"), n, replace = TRUE)
binaryOther <- sample(c("no", "yes"), n, replace = TRUE)
binaryOther[1:2] <- c("no", "yes")
code <- function(labels) as.double(labels == "yes")
for (family in c("probit", "logistic")) {
  sampler <- build(factor(binaryLabels), family = family)
  expect_identical(sampler$data@response.levels, c("no", "yes"))
  fresh <- build(factor(binaryOther), family = family)

  # whatever the new factor's level order, "no" is 0 and "yes" is 1
  for (response in list(
    factor(binaryOther),
    factor(binaryOther, levels = c("yes", "no")),
    factor(binaryOther, levels = c("yes", "no", "unused")),
    binaryOther
  )) {
    sampler$setResponse(response)
    expect_identical(sampler$data@y, code(binaryOther), info = family)
  }
  sampler$setResponse(factor(binaryOther, levels = c("yes", "no")))
  expect_identical(sampler$run(2L, 3L), fresh$run(2L, 3L), info = family)

  # a releveled factor does not flip the model
  flipped <- factor(binaryOther, levels = c("yes", "no"))
  expect_identical(as.double(flipped), 2 - code(binaryOther))
  sampler$setResponse(flipped)
  expect_identical(sampler$data@y, code(binaryOther), info = family)

  # a label the sampler lacks, and a third level, are refused by name
  expect_error(
    sampler$setResponse(factor(replace(binaryOther, 3L, "maybe"))),
    "'maybe'",
    fixed = TRUE,
    info = family
  )
  expect_error(
    sampler$setResponse(replace(binaryOther, 3L, "maybe")),
    "not among them",
    fixed = TRUE,
    info = family
  )
  expect_identical(sampler$data@y, code(binaryOther), info = family)

  # 0 and 1 as before, and a logical
  sampler$setResponse(code(binaryLabels))
  expect_identical(sampler$data@y, code(binaryLabels), info = family)
  sampler$setResponse(binaryLabels == "yes")
  expect_identical(sampler$data@y, code(binaryLabels), info = family)
}

# a sampler whose response was written as numbers records no levels and keeps
# refusing a factor, naming 0 and 1; a factor of one level in use is no 0/1
numeric <- build(code(binaryLabels), family = "probit")
expect_null(numeric$data@response.levels)
for (response in list(
  factor(binaryOther),
  binaryOther,
  factor(rep("1", n))
)) {
  expect_error(numeric$setResponse(response), "give 0 and 1", fixed = TRUE)
}
numeric$setResponse(code(binaryOther))
expect_identical(numeric$data@y, code(binaryOther))
