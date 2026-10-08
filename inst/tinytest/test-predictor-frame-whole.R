# every route that takes new predictors - setPredictor, setTestPredictor,
# setTestPredictorAndOffset and predict - takes a data frame, coded by label,
# and a numeric matrix only where every predictor column is numeric, refused
# by name on a design with a factor column (dec-B359, dec-B360)

set.seed(5)
n <- 60L
frame <- data.frame(
  a = runif(n),
  f = factor(sample(c("u", "v", "w"), n, replace = TRUE)),
  b = runif(n)
)
frame$y <- frame$a + as.integer(frame$f) / 3 + rnorm(n, 0, 0.2)
control <- dbartsControl(
  n.trees = 10L,
  n.chains = 1L,
  n.threads = 1L,
  n.samples = 5L,
  n.burn = 5L,
  updateState = FALSE
)
make <- function(test = NULL) {
  dbarts(y ~ a + f + b, frame, test = test, control = control)
}
predictors <- frame[c("a", "f", "b")]

# setPredictor, no column: a data frame coded by label
sampler <- make()
invisible(sampler$run())
levelsBefore <- attr(sampler$data@x, "factor.levels")
replacement <- predictors[rev(seq_len(n)), ]
sampler$setPredictor(replacement, forceUpdate = TRUE)
expect_equal(
  unname(sampler$data@x[, "a"]),
  replacement$a
)
expect_equal(
  unname(sampler$data@x[, "f"]),
  as.double(match(as.character(replacement$f), levelsBefore[[2L]]) - 1L)
)
expect_identical(attr(sampler$data@x, "factor.levels"), levelsBefore)
# columns found by name, in any order
shuffled <- replacement[c("b", "f", "a")]
sampler$setPredictor(shuffled, forceUpdate = TRUE)
expect_equal(unname(sampler$data@x[, "b"]), replacement$b)
# a level the rows miss keeps its declared place
sampler$setPredictor(
  within(predictors, f <- factor(rep("u", n), levels = c("u", "v", "w"))),
  forceUpdate = TRUE
)
expect_identical(attr(sampler$data@x, "factor.levels"), levelsBefore)
# a label the column does not declare, a number for a factor, and a label for
# a number are each refused by name
expect_error(
  sampler$setPredictor(
    within(predictors, f <- factor(rep("z", n))),
    forceUpdate = TRUE
  ),
  "'f'"
)
expect_error(
  sampler$setPredictor(
    within(predictors, f <- as.double(f)),
    forceUpdate = TRUE
  ),
  "'f' is categorical"
)
expect_error(
  sampler$setPredictor(
    within(predictors, a <- as.character(a)),
    forceUpdate = TRUE
  ),
  "'a' is numeric in the sampler"
)
# a matrix is refused on a design with a factor column
expect_error(
  sampler$setPredictor(as.matrix(sampler$data@x), forceUpdate = TRUE),
  "'x' is a numeric matrix, but the predictor 'f' is a factor",
  fixed = TRUE
)
# the sampler still runs and copies
expect_true(all(is.finite(sampler$run(0L, 2L)$train)))
expect_true(all(is.finite(sampler$copy()$run(0L, 2L)$train)))

# an all-numeric design takes a matrix, and a data frame alike
numeric.sampler <- dbarts(
  frame$y ~ a + b,
  frame,
  control = control
)
invisible(numeric.sampler$run())
numeric.matrix <- as.matrix(frame[c("a", "b")])[rev(seq_len(n)), ]
numeric.sampler$setPredictor(numeric.matrix, forceUpdate = TRUE)
expect_equal(
  unname(numeric.sampler$data@x[, "a"]),
  unname(numeric.matrix[, "a"])
)
numeric.sampler$setPredictor(frame[c("b", "a")], forceUpdate = TRUE)
expect_equal(unname(numeric.sampler$data@x[, "a"]), frame$a)

# setTestPredictor and setTestPredictorAndOffset: the same rule
test.rows <- predictors[1:10, ]
test.sampler <- make(test = test.rows)
test.sampler$setTestPredictor(test.rows[10:1, ])
expect_equal(unname(test.sampler$data@x.test[, "a"]), rev(test.rows$a))
expect_error(
  test.sampler$setTestPredictor(as.matrix(test.sampler$data@x)[1:10, ]),
  "'x.test' is a numeric matrix, but the predictor 'f' is a factor",
  fixed = TRUE
)
expect_error(
  test.sampler$setTestPredictorAndOffset(
    as.matrix(test.sampler$data@x)[1:10, ],
    NULL
  ),
  "the predictor 'f' is a factor"
)
test.sampler$setTestPredictorAndOffset(test.rows, NULL)
expect_equal(nrow(test.sampler$data@x.test), 10L)
expect_true(all(is.finite(test.sampler$run(0L, 2L)$test)))

# predict: a matrix of codes read 2 as the second level; now refused by name
fit <- suppressWarnings(bart(
  y ~ a + f + b,
  frame,
  n.trees = 10L,
  n.burn = 5L,
  n.samples = 5L,
  n.chains = 1L,
  n.threads = 1L,
  keepTrees = TRUE,
  verbose = FALSE
))
expect_error(
  predict(fit, cbind(a = 0.5, f = 2, b = 0.5)),
  "'newdata' is a numeric matrix, but the predictor 'f' is a factor",
  fixed = TRUE
)
expect_equal(
  dim(predict(fit, data.frame(a = 0.5, f = "v", b = 0.5))),
  c(5L, 1L)
)
# the numeric design's fits keep taking matrices
fit.numeric <- suppressWarnings(bart(
  frame[c("a", "b")],
  frame$y,
  n.trees = 10L,
  n.burn = 5L,
  n.samples = 5L,
  n.chains = 1L,
  n.threads = 1L,
  keepTrees = TRUE,
  verbose = FALSE
))
expect_equal(
  dim(predict(fit.numeric, cbind(a = 0.5, b = 0.5))),
  c(5L, 1L)
)

# the sampler's own predict method and bart(test = ) on a design with a factor
# column refuse a matrix by name as well
sampler.fit <- make()
invisible(sampler.fit$run())
expect_error(
  sampler.fit$predict(as.matrix(sampler.fit$data@x)[1:5, ]),
  "'x.test' is a numeric matrix, but the predictor 'f' is a factor",
  fixed = TRUE
)
expect_equal(length(sampler.fit$predict(predictors[1:5, ])), 5L)
expect_error(
  dbarts(
    y ~ a + f + b,
    frame,
    test = as.matrix(sampler.fit$data@x)[1:5, ],
    control = control
  ),
  "the predictor 'f' is a factor"
)
expect_error(
  suppressWarnings(bart(
    predictors,
    frame$y,
    test = as.matrix(sampler.fit$data@x)[1:5, ],
    n.trees = 5L,
    n.burn = 2L,
    n.samples = 2L,
    n.chains = 1L,
    n.threads = 1L,
    verbose = FALSE
  )),
  "the predictor 'f' is a factor"
)
