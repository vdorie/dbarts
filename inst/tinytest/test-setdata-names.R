# setData and setPredictor read data by name when it carries names, and the names
# are the sampler's or the data are refused (dec-B478), as predict and
# setTestPredictor refuse names that are not the training data's: the sampler's
# own names in another order are matched by name, a name that is not the
# sampler's is refused with the sampler left as it was, and data with no column
# names are read by position with no message. Every sampler is seeded before
# creation: its chain seeds are drawn then.

n.nm <- 60L
set.seed(1L)
frame.nm <- data.frame(a = runif(n.nm), b = runif(n.nm))
frame.nm$y <- frame.nm$a + 2 * frame.nm$b + rnorm(n.nm, 0, 0.1)
control.nm <- dbarts::dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 10L,
  n.samples = 5L,
  n.burn = 0L,
  keepTrees = TRUE,
  updateState = TRUE,
  verbose = FALSE
)
samplerOf <- function(data = frame.nm) {
  set.seed(4L)
  sampler <- dbarts::dbarts(y ~ a + b, data, control = control.nm)
  set.seed(5L)
  invisible(sampler$run(10L, 5L))
  sampler$storeState()
  sampler
}
untouched <- function(sampler) list(sampler$data, sampler$state)
# the warnings a call raises
warningsOf <- function(expr) {
  seen <- character()
  withCallingHandlers(
    expr,
    warning = function(w) {
      seen <<- c(seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  seen
}
refusal <- function(
  given,
  sampler,
  caller = "setData",
  what = "the data's columns"
) {
  paste0(
    "^\\$",
    caller,
    ": names that are not the sampler's: ",
    what,
    " are ",
    given,
    " and the sampler's are ",
    sampler,
    "; columns are matched by name, so name them as the sampler does, or give ",
    "a matrix with no column names to read them by position$"
  )
}
refusalSetPredictor <- function(given, sampler) {
  refusal(given, sampler, "setPredictor", "the columns given")
}
renamed.nm <- data.frame(c = frame.nm$a, e = frame.nm$b, y = frame.nm$y)
other.nm <- as.matrix(renamed.nm[c("c", "e")])
sameOrder.nm <- as.matrix(frame.nm[c("a", "b")])
swapped.nm <- as.matrix(frame.nm[c("b", "a")])

# setData: names that are not the sampler's are refused, every one of them or
# one, and the sampler is left as it was
sampler <- samplerOf()
held <- untouched(sampler)
expect_error(
  sampler$setData(dbarts::dbartsData(y ~ c + e, renamed.nm)),
  pattern = refusal("'c' and 'e'", "'a' and 'b'")
)
expect_error(
  sampler$setData(dbarts::dbartsData(
    y ~ a + e,
    data.frame(a = frame.nm$a, e = frame.nm$b, y = frame.nm$y)
  )),
  pattern = refusal("'a' and 'e'", "'a' and 'b'")
)
# a matrix that carries names is a frame in this: its names are checked
expect_error(
  sampler$setData(dbarts::dbartsData(other.nm, frame.nm$y)),
  pattern = refusal("'c' and 'e'", "'a' and 'b'")
)
# test rows carry the same names, and so does the sparse source the bridge
# refuses for another reason
expect_error(
  sampler$setData(dbarts::dbartsData(
    Matrix::Matrix(other.nm, sparse = TRUE),
    frame.nm$y
  )),
  pattern = refusal("'c' and 'e'", "'a' and 'b'")
)
expect_identical(untouched(sampler), held)

# the sampler's names in another order are matched by name, whatever the carrier
for (carrier in c("formula", "matrix", "frame")) {
  sampler <- samplerOf()
  twin <- samplerOf()
  taken <- switch(
    carrier,
    formula = dbarts::dbartsData(y ~ b + a, frame.nm),
    matrix = dbarts::dbartsData(swapped.nm, frame.nm$y),
    frame = dbarts::dbartsData(frame.nm[c("b", "a")], frame.nm$y)
  )
  expect_identical(
    warningsOf(sampler$setData(taken)),
    character(),
    info = carrier
  )
  twin$setData(dbarts::dbartsData(y ~ a + b, frame.nm))
  expect_identical(
    as.matrix(sampler$data@x),
    as.matrix(twin$data@x),
    info = carrier
  )
  expect_identical(colnames(sampler$data@x), c("a", "b"), info = carrier)
  set.seed(9L)
  fits <- sampler$run(0L, 3L)$train
  set.seed(9L)
  expect_identical(fits, twin$run(0L, 3L)$train, info = carrier)
}

# data with no column names are read by position, with no message
sampler <- samplerOf()
twin <- samplerOf()
unnamed <- unname(sameOrder.nm)
expect_identical(
  warningsOf(sampler$setData(dbarts::dbartsData(unnamed, frame.nm$y))),
  character()
)
twin$setData(dbarts::dbartsData(sameOrder.nm, frame.nm$y))
expect_identical(unname(sampler$data@x), unname(twin$data@x))
# a sampler with no names takes data that carry some, by position, and its
# record takes them
bare <- dbarts::dbarts(unnamed, frame.nm$y, control = control.nm)
expect_null(colnames(bare$data@x))
expect_identical(
  warningsOf(bare$setData(dbarts::dbartsData(y ~ c + e, renamed.nm))),
  character()
)
expect_identical(colnames(bare$data@x), c("c", "e"))

# the names of a sampler that repeats one: the same names are taken, other
# ones are refused as names that are not the sampler's, and the same ones in
# another order cannot be told apart by name
repeated <- cbind(x = frame.nm$a, b = frame.nm$b, x = frame.nm$a^2)
sampler <- dbarts::dbarts(repeated, frame.nm$y, control = control.nm)
expect_identical(
  warningsOf(sampler$setData(dbarts::dbartsData(repeated, frame.nm$y))),
  character()
)
expect_error(
  sampler$setData(dbarts::dbartsData(
    cbind(x = frame.nm$a, b = frame.nm$b, y = frame.nm$a^2),
    frame.nm$y
  )),
  pattern = refusal("'x', 'b' and 'y'", "'x', 'b' and 'x'")
)
expect_error(
  sampler$setData(dbarts::dbartsData(
    repeated[, c(1L, 3L, 2L)],
    frame.nm$y
  )),
  pattern = paste0(
    "^\\$setData: the data's columns cannot be matched by name, a name ",
    "repeating: the column 'x' is in the place of the sampler's 'b'; give ",
    "them in the sampler's order$"
  )
)

# a different number of columns is refused as it was
sampler <- samplerOf()
expect_error(
  sampler$setData(dbarts::dbartsData(frame.nm[c("a", "b", "y")], frame.nm$y)),
  pattern = "requires the same predictors"
)

# setPredictor with the whole of the predictors: a frame under other names is
# refused, one under the sampler's names is matched by name in any order
sampler <- samplerOf()
held <- untouched(sampler)
new.nm <- data.frame(a = rev(frame.nm$a), b = rev(frame.nm$b))
expect_error(
  sampler$setPredictor(renamed.nm[c("c", "e")]),
  pattern = refusalSetPredictor("'c' and 'e'", "'a' and 'b'")
)
expect_error(
  sampler$setPredictor(data.frame(a = new.nm$a, e = new.nm$b)),
  pattern = refusalSetPredictor("'a' and 'e'", "'a' and 'b'")
)
expect_error(
  sampler$setPredictor(other.nm),
  pattern = refusalSetPredictor("'c' and 'e'", "'a' and 'b'")
)
expect_error(
  sampler$setPredictor(Matrix::Matrix(other.nm, sparse = TRUE)),
  pattern = refusalSetPredictor("'c' and 'e'", "'a' and 'b'")
)
expect_identical(untouched(sampler), held)
for (carrier in c("frame", "matrix", "sparse")) {
  sampler <- samplerOf()
  twin <- samplerOf()
  given <- switch(
    carrier,
    frame = new.nm[c("b", "a")],
    matrix = as.matrix(new.nm[c("b", "a")]),
    sparse = Matrix::Matrix(as.matrix(new.nm[c("b", "a")]), sparse = TRUE)
  )
  expect_identical(
    warningsOf(invisible(sampler$setPredictor(given))),
    character(),
    info = carrier
  )
  invisible(twin$setPredictor(as.matrix(new.nm)))
  expect_identical(
    unname(as.matrix(sampler$data@x)),
    unname(as.matrix(twin$data@x)),
    info = carrier
  )
  expect_identical(colnames(sampler$data@x), c("a", "b"), info = carrier)
}
# a matrix with no names is read by position, with no message
sampler <- samplerOf()
expect_identical(
  warningsOf(invisible(sampler$setPredictor(unname(as.matrix(new.nm))))),
  character()
)
expect_equal(
  as.numeric(as.matrix(sampler$data@x)),
  as.numeric(as.matrix(new.nm))
)
# a sampler with no names takes a frame by position
bare <- dbarts::dbarts(unnamed, frame.nm$y, control = control.nm)
expect_identical(
  warningsOf(invisible(bare$setPredictor(renamed.nm[c("c", "e")]))),
  character()
)
expect_equal(as.numeric(as.matrix(bare$data@x)), as.numeric(other.nm))
# a column update has no names to check, and a count of columns is refused
sampler <- samplerOf()
invisible(sampler$setPredictor(rev(frame.nm$a), 1L))
expect_equal(as.numeric(sampler$data@x[, 1L]), rev(frame.nm$a))
expect_error(
  sampler$setPredictor(data.frame(a = new.nm$a, b = new.nm$b, a2 = new.nm$a)),
  pattern = "dimension of x must be equal to 2"
)
