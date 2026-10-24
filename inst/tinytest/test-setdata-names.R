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

# three columns, handed in each cyclic order: a permutation of two is its own
# inverse, so only three columns tell the match by name from its inverse. The
# columns sit at different places (means 0.5, 10.5 and 20.5), so a column in
# another's place shows
three.nm <- data.frame(
  a = frame.nm$a,
  b = 10 + frame.nm$b,
  c = 20 + frame.nm$a * frame.nm$b
)
three.nm$y <- three.nm$a + three.nm$c / 20 + rnorm(n.nm, 0, 0.1)
threeOf <- function() {
  set.seed(4L)
  sampler <- dbarts::dbarts(y ~ a + b + c, three.nm, control = control.nm)
  set.seed(5L)
  invisible(sampler$run(10L, 5L))
  sampler$storeState()
  sampler
}
columnsOf <- function(sampler) {
  unclass(as.matrix(sampler$data@x))[, c("a", "b", "c")]
}
original <- as.matrix(three.nm[c("a", "b", "c")])
for (order in list(c("c", "a", "b"), c("b", "c", "a"), c("c", "b", "a"))) {
  info <- paste(order, collapse = "")
  # setPredictor: a frame, a named matrix and a sparse matrix
  carriers <- list(
    frame = three.nm[order],
    matrix = as.matrix(three.nm[order]),
    sparse = Matrix::Matrix(as.matrix(three.nm[order]), sparse = TRUE)
  )
  for (carrier in names(carriers)) {
    sampler <- threeOf()
    sampler$setPredictor(
      carriers[[carrier]],
      forceUpdate = TRUE
    )
    expect_equal(
      unname(columnsOf(sampler)),
      unname(original),
      info = paste(carrier, info)
    )
    expect_identical(
      colnames(sampler$data@x),
      c("a", "b", "c"),
      info = paste(carrier, info)
    )
  }
  # setData: a formula, a matrix and a frame
  for (carrier in c("formula", "matrix", "frame")) {
    sampler <- threeOf()
    sampler$setData(switch(
      carrier,
      formula = dbarts::dbartsData(
        stats::as.formula(paste("y ~", paste(order, collapse = " + "))),
        three.nm
      ),
      matrix = dbarts::dbartsData(as.matrix(three.nm[order]), three.nm$y),
      frame = dbarts::dbartsData(three.nm[order], three.nm$y)
    ))
    expect_equal(
      unname(columnsOf(sampler)),
      unname(original),
      info = paste("setData", carrier, info)
    )
  }
}

# a sampler made from a matrix with no column names and a table of levels for
# its second column names no column of the data's, so a number in the factor's
# place is named by the data's column and the place
codes.nm <- unname(cbind(frame.nm$a, as.double(rep_len(0:2, n.nm))))
attr(codes.nm, "varTypes") <- c(0L, 1L)
attr(codes.nm, "factor.levels") <- list(NULL, c("u", "v", "w"))
bare.factor <- dbarts::dbarts(codes.nm, frame.nm$y, control = control.nm)
expect_null(colnames(bare.factor$data@x))
expect_error(
  bare.factor$setData(dbarts::dbartsData(
    data.frame(a = frame.nm$a, b = frame.nm$b),
    frame.nm$y
  )),
  pattern = paste0(
    "^\\$setData: the data's column 'b', in the place of factor 2, holds ",
    "numbers; give it as a factor of the sampler's labels$"
  )
)
