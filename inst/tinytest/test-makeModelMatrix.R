getTestDataFrame <- function() {
  set.seed(42)
  df <- data.frame(
    iv = seq.int(10L),
    rv = runif(10),
    f = factor(sample(3L, 10L, TRUE), labels = c("a", "b", "c"))
  )
  df[[4L]] <- matrix(
    rbinom(20L, 3L, 0.5),
    10L,
    dimnames = list(NULL, c("a", "b"))
  )
  df[[5L]] <- matrix(rnorm(20L), 10L, dimnames = list(NULL, c("a", "b")))
  names(df) <- c(names(df)[1L:3L], "im", "rm")
  df
}

# test that make model matrix works on default
df <- getTestDataFrame()
mm <- dbarts::makeModelMatrixFromDataFrame(df)

expect_equal(ncol(mm), 9L)
expect_equal(
  colnames(mm),
  c("iv", "rv", "f.a", "f.b", "f.c", "im.a", "im.b", "rm.a", "rm.b")
)
expect_equal(mm[, "iv"], df$iv)
expect_equal(mm[, "rv"], df$rv)
expect_equal(mm[, "f.a"], ifelse(df$f == "a", 1L, 0L))
expect_equal(mm[, "f.b"], ifelse(df$f == "b", 1L, 0L))
expect_equal(mm[, "f.c"], ifelse(df$f == "c", 1L, 0L))
expect_equal(mm[, "im.a"], df$im[, "a"])
expect_equal(mm[, "im.b"], df$im[, "b"])
expect_equal(mm[, "rm.a"], df$rm[, "a"])
expect_equal(mm[, "rm.b"], df$rm[, "b"])

rm(mm, df)

# test that make model matrix handles empty names
df <- getTestDataFrame()
names(df) <- NULL
mm <- dbarts::makeModelMatrixFromDataFrame(df)
expect_equal(
  colnames(mm),
  c("", "", "a", "b", "c", "a", "b", "a", "b")
)

df <- getTestDataFrame()
df$f <- as.factor(as.integer(df$f))
mm <- dbarts::makeModelMatrixFromDataFrame(df)
expect_equal(
  colnames(mm),
  c("iv", "rv", "f.1", "f.2", "f.3", "im.a", "im.b", "rm.a", "rm.b")
)

df <- getTestDataFrame()
colnames(df$im) <- NULL
mm <- dbarts::makeModelMatrixFromDataFrame(df)
expect_equal(
  colnames(mm),
  c("iv", "rv", "f.a", "f.b", "f.c", "im.1", "im.2", "rm.a", "rm.b")
)

df <- getTestDataFrame()
colnames(df$rm) <- NULL
mm <- dbarts::makeModelMatrixFromDataFrame(df)
expect_equal(
  colnames(mm),
  c("iv", "rv", "f.a", "f.b", "f.c", "im.a", "im.b", "rm.1", "rm.2")
)

rm(mm, df)


# test that make model matrix drops useless columns
df <- getTestDataFrame()
df$iv <- rep(1L, 10L)
mm <- dbarts::makeModelMatrixFromDataFrame(df)
expect_equal(ncol(mm), 8L)
expect_equal(attr(mm, "drop")$iv, TRUE)
expect_equal(
  colnames(mm),
  c("rv", "f.a", "f.b", "f.c", "im.a", "im.b", "rm.a", "rm.b")
)

df <- getTestDataFrame()
df$rv <- rep(pi, 10L)
mm <- dbarts::makeModelMatrixFromDataFrame(df)
expect_equal(ncol(mm), 8L)
expect_equal(attr(mm, "drop")$rv, TRUE)
expect_equal(
  colnames(mm),
  c("iv", "f.a", "f.b", "f.c", "im.a", "im.b", "rm.a", "rm.b")
)

df <- getTestDataFrame()
## creates a factor with one unused level
df$f <- factor(
  rep(seq.int(3L), c(5L, 5L, 1L)),
  labels = c("a", "b", "c")
)[1L:10L]
mm <- dbarts::makeModelMatrixFromDataFrame(df)
expect_equal(ncol(mm), 7L)
expect_equal(attr(mm, "drop")$f, c(5L, 5L, 0L))
expect_equal(
  colnames(mm),
  c("iv", "rv", "f.b", "im.a", "im.b", "rm.a", "rm.b")
)

df <- getTestDataFrame()
df$im[, 1L] <- rep(1L, 10L)
mm <- dbarts::makeModelMatrixFromDataFrame(df)
expect_equal(ncol(mm), 8L)
expect_equal(attr(mm, "drop")$im, c(TRUE, FALSE))
expect_equal(
  colnames(mm),
  c("iv", "rv", "f.a", "f.b", "f.c", "im.b", "rm.a", "rm.b")
)
expect_equal(as.double(mm[, 7L:8L]), as.double(df$rm))

df <- getTestDataFrame()
df$im[, 1L] <- rep(1L, 10L)
df$im[, 2L] <- rep(2L, 10L)
mm <- dbarts::makeModelMatrixFromDataFrame(df)
expect_equal(ncol(mm), 7L)
expect_equal(attr(mm, "drop")$im, c(TRUE, TRUE))
expect_equal(
  colnames(mm),
  c("iv", "rv", "f.a", "f.b", "f.c", "rm.a", "rm.b")
)

df <- getTestDataFrame()
df$rm[, 2L] <- rep(pi, 10L)
mm <- dbarts::makeModelMatrixFromDataFrame(df)
expect_equal(ncol(mm), 8L)
expect_equal(attr(mm, "drop")$rm, c(FALSE, TRUE))
expect_equal(
  colnames(mm),
  c("iv", "rv", "f.a", "f.b", "f.c", "im.a", "im.b", "rm.a")
)
expect_equal(as.integer(mm[, 6L:7L]), as.integer(df$im))

rm(mm, df)


# test that make model matrix doesn't drop useless columns when drop = FALSE
df <- getTestDataFrame()
df$iv <- rep(1L, 10L)
mm <- dbarts::makeModelMatrixFromDataFrame(df, FALSE)
expect_equal(ncol(mm), 9L)
expect_equal(mm[, "iv"], df$iv)

df <- getTestDataFrame()
df$rv <- rep(pi, 10L)
mm <- dbarts::makeModelMatrixFromDataFrame(df, FALSE)
expect_equal(mm[, "rv"], df$rv)

df <- getTestDataFrame()
df$f <- factor(
  rep(seq.int(3L), c(5L, 5L, 1L)),
  labels = c("a", "b", "c")
)[1L:10L]
mm <- dbarts::makeModelMatrixFromDataFrame(df, FALSE)
expect_equal(mm[, "f.a"], c(rep(1L, 5), rep(0L, 5)))
expect_equal(mm[, "f.b"], c(rep(0L, 5), rep(1L, 5)))
expect_equal(mm[, "f.c"], rep(0L, 10))

df <- getTestDataFrame()
df$im[, 1L] <- rep(1L, 10L)
mm <- dbarts::makeModelMatrixFromDataFrame(df, FALSE)
expect_equal(as.integer(mm[, 6L:7L]), as.integer(df$im))

df <- getTestDataFrame()
df$rm[, 2L] <- rep(pi, 10L)
mm <- dbarts::makeModelMatrixFromDataFrame(df, FALSE)
expect_equal(as.double(mm[, 8L:9L]), as.double(df$rm))

rm(mm, df)


# test that make model matrix respects drop argument when a list
df <- getTestDataFrame()
drop <- list(
  TRUE,
  FALSE,
  as.integer(table(df$f)),
  c(FALSE, FALSE),
  c(FALSE, FALSE)
)
names(drop) <- names(df)

mm <- dbarts::makeModelMatrixFromDataFrame(df, drop)
expect_equal(ncol(mm), 8L)
expect_equal(
  colnames(mm),
  c("rv", "f.a", "f.b", "f.c", "im.a", "im.b", "rm.a", "rm.b")
)

drop$iv <- FALSE
drop$rv <- TRUE
mm <- dbarts::makeModelMatrixFromDataFrame(df, drop)
expect_equal(ncol(mm), 8L)
expect_equal(
  colnames(mm),
  c("iv", "f.a", "f.b", "f.c", "im.a", "im.b", "rm.a", "rm.b")
)

drop$rv <- FALSE
drop$f <- c(1L, 0L, 1L)
mm <- dbarts::makeModelMatrixFromDataFrame(df, drop)
expect_equal(ncol(mm), 7L)
expect_equal(
  colnames(mm),
  c("iv", "rv", "f.c", "im.a", "im.b", "rm.a", "rm.b")
)
expect_equal(mm[, "f.c"], ifelse(df$f == "c", 1L, 0L))

drop$f <- as.integer(table(df$f))
drop$im <- c(FALSE, TRUE)
mm <- dbarts::makeModelMatrixFromDataFrame(df, drop)
expect_equal(ncol(mm), 8L)
expect_equal(
  colnames(mm),
  c("iv", "rv", "f.a", "f.b", "f.c", "im.a", "rm.a", "rm.b")
)
expect_equal(as.integer(mm[, "im.a"]), as.integer(df$im[, "a"]))

drop$im <- c(FALSE, FALSE)
drop$rm <- c(TRUE, TRUE)
mm <- dbarts::makeModelMatrixFromDataFrame(df, drop)
expect_equal(ncol(mm), 7L)
expect_equal(
  colnames(mm),
  c("iv", "rv", "f.a", "f.b", "f.c", "im.a", "im.b")
)

rm(mm, drop, df)

rm(getTestDataFrame)

# test that make model matrix handles character vectors correctly
n <- 1000L
if (getRversion() >= "3.6.0") {
  suppressWarnings(set.seed(
    0L,
    kind = "Mersenne-Twister",
    normal.kind = "Inversion",
    sample.kind = "Rounding"
  ))
} else {
  set.seed(0L, kind = "Mersenne-Twister", normal.kind = "Inversion")
}
mf <- data.frame(
  x1 = runif(n),
  x2 = c(rep.int(0L, n - 1L), 1L),
  x3 = factor(sample(letters[1L:5L], n, TRUE)),
  x4 = sample(letters[1L:5L], n, TRUE),
  x5 = c("a", rep("b", n - 1L)),
  x6 = c("a", rep("b", n - 2L), "c")
)

mm <- dbarts::makeModelMatrixFromDataFrame(mf)


drop <- attr(mm, "drop")

expect_true(all(!is.null(drop)))
expect_true(all(sapply(drop[sapply(drop, is.numeric)], sum) == n))

factorCols <- which(sapply(mf, function(col) {
  is.factor(col) || is.character(col)
}))

for (col in factorCols) {
  col.table <- table(mf[, col])
  col.name <- colnames(mf)[col]
  col.nvals <- length(col.table)

  expect_true(
    sum(grepl(paste0("^", col.name, "\\."), colnames(mm))) ==
      (if (col.nvals > 2L) col.nvals else col.nvals - 1L)
  )
  expect_true(all(drop[[col.name]] == col.table))
}

rm(col.nvals, col.name, col.table, col)
rm(factorCols, drop, mm, mf, n)

# makeind is the BayesTree-compatible spelling of makeModelMatrixFromDataFrame
# with drop = TRUE; 'all' is documented as not implemented and ignored
df <- data.frame(
  a = c(1.5, 2.5, 3.5, 4.5),
  f = factor(c("x", "y", "z", "x"))
)
expect_identical(
  dbarts::makeind(df),
  dbarts::makeModelMatrixFromDataFrame(df, TRUE)
)
expect_identical(dbarts::makeind(df, all = FALSE), dbarts::makeind(df))

rm(df)


# makeTestModelMatrix aligns a test frame to an existing dbartsData design
set.seed(42)
df <- data.frame(
  a = rnorm(20L),
  f = factor(sample(c("x", "y", "z"), 20L, TRUE)),
  y = rnorm(20L)
)
data <- dbarts::dbartsData(y ~ a + f, df)

newdata <- data.frame(
  a = c(0, 0, 0),
  f = factor(c("z", "y", "x"), levels = levels(df$f))
)
mm <- dbarts::makeTestModelMatrix(data, newdata)

expect_equal(dim(mm), c(3L, ncol(data@x)))
expect_equal(colnames(mm), colnames(data@x))
expect_equal(mm[, "a"], newdata$a)
# categorical predictors carry the training level coding, not the test frame's
expect_equal(mm[, "f"], c(2, 1, 0))

# a test frame in a different column order is matched by name
expect_equal(
  dbarts::makeTestModelMatrix(data, newdata[, c("f", "a")]),
  mm
)

expect_null(dbarts::makeTestModelMatrix(data, NULL))

expect_error(
  dbarts::makeTestModelMatrix(data, data.frame(a = rnorm(3L))),
  pattern = "missing variable required by the model"
)
expect_error(
  dbarts::makeTestModelMatrix(
    data,
    data.frame(a = rnorm(3L), f = factor(c("x", "y", "w")))
  ),
  pattern = "levels not present in the training data"
)

rm(mm, newdata, data, df)


# A missing factor level expands to a missing indicator value on every column
# of that factor's block - the same value the unexpanded categorical coding
# gives it - so the frame's missingness survives the expansion rather than
# reading as "none of the levels".
set.seed(42)
df <- data.frame(
  f = factor(c("a", "b", "c", NA, "a", "b", NA, "c", "a", "b")),
  z = rnorm(10L)
)
mm <- dbarts::makeModelMatrixFromDataFrame(df)
expect_equal(colnames(mm), c("f.a", "f.b", "f.c", "z"))
expect_equal(which(is.na(mm[, "f.a"])), c(4L, 7L))
expect_equal(mm[, "f.a"], ifelse(is.na(df$f), NA_real_, (df$f == "a") * 1))
expect_equal(mm[, "f.b"], ifelse(is.na(df$f), NA_real_, (df$f == "b") * 1))
expect_equal(mm[, "f.c"], ifelse(is.na(df$f), NA_real_, (df$f == "c") * 1))
expect_equal(mm[, "z"], df$z)

# the level counts that pick the retained columns count observed levels only
expect_equal(attr(mm, "drop")$f, c(3L, 3L, 2L))

# a two-level factor keeps its single indicator, missing where the code is
mm <- dbarts::makeModelMatrixFromDataFrame(
  data.frame(g = factor(c("a", "b", NA, "a", "b", "a")))
)
expect_equal(colnames(mm), "g.b")
expect_equal(as.vector(mm), c(0, 1, NA, 0, 1, 0))

# an entirely missing factor observes no level and contributes no column
mm <- dbarts::makeModelMatrixFromDataFrame(
  data.frame(g = factor(c(NA, NA, NA), levels = c("a", "b")), z = rnorm(3L))
)
expect_equal(colnames(mm), "z")

# the expansion drives the same missingness policy the unexpanded coding does:
# bart() rejects it by name, and "incorporate" models it
set.seed(6)
n <- 40L
df <- data.frame(
  f = factor(sample(c("a", "b", "c", NA), n, TRUE)),
  z = rnorm(n)
)
expect_error(
  dbarts::bart(
    df,
    rnorm(n),
    ndpost = 10L,
    nskip = 5L,
    ntree = 5L,
    nthread = 1L,
    verbose = FALSE
  ),
  pattern = "predictors contain missing values"
)
data <- dbarts::dbartsData(
  rnorm(n) ~ f + z,
  df,
  factors = "indicators",
  missing = "incorporate"
)
expect_true(anyNA(data@x))
expect_equal(ncol(data@x), 4L)

# a factor code outside its own level table is refused by column, rather than
# counted into the per-level table it does not index
df <- data.frame(z = rnorm(4L))
df$f <- structure(c(1L, 2L, 7L, 1L), levels = c("a", "b"), class = "factor")
expect_error(
  dbarts::makeModelMatrixFromDataFrame(df),
  pattern = "factor column 'f' has a level code outside its level table"
)
names(df) <- NULL
expect_error(
  dbarts::makeModelMatrixFromDataFrame(as.data.frame(df)),
  pattern = "factor column 2 has a level code outside its level table"
)

# a replayed drop pattern is sized by the table the training column carried,
# and is indexed positionally: a column whose table is longer is refused at
# the entry point rather than read past the end of the pattern. This is the
# builder's own guard, independent of the level comparison validateXTest makes
# on the same route.
set.seed(2718)
train <- data.frame(z = rnorm(20L), f = factor(rep_len(c("a", "b", "c"), 20L)))
drop <- attr(dbarts::makeModelMatrixFromDataFrame(train, TRUE), "drop")
test <- data.frame(
  z = rnorm(5L),
  f = factor(c("a", "b", "c", "a", "b"), levels = c("a", "b", "c", "d"))
)
expect_error(
  dbarts::makeModelMatrixFromDataFrame(test, drop),
  pattern = "factor column 'f' has 4 levels but its drop pattern was built for 3"
)
names(test) <- NULL
expect_error(
  dbarts::makeModelMatrixFromDataFrame(as.data.frame(test), drop),
  pattern = "factor column 2 has 4 levels but its drop pattern was built for 3"
)
# the walk off the end is unbounded in the number of undeclared levels
test <- data.frame(
  z = rnorm(5L),
  f = factor(
    c("a", "b", "c", "a", "b"),
    levels = c("a", "b", "c", paste0("z", 1:5000))
  )
)
expect_error(
  dbarts::makeModelMatrixFromDataFrame(test, drop),
  pattern = "factor column 'f' has 5003 levels"
)
# a matrix column's pattern is one flag per training column, indexed the same way
train <- data.frame(z = rnorm(20L))
train$m <- matrix(rnorm(60L), 20L, 3L, dimnames = list(NULL, c("p", "q", "r")))
drop <- attr(dbarts::makeModelMatrixFromDataFrame(train, TRUE), "drop")
test <- data.frame(z = rnorm(5L))
test$m <- matrix(
  rnorm(20L),
  5L,
  4L,
  dimnames = list(NULL, c("p", "q", "r", "s"))
)
expect_error(
  dbarts::makeModelMatrixFromDataFrame(test, drop),
  pattern = "matrix column 'm' has 4 columns but its drop pattern was built for 3"
)
# the shorter direction indexes in bounds and stays accepted: a training level
# the data never observed contributes no column, so a test frame that lacks it
# expands to exactly the training columns
train <- data.frame(
  f = factor(c("a", "b", "c", "a", "b", "c"), levels = c("a", "b", "c", "d")),
  z = rnorm(6L)
)
drop <- attr(dbarts::makeModelMatrixFromDataFrame(train, TRUE), "drop")
expect_equal(drop$f, c(2L, 2L, 2L, 0L))
test <- data.frame(
  f = factor(c("a", "c"), levels = c("a", "b", "c")),
  z = rnorm(2L)
)
expect_equal(
  colnames(dbarts::makeModelMatrixFromDataFrame(test, drop)),
  c("f.a", "f.b", "f.c", "z")
)

rm(train, test, drop)

# a frame with no columns names its own emptiness
expect_error(
  dbarts::makeModelMatrixFromDataFrame(data.frame()),
  pattern = "at least one column"
)

rm(mm, data, df, n)
