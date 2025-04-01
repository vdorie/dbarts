getTestDataFrame <- function() {
  set.seed(42)
  df <- data.frame(iv = seq.int(10L), rv = runif(10),
                   f = factor(sample(3L, 10L, TRUE), labels = c("a", "b", "c")))
  df[[4L]] <- matrix(rbinom(20L, 3L, 0.5), 10L, dimnames = list(NULL, c("a", "b")))
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
expect_equal(mm[,"iv"], df$iv)
expect_equal(mm[,"rv"], df$rv)
expect_equal(mm[,"f.a"], ifelse(df$f == "a", 1L, 0L))
expect_equal(mm[,"f.b"], ifelse(df$f == "b", 1L, 0L))
expect_equal(mm[,"f.c"], ifelse(df$f == "c", 1L, 0L))
expect_equal(mm[,"im.a"], df$im[,"a"])
expect_equal(mm[,"im.b"], df$im[,"b"])
expect_equal(mm[,"rm.a"], df$rm[,"a"])
expect_equal(mm[,"rm.b"], df$rm[,"b"])

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
df$im[,1L] <- rep(1L, 10L)
mm <- dbarts::makeModelMatrixFromDataFrame(df)
expect_equal(ncol(mm), 8L)
expect_equal(attr(mm, "drop")$im, c(TRUE, FALSE))
expect_equal(
  colnames(mm),
  c("iv", "rv", "f.a", "f.b", "f.c", "im.b", "rm.a", "rm.b")
)
expect_equal(as.double(mm[,7L:8L]), as.double(df$rm))

df <- getTestDataFrame()
df$im[,1L] <- rep(1L, 10L); df$im[,2L] <- rep(2L, 10L)
mm <- dbarts::makeModelMatrixFromDataFrame(df)
expect_equal(ncol(mm), 7L)
expect_equal(attr(mm, "drop")$im, c(TRUE, TRUE))
expect_equal(
  colnames(mm),
  c("iv", "rv", "f.a", "f.b", "f.c", "rm.a", "rm.b")
)

df <- getTestDataFrame()
df$rm[,2L] <- rep(pi, 10L)
mm <- dbarts::makeModelMatrixFromDataFrame(df)
expect_equal(ncol(mm), 8L)
expect_equal(attr(mm, "drop")$rm, c(FALSE, TRUE))
expect_equal(
  colnames(mm),
  c("iv", "rv", "f.a", "f.b", "f.c", "im.a", "im.b", "rm.a")
)
expect_equal(as.integer(mm[,6L:7L]), as.integer(df$im))

rm(mm, df)


# test that make model matrix doesn't drop useless columns when drop = FALSE
df <- getTestDataFrame()
df$iv <- rep(1L, 10L)
mm <- dbarts::makeModelMatrixFromDataFrame(df, FALSE)
expect_equal(ncol(mm), 9L)
expect_equal(mm[,"iv"], df$iv)

df <- getTestDataFrame()
df$rv <- rep(pi, 10L)
mm <- dbarts::makeModelMatrixFromDataFrame(df, FALSE)
expect_equal(mm[,"rv"], df$rv)

df <- getTestDataFrame()
df$f <- factor(
  rep(seq.int(3L), c(5L, 5L, 1L)),
  labels = c("a", "b", "c")
)[1L:10L]
mm <- dbarts::makeModelMatrixFromDataFrame(df, FALSE)
expect_equal(mm[,"f.a"], c(rep(1L, 5), rep(0L, 5)))
expect_equal(mm[,"f.b"], c(rep(0L, 5), rep(1L, 5)))
expect_equal(mm[,"f.c"], rep(0L, 10))

df <- getTestDataFrame()
df$im[,1L] <- rep(1L, 10L)
mm <- dbarts::makeModelMatrixFromDataFrame(df, FALSE)
expect_equal(as.integer(mm[,6L:7L]), as.integer(df$im))

df <- getTestDataFrame()
df$rm[,2L] <- rep(pi, 10L)
mm <- dbarts::makeModelMatrixFromDataFrame(df, FALSE)
expect_equal(as.double(mm[,8L:9L]), as.double(df$rm))

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

drop$iv <- FALSE; drop$rv <- TRUE
mm <- dbarts::makeModelMatrixFromDataFrame(df, drop)
expect_equal(ncol(mm), 8L)
expect_equal(
  colnames(mm),
  c("iv", "f.a", "f.b", "f.c", "im.a", "im.b", "rm.a", "rm.b")
)

drop$rv <- FALSE; drop$f <- c(1L, 0L, 1L)
mm <- dbarts::makeModelMatrixFromDataFrame(df, drop)
expect_equal(ncol(mm), 7L)
expect_equal(
  colnames(mm),
  c("iv", "rv", "f.c", "im.a", "im.b", "rm.a", "rm.b")
)
expect_equal(mm[,"f.c"], ifelse(df$f == "c", 1L, 0L))

drop$f <- as.integer(table(df$f)); drop$im <- c(FALSE, TRUE)
mm <- dbarts::makeModelMatrixFromDataFrame(df, drop)
expect_equal(ncol(mm), 8L)
expect_equal(
  colnames(mm),
  c("iv", "rv", "f.a", "f.b", "f.c", "im.a", "rm.a", "rm.b")
)
expect_equal(as.integer(mm[,"im.a"]), as.integer(df$im[,"a"]))

drop$im <- c(FALSE, FALSE); drop$rm <- c(TRUE, TRUE)
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

factorCols <- which(sapply(mf, function(col) is.factor(col) || is.character(col)))

for (col in factorCols) {
  col.table <- table(mf[,col])
  col.name  <- colnames(mf)[col]
  col.nvals <- length(col.table)
  
  expect_true(
    sum(grepl(paste0("^", col.name, "\\."),
    colnames(mm))) == (if (col.nvals > 2L) col.nvals else col.nvals - 1L)
  )
  expect_true(all(drop[[col.name]] == col.table))
}

rm(col.nvals, col.name, col.table, col)
rm(factorCols, drop, mm, mf, n)

