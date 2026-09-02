# A data frame's factor columns reach the engine as their own integer codes,
# and a plain matrix's reach it as doubles. The two channels must produce the
# same fit: a dbartsData built before the code channel existed stores its
# factor columns as a double matrix with varTypes and a level table, which is
# exactly what the matrix arm below is, so this is also the pin that such an
# object still loads and fits as it did.

set.seed(99L)
n.chan <- 60L
f.chan <- factor(
  sample(c("a", "b", "c"), n.chan, replace = TRUE),
  levels = c("a", "b", "c", "d") # "d" is declared but never observed
)
f.chan[c(3L, 17L)] <- NA
o.chan <- ordered(
  sample(seq_len(4L), n.chan, replace = TRUE),
  levels = seq_len(4L)
)
o.chan[c(5L, 41L)] <- NA
z.chan <- rnorm(n.chan)
y.chan <- rnorm(n.chan)
df.chan <- data.frame(f = f.chan, o = o.chan, z = z.chan)

n.test.chan <- 10L
f.test.chan <- factor(
  rep_len(c("a", "b"), n.test.chan),
  levels = levels(f.chan)
)
o.test.chan <- ordered(rep_len(c(1L, 3L), n.test.chan), levels = seq_len(4L))
df.test.chan <- data.frame(
  f = f.test.chan,
  o = o.test.chan,
  z = rnorm(n.test.chan)
)

# the same data as the double matrix an older data object stored: 0-based
# codes, the kinds on varTypes, and the level tables the declared counts come
# from
asCodes <- function(column) as.double(as.integer(column)) - 1
x.chan <- cbind(asCodes(f.chan), asCodes(o.chan), z.chan)
colnames(x.chan) <- c("f", "o", "z")
attr(x.chan, "varTypes") <- c(1L, 2L, 0L)
attr(x.chan, "factor.levels") <- list(levels(f.chan), levels(o.chan), NULL)
x.test.chan <- cbind(
  asCodes(f.test.chan),
  asCodes(o.test.chan),
  df.test.chan$z
)
colnames(x.test.chan) <- c("f", "o", "z")
attr(x.test.chan, "varTypes") <- c(1L, 2L, 0L)

# the data objects agree on the kinds and on the level tables the grids come
# from, which is what makes the fits below comparable at all
data.frame.chan <- dbarts::dbartsData(df.chan, y.chan, df.test.chan)
data.matrix.chan <- dbarts::dbartsData(x.chan, y.chan, x.test.chan)
expect_equal(data.frame.chan@varTypes, data.matrix.chan@varTypes)
expect_equal(data.frame.chan@varTypes, c(1L, 2L, 0L))

control.chan <- dbarts::dbartsControl(
  n.trees = 5L,
  n.chains = 1L,
  n.threads = 1L,
  n.burn = 0L,
  updateState = FALSE
)

set.seed(7L)
sampler.frame.chan <- dbarts::dbarts(
  df.chan,
  y.chan,
  df.test.chan,
  control = control.chan
)
fit.frame.chan <- sampler.frame.chan$run(10L, 20L)

set.seed(7L)
sampler.matrix.chan <- dbarts::dbarts(
  x.chan,
  y.chan,
  x.test.chan,
  control = control.chan
)
fit.matrix.chan <- sampler.matrix.chan$run(10L, 20L)

# bitwise: the channel changes how the values cross the bridge, nothing about
# what they are
expect_identical(fit.frame.chan$train, fit.matrix.chan$train)
expect_identical(fit.frame.chan$test, fit.matrix.chan$test)
expect_identical(fit.frame.chan$sigma, fit.matrix.chan$sigma)
expect_identical(fit.frame.chan$varcount, fit.matrix.chan$varcount)

# the NA-bearing factor columns are flagged as missing on both channels, which
# is what makes their rules draw a missing direction: a column that lost the
# flag would consume no such draw and shift every draw after it. The identical
# fits above are the evidence; this asserts the data reached the samplers with
# its NAs intact.
expect_true(anyNA(data.frame.chan@x$dense[[1L]]))
expect_true(anyNA(data.matrix.chan@x[, 1L]))

# and a level declared but never observed still counts, on both channels: the
# test rows may carry it
expect_equal(nlevels(f.chan), 4L)

# A refused whole-data replacement leaves the sampler exactly as it was: the
# level check runs on both sides before either is installed, so a sampler that
# saw the refusal draws bitwise identically to one that never saw the call.
x.bad.chan <- x.chan
x.bad.chan[1L, 1L] <- 4 # one level past the four "f" declares
attr(x.bad.chan, "varTypes") <- c(1L, 2L, 0L)
attr(x.bad.chan, "factor.levels") <- attr(x.chan, "factor.levels")
data.bad.chan <- dbarts::dbartsData(x.bad.chan, y.chan)

set.seed(11L)
sampler.kept.chan <- dbarts::dbarts(df.chan, y.chan, control = control.chan)
set.seed(11L)
sampler.refused.chan <- dbarts::dbarts(df.chan, y.chan, control = control.chan)
expect_error(
  sampler.refused.chan$setData(data.bad.chan),
  pattern = "categorical predictor values must be existing category codes"
)
expect_identical(
  sampler.refused.chan$run(10L, 20L)$train,
  sampler.kept.chan$run(10L, 20L)$train
)

rm(
  n.chan,
  f.chan,
  o.chan,
  z.chan,
  y.chan,
  df.chan,
  n.test.chan,
  f.test.chan,
  o.test.chan,
  df.test.chan,
  asCodes,
  x.chan,
  x.test.chan,
  data.frame.chan,
  data.matrix.chan,
  control.chan,
  sampler.frame.chan,
  fit.frame.chan,
  sampler.matrix.chan,
  fit.matrix.chan,
  x.bad.chan,
  data.bad.chan,
  sampler.kept.chan,
  sampler.refused.chan
)
