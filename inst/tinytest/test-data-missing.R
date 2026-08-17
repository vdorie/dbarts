# missing = "incorporate": NAs in the predictors are modeled, every split
# rule learning a direction for missing values; "error" rejects them

set.seed(42)
n <- 400L
x1 <- runif(n)
x2 <- runif(n)
isMissing <- runif(n) < 0.3
x1[isMissing] <- NA_real_
g <- factor(sample(c("a", "b", "c"), n, replace = TRUE))
g[runif(n) < 0.15] <- NA
y <- 2 * isMissing + x2 + 0.5 * (!is.na(g) & g == "b") + rnorm(n, 0, 0.25)
df <- data.frame(x1 = x1, x2 = x2, g = g, y = y)

# ingestion keeps incomplete rows (previous versions na.omit-dropped them
# for formula inputs) and codes factor NAs alongside level codes
data.mia <- dbartsData(y ~ x1 + x2 + g, df)
expect_equal(length(data.mia@y), n)
expect_equal(nrow(data.mia@x), n)
expect_true(anyNA(data.mia@x[, "x1"]) && anyNA(data.mia@x[, "g"]))
expect_equal(data.mia@missing, "incorporate")

# the error escape rejects incomplete predictors; the response side always
# must be complete; a column of nothing but NA has no observed values
expect_error(
  dbartsData(y ~ x1 + x2 + g, df, missing = "error"),
  pattern = "missing values"
)
expect_inherits(dbartsData(y ~ x2, df, missing = "error"), "dbartsData")
df.badY <- df
df.badY$y[1L] <- NA
expect_error(
  dbartsData(y ~ x1 + x2 + g, df.badY),
  pattern = "response contains missing"
)
df.allNA <- df
df.allNA$x2 <- NA_real_
expect_error(
  dbartsData(y ~ x1 + x2 + g, df.allNA),
  pattern = "entirely missing"
)

# end to end: the missingness signal is recovered, and test rows with NAs
# route through the learned directions
test.df <- df[seq_len(20L), c("x1", "x2", "g")]
control <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 50L,
  updateState = FALSE
)
sampler <- dbarts(y ~ x1 + x2 + g, df, test = test.df, control = control)
samples <- sampler$run(200L, 200L)
fits <- rowMeans(samples$train)
expect_true(mean(fits[isMissing]) - mean(fits[!isMissing]) > 1)
expect_equal(dim(samples$test), c(20L, 200L))
expect_true(!anyNA(samples$test))

# getTrees decodes each rule's missing route into a 'missing' column: L/R
# on columns that contain NAs (x1, g), NA on complete columns and leaves
control.keep <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 25L,
  n.samples = 5L,
  keepTrees = TRUE,
  updateState = FALSE
)
sampler.keep <- dbarts(y ~ x1 + x2 + g, df, control = control.keep)
invisible(sampler.keep$run(100L, 5L))
trees <- sampler.keep$getTrees()
expect_true("missing" %in% names(trees))
isRule <- trees$var > 0L
expect_true(all(is.na(trees$missing[!isRule])))
onMissingColumn <- isRule & trees$var %in% c(1L, 3L)
expect_true(any(onMissingColumn))
expect_true(all(trees$missing[onMissingColumn] %in% c("L", "R")))
expect_true(all(is.na(trees$missing[isRule & trees$var == 2L])))

# at a root split on x1, the left child's n counts the rows the rule and
# its missing route send left
roots <- which(!duplicated(trees[c("sample", "tree")]) & trees$var == 1L)
codes.x1 <- df$x1
for (i in roots) {
  goesLeft <- ifelse(
    is.na(codes.x1),
    trees$missing[i] == "L",
    codes.x1 <= trees$value[i]
  )
  expect_equal(trees$n[i + 1L], sum(goesLeft))
}

pdf(NULL)
expect_silent(sampler.keep$plotTree(1L))
dev.off()

# predict accepts NAs in new data
predictions <- sampler.keep$predict(test.df)
expect_equal(dim(predictions), c(20L, 5L))
expect_true(!anyNA(predictions))

# a sampler built with missing = "error" refuses NA updates
sampler.strict <- dbarts(y ~ x2, df, control = control, missing = "error")
x2.bad <- df$x2
x2.bad[1L] <- NA
expect_error(
  sampler.strict$setPredictor(x2.bad, "x2"),
  pattern = "missing = \"error\""
)

# on an incorporate sampler setPredictor may MOVE a column's NA pattern:
# install a replacement whose missingness carries the response signal, and the
# refit picks the new pattern up (it was unusable before), fits stay finite,
# and the installed column carries the moved NAs
set.seed(11)
m1.mv <- runif(n) < 0.3 # initial NA pattern, unrelated to the signal
m2.mv <- runif(n) < 0.3 # moved-in NA pattern, carries the signal
y.mv <- 3 * m2.mv + x2 + rnorm(n, 0, 0.25)
x1.mv <- runif(n)
x1.mv[m1.mv] <- NA_real_
x1.moved <- runif(n)
x1.moved[m2.mv] <- NA_real_
df.mv <- data.frame(x1 = x1.mv, x2 = x2, y = y.mv)
control.mv <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 50L,
  updateState = FALSE,
  seed = 3L
)
sampler.mv <- dbarts(y.mv ~ x1 + x2, df.mv, control = control.mv)
fits.before <- rowMeans(sampler.mv$run(200L, 200L)$train)
sampler.mv$setPredictor(x1.moved, "x1")
run.after <- sampler.mv$run(200L, 200L)
fits.after <- rowMeans(run.after$train)
expect_true(mean(fits.before[m2.mv]) - mean(fits.before[!m2.mv]) < 1)
expect_true(mean(fits.after[m2.mv]) - mean(fits.after[!m2.mv]) > 2)
expect_true(all(is.finite(run.after$train)))
expect_equal(which(is.na(sampler.mv$data@x[, "x1"])), which(m2.mv))

source(
  system.file("common", "stateContinuation.R", package = "dbarts"),
  local = TRUE
)
# state serialization carries the missing directions: a restored sampler
# reproduces the model
control.state <- dbartsControl(
  n.chains = 2L,
  n.threads = 1L,
  n.trees = 10L,
  n.samples = 5L,
  updateState = FALSE
)
sampler.state <- dbarts(y ~ x1 + x2 + g, df, control = control.state)
invisible(sampler.state$run(30L, 2L))
sampler.state$storeState()
sampler.restored <- dbarts(y ~ x1 + x2 + g, df, control = control.state)
sampler.restored$setState(sampler.state$state)
sampler.restored$storeState()
statesAgree(sampler.restored$state, sampler.state$state)

# xbart's folds gather the reserved codes through the data handle
xval <- xbart(
  y ~ x1 + x2 + g,
  df,
  method = "k-fold",
  n.test = 5,
  n.reps = 2L,
  n.samples = 6L,
  n.burn = c(5L, 3L, 1L),
  n.trees = 5L,
  n.threads = 1L
)
expect_true(!anyNA(xval))

rm(
  sampler,
  sampler.keep,
  sampler.strict,
  sampler.state,
  sampler.restored,
  trees,
  samples,
  predictions,
  xval,
  df,
  test.df,
  df.badY,
  df.allNA,
  data.mia,
  x1,
  x2,
  g,
  y,
  isMissing,
  codes.x1,
  roots,
  fits,
  control,
  control.keep,
  control.state,
  onMissingColumn,
  isRule,
  n,
  sampler.mv,
  control.mv,
  df.mv,
  m1.mv,
  m2.mv,
  y.mv,
  x1.mv,
  x1.moved,
  fits.before,
  fits.after,
  run.after,
  x2.bad
)
