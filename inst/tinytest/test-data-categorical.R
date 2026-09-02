# factors = "categorical": data frame ingestion without dummy expansion,
# split by category subset

set.seed(99)
n <- 400L
df <- data.frame(
  g = factor(sample(c("red", "green", "blue"), n, replace = TRUE)),
  o = factor(
    sample(c("low", "mid", "high"), n, replace = TRUE),
    levels = c("low", "mid", "high"),
    ordered = TRUE
  ),
  z = runif(n),
  b = sample(c(TRUE, FALSE), n, replace = TRUE),
  s = sample(c("cat", "dog"), n, replace = TRUE),
  stringsAsFactors = FALSE
)
groupMeans <- c(blue = 2, green = -1, red = 0)
df$y <- groupMeans[as.character(df$g)] + df$z + rnorm(n, 0, 0.4)

# structure: one column per predictor, coded and typed
data.cat <- dbartsData(y ~ g + o + z + b + s, df, factors = "categorical")
expect_equal(ncol(data.cat@x), 5L)
expect_equal(colnames(data.cat@x), c("g", "o", "z", "b", "s"))
expect_equal(data.cat@varTypes, c(1L, 2L, 0L, 0L, 1L))
expect_equal(as.double(data.cat@x[, "g"]), as.double(as.integer(df$g) - 1L))
expect_equal(as.double(data.cat@x[, "o"]), as.double(as.integer(df$o) - 1L))
expect_equal(attr(data.cat@x, "factor.levels")[[1L]], c("blue", "green", "red"))
expect_equal(attr(data.cat@x, "factor.levels")[[5L]], c("cat", "dog"))
expect_null(attr(data.cat@x, "factor.levels")[[3L]])

# categorical is the default; "indicators" expands as previous versions did
data.def <- dbartsData(y ~ g + o + z + b + s, df)
expect_equal(data.def@varTypes, data.cat@varTypes)
data.ind <- dbartsData(y ~ g + o + z + b + s, df, factors = "indicators")
expect_true(ncol(data.ind@x) > 5L)
expect_true(all(data.ind@varTypes == 0L))

# unordered factors cap at the code type's limit (see
# test-data-categorical-wide.R for levels past 53); ordered factors are
# ordinal and never cap
df.wide <- data.frame(f = factor(paste0("l", seq_len(54L))), y = rnorm(54L))
expect_inherits(
  dbartsData(y ~ f, df.wide, factors = "categorical"),
  "dbartsData"
)
df.wide$f <- as.ordered(df.wide$f)
expect_inherits(
  dbartsData(y ~ f, df.wide, factors = "categorical"),
  "dbartsData"
)

# test data code against the training level tables, whatever their order or
# subset, and unseen levels are an error
df.test <- df[seq_len(20L), c("s", "b", "z", "o", "g")]
df.test$g <- factor(as.character(df.test$g), levels = c("red", "blue", "green"))
data.cat2 <- dbartsData(
  y ~ g + o + z + b + s,
  df,
  test = df.test,
  factors = "categorical"
)
expect_equal(
  as.double(data.cat2@x.test[, "g"]),
  as.double(data.cat@x[seq_len(20L), "g"])
)
expect_equal(
  unname(data.cat2@x.test[, c("o", "z", "b", "s")]),
  unname(data.cat@x[seq_len(20L), c("o", "z", "b", "s")])
)
df.bad <- df.test
levels(df.bad$g) <- c("red", "blue", "purple")
expect_error(
  dbartsData(y ~ g + o + z + b + s, df, test = df.bad, factors = "categorical"),
  pattern = "levels not present"
)

# end to end: group means recovered through subset splits, and factor test
# data route through the same coding
control <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 75L,
  updateState = FALSE
)
sampler <- dbarts(
  y ~ g + o + z + b + s,
  df,
  test = df.test,
  control = control,
  factors = "categorical"
)
samples <- sampler$run(150L, 300L)
fitMeans <- rowMeans(samples$train) - df$z
recovered <- tapply(fitMeans, df$g, mean)
expect_true(max(abs(recovered - groupMeans[names(recovered)])) < 0.3)
expect_equal(dim(samples$test), c(20L, 300L))

# a later data.frame through setTestPredictor recodes identically
sampler$setTestPredictor(df.test)
expect_equal(sampler$data@x.test[, "g"], data.cat2@x.test[, "g"])

# getTrees keeps the raw direction mask in 'value' and decodes it in a
# 'directions' column: one L/R per level in level order, bit k - 1 of the
# mask set sending level k right; ordinal rules and leaves are NA
control.keep <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 25L,
  n.samples = 5L,
  keepTrees = TRUE,
  updateState = FALSE
)
sampler.keep <- dbarts(
  y ~ g + o + z + b + s,
  df,
  control = control.keep,
  factors = "categorical"
)
invisible(sampler.keep$run(50L, 5L))
trees <- sampler.keep$getTrees()
expect_true("directions" %in% names(trees))

varTypes <- sampler.keep$data@varTypes
isCategoricalRule <- trees$var > 0L & varTypes[pmax(trees$var, 1L)] == 1L
expect_true(any(isCategoricalRule))
expect_true(all(is.na(trees$directions[!isCategoricalRule])))
expect_true(!anyNA(trees$directions[isCategoricalRule]))

numLevels <- lengths(attr(sampler.keep$data@x, "factor.levels"))
expect_equal(
  nchar(trees$directions[isCategoricalRule]),
  unname(numLevels[trees$var[isCategoricalRule]])
)
# a categorical rule reports its split in directions; value is NA
expect_true(all(is.na(trees$value[isCategoricalRule])))
# masks put at least one level on each side
expect_true(all(
  grepl("L", trees$directions[isCategoricalRule]) &
    grepl("R", trees$directions[isCategoricalRule])
))

# the decode matches what the engine does: at a categorical root split the
# left child's n counts the training observations whose level decodes to L
categoricalRoots <- which(
  isCategoricalRule &
    !duplicated(trees[c("sample", "tree")])
)
for (i in categoricalRoots) {
  goesLeft <- strsplit(trees$directions[i], "")[[1L]] == "L"
  codes <- sampler.keep$data@x[, trees$var[i]]
  expect_equal(trees$n[i + 1L], sum(goesLeft[codes + 1L]))
}

# plotTree labels categorical rules with the left-branch level set
pdf(NULL)
expect_silent(sampler.keep$plotTree(1L))
dev.off()

# a hand-typed column of either factor kind is bounded before ingestion: the
# store derives a factor column's level count by casting its observed maximum
# to an unsigned count whichever kind it is, so a value outside the code range
# is refused at creation rather than converted
n.codes <- 40L
x.codes <- matrix(as.double(rep(c(0, 1, 2, 3), 10L)), ncol = 1L)
x.codes[1L] <- 1e13
y.codes <- rnorm(n.codes)
control.codes <- dbartsControl(
  n.trees = 1L,
  n.chains = 1L,
  n.threads = 1L,
  updateState = FALSE
)
attr(x.codes, "varTypes") <- 2L # ordered factor
expect_error(
  dbarts(x.codes, y.codes, control = control.codes),
  pattern = "ordered factor predictors must hold integer level codes"
)
attr(x.codes, "varTypes") <- 1L # categorical
expect_error(
  dbarts(x.codes, y.codes, control = control.codes),
  pattern = "categorical predictors must hold integer category codes"
)
# and a column whose codes ARE in range is accepted under either kind
x.codes[1L] <- 3
attr(x.codes, "varTypes") <- 2L # ordered factor
expect_inherits(
  dbarts(x.codes, y.codes, control = control.codes),
  "dbartsSampler"
)
