# factors = "categorical": data frame ingestion without dummy expansion,
# split by category subset on the bartcore engine

set.seed(99)
n <- 400L
df <- data.frame(
  g = factor(sample(c("red", "green", "blue"), n, replace = TRUE)),
  o = factor(sample(c("low", "mid", "high"), n, replace = TRUE),
             levels = c("low", "mid", "high"), ordered = TRUE),
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
expect_equal(data.cat@varTypes, c(1L, 0L, 0L, 0L, 1L))
expect_equal(as.double(data.cat@x[, "g"]), as.double(as.integer(df$g) - 1L))
expect_equal(as.double(data.cat@x[, "o"]), as.double(as.integer(df$o) - 1L))
expect_equal(attr(data.cat@x, "factor.levels")[[1L]],
             c("blue", "green", "red"))
expect_equal(attr(data.cat@x, "factor.levels")[[5L]], c("cat", "dog"))
expect_null(attr(data.cat@x, "factor.levels")[[3L]])

# categorical is the default; "indicators" expands as previous versions did
data.def <- dbartsData(y ~ g + o + z + b + s, df)
expect_equal(data.def@varTypes, data.cat@varTypes)
data.ind <- dbartsData(y ~ g + o + z + b + s, df, factors = "indicators")
expect_true(ncol(data.ind@x) > 5L)
expect_true(all(data.ind@varTypes == 0L))

# unordered factors cap at 53 levels; ordered factors are ordinal and do not
df.wide <- data.frame(f = factor(paste0("l", seq_len(54L))), y = rnorm(54L))
expect_error(dbartsData(y ~ f, df.wide, factors = "categorical"),
             pattern = "more than 53 levels")
df.wide$f <- as.ordered(df.wide$f)
expect_inherits(dbartsData(y ~ f, df.wide, factors = "categorical"),
                "dbartsData")

# test data code against the training level tables, whatever their order or
# subset, and unseen levels are an error
df.test <- df[seq_len(20L), c("s", "b", "z", "o", "g")]
df.test$g <- factor(as.character(df.test$g), levels = c("red", "blue", "green"))
data.cat2 <- dbartsData(y ~ g + o + z + b + s, df, test = df.test,
                        factors = "categorical")
expect_equal(as.double(data.cat2@x.test[, "g"]),
             as.double(data.cat@x[seq_len(20L), "g"]))
expect_equal(unname(data.cat2@x.test[, c("o", "z", "b", "s")]),
             unname(data.cat@x[seq_len(20L), c("o", "z", "b", "s")]))
df.bad <- df.test
levels(df.bad$g) <- c("red", "blue", "purple")
expect_error(dbartsData(y ~ g + o + z + b + s, df, test = df.bad,
                        factors = "categorical"),
             pattern = "levels not present")

# the classic engine has no subset splits
expect_error(dbarts(y ~ g + o + z + b + s, df, factors = "categorical",
                    control = dbartsControl(engine = "classic")),
             pattern = "categorical predictors require")

# end to end on the bartcore engine: group means recovered through subset
# splits, and factor test data route through the same coding
control <- dbartsControl(engine = "bartcore", n.chains = 1L, n.threads = 1L,
                         n.trees = 75L, updateState = FALSE)
sampler <- dbarts(y ~ g + o + z + b + s, df, test = df.test,
                  control = control, factors = "categorical")
samples <- sampler$run(150L, 300L)
fitMeans <- rowMeans(samples$train) - df$z
recovered <- tapply(fitMeans, df$g, mean)
expect_true(max(abs(recovered - groupMeans[names(recovered)])) < 0.3)
expect_equal(dim(samples$test), c(20L, 300L))

# a later data.frame through setTestPredictor recodes identically
sampler$setTestPredictor(df.test)
expect_equal(sampler$data@x.test[, "g"], data.cat2@x.test[, "g"])
