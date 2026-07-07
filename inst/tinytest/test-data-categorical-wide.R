# wide categorical predictors: up to 63 levels keep inline rule masks (inline
# in the flattened node too), wider columns pool their mask words per tree and
# the flattened format references them through a side channel; getTrees reports
# every categorical rule decoded (docs/design/pooled-masks.md)

set.seed(4247)
n <- 700L
levels.g <- paste0("g", sprintf("%02d", seq_len(70L)))
levels.h <- paste0("h", sprintf("%02d", seq_len(60L)))
df <- data.frame(
  g = factor(sample(levels.g, n, replace = TRUE), levels = levels.g),
  h = factor(sample(levels.h, n, replace = TRUE), levels = levels.h),
  z = runif(n)
)
signal.g <- ifelse(seq_len(70L) %% 3L == 0L, 2, 0)
df$y <- signal.g[as.integer(df$g)] + 0.5 * df$z + rnorm(n, 0, 0.25)

# ingestion types both wide factors categorical with their full level tables
data.wide <- dbartsData(y ~ g + h + z, df)
expect_equal(data.wide@varTypes, c(1L, 1L, 0L))
expect_equal(length(attr(data.wide@x, "factor.levels")[[1L]]), 70L)

control <- dbartsControl(
  n.chains = 2L,
  n.threads = 1L,
  n.trees = 50L,
  n.samples = 60L,
  keepTrees = TRUE,
  updateState = FALSE
)
sampler <- dbarts(
  y ~ g + h + z,
  df,
  test = df[seq_len(50L), ],
  control = control,
  factors = "categorical"
)
samples <- sampler$run(400L, 60L)

# subset splits over 70 levels recover the signal
fits <- rowMeans(samples$train[,, 1L]) - 0.5 * df$z
contrast <- mean(fits[as.integer(df$g) %% 3L == 0L]) -
  mean(fits[as.integer(df$g) %% 3L != 0L])
expect_true(abs(contrast - 2) < 0.4)

# getTrees: a wide rule's value is NA and its directions string covers the
# declared levels; both wide tiers decode, one L/R per level
trees <- sampler$getTrees()
expect_true("directions" %in% names(trees))
isRule.g <- trees$var == 1L
isRule.h <- trees$var == 2L
expect_true(any(isRule.g) && any(isRule.h))
expect_true(all(is.na(trees$value[isRule.g])))
expect_true(all(is.na(trees$value[isRule.h])))
expect_true(all(nchar(trees$directions[isRule.g]) == 70L))
expect_true(all(nchar(trees$directions[isRule.h]) == 60L))
expect_true(all(is.na(trees$directions[trees$var %in% c(-1L, 3L)])))
# every mask leaves at least one observed level on each side
expect_true(all(
  grepl("L", trees$directions[isRule.g]) &
    grepl("R", trees$directions[isRule.g])
))

# the decode matches the engine's routing: at a categorical root split the
# left child's n counts the training rows whose level decodes to L
roots.g <- which(isRule.g & !duplicated(trees[c("chain", "sample", "tree")]))
for (i in roots.g[seq_len(min(10L, length(roots.g)))]) {
  goesLeft <- strsplit(trees$directions[i], "")[[1L]] == "L"
  codes <- sampler$data@x[, 1L]
  expect_equal(trees$n[i + 1L], sum(goesLeft[codes + 1L]))
}

# saved-tree replay reproduces the recorded test fits
predictions <- sampler$predict(df[seq_len(50L), c("g", "h", "z")])
expect_equal(dim(predictions), dim(samples$test))
expect_true(max(abs(predictions - samples$test)) < 1e-10)

source(system.file("common", "stateContinuation.R", package = "dbarts"))
# the state round-trips, mask channels included
invisible(sampler$storeState())
state <- sampler$state
expect_true(all(sapply(state, function(chain) {
  "tree.masks" %in% names(chain) && length(chain[["tree.masks"]]) > 0L
})))
sampler2 <- dbarts(
  y ~ g + h + z,
  df,
  test = df[seq_len(50L), ],
  control = control,
  factors = "categorical"
)
sampler2$setState(state)
sampler2$storeState()
expect_true(statesAgree(sampler2$state, state))
for (ci in seq_along(state)) {
  expect_identical(sampler2$state[[ci]]$tree.masks, state[[ci]]$tree.masks)
}

pdf(NULL)
expect_silent(sampler$plotTree(1L, chainNum = 1L))
dev.off()

# missing values compose: the NA pseudo-category routes like any level
df.na <- df
df.na$g[seq_len(70L)] <- NA
sampler.na <- dbarts(
  y ~ g + h + z,
  df.na,
  control = control,
  factors = "categorical"
)
samples.na <- sampler.na$run(50L, 10L)
expect_true(all(is.finite(samples.na$sigma)))
trees.na <- sampler.na$getTrees()
expect_true("missing" %in% names(trees.na))
naRules.g <- trees.na$var == 1L & !is.na(trees.na$missing)
expect_true(any(naRules.g))
expect_true(all(trees.na$missing[naRules.g] %in% c("L", "R")))

rm(
  sampler,
  sampler2,
  sampler.na,
  samples,
  samples.na,
  state,
  trees,
  trees.na,
  predictions,
  fits,
  contrast,
  roots.g,
  isRule.g,
  isRule.h,
  naRules.g,
  data.wide,
  control,
  df,
  df.na,
  levels.g,
  levels.h,
  signal.g,
  n,
  i,
  goesLeft,
  codes
)
