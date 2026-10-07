# updatePredictorPerObservationJointly given a factor column: labels are matched to the column's
# levels by name, as setPredictor by column matches them; numbers are the codes from 0 the
# samplers hold the column in.

lv <- c(f = list(letters[1:4]), o = list(c("lo", "mid", "hi")))

makeData <- function(n, seed) {
  set.seed(seed)
  d <- data.frame(
    x1 = runif(n),
    f = factor(sample(lv$f, n, TRUE), levels = lv$f),
    o = factor(sample(lv$o, n, TRUE), levels = lv$o, ordered = TRUE)
  )
  d$y <- d$x1 + as.integer(d$f) + rnorm(n, sd = 0.3)
  d
}

# a sampler on formula; `form` fixes the column order, so two samplers can hold f and o at
# different positions
makeSampler <- function(d, seed, trees = 5L, form = y ~ x1 + f + o, burn = 0L) {
  ctl <- dbarts::dbartsControl(
    n.trees = trees,
    n.chains = 1L,
    n.threads = 1L,
    updateState = FALSE,
    verbose = FALSE
  )
  s <- dbarts::dbarts(form, d, control = ctl, seed = seed)
  if (burn > 0L) {
    invisible(s$run(burn, 1L))
  }
  s
}

# what a sampler holds, byte for byte
snapshot <- function(s) {
  s$storeState()
  serialize(list(unclass(s$state), s$data@x), NULL)
}

# the next level of a column, as labels
nextLabels <- function(x, levels) {
  levels[(match(as.character(x), levels) %% length(levels)) + 1L]
}

# every sampler's column, level by level: rows installed hold the label given, the others what
# they held before
expectInstalled <- function(
  samplers,
  name,
  before,
  labels,
  levels,
  mask,
  info
) {
  for (i in seq_along(samplers)) {
    held <- levels[samplers[[i]]$data@x[, name] + 1L]
    want <- ifelse(mask, as.character(labels), before[[i]])
    expect_identical(held, want, info = paste(info, "sampler", i))
  }
}

heldLabels <- function(s, name, levels) levels[s$data@x[, name] + 1L]

d <- makeData(200L, 1L)
n <- nrow(d)
set.seed(11)
newF <- sample(lv$f, n, TRUE)
newO <- sample(lv$o, n, TRUE)

kinds <- function(labels, levels, name) {
  ordered <- name == "o"
  list(
    "every level present" = factor(labels, levels = levels, ordered = ordered),
    "no row at the last level" = factor(
      ifelse(labels == levels[length(levels)], levels[1L], labels),
      levels = levels,
      ordered = ordered
    ),
    "levels declared in reverse" = factor(
      labels,
      levels = rev(levels),
      ordered = ordered
    ),
    "unused levels dropped" = droplevels(factor(
      ifelse(labels == levels[1L], levels[2L], labels),
      levels = levels,
      ordered = ordered
    )),
    "character" = labels,
    "sparseFactor" = dbarts::sparseFactor(factor(labels, levels = levels))
  )
}

# no warning from any accepted call
countWarn <- function(expr) {
  k <- 0L
  value <- withCallingHandlers(
    expr,
    warning = function(w) {
      k <<- k + 1L
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = k)
}

# ---- labels install at their own level, one sampler and two, f and o ------------------------
for (name in c("f", "o")) {
  levels <- lv[[name]]
  labels <- if (name == "f") newF else newO
  for (kind in names(kinds(labels, levels, name))) {
    given <- kinds(labels, levels, name)[[kind]]
    wantLabels <- as.character(given)
    one <- makeSampler(d, 2L)
    two <- list(
      makeSampler(d, 3L),
      makeSampler(d, 4L, form = y ~ o + f + x1)
    )
    info <- paste(name, kind)

    before <- list(heldLabels(one, name, levels))
    r <- countWarn(updatePredictorPerObservationJointly(list(one), given, name))
    expect_equal(r$warnings, 0L, info = info)
    expect_equal(length(r$value), n, info = info)
    expect_true(all(r$value), info = info)
    expectInstalled(list(one), name, before, wantLabels, levels, r$value, info)
    expect_identical(heldLabels(one, name, levels), wantLabels, info = info)
    one$run(0L, 2L)
    expect_identical(heldLabels(one, name, levels), wantLabels, info = info)

    before <- lapply(two, heldLabels, name, levels)
    r <- countWarn(updatePredictorPerObservationJointly(two, given, name))
    expect_equal(r$warnings, 0L, info = info)
    expectInstalled(two, name, before, wantLabels, levels, r$value, info)
    for (s in two) {
      s$run(0L, 2L)
    }
    expectInstalled(
      two,
      name,
      list(wantLabels, wantLabels),
      wantLabels,
      levels,
      r$value,
      info
    )
  }
}

# ---- identity with the column form, on a fixture where rows are declined --------------------
dSmall <- makeData(40L, 5L)
declined <- 0L
for (seed in 1:6) {
  a <- makeSampler(dSmall, seed, trees = 20L, burn = 50L)
  b <- makeSampler(dSmall, seed, trees = 20L, burn = 50L)
  given <- factor(
    sample(lv$f, 40L, TRUE),
    levels = rev(lv$f)
  )
  ma <- updatePredictorPerObservationJointly(list(a), given, "f")
  mb <- b$setPredictor(given, "f", forceUpdate = "partial")
  declined <- declined + sum(!ma)
  expect_identical(ma, mb)
  expect_identical(snapshot(a), snapshot(b))
  expect_identical(a$run(0L, 3L)$train, b$run(0L, 3L)$train)
}
expect_true(declined > 0L)

# ---- refusals leave every sampler as it was -------------------------------------------------
twinsOf <- function(seed) {
  list(
    makeSampler(d, seed),
    makeSampler(d, seed + 100L, form = y ~ o + f + x1)
  )
}
expectRefused <- function(given, name, message, info) {
  live <- twinsOf(7L)
  untouched <- twinsOf(7L)
  before <- lapply(live, snapshot)
  seedBefore <- get(".Random.seed", globalenv())
  expect_error(
    updatePredictorPerObservationJointly(live, given, name),
    message,
    fixed = TRUE,
    info = info
  )
  expect_identical(get(".Random.seed", globalenv()), seedBefore, info = info)
  for (i in 1:2) {
    expect_identical(snapshot(live[[i]]), before[[i]], info = info)
    expect_identical(snapshot(untouched[[i]]), before[[i]], info = info)
    expect_identical(
      live[[i]]$run(0L, 3L)$train,
      untouched[[i]]$run(0L, 3L)$train,
      info = info
    )
  }
}

missingLabel <- newF
missingLabel[1L] <- NA
expectRefused(
  c(newF[-1L], "z"),
  "f",
  "column 'f' has label 'z' not among its training levels",
  "unknown label, character"
)
expectRefused(
  factor(c(newF[-1L], "z")),
  "f",
  "column 'f' has label 'z' not among its training levels",
  "unknown label, factor"
)
expectRefused(
  factor(missingLabel, levels = lv$f),
  "f",
  "column 'f' has missing values, which its training values do not",
  "missing, factor"
)
expectRefused(
  missingLabel,
  "f",
  "column 'f' has missing values, which its training values do not",
  "missing, character"
)
expectRefused(
  rep(c(TRUE, FALSE), n / 2L),
  "f",
  "column 'f' is categorical; give its values as a factor or character vector of its labels, or as numbers for its codes from 0",
  "logical"
)
expectRefused(
  rep(c(TRUE, FALSE), n / 2L),
  "o",
  "column 'o' is categorical; give its values",
  "logical, ordered"
)

# ---- a column that holds a missing value takes labels with another --------------------------
dNA <- d
dNA$f[3L] <- NA
sNA <- makeSampler(dNA, 8L)
withMissing <- factor(newF, levels = lv$f)
withMissing[5L] <- NA
r <- updatePredictorPerObservationJointly(list(sNA), withMissing, "f")
expect_true(all(r))
expect_true(is.na(sNA$data@x[5L, "f"]))
expect_identical(
  heldLabels(sNA, "f", lv$f)[-5L],
  newF[-5L]
)

# ---- numbers are codes from 0, with or without a missing one --------------------------------
codes <- as.integer(factor(newF, levels = lv$f)) - 1L
for (given in list(codes, as.double(codes))) {
  s <- makeSampler(d, 9L)
  r <- updatePredictorPerObservationJointly(list(s), given, "f")
  expect_identical(s$data@x[, "f"][r], as.double(given)[r])
  expect_true(all(r))
}
s <- makeSampler(d, 9L)
withNA <- as.double(codes)
withNA[2L] <- NA
r <- updatePredictorPerObservationJointly(list(s), withNA, "f")
expect_true(is.na(s$data@x[2L, "f"]))
s <- makeSampler(d, 9L)
expect_error(
  updatePredictorPerObservationJointly(list(s), rep(4, n), "f"),
  "existing category codes"
)
expect_error(
  updatePredictorPerObservationJointly(list(s), rep(0.5, n), "f"),
  "existing category codes"
)
s <- makeSampler(d, 9L)
expect_error(
  updatePredictorPerObservationJointly(list(s), rep(3, n), "o"),
  "existing level codes"
)

# ---- numerals as labels ---------------------------------------------------------------------
dNum <- d
dNum$f <- factor(as.character(as.integer(d$f)), levels = as.character(1:4))
labelsNum <- as.character(sample(1:4, n, TRUE))
s <- makeSampler(dNum, 10L)
before <- list(heldLabels(s, "f", as.character(1:4)))
r <- updatePredictorPerObservationJointly(list(s), labelsNum, "f")
expect_identical(
  heldLabels(s, "f", as.character(1:4)),
  ifelse(r, labelsNum, before[[1L]])
)
s <- makeSampler(dNum, 10L)
r <- updatePredictorPerObservationJointly(list(s), rep(1L, n), "f")
expect_true(all(s$data@x[r, "f"] == 1))
expect_identical(
  heldLabels(s, "f", as.character(1:4))[r],
  rep("2", sum(r))
)

# ---- samplers that differ -------------------------------------------------------------------
dRev <- d
dRev$f <- factor(as.character(d$f), levels = rev(lv$f))
a <- makeSampler(d, 12L)
b <- makeSampler(dRev, 13L)
snapA <- snapshot(a)
snapB <- snapshot(b)
expect_error(
  updatePredictorPerObservationJointly(list(a, b), newF, "f"),
  "column 'f' has other levels in sampler 2 than in sampler 1, so its labels cannot be installed in both; give numbers",
  fixed = TRUE
)
expect_identical(snapshot(a), snapA)
expect_identical(snapshot(b), snapB)
r <- updatePredictorPerObservationJointly(list(a, b), codes, "f")
expect_true(all(r))
expect_identical(a$data@x[, "f"], b$data@x[, "f"])

dNumF <- d
dNumF$f <- as.numeric(d$f)
a <- makeSampler(d, 14L)
b <- makeSampler(dNumF, 15L)
snapA <- snapshot(a)
snapB <- snapshot(b)
expect_error(
  updatePredictorPerObservationJointly(list(a, b), newF, "f"),
  "column 'f' is categorical in sampler 1 and not in sampler 2, so its labels cannot be installed in both; give numbers",
  fixed = TRUE
)
expect_error(
  updatePredictorPerObservationJointly(list(b, a), newF, "f"),
  "column 'f' is categorical in sampler 2 and not in sampler 1",
  fixed = TRUE
)
expect_identical(snapshot(a), snapA)
expect_identical(snapshot(b), snapB)
r <- updatePredictorPerObservationJointly(list(a, b), codes, "f")
expect_true(all(r))

# ---- a numeric column takes numbers as before -----------------------------------------------
s <- makeSampler(d, 16L)
newX <- runif(n)
r <- updatePredictorPerObservationJointly(list(s), newX, "x1")
expect_identical(s$data@x[r, "x1"], newX[r])
