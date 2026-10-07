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
  d$y <- d$x1 + 2 * as.integer(d$f) + 3 * as.integer(d$o) + rnorm(n, sd = 0.3)
  d
}

# a sampler on formula, swept so its trees split on the factor columns (plain() has no split); `form` fixes the column order, so two samplers can hold f and o at
# different positions
makeSampler <- function(
  d,
  seed,
  trees = 20L,
  form = y ~ x1 + f + o,
  burn = 40L
) {
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

plain <- function(d, seed, form = y ~ x1 + f + o) {
  makeSampler(d, seed, trees = 5L, form = form, burn = 0L)
}

# the position in the design of a column, and whether some tree splits on it
splitsOn <- function(s, name) {
  any(s$getTrees()$var == match(name, colnames(s$data@x)))
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

# ---- labels install at their own level, one sampler and several, f and o -----------------------
# Each case is read on the engine: a twin given the same labels through setPredictor "partial"
# (one sampler) or the codes from 0 (several, which a shared scan cannot take separately) must
# return the same mask and hold the same stored state and design, and its next draws must be
# the same, which they cannot be if the engine was handed other codes than the design holds.
forms <- list(y ~ x1 + f + o, y ~ o + f + x1, y ~ f + x1 + o)
makeSet <- function(k) {
  lapply(seq_len(k), function(i) makeSampler(d, 1L + i, form = forms[[i]]))
}
anyDeclined <- FALSE
for (name in c("f", "o")) {
  levels <- lv[[name]]
  labels <- if (name == "f") newF else newO
  for (kind in names(kinds(labels, levels, name))) {
    for (k in 1:3) {
      given <- kinds(labels, levels, name)[[kind]]
      wantLabels <- as.character(given)
      codes <- match(wantLabels, levels) - 1
      live <- makeSet(k)
      twin <- makeSet(k)
      info <- paste(name, kind, k, "samplers")
      expect_true(splitsOn(live[[1L]], name), info = info)

      before <- lapply(live, heldLabels, name, levels)
      r <- countWarn(updatePredictorPerObservationJointly(live, given, name))
      expect_equal(r$warnings, 0L, info = info)
      mask <- r$value
      expect_equal(length(mask), n, info = info)
      anyDeclined <- anyDeclined || !all(mask)
      expectInstalled(live, name, before, wantLabels, levels, mask, info)

      maskTwin <- if (k == 1L) {
        twin[[1L]]$setPredictor(given, name, forceUpdate = "partial")
      } else {
        updatePredictorPerObservationJointly(twin, codes, name)
      }
      expect_identical(mask, maskTwin, info = info)
      for (i in seq_len(k)) {
        expect_identical(snapshot(live[[i]]), snapshot(twin[[i]]), info = info)
        expect_identical(
          live[[i]]$run(0L, 2L)$train,
          twin[[i]]$run(0L, 2L)$train,
          info = paste(info, "next draws, sampler", i)
        )
      }
    }
  }
}
expect_true(anyDeclined)

# ---- a copy and a reload of an updated sampler continue as the twin's does -------------------
for (name in c("f", "o")) {
  levels <- lv[[name]]
  given <- factor(
    if (name == "f") newF else newO,
    levels = rev(levels)
  )
  live <- makeSampler(d, 2L)
  twin <- makeSampler(d, 2L)
  updatePredictorPerObservationJointly(list(live), given, name)
  twin$setPredictor(given, name, forceUpdate = "partial")
  copyLive <- live$copy()
  copyTwin <- twin$copy()
  expect_identical(copyLive$data@x[, name], twin$data@x[, name], info = name)
  expect_identical(
    copyLive$run(0L, 2L)$train,
    copyTwin$run(0L, 2L)$train,
    info = paste(name, "copy")
  )
  live$storeState()
  twin$storeState()
  file <- tempfile(fileext = ".rds")
  saveRDS(live, file)
  reloadLive <- readRDS(file)
  saveRDS(twin, file)
  reloadTwin <- readRDS(file)
  unlink(file)
  expect_identical(reloadLive$data@x[, name], twin$data@x[, name], info = name)
  expect_identical(
    reloadLive$run(0L, 2L)$train,
    reloadTwin$run(0L, 2L)$train,
    info = paste(name, "reload")
  )
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
expectRefused <- function(given, name, message, info, make = twinsOf) {
  live <- make(7L)
  untouched <- make(7L)
  before <- lapply(live, snapshot)
  expect_error(
    updatePredictorPerObservationJointly(live, given, name),
    message,
    fixed = TRUE,
    info = info
  )
  for (i in seq_along(live)) {
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

# ---- a refusal a later sampler raises names it, and touches nothing -------------------------
# samplers 1 and 2 hold a missing value in f, sampler 3 holds none
dHole <- d
dHole$f[3L] <- NA
threeWithHole <- function(seed) {
  list(
    makeSampler(dHole, seed),
    makeSampler(dHole, seed + 1L, form = y ~ o + f + x1),
    makeSampler(d, seed + 2L, form = y ~ f + x1 + o)
  )
}
expectRefused(
  factor(missingLabel, levels = lv$f),
  "f",
  "column 'f' has missing values, which its training values do not (sampler 3)",
  "missing value refused by the third sampler",
  make = threeWithHole
)
# the same labels with a missing value are taken when every sampler holds one
allHole <- list(
  plain(dHole, 8L),
  plain(dHole, 9L, form = y ~ o + f + x1),
  plain(dHole, 10L, form = y ~ f + x1 + o)
)
r <- updatePredictorPerObservationJointly(
  allHole,
  factor(missingLabel, levels = lv$f),
  "f"
)
expect_true(all(r))
for (s in allHole) {
  expect_true(is.na(s$data@x[1L, "f"]))
  expect_identical(heldLabels(s, "f", lv$f)[-1L], newF[-1L])
}

# ---- a column that holds a missing value takes labels with another --------------------------
sNA <- plain(dHole, 8L)
withMissing <- factor(newF, levels = lv$f)
withMissing[5L] <- NA
r <- updatePredictorPerObservationJointly(list(sNA), withMissing, "f")
expect_true(all(r))
expect_true(is.na(sNA$data@x[5L, "f"]))
expect_identical(heldLabels(sNA, "f", lv$f)[-5L], newF[-5L])

# ---- numbers are codes from 0, with or without a missing one --------------------------------
codes <- as.integer(factor(newF, levels = lv$f)) - 1L
for (given in list(codes, as.double(codes))) {
  s <- plain(d, 9L)
  r <- updatePredictorPerObservationJointly(list(s), given, "f")
  expect_true(all(r))
  expect_identical(s$data@x[, "f"], as.double(given))
}
s <- plain(d, 9L)
withNA <- as.double(codes)
withNA[2L] <- NA
r <- updatePredictorPerObservationJointly(list(s), withNA, "f")
expect_true(is.na(s$data@x[2L, "f"]))
# a missing number installs, where a missing label in a column with none is refused
expect_true(r[2L])
s <- plain(d, 9L)
expect_error(
  updatePredictorPerObservationJointly(list(s), rep(4, n), "f"),
  "existing category codes"
)
expect_error(
  updatePredictorPerObservationJointly(list(s), rep(0.5, n), "f"),
  "existing category codes"
)
s <- plain(d, 9L)
expect_error(
  updatePredictorPerObservationJointly(list(s), rep(3, n), "o"),
  "existing level codes"
)

# ---- numerals as labels ---------------------------------------------------------------------
dNum <- d
dNum$f <- factor(as.character(as.integer(d$f)), levels = as.character(1:4))
labelsNum <- as.character(sample(1:4, n, TRUE))
s <- plain(dNum, 10L)
r <- updatePredictorPerObservationJointly(list(s), labelsNum, "f")
expect_identical(heldLabels(s, "f", as.character(1:4)), labelsNum)
s <- plain(dNum, 10L)
r <- updatePredictorPerObservationJointly(list(s), rep(1L, n), "f")
expect_true(all(s$data@x[, "f"] == 1))
expect_identical(heldLabels(s, "f", as.character(1:4)), rep("2", n))

# ---- samplers that hold the column with different levels -----------------------------------
# one vector would be a different level in each, so neither labels nor numbers are taken
dRev <- d
dRev$f <- factor(as.character(d$f), levels = rev(lv$f))
dSub <- d
dSub$f <- droplevels(factor(ifelse(d$f == "d", "a", as.character(d$f))))
differing <- function(second, third = d) {
  list(
    plain(d, 12L),
    plain(second, 13L, form = y ~ o + f + x1),
    plain(third, 14L, form = y ~ f + x1 + o)
  )
}
refusedNothingTouched <- function(samplers, given, name, message, info) {
  snaps <- lapply(samplers, snapshot)
  expect_error(
    updatePredictorPerObservationJointly(samplers, given, name),
    message,
    fixed = TRUE,
    info = info
  )
  for (i in seq_along(samplers)) {
    expect_identical(snapshot(samplers[[i]]), snaps[[i]], info = info)
  }
}
tail2 <- paste0(
  "), so one value would be a different level in each; update them in ",
  "separate calls, or create them with the same levels in the same order"
)
codes <- as.integer(factor(newF, levels = lv$f)) - 1L
for (given in list(newF, factor(newF), codes, as.double(codes), rep(TRUE, n))) {
  refusedNothingTouched(
    differing(dRev),
    given,
    "f",
    paste0(
      "the samplers hold column 'f' with different levels (sampler 2 differs from sampler 1",
      tail2
    ),
    paste("reversed second,", class(given)[1L])
  )
  refusedNothingTouched(
    differing(d, dRev),
    given,
    "f",
    paste0(
      "the samplers hold column 'f' with different levels (sampler 3 differs from sampler 1",
      tail2
    ),
    paste("reversed third,", class(given)[1L])
  )
}
refusedNothingTouched(
  differing(dSub),
  newF,
  "f",
  "(sampler 2 differs from sampler 1",
  "a subset of the levels"
)
# the same table in every sampler takes labels and numbers
same <- differing(d, d)
r <- updatePredictorPerObservationJointly(same, codes, "f")
expect_true(all(r))
expect_identical(same[[1L]]$data@x[, "f"], same[[3L]]$data@x[, "f"])

# ---- a column held as a factor in one sampler and as a number in another --------------------
dNumF <- d
dNumF$f <- as.numeric(d$f)
a <- plain(d, 14L)
b <- plain(dNumF, 15L)
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

# ---- a one-column character matrix is refused as a matrix -----------------------------------
expectRefused(
  matrix(newF, ncol = 1L),
  "f",
  "column 'f' is categorical; give its values as a factor or character vector of its labels, or as numbers for its codes from 0",
  "character matrix"
)

# ---- a numeric column takes numbers, and refuses labels -------------------------------------
s <- plain(d, 16L)
newX <- runif(n)
r <- updatePredictorPerObservationJointly(list(s), newX, "x1")
expect_identical(s$data@x[r, "x1"], newX[r])
r <- updatePredictorPerObservationJointly(list(s), as.character(newX), "x1")
expect_equal(s$data@x[r, "x1"], newX[r])
for (given in list(
  factor(newF),
  newF,
  dbarts::sparseFactor(factor(newF))
)) {
  s <- plain(d, 16L)
  snap <- snapshot(s)
  expect_error(
    updatePredictorPerObservationJointly(list(s), given, "x1"),
    "column 'x1' is numeric and cannot take labels",
    fixed = TRUE,
    info = class(given)[1L]
  )
  expect_identical(snapshot(s), snap)
}
