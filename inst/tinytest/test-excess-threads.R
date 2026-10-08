source(
  system.file("common", "friedmanData.R", package = "dbarts"),
  local = TRUE
)

# dec-B115: within-chain threading does not ship, so dbartsControl's own
# default caps the core probe at n.chains rather than handing out a budget
# tree sampling can never use past one thread per chain.
cores <- dbarts::guessNumCores()
expect_equal(dbarts::dbartsControl(n.chains = 1L)@n.threads, min(cores, 1L))
expect_equal(dbarts::dbartsControl(n.chains = 2L)@n.threads, min(cores, 2L))
# a bare new() keeps the prototype's conservative 1L, untouched by the
# constructor's capping arithmetic
expect_equal(methods::new("dbartsControl")@n.threads, 1L)

# a caller-supplied budget above n.chains still warns, once per fit, naming
# both numbers - dbarts() and bart() share the one site (bart() forwards a
# built control into dbarts()), so pinning dbarts() covers both doors
excessControl <- dbarts::dbartsControl(
  n.chains = 4L,
  n.threads = 8L,
  n.trees = 3L,
  n.samples = 2L,
  n.burn = 1L,
  verbose = FALSE,
  updateState = FALSE
)
expect_warning(
  dbarts::dbarts(y ~ x, testData, control = excessControl),
  pattern = "n.threads (8) exceeds n.chains (4)",
  fixed = TRUE
)

# at or below n.chains - the shape the default itself produces - it is
# silent
evenControl <- dbarts::dbartsControl(
  n.chains = 4L,
  n.threads = 4L,
  n.trees = 3L,
  n.samples = 2L,
  n.burn = 1L,
  verbose = FALSE,
  updateState = FALSE
)
expect_silent(dbarts::dbarts(y ~ x, testData, control = evenControl))

# bartBT's own arguments are the 0.9-x names, so its warning names what its
# caller typed
expect_warning(
  dbarts::bartBT(
    testData$x,
    testData$y,
    ntree = 3L,
    ndpost = 2L,
    nskip = 1L,
    nthread = 2L,
    verbose = FALSE
  ),
  pattern = "nthread (2) exceeds nchain (1)",
  fixed = TRUE
)

# the re-issued class vector is exactly dbarts()'s own, with no doubled tail
# a hurdle fit is two samplers but one fit: the warning is raised once
hurdleY <- pmax(testData$y - 15, 0)
nExcess <- 0L
withCallingHandlers(
  dbarts::bart(
    testData$x,
    hurdleY,
    family = "hurdle.lognormal",
    n.trees = 3L,
    n.samples = 2L,
    n.burn = 1L,
    n.chains = 1L,
    n.threads = 3L,
    verbose = FALSE
  ),
  warning = function(w) {
    if (grepl("exceeds n.chains", conditionMessage(w), fixed = TRUE)) {
      nExcess <<- nExcess + 1L
    }
    invokeRestart("muffleWarning")
  }
)
expect_equal(nExcess, 1L)

# dec-B306: where the cores cannot be counted, guessNumCores() is NA and an
# n.threads that is not stated is one, with no message, at every door that
# defaults it; a stated NA is still refused. guessNumCores is replaced in the
# namespace for the length of the block, in process.
ns <- asNamespace("dbarts")
realGuess <- get("guessNumCores", ns)
replaceGuess <- function(value) {
  unlockBinding("guessNumCores", ns)
  assign("guessNumCores", value, ns)
  lockBinding("guessNumCores", ns)
}
warnState <- dbarts:::onceWarnState
onceKeys <- c(dbarts:::frontDoorDefaultsKey, "tombstone.rbart_vi")
onceOnEntry <- mget(onceKeys, warnState, ifnotfound = list(NULL))
for (key in onceKeys) {
  warnState[[key]] <- TRUE
}
groupData <- data.frame(
  y = testData$y,
  x = testData$x[, 1L],
  g = rep_len(1:3, length(testData$y))
)
replaceGuess(function(logical = FALSE) NA_integer_)
noCores <- tryCatch(
  {
    expect_identical(dbarts::guessNumCores(), NA_integer_)
    expect_identical(dbarts::dbartsControl()@n.threads, 1L)
    expect_identical(dbarts::dbartsControl(n.chains = 3L)@n.threads, 1L)
    expect_silent(fit <- dbarts::bart(
      y ~ x,
      groupData,
      n.trees = 3L,
      n.samples = 2L,
      n.burn = 1L,
      n.chains = 2L,
      verbose = FALSE
    ))
    expect_silent(dbarts::bartBT(
      testData$x,
      testData$y,
      ntree = 3L,
      ndpost = 2L,
      nskip = 1L,
      verbose = FALSE
    ))
    expect_silent(dbarts::rbart_vi(
      y ~ x,
      groupData,
      group.by = g,
      n.trees = 3L,
      n.samples = 2L,
      n.burn = 1L,
      n.chains = 2L,
      n.thin = 1L,
      verbose = FALSE
    ))
    expect_silent(dbarts::xbart(
      y ~ x,
      groupData,
      n.trees = 3L,
      n.reps = 1L,
      n.samples = 2L,
      n.burn = c(2L, 1L),
      verbose = FALSE
    ))
    # a stated NA is refused, at the control and at the doors that read it
    message <- paste0(
      "'n.threads' must be a positive integer, not NA; guessNumCores() ",
      "returns NA when it cannot count this system's cores, and 'n.threads' ",
      "is then one unless it is given a count"
    )
    expect_error(
      dbarts::dbartsControl(n.threads = NA_integer_),
      message,
      fixed = TRUE
    )
    expect_error(
      dbarts::bart(y ~ x, groupData, n.threads = NA_integer_),
      message,
      fixed = TRUE
    )
    expect_error(
      dbarts::xbart(y ~ x, groupData, n.threads = NA_integer_),
      message,
      fixed = TRUE
    )
    expect_error(
      dbarts::rbart_vi(y ~ x, groupData, group.by = g, n.threads = NA_integer_),
      message,
      fixed = TRUE
    )
    TRUE
  },
  finally = {
    replaceGuess(realGuess)
    for (key in onceKeys) {
      warnState[[key]] <- onceOnEntry[[key]]
    }
  }
)
expect_identical(get("guessNumCores", ns), realGuess)
rm(ns, realGuess, replaceGuess, groupData, noCores, warnState, onceKeys, onceOnEntry, key)

rm(cores, excessControl, evenControl, testData, hurdleY, nExcess)
