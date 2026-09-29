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
  class = "dbartsExcessThreadsWarning",
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
  class = "dbartsExcessThreadsWarning",
  fixed = TRUE
)

# the re-issued class vector is exactly dbarts()'s own, with no doubled tail
seenClass <- NULL
withCallingHandlers(
  dbarts::bartBT(
    testData$x,
    testData$y,
    ntree = 3L,
    ndpost = 2L,
    nskip = 1L,
    nthread = 2L,
    verbose = FALSE
  ),
  warning = function(w) {
    seenClass <<- class(w)
    invokeRestart("muffleWarning")
  }
)
expect_identical(
  seenClass,
  c("dbartsExcessThreadsWarning", "dbartsWarning", "warning", "condition")
)

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
  dbartsExcessThreadsWarning = function(w) {
    nExcess <<- nExcess + 1L
    invokeRestart("muffleWarning")
  }
)
expect_equal(nExcess, 1L)

rm(cores, excessControl, evenControl, testData, seenClass, hurdleY, nExcess)
