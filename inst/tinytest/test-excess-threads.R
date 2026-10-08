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
# what an expression said: messages and warnings counted, printed lines kept
said <- function(expr) {
  value <- NULL
  nMessages <- 0L
  nWarnings <- 0L
  lines <- capture.output(
    withCallingHandlers(
      value <- expr,
      message = function(m) {
        nMessages <<- nMessages + 1L
        invokeRestart("muffleMessage")
      },
      warning = function(w) {
        nWarnings <<- nWarnings + 1L
        invokeRestart("muffleWarning")
      }
    )
  )
  list(value = value, said = c(nMessages, nWarnings, length(lines)))
}
nothing <- c(0L, 0L, 0L)
replaceGuess(function(logical = FALSE) NA_integer_)
noCores <- tryCatch(
  {
    expect_identical(dbarts::guessNumCores(), NA_integer_)
    out <- said(dbarts::dbartsControl())
    expect_identical(out$value@n.threads, 1L)
    expect_identical(out$said, nothing)
    out <- said(dbarts::dbartsControl(n.chains = 3L))
    expect_identical(out$value@n.threads, 1L)
    expect_identical(out$said, nothing)
    out <- said(dbarts::bart(
      y ~ x,
      groupData,
      n.trees = 3L,
      n.samples = 2L,
      n.burn = 1L,
      n.chains = 2L,
      keepSampler = TRUE,
      verbose = FALSE
    ))
    expect_identical(out$value$fit$control@n.threads, 1L)
    expect_identical(out$said, nothing)
    out <- said(dbarts::rbart_vi(
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
    expect_identical(out$said, nothing)
    out <- said(dbarts::xbart(
      y ~ x,
      groupData,
      n.trees = 3L,
      n.reps = 1L,
      n.samples = 2L,
      n.burn = c(2L, 1L),
      verbose = FALSE
    ))
    expect_identical(out$said, nothing)
    # a wrapper that hands on its dots leaves the count unstated; one that
    # states a count of its own, NA included, is refused
    wrapControl <- function(...) dbarts::dbartsControl(...)
    wrapBart <- function(...) {
      dbarts::bart(
        y ~ x,
        groupData,
        n.trees = 3L,
        n.samples = 2L,
        n.burn = 1L,
        n.chains = 1L,
        keepSampler = TRUE,
        verbose = FALSE,
        ...
      )
    }
    expect_identical(wrapControl()@n.threads, 1L)
    expect_identical(wrapBart()$fit$control@n.threads, 1L)
    expect_error(wrapControl(n.threads = NA_integer_), "not NA")
    expect_error(wrapBart(n.threads = NA_integer_), "not NA")
    # a control given and no count stated: the control's own count, or the
    # default when its slot is the fresh control's, is one and nothing is said
    controlFit <- function(control, ...) {
      said(dbarts::bart(
        y ~ x,
        groupData,
        control = control,
        n.samples = 2L,
        keepSampler = TRUE,
        verbose = FALSE,
        ...
      ))
    }
    heldControl <- dbarts::dbartsControl()
    shapes <- list(
      controlFit(dbarts::dbartsControl()),
      controlFit(dbarts::dbartsControl(n.trees = 3L, n.burn = 1L)),
      controlFit(dbarts::dbartsControl(n.chains = 2L, n.burn = 1L)),
      controlFit(dbarts::dbartsControl(n.burn = 1L), n.chains = 2L),
      controlFit(methods::new("dbartsControl"), n.burn = 1L),
      controlFit(heldControl, n.burn = 1L, n.trees = 3L)
    )
    for (shape in shapes) {
      expect_identical(shape$value$fit$control@n.threads, 1L)
      expect_identical(shape$said, nothing)
    }
    multinomial <- said(dbarts::bart(
      factor(rep_len(c("a", "b", "c"), nrow(groupData))) ~ x,
      groupData,
      control = dbarts::dbartsControl(n.burn = 1L),
      n.trees = 3L,
      n.samples = 2L,
      family = "multinomial",
      verbose = FALSE
    ))
    expect_identical(multinomial$said, nothing)
    # a stated NA is refused, at the control and at the doors that read it,
    # in words that tell its writer what to do
    message <- paste0(
      "'n.threads' must be a positive integer, not NA; leave it out to take ",
      "the default"
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
    # a control whose count was set to NA after it was built has a slot to
    # name, not an argument to leave out
    badControl <- dbarts::dbartsControl(n.trees = 3L, n.chains = 1L)
    badControl@n.threads <- NA_integer_
    slotMessage <- "a control's 'n.threads' slot must be a positive integer"
    expect_error(
      dbarts::dbarts(y ~ x, groupData, control = badControl),
      slotMessage,
      fixed = TRUE
    )
    expect_error(
      dbarts::bart(y ~ x, groupData, control = badControl),
      slotMessage,
      fixed = TRUE
    )
    expect_error(
      dbarts::bartBT(testData$x, testData$y, nthread = NA_integer_),
      "'nthread' must be a positive integer, not NA",
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
rm(
  ns,
  realGuess,
  replaceGuess,
  said,
  nothing,
  groupData,
  noCores,
  warnState,
  onceKeys,
  onceOnEntry,
  key
)

rm(cores, excessControl, evenControl, testData, hurdleY, nExcess)
