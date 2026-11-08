# rbart_vi says each creation warning once per fit, whatever its chain and
# thread counts: its checking sampler says it, and the samplers of its chains
# say nothing the checking sampler said. One row of positive weight leaves the
# likelihood without two distinct values of the response, which a sampler's
# creation says, and its warm-up repeats where the rows are set again.

# the port's deprecation warning is pinned in test-rbart-port.R
onceState <- dbarts:::onceWarnState
onceState[["tombstone.rbart_vi"]] <- TRUE
rm(onceState)

set.seed(20261108)
n <- 60L
d <- data.frame(y = rnorm(n), x = rnorm(n), g = factor(rep(1:3, n / 3L)))
d$w <- c(1, rep(0, n - 1L))
# the grouping, read from the calling frame as rbart_vi reads it
g <- d$g

said <- function(n.chains, n.threads) {
  texts <- character(0L)
  withCallingHandlers(
    suppressMessages(dbarts::rbart_vi(
      y ~ x,
      d,
      group.by = g,
      weights = w,
      n.samples = 4L,
      n.burn = 4L,
      n.thin = 1L,
      n.trees = 5L,
      n.chains = n.chains,
      n.threads = n.threads,
      verbose = FALSE
    )),
    warning = function(w) {
      texts <<- c(texts, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  texts
}

once <- said(1L, 1L)
expect_identical(length(once), 1L)
expect_true(grepl("taken from all 60", once, fixed = TRUE))
expect_false(any(grepl("said once per sampler", once, fixed = TRUE)))
for (chains in c(3L, 4L)) {
  expect_identical(said(chains, 1L), once, info = paste(chains, "chains"))
}
for (threads in 1:2) {
  expect_identical(said(3L, threads), once, info = paste(threads, "threads"))
}

# a creation warning of another kind, the thread budget's, is said once too,
# beside the scale fallback
# the budget's warning is said once a session, so it is cleared first
onceState <- dbarts:::onceWarnState
onceState[[dbarts:::excessThreadsKey]] <- NULL
both <- said(1L, 2L)
expect_identical(length(both), 2L)
expect_identical(sum(grepl("exceeds n.chains", both, fixed = TRUE)), 1L)
expect_identical(sum(grepl("taken from all 60", both, fixed = TRUE)), 1L)
