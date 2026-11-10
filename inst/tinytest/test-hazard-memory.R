# A hazard fit refuses no row count: before its person-period rows are built
# it projects its peak memory, and when that exceeds the machine's physical
# memory it warns once per session (dbartsMemoryWarning) and runs. The
# unexported option dbarts.physicalMemory stands in for the machine's figure,
# NA for one that cannot be read.

memoryKey <- "hazardMemory"
onceState <- dbarts:::onceWarnState
resetKey <- function() onceState[[memoryKey]] <- NULL
oldOptions <- options(dbarts.physicalMemory = NULL)

# ---- the projection's arithmetic ----
hazardFitBytes <- dbarts:::hazardFitBytes
hazardRunSettings <- dbarts:::hazardRunSettings
defaults <- hazardRunSettings(
  dbartsControl(n.chains = 4L, n.threads = 1L, n.samples = 500L),
  75L
)
# bart()'s defaults: 200 + 34 * 3 shared, 20 + 8 * 75 a chain, and 16 bytes a
# kept draw at the peak, 4 x 500 of them a training row
expect_equal(hazardFitBytes(1, 0, 3L, defaults), 34782)
expect_equal(
  hazardFitBytes(1e7, 0, 3L, defaults),
  1e7 * (200 + 34 * 3 + 4 * (20 + 8 * 75) + 16 * 4 * 500)
)
# a held-out row: 30 + 36 a predictor column and every kept draw
expect_equal(hazardFitBytes(0, 1, 3L, defaults), 30 + 36 * 3 + 16 * 2000)
expect_equal(
  hazardFitBytes(2, 5, 10L, defaults),
  2 *
    hazardFitBytes(1, 0, 10L, defaults) +
    5 * hazardFitBytes(0, 1, 10L, defaults)
)
# keepTrainingFits = FALSE drops the training rows' kept draws only
noTraining <- hazardRunSettings(
  dbartsControl(
    n.chains = 4L,
    n.threads = 1L,
    n.samples = 500L,
    keepTrainingFits = FALSE
  ),
  75L
)
expect_equal(hazardFitBytes(1, 0, 3L, noTraining), 2782)
expect_equal(hazardFitBytes(0, 1, 3L, noTraining), 30 + 36 * 3 + 16 * 2000)
# keepFits = FALSE drops every row's
noFits <- hazardRunSettings(
  dbartsControl(
    n.chains = 4L,
    n.threads = 1L,
    n.samples = 500L,
    keepFits = FALSE
  ),
  75L
)
expect_equal(hazardFitBytes(1, 1, 3L, noFits), 2782 + 30 + 36 * 3)
# a row count past the integer range is projected in double precision
expect_equal(hazardFitBytes(3e9, 0, 3L, defaults), 3e9 * 34782)

# ---- the machine's memory, read in C ----
options(dbarts.physicalMemory = NULL)
machine <- dbarts:::physicalMemory()
expect_true(is.double(machine) && length(machine) == 1L)
if (Sys.info()[["sysname"]] %in% c("Darwin", "Linux", "Windows")) {
  expect_true(is.finite(machine) && machine > 0 && machine == round(machine))
}

# ---- the warning, once per session, and the fit runs ----
set.seed(3L)
numSubjects <- 60L
x <- matrix(rnorm(numSubjects * 2L), numSubjects, 2L)
time <- rexp(numSubjects)
status <- rbinom(numSubjects, 1L, 0.7)
xTest <- matrix(rnorm(10L * 2L), 10L, 2L)
fitHazard <- function(...) {
  dbarts::bart(
    x,
    cbind(time, status),
    family = "hazard",
    n.trees = 5L,
    n.samples = 10L,
    n.burn = 5L,
    n.chains = 2L,
    n.threads = 1L,
    verbose = FALSE,
    ...
  )
}
collect <- function(expr) {
  memory <- list()
  value <- withCallingHandlers(expr, dbartsMemoryWarning = function(w) {
    memory[[length(memory) + 1L]] <<- w
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = memory)
}

# a quarter of a megabyte, which every fit below exceeds
options(dbarts.physicalMemory = 2^18)
resetKey()
first <- collect(fitHazard(test = xTest))
expect_equal(length(first$warnings), 1L)
expect_inherits(first$warnings[[1L]], "dbartsWarning")
# no refusal: the fit runs whatever the projection says
expect_inherits(first$value, "bart")
rows <- length(first$value$y)
testRows <- 10L * length(first$value$periods)
# distinct times: subject i of the sorted times is at risk i periods
expect_equal(rows, numSubjects * (numSubjects + 1L) / 2L)
projected <- hazardFitBytes(
  rows,
  testRows,
  2L,
  hazardRunSettings(
    dbartsControl(n.chains = 2L, n.threads = 1L, n.samples = 10L),
    5L
  )
)
expect_equal(
  conditionMessage(first$warnings[[1L]]),
  paste0(
    "the hazard fit expands to ",
    rows,
    " person-period rows and ",
    testRows,
    " test rows; at 2 chains keeping 10 draws a chain of the training fits ",
    "it needs about ",
    dbarts:::formatBytes(projected),
    " at its peak, and this machine has 0.25 MB; coarsen the grid with ",
    "hazard(breaks = ), set keepTrainingFits = FALSE, or run fewer chains. ",
    "This warning does not stop the fit and is shown once per session"
  )
)
# once per session
expect_equal(length(collect(fitHazard(test = xTest))$warnings), 0L)
# the trigger is the projection itself: a machine of exactly that memory is
# not warned, and one a byte smaller is
options(dbarts.physicalMemory = projected)
resetKey()
expect_equal(length(collect(fitHazard(test = xTest))$warnings), 0L)
options(dbarts.physicalMemory = projected - 1)
expect_equal(length(collect(fitHazard(test = xTest))$warnings), 1L)
# and the fit is the one an unwarned session makes
options(dbarts.physicalMemory = NA)
unwarned <- fitHazard(test = xTest, seed = 2L)
options(dbarts.physicalMemory = 2^18)
resetKey()
warned <- collect(fitHazard(test = xTest, seed = 2L))
expect_equal(length(warned$warnings), 1L)
expect_identical(warned$value$yhat.train, unwarned$yhat.train)
expect_identical(warned$value$yhat.test, unwarned$yhat.test)

# the formula door warns from the same projection, and names only the levers
# the run has: no kept training fits, one chain
frame <- data.frame(x1 = x[, 1L], x2 = x[, 2L], time = time, status = status)
resetKey()
formulaFit <- collect(dbarts::dbarts(
  survival::Surv(time, status) ~ x1 + x2,
  frame,
  frame[1:4, ],
  family = "hazard",
  n.trees = 5L,
  control = dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.samples = 10L,
    keepTrainingFits = FALSE
  )
))
formulaWarnings <- formulaFit$warnings
expect_equal(length(formulaWarnings), 1L)
expect_true(grepl(
  paste0(
    " and ",
    4L * length(formulaFit$value$model@hazard.periods),
    " test rows; at 1 chain keeping 10 draws a chain of the test fits it ",
    "needs about .*; coarsen the grid with hazard\\(breaks = \\)\\. This"
  ),
  conditionMessage(formulaWarnings[[1L]])
))

# no warning where the machine's memory cannot be read, nor below it
for (memory in list(NA_real_, 1e15)) {
  options(dbarts.physicalMemory = memory)
  resetKey()
  expect_equal(length(collect(fitHazard())$warnings), 0L)
  expect_null(onceState[[memoryKey]])
}

# ---- $setData warns from the sampler's own projection ----
form <- survival::Surv(time, status) ~ x1 + x2
newFrame <- frame[1:30, ]
newTest <- frame[1:3, c("x1", "x2")]
# one warning a session whichever door gives it: a sampler created over its
# memory warns, and its $setData then does not
options(dbarts.physicalMemory = 2^18)
resetKey()
created <- collect(dbarts::dbarts(
  form,
  frame,
  family = "hazard",
  n.trees = 5L,
  control = dbartsControl(n.chains = 2L, n.threads = 1L, n.samples = 10L)
))
expect_equal(length(created$warnings), 1L)
sampler <- created$value
expect_equal(
  length(collect(sampler$copy()$setData(dbartsData(form, newFrame)))$warnings),
  0L
)
# with the key unspent, $setData warns, counting its test subjects at every
# period, and names only what an existing sampler can do
resetKey()
swapped <- sampler$copy()
swapWarnings <- collect(swapped$setData(dbartsData(
  form,
  newFrame,
  newTest
)))$warnings
expect_equal(length(swapWarnings), 1L)
swapRows <- nrow(swapped$data@x)
swapTestRows <- 3L * length(sampler$model@hazard.periods)
expect_equal(nrow(swapped$data@x.test), swapTestRows)
expect_equal(
  conditionMessage(swapWarnings[[1L]]),
  paste0(
    "$setData: the hazard fit expands to ",
    swapRows,
    " person-period rows and ",
    swapTestRows,
    " test rows on the sampler's periods; at 2 chains keeping 10 ",
    "draws a chain of the training fits it needs about ",
    dbarts:::formatBytes(hazardFitBytes(
      swapRows,
      swapTestRows,
      2L,
      hazardRunSettings(sampler$control, 5L)
    )),
    " at its peak, and this machine has 0.25 MB; hand over fewer ",
    "subjects, set keepTrainingFits = FALSE, or make a new sampler with ",
    "hazard(breaks = ) or fewer chains. This warning does not stop the fit ",
    "and is shown once per session"
  )
)
# what it names can be done there: keepTrainingFits is taken, the chain count
# is not
expect_silent(swapped$copy()$setControl(dbartsControl(
  n.chains = 2L,
  n.threads = 1L,
  n.samples = 10L,
  keepTrainingFits = FALSE
)))
expect_error(
  swapped$copy()$setControl(dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.samples = 10L
  )),
  "n.chains"
)
# the data are taken
expect_equal(
  swapRows,
  sum(
    findInterval(
      newFrame$time,
      sampler$model@hazard.periods,
      left.open = TRUE
    ) +
      1L
  )
)

# ---- max.rows is gone, refused as any unknown hazard() argument is ----
unknown <- function(expr) tryCatch(expr, error = conditionMessage)
expect_true(startsWith(
  unknown(dbartsFamilies$hazard(max.rows = 1e7)),
  "unused argument"
))
expect_true(startsWith(
  unknown(dbartsFamilies$hazard(cap = 1e7)),
  "unused argument"
))
# inside a fitting function's 'family', where hazard() resolves by bare name
expect_identical(
  unknown(dbarts::bart(
    x,
    cbind(time, status),
    family = hazard(max.rows = 10),
    verbose = FALSE
  )),
  "unused argument (max.rows = 10)"
)
expect_identical(
  unknown(dbarts::bart(
    x,
    cbind(time, status),
    family = hazard(cap = 10),
    verbose = FALSE
  )),
  "unused argument (cap = 10)"
)
expect_identical(names(dbartsFamilies$hazard()@settings), "breaks")

options(oldOptions)
resetKey()
