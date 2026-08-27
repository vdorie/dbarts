# family = "auto" response-type detection and categorical-response routing:
# the formula path no longer hard-errors on a factor response, auto detects the
# response type at every entry point, logical and character responses coerce,
# and a family that cannot fit a categorical response meets an informative
# error instead of silently fitting the integer level codes.

set.seed(101)
n <- 120L
x <- matrix(runif(n * 3L), n, 3L)
colnames(x) <- c("x1", "x2", "x3")
lin <- 2 * (x[, 1L] - 0.5) + (x[, 2L] - x[, 3L])
yb <- rbinom(n, 1L, plogis(lin))
yf2 <- factor(ifelse(yb == 1L, "yes", "no"), levels = c("no", "yes"))
ylog <- yb == 1L
ychar <- as.character(yf2)
y3 <- factor(sample(c("a", "b", "c"), n, replace = TRUE))
df <- data.frame(
  x,
  yf2 = yf2,
  y3 = y3,
  ychar = ychar,
  stringsAsFactors = FALSE
)

controlFactorResponse <- function() {
  dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 25L,
    updateState = FALSE
  )
}

# --- dbarts formula + 2-level factor fits probit, matching x/y bit for
# bit at the same seed (formerly "range not meaningful for factors") ---
set.seed(7)
s.form <- suppressMessages(dbarts(
  yf2 ~ x1 + x2 + x3,
  df,
  control = controlFactorResponse()
))
expect_equal(s.form$model@family, "probit")
expect_true(s.form$control@binary)
set.seed(7)
s.xy <- suppressMessages(dbarts(x, yf2, control = controlFactorResponse()))
expect_equal(s.xy$model@family, "probit")
expect_identical(s.form$run(30L, 30L)$train, s.xy$run(30L, 30L)$train)

# the verdict names the detected type and resolved family
expect_message(
  dbarts(x, yf2, control = controlFactorResponse()),
  "2-level factor response detected, fitting family = \"probit\""
)

# --- logical response == 0/1-numeric probit, bit for bit ---
set.seed(8)
s.log <- suppressMessages(dbarts(x, ylog, control = controlFactorResponse()))
set.seed(8)
s.num <- dbarts(x, as.double(yb), control = controlFactorResponse())
expect_equal(s.log$model@family, "probit")
expect_equal(s.num$model@family, "probit")
expect_identical(s.log$run(30L, 30L)$train, s.num$run(30L, 30L)$train)

# --- character response routes as a factor (probit here, 2 levels) ---
set.seed(7)
s.char <- suppressMessages(dbarts(x, ychar, control = controlFactorResponse()))
expect_equal(s.char$model@family, "probit")
set.seed(7)
s.charForm <- suppressMessages(
  dbarts(ychar ~ x1 + x2 + x3, df, control = controlFactorResponse())
)
expect_identical(s.char$run(30L, 30L)$train, s.charForm$run(30L, 30L)$train)

# --- 2-level ordered factor: binary, so probit (a 3+-level ordered factor
# --- auto-dispatches to ordinal instead; test-ordinal.R covers that) ---
expect_message(
  dbarts(x, ordered(yf2), control = controlFactorResponse()),
  "2-level ordered factor response detected, fitting family = \"probit\""
)

# --- an explicit family that contradicts a factor response errors ---
expect_error(
  dbarts(x, yf2, family = "gaussian", control = controlFactorResponse()),
  "cannot fit a factor response"
)

# --- a 3+-level factor errors informatively in the single-forest
# entry points (never a silent gaussian on the integer level codes) ---
expect_error(dbarts(x, y3, control = controlFactorResponse()), "multinomial")
expect_error(
  dbarts(y3 ~ x1 + x2 + x3, df, control = controlFactorResponse()),
  "multinomial"
)
expect_error(
  bart(
    x,
    y3,
    ndpost = 5L,
    nskip = 5L,
    ntree = 10L,
    nchain = 1L,
    nthread = 1L,
    verbose = FALSE
  ),
  "multinomial"
)
# xbart and rbart_vi resolve an EXPLICIT family, so they refuse under their
# own caller-naming text rather than the family = "auto" arm's
expect_error(
  xbart(
    x,
    y3,
    n.samples = 5L,
    n.reps = 1L,
    n.trees = 10L,
    k = 2,
    n.threads = 1L
  ),
  "xbart does not fit a 3-level factor response"
)
expect_error(
  rbart_vi(
    y3 ~ x1 + x2 + x3,
    df,
    group.by = rep_len(seq_len(4L), n),
    n.samples = 5L,
    n.burn = 0L,
    n.trees = 10L,
    n.chains = 1L,
    n.threads = 1L,
    verbose = FALSE
  ),
  "rbart_vi does not fit a 3-level factor response"
)

# --- bart2: formula + 2-level factor probit matches x/y bit for bit ---
b.args <- list(
  n.trees = 20L,
  n.chains = 1L,
  n.threads = 1L,
  n.burn = 15L,
  n.samples = 15L,
  verbose = FALSE,
  keepTrees = TRUE
)
set.seed(9)
b.form <- suppressMessages(
  do.call(bart2, c(list(yf2 ~ x1 + x2 + x3, data = df), b.args))
)
set.seed(9)
b.xy <- suppressMessages(do.call(bart2, c(list(x, yf2), b.args)))
expect_equal(b.form$family, "probit")
expect_identical(b.form$yhat.train, b.xy$yhat.train)

# --- bart2: auto 3-level factor -> multinomial == explicit, with verdict ---
m.args <- list(
  n.trees = 20L,
  n.chains = 1L,
  n.threads = 1L,
  n.burn = 12L,
  n.samples = 12L,
  verbose = FALSE
)
set.seed(11)
expect_message(
  m.auto <- do.call(bart2, c(list(x, y3), m.args)),
  "3-level factor response detected, fitting family = \"multinomial\""
)
set.seed(11)
m.exp <- do.call(bart2, c(list(x, y3, family = "multinomial"), m.args))
expect_inherits(m.auto, "bartMultinomial")
expect_identical(m.auto$yhat.train, m.exp$yhat.train)

# --- rbart_vi: 2-level factor fits probit (== x/y), 3-level errors ---
g <- sample(1:6, n, replace = TRUE)
r.args <- list(
  group.by = g,
  n.samples = 12L,
  n.burn = 12L,
  n.trees = 20L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE,
  seed = 22L
)
r.form <- suppressMessages(
  do.call(rbart_vi, c(list(yf2 ~ x1 + x2 + x3, df), r.args))
)
r.xy <- suppressMessages(do.call(rbart_vi, c(list(x, yf2), r.args)))
expect_true(r.form$fit[[1L]]$control@binary)
expect_identical(r.form$yhat.train, r.xy$yhat.train)
expect_error(
  do.call(rbart_vi, c(list(y3 ~ x1 + x2 + x3, df), r.args)),
  "multinomial"
)
