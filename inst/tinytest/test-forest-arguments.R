# forest()'s coefficient and size arguments. A forest's coefficient law is one
# argument, 'amplitude': fixed() holds the coefficient at the value its shape
# gives it and nothing stated draws it. A forest's 'sd' is one unnamed number,
# judged one way at creation, through $setLeafPrior(forests = ) and in the
# front door's normal(sd = ).

forest <- dbartsForests$forest

set.seed(23)
n <- 150L
x <- matrix(runif(n * 3L), n, 3L)
colnames(x) <- paste0("x", 1:3)
z <- rbinom(n, 1L, 0.5)
y <- 2 * sin(pi * x[, 1L]) + z * (1 + x[, 2L]) + rnorm(n, sd = 0.3)

argumentControl <- function() {
  dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 10L,
    n.samples = 5L,
    updateState = FALSE,
    seed = 23L
  )
}
twoForests <- function(...) {
  dbarts(
    x,
    y,
    forests = list(forest(...), forest(basis = ~ factor(z), ...)),
    control = argumentControl()
  )
}

# --- fixed() holds the coefficient ---
held <- twoForests(amplitude = dbartsForests$fixed())
drawn <- twoForests()
heldRun <- held$run(0L, 50L)
drawnRun <- drawn$run(0L, 50L)
expect_equal(held$getForestAmplitudes()[, 1L], c(1, 0, 1))
expect_false(isTRUE(all.equal(drawn$getForestAmplitudes()[, 1L], c(1, 0, 1))))
expect_false(isTRUE(all.equal(heldRun$train, drawnRun$train)))
expect_identical(
  attr(held$control, "bartcore.forests")$params[[2L]][8L],
  0
)
expect_identical(
  attr(drawn$control, "bartcore.forests")$params[[2L]][8L],
  1
)

# --- what 'amplitude' takes (a wrapper's dots hand over values, so the
# constructor is reached through its list there) ---
for (accepted in list(
  quote(fixed),
  quote(fixed()),
  quote(fixed(1)),
  quote(fixed(1L)),
  quote(NULL)
)) {
  expect_silent(
    eval(bquote(dbarts(
      x,
      y,
      forests = list(
        forest(),
        forest(basis = ~ factor(z), amplitude = .(accepted))
      ),
      control = argumentControl()
    )))
  )
}
expect_error(
  twoForests(amplitude = dbartsForests$fixed(2)),
  "a held coefficient takes the value its forest's shape gives it"
)
expect_error(
  twoForests(amplitude = dbartsForests$fixed(c(1, 1))),
  "'value' must be a single positive number"
)
expect_error(
  twoForests(amplitude = dbartsForests$fixed(2)),
  "'amplitude = fixed(2)'",
  fixed = TRUE
)
expect_error(
  twoForests(amplitude = dbartsForests$fixed(TRUE)),
  "invalid object for slot"
)
for (refused in list(quote(dbartsPriors$normal()), "fixed", FALSE)) {
  expect_error(
    eval(bquote(dbarts(
      x,
      y,
      forests = list(
        forest(),
        forest(basis = ~ factor(z), amplitude = .(refused))
      ),
      control = argumentControl()
    ))),
    "a forest's 'amplitude' must be fixed(), which holds its coefficient",
    fixed = TRUE
  )
}
# a model of one forest has no coefficient law
expect_error(
  dbarts(
    x,
    y,
    forests = list(forest(amplitude = fixed())),
    control = argumentControl()
  ),
  "'amplitude' is the law of the coefficient that a model of several forests"
)

# --- update.amplitude is gone, at both doors ---
expect_error(
  dbarts(
    x,
    y,
    forests = list(
      forest(),
      forest(basis = ~ factor(z), update.amplitude = FALSE)
    ),
    control = argumentControl()
  ),
  "unused argument (update.amplitude = FALSE)",
  fixed = TRUE
)
dataFrame <- data.frame(x, y = y, z = z)
expect_error(
  bart(
    y ~ x1 + x2 + forest(x1 + x2, basis = ~z, update.amplitude = FALSE),
    dataFrame,
    n.trees = 5L,
    n.samples = 5L,
    n.burn = 5L,
    n.chains = 1L,
    n.threads = 1L,
    verbose = FALSE
  ),
  "unused argument (update.amplitude = FALSE)",
  fixed = TRUE
)
formulaHeld <- bart(
  y ~ x1 + x2 + forest(x1 + x2, basis = ~z, amplitude = fixed()),
  dataFrame,
  n.trees = 5L,
  n.samples = 5L,
  n.burn = 5L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)
expect_true(inherits(formulaHeld, "bart"))

# --- fixed() is the constructor wherever the call is built and whatever the
# caller has bound to the name ---
built <- as.call(list(
  dbarts::dbarts,
  x,
  y,
  forests = quote(list(
    forest(amplitude = fixed()),
    forest(basis = ~ factor(z), amplitude = fixed())
  )),
  control = argumentControl()
))
bareEnv <- list2env(list(z = z), parent = baseenv())
bareHeld <- eval(built, bareEnv)
bareHeld$run(0L, 5L)
expect_equal(bareHeld$getForestAmplitudes()[, 1L], c(1, 0, 1))
shadowed <- local({
  fixed <- TRUE
  dbarts(
    x,
    y,
    forests = list(
      forest(amplitude = fixed()),
      forest(basis = ~ factor(z), amplitude = fixed())
    ),
    control = argumentControl()
  )
})
shadowed$run(0L, 5L)
expect_equal(shadowed$getForestAmplitudes()[, 1L], c(1, 0, 1))
expect_identical(dbartsForests$fixed, dbartsPriors$fixed)

# --- one verdict on a stated sd at the three places ---
sdValues <- list(
  list("2", FALSE),
  list(TRUE, FALSE),
  list(factor("2"), FALSE),
  list(list(2), FALSE),
  list(matrix(2), FALSE),
  list(as.Date("2026-01-01"), FALSE),
  list(c(a = 1), FALSE),
  list(c(1, 2), FALSE),
  list(NA, FALSE),
  list(NaN, FALSE),
  list(Inf, FALSE),
  list(0, FALSE),
  list(-1, FALSE),
  list(2, TRUE),
  list(2L, TRUE)
)
writer <- twoForests()
for (entry in sdValues) {
  value <- entry[[1L]]
  label <- paste(deparse(value), collapse = "")
  verdicts <- c(
    creation = !inherits(
      try(twoForests(sd = value), silent = TRUE),
      "try-error"
    ),
    writer = !inherits(
      try(
        writer$setLeafPrior(forests = list(forest(), forest(sd = value))),
        silent = TRUE
      ),
      "try-error"
    ),
    front = !inherits(
      try(dbartsPriors$normal(sd = value), silent = TRUE),
      "try-error"
    )
  )
  expect_identical(
    unname(verdicts),
    rep(entry[[2L]], 3L),
    info = label
  )
}

# --- the text for each ---
sdMessage <- function(value, place = c("creation", "writer", "front")) {
  place <- match.arg(place)
  conditionMessage(tryCatch(
    switch(
      place,
      creation = twoForests(sd = value),
      writer = writer$setLeafPrior(
        forests = list(forest(), forest(sd = value))
      ),
      front = dbartsPriors$normal(sd = value)
    ),
    error = identity
  ))
}
for (place in c("creation", "writer")) {
  expect_match(sdMessage("2", place), "must be a number, not a string")
  expect_match(sdMessage(TRUE, place), "must be a number, not a logical")
  expect_match(sdMessage(factor("2"), place), "must be a number, not a factor")
  expect_match(
    sdMessage(as.Date("2026-01-01"), place),
    "must be a number, not a Date"
  )
  expect_match(sdMessage(list(2), place), "must be a number, not a list")
  expect_match(sdMessage(matrix(2), place), "must be a number, not a matrix")
  expect_match(
    sdMessage(c(dose = 1), place),
    "forest 'sd' must not be named (\"dose\"): it is one number",
    fixed = TRUE
  )
  expect_match(
    sdMessage(c(1, 2), place),
    "forest 'sd' must be a single number, not 2: a forest states one sd",
    fixed = TRUE
  )
  expect_match(
    sdMessage(NA_real_, place),
    "forest 'sd' must not be NA; leave it out for the default",
    fixed = TRUE
  )
  expect_match(sdMessage(NaN, place), "forest 'sd' must not be NA")
  for (bad in list(0, -1, Inf)) {
    expect_match(
      sdMessage(bad, place),
      "forest 'sd' must be positive and finite"
    )
  }
}
expect_match(
  sdMessage(TRUE, "front"),
  "'sd' must be a number or invchi(), not a logical",
  fixed = TRUE
)
expect_match(
  sdMessage(c(a = 1), "front"),
  "'sd' must not be named (\"a\")",
  fixed = TRUE
)
# the front door keeps the texts it has for what it already refused
expect_match(sdMessage("2", "front"), "unlike 'k' it takes no string form")
expect_match(sdMessage(c(1, 2), "front"), "'sd' must be a single number")
expect_match(sdMessage(0, "front"), "'sd' must be positive")
