# forest()'s coefficient and size arguments. A forest's coefficient law is one
# argument, 'amplitude': fixed() holds the coefficient and nothing stated draws
# it. A forest's 'sd' is one unnamed number, judged one way at creation,
# through $setLeafPrior(forests = ) and in the front door's normal(sd = ).

forest <- dbartsForests$forest

set.seed(23)
n <- 150L
x <- matrix(runif(n * 3L), n, 3L)
colnames(x) <- paste0("x", 1:3)
z <- rbinom(n, 1L, 0.5)
y <- 2 * sin(pi * x[, 1L]) + z * (1 + x[, 2L]) + rnorm(n, sd = 0.3)
dataFrame <- data.frame(x, y = y, z = z)

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
heldParams <- function(sampler, forest) {
  attr(sampler$control, "bartcore.forests")$params[[forest]][8L]
}
termControl <- list(
  n.trees = 5L,
  n.samples = 5L,
  n.burn = 5L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)

# --- fixed() holds the coefficient ---
held <- twoForests(amplitude = dbartsPriors$fixed())
drawn <- twoForests()
heldRun <- held$run(0L, 50L)
drawnRun <- drawn$run(0L, 50L)
expect_equal(held$getForestAmplitudes()[, 1L], c(1, 0, 1))
expect_false(isTRUE(all.equal(drawn$getForestAmplitudes()[, 1L], c(1, 0, 1))))
expect_false(isTRUE(all.equal(heldRun$train, drawnRun$train)))
expect_identical(heldParams(held, 2L), 0)
expect_identical(heldParams(drawn, 2L), 1)

# the shapes a hold is supported on, each with the values it holds: a forest
# with no basis at 1, a forest on a factor or on the two columns of a
# complement pair at 0 for the first column and 1 for the other
for (basis in list(quote(~ factor(z)), quote(~ cbind(1 - z, z)))) {
  sampler <- eval(bquote(dbarts(
    x,
    y,
    forests = list(
      forest(amplitude = fixed()),
      forest(basis = .(basis), amplitude = fixed())
    ),
    control = argumentControl()
  )))
  sampler$run(0L, 10L)
  expect_equal(
    sampler$getForestAmplitudes()[, 1L],
    c(1, 0, 1),
    info = deparse(basis)
  )
}
# a held forest beside a drawn one, and the reverse
plainHeld <- dbarts(
  x,
  y,
  forests = list(forest(amplitude = fixed()), forest(basis = ~ factor(z))),
  control = argumentControl()
)
plainHeld$run(0L, 10L)
expect_equal(plainHeld$getForestAmplitudes()[1L, 1L], 1)
expect_false(isTRUE(all.equal(
  plainHeld$getForestAmplitudes()[2:3, 1L],
  c(0, 1)
)))
factorHeld <- dbarts(
  x,
  y,
  forests = list(forest(), forest(basis = ~ factor(z), amplitude = fixed())),
  control = argumentControl()
)
factorHeld$run(0L, 10L)
expect_equal(factorHeld$getForestAmplitudes()[2:3, 1L], c(0, 1))
termHeld <- dbarts(
  y ~ x1 + x2 + forest(x1 + x2, basis = ~ factor(z), amplitude = fixed()),
  dataFrame,
  control = argumentControl()
)
termHeld$run(0L, 10L)
expect_equal(termHeld$getForestAmplitudes()[2:3, 1L], c(0, 1))

# --- a hold on one numeric column would multiply the forest by zero ---
oneColumn <- "amplitude = fixed\\(\\) on a basis of one numeric column is not supported yet"
expect_error(
  dbarts(
    x,
    y,
    forests = list(forest(), forest(basis = ~z, amplitude = fixed())),
    control = argumentControl()
  ),
  oneColumn
)
expect_error(
  dbarts(
    x,
    y,
    forests = list(forest(), forest(basis = z, amplitude = fixed())),
    control = argumentControl()
  ),
  "forest 2: amplitude = fixed() on a basis of one numeric column",
  fixed = TRUE
)
expect_error(
  dbartsSpec(
    dbartsData(x, y),
    argumentControl(),
    forests = list(forest(), forest(basis = ~z, amplitude = fixed()))
  ),
  oneColumn
)
expect_error(
  dbarts(
    y ~ x1 + x2 + forest(x1 + x2, basis = ~z, amplitude = fixed()),
    dataFrame,
    control = argumentControl()
  ),
  oneColumn
)
expect_error(
  do.call(
    bart,
    c(
      list(y ~ x1 + x2 + forest(x1 + x2, basis = ~z, amplitude = fixed())),
      list(dataFrame),
      termControl
    )
  ),
  oneColumn
)
# the same column drawn is accepted
expect_silent(dbarts(
  x,
  y,
  forests = list(forest(), forest(basis = ~z)),
  control = argumentControl()
))

# --- what 'amplitude' takes, and what each spelling holds (a wrapper's dots
# hand over values, so the constructor is reached through dbartsPriors there) ---
for (entry in list(
  list(quote(fixed), 0),
  list(quote(fixed()), 0),
  list(quote(fixed(1)), 0),
  list(quote(fixed(1L)), 0),
  list(quote(NULL), 1)
)) {
  sampler <- eval(bquote(dbarts(
    x,
    y,
    forests = list(
      forest(),
      forest(basis = ~ factor(z), amplitude = .(entry[[1L]]))
    ),
    control = argumentControl()
  )))
  expect_identical(
    heldParams(sampler, 2L),
    entry[[2L]],
    info = deparse(entry[[1L]])
  )
}
expect_error(
  twoForests(amplitude = dbartsPriors$fixed(2)),
  "a held coefficient takes the value its forest's shape gives it"
)
expect_error(
  twoForests(amplitude = dbartsPriors$fixed(2)),
  "'amplitude = fixed(2)'",
  fixed = TRUE
)
for (bad in list(c(1, 1), NA_real_, NaN)) {
  expect_error(
    twoForests(amplitude = dbartsPriors$fixed(bad)),
    "'value' must be a single positive number"
  )
}
expect_error(twoForests(amplitude = dbartsPriors$fixed(TRUE)))
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
# a model of one forest has no coefficient to hold, and no forest size
expect_error(
  dbarts(
    x,
    y,
    forests = list(forest(amplitude = fixed())),
    control = argumentControl()
  ),
  "this model has one forest, which has no coefficient to hold; 'amplitude' needs a model of several forests",
  fixed = TRUE
)
expect_error(
  dbarts(x, y, forests = list(forest(sd = 2)), control = argumentControl()),
  "this model has one forest, so its size is the fitting function's leaf.prior = normal(sd = ), not forest(sd = )",
  fixed = TRUE
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
expect_error(
  do.call(
    bart,
    c(
      list(
        y ~ x1 + x2 + forest(x1 + x2, basis = ~z, update.amplitude = FALSE),
        dataFrame
      ),
      termControl
    )
  ),
  "unused argument (update.amplitude = FALSE)",
  fixed = TRUE
)

# --- fixed() is the constructor wherever the call is built and whatever the
# caller has bound to the name; in the argument only ---
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
# the name is the constructor nowhere else: a column called 'fixed' is a column
frame <- data.frame(dataFrame, fixed = z)
columnBasis <- dbarts(
  y ~ x1 + x2 + forest(x1 + x2, basis = ~ factor(fixed), amplitude = fixed()),
  frame,
  control = argumentControl()
)
columnBasis$run(0L, 5L)
expect_equal(columnBasis$getForestAmplitudes()[2:3, 1L], c(0, 1))
expect_error(
  dbarts(x, y, interactions = fixed(), control = argumentControl()),
  "fixed"
)

# --- one verdict on a stated sd at the three places ---
sdValues <- list(
  list("2", FALSE),
  list(TRUE, FALSE),
  list(factor("2"), FALSE),
  list(list(2), FALSE),
  list(matrix(2), FALSE),
  list(as.Date("2026-01-01"), FALSE),
  list(Sys.time(), FALSE),
  list(as.difftime(2, units = "secs"), FALSE),
  list(2 + 0i, FALSE),
  list(as.raw(2), FALSE),
  list(c(a = 1), FALSE),
  list(c(a = -1), FALSE),
  list(c(a = 0), FALSE),
  list(numeric(0), FALSE),
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
  tryCatch(
    {
      switch(
        place,
        creation = twoForests(sd = value),
        writer = writer$setLeafPrior(
          forests = list(forest(), forest(sd = value))
        ),
        front = dbartsPriors$normal(sd = value)
      )
      "accepted"
    },
    error = conditionMessage
  )
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
    "forest 'sd' must be a single number, not a vector of length 2: a forest states one sd",
    fixed = TRUE
  )
  expect_match(
    sdMessage(numeric(0), place),
    "not a vector of length 0",
    fixed = TRUE
  )
  expect_match(
    sdMessage(Sys.time(), place),
    "must be a number, not a date-time"
  )
  expect_match(
    sdMessage(as.difftime(2, units = "secs"), place),
    "must be a number, not a time difference"
  )
  expect_match(sdMessage(2 + 0i, place), "must be a number, not a complex")
  expect_match(sdMessage(as.raw(2), place), "must be a number, not a raw")
  expect_match(
    sdMessage(new.env(), place),
    "must be a number, not an environment"
  )
  expect_match(sdMessage(~x, place), "must be a number, not a formula")
  expect_match(sdMessage(mean, place), "must be a number, not a function")
  expect_match(
    sdMessage(dbartsPriors$invchi(3, 1), place),
    "must be a number, not invchi(): a law on a forest's sd is not supported yet",
    fixed = TRUE
  )
  # a named number is refused as named whatever its value
  for (bad in list(c(a = -1), c(a = 0))) {
    expect_match(sdMessage(bad, place), "must not be named")
  }
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
  "'sd' must not be named (\"a\"): it is one number",
  fixed = TRUE
)
expect_match(sdMessage(c(a = -1), "front"), "must not be named")
expect_match(sdMessage(2 + 0i, "front"), "not a complex")
expect_match(sdMessage(as.raw(2), "front"), "not a raw")
expect_match(sdMessage(Sys.time(), "front"), "not a date-time")
# what the front door documents it takes
expect_identical(sdMessage(dbartsPriors$invchi(3, 1), "front"), "accepted")
# the front door keeps the texts it has for what it already refused
expect_match(sdMessage("2", "front"), "unlike 'k' it takes no string form")
expect_match(sdMessage(c(1, 2), "front"), "'sd' must be a single number")
expect_match(sdMessage(0, "front"), "'sd' must be positive")
