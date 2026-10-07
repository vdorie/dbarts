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
for (basis in list(quote(~ factor(z)), quote(I(1 - z) + z), cbind(1 - z, z))) {
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
# a held FIRST forest of an all-multiplied model is refused as well; no forest
# of that model is without a basis, so the count and the tree prior are
# stated on a forest and the control names no count
expect_error(
  dbarts(
    x,
    y,
    forests = list(
      forest(
        basis = ~z,
        amplitude = fixed(),
        n.trees = 10L,
        base = 0.95,
        power = 2
      ),
      forest(basis = ~ factor(z))
    ),
    control = dbartsControl(
      n.chains = 1L,
      n.threads = 1L,
      n.samples = 5L,
      updateState = FALSE,
      seed = 23L
    )
  ),
  "forest 1: amplitude = fixed() on a basis of one numeric column",
  fixed = TRUE
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
  "'amplitude = fixed(2)': a held coefficient is 1 for a forest with no basis, and 0 for the first level of a factor and 1 for the others; fixed() takes no other value here. Write fixed(), and state the forest's size with 'sd'",
  fixed = TRUE
)
for (bad in list(c(1, 1), NA_real_, NaN)) {
  expect_error(
    twoForests(amplitude = dbartsPriors$fixed(bad)),
    "'value' must be a single positive number"
  )
}
expect_error(twoForests(amplitude = dbartsPriors$fixed(TRUE)), ".")
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

# --- fixed resolves as interactions and blocks do: by the vocabulary the
# argument is evaluated in, a call being the constructor and a bare name the
# caller has bound being the caller's value ---
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
heldBy <- function(forests) {
  sampler <- eval(
    bquote(dbarts(x, y, forests = .(forests), control = argumentControl())),
    parent.frame()
  )
  heldParams(sampler, 2L)
}
# a call is the constructor whatever the caller has bound
local({
  fixed <- TRUE
  expect_identical(
    heldBy(quote(list(
      forest(),
      forest(basis = ~ factor(z), amplitude = fixed())
    ))),
    0
  )
})
# a bare name the caller has bound is the caller's value: NULL draws, anything
# else is refused as any value that is not fixed() is
local({
  fixed <- NULL
  expect_identical(
    heldBy(quote(list(
      forest(),
      forest(basis = ~ factor(z), amplitude = fixed)
    ))),
    1
  )
  expect_identical(
    dbarts(
      y ~ x1 + x2 + forest(x1 + x2, basis = ~ factor(z), amplitude = fixed),
      dataFrame,
      control = argumentControl()
    )$control |>
      attr("bartcore.forests") |>
      (\(info) info$params[[2L]][8L])(),
    1
  )
})
local({
  fixed <- FALSE
  expect_error(
    heldBy(quote(list(
      forest(),
      forest(basis = ~ factor(z), amplitude = fixed)
    ))),
    "a forest's 'amplitude' must be fixed()",
    fixed = TRUE
  )
  expect_error(
    dbarts(
      y ~ x1 + x2 + forest(x1 + x2, basis = ~ factor(z), amplitude = fixed),
      dataFrame,
      control = argumentControl()
    ),
    "a forest's 'amplitude' must be fixed()",
    fixed = TRUE
  )
})
# a bare fixed with nothing bound is the constructor, nested or positional
for (hold in c(TRUE, FALSE)) {
  expect_identical(
    heldBy(quote(list(
      forest(),
      forest(basis = ~ factor(z), amplitude = if (hold) fixed)
    ))),
    if (hold) 0 else 1
  )
}
expect_identical(
  heldBy(quote(list(
    forest(),
    forest(basis = ~ factor(z), amplitude = (fixed))
  ))),
  0
)
# and not by position: only a forest's predictors are given unnamed
expect_error(
  heldBy(quote(list(
    forest(),
    forest(NULL, ~ factor(z), NULL, NULL, NULL, NULL, fixed())
  ))),
  "forest() takes one unnamed argument",
  fixed = TRUE
)
# inside a function literal
expect_identical(
  heldBy(quote(lapply(1:2, function(i) {
    if (i == 1L) forest() else forest(basis = ~ factor(z), amplitude = fixed())
  }))),
  0
)
# an argument named 'amplitude' of another function is that function's
local({
  wave <- function(amplitude) amplitude
  fixed <- 2
  sampler <- dbarts(
    x,
    y,
    forests = list(
      forest(),
      forest(basis = ~ factor(z), sd = wave(amplitude = fixed))
    ),
    control = argumentControl()
  )
  expect_identical(
    attr(sampler$control, "bartcore.forests")$params[[2L]][4L],
    2
  )
})
# forests forwarded through a wrapper's dots resolve it as they resolve
# interactions(), one wrapper deep and two
wrap <- function(...) dbarts(x, y, ..., control = argumentControl())
wrapTwice <- function(...) wrap(...)
for (wrapper in list(wrap, wrapTwice)) {
  sampler <- wrapper(
    forests = list(
      forest(interactions = interactions(max.order = 2L)),
      forest(basis = ~ factor(z), amplitude = fixed())
    )
  )
  expect_identical(heldParams(sampler, 2L), 0)
}
# a constructor missing where it is forced points at where it lives
expect_error(
  twoForests(amplitude = fixed()),
  "dbartsPriors$fixed",
  fixed = TRUE
)
# the hint is for fixed alone: another constructor reads as R's own message
expect_identical(
  tryCatch(twoForests(sd = chisq()), error = conditionMessage),
  "could not find function \"chisq\""
)
# the other forest arguments do not take fixed
for (name in c("interactions", "blocks", "monotone")) {
  expect_match(
    tryCatch(
      eval(bquote(
        dbarts(
          x,
          y,
          ..(setNames(list(quote(fixed())), name)),
          control = argumentControl()
        ),
        splice = TRUE
      )),
      error = conditionMessage
    ),
    "could not find function \"fixed\"",
    fixed = TRUE,
    info = name
  )
}
# the writer resolves it too, and refuses it as fixed at creation
writerHeld <- twoForests()
expect_error(
  writerHeld$setLeafPrior(
    forests = list(forest(), forest(amplitude = fixed()))
  ),
  "'amplitude' is fixed at creation"
)
expect_error(
  writerHeld$setLeafPrior(forests = list(forest(), forest(amplitude = fixed))),
  "'amplitude' is fixed at creation"
)

# --- a held forest keeps the width it was created with ---
heldSwap <- function() {
  dbarts(
    x,
    y,
    forests = list(
      forest(amplitude = fixed()),
      forest(basis = ~ factor(z), amplitude = fixed())
    ),
    control = argumentControl()
  )
}
swapped <- heldSwap()
twin <- heldSwap()
expect_error(
  swapped$setForestBasis(2L, x[, 1L]),
  "$setForestBasis cannot change the width of forest 2's basis (2 to 1): its coefficients are held (amplitude = fixed()), and the held value is defined for that width only; make a new sampler",
  fixed = TRUE
)
expect_identical(swapped$data@bases, twin$data@bases)
# the stored description too; the record of a basis written as code holds
# the place its formula was written, which the two builds do not share
describedBy <- function(sampler) {
  described <- attr(sampler$control, "bartcore.forests")
  described$basisTerms <- lapply(described$basisTerms, function(record) {
    record[c("label", "levels", "xlev", "rows")]
  })
  described
}
expect_identical(describedBy(swapped), describedBy(twin))
expect_identical(describedBy(swapped)$basisTerms[[2L]]$label, "factor(z)")
expect_identical(swapped$run(0L, 5L), twin$run(0L, 5L))
expect_equal(swapped$getForestAmplitudes()[, 1L], c(1, 0, 1))
expect_error(
  swapped$setForestBasis(2L, ~z),
  "$setForestBasis cannot change the width of forest 2's basis (2 to 1)",
  fixed = TRUE
)
expect_identical(swapped$data@bases, twin$data@bases)
expect_identical(describedBy(swapped), describedBy(twin))
# a held forest created without a basis takes none
plainSwap <- heldSwap()
expect_error(
  plainSwap$setForestBasis(1L, x[, 1L]),
  "$setForestBasis cannot change the width of forest 1's basis (0 to 1): its coefficients are held (amplitude = fixed()), and the held value is defined for that width only; make a new sampler",
  fixed = TRUE
)
plainTwin <- heldSwap()
expect_null(plainSwap$data@bases[[1L]])
expect_identical(plainSwap$data@bases, plainTwin$data@bases)
expect_identical(describedBy(plainSwap), describedBy(plainTwin))
expect_identical(plainSwap$run(0L, 5L), plainTwin$run(0L, 5L))
expect_equal(plainSwap$getForestAmplitudes()[, 1L], c(1, 0, 1))
# a drawn forest swaps as before, and a held one to a factor or two columns
drawnSwap <- twoForests()
drawnSwap$setForestBasis(2L, x[, 1L])
expect_equal(ncol(drawnSwap$data@bases[[2L]]), 1L)
swapped$setForestBasis(2L, factor(z))
swapped$setForestBasis(2L, cbind(x[, 1L], x[, 2L]))
swapped$setForestBasis(2L, factor(z))
expect_equal(swapped$run(0L, 3L)$train |> dim() |> length(), 2L)
# a held forest copies, saves and reloads
copied <- swapped$copy()
copied$run(0L, 3L)
expect_equal(copied$getForestAmplitudes()[, 1L], c(1, 0, 1))
savedFile <- tempfile(fileext = ".rds")
swapped$storeState()
saveRDS(swapped, savedFile)
reloaded <- readRDS(savedFile)
reloaded$run(0L, 3L)
expect_equal(reloaded$getForestAmplitudes()[, 1L], c(1, 0, 1))
unlink(savedFile)

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
  list(c(a = 1), c(FALSE, FALSE, TRUE)),
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
    rep_len(entry[[2L]], 3L),
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
    "forest 'sd' must not be named (\"dose\"): it is one number, for every column of a basis; drop the name with unname()",
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
# a number indexed out of a named vector is a number: the name is dropped
expect_identical(sdMessage(c(a = 1), "front"), "accepted")
expect_identical(dbartsPriors$normal(sd = c(a = 1.5))@prior.sd, 1.5)
expect_match(sdMessage(c(a = -1), "front"), "'sd' must be positive")
frontSampler <- dbarts(x, y, control = argumentControl())
frontSampler$setLeafPrior(dbartsPriors$normal(sd = c(a = 0.7)))
expect_equal(
  unname(frontSampler$getLeafPrior()$k.scale / frontSampler$getK()[1L]),
  0.7
)
expect_match(sdMessage(2 + 0i, "front"), "not a complex")
expect_match(sdMessage(as.raw(2), "front"), "not a raw")
expect_match(sdMessage(Sys.time(), "front"), "not a date-time")
# what the front door documents it takes
expect_identical(sdMessage(dbartsPriors$invchi(3, 1), "front"), "accepted")
# the front door keeps the texts it has for what it already refused
expect_match(sdMessage("2", "front"), "unlike 'k' it takes no string form")
expect_match(sdMessage(c(1, 2), "front"), "'sd' must be a single number")
expect_match(sdMessage(0, "front"), "'sd' must be positive")
