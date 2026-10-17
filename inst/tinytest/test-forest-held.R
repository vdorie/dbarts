# A held coefficient is 1 and what it multiplies is the basis (dec-B354): 1
# for a forest with no basis and 0 and 1 for the two levels of a factor,
# wherever the forest stands and at every door, and the forest's sd is then
# the size of what it contributes, exactly. A held forest with no basis that
# states no sd has the size bart() gives its one forest under a continuous
# response (dec-B290) and one unit of the link's error otherwise (dec-B293).
# A held numeric column is held at 1, the column as given (dec-B354).

forest <- dbartsForests$forest
fixed <- dbartsPriors$fixed

set.seed(20261017)
n <- 150L
x1 <- runif(n)
x2 <- runif(n)
z <- as.double(seq_len(n) %in% sample(n, 70L))
zf <- factor(z)
zl <- z == 1
dose <- rnorm(n, 0.05, 0.02)
off <- rnorm(n, 0, 0.1)
w <- rep(1, n)
w[sample(n, 10L)] <- 0
y <- 3 + 2 * sin(3 * x1) + z * (1 + x2) + 20 * dose + rnorm(n, sd = 0.6)
yb <- as.double(y > median(y))
x <- cbind(x1 = x1, x2 = x2)
d <- data.frame(y, yb, x1, x2, zf, zl, dose)

heldControl <- function(...) {
  dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.samples = 10L,
    n.burn = 0L,
    updateState = FALSE,
    seed = 41L,
    ...
  )
}

# --- by kind, at every door and position ---
# A shape is the forests of a model in order: "none" (no basis), "factor",
# "logical" and "dose", a trailing "*" holding the coefficient.
basisValue <- list(factor = zf, logical = zl, dose = dose)
basisCode <- c(factor = "zf", logical = "zl", dose = "dose")
kindOf <- function(entry) sub("[*]$", "", entry)
heldIn <- function(entry) endsWith(entry, "*")
valueForest <- function(entry) {
  kind <- kindOf(entry)
  arguments <- list()
  if (kind != "none") {
    arguments$basis <- basisValue[[kind]]
  }
  if (heldIn(entry)) {
    arguments$amplitude <- fixed()
  }
  do.call(forest, arguments)
}
heldDoors <- list(
  term = function(shape) {
    terms <- vapply(
      shape,
      function(entry) {
        kind <- kindOf(entry)
        sprintf(
          "forest(x1 + x2%s%s)",
          if (kind == "none") "" else paste0(", basis = ", basisCode[[kind]]),
          if (heldIn(entry)) ", amplitude = fixed()" else ""
        )
      },
      ""
    )
    formula <- as.formula(paste("y ~", paste(terms, collapse = " + ")))
    dbarts(formula, d, control = heldControl())
  },
  list.formula = function(shape) {
    dbarts(
      y ~ x1 + x2,
      d,
      forests = lapply(shape, valueForest),
      control = heldControl()
    )
  },
  list.matrix = function(shape) {
    dbarts(x, y, forests = lapply(shape, valueForest), control = heldControl())
  },
  spec = function(shape) {
    spec <- dbartsSpec(
      dbartsData(x, y),
      heldControl(),
      forests = lapply(shape, valueForest)
    )
    new("dbartsSampler", spec$control, spec$model, spec$data)
  },
  data = function(shape) {
    bases <- lapply(shape, function(entry) basisValue[[kindOf(entry)]])
    held <- lapply(shape, function(entry) {
      if (heldIn(entry)) forest(amplitude = fixed()) else forest()
    })
    dbarts(
      dbartsData(x, y, bases = bases),
      forests = held,
      control = heldControl()
    )
  }
)
# the shapes, each read at every door that writes it: a formula puts the
# forest with no basis first wherever its term stands, so a shape with that
# forest elsewhere is a list's
shapes <- list(
  c("none*", "factor", "dose"),
  c("none", "factor*", "dose"),
  c("none", "dose", "factor*"),
  c("none", "dose", "logical*"),
  c("factor*", "dose"),
  c("factor*", "none"),
  c("dose", "none*", "factor*"),
  c("none*", "factor*"),
  c("none", "dose*", "factor"),
  c("none", "factor", "dose*"),
  c("dose*", "none")
)
heldBlock <- function(kind) if (kind %in% c("none", "dose")) 1 else c(0, 1)
for (shape in shapes) {
  plainFirst <- !("none" %in% kindOf(shape[-1L]))
  at <- if (plainFirst) names(heldDoors) else setdiff(names(heldDoors), "term")
  for (door in at) {
    info <- paste(paste(shape, collapse = " "), "at", door)
    sampler <- heldDoors[[door]](shape)
    for (f in seq_along(shape)) {
      prior <- sampler$getLeafPrior(f)
      kind <- kindOf(shape[[f]])
      expect_identical(
        prior$multiplier,
        if (kind %in% c("factor", "logical")) {
          "factor"
        } else if (kind == "dose") {
          "numeric"
        } else {
          "none"
        },
        info = info
      )
      expect_identical(
        is.null(prior[["amplitude"]]),
        !heldIn(shape[[f]]),
        info = info
      )
    }
    heldAt <- which(heldIn(shape))
    before <- lapply(heldAt, function(f) c(sampler$getForestAmplitudes(f)))
    invisible(sampler$run(0L, 50L))
    after <- lapply(heldAt, function(f) c(sampler$getForestAmplitudes(f)))
    expected <- lapply(heldAt, function(f) heldBlock(kindOf(shape[[f]])))
    expect_identical(before, expected, info = info)
    expect_identical(after, expected, info = info)
  }
}

# --- exact: a held forest's own sd is its sd ---
unitOf <- function(sampler) {
  attr(sampler$control, "bartcore.forests", exact = TRUE)$unit
}
families <- list(
  gaussian = list(response = y, family = "gaussian", unit = NA),
  probit = list(response = yb, family = binomial(link = "probit"), unit = 1),
  logistic = list(response = yb, family = "binomial", unit = pi / sqrt(3))
)
for (name in names(families)) {
  family <- families[[name]]
  for (sds in list(c(0.3, 0.6), c(NA, NA))) {
    statement <- function(index, ...) {
      if (is.na(sds[index])) forest(...) else forest(..., sd = sds[index])
    }
    sampler <- dbarts(
      x,
      family$response,
      family = family$family,
      forests = list(
        statement(1L, amplitude = fixed()),
        statement(2L, basis = zf, amplitude = fixed())
      ),
      control = heldControl()
    )
    info <- paste(name, if (is.na(sds[1L])) "default" else "stated")
    for (f in 1:2) {
      prior <- sampler$getLeafPrior(f)
      expect_equal(prior$k.scale, prior$sd, tolerance = 1e-12, info = info)
      if (!is.na(sds[f])) {
        expect_identical(prior$sd, sds[f], info = info)
      }
      expect_identical(
        prior$prior.sd.of,
        if (f == 1L) "forest" else "level difference",
        info = info
      )
    }
    if (is.na(sds[1L])) {
      # a held factor's default is the unit, and under a latent family so is
      # a held forest with no basis
      expect_identical(
        sampler$getLeafPrior(2L)$sd,
        unitOf(sampler),
        info = info
      )
      if (name != "gaussian") {
        expect_identical(unitOf(sampler), family$unit, info = info)
        expect_identical(sampler$getLeafPrior(1L)$sd, family$unit, info = info)
      }
    }
  }
}
# under gaussian a held forest with no basis that states no sd has the size
# bart() gives its one forest, k.scale over k, on the same response; with ten
# rows of weight zero far out too, both reading the rows in the likelihood
outlying <- y
outlying[w == 0] <- 60
for (case in list(
  list(response = y, weights = NULL),
  list(response = outlying, weights = w)
)) {
  held <- dbarts(
    x,
    case$response,
    weights = case$weights,
    forests = list(forest(amplitude = fixed()), forest(basis = zf)),
    control = heldControl()
  )
  single <- dbarts(
    x,
    case$response,
    weights = case$weights,
    control = heldControl()
  )
  onePrior <- single$getLeafPrior()
  expect_equal(
    held$getLeafPrior(1L)$sd,
    onePrior$k.scale / single$getK()[[1L]],
    tolerance = 1e-12
  )
  expect_false(held$getLeafPrior(1L)$sd.stated)
}
# the default is taken from the anchor retained at creation, so neither a new
# response with its scale kept nor a copy made after it moves it
held <- dbarts(
  x,
  y,
  forests = list(forest(amplitude = fixed()), forest(basis = zf)),
  control = heldControl()
)
created <- held$getLeafPrior(1L)
held$setResponse(y + 3 * x2^4, updateScale = FALSE)
expect_identical(held$getLeafPrior(1L)$sd, created$sd)
copied <- held$copy()
expect_identical(copied$getLeafPrior(1L)$sd, created$sd)
expect_identical(copied$getLeafPrior(1L)$k.scale, created$k.scale)
# a write of an sd to a held factor keeps it held, and the sd is exact
held <- dbarts(
  x,
  y,
  forests = list(forest(), forest(basis = zf, amplitude = fixed())),
  control = heldControl()
)
held$setLeafPrior(forests = list(forest(), forest(sd = 0.4)))
written <- held$getLeafPrior(2L)
expect_identical(written$sd, 0.4)
expect_equal(written$k.scale, 0.4, tolerance = 1e-12)
expect_inherits(written[["amplitude"]], "dbartsFixedPrior")
# a held factor's default is one unit at two, three and four forests
for (numForests in 2:4) {
  forests <- c(
    list(forest(), forest(basis = zf, amplitude = fixed())),
    rep(list(forest(basis = zl)), numForests - 2L)
  )
  sampler <- dbarts(x, y, forests = forests, control = heldControl())
  expect_identical(
    sampler$getLeafPrior(2L)$sd,
    unitOf(sampler),
    info = numForests
  )
  expect_equal(
    sampler$getLeafPrior(2L)$k.scale,
    unitOf(sampler),
    tolerance = 1e-12,
    info = numForests
  )
}

# --- what a fit records: a 1 for each held forest, none for a drawn one ---
fitOf <- function(formula) {
  bart(
    formula,
    d,
    n.samples = 5L,
    n.burn = 0L,
    n.chains = 1L,
    n.threads = 1L,
    verbose = FALSE
  )
}
expect_identical(
  fitOf(
    y ~ forest(x1 + x2, amplitude = fixed()) + forest(x1, basis = zf)
  )$fixed$amplitude,
  c(forest1 = 1)
)
expect_identical(
  fitOf(
    y ~ forest(x1 + x2) +
      forest(x1, basis = dose) +
      forest(x1, basis = zf, amplitude = fixed())
  )$fixed$amplitude,
  c(forest3 = 1)
)
expect_identical(
  fitOf(
    y ~ forest(x1, basis = zf, amplitude = fixed()) + forest(x2, basis = dose)
  )$fixed$amplitude,
  c(forest1 = 1)
)
expect_null(
  fitOf(y ~ forest(x1 + x2) + forest(x1, basis = zf))$fixed$amplitude
)

# --- the size, from the engine's prior draws: the held forest's term at a
# row, the second level's for the factor, has sd its stated sd ---
priorSize <- function(sampler, f, row, draws) {
  scale <- sampler$getLeafPrior(f)$response.scale
  values <- vapply(
    seq_len(draws),
    function(draw) {
      sampler$sampleTreesFromPrior(updateState = FALSE)
      sampler$sampleLeafParametersFromPrior(updateState = FALSE)
      scale * sampler$getForestFits(f)[row, 1L]
    },
    0
  )
  sd(values)
}
treatedRow <- which(z == 1)[1L]
sizeCheck <- function(sds, draws, tolerance) {
  for (name in names(families)) {
    family <- families[[name]]
    sampler <- dbarts(
      x,
      family$response,
      family = family$family,
      forests = list(
        forest(amplitude = fixed(), sd = sds[1L]),
        forest(basis = zf, amplitude = fixed(), sd = sds[2L])
      ),
      control = heldControl(n.trees = 20L)
    )
    set.seed(5)
    for (f in 1:2) {
      measured <- priorSize(sampler, f, treatedRow, draws)
      expect_true(
        abs(measured / sds[f] - 1) < tolerance,
        info = paste(name, f, measured)
      )
    }
  }
}
sizeCheck(c(0.8, 0.35), 1500L, 0.1)
if (at_home()) {
  for (sds in list(c(0.2, 1.5), c(1.1, 0.5), c(2.5, 0.15))) {
    sizeCheck(sds, 5000L, 0.04)
  }
}

# --- refused ---
for (call in list(
  quote(dbarts(
    y ~ forest(x1 + x2) +
      forest(x1, basis = zf, amplitude.prior.variance = 0.5),
    d,
    control = heldControl()
  )),
  quote(dbarts(
    x,
    y,
    forests = list(
      forest(),
      forest(basis = zf, amplitude.prior.variance = 0.5)
    ),
    control = heldControl()
  )),
  quote(do.call(forest, list(basis = zf, amplitude.prior.variance = 0.5)))
)) {
  expect_error(
    eval(call),
    "unused argument (amplitude.prior.variance = 0.5)",
    fixed = TRUE
  )
}
expect_error(
  dbarts(
    x,
    y,
    forests = list(forest(), forest(basis = zf, amplitude = fixed(2))),
    control = heldControl()
  ),
  paste0(
    "'amplitude = fixed(2)': a held coefficient is 1, and what it multiplies ",
    "is the basis: nothing, a numeric column as given, or 0 and 1 for the ",
    "two levels of a factor; fixed() takes no other value. Write fixed(), ",
    "state the forest's size with 'sd', and write another coding as the ",
    "basis, as basis = z - 0.5"
  ),
  fixed = TRUE
)

# --- not moved: a drawn numeric column beside a held factor keeps its prior
numericOf <- function(holdFactor) {
  dbarts(
    x,
    y,
    forests = list(
      forest(),
      forest(basis = zf, amplitude = if (holdFactor) fixed()),
      forest(basis = dose)
    ),
    control = heldControl()
  )$getLeafPrior(3L)
}
beside <- numericOf(TRUE)
expect_identical(beside$k.scale, numericOf(FALSE)$k.scale)
expect_identical(beside$prior.sd.of, "forest, per unit of basis (prior median)")
expect_null(beside[["amplitude"]])

# --- a held numeric column: its sd exact, stated or one unit for one sd of
# the column, and the hold kept by a swap, a column of zeros among them ---
for (sd in c(0.3, NA)) {
  for (shape in list(c("none", "dose*"), c("none", "factor", "dose*"))) {
    for (door in c("list.formula", "list.matrix", "spec", "data")) {
      info <- paste(door, length(shape), if (is.na(sd)) "default" else sd)
      sampler <- if (is.na(sd)) {
        heldDoors[[door]](shape)
      } else {
        dbarts(
          x,
          y,
          forests = c(
            lapply(shape[-length(shape)], valueForest),
            list(forest(basis = dose, amplitude = fixed(), sd = sd))
          ),
          control = heldControl()
        )
      }
      f <- length(shape)
      prior <- sampler$getLeafPrior(f)
      unit <- attr(sampler$control, "bartcore.forests", exact = TRUE)$unit
      expect_equal(
        unname(prior$sd),
        if (is.na(sd)) unit / sd(dose) else sd,
        tolerance = 1e-12,
        info = info
      )
      expect_equal(
        prior$k.scale,
        unname(prior$sd),
        tolerance = 1e-12,
        info = info
      )
      expect_identical(
        prior$prior.sd.of,
        "forest, per unit of basis",
        info = info
      )
      for (column in list(dose / 10, rep(0, n))) {
        sampler$setForestBasis(f, column)
        expect_identical(
          sampler$getLeafPrior(f)$k.scale,
          prior$k.scale,
          info = info
        )
        invisible(sampler$run(0L, 5L))
        expect_identical(c(sampler$getForestAmplitudes(f)), 1, info = info)
      }
    }
  }
}
heldDose <- bart(
  y ~ forest(x1 + x2) + forest(x1, basis = dose, amplitude = fixed()),
  d,
  n.samples = 10L,
  n.burn = 0L,
  n.chains = 1L,
  n.threads = 1L,
  seed = 43L,
  verbose = FALSE
)
expect_identical(unname(heldDose$fixed$amplitude), 1)
expect_identical(names(heldDose$fixed$amplitude), "forest2")

# --- another coding is another basis: z - 0.5 held at 1 is the forest's fit
# times -0.5 and 0.5 ---
d$zc <- z - 0.5
centered <- bart(
  y ~ forest(x1 + x2) + forest(x1, basis = zc, amplitude = fixed()),
  d,
  n.samples = 10L,
  n.burn = 0L,
  n.chains = 1L,
  n.threads = 1L,
  seed = 44L,
  verbose = FALSE
)
term <- extract(centered, type = "forest", contribution = TRUE)[,, 2L]
own <- extract(centered, type = "forest")[,, 2L]
expect_equal(term, sweep(own, 2L, z - 0.5, `*`), tolerance = 1e-12)
expect_true(all(centered$glue[, 2L] == 1))

# --- a held 0/1 number and the held factor of the same column, at one stated
# sd, draw the same ---
heldPair <- function(basis) {
  dbarts(
    x,
    y,
    forests = list(
      forest(),
      forest(basis = basis, amplitude = fixed(), sd = 0.7)
    ),
    control = heldControl()
  )$run()$train
}
expect_identical(heldPair(z), heldPair(zf))
expect_error(
  heldDoors$list.matrix(c("none", "dose"))$setForestBasis(
    2L,
    dose,
    updateBasisScale = TRUE
  ),
  "unused argument (updateBasisScale = TRUE)",
  fixed = TRUE
)
