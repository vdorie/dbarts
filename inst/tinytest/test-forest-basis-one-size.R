# A forest multiplies one numeric column, or a factor, a character or a
# logical vector of two levels. Three or more levels, and several numeric
# columns however they are written, are refused by name at every door, drawn
# or held, with the form that works: a forest for each level but the first,
# or for each column. The size is judged on the rows the fit keeps, after a
# level they empty is dropped, at a creation; afterwards no forest changes
# width. A basis of several terms is still read and named before it is
# refused, so that the names stay those a later version gives its columns.

forest <- dbartsForests$forest
fixed <- dbartsPriors$fixed

set.seed(20261017)
n <- 120L
x1 <- runif(n)
x2 <- runif(n)
z <- rep_len(c(0, 1), n)
g <- sample(c("a", "b", "c"), n, replace = TRUE)
arm <- factor(g)
armOrdered <- factor(g, ordered = TRUE)
dose <- runif(n, 0.5, 3)
age <- rnorm(n, 50, 10)
y <- sin(3 * x1) + z * (0.5 + x2) + 0.2 * dose + rnorm(n, sd = 0.3)
x <- cbind(x1 = x1, x2 = x2)
d <- data.frame(y, x1, x2, z, g, arm, armOrdered, dose, age)
B <- cbind(dose, age)

sizeControl <- function() {
  dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 5L,
    n.samples = 3L,
    n.burn = 0L,
    updateState = FALSE,
    seed = 17L
  )
}
# the second forest, its basis written as code, drawn or held
codeForest <- function(code, held) {
  if (held) {
    bquote(forest(x1, basis = .(code), amplitude = fixed()))
  } else {
    bquote(forest(x1, basis = .(code)))
  }
}
# a basis written as code at the five doors that read code
codeDoors <- list(
  term.dbarts = function(code, held, ...) {
    eval(bquote(dbarts(
      y ~ x1 + x2 + .(codeForest(code, held)),
      d,
      control = sizeControl(),
      ...
    )))
  },
  term.bart = function(code, held, ...) {
    eval(bquote(bart(
      y ~ x1 + x2 + .(codeForest(code, held)),
      d,
      n.trees = 5L,
      n.samples = 3L,
      n.burn = 0L,
      n.chains = 1L,
      n.threads = 1L,
      verbose = FALSE,
      samplerOnly = TRUE,
      ...
    )))
  },
  list.formula = function(code, held, ...) {
    eval(bquote(dbarts(
      y ~ x1 + x2,
      d,
      forests = list(forest(), .(codeForest(code, held))),
      control = sizeControl(),
      ...
    )))
  },
  list.matrix = function(code, held, ...) {
    eval(bquote(dbarts(
      x,
      y,
      forests = list(forest(), .(codeForest(code, held))),
      control = sizeControl(),
      ...
    )))
  },
  spec = function(code, held, ...) {
    eval(bquote(dbartsSpec(
      dbartsData(x, y, ...),
      sizeControl(),
      forests = list(forest(), .(codeForest(code, held)))
    )))
  }
)
# a basis handed over as a value at the four doors that take one
valueForests <- function(value, held) {
  list(
    forest(),
    do.call(
      forest,
      c(list("x1", basis = value), if (held) list(amplitude = fixed()))
    )
  )
}
valueDoors <- list(
  list.formula = function(value, held, ...) {
    dbarts(
      y ~ x1 + x2,
      d,
      forests = valueForests(value, held),
      control = sizeControl(),
      ...
    )
  },
  list.matrix = function(value, held, ...) {
    dbarts(
      x,
      y,
      forests = valueForests(value, held),
      control = sizeControl(),
      ...
    )
  },
  spec = function(value, held, ...) {
    dbartsSpec(
      dbartsData(x, y, ...),
      sizeControl(),
      forests = valueForests(value, held)
    )
  },
  data = function(value, held, ...) {
    dbarts(
      dbartsData(x, y, bases = list(NULL, value), ...),
      forests = if (held) list(forest(), forest(amplitude = fixed())),
      control = sizeControl()
    )
  }
)
refusal <- function(expr) {
  tryCatch(
    {
      force(expr)
      NA_character_
    },
    error = conditionMessage
  )
}

## --- the loop: every refused object, at every door it can be written at ----
# %1$s is the basis as written, its label, and %2$s the factor's one term
levelsText <- paste0(
  "forest 2 ('%1$s') has a basis of 3 levels ('a', 'b', 'c'); in this ",
  "version a factor, a character or a logical basis has two levels. Give ",
  "each level but the first a forest of its own, as forest(x1, basis = %2$s ",
  "== \"b\") + forest(x1, basis = %2$s == \"c\")"
)
columnsText <- function(label, names, advice) {
  paste0(
    "forest 2 ('",
    label,
    "') has a basis of ",
    length(names),
    " numeric columns (",
    toString(sQuote(names, FALSE)),
    "); in this version a basis is one numeric column, or a factor, a ",
    "character or a logical vector of two levels. ",
    advice
  )
}
indicatorWords <- paste0(
  "If the columns are a factor's indicators, give the factor itself; ",
  "otherwise give each column a forest of its own"
)
refusedCode <- list(
  arm = sprintf(levelsText, "arm", "arm"),
  armOrdered = sprintf(levelsText, "armOrdered", "armOrdered"),
  g = sprintf(levelsText, "g", "g"),
  # the factor's term without the constant column a basis has none of anyway
  "0 + arm" = sprintf(levelsText, "0 + arm", "arm"),
  "arm - 1" = sprintf(levelsText, "arm - 1", "arm"),
  "arm + 0" = sprintf(levelsText, "arm + 0", "arm"),
  "-1 + arm" = sprintf(levelsText, "-1 + arm", "arm"),
  "dose + age" = columnsText(
    "dose + age",
    c("dose", "age"),
    "Give each column a forest of its own, as forest(x1, basis = dose) + forest(x1, basis = age)"
  ),
  "splines::ns(dose, 3)" = columnsText(
    "splines::ns(dose, 3)",
    paste0("splines::ns(dose, 3)", 1:3),
    indicatorWords
  ),
  "poly(dose, 2)" = columnsText(
    "poly(dose, 2)",
    paste0("poly(dose, 2)", 1:2),
    indicatorWords
  ),
  "1 + dose" = paste0(
    "forest 2 ('1 + dose') has a basis of 2 numeric columns ('(Intercept)', ",
    "'dose'): '1 +' in a basis adds a constant column, and in this version a ",
    "basis is one numeric column. Write the column alone, as basis = dose"
  ),
  B = columnsText("B", c("Bdose", "Bage"), indicatorWords)
)
for (text in names(refusedCode)) {
  code <- str2lang(text)
  for (held in c(FALSE, TRUE)) {
    for (door in names(codeDoors)) {
      info <- paste(text, door, if (held) "held" else "drawn")
      expect_identical(
        refusal(codeDoors[[door]](code, held)),
        refusedCode[[text]],
        info = info
      )
    }
  }
}
# the example a text gives is the caller's own code, and runs as written
exampleOf <- function(text) {
  sub(".*, as ", "", text)
}
for (text in c(
  "arm",
  "g",
  "0 + arm",
  "arm - 1",
  "arm + 0",
  "-1 + arm",
  "dose + age"
)) {
  example <- exampleOf(refusedCode[[text]])
  formula <- str2lang(paste("y ~ forest(x1 + x2) +", example))
  sampler <- eval(bquote(dbarts(.(formula), d, control = sizeControl())))
  expect_identical(length(sampler$data@bases), 3L, info = text)
  expect_identical(
    vapply(sampler$data@bases[2:3], NCOL, 0L),
    if (text == "dose + age") c(1L, 1L) else c(2L, 2L),
    info = text
  )
}
expect_identical(
  ncol(
    dbarts(
      y ~ x1 + x2 + forest(x1, basis = dose),
      d,
      control = sizeControl()
    )$data@bases[[2L]]
  ),
  1L
)

# a value has no code of the caller's: its text shows no forest() call and
# says what to write in words, a factor's indicator columns being told to
# hand the factor over
valueLevels <- paste0(
  "forest 2 has a basis of 3 levels ('a', 'b', 'c'); in this version a ",
  "factor, a character or a logical basis has two levels. Give each level ",
  "but the first a forest of its own, its basis a logical vector that is ",
  "TRUE on that level's rows"
)
valueColumns <- function(names = NULL) {
  paste0(
    "forest 2 has a basis of ",
    if (is.null(names)) 2L else length(names),
    " numeric columns",
    if (!is.null(names)) paste0(" (", toString(sQuote(names, FALSE)), ")"),
    "; in this version a basis is one numeric column, or a factor, a ",
    "character or a logical vector of two levels. ",
    indicatorWords
  )
}
gIndicators <- model.matrix(~ g - 1)
refusedValues <- list(
  list("a factor of three levels", arm, valueLevels),
  list("an ordered factor of three levels", armOrdered, valueLevels),
  list("a character vector of three values", g, valueLevels),
  list("two numeric columns", B, valueColumns(c("dose", "age"))),
  list("indicator columns", cbind(1 - z, z), valueColumns()),
  list("a constant beside a column", cbind(1, dose), valueColumns()),
  list(
    "three indicator columns",
    gIndicators,
    valueColumns(colnames(gIndicators))
  ),
  list("a level omitted", model.matrix(~g)[, -1L], valueColumns(c("gb", "gc")))
)
for (object in refusedValues) {
  for (held in c(FALSE, TRUE)) {
    for (door in names(valueDoors)) {
      info <- paste(object[[1L]], door, if (held) "held" else "drawn")
      said <- refusal(valueDoors[[door]](object[[2L]], held))
      expect_identical(said, object[[3L]], info = info)
      expect_false(grepl("forest(", said, fixed = TRUE), info = info)
    }
  }
}

# a forest is named in the text by its label, a list name before its code,
# at every place in the list
expect_identical(
  refusal(dbarts(
    y ~ x1 + x2,
    d,
    forests = list(forest(), treat = do.call(forest, list("x1", basis = arm))),
    control = sizeControl()
  )),
  sub("forest 2", "forest 2 ('treat')", valueLevels, fixed = TRUE)
)
expect_identical(
  refusal(dbarts(
    y ~ x1 + x2,
    d,
    forests = list(forest(), treat = forest(x1, basis = arm)),
    control = sizeControl()
  )),
  sub("('arm')", "('treat')", refusedCode$arm, fixed = TRUE)
)
# every forest's size is judged, not the second's alone: a third forest as
# a value at dbartsSpec(), which builds no data object, and as code
expect_identical(
  refusal(dbartsSpec(
    dbartsData(x, y),
    sizeControl(),
    forests = list(
      forest(),
      do.call(forest, list(basis = dose)),
      do.call(forest, list(basis = arm))
    )
  )),
  sub("forest 2", "forest 3", valueLevels, fixed = TRUE)
)
expect_identical(
  refusal(dbarts(
    y ~ forest(x1 + x2) + forest(x1, basis = dose) + forest(x1, basis = arm),
    d,
    control = sizeControl()
  )),
  sub("forest 2", "forest 3", refusedCode$arm, fixed = TRUE)
)

## --- what stays: one numeric column, or two levels ------------------------
two <- factor(z)
doseMatrix <- cbind(dose)
stayCode <- c(
  "two",
  "as.character(two)",
  "z == 1",
  "arm == \"b\"",
  "dose",
  "doseMatrix",
  "poly(dose, 1)",
  "scale(age)",
  "dose:age",
  "I(dose * age)",
  "0 + dose",
  "dose - 1"
)
for (text in stayCode) {
  code <- str2lang(text)
  blocks <- lapply(names(codeDoors), function(door) {
    codeDoors[[door]](code, FALSE)$data@bases[[2L]]
  })
  width <- if (text %in% stayCode[1:4]) 2L else 1L
  for (index in seq_along(blocks)) {
    expect_identical(ncol(blocks[[index]]), width, info = text)
    expect_identical(unname(blocks[[index]]), unname(blocks[[1L]]), info = text)
  }
}
stayValues <- list(
  two,
  as.character(two),
  z == 1,
  arm == "b",
  dose,
  matrix(dose),
  as.vector(scale(age))
)
for (value in stayValues) {
  blocks <- lapply(valueDoors, function(door) {
    made <- door(value, FALSE)
    made$data@bases[[2L]]
  })
  for (index in seq_along(blocks)) {
    expect_identical(
      unname(blocks[[index]]),
      unname(blocks[[1L]]),
      info = names(valueDoors)[index]
    )
  }
  expect_identical(ncol(blocks[[1L]]), if (is.numeric(value)) 1L else 2L)
}

## --- after the rows --------------------------------------------------------
# a factor of three levels that 'subset' leaves with two is a basis of two
# levels, as code and as a value; under a 'subset' that keeps all three it is
# refused, at both
leavesTwo <- arm != "c"
keepsAll <- seq_len(n) <= 100L
# dbartsSpec() takes a data object whose rows are already cut
for (door in setdiff(names(codeDoors), "spec")) {
  made <- codeDoors[[door]](quote(arm), FALSE, subset = leavesTwo)
  expect_identical(made$data@basis.levels[[2L]], c("a", "b"), info = door)
  expect_identical(
    refusal(codeDoors[[door]](quote(arm), FALSE, subset = keepsAll)),
    refusedCode$arm,
    info = door
  )
}
for (door in c("list.formula", "list.matrix", "data")) {
  made <- valueDoors[[door]](arm, FALSE, subset = leavesTwo)
  expect_identical(made$data@basis.levels[[2L]], c("a", "b"), info = door)
  expect_identical(
    refusal(valueDoors[[door]](arm, FALSE, subset = keepsAll)),
    valueLevels,
    info = door
  )
}
# the rows are judged before the size: two numeric columns of which 'subset'
# empties one are refused for the empty column, as the zero column it is
expect_identical(
  refusal(dbartsData(
    x,
    y,
    subset = z == 1,
    bases = list(NULL, cbind(1 - z, z))
  )),
  "a 'bases' column of all zeros contributes nothing to a forest"
)

## --- the forms the messages show -------------------------------------------
perLevel <- dbarts(
  y ~ forest(x1 + x2) +
    forest(x1, basis = arm == "b") +
    forest(x1, basis = arm == "c"),
  d,
  control = sizeControl()
)
expect_identical(
  attr(perLevel$control, "bartcore.forests")$labels,
  c("forest1", "arm == \"b\"", "arm == \"c\"")
)
expect_identical(
  perLevel$data@basis.levels,
  list(NULL, c("FALSE", "TRUE"), c("FALSE", "TRUE"))
)
perColumn <- dbarts(
  y ~ forest(x1 + x2) + forest(x1, basis = dose) + forest(x1, basis = age),
  d,
  control = sizeControl()
)
expect_identical(length(perColumn$data@bases), 3L)
expect_identical(vapply(perColumn$data@bases[2:3], ncol, 0L), c(1L, 1L))

## --- after creation: no forest changes width ---------------------------------
numericSampler <- function() {
  dbarts(
    y ~ x1 + x2 + forest(x1, basis = dose),
    d,
    control = sizeControl()
  )
}
swapped <- numericSampler()
twin <- numericSampler()
widthText <- paste0(
  "forest 2 ('dose') multiplies one numeric column, and '%s' gives it 2; in ",
  "this version a basis is one numeric column. Give the one column"
)
expect_identical(
  refusal(swapped$setForestBasis(2L, B)),
  sprintf(widthText, "basis")
)
expect_identical(
  refusal(swapped$setForestBasis(2L, ~ dose + age)),
  sprintf(widthText, "basis")
)
# an array of three dimensions is as many columns as its second
expect_identical(
  refusal(swapped$setForestBasis(2L, array(c(dose, age), c(n, 2L, 1L)))),
  sprintf(widthText, "basis")
)
expect_identical(swapped$data, twin$data)
expect_identical(swapped$run(0L, 5L), twin$run(0L, 5L))
expect_identical(swapped$getForestAmplitudes(), twin$getForestAmplitudes())
swapped$setForestBasis(2L, cbind(age))
expect_identical(unname(swapped$data@bases[[2L]]), unname(cbind(age)))
expect_identical(colnames(swapped$data@bases[[2L]]), "dose")
fit <- bart(
  y ~ x1 + x2 + forest(x1, basis = dose),
  d,
  n.trees = 5L,
  n.samples = 3L,
  n.burn = 0L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE,
  keepTrees = TRUE
)
expect_identical(
  refusal(predict(fit, d, bases = list(NULL, B))),
  sprintf(widthText, "bases")
)
expect_identical(dim(predict(fit, d, bases = list(NULL, age))), c(3L, n))

## --- edited by hand, each at a named door ------------------------------------
numericData <- dbartsData(x, y, bases = list(NULL, dose))
proto <- dbartsSpec(numericData, sizeControl())
madeFrom <- list(
  new = function(data) {
    new("dbartsSampler", proto$control, proto$model, data)
  },
  dbarts = function(data) dbarts(data, control = sizeControl()),
  dbartsSpec = function(data) dbartsSpec(data, sizeControl()),
  validObject = function(data) validObject(data)
)
twoColumns <- numericData
twoColumns@bases[[2L]] <- unname(B)
threeLevels <- dbartsData(x, y, bases = list(NULL, two))
threeLevels@bases[[2L]] <- unname(gIndicators)
threeLevels@basis.levels[[2L]] <- c("a", "b", "c")
for (door in names(madeFrom)) {
  expect_identical(
    sub(
      "^invalid class .dbartsData. object: ",
      "",
      refusal(madeFrom[[door]](twoColumns))
    ),
    valueColumns(),
    info = door
  )
  expect_identical(
    sub(
      "^invalid class .dbartsData. object: ",
      "",
      refusal(madeFrom[[door]](threeLevels))
    ),
    valueLevels,
    info = door
  )
}
# the size is the record's and never the rows': a copy of a sampler whose
# level a swap emptied, and of one swapped to a column of zeros, is made and
# runs
levelsSampler <- dbarts(
  y ~ x1 + x2 + forest(x1, basis = factor(z)),
  d,
  control = sizeControl()
)
levelsSampler$setForestBasis(2L, factor(rep(1, n), levels = 0:1))
expect_identical(dim(levelsSampler$copy()$run(0L, 2L)$train), c(n, 2L))
zeroSampler <- numericSampler()
zeroSampler$setForestBasis(2L, rep(0, n))
expect_identical(dim(zeroSampler$copy()$run(0L, 2L)$train), c(n, 2L))

## --- the names, for later (dec-B272) ---------------------------------------
# a basis of several terms is read and named by the builder, on the model
# frame and with no sampler, as lm() names the columns of ~ 0 + <basis>; the
# refusal quotes those names
builtNames <- function(text) {
  read <- dbarts:::readForestBasis(
    dbarts:::basisCode(str2lang(text), environment()),
    d,
    n
  )
  colnames(dbarts:::buildCodeBasis(read)$basis)
}
lmNames <- function(text) {
  formula <- stats::as.formula(paste("y ~ 0 +", text))
  names(stats::coef(stats::lm(formula, d)))
}
for (text in c("dose + age", "poly(dose, 2)", "splines::ns(dose, 3)", "B")) {
  names <- builtNames(text)
  expect_identical(names, lmNames(text), info = text)
  expect_true(
    grepl(toString(sQuote(names, FALSE)), refusedCode[[text]], fixed = TRUE),
    info = text
  )
}
expect_identical(builtNames("1 + dose"), c("(Intercept)", "dose"))
expect_identical(builtNames("arm"), c("arma", "armb", "armc"))
