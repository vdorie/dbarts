# What a forest multiplies is decided by the class of what was handed over
# and by nothing else, at every door, and is fixed at creation: a factor, a
# character vector or a logical vector is levels, a number is numbers
# whatever its values. Every object below is taken through all seven creation
# doors and held to the FIRST door, so that no door can be right alone.
# `expected2` in the objects' table is what the next version of the rule
# makes of an object, one size a forest, and is not asserted yet.

forest <- dbartsForests$forest
fixed <- dbartsPriors$fixed

set.seed(20261011)
n <- 150L
x1 <- runif(n)
x2 <- runif(n)
z <- as.double(seq_len(n) %in% sample(n, 84L))
zf <- factor(z)
zc <- ifelse(z == 1, "treated", "control")
g <- sample(c("a", "b", "c"), n, replace = TRUE)
arm <- factor(g)
dose <- runif(n, 0.5, 3)
age <- rnorm(n, 50, 10)
y <- sin(3 * x1) + z * (0.5 + x2) + rnorm(n, sd = 0.3)
x <- cbind(x1 = x1, x2 = x2)
d <- data.frame(y, x1, x2)

kindControl <- function() {
  dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 15L,
    n.samples = 10L,
    n.burn = 0L,
    updateState = FALSE
  )
}
valueForests <- function(B) {
  list(forest(), do.call(forest, list(basis = B)))
}
# The seven creation doors. `B` is written as code at the first three and is
# handed over as a value at the others; `more` is 'subset' and 'weights'.
doors <- list(
  term.dbarts = function(B, y, more) {
    d$y <- y
    formula <- y ~ forest(x1 + x2) + forest(x1 + x2, basis = B)
    do.call(dbarts, c(list(formula, d, control = kindControl()), more))
  },
  term.bart = function(B, y, more) {
    d$y <- y
    formula <- y ~ forest(x1 + x2) + forest(x1 + x2, basis = B)
    fixedArgs <- list(
      formula,
      d,
      n.trees = 15L,
      n.samples = 10L,
      n.burn = 0L,
      n.chains = 1L,
      n.threads = 1L,
      verbose = FALSE,
      samplerOnly = TRUE
    )
    do.call(bart, c(fixedArgs, more))
  },
  list.code = function(B, y, more) {
    d$y <- y
    forests <- list(forest(), forest(basis = B))
    fixedArgs <- list(y ~ x1 + x2, d, forests = forests)
    do.call(dbarts, c(fixedArgs, list(control = kindControl()), more))
  },
  list.value = function(B, y, more) {
    d$y <- y
    fixedArgs <- list(y ~ x1 + x2, d, forests = valueForests(B))
    do.call(dbarts, c(fixedArgs, list(control = kindControl()), more))
  },
  matrix.value = function(B, y, more) {
    fixedArgs <- list(x, y, forests = valueForests(B))
    do.call(dbarts, c(fixedArgs, list(control = kindControl()), more))
  },
  # dbartsSpec() takes a data object, whose rows are already cut: the value
  # it is handed has those rows
  spec = function(B, y, more) {
    data <- do.call(dbartsData, c(list(x, y), more))
    kept <- if (is.null(more$subset)) !is.na(y) else more$subset & !is.na(y)
    cut <- if (is.null(dim(B))) B[kept] else B[kept, , drop = FALSE]
    spec <- dbartsSpec(data, kindControl(), forests = valueForests(cut))
    new("dbartsSampler", spec$control, spec$model, spec$data)
  },
  data = function(B, y, more) {
    data <- do.call(dbartsData, c(list(x, y, bases = list(NULL, B)), more))
    dbarts(data, control = kindControl())
  },
  data.formula = function(B, y, more) {
    d$y <- y
    fixedArgs <- list(y ~ x1 + x2, d, bases = list(NULL, B))
    dbarts(do.call(dbartsData, c(fixedArgs, more)), control = kindControl())
  }
)
codeDoors <- c("term.dbarts", "term.bart", "list.code")
dataDoors <- c("data", "data.formula")
sevenDoors <- setdiff(names(doors), "data.formula")
# the doors a row the fit drops can reach: dbartsSpec() has no such row
rowDoors <- setdiff(names(doors), "spec")

# A text as `door` words it: the data object says 'bases', its own argument,
# where a forest says 'basis'.
atDoor <- function(text, door) {
  if (door %in% dataDoors) {
    gsub("'basis'", "'bases'", text, fixed = TRUE)
  } else {
    text
  }
}
# The doors whose blocks are named alike: code names a column by its term,
# a value handed to a forest keeps its own names or has none, and the data
# door's names are its own.
doorFamily <- function(door) {
  if (door %in% codeDoors) {
    "code"
  } else if (door %in% dataDoors) {
    "data"
  } else {
    "value"
  }
}

# What a door makes of `B`: the forest's kind, levels and block, the block's
# column names and ten seeded draws, or the refusal's text as the door words
# it.
made <- function(door, B, response = y, ...) {
  more <- list(...)
  tryCatch(
    {
      set.seed(7)
      sampler <- doors[[door]](B, response, more)
      block <- sampler$data@bases[[2L]]
      levels <- sampler$data@basis.levels
      set.seed(7)
      list(
        sampler = sampler,
        kind = dbarts:::forestBasisKind(block, levels[[2L]]),
        levels = levels,
        # a value's block keeps what as.matrix() left on it, so the numbers
        block = matrix(as.vector(block), nrow(block)),
        names = colnames(block),
        attributes = names(attributes(block)),
        draws = sampler$run(0L, 10L)$train
      )
    },
    error = function(e) list(refused = conditionMessage(e))
  )
}
# A call that may be refused never stops this file: its refusal's text, or
# NULL where it is taken; and its value, or the refusal's text in its place.
refusalOf <- function(expr) {
  tryCatch(
    {
      expr
      NULL
    },
    error = function(e) conditionMessage(e)
  )
}
valueOf <- function(expr) {
  tryCatch(expr, error = function(e) conditionMessage(e))
}
# `made` is a sampler and not the text of a refusal, asserted
expectSampler <- function(made, info = NA_character_) {
  isSampler <- inherits(made, "dbartsSampler")
  expect_true(isSampler, info = info)
  isSampler
}
shapeOf <- function(result) {
  result[c("kind", "levels", "block", "draws")]
}
# Two doors' results are one creation: neither refused, and one kind, one
# set of levels, one block and one sum of draws.
expectSameShape <- function(other, first, info) {
  created <- is.null(first$refused) && is.null(other$refused)
  expect_true(created, info = info)
  if (created) {
    expect_identical(shapeOf(other), shapeOf(first), info = info)
  }
}
# `B` at each of `at`, every door held to the first: one kind, one block and
# one sum of draws, with one set of names at the doors of a family, or one
# refusal in each door's own word. Returns the first door's result.
expectOneReading <- function(B, at, info, ...) {
  first <- made(at[1L], B, ...)
  named <- list()
  named[[doorFamily(at[1L])]] <- list(first$names)
  for (door in at[-1L]) {
    other <- made(door, B, ...)
    where <- paste(info, "at", door)
    if (is.null(first$refused)) {
      expectSameShape(other, first, where)
      family <- doorFamily(door)
      if (is.null(named[[family]])) {
        named[[family]] <- list(other$names)
      }
      expect_identical(list(other$names), named[[family]], info = where)
      # the levels are kept beside the block and never on it
      bare <- !identical(other$kind, "levels") ||
        all(other$attributes %in% c("dim", "dimnames"))
      expect_true(bare, info = where)
    } else if (door %in% codeDoors == at[1L] %in% codeDoors) {
      expect_identical(other$refused, atDoor(first$refused, door), info = where)
    } else {
      expect_false(is.null(other$refused), info = where)
    }
  }
  first
}

# --- the objects, each at the seven doors ---
objects <- list(
  list("a factor", zf, "levels", 2L, "stays"),
  list("an ordered factor", factor(z, ordered = TRUE), "levels", 2L, "stays"),
  list("I() of a factor", I(zf), "levels", 2L, "stays"),
  list("a character vector", zc, "levels", 2L, "stays"),
  list("a logical vector", z == 1, "levels", 2L, "stays"),
  list("I() of a logical vector", I(z == 1), "levels", 2L, "stays"),
  list(
    "a factor with a level no row has",
    factor(z, 0:2),
    "levels",
    2L,
    "stays"
  ),
  list("addNA() of a complete factor", addNA(zf), "levels", 2L, "stays"),
  list("a 0/1 integer", as.integer(z), "numeric", 1L, "stays"),
  list("a 0/1 double", z, "numeric", 1L, "stays"),
  list("I() of a 0/1 double", I(z), "numeric", 1L, "stays"),
  list("a matrix of one 0/1 column", matrix(z), "numeric", 1L, "stays"),
  list("as.numeric() of a logical", as.numeric(z == 1), "numeric", 1L, "stays"),
  list("a time series", ts(z), "numeric", 1L, "stays"),
  list("-0.5 and 0.5", z - 0.5, "numeric", 1L, "stays"),
  list("1 and 2", z + 1, "numeric", 1L, "stays"),
  list("three values", as.double(as.integer(arm)), "numeric", 1L, "stays"),
  list("dose", dose, "numeric", 1L, "stays"),
  list("a constant column", rep(1, n), "numeric", 1L, "stays"),
  list("indicator columns", cbind(1 - z, z), "numeric", 2L, "refused"),
  list("two numeric columns", cbind(dose, age), "numeric", 2L, "refused"),
  list(
    "three indicator columns",
    model.matrix(~ g - 1),
    "numeric",
    3L,
    "refused"
  ),
  list("a factor of three levels", arm, "levels", 3L, "refused"),
  list("a character vector of three values", g, "levels", 3L, "refused")
)
for (object in objects) {
  info <- object[[1L]]
  first <- expectOneReading(object[[2L]], sevenDoors, info)
  expect_null(first$refused, info = info)
  expect_identical(first$kind, object[[3L]], info = info)
  expect_identical(ncol(first$block), object[[4L]], info = info)
  expect_identical(
    lengths(first$levels),
    c(0L, if (identical(first$kind, "levels")) object[[4L]] else 0L),
    info = info
  )
  expect_true(all(first$attributes %in% c("dim", "dimnames")), info = info)
}
# the levels, in the factor's own order, a character vector's sorted values,
# FALSE then TRUE; a forest with no basis and one of numbers record none
expect_identical(made("data", zf)$levels, list(NULL, c("0", "1")))
expect_identical(made("data", zc)$levels, list(NULL, c("control", "treated")))
expect_identical(made("data", z == 1)$levels, list(NULL, c("FALSE", "TRUE")))
expect_identical(
  made("data", factor(zc, levels = c("treated", "control")))$levels,
  list(NULL, c("treated", "control"))
)
expect_identical(made("data", dose)$levels, list(NULL, NULL))
# a basis written as code records the levels its control's record holds
for (door in codeDoors) {
  record <- attr(made(door, zc)$sampler$control, "bartcore.forests")
  expect_identical(
    record$basisTerms[[2L]]$levels,
    c("control", "treated"),
    info = door
  )
}

# the content decides nothing: a 0/1 number is one numeric column at every
# door and is not the factor's fit, and a factor's indicator columns handed
# over as numbers are two numeric columns
factorDraws <- made("term.dbarts", zf)$draws
for (door in sevenDoors) {
  for (number in list(as.integer(z), z, I(z), matrix(z))) {
    result <- made(door, number)
    expect_identical(result$kind, "numeric", info = door)
    expect_identical(result$block, matrix(z), info = door)
    expect_false(identical(result$draws, factorDraws), info = door)
  }
  expect_identical(made(door, z == 1)$draws, factorDraws, info = door)
}
expect_identical(made("list.value", cbind(1 - z, z))$kind, "numeric")
expect_identical(made("data", cbind(1 - z, z))$levels, list(NULL, NULL))

# --- the objects every door refuses ---
classText <- function(what, advice = "") {
  paste0(
    "'basis' gives forest 2 ",
    what,
    "; a basis is a factor, a character vector, a logical vector, or numbers",
    advice
  )
}
termText <- function(what, column = "") {
  paste0(
    "'basis' (B) has a term that is ",
    what,
    ": a basis is a factor, a character or logical vector, or numeric",
    column
  )
}
logicalAdvice <- paste0(
  ". Give its column as a logical vector for two levels, or as numbers for ",
  "one coefficient"
)
characterAdvice <- ". Give its column as a character vector"
factorColumn <- ". Give its column as a factor"
# a matrix of several columns is told to give one of them
ofSeveral <- function(advice) {
  sub("its column", "one of its columns", advice, fixed = TRUE)
}
# a factor with dimensions: a vector is levels, a matrix is numbers or nothing
factorMatrix <- function(columns) {
  value <- factor(rep(zc, columns))
  dim(value) <- c(n, columns)
  value
}
oneValue <- " on every row the fit keeps, so its other level has no observations and the forest would be multiplied by a constant"
naLevel <- " a factor with NA as a level, and a level of a basis has a name: give that level one"
asVector <- "; write its column as a vector, as B[, 1]"
normText <- "a 'basis' row's norm is not representable; rescale the basis to moderate values"
withNaLevel <- factor(replace(zc, 3L, NA), exclude = NULL)
refusedObjects <- list(
  list(
    cbind(z == 1, z == 0),
    termText("a logical matrix"),
    classText("a logical matrix", ofSeveral(logicalAdvice))
  ),
  list(
    matrix(z == 1),
    termText("a logical matrix", asVector),
    classText("a logical matrix", logicalAdvice)
  ),
  list(
    matrix(zc),
    termText("a character matrix", asVector),
    classText("a character matrix", characterAdvice)
  ),
  list(
    cbind(zc, zc),
    termText("a character matrix"),
    classText("a character matrix", ofSeveral(characterAdvice))
  ),
  list(
    factorMatrix(1L),
    termText("a factor matrix", asVector),
    classText("a factor matrix", factorColumn)
  ),
  list(
    factorMatrix(2L),
    termText("a factor matrix"),
    classText("a factor matrix", ofSeveral(factorColumn))
  ),
  list(
    data.frame(z = z),
    "'basis' (B): invalid type (list) for variable 'B'",
    classText("a data frame", ". Give one of its columns")
  ),
  list(
    as.list(z),
    "'basis' (B): invalid type (list) for variable 'B'",
    classText("a list", ". Give one of its columns")
  ),
  list(as.Date("2026-01-01") + z, termText("a Date"), classText("a Date")),
  list(complex(real = z), termText("a complex"), classText("a complex")),
  list(
    factor(rep("a", n)),
    "a 'basis' factor must have at least two levels",
    "a 'basis' factor must have at least two levels"
  ),
  list(
    replace(dose, 3L, NA),
    "a 'basis' cannot be NA",
    "a 'basis' cannot be NA"
  ),
  list(replace(zf, 3L, NA), "a 'basis' cannot be NA", "a 'basis' cannot be NA"),
  list(
    rep(0, n),
    "a 'basis' column of all zeros contributes nothing to a forest",
    "a 'basis' column of all zeros contributes nothing to a forest"
  ),
  list(
    rep(TRUE, n),
    paste0("'basis' (B) is TRUE", oneValue),
    paste0("'basis' is TRUE", oneValue)
  ),
  # a row whose norm underflows to zero
  list(replace(rep(0, n), 3L, 1e-200), normText, normText),
  list(
    withNaLevel,
    paste0("'basis' (B) is", naLevel),
    paste0("'basis' gives forest 2", naLevel)
  )
)
for (index in seq_along(refusedObjects)) {
  object <- refusedObjects[[index]]
  info <- paste("refused object", index)
  asCode <- expectOneReading(object[[1L]], codeDoors, info)
  expect_identical(asCode$refused, object[[2L]], info = info)
  asValue <- expectOneReading(
    object[[1L]],
    setdiff(names(doors), codeDoors),
    info
  )
  expect_identical(asValue$refused, object[[3L]], info = info)
}
# a matrix of no columns, at each door that takes a value
for (door in setdiff(names(doors), codeDoors)) {
  expect_identical(
    made(door, matrix(0, n, 0L))$refused,
    atDoor("a 'basis' must have at least one column", door),
    info = door
  )
}
# the data object says 'bases', its own argument, where a forest says 'basis'
expect_error(
  dbartsData(x, y, bases = list(NULL, matrix(z == 1))),
  "'bases' gives forest 2 a logical matrix",
  fixed = TRUE
)
expect_error(
  dbartsData(x, y, bases = list(NULL, rep(0, n))),
  "a 'bases' column of all zeros contributes nothing to a forest",
  fixed = TRUE
)
# a forest is named by its label where it has one of its own
expect_error(
  dbarts(
    x,
    y,
    forests = list(mu = forest(), tau = do.call(forest, list(basis = d))),
    control = kindControl()
  ),
  "'basis' gives forest 2 ('tau') a data frame",
  fixed = TRUE
)
# and where its rows are judged, in a list at dbarts() and at dbartsSpec()
labelled <- list(
  mu = forest(),
  tau = do.call(forest, list(basis = withNaLevel))
)
labelledText <- paste0("'basis' gives forest 2 ('tau')", naLevel)
expect_error(
  dbarts(x, y, forests = labelled, control = kindControl()),
  labelledText,
  fixed = TRUE
)
expect_error(
  dbartsSpec(dbartsData(x, y), kindControl(), forests = labelled),
  labelledText,
  fixed = TRUE
)
# a value's block is named whole or not at all: columns of which one has no
# name have none at a forest door
expect_null(made("list.value", cbind(1 - z, z))$names)

# --- the rows a missing value is judged on ---
# one level of g emptied, by 'subset' and by a missing response
keep <- g != "c"
yDropped <- replace(y, !keep, NA)
dropped <- which(!keep)[1L]
kept <- which(keep)[1L]
for (complete in list(dose, zf)) {
  onDropped <- replace(complete, dropped, NA)
  onKept <- replace(complete, kept, NA)
  for (how in c("subset", "a missing response")) {
    rows <- if (how == "subset") {
      list(subset = keep)
    } else {
      list(response = yDropped)
    }
    info <- paste(class(complete), how)
    # on a row the fit drops it is not looked at: the complete column's fit
    first <- do.call(expectOneReading, c(list(onDropped, rowDoors, info), rows))
    expect_null(first$refused, info = info)
    whole <- do.call(made, c(list("term.dbarts", complete), rows))
    expectSameShape(first, whole, info)
    # on a row the fit keeps it is refused, in one text
    for (door in rowDoors) {
      refused <- do.call(made, c(list(door, onKept), rows))$refused
      expected <- atDoor("a 'basis' cannot be NA", door)
      expect_identical(refused, expected, info = door)
    }
  }
}
# a value that is not finite goes as a missing one does
expect_null(made("data", replace(dose, dropped, Inf), subset = keep)$refused)
expect_identical(
  made("data", replace(dose, kept, Inf), subset = keep)$refused,
  "a 'bases' must be finite"
)
# dbartsSpec() over a data object cut by 'subset': the value has the kept
# rows, and one of the data's full length is refused for its length
subsetData <- dbartsData(x, y, subset = keep)
expect_error(
  dbartsSpec(subsetData, kindControl(), forests = valueForests(zf)),
  "'basis' must have the same length as 'y'",
  fixed = TRUE
)
expect_identical(
  dbartsSpec(
    subsetData,
    kindControl(),
    forests = valueForests(zf[keep])
  )$data@basis.levels,
  list(NULL, c("0", "1"))
)

# --- the rows at creation ---
# a level the kept rows leave empty is no column, in every spelling: the
# fit of the same factor written as code, whose levels are recorded so
for (how in c("subset", "a missing response")) {
  rows <- if (how == "subset") {
    list(subset = keep)
  } else {
    list(response = yDropped)
  }
  for (value in list(arm, g)) {
    info <- paste(class(value), how)
    first <- do.call(expectOneReading, c(list(value, sevenDoors, info), rows))
    expect_identical(first$levels, list(NULL, c("a", "b")), info = info)
    expect_identical(ncol(first$block), 2L, info = info)
  }
}
# fewer than two levels left, a logical left with one value and a number
# left all zero are refused at every door
twoLeftWithOne <- factor(ifelse(keep, "ab", "c"))
onlyOnDropped <- as.double(!keep) * dose
oneLevel <- expectOneReading(twoLeftWithOne, sevenDoors, "", subset = keep)
expect_identical(
  oneLevel$refused,
  "a 'basis' factor must have at least two levels"
)
oneLogical <- expectOneReading(keep, codeDoors, "", subset = keep)
expect_identical(oneLogical$refused, paste0("'basis' (B) is TRUE", oneValue))
for (door in setdiff(sevenDoors, codeDoors)) {
  expect_identical(
    made(door, keep, subset = keep)$refused,
    atDoor(paste0("'basis' is TRUE", oneValue), door),
    info = door
  )
  expect_identical(
    made(door, !keep, response = yDropped)$refused,
    atDoor(paste0("'basis' is FALSE", oneValue), door),
    info = door
  )
}
allZero <- expectOneReading(onlyOnDropped, sevenDoors, "", subset = keep)
expect_identical(
  allZero$refused,
  "a 'basis' column of all zeros contributes nothing to a forest"
)
expect_identical(
  made("data", onlyOnDropped, response = yDropped)$refused,
  atDoor(allZero$refused, "data")
)
# a row of weight zero is a row the fit keeps: its level stays
zeroWeights <- as.double(z == 1)
weighted <- expectOneReading(zf, sevenDoors, "", weights = zeroWeights)
expect_identical(weighted$levels, list(NULL, c("0", "1")))
expect_identical(weighted$block, unname(cbind(1 - z, z)))
# NA as a level: refused where a kept row is at it, at every door, and
# dropped as any level no row has where none is
atNaLevel <- is.na(as.character(withNaLevel))
named <- expectOneReading(withNaLevel, sevenDoors, "", subset = !atNaLevel)
expect_identical(named$levels, list(NULL, c("control", "treated")))

# --- the record on the data object ---
# dbartsSpec() over a data object that carries bases: a declaration's basis
# replaces the block and its levels, and a list that states none leaves both
numericData <- dbartsData(x, y, bases = list(NULL, dose))
levelsData <- dbartsData(x, y, bases = list(NULL, zc))
declared <- dbartsSpec(numericData, kindControl(), forests = valueForests(zc))
expect_identical(
  declared$data@basis.levels,
  list(NULL, c("control", "treated"))
)
expect_identical(declared$data@bases, levelsData@bases)
settingsOnly <- dbartsSpec(
  levelsData,
  kindControl(),
  forests = list(forest(n.trees = 12L), forest(n.trees = 8L))
)
expect_identical(settingsOnly$data@basis.levels, levelsData@basis.levels)
expect_identical(settingsOnly$data@bases, levelsData@bases)
renumbered <- dbartsSpec(
  levelsData,
  kindControl(),
  forests = valueForests(dose)
)
expect_identical(renumbered$data@basis.levels, list(NULL, NULL))
expect_identical(renumbered$data@bases, numericData@bases)

# a data object saved before the slot existed has none: it creates at every
# door that takes a data object, runs, and reads as numbers. The class's
# validity refuses it, in R's words, and no door calls that
proto <- dbartsSpec(levelsData, kindControl())
slotless <- proto$data
attr(slotless, "basis.levels") <- NULL
expect_error(validObject(slotless), "basis.levels", fixed = TRUE)
madeFrom <- list(
  new = function(data) {
    new("dbartsSampler", proto$control, proto$model, data)
  },
  dbarts = function(data) dbarts(data, control = kindControl()),
  dbartsSpec = function(data) {
    spec <- dbartsSpec(data, kindControl())
    new("dbartsSampler", spec$control, spec$model, spec$data)
  }
)
for (door in names(madeFrom)) {
  sampler <- valueOf(madeFrom[[door]](slotless))
  if (!expectSampler(sampler, door)) {
    next
  }
  expect_identical(dim(sampler$run(0L, 3L)$train), c(n, 3L), info = door)
  expect_null(dbarts:::dataBasisLevels(sampler$data), info = door)
  # as a forest of numbers it takes numbers
  swapped <- refusalOf(sampler$setForestBasis(2L, cbind(1 - z, z)))
  expect_null(swapped, info = door)
}
# recorded levels of another count than the block's columns are refused
miscounted <- proto$data
miscounted@basis.levels <- list(NULL, c("a", "b", "c"))
countText <- paste0(
  "'basis.levels' gives forest 2 3 levels and its block in 'bases' has 2 ",
  "columns; dbartsData() and dbarts() write the two together"
)
for (door in names(madeFrom)) {
  expect_error(madeFrom[[door]](miscounted), countText, fixed = TRUE)
}
expect_error(dbartsSpec(miscounted, kindControl()), countText, fixed = TRUE)
expect_error(validObject(miscounted), countText, fixed = TRUE)
# fewer levels than columns are refused as more are
tooFew <- proto$data
tooFew@basis.levels <- list(NULL, "a")
expect_error(
  dbarts(tooFew, control = kindControl()),
  sub("3 levels", "1 level", countText, fixed = TRUE),
  fixed = TRUE
)
# and levels recorded for a forest whose block was taken off by hand
stripped <- dbartsData(x, rev(y), bases = list(NULL, zc))
stripped@bases <- NULL
strippedText <- sub(
  "3 levels and its block in 'bases' has 2",
  "2 levels and its block in 'bases' has 0",
  countText,
  fixed = TRUE
)
expect_error(
  dbarts(stripped, control = kindControl()),
  strippedText,
  fixed = TRUE
)

# --- what a fit stores ---
fitArgs <- list(
  n.trees = 15L,
  n.samples = 20L,
  n.burn = 5L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE,
  keepTrees = TRUE,
  seed = 11L
)
# bart() takes no 'forests': a list's fit is packaged from its sampler
packaged <- function(forests, ...) {
  control <- dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 15L,
    n.samples = 20L,
    n.burn = 5L,
    keepTrees = TRUE,
    updateState = FALSE,
    seed = 11L
  )
  sampler <- dbarts(x, y, forests = forests, control = control, ...)
  run <- dbarts:::runWithBurnIn(sampler, sampler$control, TRUE)
  dbarts:::packageBartResults(
    sampler,
    run$samples,
    run$burnInSigma,
    run$burnInK,
    TRUE,
    TRUE
  )
}
fitDoors <- list(
  term = function(B, ...) {
    formula <- y ~ forest(x1 + x2) + forest(x1 + x2, basis = B)
    do.call(bart, c(list(formula, d, ...), fitArgs))
  },
  list = function(B, ...) packaged(valueForests(B), ...),
  data = function(B, ...) {
    data <- dbartsData(x, y, bases = list(NULL, B), ...)
    do.call(bart, c(list(data), fitArgs))
  }
)
reversed <- factor(zc, levels = c("treated", "control"))
stored <- list(
  list(zf, factor(z)),
  list(zc, factor(zc)),
  list(z == 1, factor(z == 1)),
  # the levels recorded, in their own order, and never the sorted ones
  list(reversed, reversed)
)
for (door in names(fitDoors)) {
  for (entry in stored) {
    fit <- fitDoors[[door]](entry[[1L]])
    expect_identical(fit$bases, list(NULL, entry[[2L]]), info = door)
    # the block of an entry is the kept sampler's own
    block <- dbarts:::storedBasisBlock(fit$bases[[2L]])
    expect_identical(block, unname(fit$fit$data@bases[[2L]]), info = door)
  }
  # the levels the kept rows leave
  fit <- fitDoors[[door]](arm, subset = keep)
  expect_identical(fit$bases[[2L]], factor(g[keep]), info = door)
  # a forest of numbers stores its block, as it did
  fit <- fitDoors[[door]](dose)
  expect_identical(unname(fit$bases[[2L]]), matrix(dose), info = door)
  expect_identical(dbarts:::storedBasisBlock(fit$bases[[2L]]), fit$bases[[2L]])
}
# a forest's contribution is the product with the kept sampler's block
fit <- fitDoors$term(reversed)
block <- unname(fit$fit$data@bases[[2L]])
glue <- fit$glue[, 2:3]
contribution <- valueOf(
  unname(extract(fit, type = "forest", contribution = TRUE)[,, 2L])
)
expect_identical(dim(contribution), c(20L, n))
expect_identical(
  contribution,
  valueOf(unname((glue %*% t(block)) * extract(fit, type = "forest")[,, 2L]))
)

# --- after creation: a forest takes the kind it was created with ---
z2 <- rev(z)
zc2 <- rev(zc)
newRows <- data.frame(x1 = rev(x1), x2 = rev(x2))
forestKinds <- list(
  factor = zf,
  character = zc,
  logical = z == 1,
  number = z,
  dose = dose
)
ownValue <- list(
  factor = factor(z2),
  character = zc2,
  logical = z2 == 1,
  number = z2,
  dose = rev(dose)
)
levelsText <- function(levels) {
  paste0(
    "forest 2 multiplies a factor with the levels ",
    levels,
    ", and '%s' gives it numbers. Give a factor, a character vector or a ",
    "logical vector with those levels; if the numbers are a factor's ",
    "indicator columns, give the factor itself. To multiply numbers, create ",
    "the model again with them"
  )
}
numbersText <- function(what, advice) {
  paste0("forest 2 multiplies numbers, and '%s' gives it ", what, ". ", advice)
}
unknownText <- function(level, levels) {
  paste0(
    "'%s' gives forest 2 the level '",
    level,
    "', which it does not have: its levels are ",
    levels,
    ", fixed at creation. Give every row one of them"
  )
}
factorAdvice <- paste0(
  "Give the column the forest was created with, at these rows; to multiply ",
  "a factor's levels, create the model again with the factor"
)
logicalToNumbers <- paste0(
  "Give numbers: as.numeric() of a logical vector is its 0/1 column. To ",
  "multiply two levels, create the model again with the vector"
)
noneText <- paste0(
  "'%s' gives forest 1 a basis, and forest 1 was created without one: what ",
  "a forest multiplies is fixed at creation. Create the model again, with ",
  "the basis, to give it one"
)
quotedLevels <- c(
  factor = "'0', '1'",
  character = "'control', 'treated'",
  logical = "'FALSE', 'TRUE'"
)
firstLevel <- c(factor = "0", character = "control", logical = "FALSE")
# what a forest of `kind` makes of the value named `value`: NULL where it is
# taken, and otherwise the refusal with '%s' where the argument is named
expectedText <- function(kind, value) {
  ofLevels <- kind %in% names(quotedLevels)
  if (value %in% c("logical matrix", "factor matrix", "data frame")) {
    advice <- switch(
      value,
      "data frame" = ". Give one of its columns",
      "factor matrix" = factorColumn,
      logicalAdvice
    )
    sub("'basis'", "'%s'", classText(paste("a", value), advice), fixed = TRUE)
  } else if (value %in% c("number", "dose", "zeros", "indicators", "~ dose")) {
    if (ofLevels) levelsText(quotedLevels[[kind]])
  } else {
    given <- sub("~ ", "", value, fixed = TRUE)
    if (ofLevels && given != kind) {
      unknownText(firstLevel[[given]], quotedLevels[[kind]])
    } else if (given == "logical" && !ofLevels) {
      numbersText("a logical vector", logicalToNumbers)
    } else if (!ofLevels) {
      what <- if (given == "factor") "a factor" else "a character vector"
      numbersText(what, factorAdvice)
    }
  }
}
swapValues <- c(
  ownValue,
  list(
    zeros = rep(0, n),
    indicators = cbind(1 - z2, z2),
    "logical matrix" = matrix(z2 == 1),
    "factor matrix" = factorMatrix(1L),
    "data frame" = data.frame(z2),
    "~ dose" = ~dose,
    "~ factor" = ~ factor(z2),
    "~ character" = ~zc2,
    "~ logical" = ~ z2 == 1
  )
)
swapSampler <- function(B) {
  set.seed(7)
  sampler <- dbarts(x, y, forests = valueForests(B), control = kindControl())
  sampler$run(0L, 3L)
  sampler
}
for (kind in names(forestKinds)) {
  sampler <- swapSampler(forestKinds[[kind]])
  twin <- swapSampler(forestKinds[[kind]])
  taken <- character()
  for (value in names(swapValues)) {
    info <- paste(kind, "forest given", value)
    expected <- expectedText(kind, value)
    refused <- refusalOf(sampler$setForestBasis(2L, swapValues[[value]]))
    if (is.null(expected)) {
      expect_null(refused, info = info)
    } else {
      expect_identical(refused, sprintf(expected, "basis"), info = info)
    }
    if (is.null(refused)) {
      taken <- c(taken, value)
      # the twin takes it too, so that the two stay one sampler
      twin$setForestBasis(2L, swapValues[[value]])
    }
    # a value's block is named whole or not at all, at a swap as at creation
    if (value == "indicators" && is.null(refused)) {
      expect_null(colnames(sampler$data@bases[[2L]]), info = info)
    }
    # a refused swap has written nothing
    expect_identical(sampler$data, twin$data, info = info)
    expect_identical(
      sampler$getForestAmplitudes(),
      twin$getForestAmplitudes(),
      info = info
    )
  }
  # the forest with no basis takes none, whatever the other was created with
  expect_identical(
    refusalOf(sampler$setForestBasis(1L, dose)),
    sprintf(noneText, "basis"),
    info = kind
  )
  # a value of another length is refused in these words, before the engine
  expect_identical(
    refusalOf(sampler$setForestBasis(2L, ownValue[[kind]][-1L])),
    "'basis' must have the same length as 'y'",
    info = kind
  )
  expect_identical(sampler$data, twin$data, info = kind)
  expect_identical(sampler$run(0L, 5L), twin$run(0L, 5L), info = kind)
  # every forest takes its own kind, as a value and as a formula
  own <- c(
    factor = "factor ~ factor",
    character = "character ~ character",
    logical = "logical ~ logical",
    number = "number dose zeros indicators ~ dose",
    dose = "number dose zeros indicators ~ dose"
  )[[kind]]
  expect_identical(paste(taken, collapse = " "), own, info = kind)
}

# the same table at predict, a list's entry and the bare value alike
predictValues <- c(
  ownValue,
  list(
    zeros = rep(0, n),
    "logical matrix" = matrix(z2 == 1),
    "factor matrix" = factorMatrix(1L),
    "data frame" = data.frame(z2)
  )
)
fits <- lapply(forestKinds, fitDoors$list)
for (kind in names(forestKinds)) {
  fit <- fits[[kind]]
  for (value in names(predictValues)) {
    info <- paste(kind, "forest given", value, "at predict")
    expected <- expectedText(kind, value)
    given <- predictValues[[value]]
    inList <- refusalOf(predict(fit, newRows, bases = list(NULL, given)))
    bare <- refusalOf(predict(fit, newRows, bases = given))
    if (is.null(expected)) {
      expect_null(inList, info = info)
    } else {
      expect_identical(inList, sprintf(expected, "bases"), info = info)
    }
    expect_identical(bare, inList, info = info)
  }
  expect_identical(
    refusalOf(predict(fit, newRows, bases = list(dose, ownValue[[kind]]))),
    sprintf(noneText, "bases"),
    info = kind
  )
}
# a column of zeros predicts the fit with that forest left out, a
# prediction being linear in a numeric forest's column
atOne <- predict(fits$dose, newRows, bases = rep(1, n), type = "ev")
atTwo <- predict(fits$dose, newRows, bases = rep(2, n), type = "ev")
atZero <- valueOf(predict(fits$dose, newRows, bases = rep(0, n), type = "ev"))
expect_equal(atZero, 2 * atOne - atTwo, tolerance = 1e-12)
expect_false(isTRUE(all.equal(atZero, atOne)))

# --- which model, by name ---
# a held forest's term stays on the rows of its own level when the same
# vector arrives with its levels in another order
heldForests <- list(
  forest(),
  do.call(forest, list(basis = zc, amplitude = fixed()))
)
held <- dbarts(x, y, forests = heldForests, control = kindControl())
creationBlock <- held$data@bases[[2L]]
expect_identical(creationBlock, unname(cbind(1 - z, z)))
expect_null(refusalOf(held$setForestBasis(2L, reversed)))
expect_identical(held$data@bases[[2L]], creationBlock)
expect_identical(held$data@basis.levels, list(NULL, c("control", "treated")))

drawn <- swapSampler(zc)
before <- drawn$getForestAmplitudes()
# a vector left with one value: the forest keeps its two columns, one of
# them zero, and its coefficients, and has them back with the vector
expect_null(refusalOf(drawn$setForestBasis(2L, rep("treated", n))))
expect_identical(drawn$data@bases[[2L]], cbind(rep(0, n), rep(1, n)))
expect_identical(drawn$getForestAmplitudes(), before)
expect_null(refusalOf(drawn$setForestBasis(2L, zc)))
expect_identical(drawn$data@bases[[2L]], creationBlock)
expect_identical(drawn$getForestAmplitudes(), before)
# a row at a level the forest lacks is refused; the same level declared and
# on no row is passed by, and the recorded levels are creation's afterwards
expect_error(
  drawn$setForestBasis(2L, replace(zc, 5L, "other")),
  sprintf(unknownText("other", quotedLevels[["character"]]), "basis"),
  fixed = TRUE
)
declaredMore <- factor(zc2, levels = c("other", "treated", "control"))
expect_null(refusalOf(drawn$setForestBasis(2L, declaredMore)))
expect_identical(drawn$data@bases[[2L]], creationBlock[n:1L, ])
expect_identical(drawn$data@basis.levels, list(NULL, c("control", "treated")))
# a name is the forest's level character for character
nearLevels <- ifelse(z == 1, "treated", "a")
near <- swapSampler(nearLevels)
for (name in c("A", "a ", "tr")) {
  renamed <- replace(
    nearLevels,
    nearLevels == tolower(substr(name, 1L, 1L)),
    name
  )
  renamed <- sub("^treated$", if (name == "tr") name else "treated", renamed)
  expect_error(
    near$setForestBasis(2L, renamed),
    sprintf(unknownText(name, "'a', 'treated'"), "basis"),
    fixed = TRUE,
    info = name
  )
}
# nor is a name read as a number: "1.0" is not the level "1"
numbered <- swapSampler(zf)
expect_error(
  numbered$setForestBasis(2L, ifelse(z2 == 1, "1.0", "0.0")),
  sprintf(unknownText("0.0", quotedLevels[["factor"]]), "basis"),
  fixed = TRUE
)
# a row at an NA level is refused, and one declared and on no row passed by
expect_error(
  drawn$setForestBasis(2L, withNaLevel),
  paste0("'basis' gives forest 2", naLevel),
  fixed = TRUE
)
expect_null(refusalOf(drawn$setForestBasis(2L, addNA(factor(zc)))))
# the factor a forest was created from, an unused level among its levels,
# is taken again and gives the block of creation
withUnused <- factor(zc, levels = c("control", "never", "treated"))
fromUnused <- swapSampler(withUnused)
expect_identical(fromUnused$data@bases[[2L]], creationBlock)
expect_null(refusalOf(fromUnused$setForestBasis(2L, rev(withUnused))))
expect_null(refusalOf(fromUnused$setForestBasis(2L, withUnused)))
expect_identical(fromUnused$data@bases[[2L]], creationBlock)
# the sampler's own stored block is numbers, which a forest of levels
# refuses: a caller keeps the factor it handed over
expect_error(
  drawn$setForestBasis(2L, drawn$data@bases[[2L]]),
  sprintf(levelsText(quotedLevels[["character"]]), "basis"),
  fixed = TRUE
)
# the kind is not a position: with the forest with no basis second, it
# refuses a basis and the first takes one
basisFirst <- dbarts(
  x,
  y,
  forests = list(do.call(forest, list(basis = zc)), forest()),
  control = kindControl()
)
expect_identical(
  basisFirst$data@basis.levels,
  list(c("control", "treated"), NULL)
)
expect_error(
  basisFirst$setForestBasis(2L, zc2),
  gsub("forest 1", "forest 2", sprintf(noneText, "basis"), fixed = TRUE),
  fixed = TRUE
)
expect_null(refusalOf(basisFirst$setForestBasis(1L, reversed)))
expect_identical(basisFirst$data@bases[[1L]], creationBlock)

# --- predict, by name ---
trainRows <- data.frame(x1, x2)
fromTerm <- fitDoors$term(zf)
fromData <- fitDoors$data(zf)
predicted <- function(fit, rows, bases) {
  valueOf(unname(predict(fit, rows, bases = bases)))
}
aligned <- predicted(fromTerm, newRows, factor(z2))
expect_identical(dim(aligned), c(20L, n))
# another order of levels is the same prediction
flipped <- factor(z2, levels = c("1", "0"))
expect_identical(predicted(fromTerm, newRows, flipped), aligned)
expect_identical(predicted(fromTerm, newRows, as.character(z2)), aligned)
# a fit over a data object from the data door predicts by name as well
expect_identical(
  predicted(fromData, newRows, flipped),
  predicted(fromData, newRows, factor(z2))
)
expect_identical(dim(predicted(fromData, newRows, flipped)), c(20L, n))
# every row at one level
allTreated <- predicted(fromTerm, newRows, rep("1", n))
expect_identical(
  predicted(fromTerm, newRows, factor(rep(1, n), levels = 0:1)),
  allTreated
)
expect_identical(dim(allTreated), c(20L, n))
expect_false(identical(allTreated, aligned))
expect_error(
  predict(fromTerm, newRows, bases = ifelse(z2 == 1, "one", "0")),
  "'bases' gives forest 2 ('B') the level 'one', which it does not have: its levels are '0', '1', fixed at creation",
  fixed = TRUE
)
expect_error(
  predict(fromTerm, newRows, bases = cbind(1 - z2, z2)),
  "forest 2 ('B') multiplies a factor with the levels '0', '1', and 'bases' gives it numbers",
  fixed = TRUE
)
expect_error(
  predict(fits$number, newRows, bases = factor(z2)),
  sprintf(numbersText("a factor", factorAdvice), "bases"),
  fixed = TRUE
)
# the training rows, on a factor whose levels are not in sorted order: the
# fit's own entry handed back is the stored training fit, to the blend's
# other order of summation, and so are the creation factor and the same
# vector as a character
unsorted <- fitDoors$term(reversed)
expect_identical(levels(unsorted$bases[[2L]]), c("treated", "control"))
for (handed in list(unsorted$bases, reversed, zc)) {
  expect_equal(
    predicted(unsorted, trainRows, handed),
    unname(unsorted$yhat.train),
    tolerance = 1e-12
  )
}
# with no row of 'newdata' left the fit's own entry stands in for the rows,
# a factor among them: draws by no rows
noRows <- replace(newRows, "x1", NA_real_)
expect_identical(
  dim(valueOf(predict(fromTerm, noRows, na.action = na.omit))),
  c(20L, 0L)
)
# the kind is not a position at predict either: with the forest with no
# basis second, it refuses a basis and the first takes one by name
noneSecond <- packaged(list(do.call(forest, list(basis = zc)), forest()))
expect_identical(
  refusalOf(predict(noneSecond, newRows, bases = list(zc2, dose))),
  gsub("forest 1", "forest 2", sprintf(noneText, "bases"), fixed = TRUE)
)
expect_identical(
  predicted(noneSecond, newRows, list(rev(reversed), NULL)),
  predicted(noneSecond, newRows, zc2)
)

# --- two forests of levels: each is held to its own ---
half <- ifelse(x2 > 0.5, "upper", "lower")
half2 <- rev(half)
upperFirst <- factor(half2, levels = c("upper", "lower"))
twoOfLevels <- list(
  forest(),
  do.call(forest, list(basis = zc)),
  do.call(forest, list(basis = half))
)
three <- dbarts(x, y, forests = twoOfLevels, control = kindControl())
expect_identical(
  three$data@basis.levels,
  list(NULL, c("control", "treated"), c("lower", "upper"))
)
expect_null(refusalOf(three$setForestBasis(3L, upperFirst)))
expect_identical(
  three$data@bases[[3L]],
  cbind(half2 == "lower", half2 == "upper") + 0
)
threeFit <- packaged(twoOfLevels)
expect_identical(threeFit$bases, list(NULL, factor(zc), factor(half)))
byName <- predicted(threeFit, newRows, list(NULL, zc2, upperFirst))
expect_identical(dim(byName), c(20L, n))
expect_identical(byName, predicted(threeFit, newRows, list(NULL, zc2, half2)))

# --- the same sampler again ---
# a copy, a reload and a sampler made from the sampler's own three have its
# kind and levels: each matches by name, and refuses what the original does
again <- function(sampler) {
  file <- tempfile(fileext = ".rds")
  sampler$storeState()
  saveRDS(sampler, file)
  reloaded <- readRDS(file)
  unlink(file)
  list(
    copy = valueOf(sampler$copy()),
    reload = valueOf({
      reloaded$run(0L, 1L)
      reloaded
    }),
    new = valueOf(
      new("dbartsSampler", sampler$control, sampler$model, sampler$data)
    )
  )
}
original <- swapSampler(zc)
for (how in c("copy", "reload", "new")) {
  sampler <- again(original)[[how]]
  if (!expectSampler(sampler, how)) {
    next
  }
  expect_identical(sampler$data@basis.levels, original$data@basis.levels)
  flippedRows <- factor(zc2, levels = c("treated", "control"))
  expect_null(refusalOf(sampler$setForestBasis(2L, flippedRows)), info = how)
  expect_identical(sampler$data@bases[[2L]], creationBlock[n:1L, ], info = how)
  expect_error(
    sampler$setForestBasis(2L, cbind(1 - z, z)),
    sprintf(levelsText(quotedLevels[["character"]]), "basis"),
    fixed = TRUE,
    info = how
  )
}
# a level, or a numeric column, that a swap left with no row is in the
# sampler's data for good: the rows are judged at a creation and nowhere else
emptied <- swapSampler(zc)
expect_null(refusalOf(emptied$setForestBasis(2L, rep("treated", n))))
zeroed <- swapSampler(dose)
expect_null(refusalOf(zeroed$setForestBasis(2L, rep(0, n))))
for (how in c("copy", "reload", "new")) {
  for (swapped in list(emptied, zeroed)) {
    sampler <- again(swapped)[[how]]
    if (!expectSampler(sampler, how)) {
      next
    }
    expect_identical(sampler$data@bases, swapped$data@bases, info = how)
    expect_identical(dim(sampler$run(0L, 5L)$train), c(n, 5L), info = how)
  }
}
# a state carries no data: installed in a sampler whose factor has the same
# levels in another order it goes by position, and the recipient's levels
# and block are its own before and after
donor <- swapSampler(factor(zc, levels = c("control", "treated")))
recipient <- swapSampler(reversed)
donor$run(0L, 5L)
donor$storeState()
ownBlock <- recipient$data@bases[[2L]]
expect_identical(ownBlock, unname(cbind(z, 1 - z)))
expect_true(valueOf(recipient$setState(donor$state)))
expect_identical(
  recipient$getForestAmplitudes(),
  donor$getForestAmplitudes()
)
expect_identical(recipient$data@bases[[2L]], ownBlock)
expect_identical(
  recipient$data@basis.levels,
  list(NULL, c("treated", "control"))
)

# --- $setData on a sampler of one forest ---
plain <- function() {
  set.seed(7)
  dbarts(x, y, control = kindControl())
}
oneForest <- plain()
untouched <- plain()
expect_error(
  oneForest$setData(levelsData),
  "$setData: the data carry forest bases, and this sampler was created with one forest and none; a sampler's forests are fixed at creation. Make a new sampler from the data with dbarts(), or hand over data without 'bases'",
  fixed = TRUE
)
expect_identical(oneForest$data, untouched$data)
# data whose bases were taken off by hand and whose levels stayed are
# refused as a sampler made from them is, before anything is stored: the
# sampler's data, its levels, its state and its next draws are a twin's
expect_error(oneForest$setData(stripped), strippedText, fixed = TRUE)
expect_identical(oneForest$data, untouched$data)
expect_null(oneForest$data@basis.levels)
oneForest$storeState()
untouched$storeState()
expect_identical(oneForest$state, untouched$state)
set.seed(3)
copied <- valueOf(oneForest$copy())
set.seed(3)
twinCopied <- untouched$copy()
if (expectSampler(copied)) {
  expect_identical(copied$data, untouched$data)
  expect_identical(copied$run(0L, 5L), twinCopied$run(0L, 5L))
}
expect_identical(oneForest$run(0L, 5L), untouched$run(0L, 5L))
# the same data without its bases is taken
expect_null(refusalOf(oneForest$setData(dbartsData(x, rev(y)))))
