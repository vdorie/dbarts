# The forest constructors (dec-A67): interactions, blocks, monotone, forest,
# varianceForest and fixed are not exported. Inside the arguments that take
# them they resolve by bare name, a bare name the caller binds taking the
# caller's value; outside, dbartsForests is their exported face.

source(
  system.file("common", "friedmanData.R", package = "dbarts"),
  local = TRUE
)

x <- testData$x
y <- testData$y
z <- rep_len(c(0, 1), nrow(x))
obj <- dbartsForests$interactions(max.order = 1)
hint <- "outside the argument that takes it, write dbartsForests\\$interactions"

maxOrder <- function(sampler) attr(sampler$model, "interaction.max.order")
blockOf <- function(sampler) attr(sampler$model, "block.of.column")
forestsOf <- function(sampler) attr(sampler$control, "bartcore.forests")
varianceOf <- function(sampler) attr(sampler$control, "bartcore.variance")

# --- the list ---------------------------------------------------------------

expect_true(is.list(dbartsForests))
expect_true(all(vapply(dbartsForests, is.function, logical(1L))))
expect_equal(
  sort(names(dbartsForests)),
  sort(c(
    "interactions",
    "blocks",
    "monotone",
    "forest",
    "varianceForest",
    "fixed"
  ))
)
for (name in names(dbartsForests)) {
  expect_false(name %in% getNamespaceExports("dbarts"))
}
expect_error(dbarts::interactions, "not an exported object")

# --- (1) the bare name equals the list spelling at each door -----------------

groups <- list(1:5, 6:10)
basisForests <- list(dbartsForests$forest(), dbartsForests$forest(basis = ~z))
bare <- dbarts::dbarts(
  x,
  y,
  interactions = interactions(max.order = 1),
  blocks = blocks(groups = groups),
  variance = varianceForest(n.trees = 10L)
)
listed <- dbarts::dbarts(
  x,
  y,
  interactions = dbartsForests$interactions(max.order = 1),
  blocks = dbartsForests$blocks(groups = groups),
  variance = dbartsForests$varianceForest(n.trees = 10L)
)
expect_identical(maxOrder(bare), maxOrder(listed))
expect_identical(blockOf(bare), blockOf(listed))
expect_identical(varianceOf(bare), varianceOf(listed))
expect_identical(
  forestsOf(dbarts::dbarts(x, y, forests = list(forest(), forest(basis = ~z)))),
  forestsOf(dbarts::dbarts(x, y, forests = basisForests))
)

data <- dbarts::dbartsData(x, y)
specOf <- function(spec) {
  list(attributes(spec$model), attributes(spec$control))
}
expect_identical(
  specOf(dbarts::dbartsSpec(
    data,
    interactions = interactions(max.order = 1),
    blocks = blocks(groups = groups),
    variance = varianceForest(n.trees = 10L)
  )),
  specOf(dbarts::dbartsSpec(
    data,
    interactions = dbartsForests$interactions(max.order = 1),
    blocks = dbartsForests$blocks(groups = groups),
    variance = dbartsForests$varianceForest(n.trees = 10L)
  ))
)
expect_identical(
  forestsOf(dbarts::dbartsSpec(
    data,
    forests = list(forest(), forest(basis = ~z))
  )),
  forestsOf(dbarts::dbartsSpec(data, forests = basisForests))
)

bartArgs <- list(
  n.samples = 5L,
  n.burn = 5L,
  n.chains = 1L,
  n.trees = 5L,
  seed = 1L,
  verbose = FALSE
)
bartFit <- function(...) {
  do.call(dbarts::bart, c(list(x, y), bartArgs, list(...)))
}
expect_identical(
  dbarts::bart(
    x,
    y,
    interactions = interactions(max.order = 1),
    blocks = blocks(groups = groups),
    variance = varianceForest(n.trees = 10L),
    n.samples = 5L,
    n.burn = 5L,
    n.chains = 1L,
    n.trees = 5L,
    seed = 1L,
    verbose = FALSE
  )$yhat.train,
  bartFit(
    interactions = dbartsForests$interactions(max.order = 1),
    blocks = dbartsForests$blocks(groups = groups),
    variance = dbartsForests$varianceForest(n.trees = 10L)
  )$yhat.train
)
expect_error(
  dbarts::bart(
    x,
    factor(rep_len(c("a", "b", "c"), nrow(x))),
    family = "multinomial",
    variance = varianceForest(n.trees = 10L)
  ),
  "does not support 'variance'"
)

# --- a name the caller also binds --------------------------------------------

# (2) as it is with an attached mask of every constructor name
masked <- function(...) stop("masked")
attach(
  list(interactions = masked, blocks = masked, forest = masked),
  name = "vocabularyMask"
)
expect_equal(
  maxOrder(dbarts::dbarts(x, y, interactions = interactions(max.order = 1))),
  1L
)
expect_identical(
  forestsOf(dbarts::dbarts(
    x,
    y,
    forests = list(
      forest(),
      forest(basis = ~z, blocks = blocks(groups = groups))
    )
  )),
  forestsOf(dbarts::dbarts(
    x,
    y,
    forests = list(
      dbartsForests$forest(),
      dbartsForests$forest(
        basis = ~z,
        blocks = dbartsForests$blocks(groups = groups)
      )
    )
  ))
)
detach("vocabularyMask")

# (3) a bare name the caller binds is the caller's
interactions <- obj
expect_equal(maxOrder(dbarts::dbarts(x, y, interactions = interactions)), 1L)

# (14) at any depth
nested <- forestsOf(dbarts::dbarts(
  x,
  y,
  forests = list(forest(), forest(basis = ~z, interactions = interactions))
))
expect_equal(nested$interactions[[2L]]$max.order, 1L)

# (21) and in a formula term's knobs, which resolve in the formula's environment
df <- data.frame(y = y, x1 = x[, 1L], x2 = x[, 2L], z = z)
expect_identical(
  forestsOf(dbarts::dbarts(
    y ~ x1 + x2 + forest(x1, basis = ~z, interactions = interactions),
    df
  )),
  forestsOf(dbarts::dbarts(
    y ~ x1 +
      x2 +
      forest(x1, basis = ~z, interactions = interactions(max.order = 1)),
    df
  ))
)

# (22) a wrapper whose top environment is a namespace never sees a user's
# global, so its bare name is the constructor, called for its defaults
namespaced <- new.env(parent = asNamespace("stats"))
namespaced$x <- x
namespaced$y <- y
nsWrapper <- function() dbarts::dbarts(x, y, interactions = interactions)
environment(nsWrapper) <- namespaced
expect_error(nsWrapper(), "needs at least one")

# (23) nor does the dbarts namespace's own binding count as the caller's
inDbarts <- new.env(parent = asNamespace("dbarts"))
inDbarts$x <- x
inDbarts$y <- y
dbartsWrapper <- function() dbarts::dbarts(x, y, interactions = interactions)
environment(dbartsWrapper) <- inDbarts
expect_error(dbartsWrapper(), "needs at least one")

# (11) NULL is a value, the door default
interactions <- NULL
expect_null(maxOrder(dbarts::dbarts(x, y, interactions = interactions)))
rm(interactions)

# (15) unbound, a bare name nested in 'forests' is the constructor, refused
expect_error(
  dbarts::dbarts(
    x,
    y,
    forests = list(forest(), forest(basis = ~z, interactions = interactions))
  ),
  "see \\?dbartsForests"
)

# (12, 13) a caller's function is no value, so the bare name is the
# constructor, called for its defaults
blocks <- function(...) "unrelated"
expect_error(dbarts::dbarts(x, y, blocks = blocks), "requires 'groups'")
rm(blocks)
varianceForest <- function(...) "unrelated"
expect_identical(
  varianceOf(dbarts::dbarts(x, y, variance = varianceForest)),
  varianceOf(dbarts::dbarts(x, y, variance = dbartsForests$varianceForest()))
)
rm(varianceForest)

# (16) a caller's list named 'forest'
forest <- basisForests
expect_identical(
  forestsOf(dbarts::dbarts(x, y, forests = forest)),
  forestsOf(dbarts::dbarts(x, y, forests = basisForests))
)
rm(forest)

# --- wrappers -----------------------------------------------------------------

# (4, 5) a wrapper's formal is forced where the wrapper was called, so an
# object reaches the door and a constructor call does not
passOn <- function(interactions) {
  dbarts::dbarts(x, y, interactions = interactions)
}
expect_equal(maxOrder(passOn(obj)), 1L)
expect_error(passOn(interactions(max.order = 1)), hint)

passOnMonotone <- function(monotone) {
  attr(dbarts::dbarts(x, y, monotone = monotone)$model, "monotone")
}
directions <- c(1, rep(0, ncol(x) - 1L))
expect_equal(
  passOnMonotone(dbartsForests$monotone(directions)),
  as.integer(directions)
)
expect_error(
  passOnMonotone(monotone(directions)),
  "outside the argument that takes it, write dbartsForests\\$monotone"
)

# (6, 7) an unsupplied formal is the door default, whatever its default
passOnNull <- function(interactions = NULL) {
  dbarts::dbarts(x, y, interactions = interactions)
}
expect_null(maxOrder(passOnNull()))
expect_null(maxOrder(passOn()))

# (8) a formal's default evaluates in the wrapper's frame
passOnDefault <- function(k, ints = interactions(max.order = k)) {
  dbarts::dbarts(x, y, interactions = ints)
}
expect_error(passOnDefault(1), hint)

# (9, 10) through lapply, a named formal is forced there, and ..N is recovered
expect_error(
  lapply(
    1,
    function(i, ints) dbarts::dbarts(x, y, interactions = ints),
    ints = interactions(max.order = 1)
  ),
  hint
)
expect_equal(
  maxOrder(lapply(
    1,
    function(i, ...) dbarts::dbarts(x, y, interactions = ..1),
    interactions(max.order = 1)
  )[[1L]]),
  1L
)

# (24) a caller formal the argument does not name is never forced
blocksOnly <- function(interactions, blocks) {
  dbarts::dbarts(x, y, blocks = blocks)
}
expect_identical(
  blockOf(blocksOnly(stop("forced"), dbartsForests$blocks(groups = groups))),
  blockOf(listed)
)

# --- the prior and family vocabularies give the same hint (dec-A75) ----------

passOnPrior <- function(prior) {
  dbarts::dbarts(x, y, tree.prior = prior)
}
expect_equal(
  passOnPrior(dbartsPriors$cgm(power = 3))$model@tree.prior@power,
  3
)
expect_error(
  passOnPrior(cgm(power = 3)),
  "outside the argument that takes it, write dbartsPriors\\$cgm"
)

passOnFamily <- function(family) {
  dbarts::dbarts(x, y, family = family)
}
expect_equal(
  attr(passOnFamily(dbartsFamilies$student(df = 6))$model, "resid.df"),
  6
)
expect_error(
  passOnFamily(student(df = 6)),
  "outside the argument that takes it, write dbartsFamilies\\$student"
)

# a name from a vocabulary this site never reads gets no hint: tree.prior
# only takes dbartsPriors, so a family constructor forced there (a wrapper's
# named formal, as above) is left as R's own plain message, not a hint
# naming dbartsFamilies - which tree.prior could never have meant
studentAtTreePrior <- tryCatch(
  passOnPrior(student(df = 6)),
  error = conditionMessage
)
expect_true(grepl(
  "could not find function \"student\"",
  studentAtTreePrior,
  fixed = TRUE
))
expect_false(grepl("dbartsFamilies", studentAtTreePrior, fixed = TRUE))
expect_false(grepl("outside the argument", studentAtTreePrior, fixed = TRUE))

# nor a forest constructor: tree.prior takes neither dbartsForests
interactionsAtTreePrior <- tryCatch(
  passOnPrior(interactions(max.order = 1)),
  error = conditionMessage
)
expect_true(grepl(
  "could not find function \"interactions\"",
  interactionsAtTreePrior,
  fixed = TRUE
))
expect_false(grepl("dbartsForests", interactionsAtTreePrior, fixed = TRUE))

# --- forwarded through dots ---------------------------------------------------

# (17) forwarded through one or more wrappers' dots
viaDots <- function(...) dbarts::dbarts(x, y, ...)
viaNested <- function(...) viaDots(...)
expect_equal(maxOrder(viaDots(interactions = interactions(max.order = 1))), 1L)
expect_equal(
  maxOrder(viaNested(interactions = interactions(max.order = 1))),
  1L
)

# (18) a ..N nested inside the argument
nestedDots <- function(...) {
  dbarts::dbarts(
    x,
    y,
    forests = list(forest(), forest(basis = ~z, interactions = ..1))
  )
}
expect_equal(
  forestsOf(nestedDots(interactions(max.order = 1)))$interactions[[
    2L
  ]]$max.order,
  1L
)

# a function literal's ..N are its own dots, not the caller's
lambdaDots <- function(...) {
  dbarts::dbarts(
    x,
    y,
    forests = c(
      list(forest()),
      Map(function(i, ...) forest(basis = ~z, interactions = ..1), 1, list(obj))
    )
  )
}
expect_equal(
  forestsOf(lambdaDots())$interactions[[2L]]$max.order,
  1L
)

# (19, 20) do.call names the caller unless given an environment that is no
# frame on the stack
expect_equal(
  maxOrder(do.call(
    viaDots,
    list(interactions = quote(interactions(max.order = 1)))
  )),
  1L
)
expect_error(
  do.call(
    viaDots,
    list(interactions = quote(interactions(max.order = 1))),
    envir = new.env()
  ),
  hint
)

# (25) a failing ..N read by two sites fails once, without R's warning on
# forcing the failed promise again
twoSites <- function(...) {
  dbarts::dbarts(x, y, interactions = ..1, blocks = ..1)
}
warned <- 0L
failure <- tryCatch(
  withCallingHandlers(
    twoSites(stop("forwarded")),
    warning = function(w) {
      warned <<- warned + 1L
      invokeRestart("muffleWarning")
    }
  ),
  error = conditionMessage
)
expect_identical(failure, "forwarded")
expect_identical(warned, 0L)
