# The tombstone registry is the one list of everything kept reachable for a
# release past its removal. It has to agree with what the package actually
# exports and with what the news file tells a user, or a stub outlives its
# announcement (or, worse, is announced and not there).

registry <- dbarts:::dbartsTombstones
expect_true(length(registry) > 0L)

field <- function(name) vapply(registry, `[[`, character(1L), name)

names.t <- field("name")
kinds <- field("kind")
expires <- field("expires")

# every entry carries the four fields the registry is read by
expect_true(all(vapply(
  registry,
  function(e) {
    all(c("name", "kind", "owner", "successor", "expires") %in% names(e))
  },
  logical(1L)
)))
expect_true(all(nzchar(names.t)))
expect_true(all(
  kinds %in%
    c("function", "method", "rcMethod", "argument", "family", "behaviour")
))

# one expiry for the whole set: the release that drops them deletes the file,
# so a per-entry date would be a lie the next reader has to check
expect_identical(unique(expires), dbarts:::tombstoneExpiry)
expect_true(
  package_version(dbarts:::tombstoneExpiry) >
    package_version(as.character(packageVersion("dbarts")))
)

# --- the registry agrees with NAMESPACE ---
namespaceFile <- file.path(
  dirname(system.file("DESCRIPTION", package = "dbarts")),
  "NAMESPACE"
)
expect_true(file.exists(namespaceFile))
namespaceText <- paste(readLines(namespaceFile), collapse = "\n")

exported <- getNamespaceExports("dbarts")
for (entry in registry[kinds == "function"]) {
  expect_true(entry$name %in% exported)
  expect_true(grepl(
    paste0("\\b", entry$name, "\\b"),
    namespaceText
  ))
}

# a method tombstone is registered on its generic for its own class, and the
# registration is what makes the stub reachable at all
for (entry in registry[kinds == "method"]) {
  generic <- sub("\\.[^.]+$", "", entry$name)
  expect_true(grepl(
    paste0("S3method(", generic, ", ", entry$owner, ")"),
    namespaceText,
    fixed = TRUE
  ))
  expect_inherits(getFromNamespace(entry$name, "dbarts"), "function")
}

# a reference-class tombstone is a live method on its generator
for (entry in registry[kinds == "rcMethod"]) {
  generator <- getFromNamespace("dbartsSampler", "dbarts")
  expect_true(entry$name %in% generator$methods())
}

# an argument tombstone is reachable on its owner: either as a formal still
# carried under the old name, or through that entry point's '...'
for (entry in registry[kinds == "argument"]) {
  owner <- getFromNamespace(entry$owner, "dbarts")
  ownFormals <- names(formals(owner))
  expect_true(entry$name %in% ownFormals || "..." %in% ownFormals)
  # the successor is either a formal that still stands on this entry point
  # (a rename) or the argument of an object the setting moved onto, whose
  # own formal is named first
  target <- sub(" =.*$", "", entry$successor)
  expect_true(entry$successor %in% ownFormals || target %in% ownFormals)
}

# --- the stubs themselves ---
expect_error(dbarts::rbart_vi(), pattern = "stan4bart")
expect_error(dbarts::rbart_vi(), pattern = "results move")
fakeRbart <- structure(list(), class = "rbart")
expect_error(predict(fakeRbart), pattern = "stan4bart")
expect_error(fitted(fakeRbart), pattern = "stan4bart")
expect_error(residuals(fakeRbart), pattern = "stan4bart")
expect_error(dbarts::extract(fakeRbart), pattern = "stan4bart")

# the two thread methods are no-ops, warned once, not errors: a 0.9-x Gibbs
# loop that brackets its sweeps with them still runs
threadSampler <- dbarts::dbarts(
  matrix(rnorm(40L), 20L, 2L),
  rnorm(20L),
  control = dbarts::dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 3L,
    n.samples = 2L,
    updateState = FALSE
  )
)
warnEnv <- dbarts:::onceWarnState
warnEnv[["tombstone.startThreads"]] <- NULL
warnEnv[["tombstone.stopThreads"]] <- NULL
expect_warning(threadSampler$startThreads(), pattern = "does nothing")
expect_warning(threadSampler$stopThreads(), pattern = "does nothing")
expect_null(suppressWarnings(threadSampler$startThreads()))
expect_null(suppressWarnings(threadSampler$stopThreads()))

# --- the consolidated argument names (dec-B98) ---

# each is accepted for one release, warned about once, and MAPPED onto the
# object it now rides, so the run is the one the old spelling asked for
consolidated <- dbarts:::onceWarnState
resetConsolidatedWarning <- function(name, caller) {
  consolidated[[paste0("tombstone.consolidated.", name, ".", caller)]] <- NULL
}

set.seed(717L)
nCons <- 40L
xCons <- matrix(rnorm(nCons * 2L), nCons, 2L)
yCons <- xCons[, 1L] + rnorm(nCons)
consControl <- dbarts::dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 5L,
  n.samples = 5L,
  updateState = FALSE
)

resetConsolidatedWarning("resid.dist", "dbarts")
expect_warning(
  samplerResidDist <- dbarts::dbarts(
    xCons,
    yCons,
    control = consControl,
    resid.dist = student(df = 5)
  ),
  pattern = "family = student"
)
expect_equal(attr(samplerResidDist$model, "resid.df"), 5)
# once per session: the second call is silent and still maps
expect_silent(
  samplerResidDistAgain <- dbarts::dbarts(
    xCons,
    yCons,
    control = consControl,
    resid.dist = student(df = 5)
  )
)
expect_equal(attr(samplerResidDistAgain$model, "resid.df"), 5)

resetConsolidatedWarning("dart", "bart")
expect_warning(
  fitDart <- dbarts::bart(
    xCons,
    yCons,
    dart = TRUE,
    n.trees = 5L,
    n.samples = 5L,
    n.burn = 2L,
    n.chains = 1L,
    n.threads = 1L,
    keepSampler = TRUE,
    verbose = FALSE
  ),
  pattern = "tree.prior = dart"
)
expect_inherits(fitDart$fit$model@tree.prior, "dbartsDartPrior")

resetConsolidatedWarning("levelGibbs", "bart")
expect_warning(
  fitLevelGibbs <- dbarts::bart(
    xCons,
    yCons,
    levelGibbs = TRUE,
    n.trees = 5L,
    n.samples = 5L,
    n.burn = 2L,
    n.chains = 1L,
    n.threads = 1L,
    keepSampler = TRUE,
    verbose = FALSE
  ),
  pattern = "tree.prior = cgm"
)
expect_true(fitLevelGibbs$fit$control@levelGibbs)

# an old spelling that names a family the call cannot fit is refused rather
# than resolved one way in silence
resetConsolidatedWarning("resid.dist", "bart")
expect_error(
  suppressWarnings(dbarts::bart(
    xCons,
    as.numeric(yCons > 0),
    family = "probit",
    resid.dist = student(df = 5),
    n.trees = 5L,
    n.samples = 5L,
    n.burn = 2L,
    n.chains = 1L,
    n.threads = 1L,
    verbose = FALSE
  )),
  pattern = "student residuals require a continuous gaussian response"
)

# --- xbart's own two (front-door S3) ---

# 'control' is refused by a message naming the flat arguments its settings
# became, not by R's own "unused argument": xbart builds its own control, so
# there is nothing to honour and the tombstone adds no capability
expect_error(
  dbarts::xbart(xCons, yCons, n.reps = 1L, control = consControl),
  pattern = "'control' has left 'xbart'"
)
expect_error(
  dbarts::xbart(xCons, yCons, n.reps = 1L, control = consControl),
  pattern = "n.cuts, useQuantiles, n.thin, storage"
)
# a three-element n.burn was 0.9-x's per-replication burn-in; chains are
# never carried between replications now, so the element is named and gone
expect_error(
  dbarts::xbart(xCons, yCons, n.reps = 1L, n.burn = c(2L, 1L, 1L)),
  pattern = "per-replication burn-in"
)

# and the consolidated names are gone from the signatures they left
for (name in c("resid.dist", "dispersion", "breaks", "max.rows")) {
  expect_false(name %in% names(formals(dbarts::bart)))
  expect_false(name %in% names(formals(dbarts::dbarts)))
}
expect_false("dart" %in% names(formals(dbarts::bart)))
expect_false("dart" %in% names(formals(dbarts::xbart)))
expect_false("levelGibbs" %in% names(formals(dbarts::bart)))

# --- the registry agrees with the news file ---
# The 1.0-0 news entry is the user-facing half of this list; the assertion
# below is the same one, read off inst/NEWS.Rd. Until the news entry lands
# (it is written with the manual, in the last slice of this arc), the file
# carries no tombstone list to check and this file stops here rather than
# pinning a section that does not exist yet.
newsFile <- file.path(
  dirname(system.file("DESCRIPTION", package = "dbarts")),
  "NEWS.Rd"
)
expect_true(file.exists(newsFile))
newsText <- paste(readLines(newsFile), collapse = "\n")
version1 <- regmatches(
  newsText,
  regexpr(
    "(?s)CHANGES IN VERSION 1\\.0-0\\}\\{.*?CHANGES IN VERSION 0",
    newsText,
    perl = TRUE
  )
)
expect_equal(length(version1), 1L)
if (!grepl("tombstone", version1, ignore.case = TRUE)) {
  exit_file(
    "NEWS tombstone list is written in the manual/news slice of this arc"
  )
}
for (nm in unique(names.t[kinds %in% c("function", "method", "argument")])) {
  expect_true(grepl(nm, version1, fixed = TRUE))
}
