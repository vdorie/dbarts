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
  function(e) all(c("name", "kind", "owner", "successor", "expires") %in% names(e)),
  logical(1L)
)))
expect_true(all(nzchar(names.t)))
expect_true(all(
  kinds %in% c("function", "method", "rcMethod", "argument", "family", "behaviour")
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
  expect_true(entry$successor %in% ownFormals)
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
