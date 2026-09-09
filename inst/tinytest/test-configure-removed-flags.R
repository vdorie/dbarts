# Three 0.9-x configure options - --enable-match-bayes-tree,
# --enable-thread-safe-unload, --with-xint-size - are one-release stubs
# (dec-B107): each stops configure with a nonzero exit and a message naming
# the removal when passed explicitly, and configure with none of them still
# succeeds. configure.win never read these options, so this check is
# POSIX-configure only.
if (identical(.Platform$OS.type, "windows")) {
  exit_file(
    "configure is POSIX-only; Windows has no configure script (dec-B107, configure.win)"
  )
}

# tinytest runs against the INSTALLED package, which ships only inst/* and
# the compiled shared object - not configure, configure.ac, or the src/*.in
# templates AC_OUTPUT reads. This check can only run against a package
# source tree, which is not what dbarts installs; skip cleanly rather than
# guess at a source location (e.g. relative to getwd() or an env var), which
# would be fragile under a plain tinytest::test_package() run.
configurePath <- system.file("configure", package = "dbarts")
if (!nzchar(configurePath)) {
  exit_file(paste(
    "configure is not part of the installed package tree;",
    "this check needs a package source checkout",
    "(see docs/plans/interfaces-and-dependencies.md's Verification block",
    "for the equivalent manual shell check)"
  ))
}

# Runs configure out-of-tree (cwd is a scratch build dir; configure resolves
# its own srcdir from its real location, same as the manual verification:
# `cd $(mktemp -d) && /path/to/configure --enable-match-bayes-tree`).
runConfigure <- function(flags = character(0L)) {
  buildDir <- tempfile("dbarts-configure-stub-")
  dir.create(buildDir)
  on.exit(unlink(buildDir, recursive = TRUE), add = TRUE)
  owd <- setwd(buildDir)
  on.exit(setwd(owd), add = TRUE)
  output <- suppressWarnings(system2(
    configurePath,
    flags,
    stdout = TRUE,
    stderr = TRUE
  ))
  status <- attr(output, "status")
  list(status = if (is.null(status)) 0L else status, output = output)
}

removedFlags <- c(
  "--enable-match-bayes-tree" = "match-bayes-tree",
  "--enable-thread-safe-unload" = "thread-safe-unload",
  "--with-xint-size=32" = "xint-size"
)

for (flag in names(removedFlags)) {
  result <- runConfigure(flag)
  expect_true(result$status != 0L, info = flag)
  expect_true(
    any(grepl("removed for 1.0-0", result$output, fixed = TRUE)),
    info = flag
  )
  expect_true(
    any(grepl(removedFlags[[flag]], result$output, fixed = TRUE)),
    info = flag
  )
}

# Passing none of the three stubs still configures successfully.
noFlags <- runConfigure(character(0L))
expect_true(noFlags$status == 0L, info = paste(noFlags$output, collapse = "\n"))

rm(configurePath, runConfigure, removedFlags, flag, result, noFlags)
