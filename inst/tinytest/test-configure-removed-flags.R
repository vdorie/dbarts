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

# configure is not part of the installed package (only inst/* and the
# compiled shared object are); run_test_file/run_test_dir set the working
# directory to this file's own directory, so a source checkout's copy is two
# levels up from inst/tinytest.
configurePath <- system.file("configure", package = "dbarts")
if (!nzchar(configurePath)) {
  sourceConfigure <- file.path("..", "..", "configure")
  if (file.exists(sourceConfigure)) {
    configurePath <- normalizePath(sourceConfigure)
  }
}
if (!nzchar(configurePath)) {
  exit_file(
    "configure found in neither the installed package nor the source tree"
  )
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
