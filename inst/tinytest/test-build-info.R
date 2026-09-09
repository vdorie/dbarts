# The build report and the gate split it drives. dbarts:::buildInfo()$mode is
# the only thing that tells R whether the loaded shared object has the scalar,
# fixed-order draw path the recorded bitwise values are pinned to, so this
# file checks both halves of that contract: the report's own shape, and that
# the seeded-drift snapshot files actually key off it.

info <- dbarts:::buildInfo()

expect_true(is.list(info))
expect_identical(names(info), c("mode", "compiled.isa", "simd.level"))
expect_true(is.character(info$mode) && length(info$mode) == 1L)
expect_true(info$mode %in% c("reference", "shipped"))
expect_true(is.character(info$compiled.isa))
expect_identical(
  info$simd.level,
  .Call(dbarts:::C_dbarts_getMaxSIMDInstructionSet)
)

# Every instruction set the binary reports must be one the dispatcher knows,
# and a nonzero dispatch level cannot be reached with nothing compiled in.
expect_true(all(
  info$compiled.isa %in% c("sse2", "sse4.1", "avx", "avx2", "neon")
))
expect_true(info$simd.level == 0L || length(info$compiled.isa) > 0L)

# The four seeded-drift files and their recorded assertion counts on the
# reference build. On the shipped build each exits at its guard, which leaves
# no assertion behind at all - so the count, not a skip attribute, is what
# distinguishes a guarded run from a green one, and a guard that failed to
# fire would show up here as a nonzero count under the shipped build.
snapshotFiles <- c(
  "test-reproducibility-continuousResponse-singleThreaded.R" = 14L,
  "test-reproducibility-continuousResponse-multithreaded.R" = 3L,
  "test-reproducibility-binaryResponse.R" = 7L,
  "test-reproducibility-xbart.R" = 3L
)
expected <- if (identical(info$mode, "reference")) {
  snapshotFiles
} else {
  setNames(rep(0L, length(snapshotFiles)), names(snapshotFiles))
}

observed <- vapply(
  names(snapshotFiles),
  function(file) {
    length(tinytest::run_test_file(
      system.file("tinytest", file, package = "dbarts"),
      verbose = 0L
    ))
  },
  integer(1L)
)
expect_identical(observed, expected)

rm(info, snapshotFiles, expected, observed)
