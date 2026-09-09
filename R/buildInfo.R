# Reports which build of the shared object is loaded, for gates that are not
# valid on both. mode is "reference" when the package was configured with
# --enable-reference-build (on Windows, with the DBARTS_REFERENCE_BUILD
# environment variable set), which compiles the scalar, fixed-order draw-path
# kernels; it is "shipped" otherwise. Recorded bitwise baselines and
# seed-locked snapshot values are pinned to the reference build, so a test
# holding such a value must exit unless mode is "reference". compiled.isa
# names the instruction sets the binary carries and simd.level is the dispatch
# level chosen at load, the same integer
# C_dbarts_getMaxSIMDInstructionSet returns.
#
# Unexported on purpose: a development and test hook, not package API.
buildInfo <- function() {
  .Call(C_dbarts_buildInfo)
}
