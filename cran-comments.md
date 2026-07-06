## dbarts 1.0-0

This is a major release: the sampling engine is rewritten (C++20), with
new features documented in NEWS. Fits are statistically equivalent to
previous versions but not draw-for-draw identical, so seeded results
change.

## Breaking change for a compiled reverse dependency

dbarts 1.0-0 removes the C++ interface (dbarts/R_C_interface.hpp) whose
only CRAN consumer is stan4bart, and replaces it with a flat C API
(dbarts/dbarts.h). stan4bart 0.0-13, submitted alongside this release,
is ported to the new interface; the pair must be published together, as
the CRAN version of stan4bart (0.0-12) cannot build against dbarts
1.0-0. stan4bart is maintained by the same maintainer as dbarts, and
its 0.0-13 sources pass R CMD check against this dbarts release
(Status: OK).

## Reverse dependencies

Checked against this release on macOS (arm64), with
_R_CHECK_FORCE_SUGGESTS_=false:

- stan4bart 0.0-13 (the lockstep development version): OK
- OK: bartMan, EBcoBART, funcml, glossa, nlfh, riAFTBART, voi
- bartCause: 5 of 393 tests fail, all hardcoded expectations of seeded
  draws that this release documents as changing (same maintainer as
  dbarts; an update regenerating them will be submitted alongside)
- tidytreatment: one vignette WARNING from loading the CRAN stan4bart
  binary, i.e. the documented ABI break; resolves when stan4bart 0.0-13
  is published
- bartXViz: not checkable on the local machine (its own Fortran
  toolchain requirement); it uses dbarts from R code only
- Suggests-level dependencies (adrftools, bundle, butcher, countSTAR,
  insight, marginaleffects, MatchIt, mcmcsae, tidyAML, tmle,
  twoStageDesignTMLE, WeightIt) use dbarts conditionally and were not
  checked individually

## Test environments

- local: macOS 15 (arm64), R 4.6
- GitHub Actions: macOS (release), Windows (release), Ubuntu (devel,
  release, oldrel-1)
- clang ASAN/UBSAN (r-hub clang-asan container, R-devel): full test
  suite clean

## R CMD check results

0 errors, 0 warnings, 0 notes.
