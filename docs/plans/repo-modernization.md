# repo-modernization

agent: sonnet
rng: neutral
budget: small, recurring

## Goal

The standing hygiene bundle, unchanged in intent from the pre-rewrite
TODO: CI pins current, sanitizer coverage widened, optional tooling
migrations decided when someone wants them.

## Context / Steps

1. CI action pins: recheck on each push; r-lib/actions stays @v2 (the
   latest major); bump when GitHub deprecates a runtime.
2. Sanitizers: add a gcc-asan/valgrind variant beside the existing
   workflow; gate the sanitizers workflow on pull_request once it has
   proven stable on pushes.
3. codecov: optional; wire only if a public badge is wanted (coverage
   measured locally on demand; 88 percent at last run).
4. roxygen2 / NEWS.md migrations: optional, taste-dependent, huge
   one-time diffs; do not start without VD asking (hand-written Rd is
   codoc-clean today; pkgdown renders NEWS.md natively).

## Constraints

- Each sub-item is independent; land separately.
- Out of scope: new CI jobs beyond the sanitizer variant
  (equivalence-ci owns the equivalence job).

## Verification

- Workflows green on push after each change; for sanitizers, a known
  clean run plus one deliberately injected error caught in a scratch
  branch (then reverted).

## Landing note (2026-07-07, 56e47f6): sub-items 1-2

Pin audit: all 8 workflows checked against current majors; only
sanitizers.yaml lagged (checkout@v5 -> v7); no deprecated runtimes.
Sanitizer widening: the clang job became a fail-fast-off matrix over
the r-hub clang-asan/gcc-asan containers; a new valgrind job runs the
suite under R -d valgrind on the r-hub valgrind container (180-minute
timeout, ERROR SUMMARY grep gate), same install/test/grep shape as
the proven asan job. Triggers unchanged (push + dispatch; the
pull_request gate stays deferred until stable). Gates: YAML parses;
runtime behavior maintainer-deferred to the first push. Sub-item 3
(codecov badge) awaits VD's call; sub-item 4 remains dormant by
design.
